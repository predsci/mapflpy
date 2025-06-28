"""
Test the functionality and accuracy of the Fortran tracing module.
"""
import copy
import os

import numpy as np

import psi_io

from mapflpy.utils import save_trace_info, load_trace_info, check_shared_object

from mapflpy.run import run_tracing_subprocess_caller


# ------------------------------------------------------------------------------
# Functions required by the test calculation
# ------------------------------------------------------------------------------
def dipole_field(lon=180., lat=45., resfac=1.0):
    """Generate a 3D dipole field in spherical coordinates.

    This function is for testing B tracing --> just assumes the dipole is at the origin.

    Parameters
    ----------
    lon: float
        longitude where the dipole moment vector is pointing (degrees). Default is 180.
    lat: float
        latitude where the dipole moment vector is pointing (degrees). Default is 45.
    resfac: float
        factor to multiply the resolution of the field by for testing purposes.

    Returns
    -------
    br: ndarray
        3D numpy array of the radial component of the field, Br.
    bt: ndarray
        3D numpy array of the theta component of the field, Bt.
    bp: ndarray
        3D numpy array of the phi component of the field, Bp.
    """
    # hardcoded global resolution
    n_p = int(np.round(181*resfac))
    n_t = int(np.round(91*resfac))
    n_r = int(np.round(71*resfac))
    r0 = 1.0
    r1 = 10.0

    # uniform grid in log(r), t, p
    logr = np.linspace(np.log10(r0), np.log10(r1), n_r)
    r = 10**logr
    t = np.linspace(0, np.pi, n_t)
    p = np.linspace(0, 2*np.pi, n_p)

    # dipole moment, normalize magnitude to max Br at 1Rs
    # dipole formula: mu0/4/pi*(3\vec{r}(\vec{m}\cdot\vec{r}/r^5 - \vec{m}/r^3)
    # assume it is at the sun center --> easy to get m0 based on B
    b0 = 2.0
    m0 = b0/2.0

    # get the 3D positions of each point
    p3d, t3d, r3d = np.meshgrid(p, t, r, indexing='ij')
    x3d, y3d, z3d = s2c(r3d, t3d, p3d)

    # get the dipole orientation
    mhat_t = np.deg2rad(90.0 - lat)
    mhat_p = np.deg2rad(lon)
    mhat_x = np.cos(mhat_p)*np.sin(mhat_t)
    mhat_y = np.sin(mhat_p)*np.sin(mhat_t)
    mhat_z = np.cos(mhat_t)
    mvec = m0*np.array([mhat_x, mhat_y, mhat_z])

    mdotr = mvec[0]*x3d + mvec[1]*y3d + mvec[2]*z3d

    bx = (3*x3d*mdotr/r3d**2 - mvec[0])/r3d**3
    by = (3*y3d*mdotr/r3d**2 - mvec[1])/r3d**3
    bz = (3*z3d*mdotr/r3d**2 - mvec[2])/r3d**3

    br, bt, bp = cvtosv(bx, by, bz, t3d, p3d)

    return br, bt, bp, r, t, p


def s2c(r, t, p):
    """
    convert numpy arrays of r,t,p (spherical) to x,y,z (cartesian)
    - r, t, p are numpy arrays of any shape (must match). t and p must be in radians.
    - x, y, z are returned in the same units as r.
    """
    ct = np.cos(t)
    st = np.sin(t)
    cp = np.cos(p)
    sp = np.sin(p)
    x = r*cp*st
    y = r*sp*st
    z = r*ct
    return x, y, z


def cvtosv(vx, vy, vz, t, p):
    """
    Convert a vector field in cartesian coordinates to spherical coordinates.
    Each input is a numpy array. t and p are the theta phi locations.
    """
    st = np.sin(t)
    ct = np.cos(t)
    sp = np.sin(p)
    cp = np.cos(p)

    # rotate the vector field components
    vr = vx*st*cp + vy*st*sp + vz*ct
    vt = vx*ct*cp + vy*ct*sp - vz*st
    vp = -vx*sp + vy*cp

    return vr, vt, vp


# -----------------------------------------------------------------------------
# Test Execution
# -----------------------------------------------------------------------------
def test_shared_object():
    check_shared_object()


def test_reference_dipole_trace(save_traces=False):
    """Check a new tracing against a reference tracing for a tilted dipole.

    Here we compute a simple tilted 3D dipole field on a simple grid.

    The B field is written to HDF5 files and used to trace field lines with mapflpy.

    The traces are compared to saved reference traces from a calculation where mapflpy
    was working as expected.

    If mapflpy has been updated, one can use this same function but set save_traces to True.
    In that case, a new reference datafile will be created.
    """

    # determine where the test is being called from
    work_dir = os.getcwd()

    # where the reference data is relative to
    ref_data_dir = os.path.join(os.path.dirname(__file__), 'data')

    print(f'### Directories:')
    print(f'  work_dir: {work_dir}')
    print(f'  ref_data_dir: {ref_data_dir}')

    # get the tilted dipole field (with resfac=0.5 its not a lot of data).
    print(f'### Generating B data')
    br, bt, bp, r, t, p = dipole_field(180, 45., resfac=0.5)

    # write it out to files.
    file_dict = {}
    print(f'### Writing B data')
    for vec, data in zip(['br', 'bt', 'bp'], [br, bt, bp]):
        file_dict[f'{vec}'] = os.path.join(work_dir, f'{vec}_dipole.h5')
        print(f'  {vec}: {file_dict[f'{vec}']}')
        psi_io.wrhdf_3d(file_dict[vec], r, t, p, data)

    # build some launch points
    t = np.deg2rad(90 - np.array([45., 35., 0., -45.]))
    p = np.deg2rad(np.array([90., 160, 180, 270.]))
    r = np.ones_like(t)*1.05
    launch_pts = np.stack([r, t, p])

    print(f'### Tracing Field Lines')
    fwd = True
    bwd = True
    traces_ = run_tracing_subprocess_caller(file_dict['br'], file_dict['bt'], file_dict['bp'], launch_pts, fwd, bwd)
    trace_list = list([arr[:, ~np.isnan(arr).any(axis=0)] for arr in traces_.traces.T])

    # load the reference tracing (or save a new one if save_traces is True)
    trace_reference_file = os.path.join(ref_data_dir, 'reference_tracing_dipole.npz')
    if save_traces:
        save_trace_info(trace_reference_file, trace_list, dtype=np.float64)
        traces_ref = copy.deepcopy(trace_list)
    else:
        traces_ref, dummy = load_trace_info(trace_reference_file)

    print(f'### Comparing current tracing to reference tracing:')
    print(f'  trace_file:\n    {trace_reference_file}')
    for i, trace in enumerate(trace_list):
        trace_ref = traces_ref[i]
        npts = len(trace[0, :])
        npts_ref = len(trace_ref[0, :])
        print(f'  comparing trace {i}: npts: {npts}, ', end='')

        # check that they are the same shape
        assert npts == npts_ref
        if npts != npts_ref:
            raise Exception(f'  trace {i} has a different length! Expected {npts_ref} but got {npts}!')

        # check the relative tolerance of the endpoints
        rtol = 1e-12
        for idir in [0, 1, 2]:
            isclose = np.isclose(trace[idir, -1], trace_ref[idir, -1], rtol=rtol)
            assert isclose == True
            if not isclose:
                raise Exception(f'\n  ERROR! trace {i} endpoint is not within {rtol:.1e} of reference endpoint!'
                                f'\n       trace[:,-1]: {trace[:, -1]},     trace.shape: {trace.shape}'
                                f'\n   trace_ref[:,-1]: {trace_ref[:, -1]}, trace_ref.shape: {trace.shape}'
                                f'\n   difference: {trace[:, -1] - trace_ref[:, -1]}')

        # summarize the difference
        print(f'maxdiff: {np.max(trace - trace_ref):.7e}, ')
