# ------------------------------------------------------------------------------
# Functions required by the test calculation
# ------------------------------------------------------------------------------
import gc
import importlib
from pathlib import Path
from typing import Tuple
import json

import numpy as np
from mapflpy.tracer import Tracer
from numpy._typing import NDArray

from mapflpy.utils import combine_fwd_bwd_traces, trim_fieldline_nan_buffer


def _recreate_lps_from_defaults(r_values: NDArray[float],
                                t_values: NDArray[float],
                                p_values: NDArray[float],
                                ) -> NDArray[float]:
    """Recreate launch points from default parameters."""
    mesh = np.meshgrid(r_values, t_values, p_values, indexing='ij')
    return np.stack([m.ravel() for m in mesh])


def dipole_field(lon: float = 180.,
                 lat: float = 45.,
                 resfac: float = 1.0,
                 r0: float = 1.0,
                 r1: float = 10.0,
                 r_resolution: int = 71,
                 t_resolution: int = 91,
                 p_resolution: int = 181
                 ) -> Tuple[NDArray[float], NDArray[float], NDArray[float], NDArray[float], NDArray[float], NDArray[float]]:
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
    r0: float
        inner radius of the grid (in Rs). Default is 1.
    r1: float
        outer radius of the grid (in Rs). Default is 10.
    r_resolution: int
        number of radial points. Default is 71.
    t_resolution: int
        number of theta points. Default is 91.
    p_resolution: int
        number of phi points. Default is 181.

    Returns
    -------
    br: ndarray
        3D numpy array of the radial component of the field, Br.
    bt: ndarray
        3D numpy array of the theta component of the field, Bt.
    bp: ndarray
        3D numpy array of the phi component of the field, Bp.
    r: ndarray
        1D numpy array of the radial scale.
    t: ndarray
        1D numpy array of the theta scale.
    p: ndarray
        1D numpy array of the phi scale.
    """
    # hardcoded global resolution
    n_p = int(np.round(p_resolution * resfac))
    n_t = int(np.round(t_resolution * resfac))
    n_r = int(np.round(r_resolution * resfac))

    # uniform grid in log(r), t, p
    logr = np.linspace(np.log10(r0), np.log10(r1), n_r)
    r = 10 ** logr
    t = np.linspace(0, np.pi, n_t)
    p = np.linspace(0, 2 * np.pi, n_p)

    # dipole moment, normalize magnitude to max Br at 1Rs
    # dipole formula: mu0/4/pi*(3\vec{r}(\vec{m}\cdot\vec{r}/r^5 - \vec{m}/r^3)
    # assume it is at the sun center --> easy to get m0 based on B
    b0 = 2.0
    m0 = b0 / 2.0

    # get the 3D positions of each point
    p3d, t3d, r3d = np.meshgrid(p, t, r, indexing='ij')
    x3d, y3d, z3d = s2c(r3d, t3d, p3d)

    # get the dipole orientation
    mhat_t = np.deg2rad(90.0 - lat)
    mhat_p = np.deg2rad(lon)
    mhat_x = np.cos(mhat_p) * np.sin(mhat_t)
    mhat_y = np.sin(mhat_p) * np.sin(mhat_t)
    mhat_z = np.cos(mhat_t)
    mvec = m0 * np.array([mhat_x, mhat_y, mhat_z])

    mdotr = mvec[0] * x3d + mvec[1] * y3d + mvec[2] * z3d

    bx = (3 * x3d * mdotr / r3d ** 2 - mvec[0]) / r3d ** 3
    by = (3 * y3d * mdotr / r3d ** 2 - mvec[1]) / r3d ** 3
    bz = (3 * z3d * mdotr / r3d ** 2 - mvec[2]) / r3d ** 3

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
    x = r * cp * st
    y = r * sp * st
    z = r * ct
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
    vr = vx * st * cp + vy * st * sp + vz * ct
    vt = vx * ct * cp + vy * ct * sp - vz * st
    vp = -vx * sp + vy * cp

    return vr, vt, vp


def check_shared_object(verbose=False):
    """Check that the mapflpy_fortran shared object can found on this installation.

    It seems that if the mapflpy_forran object can't be found, the subprocess pipes don't
    always indicate exactly what happened.

    This function will look for it and if it can't find it will raise an exception.

    This *shouldn't* actually import the shared object in case we are worried about
    thread safety (usually it gets imported by the tracing subprocess pipes, see `run.py`)
    """
    # first check if find_spec can run at all without crashing.
    try:
        module_spec = importlib.util.find_spec(".fortran.mapflpy_fortran", package='mapflpy')
    except Exception as e:
        raise Exception(f'{e}')

    # if it returns none, that means it cant find the module
    if module_spec is None:
        raise Exception(f'\n### Could not find the mapflpy shared object for Fortran tracing!')
    if verbose:
        print(f'### check_shared_object: mapflpy.fortran.mapflpy_fortran information:')
        print(module_spec)
        print('')
        # show which mapflpy modules have been loaded
        import sys
        for key in sys.modules.keys():
            if key.startswith('mapflpy'):
                print(sys.modules[key])

    # return the absolute path to the shared object file
    return module_spec.origin


def read_defaults() -> dict:
    """Read the project-wide defaults JSON file."""
    data_dir_path = Path(__file__).parent / 'data'
    defaults_filename = data_dir_path / 'defaults.json'
    with open(defaults_filename, "r") as f:
        data = json.load(f)
    return data


def generate_reference_traces() -> None:
    """
    Generate and persist a reference set of magnetic field-line traces.

    This routine iterates over all mesh presets and launch-point (LP) presets
    defined by :func:`read_defaults`. For each mesh, it synthesizes a dipole
    field, configures a :class:`Tracer` with the (Br, Bt, Bp) components and
    their spherical grids (r, theta, phi), then performs tracing for each LP
    preset in the requested direction(s) (forward, backward, or both).
    Any NaN padding/buffer added by the tracer is trimmed, and the results are
    written to a single NPZ archive at:
        ``<this_file_dir>/data/reference_traces.npz``

    Saved file layout
    -----------------
    The NPZ contains a flat mapping from string keys to arrays:

      - **key format**: ``"{mesh_key}_{lp_key}_{i}"``
          * ``mesh_key`` — name of the mesh preset (e.g. "coarse", "normal", "fine")
          * ``lp_key``   — one of {"fwd", "bwd", "both"}
          * ``i``        — zero-based index of the launch point within that preset
      - **value**: a NumPy array for a single trace produced by the tracer.

    Notes
    -----
    - This function has **side effects** only (writes a file) and returns ``None``.
    - A ``ValueError`` is raised if an unknown LP direction key is encountered.

    Returns
    -------
    None
    """
    metadata = {
        'description': 'Reference magnetic field-line traces for testing the Tracer class.',
        'tracefunc': 'mapflpy.tests.utils.generate_reference_traces',
        'meshfunc': 'mapflpy.tests.utils.dipole_field',
        'defaults': read_defaults(),
    }

    # Directory where the reference file will be saved (next to this module)
    data_dir_path = Path(__file__).parent / 'data'

    # Project-wide defaults that define mesh presets and launch-point (LP) presets
    defaults = read_defaults()

    # We'll accumulate all trace arrays here and dump them in one NPZ at the end.
    # Keys are "{mesh_key}_{lp_key}_{i}" → value is a single trace array.
    output_mapping = {}

    # Iterate through each mesh resolution/preset (e.g., "coarse", "normal", "fine")
    for mesh_key, mesh_value in defaults['mesh']['params'].items():
        # Synthesize a dipole field and its spherical grids for this mesh
        br, bt, bp, r, t, p = dipole_field(
            **defaults['mesh']['base'],
            **mesh_value,
        )

        # Create and configure the tracer with each field component and its grid.
        # The Tracer expects tuples of (values, r, theta, phi) for each component.
        tracer = Tracer()
        tracer.br = br, r, t, p
        tracer.bt = bt, r, t, p
        tracer.bp = bp, r, t, p

        # Iterate through LP presets and trace in the requested direction(s)
        for lp_key, lp_value in defaults['lps']['params'].items():
            match lp_key:
                case 'fwd':
                    # Recreate launch points for this preset and trace forward
                    lps = _recreate_lps_from_defaults(**defaults['lps']['base'], **lp_value)
                    tracer.set_tracing_direction('f')
                    traces = tracer.trace(lps, buffer_size=defaults['lps']['BUFFER'])

                case 'bwd':
                    # Recreate launch points and trace backward
                    lps = _recreate_lps_from_defaults(**defaults['lps']['base'], **lp_value)
                    tracer.set_tracing_direction('b')
                    traces = tracer.trace(lps, buffer_size=defaults['lps']['BUFFER'])

                case 'both':
                    # Trace both forward and backward, then combine into a single trace
                    lps = _recreate_lps_from_defaults(**defaults['lps']['base'], **lp_value)
                    tracer.set_tracing_direction('f')
                    fwd_traces = tracer.trace(lps, buffer_size=defaults['lps']['BUFFER'])
                    tracer.set_tracing_direction('b')
                    bwd_traces = tracer.trace(lps, buffer_size=defaults['lps']['BUFFER'])
                    traces = combine_fwd_bwd_traces(fwd_traces, bwd_traces)

                case _:
                    # Guard against configuration typos / unsupported directions
                    raise ValueError(f'Unknown tracing direction: {lp_key}')

            # Trim any NaN padding the tracer may add to pre-allocated buffers
            traces_list = trim_fieldline_nan_buffer(traces)

            # Store each individual trace in the output mapping with a stable key
            for i, arr in enumerate(traces_list):
                output_mapping[f'{mesh_key}_{lp_key}_{i}'] = arr

        # Explicit cleanup between mesh presets to control memory
        del tracer
        gc.collect()

    # Persist all traces into a single NPZ. The consumer can load by key.
    # NOTE: np.savez does not take an `allow_pickle` kwarg; that's used by np.load.
    np.savez(data_dir_path / "reference_traces", __meta__=np.array(json.dumps(metadata)), **output_mapping)


def compute_fieldline_length(trace: NDArray[float]) -> float:
    """Compute the length of a field line trace.

    Parameters
    ----------
    trace: ndarray
        2D numpy array of shape (3, N) representing the field line trace in spherical coordinates.

    Returns
    -------
    length: float
        Length of the field line trace.
    """
    # Convert spherical to Cartesian coordinates for distance calculation
    cartesian_spline = np.asarray(s2c(*trace))
    spline_difference = np.diff(cartesian_spline, axis=1)  # (3, N-1): consecutive segment vectors
    seg_lengths = np.linalg.norm(spline_difference, axis=0)  # (N-1,)
    length = seg_lengths.sum()
    return length


def compute_weighted_fieldline_difference(trace_test: NDArray[float],
                                          trace_ref: NDArray[float]) -> NDArray[float]:
    """Compute the weighted difference between two field line traces.

    Parameters
    ----------
    trace_test: ndarray
        2D numpy array of shape (3, N) representing the first field line trace in spherical coordinates.
    trace_ref: ndarray
        2D numpy array of shape (3, N) representing the second field line trace in spherical coordinates.

    Returns
    -------
    weighted_difference: ndarray
        Weighted difference between the two field line traces.
    """
    # Convert spherical to Cartesian coordinates for distance calculation
    cartesian_trace_test = np.asarray(s2c(*trace_test))
    cartesian_trace_ref = np.asarray(s2c(*trace_ref))

    distances = np.linalg.norm(cartesian_trace_test - cartesian_trace_ref, axis=0)  # (N,)
    return distances / trace_ref[0, ...]


if __name__ == "__main__":
    # If this script is run directly, generate the reference traces
    print("Generating reference traces...")
    generate_reference_traces()
    print("Reference traces generated successfully.")
