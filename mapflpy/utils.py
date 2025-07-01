"""
Utility functions for mapflpy.
"""
import importlib.util
import math
import random

import numpy as np


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


def save_trace_info(filename, traces, index={}, dtype=np.float32):
    """
    Save a list of mapflpy field line traces to binary file using numpy savez .npz format.

    This is required because sets of traces will generally not all have the same length.

    WARNING: This implementation uses the fact that we can save an numpy array of objects
    when `allow_pickle` is True. If we must save/load files w/out pickling for security
    reasons, then the format must be refactored (e.g. making a giant 3D array filled
    with NaNs and then compressing it or something and saving the index dict as a json with
    the arrays separated out.

    Parameters
    ----------
    filename : string
        path/name of the .npz file to be written (e.g. test.npz).
    traces : list
        A python list of field line traces. Each trace should be a (3,N) numpy array, where
        N is the number of points for each individual traces.
    index : dict, optional
        A dictionary that contains information about each trace. This can be whatever suites you
        for the project/application but should only contain basic python types or numpy arrays.
    dtype : numpy.dtype
        The data type to save the trace data (e.g. np.float32 or np.float64). Default is np.float32.
    """
    # turn each individual trace into a 1D array with the requested datatype
    traces_ravel = []
    for trace in traces:
        traces_ravel.append(trace.ravel().astype(dtype))

    # convert the trace list to an "array" of numpy arrays
    combo_array = np.array(traces_ravel, dtype=object)

    # save the file and the index
    np.savez(filename, combo_array=combo_array, index=index, allow_pickle=True)


def load_trace_info(filename):
    """
    Load a list of mapflpy field line traces that was saved by `save_trace_info`.

    Parameters
    ----------
    filename : string
        path/name of the numpy .npz file that was saved by `save_trace_info`.

    Returns
    -------
    traces : list
        A python list of field line traces. Each trace is assumed tobe a (3,N) numpy array, where
        N is the number of points for each individual traces.
    index : dict
        A dictionary that contains information about each trace.
    """
    # load the file
    data = np.load(filename, allow_pickle=True)

    # create the trace array by reshaping each of the raveled traces
    traces = []
    for trace1d in data['combo_array']:
        nseg = int(len(trace1d)/3)
        traces.append(np.reshape(trace1d, (3, nseg)))

    # get the index back as a regular dictionary
    index = dict(data['index'].item())

    return traces, index


def shift_phi_lps(lp, phi_shift=0.0):
    """
    Shift a Fortran ordered (3,N) launch point array in longitude by phi_shift radians.

    Parameters
    ----------
    lp : any
        Launch points for fieldline tracing. Here we assume `lp[2,:]` is the phi coordinate.
    phi_shift : float
        The longitudinal shift in radians. Defualt is 0.0.

    Returns
    -------
    A copy of lp shifted in longitude.

    """
    if np.isclose(phi_shift, 0.0):
        return lp
    else:
        lp[2, :] = np.mod(lp[2, :] + phi_shift, 2*np.pi)
        return lp


def shift_phi_traces(traces, phi_shift=0.0):
    """
    Shift mapflpy traces in longitude by phi_shift radians.

    Parameters
    ----------
    traces : list of numpy.ndarrrays.
        List of mapflpy traces. It assumes that for each `trace` in `traces` that `trace[2,:]` is the phi coordinate.
    phi_shift : float
        The longitudinal shift in radians. Defualt is 0.0.

    Returns
    -------
    A copy of the traces list shifted in longitude.

    """
    if np.isclose(phi_shift, 0.0):
        return traces
    else:
        for trace in traces:
            trace[2, :] = np.mod(trace[2, :] + phi_shift, 2*np.pi)
        return traces


def fibonacci_sphere(samples=100, randomize=False):
    """
    Generate a set of N points that will evenly sample a unit sphere. Unlike a uniform grid in phi/theta,
    These points are roughly equidistant at *all* lat lon locations (think a soccer ball).

    This code was adapted by RC from various samples on the internet and used in an MIDM poster.

    Parameters
    ----------
    samples : int
        Number of points to spread out over the unit sphere.
    randomize : bool
        Option to randomize where the N points show up (breaks up the eveness), default is False.

    Returns
    -------
    p: numpy.ndarray
        1D numpy array of phi (longitude) positions [radians, 0-2pi].
    t: numpy.ndarray
        1D numpy array of theta (co-latitude) positions [radians, 0-pi].
    """
    rnd = 1.
    if randomize:
        rnd = random.random()*samples

    points = []
    offset = 2./samples
    increment = math.pi*(3. - math.sqrt(5.))

    for i in range(samples):
        pid2 = .5*math.pi
        pi2 = 2*math.pi

        y = ((i*offset) - 1) + (offset/2)
        r = math.sqrt(1 - pow(y, 2))

        phi = ((i + rnd)%samples)*increment

        x = math.cos(phi)*r
        z = math.sin(phi)*r

        r = math.sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))

        if (r == 0):
            t = 0
        else:
            t = math.acos(z/r)
            # t=pid2 - t

        if (x == 0):
            if (y >= 0):
                p = pid2
            else:
                p = -pid2
        else:
            p = math.atan2(y, x)

        if (p < 0):
            p = p + pi2

        points.append([r, t, p])

    # make it a 2D array and return r, t, p separately
    points = np.array(points)
    r = points[:, 0]
    t = points[:, 1]
    p = points[:, 2]

    return p, t


def fetch_default_launch_points(n, r=1.01):
    """Generate a default set of N launch points for mapfl.

    The N launch points will roughly uniformly sample the sphere
    at a given radius using the `fibonacci_sphere` algorithm (default 1.01).

    Parameters
    ----------
    n: int
        Number of launch points to generate.
    r: float, optional
        Radius at which to place the launch points. Defaults to 1.01.

    Returns
    -------
    launch_points: ndarray
        A 3xN array of launch points.

    """
    p, t = fibonacci_sphere(n)
    return np.array((np.full_like(t, 1.01), t, p), order='F')
