"""
Utility functions for mapflpy.
"""
import math
import random

import numpy as np


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
    increment = math.pi*(3. - math.sqrt(5.));

    for i in range(samples):
        pid2 = .5*math.pi
        pi2 = 2*math.pi

        y = ((i*offset) - 1) + (offset/2);
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
