"""
Methods to instantiate and call the mapflpy_fortran tracer using subprocess commands.
"""
import copy
import warnings
import multiprocessing as mp

import numpy as np

from .tracer import Tracer, fetch_default_launch_points, Traces
from .utils import shift_phi_traces, shift_phi_lps, check_shared_object

def run_tracing_subprocess(pipe, br, bt, bp, lp, fwd, bwd, **mapfl_params):
    """
    Run tracing in a subprocess and send results via a pipe.

    This function creates a Tracer instance from the `tracer` module using the provided
    parameters and performs forward and/or backward tracing based on the flags. The resulting
    traces are processed and sent back to the parent process through the given pipe.

    Parameters
    ----------
    pipe : multiprocessing.connection.Connection
        The pipe used to send the tracing results back to the parent process.
    br : str
        Path to hdf4 or hdf5 Br file.
    bt : str
        Path to hdf4 or hdf5 Bt file.
    bp : str
        Path to hdf4 or hdf5 Bp file.
    lp : any
        Launch points for fieldline tracing.
    fwd : bool
        Flag indicating whether to perform forward tracing.
    bwd : bool
        Flag indicating whether to perform backward tracing.
    **mapfl_params
        Additional keyword arguments to be passed to the Tracer initialization.

    Returns
    -------
    None
        The function does not return a value directly. It sends the tracing results (or `None` if an
        error occurs) through the provided pipe before closing it.
    """
    try:
        import mapflpy.fortran.mapflpy_fortran as mapflpy_fortran

        tracer = Tracer(mapflpy_fortran, br, bt, bp, lp, **mapfl_params)
        tracer.load_fields()
        if fwd and bwd:
            tracer.prep_tracing('f')
            fwd_traces = tracer.run_tracing()
            tracer.prep_tracing('b')
            bwd_traces = tracer.run_tracing()
            cmb_traces = Traces(np.concatenate([np.flip(bwd_traces.traces, axis=0)[:-1,:,:], fwd_traces.traces]), fwd_traces.end_pos, fwd_traces.traced_to_boundary & bwd_traces.traced_to_boundary)
            pipe.send(cmb_traces)
        elif fwd:
            tracer.prep_tracing('f')
            fwd_traces = tracer.run_tracing()
            pipe.send(fwd_traces)
        elif bwd:
            tracer.prep_tracing('b')
            bwd_traces = tracer.run_tracing()
            pipe.send(bwd_traces)
        else:
            warnings.warn("You must provide at least one of fwd or bwd to run the tracer.")
            pipe.send(Traces(np.array([]), np.array([]), np.array([])))
    except Exception as e:
        warnings.warn(str(e))
        pipe.send(None)
    finally:
        pipe.close()


def run_tracing_subprocess_caller(br, bt, bp, lp=None, fwd=True, bwd=True, **mapfl_params):
    """
    Run a tracing subprocess and retrieve the results.

    This function creates a new process to execute the `run_tracing_subprocess` function.
    It sets up a multiprocessing pipe for inter-process communication and returns the tracing
    results received from the subprocess.

    Parameters
    ----------
    br : str
        Path to hdf4 or hdf5 Br file.
    bt : str
        Path to hdf4 or hdf5 Bt file.
    bp : str
        Path to hdf4 or hdf5 Bp file.
    lp : any
        Launch points for fieldline tracing.
    fwd : bool, optional
        Flag indicating whether to perform forward tracing; default is True.
    bwd : bool, optional
        Flag indicating whether to perform backward tracing; default is True.
    **mapfl_params
        Additional keyword arguments to be passed to the `run_tracing_subprocess` function.

    Returns
    -------
    result : any
        The tracing results received from the subprocess via the pipe. If an error occurs or no
        result is received, `None` is returned.

    Notes
    -----
    The function uses the Python multiprocessing module to handle process creation and pipe communication.
    """
    # check that the mapflpy_fortran module exists
    check_shared_object()

    # create the pipe
    ctx = mp.get_context('fork')
    reciever, sender = ctx.Pipe()

    # call the tracing within the subprocess
    p = ctx.Process(target=run_tracing_subprocess, args=(sender, br, bt, bp, lp, fwd, bwd), kwargs=mapfl_params)
    p.start()

    # check for an exception
    try:
        result = reciever.recv()
    except Exception as e:
        print("No result received from subprocess.")
        print(f"EXCEPTION: {str(e)}")
        result = None

    # error out if the result from either case (success/failure) is still None
    if result is None:
        # check for an import error
        try:
            import mapflpy.fortran.mapflpy_fortran as mapflpy_fortran
            print(f'asdfasdfasdf')
        except Exception as e:
            raise Exception(f'\n### mapflpy.fortran.mapflpy_fortran object could not be imported!'
                            f'\n  Exception: {e}')
        # otherwise complain about the value
        raise ValueError("\n### The mapflpy trace returned no results. This shouldn't happen!"
                         f"\n  result: {result}")

    return result


def tracer_listener(pipe, br, bt, bp, lp=None, **mapfl_params):
    """
    Listen for tracer commands on a pipe and execute them.

    This function creates a Tracer instance from the `tracer` module and continuously listens for
    incoming commands on the provided pipe. Upon receiving a command, it dynamically calls the corresponding
    method on the tracer instance and sends back the result.

    Parameters
    ----------
    pipe : multiprocessing.connection.Connection
        The pipe used for receiving commands and sending back the execution results.
    br : str
        Path to hdf4 or hdf5 Br file.
    bt : str
        Path to hdf4 or hdf5 Bt file.
    bp : str
        Path to hdf4 or hdf5 Bp file.
    lp : any
        Launch points for fieldline tracing.
    **mapfl_params
        Additional keyword arguments to be passed to the Tracer initialization.

    Returns
    -------
    None
        The function runs indefinitely until a "STOP" command is received. It does not return any value.
    """
    import mapflpy.fortran.mapflpy_fortran as mapflpy_fortran
    tracer = Tracer(mapflpy_fortran, br, bt, bp, lp, **mapfl_params)
    while True:
        try:
            method, *args = pipe.recv()
            if method == "STOP":
                break
            result = getattr(tracer, method)(*args)
            if result is not None:
                pipe.send(result)
        except Exception as e:
            pipe.send(f"Error: {e}")
    pipe.close()


def inter_domain_tracing_caller(br_cor='', bt_cor='', bp_cor='',
                                br_hel='', bt_hel='', bp_hel='',
                                lp=None, maxiter=10,
                                r_interface=30.0,
                                helio_shift=0.0,
                                rtol=1e-5,
                                **mapfl_params):
    """
    Perform inter-domain tracing using two tracer processes.

    This method sets up two tracer processes (e.g., for different magnetic domains) that run concurrently.
    It coordinates the tracing between these two processes via multiprocessing pipes. Because launch points
    that start in the corona or heliosphere are handled differently, this function wraps the lower-level inter-domain
    tracing methods to trace forward and backwards from launch points in any domain, joins them together and returns
    all traces.

    Parameters
    ----------
    br_cor : str
        Path to hdf4 or hdf5 Br file (coronal domain).
    bt_cor : str
        Path to hdf4 or hdf5 Bt file (coronal domain).
    bp_cor : str
        Path to hdf4 or hdf5 Bp file (coronal domain).
    br_hel : str
        Path to hdf4 or hdf5 Br file (heliospheric domain).
    bt_hel : str
        Path to hdf4 or hdf5 Bt file (heliospheric domain).
    bp_hel : str
        Path to hdf4 or hdf5 Bp file (heliospheric domain).
    lp : any, optional
        Launch points used by the tracer. Default is None.
    maxiter : int, optional
        Maximum number of iterations for handling boundary recrossings. Default is 10.
    r_interface : float, optional
        Radius at which to connect the traces between domains. Default is 30.
    helio_shift : float, optional
        Longitudinal shift angle between the heliospheric domain and the coronal domain in RADIANS.
        This shift is ADDED to the coronal launch point phi positions. Default is 0.0.
    rtol : float, optional
        Relative tolerance for `np.isclose` for checking a trace has hit the interface boundary. Default is 1e-5.
    **mapfl_params
        Additional keyword arguments to be passed to both tracer initializations.

    Returns
    -------
    final_traces : list
        A list of numpy arrays representing the concatenated tracing results for each launch point.
    traced_to_boundary : numpy.ndarray
        A boolean array indicating whether this trace hit the inner cor or outer hel boundary on both ends.
    boundary_recross : numpy.ndarray
        A boolean array indicating whether a boundary recrossing occurred for each launch point.

    Notes
    -----
    The function uses two separate processes to avoid sharing `mapflpy_fortran` objects between domains.
    """
    # check that the mapflpy_fortran module exists
    check_shared_object()

    # setup the pipes
    cor_reciever, cor_sender = mp.Pipe()
    hel_reciever, hel_sender = mp.Pipe()

    # make sure the interface radius is set properly for both corona and helio
    mapfl_params_cor = copy.deepcopy(mapfl_params)
    mapfl_params_hel = copy.deepcopy(mapfl_params)
    mapfl_params_cor['domain_r_max_'] = r_interface
    mapfl_params_hel['domain_r_min_'] = r_interface

    # initialize the pipes
    cor_p = mp.Process(target=tracer_listener, args=(cor_sender, br_cor, bt_cor, bp_cor), kwargs=mapfl_params_cor)
    hel_p = mp.Process(target=tracer_listener, args=(hel_sender, br_hel, bt_hel, bp_hel), kwargs=mapfl_params_hel)
    cor_p.start()
    hel_p.start()

    # double check that they are indeed separated
    cor_reciever.send(('get_mapflpy_id',))
    hel_reciever.send(('get_mapflpy_id',))
    assert cor_reciever.recv() != hel_reciever.recv(), 'mapfl object is shared'

    # load default launch points if necessary
    if lp is None:
        lp = fetch_default_launch_points(128)

    # prepare the final arrays
    n_lp = lp.shape[1]
    final_traces = [None]*n_lp
    boundary_recross = np.array([False]*n_lp)
    traced_to_boundary = np.array([False]*n_lp)

    # determine which indexes of the launch points are coronal and which are heliospheric
    inds_coronal = np.where(lp[0, :] <= r_interface)[0]
    inds_helio = np.where(lp[0, :] > r_interface)[0]

    # separate the launch points by domain (length 0 arrays are OK here)
    cor_lp = lp[:, inds_coronal]
    hel_lp = lp[:, inds_helio]

    # CORONAL INTERDOMAIN TRACE
    if len(inds_coronal) > 0:
        # trace coronal launch points forward from the coronal domain
        traces_cor_fwd, bndry_cor_fwd, recross_cor_fwd = inter_domain_tracing_from_cor(cor_reciever, hel_reciever,
                                                                                       direction='f',
                                                                                       lp=cor_lp, maxiter=maxiter,
                                                                                       r_interface=r_interface,
                                                                                       helio_shift=helio_shift)

        # trace coronal launch points backwards from the coronal domain
        traces_cor_bwd, bndry_cor_bwd, recross_cor_bwd = inter_domain_tracing_from_cor(cor_reciever, hel_reciever,
                                                                                       direction='b',
                                                                                       lp=cor_lp, maxiter=maxiter,
                                                                                       r_interface=r_interface,
                                                                                       helio_shift=helio_shift)

        # join the traces, flipping the backwards trace after dropping its first point (first point is the starting point)
        traces_cor = [np.concatenate((np.flip(traces_cor_bwd[i][:, 1:], axis=1), traces_cor_fwd[i]), axis=1)
                      for i in range(len(traces_cor_bwd))]

        # combine the tracing flags
        recross_cor = np.logical_and(recross_cor_bwd, recross_cor_fwd)
        bndry_cor = np.logical_and(bndry_cor_bwd, bndry_cor_fwd)

        # populate the coronal launch points
        for i, trace in zip(inds_coronal, traces_cor):
            final_traces[i] = trace
        boundary_recross[inds_coronal] = recross_cor
        traced_to_boundary[inds_coronal] = bndry_cor

    # HELIOSPHERIC INTERDOMAIN TRACE
    if len(inds_helio) > 0:
        # trace heliospheric launch points forward from the heliospheric domain
        traces_hel_fwd, bndry_hel_fwd, recross_hel_fwd = inter_domain_tracing_from_hel(cor_reciever, hel_reciever,
                                                                                       direction='f',
                                                                                       lp=hel_lp, maxiter=maxiter,
                                                                                       r_interface=r_interface,
                                                                                       helio_shift=helio_shift,
                                                                                       rtol=rtol)

        # trace heliospheric launch points backwards from the heliospheric domain
        traces_hel_bwd, bndry_hel_bwd, recross_hel_bwd = inter_domain_tracing_from_hel(cor_reciever, hel_reciever,
                                                                                       direction='b',
                                                                                       lp=hel_lp, maxiter=maxiter,
                                                                                       r_interface=r_interface,
                                                                                       helio_shift=helio_shift,
                                                                                       rtol=rtol)

        # join the traces, flipping the backwards trace after dropping its first point (first point is the starting point)
        traces_hel = [np.concatenate((np.flip(traces_hel_bwd[i][:, 1:], axis=1), traces_hel_fwd[i]), axis=1)
                      for i in range(len(traces_hel_bwd))]

        # combine the tracing flags
        recross_hel = np.logical_and(recross_hel_bwd, recross_hel_fwd)
        bndry_hel = np.logical_and(bndry_hel_bwd, bndry_hel_fwd)

        # populate the heliospheric launch points
        for i, trace in zip(inds_helio, traces_hel):
            final_traces[i] = trace
        boundary_recross[inds_helio] = recross_hel
        traced_to_boundary[inds_helio] = bndry_hel

    # kill the mapflpy_fortran subprocess pipes.
    cor_reciever.send(('STOP',))
    hel_reciever.send(('STOP',))
    cor_reciever.close()
    hel_reciever.close()

    return final_traces, traced_to_boundary, boundary_recross


def inter_domain_tracing_from_cor(cor_reciever, hel_reciever, direction='f', lp=None,
                                  maxiter=10, r_interface=30.0, helio_shift=0.0, rtol=1e-5):
    """
    Perform inter-domain coronal and heliospheric tracing for CORONAL launch points in the specified direction.

    This method receives two tracer processes (one for each domain) that are run concurrently. The function initiates
    tracing in one process, checks for boundary recrossings, and if necessary, alternates the tracing between
    the two processes until the tracing endpoints no longer cross a defined boundary or the maximum number of
    iterations is reached.

    Parameters
    ----------
    cor_reciever : multiprocessing.connection.Connection
        The coronal domain pipe that does the mapfl tracing (see `tracer_listener`).
    hel_reciever : multiprocessing.connection.Connection
        The heliospheric domain pipe that does the mapfl tracing (see `tracer_listener`).
    direction : str
        The direction of the mapfl tracings. This must be either 'f' or 'b' (forwards or backwards). Default is 'f'.
    lp : any
        Launch points for fieldline tracing.
    maxiter : int, optional
        Maximum number of iterations for handling boundary recrossings. Default is 10.
    r_interface : float, optional
        Radius at which to connect the traces between domains. Default is 30.
    helio_shift : float, optional
        Longitudinal shift angle between the heliospheric domain and the coronal domain in RADIANS.
        This shift is ADDED to the coronal launch point phi positions. Default is 0.0.
    rtol : float, optional
        Relative tolerance for `np.isclose` for checking a trace has hit the interface boundary. Default is 1e-5.

    Returns
    -------
    final_traces : list
        A list of numpy arrays representing the concatenated tracing results for each launch point.
    traced_to_boundary : numpy.ndarray
        A boolean array indicating whether this trace hit the inner cor or outer hel boundary.
    boundary_recross : numpy.ndarray
        A boolean array indicating whether a boundary recrossing occurred for each launch point.

    """
    # set the launch points, make a copy so input lp doesn't get vaporized on succesive traces
    cor_lps = copy.deepcopy(lp)
    cor_reciever.send(('set_launch_points', cor_lps))

    # trace in the corona going the specified direction
    cor_reciever.send(('prep_tracing', direction,))
    cor_reciever.send(('run_tracing',))

    # get the raw coronal traces back and convert to the trimmed lists
    traces_ = cor_reciever.recv()
    final_traces = list([arr[:, ~np.isnan(arr).any(axis=0)] for arr in traces_.traces.T])

    # get the end positions of the traces in r,t,p, shape: (3,n)
    radial_end_pos = np.copy(traces_.end_pos)

    # determine which traces hit the interface boundary
    midboundary_mask = np.isclose(radial_end_pos[0, :], r_interface, rtol=rtol)

    # initialize the array for checking that you went back through
    boundary_recross = np.full_like(midboundary_mask, False)

    # check that you hit a boundary to end the trace
    traced_to_boundary = traces_.traced_to_boundary

    while np.any(midboundary_mask) and maxiter > 0:
        # set the new heliospheric launchpoints using any that hit the interface from the corona
        # these must also be shifted FORWARD in phi by the helio shift value
        hel_lps = shift_phi_lps(radial_end_pos[:, midboundary_mask], helio_shift)
        hel_reciever.send(('set_launch_points', hel_lps))

        # trace forwards in the helio run
        hel_reciever.send(('prep_tracing', 'f',))
        hel_reciever.send(('run_tracing',))

        # return and digest the helio traces
        traces_ = hel_reciever.recv()
        temp_traces = list([arr[:, ~np.isnan(arr).any(axis=0)] for arr in traces_.traces.T])

        # shift these traces BACK to the coronal/carrington frame
        temp_traces = shift_phi_traces(temp_traces, -helio_shift)

        # add this trace segment neglecting the first point since it duplicates the last point of the previous segment.
        for i, trace in zip(np.where(midboundary_mask)[0], temp_traces):
            final_traces[i] = np.concatenate([final_traces[i], trace[:, 1:]], axis=1)

        # check that you hit a boundary
        traced_to_boundary[midboundary_mask] = traces_.traced_to_boundary

        # update the radial end positions (SHIFTED BACK!)
        radial_end_pos[:, midboundary_mask] = shift_phi_lps(np.copy(traces_.end_pos), -helio_shift)

        # update the flag for traces that hit the interface
        midboundary_mask = np.isclose(radial_end_pos[0, :], r_interface, rtol=rtol)

        # now trace through the corona BACKWARDS
        if np.any(midboundary_mask):
            boundary_recross |= midboundary_mask

            # take the subset of launch points and trace.
            cor_lps = radial_end_pos[:, midboundary_mask]
            cor_reciever.send(('set_launch_points', cor_lps))
            cor_reciever.send(('prep_tracing', 'b',))
            cor_reciever.send(('run_tracing',))

            # digest the traces and add these segments, again removing the first point of the new segment.
            traces_ = cor_reciever.recv()
            temp_traces = list([arr[:, ~np.isnan(arr).any(axis=0)] for arr in traces_.traces.T])
            for i, trace in zip(np.where(midboundary_mask)[0], temp_traces):
                final_traces[i] = np.concatenate([final_traces[i], trace[:, 1:]], axis=1)

            # check the trace, update the end positions and the midboundary flag, continue the loop
            traced_to_boundary[midboundary_mask] = traces_.traced_to_boundary
            radial_end_pos[:, midboundary_mask] = traces_.end_pos
            midboundary_mask = np.isclose(radial_end_pos[0, :], r_interface, rtol=rtol)
            boundary_recross |= midboundary_mask

        # if no more work to be done, break the loop
        else:
            break
        maxiter -= 1

    # return the final traces and tracing checks
    return final_traces, traced_to_boundary, boundary_recross


def inter_domain_tracing_from_hel(cor_reciever, hel_reciever, direction='f', lp=None,
                                  maxiter=10, r_interface=30.0, helio_shift=0.0, rtol=1e-5):
    """
    Perform inter-domain coronal and heliospheric tracing for HELIOSPHERIC launch points in the specified direction.

    This method receives two tracer processes (one for each domain) that are run concurrently. The function initiates
    tracing in one process, checks for boundary recrossings, and if necessary, alternates the tracing between
    the two processes until the tracing endpoints no longer cross a defined boundary or the maximum number of
    iterations is reached.

    Parameters
    ----------
    cor_reciever : multiprocessing.connection.Connection
        The coronal domain pipe that does the mapfl tracing (see `tracer_listener`).
    hel_reciever : multiprocessing.connection.Connection
        The heliospheric domain pipe that does the mapfl tracing (see `tracer_listener`).
    direction : str
        The direction of the mapfl tracings. This must be either 'f' or 'b' (forwards or backwards). Default is 'f'.
    lp : any
        Launch points for fieldline tracing.
    maxiter : int, optional
        Maximum number of iterations for handling boundary recrossings. Default is 10.
    r_interface : float, optional
        Radius at which to connect the traces between domains. Default is 30.
    helio_shift : float, optional
        Longitudinal shift angle between the heliospheric domain and the coronal domain in RADIANS.
        This shift is ADDED to the coronal launch point phi positions. Default is 0.0.
    rtol : float, optional
        Relative tolerance for `np.isclose` for checking a trace has hit the interface boundary. Default is 1e-5.

    Returns
    -------
    final_traces : list
        A list of numpy arrays representing the concatenated tracing results for each launch point.
    traced_to_boundary : numpy.ndarray
        A boolean array indicating whether this trace hit the inner cor or outer hel boundary.
    boundary_recross : numpy.ndarray
        A boolean array indicating whether a boundary recrossing occurred for each launch point.

    """
    # copy the launch points and shift thim accordingly
    hel_lps = copy.deepcopy(lp)
    hel_lps = shift_phi_lps(hel_lps, helio_shift)

    # set the launch points
    hel_reciever.send(('set_launch_points', hel_lps))

    # trace in the heliosphere going the specified direction
    hel_reciever.send(('prep_tracing', direction,))
    hel_reciever.send(('run_tracing',))

    # get the raw helio traces back and convert to the trimmed lists
    traces_ = hel_reciever.recv()
    final_traces = list([arr[:, ~np.isnan(arr).any(axis=0)] for arr in traces_.traces.T])

    # shift these traces BACK to the coronal/carrington frame
    final_traces = shift_phi_traces(final_traces, -helio_shift)

    # get and shift the end positions of the traces in r,t,p, shape: (3,n)
    radial_end_pos = shift_phi_lps(np.copy(traces_.end_pos), -helio_shift)

    # determine which traces hit the interface boundary
    midboundary_mask = np.isclose(radial_end_pos[0, :], r_interface, rtol=rtol)

    # initialize the array for checking that you went back through
    boundary_recross = np.full_like(midboundary_mask, False)

    # check that you hit a boundary to end the trace
    traced_to_boundary = traces_.traced_to_boundary

    while np.any(midboundary_mask) and maxiter > 0:
        # set the new coronal launchpoints using any that hit the interface from the heliosphere
        cor_lps = radial_end_pos[:, midboundary_mask]
        cor_reciever.send(('set_launch_points', cor_lps))

        # trace backwards in the coronal run
        cor_reciever.send(('prep_tracing', 'b',))
        cor_reciever.send(('run_tracing',))

        # return and digest the coronal traces
        traces_ = cor_reciever.recv()
        temp_traces = list([arr[:, ~np.isnan(arr).any(axis=0)] for arr in traces_.traces.T])

        # add this trace segment neglecting the first point since it duplicates the last point of the previous segment.
        for i, trace in zip(np.where(midboundary_mask)[0], temp_traces):
            final_traces[i] = np.concatenate([final_traces[i], trace[:, 1:]], axis=1)

        # check that you hit a boundary
        traced_to_boundary[midboundary_mask] = traces_.traced_to_boundary

        # update the radial end positions
        radial_end_pos[:, midboundary_mask] = np.copy(traces_.end_pos)

        # update the flag for traces that hit the interface
        midboundary_mask = np.isclose(radial_end_pos[0, :], r_interface, rtol=rtol)

        # now trace through the heliosphere FORWARDS
        if np.any(midboundary_mask):
            boundary_recross |= midboundary_mask

            # set the new heliospheric launchpoints using any that hit the interface from the corona
            # these must also be shifted FORWARD in phi by the helio shift value
            hel_lps = shift_phi_lps(radial_end_pos[:, midboundary_mask], helio_shift)

            # take the subset of launch points and trace.
            hel_reciever.send(('set_launch_points', hel_lps))
            hel_reciever.send(('prep_tracing', 'f',))
            hel_reciever.send(('run_tracing',))

            # digest the traces and add these segments, again removing the first point of the new segment.
            traces_ = hel_reciever.recv()
            temp_traces = list([arr[:, ~np.isnan(arr).any(axis=0)] for arr in traces_.traces.T])

            # shift these traces BACK to the coronal/carrington frame
            temp_traces = shift_phi_traces(temp_traces, -helio_shift)

            # add this trace segment neglecting the first point since it duplicates the last point of the previous segment.
            for i, trace in zip(np.where(midboundary_mask)[0], temp_traces):
                final_traces[i] = np.concatenate([final_traces[i], trace[:, 1:]], axis=1)

            # check the trace, update the end positions (SHIFTED BACK!) and the midboundary flag, continue the loop
            traced_to_boundary[midboundary_mask] = traces_.traced_to_boundary
            radial_end_pos[:, midboundary_mask] = shift_phi_lps(np.copy(traces_.end_pos), -helio_shift)
            midboundary_mask = np.isclose(radial_end_pos[0, :], r_interface, rtol=rtol)
            boundary_recross |= midboundary_mask

        # if no more work to be done, break the loop
        else:
            break
        maxiter -= 1

    # return the final traces and tracing checks
    return final_traces, traced_to_boundary, boundary_recross
