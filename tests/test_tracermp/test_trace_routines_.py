import pytest
from numpy.testing import assert_allclose

from mapflpy.scripts import run_forward_tracing, run_backward_tracing, run_fwdbwd_tracing, inter_domain_tracing
from tests.utils import compute_fieldline_length, compute_weighted_fieldline_difference, _BUFFER_SIZE, _TEST_TOLERANCES, _DOMAIN_RANGES
from mapflpy.utils import trim_fieldline_nan_buffer, combine_fwd_bwd_traces


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_tracing_against_reference_traces(tracermp_instance,
                                          mesh_resolution,
                                          mas_type,
                                          tracing_direction,
                                          mesh_fields_aspaths, 
                                          launch_points, 
                                          reference_traces):
    br, bt, bp = mesh_fields_aspaths
    lps = launch_points
    tracermp_instance.br = br
    tracermp_instance.bt = bt
    tracermp_instance.bp = bp

    match tracing_direction:
        case 'fwd':
            tracermp_instance.set_tracing_direction('f')
            traces = tracermp_instance.trace(lps, buffer_size=_BUFFER_SIZE)
        case 'bwd':
            tracermp_instance.set_tracing_direction('b')
            traces = tracermp_instance.trace(lps, buffer_size=_BUFFER_SIZE)
        case 'both':
            tracermp_instance.set_tracing_direction('f')
            fwd_traces = tracermp_instance.trace(lps, buffer_size=_BUFFER_SIZE)
            tracermp_instance.set_tracing_direction('b')
            bwd_traces = tracermp_instance.trace(lps, buffer_size=_BUFFER_SIZE)
            traces = combine_fwd_bwd_traces(fwd_traces, bwd_traces)
        case _:
            raise ValueError(f'Unknown launch points id: {tracing_direction}')

    traces_trimmed = trim_fieldline_nan_buffer(traces)
    for i, arr in enumerate(traces_trimmed):
        wdist = compute_weighted_fieldline_difference(arr, reference_traces[f'{mesh_resolution}_{mas_type}_{tracing_direction}_{i}'])
        assert_allclose(wdist, 0, atol=_TEST_TOLERANCES['atol_exact'])
        

@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_old_mas_traces_against_new_mas_traces(tracermp_instance,
                                               mesh_resolution,
                                               tracing_direction,
                                               old_and_new_mas,
                                               launch_points,
                                               reference_traces):
    mas_old, mas_new = old_and_new_mas
    tracermp_instance.load_fields(*mas_old)

    match tracing_direction:
        case 'fwd':
            tracermp_instance.set_tracing_direction('f')
            traces = tracermp_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
        case 'bwd':
            tracermp_instance.set_tracing_direction('b')
            traces = tracermp_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
        case 'both':
            tracermp_instance.set_tracing_direction('f')
            fwd_traces = tracermp_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
            tracermp_instance.set_tracing_direction('b')
            bwd_traces = tracermp_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
            traces = combine_fwd_bwd_traces(fwd_traces, bwd_traces)
        case _:
            raise ValueError(f'Unknown launch points id: {tracing_direction}')

    traces_trimmed = trim_fieldline_nan_buffer(traces)
    for i, arr in enumerate(traces_trimmed):
        wdist = compute_weighted_fieldline_difference(arr, reference_traces[f'{mesh_resolution}_new_mas_{tracing_direction}_{i}'])
        assert_allclose(wdist, 0, atol=_TEST_TOLERANCES['atol_fuzzy'])


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_new_mas_traces_against_old_mas_traces(tracermp_instance,
                                               mesh_resolution,
                                               tracing_direction,
                                               old_and_new_mas,
                                               launch_points,
                                               reference_traces):
    mas_old, mas_new = old_and_new_mas
    tracermp_instance.load_fields(*mas_new)

    match tracing_direction:
        case 'fwd':
            tracermp_instance.set_tracing_direction('f')
            traces = tracermp_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
        case 'bwd':
            tracermp_instance.set_tracing_direction('b')
            traces = tracermp_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
        case 'both':
            tracermp_instance.set_tracing_direction('f')
            fwd_traces = tracermp_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
            tracermp_instance.set_tracing_direction('b')
            bwd_traces = tracermp_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
            traces = combine_fwd_bwd_traces(fwd_traces, bwd_traces)
        case _:
            raise ValueError(f'Unknown launch points id: {tracing_direction}')

    traces_trimmed = trim_fieldline_nan_buffer(traces)
    for i, arr in enumerate(traces_trimmed):
        wdist = compute_weighted_fieldline_difference(arr, reference_traces[f'{mesh_resolution}_old_mas_{tracing_direction}_{i}'])
        assert_allclose(wdist, 0, atol=_TEST_TOLERANCES['atol_fuzzy'])


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_tracing_scripts_against_reference_traces(mesh_resolution,
                                                  mas_type,
                                                  tracing_direction,
                                                  mesh_fields_aspaths,
                                                  launch_points,
                                                  reference_traces):
    br, bt, bp = mesh_fields_aspaths
    lps = launch_points

    match tracing_direction:
        case 'fwd':
            traces = run_forward_tracing(br, bt, bp, lps, buffer_size=_BUFFER_SIZE)
        case 'bwd':
            traces = run_backward_tracing(br, bt, bp, lps, buffer_size=_BUFFER_SIZE)
        case 'both':
            traces = run_fwdbwd_tracing(br, bt, bp, lps, buffer_size=_BUFFER_SIZE)
        case _:
            raise ValueError(f'Unknown launch points id: {tracing_direction}')

    traces_trimmed = trim_fieldline_nan_buffer(traces)
    for i, arr in enumerate(traces_trimmed):
        wdist = compute_weighted_fieldline_difference(arr, reference_traces[f'{mesh_resolution}_{mas_type}_{tracing_direction}_{i}'])
        assert_allclose(wdist, 0, atol=_TEST_TOLERANCES['atol_exact'])


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_interdomain_tracing_against_reference_traces(mesh_resolution,
                                                      mas_type,
                                                      tracing_direction,
                                                      interdomain_files,
                                                      launch_points,
                                                      reference_traces):
    """
    This test compares an interdomain trace where the domain has been split at r_interface to the reference traces
    where there was no split-domain. These traces can differ more because exactly where the interface lies along
    a reference field line segment can vary, which is effectively like seeding a part (or parts) of the trace
    with a slightly different start location. This means the traces *will not* be the same length and we must
    use the fuzzy tolerances because the errors are related to the mesh and discretization of B (not the tracer itself).
    """
    br_cor, bt_cor, bp_cor, br_hel, bt_hel, bp_hel = interdomain_files
    lps = launch_points
    buffer = _BUFFER_SIZE

    assert _DOMAIN_RANGES['cor']['r1'] == _DOMAIN_RANGES['hel']['r0'], \
        "Inconsistent domain interface radii."
    r_interface = _DOMAIN_RANGES['cor']['r1']

    if tracing_direction != 'both':
        pytest.skip("Interdomain tracing only implemented for both directions.")
    else:
        traces, *_ = inter_domain_tracing(br_cor,
                                          bt_cor,
                                          bp_cor,
                                          br_hel,
                                          bt_hel,
                                          bp_hel,
                                          lps,
                                          r_interface=r_interface,
                                          buffer_size=buffer)
    for i, arr in enumerate(traces):
        reference_trace = reference_traces[f'{mesh_resolution}_{mas_type}_{tracing_direction}_{i}']
        # compare the distance of the first and last points (footprints)
        wdist = compute_weighted_fieldline_difference(arr[:, [0, -1]], reference_trace[:, [0, -1]])
        assert_allclose(wdist, 0, atol=_TEST_TOLERANCES['atol_fuzzy'])
        # compare the lengths of the traces
        len_test = compute_fieldline_length(arr)
        len_ref = compute_fieldline_length(reference_trace)
        assert_allclose(len_test, len_ref, rtol=_TEST_TOLERANCES['rtol_fuzzy'])
