import pytest
from numpy.testing import assert_allclose

from tests.utils import compute_weighted_fieldline_difference, _BUFFER_SIZE, _TEST_TOLERANCES
from mapflpy.utils import trim_fieldline_nan_buffer, combine_fwd_bwd_traces


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_tracing_against_reference_traces(tracer_instance,
                                          mesh_resolution,
                                          mas_type,
                                          tracing_direction,
                                          mesh_fields_asarray,
                                          launch_points,
                                          reference_traces):
    br, bt, bp, r, t, p = mesh_fields_asarray
    lps = launch_points
    tracer_instance.br = br, r, t, p
    tracer_instance.bt = bt, r, t, p
    tracer_instance.bp = bp, r, t, p

    match tracing_direction:
        case 'fwd':
            tracer_instance.set_tracing_direction('f')
            traces = tracer_instance.trace(lps, buffer_size=_BUFFER_SIZE)
        case 'bwd':
            tracer_instance.set_tracing_direction('b')
            traces = tracer_instance.trace(lps, buffer_size=_BUFFER_SIZE)
        case 'both':
            tracer_instance.set_tracing_direction('f')
            fwd_traces = tracer_instance.trace(lps, buffer_size=_BUFFER_SIZE)
            tracer_instance.set_tracing_direction('b')
            bwd_traces = tracer_instance.trace(lps, buffer_size=_BUFFER_SIZE)
            traces = combine_fwd_bwd_traces(fwd_traces, bwd_traces)
        case _:
            raise ValueError(f'Unknown launch points id: {tracing_direction}')

    traces_trimmed = trim_fieldline_nan_buffer(traces)
    for i, arr in enumerate(traces_trimmed):
        wdist = compute_weighted_fieldline_difference(arr, reference_traces[f'{mesh_resolution}_{mas_type}_{tracing_direction}_{i}'])
        assert_allclose(wdist, 0, atol=_TEST_TOLERANCES['atol_exact'])


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_old_mas_traces_against_new_mas_traces(tracer_instance,
                                               mesh_resolution,
                                               tracing_direction,
                                               old_and_new_mas,
                                               launch_points,
                                               reference_traces):
    mas_old, mas_new = old_and_new_mas
    tracer_instance.load_fields(*mas_old)

    match tracing_direction:
        case 'fwd':
            tracer_instance.set_tracing_direction('f')
            traces = tracer_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
        case 'bwd':
            tracer_instance.set_tracing_direction('b')
            traces = tracer_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
        case 'both':
            tracer_instance.set_tracing_direction('f')
            fwd_traces = tracer_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
            tracer_instance.set_tracing_direction('b')
            bwd_traces = tracer_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
            traces = combine_fwd_bwd_traces(fwd_traces, bwd_traces)
        case _:
            raise ValueError(f'Unknown launch points id: {tracing_direction}')

    traces_trimmed = trim_fieldline_nan_buffer(traces)
    for i, arr in enumerate(traces_trimmed):
        wdist = compute_weighted_fieldline_difference(arr, reference_traces[f'{mesh_resolution}_new_mas_{tracing_direction}_{i}'])
        assert_allclose(wdist, 0, atol=_TEST_TOLERANCES['atol_fuzzy'])


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_new_mas_traces_against_old_mas_traces(tracer_instance,
                                               mesh_resolution,
                                               tracing_direction,
                                               old_and_new_mas,
                                               launch_points,
                                               reference_traces):
    mas_old, mas_new = old_and_new_mas
    tracer_instance.load_fields(*mas_new)

    match tracing_direction:
        case 'fwd':
            tracer_instance.set_tracing_direction('f')
            traces = tracer_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
        case 'bwd':
            tracer_instance.set_tracing_direction('b')
            traces = tracer_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
        case 'both':
            tracer_instance.set_tracing_direction('f')
            fwd_traces = tracer_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
            tracer_instance.set_tracing_direction('b')
            bwd_traces = tracer_instance.trace(launch_points, buffer_size=_BUFFER_SIZE)
            traces = combine_fwd_bwd_traces(fwd_traces, bwd_traces)
        case _:
            raise ValueError(f'Unknown launch points id: {tracing_direction}')

    traces_trimmed = trim_fieldline_nan_buffer(traces)
    for i, arr in enumerate(traces_trimmed):
        wdist = compute_weighted_fieldline_difference(arr, reference_traces[f'{mesh_resolution}_old_mas_{tracing_direction}_{i}'])
        assert_allclose(wdist, 0, atol=_TEST_TOLERANCES['atol_fuzzy'])