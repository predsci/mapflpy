import pytest
from numpy.testing import assert_allclose

from tests.utils import compute_weighted_fieldline_difference
from mapflpy.utils import trim_fieldline_nan_buffer, combine_fwd_bwd_traces


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_tracing_against_reference_traces(tracer_instance, mesh_fields_asarray, launch_points, default_params, reference_traces):
    mesh_id, (br, bt, bp, r, t, p) = mesh_fields_asarray
    lps_id, lps = launch_points
    tracer_instance.br = br, r, t, p
    tracer_instance.bt = bt, r, t, p
    tracer_instance.bp = bp, r, t, p

    match lps_id:
        case 'fwd':
            tracer_instance.set_tracing_direction('f')
            traces = tracer_instance.trace(lps, buffer_size=default_params['lps']['BUFFER'])
        case 'bwd':
            tracer_instance.set_tracing_direction('b')
            traces = tracer_instance.trace(lps, buffer_size=default_params['lps']['BUFFER'])
        case 'both':
            tracer_instance.set_tracing_direction('f')
            fwd_traces = tracer_instance.trace(lps, buffer_size=default_params['lps']['BUFFER'])
            tracer_instance.set_tracing_direction('b')
            bwd_traces = tracer_instance.trace(lps, buffer_size=default_params['lps']['BUFFER'])
            traces = combine_fwd_bwd_traces(fwd_traces, bwd_traces)
        case _:
            raise ValueError(f'Unknown launch points id: {lps_id}')

    traces_trimmed = trim_fieldline_nan_buffer(traces)
    for i, arr in enumerate(traces_trimmed):
        wdist = compute_weighted_fieldline_difference(arr, reference_traces[f'{mesh_id}_{lps_id}_{i}'])
        assert_allclose(wdist, 0, atol=default_params['_testing']['tolerances']['atol_exact'])
