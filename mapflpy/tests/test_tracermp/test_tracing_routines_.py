import pytest
from numpy.testing import assert_allclose

from mapflpy.scripts import run_foward_tracing, run_backward_tracing, run_fwdbwd_tracing, inter_domain_tracing
from mapflpy.tests.utils import compute_fieldline_length
from mapflpy.utils import trim_fieldline_nan_buffer, combine_fwd_bwd_traces


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_tracing_against_reference_traces(tracermp_instance, mesh_fields_aspaths, launch_points, default_lps_params, reference_traces):
    mesh_id, (br, bt, bp) = mesh_fields_aspaths
    lps_id, lps = launch_points
    tracermp_instance.br = br
    tracermp_instance.bt = bt
    tracermp_instance.bp = bp

    match lps_id:
        case 'fwd':
            tracermp_instance.set_tracing_direction('f')
            traces = tracermp_instance.trace(lps, buffer_size=default_lps_params['BUFFER'])
        case 'bwd':
            tracermp_instance.set_tracing_direction('b')
            traces = tracermp_instance.trace(lps, buffer_size=default_lps_params['BUFFER'])
        case 'both':
            tracermp_instance.set_tracing_direction('f')
            fwd_traces = tracermp_instance.trace(lps, buffer_size=default_lps_params['BUFFER'])
            tracermp_instance.set_tracing_direction('b')
            bwd_traces = tracermp_instance.trace(lps, buffer_size=default_lps_params['BUFFER'])
            traces = combine_fwd_bwd_traces(fwd_traces, bwd_traces)
        case _:
            raise ValueError(f'Unknown launch points id: {lps_id}')

    traces_trimmed = trim_fieldline_nan_buffer(traces)
    for i, arr in enumerate(traces_trimmed):
        assert_allclose(arr, reference_traces[f'{mesh_id}_{lps_id}_{i}'])


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_tracing_scripts_against_reference_traces(mesh_fields_aspaths, launch_points, default_lps_params, reference_traces):
    mesh_id, (br, bt, bp) = mesh_fields_aspaths
    lps_id, lps = launch_points

    match lps_id:
        case 'fwd':
            traces = run_foward_tracing(br, bt, bp, lps, buffer_size=default_lps_params['BUFFER'])
        case 'bwd':
            traces = run_backward_tracing(br, bt, bp, lps, buffer_size=default_lps_params['BUFFER'])
        case 'both':
            traces = run_fwdbwd_tracing(br, bt, bp, lps, buffer_size=default_lps_params['BUFFER'])
        case _:
            raise ValueError(f'Unknown launch points id: {lps_id}')

    traces_trimmed = trim_fieldline_nan_buffer(traces)
    for i, arr in enumerate(traces_trimmed):
        assert_allclose(arr, reference_traces[f'{mesh_id}_{lps_id}_{i}'])


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_interdomain_tracing_against_reference_traces(interdomain_files, launch_points, default_params, reference_traces):
    mesh_id, (br_cor, bt_cor, bp_cor, br_hel, bt_hel, bp_hel) = interdomain_files
    lps_id, lps = launch_points
    buffer = default_params['lps']['BUFFER']

    assert default_params['_domains']['cor']['r1'] == default_params['_domains']['hel']['r0'], "Inconsistent domain interface radii."
    r_interface = default_params['_domains']['cor']['r1']

    if lps_id != 'both':
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
        len_test = compute_fieldline_length(arr)
        len_ref = compute_fieldline_length(reference_traces[f'{mesh_id}_{lps_id}_{i}'])
        assert_allclose(arr[:, [0, -1]], reference_traces[f'{mesh_id}_{lps_id}_{i}'][:, [0, -1]], atol=1e-3)
        assert_allclose(len_test, len_ref, atol=1e-3)
