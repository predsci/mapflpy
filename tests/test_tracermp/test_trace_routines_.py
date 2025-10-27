import pytest
from numpy.testing import assert_allclose

from mapflpy.scripts import run_forward_tracing, run_backward_tracing, run_fwdbwd_tracing, inter_domain_tracing
from tests.utils import compute_fieldline_length, compute_weighted_fieldline_difference
from mapflpy.utils import trim_fieldline_nan_buffer, combine_fwd_bwd_traces


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_tracing_against_reference_traces(tracermp_instance, mesh_fields_aspaths, launch_points, default_params, reference_traces):
    mesh_id, (br, bt, bp) = mesh_fields_aspaths
    lps_id, lps = launch_points
    tracermp_instance.br = br
    tracermp_instance.bt = bt
    tracermp_instance.bp = bp

    match lps_id:
        case 'fwd':
            tracermp_instance.set_tracing_direction('f')
            traces = tracermp_instance.trace(lps, buffer_size=default_params['lps']['BUFFER'])
        case 'bwd':
            tracermp_instance.set_tracing_direction('b')
            traces = tracermp_instance.trace(lps, buffer_size=default_params['lps']['BUFFER'])
        case 'both':
            tracermp_instance.set_tracing_direction('f')
            fwd_traces = tracermp_instance.trace(lps, buffer_size=default_params['lps']['BUFFER'])
            tracermp_instance.set_tracing_direction('b')
            bwd_traces = tracermp_instance.trace(lps, buffer_size=default_params['lps']['BUFFER'])
            traces = combine_fwd_bwd_traces(fwd_traces, bwd_traces)
        case _:
            raise ValueError(f'Unknown launch points id: {lps_id}')

    traces_trimmed = trim_fieldline_nan_buffer(traces)
    for i, arr in enumerate(traces_trimmed):
        wdist = compute_weighted_fieldline_difference(arr, reference_traces[f'{mesh_id}_{lps_id}_{i}'])
        assert_allclose(wdist, 0, atol=default_params['_testing']['tolerances']['atol_exact'])


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_tracing_scripts_against_reference_traces(mesh_fields_aspaths, launch_points, default_params, reference_traces):
    mesh_id, (br, bt, bp) = mesh_fields_aspaths
    lps_id, lps = launch_points

    match lps_id:
        case 'fwd':
            traces = run_forward_tracing(br, bt, bp, lps, buffer_size=default_params['lps']['BUFFER'])
        case 'bwd':
            traces = run_backward_tracing(br, bt, bp, lps, buffer_size=default_params['lps']['BUFFER'])
        case 'both':
            traces = run_fwdbwd_tracing(br, bt, bp, lps, buffer_size=default_params['lps']['BUFFER'])
        case _:
            raise ValueError(f'Unknown launch points id: {lps_id}')

    traces_trimmed = trim_fieldline_nan_buffer(traces)
    for i, arr in enumerate(traces_trimmed):
        wdist = compute_weighted_fieldline_difference(arr, reference_traces[f'{mesh_id}_{lps_id}_{i}'])
        assert_allclose(wdist, 0, atol=default_params['_testing']['tolerances']['atol_exact'])

@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_interdomain_tracing_against_reference_traces(interdomain_files, launch_points, default_params, reference_traces):
    """
    This test compares an interdomain trace where the domain has been split at r_interface to the reference traces
    where there was no split-domain. These traces can differ more because exactly where the interface lies along
    a reference field line segment can vary, which is effectively like seeding a part (or parts) of the trace
    with a slightly different start location. This means the traces *will not* be the same length and we must
    use the fuzzy tolerances because the errors are related to the mesh and discretization of B (not the tracer itself).
    """
    mesh_id, (br_cor, bt_cor, bp_cor, br_hel, bt_hel, bp_hel) = interdomain_files
    lps_id, lps = launch_points
    buffer = default_params['lps']['BUFFER']

    assert default_params["_testing"]["domain_ranges"]['cor']['r1'] == default_params["_testing"]["domain_ranges"]['hel']['r0'], \
        "Inconsistent domain interface radii."
    r_interface = default_params["_testing"]["domain_ranges"]['cor']['r1']

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
        # compare the distance of the first and last points (footprints)
        wdist = compute_weighted_fieldline_difference(arr[:, [0, -1]], reference_traces[f'{mesh_id}_{lps_id}_{i}'][:, [0, -1]])
        assert_allclose(wdist, 0, atol=default_params["_testing"]["tolerances"]['atol_fuzzy'])
        # compare the lengths of the traces
        len_test = compute_fieldline_length(arr)
        len_ref = compute_fieldline_length(reference_traces[f'{mesh_id}_{lps_id}_{i}'])
        assert_allclose(len_test, len_ref, rtol=default_params["_testing"]["tolerances"]['rtol_fuzzy'])
