import numpy as np
import pytest
from numpy.testing import assert_allclose


@pytest.mark.parametrize("direction", ["f", "b"])
def test_setting_direction(tracer_instance, direction):
    tracer_instance.set_tracing_direction(direction)
    assert tracer_instance['trace_fwd_'] == (direction == 'f')
    assert tracer_instance['trace_bwd_'] == (direction == 'b')


@pytest.mark.parametrize("dimension", ["br", "bt", "bp"])
def test_setting_fields_with_arrays(tracer_instance, default_fields_asarrays, dimension):
    br_in, bt_in, bp_in, r_in, t_in, p_in = default_fields_asarrays
    match dimension:
        case "br":
            bx_in = br_in
            tracer_instance.br = bx_in, r_in, t_in, p_in
            b_out, r_out, t_out, p_out = tracer_instance.br
        case "bt":
            bx_in = bt_in
            tracer_instance.bt = bx_in, r_in, t_in, p_in
            b_out, r_out, t_out, p_out = tracer_instance.bt
        case "bp":
            bx_in = bp_in
            tracer_instance.bp = bx_in, r_in, t_in, p_in
            b_out, r_out, t_out, p_out = tracer_instance.bp
        case _:
            raise ValueError(f"Unknown dimension: {dimension}")
    for arr1, arr2 in zip((b_out.T, r_out, t_out, p_out), (bx_in, r_in, t_in, p_in)):
        assert_allclose(arr1, arr2)


@pytest.mark.parametrize("dimension", ["br", "bt", "bp"])
def test_setting_fields_with_filepaths(tracer_instance, default_fields_aspaths, default_fields_asarrays, dimension):
    br_filepath, bt_filepath, bp_filepath = default_fields_aspaths
    br_in, bt_in, bp_in, r_in, t_in, p_in = default_fields_asarrays
    match dimension:
        case "br":
            bx_in = br_in
            tracer_instance.br = br_filepath
            b_out, r_out, t_out, p_out = tracer_instance.br
        case "bt":
            bx_in = bt_in
            tracer_instance.bt = bt_filepath
            b_out, r_out, t_out, p_out = tracer_instance.bt
        case "bp":
            bx_in = bp_in
            tracer_instance.bp = bp_filepath
            b_out, r_out, t_out, p_out = tracer_instance.bp
        case _:
            raise ValueError(f"Unknown dimension: {dimension}")
    for arr1, arr2 in zip((b_out.T, r_out, t_out, p_out), (bx_in, r_in, t_in, p_in)):
        assert_allclose(arr1, arr2)


def test_setting_br_with_valid_array_shapes(tracer_instance, default_fields_asarrays):
    br_in, *_ = default_fields_asarrays
    r_len, t_len, p_len = br_in.shape
    tracer_instance.br = br_in, np.full(p_len, 1), np.full(t_len, 1), np.full(r_len, 1)
    tracer_instance.br = br_in, np.full(r_len, 1), np.full(t_len, 1), np.full(p_len, 1)


def test_setting_br_with_invalid_array_shapes(tracer_instance, default_fields_asarrays):
    br_in, _, _, r_in, t_in, p_in = default_fields_asarrays
    with pytest.raises(ValueError):
        tracer_instance.br = br_in, np.full(1, 1), np.full(1, 1), np.full(1, 1)


def test_setting_br_with_invalid_filepaths(tracer_instance):
    with pytest.raises(FileNotFoundError):
        tracer_instance.br = "non_existent_file.h5"


def test_chainmap_default_overloading(tracer_instance):
    with pytest.raises(KeyError):
        _ = tracer_instance['dsmult']

    dsmult = tracer_instance['dsmult_']
    tracer_instance['dsmult_'] = 123

    assert tracer_instance['dsmult_'] == 123
    assert tracer_instance._mapfl_params.maps[0]['dsmult_'] == 123
    assert tracer_instance._mapfl_params.maps[1]['dsmult_'] == dsmult

    tracer_instance.clear()
    assert tracer_instance['dsmult_'] == dsmult
    assert tracer_instance._mapfl_params.maps[1]['dsmult_'] == dsmult
    with pytest.raises(KeyError):
        _ = tracer_instance._mapfl_params.maps[0]['dsmult_']

