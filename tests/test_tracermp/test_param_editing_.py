import pytest
from numpy.testing import assert_allclose


@pytest.mark.parametrize("direction", ["f", "b"])
def test_params_setting_direction(tracermp_instance, direction):
    tracermp_instance.set_tracing_direction(direction)
    assert tracermp_instance['trace_fwd_'] == (direction == 'f')
    assert tracermp_instance['trace_bwd_'] == (direction == 'b')


@pytest.mark.parametrize("dimension", ["br", "bt", "bp"])
def test_setting_fields_with_arrays(tracermp_instance, default_fields_asarrays, dimension):
    br_in, bt_in, bp_in, r_in, t_in, p_in = default_fields_asarrays
    match dimension:
        case "br":
            with pytest.raises(ValueError):
                tracermp_instance.br = br_in, r_in, t_in, p_in
        case "bt":
            with pytest.raises(ValueError):
                tracermp_instance.bt = bt_in, r_in, t_in, p_in
        case "bp":
            with pytest.raises(ValueError):
                tracermp_instance.bt = bt_in, r_in, t_in, p_in
        case _:
            raise ValueError(f"Unknown dimension: {dimension}")


@pytest.mark.parametrize("dimension", ["br", "bt", "bp"])
def test_setting_fields_with_filepaths(tracermp_instance, default_fields_aspaths, dimension):
    br_filepath, bt_filepath, bp_filepath = default_fields_aspaths
    match dimension:
        case "br":
            tracermp_instance.br = br_filepath
            assert tracermp_instance.br == br_filepath
            assert tracermp_instance._magnetic_fields[dimension] == br_filepath
        case "bt":
            tracermp_instance.bt = bt_filepath
            assert tracermp_instance.bt == bt_filepath
            assert tracermp_instance._magnetic_fields[dimension] == bt_filepath
        case "bp":
            tracermp_instance.bp = bp_filepath
            assert tracermp_instance.bp == bp_filepath
            assert tracermp_instance._magnetic_fields[dimension] == bp_filepath
        case _:
            raise ValueError(f"Unknown dimension: {dimension}")


def test_setting_br_with_invalid_filepaths(tracermp_instance):
    with pytest.raises(FileNotFoundError):
        tracermp_instance.br = "non_existent_file.h5"


def test_chainmap_default_overloading(tracermp_instance):
    with pytest.raises(KeyError):
        _ = tracermp_instance['dsmult']

    dsmult = tracermp_instance['dsmult_']
    tracermp_instance['dsmult_'] = 123

    assert tracermp_instance['dsmult_'] == 123
    assert tracermp_instance._mapfl_params.maps[0]['dsmult_'] == 123
    assert tracermp_instance._mapfl_params.maps[1]['dsmult_'] == dsmult

    tracermp_instance.clear()
    assert tracermp_instance['dsmult_'] == dsmult
    assert tracermp_instance._mapfl_params.maps[1]['dsmult_'] == dsmult
    with pytest.raises(KeyError):
        _ = tracermp_instance._mapfl_params.maps[0]['dsmult_']