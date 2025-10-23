import pytest

from mapflpy.globals import Traces


def test_run_without_setting_fields(tracermp_instance):
    with pytest.raises(TypeError):
        tracermp_instance.run()



def test_run_with_setting_fields(tracermp_instance, default_fields_aspaths):
    br_in, bt_in, bp_in = default_fields_aspaths
    tracermp_instance.br = br_in
    tracermp_instance.bt = bt_in
    tracermp_instance.bp = bp_in
    tracermp_instance.run()
    assert tracermp_instance.stale is False


def test_trace_without_fieldline_validation(tracermp_instance, default_fields_aspaths):
    br_in, bt_in, bp_in = default_fields_aspaths
    tracermp_instance.br = br_in
    tracermp_instance.bt = bt_in
    tracermp_instance.bp = bp_in
    tracermp_instance.run()
    response = tracermp_instance.trace(launch_points=[1, 1, 1], buffer_size=10)
    assert isinstance(response, Traces)