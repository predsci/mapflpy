import pytest
from numpy.testing import assert_allclose

from mapflpy.globals import Traces


def test_run_without_setting_fields(tracer_instance):
    with pytest.raises(TypeError):
        tracer_instance.run()


def test_run_with_setting_fields(tracer_instance, default_fields_asarrays):
    br_in, bt_in, bp_in, r_in, t_in, p_in = default_fields_asarrays
    tracer_instance.br = br_in, r_in, t_in, p_in
    tracer_instance.bt = bt_in, r_in, t_in, p_in
    tracer_instance.bp = bp_in, r_in, t_in, p_in
    tracer_instance.run()
    assert tracer_instance.stale is False


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_trace_without_fieldline_validation(tracer_instance, default_fields_asarrays):
    br_in, bt_in, bp_in, r_in, t_in, p_in = default_fields_asarrays
    tracer_instance.br = br_in, r_in, t_in, p_in
    tracer_instance.bt = bt_in, r_in, t_in, p_in
    tracer_instance.bp = bp_in, r_in, t_in, p_in
    tracer_instance.run()
    response = tracer_instance.trace(launch_points=[1, 1, 1], buffer_size=10)
    assert isinstance(response, Traces)


