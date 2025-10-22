import pytest


def test_tracer_module_import(tracer_cls):
    assert hasattr(tracer_cls, 'TracerMP'), "TracerMP class not found in mapflpy.tracer"


def test_class_instantiation(TracerMP):
    with TracerMP() as t:
        assert isinstance(t, TracerMP)


def test_multi_instance(TracerMP):
    with TracerMP() as t1, TracerMP() as t2:
        assert isinstance(t1, TracerMP)
        assert isinstance(t2, TracerMP)


def test_mapflpy_fortran_uniqueness(tracermp_pair):
    t1, t2 = tracermp_pair
    assert t1.mapfl_id != t2.mapfl_id, "TracerMP instances should use distinct mapfl instances"


def test_initial_state(tracermp_instance, defaults):
    DEFAULT_PARAMS, _, MAGNETIC_FIELD_PATHS = defaults
    for key in ('mapfl_id', '_parent', '_child', '_process'):
        assert hasattr(tracermp_instance, key)
        assert getattr(tracermp_instance, key) is not None
    assert tracermp_instance._stale is True
    assert dict(tracermp_instance._mapfl_params) == dict(DEFAULT_PARAMS)
    assert tracermp_instance._magnetic_fields == MAGNETIC_FIELD_PATHS
    assert float(tracermp_instance._timeout) == pytest.approx(30.0)