import gc

import pytest


def test_tracer_module_import(tracer_cls):
    assert hasattr(tracer_cls, 'Tracer'), "Tracer class not found in mapflpy.tracer"


def test_class_instantiation(Tracer):
    t = Tracer()
    assert isinstance(t, Tracer)
    del t
    gc.collect()


def test_instance_enforcement(Tracer):
    t1 = Tracer()
    try:
        with pytest.raises(RuntimeError):
            _ = Tracer()
    finally:
        del t1
        gc.collect()


def test_instance_enforcement_with_removal(Tracer):
    t1 = Tracer()
    del t1
    gc.collect()
    t2 = Tracer()
    assert isinstance(t2, Tracer)
    del t2
    gc.collect()


def test_tracer_and_tracermp_coexistence(Tracer, TracerMP):
    t = Tracer()
    try:
        with TracerMP() as tmp:
            assert isinstance(t, Tracer)
            assert isinstance(tmp, TracerMP)
    finally:
        del t
        gc.collect()


def test_tracermp_then_tracer_have_independent_instances(Tracer, TracerMP):
    with TracerMP() as tmp:
        t = Tracer()
        try:
            assert t.mapfl_id != tmp.mapfl_id
        finally:
            del t
            gc.collect()


def test_tracer_then_tracermp_have_independent_instances(Tracer, TracerMP):
    t = Tracer()
    try:
        with TracerMP() as tmp:
            assert t.mapfl_id != tmp.mapfl_id
    finally:
        del t
        gc.collect()


def test_initial_state(tracer_instance, defaults):
    DEFAULT_PARAMS, DEFAULT_FIELDS, _ = defaults
    merged = {**DEFAULT_PARAMS, **DEFAULT_FIELDS}
    assert tracer_instance.mapfl_id is not None
    assert tracer_instance._stale is True
    assert dict(tracer_instance._mapfl_params) == merged
