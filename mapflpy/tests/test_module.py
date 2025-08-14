# tests/test_tracer.py
import gc
import json

import pytest


def test_shared_object():
    from mapflpy.tests.utils import check_shared_object
    location = check_shared_object()
    assert location and isinstance(location, str)


def test_reference_tracers_meta(reference_traces, default_params):
    tracer_metadata = json.loads(reference_traces["__meta__"].item())
    tracer_defaults = tracer_metadata["defaults"]
    for key, value in default_params.items():
        if isinstance(value, dict):
            for subkey, subvalue in value.items():
                assert tracer_defaults[key][subkey] == subvalue
        else:
            assert tracer_defaults[key] == value
