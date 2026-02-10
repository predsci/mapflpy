# tests/test_tracer.py
import json


def test_version():
    from mapflpy import __version__
    assert isinstance(__version__, str)
    assert len(__version__) > 0
    assert __version__ != "0+unknown"


def test_shared_object():
    from tests.utils import check_shared_object
    location = check_shared_object()
    assert location and isinstance(location, str)
