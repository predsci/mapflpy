# tests/conftest.py
import gc
import shutil
import multiprocessing as mp
from _weakrefset import WeakSet
from pathlib import Path

import numpy as np
import pytest

from psi_io import wrhdf_3d


@pytest.fixture(autouse=True)
def reset_singleton(monkeypatch):
    # Ensure a clean singleton set before each test
    from mapflpy.tracer import Tracer
    monkeypatch.setattr(Tracer, "_instances", WeakSet())
    yield
    # Also clean after, in case a test leaves junk around
    monkeypatch.setattr(Tracer, "_instances", WeakSet())
    gc.collect()


@pytest.fixture(scope="session", autouse=True)
def _mp_start_method():
    try:
        mp.set_start_method("spawn")
    except RuntimeError:
        # already set by another test session
        pass
    yield


@pytest.fixture(scope="session")
def dipole_field_factory():
    from tests.utils import dipole_field
    return dipole_field


@pytest.fixture(scope="session")
def launch_point_factory():
    from tests.utils import _recreate_lps_from_defaults
    return _recreate_lps_from_defaults


@pytest.fixture(scope="session")
def default_params():
    from tests.utils import read_defaults
    return read_defaults()


@pytest.fixture(scope="session")
def default_mesh_params(default_params):
    return default_params["mesh"]


@pytest.fixture(scope="session")
def default_lps_params(default_params):
    return default_params["lps"]


@pytest.fixture(scope="session")
def reference_traces():
    ref_path = Path(__file__).parent / "data" / "reference_traces.npz"
    ref_traces = np.load(ref_path, allow_pickle=True)
    return ref_traces


@pytest.fixture(scope="session")
def _default_fields_cached(default_mesh_params):
    from tests.utils import dipole_field
    return dipole_field(**default_mesh_params["base"], **default_mesh_params["params"]["normal"])


@pytest.fixture(scope="session", params=["coarse", "normal", "fine"], ids=lambda x: x)
def _mesh_fields_cached(tmp_path_factory, default_mesh_params, dipole_field_factory, request):
    level = request.param
    br, bt, bp, r, t, p = dipole_field_factory(
        **default_mesh_params["base"],
        **default_mesh_params["params"][level],
    )
    data_dir = tmp_path_factory.mktemp(f"{level}_magnetic_field_files")
    for dim, data in zip(['br', 'bt', 'bp'], [br, bt, bp]):
        filepath = data_dir / f"{dim}.h5"
        wrhdf_3d(str(filepath), r, t, p, data)
    yield level, tuple(data_dir / f"{dim}.h5" for dim in ['br', 'bt', 'bp']), (br, bt, bp, r, t, p)
    shutil.rmtree(data_dir, ignore_errors=True)


@pytest.fixture(scope="session", params=["fwd", "bwd", "both"], ids=lambda x: x)
def _launch_points_cached(default_lps_params, launch_point_factory, request):
    level = request.param
    lps = launch_point_factory(
        **default_lps_params['base'],
        **default_lps_params['params'][level]
    )
    return level, lps


@pytest.fixture(scope="session", params=["coarse", "normal", "fine"], ids=lambda x: x)
def interdomain_files(tmp_path_factory, default_params, dipole_field_factory, request):
    level = request.param
    base_params = {
        **default_params["mesh"]["base"],
        **default_params["mesh"]["params"][level]
    }
    cor_params = {**base_params, **default_params["_testing"]["domain_ranges"]["cor"]}
    hel_params = {**base_params, **default_params["_testing"]["domain_ranges"]["hel"]}
    
    br_cor, bt_cor, bp_cor, r_cor, t_cor, p_cor = dipole_field_factory(**cor_params)
    br_hel, bt_hel, bp_hel, r_hel, t_hel, p_hel = dipole_field_factory(**hel_params)
    
    data_dir = tmp_path_factory.mktemp(f"{level}_magnetic_field_files")
    for dim, data in zip(['br_cor', 'bt_cor', 'bp_cor'], [br_cor, bt_cor, bp_cor]):
        filepath = data_dir / f"{dim}.h5"
        wrhdf_3d(str(filepath), r_cor, t_cor, p_cor, data)
    for dim, data in zip(['br_hel', 'bt_hel', 'bp_hel'], [br_hel, bt_hel, bp_hel]):
        filepath = data_dir / f"{dim}.h5"
        wrhdf_3d(str(filepath), r_hel, t_hel, p_hel, data)
    yield level, tuple(data_dir / f"{dim}_{dom}.h5" for dom in ['cor', 'hel'] for dim in ['br', 'bt', 'bp'])
    shutil.rmtree(data_dir, ignore_errors=True)


@pytest.fixture
def mesh_fields_asarray(_mesh_fields_cached):
    level, _, fields = _mesh_fields_cached
    # Return fresh copies of arrays for each test
    copied = tuple(np.copy(a) if hasattr(a, "dtype") else a for a in fields)
    return level, copied


@pytest.fixture
def mesh_fields_aspaths(_mesh_fields_cached):
    level, paths, _ = _mesh_fields_cached
    return level, tuple(str(p) for p in paths)


@pytest.fixture
def launch_points(_launch_points_cached):
    level, lps = _launch_points_cached
    # Return fresh copies of arrays for each test
    copied = np.copy(lps)
    return level, copied


@pytest.fixture
def default_fields_asarrays(_default_fields_cached):
    return tuple(np.copy(arr) for arr in _default_fields_cached)


@pytest.fixture(scope="session")
def _default_datadir(tmp_path_factory, _default_fields_cached):
    """
    Creates a directory with test files before any tests run.
    Cleans up after the whole test suite is done.
    """
    # Create a unique temporary directory for this session
    data_dir = tmp_path_factory.mktemp("magnetic_field_files")
    br, bt, bp, r, t, p = tuple(np.copy(arr) for arr in _default_fields_cached)
    for dim, data in zip(['br', 'bt', 'bp'], [br, bt, bp]):
        filepath = data_dir / f"{dim}.h5"
        wrhdf_3d(str(filepath), r, t, p, data)
    yield data_dir
    shutil.rmtree(data_dir, ignore_errors=True)


@pytest.fixture(scope="session")
def default_fields_aspaths(_default_datadir):
    br_path = str(_default_datadir / "br.h5")
    bt_path = str(_default_datadir / "bt.h5")
    bp_path = str(_default_datadir / "bp.h5")
    return br_path, bt_path, bp_path


@pytest.fixture
def tracer_cls():
    from mapflpy import tracer
    return tracer


@pytest.fixture
def Tracer(tracer_cls):
    return tracer_cls.Tracer


@pytest.fixture
def TracerMP(tracer_cls):
    return tracer_cls.TracerMP


@pytest.fixture
def defaults():
    from mapflpy.globals import DEFAULT_PARAMS, DEFAULT_FIELDS, MAGNETIC_FIELD_PATHS
    return DEFAULT_PARAMS, DEFAULT_FIELDS, MAGNETIC_FIELD_PATHS


@pytest.fixture
def tracer_instance(Tracer):
    # Create, yield, then cleanup to release the WeakSet singleton
    t = Tracer()
    try:
        yield t
    finally:
        del t
        gc.collect()


@pytest.fixture
def tracermp_instance(TracerMP):
    # Use context manager to ensure subprocess is closed
    with TracerMP() as t:
        yield t


@pytest.fixture
def tracermp_pair(TracerMP):
    with TracerMP() as t1, TracerMP() as t2:
        yield t1, t2
