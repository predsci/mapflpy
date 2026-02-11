# tests/conftest.py
import gc
import shutil
import multiprocessing as mp
from _weakrefset import WeakSet
from pathlib import Path
from types import MappingProxyType

import numpy as np
import pytest

from psi_io import wrhdf_3d

from tests.utils import (_MESH_RESOLUTION, _TRACING_DIRECTION, _OLD_MAS, _BASE_MESH_PARAMS,
                         _BASE_LPS_PARAMS, _MESH_RESO_MAPPING, _TRACING_DIR_MAPPING,
                         _OLD_MAS_MAPPING, _DOMAIN_RANGES)


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
def reference_traces():
    ref_path = Path(__file__).parent / "data" / "reference_traces.npz"
    ref_traces = np.load(ref_path, allow_pickle=True)
    return ref_traces


@pytest.fixture(scope="session",
                params=_MESH_RESOLUTION)
def mesh_resolution(request):
    return request.param


@pytest.fixture(scope="session",
                params=_TRACING_DIRECTION)
def tracing_direction(request):
    return request.param


@pytest.fixture(scope="session",
                params=_OLD_MAS)
def mas_type(request):
    return request.param


@pytest.fixture(scope="session")
def tmp_data_dir(tmp_path_factory):
    # One temp directory for the whole test session
    return tmp_path_factory.mktemp("test_data_dir")


@pytest.fixture(scope="session")
def _default_fields_cached(dipole_field_factory, mas_type):
    return dipole_field_factory(**_BASE_MESH_PARAMS, **_MESH_RESO_MAPPING["normal"], **_OLD_MAS_MAPPING[mas_type])


@pytest.fixture(scope="session")
def _mesh_fields_cached(tmp_data_dir, dipole_field_factory,
                        mesh_resolution, mas_type):
    br, bt, bp, r, t, p = dipole_field_factory(
        **_BASE_MESH_PARAMS,
        **_MESH_RESO_MAPPING[mesh_resolution],
        **_OLD_MAS_MAPPING[mas_type]
    )
    data_dir = tmp_data_dir / mesh_resolution / mas_type
    data_dir.mkdir(parents=True, exist_ok=True)
    for dim, data in zip(['br', 'bt', 'bp'], [br, bt, bp]):
        filepath = data_dir / f"{dim}.h5"
        if not filepath.exists():
            wrhdf_3d(str(filepath), r, t, p, data)
    return tuple(data_dir / f"{dim}.h5" for dim in ['br', 'bt', 'bp']), (br, bt, bp, r, t, p)


@pytest.fixture(scope="session")
def _launch_points_cached(launch_point_factory, tracing_direction):
    lps = launch_point_factory(
        **_BASE_LPS_PARAMS,
        **_TRACING_DIR_MAPPING[tracing_direction],
    )
    return lps


@pytest.fixture(scope="session")
def interdomain_files(tmp_data_dir, dipole_field_factory,
                      mesh_resolution, mas_type):
    base_params = {
        **_BASE_MESH_PARAMS,
        **_MESH_RESO_MAPPING[mesh_resolution],
        **_OLD_MAS_MAPPING[mas_type]
    }
    cor_params = {**base_params, **_DOMAIN_RANGES['cor']}
    hel_params = {**base_params, **_DOMAIN_RANGES['hel']}
    
    br_cor, bt_cor, bp_cor, r_cor, t_cor, p_cor = dipole_field_factory(**cor_params)
    br_hel, bt_hel, bp_hel, r_hel, t_hel, p_hel = dipole_field_factory(**hel_params)
    
    data_dir = tmp_data_dir / mesh_resolution / mas_type
    data_dir.mkdir(parents=True, exist_ok=True)
    for dim, data in zip(['br_cor', 'bt_cor', 'bp_cor'], [br_cor, bt_cor, bp_cor]):
        filepath = data_dir / f"{dim}.h5"
        if not filepath.exists():
            wrhdf_3d(str(filepath), r_cor, t_cor, p_cor, data)
    for dim, data in zip(['br_hel', 'bt_hel', 'bp_hel'], [br_hel, bt_hel, bp_hel]):
        filepath = data_dir / f"{dim}.h5"
        if not filepath.exists():
            wrhdf_3d(str(filepath), r_hel, t_hel, p_hel, data)
    return tuple(data_dir / f"{dim}_{dom}.h5" for dom in ['cor', 'hel'] for dim in ['br', 'bt', 'bp'])


@pytest.fixture(scope="session")
def old_and_new_mas(tmp_data_dir, dipole_field_factory, mesh_resolution):
    data_dir_old_mas = tmp_data_dir / mesh_resolution /'old_mas'
    data_dir_new_mas = tmp_data_dir / mesh_resolution /'new_mas'
    if any(not (data_dir_old_mas / f'{dim}.h5').exists() for dim in ('br', 'bt', 'bp')):
        br_old, bt_old, bp_old, r_old, t_old, p_old = dipole_field_factory(
            **_BASE_MESH_PARAMS,
            **_MESH_RESO_MAPPING[mesh_resolution],
            **_OLD_MAS_MAPPING['old_mas']
        )
        data_dir_old_mas.mkdir(parents=True, exist_ok=True)
        for dim, data in zip(['br', 'bt', 'bp'], [br_old, bt_old, bp_old]):
            fp = data_dir_old_mas / f"{dim}.h5"
            wrhdf_3d(str(fp), r_old, t_old, p_old, data)

    if any(not (data_dir_new_mas / f'{dim}.h5').exists() for dim in ('br', 'bt', 'bp')):
        br_new, bt_new, bp_new, r_new, t_new, p_new = dipole_field_factory(
            **_BASE_MESH_PARAMS,
            **_MESH_RESO_MAPPING[mesh_resolution],
            **_OLD_MAS_MAPPING['new_mas']
        )
        data_dir_new_mas.mkdir(parents=True, exist_ok=True)
        for dim, data in zip(['br', 'bt', 'bp'], [br_new, bt_new, bp_new]):
            fp = data_dir_new_mas / f"{dim}.h5"
            wrhdf_3d(str(fp), r_new, t_new, p_new, data)

    return (tuple(data_dir_old_mas / f"{dim}.h5" for dim in ['br', 'bt', 'bp']),
           tuple(data_dir_new_mas / f"{dim}.h5" for dim in ['br', 'bt', 'bp']))


@pytest.fixture
def mesh_fields_asarray(_mesh_fields_cached):
    _, fields = _mesh_fields_cached
    # Return fresh copies of arrays for each test
    copied = tuple(np.copy(a) if hasattr(a, "dtype") else a for a in fields)
    return copied


@pytest.fixture
def mesh_fields_aspaths(_mesh_fields_cached):
    paths, _= _mesh_fields_cached
    return tuple(str(p) for p in paths)


@pytest.fixture
def launch_points(_launch_points_cached):
    lps = _launch_points_cached
    # Return fresh copies of arrays for each test
    copied = np.copy(lps)
    return copied


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
    return data_dir


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
