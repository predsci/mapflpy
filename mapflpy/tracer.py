"""
Core class for managing tracing with mapflpy_fortran.
"""
import os
from collections import namedtuple, ChainMap
from types import MappingProxyType
from typing import Tuple, Literal, Iterable, Mapping

import numpy as np

from psi_io import read_hdf_by_value
from .utils import fetch_default_launch_points

# ------------------------------------------------------------------------------
# The following MAPFL_PARAMS are modeled after the mapfl.in file. However, the
# following key-value pairs are not exhaustive; rather, they are a *working*
# configuration that will suite most general uses of the Tracer object.
# ------------------------------------------------------------------------------

BFiles = namedtuple("BFiles", ('br', 'bt', 'bp'))
Traces = namedtuple('Traces', ['traces', 'end_pos', 'traced_to_boundary'])

DIRECTION = Literal['f', 'b']

DEFAULT_BUFFER_SIZE = 2000

FIELD_DEFAULTS = MappingProxyType({
    'br': None,
    'br_r': None,
    'br_nr': None,
    'br_t': None,
    'br_nt': None,
    'br_p': None,
    'br_np': None,
    'bt': None,
    'bt_r': None,
    'bt_nr': None,
    'bt_t': None,
    'bt_nt': None,
    'bt_p': None,
    'bt_np': None,
    'bp': None,
    'bp_r': None,
    'bp_nr': None,
    'bp_t': None,
    'bp_nt': None,
    'bp_p': None,
    'bp_np': None,
})




class Tracer:
    """Fieldline tracing wrapper class for mapflpy.

    Attributes
    ----------
    _mapfl : Python extension module
        The currently imported instance of the mapflpy_fortran.mapfl Fortran interface.
    _mapfl_params : dict
        Dictionary of key-value pairs to instantiate the mapflpy_fortran object. The
        keys used to instantiate mapflpy_fortran must have a trailing underscore.
    _launch_points : numpy.ndarray
        Launch points for fieldline tracing.
    _bfiles : Tuple[str, str, str])
        A tuple containing the paths to the 3D datafiles of the B field vector components Br, Bt, Bp.
    """
    _default_params = MappingProxyType({
        'integrate_along_fl_': False,
        'scalar_input_file_': '',
        'verbose_': False,
        'cubic_': False,
        'trace_fwd_': False,
        'trace_bwd_': False,
        'trace_3d_': False,
        'trace_slice_': False,
        'compute_ch_map_': False,
        'debug_level_': 0,
        'use_analytic_function_': False,
        'function_params_file_': ' ',
        'domain_r_min_': 1,
        'domain_r_max_': 300,
        'bfile_r_': 'br002.h5',
        'bfile_t_': 'bt002.h5',
        'bfile_p_': 'bp002.h5',
        'ds_variable_': True,
        'ds_over_rc_': 0.0025,
        'ds_min_': 0.00005,
        'ds_max_': 10.0,
        'ds_limit_by_local_mesh_': True,
        'ds_local_mesh_factor_': 1.0,
        'ds_lmax_': 5000.0,
        'set_ds_automatically_': True,
        'dsmult_': 1.0,
        'rffile_': ' ',
        'tffile_': ' ',
        'pffile_': ' ',
        'effile_': ' ',
        'kffile_': ' ',
        'qffile_': ' ',
        'lffile_': ' ',
        'rbfile_': ' ',
        'tbfile_': ' ',
        'pbfile_': ' ',
        'ebfile_': ' ',
        'kbfile_': ' ',
        'qbfile_': ' ',
        'lbfile_': ' ',
        'new_r_mesh_': True,
        'mesh_file_r_': ' ',
        'nrss_': 1,
        'r0_': 0.0,
        'r1_': 0.0,
        'new_t_mesh_': True,
        'mesh_file_t_': ' ',
        'ntss_': 2,
        't0_': 0.5,
        't1_': 0.6,
        'new_p_mesh_': True,
        'mesh_file_p_': ' ',
        'npss_': 2,
        'p0_': 0.1,
        'p1_': 0.2,
        'volume3d_output_file_r_': ' ',
        'volume3d_output_file_t_': ' ',
        'volume3d_output_file_p_': ' ',
        'slice_coords_are_xyz_': False,
        'trace_slice_direction_is_along_b_': True,
        'compute_q_on_slice_': False,
        'q_increment_h_': 0.0001,
        'slice_input_file_r_': 'lp_r.hdf',
        'slice_input_file_t_': 'lp_t.hdf',
        'slice_input_file_p_': 'lp_p.hdf',
        'trace_from_slice_forward_': True,
        'slice_output_file_forward_r_': 'rf.hdf',
        'slice_output_file_forward_t_': 'tf.hdf',
        'slice_output_file_forward_p_': 'pf.hdf',
        'trace_from_slice_backward_': True,
        'slice_output_file_backward_r_': 'rb.hdf',
        'slice_output_file_backward_t_': 'tb.hdf',
        'slice_output_file_backward_p_': 'pb.hdf',
        'slice_q_output_file_': 'q_output.hdf',
        'slice_length_output_file_': ' ',
        'ch_map_r_': 1.0,
        'ch_map_output_file_': 'ch_output.hdf',
        'compute_ch_map_3d_': False,
        'ch_map_3d_output_file_': ' ',
        'write_traces_to_hdf_': False,
        'write_traces_root_': 'fl',
        'write_traces_as_xyz_': True
    })

    def __init__(self, mapflpy_fortran, br, bt, bp, lp=None, **mapfl_params):
        """Constructor for Tracer object given launch points and B files.

        This method calls Tracer, a lightweight wrapper which calls mapfl. This
        requires the user to have the mapfl.so compiled and accessible (see the
        References section below for more details). The ``launch_pts``
        parameter is an ``ndarray``. that must have a 3xN  shape, where the
        columns are: the launch points in Heliographic coordinates (Rsun, rad, rad).

        Parameters
        ----------
        mapflpy_fortran: Python extension module
            The currently imported mapflpy_fortran shared object module.
        br_path : str | pathlib.Path
            Path to Br hdf4 or hdf5 file.
        bt_path : str | pathlib.Path
            Path to Bt hdf4 or hdf5 file.
        bp_path : str | pathlib.Path
            Path to Bt hdf4 or hdf5 file.
        launch_pts : numpy.ndarray, optional
            A 3xN array of seed point locations for the fieldline traces in Heliographic coordinates.
            These should be floating point values ordered in radius, theta, phi (r,t,p) where
            - r is in units of solar radii (1 is the unit-sphere)
            - t is co-latitude in radians (0-pi, 0 is the north pole, pi is the south pole).
            - p is longitude in radians (0-2pi).
            If not provided, `utils.fetch_default_launch_points` is called to generate them on a sphere.
        **mapfl_params : optional
            Keyword arguments are passes to the ``_load_params()`` function, which
            overrides the corresponding key-value pair in the MAPFL_PARAMS dict.
        """

        self._bfiles = BFiles(br, bt, bp)
        self._mapfl = mapflpy_fortran.mapfl
        self._mapfl_params = ChainMap(mapfl_params | FIELD_DEFAULTS, self._default_params)
        self._launch_points = lp

    def __getitem__(self, item):
        return self._mapfl_params[item]

    def __setitem__(self, key, value):
        self._mapfl_params[key] = value

    @property
    def mapfl(self):
        return self._mapfl

    @mapfl.setter
    def mapfl(self, mapflpy_fortran):
        self._mapfl = mapflpy_fortran.mapfl

    @property
    def bfiles(self):
        return self._bfiles

    @bfiles.setter
    def bfiles(self, paths: Tuple[str, str, str]):
        self._bfiles = BFiles(*paths)
        self._mapfl_params |= FIELD_DEFAULTS

    @property
    def params(self):
        return self._mapfl_params

    @property
    def launch_points(self):
        return self._launch_points

    @launch_points.setter
    def launch_points(self, lp: Iterable):
        self._launch_points = np.array(lp)

    def get_params(self):
        return self._mapfl_params

    def set_params(self, params: Mapping):
        self._mapfl_params.update(**params)

    def get_bfiles(self):
        return self._bfiles

    def set_bfiles(self, paths: Tuple[str, str, str]):
        self._bfiles = BFiles(*paths)
        self._mapfl_params |= FIELD_DEFAULTS

    def get_launch_points(self):
        return self._launch_points

    def set_launch_points(self, lp: Iterable):
        self._launch_points = np.array(lp)

    def get_mapflpy_fortran_id(self):
        return id(self._mapfl)

    def load_default_lps(self, n=128):
        self._launch_points = fetch_default_launch_points(n)

    def load_fields(self):
        for field in ('br', 'bt', 'bp'):
            result = read_hdf_by_value(ifile=getattr(self._bfiles, field))
            self._mapfl_params[field] = result[0].astype('float64').T
            for i, dim in enumerate(('r', 't', 'p')):
                self._mapfl_params[f"{field}_{dim}"] = result[i+1].astype('float64')
                self._mapfl_params[f"{field}_n{dim}"] = len(result[i+1])

    def prep_tracing(self, direction: DIRECTION = 'f'):
        assert all(os.path.exists(i) for i in self._bfiles), "One or more b-files cannot be found"
        if self._launch_points is None:
            self._launch_points = fetch_default_launch_points(128)
        if not isinstance(self._launch_points, np.ndarray):
            self._launch_points = np.array(self._launch_points)
        match len(self._launch_points.shape):
            case 1:
                self._launch_points = self._launch_points.reshape((3, 1), order='F')
            case 2:
                if self._launch_points.shape[0] != 3:
                    self._launch_points = self._launch_points.reshape((3, self._launch_points.shape[0]), order='F')
            case _:
                raise ValueError("Launch points must be an array of r, t, p values")
        if any((self._mapfl_params[k] is None for k in FIELD_DEFAULTS.keys())):
            self.load_fields()

        self["trace_fwd_"] = direction == 'f'
        self["trace_bwd_"] = direction == 'b'
        self.run()

    def run_tracing(self, buffer_size=DEFAULT_BUFFER_SIZE):
        traces = np.full((buffer_size, *self._launch_points.shape), np.nan, order='F').astype(np.float64)
        s1 = np.zeros(self._launch_points.shape, np.float64, order='F')
        mask = np.full((1, self._launch_points.shape[1]), False, order='F')
        for i in range(self._launch_points.shape[1]):
            trace_args = dict(
                s0=self._launch_points[:,i],
                s1=s1[:,i],
                bs0=np.zeros(3, order='F').astype(np.float64),
                bs1=np.zeros(3, order='F').astype(np.float64),
                s=np.zeros(1, order='F').astype(np.float64),
                traced_to_r_boundary=mask[:,i],
                svec=traces[:, :, i],
                svec_n=buffer_size
            )
            self.trace(trace_args)

        # return the traces without the first point since that is always a duplicated in the trace
        # NOTE: even if the trace has only one point, this will return an array with the same dimensionality/order
        return Traces(traces[:, :, :], s1, mask[0, :])

    def run(self):
        """Calls ``mapfl.run()`` with provisioned parameters."""
        self._mapfl.run(**self._mapfl_params)

    def init(self, init_params):
        """Calls ``mapfl.init()`` with provisioned init parameters.

        *NOTE: In the current state of mapflpy_fortran, this routine is not fully implemented.
        In fact, it does nothing. This method is provided with future development
        of mapflpy_fortran in mind. For now, the only way to modify the parameters used
        by mapfl is to call* ``run()``.
        """
        self._mapfl.init(**init_params)

    def trace(self, trace_args):
        """Calls ``mapfl.trace()`` with provisioned trace arguments."""
        self._mapfl.trace(**trace_args)

