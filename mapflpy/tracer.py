import os
import sys
from collections import namedtuple, ChainMap
from types import MappingProxyType
from typing import Tuple, Literal, Iterable, Mapping

import numpy as np
from psipytools.geometry.points import fibonacci_sphere

from psipytools.psihdf import read_hdf_by_value

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


# @lru_cache(maxsize=32)
def fetch_default_launch_points(n):
    p, t = fibonacci_sphere(n)
    return np.array((np.full_like(t, 1.01), t, p), order='F')


class Tracer:
    """Fieldline tracing wrapper class for Mapflpy.

    Attributes
    ----------
    _mapflpy : mapflpy.so
        Path to mapfl shared object directory.
    _mapfl_params : dict
        Dictionary of key-value pairs to instantiate the mapflpy object. The
        keys used to instantiate mapflpy must have a trailing underscore.
    _launch_points : numpy.ndarray
        Seed points for the fieldline traces in Heliographic coordinates.
    _trace_rgb : numpy.ndarray
        RGB values for each field line. Must have the same shape as the
        ``_launch_points`` attribute *i.e.* Nx3.
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
        'bfile_r_': 'br002.hdf',
        'bfile_t_': 'bt002.hdf',
        'bfile_p_': 'bp002.hdf',
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

    def __init__(self, mapflpy, br, bt, bp, lp=None, **mapfl_params):
        """Constructor for Tracer object given launch points and mag vr_files.

        This method calls Tracer, a lightweight wrapper which calls mapflpy. This
        requires the user to have the mapfl.so compiled and accessible (see the
        References section below for more details). The ``launch_pts``
        parameter can be either a path to an ``ascii`` file or an ``ndarray``.
        Both the former (table) and the latter must have a Nx3 or Nx6 shape,
        where the columns are: the launch points in Heliographic coordinates
        (Rsun, rad, rad), or *both* the launch points in Heliographic coordinates
        and the field lines' RGB values (R, G, B). When no RGB values are specified,
        a random color from the HSV colorspace (converted to RGB) will be assigned.

        Parameters
        ----------
        br_path : str | pathlib.Path
            Path to Br hdf4 or hdf5 file.
        bt_path : str | pathlib.Path
            Path to Bt hdf4 or hdf5 file.
        bp_path : str | pathlib.Path
            Path to Bt hdf4 or hdf5 file.
        launch_pts : str | pathlib.Path | numpy.ndarray
            Launch points for fieldline tracing. This can be either a path to an
            ``ascii`` file or an ``ndarray``.
        mapflpy_path : str | pathlib.Path, default = None
            Path to mapfl.so directory. When no directory is supplied, the ``Tracer``
            object will search for the mapfl shared object in the user's CORHEL_HOME
            path.
        **mapfl_params : optional
            Keyword arguments are passes to the ``_load_params()`` function, which
            overrides the corresponding key-value pair in the MAPFL_PARAMS dict.
        """

        self._bfiles = BFiles(br, bt, bp)
        self._mapfl = mapflpy.mapfl
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
    def mapfl(self, mapflpy):
        self._mapfl = mapflpy.mapfl

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

    def get_mapflpy_id(self):
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
        """Calls ``mapflpy.run()`` with provisioned parameters."""
        self._mapfl.run(**self._mapfl_params)

    def init(self, init_params):
        """Calls ``mapflpy.init()`` with provisioned init parameters.

        *NOTE: In the current state of mapflpy, this routine is not fully implemented.
        In fact, it does nothing. This method is provided with future development
        of mapflpy in mind. For now, the only way to modify the parameters used
        by mapfl is to call* ``run()``.
        """
        self._mapfl.init(**init_params)

    def trace(self, trace_args):
        """Calls ``mapflpy.trace()`` with provisioned trace arguments."""
        self._mapfl.trace(**trace_args)

    # def get_lps_rgb(self):
    #     """Getter for ``self._trace_rgb``"""
    #     return self._trace_rgb
    #
    # def fwd_trace(self):
    #     """Perform forward traces from each launch point with mapfl.
    #
    #     Returns
    #     -------
    #     tuple
    #         - list[numpy.ndarray]
    #             List of fieldline traces returned from ``_mapflpy.mapfl.trace()``.
    #             This list is of length N, where N is the number of launch points.
    #             Each element of this list is an Mx3 ndarray of the points making up
    #             the given field line (in Heliographic coordinates) starting from
    #             a launch point and tracing forward.
    #         - list[bool]
    #             Mask of length N, where N is the number of launch points, indicating
    #             which of the returned field lines do not map back to the Sun's
    #             surface, or rather ``self._mapfl_params["domain_r_min_"]``.
    #     """
    #     self._mapfl_params["trace_fwd_"] = True
    #     self._mapfl_params["trace_bwd_"] = False
    #     self.run()
    #     traces_fwd = []
    #     traces_fwd_map = np.zeros((len(self._launch_points))).astype(bool)
    #     for i, lp in enumerate(self._launch_points):
    #         trace_args = self._generate_trace_args(lp)
    #         self.trace(trace_args)
    #         svec_cleaned = trace_args['svec'][~np.isnan(trace_args['svec']).any(axis=1)]
    #         traces_fwd_map[i] = svec_cleaned[-1][0] < self._mapfl_params["domain_r_min_"] + self._mapfl_params[
    #             "ds_min_"]
    #         traces_fwd.append(svec_cleaned)
    #
    #     return traces_fwd, traces_fwd_map
    #
    # def bwd_trace(self, flip=True):
    #     """Perform backward traces from each launch point with mapfl.
    #
    #     Parameters
    #     ----------
    #     flip : bool, default = True
    #         Reverse the ordering of the points that make up each fieldline *i.e.*
    #         instead of the trace beginning at each respective launch point, the
    #         launch point will be the final point of the ndarray.
    #
    #     Returns
    #     -------
    #     tuple
    #         - list[numpy.ndarray]
    #             List of fieldline traces returned from ``_mapflpy.mapfl.trace()``.
    #             This list is of length N, where N is the number of launch points.
    #             Each element of this list is an Mx3 ndarray of the points making up
    #             the given field line (in Heliographic coordinates) starting from
    #             a launch point and tracing backward.
    #         - list[bool]
    #             Mask of length N, where N is the number of launch points, indicating
    #             which of the returned field lines do not map back to the Sun's
    #             surface, or rather ``self._mapfl_params["domain_r_min_"]``.
    #     """
    #     self._mapfl_params["trace_fwd_"] = False
    #     self._mapfl_params["trace_bwd_"] = True
    #     self.run()
    #     traces_bwd = []
    #     traces_bwd_map = np.zeros((len(self._launch_points))).astype(bool)
    #     for i, lp in enumerate(self._launch_points):
    #         trace_args = self._generate_trace_args(lp)
    #         self.trace(trace_args)
    #         svec_cleaned = trace_args['svec'][~np.isnan(trace_args['svec']).any(axis=1)]
    #         traces_bwd_map[i] = svec_cleaned[-1][0] < self._mapfl_params["domain_r_min_"] + self._mapfl_params[
    #             "ds_min_"]
    #         if flip:
    #             traces_bwd.append(np.flip(svec_cleaned, axis=0))
    #         else:
    #             traces_bwd.append(svec_cleaned)
    #     return traces_bwd, traces_bwd_map
    #
    # def fwd_bwd_trace(self):
    #     """Perform forward and backward traces from each launch point with mapfl.
    #
    #     This function independently calls ``fwd_trace()`` then ``bwd_trace()``
    #     and concatenates the (flipped) backward traces with the forward traces.
    #     The final point of the backward trace is clipped because the first point
    #     of the forward trace and the final point of the backward trace are the
    #     same *viz.* the original launch point.
    #
    #     Just as ``fwd_trace()`` and ``bwd_trace()`` return a mask indicating
    #     which field lines trace back to some "domain_r_min_", ``fwd_bwd_trace()``
    #     conjunctively combines the forward and backward mapping masks.
    #
    #     Returns
    #     -------
    #     tuple
    #         list[numpy.ndarray]
    #             List of fieldline traces returned from ``_mapflpy.mapfl.trace()``.
    #             This list is of length N, where N is the number of launch points.
    #             Each element of this list is an (F+B-1)x3 ndarray of the points
    #             making up field line (in Heliographic coordinates) where F is the
    #             number of coordinates returned from the ``fwd_trace()`` call and
    #             B is the number of points returned from the ``bwd_trace()`` call.
    #         list[bool]
    #             Mask of length N, where N is the number of launch points, indicating
    #             which of the returned field lines do not map back to the Sun's
    #             surface, or rather ``self._mapfl_params["domain_r_min_"]``.
    #     """
    #     fwd, fwd_map = self.fwd_trace()
    #     bwd, bwd_map = self.bwd_trace()
    #     mask = fwd_map & bwd_map
    #     full_traces = [np.concatenate((bwd[i][:-1], fwd[i])) for i in range(len(mask))]
    #     return full_traces, mask
    #
    # def get_rgb_colors(self):
    #     """Return an Nx3 list of RGB colors corresponding to N launch points.
    #
    #     If ``Tracer`` object was instantiated with a ``launch_pts`` parameter
    #     that contained user-defined RGB colors for each trace, this function
    #     acts as a getter for ``_trace_rgb``. If not, this function generates
    #     a list of RGB values using the psipytools ``get_random_rgb()`` function.
    #
    #     Returns
    #     -------
    #     numpy.ndarray | list[list[int]]
    #         RGB values for each launch point in `_launch_points``
    #     """
    #     return self._trace_rgb
    #     # if np.any(self._trace_rgb):
    #     #     return self._trace_rgb.astype(int)
    #     # else:
    #     #     return get_random_color_sequence(len(self._launch_points))
    #
    # def _generate_trace_args(self, lp, buffer_size=DEFAULT_BUFFER_SIZE):
    #     # These arguments are passed as keyword arguments to mapfl.trace()
    #     # This class is primarily concerned with the output of svec, the Nx3
    #     # array that holds the trace points returned from mapfl.
    #
    #     trace_args = dict(
    #         s0=np.array(lp).astype(np.float64),
    #         s1=np.zeros(3).astype(np.float64).T,
    #         bs0=np.zeros(3).astype(np.float64).T,
    #         bs1=np.zeros(3).astype(np.float64).T,
    #         s=np.zeros(1).astype(np.float64).T,
    #         traced_to_r_boundary=np.array([False]),
    #         svec=np.full((3, buffer_size), np.nan).astype(np.float64).T,
    #         svec_n=buffer_size,  # surprisingly not necessary
    #     )
    #     return trace_args
    #
    # def _load_mapflpy(self, mapflpy_path):
    #     # mapflpy_path can be explicitly passed to this function if one wishes
    #     # to define a particular version of mapflpy, or if MAPFLPY_DIR has not been
    #     # added to the user's environment.
    #     # Alternatively, one can add MAPFLPY_DIR=path/to/mapfl.so/dir to their
    #     # .bashrc or .env (as outlined in the ``examples`` directory).
    #
    #     match mapflpy_path:
    #         case None:
    #             mapflpy_path = os.getenv("MAPFLPY")
    #             assert mapflpy_path, "No MAPFLPY found in environ"
    #             sys.path.append(mapflpy_path)
    #             import mapflpy
    #         case str():
    #             sys.path.append(mapflpy_path)
    #             import mapflpy
    #         case _:
    #             mapflpy = mapflpy_path
    #
    #     self.mapfl = mapflpy.mapfl
