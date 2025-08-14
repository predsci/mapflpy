from collections import namedtuple
from enum import IntEnum
from os import PathLike
from pathlib import Path
from types import MappingProxyType
from typing import Literal, Tuple

from numpy import float64, float32
from numpy._typing import NDArray


class Polarity(IntEnum):
    """
    Enumeration of magnetic fieldline polarity classifications.

    This enumeration is used to label the connectivity of magnetic
    fieldlines based on their endpoints. The values encode both the
    type of connectivity (open, closed, or invalid) and, for open
    lines, the sign of the radial magnetic field at the footpoint.

    Members
    -------
    R0_R1_NEG : int (-2)
        Open fieldline connecting from the inner boundary (R0)
        to the outer boundary (R1) with a **negative** radial
        magnetic field (Br < 0) at the inner boundary footpoint.

    R0_R0 : int (-1)
        Closed fieldline with both endpoints anchored at the
        inner boundary (R0).

    ERROR : int (0)
        Indicates an undefined or unclassified trace, typically used
        when the endpoints do not terminate at a boundary.

    R1_R1 : int (1)
        Closed fieldline with both endpoints anchored at the
        outer boundary (R1) *e.g.* "Disconnected" trace.

    R0_R1_POS : int (2)
        Open fieldline connecting from the inner boundary (R0)
        to the outer boundary (R1) with a **positive** radial
        magnetic field (Br > 0) at the inner boundary footpoint.

    Notes
    -----
    The sign convention assumes that the radial magnetic field (Br)
    is evaluated at the **inner boundary** footpoint of the fieldline.
    """
    R0_R1_NEG = -2
    R0_R0 = -1
    ERROR = 0
    R1_R1 = 1
    R0_R1_POS = 2


# ------------------------------------------------------------------------------
# Named tuple for storing trace information.
# This is used to return the results of tracing operations.
# It contains the geometry of the traces, their start and end positions,
# and whether they were traced to a boundary.
#------------------------------------------------------------------------------
Traces = namedtuple('Traces', ['geometry', 'start_pos', 'end_pos', 'traced_to_boundary'])

DirectionType = Literal['f', 'b']
MagneticFieldLabelType = Literal['br', 'bt', 'bp']
PathType = str | Path | PathLike[str]
ArrayType = NDArray[float64 | float32]
MagnetifFieldArrayType = Tuple[ArrayType, ArrayType, ArrayType, ArrayType] | PathType

DEFAULT_BUFFER_SIZE = 2000                  # Default buffer size for reading HDF files
MAGNETIC_FIELD_LABEL = ('br', 'bt', 'bp')   # Magnetic field labels
DIRECTION = ('f', 'b')                      # Direction of tracing (forward or backward)


# ------------------------------------------------------------------------------
# Base dictionary for maintaining magnetic field filepaths.
# Used to pass field data within the TracerMP class.
# ------------------------------------------------------------------------------
MAGNETIC_FIELD_PATHS = MappingProxyType({
    'br': '',
    'bt': '',
    'bp': '',
})


# ------------------------------------------------------------------------------
# Base dictionary for maintaining magnetic field array/scale data.
# Used to pass field data within the Tracer class.
# ------------------------------------------------------------------------------
DEFAULT_FIELDS = MappingProxyType({
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

# ------------------------------------------------------------------------------
# The following MAPFL_PARAMS are modeled after the mapfl.in file. However, the
# following key-value pairs are not exhaustive; rather, they are a *working*
# configuration that will suit most general uses of the _Tracer class.
# ------------------------------------------------------------------------------
DEFAULT_PARAMS = MappingProxyType({
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
    'write_traces_as_xyz_': True,
})
