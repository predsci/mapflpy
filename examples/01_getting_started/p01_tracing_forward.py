"""
Trace Forward
=============

Perform forward tracing of magnetic field lines.

This example demonstrates how to use the :func:`~mapflpy.scripts.run_foward_tracing`
function to trace magnetic field lines forward from a set of default starting points.
It also shows how to load magnetic field data files and visualize the traced field lines in 3D.
"""
import os
if 'SPHINX_GALLERY_BUILD' not in os.environ:
    import matplotlib
    matplotlib.use('TkAgg')

from mapflpy.scripts import run_foward_tracing
from mapflpy.utils import plot_traces
from mapflpy.data import fetch_cor_magfiles
import matplotlib.pyplot as plt

# %%
# Load in the magnetic field files
magnetic_field_files = fetch_cor_magfiles()

# %%
# The :func:`~mapflpy.utils._fetch_magnetic_field_files` function returns a tuple of file paths
# to the example magnetic field data files for demonstration purposes. These file paths can be
# ammended to point to user-specific magnetic field data files.
for filepath in magnetic_field_files:
    print(f"Magnetic field file: {filepath}")

# %%
# Run tracing and print a summary
traces = run_foward_tracing(*magnetic_field_files, context='fork')
print("geometry shape:", traces.geometry.shape)  # stdout should display

# %%
# Plot
ax = plt.figure().add_subplot(projection='3d')
plot_traces(traces, ax=ax)
plt.show()