# mapflpy :: IN DEVELOPMENT
***!!! BETA: Please consult Cooper or Ryder before using package. !!!***

---

---

Python extension for tracing field lines using the Fortran tracer from 
[MapFL](https://github.com/predsci/MapFL).

The goal of mapflpy is to provide fast and accurate tracing capabilities for 
spherical vector fields inside a convenient Python interface.

`mapflpy` is designed to work natively with the staggered meshes produced by 
Predictive Science Inc.'s codes for simulating the solar corona, and 
inner heliosphere (e.g. [MAS](https://www.predsci.com/mas) or 
[POT3D](https://github.com/predsci/POT3D)), but it should be generally compatible
with any global vector field that can be described on a rectilinear grid in 
spherical coordinates.

Cooper can provide a wheel for those interested in using this package while it
is in production.

### Disclaimer

This package is currently in a pre-release state as we make an attempt to 
simplify the process for end-users to get working mapflpy libraries on
their system. The routines, organization, and interfaces are subject to change.

Automated builds, basic documentation, and an initial release are coming soon!
Thank you for your patience.

---

---

## I. Installation

---

Install from pre-built wheels (eventually...)
```bash
# pip install mapflpy  # this will work when mapflpy is on PyPi
```

Install from source:
- Right now the alternative of installing from source with pip isn't perfect yet 
(meson linking vagaries).

**Note**: Because mapflpy includes a Python extension built from cross-compiled
Fortran code, the pre-built wheels are platform specific. If necessary, please 
see the instructions for building mapflpy from source.


### A. Building
If you really want to build a wheel yourself, it *should* work with the 
`mapflpy-dev` environment prescribed in `environment.yml`. Have conda set that
up and then from within that conda environment type:
```bash
python -m build
```
This will build a `.whl` wheel file specific to your OS and CPU architecture. 
Then install the wheel, e.g.:
```bash
pip install <WHEEL_FILE>
```
Unless you've used a "repair tool" to consolidate the shared libraries to the
wheel (e.g. 
[delocate](https://github.com/matthew-brett/delocate) or 
[auditwheel](https://github.com/pypa/auditwheel))
then this might only work *within* the conda environment that built the package.

I don't recommend playing with these build tools manually to build for many
systems, and there are several *idiosyncrasies* with the current setup. 
We will likely automate this process using github actions (i.e. talk to Cooper).

### B. Testing
The automated tests will confirm if tracing is working properly or not. Here
we use pytest. Run this after installing mapflpy
```bash
pytest -rA -v --pyargs mapflpy.tests
```

## II. Usage

---
`mapflpy` centers around the `Tracer` and `TracerMP` classes (both of which inherit
from the base `_Tracer` class interface). 
The `Tracer` class is for single-threaded tracing, while the `TracerMP` class is for 
multi-processed tracing using python's`multiprocessing` module.

The `mapflpy_fortran` shared-object is built from `mapfl` Fortran source code using
F2PY and Meson. 
The Fortran routines are not intended to be called directly by end-users, but rather 
through the `Tracer` and `TracerMP` classes.
`mapfl` itself is not designed for object-oriented programming; rather, it is generally 
called from the command line using a `.in` file to specify input arguments, and relies 
on global variables to hold state information.

This package wraps the Fortran routines in a "pythonic" interface, allowing for a
more flexible approach to tracing fieldlines. 
Until there is time (...) to refactor the Fortran code into a more modular design, 
these wrapper classes are the most feasible way to allow a broader audience to use 
the Fortran tracer.

### A. Tracer vs. TracerMP

The `Tracer` class directly imports the `mapflpy_fortran` cross-compiled Fortran
module.
Since imports in Python are singletons, only one instance of the Fortran module can
exist within a process at a time.
As such, the `Tracer` class also enforces a singleton pattern (due to the fact that 
the Fortran routines rely on global state variables, and multiple instances of the 
`Tracer` class would lead to conflicts and race-conditions).

The `TracerMP` class, on the other hand, spawns multiple processes using the
`multiprocessing` module. 
Each instance of the `TracerMP` class creates and links to a distinct subprocess, 
which imports its own instance of the Fortran module. 
A two-way pipe is used to communicate between the main process and the subprocess.
*This allows multiple instances of the Fortran routines to be used simultaneously,
without conflicts (at the cost of inter-process communication overhead).*
In an attempt to limit this overhead, magnetic field data is passed as filepaths.
These filepaths are passed to the respective subprocess over the pipe, then
loaded into memory by the subprocess itself.

### B. Tracer as Dictionary

The base `_Tracer` class implements the `MutableMapping` interface, allowing for
`Tracer` and `TracerMP` instances to behave like dictionaries. Under the hood,
the "mapfl.in" parameters are stored in a `ChainMap` which contains a set of default
parameters, as well as any user-specified changes (or additions) to the parameters.

```python
from mapflpy.tracer import Tracer

tracer = Tracer()
current_params = dict(tracer)                 # get current parameters as a dictionary

tracer['verbose_'] = True                     # set verbose_ parameter to True
tracer.update(ds_min_=0.00001, ds_max_=10.1)  # update multiple parameters at once
updated_params = dict(tracer)                 # get updated parameters as a dictionary

tracer.clear()                                # reset parameters to defaults
```

The one exception to this is the magnetic field data itself. Because this data is
handled differently in the `Tracer` and `TracerMP` classes, the magnetic field data
should be set using the `br`, `bt`, and `bp` properties (or the `set_field_data` method).
`TracerMP` must be set with filepaths, whereas `Tracer` can be set with numpy arrays **or**
filepaths *i.e.*

```python
from mapflpy.tracer import Tracer
from psi_io import read_hdf_by_value

# load magnetic field data from HDF files
# read_hdf_by_value returns the data array followed by any scale arrays
# e.g. values, r_scale, t_scale, p_scale
br, *br_scales = read_hdf_by_value(ifile="br_file.h5")
bt, *bt_scales = read_hdf_by_value(ifile="bt_file.h5")
bp = "bp_file.h5"

tracer = Tracer()
tracer.br = br, *br_scales 
tracer.bt = bt, *bt_scales
tracer.bp = bp  # can be a filepath
```

or

```python
from mapflpy.tracer import TracerMP

tracer_mp = TracerMP()
tracer_mp.br = "br_file.h5"
tracer_mp.bt = "bt_file.h5"
tracer_mp.bp = "bp_file.h5"
```



### C. Tracing Fieldines
Once a `Tracer` or `TracerMP` instance has been created, and the magnetic field
data has been set, fieldlines can be traced using the `trace` method. 
***NOTE: prior to tracing fieldlines the `run` method must be called viz. to
populate any changes made to the input params or magnetic field data.***

With that said, the current "staleness" of the tracer instance can be checked using
the `stale` property *i.e.* whether there have been any changes to the input parameters
or magnetic field data since the last time the `run` method was called.

```python
from mapflpy.tracer import Tracer

tracer = Tracer()
print(tracer.stale)         # True, since nothing has been run yet
tracer.run()                # run the tracer to initialize
print(tracer.stale)         # False, since tracer is up-to-date
tracer['ds_min_'] = 0.0001  # change a parameter
print(tracer.stale)         # True, since a parameter has changed
```

NOTE: If `trace` is called while the tracer is stale, `run` will be called automatically.

```python
from mapflpy.tracer import Tracer
import numpy as np

lps = [
    np.full(10, 1.01),              # r launch points [R_sun]
    np.linspace(0, np.pi, 10),      # theta launch points [rad]
    np.zeros(10)                    # phi launch points [rad]
]

tracer = Tracer()
tracer.load_fields( ... )            # load magnetic field data
tracer.trace(lps, buffer_size=1000)  # will call run() if stale
```

Traces can be performed "forward" or "backward" by calling the `set_trace_direction`
method with either `'f'` or `'b'`. 
*Two separate calls must be made to trace in
both directions.* 
The resulting trace geometry must then be combined manually.
With that said, a utility function `combine_fwd_bwd_traces` is provided to help with 
this common use case.

```python
from mapflpy.tracer import Tracer
from mapflpy.utils import combine_fwd_bwd_traces

tracer = Tracer()
lps = [ ... ]                           # launch points
tracer.load_fields( ... )               # load magnetic field data
tracer.set_tracing_direction('f')       # set to forward tracing
fwd_traces = tracer.trace(lps)          
tracer.set_tracing_direction('b')       # set to backward tracing
bwd_traces = tracer.trace(lps)

# combine the two traces into one, correctly ordered
combined = combine_fwd_bwd_traces(fwd_traces, bwd_traces)
```

Example of iterating through a series of states within a time-dependent run:

```python
from mapflpy.tracer import Tracer
import numpy as np

traces =[]
states = range(10, 20)

lps = [
    np.full(10, 1.01),              # r launch points [R_sun]
    np.linspace(0, np.pi, 10),      # theta launch points [rad]
    np.zeros(10)                    # phi launch points [rad]
]

tracer = Tracer()
for state in states:
    tracer.load_fields(
        br=f"br_0000{state}.h5",
        bt=f"bt_0000{state}.h5",
        bp=f"bp_0000{state}.h5"
    )            
    traces.append(tracer.trace(lps, buffer_size=1000))  # will call run() if stale
```

The result of `trace` is a `Traces` object – a named-tuple-like container for the
traced fieldlines. This structure contains the following attributes:
- `geometry` : an Nx3xM array of the traced fieldline coordinates, *i.e.* the 
radial-theta-phi coordinates of the traces where N is the buffer size, and M is 
the number of fieldlines. 
(Note, to preserve a homogeneous array, fieldlines shorter than N are NaN padded).
- `start_pos` : an Mx3 array of the starting positions of each fieldline.
- `end_pos` : an Mx3 array of the ending positions of each fieldline.
- `traced_to_boundary` : a boolean array of length M indicating whether each fieldline
traced to a boundary (True) or was terminated early due to step-size constraints
(False).


### D. Using Scripts

Several scripts are provided in the `scripts/` directory.
These standalone functions are used to perform common "one-off" tracing tasks
(similar, in many respects, to how `mapfl` itself is used from the command line).

Any additional keyword arguments (see function signature) provided within the calls 
to these scripts are passed to the instantiation of the `TracerMP` class, *i.e.* 
arguments used to set the mapfl parameters.

```python
from mapflpy.scripts import run_forward_tracing

lps = [...]
bfiles = {
    'br': 'br_file.h5',
    'bt': 'bt_file.h5',
    'bp': 'bp_file.h5'
}
traces = run_forward_tracing(
    **bfiles,
    launch_points=lps,  # <-- passed to trace() method
    buffer_size=1000,  # <-- passed to trace() method
    domain_r_max_=100  # <-- example of passing mapfl parameter
)
```

or

```python
from mapflpy.scripts import run_fwdbwd_tracing

lps = [ ... ]
bfiles = {
    'br': 'br_file.h5',
    'bt': 'bt_file.h5',
    'bp': 'bp_file.h5'
}
traces = run_fwdbwd_tracing(
    **bfiles,
    launch_points=lps,  # <-- passed to trace() method
    buffer_size=1000,   # <-- passed to trace() method
)
```

This module also contains a script for performing interdomain tracing between
two different magnetic field domains (e.g. coronal to heliospheric). A more 
comprehensive explanation of the function signature can be found in the docstring for
the `run_interdomain_tracing` function. 

A few general notes about this function:
- requires 6 magnetic field files: 3 for the inner domain, and 3 for the outer domain.
- the `r_interface` parameter must be specified to indicate the radial location of the
recross boundary between the two domains.
- the `helio_shift` parameter can be used to account for the longitudinal shift angle 
between the heliospheric domain and the coronal domain.

```python
from mapflpy.scripts import inter_domain_tracing
from math import pi
lps = [ ... ]
coronal_bfiles = {
    'br_cor': 'br_coronal.h5',
    'bt_cor': 'bt_coronal.h5',
    'bp_cor': 'bp_coronal.h5'
}
heliospheric_bfiles = {
    'br_hel': 'br_helio.h5',
    'bt_hel': 'bt_helio.h5',
    'bp_hel': 'bp_helio.h5'
}
traces = inter_domain_tracing(
    **coronal_bfiles,
    **heliospheric_bfiles,
    launch_points=lps,      
    r_interface=30.0,         # <-- radial location of domain interface [R_sun]
    helio_shift=pi/6,         # <-- longitudinal shift between domains [rad]
    rtol_=1e-6,               # <-- relative tolerance used to determine when
                              #     a fieldline has crossed the interface boundary
)
```

### E. Utilities

A few utility functions are provided in the `mapflpy.utils` module.
`get_fieldline_endpoints`, `get_fieldline_npoints`, and `trim_fieldline_nan_buffer`
can be used on raw fieldline geometry (numpy arrays) or on `Traces` objects to extract
information about the traced fieldlines, or (in the case of `trim_fieldline_nan_buffer`)
to remove the NaN padding from fieldlines – returning a list of heterogeneously sized
fieldlines.

Lastly, `get_fieldline_polarity` can be used to determine and classify fieldlines as: 
- "open" (with negative polarity)
- "closed"
- "undefined" (e.g. if the fieldline did not trace to a boundary)
- "disconnected"
- "open" (with positive polarity)

The result of this function is an array of type `Polarity` (an `IntEnum` class) which
provides a mapping between integer codes and string labels for each fieldline.

## III. Requirements

---

This package requires the python HDF5 interface, `h5py`, to work with `.h5` files.
It also requires the `psi-io` python package, which has reader/writers for 
PSI-style HDF data files.

If you are working with HDF4 `.hdf` files then you must also have the optional
`pyhdf` HDF4 python interface installed.

Because these packages require the underlying C HDF libraries to be installed, 
we generally recommend using conda to install these dependencies. This often
makes life much easier than installing from source:

for HDF5 and psi-io
```bash
conda install h5py
pip install psi-io
```


for HDF4
```bash
conda install pyhdf
```

To isolate these to a specific environment, see `environment.yml`

