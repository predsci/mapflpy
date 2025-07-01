# mapflpy

---

Python extension for tracing field lines using the Fortran tracer from 
[MapFL](https://github.com/predsci/MapFL).

The goal of mapflpy is to provide fast and accurate tracing capabilities for 
spherical vector fields inside a convenient Python interface.

mapflpy is designed to work natively with the staggered meshes produced by 
Predictive Science Inc.'s codes for simulating the solar corona, and 
inner heliosphere (e.g. [MAS](https://www.predsci.com/mas) or 
[POT3D](https://github.com/predsci/POT3D)), but it should be generally compatible
with any global vector field that can be described on a rectilinear grid in 
spherical coordinates.

# NOT READY YET!!!
I wanted to put the current state of the mapflpy package up on github to test
building platform dependent wheels using github actions.

Cooper can provide a wheel in the meantime (or see "Building" below) but the
installation and usage instructions are not done yet!

## Installation (not ready yet)

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


### Building
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

### Testing
The automated tests will confirm if tracing is working properly or not. Here
we use pytest. Run this after installing mapflpy
```bash
pytest -rA -v --pyargs mapflpy.tests
```

## Usage (not ready yet)

---
The main methods that call the tracer are in `mapflpy.run`.

### Basic Tracing
```python
# Code goes here
```

### Modifying Trace Arguments

## Requirements

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

## Disclaimer

---

This package is currently in a pre-release state as we make an attempt to 
simplify the process for end-users to get working mapflpy libraries on
their system. The routines, organization, and interfaces are subject to change!

Automated builds, basic documentation, and an initial release is coming soon!
