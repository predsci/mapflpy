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

## Installation

---

```bash
pip install mapflpy
```

**Note**: Because mapflpy includes a Python extension built from cross-compiled
Fortran code, the pre-built wheels are platform specific. If necessary, please 
see the instructions for building mapflpy from source. 

## Usage

---

### Basic Tracing
```python
# Code goes here
```

## Requirements

---

This package requires the python HDF5 interface, `h5py`, to work with `.h5` files. 

If you are working with HDF4 `.hdf` files then you must also have the optional
`pyhdf` HDF4 python interface installed.

Because these packages require the underlying C HDF libraries to be installed, 
we generally recommend using conda to install these dependencies. This often
makes life much easier than installing from source:

for HDF5
```bash
conda install h5py
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
their system.

Automated tests and basic documentation are coming soon!
