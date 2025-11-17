.. raw:: html

   <p style="display: flex; align-items: center; justify-content: center; gap: 0.5rem;">
     <a href="https://predsci.com">
       <img
         src="https://www.predsci.com/corona/apr2024eclipse/images/psi_logo.png"
         alt="Predictive Science Inc. logo"
         width="50"
       />
     </a>
     <span style="font-size: 2em; font-weight: bold;">mapflpy</span>
   </p>

.. image:: https://img.shields.io/pypi/v/mapflpy?logo=pypi&logoColor=white
   :target: https://pypi.org/project/mapflpy
   :align: center

.. image:: https://img.shields.io/pypi/pyversions/mapflpy.svg?logo=python&label=python&logoColor=white
   :target: https://pypi.org/project/mapflpy
   :align: center

.. image:: https://img.shields.io/pypi/l/mapflpy?logo=apache&logoColor=white
   :target: https://opensource.org/license/apache-2-0/
   :align: center

.. image:: https://img.shields.io/librariesio/github/predsci/mapflpy?logo=Libraries.io&logoColor=white
   :target: https://github.com/predsci/mapflpy/blob/main/pyproject.toml
   :alt: Libraries.io dependency status for GitHub repo
   :align: center

.. image:: https://github.com/predsci/mapflpy/actions/workflows/publish.yml/badge.svg
   :target: https://github.com/predsci/mapflpy/actions/workflows/publish.yml
   :align: center

.. image:: https://github.com/predsci/mapflpy/actions/workflows/docs.yml/badge.svg
   :target: https://predsci.com/doc/mapflpy
   :align: center

----

Python extension for tracing field lines using the
`MapFL <https://github.com/predsci/MapFL>`_ Fortran tracer.

The goal of **mapflpy** is to provide fast and accurate tracing capabilities for
spherical vector fields inside a convenient Python interface. **mapflpy** is
designed to work natively with the staggered meshes produced by
Predictive Science Inc.'s codes for simulating the solar corona, and 
inner heliosphere (*e.g.* `MAS <https://www.predsci.com/mas>`_ or
`POT3D <https://github.com/predsci/POT3D>`_), but it should be generally compatible
with any global vector field that can be described on a rectilinear grid in 
spherical coordinates.

Resources
---------

- **Predictive Science Inc.:** https://predsci.com
- **Source Code:** https://github.com/predsci/mapflpy
- **Documentation:** https://predsci.com/doc/mapflpy
- **Distribution:** https://pypi.org/project/mapflpy

