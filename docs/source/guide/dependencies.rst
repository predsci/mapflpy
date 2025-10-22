Dependencies
============

Required Dependencies
----------------------
- `Python >= 3.11 <https://www.python.org/>`_
- `numpy <https://numpy.org/>`_
- `h5py <https://www.h5py.org/>`_
- `gfortran <https://gcc.gnu.org/fortran/>`_
- `psi-io <https://pypi.org/project/psi-io/>`_

Optional Dependencies
----------------------
- `pyhdf <https://github.com/fhs/pyhdf>`_
- `matplotlib <https://matplotlib.org/>`_

Sample Conda Environment
------------------------

.. code-block:: yaml

    name: mapflpy-env
    channels:
      - conda-forge
    dependencies:
      - python>=3.11
      - numpy
      - h5py
      - gfortran
      - pyhdf
      - matplotlib
      - pip
      - pip:
        - psi-io
        - mapflpy @ git+ssh://git@github.com/predsci/mapflpy.git

This environment can be created using the following command:
.. code-block:: bash

    conda env create -f environment.yml
