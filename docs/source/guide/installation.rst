.. _installation:

Installation
============

.. attention::

    ``mapflpy`` is not yet available on `PyPI <https://pypi.org/>`_. For now, you can install it from Predictive Science
    Inc.'s `GitHub <https://github.com/predsci/>`_ repository.

    We highly recommend using a virtual environment to manage your Python packages and avoid conflicts with other
    projects. For the best results, we recommend using ``conda`` – *via* Miniforge (preferred), Miniconda, or Anaconda
    – to create and manage your virtual environments.

To get started with **mapflpy**, you can pip install it directly from the GitHub repository using the following commands
(depending on your preferred method of access):

.. code-block:: bash

    pip install git+ssh://git@github.com/predsci/mapflpy.git

or

.. code-block:: bash

    pip install git+https://github.com/predsci/mapflpy.git

**Ensure that you have** ``git`` **installed and configured on your system to use either of the above methods.**

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
- `pooch <https://www.fatiando.org/pooch/>`_

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
      - pooch
      - pip
      - pip:
        - psi-io
        - mapflpy @ git+ssh://git@github.com/predsci/mapflpy.git

This environment can be written to a file, named `environment.yml`, and created using the following command:

.. code-block:: bash

    conda env create -f environment.yml