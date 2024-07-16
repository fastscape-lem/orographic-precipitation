Linear Theory of Orographic Precipitation
=========================================

|Build Status| |Notebooks status| |Coverage| |License|

A `Python`_ framework that implements the Linear Theory of Orographic Precipitation
following `Smith & Barstad (2004)`_.

.. |Build Status| image:: https://github.com/EstebanAce/orographic-precipitation/workflows/test/badge.svg?branch=master
   :target: https://github.com/EstebanAce/orographic-precipitation/actions
   :alt: Build Status
.. |Notebooks status| image:: https://github.com/EstebanAce/orographic-precipitation/workflows/Test%20notebooks/badge.svg
   :target: https://github.com/EstebanAce/orographic-precipitation/actions
   :alt: Notebooks status
.. |Coverage| image:: https://img.shields.io/coveralls/github/rlange2/orographic-precipitation/master
   :target: https://coveralls.io/github/rlange2/orographic-precipitation?branch=master
   :alt: Coverage Status
.. |License| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT

.. _`Python`: https://www.python.org
.. _`Smith & Barstad (2004)`: https://journals.ametsoc.org/doi/full/10.1175/1520-0469%282004%29061%3C1377%3AALTOOP%3E2.0.CO%3B2

The model includes airflow dynamics, condensed water advection, and downslope
evaporation. It consists of two vertically-integrated steady-state advection
equations describing: (i) the cloud water density and (ii) the hydrometeor
density. Solving these equations using Fourier transform techniques,
derives a single formula relating terrain and precipitation.

Please refer to the original manuscript of `Smith & Barstad (2004)`_ to understand
the model physics and limitations.

Installation
------------

Required dependencies:

* Python 3.9 or later
* `numpy`_

Optional dependencies (required for fastscape extension):

* `fastscape`_
* `xarray-simlab`_

.. _`numpy`: https://numpy.org
.. _`fastscape`: https://github.com/fastscape-lem/fastscape
.. _`xarray-simlab`: https://github.com/benbovy/xarray-simlab

To install the package simply go to terminal and use pip to install the latest version from our GitHub repository:

.. code:: shell

    pip install git+https://github.com/fastscape-lem/orographic-precipitation.git

Usage
-----
Some examples are shown in the `notebooks`_ folder (Jupyter Notebooks).

.. _`notebooks`: https://github.com/fastscape-lem/orographic-precipitation/tree/master/orographic_precipitation/notebooks

Acknowledgement
---------------

This project is supported by the `Earth Surface Process Modelling`_ group of
the German Research Centre for Geosciences in Potsdam, Germany.

.. _`Earth Surface Process Modelling`: http://www.gfz-potsdam.de/en/section/earth-surface-process-modelling/
