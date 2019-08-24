Welcome to GeLaTo
=================

|build-status| |docs| |binder|


**GeLaTo** is a python library for the Generalized Locally Toeplitz theory for IsoGeometric Analysis.

Install
*******

From PyPi
^^^^^^^^^

Simply run, for a local installation::

  pip3 install --user gelato 

or::

  pip3 install gelato 

for a global installation.

From sources
^^^^^^^^^^^^

* **Standard mode**::

    python3 -m pip install .

* **Development mode**::

    python3 -m pip install --user -e .

Examples
********

   * `1- Laplace 2D <http://nbviewer.jupyter.org/github/pyccel/gelato/blob/master/notebooks/Laplace_2d.ipynb>`_

   * `2- Vibration 1D <http://nbviewer.jupyter.org/github/pyccel/gelato/blob/master/notebooks/Vibration_1d.ipynb>`_


.. |build-status| image:: https://travis-ci.org/pyccel/GeLaTo.svg?branch=master
    :alt: build status
    :scale: 100%
    :target: https://travis-ci.org/pyccel/GeLaTo

.. |docs| image:: https://readthedocs.org/projects/gelato/badge/?version=latest
    :target: http://gelato.readthedocs.io/en/latest/?badge=latest
    :scale: 100%
    :alt: Documentation Status

.. |binder| image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/pyccel/gelato/master
