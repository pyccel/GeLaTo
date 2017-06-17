Introduction
============

A python library for the Generalized Locally Toeplitz theory for IsoGeometric Analysis.

Installation
************

Installation is done using::

  sudo python setup.py install

If you are a member of the **CLAPP** framework, you may need to add the **CLAPP_DIR** as a *prefix*::

  sudo python setup.py install --prefix=$CLAPP_DIR

Dependencies
^^^^^^^^^^^^

In order to export notebooks in a **rest** format, you need to install **pandoc**::

  sudo apt install pandoc

Tutorials
*********

Different notebooks have been implemented, where we explain step by step how to use **glt** for diffrent problems. To run them, you may need to install **jupyter**::

  sudo apt-get install ipython-notebook
  sudo -H pip install jupyter

then run it using::

  jupyter notebook

Converting notebooks to rst
^^^^^^^^^^^^^^^^^^^^^^^^^^^

run::

  jupyter nbconvert --to=rst --execute getting_started_1.ipynb
