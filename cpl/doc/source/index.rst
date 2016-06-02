.. cpl documentation master file, created by
   sphinx-quickstart on Tue May 17 17:27:47 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to cpl's documentation!
===============================

First Steps with Sphinx
=======================

This document is meant to give a tutorial-like overview of all common tasks
while using Sphinx.

The green arrows designate "more info" links leading to advanced sections about
the described task.


Install Sphinx
--------------

Install Sphinx, either from a distribution package or from
`PyPI <https://pypi.python.org/pypi/Sphinx>`_ with ::

   $ pip install Sphinx


.. toctree::
    :maxdepth: 2

.. autoclass:: cplpy.CPL
    :members:

.. doxygennamespace:: CPL
    :members:
    :project: cpp_cpl

.. f:autofunction:: CPL_init
.. f:autofunction:: cpl_setup_md
.. f:autofunction:: cpl_setup_cfd
.. f:autofunction:: cpl_proc_extents
.. f:autofunction:: cpl_my_proc_extents
.. f:autofunction:: cpl_proc_portion
.. f:autofunction:: cpl_my_proc_portion
.. f:autofunction:: cpl_map_coord2cell
.. f:autofunction:: cpl_map_cell2coord
.. f:autofunction:: cpl_get_no_cells
.. f:autofunction:: cpl_map_glob2loc_cell
.. f:autofunction:: cpl_get_olap_limits
.. f:autofunction:: cpl_get_cnst_limits
.. f:autofunction:: cpl_map_cfd2md_coord
.. f:autofunction:: cpl_map_md2cfd_coord
.. f:autofunction:: cpl_overlap


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

