.. cpl documentation master file, created by
   sphinx-quickstart on Tue May 17 17:27:47 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to CPL library's documentation
======================================

CPL Library is a communications and topology management system for coupling any continuum fluid dynamics (CFD) solver to any molecular dynamics (MD) code. CPL Library is free and open-source, please see `CPL <http://www.cpl-library.org/download.shtml>`_ for more information.

First Steps with CPL
====================

This document list the interface for CPL library in python, C/C++ and Fortran.

Install CPL library
-------------------

Install CPL library by following instructions on our website
`CPL <http://www.cpl-library.org/download.shtml>`_ or clone with ::

   $ git clone https://github.com/Crompulence/cpl-library.git


.. toctree::
    :maxdepth: 2

CPL python Bindings
-------------------

.. autoclass:: cplpy.CPL
    :members:

CPL C++ Bindings
----------------

.. doxygennamespace:: CPL
    :members:
    :project: cpp_cpl

CPL Fortran Bindings
--------------------

.. f:autofunction:: CPL_init
.. f:autofunction:: cpl_setup_md
.. f:autofunction:: cpl_setup_cfd
.. f:autofunction:: cpl_send
.. f:autofunction:: cpl_recv
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
.. f:autofunction:: cpl_get


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

