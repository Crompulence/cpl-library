
Fortran Bindings
==================



=================
Main Functions
=================
In general usage, the aim of CPL library is to provide a minimal set of functions
to facilitate coupling. These include an initialisation, setup used on both CFD
and MD side then a send and recieve command. 
There is also a helper function to provide 

----------------------
CPL init
----------------------

.. f:autosubroutine:: cpl_init

----------------------
CPL Setup MD
----------------------

.. f:autosubroutine:: cpl_setup_md

----------------------
CPL Setup CFD
----------------------

.. f:autosubroutine:: cpl_setup_cfd

----------------------
CPL Recv
----------------------
Note two interfaces exist, CPL_recv_full and CPL_recv_min

.. f:autosubroutine:: cpl_recv_full

----------------------
CPL Send
----------------------
Note two interfaces exist, CPL_send_full and CPL_send_min  

.. f:autosubroutine:: cpl_send_full
  
----------------------
CPL Get Arrays
----------------------

.. f:autosubroutine:: CPL_get_arrays

=================
Full Listing
=================
The complete list of module functions is included here.

----------------------
Coupler Module
----------------------
Contains all the setup subroutines and variables
for a coupled run

.. f:automodule:: coupler_module
----------------------
Coupler Methods Module
----------------------
All the subroutines and functions which would be called
during a run, such as CPL_Send and CPL_recv

.. f:automodule:: coupler