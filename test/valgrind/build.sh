#!/bin/bash

# Build md side
mpif90 array_stuff.f90 md_recvsend_cells.f90  -I$CPL_PATH/include  -L$CPL_PATH/lib  -Wl,-rpath=$CPL_PATH/lib/ -lcpl  -o ./md
# Build cfd side
mpif90 array_stuff.f90 cfd_sendrecv_cells.f90  -I$CPL_PATH/include  -L$CPL_PATH/lib  -Wl,-rpath=$CPL_PATH/lib/ -lcpl  -o ./cfd
