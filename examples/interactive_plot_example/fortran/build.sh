#!/bin/bash

# Build md side
mpif90 md_send_cells.f90  -I$CPL_PATH/include  -L$CPL_PATH/lib  -Wl,-rpath=$CPL_PATH/lib/ -lcpl  -o ./md