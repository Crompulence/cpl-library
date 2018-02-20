#!/bin/bash

# Build md side
mpif90 md_send_cells.f90  -I$CPL_PATH/include/cpl  -L$CPL_PATH/lib  -Wl,-rpath=$CPL_PATH/lib/  -o ./md -lcpl
