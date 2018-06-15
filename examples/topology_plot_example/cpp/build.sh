#!/bin/bash

# Build md side
#mpic++ md_send_cells.cpp -std=c++11 -I$CPL_PATH/include  -L$CPL_PATH/lib  -Wl,-rpath=$CPL_PATH/lib/ -lcpl  -o ./md
cplc++ md_send_cells.cpp -o ./md
