#!/bin/bash

# Build md side
mpic++ md_recv_cells.cpp -std=c++11 -I$CPL_PATH/include  -L$CPL_PATH/lib  -Wl,-rpath=$CPL_PATH/lib/ -lcpl  -o ./md
# Build cfd side
mpic++ cfd_send_cells.cpp  -std=c++11 -I$CPL_PATH/include  -L$CPL_PATH/lib  -Wl,-rpath=$CPL_PATH/lib/ -lcpl  -o ./cfd
