#!/bin/bash
WEBSITE_DIR=$CPL_PATH/website
WEBSITE_DEMO_DIR=$WEBSITE_DIR/user-guide/demo-codes/
EXAMPLES_DIR=$CPL_PATH/examples

#Start with sendrecv_globcell
SENDRECV_DIR=$EXAMPLES_DIR/sendrecv_globcell
FORT_DIR=$SENDRECV_DIR/fortran
CPP_DIR=$SENDRECV_DIR/cpp
PYTHON_DIR=$SENDRECV_DIR/python
source-highlight --src-lang fortran --out-format html --input $FORT_DIR/cfd_send_cells.f90 --output $WEBSITE_DEMO_DIR/cfd_sendrecv_fort.html
source-highlight --src-lang fortran --out-format html --input $FORT_DIR/md_send_cells.f90 --output $WEBSITE_DEMO_DIR/md_sendrecv_fort.html
source-highlight --src-lang C --out-format html --input $CPP_DIR/cfd_send_cells.cpp --output $WEBSITE_DEMO_DIR/cfd_sendrecv_cpp.html
source-highlight --src-lang C --out-format html --input $CPP_DIR/md_recv_cells.cpp --output $WEBSITE_DEMO_DIR/md_sendrecv_cpp.html
source-highlight --src-lang python --out-format html --input $PYTHON_DIR/cfd_send_cells.py --output $WEBSITE_DEMO_DIR/cfd_sendrecv_python.html
source-highlight --src-lang python --out-format html --input $PYTHON_DIR/md_recv_cells.py --output $WEBSITE_DEMO_DIR/md_sendrecv_python.html

#generate topology example
SENDRECV_DIR=$EXAMPLES_DIR/topology_plot_example
FORT_DIR=$SENDRECV_DIR/fortran
CPP_DIR=$SENDRECV_DIR/cpp
PYTHON_DIR=$SENDRECV_DIR/python
source-highlight --src-lang fortran --out-format html --input $FORT_DIR/md_send_cells.f90 --output $WEBSITE_DEMO_DIR/md_send_topo_fortran.html
source-highlight --src-lang C --out-format html --input $CPP_DIR/md_send_cells.cpp --output $WEBSITE_DEMO_DIR/md_send_topo_cpp.html
source-highlight --src-lang python --out-format html --input $PYTHON_DIR/CFD_recv_and_plot.py --output $WEBSITE_DEMO_DIR/cfd_recv_topo_python.html