#!/bin/bash
EXAMPLES_DIR=$CPL_PATH/examples/sendrecv_globcell
WEBSITE_DIR=$CPL_PATH/website
WEBSITE_DEMO_DIR=$WEBSITE_DIR/user-guide/demo-codes/
FORT_DIR=$EXAMPLES_DIR/fortran
CPP_DIR=$EXAMPLES_DIR/cpp
PYTHON_DIR=$EXAMPLES_DIR/python
source-highlight --src-lang fortran --out-format html --input $FORT_DIR/cfd_send_cells.f90 --output $WEBSITE_DEMO_DIR/cfd_sendrecv_fort.html 
source-highlight --src-lang fortran --out-format html --input $FORT_DIR/cfd_send_cells.f90 --output $WEBSITE_DEMO_DIR/cfd_sendrecv_fort.html 
source-highlight --src-lang C --out-format html --input $CPP_DIR/cfd_send_cells.cpp --output $WEBSITE_DEMO_DIR/cfd_sendrecv_cpp.html 
source-highlight --src-lang C --out-format html --input $CPP_DIR/md_recv_cells.cpp --output $WEBSITE_DEMO_DIR/md_sendrecv_cpp.html 
source-highlight --src-lang python --out-format html --input $PYTHON_DIR/md_recv_cells.py --output $WEBSITE_DEMO_DIR/md_sendrecv_python.html 
source-highlight --src-lang python --out-format html --input $PYTHON_DIR/md_recv_cells.py --output $WEBSITE_DEMO_DIR/md_sendrecv_python.html 
