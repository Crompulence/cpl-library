#!/bin/bash
if [ -z "$CPL_PATH" ]
then
    CPL_PATH=../
fi

WEBSITE_DIR=$CPL_PATH/website
WEBSITE_DEMO_DIR=$WEBSITE_DIR/user-guide/demo-codes/
EXAMPLES_DIR=$CPL_PATH/examples

#Start with init examples
INIT_DIR=$EXAMPLES_DIR/cpl_init
source-highlight --src-lang fortran --out-format html --input $INIT_DIR/md_init.f90 --output $WEBSITE_DEMO_DIR/md_init_fortran.html
source-highlight --src-lang C --out-format html --input $INIT_DIR/cfd_init.cpp --output $WEBSITE_DEMO_DIR/cfd_init_cpp.html
source-highlight --src-lang python --out-format html --input $INIT_DIR/md_init.py --output $WEBSITE_DEMO_DIR/md_init_python.html

#then do sendrecv_globcell
SENDRECV_DIR=$EXAMPLES_DIR/sendrecv_globcell
FORT_DIR=$SENDRECV_DIR/fortran
CPP_DIR=$SENDRECV_DIR/cpp
PYTHON_DIR=$SENDRECV_DIR/python
source-highlight --src-lang fortran --out-format html --input $FORT_DIR/cfd_send_cells.f90 --output $WEBSITE_DEMO_DIR/cfd_send_cells_fortran.html
source-highlight --src-lang fortran --out-format html --input $FORT_DIR/md_recv_cells.f90 --output $WEBSITE_DEMO_DIR/md_recv_cells_fortran.html
source-highlight --src-lang C --out-format html --input $CPP_DIR/cfd_send_cells.cpp --output $WEBSITE_DEMO_DIR/cfd_send_cells_cpp.html
source-highlight --src-lang C --out-format html --input $CPP_DIR/md_recv_cells.cpp --output $WEBSITE_DEMO_DIR/md_recv_cells_cpp.html
source-highlight --src-lang python --out-format html --input $PYTHON_DIR/cfd_send_cells.py --output $WEBSITE_DEMO_DIR/cfd_send_cells_python.html
source-highlight --src-lang python --out-format html --input $PYTHON_DIR/md_recv_cells.py --output $WEBSITE_DEMO_DIR/md_recv_cells_python.html

#next generate topology example
TOPO_DIR=$EXAMPLES_DIR/topology_plot_example
FORT_DIR=$TOPO_DIR/fortran
CPP_DIR=$TOPO_DIR/cpp
PYTHON_DIR=$TOPO_DIR/python
source-highlight --src-lang fortran --out-format html --input $FORT_DIR/md_send_cells.f90 --output $WEBSITE_DEMO_DIR/md_send_topo_fortran.html
source-highlight --src-lang C --out-format html --input $CPP_DIR/md_send_cells.cpp --output $WEBSITE_DEMO_DIR/md_send_topo_cpp.html
source-highlight --src-lang python --out-format html --input $PYTHON_DIR/CFD_recv_and_plot.py --output $WEBSITE_DEMO_DIR/CFD_recv_and_plot_python.html
source-highlight --src-lang python --out-format html --input $PYTHON_DIR/draw_grid.py --output $WEBSITE_DEMO_DIR/draw_grid_python.html

#after we generate interactive topology example
TOPOI_DIR=$EXAMPLES_DIR/interactive_plot_example
FORT_DIR=$TOPOI_DIR/fortran
CPP_DIR=$TOPOI_DIR/cpp
PYTHON_DIR=$TOPOI_DIR/python
source-highlight --src-lang fortran --out-format html --input $FORT_DIR/md_send_cells.f90 --output $WEBSITE_DEMO_DIR/Interactive_fortran.html
source-highlight --src-lang C --out-format html --input $CPP_DIR/md_send_cells.cpp --output $WEBSITE_DEMO_DIR/Interactive_cpp.html
source-highlight --src-lang python --out-format html --input $PYTHON_DIR/CFD_recv_and_plot_grid_interactive.py --output $WEBSITE_DEMO_DIR/Interactive_python.html

#generate MD CFD example
MDCFD_DIR=$EXAMPLES_DIR/MD_CFD
source-highlight --src-lang python --out-format html --input $MDCFD_DIR/cfd_oo.py --output $WEBSITE_DEMO_DIR/cfd_oo.html
source-highlight --src-lang python --out-format html --input $MDCFD_DIR/md_oo.py --output $WEBSITE_DEMO_DIR/md_oo.html
source-highlight --src-lang python --out-format html --input $MDCFD_DIR/cfd_cpl.py --output $WEBSITE_DEMO_DIR/cfd_cpl.html
source-highlight --src-lang python --out-format html --input $MDCFD_DIR/md_cpl.py --output $WEBSITE_DEMO_DIR/md_cpl.html


#Generate minimal mocks examples
MOCK_DIR=$EXAMPLES_DIR/minimal_send_recv_mocks/
source-highlight --src-lang python --out-format html --input $MOCK_DIR/minimal_CFD.py --output $WEBSITE_DEMO_DIR/minimal_CFD.html
source-highlight --src-lang python --out-format html --input $MOCK_DIR/minimal_MD.py --output $WEBSITE_DEMO_DIR/minimal_MD.html

