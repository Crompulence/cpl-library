

/home/es205/installed/openmpi-2.0.1_gcc/bin/mpif90 -fbacktrace -fbacktrace -g -ggdb -O0 -W -Wsurprising -Wampersand -Waliasing -fbounds-check -ffpe-trap=zero,overflow,underflow -fno-underscoring ./minimal_MD.f90  -I${CPL_PATH}/include/cpl -L${CPL_PATH}/lib -Wl,-rpath=${CPL_PATH}/lib/ -o ./f_MD -lcpl

/home/es205/installed/openmpi-2.0.1_gcc/bin/mpif90 -fbacktrace -fbacktrace -g -ggdb -O0 -W -Wsurprising -Wampersand -Waliasing -fbounds-check -ffpe-trap=zero,overflow,underflow -fno-underscoring ./minimal_CFD.f90  -I${CPL_PATH}/include/cpl -L${CPL_PATH}/lib -Wl,-rpath=${CPL_PATH}/lib/ -o ./f_CFD -lcpl



/home/es205/installed/openmpi-2.0.1_gcc/bin/mpic++ -g -O0 -std=c++11 ./minimal_MD.cpp  -I${CPL_PATH}/include/cpl -L${CPL_PATH}/lib -Wl,-rpath=${CPL_PATH}/lib/ -o ./c_MD -lcpl

/home/es205/installed/openmpi-2.0.1_gcc/bin/mpic++ -g -O0 -W -std=c++11 ./minimal_CFD.cpp  -I${CPL_PATH}/include/cpl -L${CPL_PATH}/lib -Wl,-rpath=${CPL_PATH}/lib/ -o ./c_CFD -lcpl

/home/es205/installed/openmpi-2.0.1_gcc/bin/mpiexec -n 1 ./f_MD : -n 1 ./f_CFD 

