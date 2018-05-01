#!/bin/sh
set -e
#Check location of build
if [ $# -lt 2 ]; then
    PWD=$(pwd)
    echo "No Build dir specified, attempting to build in " $PWD
    BUILD_DIR=$PWD
else
    BUILD_DIR=$2
    echo "Build dir specified as " $BUILD_DIR
fi
INSTALL_DIR=$BUILD_DIR/$1
echo $INSTALL_DIR
#If file exist then use this
if [ -f $INSTALL_DIR/bin/mpirun ]; then
    echo "Using cached MPI in "$INSTALL_DIR
else
    if [ $# -lt 3 ]; then
        GCC_VERSION=5
    else
        GCC_VERSION=$3
    fi
    mkdir -p $BUILD_DIR
    case $1 in
      mpich3) set -x;
        cd $BUILD_DIR
        wget -q --no-check-certificate http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz

        tar -xzf mpich-3.2.1.tar.gz
        cd mpich-3.2.1
        mkdir build && cd build
        ../configure --prefix=$INSTALL_DIR --disable-static \
                     --enable-fortran=yes FC=gfortran-$GCC_VERSION \
                     --enable-cxx=yes CXX=g++-$GCC_VERSION CC=gcc-$GCC_VERSION \
                     --enable-g=all --enable-error-messages=all \
                     --enable-error-checking=all

        make -j4
        make install;;

        #Previous attempt using repo version doesn't work with repo fortran
        #sudo apt-get install -q libcr0 default-jdk;
        #wget -q http://www.cebacad.net/files/mpich/ubuntu/mpich-3.1/mpich_3.1-1ubuntu_amd64.deb;
        #sudo dpkg -i ./mpich_3.1-1ubuntu_amd64.deb;
        #rm -f ./mpich_3.1-1ubuntu_amd64.deb;;

      mpich1) set -x;
        sudo apt-get install -q mpich-shmem-bin libmpich-shmem1.0-dev;;
      mpich2) set -x;
        sudo apt-get install -q mpich2 libmpich2-3 libmpich2-dev;;
      openmpi3) set -x;
        cd $BUILD_DIR
        wget -q --no-check-certificate https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.1.tar.gz
        tar -xvf openmpi-3.0.1.tar.gz
        cd openmpi-3.0.1
        mkdir build && cd build
        ../configure --prefix=$INSTALL_DIR \
                    FC=gfortran --enable-mpi-fortran=all 
        make -j4
        make install;;
      openmpi2) set -x;
        cd $BUILD_DIR
        wget -q --no-check-certificate https://www.open-mpi.org/software/ompi/v2.0/downloads/openmpi-2.0.1.tar.gz
        tar -xvf openmpi-2.0.1.tar.gz
        cd openmpi-2.0.1
        mkdir build && cd build
        ../configure --prefix=$INSTALL_DIR \
                    FC=gfortran --enable-mpi-fortran=all 
        make -j4
        make install;;
      openmpi1) set -x;
        cd $BUILD_DIR
        wget -q --no-check-certificate https://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-1.10.7.tar.gz
        tar -xvf openmpi-1.10.7.tar.gz
        cd openmpi-1.10.7
        mkdir build && cd build
        ../configure --prefix=$INSTALL_DIR \
                    FC=gfortran-$GCC_VERSION --enable-mpi-fortran=all \
                    CC=gcc-$GCC_VERSION CXX=g++-$GCC_VERSION
        make -j4
        make install;;
        #Previous attempt using repo version doesn't work with repo fortran
        #sudo apt-get install -q openmpi-bin openmpi-common libopenmpi-dev;;
      *)
        echo "Unknown MPI implementation:" $1; exit 1;;
    esac
fi
