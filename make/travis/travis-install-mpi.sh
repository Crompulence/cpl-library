#!/bin/sh
set -e
if [ $# -lt 2 ]; then
    PWD=$(pwd)
    echo "No Build dir specified, attempting to build in " $PWD
    BUILD_DIR=$PWD
else
    BUILD_DIR=$2
    echo "Build dir specified as " $BUILD_DIR
fi
mkdir -p $BUILD_DIR
case $1 in
  mpich3) set -x;
    cd $BUILD_DIR
    wget -q --no-check-certificate http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
    tar -xzf mpich-3.2.tar.gz
    cd mpich-3.2
    mkdir build && cd build
    ../configure --prefix=$BUILD_DIR/mpich --disable-static \
                 --enable-fortran=option FC=gfortran \
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
  openmpi) set -x;
    cd $BUILD_DIR
    wget -q --no-check-certificate https://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-1.10.2.tar.gz
    tar -xvf openmpi-1.10.2.tar.gz
    cd openmpi-1.10.2
    mkdir build && cd build
    ../configure --prefix=$BUILD_DIR/open-mpi \
                FC=gfortran --enable-mpi-fortran=all 
    make -j4
    make install;;
    #Previous attempt using repo version doesn't work with repo fortran
    #sudo apt-get install -q openmpi-bin openmpi-common libopenmpi-dev;;
  *)
    echo "Unknown MPI implementation:" $1; exit 1;;
esac
