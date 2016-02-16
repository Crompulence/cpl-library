#!/bin/sh
set -e
sudo apt-get update -q
case $1 in
  gfortran) set -x;
    sudo apt-get install gfortran gfortran-5 -y;;
  ifort) set -x;
    echo "Intel is not a supported fortran compiler for tests:" $1; exit 1;;
  *)
    echo "Unknown fortran compiler:" $1; exit 1;;
esac






