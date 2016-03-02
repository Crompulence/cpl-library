#!/bin/sh
set -e
# Note bug report #5221 required apt-get update for gfortran
sudo apt-get update -q
case $1 in
  gfortran) set -x;
    sudo apt-get install gfortran -y;;
  gfortran-4.8) set -x;
    sudo apt-get install gfortran-4.8 -y;;
  gfortran-4.9) set -x;
    sudo apt-get install gfortran-4.9 -y;;
  gfortran-5) set -x;
    sudo apt-get install gfortran-5 -y;;
  ifort) set -x;
    echo "Intel is not a supported fortran compiler for tests:" $1; exit 1;;
  *)
    echo "Unknown fortran compiler:" $1; exit 1;;
esac






