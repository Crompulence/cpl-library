#!/bin/sh
set -e
sudo apt-get update
case $1 in
  gfortran) set -x;
    sudo apt-get install gfortran;;
  ifort) set -x;
    echo "Intel not supported fortran compiler:" $1; exit 1;;
  *)
    echo "Unknown fortran compiler:" $1; exit 1;;
esac






