#!/bin/sh
set -e

export MPICC=$1/bin/mpicc
sudo pip install mpi4py
