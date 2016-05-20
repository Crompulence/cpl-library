#!/bin/sh
set -e

env MPICC=$1/bin/mpicc
sudo pip install mpi4py
