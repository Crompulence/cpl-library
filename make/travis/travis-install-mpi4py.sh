#!/bin/sh
set -e

export MPICC=$1/bin/mpicc
pip install mpi4py
