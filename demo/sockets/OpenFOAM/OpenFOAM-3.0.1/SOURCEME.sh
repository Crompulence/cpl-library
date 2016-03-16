#!/bin/bash
#~~~
#    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________
#     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________
#      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________
#       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________
#        _\/\\\_____________\/\\\/////////____\/\\\_____________
#         _\//\\\____________\/\\\_____________\/\\\_____________
#          __\///\\\__________\/\\\_____________\/\\\_____________
#           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_
#            _______\/////////__\///______________\///////////////__
#~~~
#

# Environment variable for install directory
export FOAM_CPL_VERSION=3.0.x
export FOAM_INST_DIR=`pwd`
# Source the other environment variables
foamDotFile=$FOAM_INST_DIR/OpenFOAM-$FOAM_CPL_VERSION/etc/bashrc
echo $foamDotFile
if [ -f foamDotFile ]
then
    source $foamDotFile
else
    echo "ERROR:"
    echo "   Configuration file 'OpenFOAM-$FOAM_CPL_VERSION/etc/bashrc' not found."
    exit 1 
fi

echo ""
echo "FOAM_MPI environment variable is now: " $FOAM_MPI
export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/$FOAM_MPI

# CPL environment variables
echo ""
echo "New environment variables: "
echo ""

export FOAM_CPL_SOCKET=$FOAM_INST_DIR/cpl-socket
echo "    FOAM_CPL_SOCKET = " $FOAM_CPL_SOCKET

export FOAM_CPL_SOCKET_SRC=$FOAM_CPL_SOCKET/src
echo "    FOAM_CPL_SOCKET_SRC = " $FOAM_CPL_SOCKET_SRC

export FOAM_CPL_SOCKET_UTILS=$FOAM_CPL_SOCKET_SRC/utils
echo "    FOAM_CPL_SOCKET_UTILS = " $FOAM_CPL_SOCKET_UTILS

export FOAM_CPL_SOCKET_LIBBIN=$FOAM_CPL_SOCKET/lib
echo "    FOAM_CPL_SOCKET_LIBBIN = " $FOAM_CPL_SOCKET_LIBBIN

export FOAM_CPL_SOCKET_BIN=$FOAM_CPL_SOCKET/bin
echo "    FOAM_CPL_SOCKET_BIN = " $FOAM_CPL_SOCKET_BIN

echo ""

# Update paths
export PATH=$PATH:$FOAM_CPL_SOCKET_BIN
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FOAM_CPL_SOCKET_LIBBIN
echo "PATH updated to:"
echo "   $PATH:"$FOAM_CPL_SOCKET_BIN
echo
echo "LD_LIBRARY_PATH updated to: "
echo "   $LD_LIBRARY_PATH:"$FOAM_CPL_SOCKET_LIBBIN
