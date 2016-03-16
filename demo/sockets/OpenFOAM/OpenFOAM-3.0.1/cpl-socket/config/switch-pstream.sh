#!/bin/bash

# This script is meant to replace libPstream.so library in the OpenFOAM 
# installation with a copy of libCPLPstream.lib to avoid the problem 
# of MPI_Init() being called 2 times when running coupled simulations. 
# It seems that this step is not always necessary since the error does 
# not always appear if the correct pstream library (libCPLPstream.so)
# is picked at running time.

# Using this script once it performs explained above. Using it again,
# will restore the original libPstream.so library to allow standalone
# use of OpenFOAM.

cpl_socket=cpl-socket/lib/libCPLSocket.so
cpl_pstream=cpl-socket/lib/libCPLPstream.so
if [ -f $cpl_socket ]; then
    pstream=`ldd $cpl_socket | grep libPstream.so | cut -f 3 -d ' '`
    pstream_bak=$( dirname "${pstream}" )/libPstream.so.bak
    if [ -f $pstream_bak ]; then
        mv $pstream_bak $pstream
        echo "Info:"
        echo "   libPstream.so restored from libPstream.so.bak. OpenFOAM is ready for standalone simulations!"
    else
        mv $pstream $pstream_bak
        cp $cpl_pstream $pstream
        echo "Info:"
        echo "   libCPLPstream.so moved to libPstream.so. OpenFOAM is ready for coupled simulations!" 
    fi
else
    echo "Error:"
    echo "   cpl-socket/lib/libCPLSocket.so does not exists! Have you compile the socket?"
    exit 1
fi
