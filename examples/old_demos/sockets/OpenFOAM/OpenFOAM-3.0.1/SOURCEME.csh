#!/bin/csh
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
setenv FOAM_CPL_VERSION 3.0.1
setenv FOAM_INST_DIR `pwd`
# Source the other environment variables
set foamDotFile=$FOAM_INST_DIR/OpenFOAM-$FOAM_CPL_VERSION/etc/cshrc
if( -f $foamDotFile ) then
    source $foamDotFile
else
    echo "ERROR:"
    echo "   Configuration file 'OpenFOAM-$FOAM_CPL_VERSION/etc/cshrc' not found."
    return 1
endif

echo ""
echo "FOAM_MPI environment variable is now: " $FOAM_MPI
setenv FOAM_MPI_LIBBIN $FOAM_LIBBIN/$FOAM_MPI

# CPL environment variables
echo ""
echo "New environment variables: "
echo ""

setenv FOAM_CPL_SOCKET $FOAM_INST_DIR/cpl-socket
echo "    FOAM_CPL_SOCKET = " $FOAM_CPL_SOCKET

setenv FOAM_CPL_SOCKET_SRC $FOAM_CPL_SOCKET/src
echo "    FOAM_CPL_SOCKET_SRC = " $FOAM_CPL_SOCKET_SRC

setenv FOAM_CPL_SOCKET_UTILS $FOAM_CPL_SOCKET_SRC/utils
echo "    FOAM_CPL_SOCKET_UTILS = " $FOAM_CPL_SOCKET_UTILS

setenv FOAM_CPL_SOCKET_LIBBIN $FOAM_CPL_SOCKET/lib
echo "    FOAM_CPL_SOCKET_LIBBIN = " $FOAM_CPL_SOCKET_LIBBIN

setenv FOAM_CPL_SOCKET_BIN $FOAM_CPL_SOCKET/bin
echo "    FOAM_CPL_SOCKET_BIN = " $FOAM_CPL_SOCKET_BIN

echo ""

# Update paths
setenv PATH $PATH\:$FOAM_CPL_SOCKET_BIN
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$FOAM_CPL_SOCKET_LIBBIN
echo "PATH updated to:"
echo "   $PATH\:$FOAM_CPL_SOCKET_BIN"
echo
echo "LD_LIBRARY_PATH updated to:"
echo "   $LD_LIBRARY_PATH\:$FOAM_CPL_SOCKET_LIBBIN"
