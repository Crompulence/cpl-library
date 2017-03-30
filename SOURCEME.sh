#!/bin/bash
if [ -z "$CPL_PATH" ] 
then
    CPL_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    echo "CPL_PATH variable is not set. Trying "${CPL_PATH}
fi
export CPL_PATH
export CPL_BIN_PATH="$CPL_PATH/bin"
PATH=${CPL_BIN_PATH}:$PATH 
export CPL_LIBRARY_PATH="$CPL_PATH/lib"
if [ -z "$PYTHONPATH" ] 
then
    export PYTHONPATH="$CPL_PATH/src/bindings/python"
else
    export PYTHONPATH="$CPL_PATH/src/bindings/python:$PYTHONPATH"
fi
export CPLPY_PATH="$CPL_PATH/src/bindings/python"

# Third-party libs
export CPL_THIRD_PARTY="$CPL_PATH/3rd-party"
export CPL_THIRD_PARTY_INC="$CPL_THIRD_PARTY/include"
export CPL_THIRD_PARTY_LIB="$CPL_THIRD_PARTY/lib"
