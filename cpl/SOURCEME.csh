#!/bin/csh
if ( ! $?CPL_PATH ) then
    echo "CPL_PATH variable is not set."
else
    setenv CPL_LIBRARY_PATH $CPL_PATH/lib
    if ( ! $?PYTHONPATH ) then
        setenv PYTHONPATH $CPL_PATH/src/bindings/python
    else
        setenv  PYTHONPATH $CPL_PATH/src/bindings/python:$PYTHONPATH
    endif
    setenv CPLPY_PATH $CPL_PATH/src/bindings/python
endif
