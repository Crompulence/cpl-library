#!/bin/bash
PWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [ -z "$CPL_PATH" ] 
then
    CPL_PATH=${PWD}
    echo "CPL_PATH variable is not set. Trying "${CPL_PATH}
else
    echo "CPL_PATH is already set to "${CPL_PATH}
	if [ "$CPL_PATH" == "$PWD" ]; 
	then
		echo "No need to update, set CPL_PATH manually if another location than $PWD preferred"
	else
		read -r -p "Do you want to replace with ${PWD}? [y/N] " response
		case "$response" in
			[yY][eE][sS]|[yY]) 
				CPL_PATH=${PWD}
				;;
			*)
				;;
		esac
	fi
fi
export CPL_PATH
export CPL_BIN_PATH="$CPL_PATH/bin"
PATH=${CPL_BIN_PATH}:$PATH 
export CPL_LIBRARY_PATH="$CPL_PATH/lib"
export CPL_INCLUDE_PATH="$CPL_PATH/include"
if [ -z "$PYTHONPATH" ] 
then
    export PYTHONPATH="$CPL_PATH/src/bindings/python"
else
    export PYTHONPATH="$CPL_PATH/src/bindings/python:$PYTHONPATH"
fi
#Add utilities
export PYTHONPATH="$CPL_PATH/utils:$PYTHONPATH"

export CPLPY_PATH="$CPL_PATH/src/bindings/python"

# Third-party libs
export CPL_THIRD_PARTY="$CPL_PATH/3rd-party"
export CPL_THIRD_PARTY_INC="$CPL_THIRD_PARTY/include"
export CPL_THIRD_PARTY_LIB="$CPL_THIRD_PARTY/lib"
