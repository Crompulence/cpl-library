
FILENAME = run_CPLForce_drag_models

all:
	mpic++ ${FILENAME}.cpp -std=c++11 -I${CPL_PATH}/include/cpl  -L${CPL_PATH}/lib  -Wl,-rpath=${CPL_PATH}/lib/ -lcpl  -o ${FILENAME}
clean:
	rm -f ${FILENAME}

##Add CPL library includes
#CPL_INC = -I$(CPL_PATH)/include/cpl
#CPL_LIB = -L$(CPL_LIBRARY_PATH) -Wl,-rpath=$(CPL_LIBRARY_PATH) -lcpl
#LIBS = ${CPL_LIB}
#INCS = ${CPL_INC}

## Flags passed to the C++ compiler.
#CXX = mpic++
#CXXFLAGS += -g -Wall -std=c++11 -Wextra -pthread

## House-keeping build targets.
#all : ./run_CPLForce_drag_models.cpp
#	${CXX} ${CXXFLAGS} ${LIBS} $(INCS) 
#clean :
#	rm -f ./run_CPLForce_drag_models



