EXE:= sc
COMP:= g++
LIBS:= -std=c++11
SRC:= sphere_cube.cpp

default: all

all: ${SRC}
	${COMP} ${LIBS}  ${SRC} -o $(EXE)
lib: ${SRC}
	${COMP} ${LIBS} -shared -fPIC ${SRC} -o overlaplib.so
clean:
	rm sc
