
# Fortran/C++ compilers and linker
F90 = mpif90
CPP = mpic++
LINK = mpic++

# Base flags for Fortran/C++/Linker
FFLAGS = -fno-underscoring
CFLAGS = -std=c++14
LFLAGS = -lgfortran -lmpifort

# Compiler's flag for specifying where to place compiled Fortran modules
FMODKEY = -J

# Linker's flag indicating that objects are to be linked into a shared lib
LSHAREDLIB = -shared

# Build type, optimised or debug
ifndef BUILD
   BUILD = opt
endif
ifeq ($(strip $(BUILD)),debug)
    FFLAGS += -O0 -g -ggdb -Wall
    CFLAGS += -O0 -g -ggdb -Wall
else
    FFLAGS += -O3
    CFLAGS += -O3
endif
