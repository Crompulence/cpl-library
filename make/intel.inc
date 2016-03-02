
# Fortran/C++ compilers and linker
F90 = mpif90
CPP = mpic++
LINK = mpic++

INTELDIR = /opt/intel/composer_xe_2013_sp1.2.144/compiler/lib/intel64/ # Place the path to the Intel  

# Base flags for Fortran/C++/Linker
FFLAGS = #-cxxlib
CFLAGS = -std=c++14
LFLAGS = -L$(INTELDIR) -lifport -lifcore -limf -lsvml -lintlc

# Compiler's flag for specifying where to place compiled Fortran modules
FMODKEY = -module 

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
