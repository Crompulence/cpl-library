
# Fortran and C++ compilers
F90 = mpif90
CPP = mpic++

# Base flags
FFLAGS = -fno-underscoring
CFLAGS = -std=c++14

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
