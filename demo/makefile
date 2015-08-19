
# Get compilers and base flags from platform include file 
MAKEINCPATH= ../make
include $(MAKEINCPATH)/platform.inc

# Executable paths (to be built)
MDEXE = run/md/mdProgram
CFDEXE = run/cfd/cfdProgram

# Include header/module and library flags
PWD = $(shell pwd)
CPL_INCDIR = $(PWD)/../cpl/include
CPL_LIBDIR = $(PWD)/../cpl/lib
CPL_INC = -I$(CPL_INCDIR)
CPL_LIB = -L$(CPL_LIBDIR) -Wl,-rpath=$(CPL_LIBDIR) -lcpl 
FORTRAN_LIBS = -lgfortran -lmpifort
C_LIBS = -lc -lstdc++

# Targets
default: cpp-fortran 
cpp-cpp: md_cpp cfd_cpp
cpp-fortran: md_cpp cfd_fortran
fortran-cpp: md_fortran cfd_cpp
fortran-fortran: md_fortran cfd_fortran
md_cpp:
	$(CPP) $(CFLAGS) md/md.cpp $(CPL_INC) $(CPL_LIB) $(FORTRAN_LIBS) -o $(MDEXE)
md_fortran:
	$(F90) $(FFLAGS) md/md.f90 $(CPL_INC) $(CPL_LIB) $(C_LIBS) -o $(MDEXE) 
cfd_cpp:
	$(CPP) $(CFLAGS) cfd/cfd.cpp $(CPL_INC) $(CPL_LIB) $(FORTRAN_LIBS) -o $(CFDEXE) 
cfd_fortran:
	$(F90) $(FFLAGS) cfd/cfd.f90 $(CPL_INC) $(CPL_LIB) $(C_LIBS) -o $(CFDEXE) 
clean:
	rm $(MDEXE) $(CFDEXE)
