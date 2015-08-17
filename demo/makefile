# Executable paths (to be built)
MDEXE = run/md/mdProgram
CFDEXE = run/cfd/cfdProgram

# MPI compiler wrappers
CPP = mpic++
F90 = mpif90

# Include header/module and library flags
CPL_INCDIR = ../cpl/include
CPL_LIBDIR = ../cpl/lib
CPL_LIB = -lcpl
FORTRAN_LIBS = -lgfortran -lmpifort
INCS = -I$(CPL_INCDIR)
LIBS = -L$(CPL_LIBDIR) $(CPL_LIB)

# Targets
default: cpp-fortran 
cpp-cpp: md_cpp cfd_cpp
cpp-fortran: md_cpp cfd_fortran
fortran-cpp: md_fortran cfd_cpp
fortran-fortran: md_fortran cfd_fortran
md_cpp:
	$(CPP) md/md.cpp $(INCS) $(LIBS) $(FORTRAN_LIBS) -o $(MDEXE)
md_fortran:
	$(F90) md/md.f90 $(INCS) $(LIBS) -o $(MDEXE) 
cfd_cpp:
	$(CPP) cfd/cfd.cpp $(INCS) $(LIBS) $(FORTRAN_LIBS) -o $(CFDEXE) 
cfd_fortran:
	$(F90) cfd/cfd.f90 $(INCS) $(LIBS) -o $(CFDEXE) 
clean:
	rm $(MDEXE) $(CFDEXE)
