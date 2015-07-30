#
#    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________
#     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________
#      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________
#       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________
#        _\/\\\_____________\/\\\/////////____\/\\\_____________
#         _\//\\\____________\/\\\_____________\/\\\_____________
#          __\///\\\__________\/\\\_____________\/\\\_____________
#           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_
#            _______\/////////__\///______________\///////////////__
#


#                            Includes

# Platform include file location, imports C-language and Fortran compilers
MAKEINCPATH= ../make
include $(MAKEINCPATH)/platform.inc


#                          Definitions

# Directories
objdir = obj
srcdir = src
libdir = lib
includedir = include
libsrcdir = $(srcdir)/libsrc
binddir = $(srcdir)/bindings
utilsdir = $(srcdir)/utils
fbinddir = $(binddir)/fortran
cbinddir = $(binddir)/c
cppbinddir = $(binddir)/cpp

# Targets 
lib = $(libdir)/libcpl.a

# Source files, headers and objects
libsrc = CPL_module.f90 CPL_methods.f90
libsrcfiles = $(addprefix $(libsrcdir)/, $(libsrc))
libobjfiles = $(addprefix $(objdir)/, $(libsrc:.f90=.o))
libmodfiles = $(addprefix $(objdir)/, *.mod)

fbindsrc = $(fbinddir)/CPL.f90
fbindobj = $(objdir)/CPL.o

cbindheader = $(cbinddir)/CPLC.h
cbindsrc = $(cbinddir)/CPLC.f90
cbindobj = $(objdir)/CPLC.o

cppbindheader = $(cppbinddir)/CPLCPP.h $(cppbinddir)/CPL.h
cppbindsrc = $(cppbinddir)/CPLCPP.cpp
cppbindobj = $(objdir)/CPLCPP.o

utilssrc = CPL_ndArray.cpp CPL_cartCreate.cpp
utilsheaders = $(utilssrc:.cpp=.h) 
utilssrcfiles = $(addprefix $(utilsdir)/, $(utilssrc))
utilsheaderfiles = $(addprefix $(utilsdir)/, $(utilsheaders))
utilsobjfiles = $(addprefix $(objdir)/, $(utilssrc:.cpp=.o))



#                          Flags 
FFLAGS += -J$(objdir) -fPIC
CFLAGS += -fPIC



#                         Targets

# Default: make both the fortran and the c libraries
all: fortran c cpp utilities

# This makefile requires an argument, help if none provided
help:
	@echo "Target must be specified as either:"
	@echo ""
	@echo "     all"
	@echo "     fortran"
	@echo "     c"
	@echo "     cpp"
	@echo "     utilities"
	@echo ""

# Declare phony targets
.phony.: fortran c cpp utilities
fortran: $(lib)
c: $(lib) $(cbindobj)
cpp: $(lib) $(cbindobj) $(cppbindobj) $(utils)

# Fortran bindings
$(lib): $(objdir) $(libdir) $(includedir) $(libobjfiles) $(fbindsrc)
	$(F90) $(FFLAGS) -c $(fbindsrc) -o $(fbindobj)
	@ar cr $(lib) $(libobjfiles)
	@cp $(libmodfiles) $(includedir)

# C bindings: create the lib objects first, overwrite lib including CPLC
$(cbindobj): $(lib)
	$(F90) $(FFLAGS) -c $(cbindsrc) -o $(cbindobj)
	@ar cr $(lib) $(libobjfiles) $(cbindobj)
	@cp $(cbindheader) $(includedir)

# C++ bindings: create lib and c bindings first, overwrite lib including CPLCPP
$(cppbindobj): $(lib) $(cbindobj)
	$(CPP) $(CFLAGS) -I$(cbinddir) -c $(cppbindsrc) -o $(cppbindobj)
	@ar cr $(lib) $(libobjfiles) $(cbindobj) $(cppbindobj)
	@cp $(cppbindheader) $(includedir)

# Utilities
utilities: $(lib) $(utilsobjfiles)
	@ar cr $(lib) $(utilsobjfiles) 
	@cp $(utilsheaderfiles) $(includedir)

# Targets
$(objdir):
	mkdir -p $(objdir)
$(libdir):
	mkdir -p $(libdir)
$(includedir):
	mkdir -p $(includedir)

# Compilation rule for library object files (written in Fortran)
$(libobjfiles): $(objdir)/%.o : $(libsrcdir)/%.f90
	$(F90) $(FFLAGS) -c $< -o $@

$(utilsobjfiles): $(objdir)/%.o : $(utilsdir)/%.cpp
	$(CPP) $(CFLAGS) -c $< -o $@

# Clean
clean:
	rm -rf $(objdir) $(libdir) $(includedir)
