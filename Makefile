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

# Get compilers and base flags from platform include file 
MAKEINCPATH= ./make
include $(MAKEINCPATH)/platform.inc


#                          Definitions

# Directories
testdir = ./test
objdir = obj
srcdir = src
libdir = lib
includedir = include
coresrcdir = $(srcdir)/core
binddir = $(srcdir)/bindings
utilsdir = $(srcdir)/utils
fbinddir = $(binddir)/fortran
cbinddir = $(binddir)/c
cppbinddir = $(binddir)/cpp

# Targets 
lib = $(libdir)/libcpl.so

# Source files, headers and objects
coresrc = CPL_module.f90 CPL_methods.f90 
coresrcfiles = $(addprefix $(coresrcdir)/, $(coresrc))
coreobjfiles = $(addprefix $(objdir)/, $(coresrc:.f90=.o))

fbindsrc = CPL.f90
fbindsrcfile = $(addprefix $(fbinddir)/, $(fbindsrc))
fbindobjfile = $(addprefix $(objdir)/, $(fbindsrc:.f90=.o))

cbindhdrfile = $(cbinddir)/CPLC.h
cbindsrcfile = $(cbinddir)/CPLC.f90
cbindobjfile = $(objdir)/CPLC.o

cppbindsrc = CPLCPP.cpp
cppbindhdr = $(cppbindsrc:.cpp=.h) CPL.h
cppbindsrcfiles = $(addprefix $(cppbinddir)/, $(cppbindsrc))
cppbindhdrfiles = $(addprefix $(cppbinddir)/, $(cppbindhdr))
cppbindobjfiles = $(addprefix $(objdir)/, $(cppbindsrc:.cpp=.o))

utilssrc = CPL_ndArray.cpp CPL_cartCreate.cpp CPL_vector3D.cpp \
#           CPL_usherBase.cpp
utilsextrahdr = #CPL_usherExceptions.h
utilshdr = $(utilssrc:.cpp=.h) $(utilsextrahdr)
utilssrcfiles = $(addprefix $(utilsdir)/, $(utilssrc))
utilshdrfiles = $(addprefix $(utilsdir)/, $(utilshdr))
utilsobjfiles = $(addprefix $(objdir)/, $(utilssrc:.cpp=.o))

allobjfiles = $(coreobjfiles) $(fbindobjfile) $(cbindobjfile) \
              $(cppbindobjfiles) $(utilsobjfiles)

# Define flags and compilers/linkers specific to this makefile
FFLAGS += $(FMODKEY)$(includedir) -fPIC
CFLAGS += -fPIC

#                         Targets

# Default: make both the fortran and the c libraries
default: core fortran c cpp utilities link

debug: core fortran c cpp utilities link
    FFLAGS = -fPIC -O0 -fbacktrace -fbounds-check $(FMODKEY)$(includedir)

# Declare phony targets
.phony.: fortran c cpp utilities
fortran: core $(fbindobjfile) 
c: core $(fbindobjfile) $(cbindobjfile)
cpp: core $(fbindobjfile) $(cbindobjfile) $(cppbindobjfiles) $(utilsobjfiles)

# Fortran bindings
$(fbindobjfile): core $(fbindsrcfile)
	$(F90) $(FFLAGS) -c $(fbindsrcfile) -o $(fbindobjfile)

# C bindings: create the lib objects first, overwrite lib including CPLC
$(cbindobjfile): core $(cbindsrcfile)
	$(F90) $(FFLAGS) -c $(cbindsrcfile) -o $(cbindobjfile)
	@cp $(cbindhdrfile) $(includedir)

# C++ bindings: create lib and c bindings first, overwrite lib including CPLCPP
$(cppbindobjfiles): core $(cbindobjfile) 
	$(CPP) $(CFLAGS) -I$(cbinddir) -c $(cppbindsrcfiles) -o $(cppbindobjfiles)
	@cp $(cppbindhdrfiles) $(includedir)

# Utilities
utilities: core $(utilsobjfiles) 
	@cp $(utilshdrfiles) $(includedir)

# Directory rules 
$(objdir):
	mkdir -p $(objdir)
$(libdir):
	mkdir -p $(libdir)
$(includedir):
	mkdir -p $(includedir)

# Compilation rules for library object files (written in Fortran)
core: $(objdir) $(libdir) $(includedir) $(coreobjfiles)
$(coreobjfiles): $(objdir)/%.o : $(coresrcdir)/%.f90
	$(F90) $(FFLAGS) -c $< -o $@

# Utils compilation rules
$(utilsobjfiles): $(objdir)/%.o : $(utilsdir)/%.cpp
	$(CPP) $(CFLAGS) -I$(utilsdir) -c $< -o $@

# Link static lib to dynamic (shared) library
link: $(objdir) $(libobjfiles) $(utilsobjfiles)
	$(LINK) $(LSHAREDLIB) -o $(lib) $(allobjfiles) $(LFLAGS) 

test-all:
	py.test -v -s $(testdir)
	
test-mapping:
	py.test -v -s $(testdir)/mapping

test-initialisation:
	py.test -v -s $(testdir)/initialisation

test-examples:
	./examples/sendrecv_globcell/test_all.sh
	./examples/sendrecv_globcell/test_all_port.sh

test-valgrind:
	./test/valgrind/debug_all.sh

webdocs-api:
	bash ./utils/update-webdocs-api.sh

webdocs-examples:
	bash ./utils/update-webdocs-examples.sh

webdocs-all:
	bash ./utils/update-webdocs-api.sh
	bash ./utils/update-webdocs-examples.sh

# Clean
clean:
	rm -rf $(objdir) $(libdir) $(includedir) ./*.mod
