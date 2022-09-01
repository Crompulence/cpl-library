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
PREFIX = /usr/local

#                          Definitions

# Directories
testdir = ./test
objdir = obj
srcdir = src
libdir = lib
bindir = bin

#TODO: Add here 3partylibs dir for json-fortran and fdict
includedir = include/cpl
coresrcdir = $(srcdir)/core
binddir = $(srcdir)/bindings
utilsdir = $(srcdir)/utils
fbinddir = $(binddir)/fortran
cbinddir = $(binddir)/c
cppbinddir = $(binddir)/cpp

# Targets 
lib = $(libdir)/libcpl.so

# Check for building json support
# TODO: Maybe move it to an inc file.
3rdpartybuild = false 
ifdef json-support
	3rdpartybuild = true
	io_src = CPL_io.f90 CPL_write.f90
	# Conditional code preprocessor macros. Could have been previously defined.
	BUILDPPROCMACROS += JSON_SUPPORT
	LFLAGS += -ljsonfortran
else
	io_src = CPL_write.f90
endif
	
# Source files, headers and objects
coresrc = CPL_module.f90 CPL_methods.f90 commondefs.f90 $(io_src)

coresrcfiles = $(addprefix $(coresrcdir)/, $(coresrc))
coreobjfiles = $(addprefix $(objdir)/, $(coresrc:.f90=.o))

fbindsrc = CPL.f90
fbindsrcfile = $(addprefix $(fbinddir)/, $(fbindsrc))
fbindobjfile = $(addprefix $(objdir)/, $(fbindsrc:.f90=.o))

cbindhdrfile = $(cbinddir)/CPLC.h
cbindsrcFfile = $(cbinddir)/CPLC.f90  
cbindobjFfile = $(objdir)/CPLF.o
cbindsrcCfile = $(cbinddir)/CPLC.c
cbindobjCfile = $(objdir)/CPLC.o

cppbindsrc = cpl.cpp
cppbindhdr = $(cppbindsrc:.cpp=.h)
cppbindsrcfiles = $(addprefix $(cppbinddir)/, $(cppbindsrc))
cppbindhdrfiles = $(addprefix $(cppbinddir)/, $(cppbindhdr))
cppbindobjfiles = $(addprefix $(objdir)/, $(cppbindsrc:.cpp=.o))

utilssrc = CPL_ndArray.cpp CPL_cartCreate.cpp CPL_vector3D.cpp CPL_force.cpp CPL_field.cpp  TransmittingField.cpp CPLSocket.cpp #  CPL_usherBase.cpp
utilsextrahdr = CPL_misclib.h PoolElement.h #CPL_usherExceptions.h
utilshdr = $(utilssrc:.cpp=.h) $(utilsextrahdr)
utilssrcfiles = $(addprefix $(utilsdir)/, $(utilssrc))
utilshdrfiles = $(addprefix $(utilsdir)/, $(utilshdr))
utilsobjfiles = $(addprefix $(objdir)/, $(utilssrc:.cpp=.o))

allobjfiles = $(coreobjfiles) $(fbindobjfile) $(cbindobjFfile) $(cbindobjCfile) \
              $(cppbindobjfiles) $(utilsobjfiles)

# Define flags and compilers/linkers specific to this makefile
ifdef DEBUG
	FFLAGS = -O0 -fbacktrace -Wall -fbounds-check  -cpp
	CFLAGS = -std=c++11 -O0 -g3 -fvar-tracking -msse2
endif
FFLAGS += $(FMODKEY)$(includedir) -fPIC
CFLAGS += -fPIC

#                         Targets

# Default: make both the fortran and the c libraries
default: core fortran copyutilities c cpp utilities link updatebin

#debug: core fortran copyutilities c cpp utilities link
#	FFLAGS = -O0 -fbacktrace -Wall -fbounds-check $(FMODKEY)$(includedir) -fPIC 

# Declare phony targets
.phony.: fortran c cpp utilities
fortran: core $(fbindobjfile) 
c: core $(fbindobjfile) $(cbindobjFfile) $(cbindobjCfile)
cpp: core $(fbindobjfile) $(cbindobjFfile) $(cbindobjCfile) $(cppbindobjfiles) $(utilsobjfiles)

# Fortran bindings
$(fbindobjfile): core $(fbindsrcfile)
ifdef BUILDPPROCMACROS
	$(F90) $(FFLAGS) -D$(BUILDPPROCMACROS) -c $(fbindsrcfile) -o $(fbindobjfile)
else
	$(F90) $(FFLAGS) -c $(fbindsrcfile) -o $(fbindobjfile)
endif

# C bindings: create the lib objects first, overwrite lib including CPLC
$(cbindobjFfile): core $(cbindsrcFfile)
ifdef BUILDPPROCMACROS
	$(F90) $(FFLAGS) -D$(BUILDPPROCMACROS) -c $(cbindsrcFfile) -o $(cbindobjFfile)
else
	$(F90) $(FFLAGS) -c $(cbindsrcFfile) -o $(cbindobjFfile)
endif
	@cp $(cbindhdrfile) $(includedir)

$(cbindobjCfile): core $(cbindsrcCfile)
ifdef BUILDPPROCMACROS
	$(CPP) $(CFLAGS) -D$(BUILDPPROCMACROS) -c $(cbindsrcCfile) -o $(cbindobjCfile)
else
	$(CPP) $(CFLAGS) -c $(cbindsrcCfile) -o $(cbindobjCfile)
endif
	@cp $(cbindhdrfile) $(includedir)



# C++ bindings: create lib and c bindings first, overwrite lib including CPLCPP
$(cppbindobjfiles): core $(cbindobjFfile) $(cbindobjCfile) 
ifdef BUILDPPROCMACROS
	$(CPP) $(CFLAGS) -D$(BUILDPPROCMACROS) -I$(cbinddir) -c $(cppbindsrcfiles) -o $(cppbindobjfiles)
else
	$(CPP) $(CFLAGS) -I$(cbinddir) -I$(includedir) -c $(cppbindsrcfiles) -o $(cppbindobjfiles)
endif
	@cp $(cppbindhdrfiles) $(includedir)

# Utilities
utilities: core $(utilsobjfiles) $(utilshdrfiles)
	@cp $(utilshdrfiles) $(includedir)

copyutilities: core $(utilshdrfiles)
	@cp $(utilshdrfiles) $(includedir)

# Directory rules 
$(objdir):
	mkdir -p $(objdir)
$(libdir):
	mkdir -p $(libdir)
$(includedir):
	mkdir -p $(includedir)
$(CPL_THIRD_PARTY_INC):
	mkdir -p $(CPL_THIRD_PARTY_INC)
$(CPL_THIRD_PARTY_LIB):
	mkdir -p $(CPL_THIRD_PARTY_LIB)

# Compilation rules for library object files (written in Fortran)
core: $(objdir) $(libdir) $(includedir) $(coreobjfiles)
$(coreobjfiles): $(objdir)/%.o : $(coresrcdir)/%.f90
ifeq ($(3rdpartybuild),true)
		$(F90) $(FFLAGS) -I$(CPL_THIRD_PARTY_INC) -c $< -o $@
else
		$(F90) $(FFLAGS) -c $< -o $@
endif

# Utils compilation rules
$(utilsobjfiles): $(objdir)/%.o : $(utilsdir)/%.cpp
	$(CPP) $(CFLAGS) -I$(utilsdir) -I$(includedir) -c $< -o $@

# Link static lib to dynamic (shared) library
link: $(objdir) $(libobjfiles) $(utilsobjfiles)
ifeq ($(3rdpartybuild),true)
		$(LINK) $(LSHAREDLIB) -o $(lib) $(allobjfiles) -L$(CPL_THIRD_PARTY_LIB) -Wl,-rpath $(CPL_THIRD_PARTY_LIB) $(LFLAGS) 
else
		$(LINK) $(LSHAREDLIB) -o $(lib) $(allobjfiles) $(LFLAGS) 
endif


#Replace default compilers in cplf90 and cplc++ based on 
#platform specified values for compilers
updatebin:
		@sed -i 's/FC="mpif90"/FC="$(F90)"/' bin/cplf90
		@sed -i 's/CXX="mpic++"/CXX="$(CPP)"/' bin/cplc++


linkconda: $(objdir) $(libobjfiles) $(utilsobjfiles)
		$(LINK) $(LSHAREDLIB) -o $(lib) $(allobjfiles) -Wl, -rpath $(LFLAGS) 

json-fortran: $(CPL_THIRD_PARTY_LIB) $(CPL_THIRD_PARTY_INC)
	bash $(MAKEINCPATH)/json-fortran.build

3rd-party: json-fortran

test-all: test-pytest-mapping test-pytest-initialisation test-examples test-valgrind test-gtests test_Dragmodels
	echo "Running all test"
	
test-pytest:
	pytest -v $(testdir)/pytests
	
test-pytest-mapping:
	pytest -v --fulltrace $(testdir)/pytests/mapping

test-pytest-initialisation:
	pytest -v $(testdir)/pytests/initialisation

test-examples:
	pytest -v --fulltrace $(testdir)/examples
#	./examples/sendrecv_globcell/test_all.sh
#	./examples/sendrecv_globcell/test_all_port.sh

test-valgrind:
	pytest -v  $(testdir)/valgrind

test-granular:
	cd ${testdir}/granular/suzuki/ && pytest -v test_column.py


test-gtests: CPL_force_unittest
	cd $(testdir)/gtests/ && ./CPL_force_unittest

CPL_force_unittest: 
	make -C $(testdir)/gtests

test_Dragmodels:
	cd $(testdir)/drag && pytest -v ./
	#pytest -v $(testdir)/drag

examples-coupled:
	pytest -v ./examples/coupled

webdocs-api:
	bash ./utils/update-webdocs-api.sh

webdocs-examples:
	bash ./utils/update-webdocs-examples.sh

webdocs-all:
	bash ./utils/update-webdocs-api.sh
	bash ./utils/update-webdocs-examples.sh


.PHONY: install
install:
	mkdir -p $(PREFIX)/$(libdir)
	mkdir -p $(PREFIX)/$(includedir)
	mkdir -p $(PREFIX)/$(bindir)
	cp ./$(bindir)/cplexec $(PREFIX)/$(bindir)/
	cp ./$(bindir)/cplf90 $(PREFIX)/$(bindir)/
	cp ./$(bindir)/cplc++ $(PREFIX)/$(bindir)/
	cp -r ./$(libdir)/libcpl.so $(PREFIX)/$(libdir)/
	cp -r ./$(includedir) $(PREFIX)/$(includedir)/
	ldconfig -l $(PREFIX)/$(libdir)/libcpl.so

.PHONY: uninstall
uninstall:
	rm -f $(PREFIX)/lib/libcpl.so
	rm -f $(PREFIX)/include/cpl/*
	rmdir $(PREFIX)/include/cpl
	rm $(PREFIX)/bin/cplexec
	rm $(PREFIX)/bin/cplf90
	rm $(PREFIX)/bin/cplc++

# Clean
clean:
	rm -rf $(objdir) $(libdir) $(includedir) ./*.mod

clean-all:
	rm -rf $(objdir) $(libdir) $(includedir) ./*.mod
	bash $(CPL_THIRD_PARTY)/clean.sh
