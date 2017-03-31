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

#                    Compilers and flags 
# Get compilers and base flags from platform include file 
MAKEINCPATH= $(CPL_PATH)/make
include $(MAKEINCPATH)/platform.inc

#                          Definitions

# Useful source and build definitions
EXE = lmp_cpl
SRC = main.cpp CPLSocketLAMMPS.cpp
OBJ = $(SRC:.cpp=.o)
OBJDIR = ./obj
SRCDIR = ./src
BINDIR = ./bin
CONFIGDIR = ./config
SRCFILES = $(addprefix $(SRCDIR)/, $(SRC))
OBJFILES = $(addprefix $(OBJDIR)/, $(OBJ))
LAMMPSSRCDIR = $(LAMMPS_PATH)/src
LIBLAMMPS_CPL = $(LAMMPSSRCDIR)/liblammps_cpl.a

# Includes and libs
CPL_INC = -I$(CPL_PATH)/include
CPL_LIB = -L$(CPL_LIBRARY_PATH) -Wl,-rpath=$(CPL_LIBRARY_PATH) -lcpl
LAMMPS_INC = -I$(LAMMPSSRCDIR)
LAMMPS_LIB = -L$(LAMMPSSRCDIR) -llammps_cpl
LIBS = $(CPL_LIB) $(LAMMPS_LIB) $(FORTRAN_LIB) $(MPI_LIB)

#                         Targets

.DEFAULT: 
	@make $(EXE)

help:
	@echo "                                                  "
	@echo "Make targets:                                     "
	@echo "                                                  "
	@echo "    lmp_cpl                                       "
	@echo "    clean                                         "
	@echo "                                                  "

$(EXE): $(LIBLAMMPS_CPL) $(OBJDIR) $(BINDIR) $(OBJFILES)
	$(CPP) $(CFLAGS) $(OBJFILES) $(LIBS) -o $(EXE)
	@mv $(EXE) $(BINDIR)

# LAMMPS Library
$(LIBLAMMPS_CPL):
	@bash $(CONFIGDIR)/enable-packages.sh
	@cd $(LAMMPSSRCDIR) && make mode=lib cpl

# Directories
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

# Rules
$(OBJFILES): $(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CPP) $(CFLAGS) $(LAMMPS_INC) $(CPL_INC) -c $< -o $@

test-all:
	py.test -v ./test
	
test-velocities:
	py.test -v ./test/velocityP-C

test-forces:
	py.test -v ./test/forceC-P

install-user-cpl:
	cp -R lammps_packages/USER-CPL $(LAMMPS_PATH)/src
	cp $(CONFIGDIR)/Makefile.cpl $(LAMMPS_PATH)/src/MAKE

clean:
	rm -rf $(OBJDIR) $(BINDIR)

clean-all:
	cd $(LAMMPSSRCDIR) && make no-all &> /dev/null
	rm -rf $(OBJDIR) $(BINDIR)
	rm -f $(LIBLAMMPS_CPL)
