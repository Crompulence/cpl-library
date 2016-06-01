
# Get compilers and base flags from platform include file 
MAKEINCPATH= ../../make
include $(MAKEINCPATH)/platform.inc
#F_LIBS = -lgfortran -lmpifort
#C_LIBS = -lc 

# Executable paths (to be built)
MDEXE = run/md/mdProgram
CFDEXE = run/cfd/cfdProgram

# Include header/module and library/linker flags
PWD = $(shell pwd)
CPL_INCDIR = $(PWD)/../../cpl/include
CPL_LIBDIR = $(PWD)/../../cpl/lib
CPL_INC = -I$(CPL_INCDIR)
CPL_LIB = -L$(CPL_LIBDIR) -Wl,-rpath=$(CPL_LIBDIR) -lcpl 
CPL_FLAGS = $(CPL_INC) $(CPL_LIB)
LD_FLAGS = $(CPL_FLAGS)

# Source code and associated object definitions
MD_DIR = md
MD_CPP_DIR = $(MD_DIR)/cpp
MD_F90_DIR = $(MD_DIR)/fortran
MD_CPP_SRC = main.cpp
MD_F90_SRC = main.f90
MD_CPP_SRC_FILES = $(addprefix $(MD_CPP_DIR)/, $(MD_CPP_SRC)) 
MD_F90_SRC_FILES = $(addprefix $(MD_F90_DIR)/, $(MD_F90_SRC)) 
MD_CPP_OBJ_DIR = $(MD_CPP_DIR)/obj
MD_F90_OBJ_DIR = $(MD_F90_DIR)/obj
MD_CPP_OBJ_FILES = $(addprefix $(MD_CPP_OBJ_DIR)/, $(MD_CPP_SRC:.cpp=.o))
MD_F90_OBJ_FILES = $(addprefix $(MD_F90_OBJ_DIR)/, $(MD_F90_SRC:.f90=.o))

CFD_DIR = cfd
CFD_CPP_DIR = $(CFD_DIR)/cpp
CFD_F90_DIR = $(CFD_DIR)/fortran
CFD_CPP_SRC = main.cpp
CFD_F90_SRC = main.f90
CFD_CPP_SRC_FILES = $(addprefix $(CFD_CPP_DIR)/, $(CFD_CPP_SRC)) 
CFD_F90_SRC_FILES = $(addprefix $(CFD_F90_DIR)/, $(CFD_F90_SRC)) 
CFD_CPP_OBJ_DIR = $(CFD_CPP_DIR)/obj
CFD_F90_OBJ_DIR = $(CFD_F90_DIR)/obj
CFD_CPP_OBJ_FILES = $(addprefix $(CFD_CPP_OBJ_DIR)/, $(CFD_CPP_SRC:.cpp=.o))
CFD_F90_OBJ_FILES = $(addprefix $(CFD_F90_OBJ_DIR)/, $(CFD_F90_SRC:.f90=.o))

ALL_OBJ_DIRS =  $(MD_CPP_OBJ_DIR)  $(MD_F90_OBJ_DIR) \
               $(CFD_CPP_OBJ_DIR) $(CFD_F90_OBJ_DIR)

# Targets
default: fortran-cpp
cpp-cpp: md_cpp cfd_cpp
cpp-fortran: md_cpp cfd_fortran
fortran-cpp: md_fortran cfd_cpp
fortran-fortran: md_fortran cfd_fortran

md_cpp: $(MD_CPP_OBJ_DIR) $(MD_CPP_OBJ_FILES)
	$(CPP) $(MD_CPP_OBJ_FILES) $(LD_FLAGS) $(F_LIBS) -o $(MDEXE)

md_fortran: $(MD_F90_OBJ_DIR) $(MD_F90_OBJ_FILES)
	$(F90) $(MD_F90_OBJ_FILES) $(LD_FLAGS) $(C_LIBS) -o $(MDEXE) 

cfd_cpp: $(CFD_CPP_OBJ_DIR) $(CFD_CPP_OBJ_FILES)
	$(CPP) $(CFD_CPP_OBJ_FILES) $(LD_FLAGS) $(F_LIBS) -o $(CFDEXE) 

cfd_fortran: $(CFD_F90_OBJ_DIR) $(CFD_F90_OBJ_FILES)
	$(F90) $(CFD_F90_OBJ_FILES) $(LD_FLAGS) $(C_LIBS) -o $(CFDEXE) 

clean:
	rm -rf $(MDEXE) $(CFDEXE) $(ALL_OBJ_DIRS)

# Rules
$(MD_CPP_OBJ_DIR):
	mkdir -p $(MD_CPP_OBJ_DIR)
$(MD_F90_OBJ_DIR):
	mkdir -p $(MD_F90_OBJ_DIR)
$(CFD_CPP_OBJ_DIR):
	mkdir -p $(CFD_CPP_OBJ_DIR)
$(CFD_F90_OBJ_DIR):
	mkdir -p $(CFD_F90_OBJ_DIR)

$(MD_CPP_OBJ_FILES): $(MD_CPP_OBJ_DIR)/%.o: $(MD_CPP_DIR)/%.cpp
	$(CPP) $(CFLAGS) $(CPL_INC) -c $< -o $@ 
$(MD_F90_OBJ_FILES): $(MD_F90_OBJ_DIR)/%.o: $(MD_F90_DIR)/%.f90
	$(F90) $(FFLAGS) $(CPL_INC) -c $< -o $@ 
$(CFD_CPP_OBJ_FILES): $(CFD_CPP_OBJ_DIR)/%.o: $(CFD_CPP_DIR)/%.cpp
	$(CPP) $(CFLAGS) $(CPL_INC) -c $< -o $@ 
$(CFD_F90_OBJ_FILES): $(CFD_F90_OBJ_DIR)/%.o: $(CFD_F90_DIR)/%.f90
	$(F90) $(FFLAGS) $(CPL_INC) -c $< -o $@ 
