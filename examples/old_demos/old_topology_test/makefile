# Platform include file location
MAKEINCPATH= ../../make
include $(MAKEINCPATH)/platform.inc

exe = topologyTest
cleanables = $(exe) 
test_src = test_coupler.f90
flags = -O0 -traceback -g -check all
objs = $(addprefix $(obj_dir)/, $(test_src:.f90=.o))
libdir = ../../cpl/lib/
incdir = ../../cpl/include/
lib = ../../cpl/lib/libcpl.so
c_libs = -lc -lstdc++

## Compilation rules
$(exe): $(lib)
	$(F90) -o $(exe) $(test_src) -I$(incdir) -L$(libdir) -Wl,-rpath=$(libdir) -lcpl $(c_libs)

clean:
	rm -rf $(cleanables)
