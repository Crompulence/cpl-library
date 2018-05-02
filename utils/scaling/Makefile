#!/bin/bash
#module swap PrgEnv-cray/5.2.82 PrgEnv-gnu

default: md cfd
        
#	ftn -dynamic scaling_mod.f90 scaling_with_send_recv.f90  -I$(CPL_PATH)/include/cpl  -L$(CPL_PATH)/lib  -Wl,-rpath=$(CPL_PATH)/lib/ -lcpl  -o ./a.out

md: scaling_with_send_recv.f90
	ftn -dynamic -cpp -DTYPE=1 scaling_mod.f90 scaling_with_send_recv.f90  -I$(CPL_PATH)/include/cpl  -L$(CPL_PATH)/lib  -Wl,-rpath=$(CPL_PATH)/lib/ -lcpl  -o ./md

cfd: scaling_with_send_recv.f90
	ftn -dynamic -cpp -DTYPE=0 scaling_mod.f90 scaling_with_send_recv.f90  -I$(CPL_PATH)/include/cpl  -L$(CPL_PATH)/lib  -Wl,-rpath=$(CPL_PATH)/lib/ -lcpl  -o ./cfd

clean:
	rm -f md cfd a.out *.mod *~
