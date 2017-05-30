#ifdef FIX_CLASS

FixStyle(cpl/init, fixCPLInit)

#else

#ifndef LMP_FIX_CPL_INIT_H
#define LMP_FIX_CPL_INIT_H

#include "fix.h"
#include "cpl/cpl.h"
#include "CPLSocketLAMMPS.h"
#include<memory>

class fixCPLInit : public LAMMPS_NS::Fix {

public:

    fixCPLInit(class LAMMPS_NS::LAMMPS *lammps, int narg, char **arg);
	~fixCPLInit();
    int setmask();
    void setup (int vflag); 
    //void end_of_step();
	void post_force(int vflag);
	void init();
    CPLSocketLAMMPS cplsocket;

private:

};

#endif
#endif
