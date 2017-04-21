#ifdef FIX_CLASS

FixStyle(cpl/init, fixCPLInit)

#else

#ifndef LMP_FIX_CPL_INIT_H
#define LMP_FIX_CPL_INIT_H

#include "fix.h"
#include "CPL.h"
#include "CPLSocketLAMMPS.h"
#include<memory>

class fixCPLInit : public LAMMPS_NS::Fix {

public:

    fixCPLInit(class LAMMPS_NS::LAMMPS *lammps, int narg, char **arg);
    int setmask();
    void setup (int vflag); 
    void end_of_step();
    CPLSocketLAMMPS cplsocket;

private:

};

#endif
#endif
