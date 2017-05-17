#include "fix_cpl_init.h"
#include "CPL.h"
#include<iostream>

using namespace LAMMPS_NS;

fixCPLInit::fixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    : Fix (lammps, narg, arg)
{
    class LAMMPS *lmp=lammps;
    cplsocket.initComms();
    cplsocket.initMD(lmp);

}

int fixCPLInit::setmask() {
  int mask = 0;
  mask |= LAMMPS_NS::FixConst::END_OF_STEP;
  return mask;
}


void fixCPLInit::init()
{
   cplsocket.setupFixMDtoCFD(lmp);
   cplsocket.setupFixCFDtoMD(lmp);

    //This is needed to prevent seg fault from
    //use of unallocated pointer
    cplsocket.unpackStress(lmp);

}


void fixCPLInit::setup(int vflag)
{
   end_of_step();

}



void fixCPLInit::end_of_step()
{

    // Communications
    cplsocket.recvStress();
    cplsocket.unpackStress(lmp);
    cplsocket.packVelocity(lmp);
    cplsocket.sendVelocity();

}


