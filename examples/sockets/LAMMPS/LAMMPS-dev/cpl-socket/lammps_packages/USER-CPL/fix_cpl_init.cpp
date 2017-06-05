#include<iostream>
#include "fix_cpl_init.h"

fixCPLInit::fixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    : Fix (lammps, narg, arg)
{
   class LAMMPS_NS::LAMMPS *lmp=lammps;
    cplsocket.initComms();
    cplsocket.initMD(lmp);
   //nevery = cplsocket.timestep_ratio;

}

int fixCPLInit::setmask() {
  int mask = 0;
  mask |= LAMMPS_NS::FixConst::POST_FORCE;
  return mask;
}


void fixCPLInit::init()
{
   cplsocket.setupFixMDtoCFD(lmp);
   cplsocket.setupFixCFDtoMD(lmp);

    //This is needed to prevent seg fault from
    //use of unallocated pointer
    cplsocket.unpackStress(lmp);

    //cplsocket.cplfix->updateStress()
}


void fixCPLInit::setup(int vflag)
{

  	post_force(vflag);

}



void fixCPLInit::post_force(int vflag)
{
    // Communications
    cplsocket.recvStress();
	cplsocket.cplfix->post_force(vflag);

	cplsocket.packVelocity(lmp);
    cplsocket.sendVelocity();
}

fixCPLInit::~fixCPLInit(){
	cplsocket.finalizeComms();
}

