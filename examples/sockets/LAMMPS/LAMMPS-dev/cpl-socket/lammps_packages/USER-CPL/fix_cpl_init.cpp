#include<iostream>
#include "fix_cpl_init.h"
#include "update.h"

fixCPLInit::fixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    : Fix (lammps, narg, arg)
{
   class LAMMPS_NS::LAMMPS *lmp=lammps;
    cplsocket.initComms();
    cplsocket.initMD(lmp);
    nevery = 50;//cplsocket.timestep_ratio;

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
    cplsocket.unpackBuf(lmp);
    //cplsocket.cplfix->updateBuf()
}


void fixCPLInit::setup(int vflag)
{

  	post_force(vflag);

}



void fixCPLInit::post_force(int vflag)
{
    // Communications
    if (update->ntimestep%nevery == 0)
        cplsocket.receiveBuf();
	cplsocket.cplfix->post_force(vflag);

	cplsocket.packVelocity(lmp);
    if (update->ntimestep%nevery == 0){
        cplsocket.sendVelocity();
    }
}    

fixCPLInit::~fixCPLInit(){
	cplsocket.finalizeComms();
}

