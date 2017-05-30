#include "fix_cpl_init.h"
#include<iostream>


fixCPLInit::fixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    : Fix (lammps, narg, arg)
{
   class LAMMPS_NS::LAMMPS *lmp=lammps;
   std::cout << "s1" << cplsocket.nsteps <<std::endl;
    cplsocket.initComms();
   std::cout << "s2" << cplsocket.nsteps <<std::endl;
    cplsocket.initMD(lmp);
   std::cout << "s3" << cplsocket.nsteps <<std::endl;
   nevery = cplsocket.timestep_ratio;

}

int fixCPLInit::setmask() 
{
  int mask = 0;
  //mask |= LAMMPS_NS::FixConst::END_OF_STEP;
  mask |= LAMMPS_NS::FixConst::POST_FORCE;
  return mask;
}


void fixCPLInit::init()
{
   cplsocket.setupFixMDtoCFD(lmp);
   cplsocket.setupFixCFDtoMD(lmp);
   std::cout << "LAMMPS nsteps: " << cplsocket.nsteps <<std::endl;

    //This is needed to prevent seg fault from
    //use of unallocated pointer
    cplsocket.unpackStress(lmp);

}


void fixCPLInit::setup(int vflag)
{

  	post_force(vflag);
  // end_of_step();

}



//void fixCPLInit::end_of_step()
void fixCPLInit::post_force(int vflag)
{
	static int c = 0;
    // Communications
	std::cout << "SEND: " <<  c << std::endl;
	c++;
    cplsocket.recvStress();
	cplsocket.cplfix->end_of_step();

	cplsocket.packVelocity(lmp);
    cplsocket.sendVelocity();
}

fixCPLInit::~fixCPLInit(){
	cplsocket.finalizeComms();
}

