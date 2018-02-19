/*

    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________
     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________
      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________
       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________
        _\/\\\_____________\/\\\/////////____\/\\\_____________
         _\//\\\____________\/\\\_____________\/\\\_____________
          __\///\\\__________\/\\\_____________\/\\\_____________
           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_
            _______\/////////__\///______________\///////////////__
                         C P L  -  L I B R A R Y

           Copyright (C) 2012-2015 Edward Smith & David Trevelyan

License

    This file is part of CPL-Library.

    CPL-Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CPL-Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CPL-Library.  If not, see <http://www.gnu.org/licenses/>.

Description

    "Initialiser fix" for coupled simulation with CPL-Library.

Author(s)

    Edward Smith, Eduardo Ramos Fernandez

*/



#include "fix_cpl_init.h"
#include<iostream>
#include <string.h>
#include "error.h"
#include <stdlib.h>

#include<iostream>
#include "fix_cpl_init.h"
#include "update.h"

fixCPLInit::fixCPLInit(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    : Fix (lammps, narg, arg)
{
    class LAMMPS_NS::LAMMPS *lmp=lammps;
    cplsocket.initComms();
    cplsocket.initMD(lmp);
    nevery = 20;//cplsocket.timestep_ratio;

    forcetype = std::make_shared<std::string>("Undefined");
    sendtype = std::make_shared<std::string>("velocity");
    for (int iarg=0; iarg<narg; iarg+=1){
        std::string arguments(arg[iarg]);
        if (arguments == "forcetype")
            if (iarg+1<narg)
                forcetype = std::make_shared<std::string>(arg[iarg+1]);

        if (arguments == "sendtype")
            if (iarg+1<narg)
                sendtype = std::make_shared<std::string>(arg[iarg+1]);

    }

    std::string forceType(*forcetype);
    if (forceType.compare("Undefined") == 0){
        lammps->error->all(FLERR,"Must specify forcetype on cpl/init line in LAMMPS input file");
    }  

}

int fixCPLInit::setmask() {
  int mask = 0;
  mask |= LAMMPS_NS::FixConst::POST_FORCE;
  return mask;
}


void fixCPLInit::init()
{

    cplsocket.setupFixMDtoCFD(lmp);
    cplsocket.setupFixCFDtoMD(lmp, forcetype);

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

    //std::cout << "fixCPLInit " << update->ntimestep << " " << update->ntimestep%nevery << std::endl;   

    // Recieve and unpack from CFD
    if (update->ntimestep%nevery == 0)
        cplsocket.receive();
	cplsocket.cplfix->post_force(vflag);

    //Pack and send to CFD
    std::string sendType(*sendtype);
    if (sendType.compare("velocity") == 0){
        cplsocket.packVelocity(lmp);
    } else if (sendType.compare("gran") == 0) {
        //cplsocket.packGran(lmp);
        cplsocket.pack(lmp, cplsocket.FORCE | 
                            cplsocket.VOIDRATIO);
    } else if (sendType.compare("granfull") == 0) {
        cplsocket.pack(lmp, cplsocket.VEL | 
                            cplsocket.FORCE |
                            cplsocket.FORCECOEFF | 
                            cplsocket.VOIDRATIO);
    }

    if (update->ntimestep%nevery == 0){
        cplsocket.send();
    }
}    

fixCPLInit::~fixCPLInit(){
	cplsocket.finalizeComms();
}


