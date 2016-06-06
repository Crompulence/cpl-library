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

    "Socket" class for interfacing with CPL-Library.

Author(s)

    David Trevelyan

*/
#include<iostream>
#include "socket.h"
#include "update.h"
#include "modify.h"
#include "fix_ave_chunk.h"
#include "domain.h"
#include "universe.h"
#include "input.h"
#include "comm.h"
#include "CPL.h"


void SocketLAMMPS::\
initComms (int argc, char **argv) {

    MPI_Init(&argc, &argv);

    // Split MPI_COMM_WORLD into realm communicators
    CPL::init(CPL::md_realm, realmComm);
    MPI_Comm_rank(realmComm, &rankRealm);

};

void SocketLAMMPS::\
finalizeComms() {

    MPI_Finalize();
};

void SocketLAMMPS::\
initMD (LAMMPS_NS::LAMMPS *lammps) {

    // Store my own coordinates for later
    myCoords[0] = lammps->comm->myloc[0];
    myCoords[1] = lammps->comm->myloc[1];
    myCoords[2] = lammps->comm->myloc[2];

    // Parameters for coupler_md_init
    int initialstep = lammps->update->firststep;
    double dt = lammps->update->dt;
    int *npxyz_md = lammps->comm->procgrid;
    double *globaldomain = lammps->domain->prd;
    double dummydensity = -666.0;

    // Set up new cartesian communicator with same coordinates as lammps
    // interal cartesian communicator (based on mycoords)
    MPI_Comm icomm_grid;
    int periods[3] = {0, 0, 0};
    CPL::Cart_create (lammps->world, 3, npxyz_md, periods, myCoords.data(),
                      &icomm_grid);
   
    //TODO: get the origin from LAMMPS 
    double xyz_orig[3] = {0.0 ,0.0, 0.0};
    // Call MD init, this->nsteps will be reset to correct value
    CPL::setup_md (nsteps, initialstep, dt, icomm_grid, globaldomain,
                   xyz_orig, dummydensity);

    // Setup
    timestep_ratio = CPL::get<int> ("timestep_ratio");
    getCellTopology();
    setupFixMDtoCFD(lammps);
    setupFixCFDtoMD(lammps);
    allocateBuffers();

}

void SocketLAMMPS::\
getCellTopology() {

    // Cell sizes
    dx = CPL::get<double> ("dx");
    dy = CPL::get<double> ("dy");
    dz = CPL::get<double> ("dz");
   
    // Averaging region height below and above the boundary condition plane 
    VELBC_BELOW = dy;
    VELBC_ABOVE = 0.0;

    // Cell bounds for the overlap region
    CPL::get_olap_limits(olapRegion.data());
    
    // Cell bounds for velocity BCs region
    velBCRegion = olapRegion;
    velBCRegion[3] = velBCRegion[2];
    CPL::my_proc_portion (velBCRegion.data(), velBCPortion.data());

    CPL::get_no_cells(velBCPortion.data(), velBCCells);

    // Cell bounds for the constrained region
    CPL::get_cnst_limits(cnstFRegion.data());
    CPL::my_proc_portion (cnstFRegion.data(), cnstFPortion.data());

    CPL::get_no_cells(cnstFPortion.data(), cnstFCells);

}


void SocketLAMMPS::\
allocateBuffers() {
    
    // Received stress field
    int zeroShapeStress[4] = {9, 0, 0, 0};
    int recvShape[4] = {9, cnstFCells[0], cnstFCells[1], cnstFCells[2]};
    sendStressBuff.resize(4, zeroShapeStress);
    recvStressBuff.resize(4, recvShape);

    // LAMMPS computed velocity field
    int sendShape[4] = {4, velBCCells[0], velBCCells[1], velBCCells[2]};
    int zeroShapeVel[4] = {4, 0, 0, 0};
    recvVelocityBuff.resize (4, zeroShapeVel);
    sendVelocityBuff.resize (4, sendShape);
}



void SocketLAMMPS::\
setupFixMDtoCFD(LAMMPS_NS::LAMMPS *lammps) {

    double botLeft[3];
    CPL::map_cell2coord(velBCRegion[0] , velBCRegion[2], velBCRegion[4], botLeft);
    botLeft[1] -= VELBC_BELOW;

    double topRight[3];
    CPL::map_cell2coord(velBCRegion[1] , velBCRegion[3], velBCRegion[5], topRight);
    topRight[0] += dx;
    topRight[1] += VELBC_ABOVE;
    topRight[2] += dz;

    // LAMMPS commands to set up the fix

    // CFD BC region 
    std::string cmd = "region cfdbcregion block ";
    cmd += std::to_string (botLeft[0]) + " ";
    cmd += std::to_string (topRight[0]) + " ";
    cmd += std::to_string (botLeft[1]) + " ";
    cmd += std::to_string (topRight[1]) + " ";
    cmd += std::to_string (botLeft[2]) + " ";
    cmd += std::to_string (topRight[2]) + " ";
    cmd += "units box";
    //std::cout << cmd << std::endl;
    lammps->input->one (cmd.c_str());

    // CFD BC compute chunk 2d bins in single y slice
    cmd = "compute cfdbccompute all chunk/atom bin/2d";
    cmd += " x lower " + std::to_string (dx);
    cmd += " z lower " + std::to_string (dz);
    cmd += " ids every";
    cmd += " region cfdbcregion"; 
    cmd += " units box";
    lammps->input->one (cmd.c_str());

    // CFD BC averaging fix 
    // Average values are generated every Nfreq time steps, taken
    // from the average of the Nrepeat preceeding timesteps separated
    // by Nevery. For example, consider 
    // 
    //       Nfreq = 100;   Nrepeat = 6;   Nevery = 2;
    //
    // The average here would be taken over instantaneous snapshots at
    // timesteps 90, 92, 94, 96, 98, and 100 for the first value. The
    // next output would be an average from 190, 192, 194, 196, 198 and
    // 200 (and so on).
    int Nfreq = timestep_ratio;
    int Nrepeat = 1;
    int Nevery = 1;
    cmd = "fix cfdbcfix all ave/chunk ";
    cmd += std::to_string (Nevery) + " "; 
    cmd += std::to_string (Nrepeat) + " "; 
    cmd += std::to_string (Nfreq) + " "; 
    cmd += "cfdbccompute vx vy vz ";
    cmd += "norm all ";
    lammps->input->one (cmd.c_str());

    // Work around what SEEMS to be a LAMMPS bug in allocation
    // of fix ave/chunk's internal data during construction. This forces
    // allocation.
    lammps->input->one ("run 0");
    
};

void SocketLAMMPS::\
setupFixCFDtoMD(LAMMPS_NS::LAMMPS *lammps) {

    double botLeft[3];
    CPL::map_cell2coord(cnstFRegion[0] , cnstFRegion[2], cnstFRegion[4], botLeft);

    double topRight[3];
    CPL::map_cell2coord(cnstFRegion[1], cnstFRegion[3], cnstFRegion[5], topRight);
    topRight[0] += dx;
    topRight[1] += dy;
    topRight[2] += dz;

    // Tell LAMMPS to keep track of atoms in constrained region
    std::string cmd = "region cplforceregion block ";
    cmd += std::to_string(botLeft[0]) + " " + std::to_string(topRight[0]) + " ";
    cmd += std::to_string(botLeft[1]) + " " + std::to_string(topRight[1]) + " ";
    cmd += std::to_string(botLeft[2]) + " " + std::to_string(topRight[2]) + " ";
    cmd += "units box"; 
    //std::cout << cmd << std::endl;
    lammps->input->one (cmd.c_str());
    // TODO (d.trevelyan@ic.ac.uk) every 1? , Maybe just every ratio timestep
    cmd = "group cplforcegroup dynamic all region cplforceregion every 1";
    lammps->input->one (cmd.c_str());

    // Create a FixCPLForce instance
    cmd = "fix cplforcefix all cpl/force region cplforceregion";
    lammps->input->one (cmd.c_str());
    int fixindex = lammps->modify->find_fix ("cplforcefix");
    cplfix = dynamic_cast<FixCPLForce*> 
                  (lammps->modify->fix[fixindex]);
    cplfix->updateProcPortion (cnstFPortion.data());
}

// TODO develop a custom fix so that lammps doesn't need to do 
// a global reduce (d.trevelyan@ic.ac.uk) ?
void SocketLAMMPS::\
packVelocity (const LAMMPS_NS::LAMMPS *lammps) {

    auto fixindex = lammps->modify->find_fix("cfdbcfix");
    auto fix = lammps->modify->fix[fixindex]; 
   
    int *npxyz_md = lammps->comm->procgrid;
    int nx = velBCCells[0];
    int nz = velBCCells[2];

    for (int i = 0; i != nx; ++i)
    {
        for (int k = 0; k != nz; ++k)
        {
            int icoord = myCoords[0] * nx + i;
            int kcoord = myCoords[2] * nz + k;
            int row = icoord * nz * npxyz_md[2] + kcoord;
            // Not needed - BEGIN
            double x = fix->compute_array(row, 0);
            double z = fix->compute_array(row, 1);  
            double ncount = fix->compute_array(row, 2);  
            // Not needed - END
            double vx = fix->compute_array(row, 3);  
            double vy = fix->compute_array(row, 4);  
            double vz = fix->compute_array(row, 5);  

            sendVelocityBuff(0, i, 0, k) = vx; // vx;
            sendVelocityBuff(1, i, 0, k) = vy; // vy;
            sendVelocityBuff(2, i, 0, k) = vz; // vz;
            sendVelocityBuff(3, i, 0, k) = 1.0;// dummydensity;
        }
    }
    
}

void SocketLAMMPS::\
unpackStress (const LAMMPS_NS::LAMMPS *lammps) {

    shaPtrArrayDoub stress = std::make_shared<arrayDoub> (recvStressBuff);
    cplfix->updateStress (stress);

};

void SocketLAMMPS::\
sendVelocity() {

    // Send the data to CFD
    CPL::gather(sendVelocityBuff.data(), sendVelocityBuff.shapeData(),
                velBCRegion.data(), recvVelocityBuff.data(), 
                recvVelocityBuff.shapeData());

};

void SocketLAMMPS::\
recvStress() {

    // Receive from CFD
    CPL::scatter(sendStressBuff.data(), sendStressBuff.shapeData(),
                 cnstFRegion.data(), recvStressBuff.data(),
                 recvStressBuff.shapeData());
};

