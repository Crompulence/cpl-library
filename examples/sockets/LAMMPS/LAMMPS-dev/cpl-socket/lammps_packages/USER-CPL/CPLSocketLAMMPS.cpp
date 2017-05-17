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
#include "CPLSocketLAMMPS.h"
#include "update.h"
#include "modify.h"
#include "fix_ave_chunk.h"
#include "domain.h"
#include "universe.h"
#include "input.h"
#include "comm.h"
#include "CPL.h"
#include "error.h"

void CPLSocketLAMMPS::initComms() {

    // Split MPI_COMM_WORLD into realm communicators
    CPL::init(CPL::md_realm, realmComm);
    MPI_Comm_rank(realmComm, &rankRealm);

};

void CPLSocketLAMMPS::finalizeComms() {

    CPL::finalize();
    MPI_Finalize();
};

void CPLSocketLAMMPS::initMD(LAMMPS_NS::LAMMPS *lammps) {

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
    //NOTE: Make sure set_timing is called before setup_cfd due to a unfixed bug
    CPL::set_timing(initialstep, 100, dt);
    CPL::setup_md (icomm_grid, globaldomain, xyz_orig);

    // Setup
    timestep_ratio = CPL::get<int> ("timestep_ratio");
    getCellTopology();
    allocateBuffers();


}

void CPLSocketLAMMPS::getCellTopology() {

    // Cell sizes
    dx = CPL::get<double> ("dx");
    dy = CPL::get<double> ("dy");
    dz = CPL::get<double> ("dz");
   
    // Averaging region height below and above the boundary condition plane 
	// TODO: THis is for testing!!!!
    VELBC_BELOW = 0.0;//dy/2.0;
    VELBC_ABOVE = dy;///2.0;

    // Cell bounds for the overlap region
    CPL::get_bnry_limits(velBCRegion.data());
    
    // Cell bounds for velocity BCs region
    //velBCRegion = olapRegion;
    //velBCRegion[3] = velBCRegion[2];
    CPL::my_proc_portion (velBCRegion.data(), velBCPortion.data());

    CPL::get_no_cells(velBCPortion.data(), velBCCells);

    // Cell bounds for the constrained region
    CPL::get_cnst_limits(cnstFRegion.data());
    CPL::my_proc_portion (cnstFRegion.data(), cnstFPortion.data());

    CPL::get_no_cells(cnstFPortion.data(), cnstFCells);

}


void CPLSocketLAMMPS::allocateBuffers() {
    
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



void CPLSocketLAMMPS::setupFixMDtoCFD(LAMMPS_NS::LAMMPS *lammps) {

    double botLeft[3];
    CPL::map_cell2coord(velBCRegion[0] , velBCRegion[2], velBCRegion[4], botLeft);
    botLeft[1] -= VELBC_BELOW;

    double topRight[3];
    CPL::map_cell2coord(velBCRegion[1] , velBCRegion[3], velBCRegion[5], topRight);
    topRight[0] += dx;
    topRight[1] += VELBC_ABOVE;
    topRight[2] += dz;

    // Cell sizes
    std::cout << "Cell size = " << dx <<  " " << dy << " " << dz << std::endl;

    //////////////////////////////////////////
    //This is the code to try to set region
    //////////////////////////////////////////
    int ret;
    char topRight0str[20], topRight1str[20], topRight2str[20];
    char botLeft0str[20], botLeft1str[20], botLeft2str[20];
    ret = sprintf(topRight0str, "%f6", topRight[0]);
    ret = sprintf(topRight1str, "%f6", topRight[1]);
    ret = sprintf(topRight2str, "%f6", topRight[2]);
    ret = sprintf(botLeft0str, "%f6", botLeft[0]);
    ret = sprintf(botLeft1str, "%f6", botLeft[1]);
    ret = sprintf(botLeft2str, "%f6", botLeft[2]);

    char **regionarg = new char*[10];
    regionarg[0] = (char *) "cfdbcregion";
    regionarg[1] = (char *) "block";
    regionarg[2] = (char *) botLeft0str;
    regionarg[3] = (char *) topRight0str;
    regionarg[4] = (char *) botLeft1str;
    regionarg[5] = (char *) topRight1str;
    regionarg[6] = (char *) botLeft2str;
    regionarg[7] = (char *) topRight2str;
    regionarg[8] = (char *) "units";
    regionarg[9] = (char *) "box";
    lammps->domain->add_region(10, regionarg);
    delete [] regionarg;

    // CFD BC region 
//    std::string cmd = "region cfdbcregion block ";
//    cmd += std::to_string(botLeft[0]) + " ";
//    cmd += std::to_string(topRight[0]) + " ";
//    cmd += std::to_string(botLeft[1]) + " ";
//    cmd += std::to_string(topRight[1]) + " ";
//    cmd += std::to_string(botLeft[2]) + " ";
//    cmd += std::to_string(topRight[2]) + " ";
//    cmd += "units box";
//    std::cout << cmd << std::endl;
    //lammps->input->one(cmd.c_str());
    int iregion = lammps->domain->find_region("cfdbcregion");
    if (iregion < 0) lammps->error->all(FLERR,"Fix ID for iregion cfdbcregion does not exist");
    cfdbcregion = lammps->domain->regions[iregion];

    std::cout << "setupFixMDtoCFD Region " << iregion << " " << cfdbcregion->dynamic_check()
              << " " << cfdbcregion->extent_xhi << " " 
              << " " << cfdbcregion->extent_yhi << " " 
              << " " << cfdbcregion->extent_zhi << " " 
                     << cfdbcregion->varshape << std::endl;


    //////////////////////////////////////////
    //This is the code that would set the compute
    //////////////////////////////////////////
    char dxstr[20], dystr[20], dzstr[20];
    ret = sprintf(dxstr, "%f6", dx);
    ret = sprintf(dystr, "%f6", dy);
    ret = sprintf(dzstr, "%f6", dz);

    // CFD BC compute chunk 3d bins in y slice
    char **computearg = new char*[19];
    computearg[0] = (char *) "cfdbccompute";
    computearg[1] = (char *) "all";
    computearg[2] = (char *) "chunk/atom";
    computearg[3] = (char *) "bin/3d";
    computearg[4] = (char *) "x";
    computearg[5] = (char *) "lower";
    computearg[6] = (char *) dxstr;
    computearg[7] = (char *) "y";
    computearg[8] = (char *) "lower";
    computearg[9] = (char *) dystr;
    computearg[10] = (char *) "z";
    computearg[11] = (char *) "lower";
    computearg[12] = (char *) dzstr;
    computearg[13] = (char *) "ids";
    computearg[14] = (char *) "every";
    computearg[15] = (char *) "region";
    computearg[16] = (char *) "cfdbcregion";
    computearg[17] = (char *) "units";
    computearg[18] = (char *) "box";
    lammps->modify->add_compute(19, computearg);
    //Get handle for compute
    int icompute = lammps->modify->find_compute("cfdbccompute");
    if (icompute < 0) lammps->error->all(FLERR,"Fix ID for compute cfdbccompute does not exist");
    cfdbccompute = lammps->modify->compute[icompute];
    delete [] computearg;

//    std::cout << "Region " << iregion 
//              << " " << cfdbccompute->extent_xhi << " " 
//              << " " << cfdbccompute->extent_yhi << " " 
//              << " " << cfdbccompute->extent_zhi << std::endl;

//    cmd = "compute cfdbccompute all chunk/atom bin/3d";
//    cmd += " x lower " + std::to_string(dx);
//    cmd += " y lower " + std::to_string(dy);
//    cmd += " z lower " + std::to_string(dz);
//    cmd += " ids every";
//    cmd += " region cfdbcregion"; 
//    cmd += " units box";
//    std::cout << cmd << std::endl;
    //lammps->input->one(cmd.c_str());

//    int Nfreq = 1; //timestep_ratio;
//    int Nrepeat = 1;
//    int Nevery = 1;
//    cmd = "fix cfdbcfix all ave/chunk ";
//    cmd += std::to_string(Nevery) + " "; 
//    cmd += std::to_string(Nrepeat) + " "; 
//    cmd += std::to_string(Nfreq) + " "; 
//    cmd += "cfdbccompute vx vy vz ";
//    cmd += "norm all";
//    std::cout << cmd << std::endl;
//    lammps->input->one(cmd.c_str());


    //////////////////////////////////////////
    //This is the code that would set the fix
    //////////////////////////////////////////
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
    int Nfreq = 1; //timestep_ratio;
    int Nrepeat = 1;
    int Nevery = 1;

    char Neverystr[20], Nrepeatstr[20], Nfreqstr[20];
    ret = sprintf(Neverystr, "%d", Nevery);
    ret = sprintf(Nrepeatstr, "%d", Nrepeat);
    ret = sprintf(Nfreqstr, "%d", Nfreq);

    char **fixarg = new char*[12];
    fixarg[0] = (char *) "cfdbcfix";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "ave/chunk";
    fixarg[3] = (char *) Neverystr;
    fixarg[4] = (char *) Nrepeatstr; 
    fixarg[5] = (char *) Nfreqstr; 
    fixarg[6] = (char *) "cfdbccompute";
    fixarg[7] = (char *) "vx";
    fixarg[8] = (char *) "vy";
    fixarg[9] = (char *) "vz";
    fixarg[10] = (char *) "norm";
    fixarg[11] = (char *) "all";
    lammps->modify->add_fix(12, fixarg);
    delete [] fixarg;

    //~ Set pointers for this newly-created fix
    int ifix = lammps->modify->find_fix("cfdbcfix");
    if (ifix < 0) lammps->error->all(FLERR,"Fix ID for fix cfdbcfix does not exist");
    cfdbcfix = lammps->modify->fix[ifix];


//    std::cout << "Fix " << ifix 
//              << " " << cfdbcfix->extent_xhi << " " 
//              << " " << cfdbcfix->extent_yhi << " " 
//              << " " << cfdbcfix->extent_zhi << std::endl;

    //cfdbcfix->init();
    // Work around what SEEMS to be a LAMMPS bug in allocation
    // of fix ave/chunk's internal data during construction. This forces
    // allocation.
    //lammps->input->one ("run 0");


    //int fixindex = lammps->modify->find_fix("cfdbcfix");
    //lammps->modify->fix[fixindex]->init();
    
};

void CPLSocketLAMMPS::setupFixCFDtoMD(LAMMPS_NS::LAMMPS *lammps) {

    double botLeft[3];
    CPL::map_cell2coord(cnstFRegion[0] , cnstFRegion[2], cnstFRegion[4], botLeft);

    double topRight[3];
    CPL::map_cell2coord(cnstFRegion[1], cnstFRegion[3], cnstFRegion[5], topRight);
    topRight[0] += dx;
    topRight[1] += dy;
    topRight[2] += dz;

    // Tell LAMMPS to keep track of atoms in constrained region
    int ret;
    char topRight0str[20], topRight1str[20], topRight2str[20];
    char botLeft0str[20], botLeft1str[20], botLeft2str[20];
    ret = sprintf(topRight0str, "%f6", topRight[0]);
    ret = sprintf(topRight1str, "%f6", topRight[1]);
    ret = sprintf(topRight2str, "%f6", topRight[2]);
    ret = sprintf(botLeft0str, "%f6", botLeft[0]);
    ret = sprintf(botLeft1str, "%f6", botLeft[1]);
    ret = sprintf(botLeft2str, "%f6", botLeft[2]);

    char **regionarg = new char*[10];
    regionarg[0] = (char *) "cplforceregion";
    regionarg[1] = (char *) "block";
    regionarg[2] = (char *) botLeft0str;
    regionarg[3] = (char *) topRight0str;
    regionarg[4] = (char *) botLeft1str;
    regionarg[5] = (char *) topRight1str;
    regionarg[6] = (char *) botLeft2str;
    regionarg[7] = (char *) topRight2str;
    regionarg[8] = (char *) "units";
    regionarg[9] = (char *) "box";
    lammps->domain->add_region(10, regionarg);
    delete [] regionarg;

    int iregion = lammps->domain->find_region("cplforceregion");
    if (iregion < 0) lammps->error->all(FLERR,"Fix ID for iregion cplforceregion does not exist");
    cplforceregion = lammps->domain->regions[iregion];

    std::cout << "setupFixCFDtoMD Region " << iregion << " " << cplforceregion->dynamic_check()
              << " " << cplforceregion->extent_xhi << " " 
              << " " << cplforceregion->extent_yhi << " " 
              << " " << cplforceregion->extent_zhi << " " 
                     << cplforceregion->varshape << std::endl;

//    std::string cmd = "region cplforceregion block ";
//    cmd += std::to_string(botLeft[0]) + " " + std::to_string(topRight[0]) + " ";
//    cmd += std::to_string(botLeft[1]) + " " + std::to_string(topRight[1]) + " ";
//    cmd += std::to_string(botLeft[2]) + " " + std::to_string(topRight[2]) + " ";
//    cmd += "units box"; 
//    std::cout << cmd << std::endl;
    //lammps->input->one (cmd.c_str());

    // CFD BC compute chunk 3d bins in y slice
    std::string cmd = "group cplforcegroup dynamic all region cplforceregion every 1";
    std::cout << cmd << std::endl;
    lammps->input->one (cmd.c_str());

//    char **grouparg = new char*[7];
//    grouparg[0] = (char *) "cplforcegroup";
//    grouparg[1] = (char *) "dynamic";
//    grouparg[2] = (char *) "all";
//    grouparg[3] = (char *) "region";
//    grouparg[4] = (char *) "cplforceregion";
//    grouparg[5] = (char *) "every";
//    grouparg[6] = (char *) "1";
//    lammps->input->group->assign(7, grouparg);
//    //Get handle for compute
//    int igroup = lammps->modify->group->find("cplforcegroup");
//    if (igroup < 0) lammps->error->all(FLERR,"Fix ID for group cplforcegroup does not exist");
//    cplforcegroup = lammps->modify->group[igroup];
//    delete [] grouparg;

    // Create a FixCPLForce instance
    //cmd = "fix cplforcefix all cpl/force region cplforceregion";
    //std::cout << cmd << std::endl;
    //lammps->input->one (cmd.c_str());
    char **fixarg = new char*[5];
    fixarg[0] = (char *) "cplforcefix";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "cpl/force";
    fixarg[3] = (char *) "region";
    fixarg[4] = (char *) "cplforceregion"; 
    lammps->modify->add_fix(5, fixarg);
    delete [] fixarg;

    int ifix = lammps->modify->find_fix("cplforcefix");
    if (ifix < 0) lammps->error->all(FLERR,"Fix ID for fix cplforcefix does not exist");

    //Copy coupling fix object??
    cplfix = dynamic_cast<FixCPLForce*>(lammps->modify->fix[ifix]);
    cplfix->updateProcPortion (cnstFPortion.data());


}

// TODO develop a custom fix so that lammps doesn't need to do 
// a global reduce (d.trevelyan@ic.ac.uk) ?
void CPLSocketLAMMPS::packVelocity (const LAMMPS_NS::LAMMPS *lammps) {

    //std::cout << "packVelocity\n" << std::endl;

    //auto fixindex = lammps->modify->find_fix("cfdbcfix");
    //auto fix = lammps->modify->fix[fixindex]; 
   
    int *npxyz_md = lammps->comm->procgrid;
    int nx = velBCCells[0];
    int ny = velBCCells[1];
    int nz = velBCCells[2];

    int row = -1;
    //std::cout << nx << " " << ny << " " << nz << "\n" << std::endl;
    //std::cout << "test cfdbcfix " << cfdbcfix->style << " " << cfdbcfix->compute_array(0, 0) << std::endl;
    for (int i = 0; i != nx; ++i)
    {
        for (int j = 0; j != ny; ++j)
        {
		    for (int k = 0; k != nz; ++k)
		    {
            	int icoord = myCoords[0] * nx + i;
                int jcoord = myCoords[1]*ny + j;
            	int kcoord = myCoords[2] * nz + k;
                row += 1;// icoord*ny*nz*npxyz_md[2] + jcoord*nz*npxyz_md[1] + kcoord;

                // v v v Not needed v v v 
                double x = cfdbcfix->compute_array(row, 0);
                double y = cfdbcfix->compute_array(row, 1);  
                double z = cfdbcfix->compute_array(row, 2);  
                double ncount = cfdbcfix->compute_array(row, 3);  
                // ^ ^ ^ Not needed ^ ^ ^

                double vx = cfdbcfix->compute_array(row, 4);  
                double vy = cfdbcfix->compute_array(row, 5);  
                double vz = cfdbcfix->compute_array(row, 6);  

//                std::cout << i << " " << j << " " << k << " " << x << " " << y << " " << z 
//                               << " " << vx << " " << vy << " " << vz << "\n" << std::endl;

                sendVelocityBuff(0, i, j, k) = vx; // vx;
                sendVelocityBuff(1, i, j, k) = vy; // vy;
                sendVelocityBuff(2, i, j, k) = vz; // vz;
                sendVelocityBuff(3, i, j, k) = ncount;// dummydensity;
            }
        }
    }
    
}

void CPLSocketLAMMPS::unpackStress (const LAMMPS_NS::LAMMPS *lammps) {

    shaPtrArrayDoub stress = std::make_shared<arrayDoub> (recvStressBuff);
    cplfix->updateStress (stress);

};

// Sends buffer to overlapping CFD processes.
//void CPLSocketLAMMPS::send() 
//{
//    CPL::send(sendBuf.data(), sendBuf.shapeData(), velBCRegion.data());
//}

//// Receives recvBuf from overlapping CFD processes.
//void CPLSocketLAMMPS::recv()
//{
//    CPL::recv(recvBuf.data(), recvBuf.shapeData(), cnstFRegion.data());
//}


void CPLSocketLAMMPS::sendVelocity() {

    // Send the data to CFD
    CPL::send(sendVelocityBuff.data(), sendVelocityBuff.shapeData(), velBCRegion.data());

};

void CPLSocketLAMMPS::recvStress() {

    // Receive from CFD
    CPL::recv(recvStressBuff.data(), recvStressBuff.shapeData(), cnstFRegion.data());

};

