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

    David Trevelyan, Edward Smith, Eduardo Fernandez-Ramos

*/
#include<iostream>
#include "update.h"
#include "modify.h"
#include "fix_ave_chunk.h"
#include "domain.h"
#include "universe.h"
#include "input.h"
#include "comm.h"
#include "error.h"
#include "CPLSocketLAMMPS.h"

#include "cpl/CPL_cartCreate.h"

void CPLSocketLAMMPS::initComms() {

    // Split MPI_COMM_WORLD into realm communicators
    CPL::init(CPL::md_realm, realmComm);
    MPI_Comm_rank(realmComm, &rankRealm);

};

void CPLSocketLAMMPS::finalizeComms() {

    CPL::finalize();
    //MPI_Finalize();
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
    //CPL::set_timing(initialstep, 100, dt);
    CPL::setup_md (icomm_grid, globaldomain, xyz_orig);

    // Setup
    timestep_ratio = CPL::get<int> ("timestep_ratio");
    getCellTopology();

    allocateBuffers(lammps, VEL | FORCE | FORCECOEFF | VOIDRATIO);

}

void CPLSocketLAMMPS::getCellTopology() {

    // Cell sizes
    dx = CPL::get<double> ("dx");
    dy = CPL::get<double> ("dy");
    dz = CPL::get<double> ("dz");
   
    // Averaging region height below and above the boundary condition plane 
    VELBC_BELOW = dy;
    VELBC_ABOVE = 0.0;

    // Cell bounds for velocity BCs region
    CPL::get_bnry_limits(velBCRegion.data());
    CPL::my_proc_portion (velBCRegion.data(), velBCPortion.data());
    CPL::get_no_cells(velBCPortion.data(), velBCCells);

    // Cell bounds for the constrained region
    CPL::get_cnst_limits(cnstFRegion.data());
    CPL::my_proc_portion (cnstFRegion.data(), cnstFPortion.data());
    CPL::get_no_cells(cnstFPortion.data(), cnstFCells);

}


//void CPLSocketLAMMPS::allocateBuffers() {
//    
//    // Received Buf field
//    int recvShape[4] = {9, cnstFCells[0], cnstFCells[1], cnstFCells[2]};
//    recvBuf.resize(4, recvShape);

//    // LAMMPS computed velocity field
//    int sendShape[4] = {4, velBCCells[0], velBCCells[1], velBCCells[2]};
//    sendBuf.resize (4, sendShape);

//}



void CPLSocketLAMMPS::setupFixMDtoCFD(LAMMPS_NS::LAMMPS *lammps)
{
    double botLeft[3];
    CPL::map_cell2coord(velBCRegion[0] , velBCRegion[2], velBCRegion[4], botLeft);

    double topRight[3];
    CPL::map_cell2coord(velBCRegion[1] , velBCRegion[3], velBCRegion[5], topRight);
    topRight[0] += dx;
    topRight[1] += dy;
    topRight[2] += dz;

    //////////////////////////////////////////
    //This is the code sets the region
    //////////////////////////////////////////
    int ret;
    char topRight0str[20], topRight1str[20], topRight2str[20];
    char botLeft0str[20], botLeft1str[20], botLeft2str[20];
    ret = sprintf(topRight0str, "%f", topRight[0]);
    ret = sprintf(topRight1str, "%f", topRight[1]);
    ret = sprintf(topRight2str, "%f", topRight[2]);
    ret = sprintf(botLeft0str, "%f", botLeft[0]);
    ret = sprintf(botLeft1str, "%f", botLeft[1]);
    ret = sprintf(botLeft2str, "%f", botLeft[2]);

    // CFD BC region 
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

    int iregion = lammps->domain->find_region("cfdbcregion");
    if (iregion < 0) lammps->error->all(FLERR,"Fix ID for iregion cfdbcregion does not exist");
    cfdbcregion = lammps->domain->regions[iregion];

    //////////////////////////////////////////
    //This is the code sets the compute
    //////////////////////////////////////////
    char dxstr[20], dystr[20], dzstr[20], low_y[20], hi_y[20];
    ret = sprintf(dxstr, "%f", dx);
    ret = sprintf(dystr, "%f", dy);
    ret = sprintf(dzstr, "%f", dz);
    ret = sprintf(low_y, "%f", 200.0);
    ret = sprintf(hi_y, "%f", 220.0);

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

    //////////////////////////////////////////
    //This code sets the fix
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
    //fixarg[12] = (char *) "file";
    //fixarg[13] = (char *) "cplchunk";
    lammps->modify->add_fix(12, fixarg);
    delete [] fixarg;

    //~ Set pointers for this newly-created fix
    int ifix = lammps->modify->find_fix("cfdbcfix");
    if (ifix < 0) lammps->error->all(FLERR,"Fix ID for fix cfdbcfix does not exist");
    cfdbcfix = lammps->modify->fix[ifix];
   
};

void CPLSocketLAMMPS::setupFixCFDtoMD(LAMMPS_NS::LAMMPS *lammps, 
                                      std::shared_ptr<std::string> forcetype) {

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
    char botLeft0str[20],  botLeft1str[20],  botLeft2str[20];
    //char topRight0str[sizeof(topRight[0])], topRight1str[sizeof(topRight[1])], topRight2str[sizeof(topRight[2])];
    //char botLeft0str[sizeof(botLeft[0])],  botLeft1str[sizeof(botLeft[1])],  botLeft2str[sizeof(botLeft[2])];

    ret = sprintf(topRight0str, "%f", topRight[0]);
    ret = sprintf(topRight1str, "%f", topRight[1]);
    ret = sprintf(topRight2str, "%f", topRight[2]);
    ret = sprintf(botLeft0str, "%f", botLeft[0]);
    ret = sprintf(botLeft1str, "%f", botLeft[1]);
    ret = sprintf(botLeft2str, "%f", botLeft[2]);

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

    // CFD BC compute chunk 3d bins in y slice
    std::string cmd = "group cplforcegroup dynamic all region cplforceregion every 1";
    std::cout << cmd << std::endl;
    lammps->input->one (cmd.c_str());

    // Create a FixCPLForce instance
    std::string str = *forcetype;
    char * writable = new char[str.size() + 1];
    std::copy(str.begin(), str.end(), writable);
    writable[str.size()] = '\0'; // terminating 0
    char **fixarg = new char*[7];
    fixarg[0] = (char *) "cplforcefix";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "cpl/force";
    fixarg[3] = (char *) "region";
    fixarg[4] = (char *) "cplforceregion";
    fixarg[5] = (char *) "forcetype";
    fixarg[6] = writable;
    lammps->modify->add_fix(7, fixarg);
    delete writable;
    delete [] fixarg;

    int ifix = lammps->modify->find_fix("cplforcefix");
    if (ifix < 0) lammps->error->all(FLERR,"Fix ID for fix cplforcefix does not exist");

    //Upcast Fix to child class FixCPLForce
    cplfix = dynamic_cast<FixCPLForce*>(lammps->modify->fix[ifix]);
    cplfix->updateProcPortion (cnstFPortion.data());

}

// TODO develop a custom fix so that lammps doesn't need to do 
// a global reduce (d.trevelyan@ic.ac.uk) ?
void CPLSocketLAMMPS::packVelocity(const LAMMPS_NS::LAMMPS *lammps) {

    // For some reason, this is not called as it should be as part of loop.
    // Explictly called here to ensure velocity ready to be packed.
    cfdbcfix->end_of_step();
    //std::cout << "cfdbcfix->end_of_step called in CPLSocketLAMMPS::packVelocity but nevery=" << cfdbcfix->nevery << std::endl;   

	if (CPL::is_proc_inside(velBCPortion.data())) 
    {
        int *npxyz_md = lammps->comm->procgrid;
	    int nc_velBCRegion[3];
        CPL::get_no_cells(velBCRegion.data(), nc_velBCRegion);
        int row;
	    int glob_cell[3], loc_cell[3];
        //Chosen arbitarily for now
        double V = 0.01;           
        double Vcell = 1.;
        //Chosen arbitarily for now

        for (int i = velBCPortion[0]; i <= velBCPortion[1]; i++)
        {
            for (int j = velBCPortion[2]; j <= velBCPortion[3]; j++)
            {
		        for (int k = velBCPortion[4]; k <= velBCPortion[5]; k++)
		        {
                    row = i*nc_velBCRegion[1]*nc_velBCRegion[2] + j*nc_velBCRegion[2] + k;
				     
                    double x = cfdbcfix->compute_array(row, 0);
                    double y = cfdbcfix->compute_array(row, 1);  
                    double z = cfdbcfix->compute_array(row, 2);  
                    double ncount = cfdbcfix->compute_array(row, 3);  

                    double vx = cfdbcfix->compute_array(row, 4);  
                    double vy = cfdbcfix->compute_array(row, 5);  
                    double vz = cfdbcfix->compute_array(row, 6);  

				    glob_cell[0] = i; glob_cell[1] = j; glob_cell[2] = k;
				    CPL::map_glob2loc_cell(velBCPortion.data(), glob_cell, loc_cell);

                    sendBuf(0, loc_cell[0], loc_cell[1], loc_cell[2]) = vx;
                    sendBuf(1, loc_cell[0], loc_cell[1], loc_cell[2]) = vy;
                    sendBuf(2, loc_cell[0], loc_cell[1], loc_cell[2]) = vz; 
                    sendBuf(3, loc_cell[0], loc_cell[1], loc_cell[2]) = ncount; 
                }
            }
        }
    }
}


//Pack porosity and forces
void CPLSocketLAMMPS::packGran(const LAMMPS_NS::LAMMPS *lammps) {

    // For some reason, this is not called as it should be as part of loop.
    // Explictly called here to ensure velocity ready to be packed.
    //cfdbcfix->end_of_step();
	if (CPL::is_proc_inside(velBCPortion.data())) 
    {
        int *npxyz_md = lammps->comm->procgrid;
	    int nc_velBCRegion[3];
        CPL::get_no_cells(velBCRegion.data(), nc_velBCRegion);
        int row;
	    int glob_cell[3], loc_cell[3];
        double Vcell = dx*dy*dz;

        //Downcast to CPLForceDrag type here
        CPLForceDrag& Granfxyz = dynamic_cast<CPLForceDrag&>(*cplfix->fxyz);

        //Chosen arbitarily for now
        for (int i = velBCPortion[0]; i <= velBCPortion[1]; i++) {
        for (int j = velBCPortion[2]; j <= velBCPortion[3]; j++) {
	    for (int k = velBCPortion[4]; k <= velBCPortion[5]; k++) {			                       

		    glob_cell[0] = i; glob_cell[1] = j; glob_cell[2] = k;
		    CPL::map_glob2loc_cell(velBCPortion.data(), glob_cell, loc_cell);

            sendBuf(0, loc_cell[0], loc_cell[1], loc_cell[2]) = Granfxyz.FSums(0,i,j,k);
            sendBuf(1, loc_cell[0], loc_cell[1], loc_cell[2]) = Granfxyz.FSums(1,i,j,k);
            sendBuf(2, loc_cell[0], loc_cell[1], loc_cell[2]) = Granfxyz.FSums(2,i,j,k);

            if (Granfxyz.eSums(i,j,k)/Vcell > 1.) {
//                std::cout << "cell i,j,k = " << i << " " << j << " " << k << " " << 
//                            " Particle Volume = " << Granfxyz.eSums(i,j,k) <<
//                            " Particle No. = " << Granfxyz.nSums(i,j,k) <<
//                            " Cell Volume = " << Vcell << 
//                            " Ratio = " << Granfxyz.eSums(i,j,k)/Vcell << std::endl;
                //lammps->error->all(FLERR,"Error packGran -- Particle Volume > Cell Volume");

                std::cout << "Warning, eps = 0 so set to 0.1 in CPLSocketLAMMPS::packGran" << std::endl;
                sendBuf(3, loc_cell[0], loc_cell[1], loc_cell[2]) = 0.1;
            } else {
                sendBuf(3, loc_cell[0], loc_cell[1], loc_cell[2]) = 1. - Granfxyz.eSums(i,j,k)/Vcell;
            }
            //if (Granfxyz.eSums(i,j,k) != 0.)
            //    std::cout << "eSums " << i << " " << j << " " << k << " " << Granfxyz.eSums(i,j,k) 
            //              << " " << Vcell << " " << sendBuf(3, loc_cell[0], loc_cell[1], loc_cell[2]) << std::endl;
        } } }
    }
}




//Pack general using bitflag
void CPLSocketLAMMPS::allocateBuffers(const LAMMPS_NS::LAMMPS *lammps, int sendtype) {

    // Received Buf field
    int recvShape[4] = {9, cnstFCells[0], cnstFCells[1], cnstFCells[2]};
    recvBuf.resize(4, recvShape);

    //Check what is to be packed and sent
    int packsize=0;
    if ((sendtype & VEL) == VEL){
        packsize += VELSIZE;
    }
    if ((sendtype & NBIN) == NBIN){
        packsize += NBINSIZE;
    }
    if ((sendtype & STRESS) == STRESS){
        packsize += STRESSSIZE;
        lammps->error->all(FLERR," sendtype stress not developed. Aborting.");
    }
    if ((sendtype & FORCE) == FORCE){
        packsize += FORCESIZE;
    }
    if ((sendtype & FORCECOEFF) == FORCECOEFF){
        packsize += FORCECOEFFSIZE;
    }
    if ((sendtype & VOIDRATIO) == VOIDRATIO){
        packsize += VOIDRATIOSIZE;
    }

    // LAMMPS computed velocity field
    int sendShape[4] = {packsize, velBCCells[0], velBCCells[1], velBCCells[2]};
    sendBuf.resize(4, sendShape);

    if (sendtype > 63)
        lammps->error->all(FLERR," sendtype bit flag unknown type. Aborting.");
}


//Pack general using bitflag
void CPLSocketLAMMPS::pack(const LAMMPS_NS::LAMMPS *lammps, int sendtype) {

        int *npxyz_md = lammps->comm->procgrid;
	    int nc_velBCRegion[3];
        CPL::get_no_cells(velBCRegion.data(), nc_velBCRegion);
        int row;
	    int glob_cell[3], loc_cell[3];
        double Vcell = dx*dy*dz;

        //Allocate buffers to send
        allocateBuffers(lammps, sendtype);

        //Downcast to CPLForceDrag type here
        CPLForceDrag& Granfxyz = dynamic_cast<CPLForceDrag&>(*cplfix->fxyz);

        //Chosen arbitarily for now
        for (int i = velBCPortion[0]; i <= velBCPortion[1]; i++) {
        for (int j = velBCPortion[2]; j <= velBCPortion[3]; j++) {
	    for (int k = velBCPortion[4]; k <= velBCPortion[5]; k++) {			                       

		    glob_cell[0] = i; glob_cell[1] = j; glob_cell[2] = k;
		    CPL::map_glob2loc_cell(velBCPortion.data(), glob_cell, loc_cell);
            row = i*nc_velBCRegion[1]*nc_velBCRegion[2] + j*nc_velBCRegion[2] + k;
            int npack = 0;

            //Check what is to be packed and sent
            if ((sendtype & VEL) == VEL){
                double vx = cfdbcfix->compute_array(row, 4);  
                double vy = cfdbcfix->compute_array(row, 5);  
                double vz = cfdbcfix->compute_array(row, 6);  

                sendBuf(npack+0, loc_cell[0], loc_cell[1], loc_cell[2]) = vx;
                sendBuf(npack+1, loc_cell[0], loc_cell[1], loc_cell[2]) = vy;
                sendBuf(npack+2, loc_cell[0], loc_cell[1], loc_cell[2]) = vz; 
                npack += VELSIZE;
            }
            if ((sendtype & NBIN) == NBIN){
                double ncount = cfdbcfix->compute_array(row, 3);
                sendBuf(npack, loc_cell[0], loc_cell[1], loc_cell[2]) = ncount; 
                npack += NBINSIZE;
            }
            if ((sendtype & STRESS) == STRESS){
                lammps->error->all(FLERR," sendtype stress not developed. Aborting.");
                npack += STRESSSIZE;
            }
            if ((sendtype & FORCE) == FORCE){
                sendBuf(npack+0, loc_cell[0], loc_cell[1], loc_cell[2]) = Granfxyz.FSums(0,i,j,k);
                sendBuf(npack+1, loc_cell[0], loc_cell[1], loc_cell[2]) = Granfxyz.FSums(1,i,j,k);
                sendBuf(npack+2, loc_cell[0], loc_cell[1], loc_cell[2]) = Granfxyz.FSums(2,i,j,k);
                npack += FORCESIZE;
            }
            if ((sendtype & FORCECOEFF) == FORCECOEFF){
                sendBuf(npack+0, loc_cell[0], loc_cell[1], loc_cell[2]) = Granfxyz.FcoeffSums(i,j,k);
                npack += FORCECOEFFSIZE;
            }
            if ((sendtype & VOIDRATIO) == VOIDRATIO){
                if (Granfxyz.eSums(i,j,k)/Vcell > 1.) {
                    std::cout << "Warning, eps = 0 so set to 0.1 in CPLSocketLAMMPS::packGran" << std::endl;
                    sendBuf(npack, loc_cell[0], loc_cell[1], loc_cell[2]) = 0.1;
                } else {
                    sendBuf(npack, loc_cell[0], loc_cell[1], loc_cell[2]) = 1. - Granfxyz.eSums(i,j,k)/Vcell;
                }
                npack += VOIDRATIOSIZE;
            }

//            for (int n = 0; n < 8; n++) {
//                if (i == 5 && k == 5) {
//                    std::cout << "PACK " << i << " " << j << " " << k << " " << n << " " 
//                              << loc_cell[0] << " " << loc_cell[1] << " " << loc_cell[2] << " " 
//                              << sendBuf(n, i, j, k) << std::endl;
//                }
//            }

        }}}

}
    
    



void CPLSocketLAMMPS::unpackBuf (const LAMMPS_NS::LAMMPS *lammps) {
    // Unpack buffer
    cplfix->updateBuf(recvBuf);
};


void CPLSocketLAMMPS::send() {

//    for (int i = velBCPortion[0]; i <= velBCPortion[1]; i++) {
//    for (int j = velBCPortion[2]; j <= velBCPortion[3]; j++) {
//    for (int k = velBCPortion[4]; k <= velBCPortion[5]; k++) {		
//    for (int n = 0; n < 8; n++) {
//        if (i == 5 && k == 5) {
//            std::cout << "SEND " << i << " " << j << " " << k << " " << n << " " 
//                      << sendBuf(n, i, j, k) << std::endl;
//        }
//    } } } }

    // Send the data to CFD
    CPL::send(sendBuf.data(), sendBuf.shapeData(), velBCRegion.data());
};

void CPLSocketLAMMPS::receive() {
    // Receive from CFD
    CPL::recv(recvBuf.data(), recvBuf.shapeData(), cnstFRegion.data());

//    for (int i = velBCPortion[0]; i <= velBCPortion[1]; i++) {
//    for (int j = velBCPortion[2]; j <= velBCPortion[3]; j++) {
//    int i = 14; int j = 14; 
//    for (int k = velBCPortion[4]; k <= velBCPortion[5]; k++) {			                       
//        std::cout << "recv " << i << " " << j << " " << k  << " " << 
//                  recvBuf(0,i,j,k) << " " << recvBuf(1,i,j,k) << " " << recvBuf(2,i,j,k) << " " << 
//                  recvBuf(3,i,j,k) << " " << recvBuf(4,i,j,k) << " " << recvBuf(5,i,j,k) << " " << 
//                  recvBuf(6,i,j,k) << " " << recvBuf(7,i,j,k) << " " << recvBuf(8,i,j,k) <<  std::endl;
//    }//}}
};

