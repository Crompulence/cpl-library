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
#include "fix_cpl_force.h"

namespace CPL
{

    SocketLAMMPS::SocketLAMMPS(){};

    void SocketLAMMPS::initComms (int argc, char **argv)
    {

        MPI_Init(&argc, &argv);

        // Split MPI_COMM_WORLD into realm communicators
        CPL::create_comm(CPL::md_realm, realmComm);

        MPI_Comm_rank(realmComm, &rankRealm);

    };

    void SocketLAMMPS::finalizeComms()
    {
        MPI_Finalize();
    };

    void SocketLAMMPS::initMD (LAMMPS_NS::LAMMPS *lammps)
    {

        // Store my own coordinates for later
        myCoords.assign
        (
            std::begin(lammps->comm->myloc),
            std::end(lammps->comm->myloc)
        );

        // Fortran coords count from 1 not 0
        myFCoords = myCoords;
        for (auto &c : myFCoords) c += 1;

        // Parameters for coupler_md_init
        int initialstep = lammps->update->firststep;
        double dt = lammps->update->dt;
        int *npxyz_md = lammps->comm->procgrid;
        double *globaldomain = lammps->domain->prd;
        double dummydensity = -666.0;

        // Establish coordinates of all processors
        int nprocs = lammps->comm->nprocs;
        int tempShape[2] = {3, nprocs};
        ndArray<int> allCoords(2, tempShape);
        for (int i = 0; i != npxyz_md[0]; ++i)
        {
            for (int j = 0; j != npxyz_md[1]; ++j)
            {
                for (int k = 0; k != npxyz_md[2]; ++k)
                {
                    int procno = lammps->comm->grid2proc[i][j][k];
                    allCoords(0, procno) = i;
                    allCoords(1, procno) = j;
                    allCoords(2, procno) = k;
                }
            }
        }

        // Fortran coordinate indices count from 1, not 0
        allFCoords = allCoords + 1;
        
        // Set up new cartesian communicator with same coordinates as lammps
        // interal cartesian communicator (based on mycoords)
        MPI_Comm icomm_grid;
        int periods[3] = {0, 0, 0};
        CPL::Cart_create
        (
            lammps->world,
            3,
            npxyz_md,
            periods,
            myCoords.data(),
            &icomm_grid
        );

        // Call MD init, this->nsteps will be reset to correct value
        CPL::md_init
        (
            this->nsteps,
            initialstep,
            dt,
            icomm_grid,
            allFCoords.data(),
            npxyz_md,
            globaldomain,
            dummydensity
        );

        getTimestepRatio();
        getCellTopology();
        setupFixMDtoCFD(lammps);
        setupFixCFDtoMD(lammps);

    }
 
    void SocketLAMMPS::getTimestepRatio()
    {
        timestep_ratio = CPL::get<int> ("timestep_ratio");
    }

    void SocketLAMMPS::getCellTopology()
    {

        // Minimum and maximum CPL cell indices on this processor
        extents.resize (6);
        CPL::proc_extents (myFCoords.data(), CPL::md_realm, extents.data());

        // Number of CPL cells on this processor in each direction
        ncxl = extents[1] - extents[0] + 1;
        ncyl = extents[3] - extents[2] + 1;
        nczl = extents[5] - extents[4] + 1;

        // Minimum and maximum CPL cell indices of global overlap region
        olap_limits.resize (6);
        olap_limits[0] = CPL::get<int> ("icmin_olap");
        olap_limits[1] = CPL::get<int> ("icmax_olap");
        olap_limits[2] = CPL::get<int> ("jcmin_olap");
        olap_limits[3] = CPL::get<int> ("jcmax_olap");
        olap_limits[4] = CPL::get<int> ("kcmin_olap");
        olap_limits[5] = CPL::get<int> ("kcmax_olap");

    }

    void SocketLAMMPS::setupFixMDtoCFD(LAMMPS_NS::LAMMPS *lammps)
    {

        // Get the positions of CPL cell grid nodes
        int ncx = CPL::get<int> ("ncx"); 
        int ncy = CPL::get<int> ("ncy"); 
        int ncz = CPL::get<int> ("ncz"); 
        int xyShape[2] = {ncx + 1, ncy + 1};
        int zShape[1] = {ncz + 1};
        CPL::ndArray<double> cfd_xg (CPL::get<double*> ("xg"), 2, xyShape);
        CPL::ndArray<double> cfd_yg (CPL::get<double*> ("yg"), 2, xyShape);
        CPL::ndArray<double> cfd_zg (CPL::get<double*> ("zg"), 1, zShape);

        // Assume halos were not part of CFD grid passed in during cfd_init, 
        // minimum overlap should be bottom left of the CFD domain
        // -1 for F to C indexing!
        double xmin_cfd = cfd_xg(olap_limits[0] - 1, 1 - 1);
        double ymin_cfd = cfd_yg(1 - 1, olap_limits[2] - 1);
        double zmin_cfd = cfd_zg(olap_limits[4] - 1);

        // Cell sizes
        double dx = CPL::get<double> ("dx");
        double dy = CPL::get<double> ("dy");
        double dz = CPL::get<double> ("dz");
        // + dxyz to get top vertex of top cell
        // -1 for F to C indexing!
        double xmax_cfd = cfd_xg(olap_limits[1] - 1, 1 - 1) + dx;
        //double ymax_cfd = cfd_yg(1 - 1, olap_limits[3] - 1) + dy;
        double zmax_cfd = cfd_zg(olap_limits[5] - 1) + dz;

        // Averaging region in y to be cell below bottom cfd y cell
        double ymin_cfd_halo = ymin_cfd - dy;
        double ymax_cfd_halo = ymin_cfd_halo + dy;

        // Setup CFD BC averaging region (halo only in y direction)
        double bot_l[3] = {xmin_cfd, ymin_cfd_halo, zmin_cfd};
        double top_r[3] = {xmax_cfd, ymax_cfd_halo, zmax_cfd};
        SocketLAMMPS::map_cfd_to_lammps(bot_l, lammps);
        SocketLAMMPS::map_cfd_to_lammps(top_r, lammps);

        // LAMMPS commands to set up the fix

        // CFD BC region 
        std::string cmd = "region cfdbcregion block ";
        cmd += std::to_string (bot_l[0]) + " ";
        cmd += std::to_string (top_r[0]) + " ";
        cmd += std::to_string (bot_l[1]) + " ";
        cmd += std::to_string (top_r[1]) + " ";
        cmd += std::to_string (bot_l[2]) + " ";
        cmd += std::to_string (top_r[2]) + " ";
        cmd += "units box";
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

    void SocketLAMMPS::setupFixCFDtoMD(LAMMPS_NS::LAMMPS *lammps)
    {

        // Get the positions of CPL cell grid nodes
        int ncx = CPL::get<int> ("ncx"); 
        int ncy = CPL::get<int> ("ncy"); 
        int ncz = CPL::get<int> ("ncz"); 
        int xyShape[2] = {ncx + 1, ncy + 1};
        int zShape[1] = {ncz + 1};
        CPL::ndArray<double> cfd_xg (CPL::get<double*> ("xg"), 2, xyShape);
        CPL::ndArray<double> cfd_yg (CPL::get<double*> ("yg"), 2, xyShape);
        CPL::ndArray<double> cfd_zg (CPL::get<double*> ("zg"), 1, zShape);
        
        // Get the CFD -> MD constrained/forced cell extents
        // -1 for Fortran to C indexing
        int icmin_cnst = CPL::get<int> ("icmin_cnst") - 1; 
        int jcmin_cnst = CPL::get<int> ("jcmin_cnst") - 1; 
        int kcmin_cnst = CPL::get<int> ("kcmin_cnst") - 1; 
        int icmax_cnst = CPL::get<int> ("icmax_cnst") - 1; 
        int jcmax_cnst = CPL::get<int> ("jcmax_cnst") - 1; 
        int kcmax_cnst = CPL::get<int> ("kcmax_cnst") - 1; 

        // Get the spatial edges of the box defined by icmin_cnst -> kcmax_cnst
        double xlo = cfd_xg(icmin_cnst, 0); 
        double xhi = cfd_xg(icmax_cnst + 1, 0); 
        double ylo = cfd_yg(0, jcmin_cnst); 
        double yhi = cfd_yg(0, jcmax_cnst + 1); 
        double zlo = cfd_zg(kcmin_cnst); 
        double zhi = cfd_zg(kcmax_cnst + 1); 

        // Map CFD xlo, xhi, ylo, yhi, zlo, zhi to MD coords
        double bot_l[3] = {xlo, ylo, zlo};
        double top_r[3] = {xhi, yhi, zhi};
        SocketLAMMPS::map_cfd_to_lammps(bot_l, lammps);
        SocketLAMMPS::map_cfd_to_lammps(top_r, lammps);
        xlo = bot_l[0]; ylo = bot_l[1]; zlo = bot_l[2];
        xhi = top_r[0]; yhi = top_r[1]; zhi = top_r[2];
        
        // Tell LAMMPS to keep track of atoms in constrained region
        std::string cmd = "region cplforceregion block ";
        cmd += std::to_string(xlo) + " " + std::to_string(xhi) + " ";
        cmd += std::to_string(ylo) + " " + std::to_string(yhi) + " ";
        cmd += std::to_string(zlo) + " " + std::to_string(zhi) + " ";
        cmd += "units box";
        std::cout << cmd << std::endl;
        lammps->input->one (cmd.c_str());
        // TODO (d.trevelyan@ic.ac.uk) every 1? 
        cmd = "group cplforcegroup dynamic all region cplforceregion every 1";
        lammps->input->one (cmd.c_str());
        

        // Create a FixCPLForce instance
        cmd = "fix cplforcefix all cpl/force region cplforceregion";
        lammps->input->one (cmd.c_str());

    }

    // TODO develop a custom fix so that lammps doesn't need to do 
    // a global reduce (d.trevelyan@ic.ac.uk) ?
    void SocketLAMMPS::pack (const LAMMPS_NS::LAMMPS *lammps)
    {
        // Get LAMMPS computed velocity field
        int gridShape[4] = {4, ncxl, 1, nczl};
        int zeroShape[4] = {4, 0, 0, 0};
        recvVelocity.resize (4, zeroShape);
        sendVelocity.resize (4, gridShape);

        auto fixindex = lammps->modify->find_fix("cfdbcfix");
        auto fix = lammps->modify->fix[fixindex]; 

        for (int i = 0; i != sendVelocity.shape(1); ++i)
        {
            for (int k = 0; k != sendVelocity.shape(3); ++k)
            {
                int icoord = myCoords[0]*ncxl + i;
                int kcoord = myCoords[2]*ncxl + k;
                int row = icoord*sendVelocity.shape(3) + kcoord;
                //double x = fix->compute_array(row, 0);
                //double z = fix->compute_array(row, 1);
                double vx = fix->compute_array(row, 3);  
                double vy = fix->compute_array(row, 4);  
                double vz = fix->compute_array(row, 5);  

                sendVelocity(0, i, 0, k) = vx; // vx;
                sendVelocity(1, i, 0, k) = vy; // vy;
                sendVelocity(2, i, 0, k) = vz; // vz;
                sendVelocity(3, i, 0, k) = 1.0;
            }
        }
        
    }

    void SocketLAMMPS::unpack (const LAMMPS_NS::LAMMPS *lammps)
    {
        auto fixindex = lammps->modify->find_fix ("cplforcefix");
        auto cplfix = dynamic_cast<FixCPLForce*>
        (
            lammps->modify->fix[fixindex] 
        );
        auto stress = std::make_shared<CPL::ndArray<double>> (recvStress);
        cplfix->updateStress (stress);

        // Update the processor portion
        std::vector<int> portion (6);
        CPL::proc_portion (myFCoords.data(), CPL::md_realm, olap_limits.data(), portion.data());

        // - 1 for Fortran to C indexing
        for (auto &c : portion)
        {
            c -= 1;
        } 

        cplfix->updateProcPortion (portion.data());
    };

    void SocketLAMMPS::send()
    {

        // Limits are the same as overlap limits, but only bottom y cell
        std::vector<int> limits(olap_limits);
        limits[3] = limits[2];

        // Send the data to CFD
        CPL::gather
        (
            sendVelocity.data(),
            sendVelocity.shapeData(),
            limits.data(),
            recvVelocity.data(),
            recvVelocity.shapeData()
        );

    };

    void SocketLAMMPS::receive()
    {

        int zeroShape[4] = {9, 0, 0, 0};
        int recvShape[4] = {9, ncxl, ncyl, nczl};

        sendStress.resize(4, zeroShape);
        recvStress.resize(4, recvShape);

        // Receive from CFD
        CPL::scatter
        (
            sendStress.data(),
            sendStress.shapeData(),
            olap_limits.data(),
            recvStress.data(),
            recvStress.shapeData()
        );

    };

    void SocketLAMMPS::map_cfd_to_lammps
    (
        double *coord,
        LAMMPS_NS::LAMMPS *lammps
    )
    {
        double *globaldomain = lammps->domain->prd;
        double *temp = CPL::map_cfd2md_global(coord);
        for (int c = 0; c < 3; ++c)
        {
            coord[c] = temp[c] + globaldomain[c]/2.0;
        }
    }

}
