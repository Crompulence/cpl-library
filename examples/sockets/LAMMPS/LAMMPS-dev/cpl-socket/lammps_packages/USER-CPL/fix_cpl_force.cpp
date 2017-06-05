#include<iostream>
#include<memory>
#include<fstream>
#include <cmath>

#include "fix_cpl_force.h"
#include "atom.h"
#include "universe.h"
#include "error.h"
#include "group.h"
#include "domain.h"
#include "region.h"

#include "cpl/CPL_ndArray.h"
#include "cpl/CPL_force.h"

FixCPLForce::FixCPLForce
(
    LAMMPS_NS::LAMMPS *lammps,
    int narg,
    char **arg
)
    : Fix (lammps, narg, arg)
{
   //nevery = 1;//cplsocket.timestep_ratio;
}

//NOTE: It is called from fixCPLInit initial_integrate() now.
int FixCPLForce::setmask() {
  int mask = 0;
  //mask |= LAMMPS_NS::FixConst::POST_FORCE;
  return mask;
}


/* ---------------------------------------------------------------------- */

void FixCPLForce::setup(int vflag)
{
    post_force(vflag);
}


void FixCPLForce::post_force(int vflag) {

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double dx = CPL::get<double> ("dx");     
    double dy = CPL::get<double> ("dy");     
    double dz = CPL::get<double> ("dz");     
    double dA = dx*dz;
    double dV = dx*dy*dz;
    
    char* groupstr = "cplforcegroup";
    char* regionstr = "cplforceregion";

    int cplforcegroup = group->find (groupstr);
    int groupbit = group->bitmask[cplforcegroup];

    int rid = domain->find_region (regionstr);
    auto cplforceregion = domain->regions[rid];

    // Preliminary summation, only 1 value per cell so can slice
    // away cfdStress->shape(0)
    int sumsShape[3] = {cfdStress->shape(1), cfdStress->shape(2), cfdStress->shape(3)};
    gSums.resize(3, sumsShape); // Sum of Flekkøy g weights
    nSums.resize(3, sumsShape); // Sum of number of particles
    nSums = 0.0; gSums = 0.0;

    std::string CPLForce = "test";
//    if (CPLForce == "Flekkoy"){
        auto fxyz = CPLForceFlekkoy(9, cfdStress->shape(1), 
                                       cfdStress->shape(2), 
                                       cfdStress->shape(3));
//    } else if (CPLForce == "Velocity") {
//        auto fxyz = CPLForceVelocity(3, cfdStress->shape(1), 
//                                        cfdStress->shape(2), 
//                                        cfdStress->shape(3));
//    } else if (CPLForce == "test") {
//        auto fxyz = CPLForceTest(3, cfdStress->shape(1), 
//                                    cfdStress->shape(2), 
//                                    cfdStress->shape(3));
//    } else if (CPLForce == "Granular") {
//        auto fxyz = CPLForceGranular(3, cfdStress->shape(1), 
//                                        cfdStress->shape(2), 
//                                        cfdStress->shape(3));
//    }

	double min[3]; double max[3];
    std::vector<int> cnstFPortion(6);
    std::vector<int> cnstFRegion(6);
    CPL::get_cnst_limits(cnstFRegion.data());
    CPL::my_proc_portion (cnstFRegion.data(), cnstFPortion.data());
	//MIN
	CPL::map_cell2coord(cnstFPortion[0], cnstFPortion[2], cnstFPortion[4], min);
	//MAX
	CPL::map_cell2coord(cnstFPortion[1], cnstFPortion[3], cnstFPortion[5], max);
	max[0] += dx;
	max[1] += dy;
	max[2] += dz;
    //These are global coordinates
//    min[0] = cplforceregion->extent_xlo;
//    min[1] = cplforceregion->extent_ylo;
//    min[2] = cplforceregion->extent_zlo;
//    max[0] = cplforceregion->extent_xhi;
//    max[1] = cplforceregion->extent_yhi;
//    max[2] = cplforceregion->extent_zhi;
	fxyz.set_minmax(min, max);

//    std::cout << "cells " << min[0] << " " << min[1] << " " << min[2] << " " <<
//                             max[0] << " " << max[1] << " " << max[2] << " " << 
//                             cplforceregion->extent_xlo << " " << 
//                             cplforceregion->extent_ylo << " " <<  
//                             cplforceregion->extent_zlo << " " << 
//                             cplforceregion->extent_xhi << " " <<  
//                             cplforceregion->extent_yhi << " " <<  
//                             cplforceregion->extent_zhi << " " << std::endl;

//    std::cout << "cells " << min[0] << " " << min[1] << " " << min[2] << " " <<
//                             max[0] << " " << max[1] << " " << max[2] << " " << std::endl;
//    std::cout << sumsShape[0] << " " << sumsShape[1] << " " << sumsShape[2]  << " " <<
//                cnstFPortion[0] << " " << cnstFPortion[1] << " " << cnstFPortion[2] << " " <<
//                cnstFPortion[3] << " " << cnstFPortion[4] << " " << cnstFPortion[5] << std::endl;
//    for (int i = cnstFPortion.data()[0]; i<cnstFPortion.data()[1]; ++i){ 
//        for (int j = cnstFPortion.data()[2]; j<cnstFPortion.data()[3]; ++j){ 
//            for (int k = cnstFPortion.data()[4]; k<cnstFPortion.data()[5]; ++k){
//                    std::cout << "Pre-loop nsum " << i << " " << j << " " << k << " " 
//                              << nSums(i,j,k) << std::endl;
//    }}}


    double fx=0., fy=0., fz=0.;
    double mi, xi[3], vi[3], ai[3];
    std::vector<int> cell;
    std::vector<double> fi(3);
    for (int i = 0; i < nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {

            //Get local molecule data
            mi = rmass[i];
            for (int n=0; n<3; n++){
                xi[n]=x[i][n]; 
                vi[n]=v[i][n]; 
                ai[n]=f[i][n];
            }

            // Find in which cell number (local to processor) is the particle
            // and sum all the Flekkøy weights for each cell.
            int glob_cell[3];
            CPL::map_coord2cell(xi[0], xi[1], xi[2], glob_cell);

            int loc_cell[3];
            bool validCell = CPL::map_glob2loc_cell(procPortion.data(), glob_cell, loc_cell);
            /**
            if (validCell)
            std::cout << "local: " << icell << " " << jcell << " " << kcell << 
                        " global: " << glob_cell[0] << " " << glob_cell[1] << " " << glob_cell[2] \
                     << "portion:" << procPortion[0] << " "<< procPortion[1] << " "<< procPortion[2] 
                     << " "<< procPortion[3] << " "<< procPortion[4] << " "<< procPortion[5] << " " << std::endl;
            **/
            if (! validCell) {
               std::cout << "Warning: an atom in the constrained region is within an invalid cell. \n"
                          << "This should never happen and it is likely a BUG. Report." << std::endl;
                continue;
            }

            int icell = loc_cell[0];
            int jcell = loc_cell[1];
            int kcell = loc_cell[2];

            cell = fxyz.get_cell(xi);

            if (cell[0] != icell)
                throw std::domain_error("x cell Error CPL map_coord2cell doesn't match CPL_force");
            if (cell[1] != jcell)
                throw std::domain_error("y cell Error CPL map_coord2cell doesn't match CPL_force");
            if (cell[2] != kcell)
                throw std::domain_error("z cell Error CPL map_coord2cell doesn't match CPL_force");

            double g = flekkoyGWeight (xi[1], 
                                       cplforceregion->extent_ylo, 
                                       cplforceregion->extent_yhi);
            nSums(icell, jcell, kcell) += 1.0; 
            gSums(icell, jcell, kcell) += g;

            fxyz.pre_force(xi, vi, ai);            
            //std::cout << "FLEKKOY: " << gSums(icell, jcell, kcell) << " " << cplforceregion->extent_ylo\
                << " " << cplforceregion->extent_yhi << " " << xi[1]<< std::endl;
        }
    }
    for (int i = cnstFPortion.data()[0]; i<cnstFPortion.data()[1]; ++i){ 
        for (int j = cnstFPortion.data()[2]; j<cnstFPortion.data()[3]; ++j){ 
            for (int k = cnstFPortion.data()[4]; k<cnstFPortion.data()[5]; ++k){
                if (std::abs(fxyz.nSums(i,j,k)-nSums(i,j,k)) > 1e-5){
                    std::cout << "Error in nsum " << i << " " << j << " " << k << " " 
                              << fxyz.nSums(i,j,k) << " " << nSums(i,j,k) << " end" << std::endl;
                    //throw std::domain_error("CPL fix nsum doesn't match CPL_force");
                }
                if (std::abs(fxyz.gSums(i,j,k)-gSums(i,j,k)) > 1e-5){
                    std::cout << "Error in gsum " << i << " " << j << " " << k << " " 
                              << fxyz.gSums(i,j,k) << " " << gSums(i,j,k) << std::endl;
                    //throw std::domain_error("CPL fix gsum doesn't match CPL_force");
                }
    }}}

    // Calculate force and apply
    for (int i = 0; i < nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {

            //Get local molecule data
            mi = rmass[i];
            for (int n=0; n<3; n++){
                xi[n]=x[i][n]; 
                vi[n]=v[i][n]; 
                ai[n]=f[i][n];
            }

            int glob_cell[3];
            CPL::map_coord2cell(xi[0], xi[1], xi[2], glob_cell);

            int loc_cell[3];
            bool validCell = CPL::map_glob2loc_cell(procPortion.data(), glob_cell, loc_cell);


            if (! validCell) {
//                std::cout << "Warning: an atom in the constrained region is within an invalid cell. \n"
 //                         << "This should never happen and it is likely a BUG. Report." << std::endl;
                continue;
            }

            //Get force from object
            fxyz.set_field(*cfdStress);
            fi = fxyz.get_force(xi, vi, ai);

            int icell = loc_cell[0];
            int jcell = loc_cell[1];
            int kcell = loc_cell[2];
            double n = nSums(icell, jcell, kcell);
            if (n < 1.0) {
                std::cout << "Warning: 0 particles in cell (" 
                          << icell << ", " << jcell << ", " << kcell << ")"
                          << std::endl;
            }
            else {
                double g = flekkoyGWeight (xi[1], 
                                           cplforceregion->extent_ylo,
                                           cplforceregion->extent_yhi);

                // Since the Flekkoy weight is considered only for 0 < y < L/2, for cells 
                // that are completely in y < 0 gSums(i, j, k) will be always 0.0 so can 
                // produce a NAN in the g/gSums division below.
                if (gSums(icell, jcell, kcell) > 0.0) {
                    double gdA = (g/gSums(icell, jcell, kcell)) * dA;

                    // Normal to the X-Z plane is (0, 1, 0) so (tauxy, syy, tauxy)
                    // are the only components of the stress tensor that matter.
                    fx = gdA * cfdStress->operator()(1, icell, jcell, kcell);
                    fy = gdA * cfdStress->operator()(4, icell, jcell, kcell);
                    fz = gdA * cfdStress->operator()(7, icell, jcell, kcell);

                }
            }

            if (fi[0] != fx){
                std::cout << "Error in Force fx " << i << " " << icell << " " << jcell << " " << kcell << " " 
                          << fi[0] << " " << fx << " end" << std::endl;
                throw std::domain_error("Force Error -- doesn't match CPL_force implementation");
            }
            if (fi[1] != fy){
                std::cout << "Error in fy " << i << " " << icell << " " << jcell << " " << kcell << " " 
                          << fi[1] << " " << fy << " end" << std::endl;
                throw std::domain_error("Force Error -- doesn't match CPL_force implementation");

            }
            if (fi[2] != fz){
                std::cout << "Error in fz " << i << " " << icell << " " << jcell << " " << kcell << " "
                          << fi[2] << " " << fz << " end" << std::endl;
                throw std::domain_error("Force Error -- doesn't match CPL_force implementation");

            }

            //Apply force
            f[i][0] += fi[0];
            f[i][1] += fi[1];
            f[i][2] += fi[2];

            std::cout << "Force " << i << " " << icell <<  " " << jcell << " " << kcell << " " <<
                          cfdStress->operator()(7, icell, jcell, kcell) << " " <<
                          mi << " " << fi[0] << " " << fi[1] << " " << fi[2] << " "
                          << fx << " " << fy << " " << fz << " " << vi[2] << " "  
                          << f[i][0] << " " << f[i][1] << " " << f[i][2] << std::endl;
        }
    }
}

// See Flekkøy, Wagner & Feder, 2000 Europhys. Lett. 52 271, footnote p274
double FixCPLForce::flekkoyGWeight(double y, double ymin, double ymax) {
    
    // K factor regulates how to distribute the total force across the volume.
    // 1/K represent the fraction of the constrain region volumen used.
    // Flekøy uses K = 2.
    double K = 1.0;
    // Define re-scaled coordinate 
    double L = ymax - ymin;
    double yhat = y - ymin - (1.0 - 1.0/K)*L;
    double g = 0.0;


//   if (yhat > 0.5*L) {
    if (yhat > (1/K * L)) {
        error->all(FLERR, "Position argument y to flekkoyGWeight "
                          "(y, ymin, ymax) is greater than ymax. ");
    }
    else if (yhat > 0.0) {
        g = 2.0*(1.0/(L - K*yhat) - 1.0/L - K*yhat/(L*L));
    }
    
    return g;
    
}


void FixCPLForce::updateStress (CPL::ndArray<double>& stress) {

    cfdStress = &stress;
}

void FixCPLForce::updateProcPortion (int inputPortion[]) {

    procPortion.resize(6);
    for (int i = 0; i < 6; ++i) {
        procPortion[i] = inputPortion[i];
    }
}
