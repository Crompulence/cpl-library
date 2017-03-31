#include "fix_cpl_force.h"
#include "atom.h"
#include "universe.h"
#include "error.h"
#include "group.h"
#include "domain.h"
#include "region.h"
#include<iostream>
#include<memory>
#include<fstream>
#include "CPL.h"

FixCPLForce::FixCPLForce 
(
    LAMMPS_NS::LAMMPS *lammps,
    int narg,
    char **arg
)
    : Fix (lammps, narg, arg)
{
}


int FixCPLForce::\
setmask() {
  int mask = 0;
  mask |= LAMMPS_NS::FixConst::POST_FORCE;
  return mask;
}


void FixCPLForce::\
post_force (int vflag) {

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double dx = CPL::get<double> ("dx");     
    double dy = CPL::get<double> ("dy");     
    double dz = CPL::get<double> ("dz");     
    double dA = dx*dz;
    double dV = dx*dy*dz;
    
    //TODO: Pass this info as parameters of the fix
    char* groupstr = "cplforcegroup";
    char* regionstr = "cplforceregion";

    int cplforcegroup = group->find (groupstr);
    int groupbit = group->bitmask[cplforcegroup];

    int rid = domain->find_region (regionstr);
    auto cplforceregion = domain->regions[rid];

    // Preliminary summation, only 1 value per cell so can slice
    // away cfdStress->shape(0)
    int sumsShape[3] = {
                         cfdStress->shape(1),
                         cfdStress->shape(2),
                         cfdStress->shape(3)
                       };
    CPL::ndArray<double> gSums(3, sumsShape); // Sum of Flekkøy g weights
    CPL::ndArray<double> nSums(3, sumsShape); // Sum of number of particles
    nSums = 0.0;
    gSums = 0.0;
    for (int i = 0; i < nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {
            // Find in which cell number (local to processor) is the particle
            // and sum all the Flekkoy weights for each cell.
            int glob_cell[3];
            CPL::map_coord2cell(x[i][0], x[i][1], x[i][2], glob_cell);

            int loc_cell[3];
            bool validCell = CPL::map_glob2loc_cell(procPortion.data(), glob_cell, loc_cell);
            /**
            if (validCell)
            std::cout << "local: " << icell << " " << jcell << " " << kcell << " global: " << glob_cell[0] << " " << glob_cell[1] << " " << glob_cell[2] \
            << "portion:" << procPortion[0] << " "<< procPortion[1] << " "<< procPortion[2] << " "<< procPortion[3] << " "<< procPortion[4] << " "<< procPortion[5] << " " << std::endl;
            **/
            if (! validCell) {
               // std::cout << "Warning: an atom in the constrained region is within an invalid cell. \n"
               //           << "This should never happen and it is likely a BUG. Report." << std::endl;
                continue;
            }

            int icell = loc_cell[0];
            int jcell = loc_cell[1];
            int kcell = loc_cell[2];
            nSums(icell, jcell, kcell) += 1.0; 

            double g = flekkoyGWeight (x[i][1], cplforceregion->extent_ylo, 
                                       cplforceregion->extent_yhi);
            gSums(icell, jcell, kcell) += g;
            //std::cout << "FLEKOY: " << gSums(icell, jcell, kcell) << " " << cplforceregion->extent_ylo\
                << " " << cplforceregion->extent_yhi << " " << x[i][1]<< std::endl;
        }
    }


    // Calculate force and apply
    for (int i = 0; i < nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {
            int glob_cell[3];
            CPL::map_coord2cell(x[i][0], x[i][1], x[i][2], glob_cell);

            int loc_cell[3];
            bool validCell = CPL::map_glob2loc_cell(procPortion.data(), glob_cell, loc_cell);


            if (! validCell) {
//                std::cout << "Warning: an atom in the constrained region is within an invalid cell. \n"
 //                         << "This should never happen and it is likely a BUG. Report." << std::endl;
                continue;
            }

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
                double g = flekkoyGWeight (x[i][1], cplforceregion->extent_ylo,
                                           cplforceregion->extent_yhi);


                // Since the Flekkoy weight is considered only for 0 < y < L/2, for cells 
                // that are completely in y < 0 gSums(i, j, k) will be always 0.0 so can 
                // produce a NAN in the g/gSums division below.
                if (gSums(icell, jcell, kcell) > 0.0) {
                    double gdA = (g/gSums(icell, jcell, kcell)) * dA;
                    // Normal to the X-Z plane is (0, 1, 0) so (tauxy, syy, tauzy)
                    // are the only components of the stress tensor that matter.
                    double fx = gdA * cfdStress->operator()(1, icell, jcell, kcell);
                    double fy = gdA * cfdStress->operator()(4, icell, jcell, kcell);
                    double fz = gdA * cfdStress->operator()(7, icell, jcell, kcell);
                    f[i][0] += fx;
                    // Those should be zero for parallel flow! (CHECK WITH OPENFOAM). Unless one considers pressure from
                    // the continuum domain acting on the interface.
                    f[i][1] += fy;
                    f[i][2] += fz;
                    //f[i][1] -= (g/gsum) * dA * pressure
                }
            }
        }
    }
}

// See Flekkøy, Wagner & Feder, 2000 Europhys. Lett. 52 271, footnote p274
double FixCPLForce::\
flekkoyGWeight (double y, double ymin, double ymax) {
    
    // K factor regulates how to distribute the total force across the volumen.
    // 1/K represent the fraction of the constrain region volumen used.
    // Flekøy uses K = 2.
    double K = 1;
    // Define re-scaled coordinate 
    double L = ymax - ymin;
    //double yhat = y - ymin - 0.5*L; 
    double yhat = y - ymin - (1 - 1/K)*L;
    double g = 0.0;


//   if (yhat > 0.5*L) {
    if (yhat > (1/K * L)) {
        error->all(FLERR, "Position argument y to flekkoyGWeight "
                          "(y, ymin, ymax) is greater than ymax. ");
    }
    else if (yhat > 0.0) {
        g = 2.0*(1.0/(L - K*yhat) - 1.0/L - K*yhat/(L*L));
    }
    return 1; 
   // return g;
    
}

void FixCPLForce::\
updateStress (std::shared_ptr <CPL::ndArray <double>> stress) {

    cfdStress = std::move (stress);
}

void FixCPLForce::\
updateProcPortion (int inputPortion[]) {

    procPortion.resize(6);
    for (int i = 0; i < 6; ++i) {
        procPortion[i] = inputPortion[i];
    }
}
