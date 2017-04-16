
#include "cpl_force.h"
#include <vector>
#include "CPL.h"

//Constructors
CPLForce::CPLForce(int nd, int icells, int jcells, int kcells){
    // Fields
    int fieldShape[4] = {nd, icells, jcells, kcells};
    field.resize(4, fieldShape);
};

CPLForce::CPLForce(CPL::ndArray<double> fieldin){
    field = fieldin;
};

//CPLForce::CPLForce(std::shared_ptr <CPL::ndArray <double>> fieldin){
//    field = fieldin;
//};



void CPLForce::set_field(CPL::ndArray<double> fieldin){

    field = fieldin;
    
};

CPL::ndArray<double> CPLForce::get_field(){

    return field;

};


//==========================================================
//The Flekkoy example of this would then be CPLForceFlekkoy
// as Flekkoy's constraint IS A type of CPL force
//==========================================================

//Constructor of datatype
CPLForceFlekkoy::CPLForceFlekkoy(CPL::ndArray<double> field) : CPLForce(field){

    initialisesums(field);
}

CPLForceFlekkoy::CPLForceFlekkoy(int nd, int icells, int jcells, int kcells) : CPLForce(nd, icells, jcells, kcells){

    initialisesums(field);

}

void CPLForceFlekkoy::initialisesums(CPL::ndArray<double> f){

    int sumsShape[3] = {
                            f.shape(1),
                            f.shape(2),
                            f.shape(3)
                        };

    gSums.resize(3, sumsShape); // Sum of Flekkøy g weights
    nSums.resize(3, sumsShape); // Sum of number of particles  
    nSums = 0.0;  gSums = 0.0;
}


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
void pre_force(int icell, int jcell, int kcell) {
    double g = flekkoyGWeight(x[i][1], cplforceregion->extent_ylo, 
                              cplforceregion->extent_yhi);
    nSums(icell, jcell, kcell) += 1.0; 
    gSums(icell, jcell, kcell) += g;
    std::cout << "FLEKKOY: " << gSums(icell, jcell, kcell) << " " << cplforceregion->extent_ylo\
        << " " << cplforceregion->extent_yhi << " " << x[i][1]<< std::endl;

}
/*
//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
void apply_force(int i, int icell, int jcell, int kcell, 
                 std::shared_ptr <CPL::ndArray <double>> field){
    double n = nSums(i, icell, jcell, kcell);
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
            // Normal to the X-Z plane is (0, 1, 0) so (tauxy, syy, tauxy)
            // are the only components of the stress tensor that matter.
            double fx = gdA * field->operator()(1, icell, jcell, kcell);
            double fy = gdA * field->operator()(4, icell, jcell, kcell);
            double fz = gdA * field->operator()(7, icell, jcell, kcell);

        }
    }
}

// See Flekkøy, Wagner & Feder, 2000 Europhys. Lett. 52 271, footnote p274
double FixCPLForce::flekkoyGWeight(double y, double ymin, double ymax) {
    
    // K factor regulates how to distribute the total force across the volume.
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
    
    return g;
    
}


//=====================================================
//The post force fix in LAMMPS would then become
//=====================================================

void FixCPLForce::post_force (int vflag) {

    double **x = atom->x;
    double **f = atom->f;
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

    //Constructor here for forces based on input
    if ("FORCE_TYPE" == "Flekkoy"){
        fxyz = CPLForceFlekkoy->initialise(field)
    } else ("FORCE_TYPE" == "Nieetal"){
        fxyz = CPLNieetal->initialise(field)
    }

    // First loop to sum up all quantities and prepare to apply force
    for (int i = 0; i < nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {
            // Find in which cell number (local to processor) is the particle
            // and sum all the Flekkøy weights for each cell.
            int glob_cell[3];
            CPL::map_coord2cell(x[i][0], x[i][1], x[i][2], glob_cell);

            int loc_cell[3];
            bool validCell = CPL::map_glob2loc_cell(procPortion.data(), glob_cell, loc_cell);

            //if (validCell)
            //std::cout << "local: " << icell << " " << jcell << " " << kcell << " global: " << glob_cell[0] << " " << glob_cell[1] << " " << glob_cell[2] \
            //<< "portion:" << procPortion[0] << " "<< procPortion[1] << " "<< procPortion[2] << " "<< procPortion[3] << " "<< procPortion[4] << " "<< procPortion[5] << " " << std::endl;

            if (! validCell) {
               // std::cout << "Warning: an atom in the constrained region is within an invalid cell. \n"
               //           << "This should never happen and it is likely a BUG. Report." << std::endl;
                continue;
            }

            int icell = loc_cell[0];
            int jcell = loc_cell[1];
            int kcell = loc_cell[2];

            //This is call to update pre-force 
            fxyz->pre_force(icell, jcell, kcell)

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

            // Get forces
            [fx, fy, fz] = fxyz->get_force(i, icell, jcell, kcell, field)

            f[i][0] += fx;
            f[i][1] += fy;
            f[i][2] += fz;

        }
    }
}
*/
