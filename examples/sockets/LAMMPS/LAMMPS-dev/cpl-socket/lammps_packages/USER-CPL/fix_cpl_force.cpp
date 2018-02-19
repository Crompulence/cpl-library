#include<iostream>
#include<memory>
#include<fstream>
#include <cmath>

#include "atom.h"
#include "universe.h"
#include "error.h"
#include "group.h"
#include "domain.h"
#include "region.h"

#include "fix_cpl_force.h"

#include "cpl/CPL_ndArray.h"

FixCPLForce::FixCPLForce(LAMMPS_NS::LAMMPS *lammps, int narg, char **arg)
    : Fix (lammps, narg, arg)
{
   //nevery = 50;//cplsocket.timestep_ratio;
    for (int iarg=0; iarg<narg; iarg+=1){
        std::string arguments(arg[iarg]);
        if (arguments == "forcetype")
            forcetype = std::make_shared<std::string>(arg[iarg+1]);
    }
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

//NOTE -- Not actually called post force, for some reason
// this no longer works reliabley in LAMMPS, instead call
// explicitly in CPLInit!
void FixCPLForce::post_force(int vflag) {

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *radius = atom->radius;
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
    // away cfdBuf->shape(0)
    std::string fxyzType(*forcetype);
    // std::unique_ptr<CPLForce>  fxyz(nullptr); //Moved to header
    if (fxyzType.compare("Flekkoy") == 0) {
        fxyz.reset(new CPLForceFlekkoy(9, cfdBuf->shape(1), 
                                          cfdBuf->shape(2), 
                                          cfdBuf->shape(3)));
    } else if (fxyzType.compare("test") == 0) {
        fxyz.reset(new CPLForceTest(3, cfdBuf->shape(1), 
                                       cfdBuf->shape(2), 
                                       cfdBuf->shape(3)));
    } else if (fxyzType.compare("Velocity") == 0) {
        fxyz.reset(new CPLForceVelocity(3, cfdBuf->shape(1), 
                                           cfdBuf->shape(2), 
                                           cfdBuf->shape(3)));
    } else if (fxyzType.compare("Drag") == 0) {
        fxyz.reset(new CPLForceDrag(3, cfdBuf->shape(1), 
                                       cfdBuf->shape(2), 
                                       cfdBuf->shape(3))); 
        fxyz->calc_preforce = true;

//    } else if (fxyzType.compare("Ergun") == 0) {
//        fxyz.reset(new CPLForceGranular(3, cfdBuf->shape(1), 
//                                           cfdBuf->shape(2), 
//                                           cfdBuf->shape(3))); 
    } else if (fxyzType.compare("Di_Felice") == 0) {
        fxyz.reset(new CPLForceGranular(3, cfdBuf->shape(1), 
                                           cfdBuf->shape(2), 
                                           cfdBuf->shape(3))); 
//    } else if (fxyzType.compare("Tang_et_al") == 0) {
//        fxyz.reset(new CPLForceGranular(3, cfdBuf->shape(1), 
//                                           cfdBuf->shape(2), 
//                                           cfdBuf->shape(3))); 

    } else {
        std::string cmd("CPLForce type ");
        cmd += fxyzType + " not defined";
        throw std::runtime_error(cmd);
    }

    //std::cout << "CPL Force type  " << fxyzType << fxyz->calc_preforce << std::endl;


    //Set CPLForce min/max to local processor limits using values from CPL library 
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
    fxyz->set_minmax(min, max);

    double fx=0., fy=0., fz=0.;
    double mi, radi, pot, xi[3], vi[3], ai[3];
    pot = 1.0; //Interaction Potential should be set here
    std::vector<int> cell;
    std::vector<double> fi(3);  
    //Should we reset sums here?
    fxyz->resetsums();
    if (fxyz->calc_preforce) {
        for (int i = 0; i < nlocal; ++i)
        {
            if (mask[i] & groupbit)
            {

                //Get local molecule data
                mi = rmass[i];
                radi = radius[i];
                for (int n=0; n<3; n++){
                    xi[n]=x[i][n]; 
                    vi[n]=v[i][n]; 
                    ai[n]=f[i][n];
                }

                // Sum all the weights for each cell.
                fxyz->pre_force(xi, vi, ai, mi, radi, pot);

            }
        }

    }


    // Calculate force and apply
    for (int i = 0; i < nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {

            //Get local molecule data
            mi = rmass[i];
            radi = radius[i];
            for (int n=0; n<3; n++){
                xi[n]=x[i][n]; 
                vi[n]=v[i][n]; 
                ai[n]=f[i][n];
            }

            //Get force from object
            fxyz->set_field(*cfdBuf);
            fi = fxyz->get_force(xi, vi, ai, mi, radi, pot);

            //Apply force
            f[i][0] += fi[0];
            f[i][1] += fi[1];
            f[i][2] += fi[2];

            std::cout.precision(17);
            std::cout << "Force " << i << " " << xi[2] <<  " " << vi[2] << " " << ai[2] << " " <<
                          mi << " " << fi[0] << " " << fi[1] << " " << fi[2] << " "  
                          << f[i][0] << " " << f[i][1] << " " << f[i][2] << std::endl;

        }
    }
}

void FixCPLForce::updateBuf (CPL::ndArray<double>& Buf) {
    cfdBuf = &Buf;
}

void FixCPLForce::updateProcPortion (int inputPortion[]) {

    procPortion.resize(6);
    for (int i = 0; i < 6; ++i) {
        procPortion[i] = inputPortion[i];
    }
}

