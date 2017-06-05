#ifdef FIX_CLASS

FixStyle(cpl/force, FixCPLForce)

#else

#ifndef LMP_FIX_CPL_FORCE_H
#define LMP_FIX_CPL_FORCE_H

#include<memory>

#include "fix.h"
#include "cpl/cpl.h"
#include "cpl/CPL_ndArray.h"

class FixCPLForce : public LAMMPS_NS::Fix {

public:

    FixCPLForce 
    (
        class LAMMPS_NS::LAMMPS *lammps,
        int narg,
        char **arg
    );
    int setmask();
    void setup (int vflag); 
    void post_force (int vflag);
    void updateStress (CPL::ndArray<double>& stress);
    void updateProcPortion (int inputPortion[]);

private:

	CPL::ndArray<double>* cfdStress;
    std::vector<int> procPortion;
    double flekkoyGWeight (double y, double ymin, double ymax);

    CPL::ndArray<double> gSums; // Sum of Flekkøy g weights
    CPL::ndArray<double> nSums; // Sum of number of particles
};

#endif
#endif
