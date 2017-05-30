#ifdef FIX_CLASS

FixStyle(cpl/force, FixCPLForce)

#else

#ifndef LMP_FIX_CPL_FORCE_H
#define LMP_FIX_CPL_FORCE_H

#include "fix.h"
#include "cpl/cpl.h"
#include "cpl/CPL_ndArray.h"
#include "cpl/CPL_force.h"
#include <memory>

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
    void end_of_step (); //todo add override <=== es205 17/01/17 WTF does this mean?
    void updateStress (CPL::ndArray<double>& stress);
    void updateProcPortion (int inputPortion[]);

private:

	CPL::ndArray<double>* cfdStress;
    std::vector<int> procPortion;
    double flekkoyGWeight (double y, double ymin, double ymax);

};

#endif
#endif
