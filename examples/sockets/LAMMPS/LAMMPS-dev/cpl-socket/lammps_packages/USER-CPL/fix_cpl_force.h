#ifdef FIX_CLASS

FixStyle(cpl/force, FixCPLForce)

#else

#ifndef LMP_FIX_CPL_FORCE_H
#define LMP_FIX_CPL_FORCE_H

#include<memory>

#include "fix.h"
#include "cpl/cpl.h"
#include "cpl/CPL_ndArray.h"
#include "cpl/CPL_force.h"

class FixCPLForce : public LAMMPS_NS::Fix {

public:

    FixCPLForce(class LAMMPS_NS::LAMMPS *lammps, int narg, char **arg);
    int setmask();
    void setup (int vflag); 
    void post_force (int vflag);
    void updateBuf (CPL::ndArray<double>& stress);
    void updateProcPortion (int inputPortion[]);
    std::shared_ptr<std::string> forcetype;
    std::unique_ptr<CPLForce> fxyz;

private:

	CPL::ndArray<double>* cfdBuf;
    std::vector<int> procPortion;

};

#endif
#endif
