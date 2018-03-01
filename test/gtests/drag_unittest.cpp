#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <vector>

#include "cpl.h"
#include "CPL_field.h"
#include "CPL_force.h"
#include "CPL_ndArray.h"

// The fixture for testing class Foo.
class CPL_drag_Test : public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  CPL_drag_Test() {
    // You can do set-up work for each test here.
  }

  virtual ~CPL_drag_Test() {
    // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  virtual void SetUp() {
    // Code here will be called immediately after the constructor (right
    // before each test).
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).
  }

  // Objects declared here can be used by all tests in the test case for Foo.
};

#define trplefor(Ni,Nj,Nz) for (int i = 0; i<Ni; ++i){ \
                           for (int j = 0; j<Nj; ++j){ \
                           for (int k = 0; k<Nz; ++k)


//Test for CPL::ndArray - setup and array size 
TEST_F(CPL_drag_Test, All_Drags) {

    std::shared_ptr<CPLForceDrag> Drag, Stokes, Di_Felice, BVK, Ergun;

    int nd = 9; 
    int N = 20;
    int icell = N;
    int jcell = N;
    int kcell = N;

    //Setup one particle per cell
    double r[3] = {0.0, 0.0, 0.0};
    double v[3] = {0.001, 0.0, 0.0};
    double a[3] = {0.0, 0.0, 0.0};
    std::vector<double> F, Ur;
    //double volume = (4./3.)*M_PI*pow(radius,3);
    double m=1.; double e=1.;

//    CPLForceStokes Stokes(nd, icell, jcell, kcell);
//    CPLForceDi_Felice Di_Felice(nd, icell, jcell, kcell);
//    CPLForceBVK BVK(nd, icell, jcell, kcell);
//    CPLForceErgun Ergun(nd, icell, jcell, kcell);

    Drag = std::make_shared<CPLForceDrag>(nd, icell, jcell, kcell);
    //Granular = std::make_shared<CPLForceGranular>(nd, icell, jcell, kcell);
    Stokes = std::make_shared<CPLForceStokes>(nd, icell, jcell, kcell);
    Di_Felice = std::make_shared<CPLForceDi_Felice>(nd, icell, jcell, kcell);
    BVK = std::make_shared<CPLForceBVK>(nd, icell, jcell, kcell);
    Ergun = std::make_shared<CPLForceErgun>(nd, icell, jcell, kcell);

    std::vector<std::shared_ptr<CPLForceDrag>> forces{Drag, Stokes, Di_Felice, BVK, Ergun};


    double Volfrac = 0.5/float(N);
    double s;
    int maxn = 200;
    int i = 5; int j = 5; int k = 5;
    r[0] = i/double(icell);
    r[1] = j/double(jcell);
    r[2] = k/double(kcell);

    for ( auto &f : forces ) { 
        std::cout << "==========================================" << std::endl;
        for (int n=1; n < maxn; n++) {
            f->resetsums();
            s = Volfrac*n/maxn;

            //Setup esum
            f->pre_force(r, v, a, m, s, e);

            //Get Force
            F = f->get_force(r, v, a, m, s, e);
            std::cout << "Drag Unittest: " << typeid(*f).name() << " " << s << " "
                      << F[0] << " " << F[1] << " " << F[2] << std::endl;
        }
        std::cout << "==========================================" << std::endl;

    }
    
}


