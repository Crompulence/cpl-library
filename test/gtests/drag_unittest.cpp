#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <vector>

// basic file operations
#include <iostream>
#include <fstream>


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

#define trplefor(Ni,Nj,Nz) for (int ic = 0; ic<Ni; ++ic){ \
                           for (int jc = 0; jc<Nj; ++jc){ \
                           for (int kc = 0; kc<Nz; ++kc)


//Test for CPL::ndArray - setup and array size 
TEST_F(CPL_drag_Test, All_Drags) {

    std::shared_ptr<CPLForceDrag> Drag, Stokes, Di_Felice, BVK, Tang, Ergun;

    int nd = 9; 
    int N = 20;
    int icell = N;
    int jcell = N;
    int kcell = N;
    double Lx = 1.0;
    double Ly = 1.0;
    double Lz = 1.0;
    double xorigin = 0.0;
    double yorigin = 0.0;
    double zorigin = 0.0;
    double min[3] = {xorigin, yorigin, zorigin};
    double max[3] = {Lx, Ly, Lz};
    double dx = Lx/float(icell);
    double dy = Ly/float(jcell);
    double dz = Lz/float(kcell);

    //Setup one particle per cell
    double r[3] = {0.0, 0.0, 0.0};
    double v[3] = {0.001, 0.0, 0.0};
    double a[3] = {0.0, 0.0, 0.0};
    double m=1.; double e=1.;
    double radius;
    std::vector<double> F;

    Drag = std::make_shared<CPLForceDrag>(nd, icell, jcell, kcell);
    //Granular = std::make_shared<CPLForceGranular>(nd, icell, jcell, kcell);
    Stokes = std::make_shared<CPLForceStokes>(nd, icell, jcell, kcell);
    Di_Felice = std::make_shared<CPLForceDi_Felice>(nd, icell, jcell, kcell);
    BVK = std::make_shared<CPLForceBVK>(nd, icell, jcell, kcell);
    Tang = std::make_shared<CPLForceTang>(nd, icell, jcell, kcell);
    Ergun = std::make_shared<CPLForceErgun>(nd, icell, jcell, kcell);

    std::vector<std::shared_ptr<CPLForceDrag>> forces{Drag, Stokes, Di_Felice, BVK, Tang, Ergun};
    for ( auto &f : forces ) {
        f->set_minmax(min, max);
    }

    //Choose a cell
    int i = 5; int j = 5; int k = 5;

    //Particles per cell
    int nx = 4; int ny = 4; int nz = 4;
    int Npercell = nx*ny*nz;
    double Volfrac = 0.6*(std::min(std::min(dx,dy),dz))/float(std::min(std::min(nx,ny),nz));
    int maxn = 200; //Number of sizes to loop over

    std::string fdir("./drag_output/");
    std::ofstream myfile;
    for ( auto &f : forces ) {
        std::string force_type(typeid(*f).name());
        std::string filename(force_type.substr(10));
        myfile.open(fdir+filename, std::ios::out | std::ios::trunc);
        //Write file header
        myfile << "phi" << ", " << "D" << ", "  
               <<  "v[0]" << ", " << "v[1]" << ", " << "v[2]"  << ", " 
               <<  "F[0]" << ", " << "F[1]" << ", " << "F[2]" << std::endl;
        std::cout << "==========================================" << std::endl;
        for (int n=1; n < maxn; n++) {
            f->resetsums();
            radius = Volfrac*n/double(maxn);

            //Setup esum by sticking lots of particles per cell
            for (int ip=0; ip<nx; ip++) {
            for (int jp=0; jp<ny; jp++) {
            for (int kp=0; kp<nz; kp++) {
                r[0] = Lx*i/double(icell) + dx*(0.5/double(nx)+ip/double(nx));
                r[1] = Ly*j/double(jcell) + dy*(0.5/double(ny)+jp/double(ny));
                r[2] = Lz*k/double(kcell) + dz*(0.5/double(nz)+kp/double(nz));
                f->pre_force(r, v, a, m, radius, e);
                //std::cout << "particles " << r[0] << " " << r[1] << " " << r[2] << std::endl;

            }}}

            //Get volSums
            auto volSums = f->get_internal_fields("volSums");
            trplefor(icell,jcell,kcell) {
                double phi = volSums->get_array_value(0, ic, jc, kc)/volSums->get_dV();
                if (phi != 0)
                    std::cout << "phi: " << ic << " " << jc << " " << kc << " " << phi << std::endl;
            }}}
            double phi = volSums->get_array_value(0, i, j, k)/volSums->get_dV();
            double eps = 1.0 - phi;

            //Get Force
            F = f->get_force(r, v, a, m, radius, e);
            std::cout << "Drag Unittest: " << typeid(*f).name() << " " << 2.*radius << " " 
                      << phi << " "<< F[0] << " " << F[1] << " " << F[2] << std::endl;

            myfile << phi << ", " << 2.*radius << ", "  
                   <<  v[0] << ", " << v[1] << ", " << v[2]  << ", " 
                   <<  F[0] << ", " << F[1] << ", " << F[2] << std::endl;
        }
        std::cout << "==========================================" << std::endl;
        myfile.close();
    }
    
}


