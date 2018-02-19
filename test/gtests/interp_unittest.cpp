#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <vector>

#include "cpl.h"
#include "CPL_field.h"
#include "CPL_force.h"
#include "CPL_ndArray.h"

// The fixture for testing class Foo.
class CPL_interp_Test : public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  CPL_interp_Test() {
    // You can do set-up work for each test here.
  }

  virtual ~CPL_interp_Test() {
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
TEST_F(CPL_interp_Test, test_celltonode) {
    int nd = 3;
    int icell = 3;
    int jcell = 3;
    int kcell = 3;

    CPL::ndArray<double> buf;
    int shape[4] = {nd, icell, jcell, kcell};
    buf.resize (4, shape);

    //Test define
    int n = 0;
    CPL::CPLField field(buf);
    //Slabs in x
    for (int j = 0; j < jcell; j++ ){
    for (int k = 0; k < kcell; k++ ){
        buf(n,0,j,k) = -0.5;
        buf(n,1,j,k) = 0.5;
        buf(n,2,j,k) = 1.5;
    }}
    CPL::ndArray<double> nodesx = field.celltonode(buf, n, 0, 0, 0);
    trplefor(icell-1,jcell-1,kcell-1){
        //std::cout << "nodex = " << i << " " << j << " " << k << " " 
        //                        << nodesx(n,i,j,k) << "\n"; 
        ASSERT_DOUBLE_EQ(nodesx(n,i,j,k), float(i));
    }}}

    //Slabs in y
    for (int i = 0; i < 3; i++ ){
    for (int k = 0; k < 3; k++ ){
        buf(n,i,0,k) = -0.5;
        buf(n,i,1,k) = 0.5;
        buf(n,i,2,k) = 1.5;
    }}
    CPL::ndArray<double> nodesy = field.celltonode(buf, n, 0, 0, 0);
    trplefor(icell-1,jcell-1,kcell-1){
        //std::cout << "nodey = " << i << " " << j << " " << k << " " 
        //                        << nodesy(n,i,j,k) << "\n"; 
        ASSERT_DOUBLE_EQ(nodesy(n,i,j,k), float(j));

    }}}

    //Slabs in z
    for (int i = 0; i < 3; i++ ){
    for (int j = 0; j < 3; j++ ){
        buf(n,i,j,0) = -0.5;
        buf(n,i,j,1) = 0.5;
        buf(n,i,j,2) = 1.5;
    }}
    CPL::ndArray<double> nodesz = field.celltonode(buf, n, 0, 0, 0);
    trplefor(icell-1,jcell-1,kcell-1){
        //std::cout << "nodez = " << i << " " << j << " " << k << " " 
        //                        << nodesz(n,i,j,k) << "\n"; 
        ASSERT_DOUBLE_EQ(nodesz(n,i,j,k), float(k));

    }}}

};

//Test for CPL::interpolate 
TEST_F(CPL_interp_Test, test_interpolate) {

    int nd = 3;
    int icell = 3;
    int jcell = 3;
    int kcell = 3;
    double *xi;

    CPL::ndArray<double> buf;
    int shape[4] = {nd, icell, jcell, kcell};
    buf.resize (4, shape);

    //Test define
    int n = 0;
    //Slabs in x
    for (int j = 0; j < jcell; j++ ){
    for (int k = 0; k < kcell; k++ ){
        buf(n,0,j,k) = -0.5;
        buf(n,1,j,k) = 0.5;
        buf(n,2,j,k) = 1.5;
    }}

    //Instantiate fields
    CPL::CPLField f(buf);
    CPL::CPLField field(buf);

    //Hardwire some values
    xi = new double[3];
    int order = 2;
    double rand;
    std::vector<double> zi;
    for (int j = 0; j < 1000; j++ ) {

        for (int ixyz=0; ixyz < 3; ixyz++ ){
            rand = std::rand()/float(RAND_MAX);
            //std::cout << "rand = " << rand << "\n";
            xi[ixyz] = rand*f.dxyz[ixyz];
        }

        std::vector<int> indices = {n};
        zi = f.interpolate(xi, buf, indices, order);
//        std::cout << "pos = " << xi[0] << " "  << xi[1]  << 
//                         " "  << xi[2] << " "  << zi[0] << "\n"; 

        //Get from interpolate function
        ASSERT_DOUBLE_EQ(zi[0]*f.dxyz[0], xi[0]);

        //Get directly from array function
        ASSERT_DOUBLE_EQ(field.get_array_value_interp(indices, xi)[0]*field.dxyz[0], xi[0]);

    }

    delete xi;


}


