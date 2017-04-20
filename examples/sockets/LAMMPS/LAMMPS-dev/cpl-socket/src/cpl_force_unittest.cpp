#include "cpl_force.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <vector>
#include "CPL.h"

// The fixture for testing class Foo.
class CPL_Force_Test : public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  CPL_Force_Test() {
    // You can do set-up work for each test here.
  }

  virtual ~CPL_Force_Test() {
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

TEST_F(CPL_Force_Test, test_CPL_array_size) {
    int nd = 9;
    int icell = 3;
    int jcell = 3;
    int kcell = 3;
    CPL::ndArray<double> buf;
    int shape[4] = {nd, icell, jcell, kcell};
    buf.resize (4, shape);

    //Test sizes and shapes
    ASSERT_EQ(buf.size(), nd*icell*jcell*kcell);
    ASSERT_EQ(buf.shape(0), nd);
    ASSERT_EQ(buf.shape(1), icell);
    ASSERT_EQ(buf.shape(2), jcell);
    ASSERT_EQ(buf.shape(3), kcell);

};


TEST_F(CPL_Force_Test, test_CPL_define) {
    int nd = 9;
    int icell = 3;
    int jcell = 3;
    int kcell = 3;
    CPL::ndArray<double> buf;
    int shape[4] = {nd, icell, jcell, kcell};
    buf.resize (4, shape);

    //Test define
    buf(1,1,1,1) = 5.0;
    ASSERT_DOUBLE_EQ(buf(1,1,1,1), 5.0);

    //Test redefine
    buf(1,1,1,1) = 6.0;
    ASSERT_NE(buf(1,1,1,1), 5.0);
};


TEST_F(CPL_Force_Test, test_CPL_force_constructor) {

    int nd = 9; int icell = 3; int jcell = 3; int kcell = 3;

    //Call constructor using cell numbers
    CPLForce c(nd, icell, jcell, kcell);
    CPL::ndArray<double> buf1 = c.get_field();

    //Call constructor using buf object
    CPL::ndArray<double> buf2;
    int shape[4] = {nd, icell, jcell, kcell};
    buf2.resize(4, shape);
    CPLForce d(buf2);

    //Check they give the same values
    ASSERT_EQ(buf1.size(), buf2.size());
    ASSERT_EQ(buf1.shape(0), buf2.shape(0));
    ASSERT_EQ(buf1.shape(1), buf2.shape(1));
    ASSERT_EQ(buf1.shape(2), buf2.shape(2));
    ASSERT_EQ(buf1.shape(3), buf2.shape(3));


};

TEST_F(CPL_Force_Test, test_CPL_force_get_set_field) {

    int nd = 9; int icell = 3; int jcell = 3; int kcell = 3;

    //Call constructor using cell numbers
    CPLForce c(nd, icell, jcell, kcell);
    CPL::ndArray<double> buf1 = c.get_field();
    CPLForce d(buf1);

    buf1(4,0,2,1) = 5.;
    c.set_field(buf1);
    d.set_field(buf1);
    CPL::ndArray<double> buf2 = c.get_field();
    ASSERT_DOUBLE_EQ(buf2(4,0,2,1), 5.0);
    buf2 = d.get_field();
    ASSERT_DOUBLE_EQ(buf2(4,0,2,1), 5.0);

}


TEST_F(CPL_Force_Test, test_CPL_get_cell) {

    int nd = 3; int icell = 10; int jcell = 10; int kcell = 10;

    //Call constructor using cell numbers
    CPLForce c(nd, icell, jcell, kcell);

    //Check simple example of three cells
    double r[3] = {0.75,0.25,0.35};
    std::vector<int> cell = c.get_cell(r);
    ASSERT_EQ(cell[0], 7);
    ASSERT_EQ(cell[1], 2);
    ASSERT_EQ(cell[2], 3);

    // Check edge case for 0.3 and 0.6 which 
    // because of floor give unexpected results
    r[0]=0.000000000001; r[1]=0.3; r[2]=0.6;
    cell = c.get_cell(r);
    EXPECT_EQ(cell[0], 0);
    EXPECT_EQ(cell[1], 2);
    EXPECT_EQ(cell[2], 5);

    //Check out of domain values handled correctly
    r[0]=-1.; r[1]=0.36; r[2]=0.67;
    ASSERT_THROW(c.get_cell(r), std::domain_error);
    try {
        c.get_cell(r);
    } catch (std::domain_error& ex) {
        EXPECT_STREQ("get_cell Error: Input below domain", ex.what());
    }

    r[0]=0.5; r[1]=0.36; r[2]=2.3;
    ASSERT_THROW(c.get_cell(r), std::domain_error);
    try {
        c.get_cell(r);
    } catch (std::domain_error& ex) {
        EXPECT_STREQ("get_cell Error: Input above domain", ex.what());
    }
}


TEST_F(CPL_Force_Test, test_CPL_inhereted) {

    int nd = 9; int icell = 3; int jcell = 3; int kcell = 3;

    //Call constructor using cell numbers
    CPLForceFlekkoy c(nd, icell, jcell, kcell);
    CPL::ndArray<double> buf1 = c.get_field();
    CPLForceFlekkoy d(buf1);

    buf1(4,0,2,1) = 5.;
    c.set_field(buf1);
    d.set_field(buf1);
    CPL::ndArray<double> buf2 = c.get_field();
    ASSERT_DOUBLE_EQ(buf2(4,0,2,1), 5.0);
    buf2 = d.get_field();
    ASSERT_DOUBLE_EQ(buf2(4,0,2,1), 5.0);

}

TEST_F(CPL_Force_Test, test_flekkoyGWeight) {

    int nd = 1; int icell = 8; int jcell = 8; int kcell = 8;

    //Call constructor using cell numbers
    CPLForceFlekkoy c(nd, icell, jcell, kcell);

    // Some basic tests for flekkoyGWeight here?
//    for (int i = 0; i<100; ++i){
//         ASSERT_DOUBLE_EQ(c.flekkoyGWeight(i/float(100), 0.0, 1.0), g);
//    }

    // Error test is flekkoyGWeight outside range but I can't get this working (need header from google tests in force)
    ASSERT_THROW(c.flekkoyGWeight(1.1, 0.0, 1.0), std::domain_error);
    try {
        c.flekkoyGWeight(1.1, 0.0, 1.0);
    } catch (std::domain_error& ex) {
        EXPECT_STREQ("flekkoyGWeight Error: Position argument y greater than ymin", ex.what());
    }
    ASSERT_DOUBLE_EQ(c.flekkoyGWeight(-0.1, 0.0, 1.0), 0.0);

}


#define trplefor(Ni,Nj,Nz) for (int i = 0; i<Ni; ++i){ \
                           for (int j = 0; j<Nj; ++j){ \
                           for (int k = 0; k<Nz; ++k)

TEST_F(CPL_Force_Test, test_pre_force) {

    //Call constructor using cell numbers
    int nd = 1; int icell = 8; int jcell = 8; int kcell = 8;
    CPLForceFlekkoy c(nd, icell, jcell, kcell);

    //Default is domain between 0.0 and 1.0
    double r[3] = {0.5, 0.5, 0.5};
    double v[3] = {0.0, 0.0, 0.0};
    double a[3] = {0.0, 0.0, 0.0};

    //Check division to correct locations and binning
    trplefor(icell,jcell,kcell){
        r[0] = i/float(icell);
        r[1] = j/float(jcell);
        r[2] = k/float(kcell);
        c.pre_force(r, v, a);
    } } }

    int sum = 0;
    trplefor(icell,jcell,kcell){
        //Assert single result in each cell
        ASSERT_EQ(c.nSums(i,j,k), 1);
        //Add to sum to check this all cells correct
        sum += c.nSums(i,j,k);
        // Check adjacent gSum in cells in x and z have same value
        // as only y is expected to be changing
        if (i > 0)
            ASSERT_DOUBLE_EQ(c.gSums(i-1,j,k), c.gSums(i,j,k));
        if (k > 0)
            ASSERT_DOUBLE_EQ(c.gSums(i,j,k-1), c.gSums(i,j,k));
    } } }
    //Check that sum of results is equal to Number of cells
    ASSERT_EQ(sum, icell*jcell*kcell);

    //Check that reset works
    c.resetsums();   sum = 0;
    trplefor(icell,jcell,kcell){
        sum += c.nSums(i,j,k);
    } } }
    ASSERT_EQ(sum, 0);


}

TEST_F(CPL_Force_Test, test_pre_force_varydomain) {

    //Call constructor using cell numbers
    int nd = 1; int icell = 8; int jcell = 8; int kcell = 8;
    CPLForceFlekkoy c(nd, icell, jcell, kcell);

    //Default is domain between 0.0 and 1.0
    double r[3] = {0.0, 0.0, 0.0};
    double v[3] = {0.0, 0.0, 0.0};
    double a[3] = {0.0, 0.0, 0.0};

    //Adjust limits and check
    double min[3] = {-3.0, -2.0, -5.0};
    double max[3] = {3.0, 2.0, 5.0};
    double range[3];
    for (int i = 0; i<3; ++i){
        range[i] = max[i]-min[i];
    }
    c.set_minmax(min, max);

    //Check division to correct locations
    trplefor(icell,jcell,kcell){
        r[0] = (i*range[0])/float(icell) + min[0];
        r[1] = (j*range[1])/float(jcell) + min[1];
        r[2] = (k*range[2])/float(kcell) + min[2];
        c.pre_force(r, v, a);
    } } }
    trplefor(icell,jcell,kcell){
        ASSERT_EQ(c.nSums(i,j,k), 1);
    } } }

}


TEST_F(CPL_Force_Test, test_get_force) {

    //Call constructor using cell numbers
    int nd = 9; int icell = 8; int jcell = 8; int kcell = 8;
    CPLForceFlekkoy c(nd, icell, jcell, kcell);


}





