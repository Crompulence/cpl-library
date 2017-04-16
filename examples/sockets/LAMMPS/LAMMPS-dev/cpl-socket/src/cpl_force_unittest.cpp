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

// Tests that the Foo::Bar() method returns one.
//TEST_F(CPL_Force_Test, CheckField) {
//    int ncells=10;
//    CPLForce c(1, ncells, 1, 1);
//    int size = ncells;
//    std::vector<double> array1(size);
//    // Add increasing numbers
//    for(int i=0; i<size; ++i){
//      array1[i] = i;
//    }
//    c.set_field(array1);
//    std::vector<double> array2 = c.get_field();
//    ASSERT_THAT(array2, testing::ElementsAreArray(array1));
//};


//TEST_F(CPL_Force_Test, CellCheck) {
//    int ncells=10;
//    CPLForce c(1, ncells, 1, 1);
//    int size = ncells;
//    std::vector<double> array1(size);
//    // Add increasing numbers
//    for(int i=0; i<size; ++i){
//      array1[i] = i;
//    }
//    c.set_field(array1);
//    ASSERT_THAT(c.get_cell(5), array1[5]);
//};


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






