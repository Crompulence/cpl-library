#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <vector>

#include <memory>

#include "cpl.h"
#include "CPL_field.h"
#include "CPL_force.h"
#include "CPL_ndArray.h"

// The major difference between a thing that might go wrong and a 
// thing that cannot possibly go wrong is that when a thing that 
// cannot possibly go wrong goes wrong it usually turns out to be
// impossible to get at and repair.
//                                                  Douglas Adams

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

#define threeD for (int ixyz = 0; ixyz<3; ++ixyz) \

#define trplefor_rng(si,sj,sk,Ni,Nj,Nz) for (int i = si; i<Ni; ++i){ \
                                        for (int j = sj; j<Nj; ++j){ \
                                        for (int k = sk; k<Nz; ++k)

#define trplefor(Ni,Nj,Nz) trplefor_rng(0,0,0,Ni,Nj,Nz)

///////////////////////////////////////////////////////////////////
//                                                               //
//                          CPL::ndArray                         //
//                                                               //
///////////////////////////////////////////////////////////////////
//Test for CPL::ndArray - setup and array size 
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

//Test for CPL::ndArray - set elements
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



//Test for CPL::ndArray for parallel processors
// From within testing, we cannot vary MPI topologies
// as would need to be run using (mpiexec -n X)
// it's tricky to test this, The idea of this
// class is it is entirly local to a process and
// so all numbering (and min/max) are processor local
// As a result, parallel testing is not meaningful...

//TEST_F(CPL_Force_Test, test_CPL_define_parallel) {

//    int rankRealm = 0;
//    int nd = 3;
//    int icell = 32;
//    int jcell = 32;
//    int kcell = 32;
//    CPL::ndArray<double> buf;
//    int shape[4] = {nd, icell, jcell, kcell};
//    buf.resize (4, shape);

//    //Test define
//    buf(1,1,1,1) = 5.0;
//    ASSERT_DOUBLE_EQ(buf(1,1,1,1), 5.0);

//    //Test redefine
//    buf(1,1,1,1) = 6.0;
//    ASSERT_NE(buf(1,1,1,1), 5.0);
//};






///////////////////////////////////////////////////////////////////
//                                                               //
//                            CPLField                           //
//                                                               //
///////////////////////////////////////////////////////////////////
//Test for CPL::ndArray - setup and array size 
TEST_F(CPL_Force_Test, test_CPL_field) {
    int nd = 9;
    int icell = 3;
    int jcell = 3;
    int kcell = 3;
    CPL::CPLField field(nd, icell, jcell, kcell);
    CPL::ndArray<double> buf = field.get_array();

    //Test sizes and shapes
    ASSERT_EQ(buf.size(), nd*icell*jcell*kcell);
    ASSERT_EQ(buf.shape(0), nd);
    ASSERT_EQ(buf.shape(1), icell);
    ASSERT_EQ(buf.shape(2), jcell);
    ASSERT_EQ(buf.shape(3), kcell);

    //Test pointer version of code
    CPL::ndArray<double>& buf_ptr = field.get_array_pointer();
    //std::unique_ptr<CPL::ndArray <double>> buf_ptr = field.get_array_pointer();

    ASSERT_EQ(buf_ptr.size(), nd*icell*jcell*kcell);
    ASSERT_EQ(buf_ptr.shape(0), nd);
    ASSERT_EQ(buf_ptr.shape(1), icell);
    ASSERT_EQ(buf_ptr.shape(2), jcell);
    ASSERT_EQ(buf_ptr.shape(3), kcell);

    //Check set_array and getting values using pointer
    CPL::ndArray<double> array;
    int shape[4] = {nd, icell, jcell, kcell};
    array.resize (4, shape);
    array = 0;
    trplefor(icell,jcell,kcell){
        array(0, i, j, k) = 1.0;
    } } }
    field.set_array(array);

    //The copy of the array from get_array has not changed but
    //the pointer to the value in the field object has.
    trplefor(icell,jcell,kcell){
        ASSERT_EQ(buf_ptr(0, i, j, k), array(0, i, j, k));
        ASSERT_NE(buf(0, i, j, k), array(0, i, j, k));
    } } }

};


TEST_F(CPL_Force_Test, test_CPL_field_name) {
    int nd = 9;
    int icell = 3;
    int jcell = 3;
    int kcell = 3;

    std::string fieldstr("field");
    std::string somestrng("string");
    auto somestrng_ptr = std::make_shared<std::string>("string_ptr");

    CPL::ndArray<double> buf;
    int shape[4] = {nd, icell, jcell, kcell};
    buf.resize (4, shape);

    //Check default name is okay
    std::string defaultstr("default");
    CPL::CPLField fieldA(nd, icell, jcell, kcell);
    ASSERT_EQ(defaultstr, fieldA.name);
    auto name = fieldA.field_name();
    ASSERT_EQ(defaultstr, name);

    //Check name is set correctly
    CPL::CPLField fieldB(nd, icell, jcell, kcell, "field");
    ASSERT_EQ(fieldstr, fieldB.name);
    CPL::CPLField fieldC(nd, icell, jcell, kcell, fieldstr);
    ASSERT_EQ(fieldstr, fieldC.name);
    CPL::CPLField fieldD(nd, icell, jcell, kcell, *somestrng_ptr);
    ASSERT_EQ(*somestrng_ptr, fieldD.name);

    //Check other constructors
    CPL::CPLField fieldE(buf, somestrng);
    ASSERT_EQ(somestrng, fieldE.name);
    CPL::CPLField fieldF(nd, buf, somestrng);
    ASSERT_EQ(somestrng, fieldF.name);
    ASSERT_EQ(somestrng, fieldF.field_name());

    CPL::CPLField fieldG(nd, buf, fieldstr);
    ASSERT_NE(somestrng, fieldG.name);

};


//Test for CPL::ndArray - setup and array size 
TEST_F(CPL_Force_Test, test_CPL_field_setters) {

    //Create a field class
    int nd = 3; int icell = 10; int jcell = 10; int kcell = 10;
    CPL::CPLField field(nd, icell, jcell, kcell);

    //Get the pointer to the array
    CPL::ndArray<double>& buf_ptr = field.get_array_pointer();

    //Check example of setting 
    double r[3] = {0.75, 0.25, 0.35};
    double value[3] = {1.0, 0.5, 0.25};
    field.add_to_array(r, value);
    CPL::ndArray<double> buf = field.get_array();

    ASSERT_DOUBLE_EQ(buf(0, 7, 2, 3), value[0]);
    ASSERT_DOUBLE_EQ(buf(1, 7, 2, 3), value[1]);
    ASSERT_DOUBLE_EQ(buf(2, 7, 2, 3), value[2]);

    ASSERT_DOUBLE_EQ(buf_ptr(0, 7, 2, 3), value[0]);
    ASSERT_DOUBLE_EQ(buf_ptr(1, 7, 2, 3), value[1]);
    ASSERT_DOUBLE_EQ(buf_ptr(2, 7, 2, 3), value[2]);

    //Add another set of values
    field.add_to_array(r, value);

    //Check array which we got hasn't changed
    ASSERT_DOUBLE_EQ(buf(0, 7, 2, 3), value[0]);
    ASSERT_DOUBLE_EQ(buf(1, 7, 2, 3), value[1]);
    ASSERT_DOUBLE_EQ(buf(2, 7, 2, 3), value[2]);

    //But the pointer will have
    ASSERT_DOUBLE_EQ(buf_ptr(0, 7, 2, 3), 2.*value[0]);
    ASSERT_DOUBLE_EQ(buf_ptr(1, 7, 2, 3), 2.*value[1]);
    ASSERT_DOUBLE_EQ(buf_ptr(2, 7, 2, 3), 2.*value[2]);

    // Get array again and check (this is all fairly excessive
    // testing but good to be sure memory works as expected)
    buf = field.get_array();
    ASSERT_DOUBLE_EQ(buf(0, 7, 2, 3), 2.*value[0]);
    ASSERT_DOUBLE_EQ(buf(1, 7, 2, 3), 2.*value[1]);
    ASSERT_DOUBLE_EQ(buf(2, 7, 2, 3), 2.*value[2]);

}



////Test for CPL::ndArray - setup and array size 
//TEST_F(CPL_Force_Test, test_CPL_field_setters_olap) {

//    //Create a field class
//    int nd = 3; int icell = 10; int jcell = 10; int kcell = 10;
//    CPL::ndArray<double> buf;
//    int shape[4] = {nd, icell, jcell, kcell};
//    buf.resize (4, shape);

//    //Test define
//    int n = 0;
//    //Slabs in x
//    for (int j = 0; j < jcell; j++ ){
//    for (int k = 0; k < kcell; k++ ){
//        buf(n,0,j,k) = -0.5;
//        buf(n,1,j,k) = 0.5;
//        buf(n,2,j,k) = 1.5;
//    }}
//    CPL::CPLField field(buf);

//    //Hardwire some values
//    double xi[3] = {0.75,0.25,0.35};
//    std::vector<int> indices = {0};
//    std::vector<double> v;

//    //Get directly from array function
//    v = field.get_array_value(indices, xi);
//    ASSERT_DOUBLE_EQ(v[0]*field.dxyz[0], xi[0]);

//    //Random values
//    double rand;
//    for (int j = 0; j < 1000; j++ ) {

//        for (int ixyz=0; ixyz < 3; ixyz++ ){
//            rand = std::rand()/float(RAND_MAX);
//            xi[ixyz] = rand*field.dxyz[ixyz];
//        }

//        //Get directly from array function
//        v = field.get_array_value(indices, xi);
//        ASSERT_DOUBLE_EQ(v[0]*field.dxyz[0], xi[0]);

//    }

//}


//Test for CPL::ndArray - setup and array size 
TEST_F(CPL_Force_Test, test_CPL_field_getters) {

    //Create a field class
    int nd = 3; int icell = 10; int jcell = 10; int kcell = 10;
    CPL::CPLField field(nd, icell, jcell, kcell);

    //Check example of setting
    double r[3] = {0.75, 0.25, 0.35};
    double value[3] = {1.0, 0.5, 0.25};
    field.add_to_array(r, value);

    //Check first with elements of indices vector
    std::vector<int> cell = field.get_cell(r);
    for (int i=0; i<3; i++){
        std::vector<int> indices = {i};
        //Get cell and value from that cell
        ASSERT_DOUBLE_EQ(field.get_array_value(indices, cell[0], cell[1], cell[2])[0],
                       value[i]);
        //Get value directly
        ASSERT_DOUBLE_EQ(field.get_array_value(indices, r)[0], value[i]);
    }

    //Then pass in whole vector
    std::vector<int> indices = {0, 1, 2};
    auto v = field.get_array_value(indices, cell[0], cell[1], cell[2]);
    for (int i=0; i<3; i++){
        ASSERT_DOUBLE_EQ(v[i], value[i]);
    }

    //Get interpolated value
    CPL::CPLField newfield(nd, icell, jcell, kcell);
    trplefor(icell,jcell,kcell){
        newfield.add_to_array(0, i, j, k, value[0]);
        newfield.add_to_array(1, i, j, k, value[1]);
        newfield.add_to_array(2, i, j, k, value[2]);
    } } }
    auto v2 = newfield.get_array_value(indices, cell[0], cell[1], cell[2]);
    for (int i=0; i<3; i++){
        ASSERT_DOUBLE_EQ(v2[i], value[i]);
    }
}


///////////////////////////////////////////////////////////////////
//                                                               //
//                       CPLForceTest                            //
//                                                               //
///////////////////////////////////////////////////////////////////
//Test for CPLForce base class - constructor 
TEST_F(CPL_Force_Test, test_CPL_force_constructor) {

    int nd = 9; int icell = 3; int jcell = 3; int kcell = 3;

    //Call constructor using cell numbers
    CPLForceTest c(nd, icell, jcell, kcell);
    CPL::ndArray<double> buf1 = c.get_field();

    //Call constructor using buf object
    CPL::ndArray<double> buf2;
    int shape[4] = {nd, icell, jcell, kcell};
    buf2.resize(4, shape);
    CPLForceTest d(buf2);

    //Check they give the same values
    ASSERT_EQ(buf1.size(), buf2.size());
    ASSERT_EQ(buf1.shape(0), buf2.shape(0));
    ASSERT_EQ(buf1.shape(1), buf2.shape(1));
    ASSERT_EQ(buf1.shape(2), buf2.shape(2));
    ASSERT_EQ(buf1.shape(3), buf2.shape(3));

}

//Test for CPLForce base class - set field method
TEST_F(CPL_Force_Test, test_CPL_force_get_set_field) {

    int nd = 9; int icell = 3; int jcell = 3; int kcell = 3;

    //Call constructor using cell numbers
    CPLForceTest c(nd, icell, jcell, kcell);
    CPL::ndArray<double> buf1 = c.get_field();
    CPLForceTest d(buf1);

    buf1(4,0,2,1) = 5.;
    c.set_field(buf1);
    d.set_field(buf1);
    CPL::ndArray<double> buf2 = c.get_field();
    ASSERT_DOUBLE_EQ(buf2(4,0,2,1), 5.0);
    buf2 = d.get_field();
    ASSERT_DOUBLE_EQ(buf2(4,0,2,1), 5.0);

}

//Test for CPLForce base class - get field methods
TEST_F(CPL_Force_Test, test_CPL_get_cell) {

    int nd = 3; int icell = 10; int jcell = 10; int kcell = 10;

    //Call constructor using cell numbers
    CPLForceTest c(nd, icell, jcell, kcell);

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

    //Check near edge of domain
    r[0] = 0.9;
    r[1] = 0.5;
    r[2] = 0.5;
    cell = c.get_cell(r);
    ASSERT_EQ(cell[0], 9);
    ASSERT_EQ(cell[1], 5);
    ASSERT_EQ(cell[2], 5);

    //Check out of domain values handled correctly
//    r[0]=-1.; r[1]=0.36; r[2]=0.67;
//    ASSERT_THROW(c.get_cell(r), std::domain_error);
//    try {
//        c.get_cell(r);
//    } catch (std::domain_error& ex) {
//        EXPECT_STREQ("get_cell Error: Input below domain", ex.what());
//    }

//    r[0]=0.5; r[1]=0.36; r[2]=2.3;
//    ASSERT_THROW(c.get_cell(r), std::domain_error);
//    try {
//        c.get_cell(r);
//    } catch (std::domain_error& ex) {
//        EXPECT_STREQ("get_cell Error: Input above domain", ex.what());
//    }
}

//Test for CPLForce base class - get force 
TEST_F(CPL_Force_Test, test_CPL_Force_get_force) {

    //Setup a array which is 1 everywhere
    int nd = 3; int icell = 3; int jcell = 3; int kcell = 3;
    CPL::ndArray<double> array;
    int shape[4] = {nd, icell, jcell, kcell};
    array.resize (4, shape);
    array = 0;
    trplefor(icell,jcell,kcell){
        array(0, i, j, k) = 1.0;
    } } }

    //Create force array object
    CPLForceTest fxyz(array);

    std::vector<double> f(3);
    double r[3];
    double v[3] = {0.5, 0.5, 0.5};
    double a[3] = {0.5, 0.5, 0.5};
    double m=1.; double s=1.; double e=1.;

    //Check array is uniform
    trplefor(icell,jcell,kcell){
        r[0] = i/float(icell);
        r[1] = j/float(jcell);
        r[2] = k/float(kcell);
        f = fxyz.get_force(r, v, a, m, s, e);
        ASSERT_DOUBLE_EQ(f[0], 1.0);
        ASSERT_DOUBLE_EQ(f[1], 0.0);
        ASSERT_DOUBLE_EQ(f[2], 0.0);
    } } }

    //Set array to new value
    trplefor(icell,jcell,kcell){
        array(0, i, j, k) = 0.0;
        array(1, i, j, k) = 6.0;
        array(2, i, j, k) = 2.0;
    } } }

    //Update force array object
    fxyz.set_field(array);

    //Check array is uniform
    trplefor(icell,jcell,kcell){
        r[0] = i/float(icell);
        r[1] = j/float(jcell);
        r[2] = k/float(kcell);
        f = fxyz.get_force(r, v, a, m, s, e);
        ASSERT_DOUBLE_EQ(f[0], 0.0);
        ASSERT_DOUBLE_EQ(f[1], 6.0);
        ASSERT_DOUBLE_EQ(f[2], 2.0);
    } } }


}

//Test for CPLForce base class - constructor 
TEST_F(CPL_Force_Test, test_CPL_force_internalfields) {

    int nd = 9; int icell = 3; int jcell = 3; int kcell = 3;
    CPLForceTest c(nd, icell, jcell, kcell);

    //Get nullptr for field which doesn't exist
    std::string test("not there!");
    auto field = c.get_internal_fields(test);
    EXPECT_TRUE(field == nullptr);

    //Get otherfield internal to CPLForceTest
    std::string name("otherfield");
    auto fieldb = c.get_internal_fields(name);
    EXPECT_TRUE(fieldb != nullptr);
    ASSERT_EQ(name, fieldb->field_name());

    //Check array is same size as array pointer 
    CPL::ndArray<double>& buf_ptr = fieldb->get_array_pointer();
    ASSERT_EQ(buf_ptr.shape(1), icell);
    ASSERT_EQ(buf_ptr.shape(2), jcell);
    ASSERT_EQ(buf_ptr.shape(3), kcell);
    ASSERT_NE(buf_ptr.shape(0), nd);



    //Check array is same size as array pointer 
    CPL::ndArray<double> buf = c.get_field();
    CPLForceTest d(buf);
    //CPLForceTest d(nd, buf);
    auto fieldc = d.get_internal_fields(name);
    CPL::ndArray<double>& buf_ptrc = fieldc->get_array_pointer();
    ASSERT_EQ(buf_ptrc.shape(1), icell);
    ASSERT_EQ(buf_ptrc.shape(2), jcell);
    ASSERT_EQ(buf_ptrc.shape(3), kcell);
    ASSERT_NE(buf_ptrc.shape(0), nd);

};

///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForceVelocity                           //
//                                                               //
///////////////////////////////////////////////////////////////////


//Test for CPLForceVelocity - pre force calculations
TEST_F(CPL_Force_Test, test_velocity_pre_force) {

    //Call constructor using cell numbers
    int nd = 3; int icell = 8; int jcell = 8; int kcell = 8;
    CPLForceVelocity c(nd, icell, jcell, kcell);

    //Default is domain between 0.0 and 1.0
    double r[3] = {0.5, 0.5, 0.5};
    double v[3] = {0.1, 0.2, 0.3};
    double a[3] = {0.0, 0.0, 0.0};
    double m=1.;
    double s=1.;
    double e=1.;

    //Check division to correct locations and binning
    trplefor(icell,jcell,kcell){
        r[0] = i/float(icell);
        r[1] = j/float(jcell);
        r[2] = k/float(kcell);
        c.pre_force(r, v, a, m, s, e);
    } } }

    int sum = 0; double vsum[3] = {0.0,0.0,0.0};
    trplefor(icell,jcell,kcell){

        //Assert single result in each cell and sum to check cells
        ASSERT_EQ(c.nSums(i,j,k), 1);
        sum += c.nSums(i,j,k);

        for (int n=1; n<3 ; n++){
            ASSERT_DOUBLE_EQ(c.vSums(n,i,j,k), (n+1)*0.1);
            vsum[n] += c.vSums(n,i,j,k);

            // Check adjacent vSum in cells in x and z have same value
            if (i > 0)
                ASSERT_DOUBLE_EQ(c.vSums(n,i-1,j,k), c.vSums(n,i,j,k));
            if (j > 0)
                ASSERT_DOUBLE_EQ(c.vSums(n,i,j-1,k), c.vSums(n,i,j,k));
            if (k > 0)
                ASSERT_DOUBLE_EQ(c.vSums(n,i,j,k-1), c.vSums(n,i,j,k));
        }
    } } }
    //Check that sum of results is equal to Number of cells
    ASSERT_EQ(sum, icell*jcell*kcell);
    for (int n=1; n<3 ; n++){
       //Note these are not equal to double precision
       ASSERT_FLOAT_EQ(vsum[n], (n+1)*0.1*icell*jcell*kcell);
       EXPECT_LT(vsum[n]-(n+1)*0.1*icell*jcell*kcell,1e-8);
    }
    
    //Check that reset works
    c.resetsums();   sum = 0;
    vsum[0]=0.0; vsum[1]=0.0; vsum[2]=0.0;
    trplefor(icell,jcell,kcell){
        sum += c.nSums(i,j,k);
        for (int n=1; n<3 ; n++)
            vsum[n] += c.vSums(n,i,j,k);
    } } }
    ASSERT_EQ(sum, 0);
    for (int n=1; n<3 ; n++){
       ASSERT_DOUBLE_EQ(vsum[n], 0.0);
    }
}


//Test for CPLForceFlekkoy - get force
TEST_F(CPL_Force_Test, test_velocity_get_force) {

    //Call constructor using cell numbers
    int nd = 3; int icell = 8; int jcell = 8; int kcell = 8;

    //Setup a array which is 1 everywhere
    CPL::ndArray<double> array;
    int shape[4] = {nd, icell, jcell, kcell};
    array.resize (4, shape);
    trplefor(icell,jcell,kcell){
        array(0, i, j, k) = 1.0;
        array(1, i, j, k) = 2.0;
        array(2, i, j, k) = 3.0;
    } } }

    //Create force array object
    CPLForceVelocity fxyz(array);

    //Setup one particle per cell
    //and as zero no change
    double r[3] = {0.0, 0.0, 0.0};
    double v[3] = {0.0, 0.0, 0.0};
    double a[3] = {0.0, 0.0, 0.0};
    double m=1.; double s=1.; double e=1.;

    trplefor(icell,jcell,kcell){
        r[0] = i/float(icell);
        r[1] = j/float(jcell);
        r[2] = k/float(kcell);
        fxyz.pre_force(r, v, a, m, s, e);
    } } }

    std::vector<double> f(3);
    std::vector<int> cell(3);
    trplefor(icell,jcell,kcell){
        r[0] = i/float(icell);
        r[1] = j/float(jcell);
        r[2] = k/float(kcell);
        f = fxyz.get_force(r, v, a, m, s, e);
        ASSERT_DOUBLE_EQ(f[0], 1.0);
        ASSERT_DOUBLE_EQ(f[1], 2.0);
        ASSERT_DOUBLE_EQ(f[2], 3.0);
    } } }

    //Reset sums and check with velocity 
    //already at mean value of non-zero
    fxyz.resetsums();
    trplefor(icell,jcell,kcell){
        r[0] = i/float(icell);
        r[1] = j/float(jcell);
        r[2] = k/float(kcell);
        v[0] = 1.0; v[1] = 2.0; v[2] = 3.0;
        fxyz.pre_force(r, v, a, m, s, e);
    } } }

    trplefor(icell,jcell,kcell){
        r[0] = i/float(icell);
        r[1] = j/float(jcell);
        r[2] = k/float(kcell);
        f = fxyz.get_force(r, v, a, m, s, e);
        ASSERT_DOUBLE_EQ(f[0], 0.0);
        ASSERT_DOUBLE_EQ(f[1], 0.0);
        ASSERT_DOUBLE_EQ(f[2], 0.0);
    } } }

}



///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForceFlekkoy                            //
//                                                               //
///////////////////////////////////////////////////////////////////

//Test for CPLForceFlekkoy - constructor and fields
TEST_F(CPL_Force_Test, test_flekkoy_CPL_inhereted) {

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

//Test for CPLForceFlekkoy - flekkoyGWeight method
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

//Test for CPLForceFlekkoy - pre force calculations
TEST_F(CPL_Force_Test, test_flekkoy_pre_force) {

    //Call constructor using cell numbers
    int nd = 1; int icell = 8; int jcell = 8; int kcell = 8;

    CPLForceFlekkoy c(nd, icell, jcell, kcell);

    //Default is domain between 0.0 and 1.0
    double r[3] = {0.5, 0.5, 0.5};
    double v[3] = {0.0, 0.0, 0.0};
    double a[3] = {0.0, 0.0, 0.0};
    double m=1.; double s=1.; double e=1.;

    //Check division to correct locations and binning
    trplefor(icell,jcell,kcell){
        r[0] = i/float(icell);
        r[1] = j/float(jcell);
        r[2] = k/float(kcell);
        c.pre_force(r, v, a, m, s, e);
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

//Test for CPLForceFlekkoy - pre force calculations with non-uniform domain
TEST_F(CPL_Force_Test, test_flekkoy_pre_force_varydomain) {

    //Call constructor using cell numbers
    int nd = 1; int icell = 8; int jcell = 8; int kcell = 8;
    CPLForceFlekkoy c(nd, icell, jcell, kcell);

    //Default is domain between 0.0 and 1.0
    double r[3] = {0.0, 0.0, 0.0};
    double v[3] = {0.0, 0.0, 0.0};
    double a[3] = {0.0, 0.0, 0.0};
    double m=1.; double s=1.; double e=1.;

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
        c.pre_force(r, v, a, m, s, e);
    } } }
    trplefor(icell,jcell,kcell){
        ASSERT_EQ(c.nSums(i,j,k), 1);
    } } }

}

//Test for CPLForceFlekkoy - get force
TEST_F(CPL_Force_Test, test_flekkoy_get_force) {

    //Call constructor using cell numbers
    int nd = 9; int icell = 8; int jcell = 8; int kcell = 8;
    double m=1.; double s=1.; double e=1.;

    //Setup a field which is 1 everywhere
    CPL::ndArray<double> field;
    int shape[4] = {nd, icell, jcell, kcell};
    field.resize (4, shape);
    trplefor(icell,jcell,kcell){
        field(1, i, j, k) = 1.0;
    } } }

    //Create force field object
    CPLForceFlekkoy fxyz(field);

    //Setup one particle per cell
    double r[3] = {0.0, 0.0, 0.0};
    double v[3] = {0.0, 0.0, 0.0};
    double a[3] = {0.0, 0.0, 0.0};
    trplefor(icell,jcell,kcell){
        r[0] = i/float(icell);
        r[1] = j/float(jcell);
        r[2] = k/float(kcell);
        fxyz.pre_force(r, v, a, m, s, e);
    } } }

    std::vector<double> f(3);
    std::vector<double> dA = fxyz.get_dA(); 
    trplefor(icell,jcell,kcell){
        r[0] = i/float(icell);
        r[1] = j/float(jcell);
        r[2] = k/float(kcell);
        f = fxyz.get_force(r, v, a, m, s, e);
        if (fxyz.gSums(i,j,k) != 0.){
            ASSERT_DOUBLE_EQ(f[0], dA[1]);
        } else {
            ASSERT_DOUBLE_EQ(f[0], 0.0);
           
        }
        ASSERT_DOUBLE_EQ(f[1], 0.0);
        ASSERT_DOUBLE_EQ(f[2], 0.0);
    } } }

    //Change field to gradient
    //Set gsums and nsums to iunity throughout
    fxyz.resetsums();
    trplefor(icell,jcell,kcell){
        field(1, i, j, k) = float(i);
        field(4, i, j, k) = float(j);
        field(7, i, j, k) = float(k);
        fxyz.nSums(i,j,k) = 1.0;
        fxyz.gSums(i,j,k) = 1.0;
    } } }
    //Update force field object
    fxyz.set_field(field);

    //check force at a range of random positions  
    double g;
    std::vector<int> cell;
    for (int i = 0; i<10000; ++i){
        double x = std::rand()/float(RAND_MAX);
        double y = std::rand()/float(RAND_MAX);
        double z = std::rand()/float(RAND_MAX);
        r[0] = x; r[1] = y; r[2] = z;
        f = fxyz.get_force(r, v, a, m, s, e);
        g = fxyz.flekkoyGWeight(r[1], 0.0, 1.0);
        cell = fxyz.get_cell(r);
        ASSERT_DOUBLE_EQ(f[0], g*float(cell[0])*dA[1]);
        ASSERT_DOUBLE_EQ(f[1], g*float(cell[1])*dA[1]);
        ASSERT_DOUBLE_EQ(f[2], g*float(cell[2])*dA[1]);
    }
}




///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForceDrag                               //
//                                                               //
///////////////////////////////////////////////////////////////////

//Test for CPLForceDrag - constructor and fields
TEST_F(CPL_Force_Test, test_CPLForce_Drag) {

    int nd = 3; int icell = 3; int jcell = 3; int kcell = 3;

    //Call constructor using cell numbers
    CPLForceDrag c(nd, icell, jcell, kcell);
    CPL::ndArray<double> buf1 = c.get_field();
    CPLForceDrag d(buf1);

    buf1(2,0,2,1) = 5.;
    c.set_field(buf1);
    d.set_field(buf1);
    CPL::ndArray<double> buf2 = c.get_field();
    ASSERT_DOUBLE_EQ(buf2(2,0,2,1), 5.0);
    buf2 = d.get_field();
    ASSERT_DOUBLE_EQ(buf2(2,0,2,1), 5.0);

}


//Test for CPLForceDrag - check initial esum and Fsum arrays
TEST_F(CPL_Force_Test, test_CPLForce_Drag_initial_eSumsFsum) {

    int nd = 3; int icell = 3; int jcell = 3; int kcell = 3;

    //Call constructor using cell numbers
    CPLForceDrag c(nd, icell, jcell, kcell);

    CPL::ndArray<double> ebuf = c.eSums->get_array();
    CPL::ndArray<double> Fbuf = c.FSums->get_array();

    //Try to get porosity from force type
    trplefor(icell,jcell,kcell){
        ASSERT_DOUBLE_EQ(ebuf(0,i,j,k), 0.0);
        ASSERT_DOUBLE_EQ(Fbuf(0,i,j,k), 0.0);
        ASSERT_DOUBLE_EQ(Fbuf(1,i,j,k), 0.0);
        ASSERT_DOUBLE_EQ(Fbuf(2,i,j,k), 0.0);
    } } }

}



//Test for CPLForceDrag - check sum of esum and Fsum arrays
TEST_F(CPL_Force_Test, test_CPLForce_Drag_check_eSumsFsum) {

    int nd = 9; int icell = 3; int jcell = 3; int kcell = 3;

    //Call constructor using cell numbers
    CPLForceDrag c(nd, icell, jcell, kcell);

    //Setup one particle per cell
    double r[3] = {0.0, 0.0, 0.0};
    double v[3] = {0.0, 0.0, 0.0};
    double a[3] = {0.0, 0.0, 0.0};
    std::vector<double> F, Ur;
    double radius = 0.001;
    double volume = (4./3.)*M_PI*pow(radius,3);
    double m=1.; double s=radius; double e=1.;
    //For test force, drag coefficient is a constant for all input
    double Cd = c.drag_coefficient(r, s*2., Ur);
    //Test with PCM (volume assigned to cell based on centre)
    //c.use_overlap = false;

    //Setup esum
    trplefor(icell,jcell,kcell){
        r[0] = i/double(icell);
        r[1] = j/double(jcell);
        r[2] = k/double(kcell);
        c.pre_force(r, v, a, m, s, e);
    } } }

    //Setup Fsum for u=0
    trplefor(icell,jcell,kcell){
        r[0] = i/double(icell);
        r[1] = j/double(jcell);
        r[2] = k/double(kcell);
        v[0] = i/double(icell);
        v[1] = j/double(jcell);
        v[2] = k/double(kcell);
        F = c.get_force(r, v, a, m, s, e);
    } } }

    //Check values of esum & Fsum
    CPL::ndArray<double> ebuf = c.eSums->get_array();
    CPL::ndArray<double> Fbuf = c.FSums->get_array();
    trplefor(icell,jcell,kcell){
        //Since introducing overlap mode, need assert near
        //ASSERT_NEAR(c.eSums(i,j,k), volume, 1e-14);
        ASSERT_DOUBLE_EQ(ebuf(0,i,j,k), volume);
        ASSERT_DOUBLE_EQ(Fbuf(0,i,j,k), -Cd*i/double(icell));
        ASSERT_DOUBLE_EQ(Fbuf(1,i,j,k), -Cd*j/double(jcell));
        ASSERT_DOUBLE_EQ(Fbuf(2,i,j,k), -Cd*k/double(kcell));
    } } }


    //Setup Fsum for u=0, v=1 and w=0
    double VCFD = 1.0;
    int shape[4] = {nd, icell, jcell, kcell};
    CPL::ndArray<double> field;
    field.resize (4, shape);
    trplefor(icell,jcell,kcell){
        field(1, i, j, k) = VCFD;
    } } }

    //Force here should then be based on velocity of U=1.0
    CPLForceDrag d(field);
    trplefor(icell,jcell,kcell){
        r[0] = i/double(icell);
        r[1] = j/double(jcell);
        r[2] = k/double(kcell);
        v[0] = i/double(icell);
        v[1] = j/double(jcell);
        v[2] = k/double(kcell);
        F = d.get_force(r, v, a, m, s, e);
    } } }


    //Check values of esum & Fsum
    ebuf = d.eSums->get_array();
    Fbuf = d.FSums->get_array();
    trplefor(icell,jcell,kcell){
        ASSERT_DOUBLE_EQ(ebuf(0,i,j,k), 0.);
        ASSERT_NEAR(Fbuf(0,i,j,k), Cd*(0.-i/double(icell)),1e-13);
        ASSERT_NEAR(Fbuf(1,i,j,k), Cd*(VCFD-j/double(jcell)),1e-13);
        ASSERT_NEAR(Fbuf(2,i,j,k), Cd*(0.-k/double(kcell)),1e-13);
    } } }

    //Check force based on pressure gradient
    //Setup field with dPdz = 1.0
    double dPdz = 1.0;
    CPL::ndArray<double> newfield;
    newfield.resize (4, shape);
    trplefor(icell,jcell,kcell){
        newfield(5, i, j, k) = dPdz;
    } } }

    CPLForceDrag f(newfield);
    trplefor(icell,jcell,kcell){
        r[0] = i/double(icell);
        r[1] = j/double(jcell);
        r[2] = k/double(kcell);
        v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;
        F = f.get_force(r, v, a, m, s, e);
    } } }

    //Check values of esum & Fsum
    Fbuf = f.FSums->get_array();
    trplefor(icell,jcell,kcell){
        ASSERT_DOUBLE_EQ(Fbuf(0,i,j,k), 0.0);
        ASSERT_DOUBLE_EQ(Fbuf(1,i,j,k), 0.0);
        ASSERT_DOUBLE_EQ(Fbuf(2,i,j,k), -volume*dPdz);

    } } }


}


//Test for CPL::field - overlap fraction of particle
TEST_F(CPL_Force_Test, test_CPLForce_Drag_overlap) {
    int nd = 9;
    int icell = 3;
    int jcell = 3;
    int kcell = 3;
    CPL::CPLField field(nd, icell, jcell, kcell);
    CPL::ndArray<double> buf = field.get_array();

    //Particle in centre of cell
    double scx = 0.0;double scy=0.0; double scz=0.0; double sr = 1;
    double xb = -1.0; double yb=-1.0; double zb=-1.0;
    double xt =  1.0; double yt= 1.0; double zt= 1.0;
    double Vs = (4./3.)*M_PI*pow(sr,3);
    //double Vc = (xt-xb)*(yt-yb)*(zt-zb);

    //Particle at centre of clel
    double result = field.sphere_cube_overlap(scx, scy, scz, sr, xb, yb, zb, xt, yt, zt);
    ASSERT_NEAR(result, Vs, 1e-14);

    //Particle at face of cell
    scx = 1.0;
    result = field.sphere_cube_overlap(scx, scy, scz, sr, xb, yb, zb, xt, yt, zt);
    ASSERT_NEAR(result, 0.5*Vs, 1e-14);

    //Particle at top edge of cell
    scx = 1.0; scy = 1.0;
    result = field.sphere_cube_overlap(scx, scy, scz, sr, xb, yb, zb, xt, yt, zt);
    ASSERT_NEAR(result, 0.25*Vs, 1e-14);

    //Particle at corner of cell
    scx = 1.0; scy = 1.0; scz = 1.0;
    result = field.sphere_cube_overlap(scx, scy, scz, sr, xb, yb, zb, xt, yt, zt);
    ASSERT_NEAR(result, 0.125*Vs, 1e-14);

};



//Test for CPLForceDrag - check sum of esum and Fsum arrays
TEST_F(CPL_Force_Test, test_CPLForce_Drag_check_overlap_field) {

    int nd = 9; int icell = 3; int jcell = 3; int kcell = 3;

    //Setup one particle per cell
    double v[3] = {0.0, 0.0, 0.0};
    double a[3] = {0.0, 0.0, 0.0};
    double radius = 0.5;
    double min[3] = {0.0, 0.0, 0.0};
    double max[3] = {3.0, 3.0, 3.0};

    double Vs = (4./3.)*M_PI*pow(radius,3);
    double m=1.; double s=radius; double e=1.;

    //Call constructor using cell numbers
    std::map <std::string, std::string> args_map
    {
        { "use_overlap", "1" },
        { "use_interpolate", "0" },
        { "Cd", "0.0001" }
    };

    CPLForceDrag c(nd, icell, jcell, kcell, args_map);
    c.set_minmax(min, max);

    //Get array pointers
    CPL::ndArray<double>& ebuf = c.eSums->get_array_pointer();

    //Particle at centre of grid cell 2
    double r[3] = {1.5, 1.5, 1.5};
    c.pre_force(r, v, a, m, s, e);
    ASSERT_NEAR(ebuf(0,1,1,1), Vs, 1e-13);

    //Particle at face of grid cell 2 and 3 in x
    ebuf(0,1,1,1) = 0.0;
    r[0] = 2.0;
    c.pre_force(r, v, a, m, s, e);
    ASSERT_NEAR(ebuf(0,1,1,1), 0.5*Vs, 1e-13);
    ASSERT_NEAR(ebuf(0,2,1,1), 0.5*Vs, 1e-13);

    //Particle at edge of grid cell 2 and 3 in x,y
    ebuf(0,1,1,1) = 0.0; ebuf(0,2,1,1) = 0.0;
    r[0] = 2.0;  r[1] = 2.0;
    c.pre_force(r, v, a, m, s, e);
    ASSERT_NEAR(ebuf(0,1,1,1), 0.25*Vs, 1e-13);
    ASSERT_NEAR(ebuf(0,2,1,1), 0.25*Vs, 1e-13);
    ASSERT_NEAR(ebuf(0,1,2,1), 0.25*Vs, 1e-13);
    ASSERT_NEAR(ebuf(0,2,2,1), 0.25*Vs, 1e-13);

    //Particle at corner of grid cell 2 and 3 in x,y,z
    ebuf(0,1,1,1) = 0.0; ebuf(0,2,1,1) = 0.0;
    ebuf(0,1,2,1) = 0.0; ebuf(0,2,2,1) = 0.0;
    r[0] = 2.0;  r[1] = 2.0;  r[2] = 2.0;
    c.pre_force(r, v, a, m, s, e);
    ASSERT_NEAR(ebuf(0,1,1,1), 0.125*Vs, 1e-13);
    ASSERT_NEAR(ebuf(0,2,1,1), 0.125*Vs, 1e-13);
    ASSERT_NEAR(ebuf(0,1,2,1), 0.125*Vs, 1e-13);
    ASSERT_NEAR(ebuf(0,2,2,1), 0.125*Vs, 1e-13);
    ASSERT_NEAR(ebuf(0,1,1,2), 0.125*Vs, 1e-13);
    ASSERT_NEAR(ebuf(0,2,1,2), 0.125*Vs, 1e-13);
    ASSERT_NEAR(ebuf(0,1,2,2), 0.125*Vs, 1e-13);
    ASSERT_NEAR(ebuf(0,2,2,2), 0.125*Vs, 1e-13);

    //Particle bigger than cell so fills whole cell
    s = 1.0;
    r[0] =1.5;  r[1] = 1.5;  r[2] = 1.5;
    ebuf(0,1,1,1) = 0.0;
    c.pre_force(r, v, a, m, s, e);
    ASSERT_NEAR(ebuf(0,1,1,1), 1.0, 1e-13);

    //Domain is from 0 to 1 so cell size is 0.1
    icell = 10; jcell = 10; kcell = 10;
    std::map <std::string, std::string> args_map_
    {
        { "use_overlap", "1" },
        { "use_interpolate", "0" },
        { "Cd", "0.00001" }
    };
    CPLForceDrag d(nd, icell, jcell, kcell, args_map_);
    CPL::ndArray<double>& dbuf = d.eSums->get_array_pointer();

    //Particle bigger than lots of cell so fills many cell
    s = 0.35;
    r[0] =0.5;  r[1] = 0.5;  r[2] = 0.5;
    d.pre_force(r, v, a, m, s, e);
    double Vc = 0.1*0.1*0.1;

    trplefor_rng(3,3,3,7,7,7){
        ASSERT_NEAR(dbuf(0,i,j,k), Vc, 1e-13);
    } } }

//    CPLForceDrag f(nd, icell, jcell, kcell, true);

//    s = 0.1;
//    r[0] =0.5;  r[1] = 0.5;  r[2] = 0.5;
//    f.pre_force(r, v, a, m, s, e);

//    for (int i=0; i<10; i++){
//        
//    }

}

///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForceGranular                            //
//                                                               //
///////////////////////////////////////////////////////////////////

//Test for CPLForceGranular - constructor and fields
TEST_F(CPL_Force_Test, test_Granular_CPL_inhereted) {

    int nd = 9; int icell = 3; int jcell = 3; int kcell = 3;

    //Call constructor using cell numbers
    CPLForceGranular c(nd, icell, jcell, kcell);
    CPL::ndArray<double> buf1 = c.get_field();
    CPLForceGranular d(buf1);

    buf1(2,0,2,1) = 5.;
    c.set_field(buf1);
    d.set_field(buf1);
    CPL::ndArray<double> buf2 = c.get_field();
    ASSERT_DOUBLE_EQ(buf2(2,0,2,1), 5.0);
    buf2 = d.get_field();
    ASSERT_DOUBLE_EQ(buf2(2,0,2,1), 5.0);

}


//Test for CPLForceGranular - test forces
TEST_F(CPL_Force_Test, test_Granular_CPL_forces) {

    int nd = 9; int icell = 9; int jcell = 9; int kcell = 9;

    //Setup Fsum for u=1, v=0 and w=0
    CPL::ndArray<double> field;
    int shape[4] = {nd, icell, jcell, kcell};
    field.resize (4, shape);
    trplefor(icell,jcell,kcell){
        field(0, i, j, k) = 1.0;
    } } }

    //Call constructor using cell numbers
    CPLForceGranular c(field);

    //Setup one particle per cell
    double r[3] = {0.0, 0.0, 0.0};
    double v[3] = {0.0, 0.0, 0.0};
    double a[3] = {0.0, 0.0, 0.0};
    std::vector<double> F;
    double radius = 0.01;
    double volume = (4./3.)*M_PI*pow(radius,3);
    double m=1.;
    double s=radius;
    double e=1.;

    //Setup esum
    trplefor(icell,jcell,kcell){
        r[0] = i/double(icell);
        r[1] = j/double(jcell);
        r[2] = k/double(kcell);
        c.pre_force(r, v, a, m, s, e);
    } } }

    trplefor(icell,jcell,kcell){
        r[0] = i/double(icell);
        r[1] = j/double(jcell);
        r[2] = k/double(kcell);
        v[0] = i/double(icell);
        v[1] = j/double(jcell);
        v[2] = k/double(kcell);
        F = c.get_force(r, v, a, m, s, e);
    } } }


    //Get array pointers
    std::string name("eSums");
    auto field_ptr = c.get_internal_fields(name);

    trplefor(icell,jcell,kcell){

        if (field_ptr != nullptr){
            ASSERT_DOUBLE_EQ(field_ptr->get_array_value(0, i, j, k), volume);
        }


        
//        std::cout << i << " " << j << " " << k
//                    << " " << c.eSums(i, j, k)
//                    << " " << c.FSums(0, i, j, k)
//                    << " " << c.FSums(1, i, j, k)
//                    << " " << c.FSums(2, i, j, k) << std::endl;
    } } }


}

