
//=====================================================
// Class of type CPLField
//    In the interest of dependency injection, a field 
//    class with dxyz, min, max, interpolation and get 
//    cell. Build with CPL_ndArray as dependency and 
//    use as input to constructor of CPLField.
//=====================================================

#ifndef CPLField_H_
#define CPLField_H_

#include <vector>
#include <iostream>
#include <utility>
#include <memory>

#include "CPL_ndArray.h"

namespace CPL{

class CPLField{

public:

    //Constructors
    CPLField(int nd, int icell, int jcell, int kcell);
    CPLField(CPL::ndArray<double> arrayin);

    //Getters and setters
    void set_array(CPL::ndArray<double> arrayin);
    void set_minmax(double min_in[], double max_in[]);
    void set_dxyz();

    CPL::ndArray<double> get_array();
    //std::unique_ptr<CPL::ndArray <double>> get_array_pointer();
    CPL::ndArray<double>& get_array_pointer();

    //Get cell values
    std::vector<int> get_cell(double r[]);
    std::vector<double> get_dA();    
    std::vector<double> interpolate(double r[]);    

    //Function to get sphere cube overlaps
    double sphere_cube_overlap(double, double, double, double,
                               double, double, double, 
                               double, double, double);

    //Destructor
    virtual ~CPLField() {}

    //Variables
    CPL::ndArray<double> array;

    double min[3], max[3], dxyz[3], dA[3], dV;

//private:

    CPL::ndArray<double> celltonode(CPL::ndArray<double> cell, int n, int ic, int jc, int kc);
    std::vector<double> interpolate(double r[], CPL::ndArray<double> cell_array, int n, int order);
};

}

#endif  // CPLField_H_
