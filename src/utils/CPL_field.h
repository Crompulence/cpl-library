
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

#include "CPL_ndArray.h"

class CPLField{

public:

    //Constructors
    CPLField(int nd, int icell, int jcell, int kcell);
    CPLField(const CPL::ndArray<double>&);

    //Getters and setters
    void set_field(CPL::ndArray<double> field);
    void set_minmax(double min_in[], double max_in[]);
    CPL::ndArray<double> get_field();

    //Get cell values
    std::vector<int> get_cell(double r[]);
    std::vector<double> interpolate(double r[]);    

    //Destructor
    virtual ~CPLField() {}

private:

    double min[3], max[3], dxyz[3], dA[3], dV;
    CPL::ndArray<double> field;

    void set_dxyz();

};

#endif  // CPLField_H_
