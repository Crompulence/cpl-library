
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

const std::string noname = std::string("noname");

class CPLField{

public:

    //Constructors
    CPLField(int, int, int, int, const std::string& name_="default");
    CPLField(CPL::ndArray<double>, const std::string& name_="default");
    CPLField(int, CPL::ndArray<double>, const std::string& name_="default");

//    CPLField(int, int, int, int, const std::string&);
//    CPLField(CPL::ndArray<double> arrayin, const std::string&);
//    CPLField(int, CPL::ndArray<double>, const std::string&);

    //Getters and setters
    void set_array(CPL::ndArray<double> arrayin);
    void set_minmax(double min_in[], double max_in[]);
    void default_minmax();
    void set_dxyz();
    void zero_array();
    std::string field_name();

    CPL::ndArray<double> get_array();
    //std::unique_ptr<CPL::ndArray <double>> get_array_pointer();
    CPL::ndArray<double>& get_array_pointer();

    //Get cell values
    std::vector<int> get_cell(const double r[]);
    std::vector<double> get_dA();
    double get_dV();

    //Function to get sphere cube overlaps
    double sphere_cube_overlap(double, double, double, double,
                               double, double, double, 
                               double, double, double);

    // functions to add values to array
    void add_to_array(int n, int i, int j, int k, double value);
    void add_to_array(const double r[], const double value[]);
    void add_to_array(const double r[], const double value);
    void add_to_array(const double r[], double s, const double value[]);
    void add_to_array(const double r[], double s, const double value);

    // functions to get value from cell i,j,k
    std::vector<double> interpolate(double r[]);
    double get_array_value(const int index, int i, int j, int k);
    std::vector<double> get_array_value(const std::vector<int> indices, int i, int j, int k);
    std::vector<double> get_array_value(const std::vector<int> indices, const double r[]);
    double get_array_value(const int index, const double r[]);
    double get_array_value(const double r[]);
    std::vector<double> get_array_value_interp(const std::vector<int> indices, const double r[]);

    //Variables
    int nd;
    CPL::ndArray<double> array;
    int icell, jcell, kcell;
    std::string name;

    double min[3], max[3], dxyz[3], dA[3], dV;

//private:

    CPL::ndArray<double> celltonode(CPL::ndArray<double> cell, int n, int ic, int jc, int kc);
    std::vector<double> interpolate(const double r[], 
                                    CPL::ndArray<double> cell_array, 
                                    const std::vector<int> indices, int order);

    //Destructor
    virtual ~CPLField() {
 //       array.clear();
    }
};

}

#endif  // CPLField_H_
