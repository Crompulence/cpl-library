
//=====================================================
// Class of type CPLField
//    In the interest of dependency injection, a field 
//    class with dxyz, min, max, interpolation and get 
//    cell. Build with CPL_ndArray as dependency and 
//    use as input to constructor of CPLField.
//=====================================================

#ifndef CPLFieldMap_H_
#define CPLFieldMap_H_

#include <vector>
#include <iostream>
#include <utility>
#include <memory>
#include <valarray>
#include <map>

#include "CPL_ndArray.h"
#include "CPL_field.h"

namespace CPL{

//const std::string noname = std::string("noname");

class CPLFieldMap: public CPLField{

public:

    //Constructors
    
    CPLFieldMap(int, int, int, int, const std::string& name_="default");
    CPLFieldMap(CPL::ndArray<double>, const std::string& name_="default");
    CPLFieldMap(int, CPL::ndArray<double>, const std::string& name_="default");
        //Getters and setters
    void set_array(CPL::ndArray<double> arrayin);
    void set_minmax(double min_in[], double max_in[]);
    void default_minmax();
    void set_defaults();
    void set_dxyz();
    void zero_array();
    std::string field_name();

        //Get cell values
    std::vector<int> get_cell(const double r[]);
    std::vector<double> get_dA();
    double get_dV();

    void set_Voro_data(std::vector<int> cunts_, int sum_sample_, int total_sample_);

    void add_to_array(const int n, const int i, const int j,
                      const int k, const double value);
    void add_to_array(const double r[], const double value[]);
    void add_to_array(const double r[], const double value);
    void add_to_array(const int index, const double r[], const double value);
    void add_to_array(const double r[], double s, const double value[]);
    void add_to_array(const double r[], double s, const double value);
    void add_to_array(const int index, const double r[], double s, const double value);

    void add_to_array_map_field(const int n, const int i, const int j,
                            const int k, const double value);
    void add_to_array_map_field(const double r[], const double value[]);
    void add_to_array_map_field(const double r[], const double value);
    void add_to_array_map_field(const int index, const double r[], const double value);    
    //std::vector<double> get_array_value_map_field(const std::vector<int> indices, int i, int j, int k);      
    CPL::ndArray<double> get_array();    
    CPL::ndArray<double>& get_array_pointer();
    //double get_array_value_map_field(const int index, int i, int j, int k);
    inline double get_array_value(const int index, int i, int j, int k){
        return array(index, i, j, k);
    }
    std::vector<double> get_array_value(const std::vector<int> indices, int i, int j, int k);
    std::vector<double> get_array_value(const std::vector<int> indices, const double r[]);
    double get_array_value(const int index, const double r[]);
    double get_array_value(const double r[]);
    
    double get_array_value_map_field(const int index, const double r[]);
    double get_array_value_map_field(const double r[]);
    std::vector<double> get_array_value_map_field(const std::vector<int> indices, const double r[]);
//Variables
    double get_array_value_map_field(const int index, int i, int j, int k);
    int nd;
    CPL::ndArray<double> array;
    int icell, jcell, kcell;
    std::string name;

    double min[3], max[3], dxyz[3], dA[3], dV;


    bool quickcalc;
    int quickcheck; 
    int facequickcheck;
    int edgequickcheck;
    int nonquickcheck;

    int sum_sample = 0;
    int total_sample = 0;

    std::vector<int> cuntsvector;
};

}

#endif  // CPLField_H_
