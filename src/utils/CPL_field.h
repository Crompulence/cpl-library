
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
#include <valarray>
#include <map>

#include "CPL_ndArray.h"

namespace CPL{

const std::string noname = std::string("noname");

class Domain{
    public:
        Domain() : bounds(6){}
        Domain(const std::valarray<double>& domain_orig,
                  const std::valarray<double>& domain_length);
        Domain(const std::valarray<double>& domain_bounds);
        std::valarray<double> bounds;
        double V;
        double lx, ly, lz;
        double Axy, Axz, Ayz;
        std::valarray<double> length();
        std::valarray<double> area();
        virtual ~Domain() {}
    protected:
        void computeAV_();
};

class Field: public Domain {

public:
    //Constructors
    Field(): Domain(), nCells(3) {}
    Field(const CPL::Domain& domain,
          const std::vector<int>& field_cells);
    double dx, dy, dz;
    double dAxy, dAxz, dAyz;
    double dV;
    std::vector<int> nCells;
    std::valarray<double> cellLength();
    std::valarray<double> cellArea();
    std::valarray<double> getCoord(std::vector<int>& cell);
    std::vector<int> getCell(std::valarray<double>& coord);
    virtual ~Field() {};

protected:
    void computedAdV_();
};

class PortionField: public Field {
    public:
        PortionField() : Field(), cellBounds(6) {}
        PortionField(const CPL::Domain& domain,
                     const std::vector<int>& field_cells,
                     const std::vector<int>& cell_bounds):
                     Field(domain, field_cells), cellBounds(cell_bounds){}
        std::vector<int> cellBounds;
        std::vector<int> getLocalCell(std::vector<int>& glob_cell);
        std::vector<int> getLocalCell(std::vector<int>& glob_cell, bool& valid_cell);

};

class DataField: public Field {
    public:
        DataField() {}
        DataField(const std::valarray<double>& domain_bounds,
                  const std::vector<int>& field_cells, int data_length);
        DataField(const std::valarray<double>& domain_bounds,
                  CPL::ndArray<double>& field_data, bool copy=true);
        bool own_memmory;
        CPL::ndArray<double>* data;
};


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
    void set_defaults();
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
    void add_to_array(const int n, const int i, const int j,
                      const int k, const double value);
//    inline void add_to_array(const int n, const int i, 
//                             const int j, const int k, 
//                             const double value){
//        std::cout << "n=" << n << " value=" << value << std::endl;
//        std::cout << "i=" << i << " value=" << value << std::endl;
//        std::cout << "j=" << j << " value=" << value << std::endl;
//        std::cout << "k=" << k << " value=" << value << std::endl;

//        array(n, i, j, k) += value; 
//    }
    void add_to_array(const double r[], const double value[]);
    void add_to_array(const double r[], const double value);
    void add_to_array(const int index, const double r[], const double value);
    void add_to_array(const double r[], double s, const double value[]);
    void add_to_array(const double r[], double s, const double value);
    void add_to_array(const int index, const double r[], double s, const double value);

    // functions to get value from cell i,j,k
    std::vector<double> interpolate(double r[]);
    inline double get_array_value(const int index, int i, int j, int k){
        return array(index, i, j, k);
    }
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


    bool quickcalc;
    int quickcheck; 
    int facequickcheck;
    int edgequickcheck;
    int nonquickcheck;

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
