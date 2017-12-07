
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
#include <valarray>
#include <map>

namespace CPL{


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

    //TODO: Ask Ed about this
    // std::vector<double> interpolate(double r[]);    
    // CPL::ndArray<double> celltonode(CPL::ndArray<double> cell, int n, int ic, int jc, int kc);
    // std::vector<double> interpolate(double r[], CPL::ndArray<double> cell_array, int n, int order);
    // std::vector<int> cell(std::valarray<double> r);
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


} //End Namespace
#endif  // CPLField_H_
