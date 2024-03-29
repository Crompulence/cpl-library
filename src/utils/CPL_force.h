/*

    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________
     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________
      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________
       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________
        _\/\\\_____________\/\\\/////////____\/\\\_____________
         _\//\\\____________\/\\\_____________\/\\\_____________
          __\///\\\__________\/\\\_____________\/\\\_____________
           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_
            _______\/////////__\///______________\///////////////__


                         C P L  -  L I B R A R Y

           Copyright (C) 2012-2015 Edward Smith & David Trevelyan

License

    This file is part of CPL-Library.

    CPL-Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CPL-Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CPL-Library.  If not, see <http://www.gnu.org/licenses/>.

Description

    "CPLForce" class for interfacing with CPL-Library.

Author(s)

    Edward Smith

*/




//=====================================================
//This would be an abstract base class of type CPLForce
//=====================================================

#ifndef CPLForce_H_
#define CPLForce_H_

#include <vector>
#include <iostream>
#include <memory>
#include <map>

#include "CPL_ndArray.h"
#include "CPL_field.h"

typedef std::map <std::string, std::string> map_strstr;

//typedef std::unique_ptr<CPL::CPLField> CPLuFieldPtr;
//typedef std::shared_ptr<CPL::CPLField> CPLsFieldPtr;

class CPLForce{

protected:

    int shapeVector[4];
    double min[3], max[3], dxyz[3], dA[3], dV;

    //CPL::CPLField* cfd_array_field;
    CPL::ndArray<double> array;
    std::shared_ptr<CPL::CPLField> cfd_array_field;
    std::vector<std::shared_ptr<CPL::CPLField>> fields;

public:

    //Counters for number of records in preforce, force and postfoce
    int Npre_force=0, Nforce=0, Npost_force=0;

    //Constructors
    CPLForce(int, int, int, int);
    CPLForce(CPL::ndArray<double>);

    //Getters and setters
    void set_field(CPL::ndArray<double>);
    void set_minmax(double min_in[], double max_in[]);
    CPL::ndArray<double> get_field();
    std::shared_ptr<CPL::CPLField> get_internal_fields(const std::string&);

    //Get cell values
    std::vector<int> get_cell(double r[]);    
    std::vector<double> get_dA();
    double get_dV();

    //Pre force collection and get force calculation
    // position, velocity, acceleration, mass, sigma, epsilon
    virtual void pre_force(double r[], double v[], double a[], 
                           double m, double s, double e) = 0;
    virtual std::vector<double> get_force(double r[], double v[], double a[], 
                                          double m, double s, double e) = 0;
    virtual void post_force(double r[], double v[], double a[], 
                            double m, double s, double e) = 0;
    virtual void resetsums() = 0;
    virtual void reset_instant() = 0;

    bool calc_preforce;
    bool calc_preforce_everytime;

    bool calc_postforce;
    bool calc_postforce_everytime;

//    void setup_fast_array();

    //Destructor
    virtual ~CPLForce() {
    //    delete cfd_array_field;
    }

private:

    void set_dxyz();

};

class CPLForceTest : public CPLForce {

public:

    //Constructors
    CPLForceTest(CPL::ndArray<double> field);
    CPLForceTest(int nd, int icell, int jcell, int kcell);

    //Pre force collection and get force calculation
    // position, velocity, acceleration, mass, sigma, epsilon
    void pre_force(double r[], double v[], double a[], 
                   double m, double s, double e);
    std::vector<double> get_force(double r[], double v[], double a[], 
                                  double m, double s, double e);
    void post_force(double r[], double v[], double a[], 
                   double m, double s, double e);

    void initialisesums(CPL::ndArray<double>);
    void build_fields_list();
    void resetsums();
    void reset_instant();

};


class CPLForceVelocity : public CPLForce {

public:

    //Constructors
    CPLForceVelocity(CPL::ndArray<double> field);
    CPLForceVelocity(int nd, int icell, int jcell, int kcell);

    CPLForceVelocity(CPL::ndArray<double>, map_strstr);
    CPLForceVelocity(int, int, int, int, map_strstr);

    //Pre force collection and get force calculation
    // position, velocity, acceleration, mass, sigma, epsilon
    void pre_force(double r[], double v[], double a[], 
                   double m, double s, double e);
    std::vector<double> get_force(double r[], double v[], double a[], 
                                  double m, double s, double e);
    void post_force(double r[], double v[], double a[], 
                   double m, double s, double e);

    void resetsums();
    void reset_instant();

    //Coefficient for strength of application
    double xi = 1.0;

private:

    std::shared_ptr<CPL::CPLField> nSums;
    std::shared_ptr<CPL::CPLField> nSums_mdt;
    std::shared_ptr<CPL::CPLField> vSums;
    std::shared_ptr<CPL::CPLField> vSums_mdt; //Previous timestep

    void unpack_arg_map(map_strstr arg_map);

    friend class CPL_Force_Test_test_velocity_pre_force_Test;

    void initialisesums(CPL::ndArray<double> f);

};


class CPLForceFlekkoy : public CPLForce {

public:

    //Constructors
    CPLForceFlekkoy(CPL::ndArray<double> field);
    CPLForceFlekkoy(int nd, int icell, int jcell, int kcell);

    //Pre force collection and get force calculation
    // position, velocity, acceleration, mass, sigma, epsilon
    void pre_force(double r[], double v[], double a[], 
                   double m, double s, double e);
    std::vector<double> get_force(double r[], double v[], double a[], 
                                  double m, double s, double e);
    void post_force(double r[], double v[], double a[], 
                   double m, double s, double e);

    //Force specific things
    double flekkoyGWeight(double y, double ymin, double ymax);

    CPL::ndArray<double> gSums;
    CPL::ndArray<double> nSums;

    void resetsums();
    void reset_instant();

private:

    friend class CPL_Force_Test_test_flekkoy_pre_force_Test;
    friend class CPL_Force_Test_test_flekkoy_pre_force_varydomain_Test;
    friend class CPL_Force_Test_test_flekkoy_get_force_Test;

    void initialisesums(CPL::ndArray<double> f);


};

// A base class for all drag related forces, includes overlap and interpolate
// flags, sums for porosity and valeus for density/drag coefficients and viscosity
class CPLForceDrag : public CPLForce {

public:

    //Constructors
    CPLForceDrag(CPL::ndArray<double>);
    CPLForceDrag(CPL::ndArray<double>, map_strstr);
    CPLForceDrag(int, int, int, int);
    CPLForceDrag(int, int, int, int, map_strstr);

    //Pre force collection and get force calculation
    // position, velocity, acceleration, mass, radius, interaction
    void pre_force(double r[], double v[], double a[], 
                   double m, double s, double e);
    std::vector<double> get_force(double r[], double v[], double a[], 
                                  double m, double s, double e);
    void post_force(double r[], double v[], double a[], 
                    double m, double s, double e);

    void unpack_CFD_array(CPL::ndArray<double> arrayin);

    //Force specific things
    void set_defaults();
    bool use_overlap;
    bool use_interpolate;
    bool use_gradP;
    bool use_divStress;

    virtual double drag_coefficient(double D, std::vector<double> Ui_v, double eps);
    double get_eps(double r[]);
    double Cd = 0.0000001;
    double mu = 0.001;
    double rho = 1e3;

    //Vector of fields
    void build_fields_list();

    //Shared pointer instead of unique as we also keep in fields list
    std::shared_ptr<CPL::CPLField> nSums;
    std::shared_ptr<CPL::CPLField> vSums;
    std::shared_ptr<CPL::CPLField> volSums;
    std::shared_ptr<CPL::CPLField> instant_volSums;
    std::shared_ptr<CPL::CPLField> FSums;
    std::shared_ptr<CPL::CPLField> FcoeffSums;

    void resetsums();
    void reset_instant();

protected:

    void unpack_default_arg_map(map_strstr arg_map, bool extra_args);
    virtual bool unpack_extra_arg_map(map_strstr arg_map);
    void initialisesums(CPL::ndArray<double> f);
    virtual void initialise_extrasums(CPL::ndArray<double> arrayin);

private:

    void unpack_arg_map(map_strstr arg_map);
};


class CPLForceGranular : public CPLForceDrag {

public:

    //Constructors
    using CPLForceDrag::CPLForceDrag;

//    void pre_force(double r[], double v[], double a[], 
//                   double m,   double s,   double e);

    //Granular functions
    double Stokes(double D, double mu, double eps);
    double Reynolds_number(double D, double U, double rho, double mu, double eps);
    double magnitude(std::vector<double> v);

protected:

    void initialisesums(CPL::ndArray<double> f);

};

class CPLForceStokes : public CPLForceGranular {

public:

    //Constructors
    using CPLForceGranular::CPLForceGranular;

    //Stokes specific functions
    double drag_coefficient(double D, std::vector<double> Ui_v, double eps) override;

};

class CPLForceDi_Felice : public CPLForceGranular {

public:

    //Constructors
    using CPLForceGranular::CPLForceGranular;

    //Di_Felice specific functions
    double drag_coefficient(double D, std::vector<double> Ui_v, double eps) override;

//private:

};

class CPLForceTang : public CPLForceGranular {

public:

    //Constructors
    using CPLForceGranular::CPLForceGranular;

    //BVK specific functions
    double drag_coefficient(double D, std::vector<double> Ui_v, double eps) override;

//private:

};


class CPLForceErgun : public CPLForceGranular {

public:

    //Constructors
    using CPLForceGranular::CPLForceGranular;

    //Ergun specific functions
    double drag_coefficient(double D, std::vector<double> Ui_v, double eps) override;

//private:

};


class CPLForceBVK : public CPLForceGranular {

public:

    //Constructors
    using CPLForceGranular::CPLForceGranular;

    //BVK specific functions
    double drag_coefficient(double D, std::vector<double> Ui_v, double eps) override;

//private:

};



//class CPLForceBVK_poly : public CPLForceBVK {

//public:

//    //Constructors
//    using CPLForceBVK::CPLForceBVK;
//    void initialise_extrasums(CPL::ndArray<double> arrayin) override;

//    //Shared pointer instead of unique as we also keep in fields list
//    std::shared_ptr<CPL::CPLField> DSums;
//    std::shared_ptr<CPL::CPLField> FcoeffSums_prev;

//    //BVK specific functions
//    double drag_coefficient(double r[], double D, std::vector<double> Ui_v) override;

////private:

//};

class CPLForceTenneti : public CPLForceGranular {

public:

    //Constructors
    using CPLForceGranular::CPLForceGranular;

    //Tenneti specific functions
    double drag_coefficient(double D, std::vector<double> Ui_v, double eps) override;

//private:

};

#endif  // CPLForce_H_
