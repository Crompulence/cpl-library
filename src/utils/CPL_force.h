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

#include "CPL_ndArray.h"
#include "CPL_field.h"

//typedef std::unique_ptr<CPL::CPLField> CPLFieldPtr;

class CPLForce{

protected:

    double min[3], max[3], dxyz[3], dA[3], dV;

    //CPL::CPLField* fieldptr;
    std::shared_ptr<CPL::CPLField> fieldptr;
    std::vector<std::shared_ptr<CPL::CPLField>> fields;

public:

    //Constructors
    CPLForce(int nd, int icell, int jcell, int kcell);
    CPLForce(CPL::ndArray<double> arrayin);

    //Getters and setters
    void set_field(CPL::ndArray<double> arrayin);
    void set_minmax(double min_in[], double max_in[]);
    CPL::ndArray<double> get_field();

    //Get cell values
    std::vector<int> get_cell(double r[]);    
    std::vector<double> get_dA();

    //Pre force collection and get force calculation
    // position, velocity, acceleration, mass, sigma, epsilon
    virtual void pre_force(double r[], double v[], double a[], 
                           double m, double s, double e) = 0;
    virtual std::vector<double> get_force(double r[], double v[], double a[], 
                                          double m, double s, double e) = 0;
    virtual void resetsums() = 0;

    bool calc_preforce = false;

    //Destructor
    virtual ~CPLForce() {
    //    delete fieldptr;
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

    void resetsums();

    bool calc_preforce = false;

};


class CPLForceVelocity : public CPLForce {

public:

    //Constructors
    CPLForceVelocity(CPL::ndArray<double> field);
    CPLForceVelocity(int nd, int icell, int jcell, int kcell);

    //Pre force collection and get force calculation
    // position, velocity, acceleration, mass, sigma, epsilon
    void pre_force(double r[], double v[], double a[], 
                   double m, double s, double e);
    std::vector<double> get_force(double r[], double v[], double a[], 
                                  double m, double s, double e);

    void resetsums();
    bool calc_preforce = true;

private:

    CPL::ndArray<double> vSums;
    CPL::ndArray<double> nSums;

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

    //Force specific things
    double flekkoyGWeight(double y, double ymin, double ymax);

    CPL::ndArray<double> gSums;
    CPL::ndArray<double> nSums;

    bool calc_preforce = true;

private:

    friend class CPL_Force_Test_test_flekkoy_pre_force_Test;
    friend class CPL_Force_Test_test_flekkoy_pre_force_varydomain_Test;
    friend class CPL_Force_Test_test_flekkoy_get_force_Test;

    void initialisesums(CPL::ndArray<double> f);
    void resetsums();

};


class CPLForceDrag : public CPLForce {

public:

    //Constructors
    CPLForceDrag(CPL::ndArray<double> field);
    CPLForceDrag(int nd, int icell, int jcell, int kcell);
    CPLForceDrag(int nd, int icell, int jcell, int kcell, bool overlap);

    //Pre force collection and get force calculation
    // position, velocity, acceleration, mass, radius, interaction
    void pre_force(double r[], double v[], double a[], 
                   double m, double s, double e);
    std::vector<double> get_force(double r[], double v[], double a[], 
                                  double m, double s, double e);
    std::vector<double> get_force_unopt(double r[], double v[], double a[], 
                                  double m, double s, double e);

    //Force specific things
    double drag_coefficient();

    bool calc_preforce = true;
    bool use_overlap = false;

//    CPL::ndArray<double> nSums;
//    CPL::ndArray<double> eSums;
//    CPL::ndArray<double> FSums;
//    CPL::ndArray<double> FcoeffSums;

    //Vector of fields
    void build_fields_list();

    //Shared instead of unique as we also keep in fields list
    std::shared_ptr<CPL::CPLField> nSums;
    std::shared_ptr<CPL::CPLField> eSums;
    std::shared_ptr<CPL::CPLField> FSums;
    std::shared_ptr<CPL::CPLField> FcoeffSums;

protected:

    void initialisesums(CPL::ndArray<double> f);
    void resetsums();
};


class CPLForceGranular : public CPLForceDrag {

public:

    //Constructors
    CPLForceGranular(CPL::ndArray<double> field);
    CPLForceGranular(int nd, int icell, int jcell, int kcell);

    //Pre force collection and get force calculation
    // position, velocity, acceleration, mass, radius, interaction
    void pre_force(double r[], double v[], double a[], 
                   double m, double s, double e);
    std::vector<double> get_force(double r[], double v[], double a[], 
                                  double m, double s, double e);

    //Force specific things
    double Reynolds_number(double D, double U, double rho, double mu, double eps);
    double porousity_exponent(double Re);
    double drag_coefficient(double Re);
    double magnitude(std::vector<double> v);

    bool calc_preforce = true;

    CPL::ndArray<double> vSums;
    CPL::ndArray<double> nSums;
    CPL::ndArray<double> eSums;
    CPL::ndArray<double> FSums;
    CPL::ndArray<double> FcoeffSums;

private:

    void initialisesums(CPL::ndArray<double> f);
    void resetsums();

};

class CPLForceBVK : public CPLForceGranular {

public:

    //Constructors
    CPLForceBVK(CPL::ndArray<double> field);
    CPLForceBVK(int nd, int icell, int jcell, int kcell);

    //Pre force collection and get force calculation
    // position, velocity, acceleration, mass, radius, interaction
    void pre_force(double r[], double v[], double a[], 
                   double m, double s, double e);
    std::vector<double> get_force(double r[], double v[], double a[], 
                                  double m, double s, double e);

    //Force specific things
    double Reynolds_number(double D, double U, double rho, double mu, double eps);
    double Stokes(double D, double U, double mu);
    double CPLForceBVK_expression(double eps, double D, double U, double rho, double mu);

    bool calc_preforce = true;

    CPL::ndArray<double> vSums;
    CPL::ndArray<double> nSums;
    CPL::ndArray<double> eSums;
    CPL::ndArray<double> FSums;
    CPL::ndArray<double> FcoeffSums;

private:

    void initialisesums(CPL::ndArray<double> f);
    void resetsums();

};

#endif  // CPLForce_H_
