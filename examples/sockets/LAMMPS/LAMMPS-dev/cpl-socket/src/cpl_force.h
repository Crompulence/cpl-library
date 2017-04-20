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
#include "CPL.h"


class CPLForce{

protected:

    double min[3], max[3], dxyz[3], dA[3], dV;
    CPL::ndArray<double> field;

public:

    //Constructors
    CPLForce(int nd, int icell, int jcell, int kcell);
    CPLForce(CPL::ndArray<double>);

    //Getters and setters
    void set_field(CPL::ndArray<double> field);
    void set_minmax(double min_in[], double max_in[]);
    CPL::ndArray<double> get_field();

    //Get cell values
    std::vector<int> get_cell(double r[]);    

    //Actual code
    void pre_force(double r[], double v[], double a[]);
    std::vector<double> get_force(double r[], double v[], double a[]);

    //Destructor
    ~CPLForce() {}

private:

    void set_dxyz();

};


class CPLForceFlekkoy : public CPLForce {

public:

    CPLForceFlekkoy(CPL::ndArray<double> field);
    CPLForceFlekkoy(int nd, int icell, int jcell, int kcell);

    void pre_force(double x[], double v[], double a[]);
    std::vector<double> get_force(double r[], double v[], double a[]);

    double flekkoyGWeight(double y, double ymin, double ymax);

    void resetsums();

private:

    friend class CPL_Force_Test_test_pre_force_Test;
    friend class CPL_Force_Test_test_pre_force_varydomain_Test;

    void initialisesums(CPL::ndArray<double> f);

    CPL::ndArray<double> gSums;
    CPL::ndArray<double> nSums;

};

#endif  // CPLForce_H_
