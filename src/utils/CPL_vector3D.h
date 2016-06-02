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

    A (very) simple implementation of a three-dimensional vector class,
    with definitions of some basic linear-algebra operations.

Author(s)
    
    David Trevelyan

*/

#ifndef CPL_VECTOR3D_H_INCLUDED
#define CPL_VECTOR3D_H_INCLUDED

#include<ostream>

namespace CPL
{

    struct Vector3D
    {

        Vector3D();
        Vector3D (double x, double y, double z);

        double x, y, z;

        Vector3D& operator+= (const Vector3D& right);   
        Vector3D& operator-= (const Vector3D& right);   
        Vector3D& operator*= (double right);   
        Vector3D& operator/= (double right);   
        Vector3D  operator+  (const Vector3D& right);   
        Vector3D  operator-  (const Vector3D& right);   
        Vector3D  operator-  ();   
        double    operator*  (const Vector3D& right);   
        Vector3D  operator*  (double right);   
        Vector3D  operator/  (double right);   

        double magnitude();

    };
    Vector3D operator* (double left, const Vector3D& right);
    std::ostream& operator<< (std::ostream& os, const Vector3D& vec);

} // end namespace CPL

#endif // CPL_VECTOR3D_H_INCLUDED
