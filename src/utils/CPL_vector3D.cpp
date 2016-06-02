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

    See CPL_vector3D.h

Author(s)
    
    David Trevelyan

*/
#include<string> 
#include<cmath> 

#include "CPL_vector3D.h"

namespace CPL
{

    Vector3D::Vector3D() 
        :
        Vector3D (0.0, 0.0, 0.0)
    {}
    Vector3D::Vector3D (double a, double b, double c) 
        :
        x (a), y(b), z(c)
    {}
    Vector3D& Vector3D::operator+= (const Vector3D& right)
    {
        x += right.x;
        y += right.y;
        z += right.z;
        return *this;
    }
    Vector3D& Vector3D::operator-= (const Vector3D& right)
    {
        x -= right.x;
        y -= right.y;
        z -= right.z;
        return *this;
    }
    Vector3D& Vector3D::operator*= (double right)
    {
        x *= right;
        y *= right;
        z *= right;
        return *this;
    }
    Vector3D& Vector3D::operator/= (double right)
    {
        x /= right;
        y /= right;
        z /= right;
        return *this;
    }
    Vector3D Vector3D::operator+ (const Vector3D& right)
    {
        Vector3D result = *this;
        result += right;
        return result;
    }
    Vector3D Vector3D::operator- (const Vector3D& right)
    {
        Vector3D result = *this;
        result += right;
        return result;
    }
    Vector3D Vector3D::operator- ()
    {
        Vector3D result {0.0, 0.0, 0.0};
        result -= *this;
        return result;
    }
    Vector3D Vector3D::operator* (double right)
    {
        Vector3D result = *this;
        result *= right;
        return result;
    }
    Vector3D operator* (double left, const Vector3D& right)
    {
        Vector3D result = right;
        result *= left;
        return result;
    }
    Vector3D Vector3D::operator/ (double right)
    {
        Vector3D result = *this;
        result /= right;
        return result;
    }
    double Vector3D::magnitude()
    {
        return std::sqrt (x*x + y*y + z*z); 
    }

    double Vector3D::operator* (const Vector3D& right)
    {
        return x*right.x + y*right.y + z*right.z;
    }
    std::ostream& operator<< (std::ostream& os, const Vector3D& vec)
    {
        os << "[" 
           << std::to_string(vec.x) << ", "
           << std::to_string(vec.y) << ", "
           << std::to_string(vec.z) << "]";
        return os;
    }

} // end namespace CPL
