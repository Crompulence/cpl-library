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

    misclib -- miscellaneous libraries
    A series of general library functions

Author(s)
    
    Edward Smith

*/

#ifndef CPLMiscLib_H_
#define CPLMiscLib_H_

/////////////////////////////////////////////////////////////////////
// A function to compare two strings in a case insenstive manner
/////////////////////////////////////////////////////////////////////
#include <locale>
#include <iostream>
#include <algorithm>

// templated version of my_equal so it could work with both char and wchar_t
template<typename charT>
struct my_equal {
    my_equal( const std::locale& loc ) : loc_(loc) {}
    bool operator()(charT ch1, charT ch2) {
        return std::toupper(ch1, loc_) == std::toupper(ch2, loc_);
    }
private:
    const std::locale& loc_;
};

// find substring (case insensitive)
template<typename T>
int string_contains( const T& str1, const T& str2, const std::locale& loc = std::locale() )
{
    typename T::const_iterator it = std::search( str1.begin(), str1.end(), 
        str2.begin(), str2.end(), my_equal<typename T::value_type>(loc) );
    if ( it != str1.end() ) return it - str1.begin();
    else return -1; // not found
}

//Character interface
template<typename T>
int string_contains( const T& str1, const char str2[], const std::locale& loc = std::locale() )
{
    std::string s(str2);
    return string_contains(str1, s, loc);
}


#include <iostream>
#include <sstream>

#include <algorithm>
#include <string> 
#include <cctype>

inline bool checktrue(const std::string s)
{
    //Check for number 1
    int i;
    std::istringstream(s) >> i;
    if (i == 1) {
        return true;
        std::cout << i << std::endl;
    }
    //Check word TRUE/true/True
    bool b;
    std::string scopy(s);
    std::transform(scopy.begin(), scopy.end(), 
                   scopy.begin(), ::tolower);
    std::istringstream(scopy) >> std::boolalpha >> b;
    if (b) return true;
   
    //Otherwise it's false
    return false;
}

#endif  // CPLMiscLib_H_
