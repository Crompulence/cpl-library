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

    A collection of exception classes that may be thrown during the
    USHER algorithm. See CPL_usherBase.cpp for examples. 

    Indentation in this source file shows the class inheritance 
    structure. 

Author(s)

    David Trevelyan

*/
#ifndef CPL_USHEREXCEPTIONS_H_INCLUDED
#define CPL_USHEREXCEPTIONS_H_INCLUDED
#include<exception>

/*
            Indentation shows the inheritance structure
*/

class UsherException: public std::exception {};



    class UsherTryAgain: public UsherException {};
    // For when the USHER algorithm fails in a normal way, to signal try again

        class UsherTryAgainMaxIterations : public UsherTryAgain
        {
        public:
            virtual const char* what() const noexcept
            {
                return "USHER attempt reached maximum step iteration.";
            }
        };

        class UsherTryAgainZeroForce : public UsherTryAgain 
        {
        public:
            virtual const char* what() const noexcept
            {
                return "USHER attempt found a location with zero force.";
            }
        };

        class UsherTryAgainLocalMinimum : public UsherTryAgain 
        {
        public:
            virtual const char* what() const noexcept
            {
                return "USHER attempt overstepped a local minimum/maximum";
            }
        };



    class UsherFailNonFatal : public UsherException {};
    // For when the Usher class can't do what it was asked, but isn't critical

        class UsherFailToInsertAll : public UsherFailNonFatal 
        {

        public:

            UsherFailToInsertAll (int targetInserted, int actuallyInserted)
                :
                UsherFailNonFatal{},
                targetInsertedParticles {targetInserted},
                actuallyInsertedParticles {actuallyInserted}
            {};

            virtual const char* what() const noexcept
            {
                std::string msg ("Usher could only insert ");
                msg += std::to_string (actuallyInsertedParticles) + " of ";
                msg += std::to_string (targetInsertedParticles) + " particles";
                return msg.c_str(); 
            }

        private:

            int targetInsertedParticles;
            int actuallyInsertedParticles;

        };



    // Sucess exception that may be useful for debugging
    class UsherSuccess : public UsherException 
    {
    public:
        virtual const char* what() const noexcept
        {
            return "USHER SUCCEEDED";
        }
    };


#endif // CPL_USHEREXCEPTIONS_H_INCLUDED
