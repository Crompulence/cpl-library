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

    The USHER algorithm finds a suitable position for inserting a simple
    spherical particle in a dense fluid. Details may be found in the following
    publication:

        Delgado-Buscalioni & Coveney, J. Chem. Phys. 119, 2 (2003).

Types

    VectorScalarPair       A Vector3D and a double, useful for storing
                           the force and potential energy that are 
                           used often in the USHER algorithm

Classes

    UsherBase              Generic Usher abstract base class. Derived classes 
                           are required to override the pure-virtual functions 
                           that require access to data that is internal to the
                           MD simulation. 

Author(s)

    David Trevelyan

*/
#ifndef CPL_USHERBASE_H_INCLUDED
#define CPL_USHERBASE_H_INCLUDED

#include<utility>
#include "CPL_usherExceptions.h"
#include "CPL_vector3D.h"


namespace CPL
{

    // Handy data-type used to return force and potential energy together 
    typedef std::pair<Vector3D, double> VectorScalarPair;

    // Usher base class, docs to come TODO (d.trevelyan@imperial.ac.uk)
    class UsherBase
    {

    public:

        UsherBase();
        UsherBase (int maxAttempts, double tol);
        UsherBase 
        (
            int maxAttempts,
            double tol,
            double maxStep,
            double overlapU,
            int maxIterations
        );
        virtual ~UsherBase();

        int attemptToInsertParticlesInRegion
        (
            int numberOfParticles, double targetPotential,
            double regionXMin, double regionYMin, double regionZMin, 
            double regionXMax, double regionYMax, double regionZMax
        );

    protected:

        virtual void insertOneAtPosition (double x, double y, double z) = 0;

        virtual Vector3D findPositionWithPotential 
        (
            double targetPotential, 
            double initialX, 
            double initialY,
            double initialZ
        );

        virtual VectorScalarPair calculateForceAndPotentialAt
        (
            Vector3D& testPosition 
        ) const = 0;

        virtual double getOverlapStepSize (double potentialEnergy) const = 0; 

        int maxAttempts;
        int maxIterations; 
        double maxStepSize;
        double stepSize;
        double overlapPotential; 
        double targetTolerance; 

    };

} // end namespace CPL

#endif // CPL_USHERBASE_H_INCLUDED
