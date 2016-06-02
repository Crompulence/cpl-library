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

    See CPL_usherBase.h

Author(s)
    
    David Trevelyan

*/

#include<cmath> 
#include<random> 
#include "CPL_usherBase.h"

namespace CPL
{

UsherBase::UsherBase()
    :
    maxAttempts (100),
    maxIterations (200),
    maxStepSize (0.4), // TODO get this value from density
    overlapPotential (1e5),
    targetTolerance (1e-3)
{
    stepSize = maxStepSize;
}

UsherBase::UsherBase (int attempts, double tol)
    :
    UsherBase()
{
    maxAttempts = attempts;
    targetTolerance = tol;
}

UsherBase::UsherBase 
(
    int attempts,
    double tol,
    double maxStep, 
    double overlapU,
    int iterations 
)
    :
    UsherBase(attempts, tol)
{
    maxStepSize = maxStep;
    stepSize = maxStepSize;
    overlapPotential = overlapU;
    maxIterations = iterations;
}

UsherBase::~UsherBase() {}


int UsherBase::attemptToInsertParticlesInRegion
(
    int numberOfParticlesToInsert,
    double targetPotential, 
    double regionXMin,
    double regionYMin,
    double regionZMin,
    double regionXMax,
    double regionYMax,
    double regionZMax
)
{
    // Random number generator (see
    // cppreference.com/w/cpp/numeric/random/uniform_real_distribution)
    std::random_device randomDevice;
    std::mt19937 randomNumberGenerator (randomDevice());
    std::uniform_real_distribution<> randomX (regionXMin, regionXMax);
    std::uniform_real_distribution<> randomY (regionYMin, regionYMax);
    std::uniform_real_distribution<> randomZ (regionZMin, regionZMax);

    int insertedParticlesCount = 0;

    for (int n = 0; n < numberOfParticlesToInsert; ++n)
    {
        for (int attempt = 0; attempt < maxAttempts; ++attempt)
        {

            try
            {
                auto startX = randomX (randomNumberGenerator);
                auto startY = randomY (randomNumberGenerator);
                auto startZ = randomZ (randomNumberGenerator);


                Vector3D r = findPositionWithPotential 
                (
                    targetPotential, startX, startY, startZ
                );

                if 
                (
                    r.x >= regionXMin && r.x <= regionXMax && 
                    r.y >= regionYMin && r.y <= regionYMax && 
                    r.z >= regionZMin && r.z <= regionZMax
                )
                {
                    insertOneAtPosition (r.x, r.y, r.z); 
                    ++insertedParticlesCount;
                    break;
                }

            }
            catch (UsherTryAgain& e)
            {
                // Debug: std::cout << e.what() << std::endl;
            }

        }  
    }

    if (insertedParticlesCount < numberOfParticlesToInsert)
    {
        throw UsherFailToInsertAll 
        (
            numberOfParticlesToInsert,
            insertedParticlesCount
        );
    }

    return insertedParticlesCount;

}


Vector3D UsherBase::findPositionWithPotential
(
    double targetU, 
    double initialX, 
    double initialY,
    double initialZ
)
{

    Vector3D searchPos {initialX, initialY, initialZ};
    auto fU = calculateForceAndPotentialAt (searchPos);
    auto f = std::get<Vector3D> (fU);
    auto fmag = f.magnitude();
    auto U = std::get<double> (fU);
    auto Uprev = U; 
    double upOrDown = std::copysign (1.0, (U - targetU));

    for (int iter = 0; iter < maxIterations; ++iter)
    {

        auto deltaU = std::abs (U - targetU);
        if (targetU != 0.0) deltaU /= std::abs (targetU); 

        if (deltaU <= targetTolerance)
        {
            return searchPos; 
        }
        else if (fmag == 0.0)
        {
            throw UsherTryAgainZeroForce{};
        }
        else if (upOrDown * (U - Uprev) > 0.0)
        {
            throw UsherTryAgainLocalMinimum{}; 
        }
        else
        {

            if (U > overlapPotential)
            {
                stepSize = getOverlapStepSize(U);
            }
            else
            {
                stepSize = std::min (maxStepSize, deltaU / f.magnitude());
            }

            searchPos += stepSize * upOrDown * (f / f.magnitude());

            Uprev = U;
            fU = calculateForceAndPotentialAt (searchPos);
            U = std::get<double> (fU);
            f = std::get<Vector3D> (fU);
            fmag = f.magnitude();
            upOrDown = std::copysign (1.0, (U - targetU));

        } 

    }

    throw UsherTryAgainMaxIterations{};

}

} // end namespace CPL
