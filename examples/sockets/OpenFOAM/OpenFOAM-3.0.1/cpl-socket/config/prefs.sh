#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     config/example/prefs.sh
#
# Description
#     Preset variables for the OpenFOAM configuration - POSIX shell syntax.
#
#     The prefs.sh file will be sourced by the OpenFOAM etc/bashrc when it is
#     found by foamEtcFile.
#
# See Also
#     'foamEtcFile -help' or 'foamEtcFile -list' for information about the
#     paths searched
#
#------------------------------------------------------------------------------

#- Compiler location:
#    foamCompiler= system | ThirdParty (OpenFOAM)
foamCompiler=system

#- Compiler:
#    WM_COMPILER = Gcc | Gcc45 | Gcc46 | Gcc47 | Clang | Icc (Intel icc)
export WM_COMPILER=Gcc
unset WM_COMPILER_ARCH WM_COMPILER_LIB_ARCH

#- Architecture:
#    WM_ARCH_OPTION = 32 | 64
export WM_ARCH_OPTION=64

#- Precision:
#    WM_PRECISION_OPTION = DP | SP
export WM_PRECISION_OPTION=DP

#- Optimised, debug, profiling:
#    WM_COMPILE_OPTION = Opt | Debug | Prof
export WM_COMPILE_OPTION=Opt

#- MPI implementation:
#    WM_MPLIB = SYSTEMOPENMPI | OPENMPI | MPICH | MPICH-GM | HPMPI
#               | GAMMA | MPI | QSMPI | SGIMPI | SYSTEMMPI
# Modify the MPI_ROOT to the root folder of your mpi library i.e /opt/mpich3.2
export MPI_ROOT=$(dirname $(mpicc -show | sed -e 's/.*-I\([^ ]*\).*/\1/'))
export MPI_ARCH_INC="-I$MPI_ROOT/include"
export MPI_ARCH_LIBS="-L$MPI_ROOT/lib -Wl,-rpath -Wl,$MPI_ROOT/lib -Wl,--enable-new-dtags -lmpi"
export MPI_ARCH_FLAGS=
export WM_MPLIB=SYSTEMMPI

# ----------------------------------------------------------------- end-of-file
