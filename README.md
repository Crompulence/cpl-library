
# cpl-library [![Build Status](https://travis-ci.org/Crompulence/cpl-library.png?branch=master)](https://travis-ci.org/Crompulence/cpl-library) ![Build Status](https://img.shields.io/docker/cloud/build/cpllibrary/cpl-library)

<p align="center">
  www.cpl-library.org
</p>


CPL Library is a communications and topology management system for
coupling any continuum fluid dynamics (CFD) solver to any molecular dynamics
(MD) code written in either Python, C, C++ or Fortran. For more details about
the philosophy, range of validity and aims of this software, please
see the [wikki page](https://github.com/Crompulence/cpl-library/wiki/CPL-Library)
and our dedicated website [www.cpl-library.org](http://www.cpl-library.org)
A demonstration of CPL Library's capability is included in the "demo" folder. 
Inside the demo, "dummy" MD and CFD programs are coupled via CPL Library and 
exchange data on a virtual topology. The user is encouraged to examine the 
source code of these programs, which are an easy-to-understand example of 
how to couple other programs with CPL Library.



Contents
========
    
 1) Pre-requisites for compilation
 2) Install <br />
  a) Compiling the library <br />
  b) Compiling and running the demo programs
 3) License
 

1 ) Pre-requisites for compilation
=================================

- A C++11 compiler 
- A Fortran 2008 compiler
- An MPI library 

CPL Library was developed and tested using the GCC compiler collection and
MPICH, which are both free and open-source. 
OpenMPI is not currently supported.
For modern versions of linux this should not be a problem and the
system GCC will be new enough to support C++11.
The repository version of MPICH is usually fine to use, although be
careful that existing MPI installations (usually OpenMPI by default) 
are not used instead or mixed in.
Alternativly, to ensure compatibility you can build MPICH yourself.
Be careful to ensure use the same version of GCC and MPI for all coupled codes.


2 ) Install
==========

a) Compiling the library 
------------------------

Compiling CPL Library with GCC, as a shared library, is likely to 
work with the makefile provided in the top level directory:

    $  cd /PATH/TO/cpl-library
    $  make PLATFORM=gcc

Which uses data in make/gcc.inc to build the library.
If this doesn't work, the compilers and flags may need to be specified
manually. The makefiles provided with CPL Library import flags and variables 
are in the make/[platform].inc, where the variable [platform] is specified
by the file make PLATFORM=somename. A number of template are provided
in make (intel, gcc), and the file make/PLATFORM contains the string of the 
previous build to use by default. 
The GCC template is likely to work for most machines, but
the user may create their own version. For example, to create
make/user-include.inc:

    $  cd ./make
    $  cp gcc.inc user-include.inc
    $  vi user-include.inc                  (make any necessary changes)
    $  cd ../
    $  make PLATFORM=user-include


b) Building the demo applications
---------------------------------

Please see:

[www.cpl-library.org/user-guide/quick-start.shtml ](http://www.cpl-library.org/user-guide/quick-start.shtml)


3 ) License
==========

CPL Library is released under the GNU GPL v3 license. Details are found in
the file LICENSE that is included with the release.



4 ) Directory Structure
=========================

The directory structure of CPL library is broadly as follows:

 - src - The main source code split into:
   - core - the core Fortran of the CPL library, in just two files CPL_module and CPL_methods
   - bindings - C, C++ and Python bindings 
   - utils - The array, field and force classes which provide functionality common to coupling problems, as well as numerical libraries, etc
 - examples - A range of example applications, details here http://www.cpl-library.org/user-guide/quick-start.shtml
 - test - tests for CPL library, including:
   - pytest - Basic init, setup and send/recv topology tests
   - gtest - tests used to validate most of the utils work as expected
   - valgrind - memory leak checking of the key CPL library functionality
   - examples, lammps, drag - various other tests of coupled cases, largly moved to APP directories
 - website - the source code for the www.cpl-library.org website
 - utils - General Scripts including the Dockerfile, scaling tests and a Python GUI to design coupled runs.
 - make - a series of .inc file specifying platform specific builds (e.g. for the UK supercomputer ARCHER)
 - bin - binary folder containing excutable files cplexec, cplf90 and cplc++
 - 3rd-party - third party libraries which can be optionally included (e.g. JSON-Fortran)

On the top level is a 
 - Makefile - Used to build CPL library
 - SOURCEME.sh - Used to add CPL library to your path (avoids needing to install to system directory)

A range of folders are created by the build process

 - lib - file to store the compiled library libcpl.so
 - include - header files

In order to run a coupled OpenFOAM and LAMMPS case, the applications (APP) repositories are required along with each of these codes. The APP repositories are stored at:

 - https://github.com/Crompulence/CPL_APP_OPENFOAM-3.0.1
 - https://github.com/Crompulence/CPL_APP_LAMMPS-DEV



