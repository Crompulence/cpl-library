
# cpl-library [![Build Status](https://travis-ci.org/Crompulence/cpl-library.png?branch=master)](https://travis-ci.org/Crompulence/cpl-library)

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


2 ) Install
==========

a) Compiling the library 
------------------------

Compiling CPL Library with GCC, as a shared library, is likely to 
work with the makefile provided in the top level directory:

    $  cd /PATH/TO/cpl-library
    $  make 

If this doesn't work, the compilers and flags may need to be specified
manually. The makefiles provided with CPL Library import flags and variables 
from the make/[platform].inc, where the variable [platform] is specified
by the file make/PLATFORM. A number of template are provided
in make, and the file make/PLATFORM contains the string "gcc" 
by default to use make/gcc.inc. 
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
