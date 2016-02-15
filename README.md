
# cpl-library [![Build Status](https://travis-ci.org/Crompulence/cpl-library.png?branch=master)](https://travis-ci.org/Crompulence/cpl-library)

~~~
    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________        
     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________       
      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________      
       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________     
        _\/\\\_____________\/\\\/////////____\/\\\_____________    
         _\//\\\____________\/\\\_____________\/\\\_____________   
          __\///\\\__________\/\\\_____________\/\\\_____________  
           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_ 
            _______\/////////__\///______________\///////////////__

~~~ 

<p align="center">
  www.cpl-library.org
</p>



CPL Library is a communications and topology management system for
coupling any continuum fluid dynamics (CFD) solver to any molecular dynamics
(MD) code written in either C, C++ or Fortran. A demonstration of CPL Library's
capability is included in the "demo" folder. Inside the demo, "dummy" MD and
CFD programs are coupled via CPL Library and exchange data on a virtual
topology. The user is encouraged to examine the source code of these programs,
which are an easy-to-understand example of how to couple other programs with
CPL Library.



Contents
========
    
 1. Pre-requisites for compilation
 2. Install 
 
  a) Compiling the library
  
  b) Compiling and running the demo programs
  
 3. License


1. Pre-requisites for compilation
=================================

- A C++14 compiler 
- A Fortran 2008 compiler
- An MPI library 

CPL Library was developed and tested using the GCC compiler collection and
MPICH, which are both free and open-source. 


2. Install
==========

a) Compiling the library 
------------------------

Compiling CPL Library with GCC, as a shared library, is likely to 
work with the makefile provided in cpl/:

    $  cd cpl
    $  make 

If this doesn't work, the compilers and flags may need to be specified
manually. The makefiles provided with CPL Library import the variables 
in the file make/[platform].inc, where the variable [platform] is specified
exactly by the contents of the file make/PLATFORM. A template is provided
in make/gcc.inc, and the file make/PLATFORM contains only the string "gcc" 
by default. The GCC template is likely to work for most machines, but
the user may create their own version (make/user-include.inc, for
example):

    $  cd ./make
    $  cp gcc.inc user-include.inc
    $  vi user-include.inc                     (make any necessary changes)
    $  cd ../cpl
    $  make PLATFORM=user-include




b) Building the demo applications
---------------------------------

Please see:

        http://www.cpl-library.org/user-guide/quick-start.shtml 


3. License
==========

CPL Library is released under the GNU GPL v3 license. Details are found in
the file LICENSE that is included with the release.
