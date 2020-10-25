#!/bin/bash
echo "$1"

#sed -i '/#include.*pov/s/_[^_]*_[^_]*\.pov/_'$1'/g' CPL.pov
sed -i '/#include.*CPL_p/c\#include "CPL_p'$1'"' CPL.pov
sed -i '/#include.*CPL_v/c\#include "CPL_v'$1'"' CPL.pov
 
povray CPL.pov
mv CPL.png cpl$1.png
