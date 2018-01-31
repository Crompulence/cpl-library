#include <iostream>
#include <vector>
#include "overlap.hpp"

extern "C" double sphere_cube_overlap (double scx, double scy, double scz, double sr,
                                       double xb, double yb, double zb, 
                                       double xt, double yt, double zt);

