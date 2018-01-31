#include <iostream>
#include <vector>
#include "sphere_cube.hpp"

extern "C" double sphere_cube_overlap (double scx, double scy, double scz, double sr,
                                       double xb, double yb, double zb, 
                                       double xt, double yt, double zt)
{
    vector_t v0{xb, yb, zb};
    vector_t v1{xt, yb, zb};
    vector_t v2{xt, yt, zb};
    vector_t v3{xb, yt, zb};
    vector_t v4{xb, yb, zt};
    vector_t v5{xt, yb, zt};
    vector_t v6{xt, yt, zt};
    vector_t v7{xb, yt, zt};

    Hexahedron cube{v0, v1, v2, v3, v4, v5, v6, v7};

    vector_t c{scx, scy, scz};
    Sphere s{c, sr};

    scalar_t result = overlap(s, cube);

	//std::cout << "result:    " << result << std::endl;

    return result;
}


int main( int argc, char *argv[] )
{

    vector_t v0{-1, -1, -1};
    vector_t v1{ 1, -1, -1};
    vector_t v2{ 1,  1, -1};
    vector_t v3{-1,  1, -1};
    vector_t v4{-1, -1,  1};
    vector_t v5{ 1, -1,  1};
    vector_t v6{ 1,  1,  1};
    vector_t v7{-1,  1,  1};

    Hexahedron hex{v0, v1, v2, v3, v4, v5, v6, v7};

    vector_t c{1,1,1};
    float r = 1;
    Sphere s{c, r};

    scalar_t result = overlap(s, hex);

	std::cout << "volume hex:    " << result << std::endl;

    double scx = 1;double scy=1; double scz=1; double sr = 1;
    double xb = -1; double yb=-1; double zb=-1;
    double xt =  1; double yt= 1; double zt= 1;
    result = sphere_cube_overlap (scx, scy, scz, sr, xb, yb, zb, xt, yt, zt);

	std::cout << "volume hex:    " << result << std::endl;

}
