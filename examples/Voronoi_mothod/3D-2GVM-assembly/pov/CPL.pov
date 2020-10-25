#version 3.6;
// Right-handed coordinate system where the z-axis points upwards
camera {
	location <0.12,0.12,0.05>
	sky z
	right -0.15*x*image_width/image_height
	up 0.15*z
	look_at <0,0,0.0030>
}

// White background
background{rgb 1}

// Two lights with slightly different colors
light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}
light_source{<25,-12,12> color rgb <0.43,0.45,0.45>}

// The radius of the cylinders to be used when drawing the Voronoi cells
#declare r=0.000024;

// The radius of the particles
#declare s=0.008;

// The particle packing
union{
#include "CPL_p_50000.pov"
scale <1,1,1>
	pigment{rgb <0.92,0.65,1>} finish{reflection 0.17 specular 0.3 ambient 0.42}
}
// The computed Voronoi cells
union{
#include "CPL_v_50000.pov"
scale <1,1,1>
	pigment{rgb <0.7,0.95,1>} finish{specular 0.5 ambient 0.42}
}
