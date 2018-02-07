#include <vector>
#include <math.h> 
#include <stdexcept>
#include <assert.h>

#include "CPL_ndArray.h"
#include "CPL_force.h"
#include "CPL_field.h"
#include "overlap/sphere_cube.hpp"


///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForce base class                        //
//                                                               //
///////////////////////////////////////////////////////////////////


// -> In the interest of dependency injection, I made a field 
//    class with dxyz, min, max, interpolation and get cell. 
//    Built with CPL_ndArray as dependency and can be used as input 
//    to constructor of CPLForce here. SEE CPL_field.cpp
// 


//Constructors
//CPLForce::CPLForce(int nd, int icells, int jcells, int kcells) : 
//    std::unique_ptr<CPL::CPLField> fieldptr(new CPL::CPLField(nd, icells, jcells, kcells))
//{ 
//}

//CPLForce::CPLForce(CPL::ndArray<double> arrayin) : 
//    std::unique_ptr<CPL::CPLField> fieldptr(new CPL::CPLField(arrayin))
//{
//}


//Constructors
CPLForce::CPLForce(int nd, int icells, int jcells, int kcells){ 
    //fieldptr = std::make_shared<CPL::CPLField>(nd, icells, jcells, kcells);
    fieldptr = new CPL::CPLField(nd, icells, jcells, kcells);
}



CPLForce::CPLForce(CPL::ndArray<double> arrayin) {
    //fieldptr = std::make_shared<CPL::CPLField>(arrayin);
    fieldptr = new CPL::CPLField(arrayin);
}


//Set minimum and maximum values of field application
void CPLForce::set_minmax(double min_in[], double max_in[]){
    fieldptr->set_minmax(min_in, max_in);
};

//If either min/max change or field object, we need to recalculate dx, dy and dz
void CPLForce::set_dxyz(){
    fieldptr->set_dxyz();
}

//Get dA
std::vector<double> CPLForce::get_dA(){
    return fieldptr->get_dA();
}


//Get cell from min/max and dx
std::vector<int> CPLForce::get_cell(double r[]){
    std::vector<int> cell(3);
    cell = fieldptr->get_cell(r);
    return cell;
}



//Pre force collection of sums (should this come from LAMMPS fix chunk/atom bin/3d)
void CPLForce::pre_force(double r[], double v[], double a[], double m, double s, double e) {
//    throw std::runtime_error("CPLForce::pre_force is not defined");
}


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForce::get_force(double r[], double v[], double a[], double m, double s, double e){
//    throw std::runtime_error("CPLForce::get_force is not defined");
}


void CPLForce::set_field(CPL::ndArray<double> arrayin){
    fieldptr->set_array(arrayin);   
};

CPL::ndArray<double> CPLForce::get_field(){
    return fieldptr->get_array();
};

///////////////////////////////////////////////////////////////////
//                                                               //
//                      CPLForceTest                             //
//                                                               //
///////////////////////////////////////////////////////////////////

//Constructor using cells
CPLForceTest::CPLForceTest(int nd, int icells, int jcells, int kcells) : CPLForce(nd, icells, jcells, kcells){
}

//Constructor of datatype
CPLForceTest::CPLForceTest(CPL::ndArray<double> arrayin) : CPLForce(arrayin){
}


void CPLForceTest::resetsums(){
}

//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
void CPLForceTest::pre_force(double r[], double v[], double a[], double m, double s, double e){
//    throw std::runtime_error("CPLForceTest::pre_force is not defined");
}


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForceTest::get_force(double r[], double v[], double a[], double m, double s, double e){

    std::vector<double> f(3); 
    std::vector<int> cell = get_cell(r);
    CPL::ndArray<double> array = fieldptr->get_array();
    f[0] = array(0, cell[0], cell[1], cell[2]);
    f[1] = array(1, cell[0], cell[1], cell[2]);
    f[2] = array(2, cell[0], cell[1], cell[2]);

//    std::cout << "CPLForceTest " <<  
//              cell[0] << " " << cell[1] << " " <<  cell[2] << " " 
//              << f[0] << " " << f[1] << " " << f[2] << std::endl;

    return f;
}


///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForceVelocity                           //
//                                                               //
///////////////////////////////////////////////////////////////////

//Constructor using cells
CPLForceVelocity::CPLForceVelocity(int nd, int icells, int jcells, int kcells) : CPLForce(nd, icells, jcells, kcells){
    initialisesums(fieldptr->get_array());
}

//Constructor of datatype
CPLForceVelocity::CPLForceVelocity(CPL::ndArray<double> arrayin) : CPLForce(arrayin){
    initialisesums(arrayin);
}

void CPLForceVelocity::initialisesums(CPL::ndArray<double> arrayin){
    int vsumsShape[4] = {arrayin.shape(0), arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    vSums.resize(4, vsumsShape); // Sum of velocity
    int nsumsShape[3] = {arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    nSums.resize(3, nsumsShape); // Sum of number of particles  
    resetsums();
}

void CPLForceVelocity::resetsums(){
    nSums = 0.0;  vSums = 0.0;
}

//Pre force collection of sums (should this come from LAMMPS fix chunk/atom bin/3d)
void CPLForceVelocity::pre_force(double r[], double v[], double a[], double m, double s, double e) {

    // Find in which cell number (local to processor) is the particle
    // and sum all the velocities for each cell.
    std::vector<int> cell = get_cell(r);

    nSums(cell[0], cell[1], cell[2]) += 1.0; 
    for (int i = 0; i<3; ++i)
        vSums(i, cell[0], cell[1], cell[2]) += v[i];

}


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForceVelocity::get_force(double r[], double v[], double a[], double m, double s, double e){

    std::vector<double> f(3); 
    std::vector<int> cell = get_cell(r);
    CPL::ndArray<double> array = fieldptr->get_array();
    int N = nSums(cell[0], cell[1], cell[2]);
    if (N < 1.0) {
        std::cout << "Warning: 0 particles in cell (" 
                  << cell[0] << ", " << cell[1] << ", " << cell[2] << ")"
                  << std::endl;
        f[0]=0.0; f[1]=0.0; f[2]=0.0;
        return f;
    } else {
        for (int i=0; i<3; i++){
            f[i] = ( array(i, cell[0], cell[1], cell[2])
                    -vSums(i, cell[0], cell[1], cell[2])/N);
        }
        return f;
    }
}



///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForceFlekkoy                            //
//                                                               //
///////////////////////////////////////////////////////////////////
//The Flekkoy constraint IS A type of CPL force

//Constructor using cells
CPLForceFlekkoy::CPLForceFlekkoy(int nd, int icells, int jcells, int kcells) : CPLForce(nd, icells, jcells, kcells){
    initialisesums(fieldptr->get_array());
}

//Constructor of datatype
CPLForceFlekkoy::CPLForceFlekkoy(CPL::ndArray<double> arrayin) : CPLForce(arrayin){
    initialisesums(arrayin);
}

void CPLForceFlekkoy::initialisesums(CPL::ndArray<double> arrayin){

    int sumsShape[3] = {arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    gSums.resize(3, sumsShape); // Sum of Flekkøy g weights
    nSums.resize(3, sumsShape); // Sum of number of particles  
    resetsums();
}

void CPLForceFlekkoy::resetsums(){
    nSums = 0.0;  gSums = 0.0;
}

// See Flekkøy, Wagner & Feder, 2000 Europhys. Lett. 52 271, footnote p274
double CPLForceFlekkoy::flekkoyGWeight(double y, double ymin, double ymax) {

    if (y > ymax) {
		//std::cout << "y: " << y << " ymax: " << ymax << std::endl;
        throw std::domain_error("flekkoyGWeight Error: Position argument y greater than ymin");
	}

    // K factor regulates how to distribute the total force across the volume.
    // 1/K represent the fraction of the constrain region volume used.
    // Flekøy uses K = 2.
    double K = 1.0;

    // Define re-scaled coordinate 
    double L = ymax - ymin;
    double yhat = y - ymin - (1.0 - 1.0/K)*L;
    double g;

    if (yhat > 0.0)
        g = 2.0*(1.0/(L - K*yhat) - 1.0/L - K*yhat/(L*L));
    else
        g = 0.0;
    
    return g;

}

//Pre force collection of sums (should this come from LAMMPS fix chunk/atom bin/3d)
void CPLForceFlekkoy::pre_force(double r[], double v[], double a[], double m, double s, double e) {

    // Find in which cell number (local to processor) is the particle
    // and sum all the Flekkøy weights for each cell.
    std::vector<int> cell = get_cell(r);

    double g = flekkoyGWeight(r[1], fieldptr->min[1], fieldptr->max[1]);
    nSums(cell[0], cell[1], cell[2]) += 1.0; 
    gSums(cell[0], cell[1], cell[2]) += g;

}

//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForceFlekkoy::get_force(double r[], double v[], double a[], double m, double s, double e){

    std::vector<double> dA = fieldptr->get_dA(); 
    std::vector<double> f(3); 
    std::vector<int> cell = get_cell(r);
    CPL::ndArray<double> array = fieldptr->get_array();
    double n = nSums(cell[0], cell[1], cell[2]);

    if (n < 1.0) {
        std::cout << "Warning: 0 particles in cell (" 
                  << cell[0] << ", " << cell[1] << ", " << cell[2] << ")"
                  << std::endl;
        f[0]=0.0; f[1]=0.0; f[2]=0.0;
        return f;
    } else {

        // Since the Flekkoy weight is considered only for 0 < y < L/2, for cells 
        // that are completely in y < 0 gSums(i, j, k) will be always 0.0 so can 
        // produce a NAN in the g/gSums division below.
        double g = flekkoyGWeight(r[1], fieldptr->min[1], fieldptr->max[1]);

        if (gSums(cell[0], cell[1], cell[2]) > 0.0) {
            double gdA = (g/gSums(cell[0], cell[1], cell[2])) * dA[1];

            // Normal to the X-Z plane is (0, 1, 0) so (tauxy, syy, tauxy)
            // are the only components of the stress tensor that matter.
            f[0] = gdA * array(1, cell[0], cell[1], cell[2]);
            f[1] = gdA * array(4, cell[0], cell[1], cell[2]);
            f[2] = gdA * array(7, cell[0], cell[1], cell[2]);
            return  f;
        }
    }

}


///////////////////////////////////////////////////////////////////
//                                                               //
//                      CPLForceDrag                             //
//                                                               //
///////////////////////////////////////////////////////////////////


//Constructor using cells
CPLForceDrag::CPLForceDrag(int nd, int icells, int jcells, int kcells) : CPLForce(nd, icells, jcells, kcells){
    initialisesums(fieldptr->get_array());
}

//Constructor of datatype
CPLForceDrag::CPLForceDrag(CPL::ndArray<double> arrayin) : CPLForce(arrayin){
    initialisesums(arrayin);
}

//Constructor with optional argument overlap
CPLForceDrag::CPLForceDrag(int nd, int icells, int jcells, int kcells, bool overlap) : CPLForceDrag(nd, icells, jcells, kcells){
    use_overlap = overlap;
    initialisesums(fieldptr->get_array());
}


void CPLForceDrag::initialisesums(CPL::ndArray<double> arrayin){

    int nsumsShape[3] = {arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    nSums.resize(3, nsumsShape); // Sum of number of particles

    int esumsShape[3] = {arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    eSums.resize(3, esumsShape); // Sum of porosity of particles  
    //use_overlap = true;

    int FsumsShape[4] = {3, arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    FSums.resize(4, FsumsShape); // Sum of force on particles  

    int FcoeffSumsShape[3] = {arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    FcoeffSums.resize(3, FcoeffSumsShape); // Sum of porosity of particles  

    resetsums();
}

void CPLForceDrag::resetsums(){
    nSums = 0.0; eSums = 0.0; FSums=0.0; FcoeffSums=0.0;
}

//Arbitary constant
double CPLForceDrag::drag_coefficient() {
    return 0.0000001;
}

//Pre force collection of sums (should this come from LAMMPS fix chunk/atom bin/3d)
//void CPLForceDrag::pre_force(double r[], double v[], double a[], double m, double s, double e) {

//    // Should use field.add_volume(r, radius);
//    double radius = s;
//    double volume = (4./3.)*M_PI*pow(radius,3); 
//    std::vector<int> cell = get_cell(r);
//    nSums(cell[0], cell[1], cell[2]) += 1; 


//}


//Pre force collection of sums including overlap code to assign volumes
void CPLForceDrag::pre_force(double r[], double v[], double a[], double m, double s, double e) {

    int ip, jp, kp;
    double radius = s;
    double volume = (4./3.)*M_PI*pow(radius,3); 
    std::vector<int> cell = get_cell(r);
    double box[6];
    nSums(cell[0], cell[1], cell[2]) += 1;

    if (! use_overlap)
    {
        eSums(cell[0], cell[1], cell[2]) += volume; 
    } else {
        //Get fraction of sphere in a volume
        double dx = fieldptr->dxyz[0];
        double dy = fieldptr->dxyz[1];
        double dz = fieldptr->dxyz[2];
        int nxps = ceil(radius/dx);
        int nyps = ceil(radius/dy);
        int nzps = ceil(radius/dz);
        int i = cell[0]; int j = cell[1]; int k = cell[2];

//        std::cout << "overlap calc "  
//                  << dx << " " << dy << " " << dz << " " 
//                  << nxps << " " << nyps << " " << nzps << std::endl;

        for (int ic=-nxps; ic<nxps+1; ic++) {
        for (int jc=-nyps; jc<nyps+1; jc++) {
        for (int kc=-nzps; kc<nzps+1; kc++) {
            ip = i+ic; jp = j+jc; kp = k+kc;
            box[0] = (ip  )*dx;
            box[1] = (jp  )*dy;
            box[2] = (kp  )*dz;
            box[3] = (ip+1)*dx;
            box[4] = (jp+1)*dy;
            box[5] = (kp+1)*dz;

            //Input sphere centre, radius and 6 corners of cell
            double Vsphereinbox = fieldptr->sphere_cube_overlap(r[0], r[1], r[2], radius, 
                                                                box[0], box[1], box[2], 
                                                                box[3], box[4], box[5]);

            // Here there should be extra halo padding in esums to the 
            // size of Nxps/Nyps/Nzps which store fractions which would
            // need to be sent to adjacent processes. This would then
            // happend using the CPL_swaphalo method after pre_force has
            // been called for every particle on every process. However, 
            // this is complex so instead we simply dump overlap at the
            // edge of the current process. This error should be small
            // provided radius is not much greater than cell size (which
            // we ensure by making big rigid spheres from lots of small particles).
            if (ip < 0) ip = 0;
            if (jp < 0) jp = 0;
            if (kp < 0) kp = 0;
            if (ip >= eSums.shape(0)-1) ip = eSums.shape(0)-1;
            if (jp >= eSums.shape(1)-1) jp = eSums.shape(1)-1;
            if (kp >= eSums.shape(2)-1) kp = eSums.shape(2)-1;

            eSums(ip, jp, kp) += Vsphereinbox; 

//            if (Vsphereinbox > 1e-12) {
//                std::cout << "overlap cells "  
//                      << i << " " << j << " " << k << " " 
//                      << ip << " " << jp << " " << kp << " "
//                      << i+ic << " " << j+jc << " " << k+kc << " "
//                      << r[0] << " " << r[1] << " " << r[2] << " "
//                      << box[0] << " " << box[1] << " " << box[2] << " " 
//                      << box[3] << " " << box[4] << " " << box[5] <<  " "
//                      << Vsphereinbox << std::endl;
//            }

        }}}
    }


    //std::cout << "pre_force "  << use_overlap << " " << eSums(cell[0], cell[1], cell[2]) << std::endl;
}


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForceDrag::get_force(double r[], double v[], double a[], double m, double s, double e){

    std::vector<double> f(3), Ui(3), Ui_v(3), gradP(3), divStress(4);
    std::vector<int> cell = get_cell(r);
    //CPL::ndArray<double> array = fieldptr->get_array();
    // This is much faster, array cannot exist beyond get_force scope 
    // and fieldptr won't be deleted so no need for smart pointer?
    CPL::ndArray<double>& array = fieldptr->get_array_pointer();

    assert(array.shape(0) == 9);

    //Should use std::vector<double> Ui(3) = field.interpolate(r);
    for (int i=0; i<3; i++){
        Ui[i] = array(i, cell[0], cell[1], cell[2]);
        Ui_v[i] = Ui[i]-v[i];
        gradP[i] = array(i+3, cell[0], cell[1], cell[2]);
        divStress[i] = array(i+6, cell[0], cell[1], cell[2]);

//        std::cout << "U array "  
//                  << cell[0] << " " << cell[1] << " " << cell[2] << " " 
//                  << i << " " << array(i, cell[0], cell[1], cell[2]) << std::endl;
    }

    double Cd = drag_coefficient();
    double cellvolume = fieldptr->dV;
    double volume = (4./3.)*M_PI*pow(s,3); 
    double eps = 1.0 - eSums(cell[0], cell[1], cell[2])/cellvolume;

    // Add sum of coefficients of forces 
    // Needed if you want to split implicit/explicit terms for
    // improved numerical stability according to 
    // Xiao H., Sun J. (2011) Algorithms in a Robust Hybrid
    // CFD-DEM Solver for Particle-Laden Flows, 
    // Commun. Comput. Phys. 9, 2, 297
    FcoeffSums (cell[0], cell[1], cell[2]) += Cd;

    //Calculate force
    for (int i = 0; i < 3; ++i){
        //Just drag force here
        f[i] = Cd*Ui_v[i];
        //Include pressure and stress
        //f[i] += volume*(divStress[i]-gradP[i]);
        //Add to sum of forces
        FSums(i, cell[0], cell[1], cell[2]) += f[i];
    }

    //std::cout << "Drag Force "  
    //          << r[2] << " " << v[0] << " " << Ui[0] << " "  << v[1] << " " << Ui[1] << " " << v[2] << " " << Ui[2] << " " 
    //          << divStress[2] << " " << gradP[2] << " " << f[2] << " "  << std::endl;

    return f;
}


///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForceGranular                           //
//                                                               //
///////////////////////////////////////////////////////////////////
// General Class to use for Ergun (1952), Di Felice (1994) and BVK.

/////////////////////////////////////////////////////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

//Constructor using cells
CPLForceGranular::CPLForceGranular(int nd, int icells, int jcells, int kcells) : CPLForceDrag(nd, icells, jcells, kcells){
    initialisesums(fieldptr->get_array());
}

//Constructor of datatype
CPLForceGranular::CPLForceGranular(CPL::ndArray<double> arrayin) : CPLForceDrag(arrayin){
    initialisesums(arrayin);
}

void CPLForceGranular::initialisesums(CPL::ndArray<double> arrayin){    

    int nsumsShape[3] = {arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    nSums.resize(3, nsumsShape); // Sum of number of particles

    int esumsShape[3] = {arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    eSums.resize(3, esumsShape); // Sum of porosity of particles  

    int FsumsShape[4] = {3, arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    FSums.resize(4, FsumsShape); // Sum of force on particles  

    int vsumsShape[4] = {arrayin.shape(0), arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    vSums.resize(4, vsumsShape); // Sum of velocity
    resetsums();
}

void CPLForceGranular::resetsums(){
    nSums = 0.0; eSums = 0.0; vSums = 0.0; FSums=0.0;
}

// Reynolds Number
// It is unclear here if Reynolds No. should be based
// on the mean cell velocity or particle velocity (or mean/relative from both)
double CPLForceGranular::Reynolds_number(double D, double U, double rho, double mu, double eps) {
    return rho * D * eps * U / mu;
}

// See Equation 12 in K. D. Kafui et al. / Chemical Engineering Science 57 (2002) 2395–2410
double CPLForceGranular::porousity_exponent(double Re) {
    return 3.7 - 0.65 * exp(pow(0.5*(-1.5 - log10(Re)),2));
}

// See Equation 13 in K. D. Kafui et al. / Chemical Engineering Science 57 (2002) 2395–2410
double CPLForceGranular::drag_coefficient(double Re) {
    return pow((0.63 + 4.8/pow(Re,0.5)),2);
}

double CPLForceGranular::magnitude(std::vector<double> v){
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
} 

//Pre force collection of sums (should this come from LAMMPS fix chunk/atom bin/3d)
void CPLForceGranular::pre_force(double r[], double v[], double a[], double m, double s, double e) {

    // Find in which cell number (local to processor) is the particle
    // and sum all the velocities for each cell.
    std::vector<int> cell = get_cell(r);

    // Should use field.add_volume(r, radius);
    double radius = s;
    double volume = (4./3.)*M_PI*pow(radius,3);
    nSums(cell[0], cell[1], cell[2]) += 1.; 
    eSums(cell[0], cell[1], cell[2]) += volume; 
}


//Get force using sums collected in pre force
std::vector<double> CPLForceGranular::get_force(double r[], double v[], double a[], double m, double s, double e) {

    std::vector<double> f(3), Ui(3), Ui_v(3);
    std::vector<int> cell = get_cell(r);
    CPL::ndArray<double>& array = fieldptr->get_array_pointer();

    double radius = s;
    double d = 2.0*radius;
    double volume = (4./3.)*M_PI*pow(radius,3); 

    //Porosity e is cell volume - sum in volume
    double cellvolume = fieldptr->dV;
    double eps = 1.0 - eSums(cell[0], cell[1], cell[2])/cellvolume;
    double rho = 1.0; //m/volume;
    double mu = 1.0;

    //Should use std::vector<double> Ui(3) = field.interpolate(r);
    for (int i=0; i<3; i++){
        Ui[i] = array(i, cell[0], cell[1], cell[2]);
        Ui_v[i] = Ui[i]-v[i];
    }

    //It is unclear here if Reynolds No. should be based
    //on the mean cell velocity or particle velocity
    //double Re = rho * d * eps * magnitude(Ui_v) / mu;
    double Re = Reynolds_number(d, magnitude(Ui_v), rho, mu, eps);
    double Cd = drag_coefficient(Re);
    double xi = porousity_exponent(Re);

    //Calculate force
    if (eps < 1e-5) {
        std::cout << "Warning: 0 particles in cell (" 
                  << cell[0] << ", " << cell[1] << ", " << cell[2] << ")"
                  << std::endl;
        f[0]=0.0; f[1]=0.0; f[2]=0.0;
        return f;
    } else {
        double A = 0.125*Cd*rho*M_PI*pow(d,2)*pow(eps,2)*magnitude(Ui_v)*pow(eps,xi-1.0);
        for (int i = 0; i < 3; ++i){
            f[i] = A*(Ui[i]-v[i]);
            FSums(i, cell[0], cell[1], cell[2]) += f[i];
        }
        return f;
    }
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


///////////////////////////////////////////////////////////////////
//                          BVK                                  //
///////////////////////////////////////////////////////////////////

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

//Return Stokes force
double CPLForceBVK::Stokes(double D, double U, double mu) {
    return 9.42477796076938 * mu * D * U;
}

// Reynolds Number
// It is unclear here if Reynolds No. should be based
// on the mean cell velocity or particle velocity (or mean/relative from both)
double CPLForceBVK::Reynolds_number(double D, double U, double rho, double mu, double eps) {
    return rho * D * eps * U / mu;
}

//Calculate BVK drag force per particle
double CPLForceBVK::CPLForceBVK_expression(double eps, double D, double U, double rho, double mu) {
    double phi = 1 - eps;
    double Re = Reynolds_number(D, U, rho, mu, eps);
    double BVK = 10.0*phi/pow(eps,2)
          + pow(eps,2)*(1.0+1.5*pow(phi,0.5)) 
          + (0.413*Re/(24*pow(eps,2))) 
          * (((1./eps) + 3 * phi * eps + 8.4 * pow(Re,-0.343)) 
             /(1 + pow(10,(3*phi)) * pow(Re,-0.5) *(1+4*phi)));
    return Stokes(D, U, mu) * BVK;
}


//Get force using sums collected in pre force
std::vector<double> CPLForceBVK::get_force(double r[], double v[], double a[], double m, double s, double e) {

    std::vector<double> f(3), Ui(3), Ui_v(3);
    std::vector<int> cell = get_cell(r);
    CPL::ndArray<double>& array = fieldptr->get_array_pointer();

    double radius = s;
    double d = 2.0*radius;
    double volume = (4./3.)*M_PI*pow(radius,3); 

    //Porosity e is cell volume - sum in volume
    double cellvolume = fieldptr->dV;
    double eps = 1.0 - eSums(cell[0], cell[1], cell[2])/cellvolume;
    double rho = 1.0; //m/volume;
    double mu = 1.0;

    //Should use std::vector<double> Ui(3) = field.interpolate(r);
    for (int i=0; i<3; i++){
        Ui[i] = array(i, cell[0], cell[1], cell[2]);
        Ui_v[i] = Ui[i]-v[i];
    }
    double Re = Reynolds_number(d, magnitude(Ui_v), rho, mu, eps);

    //Calculate force
    if (eps < 1e-5) {
        std::cout << "Warning: 0 particles in cell (" 
                  << cell[0] << ", " << cell[1] << ", " << cell[2] << ")"
                  << std::endl;
        f[0]=0.0; f[1]=0.0; f[2]=0.0;
        return f;
    } else {
        double A = CPLForceBVK_expression(eps, d, magnitude(Ui_v), rho, mu);
        for (int i = 0; i < 3; ++i){
            f[i] = A*(Ui[i]-v[i]);
            FSums(i, cell[0], cell[1], cell[2]) += f[i];
        }
        return f;
    }
}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
///////////// WARNING THIS IS CURRENTLY UNTESTED ////////////////
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


//Example print statement useful to copy in...
//    std::cout << "FLEKKOY: " << cell[0] << " " << cell[1]  << " " << cell[2]  
//                << nSums(icell, jcell, kcell) << " " 
//                << min[1] << " " <<  max[1] << " " << r[1]<< std::endl;


