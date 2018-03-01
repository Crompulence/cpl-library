#include <vector>
#include <math.h> 
#include <stdexcept>
#include <assert.h>
#include <map>
#include <string> 

#include "CPL_ndArray.h"
#include "CPL_force.h"
#include "CPL_field.h"
#include "CPL_misclib.h"
#include "overlap/sphere_cube.hpp"

///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForce base class                        //
//                                                               //
///////////////////////////////////////////////////////////////////


//Constructors
CPLForce::CPLForce(int nd, int icells, int jcells, int kcells){ 
    //Shared as we keep reference in fields list
    cfd_array_field = std::make_shared<CPL::CPLField>(nd, icells, jcells, kcells);
}

CPLForce::CPLForce(CPL::ndArray<double> arrayin) {
    //Shared as we keep reference in fields list
    cfd_array_field = std::make_shared<CPL::CPLField>(arrayin);
}

//Pre force collection of sums (should this come from LAMMPS fix chunk/atom bin/3d)
void CPLForce::pre_force(double r[], double v[], double a[], 
                         double m, double s, double e) {
//    throw std::runtime_error("CPLForce::pre_force is not defined");
}


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForce::get_force(double r[], double v[], double a[], 
                                        double m, double s, double e){
//    throw std::runtime_error("CPLForce::get_force is not defined");
    std::vector<double> f = {0.0, 0.0, 0.0};
    return f;
}

void CPLForce::resetsums(){
    //General function which loops over all field classes and resets
    for ( auto &f : fields ) { 
        f->zero_array(); 
    }
}

//Set minimum and maximum values of field application
void CPLForce::set_minmax(double min_in[], double max_in[]){
    cfd_array_field->set_minmax(min_in, max_in);
    for ( auto &f : fields ) { 
        f->set_minmax(min_in, max_in);
        //std::cout << "fields " << typeid(*f).name() << std::endl;
    }
};

//If either min/max change or field object, we need to recalculate dx, dy and dz
void CPLForce::set_dxyz(){
    cfd_array_field->set_dxyz();
    for ( auto &f : fields ) { 
        f->set_dxyz();
    }
}

//Get dA
std::vector<double> CPLForce::get_dA(){
    return cfd_array_field->get_dA();
}

//Get dV
double CPLForce::get_dV(){
    return cfd_array_field->get_dV();
}


//Get cell from min/max and dx
std::vector<int> CPLForce::get_cell(double r[]){
    std::vector<int> cell(3);
    cell = cfd_array_field->get_cell(r);
    return cell;
}

//This function sets the main field obtained from the CFD 
void CPLForce::set_field(CPL::ndArray<double> arrayin){
    cfd_array_field->set_array(arrayin);   
};

//Get main CFD field
CPL::ndArray<double> CPLForce::get_field(){
    return cfd_array_field->get_array();
};

std::shared_ptr<CPL::CPLField> CPLForce::get_internal_fields(const std::string& name){
    //std::cout << "Looking for internal_fields: " << name << " in list of size " << fields.size() << std::endl;
    for ( auto &f : fields ) {
        //std::cout << "getting_internal_fields " << name << " " << f->field_name() << std::endl;
        if (f->field_name() == name){
            return f;
        }
    }
    //std::shared_ptr<CPL::CPLField> nullpointer();
    return nullptr;
}

///////////////////////////////////////////////////////////////////
//                                                               //
//                      CPLForceTest                             //
//                                                               //
///////////////////////////////////////////////////////////////////

//Constructor using cells
CPLForceTest::CPLForceTest(int nd, int icells, int jcells, int kcells)
     : CPLForce(nd, icells, jcells, kcells){
    initialisesums(cfd_array_field->get_array());
}

//Constructor of datatype
CPLForceTest::CPLForceTest(CPL::ndArray<double> arrayin) : CPLForce(arrayin){
    initialisesums(arrayin);
}

//Examples of initialise sums with another internal field
void CPLForceTest::initialisesums(CPL::ndArray<double> arrayin){
    
    auto otherfield = std::make_shared<CPL::CPLField>(1, arrayin.shape(1), 
                                                         arrayin.shape(2), 
                                                         arrayin.shape(3), 
                                                         "otherfield");
    fields.push_back(otherfield);
    resetsums();
}

void CPLForceTest::resetsums(){
    //General function which loops over all field classes and resets
    for ( auto &f : fields ) { 
        f->zero_array(); 
    }
}


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
void CPLForceTest::pre_force(double r[], double v[], double a[], 
                             double m, double s, double e){
//    throw std::runtime_error("CPLForceTest::pre_force is not defined");
}


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForceTest::get_force(double r[], double v[], double a[], 
                                            double m, double s, double e){

    std::vector<int> indices = {0,1,2};
    std::vector<double> f = cfd_array_field->get_array_value(indices, r);

//    std::cout << "CPLForceTest " << f[0] << " " << f[1] << " " << f[2] << std::endl;

    return f;
}


///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForceVelocity                           //
//                                                               //
///////////////////////////////////////////////////////////////////

//Constructor using cells
CPLForceVelocity::CPLForceVelocity(int nd, int icells, int jcells, int kcells) 
    : CPLForce(nd, icells, jcells, kcells){
    initialisesums(cfd_array_field->get_array());
}

//Constructor of datatype
CPLForceVelocity::CPLForceVelocity(CPL::ndArray<double> arrayin) 
    : CPLForce(arrayin){
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
void CPLForceVelocity::pre_force(double r[], double v[], double a[], 
                                 double m, double s, double e) {

    // Find in which cell number (local to processor) is the particle
    // and sum all the velocities for each cell.
    std::vector<int> cell = get_cell(r);

    nSums(cell[0], cell[1], cell[2]) += 1.0; 
    for (int i = 0; i<3; ++i)
        vSums(i, cell[0], cell[1], cell[2]) += v[i];

}


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForceVelocity::get_force(double r[], double v[], double a[], 
                                                double m, double s, double e){

    std::vector<double> f(3); 
    std::vector<int> cell = get_cell(r);
    CPL::ndArray<double> array = cfd_array_field->get_array();
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
CPLForceFlekkoy::CPLForceFlekkoy(int nd, int icells, int jcells, int kcells) 
    : CPLForce(nd, icells, jcells, kcells){
    initialisesums(cfd_array_field->get_array());
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
void CPLForceFlekkoy::pre_force(double r[], double v[], double a[], 
                                double m, double s, double e) {

    // Find in which cell number (local to processor) is the particle
    // and sum all the Flekkøy weights for each cell.
    std::vector<int> cell = get_cell(r);

    double g = flekkoyGWeight(r[1], cfd_array_field->min[1], cfd_array_field->max[1]);
    nSums(cell[0], cell[1], cell[2]) += 1.0; 
    gSums(cell[0], cell[1], cell[2]) += g;

}

//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForceFlekkoy::get_force(double r[], double v[], double a[], 
                                               double m, double s, double e){

    std::vector<double> dA = cfd_array_field->get_dA(); 
    std::vector<double> f = {0.0, 0.0, 0.0}; 
    std::vector<int> cell = get_cell(r);
    CPL::ndArray<double> array = cfd_array_field->get_array();
    double n = nSums(cell[0], cell[1], cell[2]);

    if (n < 1.0) {
        std::cout << "Warning: 0 particles in cell (" 
                  << cell[0] << ", " << cell[1] << ", " << cell[2] << ")"
                  << std::endl;
    } else {

        // Since the Flekkoy weight is considered only for 0 < y < L/2, for cells 
        // that are completely in y < 0 gSums(i, j, k) will be always 0.0 so can 
        // produce a NAN in the g/gSums division below.
        double g = flekkoyGWeight(r[1], cfd_array_field->min[1], cfd_array_field->max[1]);

        if (gSums(cell[0], cell[1], cell[2]) > 0.0) {
            double gdA = (g/gSums(cell[0], cell[1], cell[2])) * dA[1];

            // Normal to the X-Z plane is (0, 1, 0) so (tauxy, syy, tauxy)
            // are the only components of the stress tensor that matter.
            f[0] = gdA * array(1, cell[0], cell[1], cell[2]);
            f[1] = gdA * array(4, cell[0], cell[1], cell[2]);
            f[2] = gdA * array(7, cell[0], cell[1], cell[2]);
        }
    }

    return f;
}


///////////////////////////////////////////////////////////////////
//                                                               //
//                      CPLForceDrag                             //
//                                                               //
///////////////////////////////////////////////////////////////////
// This has the form of a further abstract base class for drag 
// based coupled examples, defining flags, appropriate internal
// fields and necessary force functions.

//Constructor using cells
CPLForceDrag::CPLForceDrag(int nd, int icells, int jcells, int kcells) 
    : CPLForce(nd, icells, jcells, kcells){
    initialisesums(cfd_array_field->get_array());
}

//Constructor of datatype
CPLForceDrag::CPLForceDrag(CPL::ndArray<double> arrayin) 
    : CPLForce(arrayin){
    initialisesums(arrayin);
}

//Constructor of datatype
CPLForceDrag::CPLForceDrag(CPL::ndArray<double> arrayin, 
                           map_strstr arg_map) 
    : CPLForce(arrayin)
{
    unpack_arg_map(arg_map);
    initialisesums(arrayin);
}

//Constructor with optional argument overlap, interpolate and drag_coeff
CPLForceDrag::CPLForceDrag(int nd, int icells, int jcells, int kcells, 
                           map_strstr arg_map)
    : CPLForce(nd, icells, jcells, kcells)
{
    unpack_arg_map(arg_map);
    initialisesums(cfd_array_field->get_array());
}


void CPLForceDrag::initialisesums(CPL::ndArray<double> arrayin){
    
    int i = arrayin.shape(1);
    int j = arrayin.shape(2);
    int k = arrayin.shape(3);
    nSums = std::make_shared<CPL::CPLField>(1, i, j, k, "nSums");
    vSums = std::make_shared<CPL::CPLField>(3, i, j, k, "vSums");
    eSums = std::make_shared<CPL::CPLField>(1, i, j, k, "eSums");
    FSums = std::make_shared<CPL::CPLField>(3, i, j, k, "FSums");
    FcoeffSums = std::make_shared<CPL::CPLField>(1, i, j, k, "FcoeffSums");

    build_fields_list();
    resetsums();
}

void CPLForceDrag::build_fields_list(){

    fields.push_back(nSums);
    fields.push_back(vSums);
    fields.push_back(eSums);
    fields.push_back(FSums);
    fields.push_back(FcoeffSums);
}

//Can be inhereted from Base Class
void CPLForceDrag::resetsums(){

    for ( auto &f : fields ) { 
        f->zero_array(); 
    }
}

// Split into protected functions for a default set of inputs
// and an extended set which, if defined, will supplement.
void CPLForceDrag::unpack_arg_map(map_strstr arg_map){
    bool extra_args = unpack_extra_arg_map(arg_map);
    unpack_default_arg_map(arg_map, extra_args);
}



void CPLForceDrag::unpack_default_arg_map(map_strstr arg_map, bool extra_args){

//   std::cout << "BEFORE "
//                << "use_overlap " << use_overlap << " " 
//                << "use_interpolate " << use_interpolate << " " 
//                << "Cd " << Cd << " " 
//                << "mu " << mu << " " 
//                << "rho " << rho << " " 
//                << "use_gradP " << use_gradP << " " 
//                << "use_divStress " << use_divStress << std::endl;


    // Iterate over the map and print out all key/value pairs.
    for (const auto& arg : arg_map)
    {
        if (string_contains(arg.first, "overlap") != -1) {
            use_overlap = checktrue(arg.second);
        } else if (string_contains(arg.first, "interpolate") != -1) {
            use_interpolate = checktrue(arg.second);
        } else if (string_contains(arg.first, "gradP")  != -1) {
            use_gradP = checktrue(arg.second);
        } else if (string_contains(arg.first, "divStress")  != -1) {
            use_divStress = checktrue(arg.second);
        } else if (string_contains(arg.first, "Cd")  != -1) {
            Cd = std::stod(arg.second);
        } else if (string_contains(arg.first, "mu")  != -1) {
            mu = std::stod(arg.second);
        } else if (string_contains(arg.first, "rho")  != -1) {
            rho = std::stod(arg.second);
        } else if (extra_args) {
            std::cout << "key: " << arg.first << 
            " for forcetype not recognised as default" << '\n';
            throw std::runtime_error("CPLForceDrag input not recognised");
        }
    }

//    std::cout << "AFTER " 
//                << "use_overlap " << use_overlap << " " 
//                << "use_interpolate " << use_interpolate << " " 
//                << "Cd " << Cd << " " 
//                << "mu " << mu << " " 
//                << "rho " << rho << " " 
//                << "use_gradP " << use_gradP << " " 
//                << "use_divStress " << use_divStress << std::endl;

}

// This is the function which should be overridden
bool CPLForceDrag::unpack_extra_arg_map(map_strstr arg_map){

    // Iterate over the map and print out all key/value pairs.
    bool extra_args = false;
    for (const auto& arg : arg_map)
    {
//        std::cout << "key: " << arg.first;
//        std::cout << " value: " << arg.second
//                  << " checktrue: " << checktrue(arg.second) << '\n';
        if (string_contains(arg.first, "some_extra_argument") != -1) {
            bool extra_args_name = checktrue(arg.second);
            extra_args = true; //Set this to disable default check
        }
    }
    return extra_args;

}

//Unpack various CFD values from array into a range of fields.
//Currently we use an array with all values "cfd_array_field".
//void CPLForceDrag::unpack_CFD_array(CPL::ndArray<double> arrayin){

//    CPL::ndArray<double> CFDarrayu, CFDarrayp, CFDarrays;
//    int arrayShape[4] = {3, arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
//    CFDarrayu.resize(4, arrayShape);
//    CFDarrayp.resize(4, arrayShape);
//    CFDarrays.resize(4, arrayShape);
//    for (int n=0; n<CFDarrayu.shape(0); n++){
//    for (int i=0; i<CFDarrayu.shape(1); i++){
//    for (int j=0; j<CFDarrayu.shape(2); j++){
//    for (int k=0; k<CFDarrayu.shape(3); k++){
//        CFDarrayu(n,i,j,k) = arrayin(n,i,j,k);
//        CFDarrayp(n,i,j,k) = arrayin(n+3,i,j,k);
//        CFDarrays(n,i,j,k) = arrayin(n+6,i,j,k);
//    }}}}
//    UCFD = std::make_shared<CPL::CPLField>(CFDarrayu, "UCFD");
//    gradPCFD = std::make_shared<CPL::CPLField>(CFDarrayp, "gradPCFD");
//    gradStressCFD = std::make_shared<CPL::CPLField>(CFDarrays, "gradStressCFD");
//}

//Arbitary constant in this case so r and D are irrelevant
double CPLForceDrag::drag_coefficient(double r[], double D, 
                                      std::vector<double> Ui_v) {
    std::cout << "DRAG Cd: " << Cd << std::endl;
    return Cd; //Use default Drag value
}

//Pre force collection of sums including overlap code to assign volumes
void CPLForceDrag::pre_force(double r[], double v[], double a[], 
                             double m, double s, double e) {

    double volume = (4./3.)*M_PI*pow(s,3); 
    nSums->add_to_array(r, 1.0);
    std::cout << "Pre_force " << r[0] << " " << r[1] << " " << r[2] << " " <<volume << std::endl;
    if (! use_overlap){
        eSums->add_to_array(r, volume);
        vSums->add_to_array(r, v);
    } else {
        eSums->add_to_array(r, s, volume);
        vSums->add_to_array(r, s, v);
    }

}


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForceDrag::get_force(double r[], double v[], double a[], 
                                            double m,   double s,   double e){


    //Define variable
    std::vector<double> f(3), Ui(3), Ui_v(3), gradP(3), divStress(4);

    //Check array is the right size
    CPL::ndArray<double>& array = cfd_array_field->get_array_pointer();
    assert(array.shape(0) == 9);

    //Get all elements of recieved field
    if (! use_interpolate){
        //Based on cell
        std::vector<int> indices = {0,1,2}; 
        Ui = cfd_array_field->get_array_value(indices, r);
        for (int &n : indices) n += 3; 
        gradP = cfd_array_field->get_array_value(indices, r);
        for (int &n : indices) n += 3; 
        divStress = cfd_array_field->get_array_value(indices, r);
    } else {
        //Or interpolate to position in space
        std::vector<int> indices = {0,1,2}; 
        Ui = cfd_array_field->get_array_value_interp(indices, r);
        for (int &n : indices) n += 3; 
        gradP = cfd_array_field->get_array_value_interp(indices, r);
        for (int &n : indices) n += 3; 
        divStress = cfd_array_field->get_array_value_interp(indices, r);
    }

    //Get uCFD - uDEM
    for (int i=0; i<3; i++){
        Ui_v[i] = Ui[i]-v[i];
    }

    //Get Diameter
    double D = 2.0*s;

    //Get drag coefficient
    double A = drag_coefficient(r, D, Ui_v);

    //Calculate force
    double volume = (4./3.)*M_PI*pow(s,3); 
    for (int i = 0; i < 3; ++i){
        //Just drag force here
        f[i] = A*Ui_v[i];
        //Include pressure
        if (use_gradP)
            f[i] += -volume*gradP[i];
        // and stress
        if (use_divStress)
            f[i] += volume*divStress[i];
        //std::cout << "cell "  <<  cell[0] << " " << cell[1] << " " << cell[2] << std::endl;
    }

    // Add sum of coefficients of forces 
    // Needed if you want to split implicit/explicit terms for
    // improved numerical stability according to 
    // Xiao H., Sun J. (2011) Algorithms in a Robust Hybrid
    // CFD-DEM Solver for Particle-Laden Flows, 
    // Commun. Comput. Phys. 9, 2, 297
    FcoeffSums->add_to_array(r, A);
    FSums->add_to_array(r, &f[0]);

//    std::cout << "Drag Force "  
//              << r[2] << " " << v[0] << " " << Ui[0] << " "  << v[1] << " " << Ui[1] << " " << v[2] << " " << Ui[2] << " " 
//              << divStress[2] << " " << gradP[2] << " " << f[2] << " "  << std::endl;

    return f;
}


///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForceGranular                           //
//                                                               //
///////////////////////////////////////////////////////////////////
// General Class to use for Ergun (1952), Di Felice (1994), BVK...

// Reynolds Number
// It is unclear here if Reynolds No. should be based
// on the mean cell velocity or particle velocity (or mean/relative from both)
double CPLForceGranular::Reynolds_number(double D, double U, double rho,
                                         double mu, double eps) {
    return rho * D * eps * U / mu;
}

double CPLForceGranular::magnitude(std::vector<double> v){
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
} 

double CPLForceGranular::get_eps(double r[]){
    //Porosity e is cell volume - sum in volume
    double eps = 1.0 - eSums->get_array_value(r)/eSums->get_dV();
    std::cout << "get eps " << eSums->get_array_value(r) << " " << 
                       eSums->get_dV() << " " << eps << std::endl;
        
    if (eps < 1e-5) {
        std::vector<int> cell = get_cell(r);
        std::cout << "Warning: 0 particles in cell (" 
                  << cell[0] << ", " << cell[1] << ", " << cell[2] << ")"
                  << std::endl;
        return 0.0;
    } else {
        return eps;
    }
}

////This combined all the above
//double CPLForceGranular::drag_coefficient(double r[], double D, 
//                                        std::vector<double> Ui_v) {
//    double eps = get_eps(r);
//    std::cout  << "Granular: " << mu*(1.-eps)*(150*(1-eps))/(pow(D, 2.)*eps) << " " << eps << " " << mu << " " << D  << std::endl;
//    return mu*(1.-eps)*(150*(1-eps))/(pow(D, 2.)*eps) ;

//}

///////////////////////////////////////////////////////////////////
//                      Stokes                                  //
///////////////////////////////////////////////////////////////////


//This combined all the above
double CPLForceStokes::drag_coefficient(double r[], double D, 
                                        std::vector<double> Ui_v) {

    std::cout  << "Stokes: " << 9.42477796076938 * mu * D << std::endl;
    return 9.42477796076938 * mu * D;

}

///////////////////////////////////////////////////////////////////
//                      Di_Felice                                //
///////////////////////////////////////////////////////////////////


// See Equation 12 in K. D. Kafui et al. / Chemical Engineering Science 57 (2002) 2395–2410
double CPLForceDi_Felice::porousity_exponent(double Re) {
    return 3.7 - 0.65 * exp(pow(0.5*(-1.5 - log10(Re)),2));
}

// See Equation 13 in K. D. Kafui et al. / Chemical Engineering Science 57 (2002) 2395–2410
double CPLForceDi_Felice::drag_coefficient_Re(double Re) {
    return pow((0.63 + 4.8/pow(Re,0.5)),2);
}


//This combined all the above
double CPLForceDi_Felice::drag_coefficient(double r[], double D, 
                                          std::vector<double> Ui_v) {

    double Re, A, xi, U, eps;

    eps = get_eps(r);
    U = CPLForceGranular::magnitude(Ui_v);
    if (U > 1e-8) {
        Re = CPLForceGranular::Reynolds_number(D, U, rho, mu, eps);
        A = drag_coefficient_Re(Re);
        xi = porousity_exponent(Re);
    } else {
        return 0;
    }

    std::cout  << "Di_Felice: " << eps << " " << Re << " " << A << " " << 
               xi << " " << 0.125*A*rho*M_PI*pow(D,2)*pow(eps,2)*U*pow(eps,xi-1.0) << std::endl;
    if (eps == 0.0) {
        return 0.0;
    } else {
        return 0.125*A*rho*M_PI*pow(D,2)*pow(eps,2)*U*pow(eps,xi-1.0);
    }

}

///////////////////////////////////////////////////////////////////
//                          BVK                                  //
///////////////////////////////////////////////////////////////////

//Return Stokes force
double CPLForceBVK::Stokes(double D, double U, double mu) {
    return 9.42477796076938 * mu * D * U;
}
//Calculate BVK drag force per particle
double CPLForceBVK::CPLForceBVK_expression(double eps, double D, double U, 
                                           double rho, double mu) 
{
    double phi = 1 - eps;
    double Re = CPLForceGranular::Reynolds_number(D, U, rho, mu, eps);
    double BVK = 10.0*phi/pow(eps,2)
          + pow(eps,2)*(1.0+1.5*pow(phi,0.5)) 
          + (0.413*Re/(24*pow(eps,2))) 
          * (((1./eps) + 3 * phi * eps + 8.4 * pow(Re,-0.343)) 
             /(1 + pow(10,(3*phi)) * pow(Re,-0.5) *(1+4*phi)));
    std::cout  << "BVK: " << Stokes(D, U, mu) * BVK << std::endl;
    return Stokes(D, U, mu) * BVK;
}


double CPLForceBVK::drag_coefficient(double r[], double D, 
                                     std::vector<double> Ui_v) {

    double eps = CPLForceGranular::get_eps(r);
    if (eps == 0.0) {
        return 0.0;
    } else {
        return CPLForceBVK_expression(eps, D, 
                                      CPLForceGranular::magnitude(Ui_v), 
                                      rho, mu);
    }

}

///////////////////////////////////////////////////////////////////
//                          Ergun                                //
///////////////////////////////////////////////////////////////////


double CPLForceErgun::drag_coefficient(double r[], double D, 
                                       std::vector<double> Ui_v) {
    double eps = CPLForceGranular::get_eps(r);
    if (eps == 0.0) {
        return 0.0;
    } else {
        std::cout  << "Ergun: " << 150.0*eps*(mu/rho)*rho/(pow(eps*D, 2.0)) + 1.75*rho/(eps*D) << std::endl;
        return 150.0*eps*(mu/rho)*rho/(pow(eps*D, 2.0)) + 1.75*rho/(eps*D);
    }
}

//Example print statement useful to copy in...
//    std::cout << "FLEKKOY: " << cell[0] << " " << cell[1]  << " " << cell[2]  
//                << nSums(icell, jcell, kcell) << " " 
//                << min[1] << " " <<  max[1] << " " << r[1]<< std::endl;


