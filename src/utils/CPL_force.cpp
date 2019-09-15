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
    calc_preforce = false;
    calc_preforce_everytime = false;
    calc_postforce = false;
    calc_postforce_everytime = false;
}


CPLForce::CPLForce(CPL::ndArray<double> arrayin) {
    //Shared as we keep reference in fields list
    cfd_array_field = std::make_shared<CPL::CPLField>(arrayin);
    calc_preforce = false;
    calc_preforce_everytime = false;
    calc_postforce = false;
    calc_postforce_everytime = false;

}

//Pre force collection of sums (for things not possible with LAMMPS fix chunk/atom bin/3d)
void CPLForce::pre_force(double r[], double v[], double a[], 
                         double m, double s, double e) {
//    throw std::runtime_error("CPLForce::pre_force is not defined");
}


//Calculate force on a molecule
std::vector<double> CPLForce::get_force(double r[], double v[], double a[], 
                                        double m, double s, double e){
//    throw std::runtime_error("CPLForce::get_force is not defined");
    std::vector<double> f = {0.0, 0.0, 0.0};
    return f;
}

//Post force collection of sums  (for things not possible with LAMMPS fix chunk/atom bin/3d)
void CPLForce::post_force(double r[], double v[], double a[], 
                         double m, double s, double e) {
//    throw std::runtime_error("CPLForce::post_force is not defined");
}

void CPLForce::resetsums(){

    //Reset all counters
    Npre_force=0; Nforce=0; Npost_force=0;
    //General function which loops over all field classes and resets
    for ( auto &f : fields ) { 
        f->zero_array(); 
    }
}

//Reset any fields which are only kept for
//a single timestep
void CPLForce::reset_instant(){
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


//void CPLForce::setup_fast_array()
//{

//    int nc = cfd_array_field->array.shapeData()[0];
//    int ic = cfd_array_field->array.shapeData()[1];
//    int jc = cfd_array_field->array.shapeData()[2];
//    int kc = cfd_array_field->array.shapeData()[3];

//    //if (nc*ic*jc*kc < 1000000){
//        double *fastarray = new double[9][32][32][32];
//    //}

//    for (int i0=0; i0<nc; i0++) {
//    for (int i1=0; i1<ic; i1++) {
//    for (int i2=0; i2<jc; i2++) {
//    for (int i3=0; i3<kc; i3++) {
//        fastarray[i0][i1][i2][i3] = cfd_array_field->get_array_value(i0, i1, i2, i3);
//    }}}}

//    delete[] fastarray;

//}


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
    calc_preforce = false;
    calc_preforce_everytime = false;
    calc_postforce = false;
    calc_postforce_everytime = false;
    auto otherfield = std::make_shared<CPL::CPLField>(1, arrayin.shape(1), 
                                                         arrayin.shape(2), 
                                                         arrayin.shape(3), 
                                                         "otherfield");
    fields.push_back(otherfield);
    resetsums();
}

void CPLForceTest::resetsums(){
    //Reset all counters
    Npre_force=0; Nforce=0; Npost_force=0;
    //General function which loops over all field classes and resets
    for ( auto &f : fields ) { 
        f->zero_array(); 
    }
}

//Reset any fields which are only kept for
//a single timestep
void CPLForceTest::reset_instant(){
}


//Pre force collection of sums 
void CPLForceTest::pre_force(double r[], double v[], double a[], 
                             double m, double s, double e){
//    throw std::runtime_error("CPLForceTest::pre_force is not defined");
}


//Calculate force on a molecule
std::vector<double> CPLForceTest::get_force(double r[], double v[], double a[], 
                                            double m, double s, double e){

    std::vector<int> indices = {0,1,2};
    std::vector<double> f = cfd_array_field->get_array_value(indices, r);

    //std::cout << "CPLForceTest " << r[0] << " " << r[1] << " " << r[2] << " " << f[0] << " " << f[1] << " " << f[2] << std::endl;

    return f;
}

//Post force collection of sums  (for things not possible with LAMMPS fix chunk/atom bin/3d)
void CPLForceTest::post_force(double r[], double v[], double a[], 
                         double m, double s, double e) {
//    throw std::runtime_error("CPLForce::post_force is not defined");
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

//Constructor with optional argument overlap, interpolate and drag_coeff
CPLForceVelocity::CPLForceVelocity(int nd, int icells, int jcells, int kcells, 
                           map_strstr arg_map)
    : CPLForce(nd, icells, jcells, kcells)
{
    unpack_arg_map(arg_map);
    initialisesums(cfd_array_field->get_array());
}

//Constructor of datatype
CPLForceVelocity::CPLForceVelocity(CPL::ndArray<double> arrayin, 
                           map_strstr arg_map) 
    : CPLForce(arrayin)
{
    unpack_arg_map(arg_map);
    initialisesums(arrayin);
}


void CPLForceVelocity::unpack_arg_map(map_strstr arg_map){

    // Iterate over the map and print out all key/value pairs.
    for (const auto& arg : arg_map)
    {
        if (string_contains(arg.first, "xi") != -1) {
            xi = std::stod(arg.second);
        } else {
            std::cout << "key: " << arg.first << 
            " for CPLForceVelocity not recognised" << '\n';
            throw std::runtime_error("CPLForceVelocity input not recognised");
        }
    }
}



void CPLForceVelocity::initialisesums(CPL::ndArray<double> arrayin){
    calc_preforce = true;
    calc_preforce_everytime = true;
    calc_postforce = false;
    calc_postforce_everytime = false;
    int n = arrayin.shape(0);
    int i = arrayin.shape(1);
    int j = arrayin.shape(2);
    int k = arrayin.shape(3);
    shapeVector[0] = n;
    shapeVector[1] = i;
    shapeVector[2] = j;
    shapeVector[3] = k;
    nSums = std::make_shared<CPL::CPLField>(1, i, j, k, "nSums");
    nSums_mdt = std::make_shared<CPL::CPLField>(1, i, j, k, "nSums_mdt");
    vSums = std::make_shared<CPL::CPLField>(3, i, j, k, "vSums");
    vSums_mdt = std::make_shared<CPL::CPLField>(3, i, j, k, "vSums_mdt");

    //Get pointer to recved fields
    array = cfd_array_field->get_array_pointer();
    fields.push_back(nSums);
    fields.push_back(vSums);
    resetsums();

//    int vsumsShape[4] = {arrayin.shape(0), arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
//    vSums.resize(4, vsumsShape); // Sum of velocity
//    int nsumsShape[3] = {arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
//    nSums.resize(3, nsumsShape); // Sum of number of particles  
//    resetsums();

}

void CPLForceVelocity::resetsums(){
    //Deep copy velocity to t + dt array
    nSums_mdt = std::make_shared<CPL::CPLField>(*(nSums.get()));
    //nSums_mdt = nSums;
    vSums_mdt = std::make_shared<CPL::CPLField>(*(vSums.get()));
    //vSums_mdt = vSums;
    //Reset all counters
    Npre_force=0; Nforce=0; Npost_force=0;
    //General function which loops over all field classes and resets
    for ( auto &f : fields ) { 
        f->zero_array(); 
    }
}

//Reset any fields which are only kept for
//a single timestep
void CPLForceVelocity::reset_instant(){
}


//Pre force collection of sums (should this come from LAMMPS fix chunk/atom bin/3d)
void CPLForceVelocity::pre_force(double r[], double v[], double a[], 
                                 double m, double s, double e) {

    // Find in which cell number (local to processor) is the particle
    // and sum all the velocities for each cell.
    //std::vector<int> cell = get_cell(r);

//    if (e == 555){
//        std::cout << "CPLForceVelocity::pre_force " << m << " " << 
//                    r[0] << " " << r[1] << " " << r[2] << " " << 
//                    v[0] << " " << v[1] << " " << v[2] <<  std::endl;
//    }

    nSums->add_to_array(r, 1.0);
    vSums->add_to_array(0, r, v[0]);
    vSums->add_to_array(1, r, v[1]);
    vSums->add_to_array(2, r, v[2]);

    //nSums(cell[0], cell[1], cell[2]) = nSums(cell[0], cell[1], cell[2]) + 1.0; 
    //for (int i = 0; i<3; ++i)
    //    vSums->add_to_array(r, v[i]);
        //vSums(i, cell[0], cell[1], cell[2]) = vSums(i, cell[0], cell[1], cell[2]) + v[i];

}


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForceVelocity::get_force(double r[], double v[], double a[], 
                                                double m, double s, double e){

    std::vector<double> f(3), UCFD(3), vsum(3), vsum_mdt(3); 
    std::vector<int> cell = get_cell(r);
    //CPL::ndArray<double> array = cfd_array_field->get_array();

    UCFD[0] = cfd_array_field->get_array_value(0, cell[0], cell[1], cell[2]);
    UCFD[1] = cfd_array_field->get_array_value(1, cell[0], cell[1], cell[2]);
    UCFD[2] = cfd_array_field->get_array_value(2, cell[0], cell[1], cell[2]);

    //Take values from this timestep
    int N = nSums->get_array_value(0, cell[0], cell[1], cell[2]);
    vsum[0] = vSums->get_array_value(0, cell[0], cell[1], cell[2]);
    vsum[1] = vSums->get_array_value(1, cell[0], cell[1], cell[2]);
    vsum[2] = vSums->get_array_value(2, cell[0], cell[1], cell[2]);

    //Or values from previous timestep
    int N_mdt = nSums_mdt->get_array_value(0, cell[0], cell[1], cell[2]);
    vsum_mdt[0] = vSums_mdt->get_array_value(0, cell[0], cell[1], cell[2]);
    vsum_mdt[1] = vSums_mdt->get_array_value(1, cell[0], cell[1], cell[2]);
    vsum_mdt[2] = vSums_mdt->get_array_value(2, cell[0], cell[1], cell[2]);

//    //std::cout << "CPLForceVelocity::get_force " << m << " " << 
//    //            r[0] << " " << r[1] << " " << r[2] << " " << v[0] <<  std::endl;


//    std::cout << "UCFD[0] " << UCFD[0]  << " UCFD[1] " << UCFD[1] << " UCFD[2] " << UCFD[2] 
//              << " vsum[0] " << vsum[0] << " vsum[1] " << vsum[1] << " vsum[2] " << vsum[2]  
//              << " vsum_mdt[0] " << vsum_mdt[0] << " vsum_mdt[1] " 
//              << vsum_mdt[1] << " vsum_mdt[2] " << vsum_mdt[2] << std::endl;

    if (N_mdt < 1.0) {
//        std::cout << "Warning: 0 particles in cell (" 
//                  << cell[0] << ", " << cell[1] << ", " << cell[2] << ")"
//                  << std::endl;
        f[0]=0.0; f[1]=0.0; f[2]=0.0;
    } else {
        for (int i=0; i<3; i++){
            //f[i] = xi*(UCFD[i] - vsum[i]/N);
            f[i] = xi*(UCFD[i] - vsum_mdt[i]/N_mdt);
        }
    }

    //std::cout << "UCFD[0] " << UCFD[0]  << " vsum[0]/N " << vsum[0]/N << " f[0] = " << f[0] << std::endl;
    return f;
}

//Post force collection of sums  (for things not possible with LAMMPS fix chunk/atom bin/3d)
void CPLForceVelocity::post_force(double r[], double v[], double a[], 
                                  double m, double s, double e) {
//    throw std::runtime_error("CPLForce::post_force is not defined");
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
    calc_preforce = true;
    calc_preforce_everytime = true;
    calc_postforce = false;
    calc_postforce_everytime = false;
    int sumsShape[3] = {arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    gSums.resize(3, sumsShape); // Sum of Flekkøy g weights
    nSums.resize(3, sumsShape); // Sum of number of particles  
    resetsums();
}

void CPLForceFlekkoy::resetsums(){
    //Reset all counters
    Npre_force=0; Nforce=0; Npost_force=0;
    //Reset sums
    nSums = 0.0;  gSums = 0.0;
}

//Reset any fields which are only kept for
//a single timestep
void CPLForceFlekkoy::reset_instant(){
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

//Post force collection of sums  (for things not possible with LAMMPS fix chunk/atom bin/3d)
void CPLForceFlekkoy::post_force(double r[], double v[], double a[], 
                         double m, double s, double e) {
//    throw std::runtime_error("CPLForce::post_force is not defined");
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
    : CPLForce(nd, icells, jcells, kcells)
{
    set_defaults();
    initialisesums(cfd_array_field->get_array());
    initialise_extrasums(cfd_array_field->get_array());

}

//Constructor of datatype
CPLForceDrag::CPLForceDrag(CPL::ndArray<double> arrayin) 
    : CPLForce(arrayin)
{
    set_defaults();
    initialisesums(arrayin);
    initialise_extrasums(arrayin);
}

//Constructor of datatype
CPLForceDrag::CPLForceDrag(CPL::ndArray<double> arrayin, 
                           map_strstr arg_map) 
    : CPLForce(arrayin)
{
    set_defaults();
    unpack_arg_map(arg_map);
    initialisesums(arrayin);
    initialise_extrasums(arrayin);
}

//Constructor with optional argument overlap, interpolate and drag_coeff
CPLForceDrag::CPLForceDrag(int nd, int icells, int jcells, int kcells, 
                           map_strstr arg_map)
    : CPLForce(nd, icells, jcells, kcells)
{
    set_defaults();
    unpack_arg_map(arg_map);
    initialisesums(cfd_array_field->get_array());
    initialise_extrasums(cfd_array_field->get_array());
}

void CPLForceDrag::set_defaults(){

    use_overlap = true;
    use_interpolate = false;
    use_gradP = true;
    use_divStress = false;
    calc_preforce_everytime = false;
    calc_postforce_everytime = false;

}


void CPLForceDrag::initialisesums(CPL::ndArray<double> arrayin){

    //Default values
    calc_preforce = true;
    calc_postforce = false;
    int n = arrayin.shape(0);
    int i = arrayin.shape(1);
    int j = arrayin.shape(2);
    int k = arrayin.shape(3);
    shapeVector[0] = n;
    shapeVector[1] = i;
    shapeVector[2] = j;
    shapeVector[3] = k;
    nSums = std::make_shared<CPL::CPLField>(1, i, j, k, "nSums");
    vSums = std::make_shared<CPL::CPLField>(3, i, j, k, "vSums");
    volSums = std::make_shared<CPL::CPLField>(1, i, j, k, "volSums");
    instant_volSums = std::make_shared<CPL::CPLField>(1, i, j, k, "instant_volSums");
    FSums = std::make_shared<CPL::CPLField>(3, i, j, k, "FSums");
    FcoeffSums = std::make_shared<CPL::CPLField>(1, i, j, k, "FcoeffSums");

    //Get pointer to recved fields
    array = cfd_array_field->get_array_pointer();
    build_fields_list();
    resetsums();
}

void CPLForceDrag::initialise_extrasums(CPL::ndArray<double> arrayin){
    //Default empty, inherit and use form:
    // 1) Define new std::make_shared<CPL::CPLField> objects
    // 2) Set them to zero (or whatever initial value you want)
    // 3) Add them to fields array with fields.push_back
}


void CPLForceDrag::build_fields_list(){

    fields.push_back(nSums);
    fields.push_back(vSums);
    fields.push_back(volSums);
    fields.push_back(instant_volSums);
    fields.push_back(FSums);
    fields.push_back(FcoeffSums);
}

//Can be inhereted from Base Class
void CPLForceDrag::resetsums(){
    //Reset all counters
    Npre_force=0; Nforce=0; Npost_force=0;
    //Reset sums
    for ( auto &f : fields ) { 
        f->zero_array(); 
    }
}

//Reset any fields which are only kept for
//a single timestep
void CPLForceDrag::reset_instant(){
    instant_volSums->zero_array();  
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
        } else if (string_contains(arg.first, "preforce_everytime")  != -1) {
            calc_preforce_everytime = checktrue(arg.second);
        } else if (string_contains(arg.first, "postforce_everytime")  != -1) {
            calc_postforce_everytime = checktrue(arg.second);
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
            //bool extra_args_name = checktrue(arg.second);
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
double CPLForceDrag::drag_coefficient(double D, std::vector<double> Ui_v, double eps) {
    //std::cout << "DRAG Cd: " << Cd << std::endl;
    return Cd; //Use default Drag value
}

double CPLForceDrag::get_eps(double r[]){

    double eps;
    eps = 1.0 - instant_volSums->get_array_value(r)/volSums->get_dV();

//    std::cout << "get eps " << volSums->get_array_value(r) << " "
//                     << instant_volSums->get_array_value(r) << " " << Npre_force << " " <<
//                       volSums->get_dV() << " " << eps << std::endl;
        
    if (eps < 1e-5) {
        #if DEBUG
            std::vector<int> cell = get_cell(r);
            std::cout << "Warning: 0 particles in cell (" 
                      << cell[0] << ", " << cell[1] << ", " << cell[2] << ")"
                      << " eps = " << eps << " volSums = " << volSums->get_array_value(r)
                      << " dV " << volSums->get_dV() << " Npreforce = " << Npre_force << std::endl;
        #endif
        return 0.0;
    } else {
        return eps;
    }
}

//Pre force collection of sums including overlap code to assign volumes
void CPLForceDrag::pre_force(double r[], double v[], double a[], 
                             double m,   double s,   double e) {

//    throw std::runtime_error("CPLForceDrag::pre_force no porosity/velocity needed so don't need to call this");
    double volume = (4./3.)*M_PI*pow(s,3); 
    double v_vol[]= {v[0]*volume, v[1]*volume, v[2]*volume};
    nSums->add_to_array(r, 1.0);
#if DEBUG
    std::cout << "Pre_force " << use_overlap << " " 
              << r[0] << " " << r[1] << " " << r[2] << " " <<volume << std::endl;
#endif
    if (! use_overlap){
        instant_volSums->add_to_array(r, volume);
        volSums->add_to_array(r, volume);   
        vSums->add_to_array(r, v_vol);
    } else {
        instant_volSums->add_to_array(r, s, volume);
        volSums->add_to_array(r, s, volume);
        vSums->add_to_array(r, s, v_vol);
    }
}


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForceDrag::get_force(double r[], double v[], double a[], 
                                            double m,   double s,   double e)
{

    //Define variable
    std::vector<int> cell(3);
    std::vector<double> Avi(3), Ui_v(3), fi(3), Ui(3), gradP(3), divStress(3);

    //Get all elements of recieved field
    if (! use_interpolate){

        //Based on cell
        cell = cfd_array_field->get_cell(r);
        Ui[0] = cfd_array_field->get_array_value(0, cell[0], cell[1], cell[2]);
        Ui[1] = cfd_array_field->get_array_value(1, cell[0], cell[1], cell[2]);
        Ui[2] = cfd_array_field->get_array_value(2, cell[0], cell[1], cell[2]);
        gradP[0] = cfd_array_field->get_array_value(3, cell[0], cell[1], cell[2]);
        gradP[1] = cfd_array_field->get_array_value(4, cell[0], cell[1], cell[2]);
        gradP[2] = cfd_array_field->get_array_value(5, cell[0], cell[1], cell[2]);
        divStress[0] = cfd_array_field->get_array_value(6, cell[0], cell[1], cell[2]);
        divStress[1] = cfd_array_field->get_array_value(7, cell[0], cell[1], cell[2]);
        divStress[2] = cfd_array_field->get_array_value(8, cell[0], cell[1], cell[2]);

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
    Ui_v[0] = Ui[0]-v[0];
    Ui_v[1] = Ui[1]-v[1];
    Ui_v[2] = Ui[2]-v[2];

    //Get Diameter, drag coefficient and volume
    double D = 2.0*s;
    double volume = (4./3.)*M_PI*pow(s,3);

    // Get drag coefficients from various drag models. These are function of
    // the instantaneous porosity.
    double eps = CPLForceDrag::get_eps(r);
    double A = drag_coefficient(D, Ui_v, eps);

    //Just drag force here.
    fi[0] = eps*A*Ui_v[0];
    fi[1] = eps*A*Ui_v[1];
    fi[2] = eps*A*Ui_v[2];

    //Include pressure
    if (use_gradP)
    {
        fi[0] += -volume*gradP[0];
        fi[1] += -volume*gradP[1];
        fi[2] += -volume*gradP[2];
    }
    // and stress
    if (use_divStress)
    {
        fi[0] += volume*divStress[0];
        fi[1] += volume*divStress[1];
        fi[2] += volume*divStress[2];
    }


#if DEBUG
    std::cout << "CPLForceDrag::get_force "  
              <<  cell[0] << " " << cell[1] << " " << cell[2]
              << " " <<  volume << " " <<  A << " " 
              << v[1] << " " << Ui[1] << " " << fi[1] << std::endl;
#endif
    
    //We split A*(vi - u_CFD) into implicit and explicit part following  Heng
    //Xiao and Jin Sun (2011) Commun. Comput. Phys. Vol. 9, No. 2, pp. 297-323
    //Define Avi=A*v[i] which is explicit part of the force based on molecular
    //velocity and ForceCoeff is passed so implict A*u_CFD can be applied in
    //SediFOAM. These have a eps multiplier incorporated as this is the
    //instantaneous porosity, and we need to collect/pass the produce as
    //<eps(t)>*<A(t)*vi(t) > not equal to <eps(t)*A(t)*vi(t)>. In
    //CPLSocketFOAM, this is divided through by the time averages porosity.
    Avi[0] = eps*A*v[0];
    Avi[1] = eps*A*v[1];
    Avi[2] = eps*A*v[2];
    A = eps*A;
    if (! use_interpolate){
        FcoeffSums->add_to_array(0, cell[0], cell[1], cell[2], A);
        FSums->add_to_array(0, cell[0], cell[1], cell[2], Avi[0]);
        FSums->add_to_array(1, cell[0], cell[1], cell[2], Avi[1]);
        FSums->add_to_array(2, cell[0], cell[1], cell[2], Avi[2]);
    } else {
        FcoeffSums->add_to_array(r, A);
        FSums->add_to_array(r, Avi.data());
        // FSums->add_to_array(r, Avi);
    }

//    std::cout << "Drag Force "  << A << " " 
//              << r[2] << " " << v[0] << " " << Ui[0] << " "  << v[1] << " " << Ui[1] << " " << v[2] << " " << Ui[2] << " " 
//              << gradP[0]  << " " << gradP[1] << " " << gradP[2]  << " "
//              << divStress[0] << " " << divStress[1] << " " << divStress[2] << " "
//              << " " << fi[0] <<" " << fi[1] <<" " << fi[2] << " " << fi.size() << std::endl;

    return fi;
}


//Post force collection of sums including overlap code to assign volumes
void CPLForceDrag::post_force(double r[], double v[], double a[], 
                             double m, double s, double e) {


    throw std::runtime_error("CPLForceDrag::post_force not needed, use value from pre force");
//    double volume = (4./3.)*M_PI*pow(s,3); 
//    double v_vol[]= {v[0]*volume, v[1]*volume, v[2]*volume};
//    nSums->add_to_array(r, 1.0);

//    std::cout << "Post_force " << use_overlap << " " << r[0] << " " << r[1] << " " << r[2] 
//                                              << " " << v[0] << " " << v[1] << " " << v[2] 
//                                              << " " <<volume << std::endl;
//    if (! use_overlap){
//        volSums->add_to_array(r, volume);           
//        vSums->add_to_array(r, v_vol);
//    } else {
//        volSums->add_to_array(r, s, volume);
//        vSums->add_to_array(r, s, v_vol);
//    }
}

///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForceGranular                           //
//                                                               //
///////////////////////////////////////////////////////////////////
// General Class to use for Ergun (1952), Di Felice (1994), BVK...


void CPLForceGranular::initialisesums(CPL::ndArray<double> arrayin){

    //Default values
    calc_preforce = true;
    calc_postforce = false;
    CPLForceDrag::initialisesums(arrayin);

}

//Pre force collection of sums including overlap code to assign volumes
//void CPLForceGranular::pre_force(double r[], double v[], double a[], 
//                                 double m,   double s,   double e) {

//    double volume = (4./3.)*M_PI*pow(s,3); 
//    double v_vol[]= {v[0]*volume, v[1]*volume, v[2]*volume};
//    nSums->add_to_array(r, 1.0);
//#if DEBUG
//    std::cout << "Pre_force " << use_overlap << " " << 
//         r[0] << " " << r[1] << " " << r[2] << " " <<volume << std::endl;
//#endif
//    if (! use_overlap){
//        instant_volSums->add_to_array(r, volume);
//        volSums->add_to_array(r, volume);           
//        vSums->add_to_array(r, v_vol);
//    } else {
//        instant_volSums->add_to_array(r, s, volume);
//        volSums->add_to_array(r, s, volume);
//        vSums->add_to_array(r, s, v_vol);
//    }
//}


//Return Stokes force
double CPLForceGranular::Stokes(double D, double mu) {
    return 9.42477796076938 * mu * D;
}

// Reynolds Number
// It is unclear here if Reynolds No. should be based
// on the mean cell velocity or particle velocity (or mean/relative from both)
double CPLForceGranular::Reynolds_number(double D, double U, double rho,
                                         double mu, double eps) {
    return rho * D * eps * U / mu;
}

double CPLForceGranular::magnitude(std::vector<double> v){
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
} 



///////////////////////////////////////////////////////////////////
//                      Stokes                                  //
///////////////////////////////////////////////////////////////////
double CPLForceStokes::drag_coefficient(double D, std::vector<double> Ui_v, double eps) 
{
    double Bi = 3. * M_PI * mu * D * eps; 
    return Bi;
}

///////////////////////////////////////////////////////////////////
//                      Di_Felice                                //
///////////////////////////////////////////////////////////////////
double CPLForceDi_Felice::drag_coefficient(double D, std::vector<double> Ui_v, double eps) 
{

    double MIN_REL_VELOCITY = 1e-8;
    double U = CPLForceGranular::magnitude(Ui_v);
    if ((U < MIN_REL_VELOCITY) || (eps == 0.0))
    {
        return 0.0;
    } 
    else 
    {
        double Re = CPLForceGranular::Reynolds_number(D, U, rho, mu, eps);
        double xi = 3.7 - 0.65 * exp(-0.5*pow((1.5-log10(Re)),2));
        double Cdi = pow((0.63 + 4.8*pow(Re,-0.5)),2);
        double Bi = 0.125 * Cdi * rho * M_PI * pow(D,2) * pow(eps,2.-xi) * abs(U);
        return Bi;
    }
}

///////////////////////////////////////////////////////////////////
//                          Tang                                 //
///////////////////////////////////////////////////////////////////

double CPLForceTang::drag_coefficient(double D, std::vector<double> Ui_v, double eps) 
{

    double Re;
    double U = CPLForceGranular::magnitude(Ui_v);
    if (U > 1e-8) {
        Re = CPLForceGranular::Reynolds_number(D, U, rho, mu, 1.0);
    } else {
        return 0;
    }

    if (eps == 0.0) {
        return 0.0;
    } else {
        double phi = 1. - eps;
        double Tang = 10 * phi / pow(eps,2) 
                      + pow(eps,2) * (1 + 1.5 * pow(phi,0.5)) 
                      + Re * (  0.11 * phi * (1 + phi) 
                              - 0.00456 / pow(eps,4) 
                              + (0.169 * eps + 0.0644/ pow(eps,4))
                              * pow(Re,-0.343));
        return CPLForceGranular::Stokes(D, mu) * Tang;
    }

}


///////////////////////////////////////////////////////////////////
//                          Ergun                                //
///////////////////////////////////////////////////////////////////

double CPLForceErgun::drag_coefficient(double D, std::vector<double> Ui_v, double eps) 
{
//     double phi = 1 - eps;
//     double U = CPLForceGranular::magnitude(Ui_v);
//     if (eps == 0.0) {
//         return 0.0;
//     } else {
//         double Re = CPLForceGranular::Reynolds_number(D, U, rho, mu, eps);
//         double Ergun = 0.0555555555555555 * (  150. * (phi / pow(eps,1))
//                                              + 1.75 * ( Re / pow(eps,0)));

// //        std::cout << "Ergun: " << phi << " " << D  << " " 
// //                  << U << " " << Re << " "
// //                  <<  150. * (phi / pow(eps,2)) << " " 
// //                  << 1.75 * ( Re / pow(eps,2)) << " "
// //                  << Ergun << " "
// //                  <<  CPLForceGranular::Stokes(D, mu) << " "
// //                  << Ergun * CPLForceGranular::Stokes(D, mu) << std::endl;

//         //Form matching Chris Knights code (and Gupta 2015 thesis, normalised by stokes)
//         return Ergun * CPLForceGranular::Stokes(D, mu);

//         //OpenFOAM Form
//         //return 150.0*    mu/(pow(phi*D, 2.0)) + 1.75*rho*Ur/(phi*D);

//         //Form from paper by ...
//         //return 150.0*eps*mu/(pow(eps*D, 2.0)) + 1.75*rho/(eps*D);


//     }

    double MIN_REL_VELOCITY = 1e-8;
    double U = CPLForceGranular::magnitude(Ui_v);
    if ((U < MIN_REL_VELOCITY) || (eps == 0.0))
    {
        return 0.0;
    } 
    else 
    {
        double Bi = (150. * M_PI * mu * D / 6.) * ((1. - eps) / eps) +
                    (1.75 * M_PI * rho * pow(D,2) / 6.) * abs(U);
        return Bi;
    }
}


///////////////////////////////////////////////////////////////////
//                          BVK                                  //
///////////////////////////////////////////////////////////////////

//Calculate BVK drag force per particle
double CPLForceBVK::drag_coefficient(double D, std::vector<double> Ui_v, double eps) 
{

    double Re;
    double U = CPLForceGranular::magnitude(Ui_v);
    if (U > 1e-8) {
        Re = CPLForceGranular::Reynolds_number(D, U, rho, mu, 1.0);
    } else {
        return 0;
    }

    if (eps == 0.0) {
        return 0.0;
    } else {
        double phi = 1 - eps;
        double BVK = 10.0*phi/pow(eps,2)
              + pow(eps,2)*(1.0+1.5*pow(phi,0.5)) 
              + (0.413*Re/(24*pow(eps,2))) 
              * (((1./eps) + 3 * phi * eps + 8.4 * pow(Re,-0.343)) 
                 /(1 + pow(10,(3*phi)) * pow(Re,-0.5) *(1+4*phi)));
    //    std::cout  << "BVK: " << Stokes(D, U, mu) * BVK << std::endl;
        return CPLForceGranular::Stokes(D, mu) * BVK;
    }

}


///////////////////////////////////////////////////////////////////
//                          BVK Poly                             //
///////////////////////////////////////////////////////////////////


//void CPLForceBVK_poly::initialise_extrasums(CPL::ndArray<double> arrayin) 
//{
//    int i = arrayin.shape(1);
//    int j = arrayin.shape(2);
//    int k = arrayin.shape(3);

//    //Set by input arguments
//    int histbins = 5;
//    double Dmin = 0.0;
//    double Dmax = 10.0;
//    double binwidth = (Dmax - Dmin)/float(histbins);
//    Dhist = std::make_shared<CPL::CPLField>(histbins, i, j, k, "Dhist");
//    Dhist->zero_array();
//    fields.push_back(Dhist);

//}

////Can be inhereted from Base Class
//void CPLForceBVK_poly::resetsums(){

//    //Reset all sums as normal
//    CPLForceDrag::resetsums();
//}


////Pre force collection of sums including overlap code to assign volumes
//int CPLForceBVK_poly::get_hist_bin(double value){
//    return int(std::floor(D/binwidth));
//}


////Pre force collection of sums including overlap code to assign volumes
//void CPLForceBVK_poly::pre_force(double r[], double v[], double a[], 
//                                 double m, double s, double e) 
//{

//    //Set range of diameters based on min/max D
//    if (first_time){
//        Dmin = std::min(Dmin, 2.*s);
//        Dmax = std::max(Dmax, 2.*s);
//        binwidth = (Dmax - Dmin)/float(histbins);
//        Dmean += 2.*s;
//        Nsum += 1
//    }

//    //Call existing pre-force
//    CPLForceDrag::pre_force(r, v, a, m, s, e);

//    //Add to histogram per bin
//    int bin = get_hist_bin(double value);
//    if (! use_overlap){
//        DSums->add_to_array(bin, r, D);
//    } else {
//        DSums->add_to_array(bin, r, s, D);
//    }

//}

////Calculate BVK drag force per particle
//double CPLForceBVK_poly::drag_coefficient(double r[], double D, 
//                                          std::vector<double> Ui_v) 
//{

//    //Pre force has setup everything, turn of first time
//    if (first_time){
//        Dmean = Dmean/Nsum;
//        mean_drag = CPLForceBVK::drag_coefficient(r, Dmean, Ui_v);
//        first_time = false;
//    }

//    double Re, y1;
//    double U = CPLForceGranular::magnitude(Ui_v);
//    if (U > 1e-8) {
//        Re = CPLForceGranular::Reynolds_number(D, U, rho, mu, 1.0);
//    } else {
//        return 0;
//    }

//    double eps = CPLForceDrag::get_eps(r);
//    if (eps == 0.0) {
//        return 0.0;
//    } else {
//        double phi = 1 - eps;
//        if (sumcells){
//            Dhist = sum(DSums,(1,2,3))
//            y1 = Dhist / Dmean;
//        } else {
//            int bin = get_hist_bin(double value);
//            y1 = DSums(bin, cell[0], cell[1], cell[2]) / sum(DSums,0);
//        }
//        Fp = mean_drag * (eps * y1 + phi * y1**2 + 0.064 * eps * y1**3);
//        return np.dot(Fp, bin_sizes);

//    }

//}


//Example print statement useful to copy in...
//    std::cout << "FLEKKOY: " << cell[0] << " " << cell[1]  << " " << cell[2]  
//                << nSums(icell, jcell, kcell) << " " 
//                << min[1] << " " <<  max[1] << " " << r[1]<< std::endl;


