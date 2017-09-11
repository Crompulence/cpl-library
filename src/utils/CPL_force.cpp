#include <vector>
#include <math.h> 
#include <stdexcept>

#include "CPL_ndArray.h"
#include "CPL_force.h"
#include "CPL_field.h"


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
CPLForce::CPLForce(int nd, int icells, int jcells, int kcells){
    //fieldptr = std::make_shared<CPL::CPLField>(nd, icells, jcells, kcells);
    fieldptr = new CPL::CPLField(nd, icells, jcells, kcells);

};

CPLForce::CPLForce(CPL::ndArray<double> arrayin){
    //fieldptr = std::make_shared<CPL::CPLField>(arrayin);
    fieldptr = new CPL::CPLField(arrayin);
};


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


void CPLForceDrag::initialisesums(CPL::ndArray<double> arrayin){

    int nsumsShape[3] = {arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    nSums.resize(3, nsumsShape); // Sum of number of particles

    int esumsShape[3] = {arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    eSums.resize(3, esumsShape); // Sum of porosity of particles  

    int FsumsShape[4] = {3, arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    FSums.resize(4, FsumsShape); // Sum of force on particles  
    resetsums();
}

void CPLForceDrag::resetsums(){
    nSums = 0.0; eSums = 0.0; FSums=0.0;
}

//Arbitary constant
double CPLForceDrag::drag_coefficient() {
    return 0.00001;
}

//Pre force collection of sums (should this come from LAMMPS fix chunk/atom bin/3d)
void CPLForceDrag::pre_force(double r[], double v[], double a[], double m, double s, double e) {

    // Should use field.add_volume(r, radius);
    double radius = s;
    double volume = (4./3.)*M_PI*pow(radius,3); 
    std::vector<int> cell = get_cell(r);
    nSums(cell[0], cell[1], cell[2]) += 1; 
    eSums(cell[0], cell[1], cell[2]) += volume; 

}

//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForceDrag::get_force(double r[], double v[], double a[], double m, double s, double e){

    std::vector<double> f(3), Ui(3), Ui_v(3);
    std::vector<int> cell = get_cell(r);
    CPL::ndArray<double> array = fieldptr->get_array();

    //Should use std::vector<double> Ui(3) = field.interpolate(r);
    for (int i=0; i<3; i++){
        Ui[i] = array(i, cell[0], cell[1], cell[2]);
        Ui_v[i] = Ui[i]-v[i];
    }

    double Cd = drag_coefficient();

    //Calculate force
    for (int i = 0; i < 3; ++i){
        f[i] = Cd*Ui_v[i];
        FSums(i, cell[0], cell[1], cell[2]) += f[i];
    }

//    std::cout << "Drag Force "  
//              << r[2] << " " << v[2] << " "
//              << f[2] << " " << Ui[2] << std::endl;

    return f;
}

///////////////////////////////////////////////////////////////////
//                                                               //
//                    CPLForceGranular                           //
//                                                               //
///////////////////////////////////////////////////////////////////



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


//Pre force collection of sums (can this come from LAMMPS fix chunk/atom bin/3d)
std::vector<double> CPLForceGranular::get_force(double r[], double v[], double a[], double m, double s, double e) {

    std::vector<double> f(3), Ui(3), Ui_v(3);
    std::vector<int> cell = get_cell(r);
    CPL::ndArray<double> array = fieldptr->get_array();

    double radius = s;
    double volume = (4./3.)*M_PI*pow(radius,3); 

    //Porosity e is cell volume - sum in volume
    double cellvolume = fieldptr->dV;
    double eps = 1.0 - eSums(cell[0], cell[1], cell[2])/cellvolume;
    double rho = 1.0; //m/volume;
    double mu = 1.0;
    double d = 2.0*radius;
    //Should use std::vector<double> Ui(3) = field.interpolate(r);
    for (int i=0; i<3; i++){
        Ui[i] = array(i, cell[0], cell[1], cell[2]);
        Ui_v[i] = Ui[i]-v[i];
    }

    //It is unclear here if Reynolds No. should be based
    //on the mean cell velocity or particle velocity
    double Re = rho * d * eps * magnitude(Ui_v) / mu;
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



//Example print statement useful to copy in...
//    std::cout << "FLEKKOY: " << cell[0] << " " << cell[1]  << " " << cell[2]  
//                << nSums(icell, jcell, kcell) << " " 
//                << min[1] << " " <<  max[1] << " " << r[1]<< std::endl;


/*

//=====================================================
//The post force fix in LAMMPS would then become
//=====================================================

void FixCPLForce::post_force (int vflag) {

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double dx = CPL::get<double> ("dx");     
    double dy = CPL::get<double> ("dy");     
    double dz = CPL::get<double> ("dz");     
    double dA = dx*dz;
    double dV = dx*dy*dz;
    
    char* groupstr = "cplforcegroup";
    char* regionstr = "cplforceregion";

    int cplforcegroup = group->find (groupstr);
    int groupbit = group->bitmask[cplforcegroup];

    int rid = domain->find_region (regionstr);
    auto cplforceregion = domain->regions[rid];

    //Constructor here for forces based on input
    if ("FORCE_TYPE" == "Flekkoy"){
        fxyz = CPLForceFlekkoy->initialise(field)
    } else ("FORCE_TYPE" == "Nieetal"){
        fxyz = CPLNieetal->initialise(field)
    }

    // First loop to sum up all quantities and prepare to apply force
    for (int i = 0; i < nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {
            // Find in which cell number (local to processor) is the particle
            // and sum all the Flekkøy weights for each cell.
            int glob_cell[3];
            CPL::map_coord2cell(x[i][0], x[i][1], x[i][2], glob_cell);

            int loc_cell[3];
            bool validCell = CPL::map_glob2loc_cell(procPortion.data(), glob_cell, loc_cell);

            //if (validCell)
            //std::cout << "local: " << icell << " " << jcell << " " << kcell << " global: " << glob_cell[0] << " " << glob_cell[1] << " " << glob_cell[2] \
            //<< "portion:" << procPortion[0] << " "<< procPortion[1] << " "<< procPortion[2] << " "<< procPortion[3] << " "<< procPortion[4] << " "<< procPortion[5] << " " << std::endl;

            if (! validCell) {
               // std::cout << "Warning: an atom in the constrained region is within an invalid cell. \n"
               //           << "This should never happen and it is likely a BUG. Report." << std::endl;
                continue;
            }

            int icell = loc_cell[0];
            int jcell = loc_cell[1];
            int kcell = loc_cell[2];

            //This is call to set pre-force 
            fxyz->pre_force(icell, jcell, kcell)

        }
    }


    // Calculate force and apply
    for (int i = 0; i < nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {
            int glob_cell[3];
            CPL::map_coord2cell(x[i][0], x[i][1], x[i][2], glob_cell);

            int loc_cell[3];
            bool validCell = CPL::map_glob2loc_cell(procPortion.data(), glob_cell, loc_cell);


            if (! validCell) {
//                std::cout << "Warning: an atom in the constrained region is within an invalid cell. \n"
 //                         << "This should never happen and it is likely a BUG. Report." << std::endl;
                continue;
            }

            int icell = loc_cell[0];
            int jcell = loc_cell[1];
            int kcell = loc_cell[2];

            // Get forces
            [fx, fy, fz] = fxyz->get_force(i, icell, jcell, kcell, field)

            f[i][0] += fx;
            f[i][1] += fy;
            f[i][2] += fz;

        }
    }
}
*/
