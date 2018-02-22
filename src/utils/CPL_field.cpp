#include <vector>
#include <math.h> 
#include <stdexcept>

#include "CPL_ndArray.h"
#include "CPL_field.h"
#include "overlap/overlap.hpp"
#include "interpolate/interp_ho.h"

using namespace CPL;


//Constructor based on specified size
CPLField::CPLField(int nd, int icells, int jcells, int kcells){
    // Fields
    int arrayShape[4] = {nd, icells, jcells, kcells};
    array.resize(4, arrayShape);
    for (int i = 0; i < 3; ++i){
        min[i] = 0.0;
        max[i] = 1.0;
    }
    set_dxyz();
}

//Constructor which uses arrayin and sets size and initial value to this
CPLField::CPLField(CPL::ndArray<double> arrayin){
    array = arrayin;
    for (int i = 0; i < 3; ++i){
        min[i] = 0.0;
        max[i] = 1.0;
    }
    set_dxyz();
}


//Constructor which uses size of array in x, y and z
CPLField::CPLField(int nd, CPL::ndArray<double> arrayin){

    int arrayShape[4] = {nd, arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    array.resize(4, arrayShape);

    for (int i = 0; i < 3; ++i){
        min[i] = 0.0;
        max[i] = 1.0;
    }
    set_dxyz();
}


///////////////////////////////////////////////////
/*                                                /
// .----..----..---.  .---. .----..----.  .----.  /
/  { {__  | {_ {_   _}{_   _}| {_  | {}  }{ {__   /
/  .-._} }| {__  | |    | |  | {__ | .-. \.-._} } /
/  `----' `----' `-'    `-'  `----'`-' `-'`----'  /
/                                                 /
*//////////////////////////////////////////////////


//Set minimum and maximum values of field application
void CPLField::set_minmax(double min_in[], double max_in[]){
    for (int i = 0; i < 3; ++i){
        min[i] = min_in[i];
        max[i] = max_in[i];
    }
    set_dxyz();
}

//If either min/max change or field object, we need to recalculate dx, dy and dz
void CPLField::set_dxyz(){
    for (int i = 0; i < 3; ++i){
        dxyz[i] = (max[i] - min[i])/array.shape(i+1);
    }

    dA[0] = dxyz[1]*dxyz[2];
    dA[1] = dxyz[0]*dxyz[2];
    dA[2] = dxyz[0]*dxyz[1];
    dV = dxyz[0]*dxyz[1]*dxyz[2];
}


void CPLField::set_array(CPL::ndArray<double> arrayin){
    array = arrayin;
    set_dxyz();
}

void CPLField::zero_array(){
    array = 0.0;
}

// A function which uses particle position to determine
// which cell it is in and adds the appropriate value
// to that cell, e.g. 1 for Nsums, velocity for vsum, etc

//Add value to a particular cell
void CPLField::add_to_array(int n, int i, int j, int k, double value){
    array(n, i, j, k) += value;
}

//Just add to cell based on where centre of particle falls
//assuming that array of values is the same size as array
void CPLField::add_to_array(const double r[], const double value[]){
    std::vector<int> cell = get_cell(r);
    for (int n=0; n < array.shape(0); n++)
        add_to_array(n, cell[0], cell[1], cell[2], value[n]);
}

//If just one value, then doesn't need to be a reference to array
void CPLField::add_to_array(const double r[], const double value_){
    std::vector<int> cell = get_cell(r);
    add_to_array(0, cell[0], cell[1], cell[2], value_);
}
   

// If the radius of the particle is included, add a fraction
//Add fraction of value to cell based on fraction of sphere inside
void CPLField::add_to_array(const double r[], double s, const double value[]){

    int ip, jp, kp;
    double volume = (4./3.)*M_PI*pow(s,3); 
    std::vector<int> cell = get_cell(r);
    double box[6];

    //Get fraction of sphere in a volume
    double dx = dxyz[0];
    double dy = dxyz[1];
    double dz = dxyz[2];
    int nxps = ceil(s/dx);
    int nyps = ceil(s/dy);
    int nzps = ceil(s/dz);
    int i = cell[0]; int j = cell[1]; int k = cell[2];

//    std::cout << "overlap calc "  
//              << dx << " " << dy << " " << dz << " " 
//              << i << " " << j << " " << k << " " 
//              << nxps << " " << nyps << " " << nzps << std::endl;

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
        double Vsphereinbox = sphere_cube_overlap(r[0], r[1], r[2], s, 
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
        if (ip > array.shape(1)-1) ip = array.shape(1)-1;
        if (jp > array.shape(2)-1) jp = array.shape(2)-1;
        if (kp > array.shape(3)-1) kp = array.shape(3)-1;

        double frac = Vsphereinbox/volume;
        for (int n=0; n < array.shape(0); n++)
            array(n, ip, jp, kp) += frac*value[n]; 

//        if (Vsphereinbox > 1e-12) {
//            std::cout << "overlap cells "  
//                  << i << " " << j << " " << k << " " 
//                  << ip << " " << jp << " " << kp << " "
//                  << i+ic << " " << j+jc << " " << k+kc << " "
//                  << r[0] << " " << r[1] << " " << r[2] << " "
//                  << box[0] << " " << box[1] << " " << box[2] << " " 
//                  << box[3] << " " << box[4] << " " << box[5] <<  " "
//                  << Vsphereinbox << " " << array.shape(0) 
//                  << " " << value[0] << " " << frac << " " 
//                  << array(0, ip, jp, kp) << std::endl;
//        }

    }}}
}


//If just one value, then doesn't need to be a reference to array
void CPLField::add_to_array(const double r[], double s, const double value_){
    double value[1] = {value_};
    add_to_array(r, s, value);
}




// Use the block of 26 cells around a cell specified by "ic, jc, kc" to get the nodes of that cell
CPL::ndArray<double> CPLField::celltonode(CPL::ndArray<double> cell, int n, int ic, int jc, int kc)
{

    CPL::ndArray<double> node;
    int arrayShape[4] = {1, 2, 2, 2};
    node.resize(4, arrayShape);

//        std::cout << "cell " <<  ic << " " << jc << " " << kc << " " << cell(n, ic  , jc  , kc  ) << "\n";

    for (int i = 0; i < 2; i++ ){
    for (int j = 0; j < 2; j++ ){
    for (int k = 0; k < 2; k++ ){
        node(0,i,j,k) = 0.125*(   cell(n, ic+i  , jc+j  , kc+k  )
                                + cell(n, ic+i+1, jc+j  , kc+k  )
                                + cell(n, ic+i  , jc+j+1, kc+k  )
                                + cell(n, ic+i  , jc+j  , kc+k+1)
                                + cell(n ,ic+i  , jc+j+1, kc+k+1)
                                + cell(n, ic+i+1, jc+j  , kc+k+1)
                                + cell(n, ic+i+1, jc+j+1, kc+k  )
                                + cell(n, ic+i+1, jc+j+1, kc+k+1));

//        std::cout << "NODES " <<  i << " " << j << " " << k << " " << cell(n, ic+i  , jc+j  , kc+k  ) << " " << node(0,i,j,k) << "\n";
    }}}

    return node;

}


double CPLField::sphere_cube_overlap(double scx, double scy, double scz, double sr,
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
    

///////////////////////////////////////////////////
/*                                                /
/  .---. .----..---.  .---. .----..----.  .----.  /
/ /   __}| {_ {_   _}{_   _}| {_  | {}  }{ {__    /
/ \  {_ }| {__  | |    | |  | {__ | .-. \.-._} }  /
/  `---' `----' `-'    `-'  `----'`-' `-'`----'   /
/                                                 /
*//////////////////////////////////////////////////

//Get dA
std::vector<double> CPLField::get_dA(){

    std::vector<double> retrnd_dA(3);

    retrnd_dA[0] = dA[0];
    retrnd_dA[1] = dA[1];
    retrnd_dA[2] = dA[2];

    return retrnd_dA;
}


//Get cell from min/max and dx
std::vector<int> CPLField::get_cell(const double r[]){
    std::vector<int> cell(3);
    for (int i = 0; i < 3; ++i){
        cell[i] = floor((r[i]-min[i])/dxyz[i]);
        //Check cell is within the domain
        if (cell[i] > floor((max[i]-min[i])/dxyz[i]))
            cell[i] = floor((max[i]-min[i])/dxyz[i]);
            //throw std::domain_error("get_cell Error: Input above domain");

        if (cell[i] < 0)
            cell[i] = 0;
            //throw std::domain_error("get_cell Error: Input below domain");
    }
    return cell;
}


//This is terrible for efficiency, it creates a copy of the array!
CPL::ndArray<double> CPLField::get_array(){
    return array;
}

CPL::ndArray<double>& CPLField::get_array_pointer(){
    return array;
}

// A function which gets value n in cell i,j,k
std::vector<double> CPLField::get_array_value(const std::vector<int> indices, int i, int j, int k){
    std::vector<double> v;
    for ( auto &n : indices ) { 
        v.push_back(array(n, i, j, k));
    }
    return v;
}

// wrapper for single index
double CPLField::get_array_value(const int index, int i, int j, int k){
    return array(index, i, j, k);
}

// A function which gets value using position of molecule, no interpolation
std::vector<double> CPLField::get_array_value(const std::vector<int> indices, const double r[]){
    std::vector<int> cell = get_cell(r);
    std::vector<double> v;
    for ( auto &n : indices ) { 
        v.push_back(array(n, cell[0], cell[1], cell[2]));
    }
    return v;
}


// A function which gets interpolated value in array cell using position of molecule
std::vector<double> CPLField::get_array_value_interp(const std::vector<int> indices, const double r[]){
    int order = 2;
    std::vector<double> out = interpolate(r, array, indices, order);
    return out;
}

//Assume 3D
std::vector<double> CPLField::interpolate(const double r[], 
                                          CPL::ndArray<double> cell_array, 
                                          const std::vector<int> indices, int order)
{
    int i, m, nd;
    int *n_1d;
    double *a;
    double *b;
    double *zd;
    double *zi;

    //Setup cell based on particle
    auto cell = get_cell(r);
    int ic = cell[0];
    int jc = cell[1];
    int kc = cell[2];

    //Assume 3D and set order in each direction
    m = 3;
    n_1d = new int[m];
    for ( i = 0; i < m; i++ )
        n_1d[i] = order;

    //If we always set to 0 and 1 then no scaling applied
    a = new double[m];
    b = new double[m];
    a[0] = (ic  )*dxyz[0];
    a[1] = (jc  )*dxyz[1];
    a[2] = (kc  )*dxyz[2];
    b[0] = (ic+1)*dxyz[0];
    b[1] = (jc+1)*dxyz[1];
    b[2] = (kc+1)*dxyz[2];

    //Setup grid and get nodes
    assert(cell_array.shape(1) > ic+2);
    assert(cell_array.shape(2) > jc+2);
    assert(cell_array.shape(3) > kc+2);

    std::vector<double> v;
    for ( auto &n : indices ) {
        //We should loop over many n here and return interp for all elements 
        CPL::ndArray<double> node = celltonode(cell_array, n, ic, jc, kc);
        zd = node.data();

        double rcopy[3] = {r[0], r[1], r[2]};
        nd = lagrange_interp_nd_size(m, n_1d);
        zi = lagrange_interp_nd_value(m, n_1d, a, b, nd, zd, 1, rcopy);
        v.push_back(zi[0]);  //Push each new element of vec
    }

    delete a;
    delete b;
    delete n_1d;
    delete zi;

    // Vector of all solutions
    return v;


}



