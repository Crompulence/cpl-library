#include <vector>
#include <math.h> 
#include <stdexcept>

#include "CPL_ndArray.h"
#include "CPL_field.h"
#include "cpl.h"

//#include "interpolate/lagrange_interp_nd.hpp"
int lagrange_interp_nd_size ( int, int* );

using namespace CPL;

Domain::Domain(const std::valarray<double>& domain_orig,
               const std::valarray<double>& domain_length) : bounds(6){
    bounds[0] = domain_orig[0];
    bounds[1] = domain_orig[0] + domain_length[0];
    bounds[2] = domain_orig[1];
    bounds[3] = domain_orig[1] + domain_length[1];
    bounds[4] = domain_orig[2];
    bounds[5] = domain_orig[2] + domain_length[2];
    lx = domain_length[0];
    ly = domain_length[1];
    lz = domain_length[2];
    computeAV_();
}

Domain::Domain(const std::valarray<double>& domain_bounds) {
    lx = domain_bounds[0] - domain_bounds[1];
    ly = domain_bounds[2] - domain_bounds[3];
    lz = domain_bounds[4] - domain_bounds[5];
    bounds = domain_bounds;
    computeAV_();
}

void Domain::computeAV_() {
    Axy = lx * ly;
    Axz = lx * lz;
    Ayz = ly * lz;
    V = lx * ly * lz;
}


std::valarray<double> Domain::length() {
    return std::valarray<double>({lx, ly, lz});
}

std::valarray<double> Domain::area() {
    return std::valarray<double>({Axy, Axz, Ayz});
}


Field::Field(const CPL::Domain& domain,
             const std::vector<int>& field_cells):
             Domain(domain), nCells(field_cells) {

    computedAdV_();
}

//TODO: Add const to those when changing interface in cpllib
std::vector<int> PortionField::getLocalCell(std::vector<int>& glob_cell) {
    std::vector<int> loc_cell(3);
    CPL::map_glob2loc_cell(cellBounds.data(), glob_cell.data(), 
                           loc_cell.data());
    return loc_cell;
    
}


std::valarray<double> Field::getCoord(std::vector<int>& cell) {
}

std::vector<int> Field::getCell(std::valarray<double>& coord) {
}


std::valarray<double> Field::cellLength() {
    return std::valarray<double>({dx, dy, dz});
}


std::valarray<double> Field::cellArea() {
    return std::valarray<double>({dAxy, dAxz, dAyz});
}

void Field::computedAdV_() {
    dx = (bounds[1] - bounds[0]) / nCells[0];//data->shape(1);
    dy = (bounds[3] - bounds[2]) / nCells[1];//data->shape(2);
    dz = (bounds[5] - bounds[4]) / nCells[2];//data->shape(3);
    dAxy = dx * dy;
    dAxz = dx * dz;
    dAyz = dy * dz;
    dV = dx * dy * dz;
}





// Field:: Field(const std::valarray<double>& domain_bounds,
//               CPL::ndArray<double>& field_data, bool copy) : Domain(domain_bounds){
//     if (copy) {
//         data = new ndArray<double>(field_data);
//         own_memmory = true;
//     }
//     else {
//         data = &field_data;
//         own_memmory = false;
//     }
// }
//
// Field::~Field() {
//     if (own_memmory)
//         delete data;
// }

// //Get cell from min/max and dx
// std::vector<int> Field::cell(double r[]){
//     std::vector<int> cell(3);
//     for (int i = 0; i < 3; ++i){
//         cell[i] = floor((r[i]-min[i])/dxyz[i]);
//         //Check cell is within the domain
//         if (cell[i] > floor((max[i]-min[i])/dxyz[i]))
//             cell[i] = floor((max[i]-min[i])/dxyz[i]);
//             //throw std::domain_error("get_cell Error: Input above domain");
//
//         if (cell[i] < 0)
//             cell[i] = 0;
//             //throw std::domain_error("get_cell Error: Input below domain");
//     }
//     return cell;
// }
//

//Get interpolate of field
//std::vector<double> CPLField::interpolate(double r[]){

//    int interp[3];
//    std::vector<int> cell = get_cell(r);
//    std::vector<double> val(array.shape(0));
//    std::vector<double> rc(3);
//    double f,g,h;


//    //Check if at edge of domain and use one sided if needed
//    for (int n = 0; n < 3; n++){
//        if (cell[n]+1>array.shape(n)){
//            interp[n] = -1;
//        } else if (cell[0]-1<0){  
//            interp[n] = 1;
//        } else {
//            interp[n] = 0;
//        }
//        // Get position relative to cell minimum
//        rc[n] = r[n] - cell[n]*dxyz[n];
//    }

//    //Interpolate based on both adjacent cells
//    for (int n = 0; n<array.shape(0); n++){
//        if (interp[0] == 0){
//            f =  (   array(n, cell[0]+1,cell[1],  cell[2]  )
//                  -2*array(n, cell[0],  cell[1],  cell[2]  )  
//                   + array(n, cell[0]-1,cell[1],  cell[2]  ))/2.*dxyz[0];
//        } else if (interp[0] == 1){
//            f =  (   array(n, cell[0]+1,cell[1],  cell[2]  )
//                    -array(n, cell[0],  cell[1],  cell[2]  ))/dxyz[0];
//        } else if (interp[0] == -1){
//            f =  (   array(n, cell[0],  cell[1],  cell[2]  )
//                    -array(n, cell[0]-1,cell[1],  cell[2]  ))/dxyz[0];
//        }

//        if (interp[1] == 0){
//            g =  (   array(n, cell[0],  cell[1]+1,cell[2]  )
//                  -2*array(n, cell[0],  cell[1],  cell[2]  )  
//                   + array(n, cell[0],  cell[1]-1,cell[2]  ))/2.*dxyz[1];
//        } else if (interp[1] == 1){
//            g =  (   array(n, cell[0],  cell[1]+1,cell[2]  )
//                  -  array(n, cell[0],  cell[1],  cell[2]  ))/dxyz[1];
//        } else if (interp[1] == -1){
//            g =  (   array(n, cell[0],  cell[1]  ,cell[2]  )
//                   - array(n, cell[0],  cell[1]-1,cell[2]  ))/dxyz[1];
//        }

//        if (interp[2] == 0){
//            h =  (   array(n, cell[0],  cell[1],  cell[2]+1)
//                  -2*array(n, cell[0],  cell[1],  cell[2]  ) 
//                   + array(n, cell[0],  cell[1],  cell[2]-1))/2.*dxyz[2];
//        } else if (interp[2] == 1){
//            h =  (   array(n, cell[0],  cell[1],  cell[2]+1)
//                   - array(n, cell[0],  cell[1],  cell[2]  ))/dxyz[2];
//        } else if (interp[2] == -1){
//            h =  (   array(n, cell[0],  cell[1],  cell[2]  )
//                  -2*array(n, cell[0],  cell[1],  cell[2]-1))/dxyz[2];
//        }

//        val[n] = f*rc[0] + g*rc[1] + h*rc[2];
//    }

//    return val;
//}


// Use the block of 26 cells around a specified cell by "ic, jc, kc" into the nodes on that cell
// CPL::ndArray<double> CPLField::celltonode(CPL::ndArray<double> cell, int n, int ic, int jc, int kc)
// {
//
//     CPL::ndArray<double> node;
//     int arrayShape[4] = {1, 2, 2, 2};
//     node.resize(4, arrayShape);
//
//     for (int i = 0; i < 2; i++ ){
//     for (int j = 0; j < 2; j++ ){
//     for (int k = 0; k < 2; k++ ){
//         node(0,i,j,k) = 0.125*(   cell(n, ic+i  , jc+j  , kc+k  )
//                                 + cell(n, ic+i+1, jc+j  , kc+k  )
//                                 + cell(n, ic+i  , jc+j+1, kc+k  )
//                                 + cell(n, ic+i  , jc+j  , kc+k+1)
//                                 + cell(n ,ic+i  , jc+j+1, kc+k+1)
//                                 + cell(n, ic+i+1, jc+j  , kc+k+1)
//                                 + cell(n, ic+i+1, jc+j+1, kc+k  )
//                                 + cell(n, ic+i+1, jc+j+1, kc+k+1));
//
//         //cout << "NODES " <<  i << " " << j << " " << k << " " << node(0,i,j,k) << "\n";
//     }}}
//
//     return node;
//
// }
//


//Assume 3D
// std::vector<double> CPLField::interpolate(double r[], CPL::ndArray<double> cell_array, int n, int order)
// {
//     double *a;
//     double *b;
//     int i;
//     int j;
//     int m;
//     int *n_1d;
//     int nd;
//     int ni;
//     int seed;
//     double *xd;
//     double *xi;
//     double *zd;
//     double *ze;
//     double *zi;
//
//     //Assume 3D and set order to input
//     m = 3;
//     n_1d = new int[m];
//     for ( i = 0; i < m; i++ )
//     {
//         n_1d[i] = order;
//     }
//
//     auto cell = get_cell(r);
//     int ic = cell[0];
//     int jc = cell[1];
//     int kc = cell[2];
//
//     CPL::ndArray<double> node = celltonode(cell_array, n, ic, jc, kc);
//     zd = node.data();
//
//     //Setup grid and size
//     //nd = lagrange_interp_nd_size(m, n_1d);
//     //xd = lagrange_interp_nd_grid(m, n_1d, a, b, nd);
//     //zi = lagrange_interp_nd_value(m, n_1d, a, b, nd, zd, 1, r);
//
//     //std::vector<double> vec = {zi[0], zi[1], zi[2]};
//     std::vector<double> vec = {0., 0., 0.};
//     return vec;
//
// }
//


