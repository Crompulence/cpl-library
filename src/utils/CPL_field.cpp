#include <vector>
#include <math.h> 
#include <stdexcept>

#include "CPL_ndArray.h"
#include "CPL_field.h"
#include "overlap/overlap.hpp"

//#include "interpolate/lagrange_interp_nd.hpp"
int lagrange_interp_nd_size ( int, int* );

using namespace CPL;

//Constructors
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

//Should this be a std::shared_ptr <CPL::ndArray <double>> arrayin??
//Surely a unique pointer is better is we have to use pointers at all
CPLField::CPLField(CPL::ndArray<double> arrayin){
    array = arrayin;
    for (int i = 0; i < 3; ++i){
        min[i] = 0.0;
        max[i] = 1.0;
    }
    set_dxyz();
}

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


//Get dA
std::vector<double> CPLField::get_dA(){

    std::vector<double> retrnd_dA(3);

    retrnd_dA[0] = dA[0];
    retrnd_dA[1] = dA[1];
    retrnd_dA[2] = dA[2];

    return retrnd_dA;
}


//Get cell from min/max and dx
std::vector<int> CPLField::get_cell(double r[]){
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



void CPLField::set_array(CPL::ndArray<double> arrayin){
    array = arrayin;
    set_dxyz();
}

//This is terrible for efficiency, it creates a copy of the array!
CPL::ndArray<double> CPLField::get_array(){
    return array;
}

CPL::ndArray<double>& CPLField::get_array_pointer(){
    return array;
}


//std::unique_ptr<CPL::ndArray<double>> CPLField::get_array_pointer(){
//    return std::unique_ptr<CPL::ndArray<double>> array;
//}


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
CPL::ndArray<double> CPLField::celltonode(CPL::ndArray<double> cell, int n, int ic, int jc, int kc)
{

    CPL::ndArray<double> node;
    int arrayShape[4] = {1, 2, 2, 2};
    node.resize(4, arrayShape);

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

        //cout << "NODES " <<  i << " " << j << " " << k << " " << node(0,i,j,k) << "\n";
    }}}

    return node;

}



//Assume 3D
std::vector<double> CPLField::interpolate(double r[], CPL::ndArray<double> cell_array, int n, int order)
{
    double *a;
    double *b;
    int i;
    int j;
    int m;
    int *n_1d;
    int nd;
    int ni;
    int seed;
    double *xd;
    double *xi;
    double *zd;
    double *ze;
    double *zi;

    //Assume 3D and set order to input
    m = 3;
    n_1d = new int[m];
    for ( i = 0; i < m; i++ )
    {
        n_1d[i] = order;
    }

    auto cell = get_cell(r);
    int ic = cell[0];
    int jc = cell[1];
    int kc = cell[2];

    CPL::ndArray<double> node = celltonode(cell_array, n, ic, jc, kc);
    zd = node.data();

    //Setup grid and size
    //nd = lagrange_interp_nd_size(m, n_1d);
    //xd = lagrange_interp_nd_grid(m, n_1d, a, b, nd);
    //zi = lagrange_interp_nd_value(m, n_1d, a, b, nd, zd, 1, r);

    //std::vector<double> vec = {zi[0], zi[1], zi[2]};
    std::vector<double> vec = {0., 0., 0.};
    return vec;

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



