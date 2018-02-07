# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
#include "CPL_ndArray.h"

using namespace std;

# include "interp_ho.h"

CPL::ndArray<double> celltonode(CPL::ndArray<double> cell, int n, int ic, int jc, int kc)
{

    CPL::ndArray<double> node;
    int arrayShape[4] = {1, 2, 2, 2};
    node.resize(4, arrayShape);

    for (int i = 0; i < 2; i++ ){
    for (int j = 0; j < 2; j++ ){
    for (int k = 0; k < 2; k++ ){
//            ip = i+ic; jp = j+jc; kp = k+kc;
//            if (ip < 0) ip = 0;
//            if (jp < 0) jp = 0;
//            if (kp < 0) kp = 0;
//            if (ip >= cell.shape(0)-1) ip = cell.shape(0)-1;
//            if (jp >= cell.shape(1)-1) jp = cell.shape(1)-1;
//            if (kp >= cell.shape(2)-1) kp = cell.shape(2)-1;

        node(0,i,j,k) = 0.125*(   cell(n, ic+i  , jc+j  , kc+k  )
                                + cell(n, ic+i+1, jc+j  , kc+k  )
                                + cell(n, ic+i  , jc+j+1, kc+k  )
                                + cell(n, ic+i  , jc+j  , kc+k+1)
                                + cell(n ,ic+i  , jc+j+1, kc+k+1)
                                + cell(n, ic+i+1, jc+j  , kc+k+1)
                                + cell(n, ic+i+1, jc+j+1, kc+k  )
                                + cell(n, ic+i+1, jc+j+1, kc+k+1));

        cout << "NODES " <<  i << " " << j << " " << k << " " << node(0,i,j,k) << "\n";
    }}}

    return node;

}

//Assume 3D
std::vector<double> interpolate(double r[], int ic, int jc, int kc, CPL::ndArray<double> cell_array, int order)
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

    CPL::ndArray<double> node = celltonode(cell_array, 0, ic, jc, kc);
    zd = node.data();

    //Assume 3D and set order to input
    m = 3;
    n_1d = new int[m];
    for ( i = 0; i < m; i++ )
    {
        n_1d[i] = order;
    }

    //Setup grid and size
    nd = lagrange_interp_nd_size(m, n_1d);
    xd = lagrange_interp_nd_grid(m, n_1d, a, b, nd);
    zi = lagrange_interp_nd_value(m, n_1d, a, b, nd, zd, 1, r);

    std::vector<double> vec = {zi[0], zi[1], zi[2]};

    return vec;

}


//****************************************************************************80
//
void test03 ( )
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

    //Assume 3D and order 2
    m = 3;
    n_1d = new int[m];
    for ( i = 0; i < m; i++ )
    {
        n_1d[i] = 2;
    }

    //If we always set to 0 and 1 then no scaling applied
    a = new double[m];
    b = new double[m];
    for ( i = 0; i < m; i++ )
    {
        a[i] = 0.0;
        b[i] = 1.0;
    }

    //From block of 27 cells, get nodes
    CPL::ndArray<double> cell;
    int arrayShape[4] = {2, 9, 9, 9};
    cell.resize(4, arrayShape);
    int n = 0;
    //Assume slabs in x
    for (int j = 0; j < 3; j++ ){
    for (int k = 0; k < 3; k++ ){
        cell(n,0,j,k) = -0.5;
        cell(n,1,j,k) = 0.5;
        cell(n,2,j,k) = 1.5;
    }}

    //Assume slabs in y
    for (int i = 0; i < 3; i++ ){
    for (int k = 0; k < 3; k++ ){
        cell(n,i,0,k) = -0.5;
        cell(n,i,1,k) = 0.5;
        cell(n,i,2,k) = 1.5;
    }}

    //Assume slabs in z
    for (int i = 0; i < 3; i++ ){
    for (int j = 0; j < 3; j++ ){
        cell(n,i,j,0) = -0.5;
        cell(n,i,j,1) = 0.5;
        cell(n,i,j,2) = 1.5;
    }}

    //Setup nodes
    CPL::ndArray<double> node = celltonode(cell, n, 0, 0, 0);
    zd = node.data();

    //If we always set to 0 and 1 then no scaling applied
    a = new double[m];
    b = new double[m];
    for ( i = 0; i < m; i++ )
    {
        a[i] = 0.0;
        b[i] = 1.0;
    }

    //Setup grid and size
    nd = lagrange_interp_nd_size(m, n_1d);
    xd = lagrange_interp_nd_grid(m, n_1d, a, b, nd);
    for ( j = 0; j < nd; j++ )
        cout << "xd, zd = " << j << " " << xd[j*m] <<  
                                    " " << xd[j*m+1] << 
                                    " " << xd[j*m+2] << 
                                    " " << zd[j] << "\n";

    //Hardwire some values
    ni = 2;  
    xi = new double[m*ni];
    xi[0] = 0.8; xi[1] = 0.2; xi[2] = 0.1;
    xi[3] = 0.2; xi[4] = 0.9; xi[5] = 0.8;
    zi = lagrange_interp_nd_value(m, n_1d, a, b, nd, zd, ni, xi);
    for ( j = 0; j < ni; j++ )
        cout << "Out set = " << j << " " << xi[j*m] <<  
                                     " " << xi[j*m+1] << 
                                     " " << xi[j*m+2] << 
                                     " " << zi[j] << "\n";

    //Get random sample points and get values
    ni = 100;
    seed = 123456789;
    //Interpolate using random points
    xi = r8mat_uniform_01_new(m, ni, seed);
    zi = lagrange_interp_nd_value(m, n_1d, a, b, nd, zd, ni, xi);

    //Plot
    for ( j = 0; j < ni; j++ ) {
        cout << "Out rand = " << j << " " << xi[j*m] <<    
                                      " " << xi[j*m+1] <<   
                                      " " << xi[j*m+2] <<   
                                      " " << zi[j] << "\n";
        //z assertion
        if (std::abs(xi[j*m+2]-zi[j]) > 1e-6)
            cout << "Error = " << j << " " << xi[j*m+2]-zi[j] << "\n";
    }


    delete [] a;
    delete [] b;
    delete [] n_1d;
//    delete [] xd;
    delete [] xi;
    delete [] zi;

    return;
}


int main ( )

//****************************************************************************80
{

    test03 ( );

    return 0;
}
//*********
