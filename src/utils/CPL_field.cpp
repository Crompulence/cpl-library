#include <vector>
#include <math.h> 
#include <stdexcept>

#include "CPL_ndArray.h"
#include "CPL_field.h"


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

//Get cell from min/max and dx
std::vector<int> CPLField::get_cell(double r[]){
    std::vector<int> cell(3);
    for (int i = 0; i < 3; ++i){
        cell[i] = floor((r[i]-min[i])/dxyz[i]);
        //Check cell is within the domain
        if (cell[i] > floor((max[i]-min[i])/dxyz[i]))
            throw std::domain_error("get_cell Error: Input above domain");

        if (cell[i] < 0)
            throw std::domain_error("get_cell Error: Input below domain");
    }
    return cell;
}



void CPLField::set_array(CPL::ndArray<double> arrayin){
    array = arrayin;
}

CPL::ndArray<double> CPLField::get_array(){
    return array;
}

//Get interpolate of field
std::vector<double> CPLField::interpolate(double r[]){

    int interp[3];
    std::vector<int> cell = get_cell(r);
    std::vector<double> val(array.shape(0));
    std::vector<double> rc(3);
    double f,g,h;


    //Check if at edge of domain and use one sided if needed
    for (int n = 0; n < 3; n++){
        if (cell[n]+1>array.shape(n)){
            interp[n] = -1;
        } else if (cell[0]-1<0){  
            interp[n] = 1;
        } else {
            interp[n] = 0;
        }
        // Get position relative to cell minimum
        rc[n] = r[n] - cell[n]*dxyz[n];
    }

    //Interpolate based on both adjacent cells
    for (int n = 0; n<array.shape(0); n++){
        if (interp[0] == 0){
            f =  (   array(n, cell[0]+1,cell[1],  cell[2]  )
                  -2*array(n, cell[0],  cell[1],  cell[2]  )  
                   + array(n, cell[0]-1,cell[1],  cell[2]  ))/2.*dxyz[0];
        } else if (interp[0] == 1){
            f =  (   array(n, cell[0]+1,cell[1],  cell[2]  )
                    -array(n, cell[0],  cell[1],  cell[2]  ))/dxyz[0];
        } else if (interp[0] == -1){
            f =  (   array(n, cell[0],  cell[1],  cell[2]  )
                    -array(n, cell[0]-1,cell[1],  cell[2]  ))/dxyz[0];
        }

        if (interp[1] == 0){
            g =  (   array(n, cell[0],  cell[1]+1,cell[2]  )
                  -2*array(n, cell[0],  cell[1],  cell[2]  )  
                   + array(n, cell[0],  cell[1]-1,cell[2]  ))/2.*dxyz[1];
        } else if (interp[1] == 1){
            g =  (   array(n, cell[0],  cell[1]+1,cell[2]  )
                  -  array(n, cell[0],  cell[1],  cell[2]  ))/dxyz[1];
        } else if (interp[1] == -1){
            g =  (   array(n, cell[0],  cell[1]  ,cell[2]  )
                   - array(n, cell[0],  cell[1]-1,cell[2]  ))/dxyz[1];
        }

        if (interp[2] == 0){
            h =  (   array(n, cell[0],  cell[1],  cell[2]+1)
                  -2*array(n, cell[0],  cell[1],  cell[2]  ) 
                   + array(n, cell[0],  cell[1],  cell[2]-1))/2.*dxyz[2];
        } else if (interp[2] == 1){
            h =  (   array(n, cell[0],  cell[1],  cell[2]+1)
                   - array(n, cell[0],  cell[1],  cell[2]  ))/dxyz[2];
        } else if (interp[2] == -1){
            h =  (   array(n, cell[0],  cell[1],  cell[2]  )
                  -2*array(n, cell[0],  cell[1],  cell[2]-1))/dxyz[2];
        }

        val[n] = f*rc[0] + g*rc[1] + h*rc[2];
    }

    return val;
}



