#include <vector>
#include <math.h> 
#include <stdexcept>

#include "CPL_ndArray.h"
#include "CPL_field_map.h"
#include "overlap/overlap.hpp"
#include "interpolate/interp_ho.h"
#include "cpl.h"

using namespace CPL;
//Constructor based on specified size

CPLFieldMap::CPLFieldMap(int nd_, int icell_, 
                   int jcell_, int kcell_, 
                   const std::string& name_):CPLField(nd_,icell_,jcell_,kcell_,name_),nd(nd_), icell(icell_), jcell(jcell_), kcell(kcell_), name(name_)
{
    int arrayShape[4] = {nd_, icell_, jcell_, kcell_};
    array.resize(4, arrayShape);
    set_defaults();
    default_minmax();
    set_dxyz();

}


//Constructor which uses arrayin and sets size and initial value to this
CPLFieldMap::CPLFieldMap(CPL::ndArray<double> arrayin, const std::string& name_):CPLField(arrayin,name_),array(arrayin), name(name_)
{
     nd = arrayin.shape(0);
    icell = arrayin.shape(1);
    jcell = arrayin.shape(2);
    kcell = arrayin.shape(3);
    set_defaults();
    default_minmax();
    set_dxyz();
}

//Constructor which uses size of array in x, y and z
CPLFieldMap::CPLFieldMap(int nd_, CPL::ndArray<double> arrayin, const std::string& name_):CPLField(nd_,arrayin,name_)
, nd(nd_), array(arrayin), name(name_)
{
    icell = array.shape(1);
    jcell = array.shape(2);
    kcell = array.shape(3);
    int arrayShape[4] = {nd, arrayin.shape(1), arrayin.shape(2), arrayin.shape(3)};
    array.resize(4, arrayShape);
    set_defaults();
    default_minmax();
    set_dxyz();

}

void CPLFieldMap::set_Voro_data(std::vector<int> cunts_, int sum_sample_, int total_sample_){

    cuntsvector = cunts_;
    sum_sample = sum_sample_;
    total_sample = total_sample_;
}

void CPLFieldMap::add_to_array_map_field(const int n, const int i,
                            const int j, const int k, const double value){
        if (cuntsvector.size()%4!=0)
            std::cout<<"Wrong cuntsvector size !!!"<<std::endl;
        else if(sum_sample == 0)
            std::cout<<"More sampling points needed !!!"<<std::endl;
        else{
            int i1, i2, i3, cunts_temp;
            auto it = cuntsvector.begin();
            for (int iti = 0; iti < cuntsvector.size()/4; iti++){
                i1 = *(it + iti*4);
                i2 = *(it + iti*4 + 1);
                i3 = *(it + iti*4 + 2);
                cunts_temp = *(it + iti*4 + 3);
                array(n, i1, i2, i3) = array(n, i1, i2, i3) + (double)cunts_temp/sum_sample*value;
            }
            
        }
}

void CPLFieldMap::add_to_array_map_field(const double r[], const double value[]){
    std::vector<int> cell = get_cell(r);
    for (int n=0; n < array.shape(0); ++n)
        add_to_array_map_field(n, cell[0], cell[1], cell[2], value[n]);
}

//If just one value, then doesn't need to be a reference to array
void CPLFieldMap::add_to_array_map_field(const double r[], const double value){
    std::vector<int> cell = get_cell(r);
    add_to_array_map_field(0, cell[0], cell[1], cell[2], value);
}

//If just one value, then doesn't need to be a reference to array
void CPLFieldMap::add_to_array_map_field(const int index, const double r[], const double value){
    std::vector<int> cell = get_cell(r);
    add_to_array_map_field(index, cell[0], cell[1], cell[2], value);
}

double CPLFieldMap::get_array_value_map_field(const int index, int i, int j, int k){ 
        double VORO_cell_rev = 0;
        if (cuntsvector.size()%4!=0)
            std::cout<<"Wrong cuntsvector size !!!"<<std::endl;
        else if(sum_sample == 0)
            std::cout<<"More sampling points needed !!!"<<std::endl;
        else{
            int i1, i2, i3, cunts_temp;
            auto it = cuntsvector.begin();
            for (int iti = 0; iti < cuntsvector.size()/4; iti++){
                i1 = *(it + iti*4);
                i2 = *(it + iti*4 + 1);
                i3 = *(it + iti*4 + 2);
                cunts_temp = *(it + iti*4 + 3);               
                VORO_cell_rev = VORO_cell_rev + (double)cunts_temp/sum_sample*array(index,i1,i2,i3);
            }
                        
        }
        return VORO_cell_rev;
        
}
/*
double CPLFieldMap::get_array_value_map_field(const int index, int i, int j, int k){
        return array(index, i, j, k);
}
*/
double CPLFieldMap::get_array_value_map_field(const int index, const double r[]){
    std::vector<int> cell = get_cell(r);
    return get_array_value_map_field(index, cell[0], cell[1], cell[2]);
}

double CPLFieldMap::get_array_value_map_field(const double r[]){
    return CPLFieldMap::get_array_value_map_field(0, r);
}
std::vector<double> CPLFieldMap::get_array_value_map_field(const std::vector<int> indices, const double r[]){
    std::vector<int> cell = get_cell(r);
    std::vector<double> v;
    for ( auto &n : indices ) { 
        v.push_back(array(n, cell[0], cell[1], cell[2]));
    }
    return v;
}
CPL::ndArray<double> CPLFieldMap::get_array(){
    return array;
}

CPL::ndArray<double>& CPLFieldMap::get_array_pointer(){
    return array;
}
std::string CPLFieldMap::field_name() {
    return name;
}



void CPLFieldMap::set_defaults(){
    quickcalc = true;
}

void CPLFieldMap::default_minmax(){
    for (int i = 0; i < 3; ++i){
        min[i] = 0.0;
        max[i] = 1.0;
    } 
    set_minmax(min, max);
}

//Set minimum and maximum values of field application
void CPLFieldMap::set_minmax(double min_in[], double max_in[]){
    for (int i = 0; i < 3; ++i){
        min[i] = min_in[i];
        max[i] = max_in[i];
    }
    set_dxyz();
}

//If either min/max change or field object, we need to recalculate dx, dy and dz
void CPLFieldMap::set_dxyz(){
    for (int i = 0; i < 3; ++i){
        dxyz[i] = (max[i] - min[i])/array.shape(i+1);
    }

    dA[0] = dxyz[1]*dxyz[2];
    dA[1] = dxyz[0]*dxyz[2];
    dA[2] = dxyz[0]*dxyz[1];
    dV = dxyz[0]*dxyz[1]*dxyz[2];
}


void CPLFieldMap::set_array(CPL::ndArray<double> arrayin){
    array = arrayin;
    set_dxyz();
//    for (int n = 0; n < nd; ++n){
//    for (int i = 0; i < icell; ++i){
//    for (int j = 0; j < jcell; ++j){
//    for (int k = 0; k < kcell; ++k){
//        std::cout << array(n,i,j,k) << std::endl;
//    }}}}
}

void CPLFieldMap::zero_array(){
    //array = 0.0;
    array.zero();
}

//Get dA
std::vector<double> CPLFieldMap::get_dA(){

    std::vector<double> retrnd_dA(3);

    retrnd_dA[0] = dA[0];
    retrnd_dA[1] = dA[1];
    retrnd_dA[2] = dA[2];

    return retrnd_dA;
}

//Get dV
double CPLFieldMap::get_dV(){

    return dV;
}


//Get cell from min/max and dx
std::vector<int> CPLFieldMap::get_cell(const double r[]){
    std::vector<int> cell(3);
    for (int i = 0; i < 3; ++i){
        cell[i] = floor((r[i]-min[i])/dxyz[i]);
//        std::cout << "cell:    " << i << " " << cell[i] << " " << r[i] << " " << min[i] << " " << max[i] << " " << dxyz[i] << std::endl;
        //Check cell is within the domain
//        if (cell[i] > floor((max[i]-min[i])/dxyz[i])){
        if (r[i] > max[i]){
            std::cout << "cell above: " << i << " " << cell[i] << " " << r[i] 
                      << " " << max[i] << " " << dxyz[i] << std::endl;
            if (i == 0) cell[0] = icell-1;
            if (i == 1) cell[1] = jcell-1;
            if (i == 2) cell[2] = kcell-1;
            //cell[i] = int(floor((max[i]-min[i])/dxyz[i]));
            //throw std::domain_error("get_cell Error: Input above domain");
        }
        if (cell[i] < 0){
            std::cout << "cell below: " << i << " " << cell[i] << " " << r[i] 
                      << " " << min[i] << " " << dxyz[i] << std::endl;
            cell[i] = 0;
            //throw std::domain_error("get_cell Error: Input below domain");
        }
    }
    return cell;
}

void CPLFieldMap::add_to_array(const int n, const int i, 
                            const int j, const int k, const double value){
//        std::cout << "n=" << n << " value=" << value << std::endl;
//        std::cout << "i=" << i << " value=" << value << std::endl;
//        std::cout << "j=" << j << " value=" << value << std::endl;
//        std::cout << "k=" << k << " value=" << value << std::endl;
    array(n, i, j, k) = array(n, i, j, k) + value; 
}

//Just add to cell based on where centre of particle falls
//assuming that array of values is the same size as array
void CPLFieldMap::add_to_array(const double r[], const double value[]){
    std::vector<int> cell = get_cell(r);
    for (int n=0; n < array.shape(0); n++)
        add_to_array(n, cell[0], cell[1], cell[2], value[n]);
}

//If just one value, then doesn't need to be a reference to array
void CPLFieldMap::add_to_array(const double r[], const double value){
    std::vector<int> cell = get_cell(r);
    add_to_array(0, cell[0], cell[1], cell[2], value);
}

//If just one value, then doesn't need to be a reference to array
void CPLFieldMap::add_to_array(int index, const double r[], const double value){
    std::vector<int> cell = get_cell(r);
    add_to_array(index, cell[0], cell[1], cell[2], value);
}

// If the radius of the particle is included, add a fraction
//Add fraction of value to cell based on fraction of sphere inside
void CPLFieldMap::add_to_array(const double r[], double s, const double value[]){

    int ip, jp, kp;
    double volume = (4./3.)*M_PI*pow(s,3); 
    std::vector<int> cell = get_cell(r);
    double box[6];

    //Get fraction of sphere in a volume
    double dx = dxyz[0];
    double dy = dxyz[1];
    double dz = dxyz[2];
    //Ratio of radius to cell size
    int nxps = ceil(s/dx);
    int nyps = ceil(s/dy);
    int nzps = ceil(s/dz);
    //Min and max values
    int minc[3], maxc[3];
    //Cell indices
    int i = cell[0]; int j = cell[1]; int k = cell[2];

    //std::cout << "quickcalc " << quickcalc << std::endl;

    // We can check to see if particle is smaller than box
    // If it is then we can get add volume directly or use
    // spherical cap. If not then we need full overlap check
    if ((nxps == 1) & (nyps == 1) & (nzps == 1) & quickcalc){
        // Check if far enough inside that no overlap is possible
        box[0] = (i  )*dx; box[3] = (i+1)*dx;
        box[1] = (j  )*dy; box[4] = (j+1)*dy;
        box[2] = (k  )*dz; box[5] = (k+1)*dz;

        int sumbotfaces = (box[0] > r[0]-s) + (box[1] > r[1]-s) + (box[2] > r[2]-s);
        int sumtopfaces = (r[0]+s > box[3]) + (r[1]+s > box[4]) + (r[2]+s > box[5]);

//        std::cout << "overlap " << " x= " <<
//                    box[0]<< " < " <<r[0]-s << " to " <<
//                    r[0]+s<< " < " <<box[3] << " y= " <<
//                    box[1]<< " < " <<r[1]-s << " to " << 
//                    r[1]+s<< " < " <<box[4] << " z= " <<
//                    box[2]<< " < " <<r[2]-s << " to " << 
//                    r[2]+s<< " < " <<box[5] << " " << 
//                    (box[0] > r[0]-s) << " " <<
//                    (box[1] > r[1]-s) << " " <<
//                    (box[2] > r[2]-s) << " " <<
//                    (r[0]+s > box[3]) << " " <<
//                    (r[1]+s > box[4]) << " " <<
//                    (r[2]+s > box[5]) << " " <<
//                    sumbotfaces << " " << sumtopfaces << std::endl;

        if ((sumbotfaces + sumtopfaces) == 0) {
        //        if ((box[0] < r[0]-s) & 
        //            (r[0]+s < box[3]) &
        //            (box[1] < r[1]-s) & 
        //            (r[1]+s < box[4]) &
        //            (box[2] < r[2]-s) & 
        //            (r[2]+s < box[5])) {
            quickcheck+= 1;
//            std::cout << "overlap calc quickcheck ok"  << " x= " <<
//                        box[0]<< " < " <<r[0]-s << " to " <<
//                        r[0]+s<< " < " <<box[3] << " y= " <<
//                        box[1]<< " < " <<r[1]-s << " to " << 
//                        r[1]+s<< " < " <<box[4] << " z= " <<
//                        box[2]<< " < " <<r[2]-s << " to " << 
//                        r[2]+s<< " < " <<box[5] << std::endl;
            for (int n=0; n < array.shape(0); n++)
                add_to_array(n, i, j, k, value[n]);
            return;

        // Check if only one face is crossed -- Spherical cap case
        // if either sumbotfaces or sumtopfaces > 1
        } else if ((sumbotfaces + sumtopfaces) == 1) {
            facequickcheck+= 1;
            double h; 
            std::vector<int> shift = {0, 0, 0};
            //Check all three bottom faces
            if (sumbotfaces == 1) {
                for (int ixyz=0; ixyz<3; ixyz++){
                    if (box[ixyz] > r[ixyz]-s){
                        h = box[ixyz]-(r[ixyz]-s);
                        shift[ixyz] = -1; 
//                        std::cout << "sumbotfaces " <<  ixyz << " " <<
//                                    box[ixyz]<< " < " << r[ixyz]-s << " to " <<
//                                    r[ixyz]+s<< " < " <<box[ixyz] << h << " " << std::endl;
                        break;
                    }
                }
            //Check all three top faces
            } else if (sumtopfaces == 1) {
                for (int ixyz=0; ixyz<3; ixyz++){
                    if (box[ixyz+3] < r[ixyz]+s){
                        h = (r[ixyz]+s)-box[ixyz+3];
                        shift[ixyz] = 1; 
//                        std::cout << "sumtopfaces " <<  ixyz << " " <<
//                                    box[ixyz]<< " < " << r[ixyz]-s << " to " <<
//                                    r[ixyz]+s<< " < " << box[ixyz+3] << " " << h << " " << std::endl;
                        break;
                    }
                }
            }



            // If only just out, skip this as it causes a 
            // problem when particles on box edges
            if (h/s < 1e-9) {
                std::cout << "particle on cell edge " << " x= " <<
                            box[0]<< " < " <<r[0]-s << " to " <<
                            r[0]+s<< " < " <<box[3] << " y= " <<
                            box[1]<< " < " <<r[1]-s << " to " << 
                            r[1]+s<< " < " <<box[4] << " z= " <<
                            box[2]<< " < " <<r[2]-s << " to " << 
                            r[2]+s<< " < " <<box[5] << " " << h << " "
                            << s << " " <<  h/s << std::endl;
                for (int n=0; n < array.shape(0); n++)
                    add_to_array(n, i, j, k, value[n]);
                return;
            }


            double Vcap = (1./3.)*M_PI*pow(h,2)*(3.*s - h);
            double Vremain = volume - Vcap;
            double frac = Vremain/volume;
//            std::cout << "face overlap " << " x= " <<
//                        box[0]<< " < " <<r[0]-s << " to " <<
//                        r[0]+s<< " < " <<box[3] << " y= " <<
//                        box[1]<< " < " <<r[1]-s << " to " << 
//                        r[1]+s<< " < " <<box[4] << " z= " <<
//                        box[2]<< " < " <<r[2]-s << " to " << 
//                        r[2]+s<< " < " <<box[5] << " " << h << " "
//                        << Vcap << " " <<  frac << std::endl;
            if (h > 2.*s){
                std::cout << "face overlap " << " x= " <<
                            box[0]<< " < " <<r[0]-s << " to " <<
                            r[0]+s<< " < " <<box[3] << " y= " <<
                            box[1]<< " < " <<r[1]-s << " to " << 
                            r[1]+s<< " < " <<box[4] << " z= " <<
                            box[2]<< " < " <<r[2]-s << " to " << 
                            r[2]+s<< " < " <<box[5] << " " << h << " "
                            << Vcap << " " <<  frac << std::endl;
                assert(h < 2.*s);
            }

            assert(Vremain < volume);
            //Add to the cell that the particles centre is inside
            for (int n=0; n < array.shape(0); n++)
                add_to_array(n, i, j, k, frac*value[n]);
            //Add spherical cap to adjacent cell
            ip = i+shift[0]; jp = j+shift[1]; kp = k+shift[2];
            if (ip < 0) ip = 0;
            if (jp < 0) jp = 0;
            if (kp < 0) kp = 0;
            if (ip > array.shape(1)-1) ip = array.shape(1)-1;
            if (jp > array.shape(2)-1) jp = array.shape(2)-1;
            if (kp > array.shape(3)-1) kp = array.shape(3)-1;
            for (int n=0; n < array.shape(0); n++)
                add_to_array(n, ip, jp, kp, (1.-frac)*value[n]);
            return;

        } else if ((sumbotfaces + sumtopfaces) > 1) {
            //Edge or corner case, loop over all surrounding cells
            edgequickcheck+= 1;
            //Set iteration below
            if (box[0] > r[0]-s){minc[0] = -nxps;} else {minc[0] = 0;};
            if (box[1] > r[1]-s){minc[1] = -nyps;} else {minc[1] = 0;};
            if (box[2] > r[2]-s){minc[2] = -nzps;} else {minc[2] = 0;};
            //Set iteration above
            if (box[3] < r[0]+s){maxc[0] = nxps+1;} else {maxc[0] = 1;};
            if (box[4] < r[1]+s){maxc[1] = nyps+1;} else {maxc[1] = 1;};
            if (box[5] < r[2]+s){maxc[2] = nzps+1;} else {maxc[2] = 1;};

//            std::cout << "full overlap needed " << " x= " <<
//                        box[0]<< " < " <<r[0]-s << " to " <<
//                        r[0]+s<< " < " <<box[3] << " y= " <<
//                        box[1]<< " < " <<r[1]-s << " to " << 
//                        r[1]+s<< " < " <<box[4] << " z= " <<
//                        box[2]<< " < " <<r[2]-s << " to " << 
//                        r[2]+s<< " < " <<box[5] << " " <<
//                        minc[0] << " " << minc[1] << " " << minc[2] << " " <<
//                        maxc[0] << " " << maxc[1] << " " << maxc[2] << std::endl;


//            minc[0] = -nxps; minc[1] = -nyps; minc[2] = -nzps;
//            maxc[0] = nxps+1; maxc[1] = nyps+1; maxc[2] = nzps+1;

        }


//        std::cout << "quickcheck  " << quickcheck << 
//                     " facequickcheck "  << facequickcheck <<  
//                     " edgequickcheck " << edgequickcheck << 
//                     " nonquickcheck " << nonquickcheck << std::endl;
        
    } else {
        //Particle is bigger than cell so 
        //set iteration above and below
        nonquickcheck+= 1;
        minc[0] = -nxps; minc[1] = -nyps; minc[2] = -nzps;
        maxc[0] = nxps+1; maxc[1] = nyps+1; maxc[2] = nzps+1;
    }


//    std::cout << "overlap calc "  
//              << dx << " " << dy << " " << dz << " " 
//              << i << " " << j << " " << k << " " 
//              << nxps << " " << nyps << " " << nzps << std::endl;



    //Otherwise we have to check all possible cases
    for (int ic=minc[0]; ic<maxc[0]; ic++) {
    for (int jc=minc[1]; jc<maxc[1]; jc++) {
    for (int kc=minc[2]; kc<maxc[2]; kc++) {

        ip = i+ic; jp = j+jc; kp = k+kc;
        box[0] = (ip)*dx; box[3] = (ip+1)*dx;
        box[1] = (jp)*dy; box[4] = (jp+1)*dy;
        box[2] = (kp)*dz; box[5] = (kp+1)*dz;
        
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
            add_to_array(n, ip, jp, kp, frac*value[n]);
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
void CPLFieldMap::add_to_array(const double r[], double s, const double value_){
    double value[1] = {value_};
    add_to_array(r, s, value);
}

//If just one value, then doesn't need to be a reference to array
void CPLFieldMap::add_to_array(const int index, const double r[], 
                            double s, const double value_){
    //Allocate an empty array up to index
    double value[index];
    for (int i=0; i<index; i++)
        value[i] = 0;
    value[index] = value_;
    add_to_array(r, s, value);
}

// A function which gets value n in cell i,j,k
std::vector<double> CPLFieldMap::get_array_value(const std::vector<int> indices, int i, int j, int k){
    std::vector<double> v;
    for ( auto &n : indices ) { 
        v.push_back(array(n, i, j, k));
    }
    return v;
}

// A function which gets value using position of molecule, no interpolation
std::vector<double> CPLFieldMap::get_array_value(const std::vector<int> indices, const double r[]){
    std::vector<int> cell = get_cell(r);
    std::vector<double> v;
    for ( auto &n : indices ) { 
        v.push_back(array(n, cell[0], cell[1], cell[2]));
    }
    return v;
}

// A function which gets value using position of molecule, no interpolation
double CPLFieldMap::get_array_value(const int index, const double r[]){
    std::vector<int> cell = get_cell(r);
    return array(index, cell[0], cell[1], cell[2]);
}

double CPLFieldMap::get_array_value(const double r[]){
    return CPLFieldMap::get_array_value(0, r);
}
