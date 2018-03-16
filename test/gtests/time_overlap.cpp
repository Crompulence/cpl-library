 // basic file operations
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

#include "cpl.h"
#include "CPL_field.h"
#include "CPL_force.h"
#include "CPL_ndArray.h"

#define trplefor(Ni,Nj,Nz) for (int i = 0; i<Ni; ++i){ \
                           for (int j = 0; j<Nj; ++j){ \
                           for (int k = 0; k<Nz; ++k)

using namespace std::chrono;

//Test for CPL::ndArray - setup and array size 
int main() 
{

    bool full_overlap = false;

    int Npercell = 20;
    int nd = 3; 
    int N = 10;
    int icell = N;
    int jcell = N;
    int kcell = N;
    double Lx = 5.0;
    double Ly = 2.0;
    double Lz = 7.0;
    double xorigin = 0.0;
    double yorigin = 0.0;
    double zorigin = 0.0;
    double min[3] = {xorigin, yorigin, zorigin};
    double max[3] = {Lx, Ly, Lz};
    double dx = Lx/float(icell);
    double dy = Ly/float(jcell);
    double dz = Lz/float(kcell);

    //Setup one particle per cell
    double r[3] = {0.0, 0.0, 0.0};
    double v[3] = {0.0, 0.0, 0.0};
    double a[3] = {0.0, 0.0, 0.0};
    double m=1.; double e=1.; double s;
    double radius;
    std::vector<double> F;

    //Setup timer objects
    high_resolution_clock::time_point begin;
    high_resolution_clock::time_point end;

    //Ten cells for domain 0 to 1, cell 0.1 each
    s = 0.2*std::min(std::min(dx, dy), dz);

    ///////////////////////////////////
    //      Use overlap case
    ///////////////////////////////////
    std::map <std::string, std::string> args_map_overlap
    {
        { "use_overlap", "1" },
    };
    double overlaptime;
    if (full_overlap){
        CPLForceDrag overlap(nd, icell, jcell, kcell, args_map_overlap);
        overlap.set_minmax(min, max);
        overlap.volSums->quickcalc = false;
        begin = high_resolution_clock::now();
        trplefor(icell,jcell,kcell){

            for (int i = 0; i<Npercell; ++i){
                r[0] = i/double(icell)+0.5*dx + 2.*dx*(std::rand()/double(RAND_MAX)-0.5);
                r[1] = j/double(jcell)+0.5*dy + 2.*dy*(std::rand()/double(RAND_MAX)-0.5);
                r[2] = k/double(kcell)+0.5*dz + 2.*dz*(std::rand()/double(RAND_MAX)-0.5);
                if (r[0] < min[0]) r[0] = min[0] + 2.*s;
                if (r[1] < min[1]) r[1] = min[1] + 2.*s;
                if (r[2] < min[2]) r[2] = min[2] + 2.*s;
                if (r[0] > max[0]) r[0] = max[0] - 2.*s;
                if (r[1] > max[1]) r[1] = max[1] - 2.*s;
                if (r[2] > max[2]) r[2] = max[2] - 2.*s;
                overlap.pre_force(r, v, a, m, s, e);
            }
        } } }
        end = high_resolution_clock::now();
        overlaptime = duration_cast<microseconds>( end - begin ).count();
        std::cout << "overlap time = " << overlaptime << "e-6 s"   << std::endl;
    } else {
        overlaptime = 0.0;
    }

    ///////////////////////////////////
    //      quickcal method
    ///////////////////////////////////
    CPLForceDrag overlapquick(nd, icell, jcell, kcell, args_map_overlap);
    overlapquick.set_minmax(min, max);
    overlapquick.volSums->quickcalc = true;
    begin = high_resolution_clock::now();
    trplefor(icell,jcell,kcell){

        for (int i = 0; i<Npercell; ++i){
            r[0] = i/double(icell)+0.5*dx + 2.*dx*(std::rand()/double(RAND_MAX)-0.5);
            r[1] = j/double(jcell)+0.5*dy + 2.*dy*(std::rand()/double(RAND_MAX)-0.5);
            r[2] = k/double(kcell)+0.5*dz + 2.*dz*(std::rand()/double(RAND_MAX)-0.5);
            if (r[0] < min[0]) r[0] = min[0] + 2.*s;
            if (r[1] < min[1]) r[1] = min[1] + 2.*s;
            if (r[2] < min[2]) r[2] = min[2] + 2.*s;
            if (r[0] > max[0]) r[0] = max[0] - 2.*s;
            if (r[1] > max[1]) r[1] = max[1] - 2.*s;
            if (r[2] > max[2]) r[2] = max[2] - 2.*s;
            overlapquick.pre_force(r, v, a, m, s, e);
        }
    } } }
    end = high_resolution_clock::now();
    double overlapquicktime = duration_cast<microseconds>( end - begin ).count();
    std::cout << "overlapquick time = " << overlapquicktime << "e-6 s"   << std::endl;
    
    std::cout << "quickcheck  " << overlapquick.volSums->quickcheck << 
                 " facequickcheck "  << overlapquick.volSums->facequickcheck <<  
                 " edgenonquickcheck " << overlapquick.volSums->edgequickcheck << 
                 " cornernonquickcheck " << overlapquick.volSums->nonquickcheck << std::endl;

    ///////////////////////////////////
    //      No overlap case -- PCM
    ///////////////////////////////////
    std::map <std::string, std::string> args_map_nooverlap
    {
        { "use_overlap", "0" },
    };
    CPLForceDrag nooverlap(nd, icell, jcell, kcell, args_map_nooverlap);
    nooverlap.set_minmax(min, max);
    begin = high_resolution_clock::now();
    trplefor(icell,jcell,kcell){
        for (int n = 0; n<Npercell; ++n){
            r[0] = i/double(icell)+0.5*dx + 2.*dx*(std::rand()/double(RAND_MAX)-0.5);
            r[1] = j/double(jcell)+0.5*dy + 2.*dy*(std::rand()/double(RAND_MAX)-0.5);
            r[2] = k/double(kcell)+0.5*dz + 2.*dz*(std::rand()/double(RAND_MAX)-0.5);
            if (r[0] < min[0]) r[0] = min[0] + 2.*s;
            if (r[1] < min[1]) r[1] = min[1] + 2.*s;
            if (r[2] < min[2]) r[2] = min[2] + 2.*s;
            if (r[0] > max[0]) r[0] = max[0] - 2.*s;
            if (r[1] > max[1]) r[1] = max[1] - 2.*s;
            if (r[2] > max[2]) r[2] = max[2] - 2.*s;
            nooverlap.pre_force(r, v, a, m, s, e);
        }
    } } }
    end = high_resolution_clock::now();
    double nooverlaptime = duration_cast<microseconds>( end - begin ).count();
    std::cout << "no overlap time = " << nooverlaptime << "e-6 s"   << std::endl;
    std::cout << "ratio overlap/no_overlap " << overlaptime/nooverlaptime <<  " quick_overlap/no_overlap " << overlapquicktime/nooverlaptime << std::endl;

    return 0;
}






