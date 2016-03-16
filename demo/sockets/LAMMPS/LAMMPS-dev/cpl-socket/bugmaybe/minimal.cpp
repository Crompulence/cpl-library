// STL
#include <fstream>
#include <string>
#include <memory>
#include <iostream>

// MPI
#include "mpi.h"

// LAMMPS
#include "lammps.h"
#include "input.h"
#include "modify.h"
#include "fix_ave_chunk.h"

int main(int narg, char **arg)
{

    MPI_Init(&narg, &arg);

    // Create lammps instance on socket realm communicator
    auto lmp = std::make_unique <LAMMPS_NS::LAMMPS>
    (
        1, arg, MPI_COMM_WORLD 
    );

    // Open LAMMPS input script and perform setup line-by-line
    std::ifstream inputFile (arg[1]);
    if (!inputFile.is_open())
    {
        std::cerr << "ERROR: Could not open LAMMPS setup script "
                  << arg[1] << std::endl;
        MPI_Abort (MPI_COMM_WORLD, 1);
    }
    std::string line;
    while (!inputFile.eof())
    {
        std::getline (inputFile, line);
        lmp->input->one (line.c_str());
    }



    // This works
    std::string cmd = "fix testfix all ave/atom 1 1 1 vx vy vz"; 
    lmp->input->one (cmd.c_str());
    auto testfix = lmp->modify->find_fix("testfix");
    auto testx = lmp->modify->fix[testfix]->memory_usage();
    std::cout << "Memory usage of test fix (no chunk) = " << testx << std::endl;



    // This doesn't work 
    line = "run " + std::to_string(10);
    lmp->input->one (line.c_str());
    std::string cmd1 = "compute testchunk all chunk/atom bin/3d ";
                cmd1 += " x lower 1.0 y lower 1.0 z lower 1.0";
    std::string cmd2 = "fix testfixchunk all ave/chunk 1 1 1";
                cmd2 += " testchunk vx vy vz";
    lmp->input->one (cmd1.c_str());
    lmp->input->one (cmd2.c_str());
    int testfixchunk = lmp->modify->find_fix("testfixchunk");
    //LAMMPS_NS::FixAveChunk* fix = dynamic_cast<LAMMPS_NS::FixAveChunk*>(lmp->modify->fix[testfixchunk]);
    //auto fix = dynamic_cast<LAMMPS_NS::FixAveChunk*>(lmp->modify->fix[testfixchunk]);
    auto fix = lmp->modify->fix[testfixchunk];
    fix->setup(0);
    std::cout << "Custom fix memory usage (chunk): " << fix->memory_usage() << std::endl;
    std::cout << "Test fix chunk = " << testfixchunk << std::endl;
    auto testxchunk = lmp->modify->fix[testfixchunk]->memory_usage();
    testxchunk = lmp->modify->fix[testfixchunk]->memory_usage();
    std::cout << "Memory usage of test fix (chunk) = " << testxchunk << std::endl;


    // Tell lammps to run for cpl.nsteps/cpl.timestep_ratio
    line = "run " + std::to_string(10);
    lmp->input->one (line.c_str());
    testxchunk = lmp->modify->fix[testfixchunk]->memory_usage();
    std::cout << "Memory usage of test fix (chunk) = " << testxchunk << std::endl;
    for (int i = 0; i < 10 ; ++i)
    {
        testxchunk = lmp->modify->fix[testfixchunk]->memory_usage();
        std::cout << "Memory usage of test fix (chunk) = " << testxchunk << std::endl;
    }

    MPI_Finalize();

}
