#include "socket.h"

#include<iostream>

Socket::Socket (int argc, char** argv)
{
    MPI_Init (&argc, &argv);
    CPL::create_comm (CPL::md_realm, realm_communicator_);
}
Socket::~Socket() 
{
    MPI_Finalize();
}

void Socket::InitialiseMD (/* Provide simulation topology parameters here */)
{
    std::cout << "Initialising MD Socket" << std::endl;

    // Get these from CPL Library
    nsteps_ = 100;
    timestep_ratio_ = 50;
}

void Socket::Pack (/* Provide simulation data to be packaged here */) 
{
    std::cout << "MD Socket packaging data" << std::endl;
}

void Socket::Unpack (/* Provide simulation data into which received data is to
                        be unpackaged here */)
{
    std::cout << "MD Socket unpackaging data" << std::endl;
}
