#ifndef SOCKET_H_INCLUDED
#define SOCKET_H_INCLUDED

#include "mpi.h"
#include "CPL.h"

class Socket : public CPL::SocketMD<double>
{

public:

    Socket (int argc, char** argv);
    ~Socket();

    void InitialiseMD (/* Provide simulation topology parameters here */);
    void Pack (/* Provide simulation data to be packaged here */);
    void Unpack (/* Provide simulation data into which received data is to be
                    unpackaged here */);
    
    int nsteps() {return nsteps_;}
    int timestep_ratio() {return timestep_ratio_;}

private:

    MPI_Comm realm_communicator_;

    int nsteps_;
    int timestep_ratio_;

};

#endif // SOCKET_H_INCLUDED
