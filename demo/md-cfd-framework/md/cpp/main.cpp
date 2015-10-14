#include "socket.h"

int main (int argc, char** argv)
{
    Socket socket (argc, argv);
    
    socket.InitialiseMD (/* your simulation's topology data here */);    

    for (int iter = 0; iter < socket.nsteps(); iter += socket.timestep_ratio())
    {
        socket.Pack (/* simulation data parameters */);
        socket.Receive();
        socket.Send();
        socket.Unpack (/* simulation data parameters */);

        // Continue simulation for socket.timestep_ratio 
    }
     
}
