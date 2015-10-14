#ifndef SOCKET_H_INCLUDED
#define SOCKET_H_INCLUDED

#include "CPL.h"

class Socket : public CPL::SocketMD<double>
{
    Socket();
    ~Socket();
};

#endif // SOCKET_H_INCLUDED
