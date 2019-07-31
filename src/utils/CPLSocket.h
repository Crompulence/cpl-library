/*

    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________
     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________
      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________
       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________
        _\/\\\_____________\/\\\/////////____\/\\\_____________
         _\//\\\____________\/\\\_____________\/\\\_____________
          __\///\\\__________\/\\\_____________\/\\\_____________
           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_
            _______\/////////__\///______________\///////////////__


                         C P L  -  L I B R A R Y

           Copyright (C) 2012-2015 Edward Smith & David Trevelyan

License

    This file is part of CPL-Library.

    CPL-Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CPL-Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CPL-Library.  If not, see <http://www.gnu.org/licenses/>.

Description

    "CPLSocketLAMMPS" class for interfacing with CPL-Library.

Author(s)

    Eduardo Ramos
*/

#ifndef CPL_SOCKET_H_INCLUDED
#define CPL_SOCKET_H_INCLUDED

#include <vector>
#include <valarray>
#include "mpi.h"
#include "CPL_ndArray.h"
#include "CPL_field.h"
#include "cpl.h"
#include "TransmittingField.h"
#include <map>
#include <iostream>

typedef std::map<std::string, double> RuntimeInfoT;

class CPLSocket {

public:
    
    // Construct from no arguments
    CPLSocket(int realm_type) : noRealmProcs(0), myProcCoords(3), olapRegion(), bcRegion(),
                  cnstRegion(), bcPortionRegion(), cnstPortionRegion(),
                  procGrid(3), realmType(realm_type), cfdCells(3),
                  sendBuffAllocated(false), recvBuffAllocated(false),
                  sendEnabled(true), recvEnabled(true), paramFileLoaded(false),
                  field_pool_recv(NULL), field_pool_send(NULL) {}
    CPLSocket() : CPLSocket(-1) {}
    virtual ~CPLSocket(){}
    // Total number of timesteps and timestep ratio
    int realmType;
    int noRealmProcs;
    int nSteps;
    int timestepRatio;
    int initialStep;
    double dt;
    CPL::IncomingFieldPool* field_pool_recv;
    CPL::OutgoingFieldPool* field_pool_send;
    RuntimeInfoT outgoingRuntimeInfo;
    RuntimeInfoT incomingRuntimeInfo;
    
    // Initialisation routines 
    void initComms (bool load_param_file=false);
    virtual void init ();

    // TODO: Change to send/recv with a list of Incomming/Outgoing fields list
    //       Change them to non-abstract
    void receive(CPL::IncomingFieldPool& field_pool);
    void send(CPL::OutgoingFieldPool& field_pool);

    // Useful information for main level program
    const MPI_Comm realmCommunicator() {return realmComm;}
    const MPI_Comm cartCommunicator() {return cartComm;}
    bool isRootProcess() {return (rankRealm == 0);}

    bool isOlapRegion(){return CPL::overlap();};
    //TODO: Implement this
    bool isBcRegion(){return false;};
    bool isCnstRegion(){return false;};

    // Clean up MPI/CPL communicators
    void finalizeComms();

    //Grid and physical info associated with communication regions
    CPL::PortionField bcRegion;
    CPL::PortionField cnstRegion;
    CPL::PortionField olapRegion;
    // Portions of the regions in the local processor
    CPL::PortionField bcPortionRegion;
    CPL::PortionField cnstPortionRegion;


    void communicate();
    // Cartesian coordinates of the processor
    std::vector<int> myProcCoords;

    // Communicators for use with CPL_Library
    MPI_Comm realmComm;
    MPI_Comm cartComm;
    std::vector<int> procGrid;

    // CFD grid
    double dx, dy, dz;
    std::vector<int> cfdCells;
    // Physical info about the realm 
    CPL::Domain realmDomain;

    // Rank of this processor in realm and cartComm
    int rankRealm;
    int rankCart;

    // Data to be sent/received with CPL-Library
    CPL::ndArray<double> sendBuff;
    CPL::ndArray<double> recvBuff;

    // Configuration particular to child implementation
    // of boundary conditions and constrains
    virtual void configureBc(int mode) {};
    virtual void configureCnst(int mode) {};

    bool sendBuffAllocated, recvBuffAllocated;
    bool sendEnabled, recvEnabled;

    // void allocateBuffers(const CPL::OutgoingFieldPool& field_pool_send, 
    //                      const CPL::IncomingFieldPool& field_pool_recv);
    void setOutgoingFieldPool(CPL::OutgoingFieldPool& send_pool);
    void setIncomingFieldPool(CPL::IncomingFieldPool& recv_pool);
    
    bool paramFileLoaded;
    void loadParamFile(std::string fname="config.cpl");
    void printRuntimeInfo();

private:
    void _compute_region_runtime(const std::vector<double>& procs_runtimes, double& min, double& max, double& avg);
    void _collect_region_runtime(double proc_time, std::vector<double>& procs_runtimes);



protected:
    // Internal routines
    virtual void setTimingInfo()=0;
    virtual void setCartCommInfo()=0;
    virtual void setRealmDomainInfo()=0;
    void getTopology_();
    void getTiming_();
    void getPortionField_(std::vector<int>& region_limits, 
                                 CPL::PortionField& field);
};

#endif // CPL_SOCKET_H_INCLUDED
