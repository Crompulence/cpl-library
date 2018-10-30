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

    "CPLSocket" class for interfacing with CPL-Library.

Author(s)

    Eduardo Ramos
*/

#include "CPLSocket.h"
#include "CPL_cartCreate.h"
#include "cpl.h"
#include <iostream>

void CPLSocket::initComms(bool load_param_file) {

    // Split MPI_COMM_WORLD into realm communicators
    CPL::init(realmType, realmComm);
    MPI_Comm_rank(realmComm, &rankRealm);

};

void CPLSocket::loadParamFile(std::string fname) {
    if (!paramFileLoaded) {
#ifdef JSON_SUPPORT
        CPL::load_param_file(fname);
        paramFileLoaded = true;
#else
        paramFileLoaded = false;
#endif

    }
}

void CPLSocket::getPortionField_(std::vector<int>& region_limits, 
                                 CPL::PortionField& field) {
    field.cellBounds = region_limits;
    // Only processors inside the overlap region get the portion info
	if (CPL::is_proc_inside(region_limits.data())) {
        std::valarray<double> bot_left(3);
        CPL::map_cell2coord(region_limits[0], region_limits[2], 
                            region_limits[4], &bot_left[0]);
        std::valarray<double> top_right(3);
        CPL::map_cell2coord(region_limits[1], region_limits[3],
                            region_limits[5], &top_right[0]);
        top_right += std::valarray<double>({dx, dy, dz});
        std::vector<int> field_cells(3);
        CPL::get_no_cells(region_limits.data(), field_cells.data()); 
        field = CPL::PortionField(CPL::Domain(bot_left, top_right-bot_left),
                                  field_cells, region_limits);
    } 
}

void CPLSocket::finalizeComms() {
    CPL::finalize();
};

void CPLSocket::getTopology_() {

    // Cell sizes
    //TODO: Convert into a Field along with the mesh
    dx = CPL::get<double> ("dx");
    dy = CPL::get<double> ("dy");
    dz = CPL::get<double> ("dz");
   
    // Cell bounds for the overlap region
    std::vector<int> region_limits(6);
	CPL::get_olap_limits(region_limits.data());
    getPortionField_(region_limits, olapRegion);
    CPL::get_bnry_limits(region_limits.data());
    getPortionField_(region_limits, bcRegion);
    CPL::get_cnst_limits(region_limits.data());
    getPortionField_(region_limits, cnstRegion);
    CPL::my_proc_portion (bcRegion.cellBounds.data(), region_limits.data());
    getPortionField_(region_limits, bcPortionRegion);
    CPL::my_proc_portion (cnstRegion.cellBounds.data(),region_limits.data()); 
    getPortionField_(region_limits, cnstPortionRegion);
    
}


void CPLSocket::allocateBuffers(const CPL::OutgoingFieldPool& field_pool_send,
                                const CPL::IncomingFieldPool& field_pool_recv) {
    allocateSendBuffer(field_pool_send);
    allocateRecvBuffer(field_pool_recv);
}

void CPLSocket::allocateSendBuffer(const CPL::OutgoingFieldPool& field_pool_send) {
    field_pool_send.allocateBuffer(sendBuff);
    sendBuffAllocated = true;
}

void CPLSocket::allocateRecvBuffer(const CPL::IncomingFieldPool& field_pool_recv) {
    field_pool_recv.allocateBuffer(recvBuff);
    recvBuffAllocated = true;
}

void CPLSocket::communicate(CPL::OutgoingFieldPool& field_pool_send,
                            CPL::IncomingFieldPool& field_pool_recv) {

    bool send_cond = sendBuffAllocated && sendEnabled;
    bool recv_cond = recvBuffAllocated && recvEnabled;
    if (send_cond) {
        field_pool_send.updateAll();
        field_pool_send.packAll();
    }
    if (realmType == CPL::md_realm) {
        if (recv_cond) {
            receive(field_pool_recv);
        }
        if (send_cond) {
            send(field_pool_send);
        }
    }
    else {
        if (send_cond)
            send(field_pool_send);
        if (recv_cond)
            receive(field_pool_recv);
    }
    if (recv_cond) {
        field_pool_recv.unpackAll();
        field_pool_recv.updateAll();
    }
}



void CPLSocket::init() {
    setTimingInfo();
    if (realmType == CPL::md_realm) 
        CPL::set_timing(initialStep, 0, 1.0);
    else {
        CPL::set_timing(0.0, nSteps, dt);
    }
    setCartCommInfo();
    int periods[3] = {0, 0, 0};
    CPL::Cart_create (realmComm, 3, procGrid.data(), periods, 
                      myProcCoords.data(), &cartComm);
    MPI_Comm_rank (cartComm, &rankCart);
    setRealmDomainInfo();
    //TODO: Maybe not a good idea to force origin to (0,0,0)
    double orig[3] = {0.0 , 0.0, 0.0};
    if (realmType == CPL::md_realm) 
        CPL::setup_md (cartComm, &(realmDomain.length()[0]), orig);
    else
        //TODO: The grid should be specified in the configuration file and decouple
        //from the setup_cfd.
        CPL::setup_cfd (cartComm, &(realmDomain.length()[0]), orig, cfdCells.data());
    getTiming_();
    getTopology_();
}


void CPLSocket::getTiming_() {
    timestepRatio = CPL::get<int> ("timestep_ratio");
	nSteps = CPL::get<int> ("nsteps_coupled") * timestepRatio;
}


void CPLSocket::send(CPL::OutgoingFieldPool& field_pool) {
    CPL::send(sendBuff.data(), sendBuff.shapeData(), field_pool.field.cellBounds.data());
};

void CPLSocket::receive(CPL::IncomingFieldPool& field_pool) {
    CPL::recv(recvBuff.data(), recvBuff.shapeData(), field_pool.field.cellBounds.data());
};






