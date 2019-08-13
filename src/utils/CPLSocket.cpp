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
#include <chrono>

void CPLSocket::initComms(bool load_param_file) {

    // Split MPI_COMM_WORLD into realm communicators
    CPL::init(realmType, realmComm);
    MPI_Comm_rank(realmComm, &rankRealm);
    MPI_Comm_size(realmComm, &noRealmProcs);
};

void CPLSocket::loadParamFile(std::string fname) {
    if (!paramFileLoaded) {
        CPL::load_param_file(fname);
        paramFileLoaded = true;
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

void CPLSocket::setOutgoingFieldPool(CPL::OutgoingFieldPool& send_pool) {
    field_pool_send = &send_pool;
    field_pool_send->allocateBuffer(sendBuff);
    sendBuffAllocated = true;
    outgoingRuntimeInfo["setup"] = 0.0;
    outgoingRuntimeInfo["update"] = 0.0;
    outgoingRuntimeInfo["pack"] = 0.0;
    outgoingRuntimeInfo["send"] = 0.0;
    outgoingRuntimeInfo["total"] = 0.0;
}

void CPLSocket::setIncomingFieldPool(CPL::IncomingFieldPool& recv_pool) {
    field_pool_recv = &recv_pool;
    field_pool_recv->allocateBuffer(recvBuff);
    recvBuffAllocated = true;
    outgoingRuntimeInfo["setup"] = 0.0;
    outgoingRuntimeInfo["update"] = 0.0;
    outgoingRuntimeInfo["unpack"] = 0.0;
    outgoingRuntimeInfo["receive"] = 0.0;
    outgoingRuntimeInfo["total"] = 0.0;
}

#define  get_delta(t1, t2) \
    std::chrono::duration<double>(t2-t1).count()

#define tic() \
    std::chrono::high_resolution_clock::now()

void CPLSocket::communicate() {
    bool send_cond = sendBuffAllocated && sendEnabled;
    bool recv_cond = recvBuffAllocated && recvEnabled;
    RuntimeInfoT& ort = outgoingRuntimeInfo;
    RuntimeInfoT& irt = incomingRuntimeInfo;
    std::chrono::time_point<std::chrono::high_resolution_clock> t1, t2;
    if (send_cond) {
        t1 = tic(); 
        field_pool_send->updateAll();
        t2 = tic(); 
        ort["update"] += get_delta(t1, t2);
        t1 = tic(); 
        field_pool_send->packAll();
        t2 = tic(); 
        ort["pack"] += get_delta(t1, t2);
    }
    if (realmType == CPL::md_realm) {
        if (recv_cond) {
            t1 = tic(); 
            receive(*field_pool_recv);
            t2 = tic(); 
            irt["receive"] += get_delta(t1, t2);
        }
        if (send_cond) {
            t1 = tic(); 
            send(*field_pool_send);
            t2 = tic(); 
            ort["send"] += get_delta(t1, t2);
        }

    }
    else {
        if (send_cond) {
            t1 = tic();
            send(*field_pool_send);
            t2 = tic(); 
            ort["send"] += get_delta(t1, t2);
        }
        if (recv_cond) {
            t1 = tic();
            receive(*field_pool_recv);
            t2 = tic();
            irt["receive"] += get_delta(t1, t2);
        }
    }
    if (recv_cond) {
        t1 = tic();
        field_pool_recv->unpackAll();
        t2 = tic(); 
        irt["unpack"] += get_delta(t1, t2);
        t1 = tic();
        field_pool_recv->updateAll();
        t2 = tic();
        irt["update"] += get_delta(t1, t2);
    }
    ort["total"] = ort["update"] + ort["pack"] + ort["send"];
    irt["total"] = irt["update"] + irt["unpack"] + irt["receive"];
}

void CPLSocket::_compute_region_runtime(const std::vector<double>& procs_runtimes,
                                        double& min, double& max, double& avg) {
    if (isRootProcess()){
        max = -1;
        min = 999999;
        avg = 0.0;
        int olap_procs = 0;
        for (int i=0; i < noRealmProcs; i++) {
            double val = procs_runtimes[i];
            if (!(val < 0)) {
               if (val > max) max = val;
               if (val < min) min = val;
               avg += val;
               olap_procs += 1;
            }
        }
        avg /= olap_procs;
    }
}
void CPLSocket::_collect_region_runtime(double proc_time, std::vector<double>& procs_runtimes) {
    double time = -1;
    if (isOlapRegion())
        time = proc_time;
    double* values = NULL;
    if (isRootProcess())
        values = procs_runtimes.data();
    MPI_Gather(&time, 1, MPI_DOUBLE, values, 1, MPI_DOUBLE, 0, realmComm);
}

void CPLSocket::printRuntimeInfo() {
    RuntimeInfoT& ort = outgoingRuntimeInfo;
    RuntimeInfoT& irt = incomingRuntimeInfo;
    double irt_avg_update, irt_min_update, irt_max_update,
           ort_avg_update, ort_min_update, ort_max_update;
    double irt_avg_unpack, irt_min_unpack, irt_max_unpack,
           ort_avg_pack, ort_min_pack, ort_max_pack;
    double irt_avg_receive, irt_min_receive, irt_max_receive,
           ort_avg_send, ort_min_send, ort_max_send;
    double irt_avg_total, irt_min_total, irt_max_total,
           ort_avg_total, ort_min_total, ort_max_total;
    std::string realm_name, outgoing_name, incoming_name = "";
    std::vector<double> irt_procs_runtime(noRealmProcs), ort_procs_runtime(noRealmProcs);

    // Collect and compute times in order: update, pack/unpack, total 
    // Update
    _collect_region_runtime(ort["update"], ort_procs_runtime);
    _collect_region_runtime(irt["update"], irt_procs_runtime);
    _compute_region_runtime(ort_procs_runtime, ort_min_update, ort_max_update, ort_avg_update);
    _compute_region_runtime(irt_procs_runtime, irt_min_update, irt_max_update, irt_avg_update);
    // Pack/unpack 
    _collect_region_runtime(ort["pack"], ort_procs_runtime);
    _collect_region_runtime(irt["unpack"], irt_procs_runtime);
    _compute_region_runtime(ort_procs_runtime, ort_min_pack, ort_max_pack, ort_avg_pack);
    _compute_region_runtime(irt_procs_runtime, irt_min_unpack, irt_max_unpack, irt_avg_unpack);
    // Send/receive 
    _collect_region_runtime(ort["send"], ort_procs_runtime);
    _collect_region_runtime(irt["receive"], irt_procs_runtime);
    _compute_region_runtime(ort_procs_runtime, ort_min_send, ort_max_send, ort_avg_send);
    _compute_region_runtime(irt_procs_runtime, irt_min_receive, irt_max_receive, irt_avg_receive);
    // Total
    _collect_region_runtime(ort["total"], ort_procs_runtime);
    _collect_region_runtime(irt["total"], irt_procs_runtime);
    _compute_region_runtime(ort_procs_runtime, ort_min_total, ort_max_total, ort_avg_total);
    _compute_region_runtime(irt_procs_runtime, irt_min_total, irt_max_total, irt_avg_total);

    if (realmType == CPL::md_realm) {
        realm_name = "MD";
        incoming_name  = "Constrain";
        outgoing_name = "BC";
    }
    else {
        realm_name = "CFD";
        incoming_name  = "BC";
        outgoing_name = "Constrain";
    }
    if (isRootProcess()) { 
        std::stringstream str_out;
        str_out << std::endl;
        str_out << "Coupled " << realm_name << " runtime info (min, max, average):" << std::endl\
                << "   Outgoing fields(" << outgoing_name << " region):" << std::endl\
                << "      Update: " << ort_min_update << "," << ort_max_update << "," << ort_avg_update << std::endl\
                << "      Pack  : " << ort_min_pack << "," << ort_max_pack << "," << ort_avg_pack << std::endl\
                << "      Send  : " << ort_min_send << "," << ort_max_send << "," << ort_avg_send << std::endl\
                << "      Total : " << ort_min_total << "," << ort_max_total << "," << ort_avg_total << std::endl\
                << "   Incomming fields(" << incoming_name << " region):" << std::endl\
                << "      Update: " << irt_min_update << "," << irt_max_update << "," << irt_avg_update << std::endl\
                << "      Unpack  : " << irt_min_unpack << "," << irt_max_unpack << "," << irt_avg_unpack << std::endl\
                << "      Receive : " << irt_min_receive << "," << irt_max_receive << "," << irt_avg_receive << std::endl\
                << "      Total : " << irt_min_total << "," << irt_max_total << "," << irt_avg_total << std::endl;
        std::cout << str_out.str() << std::endl;
    }
    else {
        // Gather 
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
