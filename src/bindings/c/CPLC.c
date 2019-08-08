#include "CPLC.h"

void CPLC_init(int  calling_realm, MPI_Comm* returned_realm_comm) {
    int returned_realm_comm_int;
    CPLC_init_Fort(calling_realm, &returned_realm_comm_int);
    (*returned_realm_comm) = MPI_Comm_f2c(returned_realm_comm_int);
}

void CPLC_setup_cfd(MPI_Comm icomm_grid, double xyzL[], 
                    double xyz_orig[], int ncxyz[]) {

    MPI_Fint icomm_grid_int = MPI_Comm_c2f(icomm_grid);
    CPLC_setup_cfd_Fort(icomm_grid_int, xyzL, xyz_orig, ncxyz);
}

void CPLC_setup_md(MPI_Comm icomm_grid, double xyzL[], 
                   double xyz_orig[]) {
    MPI_Fint icomm_grid_int = MPI_Comm_c2f(icomm_grid);
    CPLC_setup_md_Fort(icomm_grid_int, xyzL, xyz_orig);
}
