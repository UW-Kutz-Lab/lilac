#include "utils/defs.hpp"
#include "mpi_controller.h"
int mpi_controller::get_mpi_rank(){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}
int mpi_controller::get_mpi_size(){
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

void mpi_controller::se send_blob(void* in, size_t dat_size, int target){
    MPI_Send(in, dat_size, MPI_BYTE, 
}
