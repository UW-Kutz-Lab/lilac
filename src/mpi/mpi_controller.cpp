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

//!This is currently a thin wrapper for mpi_send, but will be used as part of the worker queue with more options
void mpi_controller::send_blob(void* in, size_t size, MPI_Datatype type, int dest, int tag, MPI_Comm comm){
    send_job(mpi_job(in, size, type, dest, tag, comm));
}

//!Either sends job to destination or places in queue for sending
//!Destroys the passed job after successful send or addition to queue
void mpi_controller::send_job(mpi_job in){
    if(workers[in.dest].available){
        workers[in.dest].available=false;
        MPI_Send(in.data.get(), in.size, in.type, in.dest, in.tag, in.comm);
    }
    else{
        waiting_jobs[in.dest].push_back(std::move(in));
    }
}
mpi_controller::mpi_job::mpi_job(void* dat, size_t dat_s, MPI_Datatype dtype, int _dest, int _tag, MPI_Comm _comm):
    size(dat_s), type(dtype), dest(_dest), tag(_tag), comm(_comm){
        int dat_size;
        MPI_Type_size(dtype, &dat_size);
        std::unique_ptr<char> data_tmp(new char[dat_size*dat_s]);
        data = std::move(data_tmp);
    }
