#include "mpi_controller.h"
#include <thread>
//enough to cover most bases for now
//idk why mpi won't just let me query the sent datatype
std::map<int, MPI_Datatype> int_to_mpi = {
    {0, MPI_INT}, {1, MPI_FLOAT}, {2, MPI_CHAR}, {3, MPI_C_DOUBLE_COMPLEX},
    {4, MPI_C_FLOAT_COMPLEX}, {5, MPI_DOUBLE}};

static std::map<MPI_Datatype, int> __reverse_types(const std::map<int, MPI_Datatype>& m){
    std::map<MPI_Datatype, int> rm;
    for(const auto& pair : m){
        rm[pair.second]=pair.first;
    }
    return rm;
}

std::map<MPI_Datatype, int> mpi_to_int = __reverse_types(int_to_mpi);
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

void mpi_controller::register_handler(const mpi_controller::event_handler& in){
    event = in;
}

void mpi_controller::send_blob(void* in, size_t size, MPI_Datatype type, int dest, MPI_Comm comm){
    send_job(mpi_job(size, type, dest, comm, in));
}

//!Either places job in send queue or places in wait queue while destination is busy
//!Destroys the passed job after successful send or addition to queue
void mpi_controller::send_job(mpi_job in){
    if(workers[in.dest].available){
        workers[in.dest].available=false;
        send_queue.push_back(std::move(in));
    }
    else{
        waiting_jobs[in.dest].push_back(std::move(in));
    }
}




//!Main MPI loop-only one can exist within a program.
void mpi_controller::send_recv(){
    static int ms_sleep_for=100;
    //don't combine into one big statement to ensure compiler doesn't optimize away calls
    bool work_done = async_recv();
    work_done = handle_events() | work_done;
    work_done = async_send() | work_done;

    //If there was nothing to send or recieve, sleep for a few seconds and try again
    //This frees up the core for something else
    if(!work_done){
        ms_sleep_for = std::min(int(ms_sleep_for*1.5), 5000);
        std::this_thread::sleep_for(std::chrono::milliseconds(ms_sleep_for));
    }
    else{
        ms_sleep_for=100;
    }
}

bool mpi_controller::async_recv(){
    bool work_done=false;
    MPI_Status status;
    int flag;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    while(flag){
        work_done=true;
        std::unique_ptr<MPI_Request> req(new MPI_Request);
        auto dtype = int_to_mpi[status.MPI_TAG];

        int job_s=0;
        MPI_Get_count(&status, dtype, &job_s);

        mpi_job job(job_s, dtype, status.MPI_SOURCE, MPI_COMM_WORLD);

        MPI_Irecv(job.data.get(), job.size, job.type, job.dest, status.MPI_TAG, job.comm, req.get());
        active_recvs.insert(std::make_pair(std::move(req), std::move(job)));
        //see if any more receives can be done
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    }
    return work_done;
}

bool mpi_controller::handle_events(){
    bool work_done=false;
    //This simple adds the recieved data to a new thread, 
    //instead of actually passing to the event handler.
    //I could use a condition variable here, but just having the
    //handler poll for new data occasionally is so much easier and simpler
    //
    //It is assumed that the events don't need top be instantly handled, as opposed
    //to the async data writing operations
    for(auto iter = active_recvs.begin(); iter != active_recvs.end();){
        int flag=0;
        MPI_Status stat;
        MPI_Test(iter->first.get(), &flag, &stat);
        if(flag){
            recv_queue.push_back(std::move(iter->second));
            iter = active_recvs.erase(iter);
        }
        else{
            ++iter;
        }
    }
    return work_done;
}

//!deals with sending jobs out
bool mpi_controller::async_send(){
    bool work_done=false;
    //set up and begins the sending process
    while(!send_queue.empty()){
        work_done=true;
        mpi_job job = send_queue.pop_front();
        //This allows us to keep the same address once scope is left
        //sends are tagged with the type since I couldn't find a way
        //to actually query the type being sent
        std::unique_ptr<MPI_Request> req(new MPI_Request);
        MPI_Isend(job.data.get(), job.size, job.type, job.dest, mpi_to_int[job.type],
                job.comm, req.get());
        active_sends.insert(std::make_pair(std::move(req), std::move(job)));
    }
    //test which sends have finished and can be deleted
    for(auto iter = active_sends.begin(); iter != active_sends.end();){
        int flag;
        MPI_Status stat;
        MPI_Test(iter->first.get(), &flag, &stat);
        if(flag){
            //Flag being non-zero means that the operation has completed
            iter = active_sends.erase(iter);
        }
        else{
            ++iter;
        }
    }
    return work_done;
}


mpi_controller::mpi_job::mpi_job(size_t dat_s, MPI_Datatype dtype, int _dest,
        MPI_Comm _comm, void* dat):
    size(dat_s), type(dtype), dest(_dest), comm(_comm){
        int dat_size;
        MPI_Type_size(dtype, &dat_size);
        std::unique_ptr<char> data_tmp(new char[dat_size*dat_s+2*MPI_BSEND_OVERHEAD]);
        data = std::move(data_tmp);
        if(dat){
            char* cdat = (char*)dat;
            std::copy(cdat, cdat+dat_size*dat_s, data.get());
        }
    }
