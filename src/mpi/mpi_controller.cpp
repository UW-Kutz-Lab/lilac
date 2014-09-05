/*
Copyright (c) 2014, Sam Schetterer, Nathan Kutz, University of Washington
Authors: Sam Schetterer
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#include "mpi_controller.h"
#include <algorithm>
#include <thread>

mpi_controller::mpi_controller(mpi_controller::event_handler event_h,
        mpi_controller::work_checker work_c): event(event_h), more_work(work_c){
    has_work=true;
    event_thread = std::move(std::thread([this](){this->event_loop();}));
}

int mpi_controller::get_mpi_rank() const{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}
int mpi_controller::get_mpi_size() const{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
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



//!Function for MPI event thread
void mpi_controller::event_loop(){
    auto sender_lambda = [this](mpi_job j){
        this->send_job(std::move(j));
    };
    while(has_work){
        int ms_sleep_for=100;
        bool work = !recv_queue.empty();
        if(work){
            ms_sleep_for=100;
            event(recv_queue.pop_front(), sender_lambda);
        }
        else{
            ms_sleep_for = std::min(int(ms_sleep_for*1.5), 5000);
            std::this_thread::sleep_for(std::chrono::milliseconds(ms_sleep_for));
        }
    }
}
//!Main MPI loop-only one can exist within a program.
void mpi_controller::send_recv(){
    int ms_sleep_for=100;
    //continue while-jobs are still being sent
    //jobs are being recieved
    //processes are currently working
    //there is more work to be done
    while(!send_queue.empty() || !recv_queue.empty() || is_working() || more_work()){
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
    has_work=false;
    event_thread.join();
}

bool mpi_controller::async_recv(){
    bool work_done=false;
    MPI_Status status;
    int flag;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    while(flag){
        work_done=true;
        std::unique_ptr<MPI_Request> req(new MPI_Request);
        int job_s=0;
        MPI_Get_count(&status, MPI_DOUBLE, &job_s);

        mpi_job job(job_s, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD);

        MPI_Irecv(job.data.get(), job.size, job.type, job.dest, status.MPI_TAG, job.comm, req.get());
        active_recvs.insert(std::make_pair(std::move(req), std::move(job)));
        //see if any more receives can be done
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    }
    return work_done;
}

bool mpi_controller::handle_events(){
    bool work_done=false;
    //This simple adds the recieved data to a queue,
    //instead of actually passing to the event handler.
    //I could use a condition variable here, but just having the
    //handler poll for new data occasionally is so much easier and simpler
    //
    //It is assumed that the events don't need to be instantly handled, as opposed
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
        std::unique_ptr<MPI_Request> req(new MPI_Request);
        MPI_Isend(job.data.get(), job.size, job.type, job.dest, job.tag,
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

bool mpi_controller::is_working(){
    return std::any_of(workers.begin(), workers.end(),
            [](const worker& w){return !w.available;});
}

mpi_controller::mpi_job::mpi_job(size_t dat_s, int _dest,
        int _tag, MPI_Comm _comm, double* dat): size(dat_s), type(MPI_DOUBLE),
    dest(_dest), tag(_tag), comm(_comm){
        int dat_size;
        MPI_Type_size(type, &dat_size);
        std::unique_ptr<char> data_tmp(new char[dat_size*dat_s]);
        data = std::move(data_tmp);
        if(dat){
            char* cdat = (char*)dat;
            std::copy(cdat, cdat+dat_size*dat_s, data.get());
        }
    }
