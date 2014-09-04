#ifndef MPI_CONTROLLER_H
#define MPI_CONTROLLER_H
#include "threaded_queue.cpp"
#include "utils/defs.hpp"
#include <functional>
#include <utility>
#include <thread>
#include <deque>
#include <map>
/*!
 * This class provides an event-based controller for controlling other mpi nodes
 * It follows a master-worker type approce, where one core is constantly delegating
 * available work out to other available cpus
 *
 * A search class will register to be called when a process returns a value,
 * and also provides a utility for sending data to another process.
 *
 */

class mpi_controller{
    public:
        struct mpi_job{
            std::unique_ptr<char> data;
            size_t size;
            MPI_Datatype type;
            int dest;
            int tag;
            MPI_Comm comm;
            mpi_job(size_t dat_s, int _dest, int _tag, MPI_Comm _comm, double* dat=0);
            bool operator <(const mpi_job& job) const{
                return data.get() < job.data.get();
            }
        };       
        typedef std::function<void(mpi_job)> mpi_sender;
        typedef std::function<bool(void)> work_checker;
        typedef std::function<void(mpi_job, mpi_sender)> event_handler;
        int get_mpi_rank() const;
        int get_mpi_size() const;
        //registers a callable object as the event handler
        mpi_controller(event_handler event_h, work_checker work_c);
    private:
        struct worker{
            int id;
            bool available;
        };
        bool has_work;
        void event_loop();
        void send_recv();
        bool async_send();
        bool async_recv();
        bool handle_events();
        bool is_working();
        //for now, a single event handler
        //however, messages will be tagged in the future to dispatch
        //to different handlers
        event_handler event;
        work_checker more_work;
        std::vector<worker> workers;
        //!holds jobs for a worker while worker is busy
        std::vector<threaded_queue<mpi_job>> waiting_jobs;

        //!holds jobs that are ready to be sent
        threaded_queue<mpi_job> send_queue;

        //holds jobs that have been recieved and are waiting for processing
        threaded_queue<mpi_job> recv_queue;
        //
        //Puts job in queue and destroys passed job
        void send_job(mpi_job in);

        //sends active jobs, delets jobs which have completed sending
        //holds waiting send requests
        std::map<std::unique_ptr<MPI_Request>, mpi_job> active_sends;
        std::map<std::unique_ptr<MPI_Request>, mpi_job> active_recvs;

        //runs the event loop
        std::thread event_thread;
};
#endif
