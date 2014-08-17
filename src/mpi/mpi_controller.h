#ifndef MPI_CONTROLLER_H
#define MPI_CONTROLLER_H
#include "utils/defs.hpp"
#include <functional>
#include <utility>
#include <deque>
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
        //(data, size(data), sender_rank
        typedef std::function<void(void*, size_t, int)> event_handler;
        int get_mpi_rank();
        int get_mpi_size();
        template<class Lambda>
            void register_handler(Lambda&& in, int handle){
                //use lambda trickery so that callable objects will work here as well
                this->event = [&in](void* v, size_t s, int i){
                    std::forward<Lambda>(in)(v, s, i);
                };
            }
    void send_blob(void* in, size_t size, MPI_Datatype type, int dest, int tag, MPI_Comm comm);
    private:
        struct worker{
            int id;
            bool available;
        };
        struct mpi_job{
            std::unique_ptr<char> data;
            size_t size;
            MPI_Datatype type;
            int dest, tag;
            MPI_Comm comm;
            mpi_job(){}
            mpi_job(void* dat, size_t dat_s, MPI_Datatype dtype, int _dest, int _tag, MPI_Comm _comm);
        };
        //for now, a single event handler
        //however, messages will be tagged in the future to dispatch
        //to different handlers
        event_handler event;
        std::vector<worker> workers;
        std::vector<std::deque<mpi_job>> waiting_jobs;
        //Destroys passed job
        void send_job(mpi_job in);
};
#endif
