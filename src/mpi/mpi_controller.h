#ifndef mpi_controller.h
#define mpi_controller.h
#include "utils/ptr_passer.hpp"
#include <functional>
#include <utility>
#include <vector>
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
                }
            }
        void send_blob(void* in, size_t dat_size, MPI_Datatype type, int target);
    private:
        //for now, a single event handler
        //however, messages will be tagged in the future to dispatch
        //to different handlers
        event_handler event;
        //max size must be set by the underlying searcher when in use
        max_size=-1;
        std::vector<int> open_units;
};
#endif
