#ifndef THREADED_QUEUE_CPP
#define THREADED_QUEUE_CPP
#include <mutex>
#include <deque>
template<class T>
class threaded_queue{
    //A lock free queue would be more preformant
    ////But this is far simpler
    std::deque<T> queue;
    std::mutex data_lock;
    typedef std::lock_guard<std::mutex> lock_g;
    public:
    void push_back(T&& inval){
        lock_g lock(data_lock);
        queue.push_back(std::move(inval));
    }
    void push_back(const T& inval){
        lock_g lock(data_lock);
        queue.push_back(inval);
    }
    T&& pop_front(){
        lock_g lock(data_lock); 
        T val = std::move(queue.front());
        queue.pop_front();
        return std::move(val);
    }
    bool empty(){
        //empty is a thread-safe operation-
        //But I don't want the calling this in the middle of a pop operation
        //and then trying to pop an empty queue afterwards
        lock_g lock(data_lock);
        bool is_empty = queue.empty();
        return is_empty;
    }
};
#endif
