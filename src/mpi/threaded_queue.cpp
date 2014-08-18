#ifndef THREADED_QUEUE_CPP
#define THREADED_QUEUE_CPP
#include <mutex>
#include <deque>
template<class T>
class threaded_queue{
    std::deque<T> queue;
    std::mutex data_lock;
    public:
    void push_back(T&& inval){
        data_lock.lock();
        queue.push_back(std::move(inval));
        data_lock.unlock();
    }
    void push_back(const T& inval){
        data_lock.lock();
        queue.push_back(inval);
        data_lock.unlock();
    }
    T&& pop_front(){
        data_lock.lock();
        T val = std::move(queue.front());
        queue.pop_front();
        data_lock.unlock();
        return std::move(val);
    }
    bool empty(){
        //empty is a thread-safe operation-
        //But I don't want the calling this in the middle of a pop operation
        //and then trying to pop an empty queue afterwards
        data_lock.lock();
        bool is_empty = queue.empty();
        data_lock.unlock();
        return is_empty;
    }
};
#endif
