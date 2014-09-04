#ifndef MULTIRATE_TMPL_HPP
#define MULTIRATE_TMPL_HPP
#include "defs.hpp"
#include "integrator.h"
#include "types/float_traits.hpp"
#include "rhs/rhs_multirate.hpp"
template<class T>
class multirate_tmpl:public integrator{
    

    //integrator parameters
    MAKE_ALIGNED const static real_type a[] = {0.2, 0.3, 0.8, 8.0/9, 1, 1};
    MAKE_ALIGNED const static real_type b1 = 0.2;
    MAKE_ALIGNED const static real_type b2[] = {3.0/40, 9.0/40};
    MAKE_ALIGNED const static real_type b3[] = {44.0/45, -56.0/15, 32.0/9};
    MAKE_ALIGNED const static real_type b4[] = {19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729};
    MAKE_ALIGNED const static real_type b5[] = {9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656};
    MAKE_ALIGNED const static real_type b6[] = {35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84};
    MAKE_ALIGNED const static real_type c4[] = {5179.0/57600, 0, 7571.0/16695, 393.0/640,
        -92097.0/339200, 187.0/2100, 1.0/40};
    MAKE_ALIGNED       const static real_type c5[] = {35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84, 0};



    constexpr char max_depth=32;
    rhs_multirate* rhs_m;
    typedef real_type float_traits<T>::type;
    struct continuous_block{
        size_t len;
        size_t start_pos;
        T* blocks;
    };
    //These serve as holding spaces for the runge_kutta methods
    T* restr f0, * restr f1, * restr f3, * restr f4, * restr f5;
    T* restr f6, * restr f7;
    T* restr f_hold;
    struct scattered_dat{
        size_t len;
        //Here I use vector instead of pointer since I will be adding and removing elements very regularly
        std::vector<size_t> val_pos; //the actual position of values[i] is val_pos[i]
        std::vector<T> values;
    };

    struct integration_section{
        unsigned char cur_level; //level will be limited to at most 32 different levels
        std::vector<continuous_block> blocks;
        scattered_dat scattered;
    };

    double t_cur;
    void step(integration_section& sec){}
    void step_block(continuous_block& blk, tc, tn, inside_pos){
        
    }
    public:
        T* restr cur_sol;
};
#endif
