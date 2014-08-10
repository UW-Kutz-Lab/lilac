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
#include "gaba_list.hpp"
#include "utils/comp_funcs.hpp"
#include "utils/item_heads.hpp"
#include "types/type_register.hpp"
#include "writer/writer.h"
#include "engine/engineimp.h"
#include "controller/controller.h"
#include "rhs_type.hpp"
#include "rhs_sde.hpp"
#include "rhs_multirate.hpp"
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/SparseCore>
using namespace Eigen;
//!This class represents the c_elegans system
/*!
 * This class calculates the voltages and currents in the neural network of a c_elegans worm.
 * In addition, it also allows for the ablation of specific neurons within the system
 * to see how the system responds to the removal of neurons.
 */
class c_elegans:public rhs_type<double>, public rhs_sde, public rhs_multirate{
    int cur_ind;
    double dummy;
    constexpr static size_t num_neur = 279;
    typedef Eigen::SparseMatrix<double, RowMajor> sparse_type;
    sparse_type laplacian, AEchem_trans, AEchem_trans_full, ag_full;
    Array<double, num_neur, 1> vmean, sig, Echem, eqS, eqV;
    constexpr static size_t dim_v = num_neur;
    constexpr static size_t dim_s = dim_v;
    double beta, memV, memG, gchem, gelec, tau, EchemEx, EchemInh, ar, ad;
    std::vector<size_t> inj_nodes;
    std::vector<size_t> abl_neur;
    bool has_gone;
    bool first_round;
    public:
    virtual void update();
    std::string type() const;
    void postprocess(input& in);
    std::vector<std::string> dependencies() const;
    int dxdt(const ptr_passer x, ptr_passer dx, double dt);
    int dxdt(const ptr_passer x, val_blocks& blocks,
            ind_vals& inds, double t);
    int dwdt(const ptr_passer x, ptr_passer dx, double t);
    void mod_w(ptr_passer w, double t);
    void initial_condition(ptr_passer in, size_t len);
    //block helper functions
    inline void block_volt(const double* x, rhs_multirate::data_block<double>& block);
    inline void block_cur(const double* x, rhs_multirate::data_block<double>& block);
    inline void sing_volt(const double* x, rhs_multirate::single_data<double>& dat);
    inline void sing_cur(const double* x, rhs_multirate::single_data<double>& dat);
    double test[num_neur*2];
};


bool next_comb(std::vector<size_t>& inval, size_t max_num){
    class comb_helper{
        public:
            typedef std::vector<size_t>::iterator iter_t;
            typedef std::vector<size_t>::const_iterator c_iter_t;
            static bool do_next(iter_t cur, c_iter_t end, size_t max_num){
                if(cur == end){
                    //test to ensure valid ptrs
                    //returns zero, since no ablations usually means something else is important
                    return false;
                }
                iter_t temp = cur;
                temp++;
                //increment current element;
                (*cur)++;
                if((*cur) >= max_num){
                    (*cur) = 0;
                }
                //if is at end
                if(temp == end){
                    return (*cur) == 0;
                }
                //if the current value is now greater than
                //or equal to the next, increment
                //while loop allows bad inputs to reach valid combination states
                //v[0] < v[1] < v[2] etc;
                //also yes this has n^2 runtime. fortunately n is small, like around 4-5 at most
                //and is called fairly rarely AND is simple
                while(*cur >= *temp){
                    *cur = 0;
                    if(do_next(temp, end, max_num)){
                        return true;
                    }
                }
                return false;
            }
    };
    return comb_helper::do_next(inval.begin(), inval.end(), max_num);
}
constexpr size_t c_elegans::num_neur;
constexpr size_t c_elegans::dim_v;
constexpr size_t c_elegans::dim_s;
template class type_register<c_elegans>;
int c_elegans::dxdt(const ptr_passer x,  ptr_passer dx, double dt){
    has_gone=true;
    const double* restr v = x;
    double* restr dv = dx;
    double* restr ds = dv+dim_v;
    const double* restr s = v+dim_v;
    //map eigen arrays over input pointers
    Map<Array<double, dim_v, 1>> vmap(const_cast<double*>(v));
    Map<Array<double, dim_v, 1>> dvmap(dv);
    Map<Array<double, dim_v, 1>> smap(const_cast<double*>(s));
    Map<Array<double, dim_v, 1>> dsmap(ds);
    data_block<double> bl;
    const double* px = x;
    bl.dx = test;
    bl.len=dim_v;
    bl.start=0;
    block_volt(px, bl);
    bl.dx = bl.dx + dim_v;
    bl.start = dim_v;
    block_cur(x, bl);
    return 0;

}
//!Calculates derivative for specified values only
int c_elegans::dxdt(const ptr_passer x, val_blocks& _blocks,
        ind_vals& _inds, double t){
    const double* px = x;
    size_t num_blocks = _blocks.num_blocks;
    size_t num_ind = _inds.num_vals;
    single_data<double>* inds = _inds.vals;
    data_block<double>* blocks = _blocks.blocks;
    has_gone=true;

    for(size_t qq = 0; qq < num_blocks; qq++){
        data_block<double>& cur_block = blocks[qq];
        //see if block is calculating voltage, current, or both
        //if start is past voltage, just do current
        if(cur_block.start >= dim_v){
            block_cur(px, cur_block);
        }
        //case where block is entirely in the voltage
        else if(cur_block.start <= dim_v - cur_block.len){
            block_volt(px, cur_block);
        }
        //block spans voltage and current
        else{
            data_block<double> block;

            block.dx = cur_block.dx;
            block.start = cur_block.start;
            block.len = dim_v-cur_block.start;
            block_volt(px, block);

            block.dx = cur_block.dx + block.len;
            block.start = dim_v;
            block.len = (cur_block.len + cur_block.start) - dim_v;
            block_cur(px, block);
        }
    }
    //process singular data
    for(size_t qq = 0; qq < num_ind; qq++){
        single_data<double>& cdat = inds[qq];
        if(cdat.pos < dim_v){
            sing_volt(px, cdat);
        }
        else{
            sing_cur(px, cdat);
        }
    }
    return 0;
}
inline void c_elegans::block_volt(const double* x, rhs_multirate::data_block<double>& cur_block){
    Map<Array<double, Dynamic, 1>> sol(cur_block.dx, cur_block.len);
    Map<Array<double, dim_v, 1>> vmap(const_cast<double*>((const double*)x));
    Map<Array<double, dim_v, 1>> smap(const_cast<double*>((const double*)x+dim_v));

    size_t startv = cur_block.start;

    auto vmap_s = vmap.segment(startv, cur_block.len);
    auto AEchem_s = AEchem_trans.middleRows(startv, cur_block.len);
    auto lap_s = laplacian.middleRows(startv, cur_block.len);

    sol = (-1.0/tau)*(
            (memG*(vmap_s - memV)) +
            (gchem * (vmap_s*(AEchem_s*(smap*(1-Echem)).matrix()).array())) +
            (gelec*lap_s*vmap.matrix()).array()
            );
    //set injection, zero nodes
    for(auto s : inj_nodes){
        double amp = 2e4;
        if(s >= startv && s < startv+cur_block.len){
            sol[s-startv] += (-1.0/tau)*amp;
        }
    }
    for(auto s : abl_neur){
        if(s >= startv && s < startv+cur_block.len){
            sol[s-startv] = 0;
        }
    }
}

inline void c_elegans::sing_volt(const double* x, rhs_multirate::single_data<double>& dat){
    //test if is zero, then skip if it is
    for(auto s : abl_neur){
        if(s == dat.pos){
            dat.val = 0;
            return;
        }
    }
    Map<Array<double, dim_v, 1>> vmap(const_cast<double*>((const double*)x));
    Map<Array<double, dim_v, 1>> smap(const_cast<double*>((const double*)x+dim_v));
    double& sol_val = dat.val;
    size_t pos = dat.pos;
    auto AEchem_r = AEchem_trans.row(dat.pos);
    auto lap_r = laplacian.row(dat.pos);
    sol_val = (-1.0/tau)*(
            (memG*(x[pos] - memV)) +
            (gchem * x[pos]*(AEchem_r*((smap*(1-Echem)).matrix()))(0, 0)) +
            gelec*(lap_r*vmap.matrix())(0, 0));
    //set injection, zero nodes
    for(auto s : inj_nodes){
        double amp = 2e4;
        if(s==pos){
            sol_val += (-1.0/tau)*amp;
        }
    }

}
inline void c_elegans::block_cur(const double* x, rhs_multirate::data_block<double>& cur_block){
    size_t startv = cur_block.start-dim_v;
    Map<Array<double, Dynamic, 1>> sol(cur_block.dx, cur_block.len);
    Map<Array<double, dim_v, 1>> vmap(const_cast<double*>(x));
    Map<Array<double, dim_v, 1>> smap(const_cast<double*>(x+dim_v));
    //indices are mapping current to current, therefore must subtract votage dimensions
    auto sig_s = sig.segment(startv, cur_block.len);
    auto vmap_s = vmap.segment(startv, cur_block.len);
    auto vmean_s = vmean.segment(startv, cur_block.len);
    auto smap_s = smap.segment(startv, cur_block.len);
    sig_s = 1.0 / (1.0 + (-1*beta*(vmap_s-vmean_s)).exp());
    sol = ar*(sig_s * (1.0-smap_s)) - ad*smap_s;
    for(auto s : abl_neur){
        if(s >= startv && s < startv+cur_block.len){
            sol[s-startv] = 0;
        }
    }
}

inline void c_elegans::sing_cur(const double* x, rhs_multirate::single_data<double>& dat){
    size_t pos = dat.pos-dim_v;
    double& sol_val = dat.val;
    for(auto s : abl_neur){
        if(s==dat.pos){
            sol_val = 0;
            return;
        }
    }
    //indices are mapping current to current, therefore must subtract votage dimensions
    double sig;
    double vmap_v = x[pos];
    double vmean_v = vmean(pos);
    double smap_v = x[pos+dim_v];
    sig = 1.0 / (1.0 + std::exp(-1*beta*(vmap_v-vmean_v)));
    sol_val = ar*(sig * (1.0-smap_v)) - ad*smap_v;
}
int c_elegans::dwdt(const ptr_passer x, ptr_passer _dw, double dt){
    double* dw = _dw;
    std::fill(dw, dw+num_neur*2, 1);
    //doesn't need to set abl_neur to zero, since that is done by mod_w
    return 0;
}
void c_elegans::mod_w(ptr_passer _w, double t){
    if(!abl_neur.empty()){
        double *w = _w;
        for(auto abl : abl_neur){
            w[abl]= w[abl+num_neur] = 0;
        }
    }
}
std::string c_elegans::type() const{
    return "c_elegans";
}
template<class sp_t>
void read_mat(std::string fname, sp_t& in_mat){
    //intel compiler requries this...
    std::ifstream in_f;
    in_f.open(fname.c_str());
    std::vector<Eigen::Triplet<double, int> > intr;
    while(!in_f.eof()){
        int row, col;
        double val=0;
        in_f >> row >> col >> val;
        intr.push_back(Triplet<double, int>(row-1, col-1, val));
    }
    in_mat.setFromTriplets(intr.begin(), intr.end());
}
/*!
 * This function does the processing for the c_elegans class.
 *
 * It initializes the various matrices and reads values from the input files
 */
void c_elegans::postprocess(input& in){
    rhs_type::postprocess(in);
    if(dimension != num_neur*2){
        err("Dimension must be 558, which is double the number of neurons",
                "", "", FATAL_ERROR);
    }
    in.retrieve(beta, "beta", this);
    in.retrieve(tau, "tau", this);
    in.retrieve(gelec, "gelec", this);
    in.retrieve(gchem, "gchem", this);
    in.retrieve(memV, "memV", this);
    in.retrieve(memG, "memG", this);
    in.retrieve(EchemEx, "EchemEx", this);
    in.retrieve(EchemInh, "EchemInh", this);
    in.retrieve(ar, "ar", this);
    in.retrieve(ad, "ad", this);
    std::string ag_fname, a_fname;
    in.retrieve(ag_fname, "ag_mat", this);
    in.retrieve(a_fname, "a_mat", this);
    sparse_type a_m(num_neur, num_neur);
    ag_full.resize(num_neur, num_neur);
    laplacian.resize(num_neur, num_neur);
    read_mat(ag_fname, ag_full);
    read_mat(a_fname, a_m);
    //create transposed sparse matrix AEchem
    AEchem_trans_full.resize(num_neur, num_neur);
    AEchem_trans_full = a_m.transpose();
    AEchem_trans.resize(num_neur, num_neur);
    //do any needed fake iterations, must make more general at some point
    size_t num_comb;
    int iterations;
    in.retrieve(num_comb, "num_comb", this);
    in.retrieve(iterations, "iterations", this);
    in.retrieve(cur_ind, "!start_ind", this);
    inj_nodes.clear();
    inj_nodes.resize(2);
    inj_nodes[0]=276;
    inj_nodes[1]=278;
    abl_neur.resize(num_comb);
    for(auto& val : abl_neur){
        val = 0;
    } 
    if(abl_neur.size() != 1){
        next_comb(abl_neur, num_neur);
    }
    for(int i = 0; i < cur_ind; i++){
        for(int j = 0; j < iterations; j++){
            if(next_comb(abl_neur, num_neur)){
                char ind_str[20];//won't ever have a 20 digit index
                //handy buffer to overflow for those hacking this.
                sprintf(ind_str, "%d", (int)cur_ind);
                err(std::string("Combinations exhausted in index ") + ind_str,
                        "c_elegans::postprocess","rhs/c_elegans.cpp", FATAL_ERROR);
            }
        }
    }
    auto dat_inds =
        std::shared_ptr<writer>(new writer(true));
    dat_inds->add_data(data::create("Ablations", abl_neur.data(), abl_neur.size()), writer::OTHER);
    holder->add_writer(dat_inds);

    //write first ablation data

    //set up dummy connection to toroidal controller for now
    controller* cont;
    in.retrieve(cont, "controller", this);
    auto val = std::make_shared<variable>();
    val->setname("c_elegans_quickfix");
    val->holder = holder;
    val->parse("0.1");
    in.insert_item(val);
    cont->addvar(val);
    in.retrieve(dummy, val->name(), this);
    has_gone=true; //is true at first to allow update of zero index to occur
    first_round=true;
    update();
}



void c_elegans::update(){
    if(has_gone){
        has_gone=false;//prevent more than one update happening without call to rhs
        if(!first_round){
            if(next_comb(abl_neur, num_neur)){
                holder->write_dat();
                char ind_str[20];//won't ever have a 20 digit index
                //handy buffer to overflow for those hacking this.
                sprintf(ind_str, "%d", cur_ind);
                err(std::string("Combinations exhausted in index ") + ind_str,
                        "c_elegans::update","rhs/c_elegans.cpp", FATAL_ERROR);
            }
            auto dat_inds =
                std::shared_ptr<writer>(new writer(true));
            dat_inds->add_data(data::create("Ablations", abl_neur.data(), abl_neur.size()), writer::OTHER);
            holder->add_writer(dat_inds);

        }
        else{
            first_round=false;
        }

        Matrix<double, Dynamic, Dynamic> ag_dense(num_neur, num_neur);
        //insert code to zero it out later

        auto ag_m = ag_full;
        AEchem_trans = AEchem_trans_full;
        auto fncval = [this](int i, int j, double val)->bool{
            for(int kk = 0; kk < (int)this->abl_neur.size(); kk++){
                if((int)this->abl_neur[kk] == i || (int)this->abl_neur[kk] == j){
                    return false;
                }
            }
            return true;
        };
        AEchem_trans.prune(fncval);
        ag_m.prune(fncval);
        ag_dense = ag_m *  Matrix<double, num_neur, num_neur>::Identity();
        auto sumvals = ag_dense.colwise().sum();
        sparse_type temp(num_neur, num_neur);
        //generate the sparse diagonal matrix to build the lapacian
        std::vector<Triplet<double, int> > temp_tr;
        for(int i = 0; i < (int)num_neur; i++){
            temp_tr.push_back(Triplet<double, int>(i, i, sumvals[i]));
        }
        temp.setFromTriplets(temp_tr.begin(), temp_tr.end());
        laplacian = temp - ag_m;
        //initialize the Echem array
        for(size_t i = 0; i < num_neur; i++){
            if(GABAergic[i]){
                Echem[i] = EchemInh;
            }
            else{
                Echem[i] = EchemEx;
            }
        }
        //Initialize the sig array
        for(size_t i = 0; i < num_neur; i++){
            sig[i] = 0.5;
        }

        eqS = sig * (ar/(ar*sig + ad));
        //more initialization of temporary dense matrices
        Matrix<double, Dynamic, Dynamic> ldense(num_neur,num_neur);
        ldense = laplacian*Matrix<double, num_neur, num_neur>::Identity();
        Matrix<double, Dynamic, Dynamic> aedense(num_neur, num_neur);
        aedense= AEchem_trans*Matrix<double, num_neur, num_neur>::Identity();
        Matrix<double, Dynamic, Dynamic> C(num_neur, num_neur);
        //create the C matrix
        C= memG*Matrix<double, num_neur, num_neur>::Identity() + gelec*ldense;
        //initialize matrix to modify diagonal part of C
        Matrix<double, num_neur, 1> tmp =(gchem * aedense * eqS.matrix()); 
        for(size_t i = 0; i < num_neur; i++){
            C(i, i) += tmp(i);
        }
        Matrix<double, num_neur, 1> Ivals;
        Ivals.setZero();
        double amp=2e4;
        Ivals[276]=amp;
        Ivals[278]=amp;
        Matrix<double, num_neur, 1> b;
        //create B vector
        b= gchem*aedense*(eqS * Echem).matrix();
        for(size_t i = 0; i < num_neur; i++){
            b[i] += (memG*memV + Ivals[i]);
        }
        //calculate eqV
        eqV.matrix() = C.inverse()*b;
        vmean = eqV+(1.0/beta) * (1.0/sig - 1).log();
        for(auto val : abl_neur){
            eqV[val] = vmean [val] = eqS[val] = 0;
        };
    }
}

std::vector<std::string> c_elegans::dependencies() const{
    std::string deps[] = {"beta", "memV", "memG", "gchem", "gelec", "num_comb",
        "tau", "EchemEx", "EchemInh", "ar", "ad", "ag_mat", "a_mat", "iterations",
        "controller"};
    return make_append(deps, rhs_type::dependencies());
}

void c_elegans::initial_condition(ptr_passer in, size_t len){
    if(len != dimension){
        rhs_type::initial_condition(in, len);
    }
    double* vals = in;
    for(size_t i = 0; i < num_neur; i++){
        //vals[i] = eqV[i];
        vals[i] = vmean[i];
        vals[i+dim_v] = eqS[i];
    }
}
