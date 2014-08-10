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
#ifndef RHS_MULTIRATE_HPP
#define RHS_MULTIRATE_HPP
#include "utils/ptr_passer.hpp"
#include <vector>
class rhs_multirate{
    public:
        template<class T>
            struct data_block{
                size_t start;//starting index
                size_t len;//length of block
                T* restr dx;
            };
        template<class T>
            struct single_data{
                T& val;
                size_t pos; //original index
            };
        struct val_blocks{
            ptr_passer blocks;
            size_t num_blocks;
            val_blocks(): blocks((void*)0){
                num_blocks = 0;
            }

            val_blocks(ptr_passer p, size_t num_vals):blocks(p), num_blocks(num_vals)
            {}
        };
        struct ind_vals{
            ptr_passer vals;
            size_t num_vals;
            ind_vals():vals((void*)0){
                num_vals = 0;
            }

            ind_vals(ptr_passer p, size_t num_inds):vals(p), num_vals(num_inds)
            {}
        };
        //calculates specified derivatives in blocks and ind_vals
        virtual int dxdt(const ptr_passer x, val_blocks& blocks,
                ind_vals& inds, double t)=0;
};
#endif
