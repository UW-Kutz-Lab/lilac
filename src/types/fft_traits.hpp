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
#ifndef FFT_TRAITS_HPP
#define FFT_TRAITS_HPP
#include "float_traits.hpp"
#include <complex>
#include <array>
template<class T>
struct fft_traits{};

template<>
struct fft_traits<double>: public float_traits<double>{
    static constexpr size_t fft_seq=1;
    static constexpr bool mappable_to_complex=false;
    static constexpr bool is_double=true;
    std::complex<double> to_complex(const double& val){
        return val;
    }
};
template<>
struct fft_traits<float>: public float_traits<float>{
    static constexpr size_t fft_seq=1;
    static constexpr bool mappable_to_complex=false;
    static constexpr bool is_double=false;
    std::complex<float> to_complex(const float& val){
        return val;
    }
};
template<>
struct fft_traits<std::complex<double>>:public float_traits<std::complex<double>>{
    static constexpr size_t fft_seq=1;
    static constexpr bool mappable_to_complex=true;
    static constexpr bool is_double=true;
    //types which are mappable to complex<float> or complex<double> don't need to define a to_complex function
};
template<>
struct fft_traits<std::complex<float>>:public float_traits<std::complex<float>>{
    static constexpr size_t fft_seq=1;
    static constexpr bool mappable_to_complex=true;
    static constexpr bool is_double=false;
    //types which are mappable to complex<float> or complex<double> don't need to define a to_complex function
};
#endif
