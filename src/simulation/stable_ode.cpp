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
#include "stable_ode_tmpl.hpp"
#include "utils/comp_funcs.hpp"
#include "types/type_constructor.hpp"
#include "types/type_register.hpp"

template class type_register<stable_ode>;
//stable_ode functions

/*!
 * This is the function that iterates the ODE system forwards.
 * It applies the operator pre_integration_operations, then integrates
 * from tcur to tcur+tlen, and then applies the operator post_integration_operations.
 * \sa stable::simulate
 */
void stable_ode::iterate_system(){
    actual->iterate_system();
}

void stable_ode::postprocess(input& invals){
    invals.retrieve(inter, "integrator", this);
    type_constructor<stable_ode_tmpl>::create(&actual, inter);
    actual->setname(this->name() + "_actual");
    actual->holder = holder;
    actual->postprocess(invals);
    add_as_parent(actual);
}

stable_ode::~stable_ode(){
}
/*!
 * Performs operations that occur before the integration
 */
void stable_ode::pre_integration_operations(){
    actual->pre_integration_operations();
}
/*!
 * Performs operations after the integration
 */
void stable_ode::post_integration_operations(){

    actual->post_integration_operations();
}
std::string stable_ode::type() const{
    return "stable_ode";
}
/*!
 * Returns the dependencies of the stable_ode class
 * the stable_ode class depends on:
 *      - integrator integrator: The integrator that is being used
 *      - double int_len: The time to integrate between checking for convergence
 *
 * And has optional parameters:
 *      - double t0: The starting time, defaults to 0
 *
 */
std::vector<std::string> stable_ode::dependencies() const{
    std::string deps[] = {"integrator", "t0", "int_len"};
    return make_append(deps, stable::dependencies());
}
double stable_ode::score(){
    return actual->score();
}
const std::type_info& stable_ode::vtype() const{
    return actual->vtype();
}
double stable_ode::get_change(){
    return 0;
}
