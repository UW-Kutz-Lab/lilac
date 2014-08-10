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
/*! \file This is the .cpp file for the example rhs function, with docs for someone new to C++
 * It contains some basic test functions, as well as extra in-function
 * documentation explaining most steps and various details
 * */
#include "utils/item_heads.hpp"//various headers important to dealing with items, see header for current content
#include "comp_funcs.hpp" //includes various functions for complex variables and some template math ones
#include "types/type_register.hpp"//include file for the type_register class;
//
//
#include "rhs_type.hpp"//This file contains the template for the base class type, along with most basic headers
//The way this works is through polymorphism, which basically means that the engine sees a class of rhs, but is really accessing a class of type example_rhs
//
//Also, note that we do not need a header file. This is since nothing else actually
//needs to see the definition of the rhs class. If this isn't true for another class
//(try to avoid these situations), then just make a header
//
//example_rhs is exposed to the rest of the system since it inherits from rhs and can be instantiated
//through the input system. Nothing else actually has access to the definition/declaration of example_rhs


//!This class is a basic example of how to write an rhs class
/*!This class is an example of a function rhs class.
 * It implements all of the functions needed, and can even be called from
 * the input file. Make sure to read the source file, as it contains a large amount of documentations that is not present in the documentation
 * Notice that it inherits from rhs_type<comp>, since it is a complex value rhs.
 * An rhs that dealt with real variables would inherit from rhs_type<double>.
 *
 * \sa example_rhs.cpp, example_rhs.h, rhs, and item 
 */
class example_rhs:public rhs_type<comp>{
    //up here declare variables that are used internally by the function
    //These may be akin to various parameters, tunable or not
    
    //double represents the 64 bit floating point variables
    double val1, val2;
    
    //std::string represents the string type values
    std::string random_info;

    //comp values are the complex variable types
    //comp_funcs.h contains some optimized functions for manipulating complex variables
    comp something_useless;

    //pointers store a location in memory, and can act as arrays.
    //They offer high performance compared to other options,but also require
    //rather vigilant care since they directly manipulate memory
    //They are useful for storing temporary calculations
    comp* value_holder;

    public: //anything that isn't declared under public is only accessible through
    //functions of the class, but public members can be accessed outside of the class
    //
    //Functions of a c++ class are declared in a header file, but actually defined in a source file
    //This allows one to separate the implementation of the class from telling the engine
    //what the class is, ane makes everything work better in the long run
    //
    //PLEASE FOLLOW THIS PATTERN WITH USER DEFINIED CLASSES
    //
    //these functions delcared in the next block are utility functions
    //used by the engine to sort different parts of the simulation and
    //make everything work nicely together. In each class they must be declared exactly
    //as is here. The detailed documentation, and comments in the .cpp file,
    //have information on what these functions do

    virtual std::string type() const;
    virtual std::vector<std::string> dependencies() const;
    virtual void postprocess(input& indat);


    //This function updates data inside the class whenever a variable it references changes
    //use it to re-do any calculations that are done to initialize the class
    virtual void update();

    //This function calculates the derivative given the current function value
    //Note that instead of pointers, this class takes in an object of type
    //ptr_passer. This holds a pointer and type information to ensure that
    //the rest of the engine, which is quite generic and type-independent,
    //won't pass the wrong type of pointer to this function. Although the 
    //rest of the type-checking is probably strict enough to allow the use of
    //void*, this prevents a user from accidentally passing a bad pointer anyways.
    virtual int dxdt(const ptr_passer x, ptr_passer dx, double t);
    
    //This is a destructor, which frees all of your memory
    //if you allocate any memory, you must free it
    virtual ~example_rhs();
};
//
//****IMPORTANT*****IMPORTANT*****IMPORTANT*****IMPORTANT*****IMPORTANT*
//To have code that compiles with the engine, you need to edit the makefile in
//the rhs directory. It is the file that controls when pieces of code are built.
//Alhough it is possible to simply recompile everything every time, it is incredibly ineffecient
//and wastes and incredible amount of time.
//
//An example can be found, for this file, in the makefile in rhs.
//You pretty much just need to copy and paste a version with your file names
//instead of the example file names. Also, please don't remove the example build commands
//
//
//*************************
//IMPORTANT
//the line of code below makes it so that you can create an instance of the class
//through the input scripts. This, in a roundabout way, is the only method through
//which example_rhs is exposed to the rest of the engine-as an rhs object that was instantiated
//by the input system
template class type_register<example_rhs>;
//If you don't do this, then you can't create this in the input file
//
//That is desirable behavior, for example it is used by the template backends of the integrators


//!This function returns the variable type
/*!
 * This function returns a string that has the typename of your variable. If the variable
 * is meant to be instantiated in an input script, then this function must return the name
 * which corresponds to this class. See rhs_CNLS for a case where that differs from the class name
 * @return Name of the class type
 */
std::string example_rhs::type() const{
    return "example_rhs";
}

//!This function returns a list of the dependencies needed to initialize the class
/*!The engine ensures that these variables will themselves be initialized
* when they are used to initialize the class. An error occurs if the variables are not found
* in the input or there is a loop of dependencies, this is a depends on b, which depends on c,
* which in turn depends back on a.
* 
*
* \note
* This class depends on val1, val2, and random_info
* @return List of dependencies
*/
std::vector<std::string> example_rhs::dependencies() const{
    std::string deps[] = {"val1", "val2", "random_info"};
    //this array holds the dependencies

    //This function appenps an array of strings and a vector of strings,
    //returning a vector of the results
    return make_append(deps, rhs::dependencies());
    //return the dependencies for this class along with any dependencies the parent classes have
}

//!This function processes the input data for the class
/*!
 * This function process data that has been generated from an input file
 * All of the names returned in the dependencies file exist and have been processed
 * @param dat The data structure containing the variable names
 */
void example_rhs::postprocess(input& dat){
    //warning so nobody ever actually constructs one of these
    //please don't have warnings that always throw in your code ever
    err("Example_rhs created, mysteriously fails on some architectures seemingly from AVX instructions. Also does nothing useful", "rhs::create", "rhs/rhs.cpp", WARNING);

    rhs::postprocess(dat);//always postprocess the parent class first
    //Any class that inherits from rhs, or item_dim, has access to a variable
    //called dimension, which represents the dimension of the problem at hand
    //dimension is initialized in the item_dep, which is why the parent postprocessing
    //must be called first
    

    //There are two types of errors, warning and fatal errors.
    //Warnings print something to the screen to warn the user,
    //but do not stop the program. A fatal error causes an immediate exit
    //These are caused by the function err. Details can be found in the documentation
    //
    //Lets throw a fatal error if the dimension is less than 10
    if(dimension > 10){
        err("dimension > 10, which is required for example_rhs", "example_rhs::postprocess",
                "rhs/example_rhs.cpp", FATAL_ERROR);
    }
    
    //Lets also throw an error if the dimension is equal to 5
    if(dimension==5){
        err("dimension=5, which gets a warning from example_rhs", "example_rhs::postprocess",
                "rhs/example_rhs.cpp", WARNING);
    }
    //The map contains item*, which point to an item of arbitrary type.
    //I plan to implement type checking of the dependencies.
    //But sometime later
    //
    //To retrieve a value, you can call the retrieve function of the variable at hand
    //you just pass the name of parameter that is being retrieved to the map, which returns
    //an item*. Then this pointer is called to retrieve the value, and is passed the address
    //of val1

    dat.retrieve(val1, "val1", this);
    dat.retrieve(val2, "val2", this);
    dat.retrieve(random_info, "random_info", this);

    //now, we are going to allocate some memory to something
    //This may be useful for storing temporary calculations during the RHS.
    //Inheriting from item_dim provides access to a memory pool, memp.
    //just pass memp.create the dimension of the problem and the
    //addresses of the pointers and the alignment, allocation, and deallocation is
    //taken care of!
    //
    //memp.create(dimension, &ptr1, &ptr2, &ptr3, etc)
    //for custom alignment(standard is 32 byte)
    //memp.create(alignment, dimension, &ptr1, &prt2, etc)
    //
    memp.create(32, dimension, &value_holder);
    for(size_t i = 0; i < dimension; i++){
        value_holder[i] = Id*(double)i*val1 + (dimension-i)*val2;
    }
}
//!Update function for example_rhs
/*!
 * This function updates the values held by this class in case one of them
 * is a variable that is being controlled. We don't re-allocate memory,
 * since integers like dimension cannot be controlled. However, we must re-evaluate
 * the array value_holder. The values val1, val2, etc are controlled themselves
 * by the controller
 * 
 */
void example_rhs::update(){
    for(size_t i = 0; i < dimension; i++){
        value_holder[i] = Id*(double)i*val1 + (dimension-i)*val2;
    }
}

//!Destructor for example_rhs
/*!
 * This is the destructor for example_rhs. Whenever the class is deleted by the engine, this
 * frees resources. The destructor for the memory pool is called here, so don't worry
 * about deallocations
 */
example_rhs::~example_rhs(){
    //this->free_various_resources();
    //Seriously though try to use resources which destruct automatically
    //Through their own destructor
}

//!This is the rhs function that calculates the derivative
/*!
 * This functions calculates a derivative given the current time and function value
 * it takes in a pointer to x, which is the function value, and stores that data in dx.
 *
 * @param x The current function value
 * @param dx The array i which the derivative is stored
 * @param t The current time
 */
int example_rhs::dxdt(const ptr_passer _x, ptr_passer _dx, double t){
    //This pulls the pointer out of the ptr_passer object.
    //The conversion is implicit, but will throw a fatal error if the
    //types do not match.
    const comp* restr x = _x;
    comp* restr dx = _dx;
    //This doesn't really do anything useful
    //but showcases some of the compex variable functions
    //See detailed documentation for more
    //But is a nice example
    for(size_t i = 0; i < dimension; i++){
        dx[i] = (_real(x[i])*_imag(value_holder[i])) * (val2/val1);
    }
    return 0;
}
