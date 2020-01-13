#ifndef RK4_H
#define RK4_H
// This is an example of a header file for a C++ program to implement
// a fourth-order Runge-Kutta method.  It is written for the C++11
// or later language standards.  Supply option -std=c++11 to either the Intel
// C++ compiler (icpc) or GNU C++ (g++).
//
#include <valarray> // or, less optimally, <vector>

// Some convenient typedefs
using real_t = double;
using array_t = std::valarray<real_t>;

// Calling interface for the vector field that defines
// the ordinary differential equation.  N is the number of equations
// in the system.
using vecfield_t = void (*)(real_t, const array_t &y, array_t &dy, size_t num_eqs);

// Calling interfaces for the fourth-order RK method
extern void rk4(vecfield_t f, array_t &y, real_t tstart, real_t tend, size_t, size_t num_eqs);
#endif
