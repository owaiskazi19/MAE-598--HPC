// This is an example of a header file for a C program that implements
// the Sept. 12 assignment.  It is written for C99 or later language
// standards.  Supply option -std=c99 to either the Intel
// C compiler (icc) or GNU C (gcc).
//
#include <stdlib.h>
typedef double real_t;

// Calling interface for the vector field that defines
// the ordinary differential equation.  N is the number of equations
// in the system.
typedef void (*vecfield_t)(real_t, size_t n, const real_t y[n], real_t dy[n]);

// Calling interfaces for the fourth-order RK method
extern void rkstep(vecfield_t f, real_t t, real_t h, size_t n, real_t y[n]);
extern void rk4(vecfield_t f, size_t n, real_t y[n], real_t tstart, 
  real_t tend, size_t nsteps);
