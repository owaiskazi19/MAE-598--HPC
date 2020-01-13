#include "rk4.h"

// RKSTEP - take on Runge-Kutta step of size H of the vector field F(T,Y).
// This file illustrates some constructs that were introduced in the 1999 C
// standard and are not found in Kernighan and Ritchie's definitive book.
// You may find some of them helpful for scientific code.
// Supply option -std=c99 to the Intel and GNU compilers.
void
rkstep(vecfield_t f, real_t t, real_t h, size_t n, real_t y[n])
{
   real_t dy1[n], ytemp[n];  // C99+, not C89
   //  include other arrays as needed
 
   if(n == 0) return;

   // First Euler step
   f(t, n, y, dy1);
   for(size_t j = 0; j < n; j++)  // C99+, not C89
      ytemp[j] = y[j] + (h*0.5)*dy1[j];

   // your code here
   return;
}

// RK4 - integrate the vector field F from TSTART to TEND in NSTEPS equal
// time steps.
void
rk4(vecfield_t f, size_t n, real_t y[n], real_t tstart, real_t tend, 
  size_t nsteps)
{
  // your code here
}
