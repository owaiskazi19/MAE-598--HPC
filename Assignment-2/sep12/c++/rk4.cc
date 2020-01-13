#include "rk4.h"

// RKSTEP - take on Runge-Kutta step of size H of the vector field F(T,Y).
void
rkstep(vecfield_t f, real_t t, real_t h, array_t& y)
{
   size_t n = y.size();  // number of equations
   if(n == 0) return;

   array_t dy1(n), ytemp(n);
   f(t, y, dy1);
   for(size_t j = 0; j < n; j++)  // first partial RK step
      ytemp[j] = y[j] + (h*0.5)*dy1[j];

   // your code here

   return;
}

// RK4 - integrate the vector field F from TSTART to TEND in NSTEPS equal
// time steps.
void
rk4(vecfield_t f, array_t& y, real_t tstart, real_t tend, size_t nsteps)
{
  // your code here
}
