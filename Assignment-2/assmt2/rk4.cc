#include "rk4.h"

// RKSTEP - take on Runge-Kutta step of size H of the vector field F(T,Y).
void
rkstep(vecfield_t f, real_t t, real_t h, array_t& y)
{
   size_t n = y.size();  // number of equations
   if(n == 0) return;

   array_t dy1(n), dy2(n), dy3(n), dy4(n), ytemp(n);
   f(t, y, dy1);                  // Get k1
   for(size_t j = 0; j < n; j++)  // Update yn to calculate k2
      ytemp[j] = y[j] + (h*0.5)*dy1[j];
   
   f(t+(0.5 * h), ytemp, dy2);    // Get k2
   for(size_t j=0; j < n; j++)    // Update yn to calculate k3
      ytemp[j] = y[j] + (h*0.5)*dy2[j];
   
   f(t+(0.5 * h), ytemp, dy3);    // Get k3
   for(size_t j=0;j < n; j++)     // Update yn to calculate k4
      ytemp[j] = y[j] + h*dy3[j];
   
   f(t+h, ytemp, dy4);            // Get k4
   for(size_t j=0;j < n; j++)     // Update yn to calculate yn+1
      y[j] = y[j] + (h/6)*(dy1[j] + 2*dy2[j] + 2*dy3[j] + dy4[j]);

   return;
}

// RK4 - integrate the vector field F from TSTART to TEND in NSTEPS equal
// time steps.
void
rk4(vecfield_t f, array_t& y, real_t tstart, real_t tend, size_t nsteps)
{
  real_t  t =0;  // Iniatlize time
  real_t h = (tend - tstart)/nsteps;  // Calculate height
  for (size_t i =0 ;i<nsteps;i++){
     rkstep(f,t,h,y);
     t += h;
  }
}
