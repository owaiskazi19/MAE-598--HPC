#include "rk4.h"

// You need to include the following files for printf() and for exp(),
// respectively.
#include <stdio.h>
#include <math.h> 

// Skeleton of a sample  ODE.C file

// MYODE is the name of your actual vector field.  It must have the same
// arguments in the same order as the prototype in rk4.h.

void myode(real_t t, size_t neq, const real_t y[neq], real_t dy[neq])
{
  // You may assume that NEQ == 1 for the test problem, but in general
  // NEQ refers to the Number of EQuations in the system of ODEs
}

int main(void)
{
   real_t tstart = 0.0;
   real_t tend = 2.0;
   size_t nsteps = 16;  // so that h=1/8
   real_t y[1];

   y[0] = 0.0;
   // other variables as needed
   // and call RK4

   printf("RK4 result: %g\n", y[0]);
   printf("truth: %g\n", exp(1.0));
}

