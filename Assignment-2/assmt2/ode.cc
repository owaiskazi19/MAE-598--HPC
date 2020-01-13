#include "rk4.h"

// You need to include the following files for formatted i/o and for exp(),
// respectively.
#include <iostream>
#include <cmath> 

// Skeleton of a sample  ODE.CC  file.  Important: compile with the
// flag -std=c++11 with both g++ and Intel icpc to use more recent
// versions of the C++ language.

// MYODE is the name of your actual vector field.  It must have the same
// arguments in the same order as the prototype in rk4.h.

void myode(real_t t, const array_t& y, array_t& dy)
{
  // The number of equations is assumed to equal y.size().
  // For this test problem, of course, there is only one equation.
  dy[0] = y[0] / 2;
}

int main(void)
{
   real_t tstart = 0.0;
   real_t tend = 2.0;
   size_t nsteps = 16;  // so that h=1/8
   array_t y(1);
   using namespace std;  // otherwise we have to write std::cout, std::exp, etc.

   y[0] = 1.0;   // just in case
   // other variables as needed
   // and call RK4
   rk4( myode, y, tstart, tend, nsteps);

   cout << "RK4 result: " << y[0] << endl;
   cout << "truth: " << exp(1.0) << endl;
   return(EXIT_SUCCESS);
}
