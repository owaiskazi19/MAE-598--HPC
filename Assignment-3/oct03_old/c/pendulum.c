#include "rk4.h"
#include "pendulum.h"

// This code defines the Poincare' map for the forced damped nonlinear
// pendulum, x'' + damping*x' + sin(x) = force*cos(t).
// Entities declared STATIC are not exportable beyond this file.

// Damping and forcing parameters for the nonlinear pendulum considered
// for this exercise.  If you wanted to include an input routine in
// another file to set their values, then you would have to do one of
// two things: (1) remove the STATIC qualifier and document how and
// where the values might be changed; or (2) add an input/output routine
// to this file that sets their values.

static real_t damping = 0.2;
static real_t force = 1.66;

//----------------------------------------------------------------------------
//  The first-order system corresponding to the forced damped pendulum.
//  You may assume that N = 2 for this assignment.

static void
pendulum(real_t t, size_t n, const real_t y[n], real_t dy[n])
{
   dy[0] = y[1];   // derivative of position is velocity
   dy[1] = force*cos(t) - sin(y[0]) - damping*y[1];
   return;
}

//----------------------------------------------------------------------------
// Poincare map for the forced nonlinear pendulum.  
// Y is the 2-vector giving the angular (position, velocity) of the
// pendulum after every 2*PI time interval.
void
poincare(real_t y[2], size_t niter)
{
   const size_t NSTEPS = 256;

   for(size_t j = 0; j < niter; j++)  {
      rk4(pendulum, 2, y, 0.0, TWOPI, NSTEPS);
      // Compute the appropriate modulus operation to move y[0] (position)
      // into the interval [-PI, PI].
      // your code here:  y[0] = ...
   }
   return;
}

