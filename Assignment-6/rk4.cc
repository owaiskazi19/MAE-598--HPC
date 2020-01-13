// This code can use either valarray or vector as you wish.
// It also follows the 2011 C++ standard, not the 1998 one.

#include <valarray>
using real_t = double;                 // C++11, not C++98
using array_t = std::valarray<real_t>; // change to vector if you wish
using vecfield_t = void (*)(real_t t, const array_t &y, array_t &dy, size_t num_eqs);

// --------------------------------------------------------------------------
// The 4th-order RK integrator.  We integrate from TSTART to TEND in NSTEPS
// steps with no local error estimation or time-step adjustment.  The initial
// condition is contained in Y and the approximated solution is stored there
// on return.
void rk4(vecfield_t f, array_t &y, real_t tstart, real_t tend, size_t nsteps, size_t num_eqs) {
  size_t n = y.size();
  if (n == 0)
    return;

  array_t k1(num_eqs), k2(num_eqs), k3(num_eqs), k4(num_eqs), ytemp(num_eqs);
  real_t dt = (tend - tstart) / nsteps;

  // Your code here to run the *body* of RKSTEP in a large loop.
  // For efficiency reasons, you do not want to create and destroy
  // small arrays repeatedly, because doing so with C++ STL containers
  // is very expensive.

  for (size_t i = 0; i < nsteps; i++) {
    real_t t = tstart + i * dt;

    // Get k1 for every equation in y
    f(t, y, k1, num_eqs);
    // Update every value in preparation to calculate k2
    for (size_t j = 0; j < num_eqs; j++)
      ytemp[j] = y[j] + (dt/2) * k1[j];

    // Get k2 for every equation in ytemp
    f(t + (dt/2), ytemp, k2, num_eqs);
    // Update every value in preparation to calculate k3
    for (size_t j = 0; j < num_eqs; j++)
      ytemp[j] = y[j] + (dt/2) * k2[j];

    // Get k3 for every equation in ytemp
    f(t + (dt/2), ytemp, k3, num_eqs);
    // Update every value in preparation to calculate k4
    for (size_t j = 0; j < num_eqs; j++)
      ytemp[j] = y[j] + dt * k3[j];

    // Get k4 for every equation in ytemp
    f(t + dt, ytemp, k4, num_eqs);
    // Finally calculate the value for y(n+1)
    for (size_t j = 0; j < num_eqs; j++)
      y[j] = y[j] + (dt / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
  }

  return;
}
