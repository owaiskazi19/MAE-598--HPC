#ifndef PENDULUM_H
#define PENDULUM_H
#include "rk4.h"
#include <cmath>

const real_t PI = 3.1415926535897932;
const real_t TWOPI = 2 * PI;

class forced_damped {
public:
  static real_t damping; // map parameters, specified below.
  static real_t force;

  //  The first-order system that defines the forced, damped pendulum
  static void pendulum(real_t t, const array_t &y, array_t &dy, size_t num_eqs);

  // Member function that defines the Poincare' map
  static void poincare(array_t &y, size_t niter, size_t points_left);
};

#endif
