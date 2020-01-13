#ifndef PENDULUM_H
#define PENDULUM_H
#include <cmath>
#include "rk4.h"

const real_t PI = 3.1415926535897932;
const real_t MINUSPI = -3.1415926535897932;
const real_t TWOPI = 2*PI;
extern void poincare(array_t& y, size_t niter);
#endif
