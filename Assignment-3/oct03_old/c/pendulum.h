#ifndef PENDULUM_H
#define PENDULUM_H
#include <stdlib.h>
#include <tgmath.h>

// These must be #define's, not const real_t, for initialization purposes 
// in C99.
#define PI (3.1415926535897932)
#define MINUSPI (-3.1415926535897932)
#define TWOPI (2*PI)

// The Poincare' map for the forced damped pendulum.
extern void poincare(real_t y[2], size_t niter);
#endif
