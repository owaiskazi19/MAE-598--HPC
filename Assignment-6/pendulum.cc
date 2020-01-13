#include "pendulum.h"
#include "rk4.h"

// This code defines the Poincare' map for the forced damped nonlinear
// pendulum, x'' + damping*x' + sin(x) = force*cos(t).
// Of course, you may use the same sampled implementation as in the C
// code, because C++ is largely upwardly compatible with C.
// This version illustrates the use of a "static class" to structure the
// map and its parameters.

// There is only ever one instance of a static class.
// Static classes do not have constructors; instead, their
// data elements are initialized outside of the class as shown below.

// Any class (not just static classes) can declare a STATIC member function.
// The difference between static member functions and ordinary ones is that
// the static member functions do not have a THIS pointer passed implicitly
// among their arguments.  Consequently, static member functions cannot
// access private data elements unless they are also STATIC.

// The restrictions nevertheless work for this application because
// (1) we integrate exactly one forced, damped pendulum and
// (2) the calling interface to the RK4 procedure does not admit a function
//     containing a pointer to an arbitrary data structure.

// One potential advantage of this arrangement over the C version is that
// the parameters DAMPING and FORCE are in their own namespace.  Outside
// of the member functions, they can be accessed only with the qualified
// names FORCED_DAMPED::DAMPING and FORCED_DAMPED::FORCE and so cannot be
// confused with other variables named DAMPING and FORCE that may be defined
// elsewhere in the program.

inline real_t wrap_around(real_t angle) {
  real_t m = (angle + PI) - TWOPI * floor((angle + PI) / TWOPI);
  return fmod(TWOPI + m, TWOPI) - PI;
}

//  The first-order system that defines the forced, damped pendulum
void forced_damped::pendulum(real_t t, const array_t &y, array_t &dy, size_t num_eqs) {
  for (size_t i = 0; i < num_eqs/2; i++) {
    dy[2 * i] = y[2 * i + 1]; // derivative of position is velocity
    dy[2 * i + 1] = force * std::cos(t) - std::sin(y[2 * i]) - damping * y[2 * i + 1];
  }
  return; // optional
};

// Member function that defines the Poincare' map
void forced_damped::poincare(array_t &y, size_t niter, size_t points_left) {
  const size_t NSTEPS = 256;
  for (size_t j = 0; j < niter; j++) {
    rk4(pendulum, y, 0., TWOPI, NSTEPS, points_left * 2);
    // Your code here to wrap the position, y[0], back into
    // the interval [-PI,PI].
    // https://stackoverflow.com/questions/11980292/how-to-wrap-around-a-range
    for (size_t i = 0; i < points_left; i++)
      y[2 * i] = wrap_around(y[2 * i]);
  }
  return; // optional
};

real_t forced_damped::damping = 0.2;
real_t forced_damped::force = 1.66;
