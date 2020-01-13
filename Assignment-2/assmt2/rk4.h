// RK4.H - sample header file for C++ Runge-Kutta methods.
// You may use <valarray> or <vector> as you wish; make the appropriate
// changes below for the latter.
// The code here follows the 2011 C++ standard, not the 1998 one.
// Supply option -std=c++11 to the GNU and Intel compilers.

#include <valarray>
using real_t = double;  // C++11, not C++98
using array_t = std::valarray<real_t>;  // change to vector<> if you wish

// Calling interface for the vector field that defines
// the ordinary differential equation.  Y.SIZE() defines the number of equations
// in the system.  Ensure that DY.SIZE() equals Y.SIZE() for the derivatives.
using vecfield_t = void (*)(real_t t, const array_t& y, array_t& dy);

// Calling interfaces for the fourth-order RK method
extern void rkstep(vecfield_t f, real_t t, real_t h, array_t& y);
extern void rk4(vecfield_t f, array_t& y, real_t tstart, real_t tend,
  size_t nsteps);
