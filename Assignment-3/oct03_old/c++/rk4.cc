// This code can use either valarray or vector as you wish.
// It also follows the 2011 C++ standard, not the 1998 one.

#include <valarray>
using real_t = double;  // C++11, not C++98
using array_t = std::valarray<real_t>;  // change to vector if you wish
using vecfield_t = void (*)(real_t t, const array_t& y, array_t& dy);

// --------------------------------------------------------------------------
// The 4th-order RK integrator.  We integrate from TSTART to TEND in NSTEPS
// steps with no local error estimation or time-step adjustment.  The initial
// condition is contained in Y and the approximated solution is stored there
// on return.
void
rk4(vecfield_t f, array_t& y, real_t tstart, real_t tend, size_t nsteps)
{
    real_t h = (tend - tstart)/nsteps;
    size_t n = y.size();
    if(n == 0) return;

    array_t k1(n), k2(n), k3(n), k4(n), ytemp(n);
    real_t dt = (tend - tstart)/nsteps;

    // Your code here to run the *body* of RKSTEP in a large loop.
    // For efficiency reasons, you do not want to create and destroy
    // small arrays repeatedly, because doing so with C++ STL containers
    // is very expensive.
    return;
}
