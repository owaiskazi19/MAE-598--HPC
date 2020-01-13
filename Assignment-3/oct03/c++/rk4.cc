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
rkstep(vecfield_t f, real_t t, real_t h, array_t& y)
{
   size_t n = y.size();  // number of equations
   if(n == 0) return;

   array_t dy1(n), dy2(n), dy3(n), dy4(n), ytemp(n);
   f(t, y, dy1);                  // Get k1
   for(size_t j = 0; j < n; j++)  // Update yn to calculate k2
      ytemp[j] = y[j] + (h*0.5)*dy1[j];
   
   f(t+(0.5 * h), ytemp, dy2);    // Get k2
   for(size_t j=0; j < n; j++)    // Update yn to calculate k3
      ytemp[j] = y[j] + (h*0.5)*dy2[j];
   
   f(t+(0.5 * h), ytemp, dy3);    // Get k3
   for(size_t j=0;j < n; j++)     // Update yn to calculate k4
      ytemp[j] = y[j] + h*dy3[j];
   
   f(t+h, ytemp, dy4);            // Get k4
   for(size_t j=0;j < n; j++)     // Update yn to calculate yn+1
      y[j] = y[j] + (h/6)*(dy1[j] + 2*dy2[j] + 2*dy3[j] + dy4[j]);

   return;
}

void
rk4(vecfield_t f, array_t& y, real_t tstart, real_t tend, size_t nsteps)
{
    real_t h = (tend - tstart)/nsteps;
    size_t n = y.size();
    if(n == 0) return;
    real_t  t =0;               // Iniatialize time
    array_t k1(n), k2(n), k3(n), k4(n), ytemp(n);
    real_t dt = (tend - tstart)/nsteps;          // Calculate height
    for (size_t i =0 ;i<nsteps;i++){
     rkstep(f,t,h,y);
     t += h;
  }

    // Your code here to run the *body* of RKSTEP in a large loop.
    // For efficiency reasons, you do not want to create and destroy
    // small arrays repeatedly, because doing so with C++ STL containers
    // is very expensive.
    return;
}
