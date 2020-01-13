// Collaborated with Amitabh Das
#include "rk4.h"
#include "pendulum.h"
#include <fstream>
#include<iostream>
#include <vector>
#include <omp.h>
#include <chrono>
//  This is one possible way to represent the grid to be used for computing
//  the basin boundary.  We assume a square grid here but it's straightforward
//  to extend the structure to accommodate a rectangular one.
//  By default, the constructor initializs the grid to [-PI,PI] x [-PI,PI]
//  with a resolution of 400x400.
//  A declaration like 
//    grid box(50);
//  yields a box 50x50 box over the domain [-PI,PI] x [-PI,PI].

struct grid {
   real_t xmin;  // angular position limits
   real_t xmax;
   real_t ymin;  // angular velocity limits
   real_t ymax;
   size_t resolution;  // number of grid points in each direction
   // constructor with default arguments
   grid(size_t n=100, real_t left=-PI, real_t right=PI, real_t bottom=-PI,
     real_t top=PI) : xmin{left}, xmax{right}, ymin{bottom}, ymax{top},
       resolution{n} {};
};

//  Here is one possible way to collect all the information relevant to
//  tracking the periodic orbits.  In this format, only one representative
//  of each periodic orbit is stored in PT; you have to iterate the map
//  to check whether a candidate point comes within EPSILON of the stored orbit.
//  Since it is unlikely that high-period orbits can be found reliably with
//  the methods used here, and since we aren't interested in vast numbers
//  of different periodic points, the integer-valued components could be
//  declared INT; however, the SIZE_T declaration emphasizes that the values
//  are supposed to be nonnegative.

class fixedpoint {
   std::vector<real_t> pt;  // the list of periodic points
   real_t epsilon;   // convergence criterion
   std::vector<size_t> period;  // their respective periods
   size_t maxiter;  // maximum number of map iterations
   size_t nfixed;  // the number of fixed points in the list (i.e., 3)
   size_t maxperiod;  // of any orbit under consideration by this program
   fixedpoint(const fixedpoint&);  // no public copy constructor
   fixedpoint(const fixedpoint&&);  // no public move constructor
   fixedpoint& operator=(const fixedpoint&); // no public copy assignment
   fixedpoint& operator=(const fixedpoint&&); // no public move assignment
public:
   // default constructor
   fixedpoint() : 
     pt(), // period 2
     epsilon{1.0e-09},  period(), maxiter{400}, nfixed{0}, maxperiod{4}
     {};
   std::valarray<int> compute_basin(const grid&);
   int check_point(const array_t&, real_t, real_t, real_t, real_t);
};

// Include the following statement if you want to use the static class defined
// in pendulum.cc.  You can call the Poincare map defined there as 
//    forced_damped::poincare(y, niter);
// for appropriately defined variables Y and NITER.

//class forced_damped;

//----------------------------------------------------------------------------
// Check whether Y or its orbit contains a point that is already in the
// list of fixed points.  If so, then return an integer value (1, 2, or 3)
// indicating whether it is the first, second, or third point in FIXEDPOINT.
// In other words, if Y lies within fp.epsilon of pt[i], then return i+1,
// i = 0, 1, 2.  Otherwise, return 0.
// You can use any metric here.  The 2-norm in C/C++ can be computed with the
// standard library function HYPOT(), but the 0-norm also works as the
// maximum of the absolute values of the differences between any component.
// C99 and C++ include the functions FMAX() and FABS(), which are convenient
// for this purpose.

int
fixedpoint::check_point(const array_t& y, real_t y1, real_t y2, real_t damping, real_t force) // checks 2 value 
{
   using namespace std;  // for simplicity
   valarray<real_t> y_period(2);
   real_t point_distance = epsilon;

   size_t p_period = 0; // Step 3 - Find prime period
   y_period[0] = y1;
   y_period[1] = y2;

   while (p_period < maxperiod) {
     poincare(y_period, 1, damping, force);
     point_distance = hypot(fabs(y1 - y_period[0]), fabs(y2 - y_period[1]));
     if (point_distance <= epsilon/2)
        break;
     p_period++;
   }

    if (p_period == maxperiod)  // aperiodic within maxperiod, return 0
        return 0;

  int conv_pos = 0;
  #pragma omp critical
  {
    // If there are points in L, find the closest one
    valarray<real_t> y_check(2);

    for (size_t k = 0; k < nfixed && conv_pos == 0; k++) {
      if (period[k] != p_period + 1) continue;
      y_check[0] = pt[2*k];
      y_check[1] = pt[2*k + 1];
      for (size_t l =0; l<period[k]; l++) {
        point_distance = hypot(fabs(y1 - y_check[0]), fabs(y2 - y_check[1])); //hypot returns the square root of sum of square of arguments passed. Calculate point
        if (point_distance < epsilon) {         
            conv_pos = k + 1;              
            break;
        }
        poincare(y_check, 1, damping, force);
      }
    }
    
    // If given point not close to any point in L, add given point and its period to L

    if (conv_pos == 0) {
      pt.push_back(y_period[0]);
      pt.push_back(y_period[1]);
      period.push_back(p_period + 1);
      nfixed++;
      conv_pos = nfixed;
    }
  }

  return conv_pos;
  
}

//----------------------------------------------------------------------------
// Compute the basin boundary of each fixed point.  For the (J,K)th grid point,
// compute up to MAXITER iterations of the Poincare map and see where the
// initial condition ends up.  If it lies within EPSILON of one of the fixed
// points (which CHECK_POINT will determine), then mark the basin accordingly.

std::valarray<int> 
fixedpoint::compute_basin(const grid& box)
{
  real_t damping, force;
   std::cout<<"Enter damping value"<< std::endl;
   std::cin>>damping;
   std::cout<<"Enter force value"<<std::endl;
   std::cin>>force;

  size_t n = box.resolution;
  std::valarray<int> basin(n*n);

  // your code here
  
  real_t x_height = (box.xmax - box.xmin) / n;  //Calculate x_height
  real_t y_height = (box.ymax - box.ymin) / n;  //Calculate y_height
  

  #pragma omp parallel for
  for (size_t j = 0;j<n; j++) {
    array_t y(2*n);
    real_t xx = box.xmin + x_height * j;
    for (size_t k =0; k<n ;k++) {
      real_t yy = box.ymin + y_height * k;
      basin[j*(n-1)+k] = 0;
      y[2 * k] = xx;
      y[2 * k + 1] = yy;
    }  
    poincare(y,maxiter, damping, force); 
    for (size_t k =0;k<n;k++) {
        basin[j*(n-1) + k] = check_point(y, y[2*k], y[2*k + 1], damping, force); 
    }
   }

    std::cout << "Points found " << nfixed << std::endl;
    for (size_t l = 0; l < nfixed; l++) {
        std::cout << pt[2*l] << ", " << pt[2*l + 1] << " - " << period[l] << std::endl;
    } 
    return (basin);
}

//----------------------------------------------------------------------------
// Main program.  

int
main(int argc, char **argv)
{
    using namespace std;

    fixedpoint fp; //  Our only instance
    grid box;  // will have the default limits and resolution

    // Replace the contents of any previous file and output as binary
    ofstream outfile("basin.dat", ios::trunc | ios::binary);

    auto start = std::chrono::high_resolution_clock::now();
    valarray<int> basin(fp.compute_basin(box));
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time Span "<<duration.count() << endl;

    // Write out the results in a binary format compatible with the MATLAB
    // script wada.m for visualization.

    outfile.write(reinterpret_cast<const char*>(&box.xmin), sizeof(box.xmin));
    outfile.write(reinterpret_cast<const char*>(&box.xmax), sizeof(box.xmax));
    outfile.write(reinterpret_cast<const char*>(&box.ymin), sizeof(box.ymin));
    outfile.write(reinterpret_cast<const char*>(&box.ymax), sizeof(box.ymax));
    outfile.write(reinterpret_cast<const char*>(&box.resolution),
      sizeof(box.resolution));
    outfile.write(reinterpret_cast<const char*>(&basin[0]),
      sizeof(int)*basin.size());  // N x N basin boundary
    outfile.close();
    return(EXIT_SUCCESS);
}
