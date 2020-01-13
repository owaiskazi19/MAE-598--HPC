#include "rk4.h"
#include "pendulum.h"
#include <fstream>
#include<iostream>
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
   grid(size_t n=400, real_t left=-PI, real_t right=PI, real_t bottom=-PI,
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
   real_t pt[3][2];  // the list of periodic points
   real_t epsilon;   // convergence criterion
   size_t period[3];  // their respective periods
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
     pt{{-2.407661292643E+00, 6.773267621896E-01},  // period 1
       {-6.099150484926E-01, 1.677463819037E+00},  // period 2
       {-7.218820695570E-01, -4.620944521247E-01}}, // period 2
     epsilon{1.0e-09},  period{1, 2, 2}, maxiter{400}, nfixed{3}, maxperiod{2}
     {};
   std::valarray<int> compute_basin(const grid&);
   int check_point(const array_t&, double_t, double_t);
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
fixedpoint::check_point(const array_t& y, real_t y1, real_t y2) // checks 2 value 
{
   using namespace std;  // for simplicity
   valarray<real_t> y_check(2);

   // your code here
   for (size_t k = 0; k < nfixed; k++) {
     real_t point_distance = hypot(fabs(pt[k][0] - y1), fabs(pt[k][1]-y2)); //hypot returns the square root of sum of square of arguments passed. Calculate point
     // If point is near B then make B as B' and call poincare again.
      if (point_distance >= epsilon && period[k] == 2) {
      y_check[0] = pt[k][0];
      y_check[1] = pt[k][1];
      poincare(y_check, 1);
      point_distance = hypot(fabs(y1 - y_check[0]), fabs(y2 - y_check[1]));
    }

     if (point_distance < epsilon) {         //Checking whether point is in the fixed points
       return (k+1);                         //Returns FIXEDPOINT 
     }
   }
   return(0);  // replace this with the appropriate index
}

//----------------------------------------------------------------------------
// Compute the basin boundary of each fixed point.  For the (J,K)th grid point,
// compute up to MAXITER iterations of the Poincare map and see where the
// initial condition ends up.  If it lies within EPSILON of one of the fixed
// points (which CHECK_POINT will determine), then mark the basin accordingly.

std::valarray<int> 
fixedpoint::compute_basin(const grid& box)
{
    size_t n = box.resolution;
    std::valarray<int> basin(n*n);

    // your code here
    array_t y(2*n);
    array_t y2(2*n);
    array_t arr(n);
    int num;
    int in = 0;

    for (size_t j = 0;j<n; j++) {
      real_t xx = (2*j + 1)* PI/n - PI;           //Calculates angular position of the oscillator
      //Computes up to MAXITER iterations of poincare
      for (size_t k =0; k<n ;k++) {
        real_t yy = (2*k + 1) * PI /n - PI;      
        basin[j*(n-1)+k] = 0;
        y[2*k] = xx;                                //Assigns it to array y
        y[2*k + 1] = yy;
        arr[k] = k; // Create array to store the values which are not converged
      }
      int points_remaining = n;
      for (size_t m =0;  points_remaining >0 && m < maxiter; m+=maxiter/10) {
        in = 0;
        int count = 0;
        poincare(y,maxiter/10);                 // Poincare for maxiter/10
        for(size_t l =0; l<points_remaining ; l++) {
          std::cout<< "Poincare" << (j*n+l) << "->" << xx << "->" <<y[l]<<","<<y[l+1]<<"\n";
        
          num =  check_point(y, y[2*l], y[2*l+1]); //Pass the points to see whether they are converged
          if (num != 0) {                          
            basin[j*(n-1) + arr[l]] = num;               // Store converged value to basin
            count++;
          } 
          else 
          {
            arr[in] = arr[l];     // Update mapping for arr to map from in to what l was pointing
            y2[2*in] = y[2*l];     //Store the udpated value in new y2
            y2[2*in+1] = y[2*l+1];
            in++;
          } 
        }
        
        points_remaining = points_remaining - count;  //Decrement according to points remaining to converged    
        y.swap(y2);                      //Swap new y with old y        
      }
    }  
    return(basin);
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

    valarray<int> basin(fp.compute_basin(box));

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
