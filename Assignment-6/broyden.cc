// Collaborated with Amitabh Das
#include "pendulum.h"
#include "rk4.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <chrono>
#include <vector>
#include <functional>

using broyden_f = std::function<void(std::valarray<real_t> &)>;

//  This is one possible way to represent the grid to be used for computing
//  the basin boundary.  We assume a square grid here but it's straightforward
//  to extend the structure to accommodate a rectangular one.
//  By default, the constructor initializs the grid to [-PI,PI] x [-PI,PI]
//  with a resolution of 400x400.
//  A declaration like
//    grid box(50);
//  yields a box 50x50 box over the domain [-PI,PI] x [-PI,PI].

struct grid {
  real_t xmin; // angular position limits
  real_t xmax;
  real_t ymin; // angular velocity limits
  real_t ymax;
  size_t resolution; // number of grid points in each direction
  // constructor with default arguments
  grid(size_t n = 100, real_t left = -PI, real_t right = PI,
       real_t bottom = -PI, real_t top = PI)
      : xmin{left}, xmax{right}, ymin{bottom}, ymax{top}, resolution{n} {};
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
  std::vector<real_t> pt; // the list of periodic points
  real_t epsilon1;   // convergence criterion
  real_t epsilon2;   // convergence criterion
  real_t epsilon3;   // convergence criterion
  std::vector<size_t> period; // their respective periods
  std::vector<bool> stable; // their respective periods
  size_t maxiter;   // maximum number of map iterations
  size_t nfixed;    // the number of fixed points in the list (i.e., 3)
  size_t maxperiod; // of any orbit under consideration by this program
  fixedpoint(const fixedpoint &);             // no public copy constructor
  fixedpoint(const fixedpoint &&);            // no public move constructor
  fixedpoint &operator=(const fixedpoint &);  // no public copy assignment
  fixedpoint &operator=(const fixedpoint &&); // no public move assignment
  size_t not_conv;
  size_t exists;

  int closest_point_in_list(real_t y1, real_t y2, size_t prime_period_estimate, real_t epsilon);
  bool broyden(std::valarray<real_t> &, size_t maxiter, const broyden_f &);
  size_t estimate_prime_period(real_t y1, real_t y2);
  std::valarray<real_t> initial_jacobian(const std::valarray<real_t> &, const broyden_f &, real_t);

public:
  // default constructor
  fixedpoint() : pt(), epsilon1{0.1}, epsilon2{1.0e-09}, epsilon3{1.0e-12}, period(), maxiter{20}, nfixed{0}, maxperiod{6}, not_conv{0}, exists{0} {};

  std::valarray<int> compute_basin(const grid &);
  int check_point(const array_t &, size_t pos);
};

// Include the following statement if you want to use the static class defined
// in pendulum.cc.  You can call the Poincare map defined there as
//    forced_damped::poincare(y, niter);
// for appropriately defined variables Y and NITER.

class forced_damped;


size_t fixedpoint::estimate_prime_period(real_t y1, real_t y2) {
  std::valarray<real_t> y_period = {y1, y2};

  size_t prime_period_estimate = 0;
  real_t dist = 10, temp = 10;
  for (size_t i = 0; i < maxperiod; i++) {
    forced_damped::poincare(y_period, 1, 1);
    temp = hypot(fabs(y1 - y_period[0]), fabs(y2 - y_period[1]));
    if (temp < dist) {
      prime_period_estimate = i + 1;
      dist = temp;
      if (temp < epsilon1) break;
    }
  }

  return prime_period_estimate;
}


bool fixedpoint::broyden(std::valarray<real_t> &estimate, size_t maxiter, const broyden_f &f) {
  size_t n = 0;
  std::valarray<real_t> x_n(estimate);
  std::valarray<real_t> f_n(2);
  std::valarray<real_t> jacobian(initial_jacobian(estimate, f, TWOPI/100));

  real_t determinant = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
  std::swap(jacobian[0], jacobian[3]);
  jacobian[1] *= -1;
  jacobian[2] *= -1;
  jacobian /= 1/determinant;

  for (;n < maxiter; n++) {
    std::valarray<real_t> deltax = {0, 0}, deltaf = {0, 0};

    f_n = x_n;
    f(f_n);

    deltax[0] = -jacobian[0] * f_n[0] - jacobian[1] * f_n[1];
    deltax[1] = -jacobian[2] * f_n[0] - jacobian[3] * f_n[1];

    if (hypot(deltax[0], deltax[1]) < epsilon3) {
      estimate[0] = x_n[0];
      estimate[1] = x_n[1];
      break;
    }

    x_n += deltax;
    std::valarray<real_t> tempfn(f_n);
    f_n = x_n;
    f(f_n);
    deltaf = f_n - tempfn;

    std::valarray<real_t> numerator(2);
    numerator[0] = -jacobian[0] * deltaf[0] - jacobian[1] * deltaf[1];
    numerator[1] = -jacobian[2] * deltaf[0] - jacobian[3] * deltaf[1];
    numerator += deltax;

    real_t denominator = deltax[0] * (jacobian[0] * deltaf[0] + jacobian[1] * deltaf[1]) +
                         deltax[1] * (jacobian[2] * deltaf[0] + jacobian[3] * deltaf[1]);

    numerator /= denominator;

    std::valarray<real_t> F = {0, 0, 0, 0};
    F[0] = numerator[0] * deltax[0];
    F[1] = numerator[0] * deltax[1];
    F[2] = numerator[1] * deltax[0];
    F[3] = numerator[1] * deltax[1];

    F[0] += 1;
    F[3] += 1;

    std::valarray<real_t> temp(jacobian);
    jacobian[0] = F[0] * temp[0] + F[1] * temp[2];
    jacobian[1] = F[0] * temp[1] + F[1] * temp[3];
    jacobian[2] = F[2] * temp[0] + F[3] * temp[2];
    jacobian[3] = F[2] * temp[1] + F[3] * temp[3];
  }

  if (n >= maxiter) {
    return false;
  }

  return true;
}


std::valarray<real_t> fixedpoint::initial_jacobian(const std::valarray<real_t> &est, const broyden_f &f, real_t eta) {
  std::valarray<real_t> jacobian = {0, 0, 0, 0};
  std::valarray<real_t> col1(est), col2(est), def(est);
  col1[0] += eta;
  col2[1] += eta;

  f(def);
  f(col1);
  col1 -= def;
  col1 /= eta;
  f(col2);
  col2 -= def;
  col2 /= eta;

  jacobian[0] = col1[0];
  jacobian[2] = col1[1];
  jacobian[1] = col2[0];
  jacobian[3] = col2[1];

  return jacobian;
}


int fixedpoint::closest_point_in_list(real_t y1, real_t y2, size_t prime_period_estimate, real_t epsilon) {
  int convergence_pos = 0;
  real_t dist;
  std::valarray<real_t> y_check(2);

  for (size_t i = 0; i < nfixed && convergence_pos == 0; i++) {
    if (period[i] != prime_period_estimate) continue;
    y_check[0] = y1;
    y_check[1] = y2;
    for (size_t j = 0; j <= period[i]; j++) {
      dist = hypot(fabs(pt[2*i] - y_check[0]), fabs(pt[2*i+1] - y_check[1]));
      if (dist < epsilon) {
        convergence_pos = i+1;
        break;
      }
      forced_damped::poincare(y_check, 1, 1);
    }
  }

  return convergence_pos;
}

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

int fixedpoint::check_point(const array_t &y, size_t pos) // assumed size 2
{
  using namespace std; // for simplicity
  size_t prime_period_estimate;
  int convergence_pos;

  // Step 2, estimate a low tolerance prime period
  prime_period_estimate = estimate_prime_period(y[pos], y[pos+1]);
  // Step 3, Maybe the low tolerance prime period helps us find a nearby point
  // on the list.
  convergence_pos = closest_point_in_list(y[pos], y[pos+1], prime_period_estimate, epsilon1);
  if (convergence_pos != 0) {
    exists++;
    return convergence_pos;
  }

  // Step 4, apply broydens method to find new period and better estimate of y1
  // and y2.
  std::valarray<real_t> estimate = {y[pos], y[pos+1]};
  bool converged = broyden(estimate, 100, [prime_period_estimate](std::valarray<real_t> &arr) {
        // broyden's f(arr[a, b] = poincare(arr, prime_period_estimate) - arr)
        std::valarray<real_t> temp(arr);
        forced_damped::poincare(arr, prime_period_estimate, 1);
        arr -= temp;
  });

  if (!converged) {
    not_conv++;
    return 0;
  }
  // Estimate period again?
  prime_period_estimate = estimate_prime_period(estimate[0], estimate[1]);

  // Calculate DP^(prime_period_estimate)
  std::valarray<real_t> col1(estimate), col2(estimate), def(estimate);
  col1[0] += TWOPI/100;
  col2[1] += TWOPI/100;
  forced_damped::poincare(col1, prime_period_estimate, 1);
  forced_damped::poincare(col2, prime_period_estimate, 1);
  forced_damped::poincare(def, prime_period_estimate, 1);
  col1 -= def;
  col2 -= def;
  col1 /= (TWOPI/100);
  col2 /= (TWOPI/100);
  std::valarray<real_t> jacobian = {col1[0], col2[0], col1[1], col2[1]};

  // Step 5, does a point within epsilon2 exist?
  convergence_pos = closest_point_in_list(estimate[0], estimate[1], prime_period_estimate, epsilon2);
  if (convergence_pos != 0) {
    exists++;
    return convergence_pos;
  }

  // Step 6, Estimate the eigenvalues to get the new points
  real_t a = 1.0;
  real_t b = jacobian[3] - jacobian[0];
  real_t c = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
  // If negative?
  real_t squared = b * b - 4 * a * c;

  real_t root1 = 10, root2 = 10;

  if (squared > 0) {
    // Catastrophic cancellation?
    if (b > 0) {
      // root1 = (-b + sqrt(b * b - 4 * a * c))/(2 * a);
      root1 = (2 * c)/(-b - sqrt(squared));
      root2 = (-b - sqrt(squared))/(2 * a);
    } else {
      root1 = (-b + sqrt(squared))/(2 * a);
      // root2 = (-b - sqrt(b * b - 4 * a * c))/(2 * a);
      root2 = (2 * c)/(-b + sqrt(squared));
    }
  }

  // store estimate and stability of the point
  pt.push_back(estimate[0]);
  pt.push_back(estimate[1]);
  std::cout << "Adding (" << estimate[0] << ", " << estimate[1] << "), Roots: [";
  if (squared < 0)
    std::cout << -b << " + " << abs(squared) << "i" << " , " << -b << " - " << abs(squared) << "i";
  else
    std::cout << root1 << ", " << root2;
  std::cout << "]" << std::endl;
  period.push_back(prime_period_estimate);
  stable.push_back((squared < 0 && fabs(hypot(-b, squared)) < 1.0) || (fabs(root1) < 1.0 && fabs(root2) < 1.0));
  nfixed++;
  return nfixed;
}

//----------------------------------------------------------------------------
// Compute the basin boundary of each fixed point.  For the (J,K)th grid point,
// compute up to MAXITER iterations of the Poincare map and see where the
// initial condition ends up.  If it lies within EPSILON of one of the fixed
// points (which CHECK_POINT will determine), then mark the basin accordingly.

std::valarray<int> fixedpoint::compute_basin(const grid &box) {
  size_t n = box.resolution;
  std::valarray<int> basin(n * n);

  real_t x_step = (box.xmax - box.xmin) / n;
  real_t y_step = (box.ymax - box.ymin) / n;

  for (size_t i = 0; i < n; i++) {
    array_t y(2 * n);
    real_t x_i = box.xmin + x_step * i;

    for (size_t j = 0; j < n; j++) {
      real_t y_i = box.ymin + y_step * j;
      basin[i * n + j] = 0;
      y[2 * j] = x_i;
      y[2 * j + 1] = y_i;
    }

    forced_damped::poincare(y, maxiter, n);
    for (size_t j = 0; j < n; j++)
      basin[i*n+j] = check_point(y, 2*j);
  } // End basin calculation

  std::cout << "There are " << not_conv << " points not converged in broyden." << std::endl;
  std::cout << "There are " << exists << " points that matched in list." << std::endl;
  std::cout << "There are " << nfixed << " points discovered." << std::endl;
  for (size_t i = 0; i < nfixed; i++) {
      std::cout << "[" << (stable[i] ? "STABLE" : "UNSTABLE") << "] Period " << period[i] << ": (" << pt[2*i] << ", " << pt[2*i + 1] << ")";
      if (period[i] > 1) {
        std::cout << " -> [";
        std::valarray<real_t> t = {pt[2*i], pt[2*i+1]};
        for (size_t j = 1; j <= period[i]; j++) {
          forced_damped::poincare(t, 1, 1);
          std::cout << "(" << t[0] << ", " << t[1] << ")";
          if (j < period[i]) std::cout << ", ";
        }
        std::cout << "]";
      }
      std::cout << std::endl;
  }

  return (basin);
}

//----------------------------------------------------------------------------
// Main program.

int main(int argc, char **argv) {
  using namespace std;

  fixedpoint fp; //  Our only instance
  grid box;  // will have the default limits and resolution

  // Replace the contents of any previous file and output as binary
  ofstream outfile("basin.dat", ios::trunc | ios::binary);

  auto start = std::chrono::high_resolution_clock::now();
  valarray<int> basin(fp.compute_basin(box));
  auto stop = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  cout << duration.count() << endl;

  // Write out the results in a binary format compatible with the MATLAB
  // script wada.m for visualization.

  outfile.write(reinterpret_cast<const char *>(&box.xmin), sizeof(box.xmin));
  outfile.write(reinterpret_cast<const char *>(&box.xmax), sizeof(box.xmax));
  outfile.write(reinterpret_cast<const char *>(&box.ymin), sizeof(box.ymin));
  outfile.write(reinterpret_cast<const char *>(&box.ymax), sizeof(box.ymax));
  outfile.write(reinterpret_cast<const char *>(&box.resolution),
                sizeof(box.resolution));
  outfile.write(reinterpret_cast<const char *>(&basin[0]),
                sizeof(int) * basin.size()); // N x N basin boundary
  outfile.close();
  return (EXIT_SUCCESS);
}
