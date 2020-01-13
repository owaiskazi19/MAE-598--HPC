// Implemented in collaboration with Amitabh Das
#include "pendulum.h"
#include "rk4.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <mpi.h>

#include <cstdint>
#include <climits>
#include <csignal>

#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "what is happening here?"
#endif

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
  real_t epsilon;   // convergence criterion
  std::vector<size_t> period; // their respective periods
  size_t maxiter;   // maximum number of map iterations
  size_t nfixed;    // the number of fixed points in the list (i.e., 3)
  size_t maxperiod; // of any orbit under consideration by this program
  fixedpoint(const fixedpoint &);             // no public copy constructor
  fixedpoint(const fixedpoint &&);            // no public move constructor
  fixedpoint &operator=(const fixedpoint &);  // no public copy assignment
  fixedpoint &operator=(const fixedpoint &&); // no public move assignment

  int check_point(const array_t &, size_t pos);
  int add_to_list(const real_t, const real_t, size_t);
  void gather_results(int, int, std::valarray<int> &);
  void send_results(int, int, int, int, const std::valarray<int> &);

  public:
  // default constructor
  fixedpoint() : pt(), epsilon{1.0e-09}, period(), maxiter{400}, nfixed{0}, maxperiod{4} {};
  std::valarray<int> compute_basin(const grid &, int, int);
};

// Include the following statement if you want to use the static class defined
// in pendulum.cc.  You can call the Poincare map defined there as
//    forced_damped::poincare(y, niter);
// for appropriately defined variables Y and NITER.

class forced_damped;

int fixedpoint::add_to_list(const real_t y1, const real_t y2, size_t point_period) {
  int convergence_pos = 0;
  real_t dist;
  // Step 4
  // If there are points in L, find the closest one
  std::valarray<real_t> y_check(2);
  for (size_t i = 0; i < nfixed && convergence_pos == 0; i++) {
    if (period[i] != point_period) continue;
    y_check[0] = pt[2*i];
    y_check[1] = pt[2*i+1];
    for (size_t j = 0; j < period[i]; j++) {
      dist = hypot(fabs(y1 - y_check[0]), fabs(y2 - y_check[1]));
      if (dist < epsilon) {
        convergence_pos = i+1;
        break;
      }
      forced_damped::poincare(y_check, 1, 1);
    }
  }

  // Step 5
  // If given point not close to any point in L, add given point and its period
  // to L
  if (convergence_pos == 0) {
    pt.push_back(y1);
    pt.push_back(y2);
    period.push_back(point_period);
    nfixed++;
    convergence_pos = nfixed;
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
  valarray<real_t> y_period(2);
  real_t dist = epsilon;

  // Step 3
  // Find prime period
  size_t point_period = 0;
  y_period[0] = y[pos];
  y_period[1] = y[pos + 1];

  while (point_period < maxperiod) {
    forced_damped::poincare(y_period, 1, 1);
    dist = hypot(fabs(y[pos] - y_period[0]), fabs(y[pos+1] - y_period[1]));
    if (dist <= epsilon/2) break;
    point_period++;
  }

  // If aperiodic within maxperiod, return 0
  if (point_period == maxperiod)
    return 0;

  return add_to_list(y_period[0], y_period[1], point_period + 1);
}


void fixedpoint::gather_results(int world_size, int resolution, std::valarray<int> &basin) {
    MPI_Status status;
    int sender;
    size_t size, p_start, p_end, split_size = resolution / world_size;

    MPI_Recv(&size, 1, my_MPI_SIZE_T, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, &status);
    sender = status.MPI_SOURCE;

    p_start = sender * split_size;
    p_end = sender == world_size - 1 ? resolution : p_start + split_size;

    std::vector<size_t> period(size);
    std::vector<real_t> points(2 * size);
    std::valarray<int> part_basin((p_end - p_start) * resolution);
    std::cout << sender << " calculated for rows " << p_start << " to " << p_end << std::endl;

    MPI_Recv(&period[0], size, my_MPI_SIZE_T, sender, 14, MPI_COMM_WORLD, &status);
    MPI_Recv(&points[0], size * 2, MPI_DOUBLE, sender, 15, MPI_COMM_WORLD, &status);
    MPI_Recv(&part_basin[0], ((p_end - p_start) * resolution), MPI_DOUBLE, sender, 16, MPI_COMM_WORLD, &status);

    // From points index to fp.pt index
    std::valarray<size_t> map(size + 1);

    bool redo = false;
    for (size_t k = 0; k < size; k++) {
      size_t pt_period = period[k];
      real_t y1 = points[2 * k];
      real_t y2 = points[2 * k + 1];
      std::cout << "sub-process found point: (" << y1 << ", " << y2 << "), Period: " << pt_period;

      map[k + 1] = add_to_list(y1, y2, pt_period);
      std::cout << " <> " << k + 1 << " matched with point at index " << map[k + 1] << std::endl;
      redo |= k + 1 != map[k + 1];
    }

    if (redo) {
      std::cout << "Mismatch in point index detected. Re-doing part-basin." << std::endl;
      // Translate the returned values
      for (size_t i = 0; i < part_basin.size(); i++) {
        basin[p_start * resolution + i] = map[part_basin[i]];
      }
    } else {
      for (size_t i = 0; i < part_basin.size(); i++) {
        basin[p_start * resolution + i] = part_basin[i];
      }
    }

    return;
}


void fixedpoint::send_results(int target, int world_rank, int resolution, int world_size, const std::valarray<int> &basin) {
    size_t p_start, p_end, split_size = resolution / world_size;

    p_start = world_rank * split_size;
    p_end = world_rank == world_size - 1 ? resolution : p_start + split_size;

    MPI_Send(&nfixed, 1, my_MPI_SIZE_T, target, 11, MPI_COMM_WORLD);
    MPI_Send(&period[0], nfixed, my_MPI_SIZE_T, target, 14, MPI_COMM_WORLD);
    MPI_Send(&pt[0], nfixed * 2, MPI_DOUBLE, target, 15, MPI_COMM_WORLD);
    MPI_Send(&basin[0], ((p_end - p_start) * resolution), MPI_INT, target, 16, MPI_COMM_WORLD);
}


//----------------------------------------------------------------------------
// Compute the basin boundary of each fixed point.  For the (J,K)th grid point,
// compute up to MAXITER iterations of the Poincare map and see where the
// initial condition ends up.  If it lies within EPSILON of one of the fixed
// points (which CHECK_POINT will determine), then mark the basin accordingly.
std::valarray<int> fixedpoint::compute_basin(const grid &box, int world_rank, int world_size) {
  size_t n = box.resolution, split_size = n / world_size;

  // Start row to end row managed by a process
  size_t start = world_rank * split_size;
  size_t end = world_rank == world_size - 1 ? box.resolution : start + split_size;

  std::valarray<int> basin(n * (world_rank == 0 ? n : (end - start)));

  std::cout << world_rank << ": Working from " << start << " to " << end - 1 << std::endl;

  real_t x_step = (box.xmax - box.xmin) / n;
  real_t y_step = (box.ymax - box.ymin) / n;

  for (size_t i = 0; i < (end - start); i++) {
    array_t y(2 * n);
    real_t x_i = box.xmin + x_step * (i + start);

    for (size_t j = 0; j < n; j++) {
      real_t y_i = box.ymin + y_step * j;
      basin[i * n + j] = 0;
      y[2 * j] = x_i;
      y[2 * j + 1] = y_i;
    }

    forced_damped::poincare(y, maxiter, n);
    for (size_t j = 0; j < n; j++)
      basin[i * n + j] = check_point(y, 2 * j);
  } // End basin calculation

  if (world_rank != 0) {
    std::cout << world_rank << " sending results..." << std::endl;
    send_results(0, world_rank, n, world_size, basin);
    std::cout << world_rank << " finished sending results." << std::endl;
  } else {
    for (size_t k = 0; k < nfixed; k++) {
      std::cout << "Point: (" << pt[2*k] << ", " << pt[2*k+1] << "), Period " << period[k] << std::endl;
    }

    std::cout << world_rank << " gathering results..." << std::endl;
    for (int i = 1; i < world_size; i++)
      gather_results(world_size, n, basin);
    std::cout << world_rank << " finished gathering results." << std::endl;
  }

  return (basin);
}

void signalHandler( int signum ) {
  // cleanup and close up stuff here
  MPI_Abort(MPI_COMM_WORLD, signum);
  // terminate program
  exit(signum);
}

//----------------------------------------------------------------------------
// Main program.

int main(int argc, char **argv) {
  using namespace std;
  int world_size, world_rank;

  MPI_Init(nullptr, nullptr);

  try {
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    struct sigaction sigHandler;

    sigHandler.sa_handler = signalHandler;
    sigemptyset(&sigHandler.sa_mask);
    sigHandler.sa_flags = 0;

    sigaction(SIGABRT, &sigHandler, nullptr);
    sigaction(SIGILL, &sigHandler, nullptr);
    sigaction(SIGSEGV, &sigHandler, nullptr);
    sigaction(SIGTERM, &sigHandler, nullptr);
    sigaction(SIGINT, &sigHandler, nullptr);

    std::valarray<real_t> params(3);
    fixedpoint fp; //  Our only instance

    if (world_rank == 0) {
      std::cin >> params[0] >> params[1] >> params[2];
    }

    MPI_Bcast(&params[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    forced_damped::damping = params[0];
    forced_damped::force = params[1];
    grid box(params[2]);  // will have the default limits and resolution

    if (static_cast<size_t>(world_size) > box.resolution) {
      world_size = box.resolution;
    }

    if (static_cast<size_t>(world_rank) >= box.resolution) {
      std::cout << "Exiting process " << world_rank << std::endl;
      MPI_Finalize();
      return (EXIT_SUCCESS);
    }

    valarray<int> basin(fp.compute_basin(box, world_rank, world_size));

    if (world_rank == 0) {
      ofstream outfile("basin.dat", ios::trunc | ios::binary);
      outfile.write(reinterpret_cast<const char *>(&box.xmin), sizeof(box.xmin));
      outfile.write(reinterpret_cast<const char *>(&box.xmax), sizeof(box.xmax));
      outfile.write(reinterpret_cast<const char *>(&box.ymin), sizeof(box.ymin));
      outfile.write(reinterpret_cast<const char *>(&box.ymax), sizeof(box.ymax));
      outfile.write(reinterpret_cast<const char *>(&box.resolution),
          sizeof(box.resolution));
      outfile.write(reinterpret_cast<const char *>(&basin[0]),
          sizeof(int) * basin.size()); // N x N basin boundary
      outfile.close();
    }

  } catch (...) {
    std::cout << "Exception occured!" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return (EXIT_FAILURE);
  }

  // Finalize the MPI environment.
  MPI_Finalize();
  return (EXIT_SUCCESS);
}
