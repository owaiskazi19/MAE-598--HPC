#include "rk4.h"
#include "pendulum.h"
#include <stdio.h>

//  The maximum period of any orbit under consideration by this program.
static size_t maxperiod = 2;

//  Here is one possible way to collect all the information relevant to
//  tracking the periodic orbits.  In this format, only one representative
//  of each periodic orbit is stored in PT; you have to iterate the map
//  to check whether a candidate point comes within EPSILON of the stored orbit.

struct fixedpoint {
   real_t pt[3][2];  // the list of periodic points
   real_t epsilon;   // convergence criterion
   size_t period[3];  // their respective periods
   size_t maxiter;  // maximum number of map iterations
   size_t nfixed;  // the number of fixed points in the list (i.e., 3)
};

// The fixed points whose basins we will calculate.  One representative
// of each point is stored.  This is a static initialization; the assignments
// are done at compile and link time.

static struct fixedpoint fp = {
   {{-2.407661292643E+00,  6.773267621896E-01},  // period 1
    {-6.099150484926E-01,  1.677463819037E+00},  // period 2
    {-7.218820695570E-01, -4.620944521247E-01}  // period 2
   },  // pt
   1.0E-09,  //  required approach distance to confirm convergence
   {1, 2, 2}, // the respecitve periods
   400, // maxiter
   3 // nfixed
};

//  This is one possible way to represent the grid to be used for computing
//  the basin boundary.  We assume a square grid here but it's straightforward
//  to extend the structure to accommodate a rectangular one.

struct grid {
   real_t xmin;  // angular position limits
   real_t xmax;
   real_t ymin;  // angular velocity limits
   real_t ymax;
   size_t resolution;  // number of grid points in each direction
};

// The grid used to generate the basin.  Unfortunately, in Standard C,
// you cannot declare PI and MINUSPI as CONST quantities, which is less
// error-prone than simple macros.  Expressions like -PI (where PI is a
// CONST variable) are permitted in GNU C as an extension, but not in Intel C.
// For portability reasons, it is best to avoid language extensions.

static struct grid box = {
  MINUSPI, PI,  // xmin, xmax
  MINUSPI, PI,  // ymin, ymax
  400  // resolution
};
//----------------------------------------------------------------------------
// Check whether Y or its orbit contains a point that is already in the
// list of fixed points.  If so, then return an integer value (1, 2, or 3)
// indicating whether it is the first, second, or third point in FIXEDPOINT.
// In other words, if Y lies within fp.epsilon of pt[i], then return i+1,
// i = 0, 1, 2.  Otherwise, return 0.
// You can use any metric here.  The 2-norm in C can be computed with the
// standard library function HYPOT(), but the 0-norm also works as the
// maximum of the absolute values of the differences between any component.
// C99 includes the functions FMAX() and FABS(), which are convenient for this
// purpose.

int
check_point(real_t y[2])
{
   // your code here
   return(0);  // substitute the appropriate value
}
//----------------------------------------------------------------------------
// Compute the basin boundary of each fixed point.  For the (J,K)th grid point,
// compute up to MAXITER iterations of the Poincare map and see where the
// initial condition ends up.  If it lies within EPSILON of one of the fixed
// points (which CHECK_POINT will determine), then mark the basin accordingly.

int*
compute_basin(const struct grid *b)
{
    size_t n = b->resolution;
    int *basin = (int *) calloc(n*n, sizeof(*basin));  // initialized to 0

    if(basin == NULL) return(NULL);  // allocation failure

    // your code here
    return(basin);
}
//----------------------------------------------------------------------------
// Main program.  
// The version below simply allocates an NxN array of integers for the
// basin boundary array, and you will need to do all the subscript arithmetic
// explicitly to compute the offset needed to access the (j,k)th element.
// Alternatively, you can to convert BASIN to an int**,
// allocate N int*'s, then allocate the NxN block and precompute the
// address of each row (see Sept. 26 lecture notes for details).
// In this case, you'll need to modify the FWRITE call that writes out
// the BASIN array.
int
main(int argc, char **argv)
{
    size_t n = box.resolution;
    const char *filename = "basin.dat";
    FILE *outfp = fopen(filename, "wb");   // 'b' for binary
    int *basin;

    if(outfp == NULL) {  // cannot open output file, so halt
       perror(filename);
       return(EXIT_FAILURE);
    }

    basin = compute_basin(&box);
    if(basin == NULL) {  // heap disaster; printing a message is
       return(EXIT_FAILURE);   // unlikely to succeed
    }

    // Write out the results in a binary format compatible with the MATLAB
    // script wada.m for visualization.

    fwrite(&box, sizeof(box.xmin), 4, outfp);  // xmin through ymax
    fwrite(&fp.epsilon, sizeof(fp.epsilon), 1, outfp);
    fwrite(&box.resolution, sizeof(box.resolution), 1, outfp);
    fwrite(basin, sizeof(*basin), n*n, outfp);  // N x N basin boundary
    fclose(outfp);
    free(basin);
    return(EXIT_SUCCESS);
}
