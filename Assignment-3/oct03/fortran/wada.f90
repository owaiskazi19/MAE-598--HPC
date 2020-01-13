!  WADA - find the periodic points in a given domain for the
!  forced damped nonlinear pendulum with a fixed set of parameters.
!
      module wada
      use pendulum
      implicit none
!
!  Grid that defines the region over which we compute the basin.
!  The initializations below are the default values assigned to each
!  component when a new instance is created.
!
      type::grid
         real(WP):: xmin=-PI  ! grid limits in X and Y
         real(WP):: xmax=PI
         real(WP):: ymin=-PI
         real(WP):: ymax=PI
         integer:: resolution=400  ! grid resolution in each direction
      end type grid
!
!  Fixed points obtained with the default map parameters.
!  These values are all precomputed by another program.
!
      integer,parameter:: MAXPERIOD=2
      integer,parameter:: MAXPTS=3
      type:: fixedpoints
         real(WP)::pt(2,MAXPTS)=reshape([                    &
           -2.407661292643_WP, 0.6773267621896_WP,           &
           -0.6099150484926_WP,  1.677463819037_WP,          &
           -0.7218820695570_WP, -0.4620944521247_WP], [2, MAXPTS])
         real(WP):: epsilon=1.0e-09_WP  ! distance threshold
         integer:: period(MAXPTS)=[1,2,2]
         integer:: maxiter=400  ! max number of map iterations
         integer:: nfound=MAXPTS
      end type fixedpoints
!
!  The only instance of the FIXEDPOINTS type in this program.
!
      type(fixedpoints), save, private:: fp
      contains
!----------------------------------------------------------------------
      subroutine check_point(y,eps,idx)
!  CHECK_POINT - determine which of the periodic points, if any,
!  Y and its orbit approach.
!  Y : the 2-vector whose orbit is examined.
!  EPS : distance threshold; if the orbit approaches within EPS, then
!   we conclude that the limiting point is the IDXth one.
!  IDX := set to M if the orbit of Y comes within EPS of the Mth fixed point,
!   M=1,2,3.  Otherwise, IDX is set to 0.
!
      real(WP), intent(in):: y(2), eps
      integer, intent(out):: idx

      !  Your code here

      idx = 0  ! set to the appropriate value
      return
      end subroutine check_point
!----------------------------------------------------------------------
      subroutine compute_basin(box,basin)
!  COMPUTE_BASIN - compute the Wada basin boundary.
!
      type(grid), intent(in):: box
      integer, allocatable, intent(out):: basin(:,:)
!
!  Local variables
!
      integer:: n  ! among others that you will need
!
!  The array BASIN holds the fate of each point in the grid.
!  BASIN(J,K) is one of 1, 2, or 3 according as the initial condition
!  (x_j, y_k) tends to fixed point A, B or C (defined as fp%pt(:,1),
!  fp%pt(:,2), and fp%pt(:,3), respectively), or as 0 if the orbit tends
!  to none of them.  The initialization of the array defines each point
!  as none of them, but your code here will made updated determinations.
!
      n = box%resolution
      allocate(basin(n,n))  ! N x N array that defines the fate of each point
      basin = 0
!
!  Compute the basin boundary
!
      return
      end subroutine compute_basin
!----------------------------------------------------------------------
      end module wada
!----------------------------------------------------------------------
! MAIN.
      use wada
      implicit none
      integer, parameter:: IOUT=3  ! output unit number
      integer, allocatable:: basin(:,:)
      type(grid):: box
!
      open(IOUT, file='basin.dat', form='unformatted', access='stream',  &
        status='replace', position='rewind')
      call compute_basin(box,basin)
      write(IOUT) box, basin
      close(IOUT)
      stop
      end

