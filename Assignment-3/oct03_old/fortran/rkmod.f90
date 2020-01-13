!  RKMOD - interface and code for a fixed-step RK4 method.
!
      module rkmod
      use precision
      implicit none
!
!  Interface for vector fields that this code can integrate.
!  N is the number of equations.
!
      abstract interface
         pure subroutine vectorfield(y,n,t,dy)
         import
         integer,intent(in)::n
         real(WP),intent(in)::t,y(n)
         real(WP),intent(out)::dy(n)
         end subroutine vectorfield
      end interface
!
!  This declaration is optional but make RK4 the only exported procedure name.
!
      private::rkstep
      contains
!-----------------------------------------------------------------------
      pure subroutine rk4(f,n,y,tstart,tend,nsteps)
!  RK4 - integrate the vector field F from TSTART to TEND in NSTEPS equal steps.
!
      procedure(vectorfield)::f
      integer,intent(in)::n
      real(WP),intent(in)::tstart,tend
      real(WP),intent(inout)::y(n)
      integer,intent(in)::nsteps
!
!  Local variables
!
      integer::j
      real(WP)::h,t
!
      h=(tend-tstart)/nsteps
      do j=1,nsteps
         t=tstart+(j-1)*h
         call rkstep(f,y,n,t,h)
      enddo
      return
      end subroutine rk4
!-----------------------------------------------------------------------
       pure subroutine rkstep(f,y,n,t,h)
!  RKSTEP - take one Runge-Kutta step of size H.
!
      procedure(vectorfield)::f
      integer,intent(in)::n
      real(WP),intent(in)::t,h
      real(WP),intent(inout)::y(n)
!
!  Local variables
!
      real(WP)::dy1(n), dy2(n), dy3(n), dy4(n), ytemp(n)
!
!  Your previous code here
!
      return
      end subroutine rkstep
!----------------------------------------------------------------------
      end module rkmod
