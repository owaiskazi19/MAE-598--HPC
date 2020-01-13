!  RK4 - Module containing the fourth-order RK method and calling interfaces.
!  This code is written so that it can be compiled either as free- or fixed-form
!  source.
!  Fixed-form source is "old style": statements start in column 7 (or later)
!  and end by column 72 (or sooner).  If you need to continue a statement onto
!  the next line, place an ampersand (&) in column 6.
!  Free-form source is "new style": you can start a statment anywhere, and lines
!  can be up to 132 characters lone.  To continue a long statement onto the
!  next line, end it with an &.  Optionally, you can start the continuation
!  line with an & as well.
!  It is possible to combine the two: start in column 7 and place ampersands
!  in both column 73 on the current line and column 6 in the next line for
!  continuation.
!  AVOID TABS in Fortran source.
!  Include your name(s) in this source file when you turn it in.
!  You may modify the comments any way you wish, as long as you use a
!  consistent style.
!
      module rkmod
      implicit none
      integer, parameter:: WP=kind(1.0d0)  ! double precision, here
!
!  Define the calling interface for the vector field.
!  This code is Fortran 2003 (or later).
!
      abstract interface
         subroutine vectorfield(t, y, n, dy)
         import  ! WP and any others
         integer, intent(in):: n
         real(WP), intent(in):: t,y(n)
         real(WP), intent(out):: dy(n)
         end subroutine vectorfield
      end interface
      contains
!-----------------------------------------------------------------------
      subroutine rk4(f, y, n, tstart, tend, nsteps)
!  RK4 - integrate F from TSTART to TEND in NSTEPS steps.
!
      procedure(vectorfield):: f
      integer, intent(in):: n, nsteps
      real(WP), intent(in):: tstart, tend
      real(WP), intent(inout):: y(n)
!
!  Local variables  (H and loop counters as needed)
!
      ! your code here
      return
      end subroutine rk4
!-----------------------------------------------------------------------
       subroutine rkstep(f,y,n,t,h)
!  RKSTEP - take one Runge-Kutta step of size H.
!
      procedure(vectorfield):: f
      integer, intent(in):: n
      real(WP), intent(in):: t, h
      real(WP), intent(inout):: y(n)
!
!  Local variables
!
      real(WP):: dy1(n), ytemp(n)  ! and others as needed
!
      call f(t, y, n, dy1)
      ytemp = y + (0.5*h)*dy1  ! Euler step

      ! your code here
      return
      end subroutine rkstep
!----------------------------------------------------------------------
      end module rkmod
