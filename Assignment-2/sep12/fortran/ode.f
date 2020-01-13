!  MAIN - test program for the RK method.
!
      module ode
      use rkmod
      implicit none
      contains
!------------------------------------------------------------------------------
      subroutine testproblem(t, y, n, dy)
!  TESTPROBLEM - this routine's calling sequency must match that in the
!  abstract interface declared in rkmod.f.  The names of the dummy variables
!  need not be the same, but their type, kind, and rank must match.
!  Nere H is a
!
      integer, intent(in):: n  ! assumed 1 for testing
      real(WP), intent(in):: t, y(n)
      real(WP), intent(out):: dy(n)
!
!  your code here
!
      return
      end subroutine testproblem
!------------------------------------------------------------------------------
      end module ode
!------------------------------------------------------------------------------
      program test
!  TEST - test the RK solver
      use ode
      use, intrinsic:: iso_fortran_env
      implicit none
!
      real(WP):: tstart=0, tend=2, one=1, y(1)
      integer:: nsteps=16

      y(1) = 1.0
      call rk4(testproblem, y, 1, tstart, tend, nsteps)
      write(OUTPUT_UNIT,*) 'rk4:', y
      write(OUTPUT_UNIT,*) 'truth:', exp(one)
      stop
      end program test

