!  PENDULUM - the forced damped nonlinear pendulum.
!
      module pendulum
      use rkmod
      use,intrinsic:: ieee_arithmetic,only: ieee_rem  ! optional
      implicit none
      real(WP), parameter:: PI=3.1415926535897932_WP
      real(WP), parameter:: TWOPI=2*PI
!
!  Default map parameters
!
      real(WP), save, private:: damp=0.2_WP
      real(WP), save, private:: force=1.66_WP
      contains
!-----------------------------------------------------------------------
      pure subroutine poincare(y, niter)
!  POINCARE - 2-pi Poincare map for the forced damped pendulum.
!
      integer, intent(in):: niter
      real(WP), intent(inout):: y(2)
!
      integer:: k
      integer, parameter:: NSTEPS=256
      real(WP), parameter:: ZERO=0.0

      do k=1,nsteps
         call rk4(forced_damped, 2, y, ZERO, TWOPI, NSTEPS)
         !  Your code here to implement the appropriate modulus
         !  operation(s) on the position variable, y(1)
      enddo
      return
      end subroutine poincare
!-----------------------------------------------------------------------
      pure subroutine forced_damped(y, n, t, dy)
!  Y : 2-vector giving (position, velocity).
!  DY := corresponding time derivatives.
!
      integer, intent(in):: n
      real(WP), intent(in):: y(n), t
      real(WP), intent(out):: dy(n)
!
      if(n.ne.2) then  ! can also be written (n /= 2)
         dy = 0.0
      else
        dy(1) = y(2)
        dy(2) = force*cos(t) - sin(y(1)) - damp*y(2)
      endif
      return
      end subroutine forced_damped
!-----------------------------------------------------------------------
      end module pendulum
