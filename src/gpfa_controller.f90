module gpfa_controller
   use gpfa_dispatch, only: gpfa
   implicit none
contains

   subroutine cgpfa(cr, ci, trig, nblk, m, isign)
!      use iso_c_binding
      implicit none
      integer :: m, isign, nblk, n, i, inc
      real(8) :: trig(*), cr(nblk * m), ci(nblk * m)
!      real(8), pointer :: cr(:)
!      real(8) :: cr(2*nblk*m)
!      complex(8), target :: c(nblk*m)
!      call C_F_POINTER(C_LOC(c), cr, [2*nblk*m])
!      cr(1:2*nblk*m-1:2)=dble(c(1:nblk*m))
!      cr(2:2*nblk*m:2)=aimag(c(1:nblk*m))
!      inc=2*nblk
      inc = nblk
      do n = 1, nblk
!         i=2*n-1
         i = n
         call gpfa(cr(i:), ci(i:), trig, inc, 1, m, 1, isign)
      end do
!      c(1:nblk*m)=cmplx(cr(1:2*nblk*m-1:2),cr(2:2*nblk*m:2),kind=kind(0.0d0))
   end subroutine cgpfa

end module gpfa_controller
