module gpfa_controller
   use, intrinsic :: iso_fortran_env, only: real64
   use gpfa_dispatch, only: gpfa
   implicit none
contains

   subroutine cgpfa(cr, ci, trig, nblk, m, isign)
!      use iso_c_binding
      implicit none
      integer :: m, isign, nblk, inc
      real(real64) :: trig(*), cr(nblk * m), ci(nblk * m)
!      real(real64), pointer :: cr(:)
!      real(real64) :: cr(2*nblk*m)
!      complex(real64), target :: c(nblk*m)
!      call C_F_POINTER(C_LOC(c), cr, [2*nblk*m])
!      cr(1:2*nblk*m-1:2)=dble(c(1:nblk*m))
!      cr(2:2*nblk*m:2)=aimag(c(1:nblk*m))
!      inc=2*nblk
      inc = nblk
      call gpfa(cr, ci, trig, inc, 1, m, nblk, isign)
!      c(1:nblk*m)=cmplx(cr(1:2*nblk*m-1:2),cr(2:2*nblk*m:2),kind=real64)
   end subroutine cgpfa

end module gpfa_controller
