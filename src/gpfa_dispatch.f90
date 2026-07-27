module gpfa_dispatch
   use, intrinsic :: iso_fortran_env, only: real64
   use gpfa_radix2, only: gpfa2f
   use gpfa_radix3, only: gpfa3f
   use gpfa_radix5, only: gpfa5f
   implicit none
contains

   subroutine gpfa(a, b, trigs, inc, jump, n, lot, isign)
      implicit none
      integer :: inc, jump, n, lot, isign, nn, ifac, ll, kk, nj(3), ip, iq, ir, i
      real(real64) :: a(*), b(*), trigs(*)
!
!     decompose n into factors 2,3,5
!     ------------------------------
      nn = n
      ifac = 2
!
      do ll = 1, 3
         kk = 0
         do while (mod(nn, ifac) .eq. 0)
            kk = kk + 1
            nn = nn / ifac
         end do
         nj(ll) = kk
         ifac = ifac + ll
      end do
!
      if (nn .ne. 1) then
         write (6, 40) n
40       format(' *** warning!!!', i10, ' is not a legal value of n ***')
         return
      end if
!
      ip = nj(1)
      iq = nj(2)
      ir = nj(3)
!
!     compute the transform
!     ---------------------
      i = 1
      if (ip .gt. 0) then
         call gpfa2f(a, b, trigs, inc, jump, n, ip, lot, isign)
         i = i + 2 * (2**ip)
      end if
      if (iq .gt. 0) then
         call gpfa3f(a, b, trigs(i), inc, jump, n, iq, lot, isign)
         i = i + 2 * (3**iq)
      end if
      if (ir .gt. 0) then
         call gpfa5f(a, b, trigs(i), inc, jump, n, ir, lot, isign)
      end if
!
   end subroutine gpfa

end module gpfa_dispatch
