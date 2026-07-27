module gpfa_setup
   use constants, only: two_pi
   implicit none
contains

   subroutine setgpfa(trigs, n)
      implicit none
      integer :: n, nn, ifac, ll, kk, nj(3), ip, iq, ir, ni, irot, kink, k, i
      real(8) :: trigs(*), del, angle
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
!     compute list of rotated twiddle factors
!     ---------------------------------------
      nj(1) = 2**ip
      nj(2) = 3**iq
      nj(3) = 5**ir
!
      i = 1
!
      do ll = 1, 3
         ni = nj(ll)
         if (ni .eq. 1) cycle
!
         del = two_pi / dble(ni)
         irot = n / ni
         kink = mod(irot, ni)
         kk = 0
!
         do k = 1, ni
            angle = dble(kk) * del
            trigs(i) = cos(angle)
            trigs(i + 1) = sin(angle)
            i = i + 2
            kk = kk + kink
            if (kk .gt. ni) kk = kk - ni
         end do
      end do
   end subroutine setgpfa

end module gpfa_setup
