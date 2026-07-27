!
!  numerical constants
!
!
!  last revised: 15 January 2011
!
module numconstants
   use mpidefs
   implicit none
   logical, target :: light_up
   integer :: print_intermediate_results, global_rank
   integer, allocatable :: monen(:)
   integer, private :: nmax = 0
   real(8) :: pi
   real(8), allocatable :: bcof(:, :), fnr(:), vwh_coef(:, :, :, :)
   real(8), allocatable :: vcc_const(:, :, :), fnm1_const(:, :), fn_const(:, :), fnp1_const(:, :)
   real(8), allocatable :: tran_coef(:, :, :)
   data pi/3.1415926535897932385/
   data light_up/.false./

contains

   subroutine init(notd)
      implicit none
      integer :: notd, l, n, ierr, nbc, m, mm1, mp1, np1, nm1, nn1
      real(8) :: fnorm1, fnorm2
!
!  bcof(n,l)=((n+l)!/(n!l!))^(1/2)
!
      if (notd .le. nmax) return
      nmax = max(nmax, notd)
      nbc = 6 * notd + 6
      if (allocated(fnr)) deallocate (monen, fnr, bcof)
      allocate (monen(0:2 * notd), bcof(0:nbc, 0:nbc), fnr(0:2 * nbc), stat=ierr)
!         write(*,'('' nmax, bcof status:'',2i5)') nmax,ierr
      do n = 0, 2 * notd
         monen(n) = (-1)**n
      end do
      fnr(0) = 0.d0
      do n = 1, 2 * nbc
         fnr(n) = dsqrt(dble(n))
      end do
      bcof(0, 0) = 1.d0
      do n = 0, nbc - 1
         do l = n + 1, nbc
            bcof(n, l) = fnr(n + l) * bcof(n, l - 1) / fnr(l)
            bcof(l, n) = bcof(n, l)
         end do
         bcof(n + 1, n + 1) = fnr(n + n + 2) * fnr(n + n + 1) * bcof(n, n) / fnr(n + 1) / fnr(n + 1)
      end do
      if (allocated(vwh_coef)) deallocate (vwh_coef)
      allocate (vwh_coef(-notd:notd, 1:notd, -1:1, -1:1))
!
!  constants used for calculation of svwf functions.
!
      do n = 1, notd
         nn1 = n * (n + 1)
         np1 = n + 1
         nm1 = n - 1
         fnorm1 = -.5d0 / fnr(n + n + 1) / fnr(n) / fnr(n + 1)
         fnorm2 = -.5d0 * fnr(n + n + 1) / fnr(n) / fnr(n + 1)
         m = -n
         mp1 = m + 1
         mm1 = m - 1
         vwh_coef(m, n, 1, 1) = -fnorm1 * n * fnr(np1 + m) * fnr(np1 + mp1)
         vwh_coef(m, n, 1, -1) = fnorm1 * np1 * fnr(n - m) * fnr(nm1 - m)
         vwh_coef(m, n, -1, 1) = fnorm1 * n * fnr(np1 - m) * fnr(np1 - mm1)
         vwh_coef(m, n, -1, -1) = 0.d0
         vwh_coef(m, n, 0, 1) = fnorm1 * n * fnr(np1 + m) * fnr(np1 - m)
         vwh_coef(m, n, 0, -1) = 0.d0
         vwh_coef(m, n, 1, 0) = -fnorm2 * fnr(n - m) * fnr(np1 + m)
         vwh_coef(m, n, -1, 0) = -0.d0
         vwh_coef(m, n, 0, 0) = -fnorm2 * m
         do m = -n + 1, -1
            mp1 = m + 1
            mm1 = m - 1
            vwh_coef(m, n, 1, 1) = -fnorm1 * n * fnr(np1 + m) * fnr(np1 + mp1)
            vwh_coef(m, n, 1, -1) = fnorm1 * np1 * fnr(n - m) * fnr(nm1 - m)
            vwh_coef(m, n, -1, 1) = fnorm1 * n * fnr(np1 - m) * fnr(np1 - mm1)
            vwh_coef(m, n, -1, -1) = -fnorm1 * np1 * fnr(n + m) * fnr(nm1 + m)
            vwh_coef(m, n, 0, 1) = fnorm1 * n * fnr(np1 + m) * fnr(np1 - m)
            vwh_coef(m, n, 0, -1) = fnorm1 * np1 * fnr(n + m) * fnr(n - m)
            vwh_coef(m, n, 1, 0) = -fnorm2 * fnr(n - m) * fnr(np1 + m)
            vwh_coef(m, n, -1, 0) = -fnorm2 * fnr(n + m) * fnr(np1 - m)
            vwh_coef(m, n, 0, 0) = -fnorm2 * m
         end do
         do m = 0, n - 1
            mp1 = m + 1
            mm1 = m - 1
            vwh_coef(m, n, 1, 1) = -fnorm1 * n * fnr(np1 + m) * fnr(np1 + mp1)
            vwh_coef(m, n, 1, -1) = fnorm1 * np1 * fnr(n - m) * fnr(nm1 - m)
            vwh_coef(m, n, -1, 1) = fnorm1 * n * fnr(np1 - m) * fnr(np1 - mm1)
            vwh_coef(m, n, -1, -1) = -fnorm1 * np1 * fnr(n + m) * fnr(nm1 + m)
            vwh_coef(m, n, 0, 1) = fnorm1 * n * fnr(np1 + m) * fnr(np1 - m)
            vwh_coef(m, n, 0, -1) = fnorm1 * np1 * fnr(n + m) * fnr(n - m)
            vwh_coef(m, n, 1, 0) = -fnorm2 * fnr(n - m) * fnr(np1 + m)
            vwh_coef(m, n, -1, 0) = -fnorm2 * fnr(n + m) * fnr(np1 - m)
            vwh_coef(m, n, 0, 0) = -fnorm2 * m
         end do
         m = n
         mp1 = m + 1
         mm1 = m - 1
         vwh_coef(m, n, 1, 1) = -fnorm1 * n * fnr(np1 + m) * fnr(np1 + mp1)
         vwh_coef(m, n, 1, -1) = 0.d0
         vwh_coef(m, n, -1, 1) = fnorm1 * n * fnr(np1 - m) * fnr(np1 - mm1)
         vwh_coef(m, n, -1, -1) = -fnorm1 * np1 * fnr(n + m) * fnr(nm1 + m)
         vwh_coef(m, n, 0, 1) = fnorm1 * n * fnr(np1 + m) * fnr(np1 - m)
         vwh_coef(m, n, 0, -1) = 0.d0
         vwh_coef(m, n, 1, 0) = -0.d0
         vwh_coef(m, n, -1, 0) = -fnorm2 * fnr(n + m) * fnr(np1 - m)
         vwh_coef(m, n, 0, 0) = -fnorm2 * m
      end do
   end subroutine init

end module numconstants
