module wave_functions
   use constants
   use angular_functions, only: azimuthal_phase_factors, cartesian_to_spherical, &
                                complex_rotation_coefficients, rotation_coefficients
   use bessel_functions, only: riccati_bessel, riccati_hankel
   use coefficient_indexing, only: polarized_mode_index
   implicit none
contains

   subroutine vector_spherical_wave_functions(rpos, ri, nodr, itype, vwh, index_model, lr_to_mode)
      use numerical_tables
      implicit none
      logical, optional :: lr_to_mode
      logical :: lrtomode
      integer :: nodr, itype, n, nodrp1, nodrm1, nn1, np1, nm1, p, sp, imod, m, iadd(-nodr:nodr)
      integer, save :: nodrmax
      integer, optional :: index_model
      real(8) :: rpos(3), r, ct
      real(8) pmn(0:0, 0:(nodr + 1) * (nodr + 3))
      complex(8) :: ci, vwh(3, 2 * nodr * (nodr + 2)), ri(2), ephi, a(2), vtemp(3, 2)
      complex(8) :: a1vec(-nodr:nodr), &
                    b1vec(-nodr:nodr), z1vec(-nodr:nodr), a2vec(-nodr:nodr), &
                    b2vec(-nodr:nodr), z2vec(-nodr:nodr)
      complex(8) :: umn(-nodr - 2:nodr + 2, 0:nodr + 1, 2), hn(0:nodr + 1, 2), ephim(-nodr - 1:nodr + 1)
      data ci, nodrmax/(0.d0, 1.d0), 0/
      if (nodr .gt. nodrmax) then
         nodrmax = nodr
         call initialize_numerical_tables(nodr + 2)
      end if
      if (present(index_model)) then
         imod = index_model
      else
         imod = 1
      end if
      if (present(lr_to_mode)) then
         lrtomode = lr_to_mode
      else
         lrtomode = .false.
      end if
      call cartesian_to_spherical(rpos, r, ct, ephi)
      if (r .le. 1.d-5) then
         vwh(:, 1:2 * nodr * (nodr + 2)) = (0.d0, 0.d0)
         if (itype .eq. 3) return
         do p = 1, 2
            vwh(1, polarized_mode_index(-1, 1, p, nodr, imod)) = .5d0 * fnr(2) / fnr(3)
            vwh(2, polarized_mode_index(-1, 1, p, nodr, imod)) = -.5d0 * ci * fnr(2) / fnr(3)
            vwh(3, polarized_mode_index(0, 1, p, nodr, imod)) = 1.d0 * fnr(2) / fnr(6)
            vwh(1, polarized_mode_index(1, 1, p, nodr, imod)) = -.5d0 * fnr(2) / fnr(3)
            vwh(2, polarized_mode_index(1, 1, p, nodr, imod)) = -.5d0 * ci * fnr(2) / fnr(3)
            if (lrtomode) exit
         end do
         return
      end if
      nodrp1 = nodr + 1
      nodrm1 = nodr - 1
!
! this is now a vector operation w/ l/r form
!
      a = ri * r
      do p = 1, 2
         if (itype .eq. 1) then
            call riccati_bessel(nodrp1, a(p), hn(0, p))
         else
            call riccati_hankel(nodrp1, a(p), hn(0, p))
         end if
         hn(0:nodrp1, p) = hn(0:nodrp1, p) / a(p)
         if (a(2) .eq. a(1)) then
            hn(0:nodrp1, 2) = hn(0:nodrp1, 1)
            exit
         end if
      end do
      call rotation_coefficients(ct, 0, nodrp1, pmn)
      call azimuthal_phase_factors(ephi, nodrp1, ephim)
      umn = 0.d0
      do p = 1, 2
         umn(0, 0, p) = hn(0, p) * fnr(2)
         do n = 1, nodrp1
            nn1 = n * (n + 1)
            umn(-n:n, n, p) = fnr(2) * pmn(0, nn1 - n:nn1 + n) * ephim(-n:n) * hn(n, p)
            umn(-n - 1, n, p) = 0.d0
            umn(n + 1, n, p) = 0.d0
         end do
      end do
      do p = 1, 2
         sp = -(-1)**p
         do n = 1, nodr
            do m = -n, n
               iadd(m) = polarized_mode_index(m, n, p, nodr, imod)
            end do
            nn1 = n * (n + 1)
            np1 = n + 1
            nm1 = n - 1
            a1vec(-n:n) = vwh_coef(-n:n, n, 1, 1) * umn(-nm1:np1, np1, p) &
                          + vwh_coef(-n:n, n, 1, -1) * umn(-nm1:np1, nm1, p)
            b1vec(-n:n) = vwh_coef(-n:n, n, -1, 1) * umn(-np1:nm1, np1, p) &
                          + vwh_coef(-n:n, n, -1, -1) * umn(-np1:nm1, nm1, p)
            z1vec(-n:n) = vwh_coef(-n:n, n, 0, 1) * umn(-n:n, np1, p) &
                          + vwh_coef(-n:n, n, 0, -1) * umn(-n:n, nm1, p)
            a2vec(-n:n) = vwh_coef(-n:n, n, 1, 0) * umn(-nm1:np1, n, p)
            b2vec(-n:n) = vwh_coef(-n:n, n, -1, 0) * umn(-np1:nm1, n, p)
            z2vec(-n:n) = vwh_coef(-n:n, n, 0, 0) * umn(-n:n, n, p)
            vwh(1, iadd(-n:n)) = -0.5d0 * (a1vec(-n:n) + b1vec(-n:n)) &
                                 - sp * 0.5d0 * ci * (a2vec(-n:n) + b2vec(-n:n))
            vwh(2, iadd(-n:n)) = -0.5d0 * ci * (-a1vec(-n:n) + b1vec(-n:n)) &
                                 - sp * 0.5d0 * (a2vec(-n:n) - b2vec(-n:n))
            vwh(3, iadd(-n:n)) = -z1vec(-n:n) &
                                 - sp * ci * z2vec(-n:n)
         end do
      end do
      if (lrtomode) then
         do n = 1, nodr
            do m = -n, n
               do p = 1, 2
                  vtemp(:, p) = vwh(:, polarized_mode_index(m, n, p, nodr, imod))
               end do
               vwh(:, polarized_mode_index(m, n, 1, nodr, imod)) = (vtemp(:, 1) + vtemp(:, 2)) * 0.5d0
               vwh(:, polarized_mode_index(m, n, 2, nodr, imod)) = (vtemp(:, 1) - vtemp(:, 2)) * 0.5d0
            end do
         end do
      end if
   end subroutine vector_spherical_wave_functions

   subroutine scalar_wave_function(nodr, itype, x, y, z, ri, swf)
      use numerical_tables
      implicit none
      integer :: nodr, itype, n, m, mn
      real(8) :: x, y, z, r, ct, rho, ymn(0:nodr * (nodr + 2)), c, c0
      complex(8) :: ri, swf(0:nodr * (nodr + 2)), ephi, rri, hn(0:nodr)

      r = sqrt(x * x + y * y + z * z)
      if (r .lt. 1.d-10) then
         swf = 0.d0
         if (itype .eq. 1) swf(0) = 1.d0
         return
      end if
      ct = z / r
      if (x .eq. 0.d0 .and. y .eq. 0.d0) then
         ephi = 1.d0
         rho = 0.d0
      else
         rho = sqrt(x * x + y * y)
         ephi = cmplx(x, y, kind=kind(0.0d0)) / rho
      end if
      call rotation_coefficients(ct, 0, nodr, ymn)
      rri = r * ri
      if (itype .eq. 3) then
         hn(0) = -(0.d0, 1.d0) * exp((0.d0, 1.d0) * rri) / rri
         hn(1) = -exp((0.d0, 1.d0) * rri) * ((0.d0, 1.d0) + rri) / rri / rri
         do n = 2, nodr
            hn(n) = dble(n + n - 1) / rri * hn(n - 1) - hn(n - 2)
         end do
      else
         call riccati_bessel(nodr, rri, hn)
         hn = hn / rri
      end if

      c0 = 1.d0 / sqrt(four_pi)
      do n = 0, nodr
         c = c0 * sqrt(dble(n + n + 1))
         do m = -n, n
            mn = n * (n + 1) + m
            swf(mn) = hn(n) * ymn(mn) * ephi**m * c
         end do
      end do
   end subroutine scalar_wave_function

   subroutine reciprocal_scalar_wave_function(nodr, kx, ky, x, y, z, ri, swf)
      use numerical_tables
      implicit none
      integer :: nodr, n, m, mn
      real(8) :: kx, ky, x, y, z, k
      complex(8) :: ri, swf(0:nodr * (nodr + 2)), ephi, kz, ymn(0:nodr * (nodr + 2)), skz, c, cr
      k = sqrt(kx * kx + ky * ky)
      kz = sqrt((1.d0, 0.d0) - k * k / ri / ri)
      if (z .gt. 0.d0) then
         skz = kz
      else
         skz = -kz
      end if
      if (k .eq. 0.d0) then
         ephi = 1.d0
      else
         ephi = cmplx(kx, ky, kind=kind(0.0d0)) / k
      end if
      call complex_rotation_coefficients(skz, 0, nodr, ymn)
      c = exp((0.d0, 1.d0) * (kx * x + ky * y + ri * skz * z)) / ri / ri / kz / sqrt(four_pi)
      do n = 0, nodr
         cr = ((0.d0, -1.d0))**n * sqrt(dble(n + n + 1))
         do m = -n, n
            mn = n * (n + 1) + m
            swf(mn) = cr * c * ymn(mn) * ephi**m
         end do
      end do
   end subroutine reciprocal_scalar_wave_function
!
! inverse of a 2 X 2 complex matrix.
! March 2013
!
   pure subroutine invert_two_by_two_matrix(mat, imat)
      implicit none
      complex(8), intent(in) :: mat(2, 2)
      complex(8), intent(out) :: imat(2, 2)
      integer :: s, t
      complex(8) :: tmat(2, 2), det
      tmat = mat
      det = mat(1, 1) * mat(2, 2) - mat(2, 1) * mat(1, 2)
      do concurrent(s=1:2, t=1:2)
         imat(s, t) = (-1)**(s + t) * tmat(3 - t, 3 - s) / det
      end do
   end subroutine invert_two_by_two_matrix
!
! move between unequal size matrices.
! March 2013
!
   pure subroutine resize_mode_coefficients(nin, nout, cin, cout)
      implicit none
      integer, intent(in) :: nin, nout
      complex(8), intent(in) :: cin(0:nin + 1, nin, 2)
      complex(8), intent(out) :: cout(0:nout + 1, nout, 2)
      complex(8) :: ct(0:max(nin, nout) + 1, max(nin, nout), 2)
      ct = 0.d0
      ct(0:nin + 1, 1:nin, 1:2) = cin(0:nin + 1, 1:nin, 1:2)
      cout = 0.d0
      cout(0:nout + 1, 1:nout, 1:2) = ct(0:nout + 1, 1:nout, 1:2)
   end subroutine resize_mode_coefficients
   subroutine left_right_mode_transformation(nodr, alr, amode, lr_to_mode)
      implicit none
      logical :: lrtomode
      logical, optional :: lr_to_mode
      integer :: nodr
      complex(8) :: alr(nodr * (nodr + 2), 2), amode(nodr * (nodr + 2), 2), at(nodr * (nodr + 2), 2)
      if (present(lr_to_mode)) then
         lrtomode = lr_to_mode
      else
         lrtomode = .true.
      end if
      if (lrtomode) then
         at = alr(:, :)
         amode(:, 1) = at(:, 1) + at(:, 2)
         amode(:, 2) = at(:, 1) - at(:, 2)
      else
         at = amode(:, :)
         alr(:, 1) = .5d0 * (at(:, 1) + at(:, 2))
         alr(:, 2) = .5d0 * (at(:, 1) - at(:, 2))
      end if
   end subroutine left_right_mode_transformation

   pure subroutine reverse_azimuthal_modes(nodr, ain, aout)
      implicit none
      integer, intent(in) :: nodr
      complex(8), intent(in) :: ain(2 * nodr * (nodr + 2))
      complex(8), intent(out) :: aout(2 * nodr * (nodr + 2))
      integer :: m, n, p, mnp, mnp2, im, m1
      do concurrent(m=-nodr:nodr) local(m1, im, mnp, mnp2)
         m1 = max(abs(m), 1)
         im = (-1)**m
         do n = m1, nodr
            do p = 1, 2
               mnp = polarized_mode_index(m, n, p, nodr, 2)
               mnp2 = polarized_mode_index(-m, n, p, nodr, 2)
               aout(mnp2) = im * ain(mnp)
            end do
         end do
      end do
   end subroutine reverse_azimuthal_modes

   subroutine compose_group_filename(firststring, number, laststring, newstring)
      implicit none
      integer :: number
      character(len=256) :: firststring, laststring, newstring, sform, intfile
      if (number .lt. 10) then
         sform = '(a,i1,a,a)'
      elseif (number .lt. 100) then
         sform = '(a,i2,a,a)'
      elseif (number .lt. 1000) then
         sform = '(a,i3,a,a)'
      else
         sform = '(a,i4,a,a)'
      end if
      write (intfile, fmt=sform) trim(firststring), number, '_', trim(laststring)
      read (intfile, '(a)') newstring
   end subroutine compose_group_filename
end module wave_functions
