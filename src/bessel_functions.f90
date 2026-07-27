module bessel_functions
   use constants
   implicit none
contains

!
!  Riccati-Bessel function psi(n), complex argument
!
!
!  last revised: 15 January 2011
!
   pure subroutine riccati_bessel(n, ds, psi)
      implicit none
      integer, intent(in) :: n
      complex(8), intent(in) :: ds
      complex(8), intent(out) :: psi(0:n)
      integer :: i
      complex(8) :: chi(0:n)
      call complex_spherical_bessel(n, ds, psi, chi)
      do concurrent(i=0:n)
         psi(i) = psi(i) * ds
      end do
      return
   end subroutine riccati_bessel
!
!  Riccati-Hankel function xi(n), complex argument
!
!
!  last revised: 15 January 2011
!  March 2013
!  The condition abs(xi(i))/abs(psi(0)) << 1.d-6
!  implies an argument with large imag part, and use of xi = psi + i chi will have
!  round off problems.   Upwards recurrence is used in this case.
!
   pure subroutine riccati_hankel(n, ds, xi)
      implicit none
      integer, intent(in) :: n
      complex(8), intent(in) :: ds
      complex(8), intent(out) :: xi(0:n)
      complex(8), parameter :: ci = (0.d0, 1.d0)
      integer :: i
      complex(8) :: psi(0:n), chi(0:n), psi0
      xi(0) = -ci * exp(ci * ds)
      psi0 = sin(ds)
      if (abs(xi(0)) / abs(psi0) .lt. 1.d-6) then
         xi(1) = -exp(ci * ds) * (ci + ds) / ds
         do i = 2, n
            xi(i) = dble(i + i + 1) / ds * xi(i - 1) - xi(i - 2)
         end do
      else
         call complex_spherical_bessel(n, ds, psi, chi)
         do concurrent(i=1:n)
            xi(i) = (psi(i) + ci * chi(i)) * ds
         end do
      end if
   end subroutine riccati_hankel
!
!     ==========================================================
!     Purpose: Compute spherical Bessel functions jn(z) & yn(z)
!              for a complex argument
!     Input :  z --- Complex argument
!              n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
!     Output:  CSJ(n) --- jn(z)
!              CSY(n) --- yn(z)
!              NM --- Highest order computed
!     Routines called:
!              bessel_recurrence_start and bessel_recurrence_start_for_order
!              for computing the starting point for backward recurrence
!     ==========================================================
!
!    obtained from, and copywrited by, Jian-Ming Jin
!    http://jin.ece.uiuc.edu/
!
!
!  last revised: 15 January 2011
!
   pure subroutine complex_spherical_bessel(n, z, csj, csy)
      implicit none
      integer, intent(in) :: n
      complex(8), intent(in) :: z
      complex(8), intent(out) :: csj(0:n), csy(0:n)
      integer :: nm, k, m
      real(8) :: a0
      complex(8) :: csa, csb, cs, cf0, cf1, cf
      a0 = abs(z)
      nm = n
      if (a0 .lt. 1.0d-60) then
         csj = (0.d0, 0.d0)
         csy = (-1.d300, 0.d0)
         csy(0) = (1.d0, 0.d0)
         return
      end if
      csj = (0.d0, 0.d0)
      csj(0) = sin(z) / z
      csj(1) = (csj(0) - cos(z)) / z
      if (n .ge. 2) then
         csa = csj(0)
         csb = csj(1)
         m = bessel_recurrence_start(a0, 200)
         if (m .lt. n) then
            nm = m
         else
            m = bessel_recurrence_start_for_order(a0, n, 15)
         end if
         cf0 = 0.0d0
         cf1 = 1.0d0 - 100
         do k = m, 0, -1
            cf = (2.0d0 * k + 3.0d0) * cf1 / z - cf0
            if (k .le. nm) csj(k) = cf
            cf0 = cf1
            cf1 = cf
         end do
         if (abs(csa) .gt. abs(csb)) cs = csa / cf
         if (abs(csa) .le. abs(csb)) cs = csb / cf0
         do concurrent(k=0:min(nm, n))
            csj(k) = cs * csj(k)
         end do
      end if
      csy = (1.d200, 0.d0)
      csy(0) = -cos(z) / z
      csy(1) = (csy(0) - sin(z)) / z
      do k = 2, min(nm, n)
         if (abs(csj(k - 1)) .gt. abs(csj(k - 2))) then
            csy(k) = (csj(k) * csy(k - 1) - 1.0d0 / (z * z)) / csj(k - 1)
         else
            csy(k) = (csj(k) * csy(k - 2) - (2.0d0 * k - 1.0d0) / z**3) / csj(k - 2)
         end if
      end do
   end subroutine complex_spherical_bessel

   pure subroutine bessel_integer_complex(n, z, nmax, b)
      implicit none
      integer, intent(in) :: n
      complex(8), intent(in) :: z
      integer, intent(out) :: nmax
      complex(8), intent(out) :: b(0:n)
      complex(8) :: cbj(0:n + 1), cdj(0:n + 1), cby(0:n + 1), cdy(0:n + 1)

      if (aimag(z) .eq. 0.d0) then
         b = bessel_jn(0, n, real(z, kind=kind(0.d0)))
         nmax = n
         return
      end if

      call complex_bessel_jy(n, z, nmax, cbj, cdj, cby, cdy)
      nmax = min(n, nmax)
      b(0:nmax) = cbj(0:nmax)
   end subroutine bessel_integer_complex

   pure subroutine complex_bessel_jy(n, z, nm, cbj, cdj, cby, cdy)
      !*****************************************************************************80
      !
         !! Bessel functions, derivatives, Jn(z) and Yn(z) of complex argument.
      !
      !  Licensing:
      !
      !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
      !    they give permission to incorporate this routine into a user program
      !    provided that the copyright is acknowledged.
      !
      !  Modified:
      !
      !    03 August 2012
      !
      !  Author:
      !
      !    Shanjie Zhang, Jianming Jin
      !
      !  Reference:
      !
      !    Shanjie Zhang, Jianming Jin,
      !    Computation of Special Functions,
      !    Wiley, 1996,
      !    ISBN: 0-471-11963-6,
      !    LC: QA351.C45.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) N, the order of Jn(z) and Yn(z).
      !
      !    Input, complex ( kind = 8 ) Z, the argument of Jn(z) and Yn(z).
      !
      !    Output, integer ( kind = 4 ) NM, the highest order computed.
      !
      !    Output, complex ( kind = 8 ) CBJ(0:N), CDJ(0:N), CBY(0:N), CDY(0:N),
      !    the values of Jn(z), Jn'(z), Yn(z), Yn'(z).
      !
      implicit none
      integer(kind=4), intent(in) :: n
      complex(kind=8), intent(in) :: z
      integer(kind=4), intent(out) :: nm
      complex(kind=8), intent(out) :: cbj(0:n + 1), cby(0:n + 1), cdj(0:n + 1), cdy(0:n + 1)
      real(kind=8), parameter, dimension(4) :: a = (/ &
                                               -0.7031250000000000D-01, 0.1121520996093750D+00, &
                                               -0.5725014209747314D+00, 0.6074042001273483D+01/)
      real(kind=8) a0
      real(kind=8), parameter, dimension(4) :: a1 = (/ &
                                               0.1171875000000000D+00, -0.1441955566406250D+00, &
                                               0.6765925884246826D+00, -0.6883914268109947D+01/)
      real(kind=8), parameter, dimension(4) :: b = (/ &
                                               0.7324218750000000D-01, -0.2271080017089844D+00, &
                                               0.1727727502584457D+01, -0.2438052969955606D+02/)
      real(kind=8), parameter, dimension(4) :: b1 = (/ &
                                               -0.1025390625000000D+00, 0.2775764465332031D+00, &
                                               -0.1993531733751297D+01, 0.2724882731126854D+02/)
      complex(kind=8) cbj0
      complex(kind=8) cbj1
      complex(kind=8) cbjk
      complex(kind=8) cbs
      complex(kind=8) cby0
      complex(kind=8) cby1
      complex(kind=8) ce
      complex(kind=8) cf
      complex(kind=8) cf1
      complex(kind=8) cf2
      complex(kind=8) cp0
      complex(kind=8) cp1
      complex(kind=8) cq0
      complex(kind=8) cq1
      complex(kind=8) cs0
      complex(kind=8) csu
      complex(kind=8) csv
      complex(kind=8) ct1
      complex(kind=8) ct2
      complex(kind=8) cu
      complex(kind=8) cyy
      real(kind=8) el
      integer(kind=4) k
      integer(kind=4) m
      real(kind=8) r2p
      real(kind=8) y0

      el = 0.5772156649015329D+00
      r2p = 0.63661977236758D+00
      y0 = abs(aimag(z))
      a0 = abs(z)
      nm = n
      if (a0 < 1.0D-100) then
         do concurrent(k=0:n)
            cbj(k) = cmplx(0.0D+00, 0.0D+00, kind=8)
            cdj(k) = cmplx(0.0D+00, 0.0D+00, kind=8)
            cby(k) = -cmplx(1.0D+30, 0.0D+00, kind=8)
            cdy(k) = cmplx(1.0D+30, 0.0D+00, kind=8)
         end do
         cbj(0) = cmplx(1.0D+00, 0.0D+00, kind=8)
         cdj(1) = cmplx(0.5D+00, 0.0D+00, kind=8)
         return
      end if
      if (a0 <= 300.0D+00 .or. 80 < n) then
         if (n == 0) then
            nm = 1
         end if
         m = bessel_recurrence_start(a0, 200)
         if (m < nm) then
            nm = m
         else
            m = bessel_recurrence_start_for_order(a0, nm, 15)
         end if
         cbs = cmplx(0.0D+00, 0.0D+00, kind=8)
         csu = cmplx(0.0D+00, 0.0D+00, kind=8)
         csv = cmplx(0.0D+00, 0.0D+00, kind=8)
         cf2 = cmplx(0.0D+00, 0.0D+00, kind=8)
         cf1 = cmplx(1.0D-30, 0.0D+00, kind=8)
         do k = m, 0, -1
            cf = 2.0D+00 * (k + 1.0D+00) / z * cf1 - cf2
            if (k <= nm) then
               cbj(k) = cf
            end if
            if (k == 2 * int(k / 2) .and. k .ne. 0) then
               if (y0 <= 1.0D+00) then
                  cbs = cbs + 2.0D+00 * cf
               else
                  cbs = cbs + (-1.0D+00)**(k / 2) * 2.0D+00 * cf
               end if
               csu = csu + (-1.0D+00)**(k / 2) * cf / k
            else if (1 < k) then
               csv = csv + (-1.0D+00)**(k / 2) * k / (k * k - 1.0D+00) * cf
            end if
            cf2 = cf1
            cf1 = cf
         end do
         if (y0 <= 1.0D+00) then
            cs0 = cbs + cf
         else
            cs0 = (cbs + cf) / cos(z)
         end if
         do concurrent(k=0:nm)
            cbj(k) = cbj(k) / cs0
         end do
         ce = log(z / 2.0D+00) + el
         cby(0) = r2p * (ce * cbj(0) - 4.0D+00 * csu / cs0)
         cby(1) = r2p * (-cbj(0) / z + (ce - 1.0D+00) * cbj(1) &
                         - 4.0D+00 * csv / cs0)
      else
         ct1 = z - quarter_pi
         cp0 = cmplx(1.0D+00, 0.0D+00, kind=8)
         do k = 1, 4
            cp0 = cp0 + a(k) * z**(-2 * k)
         end do
         cq0 = -0.125D+00 / z
         do k = 1, 4
            cq0 = cq0 + b(k) * z**(-2 * k - 1)
         end do
         cu = sqrt(r2p / z)
         cbj0 = cu * (cp0 * cos(ct1) - cq0 * sin(ct1))
         cby0 = cu * (cp0 * sin(ct1) + cq0 * cos(ct1))
         cbj(0) = cbj0
         cby(0) = cby0
         ct2 = z - three_quarters_pi
         cp1 = cmplx(1.0D+00, 0.0D+00, kind=8)
         do k = 1, 4
            cp1 = cp1 + a1(k) * z**(-2 * k)
         end do
         cq1 = 0.375D+00 / z
         do k = 1, 4
            cq1 = cq1 + b1(k) * z**(-2 * k - 1)
         end do
         cbj1 = cu * (cp1 * cos(ct2) - cq1 * sin(ct2))
         cby1 = cu * (cp1 * sin(ct2) + cq1 * cos(ct2))
         cbj(1) = cbj1
         cby(1) = cby1
         do k = 2, nm
            cbjk = 2.0D+00 * (k - 1.0D+00) / z * cbj1 - cbj0
            cbj(k) = cbjk
            cbj0 = cbj1
            cbj1 = cbjk
         end do
      end if
      cdj(0) = -cbj(1)
      do concurrent(k=1:nm)
         cdj(k) = cbj(k - 1) - k / z * cbj(k)
      end do
      if (1.0D+00 < abs(cbj(0))) then
         cby(1) = (cbj(1) * cby(0) - 2.0D+00 / (pi * z)) / cbj(0)
      end if
      do k = 2, nm
         if (abs(cbj(k - 2)) <= abs(cbj(k - 1))) then
            cyy = (cbj(k) * cby(k - 1) - 2.0D+00 / (pi * z)) / cbj(k - 1)
         else
            cyy = (cbj(k) * cby(k - 2) - 4.0D+00 * (k - 1.0D+00) &
                   / (pi * z * z)) / cbj(k - 2)
         end if
         cby(k) = cyy
      end do
      cdy(0) = -cby(1)
      do concurrent(k=1:nm)
         cdy(k) = cby(k - 1) - k / z * cby(k)
      end do
      return
   end subroutine complex_bessel_jy

   pure integer function bessel_recurrence_start(x, mp)
      implicit none
      real(8), intent(in) :: x
      integer, intent(in) :: mp
      integer :: n0, n1, it, nn
      real(8) :: a0, f1, f, f0
      a0 = dabs(x)
      n0 = int(1.1 * a0) + 1
      f0 = bessel_order_envelope(n0, a0) - mp
      n1 = n0 + 5
      f1 = bessel_order_envelope(n1, a0) - mp
      do it = 1, 20
         nn = n1 - (n1 - n0) / (1.0d0 - f0 / f1)
         f = bessel_order_envelope(nn, a0) - mp
         if (abs(nn - n1) .lt. 1) exit
         n0 = n1
         f0 = f1
         n1 = nn
         f1 = f
      end do
      bessel_recurrence_start = nn
   end function bessel_recurrence_start

   pure integer function bessel_recurrence_start_for_order(x, n, mp)
      implicit none
      real(8), intent(in) :: x
      integer, intent(in) :: n, mp
      integer :: n0, n1, it, nn
      real(8) :: a0, hmp, ejn, obj, f0, f1, f
      a0 = dabs(x)
      hmp = 0.5d0 * dble(mp)
      ejn = bessel_order_envelope(n, a0)
      if (ejn .le. hmp) then
         obj = mp
         n0 = int(1.1 * a0)
      else
         obj = hmp + ejn
         n0 = n
      end if
      f0 = bessel_order_envelope(n0, a0) - obj
      n1 = n0 + 5
      f1 = bessel_order_envelope(n1, a0) - obj
      do it = 1, 20
         nn = n1 - (n1 - n0) / (1.0d0 - f0 / f1)
         f = bessel_order_envelope(nn, a0) - obj
         if (abs(nn - n1) .lt. 1) exit
         n0 = n1
         f0 = f1
         n1 = nn
         f1 = f
      end do
      bessel_recurrence_start_for_order = nn + 10
   end function bessel_recurrence_start_for_order

   pure real(8) function bessel_order_envelope(n, x)
      implicit none
      integer, intent(in) :: n
      real(8), intent(in) :: x
      integer :: order
      order = max(1, abs(n))
      bessel_order_envelope = 0.5d0 * dlog10(6.28d0 * order) - order * dlog10(1.36d0 * x / order)
   end function bessel_order_envelope
end module bessel_functions
