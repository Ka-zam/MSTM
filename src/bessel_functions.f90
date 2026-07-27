module bessel_functions
   implicit none
contains

   subroutine ricbessel(n, ds, eps, nmax, psi)
      implicit none
      integer :: n, nmax, ns, i
      real(8) :: ds, dns, sn, psi(0:n), psit, ds2, sum, eps, err
      if (int(ds) .lt. n) then
         ns = nint(ds + 4.*(ds**.3333d0) + 17)
         ns = max(n + 10, ns)
         dns = 0.d0
         do i = ns - 1, n, -1
            sn = dble(i + 1) / ds
            dns = sn - 1.d0 / (dns + sn)
         end do
         psi(n) = dns
         psi(n - 1) = dble(n) / ds - 1.d0 / (dns + dble(n) / ds)
         do i = n - 2, 1, -1
            sn = dble(i + 1) / ds
            psi(i) = sn - 1.d0 / (psi(i + 1) + sn)
         end do
         psit = dsin(ds)
         psi(0) = psit
         ds2 = ds * ds
         sum = psit * psit / ds2
         do i = 1, n
            psit = psit / (dble(i) / ds + psi(i))
            sum = sum + dble(i + i + 1) * psit * psit / ds2
            err = dabs(1.d0 - sum)
            psi(i) = psit
            if (err .lt. eps) then
               nmax = i
               return
            end if
         end do
         nmax = n
      else
         psi(0) = dsin(ds)
         psi(1) = psi(0) / ds - dcos(ds)
         do i = 1, n - 1
            sn = dble(i + i + 1) / ds
            psi(i + 1) = sn * psi(i) - psi(i - 1)
         end do
         nmax = n
      end if
   end subroutine ricbessel
!
!  ricatti-hankel function xi(n), real argument
!
!
!  last revised: 15 January 2011
!
   subroutine richankel(n, ds, xi)
      implicit none
      integer :: n, i, ns
      real(8) :: ds, dns, sn, chi0, chi1, chi2, psi, psi0, psi1
      complex(8) :: xi(0:n)
      if (int(ds) .lt. n) then
         ns = nint(ds + 4.*(ds**.3333) + 17)
         ns = max(n + 10, ns)
         dns = 0.d0
         do i = ns - 1, n, -1
            sn = dble(i + 1) / ds
            dns = sn - 1.d0 / (dns + sn)
         end do
         xi(n) = dns
         xi(n - 1) = dble(n) / ds - 1.d0 / (dns + dble(n) / ds)
         do i = n - 2, 1, -1
            sn = dble(i + 1) / ds
            xi(i) = sn - 1.d0 / (xi(i + 1) + sn)
         end do
         chi0 = -dcos(ds)
         psi = dsin(ds)
         chi1 = chi0 / ds - psi
         xi(0) = cmplx(psi, chi0, kind=kind(0.0d0))
         do i = 1, n
            chi2 = dble(i + i + 1) / ds * chi1 - chi0
            psi = psi / (dble(i) / ds + xi(i))
            xi(i) = cmplx(psi, chi1, kind=kind(0.0d0))
            chi0 = chi1
            chi1 = chi2
         end do
         return
      else
         chi0 = -dcos(ds)
         psi0 = dsin(ds)
         chi1 = chi0 / ds - psi0
         psi1 = psi0 / ds + chi0
         xi(0) = cmplx(psi0, chi0, kind=kind(0.0d0))
         xi(1) = cmplx(psi1, chi1, kind=kind(0.0d0))
         do i = 1, n - 1
            sn = dble(i + i + 1) / ds
            xi(i + 1) = sn * xi(i) - xi(i - 1)
         end do
         return
      end if
   end subroutine richankel
!
!  ricatti-bessel function psi(n), complex argument
!
!
!  last revised: 15 January 2011
!
   subroutine cricbessel(n, ds, psi)
      implicit none
      integer :: n, i
      complex(8) :: ds, psi(0:n), chi(0:n)
      call cspherebessel(n, ds, psi, chi)
      do i = 0, n
         psi(i) = psi(i) * ds
      end do
      return
   end subroutine cricbessel
!
!  ricatti-hankel function psi(n), complex argument
!
!
!  last revised: 15 January 2011
!  March 2013
!  The condition abs(xi(i))/abs(psi(0)) << 1.d-6
!  implies an argument with large imag part, and use of xi = psi + i chi will have
!  round off problems.   Upwards recurrence is used in this case.
!
   subroutine crichankel(n, ds, xi)
      implicit none
      integer :: n, i
      complex(8) :: ds, psi(0:n), chi(0:n), xi(0:n), ci, &
                    psi0
      data ci/(0.d0, 1.d0)/
      xi(0) = -ci * exp(ci * ds)
      psi0 = sin(ds)
      if (abs(xi(0)) / abs(psi0) .lt. 1.d-6) then
         xi(1) = -exp(ci * ds) * (ci + ds) / ds
         do i = 2, n
            xi(i) = dble(i + i + 1) / ds * xi(i - 1) - xi(i - 2)
         end do
      else
         call cspherebessel(n, ds, psi, chi)
         do i = 1, n
            xi(i) = (psi(i) + ci * chi(i)) * ds
         end do
      end if
   end subroutine crichankel
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
!              MSTA1 and MSTA2 for computing the starting
!              point for backward recurrence
!     ==========================================================
!
!    obtained from, and copywrited by, Jian-Ming Jin
!    http://jin.ece.uiuc.edu/
!
!
!  last revised: 15 January 2011
!
   subroutine cspherebessel(n, z, csj, csy)
      implicit none
      integer :: n, nm, k, m
      real(8) :: a0
      complex(8) :: z, csj(0:n), csy(0:n), csa, csb, cs, cf0, cf1, cf
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
         m = msta1(a0, 200)
         if (m .lt. n) then
            nm = m
         else
            m = msta2(a0, n, 15)
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
         do k = 0, min(nm, n)
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
   end subroutine cspherebessel

   subroutine ch12n(n, z, nm, chf1)
      use numconstants

      !*****************************************************************************80
      !
         !! CH12N computes Hankel functions of first and second kinds, complex argument.
      !
      !  Discussion:
      !
      !    Both the Hankel functions and their derivatives are computed.
      !
      !  Licensing:
      !
      !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
      !    they give permission to incorporate this routine into a user program
      !    provided that the copyright is acknowledged.
      !
      !  Modified:
      !
      !    26 July 2012
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
      !    Input, integer ( kind = 4 ) N, the order of the functions.
      !
      !    Input, complex ( kind = 8 ) Z, the argument.
      !
      !    Output, integer ( kind = 4 ) NM, the highest order computed.
      !
      !    Output, complex ( kind = 8 ) CHF1(0:n), CHD1(0:n), CHF2(0:n), CHD2(0:n),
      !    the values of Hn(1)(z), Hn(1)'(z), Hn(2)(z), Hn(2)'(z).
      !
      implicit none

      integer(kind=4) n

      complex(kind=8) cbi(0:n + 1)
      complex(kind=8) cbj(0:n + 1)
      complex(kind=8) cbk(0:n + 1)
      complex(kind=8) cby(0:n + 1)
      complex(kind=8) cdi(0:n + 1)
      complex(kind=8) cdj(0:n + 1)
      complex(kind=8) cdk(0:n + 1)
      complex(kind=8) cdy(0:n + 1)
      complex(kind=8) cf1
      complex(kind=8) cfac
      complex(kind=8) chf1(0:n)
      complex(kind=8) ci
      integer(kind=4) k
      integer(kind=4) nm
      complex(kind=8) z
      complex(kind=8) zi

      ci = cmplx(0.0D+00, 1.0D+00, kind=8)
      if (aimag(z) .le. 0.0D+00) then
         call cjynb(n, z, nm, cbj, cdj, cby, cdy)
         nm = min(n, nm)
         do k = 0, nm
            chf1(k) = cbj(k) + ci * cby(k)
         end do
      else
         zi = -ci * z
         call ciknb(n, zi, nm, cbi, cdi, cbk, cdk)
         cf1 = -ci
         cfac = 2.0D+00 / (pi * ci)
         nm = min(n, nm)
         do k = 0, nm
            chf1(k) = cfac * cbk(k)
            cfac = cfac * cf1
         end do
      end if
      return
   end subroutine ch12n

   subroutine ciknb(n, z, nm, cbi, cdi, cbk, cdk)
      use numconstants

      !*****************************************************************************80
      !
         !! CIKNB computes complex modified Bessel functions In(z) and Kn(z).
      !
      !  Discussion:
      !
      !    This procedure also evaluates the derivatives.
      !
      !  Licensing:
      !
      !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
      !    they give permission to incorporate this routine into a user program
      !    provided that the copyright is acknowledged.
      !
      !  Modified:
      !
      !    30 July 2012
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
      !    Input, integer ( kind = 4 ) N, the order of In(z) and Kn(z).
      !
      !    Input, complex ( kind = 8 ) Z, the argument.
      !
      !    Output, integer ( kind = 4 ) NM, the highest order computed.
      !
      !    Output, complex ( kind = 8 ) CB((0:N), CDI(0:N), CBK(0:N), CDK(0:N),
      !    the values of In(z), In'(z), Kn(z), Kn'(z).
      !
      implicit none

      integer(kind=4) n

      real(kind=8) a0
      complex(kind=8) ca0
      complex(kind=8) cbi(0:n + 1)
      complex(kind=8) cdi(0:n + 1)
      complex(kind=8) cbkl
      complex(kind=8) cbs
      complex(kind=8) cbk(0:n + 1)
      complex(kind=8) cdk(0:n + 1)
      complex(kind=8) cf
      complex(kind=8) cf0
      complex(kind=8) cf1
      complex(kind=8) cg
      complex(kind=8) cg0
      complex(kind=8) cg1
      complex(kind=8) ci
      complex(kind=8) cr
      complex(kind=8) cs0
      complex(kind=8) csk0
      real(kind=8) el
      real(kind=8) fac
      integer(kind=4) k
      integer(kind=4) k0
      integer(kind=4) l
      integer(kind=4) m
!           integer ( kind = 4 ) msta1
!           integer ( kind = 4 ) msta2
      integer(kind=4) nm
      real(kind=8) vt
      complex(kind=8) z
      complex(kind=8) z1
      el = 0.57721566490153D+00
      a0 = abs(z)
      nm = n
      if (a0 < 1.0D-100) then
         do k = 0, n
            cbi(k) = cmplx(0.0D+00, 0.0D+00, kind=8)
            cbk(k) = cmplx(1.0D+30, 0.0D+00, kind=8)
            cdi(k) = cmplx(0.0D+00, 0.0D+00, kind=8)
            cdk(k) = -cmplx(1.0D+30, 0.0D+00, kind=8)
         end do
         cbi(0) = cmplx(1.0D+00, 0.0D+00, kind=8)
         cdi(1) = cmplx(0.5D+00, 0.0D+00, kind=8)
         return
      end if
      ci = cmplx(0.0D+00, 1.0D+00, kind=8)
      if (real(z, kind=8) < 0.0D+00) then
         z1 = -z
      else
         z1 = z
      end if
      if (n == 0) then
         nm = 1
      end if
      m = msta1(a0, 200)
      if (m < nm) then
         nm = m
      else
         m = msta2(a0, nm, 15)
      end if
      cbs = 0.0D+00
      csk0 = 0.0D+00
      cf0 = 0.0D+00
      cf1 = 1.0D-100
      do k = m, 0, -1
         cf = 2.0D+00 * (k + 1.0D+00) * cf1 / z1 + cf0
         if (k <= nm) then
            cbi(k) = cf
         end if
         if (k /= 0 .and. k == 2 * int(k / 2)) then
            csk0 = csk0 + 4.0D+00 * cf / k
         end if
         cbs = cbs + 2.0D+00 * cf
         cf0 = cf1
         cf1 = cf
      end do
      cs0 = exp(z1) / (cbs - cf)
      do k = 0, nm
         cbi(k) = cs0 * cbi(k)
      end do
      if (a0 <= 9.0D+00) then
         cbk(0) = -(log(0.5D+00 * z1) + el) * cbi(0) + cs0 * csk0
         cbk(1) = (1.0D+00 / z1 - cbi(1) * cbk(0)) / cbi(0)
      else
         ca0 = sqrt(pi / (2.0D+00 * z1)) * exp(-z1)
         if (a0 < 25.0D+00) then
            k0 = 16
         else if (a0 < 80.0D+00) then
            k0 = 10
         else if (a0 < 200.0D+00) then
            k0 = 8
         else
            k0 = 6
         end if
         do l = 0, 1
            cbkl = 1.0D+00
            vt = 4.0D+00 * l
            cr = cmplx(1.0D+00, 0.0D+00, kind=8)
            do k = 1, k0
               cr = 0.125D+00 * cr &
                    * (vt - (2.0D+00 * k - 1.0D+00)**2) / (k * z1)
               cbkl = cbkl + cr
            end do
            cbk(l) = ca0 * cbkl
         end do
      end if
      cg0 = cbk(0)
      cg1 = cbk(1)
      do k = 2, nm
         cg = 2.0D+00 * (k - 1.0D+00) / z1 * cg1 + cg0
         cbk(k) = cg
         cg0 = cg1
         cg1 = cg
      end do
      if (real(z, kind=8) < 0.0D+00) then
         fac = 1.0D+00
         do k = 0, nm
            if (aimag(z) < 0.0D+00) then
               cbk(k) = fac * cbk(k) + ci * pi * cbi(k)
            else
               cbk(k) = fac * cbk(k) - ci * pi * cbi(k)
            end if
            cbi(k) = fac * cbi(k)
            fac = -fac
         end do
      end if
      cdi(0) = cbi(1)
      cdk(0) = -cbk(1)
      do k = 1, nm
         cdi(k) = cbi(k - 1) - k / z * cbi(k)
         cdk(k) = -cbk(k - 1) - k / z * cbk(k)
      end do
      return
   end subroutine ciknb

   subroutine bessel_integer_complex(n, z, nmax, b)
      use numconstants
      implicit none
      integer :: n, nmax
      complex(8) :: z, b(0:n), cbj(0:n + 1), cdj(0:n + 1), cby(0:n + 1), cdy(0:n + 1)

      call cjynb(n, z, nmax, cbj, cdj, cby, cdy)
      nmax = min(n, nmax)
      b(0:nmax) = cbj(0:nmax)
   end subroutine bessel_integer_complex

   subroutine cjynb(n, z, nm, cbj, cdj, cby, cdy)
      use numconstants

      !*****************************************************************************80
      !
         !! CJYNB: Bessel functions, derivatives, Jn(z) and Yn(z) of complex argument.
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
      integer(kind=4) n
      real(kind=8), save, dimension(4) :: a = (/ &
                                          -0.7031250000000000D-01, 0.1121520996093750D+00, &
                                          -0.5725014209747314D+00, 0.6074042001273483D+01/)
      real(kind=8) a0
      real(kind=8), save, dimension(4) :: a1 = (/ &
                                          0.1171875000000000D+00, -0.1441955566406250D+00, &
                                          0.6765925884246826D+00, -0.6883914268109947D+01/)
      real(kind=8), save, dimension(4) :: b = (/ &
                                          0.7324218750000000D-01, -0.2271080017089844D+00, &
                                          0.1727727502584457D+01, -0.2438052969955606D+02/)
      real(kind=8), save, dimension(4) :: b1 = (/ &
                                          -0.1025390625000000D+00, 0.2775764465332031D+00, &
                                          -0.1993531733751297D+01, 0.2724882731126854D+02/)
      complex(kind=8) cbj(0:n + 1)
      complex(kind=8) cbj0
      complex(kind=8) cbj1
      complex(kind=8) cbjk
      complex(kind=8) cbs
      complex(kind=8) cby(0:n + 1)
      complex(kind=8) cby0
      complex(kind=8) cby1
      complex(kind=8) cdj(0:n + 1)
      complex(kind=8) cdy(0:n + 1)
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
!           integer ( kind = 4 ) msta1
!           integer ( kind = 4 ) msta2
      integer(kind=4) nm
      real(kind=8) r2p
      real(kind=8) y0
      complex(kind=8) z

      el = 0.5772156649015329D+00
      r2p = 0.63661977236758D+00
      y0 = abs(aimag(z))
      a0 = abs(z)
      nm = n
      if (a0 < 1.0D-100) then
         do k = 0, n
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
         m = msta1(a0, 200)
         if (m < nm) then
            nm = m
         else
            m = msta2(a0, nm, 15)
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
         do k = 0, nm
            cbj(k) = cbj(k) / cs0
         end do
         ce = log(z / 2.0D+00) + el
         cby(0) = r2p * (ce * cbj(0) - 4.0D+00 * csu / cs0)
         cby(1) = r2p * (-cbj(0) / z + (ce - 1.0D+00) * cbj(1) &
                         - 4.0D+00 * csv / cs0)
      else
         ct1 = z - 0.25D+00 * pi
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
         ct2 = z - 0.75D+00 * pi
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
      do k = 1, nm
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
      do k = 1, nm
         cdy(k) = cby(k - 1) - k / z * cby(k)
      end do
      return
   end subroutine cjynb

   integer function msta1(x, mp)
      implicit none
      integer :: mp, n0, n1, it, nn
      real(8) :: x, a0, f1, f, f0
      a0 = dabs(x)
      n0 = int(1.1 * a0) + 1
      f0 = envj(n0, a0) - mp
      n1 = n0 + 5
      f1 = envj(n1, a0) - mp
      do it = 1, 20
         nn = n1 - (n1 - n0) / (1.0d0 - f0 / f1)
         f = envj(nn, a0) - mp
         if (abs(nn - n1) .lt. 1) exit
         n0 = n1
         f0 = f1
         n1 = nn
         f1 = f
      end do
      msta1 = nn
   end function msta1

   integer function msta2(x, n, mp)
      implicit none
      integer :: n, mp, n0, n1, it, nn
      real(8) :: x, a0, hmp, ejn, obj, f0, f1, f
      a0 = dabs(x)
      hmp = 0.5d0 * dble(mp)
      ejn = envj(n, a0)
      if (ejn .le. hmp) then
         obj = mp
         n0 = int(1.1 * a0)
      else
         obj = hmp + ejn
         n0 = n
      end if
      f0 = envj(n0, a0) - obj
      n1 = n0 + 5
      f1 = envj(n1, a0) - obj
      do it = 1, 20
         nn = n1 - (n1 - n0) / (1.0d0 - f0 / f1)
         f = envj(nn, a0) - obj
         if (abs(nn - n1) .lt. 1) exit
         n0 = n1
         f0 = f1
         n1 = nn
         f1 = f
      end do
      msta2 = nn + 10
   end function msta2

   real(8) function envj(n, x)
      implicit none
      integer :: n
      real(8) :: x
      n = max(1, abs(n))
      envj = 0.5d0 * dlog10(6.28d0 * n) - n * dlog10(1.36d0 * x / n)
   end function envj
end module bessel_functions
