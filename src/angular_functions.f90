module angular_functions
   use bessel_functions, only: riccati_bessel, riccati_hankel
   use coefficient_indexing, only: amnaddress, amnpaddress
   use constants
   implicit none
contains

   subroutine vcfunc(m, n, k, l, vcn)
      use numerical_tables
      implicit none
      integer :: m, n, k, l, wmax, wmin, w, mk
      real(8) :: vcn(0:n + l), t1, t2, t3, vcmax, vctest, rat
      vcn = 0.d0
      wmax = n + l
      wmin = max(abs(n - l), abs(m + k))
      vcn(wmax) = bcof(n + m, l + k) * bcof(n - m, l - k) / bcof(n + n, l + l)
      if (wmin .eq. wmax) return
      vcn(wmax - 1) = vcn(wmax) * (l * m - k * n) * fnr(2 * (l + n) - 1) / fnr(l) / fnr(n)&
     &  / fnr(n + l + m + k) / fnr(n + l - m - k)
      if (wmin .eq. wmax - 1) return
      mk = m + k
      vcmax = abs(vcn(wmax)) + abs(vcn(wmax - 1))
!
!  a downwards recurrence is used initially
!
      do w = wmax, wmin + 2, -1
         t1 = 2 * w * fnr(w + w + 1) * fnr(w + w - 1) / (fnr(w + mk) * fnr(w - mk)&
                                                  &     * fnr(n - l + w) * fnr(l - n + w) * fnr(n + l - w + 1) * fnr(n + l + w + 1))
         t2 = dble((m - k) * w * (w - 1) - mk * n * (n + 1) + mk * l * (l + 1))&
     &    / dble(2 * w * (w - 1))
         t3 = fnr(w - mk - 1) * fnr(w + mk - 1) * fnr(l - n + w - 1) * fnr(n - l + w - 1)&
             &     * fnr(n + l - w + 2) * fnr(n + l + w) / (dble(2 * (w - 1)) * fnr(2 * w - 3)&
                                                           &     * fnr(2 * w - 1))
         vcn(w - 2) = (t2 * vcn(w - 1) - vcn(w) / t1) / t3
         if (mod(wmax - w, 2) .eq. 1) then
            vctest = abs(vcn(w - 2)) + abs(vcn(w - 1))
            vcmax = max(vcmax, vctest)
            rat = vctest / vcmax
!
!  if/when the coefficients start to decrease in magnitude, an upwards recurrence takes over
!
            if (rat .lt. 0.01d0) exit
         end if
      end do
      if (w - 2 .gt. wmin) then
         wmax = w - 3
         call vcfuncuprec(m, n, k, l, wmax, vcn)
      end if
   end subroutine vcfunc
!
!  upwards VC coefficient recurrence
!
!
!  last revised: 15 January 2011
!
   subroutine vcfuncuprec(m, n, k, l, wmax, vcn)
      use numerical_tables
      implicit none
      integer :: m, n, k, l, wmax, w, mk, nl, m1, n1, l1, k1, w1, w2
      real(8) :: vcn(0:n + l), t1, t2, t3, vc1
      mk = abs(m + k)
      nl = abs(n - l)
      if (nl .ge. mk) then
         w = nl
         if (n .ge. l) then
            m1 = m
            n1 = n
            l1 = l
            k1 = k
         else
            m1 = k
            n1 = l
            k1 = m
            l1 = n
         end if
         vc1 = (-1)**(k1 + l1) * bcof(l1 + k1, w - m1 - k1) &
               * bcof(l1 - k1, w + m1 + k1) / bcof(l1 + l1, w + w + 1)
      else
         w = mk
         if (m + k .ge. 0) then
            vc1 = (-1)**(n + m) * bcof(n - l + w, l - k) * bcof(l - n + w, n - m) &
                  / bcof(w + w + 1, n + l - w)
         else
            vc1 = (-1)**(l + k) * bcof(n - l + w, l + k) * bcof(l - n + w, n + m) &
                  / bcof(w + w + 1, n + l - w)
         end if
      end if
      w1 = w
      vcn(w) = vc1
      w = w1 + 1
      mk = m + k
      w2 = min(wmax, n + l)
      if (w2 .gt. w1) then
         t1 = 2 * w * fnr(w + w + 1) * fnr(w + w - 1) / (fnr(w + mk) * fnr(w - mk) &
                                                        * fnr(n - l + w) * fnr(l - n + w) * fnr(n + l - w + 1) * fnr(n + l + w + 1))
         if (w1 .eq. 0) then
            t2 = .5 * dble(m - k)
         else
            t2 = dble((m - k) * w * (w - 1) - mk * n * (n + 1) + mk * l * (l + 1)) &
                 / dble(2 * w * (w - 1))
         end if
         vcn(w) = t1 * t2 * vcn(w1)
      end if
      do w = w1 + 2, w2
         t1 = 2 * w * fnr(w + w + 1) * fnr(w + w - 1) / (fnr(w + mk) * fnr(w - mk) &
                                                        * fnr(n - l + w) * fnr(l - n + w) * fnr(n + l - w + 1) * fnr(n + l + w + 1))
         t2 = dble((m - k) * w * (w - 1) - mk * n * (n + 1) + mk * l * (l + 1)) &
              / dble(2 * w * (w - 1))
         t3 = fnr(w - mk - 1) * fnr(w + mk - 1) * fnr(l - n + w - 1) * fnr(n - l + w - 1) &
              * fnr(n + l - w + 2) * fnr(n + l + w) / (dble(2 * (w - 1)) * fnr(2 * w - 3) &
                                                       * fnr(2 * w - 1))
         vcn(w) = t1 * (t2 * vcn(w - 1) - t3 * vcn(w - 2))
      end do
   end subroutine vcfuncuprec
!
!  Normalized associated legendre functions
!
!
!  last revised: 15 January 2011
!
   subroutine normalized_associated_legendre(cbe, mmax, nmax, dc)
      use numerical_tables
      implicit none
      integer :: nmax, mmax, m, n, im
      real(8) :: dc(-mmax:mmax, 0:nmax), cbe, sbe
      sbe = dsqrt((1.d0 + cbe) * (1.d0 - cbe))
      dc = 0.d0
      do m = 0, mmax
         dc(m, m) = (-1)**m * (0.5d0 * sbe)**m * bcof(m, m)
         if (m .eq. nmax) exit
         dc(m, m + 1) = fnr(m + m + 1) * cbe * dc(m, m)
         do n = m + 1, nmax - 1
            dc(m, n + 1) = (-fnr(n - m) * fnr(n + m) * dc(m, n - 1) + dble(n + n + 1) * cbe * dc(m, n)) &
                           / (fnr(n + 1 - m) * fnr(n + 1 + m))
         end do
      end do
      do m = 1, mmax
         im = (-1)**m
         do n = m, nmax
            dc(-m, n) = im * dc(m, n)
         end do
      end do
   end subroutine normalized_associated_legendre

!!
!  Generalized spherical functions
!
!  dc(m,n*(n+1)+k)=(-1)^(m + k)((n - k)!(n + k)!/(n - m)!/(n + m)!)^(1/2)
!  ((1 + x)/2)^((m + k)/2)((1 - x)/2)^((k - m)/2)JacobiP[n - k, k - m, k + m, x]
!
!  for |m| <= kmax, n=0,1,...nmax, |k| <= n
!
!
!  last revised: 15 January 2011
!
   subroutine rotation_coefficients(cbe, kmax, nmax, dc)
      use numerical_tables
      implicit none
      integer :: kmax, nmax, k, m, sin, n, knmax, nn1, kn, im, m1
      real(8) :: cbe, sbe, dc(-kmax:kmax, 0:nmax * (nmax + 2)), cbe2, sbe2, dk0(-nmax - 1:nmax + 1), &
                 dk01(-nmax - 1:nmax + 1), sben, dkt, fmn, dkm0, dkm1, dkn1
!if(light_up) then
!write(*,'('' rot1 '',3es13.5)') cbe
!flush(6)
!endif
      dc = 0.d0
      if (abs(cbe) .ge. 1.d0) then
         sbe = 0.d0
      else
         sbe = dsqrt(abs((1.d0 + cbe) * (1.d0 - cbe)))
      end if
      cbe2 = .5d0 * (1.d0 + cbe)
      sbe2 = .5d0 * (1.d0 - cbe)
      sin = 1
      dk0(0) = 1.d0
      sben = 1.d0
      dc(0, 0) = 1.d0
      dk01(0) = 0.d0
!if(light_up) then
!write(*,'('' rot2 '',3es13.5)') cbe,sbe
!flush(6)
!endif
      do n = 1, nmax
         knmax = min(n, kmax)
         nn1 = n * (n + 1)
         sin = -sin
         sben = sben * sbe / 2.d0
         if (sben .lt. 1.d-30) sben = 0.d0
         dk0(n) = sin * sben * bcof(n, n)
         dk0(-n) = sin * dk0(n)
         dk01(n) = 0.d0
         dk01(-n) = 0.d0
         dc(0, nn1 + n) = dk0(n)
         dc(0, nn1 - n) = dk0(-n)
!if(light_up.and.abs(cbe).gt.0.99999d0) then
!write(*,'('' rot2b '',i5,3es13.5)') n,sben
!flush(6)
!endif
         do k = -n + 1, n - 1
            kn = nn1 + k
            dkt = dk01(k)
            dk01(k) = dk0(k)
            dk0(k) = (cbe * dble(n + n - 1) * dk01(k) - fnr(n - k - 1) * fnr(n + k - 1) * dkt) &
                     / (fnr(n + k) * fnr(n - k))
            dc(0, kn) = dk0(k)
         end do
         im = 1
         do m = 1, knmax
            im = -im
            fmn = 1.d0 / fnr(n - m + 1) / fnr(n + m)
            m1 = m - 1
            dkm0 = 0.d0
            do k = -n, n
               kn = nn1 + k
               dkm1 = dkm0
               dkm0 = dc(m1, kn)
               if (k .eq. n) then
                  dkn1 = 0.d0
               else
                  dkn1 = dc(m1, kn + 1)
               end if
               dc(m, kn) = (fnr(n + k) * fnr(n - k + 1) * cbe2 * dkm1 &
                            - fnr(n - k) * fnr(n + k + 1) * sbe2 * dkn1 &
                            - dble(k) * sbe * dc(m1, kn)) * fmn
               dc(-m, nn1 - k) = dc(m, kn) * (-1)**(k) * im
            end do
         end do
      end do
!if(light_up) then
!write(*,'('' rot3 '',3es13.5)') cbe
!flush(6)
!endif
   end subroutine rotation_coefficients
!
!  Generalized spherical functions: complex argument
!
!  dc(m,n*(n+1)+k)=(-1)^(m + k)((n - k)!(n + k)!/(n - m)!/(n + m)!)^(1/2)
!  ((1 + x)/2)^((m + k)/2)((1 - x)/2)^((k - m)/2)JacobiP[n - k, k - m, k + m, x]
!
!  for |m| <= kmax, n=0,1,...nmax, |k| <= n
!
!
!  New: 08/25/2011
!
   subroutine complex_rotation_coefficients(cbe, kmax, nmax, dc, sin_beta)
      use numerical_tables
      implicit none
      integer :: kmax, nmax, k, m, in, n, knmax, nn1, kn, im, m1
      real(8) :: fmn
      complex(8) :: cbe, sbe, dc(-kmax:kmax, 0:nmax * (nmax + 2)), cbe2, sbe2, dk0(-nmax - 1:nmax + 1), &
                    dk01(-nmax - 1:nmax + 1), sben, dkt, dkm0, dkm1, dkn1
      complex(8), optional :: sin_beta
      if (present(sin_beta)) then
         sbe = sin_beta
      else
         sbe = sqrt((1.d0 + cbe) * (1.d0 - cbe))
      end if
      cbe2 = .5d0 * (1.d0 + cbe)
      sbe2 = .5d0 * (1.d0 - cbe)
      in = 1
      dk0(0) = 1.d0
      sben = 1.d0
      dc(0, 0) = 1.d0
      dk01(0) = 0.
      do n = 1, nmax
         knmax = min(n, kmax)
         nn1 = n * (n + 1)
         in = -in
         sben = sben * sbe / 2.d0
         dk0(n) = in * sben * bcof(n, n)
         dk0(-n) = in * dk0(n)
         dk01(n) = 0.
         dk01(-n) = 0.
         dc(0, nn1 + n) = dk0(n)
         dc(0, nn1 - n) = dk0(-n)
         do k = -n + 1, n - 1
            kn = nn1 + k
            dkt = dk01(k)
            dk01(k) = dk0(k)
            dk0(k) = (cbe * dble(n + n - 1) * dk01(k) - fnr(n - k - 1) * fnr(n + k - 1) * dkt) &
                     / (fnr(n + k) * fnr(n - k))
            dc(0, kn) = dk0(k)
         end do
         im = 1
         do m = 1, knmax
            im = -im
            fmn = 1.d0 / fnr(n - m + 1) / fnr(n + m)
            m1 = m - 1
            dkm0 = 0.
            do k = -n, n
               kn = nn1 + k
               dkm1 = dkm0
               dkm0 = dc(m1, kn)
               if (k .eq. n) then
                  dkn1 = 0.
               else
                  dkn1 = dc(m1, kn + 1)
               end if
               dc(m, kn) = (fnr(n + k) * fnr(n - k + 1) * cbe2 * dkm1 &
                            - fnr(n - k) * fnr(n + k + 1) * sbe2 * dkn1 &
                            - dble(k) * sbe * dc(m1, kn)) * fmn
               dc(-m, nn1 - k) = dc(m, kn) * (-1)**(k) * im
            end do
         end do
      end do
   end subroutine complex_rotation_coefficients
!
! vector spherical harmonic function
! november 2011
! april 2012: lr formulation
! 2020 : back to nm formulation.   Complex cb
!
   subroutine complexpivec(cb, nodr, pivec, icon, lr_model, azimuth_angle, index_model)
      use numerical_tables
      implicit none
      logical :: lrmod
      logical, optional :: lr_model
      integer :: nodr, n, m, p, mn, i, nn1, imod, mnp, q
      integer, optional :: icon, index_model
      real(8) :: fnm, const, alpha
      real(8), optional :: azimuth_angle
      complex(8) :: cb, pivec(2 * nodr * (nodr + 2), 2), drot(-1:1, 0:nodr * (nodr + 2)), tau(2), ci, cin, ephi
      data ci/(0.d0, 1.d0)/

      if (present(index_model)) then
         imod = index_model
      else
         imod = 1
      end if
      if (present(azimuth_angle)) then
         alpha = azimuth_angle
      else
         alpha = 0.d0
      end if
      if (present(icon)) then
         i = icon
      else
         i = 1
      end if
      if (present(lr_model)) then
         lrmod = lr_model
      else
         lrmod = .false.
      end if
      const = 0.5d0 / sqrt_two_pi
      call complex_rotation_coefficients(cb, 1, nodr, drot)
      do n = 1, nodr
         nn1 = n * (n + 1)
         fnm = sqrt(dble(n + n + 1) / 2.d0) / 4.d0 * const
         cin = 4.d0 * (-i * ci)**(n + 1)
         do m = -n, n
            ephi = exp(i * ci * dble(m) * alpha)
            mn = nn1 + m
            tau(1) = -fnm * (-drot(-1, mn) + drot(1, mn)) * ephi
            tau(2) = -fnm * (drot(-1, mn) + drot(1, mn)) * ephi
            do p = 1, 2
               mnp = amnpaddress(m, n, p, nodr, imod)
               pivec(mnp, 1) = cin * tau(p)
               pivec(mnp, 2) = i * ci * cin * tau(3 - p)
            end do
            if (lrmod) then
               do q = 1, 2
                  do p = 1, 2
                     mnp = amnpaddress(m, n, p, nodr, imod)
                     tau(p) = pivec(mnp, q)
                  end do
                  mnp = amnpaddress(m, n, 1, nodr, imod)
                  pivec(mnp, q) = tau(1) + tau(2)
                  mnp = amnpaddress(m, n, 2, nodr, imod)
                  pivec(mnp, q) = tau(1) - tau(2)
               end do
            end if
         end do
      end do
   end subroutine complexpivec
!
!  tau are the vector spherical harmonic functions, normalized
!
!
!  last revised: 15 January 2011
!
   subroutine taufunc(cb, nmax, tau)
      use numerical_tables
      implicit none
      integer :: nmax, n, m, nn1, mn
      real(8) :: drot(-1:1, 0:nmax * (nmax + 2)), tau(0:nmax + 1, nmax, 2), cb, fnm
      call rotation_coefficients(cb, 1, nmax, drot)
      do n = 1, nmax
         nn1 = n * (n + 1)
         fnm = sqrt(dble(n + n + 1) / 2.d0) / 4.d0
         do m = -n, -1
            mn = nn1 + m
            tau(n + 1, -m, 1) = -fnm * (-drot(-1, mn) + drot(1, mn))
            tau(n + 1, -m, 2) = -fnm * (drot(-1, mn) + drot(1, mn))
         end do
         do m = 0, n
            mn = nn1 + m
            tau(m, n, 1) = -fnm * (-drot(-1, mn) + drot(1, mn))
            tau(m, n, 2) = -fnm * (drot(-1, mn) + drot(1, mn))
         end do
      end do
   end subroutine taufunc
!
!  rotation of expansion coefficients amn by euler angles alpha,beta,gamma
!  idir=1: forward rotation, idir=-1, reverse rotation.
!
!
!  last revised: 15 January 2011
!
   subroutine rotate_expansion_coefficients(alpha, beta, gamma, nmax, mmax, amn, idir)
      use numerical_tables
      implicit none
      integer :: nmax, mmax, idir, k, n, m, in, kmax, ka, na, im, m1
      real(8) :: dc(-nmax - 1:nmax + 1, -nmax - 1:nmax + 1), dk0(-nmax - 1:nmax + 1), &
                 dk01(-nmax - 1:nmax + 1), sbe, cbe, sbe2, cbe2, sben, dkt, &
                 fmn, dkm0, dkm1, alpha, beta, gamma
      complex(8) :: ealpha, amn(0:nmax + 1, nmax, 2), ealpham(-nmax:nmax), &
                    amnt(2, -nmax:nmax), a, b, ci, egamma, egammam(-nmax:nmax)
      data ci/(0.d0, 1.d0)/
      call init(nmax)
      dc = 0.d0
      dk01 = 0.d0
      dk0 = 0.d0
      ealpha = exp(ci * alpha)
      egamma = exp(ci * gamma)
      cbe = cos(beta)
      sbe = sqrt((1.d0 + cbe) * (1.d0 - cbe))
      cbe2 = .5d0 * (1.d0 + cbe)
      sbe2 = .5d0 * (1.d0 - cbe)
      call azimuthal_phase_factors(ealpha, nmax, ealpham)
      call azimuthal_phase_factors(egamma, nmax, egammam)
      in = 1
      dk0(0) = 1.d0
      sben = 1.d0
      dk01(0) = 0.d0
      do n = 1, nmax
         kmax = min(n, mmax)
         do k = -kmax, kmax
            if (k .le. -1) then
               ka = n + 1
               na = -k
            else
               ka = k
               na = n
            end if
            if (idir .eq. 1) then
               amnt(1, k) = amn(ka, na, 1) * ealpham(k)
               amnt(2, k) = amn(ka, na, 2) * ealpham(k)
            else
               amnt(1, -k) = amn(ka, na, 1) * egammam(k)
               amnt(2, -k) = amn(ka, na, 2) * egammam(k)
            end if
         end do
         in = -in
         sben = sben * sbe / 2.d0
         dk0(n) = in * sben * bcof(n, n)
         dk0(-n) = in * dk0(n)
         dk01(n) = 0.d0
         dk01(-n) = 0.d0
         dc(0, n) = dk0(n)
         dc(0, -n) = dk0(-n)
         do k = -n + 1, n - 1
            dkt = dk01(k)
            dk01(k) = dk0(k)
            dk0(k) = (cbe * (n + n - 1) * dk01(k) - fnr(n - k - 1) * fnr(n + k - 1) * dkt) &
                     / (fnr(n + k) * fnr(n - k))
            dc(0, k) = dk0(k)
         end do
         im = 1
         do m = 1, kmax
            im = -im
            fmn = 1./fnr(n - m + 1) / fnr(n + m)
            m1 = m - 1
            dkm0 = 0.
            do k = -n, n
               dkm1 = dkm0
               dkm0 = dc(m1, k)
               dc(m, k) = (fnr(n + k) * fnr(n - k + 1) * cbe2 * dkm1 &
                           - fnr(n - k) * fnr(n + k + 1) * sbe2 * dc(m1, k + 1) &
                           - k * sbe * dc(m1, k)) * fmn
               dc(-m, -k) = dc(m, k) * (-1)**(k) * im
            end do
         end do
         do m = -n, n
            if (m .le. -1) then
               ka = n + 1
               na = -m
            else
               ka = m
               na = n
            end if
            a = 0.
            b = 0.
            do k = -kmax, kmax
               a = a + dc(-k, -m) * amnt(1, k)
               b = b + dc(-k, -m) * amnt(2, k)
            end do
            if (idir .eq. 1) then
               amn(ka, na, 1) = a * egammam(m)
               amn(ka, na, 2) = b * egammam(m)
            else
               amn(ka, na, 1) = a * ealpham(m)
               amn(ka, na, 2) = b * ealpham(m)
            end if
         end do
      end do
   end subroutine rotate_expansion_coefficients
!
!  regular vswf expansion coefficients for a plane wave: general case, complex cos beta
!
   subroutine genplanewavecoef(alpha, cb, nodr, pmnp0, lr_tran)
      use numerical_tables
      implicit none
      logical :: lrtran
      logical, optional :: lr_tran
      integer :: nodr, m, n, p, sp, nn1, mn
      real(8) :: alpha, fnm, ca, sa
      complex(8) :: drot(-1:1, 0:nodr * (nodr + 2)), tau(0:nodr + 1, nodr, 2), &
                    taulr(0:nodr + 1, nodr, 2), cb, ealpha, ci, cin, &
                    pmnp0(0:nodr + 1, nodr, 2, 2), ealpham(-nodr:nodr)
      data ci/(0.d0, 1.d0)/
      if (present(lr_tran)) then
         lrtran = lr_tran
      else
         lrtran = .true.
      end if
      call complex_rotation_coefficients(cb, 1, nodr, drot)
      do n = 1, nodr
         nn1 = n * (n + 1)
         fnm = sqrt(dble(n + n + 1) / 2.d0) / 4.d0
         do m = -n, -1
            mn = nn1 + m
            tau(n + 1, -m, 1) = -fnm * (-drot(-1, mn) + drot(1, mn))
            tau(n + 1, -m, 2) = -fnm * (drot(-1, mn) + drot(1, mn))
         end do
         do m = 0, n
            mn = nn1 + m
            tau(m, n, 1) = -fnm * (-drot(-1, mn) + drot(1, mn))
            tau(m, n, 2) = -fnm * (drot(-1, mn) + drot(1, mn))
         end do
      end do
      ca = cos(alpha)
      sa = sin(alpha)
      ealpha = cmplx(ca, sa, kind=kind(0.0d0))
      call azimuthal_phase_factors(ealpha, nodr, ealpham)
      if (lrtran) then
         taulr(:, :, 1) = (tau(:, :, 1) + tau(:, :, 2))*.5d0
         taulr(:, :, 2) = (tau(:, :, 1) - tau(:, :, 2))*.5d0
         do n = 1, nodr
            cin = 4.d0 * ci**(n + 1)
            do p = 1, 2
               sp = -(-1)**p
               do m = -n, -1
                  pmnp0(n + 1, -m, p, 1) = -cin * taulr(n + 1, -m, p) * ealpham(-m)
                  pmnp0(n + 1, -m, p, 2) = sp * ci * cin * taulr(n + 1, -m, p) * ealpham(-m)
               end do
               do m = 0, n
                  pmnp0(m, n, p, 1) = -cin * taulr(m, n, p) * ealpham(-m)
                  pmnp0(m, n, p, 2) = sp * ci * cin * taulr(m, n, p) * ealpham(-m)
               end do
            end do
         end do
      else
         do n = 1, nodr
            cin = 4.d0 * ci**(n + 1)
            do p = 1, 2
               do m = -n, -1
                  pmnp0(n + 1, -m, p, 1) = -cin * tau(n + 1, -m, p) * ealpham(-m)
                  pmnp0(n + 1, -m, p, 2) = ci * cin * tau(n + 1, -m, 3 - p) * ealpham(-m)
               end do
               do m = 0, n
                  pmnp0(m, n, p, 1) = -cin * tau(m, n, p) * ealpham(-m)
                  pmnp0(m, n, p, 2) = ci * cin * tau(m, n, 3 - p) * ealpham(-m)
               end do
            end do
         end do
      end if
   end subroutine genplanewavecoef

   subroutine gaussianbeamcoef(alpha, cbeta, cbeam, nodr, pmnp0, lr_tran)
      use numerical_tables
      implicit none
      logical :: lrtran
      logical, optional :: lr_tran
      integer :: nodr, m, n, p, k
      real(8) :: alpha, cbeta, cbeam, gbn
      complex(8) :: ccb, pmnp0(0:nodr + 1, nodr, 2, 2)
      if (present(lr_tran)) then
         lrtran = lr_tran
      else
         lrtran = .true.
      end if
      ccb = cbeta
      call genplanewavecoef(alpha, ccb, nodr, pmnp0, lr_tran=lrtran)
      do n = 1, nodr
         gbn = dexp(-((dble(n) + .5d0) * cbeam)**2.)
         do p = 1, 2
            do k = 1, 2
               do m = -n, -1
                  pmnp0(n + 1, -m, p, k) = pmnp0(n + 1, -m, p, k) * gbn
               end do
               do m = 0, n
                  pmnp0(m, n, p, k) = pmnp0(m, n, p, k) * gbn
               end do
            end do
         end do
      end do
   end subroutine gaussianbeamcoef
!
!  axial translation coefficients calculated by the diamond recurrence formula
!  new: 10 october 2011
!  april 2012: lr formulation
!  may 2012: new ordering scheme:
!  input:  itype : 1 or 3 (regular, outgoing)
!     r: axial translation distance (positive)
!     ri: rank 2 complex array: L and R refractive indices of medium
!     nmax, lmax: largest row and column orders.
!     ndim: dimension of ac
!  output:
!     ac:  rank 1 complex array, dimension ndim, containing the matrix elements.
!  storage scheme:   for each degree m, with ordering m=0, -1, 1, -2, 2, ..min(nmax,lmax),
!  the elements for degree m are stored
!
   subroutine axialtrancoefrecurrence(itype, r, ri, nmax, lmax, ndim, ac)
      use numerical_tables
      implicit none
      integer :: itype, nmax, lmax, n, l, m, p, nlmin, &
                 wmin, wmax, ml, m1, np1, nm1, lm1, lp1, sp
      integer :: iadd, nlmax, iadd0, iadd1, ndim
      integer :: ma, blockdim
      integer, save :: nlmax0
      real(8) :: r
      complex(8) :: ri(2), ci, z(2), xi(0:nmax + lmax, 2)
      complex(8) :: ac(ndim), act(nmax, lmax, 2), actt(2, 2)
      data ci, nlmax0/(0.d0, 1.d0), 0/
      nlmax = max(nmax, lmax)
      nlmin = min(nmax, lmax)
      if (nlmax .gt. nlmax0) then
         nlmax0 = nlmax
         call axialtrancoefinit(nlmax)
      end if

      if (r .le. 1.d-12) then
         ac = (0.d0, 0.d0)
         if (itype .ne. 1) return
         iadd0 = 0
         do ma = 0, nlmin
            m1 = max(1, ma)
            do m = -ma, ma, 2 * m1
               blockdim = (nmax - m1 + 1) * (lmax - m1 + 1) * 2
               iadd1 = iadd0 + blockdim
               act = 0.d0
               do l = m1, nlmin
                  act(l, l, 1) = 1.d0
                  act(l, l, 2) = 1.d0
               end do
               ac(iadd0 + 1:iadd1) = reshape(act(m1:nmax, m1:lmax, 1:2), (/blockdim/))
               iadd0 = iadd1
            end do
         end do
         return
      end if
      z = r * ri
      do p = 1, 2
         if (itype .eq. 1) then
            call riccati_bessel(nmax + lmax, z(p), xi(0:, p))
         else
            call riccati_hankel(nmax + lmax, z(p), xi(0:, p))
         end if
         xi(0:, p) = xi(0:, p) / z(p)
         if (z(1) .eq. z(2)) then
            xi(0:nmax + lmax, 2) = xi(0:nmax + lmax, 1)
            exit
         end if
      end do
      lm1 = lmax - 1

      iadd0 = 0
      do ma = 0, nlmin
         m1 = max(1, ma)
         lp1 = m1 + 1
         do m = -ma, ma, 2 * m1
            blockdim = 2 * (nmax - m1 + 1) * (lmax - m1 + 1)
            iadd1 = iadd0 + blockdim
            n = m1
            do l = m1, lmax
               wmin = abs(n - l)
               wmax = n + l
               iadd = iadd + 1
               ml = l * (l + 1) + m
               do p = 1, 2
                  actt(1, p) = sum(vcc_const(n, ml, wmin:wmax:2) * xi(wmin:wmax:2, p))
                  actt(2, p) = ci * sum(vcc_const(n, ml, wmin + 1:wmax - 1:2) * xi(wmin + 1:wmax - 1:2, p))
               end do
               act(n, l, 1) = actt(1, 1) + actt(2, 1)
               act(n, l, 2) = actt(1, 2) - actt(2, 2)
            end do
            l = lmax
            ml = l * (l + 1) + m
            do n = m1 + 1, nmax
               wmin = abs(n - l)
               wmax = n + l
               do p = 1, 2
                  actt(1, p) = sum(vcc_const(n, ml, wmin:wmax:2) * xi(wmin:wmax:2, p))
                  actt(2, p) = ci * sum(vcc_const(n, ml, wmin + 1:wmax - 1:2) * xi(wmin + 1:wmax - 1:2, p))
               end do
               act(n, l, 1) = actt(1, 1) + actt(2, 1)
               act(n, l, 2) = actt(1, 2) - actt(2, 2)
            end do

            if (m1 .lt. nlmin) then
               do n = m1, nmax - 1
                  np1 = n + 1
                  nm1 = n - 1
                  do p = 1, 2
                     sp = -(-1)**p
                     act(np1, m1:lmax - 1, p) = &
                        -act(n, m1 + 1:lmax, p) * fnp1_const(m, m1:lm1) &
                        + sp * (fn_const(m, m1:lm1) - fn_const(m, n)) * ci * act(n, m1:lm1, p)
                     act(np1, m1 + 1:lm1, p) = act(np1, m1 + 1:lm1, p) &
                                               + act(n, m1:lmax - 2, p) * fnm1_const(m, lp1:lm1)
                     if (n .gt. m1) then
                        act(np1, m1:lm1, p) = act(np1, m1:lm1, p) &
                                              + act(nm1, m1:lm1, p) * fnm1_const(m, n)
                     end if
                     act(np1, m1:lm1, p) = act(np1, m1:lm1, p) / fnp1_const(m, n)
                  end do
               end do
            end if
            ac(iadd0 + 1:iadd1) = reshape(act(m1:nmax, m1:lmax, 1:2), (/blockdim/))
            iadd0 = iadd1
         end do
      end do
   end subroutine axialtrancoefrecurrence
!
!  constants for translation coefficient calculation
!
   subroutine axialtrancoefinit(nmax)
      use numerical_tables
      implicit none
      integer :: nmax, m, n, l, w, n21, ml, ll1, wmin, wmax, nlmin, lp1, lm1
      real(8) :: c1, c2, vc1(0:2 * nmax), vc2(0:2 * nmax)
      complex(8) :: ci, inlw
      data ci/(0.d0, 1.d0)/
      if (allocated(vcc_const)) deallocate (vcc_const, fnm1_const, fn_const, fnp1_const)
      allocate (vcc_const(nmax, nmax * (nmax + 2), 0:2 * nmax), fnm1_const(-nmax:nmax, nmax), &
                fn_const(-nmax:nmax, nmax), fnp1_const(-nmax:nmax, nmax))
      do n = 1, nmax
         n21 = n + n + 1
         do l = 1, nmax
            c1 = fnr(n21) * fnr(l + l + 1)
            ll1 = l * (l + 1)
            call vcfunc(-1, n, 1, l, vc2)
            wmin = abs(n - l)
            wmax = n + l
            nlmin = min(l, n)
            do m = -nlmin, nlmin
               ml = ll1 + m
               c2 = -c1 * (-1)**m
               call vcfunc(-m, n, m, l, vc1)
               do w = wmin, wmax
                  inlw = ci**(n - l + w)
                  vcc_const(n, ml, w) = c2 * vc1(w) * vc2(w) * (dble(inlw) + aimag(inlw))
               end do
            end do
         end do
      end do
      fnm1_const = 0.
      fn_const = 0.
      fnp1_const = 0.
      do m = -nmax, nmax
         do l = max(1, abs(m)), nmax
            lp1 = l + 1
            lm1 = l - 1
            fnm1_const(m, l) = fnr(lm1) * fnr(lp1) * fnr(l - m) * fnr(l + m) / fnr(lm1 + l) / fnr(l + lp1) / dble(l)
            fn_const(m, l) = dble(m) / dble(l) / dble(lp1)
            fnp1_const(m, l) = fnr(l) * fnr(l + 2) * fnr(lp1 - m) * fnr(lp1 + m) / fnr(l + lp1) / fnr(l + l + 3) / dble(lp1)
         end do
      end do
   end subroutine axialtrancoefinit

   subroutine gentrancoefconstants(nodrmax)
      use numerical_tables
      implicit none
      integer :: nodrmax, v, w, wmax, wmin, n, l, m, k, m1m, mn, kl
      real(8) :: vc1(0:2 * nodrmax), vc2(0:2 * nodrmax)
      complex(8) :: ci, c, a
      data ci/(0.d0, 1.d0)/
      if (allocated(tran_coef)) deallocate (tran_coef)
      allocate (tran_coef(nodrmax * (nodrmax + 2), nodrmax * (nodrmax + 2), 0:2 * nodrmax))
      tran_coef = 0.d0
      do l = 1, nodrmax
         do n = 1, nodrmax
            wmax = n + l
            call vcfunc(-1, n, 1, l, vc2)
            c = -ci**(n - l) * fnr(n + n + 1) * fnr(l + l + 1)
            do k = -l, l
               kl = l * (l + 1) + k
               do m = -n, n
                  mn = n * (n + 1) + m
                  m1m = (-1)**m
                  v = k - m
                  call vcfunc(-m, n, k, l, vc1)
                  wmin = max(abs(v), abs(n - l))
                  do w = wmax, wmin, -1
                     a = ci**w * c * m1m * vc1(w) * vc2(w)
                     if (mod(wmax - w, 2) .eq. 0) then
                        tran_coef(mn, kl, w) = dble(a)
                     else
                        tran_coef(mn, kl, w) = aimag(a)
                     end if
                  end do
               end do
            end do
         end do
      end do
      return
   end subroutine gentrancoefconstants

   subroutine gentranmatrix(nodr_s, nodr_t, translation_vector, &
                            refractive_index, ac_matrix, vswf_type, &
                            mode_s, mode_t, index_model)
      use numerical_tables
      implicit none
      integer :: nodr_s, nodr_t, nodrmax, wmax, p, n, m, k, l, mn, kl, imodel, &
                 nblks, nblkt, w, v, wmin, itype, nmodes, nmodet, nmode, mna, kla
      integer, optional :: vswf_type, mode_s, mode_t, index_model
      integer, save :: setnodrmax
      real(8) :: r, ct, xp(3), ymn(-nodr_s - nodr_t:nodr_s + nodr_t, 0:nodr_s + nodr_t)
      real(8) :: translation_vector(3)
      complex(8) :: ri(2), ephi, rri, ephim(-nodr_s - nodr_t:nodr_s + nodr_t), &
                    hn(0:nodr_s + nodr_t, 2), a1, a2, b1, b2
      complex(8) :: ac_matrix(nodr_t * (nodr_t + 2), nodr_s * (nodr_s + 2), 1:2)
      complex(8), optional :: refractive_index(2)
      data setnodrmax/0/
      nblks = nodr_s * (nodr_s + 2)
      nblkt = nodr_t * (nodr_t + 2)
      if (present(mode_s)) then
         nmodes = mode_s
      else
         nmodes = 2
      end if
      if (present(mode_t)) then
         nmodet = mode_t
      else
         nmodet = 2
      end if
      if (present(index_model)) then
         imodel = index_model
      else
         imodel = 2
      end if
      nmode = max(nmodet, nmodes)
      xp = translation_vector
      if (present(refractive_index)) then
         ri = refractive_index
      else
         ri = (1.d0, 0.d0)
      end if
      if (present(vswf_type)) then
         itype = vswf_type
      else
         itype = 3
      end if
      nodrmax = max(nodr_s, nodr_t)
      if (nodrmax .gt. setnodrmax) then
         setnodrmax = nodrmax
         call gentrancoefconstants(nodrmax)
      end if
      wmax = nodr_s + nodr_t
      r = xp(1) * xp(1) + xp(2) * xp(2) + xp(3) * xp(3)
      if (r .eq. 0.d0) then
         ac_matrix = 0.d0
         if (itype .eq. 1) then
            do n = 1, min(nblks, nblkt)
               ac_matrix(n, n, 1:nmode) = 1.d0
            end do
         end if
      else
         r = sqrt(r)
         ct = xp(3) / r
         if (xp(1) .eq. 0.d0 .and. xp(2) .eq. 0.d0) then
            ephi = (1.d0, 0.d0)
         else
            ephi = cmplx(xp(1), xp(2), kind=kind(0.0d0)) / sqrt(xp(1) * xp(1) + xp(2) * xp(2))
         end if
         ephim(0) = 1.d0
         do m = 1, wmax
            ephim(m) = ephi * ephim(m - 1)
            ephim(-m) = conjg(ephim(m))
         end do
         call normalized_associated_legendre(ct, wmax, wmax, ymn)
         do p = 1, 2
            rri = r * ri(p)
            if (itype .eq. 3) then
               hn(0, p) = -(0.d0, 1.d0) * exp((0.d0, 1.d0) * rri) / rri
               hn(1, p) = -exp((0.d0, 1.d0) * rri) * ((0.d0, 1.d0) + rri) / rri / rri
               do n = 2, wmax
                  hn(n, p) = dble(n + n - 1) / rri * hn(n - 1, p) - hn(n - 2, p)
               end do
            else
               call riccati_bessel(wmax, rri, hn(:, p))
               hn(:, p) = hn(:, p) / rri
            end if
            if (ri(2) .eq. ri(1)) then
               hn(:, 2) = hn(:, 1)
               exit
            end if
         end do
         do n = 1, nodr_s
            do m = -n, n
               mna = amnaddress(m, n, nodr_s, imodel)
               mn = n * (n + 1) + m
               do l = 1, nodr_t
                  wmax = n + l
                  do k = -l, l
                     kla = amnaddress(k, l, nodr_t, imodel)
                     kl = l * (l + 1) + k
                     v = m - k
                     wmin = max(abs(v), abs(n - l))
                     a1 = 0.
                     a2 = 0.
                     b1 = 0.
                     b2 = 0.
                     do w = wmax, wmin, -1
                        if (mod(wmax - w, 2) .eq. 0) then
                           a1 = a1 + hn(w, 1) * ymn(v, w) * tran_coef(kl, mn, w)
                           if (nmode .eq. 1) cycle
                           a2 = a2 + hn(w, 2) * ymn(v, w) * tran_coef(kl, mn, w)
                        else
                           b1 = b1 + hn(w, 1) * ymn(v, w) * tran_coef(kl, mn, w)
                           b2 = b2 + hn(w, 2) * ymn(v, w) * tran_coef(kl, mn, w)
                        end if
                     end do
                     if (nmode .eq. 1) then
                        ac_matrix(kla, mna, 1) = ephim(v) * a1
                     else
                        ac_matrix(kla, mna, 1) = ephim(v) * (a1 + (0.d0, 1.d0) * b1)
                        ac_matrix(kla, mna, 2) = ephim(v) * (a2 - (0.d0, 1.d0) * b2)
                     end if
                  end do
               end do
            end do
         end do
      end if
   end subroutine gentranmatrix
!
!  test to determine convergence of regular vswf addition theorem for max. order lmax
!  and translation distance r w/ refractive index ri.
!
!
!  last revised: 15 January 2011
!
   subroutine tranordertest(r, ri, lmax, eps, nmax)
      use numerical_tables
      implicit none
      integer :: nmax, lmax, n, l, m, w, n21, wmin, wmax
      integer, parameter :: nlim = 200
      real(8) :: r, alnw, sum, eps
      real(8) :: vc1(0:nlim + lmax)
      complex(8) :: ri, ci, z, a, b, c
      complex(8) :: xi(0:nlim + lmax)
      data ci/(0.d0, 1.d0)/
      if (r .eq. 0.d0) then
         nmax = lmax
         return
      end if
      z = r * dble(ri)
      sum = 0.d0
      do n = 1, nlim
         call init(n + lmax)
         call riccati_bessel(n + lmax, z, xi)
         do l = 0, n + lmax
            xi(l) = xi(l) / z * ci**l
         end do
         n21 = n + n + 1
         l = lmax
         c = fnr(n21) * fnr(l + l + 1) * ci**(n - l)
         call vcfunc(-1, n, 1, l, vc1)
         wmin = abs(n - l)
         wmax = n + l
         m = 1
         a = 0.
         b = 0.
         do w = wmin, wmax
            alnw = vc1(w) * vc1(w)
            if (mod(n + l + w, 2) .eq. 0) then
               a = a + alnw * xi(w)
            else
               b = b + alnw * xi(w)
            end if
         end do
         a = c * a
         b = c * b
         sum = sum + a * conjg(a) + b * conjg(b)
         if (abs(1.d0 - sum) .lt. eps) exit
      end do
      nmax = min(n, nlim)
      nmax = max(nmax, lmax)
   end subroutine tranordertest
!
!  address for axial translation coefficient
!
!
!  last revised: 15 January 2011
!

   pure integer function atcadd(m, n, ntot)
      implicit none
      integer, intent(in) :: m, n, ntot
      atcadd = n - ntot + (max(1, m) * (1 + 2 * ntot - max(1, m))) / 2 + ntot * min(1, m)
   end function atcadd

   pure integer function atcdim(ntot, ltot)
      implicit none
      integer, intent(in) :: ntot, ltot
      integer :: nmin, nmax
      nmin = min(ntot, ltot)
      nmax = max(ntot, ltot)
      atcdim = 2 * (nmin * (1 - nmin * nmin + 3 * nmax * (2 + nmin))) / 3
   end function atcdim
!
! the offset (integer) for the ntot X ltot translation matrix for degree m
!
   pure integer function moffset(m, ntot, ltot)
      implicit none
      integer, intent(in) :: m, ntot, ltot
      if (m .eq. 0) then
         moffset = 0
      elseif (m .lt. 0) then
         moffset = 2 * (-((1 + m) * (2 + m) * (3 + 2 * m + 3 * ntot)) &
                        - 3 * ltot * (2 + ntot + m * (3 + m + 2 * ntot))) / 3
      else
         moffset = 2 * (-3 * ltot * (-1 + m)**2 + 6 * ltot * m * ntot &
                        + (-1 + m) * (m * (-4 + 2 * m - 3 * ntot) + 3 * (1 + ntot))) / 3
      end if
   end function moffset
!
! cartesian_to_spherical takes the Cartesian point (x,y,z) = xp(1), xp(2), xp(3)
! and converts to polar form: r: radius, ct: cos(theta), ep = exp(i phi)
!
!
!  last revised: 15 January 2011
!
   pure subroutine cartesian_to_spherical(xp, r, ct, ep)
      implicit none
      real(8), intent(in) :: xp(3)
      real(8), intent(out) :: r, ct
      complex(8), intent(out) :: ep
      r = xp(1) * xp(1) + xp(2) * xp(2) + xp(3) * xp(3)
      if (r .eq. 0.d0) then
         ct = 1.d0
         ep = (1.d0, 0.d0)
         return
      end if
      r = sqrt(r)
      ct = xp(3) / r
      if (xp(1) .eq. 0.d0 .and. xp(2) .eq. 0.d0) then
         ep = (1.d0, 0.d0)
      else
         ep = cmplx(xp(1), xp(2), kind=kind(0.0d0)) / sqrt(xp(1) * xp(1) + xp(2) * xp(2))
      end if
      return
   end subroutine cartesian_to_spherical

   pure subroutine cartesian_vectors_to_spherical(nt, xp, xps)
      implicit none
      integer, intent(in) :: nt
      real(8), intent(in) :: xp(3, nt)
      real(8), intent(out) :: xps(3, nt)
      integer :: i
      real(8) :: r, ct, phi
      do concurrent(i=1:nt) local(r, ct, phi)
         r = xp(1, i) * xp(1, i) + xp(2, i) * xp(2, i) + xp(3, i) * xp(3, i)
         if (r .eq. 0.d0) then
            ct = 1.d0
            phi = 0.d0
         else
            r = sqrt(r)
            ct = xp(3, i) / r
            if (xp(1, i) .eq. 0.d0 .and. xp(2, i) .eq. 0.d0) then
               phi = 0.d0
            else
               phi = datan2(xp(2, i), xp(1, i))
            end if
         end if
         xps(:, i) = (/ct, phi, r/)
      end do
      return
   end subroutine cartesian_vectors_to_spherical

!
! euler rotation of a point (x,y,z) = xp(1), xp(2), xp(3)
! November 2012
!
   subroutine euler_rotate_cartesian_vectors(xp, eulerangf, dir, xprot, num)
      implicit none
      integer :: dir, n, i
      integer, optional :: num
      real(8) :: xp(3, *), eulerangf(3), eulerang(3), cang(3), sang(3), &
                 mat1(3, 3), mat2(3, 3), mat3(3, 3), xprot(3, *), xpt(3)
      if (present(num)) then
         n = num
      else
         n = 1
      end if
      if (dir .eq. 1) then
         eulerang = eulerangf
      else
         eulerang(1:3) = -eulerangf(3:1:-1)
      end if
      cang = cos(eulerang)
      sang = sin(eulerang)
      mat1(1, :) = (/cang(1), sang(1), 0.d0/)
      mat1(2, :) = (/-sang(1), cang(1), 0.d0/)
      mat1(3, :) = (/0.d0, 0.d0, 1.d0/)
      mat2(1, :) = (/cang(2), 0.d0, -sang(2)/)
      mat2(2, :) = (/0.d0, 1.d0, 0.d0/)
      mat2(3, :) = (/sang(2), 0.d0, cang(2)/)
      mat3(1, :) = (/cang(3), sang(3), 0.d0/)
      mat3(2, :) = (/-sang(3), cang(3), 0.d0/)
      mat3(3, :) = (/0.d0, 0.d0, 1.d0/)
      do i = 1, n
         xpt = xp(:, i)
         xpt = matmul(mat1, xpt)
         xpt = matmul(mat2, xpt)
         xpt = matmul(mat3, xpt)
         xprot(:, i) = xpt
      end do
   end subroutine euler_rotate_cartesian_vectors
!
! azimuthal_phase_factors returns the complex array epm(m) = exp(i m phi) for
! m=-nodr,nodr.   ep =exp(i phi), and epm is dimensioned epm(-nd:nd)
!
!
!  last revised: 15 January 2011
!
   subroutine azimuthal_phase_factors(ep, nodr, epm)
      implicit none
      integer :: nodr, m
      complex(8) :: ep, epm(-nodr:nodr)
      epm(0) = (1.d0, 0.d0)
      do m = 1, nodr
         epm(m) = ep * epm(m - 1)
         epm(-m) = conjg(epm(m))
      end do
      return
   end subroutine azimuthal_phase_factors
!
!  test to determine max order of vswf expansion of a plane wave at distance r
!
!
!  last revised: 15 January 2011
!
   subroutine planewavetruncationorder(r, rimedium, eps, nodr)
      implicit none
      integer :: nodr, n1, n
      real(8) :: r, eps, err
      complex(8), allocatable :: jn(:)
      complex(8) :: sum, ci, eir, rimedium(2), rri, rib
      data ci/(0.d0, 1.d0)/
      rib = 2.d0 / (1.d0 / rimedium(1) + 1.d0 / rimedium(2))
      n1 = max(10, int(3.*r + 1))
      allocate (jn(0:n1))
      rri = r * rib
      call riccati_bessel(n1, rri, jn)
      jn(0:n1) = jn(0:n1) / rri
      eir = exp(-ci * rri)
      sum = jn(0) * eir
      do n = 1, n1
         sum = sum + ci**n * dble(n + n + 1) * jn(n) * eir
         err = abs(1.d0 - sum)
         if (err .lt. eps) then
            nodr = n
            deallocate (jn)
            return
         end if
      end do
      nodr = n1
      deallocate (jn)
   end subroutine planewavetruncationorder
end module angular_functions
