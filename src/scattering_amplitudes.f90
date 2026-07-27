module scattering_amplitudes
   use constants
   use fft_translation
   use intrinsics
   use mie
   use mpidefs
   use numerical_tables
   use periodic_lattice_subroutines
   use specialfuncs
   use spheredata
   use surface_subroutines
   implicit none
   private
   public :: amplitude_to_scattering_matrix, common_origin_amplitude_matrix, &
             common_origin_scattering_matrix, evaluate_fixed_orientation_scattering_matrix, &
             fixed_orientation_scattering_matrix_expansion, multiple_origin_amplitude_matrix, &
             multiple_origin_scattering_matrix, &
             numerical_scattering_matrix_azimuthal_average_multiple_origin, &
             numerical_scattering_matrix_azimuthal_average_single_origin, periodic_lattice_scattering, &
             s11_expansion
contains

   subroutine periodic_lattice_scattering(amnp, qsca, scat_mat, krho_vec, num_dirs, dry_run)
      implicit none
      logical :: prop, calcsmat, dryrun
      logical, optional :: dry_run
      integer :: nmax, dir, n, i, q, p, ix, iy, nang, istrt
      integer, optional :: num_dirs(2)
      real(8) :: qsca(2, 2), k0x, k0y, wx, wy, targetz, kx, ky, krho, phi, smscale
      real(8), optional :: scat_mat(32, *), krho_vec(2, *)
      complex(8) :: amnp(*), ri, s, kz, sa(4)

      if (present(dry_run)) then
         dryrun = dry_run
      else
         dryrun = .false.
      end if
      calcsmat = present(scat_mat)
      k0x = incident_lateral_vector(1)
      k0y = incident_lateral_vector(2)
      wx = cell_width(1)
      wy = cell_width(2)
      nmax = ceiling(maxval(cell_width(1:2)) / pi)
      if (present(num_dirs)) num_dirs = 0
      smscale = 1.d0 / cross_section_radius**2 * 16.d0 * four_pi / wx / wy
      do dir = 1, 2
         if (dir .eq. 1) then
            istrt = 17
         else
            istrt = 1
         end if
         if (dir .eq. 1) then
            targetz = top_boundary
            ri = layer_ref_index(number_plane_boundaries)
         else
            targetz = bot_boundary
            ri = layer_ref_index(0)
         end if
         kx = k0x
         ky = k0y
         krho = sqrt(kx * kx + ky * ky)
         s = krho
         kz = ri * sqrt(1.d0 - s * s / ri / ri)
         if (kx .eq. 0.d0 .and. ky .eq. 0.d0) then
            phi = 0.d0
         else
            phi = datan2(ky, kx)
         end if
         nang = 1
         if (.not. dryrun) then
            call multiple_origin_amplitude_matrix(amnp, s, phi, targetz, dir, sa)
            qsca(1, dir) = kz * (abs(sa(2))**2 + abs(sa(4))**2)
            qsca(2, dir) = kz * (abs(sa(1))**2 + abs(sa(3))**2)
            if (calcsmat) then
               sa = sa * sqrt(kz)
               call amplitude_to_scattering_matrix(sa, scat_mat(istrt:istrt + 15, nang))
               scat_mat(istrt:istrt + 15, nang) = scat_mat(istrt:istrt + 15, nang) * smscale
               if (present(krho_vec)) then
                  krho_vec(1, nang) = kx
                  krho_vec(2, nang) = ky
               end if
            end if
         end if

         do n = 1, nmax
            prop = .false.
            do i = 0, 8 * n - 1
               q = i / (2 * n)
               p = i - 2 * q * n
               if (q .eq. 0) then
                  ix = n
                  iy = -n + p
               elseif (q .eq. 1) then
                  ix = n - p
                  iy = n
               elseif (q .eq. 2) then
                  ix = -n
                  iy = n - p
               else
                  ix = -n + p
                  iy = -n
               end if
               kx = two_pi * dble(ix) / wx + k0x
               ky = two_pi * dble(iy) / wy + k0y
               krho = sqrt(kx * kx + ky * ky)
               if (krho .le. abs(ri)) then
!                  if(krho.le.1.d0) then
                  prop = .true.
                  s = krho
                  kz = ri * sqrt(1.d0 - s * s / ri / ri)
!                     kz=sqrt(1.d0-s*s)
                  if (kx .eq. 0.d0 .and. ky .eq. 0.d0) then
                     phi = 0.d0
                  else
                     phi = datan2(ky, kx)
                  end if
                  nang = nang + 1
                  if (.not. dryrun) then
                     call multiple_origin_amplitude_matrix(amnp, s, phi, targetz, dir, sa)
                     qsca(1, dir) = qsca(1, dir) + kz * (abs(sa(2))**2 + abs(sa(4))**2)
                     qsca(2, dir) = qsca(2, dir) + kz * (abs(sa(1))**2 + abs(sa(3))**2)
                     if (calcsmat) then
                        sa = sa * sqrt(kz)
                        call amplitude_to_scattering_matrix(sa, scat_mat(istrt:istrt + 15, nang))
                        scat_mat(istrt:istrt + 15, nang) = scat_mat(istrt:istrt + 15, nang) * smscale
                        if (present(krho_vec)) then
                           krho_vec(1, nang) = kx
                           krho_vec(2, nang) = ky
                        end if
                     end if
                  end if
               end if
            end do
            if (.not. prop) exit
         end do
         if (present(num_dirs)) num_dirs(3 - dir) = nang
      end do
      qsca = qsca / cross_section_radius**2 * 16.d0 * four_pi / wx / wy
   end subroutine periodic_lattice_scattering

   subroutine multiple_origin_amplitude_matrix(amnp, s, phi, targetz, dir, sa)
      implicit none
      integer :: p, i, dir
      real(8) :: phi, targetz
      complex(8) :: amnp(number_eqns, 2), sa(4), sat(4), s
      complex(8), allocatable :: pmnpi(:, :), amnpi(:, :)

      sa = 0.d0
      do i = 1, number_spheres
         if (host_sphere(i) .ne. 0) cycle
         allocate (pmnpi(sphere_block(i), 2), amnpi(sphere_block(i), 2))
         call layer_vector_spherical_harmonics(s, phi, targetz, dir, sphere_position(:, i), sphere_order(i), pmnpi)
         amnpi(:, :) = amnp(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i), :)
         do p = 1, 2
            call left_right_mode_transformation(sphere_order(i), amnpi(:, p), amnpi(:, p))
         end do
         sat(1) = sum(pmnpi(:, 2) * amnpi(:, 2)) * 0.5d0
         sat(2) = sum(pmnpi(:, 1) * amnpi(:, 1)) * 0.5d0
         sat(3) = -sum(pmnpi(:, 1) * amnpi(:, 2)) * 0.5d0
         sat(4) = -sum(pmnpi(:, 2) * amnpi(:, 1)) * 0.5d0
         sa(:) = sa(:) + sat(:)
         deallocate (pmnpi, amnpi)
      end do
   end subroutine multiple_origin_amplitude_matrix

   subroutine common_origin_amplitude_matrix(amnp, s, phi, targetz, dir, sa)
      implicit none
      integer :: dir
      real(8) :: phi, targetz
      complex(8) :: amnp(2 * t_matrix_order * (t_matrix_order + 2), 2), sa(4), s, pmnp(2 * t_matrix_order * (t_matrix_order + 2), 2)
      call layer_vector_spherical_harmonics(s, phi, targetz, dir, cluster_origin, t_matrix_order, pmnp)
      sa(1) = 0.5d0 * sum(pmnp(:, 2) * amnp(:, 2))
      sa(2) = 0.5d0 * sum(pmnp(:, 1) * amnp(:, 1))
      sa(3) = -0.5d0 * sum(pmnp(:, 1) * amnp(:, 2))
      sa(4) = -0.5d0 * sum(pmnp(:, 2) * amnp(:, 1))
   end subroutine common_origin_amplitude_matrix

   pure subroutine amplitude_to_scattering_matrix(sa, sm)
      implicit none
      complex(8), intent(in) :: sa(4)
      real(8), intent(out) :: sm(4, 4)
      integer :: i, j
      complex(8) :: sp(4, 4)
      do concurrent(i=1:4, j=1:4)
         sp(i, j) = sa(i) * conjg(sa(j))
      end do
      sm(1, 1) = sp(1, 1) + sp(2, 2) + sp(3, 3) + sp(4, 4)
      sm(1, 2) = -sp(1, 1) + sp(2, 2) - sp(3, 3) + sp(4, 4)
      sm(2, 1) = -sp(1, 1) + sp(2, 2) + sp(3, 3) - sp(4, 4)
      sm(2, 2) = sp(1, 1) + sp(2, 2) - sp(3, 3) - sp(4, 4)
      sm(3, 3) = 2.*(sp(1, 2) + sp(3, 4))
      sm(3, 4) = 2.*aimag(sp(2, 1) + sp(4, 3))
      sm(4, 3) = 2.*aimag(sp(1, 2) - sp(3, 4))
      sm(4, 4) = 2.*(sp(1, 2) - sp(3, 4))
      sm(1, 3) = 2.*(sp(2, 3) + sp(1, 4))
      sm(3, 1) = 2.*(sp(2, 4) + sp(1, 3))
      sm(1, 4) = -2.*aimag(sp(2, 3) - sp(1, 4))
      sm(4, 1) = -2.*aimag(sp(4, 2) + sp(1, 3))
      sm(2, 3) = -2.*(sp(2, 3) - sp(1, 4))
      sm(3, 2) = -2.*(sp(2, 4) - sp(1, 3))
      sm(2, 4) = -2.*aimag(sp(2, 3) + sp(1, 4))
      sm(4, 2) = -2.*aimag(sp(4, 2) - sp(1, 3))
   end subroutine amplitude_to_scattering_matrix

   subroutine multiple_origin_scattering_matrix(amnp, ct, phi, csca, sa, sm, rotate_plane, s11_only)
      implicit none
      logical :: rotate, s11only
      logical, optional :: rotate_plane, s11_only
      integer :: dir, nelem
      real(8) :: ct, phi, sm(*), csca, targetz
      complex(8) :: amnp(number_eqns, 2), sa(4), ri, s, amnpt(number_eqns, 2)
      if (present(rotate_plane)) then
         rotate = rotate_plane
      else
         rotate = .false.
      end if
      if (present(s11_only)) then
         s11only = s11_only
      else
         s11only = .false.
      end if
      if (s11only) then
         nelem = 2
      else
         nelem = 16
      end if

      sa = 0.d0
      if (ct .le. 0.d0) then
         ri = layer_ref_index(0)
         dir = 2
         targetz = bot_boundary
      else
         ri = layer_ref_index(number_plane_boundaries)
         dir = 1
         targetz = top_boundary
      end if
      if (rotate) then
         amnpt(:, 1) = amnp(:, 1) * cos(phi) + amnp(:, 2) * sin(phi)
         amnpt(:, 2) = amnp(:, 1) * sin(phi) - amnp(:, 2) * cos(phi)
      else
         amnpt = amnp
      end if
      s = dble(ri) * sqrt((1.d0 - ct) * (1.d0 + ct))
      call multiple_origin_amplitude_matrix(amnpt, s, phi, targetz, dir, sa)
!         sa=sa*ri*ct/sqrt(csca/16.d0/pi)
! 10-22 patch
      sa = sa * ri * ri * ct / sqrt(csca / 16.d0 / pi)
      if (s11only) then
         sm(1) = abs(sa(2))**2.+abs(sa(4))**2.
         sm(2) = abs(sa(1))**2.+abs(sa(3))**2.
      else
         call amplitude_to_scattering_matrix(sa, sm)
      end if
! 10-22 patch
      sm(1:nelem) = sm(1:nelem) / dble(ri)
   end subroutine multiple_origin_scattering_matrix
!
!  scattering amplitude sa and matrix sm calculation
!
!  original: 15 January 2011
!  revised: 21 February 2011: S11 normalization changed
!  april 2013: moved things around to try to get it to work.
!
   subroutine common_origin_scattering_matrix(amn0, nodrt, ct, phi, sa, sm, rotate_plane, normalize_s11, &
                                              s11_only)
      implicit none
      logical, optional :: rotate_plane, normalize_s11, s11_only
      logical :: rotate, norms11, s11only
      integer :: nodrt, m, n, p, m1, n1, nelem
      real(8) :: ct, phi, sm(*), cphi, sphi, qsca, tau(0:nodrt + 1, nodrt, 2)
      complex(8) :: amn0(0:nodrt + 1, nodrt, 2, 2), sa(4), ephi, ephim(-nodrt:nodrt), &
                    ci, cin, a, b
      data ci/(0.d0, 1.d0)/
      if (present(rotate_plane)) then
         rotate = rotate_plane
      else
         rotate = .false.
      end if
      if (present(normalize_s11)) then
         norms11 = normalize_s11
      else
         norms11 = .true.
      end if
      if (present(s11_only)) then
         s11only = s11_only
      else
         s11only = .false.
      end if
      if (s11only) then
         nelem = 2
      else
         nelem = 16
      end if
      call vector_spherical_harmonics(ct, nodrt, tau)
      cphi = cos(phi)
      sphi = sin(phi)
      ephi = cmplx(cphi, sphi, kind=kind(0.0d0))
      call azimuthal_phase_factors(ephi, nodrt, ephim)
      sa = (0.d0, 0.d0)
      qsca = 0.d0
      do n = 1, nodrt
         cin = (-ci)**n
         do m = -n, n
            if (m .le. -1) then
               m1 = n + 1
               n1 = -m
            else
               m1 = m
               n1 = n
            end if
            do p = 1, 2
               qsca = qsca + amn0(m1, n1, p, 1) * conjg(amn0(m1, n1, p, 1)) &
                      + amn0(m1, n1, p, 2) * conjg(amn0(m1, n1, p, 2))
               if (rotate) then
                  a = amn0(m1, n1, p, 1) * cphi + amn0(m1, n1, p, 2) * sphi
                  b = amn0(m1, n1, p, 1) * sphi - amn0(m1, n1, p, 2) * cphi
               else
                  a = amn0(m1, n1, p, 1)
                  b = -amn0(m1, n1, p, 2)
               end if
               sa(1) = sa(1) + cin * tau(m1, n1, 3 - p) * b * ephim(m)
               sa(2) = sa(2) + ci * cin * tau(m1, n1, p) * a * ephim(m)
               sa(3) = sa(3) + ci * cin * tau(m1, n1, p) * b * ephim(m)
               sa(4) = sa(4) + cin * tau(m1, n1, 3 - p) * a * ephim(m)
            end do
         end do
      end do
      qsca = qsca * 2.d0
      if (.not. norms11) qsca = 1.d0 / pi
      sa = sa * 4.d0 / sqrt(2.d0 * qsca)
      if (s11only) then
         sm(1) = abs(sa(2))**2.+abs(sa(4))**2.
         sm(2) = abs(sa(1))**2.+abs(sa(3))**2.
      else
         call amplitude_to_scattering_matrix(sa, sm)
      end if
! patch 10-22
      sm(1:nelem) = sm(1:nelem) / four_pi / dble(layer_ref_index(0))
   end subroutine common_origin_scattering_matrix

   subroutine numerical_scattering_matrix_azimuthal_average_single_origin(amn0, nodrt, ct, sm, &
                                                                          rotate_plane, normalize_s11, &
                                                                          number_angles, s11_only)
      implicit none
      logical :: rotate, norms11, s11only
      logical, optional :: rotate_plane, normalize_s11, s11_only
      integer :: i, numang, nodrt, nelem
      integer, optional :: number_angles
      real(8) :: ct, phi, sm(*), smt(16)
      complex(8) :: amn0(0:nodrt + 1, nodrt, 2, 2), sa(4)
      if (present(rotate_plane)) then
         rotate = rotate_plane
      else
         rotate = .false.
      end if
      if (present(normalize_s11)) then
         norms11 = normalize_s11
      else
         norms11 = .true.
      end if
      if (present(s11_only)) then
         s11only = s11_only
      else
         s11only = .false.
      end if
      if (present(number_angles)) then
         numang = number_angles
      else
         numang = 2 * nodrt + 2
      end if
      if (s11only) then
         nelem = 2
      else
         nelem = 16
      end if
      sm(1:nelem) = 0.
      do i = 1, numang
         phi = two_pi * dble(i - 1) / dble(numang)
         call common_origin_scattering_matrix(amn0, nodrt, ct, phi, sa, smt, &
                                              rotate_plane=rotate, normalize_s11=norms11, s11_only=s11only)
         sm(1:nelem) = sm(1:nelem) + smt(1:nelem)
      end do
      sm(1:nelem) = sm(1:nelem) / dble(numang)
   end subroutine numerical_scattering_matrix_azimuthal_average_single_origin

   subroutine numerical_scattering_matrix_azimuthal_average_multiple_origin(amnp, ct, sm, &
                                                                            number_angles, rotate_plane, s11_only)
      implicit none
      logical :: s11only, rotate
      logical, optional :: s11_only, rotate_plane
      integer :: i, numang, nelem
      integer, optional :: number_angles
      real(8) :: ct, phi, sm(*), smt(16), csca
      complex(8) :: amnp(number_eqns, 2), sa(4)
      if (present(number_angles)) then
         numang = number_angles
      else
         numang = 2 * t_matrix_order + 2
      end if
      if (present(s11_only)) then
         s11only = s11_only
      else
         s11only = .false.
      end if
      if (present(rotate_plane)) then
         rotate = rotate_plane
      else
         rotate = .false.
      end if
      if (s11only) then
         nelem = 2
      else
         nelem = 16
      end if
      sm(1:nelem) = 0.
      csca = two_pi
      do i = 1, numang
         phi = two_pi * dble(i - 1) / dble(numang)
         call multiple_origin_scattering_matrix(amnp, ct, phi, csca, sa, smt, &
                                                rotate_plane=rotate, s11_only=s11only)
         sm(1:nelem) = sm(1:nelem) + smt(1:nelem)
      end do
      sm(1:nelem) = sm(1:nelem) / dble(numang)
   end subroutine numerical_scattering_matrix_azimuthal_average_multiple_origin

!   c                                                                               c
!   c  subroutine scatexp(amn0,nodrt,nodrg,gmn) computes the expansion coefficients c
!   c  for the spherical harmonic expansion of the scattering phase function from   c
!   c  the scattering coefficients amn0.  For a complete expansion, the max. order  c
!   c  of the phase function expansion (nodrg) will be 2*nodrt, where nodrt is      c
!   c  the max. order of the scattered field expansion.   In this code nodrg is     c
!   c  typically set to 1, so that the subroutine returns the first moments         c
!   c  of the phase function; gmn(1) and gmn(2).                                    c
!   c                                                                               c
!   c  The expansion coefficients are normalized so that gmn(0)=1                   c
!   c                                                                               c
!   c  gmn(1)/3 is the asymmetry parameter.                                         c
!   c                                                                               c
   subroutine s11_expansion(amn0, nodrt, mmax, nodrg, gmn)
      implicit none
      integer :: nodrt, m, n, p, ma, na, mmax, nodrg, w, w1, w2, u, uw, ww1, &
                 l1, l2, ka, la, k, l, q, ik
      real(8) :: vc1(0:nodrt * 2 + 1), vc2(0:nodrt * 2 + 1), g0
      complex(8) :: amn0(0:nodrt + 1, nodrt, 2, 2), gmn(0:nodrg * (nodrg + 3) / 2), &
                    a(2, 2), c, c2
      gmn = (0.d0, 0.d0)
      do n = 1, nodrt
         l1 = max(1, n - nodrg)
         l2 = min(nodrt, n + nodrg)
         do l = l1, l2
            c = sqrt(dble((n + n + 1) * (l + l + 1))) &
                * cmplx(0.d0, 1.d0, kind=kind(0.0d0))**(l - n)
            w2 = min(n + l, nodrg)
            call vector_coupling_coefficients(-1, l, 1, n, vc2)
            do m = -n, n
               if (m .le. -1) then
                  ma = n + 1
                  na = -m
               else
                  ma = m
                  na = n
               end if
               do k = -l, min(l, m)
                  if (k .le. -1) then
                     ka = l + 1
                     la = -k
                  else
                     ka = k
                     la = l
                  end if
                  u = m - k
                  if (u .le. mmax) then
                     ik = (-1)**k
                     c2 = ik * c
                     do p = 1, 2
                        do q = 1, 2
                           a(p, q) = c2 * (amn0(ma, na, p, 1) * conjg(amn0(ka, la, q, 1)) &
                                           + amn0(ma, na, p, 2) * conjg(amn0(ka, la, q, 2)))
                        end do
                     end do
                     w1 = max(abs(n - l), abs(u))
                     w2 = min(n + l, nodrg)
                     call vector_coupling_coefficients(-k, l, m, n, vc1)
                     do w = w1, w2
                        uw = (w * (w + 1)) / 2 + u
                        do p = 1, 2
                           if (mod(n + l + w, 2) .eq. 0) then
                              q = p
                           else
                              q = 3 - p
                           end if
                           gmn(uw) = gmn(uw) - vc1(w) * vc2(w) * a(p, q)
                        end do
                     end do
                  end if
               end do
            end do
         end do
      end do
      g0 = dble(gmn(0))
      gmn(0) = 1.d0
      do w = 1, nodrg
         ww1 = (w * (w + 1)) / 2
         gmn(ww1) = cmplx(dble(gmn(ww1)), 0.d0, kind=kind(0.0d0)) / g0
         do u = 1, min(mmax, w)
            uw = ww1 + u
            gmn(uw) = (-1)**u * 2.d0 * gmn(uw) / g0
         end do
      end do
   end subroutine s11_expansion
!
!  calculate azimuth--averaged scattering matrix from expansion, for cos(theta) = ct
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!  this is currently not used in v. 3.0
!
   subroutine evaluate_fixed_orientation_scattering_matrix(ntot, s00, s02, sp22, sm22, ct, sm, normalize_s11)
      logical :: norms11
      logical, optional :: normalize_s11
      integer :: ntot, w, ww1
      real(8) :: s00(4, 4, 0:ntot * 2), s02(4, 4, 0:ntot * 2), sp22(4, 4, 0:ntot * 2), sm22(4, 4, 0:ntot * 2), &
                 sm(4, 4), dc(-2:2, 0:2 * ntot * (2 * ntot + 2)), ct, temp
      if (present(normalize_s11)) then
         norms11 = normalize_s11
      else
         norms11 = .true.
      end if
!if(light_up) then
!write(*,'('' foc1 '',3es13.5)') ct
!flush(6)
!endif
      call rotation_coefficients(ct, 2, 2 * ntot, dc)

      sm = 0.d0
      do w = 0, 2 * ntot
         ww1 = w * (w + 1)
         sm(:, :) = sm(:, :) + s00(:, :, w) * dc(0, ww1) + s02(:, :, w) * dc(0, ww1 + 2) &
                    + sp22(:, :, w) * dc(2, ww1 + 2) + sm22(:, :, w) * dc(-2, ww1 + 2)
      end do
!if(light_up) then
!write(*,'('' foc2 '',3es13.5)') ct,s00(1,1,0),sm(1,1)
!flush(6)
!endif
      if (norms11) then
         if (abs(s00(1, 1, 0)) .lt. 1.d-10) then
            sm = 0.d0
         else
            sm = sm / s00(1, 1, 0)
         end if
      end if
!
!  a patch
!
      sm(3, 1) = -sm(3, 1)
      sm(1, 3) = -sm(1, 3)
      sm(4, 3) = -sm(4, 3)
      temp = sm(4, 1)
      sm(4, 1) = sm(4, 2)
      sm(4, 2) = temp

      sm(1, 2) = (sm(1, 2) + sm(2, 1)) / 2.d0
      sm(2, 1) = sm(1, 2)
      sm(3, 4) = (sm(3, 4) - sm(4, 3)) / 2.d0
      sm(4, 3) = -sm(3, 4)

!         do i=1,4
!            do j=1,4
!               if(i.ne.1.or.j.ne.1) then
!                  sm(i,j)=sm(i,j)/sm(1,1)
!               endif
!            enddo
!         enddo
!if(light_up) then
!write(*,'('' foc3 '',i3)') mstm_global_rank
!flush(6)
!endif
   end subroutine evaluate_fixed_orientation_scattering_matrix
!
!  determine the generalized spherical function expansion for the azimuth-averaged scattering matrix
!  corresponding to the target-based scattering field expansion of amnp.
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: fixed flush call.
!

   subroutine fixed_orientation_scattering_matrix_expansion(ntot, amnp, s00, s02, sp22, sm22, mpi_comm)
      integer :: ntot, n, p, m, l, wmin, wmax, m1m, q, m1mq, m1mnpl, w, m1w, i, wtot, j, &
                 rank, numprocs, nsend, runprintunit, mpicomm, task
      integer, optional :: mpi_comm
      real(8) :: s00(4, 4, 0:ntot * 2), s02(4, 4, 0:ntot * 2), sp22(4, 4, 0:ntot * 2), sm22(4, 4, 0:ntot * 2), &
                 cm1p1(0:ntot * 2), cm1m1(0:ntot * 2), cmmpm(0:ntot * 2), cmmm2pm(0:ntot * 2), &
                 cmmp2pm(0:ntot * 2)
      complex(8) :: amnp(0:ntot + 1, ntot, 2, 2), a1(-ntot - 2:ntot + 2, ntot, 2), a2(-ntot - 2:ntot + 2, ntot, 2), &
                    fnl, a1122, a2112, a1p2, a1m2
      complex(8), parameter :: ci = (0.d0, 1.d0)
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call init(2 * ntot)
      runprintunit = run_print_unit
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      a1 = (0.d0, 0.d0)
      a2 = (0.d0, 0.d0)
      do n = 1, ntot
         do p = 1, 2
            do m = -n, -1
               a1(m, n, p) = amnp(n + 1, -m, p, 1)
               a2(m, n, p) = amnp(n + 1, -m, p, 2)
            end do
            do m = 0, n
               a1(m, n, p) = amnp(m, n, p, 1)
               a2(m, n, p) = amnp(m, n, p, 2)
            end do
         end do
      end do
      s00 = 0.d0
      s02 = 0.d0
      sp22 = 0.d0
      sm22 = 0.d0
      wtot = ntot + ntot
      task = 0
      do n = 1, ntot
         do l = 1, ntot
            wmin = abs(n - l)
            wmax = n + l
            fnl = sqrt(dble((n + n + 1) * (l + l + 1))) * ci**(l - n)
            call vector_coupling_coefficients(-1, n, 1, l, cm1p1)
            call vector_coupling_coefficients(-1, n, -1, l, cm1m1)
            do m = -min(n, l + 2), min(n, l + 2)
               m1m = (-1)**m
               if (abs(m) .le. l) then
                  call vector_coupling_coefficients(-m, n, m, l, cmmpm)
               else
                  cmmpm = 0.d0
               end if
               if (abs(-2 + m) .le. l) then
                  call vector_coupling_coefficients(-m, n, -2 + m, l, cmmm2pm)
               else
                  cmmm2pm = 0.d0
               end if
               if (abs(2 + m) .le. l) then
                  call vector_coupling_coefficients(-m, n, 2 + m, l, cmmp2pm)
               else
                  cmmp2pm = 0.d0
               end if
               do p = 1, 2
                  do q = 1, 2
                     m1mq = (-1)**(m + q)
                     m1mnpl = (-1)**(m + n + p + l)
                     a1122 = (a1(m, n, p) * conjg(a1(m, l, q)) + a2(m, n, p) * conjg(a2(m, l, q)))
                     a2112 = (a2(m, n, p) * conjg(a1(m, l, q)) - a1(m, n, p) * conjg(a2(m, l, q)))
                     a1p2 = (a1(m, n, p) + ci * a2(m, n, p)) * conjg(a1(m - 2, l, q) - ci * a2(m - 2, l, q))
                     a1m2 = (a1(m, n, p) - ci * a2(m, n, p)) * conjg(a1(m + 2, l, q) + ci * a2(m + 2, l, q))
                     do w = wmin, wmax
                        task = task + 1
                        if (mod(task - 1, numprocs) .ne. rank) cycle
                        m1w = (-1)**w
                        if (mod(n + l + w + p + q, 2) .eq. 0) then
                           s00(1, 1, w) = s00(1, 1, w) - (m1m * fnl * a1122 * cm1p1(w) * cmmpm(w)) / 2.
                        else
                           s00(4, 4, w) = s00(4, 4, w) + (ci / 2.*m1m * fnl * a2112 * cm1p1(w) * cmmpm(w))
                        end if
                        if (w .lt. 2) cycle
                        if (mod(n + l + w + p + q, 2) .eq. 0) then
                           s02(2, 1, w) = s02(2, 1, w) - (m1mq * a1122 * fnl * cm1m1(w) * cmmpm(w)) / 2.
                           s02(1, 2, w) = s02(1, 2, w) - (m1m * a1p2 * fnl * cm1p1(w) * cmmm2pm(w)) / 4.
                           s02(1, 2, w) = s02(1, 2, w) - (m1m * a1m2 * fnl * cm1p1(w) * cmmp2pm(w)) / 4.
                        else
                           s02(3, 4, w) = s02(3, 4, w) + aimag(-ci / 2.*m1mq * a2112 * fnl * cm1m1(w) * cmmpm(w))
                           s02(4, 3, w) = s02(4, 3, w) + aimag(m1m * a1p2 * fnl * cm1p1(w) * cmmm2pm(w)) / 4.
                           s02(4, 3, w) = s02(4, 3, w) - aimag(m1m * a1m2 * fnl * cm1p1(w) * cmmp2pm(w)) / 4.
                        end if
                        sm22(2, 2, w) = sm22(2, 2, w) - (m1mnpl * m1w * a1p2 * fnl * cm1m1(w) * cmmm2pm(w)) / 8.
                        sp22(2, 2, w) = sp22(2, 2, w) - (m1mq * a1p2 * fnl * cm1m1(w) * cmmm2pm(w)) / 8.
                        sm22(2, 2, w) = sm22(2, 2, w) - (m1mq * a1m2 * fnl * cm1m1(w) * cmmp2pm(w)) / 8.
                        sp22(2, 2, w) = sp22(2, 2, w) - (m1mnpl * m1w * a1m2 * fnl * cm1m1(w) * cmmp2pm(w)) / 8.
                     end do
                  end do
               end do
            end do
         end do
      end do
      call mstm_mpi(mpi_command='barrier')
      nsend = 4 * 4 * (2 * ntot + 1)
      call mstm_mpi(mpi_command='allreduce', mpi_recv_buf_dp=s00, &
                    mpi_number=nsend, mpi_operation=mstm_mpi_sum, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='allreduce', mpi_recv_buf_dp=s02, &
                    mpi_number=nsend, mpi_operation=mstm_mpi_sum, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='allreduce', mpi_recv_buf_dp=sp22, &
                    mpi_number=nsend, mpi_operation=mstm_mpi_sum, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='allreduce', mpi_recv_buf_dp=sm22, &
                    mpi_number=nsend, mpi_operation=mstm_mpi_sum, mpi_comm=mpicomm)
!
!  a patch
!
      do i = 3, 4
         do j = 1, i
            s00(j, i, 0:wtot) = -s00(j, i, 0:wtot)
            s02(j, i, 0:wtot) = -s02(j, i, 0:wtot)
            sm22(j, i, 0:wtot) = -sm22(j, i, 0:wtot)
            sp22(j, i, 0:wtot) = -sp22(j, i, 0:wtot)
         end do
      end do
      sm22(3, 3, :) = -sm22(2, 2, :)
      sp22(3, 3, :) = sp22(2, 2, :)
!
! patch 10-22
      s00 = 0.5d0 * s00 / dble(layer_ref_index(0))
      s02 = 0.5d0 * s02 / dble(layer_ref_index(0))
      sm22 = 0.5d0 * sm22 / dble(layer_ref_index(0))
      sp22 = 0.5d0 * sp22 / dble(layer_ref_index(0))
!         deallocate(nlindex,nlnum)

   end subroutine fixed_orientation_scattering_matrix_expansion
end module scattering_amplitudes
