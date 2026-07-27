module effective_medium_analysis
   use constants
   use input_state
   implicit none
contains

   subroutine effective_extinction_coefficient_ratio(scacoef, abscoef, srat, arat)
      implicit none
      integer :: i
      real(8) :: miesca, mieabs, scaq, absq, h, scacoef, abscoef, root0, root, func, dfunc, srat, arat, across

!         miesca=(mean_qext_mie-mean_qabs_mie)*sphere_volume_fraction*3.d0/4.d0/length_scale_factor
      mieabs = (mean_qabs_mie) * sphere_volume_fraction * 3.d0 / 4.d0 / length_scale_factor
      miesca = (mean_qext_mie) * sphere_volume_fraction * 3.d0 / 4.d0 / length_scale_factor
      if (target_shape .le. 1) then
!            scaq=1.d0-0.5d0*(-sum(dif_boundary_sca(:,0))+sum(dif_boundary_sca(:,number_plane_boundaries+1)))
         scaq = 1.d0 - 0.5d0 * (-sum(dif_boundary_sca(:, 0)) + sum(dif_boundary_sca(:, 1))) - q_eff_tot(2, 1)
         absq = 1.d0 - q_eff_tot(2, 1)
         scaq = max(1.d-5, scaq)
         absq = max(1.d-5, absq)
         if (target_shape .eq. 0) then
            across = 4.d0 * product(target_dimensions(1:2)) * length_scale_factor**2
         elseif (target_shape .eq. 1) then
            across = pi * (target_dimensions(1) * length_scale_factor)**2
         end if
         h = four_pi_over_three * dble(number_spheres) * length_scale_factor**3 / (across * sphere_volume_fraction)
         scacoef = -dlog(scaq) / h
         if (abs(mieabs) .lt. 1.d-7) then
            abscoef = 0.d0
         else
            abscoef = -dlog(absq) / h
         end if
      else
!            scaq=0.5d0*(-sum(dif_boundary_sca(:,0))+sum(dif_boundary_sca(:,number_plane_boundaries+1)))
         scaq = 0.5d0 * (-sum(dif_boundary_sca(:, 0)) + sum(dif_boundary_sca(:, 1))) + q_eff_tot(2, 1)
         absq = q_eff_tot(2, 1)
         h = length_scale_factor * (target_dimensions(1) - 1.d0)**3 / target_dimensions(1)**2
         root0 = scaq
         do i = 1, 100
            root = root0
            func = 1.d0 - (1.d0 - exp(-2.d0 * root) * (1.d0 + 2.d0 * root)) / (2.d0 * root * root) - scaq
            dfunc = exp(-2.d0 * root) * (-1.d0 + exp(2.d0 * root) - 2.d0 * root * (1.d0 + root)) / root**3
            root0 = root - func / dfunc
            if (abs(1.d0 - root / root0) .lt. 1.d-6) exit
         end do
         scacoef = root / h
         if (abs(mieabs) .lt. 1.d-7) then
            abscoef = 0.d0
         else
            root0 = absq
            do i = 1, 100
               root = root0
               func = 1.d0 - (1.d0 - exp(-2.d0 * root) * (1.d0 + 2.d0 * root)) / (2.d0 * root * root) - absq
               dfunc = exp(-2.d0 * root) * (-1.d0 + exp(2.d0 * root) - 2.d0 * root * (1.d0 + root)) / root**3
               root0 = root - func / dfunc
               if (abs(1.d0 - root / root0) .lt. 1.d-6) exit
            end do
            abscoef = root / h
         end if
      end if
      srat = scacoef / miesca
      if (abs(mieabs) .lt. 1.d-7) then
         arat = 1.d0
      else
         arat = abscoef / mieabs
      end if
   end subroutine effective_extinction_coefficient_ratio

   subroutine effective_refractive_index(ndat, edat, d, rieff, &
                                         e0)
      use mpidefs
      implicit none
      integer :: ndat, i, rank
      real(8) :: d, xdat(ndat), phase(ndat), amplitude(ndat), &
                 phaseslope, phaseintercept, &
                 ampslope, ampintercept, oldphase, &
                 newphase
      complex(8) :: edat(ndat), rieff, e0
      call mstm_mpi(mpi_command='rank', mpi_rank=rank)

      oldphase = -1.d10
      do i = 1, ndat
         xdat(i) = d * dble(i - 1)
         phase(i) = datan2(aimag(edat(i)), dble(edat(i)))
         amplitude(i) = dlog(abs(edat(i)))
      end do
      oldphase = phase(1)
      do i = 2, ndat
         newphase = phase(i)
         do while (abs(newphase - oldphase) .gt. pi)
            if (newphase .gt. oldphase) then
               newphase = newphase - two_pi
            else
               newphase = newphase + two_pi
            end if
         end do
         phase(i) = newphase
         oldphase = newphase
      end do
      call linear_regression(ndat, phase, xdat, &
                             phaseslope, phaseintercept)
      call linear_regression(ndat, amplitude, xdat, &
                             ampslope, ampintercept)
      rieff = cmplx(phaseslope, -ampslope, kind=kind(0.0d0))
      e0 = dexp(ampintercept) * exp((0.d0, 1.d0) * phaseintercept)
   end subroutine effective_refractive_index

   pure subroutine linear_regression(ndat, fdat, xdat, a, b)
      implicit none
      integer, intent(in) :: ndat
      real(8), intent(in) :: fdat(ndat), xdat(ndat)
      real(8), intent(out) :: a, b
      real(8) :: xbar, x2bar, fbar, xfbar
      fbar = sum(fdat(1:ndat)) / dble(ndat)
      xbar = sum(xdat(1:ndat)) / dble(ndat)
      xfbar = sum(xdat(1:ndat) * fdat(1:ndat)) / dble(ndat)
      x2bar = sum(xdat(1:ndat) * xdat(1:ndat)) / dble(ndat)
      a = (xfbar - xbar * fbar) / (x2bar - xbar * xbar)
      b = (fbar * x2bar - xbar * xfbar) / (x2bar - xbar * xbar)
   end subroutine linear_regression

   subroutine effective_ref_index_fit(anp, rifit, xfit, info)
      use levenberg_marquardt
      implicit none
      integer :: info, n0, m0, fitorder
      real(8) :: parm(3), vec(4 * t_matrix_order), xfit, rii, rir, fv
      complex(8) :: anp(2, t_matrix_order), rifit, ri1(2), rim(2)

      ri1 = 1.d0
      rim = layer_ref_index(0)
      fitorder = min(80, t_matrix_order)
!         fitorder=2*max_mie_order
      m0 = 4 * fitorder
      if (fit_for_radius) then
         n0 = 3
      else
         n0 = 2
      end if
      if (random_configuration_host_model .eq. 1) then
         xfit = target_dimensions(1) * length_scale_factor
      elseif (random_configuration_host_model .eq. 2) then
         xfit = vol_radius / (sphere_volume_fraction)**0.33333
      end if
      fv = (vol_radius / (target_dimensions(1) * length_scale_factor))**3.
      if (random_configuration_host) then
         rifit = rim(1)
      else
         rii = mean_qext_mie * area_mean_radius**2 * 3.d0 * dble(number_spheres) / (4.d0 * xfit**3) / 2.d0
         rir = (1.d0 - fv) + fv * dble(ref_index_scale_factor)
         rifit = cmplx(rir, rii / 2.d0, kind=kind(0.0d0))
      end if
      parm = (/dble(rifit), aimag(rifit), xfit/)
      effective_fit_order = fitorder
      effective_fit_radius = xfit
      if (allocated(effective_fit_coefficients)) deallocate (effective_fit_coefficients)
      allocate (effective_fit_coefficients(2, fitorder))
      effective_fit_coefficients = anp(:, 1:fitorder)
      call lmdif1(effective_ref_index_residual, m0, n0, parm, vec, 1.d-6, info)
      deallocate (effective_fit_coefficients)
      rifit = cmplx(parm(1), parm(2), kind=kind(0.0d0))
      xfit = parm(3)
   end subroutine effective_ref_index_fit

   subroutine effective_ref_index_residual(ndat, nparm, xparm, fdat, iflag)
      implicit none
      integer :: nparm, ndat, iflag, n, i
      real(8) :: xparm(nparm), fdat(ndat), xsp, qext, qabs, qsca
      complex(8) :: ri(2), anpmie(2, 2, effective_fit_order), a(2), ri0(2)
      ri0(:) = layer_ref_index(0)
      ri(:) = cmplx(xparm(1), xparm(2), kind=kind(0.0d0))
      if (nparm .eq. 3) then
         xsp = xparm(3)
      else
         xsp = effective_fit_radius
      end if
      if (random_orientation) then
         call optically_active_mie_coefficients(xsp, ri, effective_fit_order, 0.d0, qext, qsca, qabs, &
                                                anp_mie=anpmie, ri_medium=ri0)
      elseif (effective_medium_simulation) then
         call optically_active_mie_coefficients(xsp, ri0, effective_fit_order, 0.d0, qext, qsca, qabs, &
                                                ri_medium=ri, anp_eff_mie=anpmie)
      else
         call optically_active_mie_coefficients(xsp, ri, effective_fit_order, 0.d0, &
                                                qext, qsca, qabs, anp_mie=anpmie)
      end if
      i = 1
      do n = 1, effective_fit_order
         a(1) = anpmie(1, 1, n) + anpmie(2, 1, n)
         a(2) = anpmie(1, 1, n) - anpmie(2, 1, n)
         fdat(i) = dble(effective_fit_coefficients(1, n) - a(1))
         fdat(i + 1) = aimag(effective_fit_coefficients(1, n) - a(1))
         fdat(i + 2) = dble(effective_fit_coefficients(2, n) - a(2))
         fdat(i + 3) = aimag(effective_fit_coefficients(2, n) - a(2))
         i = i + 4
      end do
   end subroutine effective_ref_index_residual

   subroutine diffuse_scattering_effective_ref_index(a, qe, extrat)
      implicit none
      real(8) :: a, extrat, tvol, rc, ds, dela, qe, delextrat, extrat0
      call calculate_target_volume(target_dimensions, tvol)
      tvol = tvol * length_scale_factor**3
      rc = (tvol * 3.d0 / 4.d0 / pi)**(1.d0 / 3.d0)
      qe = (dif_csca_ratio(1) * q_eff_tot(3, 1) + q_eff_tot(2, 1)) * cross_section_radius**2 / rc**2
      ds = (dif_csca_ratio(1) * q_eff_tot(3, 1) + q_eff_tot(2, 1)) * cross_section_radius**2 * pi / tvol
!         a=ds/2.d0
!         dela=1.d0
!         do while(abs(dela).gt.1.d-10)
!            dela=(3.*a*(1. + 2.*a*rc)-a*exp(2.*a*rc)*(3. + 2.*a**2*rc**2*(-3. + 2.*ds*rc))) &
!               /(6. - 6.*exp(2*a*rc) + 12.*a*rc*(1. + a*rc))
!            a=a+dela
!         enddo
!         extrat=2.d0*a/((mean_qext_mie)*pi*area_mean_radius**2*dble(number_spheres)/tvol)
! patch 12/2023
!
      a = ds * rc
      dela = 1.d0
      delextrat = 1.d0
      extrat0 = 0.d0
      do while (abs(delextrat) .gt. 1.d-5)
         dela = -(a * (-1.d0 - 2.d0 * a + exp(2.d0 * a) * (1.d0 + 2.d0 * a * a * (-1.d0 + qe)))) &
                / (2.d0 + 4.d0 * a * (1.d0 + a) - 2.d0 * exp(2.d0 * a))
         a = a + dela
         extrat = a / ((mean_qext_mie) * pi * area_mean_radius**2 * dble(number_spheres) / tvol) / rc
         delextrat = extrat - extrat0
         extrat0 = extrat
      end do
      a = a / 2.d0 / rc
   end subroutine diffuse_scattering_effective_ref_index
end module effective_medium_analysis
