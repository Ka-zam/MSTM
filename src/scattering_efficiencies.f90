module scattering_efficiencies
   use, intrinsic :: iso_fortran_env, only: real64
   use constants
   use fft_translation
   use mie
   use parallel_runtime, only: mpi_comm_world, mstm_global_rank, parallel_allreduce_sum, parallel_rank, &
                               parallel_receive, parallel_send, parallel_size
   use numerical_tables
   use periodic_lattice_operations
   use bessel_functions, only: riccati_bessel, riccati_hankel
   use quadrature, only: integrate_gauss_kronrod_adaptive
   use wave_functions, only: left_right_mode_transformation
   use sphere_data
   use surface
   use scattering_amplitudes, only: common_origin_amplitude_matrix, multiple_origin_amplitude_matrix, &
                                    numerical_scattering_matrix_azimuthal_average_multiple_origin, &
                                    numerical_scattering_matrix_azimuthal_average_single_origin
   use scattering_interactions, only: sphere_interaction
   implicit none
   private
   public :: boundary_extinction, boundary_scattering, common_origin_hemispherical_scattering, &
             configuration_efficiency_factors, extinction_theorem, hemispherical_integrand, &
             hemispherical_scattering, sphere_efficiency_factors, total_efficiency_factors, &
             waveguide_mode_scattering
   logical :: hemispherical_single_origin
   complex(real64), allocatable :: hemispherical_amnp(:)
contains

!
!  general efficiency factor calcuation
!  L/R formulation
!  April 2012
!  march 2013: something is always changing in this one
!
   subroutine sphere_efficiency_factors(ri0, nodr, npol, xsp, anp, gnpinc, gnp, qe, qs, qa)
      implicit none
      integer :: nodr, npol, m, n, p, p1, p2, k, ma, na, s, t, ss, st
      real(real64) :: qe(2 * npol - 1), qa(2 * npol - 1), qs(2 * npol - 1), const, xsp, qi(2 * npol - 1)
      complex(real64) :: anp(0:nodr + 1, nodr, 2, npol), gnp(0:nodr + 1, nodr, 2, npol), &
                         gnpinc(0:nodr + 1, nodr, 2, npol), &
                         psi(0:nodr, 2), xi(0:nodr, 2), psip(2), xip(2), xri(2), &
                         rib, aamat(2, 2), ggmat(2, 2), agmat(2, 2), &
                         gamat(2, 2), cterm, ri0(2)
      qe = 0.d0
      qa = 0.d0
      qs = 0.d0
      qi = 0.d0
      xri = xsp * ri0
      rib = 2.d0 / (1.d0 / ri0(1) + 1.d0 / ri0(2))
      do p = 1, 2
         call riccati_bessel(nodr, xri(p), psi(0, p))
         call riccati_hankel(nodr, xri(p), xi(0, p))
      end do
      do n = 1, nodr
         do s = 1, 2
            psip(s) = psi(n - 1, s) - n * psi(n, s) / xri(s)
            xip(s) = xi(n - 1, s) - n * xi(n, s) / xri(s)
         end do
         do s = 1, 2
            ss = (-1)**s
            do t = 1, 2
               st = (-1)**t
               cterm = cmplx(0.d0, 1.d0, kind=real64) * conjg(1./ri0(t)) / ri0(s)
               aamat(s, t) = cterm * (xip(s) * conjg(xi(n, t)) &
                                      - ss * st * xi(n, s) * conjg(xip(t)))
               ggmat(s, t) = cterm * (psip(s) * conjg(psi(n, t)) &
                                      - ss * st * psi(n, s) * conjg(psip(t)))
               agmat(s, t) = cterm * (xip(s) * conjg(psi(n, t)) &
                                      - ss * st * xi(n, s) * conjg(psip(t)))
               gamat(s, t) = cterm * (psip(s) * conjg(xi(n, t)) &
                                      - ss * st * psi(n, s) * conjg(xip(t)))
            end do
         end do
         do m = -n, n
            if (m .le. -1) then
               ma = n + 1
               na = -m
            else
               ma = m
               na = n
            end if
            do p1 = 1, npol
               do p2 = 1, npol
                  if (p1 .eq. 1 .and. p2 .eq. 1) then
                     k = 1
                     const = 1.d0
                  elseif (p1 .eq. 2 .and. p2 .eq. 2) then
                     k = 2
                     const = 1.d0
                  else
                     k = 3
                     const = .5d0
                  end if
                  do s = 1, 2
!if(k.eq.1) write(*,'(3i3,4es12.4)') m,n,s,anp(ma,na,s,p1),gnpinc(ma,na,s,p2)
                     do t = 1, 2
!                           qi(k)=qi(k)+const*ggmat(s,t)*(gnp(ma,na,s,p1)-gnpinc(ma,na,s,p1)) &
!                             *conjg(rib*(gnp(ma,na,t,p2)-gnpinc(ma,na,t,p2)))
!                           qi(k)=qi(k)+const*agmat(s,t)*anp(ma,na,s,p1)*conjg(gnp(ma,na,t,p2)-gnpinc(ma,na,t,p2)) &
!                               *(conjg(rib)+(-1)**(s+t)*rib)
                        qi(k) = qi(k) + const * ggmat(s, t) * (gnpinc(ma, na, s, p1)) &
                                * conjg(rib * (gnpinc(ma, na, t, p2)))
                        qa(k) = qa(k) + const &
                                * (aamat(s, t) * anp(ma, na, s, p1) * conjg(rib * anp(ma, na, t, p2)) &
                                   + ggmat(s, t) * gnp(ma, na, s, p1) * conjg(rib * gnp(ma, na, t, p2)) &
                                   + agmat(s, t) * anp(ma, na, s, p1) * conjg(gnp(ma, na, t, p2)) &
                                   * (conjg(rib) + (-1)**(s + t) * rib))
                        qs(k) = qs(k) - const &
                                * aamat(s, t) * anp(ma, na, s, p1) * conjg(rib * anp(ma, na, t, p2))
                        qe(k) = qe(k) + const &
                                * (agmat(s, t) * anp(ma, na, s, p1) * conjg(gnpinc(ma, na, t, p2)) &
                                   * (conjg(rib) + (-1)**(s + t) * rib))
!if(k.eq.1) write(*,'(2i3,4es12.4)') s,t,agmat(s,t),(conjg(rib)+(-1)**(s+t)*rib)
                     end do
                  end do
               end do
            end do
         end do
      end do
      do k = 1, 2 * npol - 1
         qi(k) = qi(k) * 2./xsp / xsp
         qe(k) = qe(k) * 2./xsp / xsp
         qa(k) = qa(k) * 2./xsp / xsp
         qs(k) = (qs(k)) * 2./xsp / xsp
! 10-22: qs replaced with qi.
         qs(k) = qi(k)
!
! qi is the absorption of the incident field within the sphere volume, not zero only for
! host sphere = 0 and dissipative external.   corrects qe so that qe=qs+qa
!
!            qe(k)=qe(k)+qi(k)
!write(*,'(i5,4e13.5)') k,qe(k)+qi(k),qa(k)+qs(k)
      end do

   end subroutine sphere_efficiency_factors
!
! calling routine for efficiency calculation
! april 2012: lr formulation
! february 2013:  number of rhs, mpi comm options added.
!
   subroutine configuration_efficiency_factors(nsphere, npol, amnp, gmnp0, qeff, &
                                               mpi_comm)
      implicit none
      integer :: nsphere, i, nblk, noff, neqns, nodr, &
                 npol, mpicomm, p, &
                 b11, b12, rank, numprocs
      integer, optional :: mpi_comm
      real(real64) :: qeff(3, 2 * npol - 1, nsphere), qe(2 * npol - 1), qa(2 * npol - 1), qs(2 * npol - 1)
      real(real64), allocatable :: qeffi(:, :)
      complex(real64) :: amnp(number_eqns, npol), gmnp0(number_eqns, npol), ri0(2)
      complex(real64), allocatable :: amnpi(:, :, :, :), gmnpi(:, :, :, :), fmnpi(:, :, :, :), &
                                      gmnp(:, :)
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      neqns = number_eqns
      allocate (gmnp(neqns, npol), qeffi(3, 2 * npol - 1))
      gmnp = 0.d0
      if (light_up) then
         write (*, '('' qe1 '',i3)') mstm_global_rank
         flush (6)
      end if
      do p = 1, npol
         if (number_host_spheres .eq. 0) then
            noff = 0
            do i = 1, number_spheres
               nodr = sphere_order(i)
               nblk = 2 * nodr * (nodr + 2)
               b11 = noff + 1
               b12 = b11 + nblk - 1
               call apply_single_sphere_mie_coefficients(i, nodr, amnp(b11:b12, p), gmnp(b11:b12, p), 'i')
               noff = noff + nblk * number_field_expansions(i)
            end do
         else
            call sphere_interaction(neqns, 1, amnp(:, p), gmnp(:, p), &
                                    mpi_comm=mpicomm, mie_mult=(/.false./), &
                                    store_matrix_option=.false.)
            gmnp(:, p) = gmnp(:, p) + gmnp0(:, p)
         end if
      end do
      if (light_up) then
         write (*, '('' qe2 '',i3)') mstm_global_rank
         flush (6)
      end if
      qeff(1:3, 1:2 * npol - 1, 1:nsphere) = 0.d0
      noff = 0
      do i = 1, number_spheres
         nodr = sphere_order(i)
         nblk = 2 * nodr * (nodr + 2)
         allocate (amnpi(0:nodr + 1, nodr, 2, npol), &
                   gmnpi(0:nodr + 1, nodr, 2, npol), &
                   fmnpi(0:nodr + 1, nodr, 2, npol))
         call exterior_refractive_index(i, ri0)
!            ri0=sphere_ref_index(:,host_sphere(i))
         b11 = noff + 1
         b12 = b11 + nblk - 1
         do p = 1, npol
            amnpi(0:nodr + 1, 1:nodr, 1:2, p) = reshape(amnp(b11:b12, p), &
                                                        (/nodr + 2, nodr, 2/))
            fmnpi(0:nodr + 1, 1:nodr, 1:2, p) = reshape(gmnp0(b11:b12, p), &
                                                        (/nodr + 2, nodr, 2/))
            gmnpi(0:nodr + 1, 1:nodr, 1:2, p) = reshape(gmnp(b11:b12, p), &
                                                        (/nodr + 2, nodr, 2/))
         end do
         call sphere_efficiency_factors(ri0, nodr, npol, sphere_radius(i), amnpi, fmnpi, gmnpi, &
                                        qe, qs, qa)
         qeffi(1, :) = qe(:)
         qeffi(3, :) = qs(:)
         qeffi(2, :) = qa(:)
         noff = noff + nblk * number_field_expansions(i)
         deallocate (gmnpi, amnpi, fmnpi)
         if (npol .eq. 1) then
            qeff(:, 1, i) = qeffi(:, 1)
         else
            qeff(:, 1, i) = .5 * (qeffi(:, 1) + qeffi(:, 2))
            qeff(:, 2, i) = qeffi(:, 1)
            qeff(:, 3, i) = qeffi(:, 2)
         end if
      end do
      deallocate (gmnp, qeffi)
      if (light_up) then
         write (*, '('' qe3 '',i3)') mstm_global_rank
         flush (6)
      end if
   end subroutine configuration_efficiency_factors

   subroutine total_efficiency_factors(nsphere, nrow, xgeff, qeffp, qabsvol, qefftot)
      implicit none
      integer :: nsphere, nrow, i, j
      real(real64) :: qeffp(3, nrow, nsphere), qefftot(3, nrow), qabsvol(nrow, nsphere), xgeff
      do i = 1, nsphere
         qabsvol(:, i) = qeffp(2, :, i)
         do j = 1, nsphere
            if (host_sphere(j) .eq. i) then
               qabsvol(:, i) = qabsvol(:, i) - qeffp(2, :, j) &
                               * sphere_radius(j)**2 &
                               / sphere_radius(i)**2
            end if
         end do
      end do
!         do i=1,nsphere
!            if(aimag(sphere_ref_index(1,i)).eq.0.d0 &
!             .and.aimag(sphere_ref_index(2,i)).eq.0.d0) then
!               qabsvol(:,i)=0.d0
!            endif
!         enddo
      qeffp(2, :, :) = 0.d0
      do i = 1, nsphere
         qeffp(2, :, i) = qabsvol(:, i)
         do j = 1, nsphere
            if (host_sphere(j) .eq. i) then
               qeffp(2, :, i) = qeffp(2, :, i) + qabsvol(:, j) &
                                * sphere_radius(j)**2 &
                                / sphere_radius(i)**2
            end if
         end do
      end do
      qefftot = 0.
      do i = 1, nsphere
         if (host_sphere(i) .eq. 0) then
            qefftot(:, :) = qefftot(:, :) + qeffp(:, :, i) * sphere_radius(i)**2 / xgeff**2
         end if
      end do
! 10--22 : qsca=qext + qinc-qabs
      qefftot(3, :) = qefftot(1, :) + qefftot(3, :) - qefftot(2, :)

   end subroutine total_efficiency_factors

   subroutine waveguide_mode_scattering(amnp, qsevan, mpi_comm)
      implicit none
      integer :: i, p, pole, mpicomm
      integer, optional :: mpi_comm
      real(real64) :: qsevan(2)
      complex(real64) :: amnp(number_eqns, 2)
      complex(real64), allocatable :: amn(:), gmn(:)
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      allocate (amn(number_eqns), gmn(number_eqns))
      pole_integration = .true.
      qsevan = 0.
      do pole = 1, number_singular_points
         pole_integration_s = singular_points(pole)
         recalculate_surface_matrix = .true.
         do p = 1, 2
            amn(:) = amnp(:, p)
            gmn = 0.d0
            call sphere_interaction(number_eqns, 1, amn, gmn, &
                                    mie_mult=(/.false./), initial_run=.true., skip_external_translation=.true., &
                                    mpi_comm=mpicomm)
            do i = 1, number_spheres
               if (host_sphere(i) .ne. 0) cycle
               call left_right_mode_transformation(sphere_order(i), &
                                                   amn(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i)), &
                                                   amn(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i)))
               call left_right_mode_transformation(sphere_order(i), &
                                                   gmn(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i)), &
                                                   gmn(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i)))
               qsevan(p) = qsevan(p) + sum(dble(conjg(amn(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i))) &
                                                * gmn(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i)))) &
                           / layer_ref_index(sphere_layer(i))
            end do
         end do
      end do
      qsevan = qsevan * 2.d0 / cross_section_radius**2
      deallocate (gmn, amn)
      pole_integration = .false.
   end subroutine waveguide_mode_scattering

   subroutine boundary_scattering(amn, qbsca, mpi_comm)
      implicit none
      integer :: i, j, rank, numprocs, mpicomm, task, nsend, k
      integer, optional :: mpi_comm
      real(real64) :: qbsca(2, 0:1), qt(2), targetz, qt2(2, 0:1)
      complex(real64) :: amn(number_eqns, 2)
      complex(real64), allocatable :: amn1(:, :), amn2(:, :)
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call parallel_size(mpi_size=numprocs, &
                         mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, &
                         mpi_comm=mpicomm)

      qbsca = 0.d0
      task = 0
      do k = 0, 1
         if (k .eq. 0) then
            targetz = bot_boundary
         elseif (k .eq. 1) then
            targetz = top_boundary
         end if
         if (k .eq. 0 .and. aimag(layer_ref_index(0)) .gt. 1.d-6) then
            qbsca(:, k) = 0.d0
         elseif (k .eq. 1 .and. &
                 aimag(layer_ref_index(number_plane_boundaries)) .gt. 1.d-6) then
            qbsca(:, k) = 0.d0
         else
            do i = 1, number_spheres
               if (host_sphere(i) .ne. 0) cycle
               allocate (amn1(sphere_block(i), 2))
               amn1 = amn(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i), 1:2)
               do j = 1, number_spheres
                  if (host_sphere(j) .ne. 0) cycle
                  task = task + 1
                  if (mod(task, numprocs) .ne. rank) cycle
                  allocate (amn2(sphere_block(j), 2))
                  amn2 = amn(sphere_offset(j) + 1:sphere_offset(j) + sphere_block(j), 1:2)
                  call sphere_boundary_scattering(sphere_order(i), sphere_position(:, i), &
                                                  amn1, sphere_order(j), sphere_position(:, j), amn2, targetz, qt)
                  qbsca(:, k) = qbsca(:, k) + qt
                  deallocate (amn2)
               end do
               deallocate (amn1)
            end do
         end if
      end do
      qbsca = qbsca / cross_section_radius**2
      if (numprocs .gt. 1) then
         nsend = 4
         qt2 = qbsca
         call parallel_allreduce_sum(mpi_number=nsend, &
                                     send_buffer=qt2, receive_buffer=qbsca, &
                                     mpi_comm=mpicomm)
      end if
   end subroutine boundary_scattering

   subroutine boundary_extinction(amn, alpha, sinc, dir, qbext, common_origin)
      implicit none
      logical :: comorg
      logical, optional :: common_origin
      integer :: k, dir
      real(real64) :: qt(2), targetz, qbext(2, 0:1), alpha, sinc
      complex(real64) :: amn(*)
      if (present(common_origin)) then
         comorg = common_origin
      else
         comorg = .false.
      end if
      qbext = 0.d0
      do k = 0, 1
         if (k .eq. 0) then
            targetz = bot_boundary
         elseif (k .eq. 1) then
            targetz = top_boundary
         end if
         if (k .eq. 0 .and. aimag(layer_ref_index(0)) .gt. 1.d-6 .and. dir .eq. 2) then
            qbext(:, k) = 0.d0
         elseif (k .eq. 1 .and. &
                 aimag(layer_ref_index(number_plane_boundaries)) .gt. 1.d-6 .and. dir .eq. 1) then
            qbext(:, k) = 0.d0
         else
            call extinction_theorem(amn, sinc, dir, alpha, targetz, qt, common_origin=comorg)
            qbext(:, k) = qt
         end if
      end do
   end subroutine boundary_extinction

   subroutine common_origin_hemispherical_scattering(amn, qbsca)
      implicit none
      integer :: k
      real(real64) :: qbsca(2, 0:1), targetz
      complex(real64) :: amn(2 * t_matrix_order * (t_matrix_order + 2), 2)

      qbsca = 0.d0
      do k = 0, 1
         if (k .eq. 0) then
!               targetz=min(bot_boundary,cluster_origin(3)-1.d2)
            targetz = bot_boundary
         elseif (k .eq. 1) then
            !              targetz=max(top_boundary,cluster_origin(3)+1.d2)
            targetz = top_boundary
         end if
         call sphere_boundary_scattering(t_matrix_order, (/0.d0, 0.d0, 0.d0/), &
                                         amn, t_matrix_order, (/0.d0, 0.d0, 0.d0/), amn, targetz, qbsca(:, k), lr_to_mode=.false.)
      end do
      qbsca = qbsca / cross_section_radius**2
   end subroutine common_origin_hemispherical_scattering

   subroutine hemispherical_scattering(amnp, singleorigin, numerical, qbsca, mpi_comm)
      implicit none
      logical :: singleorigin, numerical
      integer :: mpicomm, rank, numprocs, maxnumdiv, subdiv, errorcodes, p, p1, p2, ncoef
      integer, optional :: mpi_comm
      real(real64) :: qbsca(2, 2), inteps, mindiv, c0, c1
      complex(real64) :: amnp(*), cq(2, 2)
      data maxnumdiv, inteps, mindiv, c0, c1/3, 1.d-6, 1.d-4, 0.d0, 1.d0/
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call parallel_size(mpi_size=numprocs, &
                         mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, &
                         mpi_comm=mpicomm)

      if (.not. numerical) then
         if (singleorigin) then
            call common_origin_hemispherical_scattering(amnp, qbsca)
         else
            call boundary_scattering(amnp, qbsca, mpi_comm=mpicomm)
         end if
      else
         hemispherical_single_origin = singleorigin
         if (singleorigin) then
            ncoef = 4 * t_matrix_order * (t_matrix_order + 2)
         else
            ncoef = 2 * number_eqns
         end if
         if (allocated(hemispherical_amnp)) deallocate (hemispherical_amnp)
         allocate (hemispherical_amnp(ncoef))
         hemispherical_amnp = amnp(1:ncoef)
         if (numprocs .eq. 1) then
            p1 = 1
            p2 = 2
         else
            if (rank .eq. 0) then
               p1 = 1
               p2 = 1
            elseif (rank .eq. 1) then
               p1 = 2
               p2 = 2
            else
               p1 = 0
               p2 = 0
            end if
         end if
         cq = 0.d0
         do p = p1, p2
            if (p .eq. 0) cycle
            c0 = -2.d0 + dble(p)
            c1 = c0 + 1.d0
            call integrate_gauss_kronrod_adaptive(2, c0, c1, hemispherical_integrand, cq(:, p), subdiv, &
                                                  errorcodes, inteps, mindiv, maxnumdiv)
         end do
         if (numprocs .gt. 1) then
            if (rank .eq. 1) then
               call parallel_send(mpi_number=2, &
                                  send_buffer=cq(:, 2), mpi_rank=0, &
                                  mpi_comm=mpicomm)
            elseif (rank .eq. 0) then
               call parallel_receive(mpi_number=2, &
                                     receive_buffer=cq(:, 2), mpi_rank=1, &
                                     mpi_comm=mpicomm)
            end if
         end if
         qbsca = cq * 4.d0 / cross_section_radius**2.
         qbsca(:, 1) = -qbsca(:, 1)
         deallocate (hemispherical_amnp)
      end if
   end subroutine hemispherical_scattering

   subroutine hemispherical_integrand(ntot, ct, q)
      implicit none
      integer :: ntot
      real(real64) :: ct, sm(2)
      complex(real64) :: q(ntot)
      q = 0.d0
      if (hemispherical_single_origin) then
         call numerical_scattering_matrix_azimuthal_average_single_origin( &
            hemispherical_amnp, t_matrix_order, ct, sm, &
            rotate_plane=.false., normalize_s11=.false., s11_only=.true.)
         q(1:2) = sm(1:2)
      else
         call numerical_scattering_matrix_azimuthal_average_multiple_origin(hemispherical_amnp, ct, sm, &
                                                                            rotate_plane=.false., s11_only=.true.)
         q(1:2) = sm(1:2)
      end if
   end subroutine hemispherical_integrand

   subroutine extinction_theorem(amnp, sinc, sdir, alpha, targetz, qe, common_origin)
      implicit none
      logical :: comorg
      logical, optional :: common_origin
      integer :: sdir, p, tlay
      real(real64) :: sinc, alpha, qe(2), targetz, sourcez, const3(3, 2)
      complex(real64) :: amnp(*), s, sourceri, targetri, gfs(2, 2, 2), skz, tkz, sa(4), amp(2, 2)
      if (present(common_origin)) then
         comorg = common_origin
      else
         comorg = .false.
      end if
      if (sdir .eq. 1) then
         sourceri = layer_ref_index(0)
      else
         sourceri = layer_ref_index(number_plane_boundaries)
      end if
      sourcez = incident_field_boundary
      tlay = find_layer_index(targetz)
      targetri = layer_ref_index(tlay)
      s = sinc
      call layer_green_function(s, sourcez, targetz, gfs, skz, tkz, include_direct=.true.)
      do p = 1, 2
         if (comorg) then
            call common_origin_amplitude_matrix(amnp, s, alpha, targetz, p, sa)
         else
            call multiple_origin_amplitude_matrix(amnp, s, alpha, targetz, p, sa)
         end if
         amp(p, 1) = sa(2)
         amp(p, 2) = sa(1)
      end do
      const3(1, 1) = 16.d0 * (abs(tkz)**2 + s**2 / abs(targetri)**2) * dble(targetri * tkz)
      const3(2, 1) = -const3(1, 1)
      const3(3, 1) = 16.d0 * (abs(tkz)**2 - s**2 / abs(targetri)**2) * aimag(targetri * tkz)
      const3(1, 2) = 16.d0 * dble(targetri * tkz)
      const3(2, 2) = -const3(1, 2)
      const3(3, 2) = -16.d0 * aimag(targetri * tkz)
      do p = 1, 2
        qe(p) = (const3(1, p) * dble(conjg(gfs(1, sdir, p)) * amp(1, p)) + const3(2, p) * dble(conjg(gfs(2, sdir, p)) * amp(2, p)) &
                  + const3(3, p) * aimag(conjg(gfs(2, sdir, p)) * amp(1, p) - conjg(gfs(1, sdir, p)) * amp(2, p))) &
                 / incident_field_scale(p)
      end do
      qe = qe / cross_section_radius**2
   end subroutine extinction_theorem

end module scattering_efficiencies
