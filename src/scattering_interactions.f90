module scattering_interactions
   use, intrinsic :: iso_fortran_env, only: real64
   use constants
   use fft_translation, only: fft_plan
   use mie
   use parallel_runtime, only: mpi_comm_world, mstm_global_rank, parallel_allreduce_sum, &
                               parallel_allreduce_sum_complex64_sequence, parallel_barrier, parallel_rank, parallel_size
   use numerical_tables
   use periodic_lattice_operations
   use angular_functions, only: estimate_translation_order, generate_gaussian_beam_coefficients
   use coefficient_indexing, only: mode_index
   use sphere_data
   use surface
   use translation_expansions, only: external_to_external_expansion, external_to_internal_expansion
   use translation_operator, only: translation_operator_state
   use translation_surface_interactions, only: periodic_lattice_sphere_interaction, sphere_surface_interaction
   implicit none
   private
   public :: distribute_from_common_origin, estimate_sphere_translation_orders, &
             layered_gaussian_beam_coefficients, merge_to_common_origin, phase_shift, &
             sphere_interaction, sphere_plane_wave_coefficients
contains

   subroutine sphere_interaction(neqns, nrhs, ain, aout, initial_run, &
                                 rhs_list, mpi_comm, con_tran, mie_mult, fft_option, &
                                 store_matrix_option, skip_external_translation)
      implicit none
      integer :: neqns, rank, numprocs, nrhs, nsend, &
                 mpicomm, rhs
      logical :: initrun, rhslist(nrhs), contran(nrhs), miemult(nrhs), &
                 fftopt, smopt, etopt
      logical, optional :: initial_run, rhs_list(nrhs), con_tran(nrhs), &
                           mie_mult(nrhs), fft_option, store_matrix_option, skip_external_translation
      integer, optional :: mpi_comm
      complex(real64) :: aout_t(neqns, nrhs), ain_t(neqns, nrhs), &
                         ain(neqns, nrhs), aout(neqns, nrhs), &
                         aout_t2(neqns, nrhs)
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      if (present(initial_run)) then
         initrun = initial_run
      else
         initrun = .false.
      end if
      if (present(rhs_list)) then
         rhslist = rhs_list
      else
         rhslist = .true.
      end if
      if (present(con_tran)) then
         contran = con_tran
      else
         contran = .false.
      end if
      if (present(mie_mult)) then
         miemult = mie_mult
      else
         miemult = .true.
      end if
      if (present(fft_option)) then
         fftopt = fft_option
      else
         fftopt = sphere_cluster%fft_translation_option
      end if
      if (present(store_matrix_option)) then
         smopt = store_matrix_option
      else
         smopt = .true.
      end if
      if (present(skip_external_translation)) then
         etopt = .not. skip_external_translation
      else
         etopt = .true.
      end if
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
!
! sphere-to-sphere H (i-j) interaction
!
      if (light_up) then
         write (*, '('' si1 '',i3)') mstm_global_rank
         flush (6)
      end if
      do rhs = 1, nrhs
         if (contran(rhs)) then
            ain_t(:, rhs) = conjg(ain(:, rhs))
            if (miemult(rhs)) then
               call apply_mie_coefficients(neqns, 1, -1, ain_t(:, rhs), ain_t(:, rhs))
            end if
         else
            ain_t(:, rhs) = ain(:, rhs)
         end if
      end do

      call parallel_barrier(mpi_comm=mpicomm)
      aout_t = 0.
      if (etopt) then
         if (periodic_lattice) then
!write(*,*) 'step 3'
!flush(6)
            if (light_up) then
               write (*, '('' si2a '',i3)') mstm_global_rank
               flush (6)
            end if
!write(*,'('' rank,comm,init:'',3i12,l1)') numprocs,rank,mstm_global_rank,initrun
!flush(6)
!stop
            call periodic_lattice_sphere_interaction(neqns, nrhs, ain_t, aout_t, &
                                                     store_matrix_option=smopt, initial_run=initrun, &
                                                     rhs_list=rhslist, mpi_comm=mpicomm, con_tran=contran)
!            nsend=neqns*nrhs
!            call parallel_allreduce_sum(receive_buffer=aout_t, &
!                 mpi_number=nsend,mpi_comm=mpicomm)

         else
            if (light_up) then
               write (*, '('' si2b '',i3)') mstm_global_rank
               flush (6)
            end if
            if (fftopt) then
               call fft_plan%apply(neqns, nrhs, ain_t, aout_t, &
                                   store_matrix_option=smopt, initial_run=initrun, &
                                   rhs_list=rhslist, mpi_comm=mpicomm, con_tran=contran)
            else
               call external_to_external_expansion(neqns, nrhs, ain_t, aout_t, &
                                                   store_matrix_option=smopt, initial_run=initrun, &
                                                   rhs_list=rhslist, mpi_comm=mpicomm, con_tran=contran)
            end if
         end if
      end if

      call parallel_barrier(mpi_comm=mpicomm)
      if (sphere_cluster%number_host_spheres .gt. 0) then
         aout_t2 = 0.
         call external_to_internal_expansion(neqns, nrhs, ain_t, aout_t2, &
                                             rhs_list=rhslist, mpi_comm=mpicomm, con_tran=contran)
         aout_t = aout_t + aout_t2
      end if

      if (plane_surface_present .and. (.not. periodic_lattice)) then
         aout_t2 = 0.
         call sphere_surface_interaction(neqns, nrhs, ain_t, aout_t2, &
                                         initial_run=initrun, rhs_list=rhslist, &
                                         mpi_comm=mpicomm, con_tran=contran)
         aout_t = aout_t + aout_t2
      end if

      if (numprocs .gt. 1) then
         nsend = neqns * nrhs
         call parallel_allreduce_sum(receive_buffer=aout_t, &
                                     mpi_number=nsend, mpi_comm=mpicomm)
      end if

      do rhs = 1, nrhs
         if (.not. contran(rhs)) then
            if (miemult(rhs)) then
               call apply_mie_coefficients(neqns, 1, 1, aout_t(:, rhs), aout(:, rhs))
            else
               aout(:, rhs) = aout_t(:, rhs)
            end if
         else
            aout(:, rhs) = conjg(aout_t(:, rhs))
         end if
      end do

!         if(numprocs.gt.1) then
!            nsend=neqns*nrhs
!            call parallel_allreduce_sum(receive_buffer=aout, &
!                 mpi_number=nsend,mpi_comm=mpicomm)
!         endif
      if (light_up) then
         write (*, '('' si3 '',i3)') mstm_global_rank
         flush (6)
      end if

   end subroutine sphere_interaction

!
!  determination of maximum orders for target--based expansions
!
!
!  last revised: 15 January 2011
!
   subroutine estimate_sphere_translation_orders(eps, ntran, nodrt)
      implicit none
      integer :: nodrt, ntran(*), i, host, nodrmax
      real(real64) :: r, eps, rpos0(3), xi0(3)
      complex(real64) :: ri0, riext(2)
      nodrt = 0
      nodrmax = sphere_cluster%max_mie_order
      do i = 1, sphere_cluster%number_spheres
         host = sphere_cluster%host_sphere(i)
         call exterior_refractive_index(i, riext)
         ri0 = 2.d0 / (1.d0 / riext(1) + 1.d0 / riext(2))
         if (host .eq. 0) then
            rpos0 = sphere_cluster%cluster_origin
         else
            rpos0(:) = sphere_cluster%sphere_position(:, host)
         end if
         xi0(:) = sphere_cluster%sphere_position(:, i) - rpos0(:)
         r = sqrt(dot_product(xi0, xi0))
         call estimate_translation_order(r, ri0, sphere_cluster%sphere_order(i), eps, ntran(i))
         if (host .eq. 0) nodrt = max(nodrt, ntran(i), nodrmax)
      end do
   end subroutine estimate_sphere_translation_orders
!
!  plane wave expansion coefficients at sphere origins.  uses a phase shift.
!
!
!  last revised: 15 January 2011
!
! may 2019: k is now rightmost column
   subroutine sphere_plane_wave_coefficients(alpha, sinc, dir, pmnp, excited_spheres, mpi_comm)
      implicit none
      logical :: exsphere(sphere_cluster%number_spheres)
      logical, optional :: excited_spheres(sphere_cluster%number_spheres)
      integer :: p, i, dir, mpicomm, n, m, mn, rank
      integer, optional :: mpi_comm
      real(real64) :: alpha, sinc, qext, qsca, qabs
      complex(real64) :: pmnp(sphere_cluster%number_eqns, 2), ri1(2), ri0(2), pt(2, 2)
      complex(real64), allocatable :: pmnptot(:, :), dnpeff(:, :, :), pmnp0(:, :, :)
      if (present(excited_spheres)) then
         exsphere = excited_spheres
      else
         exsphere = .true.
      end if
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      if (sphere_cluster%effective_medium_simulation) then
         allocate (dnpeff(2, 2, sphere_cluster%t_matrix_order), pmnp0(sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2), 2, 2))
         ri0 = sphere_cluster%effective_ref_index
         ri1 = layer_ref_index(0)
         call optically_active_mie_coefficients(sphere_cluster%effective_cluster_radius, ri1, sphere_cluster%t_matrix_order, &
                                                0.d0, qext, qsca, qabs, &
                                                ri_medium=ri0, dnp_eff_mie=dnpeff)
!if(rank.eq.0) then
!write(*,'(3es16.9)') sphere_cluster%effective_cluster_radius,sphere_cluster%effective_ref_index
!do n=1,2
!write(*,'(i2,4es12.4)') n,dnpeff(1,1,n)+dnpeff(2,1,n),dnpeff(1,1,n)-dnpeff(2,1,n)
!enddo
!endif
         call layer_plane_wave_coefficients(alpha, sinc, dir, (/0.d0, 0.d0, 0.d0/), sphere_cluster%t_matrix_order, &
                                            pmnp0)
         do n = 1, sphere_cluster%t_matrix_order
            do m = -n, n
               mn = mode_index(m, n, sphere_cluster%t_matrix_order, 2)
               pt = pmnp0(mn, 1:2, 1:2)
               do p = 1, 2
                  pmnp0(mn, p, :) = dnpeff(p, 1, n) * pt(1, :) + dnpeff(p, 2, n) * pt(2, :)
               end do
            end do
         end do
         call distribute_from_common_origin(sphere_cluster%t_matrix_order, pmnp0, pmnp, number_rhs=2, &
                                            vswf_type=1, mpi_comm=mpicomm)
         deallocate (dnpeff, pmnp0)
         return
      end if

      pmnp = 0.d0
      do i = 1, sphere_cluster%number_spheres
         if (sphere_cluster%host_sphere(i) .ne. 0) cycle
         if (.not. exsphere(i)) cycle
         allocate (pmnptot(sphere_cluster%sphere_block(i), 2))
         if (sphere_cluster%gaussian_beam_constant .eq. 0.d0) then
        call layer_plane_wave_coefficients(alpha, sinc, dir, sphere_cluster%sphere_position(:, i), sphere_cluster%sphere_order(i), &
                                               pmnptot)
         else
   call layered_gaussian_beam_coefficients(alpha, sinc, dir, sphere_cluster%sphere_position(:, i), sphere_cluster%sphere_order(i), &
                                                    pmnptot)
         end if
         do p = 1, 2
            pmnp(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i), p) &
               = pmnptot(1:sphere_cluster%sphere_block(i), p)
         end do
         deallocate (pmnptot)
      end do
   end subroutine sphere_plane_wave_coefficients

   subroutine layered_gaussian_beam_coefficients(alpha, sinc, sdir, rpos, nodr, pmnp, include_direct, &
                                                 include_indirect)
      implicit none
      logical :: incdir, incindir, shiftgb
      logical, optional :: include_direct, include_indirect
      integer :: p, incregion, sdir, nodr, layer, nodrgb
      real(real64) :: alpha, rpos(3), sinc, cbinc, rtran(3), r, rshft(3), zs
      complex(real64) :: pmnp(2 * nodr * (nodr + 2), 2), riinc
      complex(real64), allocatable :: pmnp0(:, :), ptvec(:)
      type(translation_operator_state) :: tranmat
      if (sdir .eq. 1) then
         riinc = layer_ref_index(0)
         incregion = 0
      else
         riinc = layer_ref_index(number_plane_boundaries)
         incregion = number_plane_boundaries
      end if
      if (present(include_direct)) then
         incdir = include_direct
      else
         incdir = .true.
      end if
      if (present(include_indirect)) then
         incindir = include_indirect
      else
         incindir = .true.
      end if
      layer = find_layer_index(rpos(3))
      cbinc = dble(sqrt((1.d0 - sinc / riinc) * (1.d0 + sinc / riinc)) * (3 - 2 * sdir))
      pmnp = 0.d0
      rtran = rpos(:) - sphere_cluster%gaussian_beam_focal_point(:)
      r = sqrt(sum(rtran**2))
      call estimate_translation_order(r, riinc, nodr, 1.d-6, nodrgb)
      nodrgb = max(nodrgb, nodr)
!write(*,'('' gb order:'',i8)') nodrgb
      nodrgb = min(nodrgb, sphere_cluster%max_t_matrix_order)
      allocate (pmnp0(2 * nodrgb * (nodrgb + 2), 2))
      pmnp = 0.d0
!write(*,*) 's1'
      call generate_gaussian_beam_coefficients(alpha, cbinc, sphere_cluster%gaussian_beam_constant, nodrgb, pmnp0)
      if (layer .eq. incregion .and. incdir) then
         call tranmat%configure(1, rtran, (/riinc, riinc/), nodrgb .ge. sphere_cluster%translation_switch_order)
         do p = 1, 2
!write(*,*) 's2', p
            call tranmat%apply(nodrgb, 2, nodr, 2, pmnp0(:, p), pmnp(:, p))
         end do
         call tranmat%clear()
      end if
      if (number_plane_boundaries .gt. 0 .and. incindir) then
         shiftgb = .false.
         if (sdir .eq. 1 .and. sphere_cluster%gaussian_beam_focal_point(3) .ge. 0.d0) then
            shiftgb = .true.
            rshft = (/0.d0, 0.d0, -1.d-5 - sphere_cluster%gaussian_beam_focal_point(3)/)
            zs = -1.d-5
         elseif (sdir .ne. 1 .and. &
                 sphere_cluster%gaussian_beam_focal_point(3) .lt. plane_boundary_position(number_plane_boundaries)) then
            shiftgb = .true.
      rshft = (/0.d0, 0.d0, plane_boundary_position(number_plane_boundaries) + 1.d-5 - sphere_cluster%gaussian_beam_focal_point(3)/)
            zs = plane_boundary_position(number_plane_boundaries) + 1.d-5
         else
            zs = sphere_cluster%gaussian_beam_focal_point(3)
         end if
         if (shiftgb) then
            allocate (ptvec(2 * nodrgb * (nodrgb + 2)))
            call tranmat%configure(1, rshft, (/riinc, riinc/), nodrgb .ge. sphere_cluster%translation_switch_order)
            do p = 1, 2
               ptvec(:) = pmnp0(:, p)
               pmnp0(:, p) = 0.d0
               call tranmat%apply(nodrgb, 2, nodrgb, 2, ptvec(:), pmnp0(:, p))
            end do
            call tranmat%clear()
            deallocate (ptvec)
         end if

!write(*,*) 's3'
!            allocate(rmat(2*nodr*(nodr+2),2*nodrgb*(nodrgb+2)))
         allocate (ptvec(2 * nodr * (nodr + 2)))
!write(*,*) 's4',nodr,nodrgb

         do p = 1, 2
            ptvec = 0.d0
            call plane_boundary_interaction(nodr, nodrgb, &
                                            rtran(1), rtran(2), zs, rpos(3), &
                                            ptvec, index_model=2, lr_transformation=.true., &
                                            make_symmetric=.false., propagating_directions_only=.true., &
                                            source_vector=pmnp0(:, p))
!               pmnp(:,p)=pmnp(:,p)+0.5d0*matmul(rmat(:,:),pmnp0(:,p))
            pmnp(:, p) = pmnp(:, p) + 0.5d0 * ptvec(:)
         end do
!write(*,*) 's5'
         deallocate (ptvec)
      end if
      deallocate (pmnp0)
   end subroutine layered_gaussian_beam_coefficients

   subroutine phase_shift(amnp, dir)
      implicit none
      integer :: dir, i
      real(real64), save :: ilv(2)
      complex(real64) :: amnp(sphere_cluster%number_eqns, 2)

      if (dir .eq. 1) then
         do i = 1, sphere_cluster%number_spheres
            if (sphere_cluster%host_sphere(i) .ne. 0) cycle
            amnp(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i), 1:2) &
               = amnp(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i), 1:2) &
                 * exp((0.d0, -1.d0) * sum(incident_lateral_vector * sphere_cluster%sphere_position(1:2, i)))
         end do
         ilv = incident_lateral_vector
!            incident_lateral_vector=0.d0
      else
!            incident_lateral_vector=ilv
         do i = 1, sphere_cluster%number_spheres
            if (sphere_cluster%host_sphere(i) .ne. 0) cycle
            amnp(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i), 1:2) &
               = amnp(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i), 1:2) &
                 * exp((0.d0, 1.d0) * sum(incident_lateral_vector * sphere_cluster%sphere_position(1:2, i)))
         end do
      end if
   end subroutine phase_shift
!
!  translation of sphere-based expansions to common target origin
!
!
!  last revised: 15 January 2011
!  april 2012: l/r formulation: amnp is in l/r basis, amnp0 is in e/m basis
!  february 2013: number rhs, mpi comm options.   This is for general sphere configurations.
!
   subroutine merge_to_common_origin(nodrt, amnp, amnp0, number_rhs, &
                                     single_sphere, origin_position, mpi_comm, merge_procs, merge_radius, &
                                     sphere_translation_list)
      implicit none
      logical :: mergeprocs, tlist(sphere_cluster%number_spheres)
      logical, optional :: merge_procs, sphere_translation_list(sphere_cluster%number_spheres)
      integer :: nodrt, i, m, n, nblk, ntrani, nrhs, task, &
                 rank, numprocs, nsend, proc, mpicomm, noi, &
                 startsphere, endsphere, rhs
      integer, optional :: number_rhs, mpi_comm, single_sphere
      real(real64) :: r0(3), mrad, ri0
      real(real64), optional :: origin_position(3), merge_radius
      complex(real64) :: amnp(sphere_cluster%number_eqns, *), amnp0(0:nodrt + 1, nodrt, 2, *), rimedium(2)
      complex(real64), allocatable :: amnpt(:, :, :)
      type(translation_operator_state) :: tranmat

      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      if (present(origin_position)) then
         r0 = origin_position
      else
         r0 = 0.d0
      end if
      if (present(merge_radius)) then
         mrad = merge_radius
      else
         mrad = 1.d10
      end if
      if (present(number_rhs)) then
         nrhs = number_rhs
      else
         nrhs = 1
      end if
      if (present(merge_procs)) then
         mergeprocs = merge_procs
      else
         mergeprocs = .true.
      end if
      if (present(single_sphere)) then
         startsphere = single_sphere
         endsphere = single_sphere
      else
         startsphere = 1
         endsphere = sphere_cluster%number_spheres
      end if
      if (present(sphere_translation_list)) then
         tlist = sphere_translation_list
      else
         tlist = .true.
      end if
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      amnp0(:, :, :, 1:nrhs) = (0.d0, 0.d0)
      task = 0
      do i = startsphere, endsphere
         if (.not. tlist(i)) cycle
         nblk = sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2) * 2
         if (sphere_cluster%host_sphere(i) .eq. 0) then
            ri0 = sqrt(sum(r0 - sphere_cluster%sphere_position(:, i))**2)
            if (ri0 .gt. mrad) cycle
            task = task + 1
            proc = mod(task, numprocs)
            if (proc .eq. rank) then
               call exterior_refractive_index(i, rimedium)
               noi = sphere_cluster%sphere_order(i)
               ntrani = min(nodrt, sphere_cluster%translation_order(i))
!                  ntrani=max(ntrani,noi)
               allocate (amnpt(0:ntrani + 1, ntrani, 2))
               amnpt = (0.d0, 0.d0)
               call tranmat%configure(1, r0 - sphere_cluster%sphere_position(:, i), rimedium, &
                                      max(noi, ntrani) .ge. sphere_cluster%translation_switch_order)
               do rhs = 1, nrhs
                  call tranmat%apply(noi, 2, ntrani, 2, &
                                     amnp(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + nblk, rhs), amnpt)
                  do n = 1, ntrani
                     do m = 0, ntrani + 1
                        amnp0(m, n, 1, rhs) = amnp0(m, n, 1, rhs) &
                                              + amnpt(m, n, 1) + amnpt(m, n, 2)
                        amnp0(m, n, 2, rhs) = amnp0(m, n, 2, rhs) &
                                              + amnpt(m, n, 1) - amnpt(m, n, 2)
                     end do
                  end do
               end do
               call tranmat%clear()
               deallocate (amnpt)
            end if
         end if
      end do
      if (numprocs .gt. 1 .and. mergeprocs) then
         nsend = 2 * nodrt * (nodrt + 2) * nrhs
         call parallel_allreduce_sum_complex64_sequence(receive_buffer=amnp0, &
                                                        mpi_number=nsend, mpi_comm=mpicomm)
      end if
   end subroutine merge_to_common_origin

   subroutine distribute_from_common_origin(nodrt, amnp0, amnp, number_rhs, &
                                            origin_position, origin_host, vswf_type, mpi_comm, merge_procs, &
                                            single_sphere, sphere_translation_list)
      implicit none
      logical :: mergeprocs, tlist(sphere_cluster%number_spheres)
      logical, optional :: merge_procs, sphere_translation_list(sphere_cluster%number_spheres)
      integer :: nodrt, i, nblk, ntrani, nrhs, task, rhs, startsphere, endsphere, &
                 rank, numprocs, nsend, proc, mpicomm, noi, vtype, ohost
      integer, optional :: number_rhs, mpi_comm, vswf_type, origin_host, single_sphere
      real(real64) :: r0(3)
      real(real64), optional :: origin_position(3)
      complex(real64) :: amnp(sphere_cluster%number_eqns, *), amnp0(0:nodrt + 1, nodrt, 2, *), rimedium(2)
      type(translation_operator_state) :: tranmat

      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      if (present(origin_position)) then
         r0 = origin_position
      else
         r0 = 0.d0
      end if
      if (present(origin_host)) then
         ohost = origin_host
      else
         ohost = 0
      end if
      if (present(number_rhs)) then
         nrhs = number_rhs
      else
         nrhs = 1
      end if
      if (present(vswf_type)) then
         vtype = vswf_type
      else
         vtype = 1
      end if
      if (present(merge_procs)) then
         mergeprocs = merge_procs
      else
         mergeprocs = .true.
      end if
      if (present(single_sphere)) then
         if (single_sphere .ne. 0) then
            startsphere = single_sphere
            endsphere = single_sphere
         else
            startsphere = 1
            endsphere = sphere_cluster%number_spheres
         end if
      else
         startsphere = 1
         endsphere = sphere_cluster%number_spheres
      end if
      if (present(sphere_translation_list)) then
         tlist = sphere_translation_list
      else
         tlist = .true.
      end if

      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      amnp(1:sphere_cluster%number_eqns, 1:nrhs) = (0.d0, 0.d0)
      task = 0
      do i = startsphere, endsphere
         if (.not. tlist(i)) cycle
         nblk = sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2) * 2
         if (sphere_cluster%host_sphere(i) .eq. ohost) then
            task = task + 1
            proc = mod(task, numprocs)
            if (proc .eq. rank) then
               call exterior_refractive_index(i, rimedium)
               noi = sphere_cluster%sphere_order(i)
               ntrani = nodrt
               call tranmat%configure(vtype, sphere_cluster%sphere_position(:, i) - r0, rimedium, &
                                      max(noi, ntrani) .ge. sphere_cluster%translation_switch_order)
               do rhs = 1, nrhs
                  call tranmat%apply(ntrani, 2, noi, 2, &
                                     amnp0(0:nodrt + 1, 1:nodrt, 1:2, rhs), &
                                     amnp(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + nblk, rhs))
               end do
               call tranmat%clear()
            end if
         end if
      end do
      if (numprocs .gt. 1 .and. mergeprocs) then
         nsend = sphere_cluster%number_eqns * nrhs
         call parallel_allreduce_sum_complex64_sequence(receive_buffer=amnp, &
                                                        mpi_number=nsend, mpi_comm=mpicomm)
      end if
   end subroutine distribute_from_common_origin
end module scattering_interactions
