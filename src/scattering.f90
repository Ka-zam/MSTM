!
! scatprops module: various subroutines for calculation of observables from the solution
!
!
!  last revised: 15 January 2011
!
module scatprops
   use mpidefs
   use intrinsics
   use numconstants
   use specialfuncs
   use surface_subroutines
   use periodic_lattice_subroutines
   use spheredata
   use mie
   use translation
   use fft_translation
   implicit none
   logical :: hemispherical_single_origin
   complex(8), allocatable :: hemispherical_amnp(:)
contains

! The general sphere interaction driver
! LR formulation
! february 2013: number of rhs now a required argument list, and mpi comm is an option
!
   subroutine sphereinteraction(neqns, nrhs, ain, aout, initial_run, &
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
      complex(8) :: aout_t(neqns, nrhs), ain_t(neqns, nrhs), &
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
         fftopt = fft_translation_option
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
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
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
               call multmiecoeffmult(neqns, 1, -1, ain_t(:, rhs), ain_t(:, rhs))
            end if
         else
            ain_t(:, rhs) = ain(:, rhs)
         end if
      end do

      call mstm_mpi(mpi_command='barrier', mpi_comm=mpicomm)
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
!            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=aout_t, &
!                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)

         else
            if (light_up) then
               write (*, '('' si2b '',i3)') mstm_global_rank
               flush (6)
            end if
            if (fftopt) then
               call fft_external_to_external_expansion(neqns, nrhs, ain_t, aout_t, &
                                                       store_matrix_option=smopt, initial_run=initrun, &
                                                       rhs_list=rhslist, mpi_comm=mpicomm, con_tran=contran)
            else
               call external_to_external_expansion(neqns, nrhs, ain_t, aout_t, &
                                                   store_matrix_option=smopt, initial_run=initrun, &
                                                   rhs_list=rhslist, mpi_comm=mpicomm, con_tran=contran)
            end if
         end if
      end if

      call mstm_mpi(mpi_command='barrier', mpi_comm=mpicomm)
      if (number_host_spheres .gt. 0) then
         aout_t2 = 0.
         call external_to_internal_expansion(neqns, nrhs, ain_t, aout_t2, &
                                             rhs_list=rhslist, mpi_comm=mpicomm, con_tran=contran)
         aout_t = aout_t + aout_t2
      end if

      if (plane_surface_present .and. (.not. periodic_lattice)) then
         aout_t2 = 0.
         call spheresurfaceinteraction(neqns, nrhs, ain_t, aout_t2, &
                                       initial_run=initrun, rhs_list=rhslist, &
                                       mpi_comm=mpicomm, con_tran=contran)
         aout_t = aout_t + aout_t2
      end if

      if (numprocs .gt. 1) then
         nsend = neqns * nrhs
         call mstm_mpi(mpi_command='allreduce', mpi_recv_buf_dc=aout_t, &
                       mpi_number=nsend, mpi_operation=mstm_mpi_sum, mpi_comm=mpicomm)
      end if

      do rhs = 1, nrhs
         if (.not. contran(rhs)) then
            if (miemult(rhs)) then
               call multmiecoeffmult(neqns, 1, 1, aout_t(:, rhs), aout(:, rhs))
            else
               aout(:, rhs) = aout_t(:, rhs)
            end if
         else
            aout(:, rhs) = conjg(aout_t(:, rhs))
         end if
      end do

!         if(numprocs.gt.1) then
!            nsend=neqns*nrhs
!            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=aout, &
!                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
!         endif
      if (light_up) then
         write (*, '('' si3 '',i3)') mstm_global_rank
         flush (6)
      end if

   end subroutine sphereinteraction

!
!  determination of maximum orders for target--based expansions
!
!
!  last revised: 15 January 2011
!
   subroutine tranorders(eps, ntran, nodrt)
      implicit none
      integer :: nodrt, ntran(*), i, host, nodrmax
      real(8) :: r, eps, rpos0(3), xi0(3)
      complex(8) :: ri0, riext(2)
      nodrt = 0
      nodrmax = max_mie_order
      do i = 1, number_spheres
         host = host_sphere(i)
         call exteriorrefindex(i, riext)
         ri0 = 2.d0 / (1.d0 / riext(1) + 1.d0 / riext(2))
         if (host .eq. 0) then
            rpos0 = cluster_origin
         else
            rpos0(:) = sphere_position(:, host)
         end if
         xi0(:) = sphere_position(:, i) - rpos0(:)
         r = sqrt(dot_product(xi0, xi0))
         call tranordertest(r, ri0, sphere_order(i), eps, ntran(i))
         if (host .eq. 0) nodrt = max(nodrt, ntran(i), nodrmax)
      end do
   end subroutine tranorders
!
!  plane wave expansion coefficients at sphere origins.  uses a phase shift.
!
!
!  last revised: 15 January 2011
!
! may 2019: k is now rightmost column
   subroutine sphereplanewavecoef(alpha, sinc, dir, pmnp, excited_spheres, mpi_comm)
      implicit none
      logical :: exsphere(number_spheres)
      logical, optional :: excited_spheres(number_spheres)
      integer :: p, i, dir, mpicomm, n, m, mn, rank
      integer, optional :: mpi_comm
      real(8) :: alpha, sinc, qext, qsca, qabs
      complex(8) :: pmnp(number_eqns, 2), ri1(2), ri0(2), pt(2, 2)
      complex(8), allocatable :: pmnptot(:, :), dnpeff(:, :, :), pmnp0(:, :, :)
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
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      if (effective_medium_simulation) then
         allocate (dnpeff(2, 2, t_matrix_order), pmnp0(t_matrix_order * (t_matrix_order + 2), 2, 2))
         ri0 = effective_ref_index
         ri1 = layer_ref_index(0)
         call mieoa(effective_cluster_radius, ri1, t_matrix_order, 0.d0, qext, qsca, qabs, &
                    ri_medium=ri0, dnp_eff_mie=dnpeff)
!if(rank.eq.0) then
!write(*,'(3es16.9)') effective_cluster_radius,effective_ref_index
!do n=1,2
!write(*,'(i2,4es12.4)') n,dnpeff(1,1,n)+dnpeff(2,1,n),dnpeff(1,1,n)-dnpeff(2,1,n)
!enddo
!endif
         call layerplanewavecoef(alpha, sinc, dir, (/0.d0, 0.d0, 0.d0/), t_matrix_order, &
                                 pmnp0)
         do n = 1, t_matrix_order
            do m = -n, n
               mn = amnaddress(m, n, t_matrix_order, 2)
               pt = pmnp0(mn, 1:2, 1:2)
               do p = 1, 2
                  pmnp0(mn, p, :) = dnpeff(p, 1, n) * pt(1, :) + dnpeff(p, 2, n) * pt(2, :)
               end do
            end do
         end do
         call distribute_from_common_origin(t_matrix_order, pmnp0, pmnp, number_rhs=2, &
                                            vswf_type=1, mpi_comm=mpicomm)
         deallocate (dnpeff, pmnp0)
         return
      end if

      pmnp = 0.d0
      do i = 1, number_spheres
         if (host_sphere(i) .ne. 0) cycle
         if (.not. exsphere(i)) cycle
         allocate (pmnptot(sphere_block(i), 2))
         if (gaussian_beam_constant .eq. 0.d0) then
            call layerplanewavecoef(alpha, sinc, dir, sphere_position(:, i), sphere_order(i), &
                                    pmnptot)
         else
            call layergaussbeamcoef(alpha, sinc, dir, sphere_position(:, i), sphere_order(i), &
                                    pmnptot)
         end if
         do p = 1, 2
            pmnp(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i), p) &
               = pmnptot(1:sphere_block(i), p)
         end do
         deallocate (pmnptot)
      end do
   end subroutine sphereplanewavecoef

   subroutine layergaussbeamcoef(alpha, sinc, sdir, rpos, nodr, pmnp, include_direct, &
                                 include_indirect)
      implicit none
      logical :: incdir, incindir, shiftgb
      logical, optional :: include_direct, include_indirect
      integer :: p, incregion, sdir, nodr, layer, nodrgb
      real(8) :: alpha, rpos(3), sinc, cbinc, rtran(3), r, rshft(3), zs
      complex(8) :: pmnp(2 * nodr * (nodr + 2), 2), riinc
      complex(8), allocatable :: pmnp0(:, :), ptvec(:)
      type(translation_data), pointer :: loc_tranmat
      type(translation_data), target :: tranmat
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
      layer = layer_id(rpos(3))
      cbinc = dble(sqrt((1.d0 - sinc / riinc) * (1.d0 + sinc / riinc)) * (3 - 2 * sdir))
      pmnp = 0.d0
      rtran = rpos(:) - gaussian_beam_focal_point(:)
      r = sqrt(sum(rtran**2))
      call tranordertest(r, riinc, nodr, 1.d-6, nodrgb)
      nodrgb = max(nodrgb, nodr)
!write(*,'('' gb order:'',i8)') nodrgb
      nodrgb = min(nodrgb, max_t_matrix_order)
      allocate (pmnp0(2 * nodrgb * (nodrgb + 2), 2))
      pmnp = 0.d0
!write(*,*) 's1'
      call gaussianbeamcoef(alpha, cbinc, gaussian_beam_constant, nodrgb, pmnp0)
      if (layer .eq. incregion .and. incdir) then
         tranmat%matrix_calculated = .false.
         tranmat%vswf_type = 1
         tranmat%translation_vector = rtran
         tranmat%refractive_index = (/riinc, riinc/)
         tranmat%rot_op = (nodrgb .ge. translation_switch_order)
         loc_tranmat => tranmat
         do p = 1, 2
!write(*,*) 's2', p
            call coefficient_translation(nodrgb, 2, nodr, 2, &
                                         pmnp0(:, p), pmnp(:, p), loc_tranmat)
         end do
         if (.not. tranmat%zero_translation) then
            if (tranmat%rot_op) then
               if (associated(tranmat%rot_mat)) deallocate (tranmat%rot_mat, tranmat%phi_mat, tranmat%z_mat)
               nullify (tranmat%rot_mat, tranmat%phi_mat, tranmat%z_mat)
               nullify (loc_tranmat)
            else
               if (associated(tranmat%gen_mat)) deallocate (tranmat%gen_mat)
               nullify (tranmat%gen_mat)
               nullify (loc_tranmat)
            end if
         end if
      end if
      if (number_plane_boundaries .gt. 0 .and. incindir) then
         shiftgb = .false.
         if (sdir .eq. 1 .and. gaussian_beam_focal_point(3) .ge. 0.d0) then
            shiftgb = .true.
            rshft = (/0.d0, 0.d0, -1.d-5 - gaussian_beam_focal_point(3)/)
            zs = -1.d-5
         elseif (sdir .ne. 1 .and. &
                 gaussian_beam_focal_point(3) .lt. plane_boundary_position(number_plane_boundaries)) then
            shiftgb = .true.
            rshft = (/0.d0, 0.d0, plane_boundary_position(number_plane_boundaries) + 1.d-5 - gaussian_beam_focal_point(3)/)
            zs = plane_boundary_position(number_plane_boundaries) + 1.d-5
         else
            zs = gaussian_beam_focal_point(3)
         end if
         if (shiftgb) then
            allocate (ptvec(2 * nodrgb * (nodrgb + 2)))
            tranmat%matrix_calculated = .false.
            tranmat%vswf_type = 1
            tranmat%translation_vector = rshft
            tranmat%refractive_index = (/riinc, riinc/)
            tranmat%rot_op = (nodrgb .ge. translation_switch_order)
            loc_tranmat => tranmat
            do p = 1, 2
               ptvec(:) = pmnp0(:, p)
               pmnp0(:, p) = 0.d0
               call coefficient_translation(nodrgb, 2, nodrgb, 2, &
                                            ptvec(:), pmnp0(:, p), loc_tranmat)
            end do
            if (.not. tranmat%zero_translation) then
               if (tranmat%rot_op) then
                  if (associated(tranmat%rot_mat)) deallocate (tranmat%rot_mat, tranmat%phi_mat, tranmat%z_mat)
                  nullify (tranmat%rot_mat, tranmat%phi_mat, tranmat%z_mat)
                  nullify (loc_tranmat)
               else
                  if (associated(tranmat%gen_mat)) deallocate (tranmat%gen_mat)
                  nullify (tranmat%gen_mat)
                  nullify (loc_tranmat)
               end if
            end if
            deallocate (ptvec)
         end if

!write(*,*) 's3'
!            allocate(rmat(2*nodr*(nodr+2),2*nodrgb*(nodrgb+2)))
         allocate (ptvec(2 * nodr * (nodr + 2)))
!write(*,*) 's4',nodr,nodrgb

         do p = 1, 2
            ptvec = 0.d0
            call plane_interaction(nodr, nodrgb, &
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
   end subroutine layergaussbeamcoef

   subroutine phase_shift(amnp, dir)
      implicit none
      integer :: dir, i
      real(8), save :: ilv(2)
      complex(8) :: amnp(number_eqns, 2)

      if (dir .eq. 1) then
         do i = 1, number_spheres
            if (host_sphere(i) .ne. 0) cycle
            amnp(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i), 1:2) &
               = amnp(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i), 1:2) &
                 * exp((0.d0, -1.d0) * sum(incident_lateral_vector * sphere_position(1:2, i)))
         end do
         ilv = incident_lateral_vector
!            incident_lateral_vector=0.d0
      else
!            incident_lateral_vector=ilv
         do i = 1, number_spheres
            if (host_sphere(i) .ne. 0) cycle
            amnp(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i), 1:2) &
               = amnp(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i), 1:2) &
                 * exp((0.d0, 1.d0) * sum(incident_lateral_vector * sphere_position(1:2, i)))
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
      logical :: mergeprocs, tlist(number_spheres)
      logical, optional :: merge_procs, sphere_translation_list(number_spheres)
      integer :: nodrt, i, m, n, nblk, ntrani, nrhs, task, &
                 rank, numprocs, nsend, proc, mpicomm, noi, &
                 startsphere, endsphere, rhs
      integer, optional :: number_rhs, mpi_comm, single_sphere
      real(8) :: r0(3), mrad, ri0
      real(8), optional :: origin_position(3), merge_radius
      complex(8) :: amnp(number_eqns, *), amnp0(0:nodrt + 1, nodrt, 2, *), rimedium(2)
      complex(8), allocatable :: amnpt(:, :, :)
      type(translation_data), pointer :: loc_tranmat
      type(translation_data), target :: tranmat

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
         endsphere = number_spheres
      end if
      if (present(sphere_translation_list)) then
         tlist = sphere_translation_list
      else
         tlist = .true.
      end if
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      amnp0(:, :, :, 1:nrhs) = (0.d0, 0.d0)
      task = 0
      do i = startsphere, endsphere
         if (.not. tlist(i)) cycle
         nblk = sphere_order(i) * (sphere_order(i) + 2) * 2
         if (host_sphere(i) .eq. 0) then
            ri0 = sqrt(sum(r0 - sphere_position(:, i))**2)
            if (ri0 .gt. mrad) cycle
            task = task + 1
            proc = mod(task, numprocs)
            if (proc .eq. rank) then
               call exteriorrefindex(i, rimedium)
               noi = sphere_order(i)
               ntrani = min(nodrt, translation_order(i))
!                  ntrani=max(ntrani,noi)
               allocate (amnpt(0:ntrani + 1, ntrani, 2))
               amnpt = (0.d0, 0.d0)
               tranmat%matrix_calculated = .false.
               tranmat%vswf_type = 1
               tranmat%translation_vector = r0 - sphere_position(:, i)
               tranmat%refractive_index = rimedium
               tranmat%rot_op = (max(noi, ntrani) .ge. translation_switch_order)
               loc_tranmat => tranmat
               do rhs = 1, nrhs
                  call coefficient_translation(noi, 2, ntrani, 2, &
                                               amnp(sphere_offset(i) + 1:sphere_offset(i) + nblk, rhs), amnpt, &
                                               loc_tranmat)
                  do n = 1, ntrani
                     do m = 0, ntrani + 1
                        amnp0(m, n, 1, rhs) = amnp0(m, n, 1, rhs) &
                                              + amnpt(m, n, 1) + amnpt(m, n, 2)
                        amnp0(m, n, 2, rhs) = amnp0(m, n, 2, rhs) &
                                              + amnpt(m, n, 1) - amnpt(m, n, 2)
                     end do
                  end do
               end do
               if (.not. tranmat%zero_translation) then
                  if (tranmat%rot_op) then
                     if (associated(tranmat%rot_mat)) deallocate (tranmat%rot_mat, tranmat%phi_mat, tranmat%z_mat)
                     nullify (tranmat%rot_mat, tranmat%phi_mat, tranmat%z_mat)
                     nullify (loc_tranmat)
                  else
                     if (associated(tranmat%gen_mat)) deallocate (tranmat%gen_mat)
                     nullify (tranmat%gen_mat)
                     nullify (loc_tranmat)
                  end if
               end if
               deallocate (amnpt)
            end if
         end if
      end do
      if (numprocs .gt. 1 .and. mergeprocs) then
         nsend = 2 * nodrt * (nodrt + 2) * nrhs
         call mstm_mpi(mpi_command='allreduce', mpi_recv_buf_dc=amnp0, &
                       mpi_number=nsend, mpi_operation=mstm_mpi_sum, mpi_comm=mpicomm)
      end if
   end subroutine merge_to_common_origin

   subroutine distribute_from_common_origin(nodrt, amnp0, amnp, number_rhs, &
                                            origin_position, origin_host, vswf_type, mpi_comm, merge_procs, &
                                            single_sphere, sphere_translation_list)
      implicit none
      logical :: mergeprocs, tlist(number_spheres)
      logical, optional :: merge_procs, sphere_translation_list(number_spheres)
      integer :: nodrt, i, nblk, ntrani, nrhs, task, rhs, startsphere, endsphere, &
                 rank, numprocs, nsend, proc, mpicomm, noi, vtype, ohost
      integer, optional :: number_rhs, mpi_comm, vswf_type, origin_host, single_sphere
      real(8) :: r0(3)
      real(8), optional :: origin_position(3)
      complex(8) :: amnp(number_eqns, *), amnp0(0:nodrt + 1, nodrt, 2, *), rimedium(2)
      type(translation_data), pointer :: loc_tranmat
      type(translation_data), target :: tranmat

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
            endsphere = number_spheres
         end if
      else
         startsphere = 1
         endsphere = number_spheres
      end if
      if (present(sphere_translation_list)) then
         tlist = sphere_translation_list
      else
         tlist = .true.
      end if

      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      amnp(1:number_eqns, 1:nrhs) = (0.d0, 0.d0)
      task = 0
      do i = startsphere, endsphere
         if (.not. tlist(i)) cycle
         nblk = sphere_order(i) * (sphere_order(i) + 2) * 2
         if (host_sphere(i) .eq. ohost) then
            task = task + 1
            proc = mod(task, numprocs)
            if (proc .eq. rank) then
               call exteriorrefindex(i, rimedium)
               noi = sphere_order(i)
               ntrani = nodrt
               tranmat%matrix_calculated = .false.
               tranmat%vswf_type = vtype
               tranmat%translation_vector = sphere_position(:, i) - r0
               tranmat%refractive_index = rimedium
               tranmat%rot_op = (max(noi, ntrani) .ge. translation_switch_order)
               loc_tranmat => tranmat
               do rhs = 1, nrhs
                  call coefficient_translation(ntrani, 2, noi, 2, &
                                               amnp0(0:nodrt + 1, 1:nodrt, 1:2, rhs), &
                                               amnp(sphere_offset(i) + 1:sphere_offset(i) + nblk, rhs), &
                                               loc_tranmat)
               end do
               if (.not. tranmat%zero_translation) then
                  if (tranmat%rot_op) then
                     if (associated(tranmat%rot_mat)) deallocate (tranmat%rot_mat, tranmat%phi_mat, tranmat%z_mat)
                     nullify (tranmat%rot_mat, tranmat%phi_mat, tranmat%z_mat)
                     nullify (loc_tranmat)
                  else
                     if (associated(tranmat%gen_mat)) deallocate (tranmat%gen_mat)
                     nullify (tranmat%gen_mat)
                     nullify (loc_tranmat)
                  end if
               end if
            end if
         end if
      end do
      if (numprocs .gt. 1 .and. mergeprocs) then
         nsend = number_eqns * nrhs
         call mstm_mpi(mpi_command='allreduce', mpi_recv_buf_dc=amnp, &
                       mpi_number=nsend, mpi_operation=mstm_mpi_sum, mpi_comm=mpicomm)
      end if
   end subroutine distribute_from_common_origin
!
!  general efficiency factor calcuation
!  L/R formulation
!  April 2012
!  march 2013: something is always changing in this one
!
   subroutine lrsphereqeff(ri0, nodr, npol, xsp, anp, gnpinc, gnp, qe, qs, qa)
      implicit none
      integer :: nodr, npol, m, n, p, p1, p2, k, ma, na, s, t, ss, st
      real(8) :: qe(2 * npol - 1), qa(2 * npol - 1), qs(2 * npol - 1), const, xsp, qi(2 * npol - 1)
      complex(8) :: anp(0:nodr + 1, nodr, 2, npol), gnp(0:nodr + 1, nodr, 2, npol), &
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
         call cricbessel(nodr, xri(p), psi(0, p))
         call crichankel(nodr, xri(p), xi(0, p))
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
               cterm = cmplx(0.d0, 1.d0, kind=kind(0.0d0)) * conjg(1./ri0(t)) / ri0(s)
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

   end subroutine lrsphereqeff
!
! calling routine for efficiency calculation
! april 2012: lr formulation
! february 2013:  number of rhs, mpi comm options added.
!
   subroutine qefficiencyfactors(nsphere, npol, amnp, gmnp0, qeff, &
                                 mpi_comm)
      implicit none
      integer :: nsphere, i, nblk, noff, neqns, nodr, &
                 npol, mpicomm, p, &
                 b11, b12, rank, numprocs
      integer, optional :: mpi_comm
      real(8) :: qeff(3, 2 * npol - 1, nsphere), qe(2 * npol - 1), qa(2 * npol - 1), qs(2 * npol - 1)
      real(8), allocatable :: qeffi(:, :)
      complex(8) :: amnp(number_eqns, npol), gmnp0(number_eqns, npol), ri0(2)
      complex(8), allocatable :: amnpi(:, :, :, :), gmnpi(:, :, :, :), fmnpi(:, :, :, :), &
                                 gmnp(:, :)
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
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
               call onemiecoeffmult(i, nodr, amnp(b11:b12, p), gmnp(b11:b12, p), 'i')
               noff = noff + nblk * number_field_expansions(i)
            end do
         else
            call sphereinteraction(neqns, 1, amnp(:, p), gmnp(:, p), &
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
         call exteriorrefindex(i, ri0)
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
         call lrsphereqeff(ri0, nodr, npol, sphere_radius(i), amnpi, fmnpi, gmnpi, &
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
   end subroutine qefficiencyfactors

   subroutine qtotcalc(nsphere, nrow, xgeff, qeffp, qabsvol, qefftot)
      implicit none
      integer :: nsphere, nrow, i, j
      real(8) :: qeffp(3, nrow, nsphere), qefftot(3, nrow), qabsvol(nrow, nsphere), xgeff
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

   end subroutine qtotcalc

   subroutine waveguide_mode_scattering(amnp, qsevan, mpi_comm)
      implicit none
      integer :: i, p, pole, mpicomm
      integer, optional :: mpi_comm
      real(8) :: qsevan(2)
      complex(8) :: amnp(number_eqns, 2)
      complex(8), allocatable :: amn(:), gmn(:)
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
            call sphereinteraction(number_eqns, 1, amn, gmn, &
                                   mie_mult=(/.false./), initial_run=.true., skip_external_translation=.true., &
                                   mpi_comm=mpicomm)
            do i = 1, number_spheres
               if (host_sphere(i) .ne. 0) cycle
               call lr_mode_transformation(sphere_order(i), amn(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i)), &
                                           amn(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i)))
               call lr_mode_transformation(sphere_order(i), gmn(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i)), &
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
      real(8) :: qbsca(2, 0:1), qt(2), targetz, qt2(2, 0:1)
      complex(8) :: amn(number_eqns, 2)
      complex(8), allocatable :: amn1(:, :), amn2(:, :)
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, &
                    mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, &
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
         call mstm_mpi(mpi_command='allreduce', mpi_number=nsend, &
                       mpi_send_buf_dp=qt2, mpi_recv_buf_dp=qbsca, &
                       mpi_operation=mstm_mpi_sum, &
                       mpi_comm=mpicomm)
      end if
   end subroutine boundary_scattering

   subroutine boundary_extinction(amn, alpha, sinc, dir, qbext, common_origin)
      implicit none
      logical :: comorg
      logical, optional :: common_origin
      integer :: k, dir
      real(8) :: qt(2), targetz, qbext(2, 0:1), alpha, sinc
      complex(8) :: amn(*)
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
      real(8) :: qbsca(2, 0:1), targetz
      complex(8) :: amn(2 * t_matrix_order * (t_matrix_order + 2), 2)

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
      real(8) :: qbsca(2, 2), inteps, mindiv, c0, c1
      complex(8) :: amnp(*), cq(2, 2)
      data maxnumdiv, inteps, mindiv, c0, c1/3, 1.d-6, 1.d-4, 0.d0, 1.d0/
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, &
                    mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, &
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
            call gkintegrate(2, c0, c1, hemispherical_integrand, cq(:, p), subdiv, &
                             errorcodes, inteps, mindiv, maxnumdiv)
         end do
         if (numprocs .gt. 1) then
            if (rank .eq. 1) then
               call mstm_mpi(mpi_command='send', mpi_number=2, &
                             mpi_send_buf_dc=cq(:, 2), mpi_rank=0, &
                             mpi_comm=mpicomm)
            elseif (rank .eq. 0) then
               call mstm_mpi(mpi_command='recv', mpi_number=2, &
                             mpi_recv_buf_dc=cq(:, 2), mpi_rank=1, &
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
      real(8) :: ct, sm(2)
      complex(8) :: q(ntot)
      q = 0.d0
      if (hemispherical_single_origin) then
         call numerical_sm_azimuthal_average_so(hemispherical_amnp, t_matrix_order, ct, sm, &
                                                rotate_plane=.false., normalize_s11=.false., s11_only=.true.)
         q(1:2) = sm(1:2)
      else
         call numerical_sm_azimuthal_average_mo(hemispherical_amnp, ct, sm, &
                                                rotate_plane=.false., s11_only=.true.)
         q(1:2) = sm(1:2)
      end if
   end subroutine hemispherical_integrand

   subroutine extinction_theorem(amnp, sinc, sdir, alpha, targetz, qe, common_origin)
      implicit none
      logical :: comorg
      logical, optional :: common_origin
      integer :: sdir, p, tlay
      real(8) :: sinc, alpha, qe(2), targetz, sourcez, const3(3, 2)
      complex(8) :: amnp(*), s, sourceri, targetri, gfs(2, 2, 2), skz, tkz, sa(4), amp(2, 2)
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
      tlay = layer_id(targetz)
      targetri = layer_ref_index(tlay)
      s = sinc
      call layer_gf(s, sourcez, targetz, gfs, skz, tkz, include_direct=.true.)
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
      smscale = 1.d0 / cross_section_radius**2 * 64.d0 * pi / wx / wy
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
               kx = 2.d0 * pi * dble(ix) / wx + k0x
               ky = 2.d0 * pi * dble(iy) / wy + k0y
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
      qsca = qsca / cross_section_radius**2 * 64.d0 * pi / wx / wy
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
         call layervsh(s, phi, targetz, dir, sphere_position(:, i), sphere_order(i), pmnpi)
         amnpi(:, :) = amnp(sphere_offset(i) + 1:sphere_offset(i) + sphere_block(i), :)
         do p = 1, 2
            call lr_mode_transformation(sphere_order(i), amnpi(:, p), amnpi(:, p))
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
      call layervsh(s, phi, targetz, dir, cluster_origin, t_matrix_order, pmnp)
      sa(1) = 0.5d0 * sum(pmnp(:, 2) * amnp(:, 2))
      sa(2) = 0.5d0 * sum(pmnp(:, 1) * amnp(:, 1))
      sa(3) = -0.5d0 * sum(pmnp(:, 1) * amnp(:, 2))
      sa(4) = -0.5d0 * sum(pmnp(:, 2) * amnp(:, 1))
   end subroutine common_origin_amplitude_matrix

   subroutine amplitude_to_scattering_matrix(sa, sm)
      implicit none
      integer :: i, j
      real(8) :: sm(4, 4)
      complex(8) :: sa(4), sp(4, 4)
      do i = 1, 4
         do j = 1, 4
            sp(i, j) = sa(i) * conjg(sa(j))
         end do
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

   subroutine multiple_origin_scatteringmatrix(amnp, ct, phi, csca, sa, sm, rotate_plane, s11_only)
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
   end subroutine multiple_origin_scatteringmatrix
!
!  scattering amplitude sa and matrix sm calculation
!
!  original: 15 January 2011
!  revised: 21 February 2011: S11 normalization changed
!  april 2013: moved things around to try to get it to work.
!
   subroutine scatteringmatrix(amn0, nodrt, ct, phi, sa, sm, rotate_plane, normalize_s11, &
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
      call taufunc(ct, nodrt, tau)
      cphi = cos(phi)
      sphi = sin(phi)
      ephi = cmplx(cphi, sphi, kind=kind(0.0d0))
      call ephicoef(ephi, nodrt, ephim)
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
      sm(1:nelem) = sm(1:nelem) / (4.d0 * pi) / dble(layer_ref_index(0))
   end subroutine scatteringmatrix

   subroutine numerical_sm_azimuthal_average_so(amn0, nodrt, ct, sm, rotate_plane, normalize_s11, &
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
         phi = 2.d0 * pi * dble(i - 1) / dble(numang)
         call scatteringmatrix(amn0, nodrt, ct, phi, sa, smt, &
                               rotate_plane=rotate, normalize_s11=norms11, s11_only=s11only)
         sm(1:nelem) = sm(1:nelem) + smt(1:nelem)
      end do
      sm(1:nelem) = sm(1:nelem) / dble(numang)
   end subroutine numerical_sm_azimuthal_average_so

   subroutine numerical_sm_azimuthal_average_mo(amnp, ct, sm, number_angles, rotate_plane, s11_only)
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
      csca = pi * 2.d0
      do i = 1, numang
         phi = 2.d0 * pi * dble(i - 1) / dble(numang)
         call multiple_origin_scatteringmatrix(amnp, ct, phi, csca, sa, smt, &
                                               rotate_plane=rotate, s11_only=s11only)
         sm(1:nelem) = sm(1:nelem) + smt(1:nelem)
      end do
      sm(1:nelem) = sm(1:nelem) / dble(numang)
   end subroutine numerical_sm_azimuthal_average_mo

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
   subroutine s11expansion(amn0, nodrt, mmax, nodrg, gmn)
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
            call vcfunc(-1, l, 1, n, vc2)
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
                     call vcfunc(-k, l, m, n, vc1)
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
   end subroutine s11expansion
!
!  calculate azimuth--averaged scattering matrix from expansion, for cos(theta) = ct
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!  this is currently not used in v. 3.0
!
   subroutine fosmcalc(ntot, s00, s02, sp22, sm22, ct, sm, normalize_s11)
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
      call rotcoef(ct, 2, 2 * ntot, dc)

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
   end subroutine fosmcalc
!
!  determine the generalized spherical function expansion for the azimuth-averaged scattering matrix
!  corresponding to the target-based scattering field expansion of amnp.
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: fixed flush call.
!

   subroutine fosmexpansion(ntot, amnp, s00, s02, sp22, sm22, mpi_comm)
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
            call vcfunc(-1, n, 1, l, cm1p1)
            call vcfunc(-1, n, -1, l, cm1m1)
            do m = -min(n, l + 2), min(n, l + 2)
               m1m = (-1)**m
               if (abs(m) .le. l) then
                  call vcfunc(-m, n, m, l, cmmpm)
               else
                  cmmpm = 0.d0
               end if
               if (abs(-2 + m) .le. l) then
                  call vcfunc(-m, n, -2 + m, l, cmmm2pm)
               else
                  cmmm2pm = 0.d0
               end if
               if (abs(2 + m) .le. l) then
                  call vcfunc(-m, n, 2 + m, l, cmmp2pm)
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

   end subroutine fosmexpansion

!
!  compute the coefficients for the GSF expansion of the random orientation
!  scattering matrix.
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!  January 2012: added computation of coherent field average
!  February 2013:  added number processors option.
!
   subroutine ranorientscatmatrix(tmatrixfile, sm, smcf, beam_width, number_processors, mpi_comm, &
                                  mean_t_matrix, override_order, keep_quiet)
      implicit none
      logical :: symmetrical, nkq
      logical, optional :: keep_quiet
      integer :: nodr, nodrw, nodr2, m, n, p, k, l, q, t, v, u, w, nblk, kl, mn, nn1, tn, &
                 lmax, ll1, tvl, ku, k1, ns, ik, ik1, m1, nu, n1s, n1e, nu1, p1, n1max, &
                 in, n1, i, kt, nodrt, nodrrhs, mnm, klm, ikm, &
                 rank, numprocs, &
                 nblkrhs, numprocscalc, orig_group, new_group, &
                 new_comm, new_rank, nblkw, wv, sizedm, sizetm, nread, mpicomm, rank0
      integer, allocatable :: windex(:), vindex(:), wvindex(:), wvnum(:), group_list(:)
      integer, optional :: number_processors, mpi_comm, override_order
      real(8) :: sm(4, 4, 0:*), fl, xv, fl2, cbeam, gbn, wvperproc, sum, &
                 time1, time2, smcf(4, 4, 0:*), qextt, qscatt
      real(8), allocatable :: vc(:)
      real(8), optional :: beam_width
      complex(8) :: ci, cin, a, tct, tcp(2)
      complex(8), allocatable :: aw(:, :, :), bw(:, :, :), cw(:), &
                                 dw(:), pp(:, :, :), bm(:, :, :), &
                                 am(:, :, :), fm(:, :, :, :, :), bmcf(:, :, :), &
                                 amcf(:, :, :), fmcf(:, :, :, :, :), awcf(:, :, :), &
                                 bwcf(:, :, :), cwcf(:), dwcf(:)
      complex(8), allocatable :: dm(:, :, :, :, :, :), dmcf(:, :, :, :, :, :)
      complex(8), optional :: mean_t_matrix(*)
      complex(4), allocatable :: tc(:, :, :, :)
      character(len=30) :: tmatrixfile
      data ci/(0.d0, 1.d0)/
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      if (present(keep_quiet)) then
         nkq = .not. keep_quiet
      else
         nkq = .true.
      end if
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank0)
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      if (present(beam_width)) then
         cbeam = beam_width
      else
         cbeam = 0.d0
      end if
      if (present(number_processors)) then
         numprocscalc = min(number_processors, numprocs)
      else
         numprocscalc = numprocs
      end if
      allocate (group_list(numprocscalc))
      group_list = (/(i, i=0, numprocscalc - 1)/)
      call mstm_mpi(mpi_command='group', mpi_group=orig_group, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='incl', mpi_group=orig_group, &
                    mpi_size=numprocscalc, mpi_new_group_list=group_list, &
                    mpi_new_group=new_group, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='create', mpi_group=new_group, &
                    mpi_new_comm=new_comm, mpi_comm=mpicomm)

      xv = cross_section_radius
      if (rank .le. numprocscalc - 1) then

         call mstm_mpi(mpi_command='rank', mpi_comm=new_comm, &
                       mpi_rank=new_rank)
         if (rank .eq. 0) time1 = mstm_mpi_wtime()
!
!  read the T matrix from the file
!
         do i = 0, numprocscalc - 1
            if (i .eq. new_rank) then
               open (3, file=tmatrixfile)
               read (3, *) nodr, nodrrhs
               if (present(override_order)) then
                  nodr = override_order
                  nodrrhs = nodr
               end if
               symmetrical = nodr .eq. nodrrhs
               nodrt = nodr
               nodr2 = nodr + nodr
               nodrw = nodr2
               nblk = nodr * (nodr + 2)
!                  nodrrhs=nodr
               nblkrhs = nodrrhs * (nodrrhs + 2)
               sizetm = 4 * nblk * nblkrhs
               allocate (tc(2, nblk, 2, nblkrhs))
               tc = (0., 0.)
               do l = 1, nodrrhs
                  if (present(mean_t_matrix)) then
                     mean_t_matrix(2 * l - 1:2 * l) = 0.d0
                  end if
                  do k = -l, l
                     kl = l * (l + 1) + k
                     klm = l * (l + 1) - k
                     do q = 1, 2
                        if (symmetrical) then
                           nread = l
                        else
                           nread = nodr
                        end if
                        do n = 1, nread
                           do m = -n, n
                              mn = n * (n + 1) + m
                              mnm = n * (n + 1) - m
                              do p = 1, 2
                                 read (3, *) tcp(p)
                                 tc(p, mn, q, kl) = tcp(p)
                              end do
                              if (present(mean_t_matrix) .and. n .eq. l .and. m .eq. k) then
                                 mean_t_matrix(2 * (l - 1) + q) = mean_t_matrix(2 * (l - 1) + q) + tcp(q)
                              end if
                              if (n .lt. l .and. symmetrical) then
                                 ikm = (-1)**(m + k)
                                 do p = 1, 2
                                    tc(q, klm, p, mnm) = tc(p, mn, q, kl) * ikm
                                 end do
                              end if
                           end do
                        end do
                     end do
                  end do
                  if (present(mean_t_matrix)) then
                     mean_t_matrix(2 * l - 1:2 * l) = mean_t_matrix(2 * l - 1:2 * l) / dble(2 * l + 1)
                  end if
               end do
!if(rank.eq.0) then
!open(30,file='meantmatrix.dat')
!do n=1,nodrt
!do p=1,2
!tct=0.
!do m=-n,n
!mn=n*(n+1)+m
!tct=tct+tc(p,mn,p,mn)
!enddo
!tct=tct/dble(n+n+1)
!write(30,'(2i5,2es13.5)') n,p,tct
!enddo
!enddo
!close(30)
!endif
               qextt = 0.d0
               qscatt = 0.d0
               do l = 1, nodrrhs
                  do k = -l, l
                     kl = l * (l + 1) + k
                     do q = 1, 2
                        do n = 1, nodr
                           do m = -n, n
                              mn = n * (n + 1) + m
                              do p = 1, 2
                                 qscatt = qscatt + cabs(tc(p, mn, q, kl))**2
                              end do
                           end do
                        end do
                        if (l .le. nodr) qextt = qextt - real(tc(q, kl, q, kl))
                     end do
                  end do
               end do
               qextt = qextt * 2./cross_section_radius**2
               qscatt = qscatt * 2./cross_section_radius**2
               close (3)
            end if
            call mstm_mpi(mpi_command='barrier', mpi_comm=new_comm)
         end do
         if (rank0 .eq. 0 .and. nkq) then
            write (run_print_unit, '('' t matrix ext, sca:'',2e13.5)') qextt, qscatt
         end if
!!
!!  send to the other processors
!!
!            if(numprocscalc.gt.1) then
!               call mstm_mpi(mpi_command='bcast',mpi_send_buf_c=tc, &
!                    mpi_number=sizetm,mpi_rank=0,mpi_comm=new_comm)
!            endif

         allocate (vc(0:4 * nodr + 2))
         allocate (aw(0:2, -1:1, 0:nodrw), bw(0:2, -1:1, 0:nodrw), cw(0:nodrw), &
                   dw(0:nodrw), pp(nodr, 2, 2), bm(2, nodr * (nodr + 2), 2), &
                   am(2, nodrrhs + 1, 2), fm(3, nodr, 2, nodr, 2), bmcf(2, nodr * (nodr + 2), 2), &
                   amcf(2, nodrrhs + 1, 2), fmcf(3, nodr, 2, nodr, 2), awcf(0:2, -1:1, 0:nodrw), &
                   bwcf(0:2, -1:1, 0:nodrw), cwcf(0:nodrw), dwcf(0:nodrw))
         allocate (dm(-nodr - 1:nodr + 1, 3, nodr, 2, nodr, 2), dmcf(-nodr - 1:nodr + 1, 3, nodr, 2, nodr, 2))
         if (rank0 .eq. 0 .and. nkq) then
            time2 = mstm_mpi_wtime() - time1
            call timewrite(run_print_unit, ' t matrix read time:', time2)
            time1 = mstm_mpi_wtime()
         end if
         dm = (0.d0, 0.d0)
         dmcf = (0.d0, 0.d0)
         sizedm = size(dm)
         call init(nodr2)
!
!  compute the GB modified T matrix
!
         do n = 1, nodrrhs
            gbn = dexp(-((dble(n) + .5d0) * cbeam)**2.)
            cin = ci**(n + 1)
            pp(n, 1, 1) = -.5d0 * cin * fnr(n + n + 1) * gbn
            pp(n, 2, 1) = -pp(n, 1, 1)
            pp(n, 1, 2) = -pp(n, 1, 1)
            pp(n, 2, 2) = pp(n, 2, 1)
         end do
         do n = 1, nodr
            nn1 = n * (n + 1)
            do m = -n, n
               mn = nn1 + m
               do p = 1, 2
                  do l = 1, nodrrhs
                     do k = -l, l
                        kl = l * (l + 1) + k
                        a = tc(p, mn, 1, kl)
                        tc(p, mn, 1, kl) = tc(p, mn, 1, kl) * pp(l, 1, 1) &
                                           + tc(p, mn, 2, kl) * pp(l, 1, 2)
                        tc(p, mn, 2, kl) = a * pp(l, 2, 1) + tc(p, mn, 2, kl) * pp(l, 2, 2)
                     end do
                  end do
               end do
            end do
         end do
!
!  determine the distribution of work load among the processors
!
         nblkw = nodr2 * (nodr2 + 2) + 1
         allocate (windex(nblkw), vindex(nblkw), wvindex(0:numprocscalc - 1), &
                   wvnum(0:numprocscalc - 1))
         w = 0
         do n = 0, nodr2
            do m = -n, n
               w = w + 1
               windex(w) = n
               vindex(w) = m
            end do
         end do
         wvperproc = dble(nblkw) / dble(numprocscalc)
         sum = 0.
         do i = 0, numprocscalc - 1
            wvindex(i) = floor(sum)
            sum = sum + wvperproc
         end do
         do i = 0, numprocscalc - 2
            wvnum(i) = wvindex(i + 1) - wvindex(i)
         end do
         wvnum(numprocscalc - 1) = nblkw - wvindex(numprocscalc - 1)
         if (rank0 .eq. 0 .and. nkq) then
            write (run_print_unit, '('' d matrix calculation, order+degree per proc.:'',f9.2)') &
               wvperproc
            flush (run_print_unit)
         end if
!
!  the big loop
!
         do wv = wvindex(rank) + 1, wvindex(rank) + wvnum(rank)
            w = windex(wv)
            v = vindex(wv)
            bm = (0.d0, 0.d0)
            bmcf = (0.d0, 0.d0)
            do n = 1, nodr
               nn1 = n * (n + 1)
               do l = max(1, abs(w - n)), min(nodrrhs, w + n)
                  am(1, l, 1) = 0.d0
                  am(1, l, 2) = 0.d0
                  am(2, l, 1) = 0.d0
                  am(2, l, 2) = 0.d0
                  amcf(1, l, 1) = 0.d0
                  amcf(1, l, 2) = 0.d0
                  amcf(2, l, 1) = 0.d0
                  amcf(2, l, 2) = 0.d0
               end do
               do t = -n, n
                  tn = nn1 + t
                  lmax = min(nodrrhs, w + n)
                  call vcfunc(v, w, -t, n, vc)
                  do l = max(1, abs(v - t), abs(n - w)), lmax
                     ll1 = l * (l + 1)
                     tvl = ll1 + t - v
                     do k = 1, 2
                        do p = 1, 2
                           am(k, l, p) = am(k, l, p) + vc(l) * tc(p, tn, k, tvl)
                           if (l .eq. n .and. v .eq. 0) then
! may 2019: orientation averaged T matrix is azimuthal independent
                              tct = 0.d0
                              do kt = -n, n
                                 tct = tct + tc(p, nn1 + kt, k, nn1 + kt)
                              end do
                              amcf(k, l, p) = amcf(k, l, p) + vc(l) * tct / dble(n + n + 1)
!                                 amcf(k,l,p)=amcf(k,l,p)+vc(l)*tcm(p,n)
!                                 amcf(k,l,p)=amcf(k,l,p)+vc(l)*tc(p,tn,k,tvl)
                           end if
                        end do
                     end do
                  end do
               end do
               do m = -n, n
                  mn = nn1 + m
                  do k = 1, 2
                     u = m - (-3 + 2 * k)
                     if (abs(u) .le. w) then
                        lmax = min(nodrrhs, w + n)
                        call vcfunc(-u, w, m, n, vc)
                        do l = max(1, abs(w - n)), lmax
                           fl = -(-1)**l * vc(l) / dble(l + l + 1)
                           do p = 1, 2
                              bm(k, mn, p) = bm(k, mn, p) + am(k, l, p) * fl
                              if (v .eq. 0) then
                                 bmcf(k, mn, p) = bmcf(k, mn, p) + amcf(k, l, p) * fl
                              end if
                           end do
                        end do
                     end if
                  end do
               end do
            end do

!               do k=1,2
!                  do n=1,nodr*(nodr+2)
!                     do p=1,2
!                        bmcf(k,n,p)=bm(k,n,p)-bmcf(k,n,p)
!                     enddo
!                  enddo
!               enddo

            do u = -min(w, nodr + 1), min(w, nodr + 1)
               do ku = 1, 3
                  if (ku .eq. 1) then
                     k = -1
                     k1 = -1
                  elseif (ku .eq. 2) then
                     k = 1
                     k1 = 1
                  else
                     k = 1
                     k1 = -1
                  end if
                  m = u + k
                  ns = max(1, abs(m))
                  ik = (k + 1) / 2 + 1
                  ik1 = (k1 + 1) / 2 + 1
                  m1 = u + k1
                  do n = ns, nodr
                     nu = n * (n + 1) + m
                     n1s = max(1, abs(m1), n - nodrw)
                     n1e = min(nodr, n + nodrw)
                     do n1 = n1s, n1e
                        cin = ci**(n - n1)
                        nu1 = n1 * (n1 + 1) + m1
                        fl = -fnr(n + n + 1) * fnr(n1 + n1 + 1) * dble(w + w + 1)
                        do p = 1, 2
                           do p1 = 1, 2
                              a = bm(ik, nu, p) * cin * fl * conjg(bm(ik1, nu1, p1))
                              dm(u, ku, n, p, n1, p1) = dm(u, ku, n, p, n1, p1) + a
                              if (v .eq. 0) then
                                 a = bmcf(ik, nu, p) * cin * fl * conjg(bmcf(ik1, nu1, p1))
                                 dmcf(u, ku, n, p, n1, p1) = dmcf(u, ku, n, p, n1, p1) + a
                              end if
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate (tc)

         call mstm_mpi(mpi_command='allreduce', mpi_recv_buf_dc=dm, &
                       mpi_number=sizedm, mpi_operation=mstm_mpi_sum, mpi_comm=new_comm)
         call mstm_mpi(mpi_command='allreduce', mpi_recv_buf_dc=dmcf, &
                       mpi_number=sizedm, mpi_operation=mstm_mpi_sum, mpi_comm=new_comm)
         if (rank0 .eq. 0 .and. nkq) then
            time2 = mstm_mpi_wtime() - time1
            call timewrite(run_print_unit, ' d matrix time:', time2)
            time1 = mstm_mpi_wtime()
         end if
!
!  compute the expansion coefficients
!
         aw = 0.d0
         bw = 0.d0
         cw = 0.d0
         dw = 0.d0
         awcf = 0.d0
         bwcf = 0.d0
         cwcf = 0.d0
         dwcf = 0.d0
         do w = 0, nodrw
            do n = 1, nodr
               n1s = max(1, abs(n - w))
               n1e = min(nodr, n + w)
               do n1 = n1s, n1e
                  do k = 1, 3
                     do p = 1, 2
                        do p1 = 1, 2
                           fm(k, n, p, n1, p1) = 0.
                           fmcf(k, n, p, n1, p1) = 0.
                        end do
                     end do
                  end do
               end do
            end do
            do u = -nodr - 1, nodr + 1
               do k = -1, 1, 2
                  m = u + k
                  ik = (k + 1) / 2 + 1
                  ns = max(1, abs(m))
                  do n = ns, nodr
                     n1max = min(w + n, nodr)
                     call vcfunc(m, n, 0, w, vc)
                     do n1 = ns, nodr
                        if ((n + n1 .lt. w) .or. (abs(n - n1) .gt. w)) cycle
                        fl = -(-1)**n * vc(n1) * fnr(w + w + 1) / fnr(n1 + n1 + 1)
                        do p = 1, 2
                           do p1 = 1, 2
                              fm(ik, n, p, n1, p1) = fm(ik, n, p, n1, p1) &
                                                     + dm(u, ik, n, p, n1, p1) * fl
                              fmcf(ik, n, p, n1, p1) = fmcf(ik, n, p, n1, p1) &
                                                       + dmcf(u, ik, n, p, n1, p1) * fl
                           end do
                        end do
                     end do
                  end do
               end do
               if (w .lt. 2) cycle
               m = u + 1
               m1 = u - 1
               ns = max(1, abs(m))
               n1s = max(1, abs(m1))
               do n = ns, nodr
                  n1max = min(w + n, nodr)
                  call vcfunc(m, n, -2, w, vc)
                  do n1 = n1s, nodr
                     if ((n + n1 .lt. w) .or. (abs(n - n1) .gt. w)) cycle
                     fl = -(-1)**n * vc(n1) * fnr(w + w + 1) / fnr(n1 + n1 + 1)
                     do p = 1, 2
                        do p1 = 1, 2
                           fm(3, n, p, n1, p1) = fm(3, n, p, n1, p1) &
                                                 + dm(u, 3, n, p, n1, p1) * fl
                           fmcf(3, n, p, n1, p1) = fmcf(3, n, p, n1, p1) &
                                                   + dmcf(u, 3, n, p, n1, p1) * fl
                        end do
                     end do
                  end do
               end do
            end do
            do n = 1, nodr
               n1s = max(1, abs(n - w))
               n1e = min(nodr, n + w)
               in = (-1)**n
               n1max = min(w + n, nodr)
               call vcfunc(1, n, 0, w, vc)
               do n1 = n1s, n1e
                  fl = 2.d0 * in * vc(n1) * fnr(w + w + 1) / fnr(n1 + n1 + 1)
                  i = mod(n + n1 - w, 2) + 1
                  do p = 1, 2
                     p1 = (2 - i) * p + (i - 1) * (3 - p)
                     do k = -1, 1, 2
                        ik = (k + 1) / 2 + 1
                        aw(0, k, w) = aw(0, k, w) + fm(ik, n, p, n1, p1) * fl
                        bw(0, k, w) = bw(0, k, w) + fm(ik, n, p, n1, 3 - p1) * fl
                        awcf(0, k, w) = awcf(0, k, w) + fmcf(ik, n, p, n1, p1) * fl
                        bwcf(0, k, w) = bwcf(0, k, w) + fmcf(ik, n, p, n1, 3 - p1) * fl
                     end do
                     bw(2, 0, w) = bw(2, 0, w) + fm(3, n, p, n1, 3 - p1) * fl
                     aw(2, 0, w) = aw(2, 0, w) + fm(3, n, p, n1, p1) * fl
                     bwcf(2, 0, w) = bwcf(2, 0, w) + fmcf(3, n, p, n1, 3 - p1) * fl
                     awcf(2, 0, w) = awcf(2, 0, w) + fmcf(3, n, p, n1, p1) * fl
                  end do
               end do
               if (w .lt. 2) cycle
               call vcfunc(1, n, -2, w, vc)
               do n1 = n1s, n1e
                  fl = 2.d0 * in * vc(n1) * fnr(w + w + 1) / fnr(n1 + n1 + 1)
                  i = mod(n + n1 - w, 2) + 1
                  do p = 1, 2
                     p1 = (2 - i) * p + (i - 1) * (3 - p)
                     do k = -1, 1, 2
                        ik = (k + 1) / 2 + 1
                        aw(2, k, w) = aw(2, k, w) + fm(ik, n, p, n1, p1) * fl * (-1)**p1
                        bw(2, k, w) = bw(2, k, w) + fm(ik, n, p, n1, 3 - p1) * fl * (-1)**(3 - p1)
                        awcf(2, k, w) = awcf(2, k, w) + fmcf(ik, n, p, n1, p1) * fl * (-1)**p1
                        bwcf(2, k, w) = bwcf(2, k, w) &
                                        + fmcf(ik, n, p, n1, 3 - p1) * fl * (-1)**(3 - p1)
                     end do
                  end do
                  fl2 = 2.*(-1)**(n1 + w) * vc(n1) * fnr(w + w + 1) / fnr(n1 + n1 + 1)
                  do p = 1, 2
                     do p1 = 1, 2
                        cw(w) = cw(w) + fm(3, n, p, n1, p1) * fl * (-1)**p1
                        dw(w) = dw(w) + fm(3, n, p, n1, p1) * fl2 * (-1)**p
                        cwcf(w) = cwcf(w) + fmcf(3, n, p, n1, p1) * fl * (-1)**p1
                        dwcf(w) = dwcf(w) + fmcf(3, n, p, n1, p1) * fl2 * (-1)**p
                     end do
                  end do
               end do
            end do
         end do
!            do w=0,nodrw
!               do k=-1,1
!                  do i=0,2
!                     aw(i,k,w)=aw(i,k,w)*2./xv/xv
!                     bw(i,k,w)=bw(i,k,w)*2./xv/xv
!                     awcf(i,k,w)=awcf(i,k,w)*2./xv/xv
!                     bwcf(i,k,w)=bwcf(i,k,w)*2./xv/xv
!                  enddo
!               enddo
!               cw(w)=cw(w)*2./xv/xv
!               dw(w)=dw(w)*2./xv/xv
!               cwcf(w)=cwcf(w)*2./xv/xv
!               dwcf(w)=dwcf(w)*2./xv/xv
!            enddo
         do n = 0, nodrw
            sm(1, 1, n) = aw(0, -1, n) + aw(0, 1, n)
            sm(1, 2, n) = aw(2, -1, n) + aw(2, 1, n)
            sm(1, 3, n) = 2.d0 * aimag(aw(2, 0, n))
            sm(1, 4, n) = aw(0, 1, n) - aw(0, -1, n)
            sm(2, 2, n) = dw(n)
            sm(2, 3, n) = aimag(dw(n))
            sm(2, 4, n) = aw(2, 1, n) - aw(2, -1, n)
            sm(3, 1, n) = -aimag(bw(2, -1, n) + bw(2, 1, n))
            sm(3, 2, n) = aimag(cw(n))
            sm(3, 3, n) = cw(n)
            sm(3, 4, n) = aimag(bw(2, -1, n) - bw(2, 1, n))
            sm(4, 1, n) = bw(0, -1, n) + bw(0, 1, n)
            sm(4, 2, n) = 2.*bw(2, 0, n)
            sm(4, 4, n) = bw(0, 1, n) - bw(0, -1, n)

            smcf(1, 1, n) = awcf(0, -1, n) + awcf(0, 1, n)
            smcf(1, 2, n) = awcf(2, -1, n) + awcf(2, 1, n)
            smcf(1, 3, n) = 2.d0 * aimag(awcf(2, 0, n))
            smcf(1, 4, n) = awcf(0, 1, n) - awcf(0, -1, n)
            smcf(2, 2, n) = dwcf(n)
            smcf(2, 3, n) = aimag(dwcf(n))
            smcf(2, 4, n) = awcf(2, 1, n) - awcf(2, -1, n)
            smcf(3, 1, n) = -aimag(bwcf(2, -1, n) + bwcf(2, 1, n))
            smcf(3, 2, n) = aimag(cwcf(n))
            smcf(3, 3, n) = cwcf(n)
            smcf(3, 4, n) = aimag(bwcf(2, -1, n) - bwcf(2, 1, n))
            smcf(4, 1, n) = bwcf(0, -1, n) + bwcf(0, 1, n)
            smcf(4, 2, n) = 2.*bwcf(2, 0, n)
            smcf(4, 4, n) = bwcf(0, 1, n) - bwcf(0, -1, n)

         end do
! patch 10-22
         sm(:, :, 0:nodrw) = sm(:, :, 0:nodrw) / dble(layer_ref_index(0)) / 2.d0
         smcf(:, :, 0:nodrw) = smcf(:, :, 0:nodrw) / dble(layer_ref_index(0)) / 2.d0
!            sm(:,:,0:nodrw)=sm(:,:,0:nodrw)*4.d0*pi
!            smcf(:,:,0:nodrw)=smcf(:,:,0:nodrw)*4.d0*pi
!
!  normalization
!
!            qsca0=sm(1,1,0)
!            do n=0,nodrw
!               do i=1,4
!                  do j=1,4
!                     sm(i,j,n)=sm(i,j,n)/qsca0
!                     smcf(i,j,n)=smcf(i,j,n)/qsca0
!                  enddo
!               enddo
!            enddo
         if (rank0 .eq. 0 .and. nkq) then
            time2 = mstm_mpi_wtime() - time1
            call timewrite(run_print_unit, ' scat matrix coef time:', time2)
         end if
         deallocate (windex, vindex, wvindex, wvnum, dm)
         deallocate (vc)
         deallocate (aw, bw, cw, dw, pp, bm, am, fm, bmcf, amcf, fmcf, awcf, &
                     bwcf, cwcf, dwcf)
      end if
      call mstm_mpi(mpi_command='barrier', mpi_comm=mpicomm)
   end subroutine ranorientscatmatrix
!
!  calculation of the RO scattering matrix from the GSF expansion
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!
   subroutine ranorienscatmatrixcalc(ct, smc, nodrexp, sm)
      implicit none
      integer :: nodrexp, n, nn0, nnp2, nnm2
      real(8) :: smc(4, 4, 0:nodrexp), sm(4, 4), dc(-2:2, 0:nodrexp * (nodrexp + 2)), &
                 ct
!
!     dc is the normalized generalized spherical function
!     dc(k,n*(n+1)+m) = ((n-k)!(n+m)!/(n+k)!/(n-m)!)^(1/2) D^k_{mn},
!     where D^k_{mn} is defined in M&M JOSA 96
!
      call rotcoef(ct, 2, nodrexp, dc)
      sm = 0.d0
      do n = 0, nodrexp
         nn0 = n * (n + 1)
         nnp2 = nn0 + 2
         nnm2 = nn0 - 2
         sm(1, 1) = sm(1, 1) + dc(0, nn0) * smc(1, 1, n)
         sm(1, 4) = sm(1, 4) + dc(0, nn0) * smc(1, 4, n)
         sm(4, 1) = sm(4, 1) + dc(0, nn0) * smc(4, 1, n)
         sm(4, 4) = sm(4, 4) + dc(0, nn0) * smc(4, 4, n)
         if (n .ge. 2) then
            sm(1, 2) = sm(1, 2) + dc(2, nn0) * smc(1, 2, n)
            sm(2, 4) = sm(2, 4) + dc(2, nn0) * smc(2, 4, n)
            sm(3, 4) = sm(3, 4) + dc(2, nn0) * smc(3, 4, n)
            sm(1, 3) = sm(1, 3) + dc(2, nn0) * smc(1, 3, n)
            sm(3, 1) = sm(1, 3) + dc(2, nn0) * smc(3, 1, n)
            sm(4, 2) = sm(4, 2) + dc(2, nn0) * smc(4, 2, n)
            sm(2, 2) = sm(2, 2) + dc(2, nnm2) * smc(2, 2, n) + dc(2, nnp2) * smc(3, 3, n)
            sm(2, 3) = sm(2, 3) + dc(2, nnp2) * smc(2, 3, n) + dc(2, nnp2) * smc(3, 2, n)
            sm(3, 3) = sm(3, 3) - dc(2, nnm2) * smc(2, 2, n) + dc(2, nnp2) * smc(3, 3, n)
            sm(3, 2) = sm(3, 2) + dc(2, nnp2) * smc(2, 3, n) - dc(2, nnp2) * smc(3, 2, n)
         end if
      end do
      sm(2, 1) = sm(1, 2)
      sm(4, 3) = -sm(3, 4)
!
!  discontiued scaling option: now done in main program
!
!            if(iscale.eq.1) then
!               do j=1,4
!                  do k=j,4
!                     if(j.ne.1.or.k.ne.1) then
!                        sm(j,k,i)=sm(j,k,i)/sm(1,1,i)
!                     endif
!                  enddo
!               enddo
!            endif
!
!    here are the VV and HH differential cross sections
!
!            gvv=.25*(sm(1,1)+sm(2,2)-2.*sm(1,2))
!            ghh=.25*(sm(1,1)+sm(2,2)+2.*sm(1,2))
!
      return
   end subroutine ranorienscatmatrixcalc

end module scatprops
