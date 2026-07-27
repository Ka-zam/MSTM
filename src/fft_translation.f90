module fft_translation
   use constants
   use gpfa_controller, only: cgpfa
   use gpfa_setup, only: setgpfa
   use mpidefs
   use intrinsics
   use numerical_tables
   use specialfuncs
   use spheredata
   use translation, only: clear_stored_translation_matrices, transform_mode_coefficients, &
                          translation_operator_state
   use mie
   implicit none
   type node_data
      integer :: number_elements
      type(linked_ilist), pointer :: members
   end type node_data
   type linked_ilist
      integer :: index
      type(linked_ilist), pointer :: next => null()
   end type linked_ilist
   type coefficient_list
      complex(8), pointer :: coefficient_vector(:, :)
   end type coefficient_list

   integer :: cell_dim(3), number_neighbor_nodes
   integer, private :: neighbor_node(3, 0:26), fft_local_host, fft_number_spheres
   integer, target :: node_order, neighbor_node_model
   integer, allocatable, private :: sphere_node(:, :)
!      real(8) :: d_cell
   real(8), private :: cell_origin(3), cell_boundary(3)
   real(8), private :: timedat(10)
   real(8), target :: cell_volume_fraction, d_cell
   complex(8), private :: host_ref_index(2)
   complex(4), allocatable, private :: cell_translation_matrix(:, :, :, :, :, :)
   type(node_data), allocatable, private :: cell_list(:, :, :), sphere_local_interaction_list(:)
   type(translation_operator_state), target, allocatable, private :: stored_local_j_mat(:), stored_local_h_mat(:)
   data fft_local_host/0/
   data neighbor_node_model/2/
   data node_order/3/
   data cell_volume_fraction/0.2d0/

contains

   subroutine clear_fft_matrix(clear_h)
      implicit none(type, external)
      logical :: clearh
      logical, optional, intent(in) :: clear_h
      integer :: local_h_matrix_count, local_j_matrix_count

      if (present(clear_h)) then
         clearh = clear_h
      else
         clearh = .false.
      end if
      local_j_matrix_count = 0
      local_h_matrix_count = 0
      if (allocated(stored_local_j_mat)) local_j_matrix_count = size(stored_local_j_mat)
      if (allocated(stored_local_h_mat)) local_h_matrix_count = size(stored_local_h_mat)
      if (light_up) then
         write (*, '('' fft cfm 1'',2i10,l)') mstm_global_rank, local_j_matrix_count, allocated(stored_local_j_mat)
         flush (6)
      end if
      call clear_stored_translation_matrices(stored_local_j_mat)
      if (light_up) then
         write (*, '('' fft cfm 2'',2i10,l)') mstm_global_rank, local_h_matrix_count, allocated(stored_local_h_mat)
         flush (6)
      end if
      call clear_stored_translation_matrices(stored_local_h_mat)
      if (clearh) then
         if (allocated(cell_translation_matrix)) deallocate (cell_translation_matrix)
      end if
      if (allocated(sphere_node)) deallocate (sphere_node)
      if (light_up) then
         write (*, '('' fft cfm 3'',i3,l)') mstm_global_rank, allocated(cell_translation_matrix)
         flush (6)
      end if

   end subroutine clear_fft_matrix

   subroutine fft_external_to_external_expansion(neqns, nrhs, ain, gout, &
                                                 store_matrix_option, initial_run, rhs_list, &
                                                 mpi_comm, con_tran)
      implicit none
      integer :: neqns, rank, numprocs, nsphere, nrhs, mpicomm, p, nsend, &
                 i, rhs, noff, groupsize, mpigroup, syncgroup, oddnumproc
      logical :: smopt, rhslist(nrhs), contran(nrhs)
      logical, save :: firstrun, inp1, inp2
      logical, optional :: store_matrix_option, initial_run, &
                           rhs_list(nrhs), con_tran(nrhs)
      integer, save :: pgroup, pcomm, synccomm1, synccomm2, p1, p2, prank
      integer, allocatable :: grouplist(:)
      integer, optional :: mpi_comm
      complex(8) :: ain(neqns, nrhs), gout(neqns, nrhs)
      complex(8), allocatable, save :: anode(:, :, :, :, :), gnode(:, :, :, :, :), &
                                       gout_loc(:, :), ain_t(:, :), gout_t(:, :)
      data firstrun/.true./
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      if (present(store_matrix_option)) then
         smopt = store_matrix_option
      else
         smopt = .true.
      end if
      if (present(initial_run)) then
         firstrun = initial_run
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
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      gout = 0.

      if (firstrun) then
         if (numprocs .gt. 1) then
            oddnumproc = mod(numprocs, 2)
            pgroup = floor(dble(2 * rank) / dble(numprocs)) + 1
            p1 = pgroup
            p2 = p1
            call mstm_mpi(mpi_command='split', &
                          mpi_color=pgroup, mpi_key=rank, &
                          mpi_new_comm=pcomm, &
                          mpi_comm=mpicomm)
            call mstm_mpi(mpi_command='rank', mpi_rank=prank, mpi_comm=pcomm)
            call mstm_mpi(mpi_command='group', mpi_group=mpigroup, mpi_comm=mpicomm)
            groupsize = numprocs / 2 + 1
            allocate (grouplist(groupsize))
            grouplist(1) = 0
            do i = 1, groupsize - 1
               grouplist(i + 1) = i + (numprocs / 2) - 1 + oddnumproc
            end do
            inp1 = .false.
            do i = 1, groupsize
               if (rank .eq. grouplist(i)) then
                  inp1 = .true.
                  exit
               end if
            end do
            call mstm_mpi(mpi_command='incl', &
                          mpi_group=mpigroup, &
                          mpi_size=groupsize, &
                          mpi_new_group_list=grouplist, &
                          mpi_new_group=syncgroup)
            call mstm_mpi(mpi_command='create', &
                          mpi_group=syncgroup, &
                          mpi_comm=mpicomm, &
                          mpi_new_comm=synccomm1)
            deallocate (grouplist)
            groupsize = groupsize + oddnumproc
            allocate (grouplist(groupsize))
            grouplist(1) = numprocs / 2 + oddnumproc
            do i = 1, groupsize - 1
               grouplist(i + 1) = i - 1
            end do
            inp2 = .false.
            do i = 1, groupsize
               if (rank .eq. grouplist(i)) then
                  inp2 = .true.
                  exit
               end if
            end do
            call mstm_mpi(mpi_command='incl', &
                          mpi_group=mpigroup, &
                          mpi_size=groupsize, &
                          mpi_new_group_list=grouplist, &
                          mpi_new_group=syncgroup)
            call mstm_mpi(mpi_command='create', &
                          mpi_group=syncgroup, &
                          mpi_comm=mpicomm, &
                          mpi_new_comm=synccomm2)
            deallocate (grouplist)
         else
            pgroup = 1
            p1 = 1
            p2 = 2
            pcomm = mpicomm
         end if
!            call node_selection(cell_volume_fraction)
         if (light_up) then
            write (*, '('' fft1 '',i3)') mstm_global_rank
            flush (6)
         end if
         call fft_translation_initialization(host_ref_index, node_order, p1, p2)
         timedat = 0.d0
      end if
!  a test to speed up 3/23

      if (firstrun) then
         if (allocated(anode)) deallocate (anode, gnode, gout_loc, ain_t, gout_t)
         allocate (anode(cell_dim(1), cell_dim(2), cell_dim(3), node_order * (node_order + 2) * 2, nrhs), &
                   gnode(cell_dim(1), cell_dim(2), cell_dim(3), node_order * (node_order + 2) * 2, nrhs), &
                   gout_loc(number_eqns, nrhs), ain_t(neqns, nrhs), gout_t(neqns, nrhs))
      end if

      if (light_up) then
         write (*, '('' fft2 '',i3)') mstm_global_rank
         flush (6)
      end if

      do rhs = 1, nrhs
         noff = 0
         do i = 1, number_spheres
            if (host_sphere(i) .ne. fft_local_host) cycle
            if (contran(rhs)) then
               call transform_mode_coefficients(sphere_order(i), 2, -1, -1, &
                                                ain(noff + 1:noff + 2 * sphere_order(i) * (sphere_order(i) + 2), rhs), &
                                                ain_t(noff + 1:noff + 2 * sphere_order(i) * (sphere_order(i) + 2), rhs))
            else
               call transform_mode_coefficients(sphere_order(i), 2, 1, 1, &
                                                ain(noff + 1:noff + 2 * sphere_order(i) * (sphere_order(i) + 2), rhs), &
                                                ain_t(noff + 1:noff + 2 * sphere_order(i) * (sphere_order(i) + 2), rhs))
            end if
            noff = noff + 2 * sphere_order(i) * (sphere_order(i) + 2) * number_field_expansions(i)
         end do
      end do

      anode = 0.d0
      gnode = 0.d0
      gout_t = 0.d0
      gout_loc = 0.d0
      if (light_up) then
         write (*, '('' fft3 '',i3)') mstm_global_rank
         flush (6)
      end if

!timedat(1)=mstm_mpi_wtime()
      call local_sphere_to_node_translation(nrhs, ain_t, anode, &
                                            store_matrix_option=smopt, initial_run=firstrun, &
                                            mpi_comm=mpicomm, local_host=fft_local_host, sphere_to_node=.true., &
                                            merge_procs=.true.)
!timedat(1)=mstm_mpi_wtime()-timedat(1)

      if (light_up) then
         write (*, '('' fft4 '',i3)') mstm_global_rank
         flush (6)
      end if

!timedat(2)=mstm_mpi_wtime()
      do p = p1, p2
         do rhs = 1, nrhs
            call fft_node_to_node_translation(anode(:, :, :, :, rhs), &
                                              cell_translation_matrix(:, :, :, :, :, p), &
                                              gnode(:, :, :, :, rhs), p, mpi_comm=pcomm)
         end do
      end do
!timedat(2)=mstm_mpi_wtime()-timedat(2)

      call mstm_mpi(mpi_command='barrier', mpi_comm=mpicomm)
      if (numprocs .gt. 1) then
         nsend = cell_dim(1) * cell_dim(2) * cell_dim(3) * node_order * (node_order + 2)
         do rhs = 1, nrhs
            if (inp1) then
               call mstm_mpi(mpi_command='bcast', &
                             mpi_send_buf_dc=gnode(:, :, :, 1:node_order * (node_order + 2), rhs), &
                             mpi_number=nsend, &
                             mpi_rank=0, &
                             mpi_comm=synccomm1)
            end if
            if (inp2) then
               call mstm_mpi(mpi_command='bcast', &
                             mpi_send_buf_dc=gnode(:, :, :, &
                                                   node_order * (node_order + 2) + 1:node_order * (node_order + 2) * 2, rhs), &
                             mpi_number=nsend, &
                             mpi_rank=0, &
                             mpi_comm=synccomm2)
            end if
         end do
      end if
      call mstm_mpi(mpi_command='barrier', mpi_comm=mpicomm)

!         if(numprocs/2.gt.1) then
!            nsend=cell_dim(1)*cell_dim(2)*cell_dim(3)*node_order*(node_order+2)*2*nrhs
!            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=gnode, &
!                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
!         endif
      if (light_up) then
         write (*, '('' fft5 '',i3)') mstm_global_rank
         flush (6)
      end if

!timedat(3)=mstm_mpi_wtime()

      call local_sphere_to_node_translation(nrhs, gout_t, gnode, &
                                            store_matrix_option=smopt, &
                                            mpi_comm=mpicomm, local_host=fft_local_host, sphere_to_node=.false., &
                                            merge_procs=.false.)

!timedat(3)=mstm_mpi_wtime()-timedat(3)
!timedat(4)=mstm_mpi_wtime()

      call local_sphere_to_sphere_expansion(nrhs, ain_t, gout_loc, &
                                            store_matrix_option=smopt, initial_run=firstrun, &
                                            mpi_comm=mpicomm, merge_procs=.false., &
                                            local_host=fft_local_host)
!timedat(4)=mstm_mpi_wtime()-timedat(4)

      gout_t = gout_t + gout_loc

      if (light_up) then
         write (*, '('' fft6 '',i3)') mstm_global_rank
         flush (6)
      end if

      do rhs = 1, nrhs
         noff = 0
         do i = 1, number_spheres
            if (host_sphere(i) .ne. fft_local_host) cycle
            if (contran(rhs)) then
               call transform_mode_coefficients(sphere_order(i), 2, -1, -1, &
                                                gout_t(noff + 1:noff + 2 * sphere_order(i) * (sphere_order(i) + 2), rhs), &
                                                gout(noff + 1:noff + 2 * sphere_order(i) * (sphere_order(i) + 2), rhs))
            else
               call transform_mode_coefficients(sphere_order(i), 2, 1, 1, &
                                                gout_t(noff + 1:noff + 2 * sphere_order(i) * (sphere_order(i) + 2), rhs), &
                                                gout(noff + 1:noff + 2 * sphere_order(i) * (sphere_order(i) + 2), rhs))
            end if
            noff = noff + 2 * sphere_order(i) * (sphere_order(i) + 2) * number_field_expansions(i)
         end do
      end do

      if (light_up) then
         write (*, '('' fft7 '',i3)') mstm_global_rank
         flush (6)
      end if
      call mstm_mpi(mpi_command='barrier', mpi_comm=mpicomm)

!         deallocate(anode,gnode,gout_loc,ain_t,gout_t)
      firstrun = .false.

!if(mstm_global_rank.eq.mstm_global_numprocs-1.or.mstm_global_rank.eq.0.or..true.) then
!write(*,'(i5,4es13.5)') mstm_global_rank,timedat(1:4)
!endif

   end subroutine fft_external_to_external_expansion

   subroutine local_sphere_to_sphere_expansion(nrhs, ain, gout, &
                                               store_matrix_option, initial_run, rhs_list, &
                                               mpi_comm, con_tran, merge_procs, local_host)
      implicit none
      integer :: rank, numprocs, nsphere, nrhs, mpicomm, rhs, &
                 i, j, npi1, npi2, npj1, npj2, noj, noi, task, proc, &
                 ndim, idim, localhost, npairs, n, nsend
      logical :: smopt, rhslist(nrhs), contran(nrhs), mergeprocs
      logical, save :: calcmat, firstrun
      logical, optional :: store_matrix_option, initial_run, &
                           rhs_list(nrhs), con_tran(nrhs), merge_procs
      integer, optional :: mpi_comm, local_host
      complex(8) :: ain(number_eqns, nrhs), gout(number_eqns, nrhs), rimedium(2)
      type(translation_operator_state), pointer :: loc_tranmat
      type(translation_operator_state), target :: tranmat
      type(linked_ilist), pointer :: llist
      data firstrun/.true./
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      if (present(store_matrix_option)) then
         smopt = store_matrix_option
      else
         smopt = .true.
      end if
      if (present(initial_run)) then
         firstrun = initial_run
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
      if (present(local_host)) then
         localhost = local_host
      else
         localhost = 0
      end if
      if (present(merge_procs)) then
         mergeprocs = merge_procs
      else
         mergeprocs = .true.
      end if
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      gout = 0.
!
!  compute offsets for scattering coefficients
!
      if (firstrun) then
         if (smopt .and. store_translation_matrix) then
            task = 0
            ndim = 0
            do i = 1, number_spheres - 1
               if (host_sphere(i) .ne. localhost) cycle
               task = task + 1
               proc = mod(task, numprocs)
               if (proc .eq. rank) then
                  ndim = ndim + sphere_local_interaction_list(i)%number_elements
               end if
            end do
            if (allocated(stored_local_h_mat)) deallocate (stored_local_h_mat)
            allocate (stored_local_h_mat(ndim))
         end if
         calcmat = .true.
         firstrun = .false.
      else
         calcmat = (.not. (smopt .and. store_translation_matrix))
      end if

      idim = 0
      task = 0
      do i = 1, number_spheres - 1
         if (host_sphere(i) .ne. localhost) cycle
         task = task + 1
         proc = mod(task, numprocs)
         if (proc .eq. rank) then
            call exterior_refractive_index(i, rimedium)
            npairs = sphere_local_interaction_list(i)%number_elements
            llist => sphere_local_interaction_list(i)%members
            do n = 1, npairs
               idim = idim + 1
               j = llist%index
               if (n .lt. npairs) llist => llist%next
               noi = sphere_order(i)
               npi1 = sphere_offset(i) + 1
               npi2 = sphere_offset(i) + sphere_block(i)
               noj = sphere_order(j)
               npj1 = sphere_offset(j) + 1
               npj2 = sphere_offset(j) + sphere_block(j)
               if (calcmat) then
                  if (smopt .and. store_translation_matrix) then
                     call stored_local_h_mat(idim)%configure(3, sphere_position(:, i) - sphere_position(:, j), &
                                                             rimedium, min(noi, noj) .ge. translation_switch_order)
                     loc_tranmat => stored_local_h_mat(idim)
                  else
                     call tranmat%configure(3, sphere_position(:, i) - sphere_position(:, j), rimedium, &
                                            min(noi, noj) .ge. translation_switch_order)
                     loc_tranmat => tranmat
                  end if
               else
                  loc_tranmat => stored_local_h_mat(idim)
               end if
               do rhs = 1, nrhs
                  if (rhslist(rhs)) then
                     call loc_tranmat%apply(noj, 2, noi, 2, ain(npj1:npj2, rhs), gout(npi1:npi2, rhs), &
                                            shift_op=.false., tran_op=contran(rhs))
                     call loc_tranmat%apply(noi, 2, noj, 2, ain(npi1:npi2, rhs), gout(npj1:npj2, rhs), &
                                            shift_op=.true., tran_op=contran(rhs))
                  end if
               end do
               if (calcmat .and. (.not. (smopt .and. store_translation_matrix))) then
                  call tranmat%clear()
               end if
            end do
         end if
      end do

      if (numprocs .gt. 1 .and. mergeprocs) then
         nsend = number_eqns * nrhs
         call mstm_mpi(mpi_command='allreduce', mpi_recv_buf_dc=gout, &
                       mpi_number=nsend, mpi_operation=mstm_mpi_sum, mpi_comm=mpicomm)
      end if
   end subroutine local_sphere_to_sphere_expansion

   subroutine local_sphere_to_node_translation(nrhs, asphere, anode, &
                                               store_matrix_option, initial_run, rhs_list, &
                                               mpi_comm, con_tran, local_host, sphere_to_node, merge_procs)
      implicit none
      integer :: neqns, rank, numprocs, nsphere, nrhs, mpicomm, &
                 i, npi1, npi2, noi, task, proc, nsend, rhs, &
                 ndim, idim, localhost, nodei(3)
      logical :: smopt, rhslist(nrhs), contran(nrhs), spheretonode, mergeprocs
      logical, save :: calcmat, firstrun
      logical, optional :: store_matrix_option, initial_run, merge_procs, &
                           rhs_list(nrhs), con_tran(nrhs), sphere_to_node
      integer, optional :: mpi_comm, local_host
      real(8) :: rtran(3)
  complex(8) :: asphere(number_eqns, nrhs), anode(cell_dim(1), cell_dim(2), cell_dim(3), node_order * (node_order + 2) * 2, nrhs), &
                    rimedium(2), anodet(node_order * (node_order + 2) * 2)
      type(translation_operator_state), pointer :: loc_tranmat
      type(translation_operator_state), target :: tranmat
      data firstrun/.true./
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      if (present(store_matrix_option)) then
         smopt = store_matrix_option
      else
         smopt = .true.
      end if
      if (present(initial_run)) then
         firstrun = initial_run
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
      if (present(local_host)) then
         localhost = local_host
      else
         localhost = 0
      end if
      if (present(sphere_to_node)) then
         spheretonode = sphere_to_node
      else
         spheretonode = .true.
      end if
      if (present(merge_procs)) then
         mergeprocs = merge_procs
      else
         mergeprocs = .true.
      end if
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      nsphere = number_spheres
      neqns = number_eqns
!
!  compute offsets for scattering coefficients
!
      if (firstrun) then
         if (smopt .and. store_translation_matrix) then
            task = 0
            ndim = 0
            do i = 1, nsphere
               if (host_sphere(i) .eq. localhost) then
                  task = task + 1
                  proc = mod(task, numprocs)
                  if (proc .eq. rank) ndim = ndim + 1
               end if
            end do
            if (allocated(stored_local_j_mat)) deallocate (stored_local_j_mat)
            allocate (stored_local_j_mat(ndim))
         end if
         calcmat = .true.
         firstrun = .false.
      else
         calcmat = (.not. (smopt .and. store_translation_matrix))
      end if

      if (spheretonode) then
         anode = 0.d0
      else
         asphere = 0.d0
      end if

      idim = 0
      task = 0
      do i = 1, nsphere
         if (host_sphere(i) .eq. localhost) then
            task = task + 1
            proc = mod(task, numprocs)
            if (proc .eq. rank) then
               idim = idim + 1
               call exterior_refractive_index(i, rimedium)
               nodei(:) = sphere_node(:, i)
               noi = sphere_order(i)
               npi1 = sphere_offset(i) + 1
               npi2 = sphere_offset(i) + sphere_block(i)
               rtran = d_cell * (dble(nodei) - .5d0) - sphere_position(:, i) + cell_boundary(:)
               if (calcmat) then
                  if (smopt .and. store_translation_matrix) then
                     call stored_local_j_mat(idim)%configure(1, rtran, rimedium, &
                                                             min(noi, node_order) .ge. translation_switch_order)
                     loc_tranmat => stored_local_j_mat(idim)
                  else
                     call tranmat%configure(1, rtran, rimedium, &
                                            min(noi, node_order) .ge. translation_switch_order)
                     loc_tranmat => tranmat
                  end if
               else
                  loc_tranmat => stored_local_j_mat(idim)
               end if
               if (spheretonode) then
                  do rhs = 1, nrhs
                     if (rhslist(rhs)) then
                        anodet = 0.d0
                        call loc_tranmat%apply(noi, 2, node_order, 2, asphere(npi1:npi2, rhs), &
                                               anodet, shift_op=.false., tran_op=contran(rhs))
                        anode(nodei(1), nodei(2), nodei(3), :, rhs) &
                           = anode(nodei(1), nodei(2), nodei(3), :, rhs) + anodet(:)
                     end if
                  end do
               else
                  do rhs = 1, nrhs
                     if (rhslist(rhs)) then
                        anodet(:) = anode(nodei(1), nodei(2), nodei(3), :, rhs)
                        call loc_tranmat%apply(node_order, 2, noi, 2, anodet, &
                                               asphere(npi1:npi2, rhs), shift_op=.true., tran_op=contran(rhs))
                     end if
                  end do
               end if

               if (calcmat .and. (.not. (smopt .and. store_translation_matrix))) then
                  call tranmat%clear()
               end if
            end if
         end if
      end do

      if (numprocs .gt. 1 .and. mergeprocs) then
         if (spheretonode) then
            nsend = cell_dim(1) * cell_dim(2) * cell_dim(3) * node_order * (node_order + 2) * 2 * nrhs
            call mstm_mpi(mpi_command='allreduce', mpi_recv_buf_dc=anode, &
                          mpi_number=nsend, mpi_operation=mstm_mpi_sum, mpi_comm=mpicomm)
         else
            nsend = number_eqns * nrhs
            call mstm_mpi(mpi_command='allreduce', mpi_recv_buf_dc=asphere, &
                          mpi_number=nsend, mpi_operation=mstm_mpi_sum, mpi_comm=mpicomm)
         end if
      end if
   end subroutine local_sphere_to_node_translation

   subroutine node_selection(fva, target_min, target_max, d_specified, local_host)
      implicit none
      logical :: dspec
      logical, optional :: d_specified
      integer :: nsphere, m, spherenode(3), n, i, node, ix, iy, iz, ir, j, icell, ncell, cell, lochost
      integer, optional :: local_host
      real(8) :: fva, r, fv, amean, svol, tvol, dd, targetmin(3), targetmax(3)
      real(8), optional :: target_min(3), target_max(3)
      type(linked_ilist), pointer :: ilist, ilist2

      if (present(target_min)) then
         targetmin = target_min
      else
         targetmin(:) = sphere_min_position
      end if
      if (present(target_max)) then
         targetmax = target_max
      else
         targetmax(:) = sphere_max_position
      end if
      if (present(d_specified)) then
         dspec = d_specified
      else
         dspec = .false.
      end if
      if (present(local_host)) then
         lochost = local_host
      else
         lochost = 0
      end if
      fft_local_host = lochost
      cell_boundary = targetmin
      if (fft_local_host .eq. 0) then
         host_ref_index = layer_ref_index(0)
      else
         host_ref_index = sphere_ref_index(:, fft_local_host)
      end if

      amean = 0.
      nsphere = 0
      svol = 0.
      do i = 1, number_spheres
         if (host_sphere(i) .eq. fft_local_host) then
            amean = amean + sphere_radius(i)
            svol = svol + sphere_radius(i)**3
            nsphere = nsphere + 1
         end if
      end do
      svol = svol**(1.d0 / 3.d0)
      amean = amean / dble(nsphere)
      fft_number_spheres = nsphere
      tvol = 1.d0
      do i = 1, 3
         dd = targetmax(i) - targetmin(i)
         dd = max(dd, amean)
         tvol = tvol * dd
      end do

      if (.not. dspec) then
         if (fva .le. 0.d0) then
            fv = four_pi_over_three * svol**3 / tvol
            fv = min(fv, 1.d0)
            fv = max(fv, .02d0)
         else
            fv = fva
         end if
         fva = fv
         d_cell = (four_pi_over_three / fv / dble(nsphere))**(1.d0 / 3.d0) * svol
      end if

!write(*,'('' host ri:'',2e13.5)') host_ref_index
!flush(6)

      cell_dim = ceiling((targetmax(:) - targetmin(:)) / d_cell)
      do m = 1, 3
         cell_dim(m) = correctn235(cell_dim(m))
      end do
      d_cell = maxval((targetmax(:) - targetmin(:)) / dble(cell_dim(:)))

!write(*,'('' dcell:'',e13.5,3i5)') d_cell,cell_dim
!flush(6)
      if (allocated(cell_list)) deallocate (cell_list)
      if (allocated(sphere_node)) deallocate (sphere_node)
      allocate (cell_list(cell_dim(1), cell_dim(2), cell_dim(3)), &
                sphere_node(3, number_spheres))
      cell_list(:, :, :)%number_elements = 0
      do iz = 1, cell_dim(3)
         do iy = 1, cell_dim(2)
            do ix = 1, cell_dim(1)
               allocate (cell_list(ix, iy, iz)%members)
            end do
         end do
      end do

!write(*,'('' d,n:'',e13.5,3i8)') d_cell,cell_dim
!flush(6)

      cell_origin(:) = d_cell * dble(cell_dim(:)) / 2.d0
      do i = 1, number_spheres
         if (host_sphere(i) .ne. fft_local_host) cycle
         do m = 1, 3
            r = sphere_position(m, i) - targetmin(m)
            node = floor(r / d_cell) + 1
            spherenode(m) = node
            spherenode(m) = min(spherenode(m), cell_dim(m))
            spherenode(m) = max(spherenode(m), 1)
         end do
         sphere_node(:, i) = spherenode
         n = cell_list(spherenode(1), spherenode(2), spherenode(3))%number_elements
!write(*,'('' i:'',5i5)') i,spherenode,n
!flush(6)
         ilist => cell_list(spherenode(1), spherenode(2), spherenode(3))%members
         do m = 1, n
            if (m .eq. n) allocate (ilist%next)
            ilist => ilist%next
         end do
         ilist%index = i
         ilist%next => null()
         cell_list(spherenode(1), spherenode(2), spherenode(3))%number_elements = n + 1
      end do

      n = -1
      do iz = -1, 1
         do iy = -1, 1
            do ix = -1, 1
               ir = ix * ix + iy * iy + iz * iz
               if (neighbor_node_model .eq. 0) then
                  if (ir .gt. 0) cycle
               elseif (neighbor_node_model .eq. 1) then
                  if (ir .gt. 1) cycle
               elseif (neighbor_node_model .eq. 2) then
                  if (ir .gt. 2) cycle
               end if
               n = n + 1
               neighbor_node(:, n) = (/ix, iy, iz/)
            end do
         end do
      end do
      number_neighbor_nodes = n

!write(*,'('' nn:'',i8)') number_neighbor_nodes
!flush(6)

      if (allocated(sphere_local_interaction_list)) deallocate (sphere_local_interaction_list)
      allocate (sphere_local_interaction_list(number_spheres))
      sphere_local_interaction_list(:)%number_elements = 0

      do i = 1, number_spheres
         if (host_sphere(i) .ne. fft_local_host) cycle
         allocate (sphere_local_interaction_list(i)%members)
         ilist2 => sphere_local_interaction_list(i)%members
         n = 0
         do cell = 0, number_neighbor_nodes
            spherenode = sphere_node(:, i) - neighbor_node(:, cell)
            if (any(spherenode .lt. (/1, 1, 1/)) .or. any(spherenode .gt. cell_dim)) cycle
            ncell = cell_list(spherenode(1), spherenode(2), spherenode(3))%number_elements
            ilist => cell_list(spherenode(1), spherenode(2), spherenode(3))%members
            do icell = 1, ncell
               j = ilist%index
               if (j .gt. i) then
                  ilist2%index = j
                  n = n + 1
                  allocate (ilist2%next)
                  ilist2 => ilist2%next
               end if
               if (icell .lt. ncell) ilist => ilist%next
            end do
         end do
         sphere_local_interaction_list(i)%number_elements = n
      end do

   end subroutine node_selection

   subroutine fft_node_to_node_translation(acoef, tranmat, gcoef, pmode, mpi_comm, tran_op)
      implicit none
      logical :: tranop
      logical, optional :: tran_op
      integer :: nblk, celldim2(3), n, l, mpicomm, numprocs, rank, ncells, task, proc, nsend, pmode, ncells8
      integer, optional :: mpi_comm
    complex(4) :: tranmat(8 * cell_dim(1) * cell_dim(2) * cell_dim(3), node_order * (node_order + 2), node_order * (node_order + 2))
      complex(8) :: acoef(cell_dim(1) * cell_dim(2) * cell_dim(3), node_order * (node_order + 2), 2), &
                    aft(8 * cell_dim(1) * cell_dim(2) * cell_dim(3)), &
                    gcoef(cell_dim(1) * cell_dim(2) * cell_dim(3), node_order * (node_order + 2), 2), &
                    gft(8 * cell_dim(1) * cell_dim(2) * cell_dim(3), node_order * (node_order + 2))

      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      if (present(tran_op)) then
         tranop = tran_op
      else
         tranop = .false.
      end if
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      nblk = node_order * (node_order + 2)
      celldim2 = 2 * cell_dim
      ncells = cell_dim(1) * cell_dim(2) * cell_dim(3)
      ncells8 = ncells * 8
      gft = 0.d0
      task = 0

      do n = 1, nblk
         task = task + 1
         proc = mod(task, numprocs)
         if (proc .eq. rank) then
            aft = 0.d0
            call fftmtx(acoef(1:ncells, n, pmode), aft(1:ncells8), 1, cell_dim, celldim2, 1)
            if (tranop) then
               do l = 1, nblk
                  gft(1:ncells8, l) = gft(1:ncells8, l) + tranmat(1:ncells8, n, l) * aft(1:ncells8)
               end do
            else
               do l = 1, nblk
                  gft(1:ncells8, l) = gft(1:ncells8, l) + tranmat(1:ncells8, l, n) * aft(1:ncells8)
               end do
            end if
         end if
      end do

      if (numprocs .gt. 1) then
         nsend = nblk * ncells8
         call mstm_mpi(mpi_command='allreduce', &
                       mpi_recv_buf_dc=gft(1:ncells8, 1:nblk), &
                       mpi_number=nsend, mpi_operation=mstm_mpi_sum, mpi_comm=mpicomm)
      end if

      task = 0
      do n = 1, nblk
         task = task + 1
         proc = mod(task, numprocs)
         if (proc .eq. rank) then
            gcoef(1:ncells, n, pmode) = 0.d0
            call fftmtx(gcoef(1:ncells, n, pmode), gft(1:ncells8, n), 1, cell_dim, celldim2, -1)
         end if
      end do
      if (numprocs .gt. 1) then
         nsend = nblk * ncells
         call mstm_mpi(mpi_command='allreduce', &
                       mpi_recv_buf_dc=gcoef(1:ncells, 1:nblk, pmode), &
                       mpi_number=nsend, mpi_operation=mstm_mpi_sum, mpi_comm=mpicomm)
      end if
      gcoef(1:ncells, 1:nblk, pmode) = gcoef(1:ncells, 1:nblk, pmode) / dble(ncells8)
   end subroutine fft_node_to_node_translation

   subroutine fft_translation_initialization(ri, nodr, p1, p2)
      implicit none
      logical :: inhole
      integer :: nodr, nx, ny, nz, is, nblk, ncells2(3), isx, isy, isz, nx1, ny1, nz1, &
                 n, i, node(3), p1, p2, p, l
      real(8) :: x, y, z, xp(3), r
      complex(8) :: hij(nodr * (nodr + 2), nodr * (nodr + 2), 2), ri(2), &
                    htemp(1:2 * cell_dim(1), 1:2 * cell_dim(2), 1:2 * cell_dim(3))

      nblk = nodr * (nodr + 2)
      ncells2 = 2 * cell_dim

      if (light_up) then
         write (*, '('' fft ti '',i3,l)') mstm_global_rank, allocated(cell_translation_matrix)
         flush (6)
      end if
      if (allocated(cell_translation_matrix)) return
      if (allocated(cell_translation_matrix)) deallocate (cell_translation_matrix)
      allocate (cell_translation_matrix(0:2 * cell_dim(1) - 1, &
                                        0:2 * cell_dim(2) - 1, 0:2 * cell_dim(3) - 1, 1:nblk, 1:nblk, p1:p2))
      cell_translation_matrix = 0.d0

      do nz = 0, cell_dim(3)
         z = d_cell * dble(nz)
         do ny = 0, cell_dim(2)
            y = d_cell * dble(ny)
            do nx = 0, cell_dim(1)
               x = d_cell * dble(nx)
               inhole = .false.
               do i = 0, number_neighbor_nodes
                  node(:) = (/nx, ny, nz/) - neighbor_node(:, i)
                  if (all(node .eq. (/0, 0, 0/))) then
                     inhole = .true.
                     exit
                  end if
               end do
               if (inhole) cycle

               r = sqrt(x * x + y * y + z * z)
               do is = 0, 7
                  isx = 1 - 2 * mod(is, 2)
                  isy = 1 - 2 * mod(int(is / 2), 2)
                  isz = 1 - 2 * mod(int(is / 4), 2)
                  nx1 = isx * nx - (isx - 1) * cell_dim(1)
                  ny1 = isy * ny - (isy - 1) * cell_dim(2)
                  nz1 = isz * nz - (isz - 1) * cell_dim(3)
                  if (nx1 .lt. ncells2(1) .and. ny1 .lt. ncells2(2) &
                      .and. nz1 .lt. ncells2(3)) then
                     xp(1) = isx * x
                     xp(2) = isy * y
                     xp(3) = isz * z
                     call generate_translation_matrix(nodr, nodr, translation_vector=xp, &
                                                      refractive_index=ri, ac_matrix=hij, vswf_type=3, &
                                                      mode_s=2, mode_t=2)
                     cell_translation_matrix(nx1, ny1, nz1, 1:nblk, 1:nblk, p1:p2) &
                        = cell_translation_matrix(nx1, ny1, nz1, 1:nblk, 1:nblk, p1:p2) &
                          + hij(1:nblk, 1:nblk, p1:p2)
                  end if
               end do
            end do
         end do
      end do
      n = 2 * nblk * nblk
      do p = p1, p2
         do n = 1, nblk
            do l = 1, nblk
               htemp(:, :, :) = cell_translation_matrix(:, :, :, l, n, p)
               call fftmtx(htemp, htemp, 1, ncells2, ncells2, 1)
               cell_translation_matrix(:, :, :, l, n, p) = htemp(:, :, :)
            end do
         end do
      end do

   end subroutine fft_translation_initialization

   subroutine fft1don3d(ain, aout, nblk, ntot1, ntot2, ntot3in, ntot3out, ndimin, &
                        ndimout, isign, looporder, trig)
      implicit none
      integer :: nblk, ntot1, ntot2, ntot3in, ntot3out, l, m, n, isign, &
                 looporder(3), &
                 triplet(3), ndimin(3), ndimout(3), ntot3, i1, i2, i3
      integer, parameter :: mxtrig = 1000
      real(8) :: trig(mxtrig), ar_temp(nblk, max(ntot3in, ntot3out)), &
                 ai_temp(nblk, max(ntot3in, ntot3out))
      complex(8) :: aout(nblk, ndimout(1), ndimout(2), ndimout(3)), &
                    ain(nblk, ndimin(1), ndimin(2), ndimin(3))
      i1 = looporder(1)
      i2 = looporder(2)
      i3 = looporder(3)
      ntot3 = max(ntot3in, ntot3out)
      do l = 1, ntot1
         triplet(i1) = l
         do m = 1, ntot2
            triplet(i2) = m
            do n = 1, ntot3in
               triplet(i3) = n
               ar_temp(1:nblk, n) = dble(ain(1:nblk, triplet(1), triplet(2), triplet(3)))
               ai_temp(1:nblk, n) = aimag(ain(1:nblk, triplet(1), triplet(2), triplet(3)))
            end do
            do n = ntot3in + 1, ntot3out
               ar_temp(1:nblk, n) = 0.d0
               ai_temp(1:nblk, n) = 0.d0
            end do
            call cgpfa(ar_temp(1:nblk, 1:ntot3), ai_temp(1:nblk, 1:ntot3), &
                       trig, nblk, ntot3, isign)
            do n = 1, ntot3out
               triplet(i3) = n
               aout(1:nblk, triplet(1), triplet(2), triplet(3)) &
                  = cmplx(ar_temp(1:nblk, n), ai_temp(1:nblk, n), kind=kind(0.0d0))
            end do
         end do
      end do
   end subroutine fft1don3d

   subroutine fftmtx(am, amf, nblk, ntot, ntot2, isign)
      implicit none
      integer :: nblk, ntot(3), ntot2(3), isign, &
                 ntotxold, ntotyold, ntotzold, &
                 nblkold
      integer, parameter :: mxtrig = 1000
      real(8) :: trig(mxtrig, 3)
      save :: trig, ntotxold, ntotyold, ntotzold, nblkold
      complex(8) :: amf(nblk, ntot2(1), ntot2(2), ntot2(3)), &
                    am(nblk, ntot(1), ntot(2), ntot(3))
      data ntotxold, ntotyold, ntotzold, nblkold/0, 0, 0, 0/
      if (ntot2(1) .ne. ntotxold .or. ntot2(2) .ne. ntotyold &
          .or. ntot2(3) .ne. ntotzold) then
         ntotxold = ntot2(1)
         ntotyold = ntot2(2)
         ntotzold = ntot2(3)
         call setgpfa(trig(:, 1), ntot2(1))
         call setgpfa(trig(:, 2), ntot2(2))
         call setgpfa(trig(:, 3), ntot2(3))
      end if
      if (isign .eq. 1) then
         call fft1don3d(am, amf, nblk, ntot(1), ntot(2), ntot(3), ntot2(3), ntot, &
                        ntot2, 1, (/1, 2, 3/), trig(:, 3))
         call fft1don3d(amf, amf, nblk, ntot2(3), ntot(1), ntot(2), ntot2(2), ntot2, &
                        ntot2, 1, (/3, 1, 2/), trig(:, 2))
         call fft1don3d(amf, amf, nblk, ntot2(3), ntot2(2), ntot(1), ntot2(1), ntot2, &
                        ntot2, 1, (/3, 2, 1/), trig(:, 1))
      else
         call fft1don3d(amf, amf, nblk, ntot2(1), ntot2(2), ntot2(3), ntot(3), ntot2, &
                        ntot2, -1, (/1, 2, 3/), trig(:, 3))
         call fft1don3d(amf, amf, nblk, ntot(3), ntot2(1), ntot2(2), ntot(2), ntot2, &
                        ntot2, -1, (/3, 1, 2/), trig(:, 2))
         call fft1don3d(amf, am, nblk, ntot(3), ntot(2), ntot2(1), ntot(1), ntot2, &
                        ntot, -1, (/3, 2, 1/), trig(:, 1))
      end if
   end subroutine fftmtx

   pure logical function checkn235(n)
      implicit none
      integer, intent(in) :: n
      integer :: nn, ifac, ll, kk
      nn = n
      ifac = 2
      do ll = 1, 3
         kk = 0
         do while (mod(nn, ifac) .eq. 0)
            kk = kk + 1
            nn = nn / ifac
         end do
         ifac = ifac + ll
      end do
      checkn235 = nn .eq. 1
   end function checkn235

   pure integer function correctn235(n)
      implicit none
      integer, intent(in) :: n
      integer :: n1
      n1 = n
      do while (.not. checkn235(n1))
         n1 = n1 + 1
      end do
      correctn235 = n1
   end function correctn235

end module fft_translation
