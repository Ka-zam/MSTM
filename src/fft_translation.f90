module fft_translation
   use, intrinsic :: iso_fortran_env, only: real32, real64
   use constants
   use gpfa_controller, only: cgpfa
   use gpfa_setup, only: setgpfa
   use parallel_runtime, only: mpi_comm_world, mstm_global_rank, parallel_allreduce_sum, parallel_barrier, &
                               parallel_broadcast, parallel_communicator_create, parallel_group, &
                               parallel_group_include, parallel_rank, parallel_size, parallel_split, parallel_wall_time
   use numerical_tables
   use angular_functions, only: generate_translation_matrix
   use sphere_data
   use translation_expansions, only: clear_stored_translation_matrices
   use translation_operator, only: transform_mode_coefficients, translation_operator_state
   use mie
   implicit none
   type node_data
      integer :: number_elements = 0
      type(linked_ilist), pointer :: members => null()
   end type node_data
   type linked_ilist
      integer :: index
      type(linked_ilist), pointer :: next => null()
   end type linked_ilist
   type coefficient_list
      complex(real64), pointer :: coefficient_vector(:, :)
   end type coefficient_list

   type, public :: fft_translation_metrics_t
      real(real64) :: initialization_time = 0.0_real64
      real(real64) :: sphere_to_node_time = 0.0_real64
      real(real64) :: node_to_node_time = 0.0_real64
      real(real64) :: node_to_sphere_time = 0.0_real64
      real(real64) :: local_interaction_time = 0.0_real64
      real(real64) :: transform_time = 0.0_real64
      integer :: interaction_calls = 0
      integer :: transform_calls = 0
   end type fft_translation_metrics_t

   type, public :: fft_translation_plan_t
      private
      integer :: cell_dim(3) = 0
      integer :: number_neighbor_nodes = 0
      integer :: neighbor_node(3, 0:26) = 0
      integer :: fft_local_host = 0
      integer :: fft_number_spheres = 0
      integer :: node_order = 3
      integer :: neighbor_node_model = 2
      integer, allocatable :: sphere_node(:, :)
      real(real64) :: cell_origin(3) = 0.0_real64
      real(real64) :: cell_boundary(3) = 0.0_real64
      real(real64) :: cell_volume_fraction = 0.2_real64
      real(real64) :: d_cell = 0.0_real64
      complex(real64) :: host_ref_index(2) = (0.0_real64, 0.0_real64)
      complex(real32), allocatable :: cell_translation_matrix(:, :, :, :, :, :)
      type(node_data), allocatable :: cell_list(:, :, :), sphere_local_interaction_list(:)
      type(translation_operator_state), allocatable :: stored_local_j_mat(:), stored_local_h_mat(:)
      complex(real64), allocatable :: anode(:, :, :, :, :), gnode(:, :, :, :, :)
      complex(real64), allocatable :: gout_local(:, :), input_work(:, :), output_work(:, :)
      logical :: first_application = .true.
      logical :: first_local_interaction = .true.
      logical :: first_node_translation = .true.
      logical :: calculate_local_matrix = .true.
      logical :: calculate_node_matrix = .true.
      logical :: in_polarization_group_1 = .false., in_polarization_group_2 = .false.
      integer :: polarization_group = 1, polarization_communicator = 0
      integer :: synchronization_communicator_1 = 0
      integer :: synchronization_communicator_2 = 0
      integer :: first_polarization = 1, last_polarization = 2, polarization_rank = 0
      real(real64) :: transform_trigonometry(1000, 3) = 0.0_real64
      integer :: previous_transform_size(3) = 0, previous_transform_batch = 0
      type(fft_translation_metrics_t) :: metrics
   contains
      procedure, public :: clear => clear_fft_matrix
      procedure, public :: configure => configure_fft_nodes
      procedure, public :: apply => fft_external_to_external_expansion
      procedure, public :: configuration => get_fft_configuration
      procedure, public :: performance_metrics => get_fft_performance_metrics
      final :: finalize_fft_plan
   end type fft_translation_plan_t

   type(fft_translation_plan_t), target, public :: fft_plan

contains

   subroutine clear_fft_matrix(self, clear_h)
      implicit none(type, external)
      class(fft_translation_plan_t), intent(inout) :: self
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
      if (allocated(self%stored_local_j_mat)) local_j_matrix_count = size(self%stored_local_j_mat)
      if (allocated(self%stored_local_h_mat)) local_h_matrix_count = size(self%stored_local_h_mat)
      if (light_up) then
         write (*, '('' fft cfm 1'',2i10,l)') mstm_global_rank, local_j_matrix_count, allocated(self%stored_local_j_mat)
         flush (6)
      end if
      call clear_stored_translation_matrices(self%stored_local_j_mat)
      if (light_up) then
         write (*, '('' fft cfm 2'',2i10,l)') mstm_global_rank, local_h_matrix_count, allocated(self%stored_local_h_mat)
         flush (6)
      end if
      call clear_stored_translation_matrices(self%stored_local_h_mat)
      if (clearh) then
         if (allocated(self%cell_translation_matrix)) deallocate (self%cell_translation_matrix)
      end if
      call clear_fft_geometry(self)
      if (allocated(self%anode)) deallocate (self%anode, self%gnode, self%gout_local, &
                                             self%input_work, self%output_work)
      self%first_application = .true.
      self%first_local_interaction = .true.
      self%first_node_translation = .true.
      self%calculate_local_matrix = .true.
      self%calculate_node_matrix = .true.
      if (light_up) then
         write (*, '('' fft cfm 3'',i3,l)') mstm_global_rank, allocated(self%cell_translation_matrix)
         flush (6)
      end if

   end subroutine clear_fft_matrix

   subroutine finalize_fft_plan(self)
      type(fft_translation_plan_t), intent(inout) :: self

      call self%clear(clear_h=.true.)
   end subroutine finalize_fft_plan

   subroutine clear_fft_geometry(self)
      class(fft_translation_plan_t), intent(inout) :: self
      integer :: ix, iy, iz, sphere

      if (allocated(self%cell_list)) then
         do iz = 1, size(self%cell_list, 3)
            do iy = 1, size(self%cell_list, 2)
               do ix = 1, size(self%cell_list, 1)
                  call clear_linked_indices(self%cell_list(ix, iy, iz)%members)
               end do
            end do
         end do
         deallocate (self%cell_list)
      end if
      if (allocated(self%sphere_local_interaction_list)) then
         do sphere = 1, size(self%sphere_local_interaction_list)
            call clear_linked_indices(self%sphere_local_interaction_list(sphere)%members)
         end do
         deallocate (self%sphere_local_interaction_list)
      end if
      if (allocated(self%sphere_node)) deallocate (self%sphere_node)
   end subroutine clear_fft_geometry

   subroutine clear_linked_indices(head)
      type(linked_ilist), pointer, intent(inout) :: head
      type(linked_ilist), pointer :: current, next

      current => head
      do while (associated(current))
         next => current%next
         deallocate (current)
         current => next
      end do
      nullify (head)
   end subroutine clear_linked_indices

   subroutine get_fft_configuration(self, cell_dimensions, volume_fraction, cell_size, &
                                    neighbor_count, expansion_order)
      class(fft_translation_plan_t), intent(in) :: self
      integer, intent(out) :: cell_dimensions(3), neighbor_count, expansion_order
      real(real64), intent(out) :: volume_fraction, cell_size

      cell_dimensions = self%cell_dim
      volume_fraction = self%cell_volume_fraction
      cell_size = self%d_cell
      neighbor_count = self%number_neighbor_nodes
      expansion_order = self%node_order
   end subroutine get_fft_configuration

   pure function get_fft_performance_metrics(self) result(metrics)
      class(fft_translation_plan_t), intent(in) :: self
      type(fft_translation_metrics_t) :: metrics

      metrics = self%metrics
   end function get_fft_performance_metrics

   subroutine fft_external_to_external_expansion(self, neqns, nrhs, ain, gout, &
                                                 store_matrix_option, initial_run, rhs_list, &
                                                 mpi_comm, con_tran)
      implicit none(type, external)
      class(fft_translation_plan_t), intent(inout), target :: self
      integer :: neqns, rank, numprocs, nsphere, nrhs, mpicomm, p, nsend, &
                 i, rhs, noff, groupsize, mpigroup, syncgroup, oddnumproc
      real(real64) :: phase_start
      logical :: smopt, rhslist(nrhs), contran(nrhs), resize_work_arrays
      logical, optional :: store_matrix_option, initial_run, &
                           rhs_list(nrhs), con_tran(nrhs)
      integer, allocatable :: grouplist(:)
      integer, optional :: mpi_comm
      complex(real64) :: ain(neqns, nrhs), gout(neqns, nrhs)
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
         self%first_application = initial_run
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
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      gout = 0.

      if (self%first_application) then
         self%metrics%initialization_time = 0.0_real64
         self%metrics%sphere_to_node_time = 0.0_real64
         self%metrics%node_to_node_time = 0.0_real64
         self%metrics%node_to_sphere_time = 0.0_real64
         self%metrics%local_interaction_time = 0.0_real64
         self%metrics%transform_time = 0.0_real64
         self%metrics%interaction_calls = 0
         self%metrics%transform_calls = 0
         if (numprocs .gt. 1) then
            oddnumproc = mod(numprocs, 2)
            self%polarization_group = floor(dble(2 * rank) / dble(numprocs)) + 1
            self%first_polarization = self%polarization_group
            self%last_polarization = self%first_polarization
            call parallel_split( &
               mpi_color=self%polarization_group, mpi_key=rank, &
               mpi_new_comm=self%polarization_communicator, &
               mpi_comm=mpicomm)
            call parallel_rank(mpi_rank=self%polarization_rank, mpi_comm=self%polarization_communicator)
            call parallel_group(mpi_group=mpigroup, mpi_comm=mpicomm)
            groupsize = numprocs / 2 + 1
            allocate (grouplist(groupsize))
            grouplist(1) = 0
            do i = 1, groupsize - 1
               grouplist(i + 1) = i + (numprocs / 2) - 1 + oddnumproc
            end do
            self%in_polarization_group_1 = .false.
            do i = 1, groupsize
               if (rank .eq. grouplist(i)) then
                  self%in_polarization_group_1 = .true.
                  exit
               end if
            end do
            call parallel_group_include( &
               mpi_group=mpigroup, &
               mpi_size=groupsize, &
               mpi_new_group_list=grouplist, &
               mpi_new_group=syncgroup)
            call parallel_communicator_create( &
               mpi_group=syncgroup, &
               mpi_comm=mpicomm, &
               mpi_new_comm=self%synchronization_communicator_1)
            deallocate (grouplist)
            groupsize = groupsize + oddnumproc
            allocate (grouplist(groupsize))
            grouplist(1) = numprocs / 2 + oddnumproc
            do i = 1, groupsize - 1
               grouplist(i + 1) = i - 1
            end do
            self%in_polarization_group_2 = .false.
            do i = 1, groupsize
               if (rank .eq. grouplist(i)) then
                  self%in_polarization_group_2 = .true.
                  exit
               end if
            end do
            call parallel_group_include( &
               mpi_group=mpigroup, &
               mpi_size=groupsize, &
               mpi_new_group_list=grouplist, &
               mpi_new_group=syncgroup)
            call parallel_communicator_create( &
               mpi_group=syncgroup, &
               mpi_comm=mpicomm, &
               mpi_new_comm=self%synchronization_communicator_2)
            deallocate (grouplist)
         else
            self%polarization_group = 1
            self%first_polarization = 1
            self%last_polarization = 2
            self%polarization_communicator = mpicomm
         end if
!            call configure_fft_nodes(self%cell_volume_fraction)
         if (light_up) then
            write (*, '('' fft1 '',i3)') mstm_global_rank
            flush (6)
         end if
         phase_start = parallel_wall_time()
 call initialize_fft_translation_matrix(self, self%host_ref_index, self%node_order, self%first_polarization, self%last_polarization)
         self%metrics%initialization_time = parallel_wall_time() - phase_start
      end if
      self%metrics%interaction_calls = self%metrics%interaction_calls + 1
!  a test to speed up 3/23

      resize_work_arrays = .not. allocated(self%anode)
      if (allocated(self%anode)) then
         resize_work_arrays = size(self%anode, 1) /= self%cell_dim(1) .or. &
                              size(self%anode, 2) /= self%cell_dim(2) .or. &
                              size(self%anode, 3) /= self%cell_dim(3) .or. &
                              size(self%anode, 4) /= self%node_order * (self%node_order + 2) * 2 .or. &
                              size(self%anode, 5) /= nrhs .or. size(self%gout_local, 1) /= sphere_cluster%number_eqns
      end if
      if (resize_work_arrays) then
         if (allocated(self%anode)) deallocate (self%anode, self%gnode, self%gout_local, self%input_work, self%output_work)
    allocate (self%anode(self%cell_dim(1), self%cell_dim(2), self%cell_dim(3), self%node_order * (self%node_order + 2) * 2, nrhs), &
              self%gnode(self%cell_dim(1), self%cell_dim(2), self%cell_dim(3), self%node_order * (self%node_order + 2) * 2, nrhs), &
                   self%gout_local(sphere_cluster%number_eqns, nrhs), self%input_work(neqns, nrhs), self%output_work(neqns, nrhs))
      end if

      if (light_up) then
         write (*, '('' fft2 '',i3)') mstm_global_rank
         flush (6)
      end if

      do rhs = 1, nrhs
         noff = 0
         do i = 1, sphere_cluster%number_spheres
            if (sphere_cluster%host_sphere(i) .ne. self%fft_local_host) cycle
            if (contran(rhs)) then
               call transform_mode_coefficients(sphere_cluster%sphere_order(i), 2, -1, -1, &
                              ain(noff + 1:noff + 2 * sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2), rhs), &
                    self%input_work(noff + 1:noff + 2 * sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2), rhs))
            else
               call transform_mode_coefficients(sphere_cluster%sphere_order(i), 2, 1, 1, &
                              ain(noff + 1:noff + 2 * sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2), rhs), &
                    self%input_work(noff + 1:noff + 2 * sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2), rhs))
            end if
 noff = noff + 2 * sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2) * sphere_cluster%number_field_expansions(i)
         end do
      end do

      self%anode = 0.d0
      self%gnode = 0.d0
      self%output_work = 0.d0
      self%gout_local = 0.d0
      if (light_up) then
         write (*, '('' fft3 '',i3)') mstm_global_rank
         flush (6)
      end if

      phase_start = parallel_wall_time()
      call local_sphere_to_node_translation(self, nrhs, self%input_work, self%anode, &
                                            store_matrix_option=smopt, initial_run=self%first_application, &
                                            mpi_comm=mpicomm, local_host=self%fft_local_host, sphere_to_node=.true., &
                                            merge_procs=.true.)
      self%metrics%sphere_to_node_time = self%metrics%sphere_to_node_time + parallel_wall_time() - phase_start

      if (light_up) then
         write (*, '('' fft4 '',i3)') mstm_global_rank
         flush (6)
      end if

      phase_start = parallel_wall_time()
      do p = self%first_polarization, self%last_polarization
         do rhs = 1, nrhs
            call fft_node_to_node_translation(self, self%anode(:, :, :, :, rhs), &
                                              self%cell_translation_matrix(:, :, :, :, :, p), &
                                              self%gnode(:, :, :, :, rhs), p, mpi_comm=self%polarization_communicator)
         end do
      end do
      self%metrics%node_to_node_time = self%metrics%node_to_node_time + parallel_wall_time() - phase_start

      call parallel_barrier(mpi_comm=mpicomm)
      if (numprocs .gt. 1) then
         nsend = self%cell_dim(1) * self%cell_dim(2) * self%cell_dim(3) * self%node_order * (self%node_order + 2)
         do rhs = 1, nrhs
            if (self%in_polarization_group_1) then
               call parallel_broadcast( &
                  send_buffer=self%gnode(:, :, :, 1:self%node_order * (self%node_order + 2), rhs), &
                  mpi_number=nsend, &
                  mpi_rank=0, &
                  mpi_comm=self%synchronization_communicator_1)
            end if
            if (self%in_polarization_group_2) then
               call parallel_broadcast( &
                  send_buffer=self%gnode(:, :, :, &
                                    self%node_order * (self%node_order + 2) + 1:self%node_order * (self%node_order + 2) * 2, rhs), &
                  mpi_number=nsend, &
                  mpi_rank=0, &
                  mpi_comm=self%synchronization_communicator_2)
            end if
         end do
      end if
      call parallel_barrier(mpi_comm=mpicomm)

!         if(numprocs/2.gt.1) then
!            nsend=self%cell_dim(1)*self%cell_dim(2)*self%cell_dim(3)*self%node_order*(self%node_order+2)*2*nrhs
!            call parallel_allreduce_sum(receive_buffer=self%gnode, &
!                 mpi_number=nsend,mpi_comm=mpicomm)
!         endif
      if (light_up) then
         write (*, '('' fft5 '',i3)') mstm_global_rank
         flush (6)
      end if

      phase_start = parallel_wall_time()
      call local_sphere_to_node_translation(self, nrhs, self%output_work, self%gnode, &
                                            store_matrix_option=smopt, &
                                            mpi_comm=mpicomm, local_host=self%fft_local_host, sphere_to_node=.false., &
                                            merge_procs=.false.)

      self%metrics%node_to_sphere_time = self%metrics%node_to_sphere_time + parallel_wall_time() - phase_start
      phase_start = parallel_wall_time()
      call local_sphere_to_sphere_expansion(self, nrhs, self%input_work, self%gout_local, &
                                            store_matrix_option=smopt, initial_run=self%first_application, &
                                            mpi_comm=mpicomm, merge_procs=.false., &
                                            local_host=self%fft_local_host)
      self%metrics%local_interaction_time = self%metrics%local_interaction_time + parallel_wall_time() - phase_start

      self%output_work = self%output_work + self%gout_local

      if (light_up) then
         write (*, '('' fft6 '',i3)') mstm_global_rank
         flush (6)
      end if

      do rhs = 1, nrhs
         noff = 0
         do i = 1, sphere_cluster%number_spheres
            if (sphere_cluster%host_sphere(i) .ne. self%fft_local_host) cycle
            if (contran(rhs)) then
               call transform_mode_coefficients(sphere_cluster%sphere_order(i), 2, -1, -1, &
                 self%output_work(noff + 1:noff + 2 * sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2), rhs), &
                               gout(noff + 1:noff + 2 * sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2), rhs))
            else
               call transform_mode_coefficients(sphere_cluster%sphere_order(i), 2, 1, 1, &
                 self%output_work(noff + 1:noff + 2 * sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2), rhs), &
                               gout(noff + 1:noff + 2 * sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2), rhs))
            end if
 noff = noff + 2 * sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2) * sphere_cluster%number_field_expansions(i)
         end do
      end do

      if (light_up) then
         write (*, '('' fft7 '',i3)') mstm_global_rank
         flush (6)
      end if
      call parallel_barrier(mpi_comm=mpicomm)

!         deallocate(self%anode,self%gnode,self%gout_local,self%input_work,self%output_work)
      self%first_application = .false.

   end subroutine fft_external_to_external_expansion

   subroutine local_sphere_to_sphere_expansion(self, nrhs, ain, gout, &
                                               store_matrix_option, initial_run, rhs_list, &
                                               mpi_comm, con_tran, merge_procs, local_host)
      implicit none(type, external)
      class(fft_translation_plan_t), intent(inout), target :: self
      integer :: rank, numprocs, nsphere, nrhs, mpicomm, rhs, &
                 i, j, npi1, npi2, npj1, npj2, noj, noi, task, proc, &
                 ndim, idim, localhost, npairs, n, nsend
      logical :: smopt, rhslist(nrhs), contran(nrhs), mergeprocs
      logical, optional :: store_matrix_option, initial_run, &
                           rhs_list(nrhs), con_tran(nrhs), merge_procs
      integer, optional :: mpi_comm, local_host
      complex(real64) :: ain(sphere_cluster%number_eqns, nrhs), gout(sphere_cluster%number_eqns, nrhs), rimedium(2)
      type(translation_operator_state), pointer :: loc_tranmat
      type(translation_operator_state), target :: tranmat
      type(linked_ilist), pointer :: llist
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
         self%first_local_interaction = initial_run
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
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      gout = 0.
!
!  compute offsets for scattering coefficients
!
      if (self%first_local_interaction) then
         if (smopt .and. sphere_cluster%store_translation_matrix) then
            task = 0
            ndim = 0
            do i = 1, sphere_cluster%number_spheres - 1
               if (sphere_cluster%host_sphere(i) .ne. localhost) cycle
               task = task + 1
               proc = mod(task, numprocs)
               if (proc .eq. rank) then
                  ndim = ndim + self%sphere_local_interaction_list(i)%number_elements
               end if
            end do
            if (allocated(self%stored_local_h_mat)) deallocate (self%stored_local_h_mat)
            allocate (self%stored_local_h_mat(ndim))
         end if
         self%calculate_local_matrix = .true.
         self%first_local_interaction = .false.
      else
         self%calculate_local_matrix = (.not. (smopt .and. sphere_cluster%store_translation_matrix))
      end if

      idim = 0
      task = 0
      do i = 1, sphere_cluster%number_spheres - 1
         if (sphere_cluster%host_sphere(i) .ne. localhost) cycle
         task = task + 1
         proc = mod(task, numprocs)
         if (proc .eq. rank) then
            call exterior_refractive_index(i, rimedium)
            npairs = self%sphere_local_interaction_list(i)%number_elements
            llist => self%sphere_local_interaction_list(i)%members
            do n = 1, npairs
               idim = idim + 1
               j = llist%index
               if (n .lt. npairs) llist => llist%next
               noi = sphere_cluster%sphere_order(i)
               npi1 = sphere_cluster%sphere_offset(i) + 1
               npi2 = sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i)
               noj = sphere_cluster%sphere_order(j)
               npj1 = sphere_cluster%sphere_offset(j) + 1
               npj2 = sphere_cluster%sphere_offset(j) + sphere_cluster%sphere_block(j)
               if (self%calculate_local_matrix) then
                  if (smopt .and. sphere_cluster%store_translation_matrix) then
      call self%stored_local_h_mat(idim)%configure(3, sphere_cluster%sphere_position(:, i) - sphere_cluster%sphere_position(:, j), &
                                                               rimedium, min(noi, noj) .ge. sphere_cluster%translation_switch_order)
                     loc_tranmat => self%stored_local_h_mat(idim)
                  else
                  call tranmat%configure(3, sphere_cluster%sphere_position(:, i) - sphere_cluster%sphere_position(:, j), rimedium, &
                                            min(noi, noj) .ge. sphere_cluster%translation_switch_order)
                     loc_tranmat => tranmat
                  end if
               else
                  loc_tranmat => self%stored_local_h_mat(idim)
               end if
               do rhs = 1, nrhs
                  if (rhslist(rhs)) then
                     call loc_tranmat%apply(noj, 2, noi, 2, ain(npj1:npj2, rhs), gout(npi1:npi2, rhs), &
                                            shift_op=.false., tran_op=contran(rhs))
                     call loc_tranmat%apply(noi, 2, noj, 2, ain(npi1:npi2, rhs), gout(npj1:npj2, rhs), &
                                            shift_op=.true., tran_op=contran(rhs))
                  end if
               end do
               if (self%calculate_local_matrix .and. (.not. (smopt .and. sphere_cluster%store_translation_matrix))) then
                  call tranmat%clear()
               end if
            end do
         end if
      end do

      if (numprocs .gt. 1 .and. mergeprocs) then
         nsend = sphere_cluster%number_eqns * nrhs
         call parallel_allreduce_sum(receive_buffer=gout, &
                                     mpi_number=nsend, mpi_comm=mpicomm)
      end if
   end subroutine local_sphere_to_sphere_expansion

   subroutine local_sphere_to_node_translation(self, nrhs, asphere, anode, &
                                               store_matrix_option, initial_run, rhs_list, &
                                               mpi_comm, con_tran, local_host, sphere_to_node, merge_procs)
      implicit none(type, external)
      class(fft_translation_plan_t), intent(inout), target :: self
      integer :: neqns, rank, numprocs, nsphere, nrhs, mpicomm, &
                 i, npi1, npi2, noi, task, proc, nsend, rhs, &
                 ndim, idim, localhost, nodei(3)
      logical :: smopt, rhslist(nrhs), contran(nrhs), spheretonode, mergeprocs
      logical, optional :: store_matrix_option, initial_run, merge_procs, &
                           rhs_list(nrhs), con_tran(nrhs), sphere_to_node
      integer, optional :: mpi_comm, local_host
      real(real64) :: rtran(3)
  complex(real64) :: asphere(sphere_cluster%number_eqns, nrhs), anode(self%cell_dim(1), self%cell_dim(2), self%cell_dim(3), self%node_order * (self%node_order + 2) * 2, nrhs), &
                         rimedium(2), anodet(self%node_order * (self%node_order + 2) * 2)
      type(translation_operator_state), pointer :: loc_tranmat
      type(translation_operator_state), target :: tranmat
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
         self%first_node_translation = initial_run
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
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      nsphere = sphere_cluster%number_spheres
      neqns = sphere_cluster%number_eqns
!
!  compute offsets for scattering coefficients
!
      if (self%first_node_translation) then
         if (smopt .and. sphere_cluster%store_translation_matrix) then
            task = 0
            ndim = 0
            do i = 1, nsphere
               if (sphere_cluster%host_sphere(i) .eq. localhost) then
                  task = task + 1
                  proc = mod(task, numprocs)
                  if (proc .eq. rank) ndim = ndim + 1
               end if
            end do
            if (allocated(self%stored_local_j_mat)) deallocate (self%stored_local_j_mat)
            allocate (self%stored_local_j_mat(ndim))
         end if
         self%calculate_node_matrix = .true.
         self%first_node_translation = .false.
      else
         self%calculate_node_matrix = (.not. (smopt .and. sphere_cluster%store_translation_matrix))
      end if

      if (spheretonode) then
         anode = 0.d0
      else
         asphere = 0.d0
      end if

      idim = 0
      task = 0
      do i = 1, nsphere
         if (sphere_cluster%host_sphere(i) .eq. localhost) then
            task = task + 1
            proc = mod(task, numprocs)
            if (proc .eq. rank) then
               idim = idim + 1
               call exterior_refractive_index(i, rimedium)
               nodei(:) = self%sphere_node(:, i)
               noi = sphere_cluster%sphere_order(i)
               npi1 = sphere_cluster%sphere_offset(i) + 1
               npi2 = sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i)
               rtran = self%d_cell * (dble(nodei) - .5d0) - sphere_cluster%sphere_position(:, i) + self%cell_boundary(:)
               if (self%calculate_node_matrix) then
                  if (smopt .and. sphere_cluster%store_translation_matrix) then
                     call self%stored_local_j_mat(idim)%configure(1, rtran, rimedium, &
                                                             min(noi, self%node_order) .ge. sphere_cluster%translation_switch_order)
                     loc_tranmat => self%stored_local_j_mat(idim)
                  else
                     call tranmat%configure(1, rtran, rimedium, &
                                            min(noi, self%node_order) .ge. sphere_cluster%translation_switch_order)
                     loc_tranmat => tranmat
                  end if
               else
                  loc_tranmat => self%stored_local_j_mat(idim)
               end if
               if (spheretonode) then
                  do rhs = 1, nrhs
                     if (rhslist(rhs)) then
                        anodet = 0.d0
                        call loc_tranmat%apply(noi, 2, self%node_order, 2, asphere(npi1:npi2, rhs), &
                                               anodet, shift_op=.false., tran_op=contran(rhs))
                        anode(nodei(1), nodei(2), nodei(3), :, rhs) &
                           = anode(nodei(1), nodei(2), nodei(3), :, rhs) + anodet(:)
                     end if
                  end do
               else
                  do rhs = 1, nrhs
                     if (rhslist(rhs)) then
                        anodet(:) = anode(nodei(1), nodei(2), nodei(3), :, rhs)
                        call loc_tranmat%apply(self%node_order, 2, noi, 2, anodet, &
                                               asphere(npi1:npi2, rhs), shift_op=.true., tran_op=contran(rhs))
                     end if
                  end do
               end if

               if (self%calculate_node_matrix .and. (.not. (smopt .and. sphere_cluster%store_translation_matrix))) then
                  call tranmat%clear()
               end if
            end if
         end if
      end do

      if (numprocs .gt. 1 .and. mergeprocs) then
         if (spheretonode) then
            nsend = self%cell_dim(1) * self%cell_dim(2) * self%cell_dim(3) * self%node_order * (self%node_order + 2) * 2 * nrhs
            call parallel_allreduce_sum(receive_buffer=anode, &
                                        mpi_number=nsend, mpi_comm=mpicomm)
         else
            nsend = sphere_cluster%number_eqns * nrhs
            call parallel_allreduce_sum(receive_buffer=asphere, &
                                        mpi_number=nsend, mpi_comm=mpicomm)
         end if
      end if
   end subroutine local_sphere_to_node_translation

   subroutine configure_fft_nodes(self, fva, target_min, target_max, d_specified, local_host, &
                                  requested_cell_size, requested_node_order, requested_neighbor_model)
      implicit none(type, external)
      class(fft_translation_plan_t), intent(inout) :: self
      logical :: dspec
      logical, optional :: d_specified
      integer :: nsphere, m, spherenode(3), n, i, node, ix, iy, iz, ir, j, icell, ncell, cell, lochost
      integer :: previous_cell_dimensions(3), previous_node_order, previous_local_host
      integer, optional :: local_host
      real(real64), intent(in) :: fva
      real(real64) :: r, fv, amean, svol, tvol, dd, targetmin(3), targetmax(3)
      real(real64), optional :: target_min(3), target_max(3)
      real(real64), optional, intent(in) :: requested_cell_size
      integer, optional, intent(in) :: requested_node_order
      integer, optional, intent(in) :: requested_neighbor_model
      complex(real64) :: previous_host_ref_index(2)
      type(linked_ilist), pointer :: ilist, ilist2

      previous_cell_dimensions = self%cell_dim
      previous_node_order = self%node_order
      previous_local_host = self%fft_local_host
      previous_host_ref_index = self%host_ref_index
      if (present(requested_neighbor_model)) self%neighbor_node_model = requested_neighbor_model

      if (present(target_min)) then
         targetmin = target_min
      else
         targetmin(:) = sphere_cluster%sphere_min_position
      end if
      if (present(target_max)) then
         targetmax = target_max
      else
         targetmax(:) = sphere_cluster%sphere_max_position
      end if
      if (present(d_specified)) then
         dspec = d_specified
      else
         dspec = .false.
      end if
      if (dspec .and. present(requested_cell_size)) self%d_cell = requested_cell_size
      if (dspec .and. self%d_cell <= 0.0_real64) dspec = .false.
      if (present(local_host)) then
         lochost = local_host
      else
         lochost = 0
      end if
      self%fft_local_host = lochost
      self%cell_boundary = targetmin
      if (self%fft_local_host .eq. 0) then
         self%host_ref_index = layer_ref_index(0)
      else
         self%host_ref_index = sphere_cluster%sphere_ref_index(:, self%fft_local_host)
      end if

      amean = 0.
      nsphere = 0
      svol = 0.
      do i = 1, sphere_cluster%number_spheres
         if (sphere_cluster%host_sphere(i) .eq. self%fft_local_host) then
            amean = amean + sphere_cluster%sphere_radius(i)
            svol = svol + sphere_cluster%sphere_radius(i)**3
            nsphere = nsphere + 1
         end if
      end do
      svol = svol**(1.d0 / 3.d0)
      amean = amean / dble(nsphere)
      self%fft_number_spheres = nsphere
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
         self%cell_volume_fraction = fv
         self%d_cell = (four_pi_over_three / fv / dble(nsphere))**(1.d0 / 3.d0) * svol
      else
         self%cell_volume_fraction = fva
      end if

!write(*,'('' host ri:'',2e13.5)') self%host_ref_index
!flush(6)

      self%cell_dim = ceiling((targetmax(:) - targetmin(:)) / self%d_cell)
      do m = 1, 3
         self%cell_dim(m) = next_supported_fft_size(self%cell_dim(m))
      end do
      self%d_cell = maxval((targetmax(:) - targetmin(:)) / dble(self%cell_dim(:)))
      if (present(requested_node_order)) then
         if (requested_node_order <= 0) then
            self%node_order = -requested_node_order + ceiling(self%d_cell)
         else
            self%node_order = requested_node_order
         end if
      end if
      if (allocated(self%cell_translation_matrix)) then
         if (any(previous_cell_dimensions /= self%cell_dim) .or. previous_node_order /= self%node_order .or. &
             previous_local_host /= self%fft_local_host .or. any(previous_host_ref_index /= self%host_ref_index)) then
            deallocate (self%cell_translation_matrix)
         end if
      end if

!write(*,'('' dcell:'',e13.5,3i5)') self%d_cell,self%cell_dim
!flush(6)
      call clear_fft_geometry(self)
      allocate (self%cell_list(self%cell_dim(1), self%cell_dim(2), self%cell_dim(3)), &
                self%sphere_node(3, sphere_cluster%number_spheres))
      self%cell_list(:, :, :)%number_elements = 0
      do iz = 1, self%cell_dim(3)
         do iy = 1, self%cell_dim(2)
            do ix = 1, self%cell_dim(1)
               allocate (self%cell_list(ix, iy, iz)%members)
            end do
         end do
      end do

!write(*,'('' d,n:'',e13.5,3i8)') self%d_cell,self%cell_dim
!flush(6)

      self%cell_origin(:) = self%d_cell * dble(self%cell_dim(:)) / 2.d0
      do i = 1, sphere_cluster%number_spheres
         if (sphere_cluster%host_sphere(i) .ne. self%fft_local_host) cycle
         do m = 1, 3
            r = sphere_cluster%sphere_position(m, i) - targetmin(m)
            node = floor(r / self%d_cell) + 1
            spherenode(m) = node
            spherenode(m) = min(spherenode(m), self%cell_dim(m))
            spherenode(m) = max(spherenode(m), 1)
         end do
         self%sphere_node(:, i) = spherenode
         n = self%cell_list(spherenode(1), spherenode(2), spherenode(3))%number_elements
!write(*,'('' i:'',5i5)') i,spherenode,n
!flush(6)
         ilist => self%cell_list(spherenode(1), spherenode(2), spherenode(3))%members
         do m = 1, n
            if (m .eq. n) allocate (ilist%next)
            ilist => ilist%next
         end do
         ilist%index = i
         ilist%next => null()
         self%cell_list(spherenode(1), spherenode(2), spherenode(3))%number_elements = n + 1
      end do

      n = -1
      do iz = -1, 1
         do iy = -1, 1
            do ix = -1, 1
               ir = ix * ix + iy * iy + iz * iz
               if (self%neighbor_node_model .eq. 0) then
                  if (ir .gt. 0) cycle
               elseif (self%neighbor_node_model .eq. 1) then
                  if (ir .gt. 1) cycle
               elseif (self%neighbor_node_model .eq. 2) then
                  if (ir .gt. 2) cycle
               end if
               n = n + 1
               self%neighbor_node(:, n) = (/ix, iy, iz/)
            end do
         end do
      end do
      self%number_neighbor_nodes = n

!write(*,'('' nn:'',i8)') self%number_neighbor_nodes
!flush(6)

      if (allocated(self%sphere_local_interaction_list)) deallocate (self%sphere_local_interaction_list)
      allocate (self%sphere_local_interaction_list(sphere_cluster%number_spheres))
      self%sphere_local_interaction_list(:)%number_elements = 0

      do i = 1, sphere_cluster%number_spheres
         if (sphere_cluster%host_sphere(i) .ne. self%fft_local_host) cycle
         allocate (self%sphere_local_interaction_list(i)%members)
         ilist2 => self%sphere_local_interaction_list(i)%members
         n = 0
         do cell = 0, self%number_neighbor_nodes
            spherenode = self%sphere_node(:, i) - self%neighbor_node(:, cell)
            if (any(spherenode .lt. (/1, 1, 1/)) .or. any(spherenode .gt. self%cell_dim)) cycle
            ncell = self%cell_list(spherenode(1), spherenode(2), spherenode(3))%number_elements
            ilist => self%cell_list(spherenode(1), spherenode(2), spherenode(3))%members
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
         self%sphere_local_interaction_list(i)%number_elements = n
      end do

   end subroutine configure_fft_nodes

   subroutine fft_node_to_node_translation(self, acoef, tranmat, gcoef, pmode, mpi_comm, tran_op)
      implicit none(type, external)
      class(fft_translation_plan_t), intent(inout) :: self
      logical :: tranop
      logical, optional :: tran_op
      integer :: nblk, celldim2(3), n, l, mpicomm, numprocs, rank, ncells, task, proc, nsend, pmode, ncells8
      integer, optional :: mpi_comm
    complex(real32) :: tranmat(8 * self%cell_dim(1) * self%cell_dim(2) * self%cell_dim(3), self%node_order * (self%node_order + 2), self%node_order * (self%node_order + 2))
     complex(real64) :: acoef(self%cell_dim(1) * self%cell_dim(2) * self%cell_dim(3), self%node_order * (self%node_order + 2), 2), &
                         aft(8 * self%cell_dim(1) * self%cell_dim(2) * self%cell_dim(3)), &
                        gcoef(self%cell_dim(1) * self%cell_dim(2) * self%cell_dim(3), self%node_order * (self%node_order + 2), 2), &
                         gft(8 * self%cell_dim(1) * self%cell_dim(2) * self%cell_dim(3), self%node_order * (self%node_order + 2))

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
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      nblk = self%node_order * (self%node_order + 2)
      celldim2 = 2 * self%cell_dim
      ncells = self%cell_dim(1) * self%cell_dim(2) * self%cell_dim(3)
      ncells8 = ncells * 8
      if (numprocs == 1) then
         call fft_node_to_node_translation_batched(self, acoef, tranmat, gcoef, pmode, tranop)
         return
      end if
      gft = 0.d0
      task = 0

      do n = 1, nblk
         task = task + 1
         proc = mod(task, numprocs)
         if (proc .eq. rank) then
            aft = 0.d0
            call transform_3d_fft(self, acoef(1:ncells, n, pmode), aft(1:ncells8), 1, self%cell_dim, celldim2, 1)
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
         call parallel_allreduce_sum( &
            receive_buffer=gft(1:ncells8, 1:nblk), &
            mpi_number=nsend, mpi_comm=mpicomm)
      end if

      task = 0
      do n = 1, nblk
         task = task + 1
         proc = mod(task, numprocs)
         if (proc .eq. rank) then
            gcoef(1:ncells, n, pmode) = 0.d0
            call transform_3d_fft(self, gcoef(1:ncells, n, pmode), gft(1:ncells8, n), 1, self%cell_dim, celldim2, -1)
         end if
      end do
      if (numprocs .gt. 1) then
         nsend = nblk * ncells
         call parallel_allreduce_sum( &
            receive_buffer=gcoef(1:ncells, 1:nblk, pmode), &
            mpi_number=nsend, mpi_comm=mpicomm)
      end if
      gcoef(1:ncells, 1:nblk, pmode) = gcoef(1:ncells, 1:nblk, pmode) / dble(ncells8)
   end subroutine fft_node_to_node_translation

   subroutine fft_node_to_node_translation_batched(self, acoef, tranmat, gcoef, pmode, tranop)
      class(fft_translation_plan_t), intent(inout) :: self
      logical, intent(in) :: tranop
      integer, intent(in) :: pmode
      integer :: cell_count, doubled_cell_count, doubled_dimensions(3), input_mode, output_mode, mode_count
      complex(real32), intent(in) :: &
         tranmat(8 * self%cell_dim(1) * self%cell_dim(2) * self%cell_dim(3), &
                 self%node_order * (self%node_order + 2), self%node_order * (self%node_order + 2))
      complex(real64), intent(in) :: &
         acoef(self%cell_dim(1) * self%cell_dim(2) * self%cell_dim(3), self%node_order * (self%node_order + 2), 2)
      complex(real64), intent(inout) :: &
         gcoef(self%cell_dim(1) * self%cell_dim(2) * self%cell_dim(3), self%node_order * (self%node_order + 2), 2)
      complex(real64), allocatable :: spatial_input(:, :), frequency_input(:, :), &
                                      frequency_output(:, :), spatial_output(:, :)

      mode_count = self%node_order * (self%node_order + 2)
      doubled_dimensions = 2 * self%cell_dim
      cell_count = product(self%cell_dim)
      doubled_cell_count = product(doubled_dimensions)
      allocate (spatial_input(mode_count, cell_count), frequency_input(mode_count, doubled_cell_count), &
                frequency_output(mode_count, doubled_cell_count), spatial_output(mode_count, cell_count))

      do concurrent(input_mode=1:mode_count)
         spatial_input(input_mode, :) = acoef(:, input_mode, pmode)
      end do
      frequency_input = 0.0_real64
      call transform_3d_fft(self, spatial_input, frequency_input, mode_count, self%cell_dim, doubled_dimensions, 1)

      frequency_output = 0.0_real64
      if (tranop) then
         do concurrent(output_mode=1:mode_count) local(input_mode)
            do input_mode = 1, mode_count
               frequency_output(output_mode, :) = frequency_output(output_mode, :) + &
                                                  tranmat(:, input_mode, output_mode) * frequency_input(input_mode, :)
            end do
         end do
      else
         do concurrent(output_mode=1:mode_count) local(input_mode)
            do input_mode = 1, mode_count
               frequency_output(output_mode, :) = frequency_output(output_mode, :) + &
                                                  tranmat(:, output_mode, input_mode) * frequency_input(input_mode, :)
            end do
         end do
      end if

      spatial_output = 0.0_real64
      call transform_3d_fft(self, spatial_output, frequency_output, mode_count, self%cell_dim, doubled_dimensions, -1)
      do concurrent(output_mode=1:mode_count)
         gcoef(:, output_mode, pmode) = spatial_output(output_mode, :) / real(doubled_cell_count, real64)
      end do
   end subroutine fft_node_to_node_translation_batched

   subroutine initialize_fft_translation_matrix(self, ri, nodr, p1, p2)
      implicit none(type, external)
      class(fft_translation_plan_t), intent(inout) :: self
      logical :: inhole
      integer :: nodr, nx, ny, nz, is, nblk, ncells2(3), isx, isy, isz, nx1, ny1, nz1, &
                 n, i, node(3), p1, p2, p, l
      real(real64) :: x, y, z, xp(3), r
      complex(real64) :: hij(nodr * (nodr + 2), nodr * (nodr + 2), 2), ri(2), &
                         htemp(1:2 * self%cell_dim(1), 1:2 * self%cell_dim(2), 1:2 * self%cell_dim(3))

      nblk = nodr * (nodr + 2)
      ncells2 = 2 * self%cell_dim

      if (light_up) then
         write (*, '('' fft ti '',i3,l)') mstm_global_rank, allocated(self%cell_translation_matrix)
         flush (6)
      end if
      if (allocated(self%cell_translation_matrix)) return
      if (allocated(self%cell_translation_matrix)) deallocate (self%cell_translation_matrix)
      allocate (self%cell_translation_matrix(0:2 * self%cell_dim(1) - 1, &
                                             0:2 * self%cell_dim(2) - 1, 0:2 * self%cell_dim(3) - 1, 1:nblk, 1:nblk, p1:p2))
      self%cell_translation_matrix = 0.d0

      do nz = 0, self%cell_dim(3)
         z = self%d_cell * dble(nz)
         do ny = 0, self%cell_dim(2)
            y = self%d_cell * dble(ny)
            do nx = 0, self%cell_dim(1)
               x = self%d_cell * dble(nx)
               inhole = .false.
               do i = 0, self%number_neighbor_nodes
                  node(:) = (/nx, ny, nz/) - self%neighbor_node(:, i)
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
                  nx1 = isx * nx - (isx - 1) * self%cell_dim(1)
                  ny1 = isy * ny - (isy - 1) * self%cell_dim(2)
                  nz1 = isz * nz - (isz - 1) * self%cell_dim(3)
                  if (nx1 .lt. ncells2(1) .and. ny1 .lt. ncells2(2) &
                      .and. nz1 .lt. ncells2(3)) then
                     xp(1) = isx * x
                     xp(2) = isy * y
                     xp(3) = isz * z
                     call generate_translation_matrix(nodr, nodr, translation_vector=xp, &
                                                      refractive_index=ri, ac_matrix=hij, vswf_type=3, &
                                                      mode_s=2, mode_t=2)
                     self%cell_translation_matrix(nx1, ny1, nz1, 1:nblk, 1:nblk, p1:p2) &
                        = self%cell_translation_matrix(nx1, ny1, nz1, 1:nblk, 1:nblk, p1:p2) &
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
               htemp(:, :, :) = self%cell_translation_matrix(:, :, :, l, n, p)
               call transform_3d_fft(self, htemp, htemp, 1, ncells2, ncells2, 1)
               self%cell_translation_matrix(:, :, :, l, n, p) = htemp(:, :, :)
            end do
         end do
      end do

   end subroutine initialize_fft_translation_matrix

   subroutine apply_fft_axis(ain, aout, nblk, ntot1, ntot2, ntot3in, ntot3out, ndimin, &
                             ndimout, isign, looporder, trig)
      implicit none
      integer :: nblk, ntot1, ntot2, ntot3in, ntot3out, l, m, n, isign, &
                 looporder(3), &
                 triplet(3), ndimin(3), ndimout(3), ntot3, i1, i2, i3
      integer, parameter :: mxtrig = 1000
      real(real64) :: trig(mxtrig), ar_temp(nblk, max(ntot3in, ntot3out)), &
                      ai_temp(nblk, max(ntot3in, ntot3out))
      complex(real64) :: aout(nblk, ndimout(1), ndimout(2), ndimout(3)), &
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
                  = cmplx(ar_temp(1:nblk, n), ai_temp(1:nblk, n), kind=real64)
            end do
         end do
      end do
   end subroutine apply_fft_axis

   subroutine transform_3d_fft(self, am, amf, nblk, ntot, ntot2, isign)
      implicit none(type, external)
      class(fft_translation_plan_t), intent(inout) :: self
      integer :: nblk, ntot(3), ntot2(3), isign
      real(real64) :: transform_start
      complex(real64) :: amf(nblk, ntot2(1), ntot2(2), ntot2(3)), &
                         am(nblk, ntot(1), ntot(2), ntot(3))
      self%metrics%transform_calls = self%metrics%transform_calls + 1
      transform_start = parallel_wall_time()
      if (ntot2(1) .ne. self%previous_transform_size(1) .or. ntot2(2) .ne. self%previous_transform_size(2) &
          .or. ntot2(3) .ne. self%previous_transform_size(3)) then
         self%previous_transform_size(1) = ntot2(1)
         self%previous_transform_size(2) = ntot2(2)
         self%previous_transform_size(3) = ntot2(3)
         call setgpfa(self%transform_trigonometry(:, 1), ntot2(1))
         call setgpfa(self%transform_trigonometry(:, 2), ntot2(2))
         call setgpfa(self%transform_trigonometry(:, 3), ntot2(3))
      end if
      if (isign .eq. 1) then
         call apply_fft_axis(am, amf, nblk, ntot(1), ntot(2), ntot(3), ntot2(3), ntot, &
                             ntot2, 1, (/1, 2, 3/), self%transform_trigonometry(:, 3))
         call apply_fft_axis(amf, amf, nblk, ntot2(3), ntot(1), ntot(2), ntot2(2), ntot2, &
                             ntot2, 1, (/3, 1, 2/), self%transform_trigonometry(:, 2))
         call apply_fft_axis(amf, amf, nblk, ntot2(3), ntot2(2), ntot(1), ntot2(1), ntot2, &
                             ntot2, 1, (/3, 2, 1/), self%transform_trigonometry(:, 1))
      else
         call apply_fft_axis(amf, amf, nblk, ntot2(1), ntot2(2), ntot2(3), ntot(3), ntot2, &
                             ntot2, -1, (/1, 2, 3/), self%transform_trigonometry(:, 3))
         call apply_fft_axis(amf, amf, nblk, ntot(3), ntot2(1), ntot2(2), ntot(2), ntot2, &
                             ntot2, -1, (/3, 1, 2/), self%transform_trigonometry(:, 2))
         call apply_fft_axis(amf, am, nblk, ntot(3), ntot(2), ntot2(1), ntot(1), ntot2, &
                             ntot, -1, (/3, 2, 1/), self%transform_trigonometry(:, 1))
      end if
      self%metrics%transform_time = self%metrics%transform_time + parallel_wall_time() - transform_start
   end subroutine transform_3d_fft

   pure logical function is_supported_fft_size(n)
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
      is_supported_fft_size = nn .eq. 1
   end function is_supported_fft_size

   pure integer function next_supported_fft_size(n)
      implicit none
      integer, intent(in) :: n
      integer :: n1
      n1 = n
      do while (.not. is_supported_fft_size(n1))
         n1 = n1 + 1
      end do
      next_supported_fft_size = n1
   end function next_supported_fft_size

end module fft_translation
