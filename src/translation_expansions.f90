module translation_expansions
   use angular_functions, only: generate_translation_matrix
   use iso_fortran_env, only: real64
   use mie, only: apply_mie_coefficients, exterior_refractive_index
   use mpidefs, only: mpi_comm_world, mstm_global_rank, mstm_mpi, mstm_mpi_sum
   use numerical_tables, only: light_up
   use periodic_lattice_subroutines, only: periodic_lattice, plane_boundary_lattice_interaction
   use spheredata, only: host_sphere, number_eqns, number_spheres, sphere_block, sphere_layer, sphere_offset, &
                         sphere_order, sphere_position, sphere_ref_index, store_translation_matrix, &
                         translation_switch_order
   use surface_subroutines, only: find_layer_index, plane_boundary_interaction, plane_surface_present
   use translation_operator, only: translation_operator_state

   implicit none(type, external)
   private

   public :: interaction_radius
   public :: clear_stored_translation_matrices
   public :: general_interaction_matrix
   public :: external_to_external_expansion
   public :: external_to_internal_expansion

   ! Non-owning geometry view. Associated arrays must outlive the view.
   type, public :: nested_sphere_geometry_view
      private
      integer :: sphere_count = 0
      integer, pointer :: hosts(:) => null(), blocks(:) => null(), offsets(:) => null(), orders(:) => null()
      real(real64), pointer :: positions(:, :) => null()
      complex(real64), pointer :: refractive_indices(:, :) => null()
   contains
      procedure, public :: configure => configure_nested_sphere_geometry_view
      procedure, public :: clear => clear_nested_sphere_geometry_view
      final :: finalize_nested_sphere_geometry_view
   end type nested_sphere_geometry_view

   real(real64), target :: interaction_radius = 1.0e10_real64
   type(translation_operator_state), target, allocatable :: stored_trans_mat(:)

contains

   subroutine configure_nested_sphere_geometry_view(self, hosts, blocks, offsets, orders, positions, refractive_indices)
      class(nested_sphere_geometry_view), intent(inout) :: self
      integer, target, intent(in) :: hosts(:)
      integer, target, intent(in) :: blocks(size(hosts)), offsets(size(hosts)), orders(size(hosts))
      real(real64), target, intent(in) :: positions(3, size(hosts))
      complex(real64), target, intent(in) :: refractive_indices(2, size(hosts))

      call self%clear()
      self%sphere_count = size(hosts)
      self%hosts => hosts
      self%blocks => blocks
      self%offsets => offsets
      self%orders => orders
      self%positions => positions
      self%refractive_indices => refractive_indices
   end subroutine configure_nested_sphere_geometry_view

   subroutine clear_nested_sphere_geometry_view(self)
      class(nested_sphere_geometry_view), intent(inout) :: self

      nullify (self%hosts, self%blocks, self%offsets, self%orders, self%positions, self%refractive_indices)
      self%sphere_count = 0
   end subroutine clear_nested_sphere_geometry_view

   subroutine finalize_nested_sphere_geometry_view(self)
      type(nested_sphere_geometry_view), intent(inout) :: self

      call self%clear()
   end subroutine finalize_nested_sphere_geometry_view

   subroutine clear_stored_translation_matrices(mat)
      implicit none(type, external)
      integer :: n
      type(translation_operator_state), allocatable, intent(inout) :: mat(:)
      if (.not. allocated(mat)) return
      n = size(mat)
      if (light_up) then
         write (*, '('' cstm 1 '',3i10)') mstm_global_rank, n
         flush (6)
      end if
      deallocate (mat)
   end subroutine clear_stored_translation_matrices

   subroutine general_interaction_matrix(matrix, mie_mult, mpi_comm)
      implicit none(type, external)
      logical :: miemult
      logical, optional, intent(in) :: mie_mult
      integer :: j, i, nbi, nbj, ilay, jlay, vtype, i1, i2, j1, j2, mpicomm, rank, numprocs, task, nsend
      integer, optional, intent(in) :: mpi_comm
      real(real64) :: rp(3), rdist
      complex(real64), intent(out) :: matrix(number_eqns, number_eqns)
      complex(real64) :: rimedium(2)
      complex(real64), allocatable :: fsmat(:, :, :), acmat(:, :), anp(:)
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      if (present(mie_mult)) then
         miemult = mie_mult
      else
         miemult = .true.
      end if
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)

      matrix = 0.0_real64
      task = 0
      do i = 1, number_spheres
         nbi = sphere_order(i) * (sphere_order(i) + 2)
         ilay = find_layer_index(sphere_position(3, i))
         do j = 1, number_spheres
            if (host_sphere(i) .ne. host_sphere(j) .and. host_sphere(j) .ne. i .and. host_sphere(i) .ne. j) cycle
            task = task + 1
            if (mod(task, numprocs) .ne. rank) cycle
            rp = sphere_position(:, i) - sphere_position(:, j)
            rdist = sqrt(sum(rp**2))
            nbj = sphere_order(j) * (sphere_order(j) + 2)
            jlay = find_layer_index(sphere_position(3, j))
            allocate (acmat(2 * nbi, 2 * nbj))
            acmat = 0.0_real64
            if (host_sphere(i) .ne. 0 .or. host_sphere(j) .ne. 0) then
               if (host_sphere(j) .eq. i) then
                  rimedium = sphere_ref_index(:, i)
                  vtype = 1
                  i1 = sphere_offset(i) + sphere_block(i) + 1
                  i2 = sphere_offset(i) + 2 * sphere_block(i)
                  j1 = sphere_offset(j) + 1
                  j2 = sphere_offset(j) + sphere_block(j)
               elseif (host_sphere(i) .eq. j) then
                  rimedium = sphere_ref_index(:, j)
                  vtype = 1
                  j1 = sphere_offset(j) + sphere_block(j) + 1
                  j2 = sphere_offset(j) + 2 * sphere_block(j)
                  i1 = sphere_offset(i) + 1
                  i2 = sphere_offset(i) + sphere_block(i)
               elseif (host_sphere(i) .eq. host_sphere(j)) then
                  if (rdist .gt. interaction_radius) cycle
                  rimedium = sphere_ref_index(:, host_sphere(i))
                  i1 = sphere_offset(i) + 1
                  i2 = sphere_offset(i) + sphere_block(i)
                  j1 = sphere_offset(j) + 1
                  j2 = sphere_offset(j) + sphere_block(j)
                  vtype = 3
               else
                  cycle
               end if
               allocate (fsmat(nbi, nbj, 2))
               call generate_translation_matrix(sphere_order(j), sphere_order(i), translation_vector=rp, &
                                                refractive_index=rimedium, ac_matrix=fsmat, vswf_type=vtype, &
                                                mode_s=2, mode_t=2, index_model=2)
               acmat(1:nbi, 1:nbj) = fsmat(1:nbi, 1:nbj, 1)
               acmat(nbi + 1:2 * nbi, nbj + 1:2 * nbj) = fsmat(1:nbi, 1:nbj, 2)
               deallocate (fsmat)
            else
               if (rdist .gt. interaction_radius) cycle
               i1 = sphere_offset(i) + 1
               i2 = sphere_offset(i) + sphere_block(i)
               j1 = sphere_offset(j) + 1
               j2 = sphere_offset(j) + sphere_block(j)
               if (periodic_lattice) then
                  if (plane_surface_present) then
                     call plane_boundary_lattice_interaction(sphere_order(i), sphere_order(j), &
                                                             rp(1), rp(2), sphere_position(3, i), sphere_position(3, j), &
                                                             acmat, include_source=.true., lr_transformation=.true., index_model=2)
                  else
                     allocate (fsmat(nbi, nbj, 2))
                     call plane_boundary_lattice_interaction(sphere_order(i), sphere_order(j), &
                                                             rp(1), rp(2), sphere_position(3, i), sphere_position(3, j), &
                                                             fsmat, include_source=.true., lr_transformation=.true., index_model=2)
                     acmat(1:nbi, 1:nbj) = acmat(1:nbi, 1:nbj) + fsmat(1:nbi, 1:nbj, 1)
                     acmat(nbi + 1:2 * nbi, nbj + 1:2 * nbj) = acmat(nbi + 1:2 * nbi, nbj + 1:2 * nbj) + fsmat(1:nbi, 1:nbj, 2)
                     deallocate (fsmat)
                  end if
               else
                  if (plane_surface_present) then
                     call plane_boundary_interaction(sphere_order(i), sphere_order(j), &
                                                     rp(1), rp(2), sphere_position(3, j), sphere_position(3, i), &
                                                     acmat, index_model=2, lr_transformation=.true., &
                                                     make_symmetric=.false.)
                  end if
                  if (ilay .eq. jlay) then
                     allocate (fsmat(nbi, nbj, 2))
                     call exterior_refractive_index(i, rimedium)
                     call generate_translation_matrix(sphere_order(j), sphere_order(i), translation_vector=rp, &
                                                      refractive_index=rimedium, ac_matrix=fsmat, vswf_type=3, &
                                                      mode_s=2, mode_t=2, index_model=2)
                     acmat(1:nbi, 1:nbj) = acmat(1:nbi, 1:nbj) + fsmat(1:nbi, 1:nbj, 1)
                     acmat(nbi + 1:2 * nbi, nbj + 1:2 * nbj) = acmat(nbi + 1:2 * nbi, nbj + 1:2 * nbj) + fsmat(1:nbi, 1:nbj, 2)
                     deallocate (fsmat)
                  end if
               end if
            end if
            matrix(i1:i2, j1:j2) = acmat(1:2 * nbi, 1:2 * nbj)
            deallocate (acmat)
         end do
      end do
      if (numprocs .gt. 1) then
         nsend = number_eqns * number_eqns
         call mstm_mpi(mpi_command='reduce', mpi_operation=mstm_mpi_sum, mpi_rank=0, &
                       mpi_recv_buf_dc=matrix, mpi_number=nsend, mpi_comm=mpicomm)
      end if
      if (rank .eq. 0) then
         if (miemult) then
            allocate (anp(number_eqns))
            do j = 1, number_eqns
               if (any(abs(matrix(1:number_eqns, j)) .ne. 0.0_real64)) then
                  call apply_mie_coefficients(number_eqns, 1, 1, matrix(1:number_eqns, j), anp)
                  matrix(1:number_eqns, j) = -anp(1:number_eqns)
               end if
               matrix(j, j) = matrix(j, j) + 1.0_real64
            end do
            deallocate (anp)
         end if
      end if
   end subroutine general_interaction_matrix
!
! outgoing translation operation:  a(i) = H(i-j) a(j).
! February 2013: number of rhs is a required argument.   mpi comm option added.
! this does not perform an allgather on output arrays.  that operation will be needed
! to use the results
!
   subroutine external_to_external_expansion(neqns, nrhs, ain, gout, &
                                             store_matrix_option, initial_run, rhs_list, &
                                             mpi_comm, con_tran)
      implicit none(type, external)
      integer, intent(in) :: neqns, nrhs
      integer :: rank, numprocs, nsphere, mpicomm, &
                 i, j, npi1, npi2, npj1, npj2, noj, noi, task, proc, &
                 nmin, ndim, idim, rhs
      logical :: smopt, rhslist(nrhs), contran(nrhs), rot, calcmat
      logical, save :: firstrun = .true.
      logical, optional, intent(in) :: store_matrix_option, initial_run, &
                                       rhs_list(nrhs), con_tran(nrhs)
      integer, optional, intent(in) :: mpi_comm
      real(real64) :: rdist
      complex(real64), intent(in) :: ain(neqns, nrhs)
      complex(real64), intent(out) :: gout(neqns, nrhs)
      complex(real64) :: rimedium(2)
      type(translation_operator_state), pointer :: loc_tranmat
      type(translation_operator_state), target :: tranmat
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      nsphere = number_spheres
      gout = 0.0_real64
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
!
!  compute offsets for scattering coefficients
!
      if (firstrun) then
         if (smopt .and. store_translation_matrix) then
            task = 0
            ndim = 0
            do i = 1, nsphere - 1
               do j = i + 1, nsphere
                  if (host_sphere(j) .eq. host_sphere(i) &
                      .and. sphere_layer(j) .eq. sphere_layer(i)) then
                     rdist = sqrt(sum((sphere_position(:, i) - sphere_position(:, j))**2))
                     if (rdist .gt. interaction_radius) cycle
                     task = task + 1
                     proc = mod(task, numprocs)
                     if (proc .eq. rank) ndim = ndim + 1
                  end if
               end do
            end do
            if (allocated(stored_trans_mat)) then
               call clear_stored_translation_matrices(stored_trans_mat)
            end if
            allocate (stored_trans_mat(ndim))
         end if
         calcmat = .true.
         firstrun = .false.
      else
         calcmat = (.not. (smopt .and. store_translation_matrix))
      end if

      idim = 0
      task = 0
      do i = 1, nsphere - 1
         do j = i + 1, nsphere
            if (host_sphere(j) .eq. host_sphere(i) &
                .and. sphere_layer(j) .eq. sphere_layer(i)) then
               rdist = sqrt(sum((sphere_position(:, i) - sphere_position(:, j))**2))
               if (rdist .gt. interaction_radius) cycle
               task = task + 1
               proc = mod(task, numprocs)
               if (proc .eq. rank) then
                  idim = idim + 1
                  call exterior_refractive_index(i, rimedium)
                  noi = sphere_order(i)
                  npi1 = sphere_offset(i) + 1
                  npi2 = sphere_offset(i) + sphere_block(i)
                  noj = sphere_order(j)
                  npj1 = sphere_offset(j) + 1
                  npj2 = sphere_offset(j) + sphere_block(j)
                  if (calcmat) then
                     nmin = min(noi, noj)
                     rot = (nmin .ge. translation_switch_order)
                     if (smopt .and. store_translation_matrix) then
                        call stored_trans_mat(idim)%configure(3, sphere_position(:, i) - sphere_position(:, j), &
                                                              rimedium, rot)
                        loc_tranmat => stored_trans_mat(idim)
                     else
                        call tranmat%configure(3, sphere_position(:, i) - sphere_position(:, j), rimedium, rot)
                        loc_tranmat => tranmat
                     end if
                  else
                     loc_tranmat => stored_trans_mat(idim)
                  end if
                  do rhs = 1, nrhs
                     if (.not. rhslist(rhs)) cycle
                     call loc_tranmat%apply(noj, 2, noi, 2, ain(npj1:npj2, rhs), gout(npi1:npi2, rhs), &
                                            shift_op=.false., tran_op=contran(rhs))
                     call loc_tranmat%apply(noi, 2, noj, 2, ain(npi1:npi2, rhs), gout(npj1:npj2, rhs), &
                                            shift_op=.true., tran_op=contran(rhs))
                  end do
                  if (calcmat .and. (.not. (smopt .and. store_translation_matrix))) then
                     call tranmat%clear()
                  end if
               end if
            end if
         end do
      end do

   end subroutine external_to_external_expansion
!
!  calculation of bmnp(i) = J(i-j) amnp(j) for i internal, host j=i, and
!  gmnp(i) = J(i-j) f(j), for host i = j.    This is the regular translation operation.   g and b
!  are returned ordered as (a,f) in the output array.
!  february 2013: number of rhs option
!  the routine does not perform an mpi reduce.
!
   subroutine external_to_internal_expansion(neqns, nrhs, ain, bout, &
                                             rhs_list, mpi_comm, con_tran, geometry)
      implicit none(type, external)
      integer, intent(in) :: neqns, nrhs
      integer :: rank, numprocs, mpicomm
      logical :: rhslist(nrhs), contran(nrhs)
      logical, optional, intent(in) :: rhs_list(nrhs), con_tran(nrhs)
      integer, optional, intent(in) :: mpi_comm
      complex(real64), intent(in) :: ain(neqns, nrhs)
      complex(real64), intent(out) :: bout(neqns, nrhs)
      type(nested_sphere_geometry_view), optional, intent(in) :: geometry
      type(nested_sphere_geometry_view) :: global_geometry

      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
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
      bout = 0.0_real64

      if (present(geometry)) then
         call apply_external_to_internal_expansion(geometry, neqns, nrhs, ain, bout, rhslist, contran, rank, numprocs)
      else
         call global_geometry%configure(host_sphere, sphere_block, sphere_offset, sphere_order, &
                                        sphere_position, sphere_ref_index)
         call apply_external_to_internal_expansion(global_geometry, neqns, nrhs, ain, bout, rhslist, contran, rank, numprocs)
      end if
   end subroutine external_to_internal_expansion

   subroutine apply_external_to_internal_expansion(geometry, neqns, nrhs, ain, bout, rhslist, contran, rank, numprocs)
      type(nested_sphere_geometry_view), intent(in) :: geometry
      integer, intent(in) :: neqns, nrhs, rank, numprocs
      logical, intent(in) :: rhslist(nrhs), contran(nrhs)
      complex(real64), intent(in) :: ain(neqns, nrhs)
      complex(real64), intent(inout) :: bout(neqns, nrhs)
      integer :: i, j, task, proc, extsurf, intsurf, ext1, ext2, int1, int2, noext, noint, rhs
      type(translation_operator_state) :: tranmat

      task = 0

      do i = 1, geometry%sphere_count - 1
         do j = i + 1, geometry%sphere_count
            if (geometry%hosts(j) .eq. i .or. geometry%hosts(i) .eq. j) then
               task = task + 1
               proc = mod(task, numprocs)
               if (proc .eq. rank) then
                  if (geometry%hosts(j) .eq. i) then
                     extsurf = j
                     intsurf = i
                  else
                     extsurf = i
                     intsurf = j
                  end if
                  noext = geometry%orders(extsurf)
                  noint = geometry%orders(intsurf)
                  call tranmat%configure(1, geometry%positions(:, intsurf) - geometry%positions(:, extsurf), &
                                         geometry%refractive_indices(:, intsurf), .true.)
                  ext1 = geometry%offsets(extsurf) + 1
                  ext2 = ext1 - 1 + geometry%blocks(extsurf)
                  int1 = geometry%offsets(intsurf) + 1 + geometry%blocks(intsurf)
                  int2 = int1 - 1 + geometry%blocks(intsurf)
                  do rhs = 1, nrhs
                     if (.not. rhslist(rhs)) cycle
                     call tranmat%apply(noext, 2, noint, 2, ain(ext1:ext2, rhs), &
                                        bout(int1:int2, rhs), shift_op=.false., tran_op=contran(rhs))
                     call tranmat%apply(noint, 2, noext, 2, ain(int1:int2, rhs), &
                                        bout(ext1:ext2, rhs), shift_op=.true., tran_op=contran(rhs))
                  end do
                  call tranmat%clear()
               end if
            end if
         end do
      end do
   end subroutine apply_external_to_internal_expansion

end module translation_expansions
