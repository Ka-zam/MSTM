module translation
   use angular_functions, only: atcdim, axialtrancoefrecurrence, cartosphere, ephicoef, gentranmatrix, moffset, rotcoef
   use coefficient_indexing, only: amnpaddress
   use iso_fortran_env, only: real64
   use mie, only: exteriorrefindex, multmiecoeffmult
   use mpidefs, only: mpi_comm_world, mstm_global_rank, mstm_mpi, mstm_mpi_sum, mstm_mpi_wtime
   use numerical_tables, only: light_up
   use periodic_lattice_subroutines, only: periodic_lattice, plane_boundary_lattice_interaction
   use spheredata, only: host_sphere, number_eqns, number_spheres, one_side_only, recalculate_surface_matrix, &
                         run_print_unit, sphere_block, sphere_layer, sphere_offset, sphere_order, sphere_position, &
                         sphere_ref_index, store_surface_matrix, store_translation_matrix, translation_switch_order
   use surface_subroutines, only: layer_id, layer_ref_index, plane_interaction, plane_surface_present
   use wave_functions, only: degree_transformation, mtransfer

   implicit none(type, external)
   private

   public :: translation_data
   public :: interaction_radius
   public :: clear_stored_trans_mat
   public :: general_interaction_matrix
   public :: periodic_lattice_sphere_interaction
   public :: spheresurfaceinteraction
   public :: external_to_external_expansion
   public :: external_to_internal_expansion
   public :: shiftcoefficient

   type translation_data
      private
      logical :: matrix_calculated = .false.
      logical :: rot_op = .false.
      logical :: zero_translation = .false.
      integer :: vswf_type = 0
      real(real64) :: translation_vector(3) = 0.0_real64
      real(real64), allocatable :: rot_mat(:, :)
      complex(real64) :: refractive_index(2) = (1.0_real64, 0.0_real64)
      complex(real64), allocatable :: phi_mat(:), z_mat(:), gen_mat(:, :, :)
   contains
      procedure, public :: configure => configure_translation
      procedure, public, pass(tranmat) :: apply => coefficient_translation
      procedure, public :: clear => clear_translation
      final :: finalize_translation
   end type translation_data
   type surface_ref_data
      logical :: symmetrical = .false.
      integer :: row_order = 0, col_order = 0
      complex(real64), allocatable :: matrix(:)
   end type surface_ref_data
   type pl_translation_data
      integer :: row_order = 0, col_order = 0
      complex(real64), allocatable :: matrix(:)
   end type pl_translation_data
   real(real64), target :: interaction_radius = 1.0e10_real64
   type(translation_data), target, allocatable :: stored_trans_mat(:)
   type(surface_ref_data), target, allocatable :: stored_ref(:)
   type(pl_translation_data), target, allocatable :: stored_plmat(:)

contains

   subroutine configure_translation(self, vswf_type, translation_vector, refractive_index, use_rotation)
      class(translation_data), intent(inout) :: self
      integer, intent(in) :: vswf_type
      real(real64), intent(in) :: translation_vector(3)
      complex(real64), intent(in) :: refractive_index(2)
      logical, intent(in) :: use_rotation

      call self%clear()
      self%vswf_type = vswf_type
      self%translation_vector = translation_vector
      self%refractive_index = refractive_index
      self%rot_op = use_rotation
   end subroutine configure_translation

   subroutine clear_translation(self)
      class(translation_data), intent(inout) :: self

      if (allocated(self%rot_mat)) deallocate (self%rot_mat)
      if (allocated(self%phi_mat)) deallocate (self%phi_mat)
      if (allocated(self%z_mat)) deallocate (self%z_mat)
      if (allocated(self%gen_mat)) deallocate (self%gen_mat)
      self%matrix_calculated = .false.
      self%zero_translation = .false.
      self%vswf_type = 0
      self%translation_vector = 0.0_real64
      self%refractive_index = (1.0_real64, 0.0_real64)
      self%rot_op = .false.
   end subroutine clear_translation

   subroutine finalize_translation(self)
      type(translation_data), intent(inout) :: self

      call self%clear()
   end subroutine finalize_translation

   subroutine clear_stored_trans_mat(mat)
      implicit none(type, external)
      integer :: n
      type(translation_data), target, allocatable, intent(inout) :: mat(:)
      if (.not. allocated(mat)) return
      n = size(mat)
      if (light_up) then
         write (*, '('' cstm 1 '',3i10)') mstm_global_rank, n
         flush (6)
      end if
      if (light_up) then
         write (*, '('' cstm 2 '',3i10)') mstm_global_rank, n
         flush (6)
      end if
      deallocate (mat)
   end subroutine clear_stored_trans_mat

   subroutine clear_stored_ref_mat(mat)
      implicit none(type, external)
      type(surface_ref_data), allocatable, intent(inout) :: mat(:)
      if (.not. allocated(mat)) return
      deallocate (mat)
   end subroutine clear_stored_ref_mat

   subroutine clear_stored_pl_mat(mat)
      implicit none(type, external)
      type(pl_translation_data), allocatable, intent(inout) :: mat(:)
      if (.not. allocated(mat)) return
      deallocate (mat)
   end subroutine clear_stored_pl_mat

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

      matrix = 0.d0
      task = 0
      do i = 1, number_spheres
         nbi = sphere_order(i) * (sphere_order(i) + 2)
         ilay = layer_id(sphere_position(3, i))
         do j = 1, number_spheres
            if (host_sphere(i) .ne. host_sphere(j) .and. host_sphere(j) .ne. i .and. host_sphere(i) .ne. j) cycle
            task = task + 1
            if (mod(task, numprocs) .ne. rank) cycle
            rp = sphere_position(:, i) - sphere_position(:, j)
            rdist = sqrt(sum(rp**2))
            nbj = sphere_order(j) * (sphere_order(j) + 2)
            jlay = layer_id(sphere_position(3, j))
            allocate (acmat(2 * nbi, 2 * nbj))
            acmat = 0.d0
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
               call gentranmatrix(sphere_order(j), sphere_order(i), translation_vector=rp, &
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
                     call plane_interaction(sphere_order(i), sphere_order(j), &
                                            rp(1), rp(2), sphere_position(3, j), sphere_position(3, i), &
                                            acmat, index_model=2, lr_transformation=.true., &
                                            make_symmetric=.false.)
                  end if
                  if (ilay .eq. jlay) then
                     allocate (fsmat(nbi, nbj, 2))
                     call exteriorrefindex(i, rimedium)
                     call gentranmatrix(sphere_order(j), sphere_order(i), translation_vector=rp, &
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
               if (any(abs(matrix(1:number_eqns, j)) .ne. 0.d0)) then
                  call multmiecoeffmult(number_eqns, 1, 1, matrix(1:number_eqns, j), anp)
                  matrix(1:number_eqns, j) = -anp(1:number_eqns)
               end if
               matrix(j, j) = matrix(j, j) + 1.d0
            end do
            deallocate (anp)
         end if
      end if
   end subroutine general_interaction_matrix

   subroutine periodic_lattice_sphere_interaction(neqns, nrhs, ain, aout, &
                                                  initial_run, rhs_list, mpi_comm, con_tran, store_matrix_option)
      implicit none(type, external)
      integer, intent(in) :: neqns, nrhs
      integer :: rank, numprocs, nmat, mpicomm, proc, i, j, rhs, &
                 i1, i2, j1, j2, task, rank0, nbi, nbj, rmatdim
      logical :: initrun, rhslist(nrhs), calcmat, contran(nrhs), smopt
      logical, optional, intent(in) :: initial_run, rhs_list(nrhs), con_tran(nrhs), store_matrix_option
      integer, save :: nmat_tot
      integer, optional, intent(in) :: mpi_comm
      real(real64) :: rp(3), time1, time2
      complex(real64), intent(inout) :: ain(neqns, nrhs)
      complex(real64), intent(out) :: aout(neqns, nrhs)
      complex(real64) :: ri
      type(pl_translation_data), target :: rmat
      type(pl_translation_data), pointer :: loc_rmat

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
      if (present(store_matrix_option)) then
         smopt = store_matrix_option
      else
         smopt = .true.
      end if
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank0)
      calcmat = initrun
      if (.not. (smopt .and. store_translation_matrix)) calcmat = .true.

      if (calcmat) then
         task = 0
         nmat = 0
         do i = 1, number_spheres
            do j = 1, number_spheres
               if ((host_sphere(i) .eq. 0) .and. (host_sphere(j) .eq. 0)) then
                  task = task + 1
                  proc = mod(task, numprocs)
                  if (proc .eq. rank) then
                     nmat = nmat + 1
                  end if
               end if
            end do
         end do
         nmat_tot = nmat
         if (smopt .and. store_translation_matrix) then
            call clear_stored_pl_mat(stored_plmat)
            allocate (stored_plmat(nmat))
         end if
      end if

      task = 0
      nmat = 0
      aout = 0.d0
      time1 = mstm_mpi_wtime()
      do i = 1, number_spheres
         i1 = sphere_offset(i) + 1
         i2 = sphere_offset(i) + sphere_block(i)
         nbi = sphere_order(i) * (sphere_order(i) + 2)
         do j = 1, number_spheres
            if ((host_sphere(i) .eq. 0) .and. (host_sphere(j) .eq. 0)) then
               task = task + 1
               proc = mod(task, numprocs)
               if (proc .eq. rank) then
                  j1 = sphere_offset(j) + 1
                  j2 = sphere_offset(j) + sphere_block(j)
                  nbj = sphere_order(j) * (sphere_order(j) + 2)
                  nmat = nmat + 1
                  if (calcmat) then
                     rp = sphere_position(:, i) - sphere_position(:, j)
                     ri = layer_ref_index(sphere_layer(i))
                     if (plane_surface_present) then
                        rmatdim = 4 * nbi * nbj
                     else
                        rmatdim = 2 * nbi * nbj
                     end if
                     if (smopt .and. store_translation_matrix) then
                        stored_plmat(nmat)%row_order = sphere_order(i)
                        stored_plmat(nmat)%col_order = sphere_order(j)
                        allocate (stored_plmat(nmat)%matrix(rmatdim))
                        loc_rmat => stored_plmat(nmat)
                     else
                        rmat%row_order = sphere_order(i)
                        rmat%col_order = sphere_order(j)
                        allocate (rmat%matrix(rmatdim))
                        loc_rmat => rmat
                     end if
                     call plane_boundary_lattice_interaction(sphere_order(i), sphere_order(j), &
                                                             rp(1), rp(2), sphere_position(3, i), sphere_position(3, j), &
                                                    loc_rmat%matrix, include_source=.true., lr_transformation=.true., index_model=2)
                     if (smopt .and. store_translation_matrix .and. rank0 .eq. 0) then
                        time2 = mstm_mpi_wtime()
                        if (time2 - time1 .ge. 15.d0) then
                           write (run_print_unit, '('' assembling pl matrix '',i5,''/'',i5)') &
                              nmat, nmat_tot
                           flush (run_print_unit)
                           time1 = time2
                        end if
                     end if
!if(mstm_global_rank.le.3) then
! write(*,'('' rank,comm:'',4i12,l2)') rank,mstm_global_rank,i,j,contran(1)
!flush(6)
!endif
                  else
                     loc_rmat => stored_plmat(nmat)
                  end if
                  do rhs = 1, nrhs
                     if (.not. rhslist(rhs)) cycle
                     if (plane_surface_present) then
                        if (.not. contran(rhs)) then
                           call pl_matrix_mult(sphere_order(i), sphere_order(j), ain(j1:j2, rhs), aout(i1:i2, rhs), &
                                               .false., pb_mat=loc_rmat%matrix)
                        else
                           call pl_matrix_mult(sphere_order(i), sphere_order(j), aout(j1:j2, rhs), ain(i1:i2, rhs), &
                                               .true., pb_mat=loc_rmat%matrix)
                        end if
                     else
                        if (.not. contran(rhs)) then
                           call pl_matrix_mult(sphere_order(i), sphere_order(j), ain(j1:j2, rhs), aout(i1:i2, rhs), &
                                               .false., fs_mat=loc_rmat%matrix)
                        else
                           call pl_matrix_mult(sphere_order(i), sphere_order(j), aout(j1:j2, rhs), ain(i1:i2, rhs), &
                                               .true., fs_mat=loc_rmat%matrix)
                        end if
                     end if
                  end do
                  if (.not. (smopt .and. store_translation_matrix)) deallocate (rmat%matrix)
               end if
            end if
         end do
      end do
      if (store_surface_matrix) recalculate_surface_matrix = .false.
   end subroutine periodic_lattice_sphere_interaction

   subroutine pl_matrix_mult(nodrt, nodrs, as, at, tran, fs_mat, pb_mat)
      implicit none(type, external)
      logical, intent(in) :: tran
      integer, intent(in) :: nodrt, nodrs
      integer :: p
      complex(real64), intent(inout) :: as(nodrs * (nodrs + 2), 2), at(nodrt * (nodrt + 2), 2)
      complex(real64), optional, intent(in) :: fs_mat(nodrt * (nodrt + 2), nodrs * (nodrs + 2), 2), &
                                               pb_mat(nodrt * (nodrt + 2), 2, nodrs * (nodrs + 2), 2)

      if (.not. tran) then
         if (present(fs_mat)) then
            do p = 1, 2
               at(:, p) = at(:, p) + matmul(fs_mat(:, :, p), as(:, p))
            end do
         else
            do p = 1, 2
               at(:, p) = at(:, p) + matmul(pb_mat(:, p, :, 1), as(:, 1)) + matmul(pb_mat(:, p, :, 2), as(:, 2))
            end do
         end if
      else
         if (present(fs_mat)) then
            do p = 1, 2
               as(:, p) = as(:, p) + matmul(at(:, p), fs_mat(:, :, p))
            end do
         else
            do p = 1, 2
               as(:, p) = as(:, p) + matmul(at(:, 1), pb_mat(:, 1, :, p)) + matmul(at(:, 2), pb_mat(:, 2, :, p))
            end do
         end if
      end if
   end subroutine pl_matrix_mult

   subroutine spheresurfaceinteraction(neqns, nrhs, ain, aout, &
                                       initial_run, rhs_list, mpi_comm, con_tran)
      implicit none(type, external)
      integer, intent(in) :: neqns, nrhs
      integer :: rank, numprocs, nmat, mpicomm, proc, i, j, rhs, &
                 i1, i2, j1, j2, task, jstart, rmatdim, rank0
      logical :: initrun, rhslist(nrhs), calcmat, contran(nrhs), rmatsymm
      logical, optional, intent(in) :: initial_run, rhs_list(nrhs), con_tran(nrhs)
      integer, save :: nmat_tot
      integer, optional, intent(in) :: mpi_comm
      real(real64) :: rp(3), time1, time2
      complex(real64), intent(inout) :: ain(neqns, nrhs)
      complex(real64), intent(out) :: aout(neqns, nrhs)
      complex(real64), allocatable :: atempi(:), atempj(:), &
                                      atempi2(:), atempj2(:)
      type(surface_ref_data), target :: rmat
      type(surface_ref_data), pointer :: loc_rmat

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
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank0)
      calcmat = initrun .and. recalculate_surface_matrix
      if (.not. store_surface_matrix) calcmat = .true.

      if (calcmat) then
         task = 0
         nmat = 0
         do i = 1, number_spheres
            if (one_side_only) then
               jstart = i
            else
               jstart = 1
            end if
            do j = jstart, number_spheres
               if ((host_sphere(i) .eq. 0) .and. (host_sphere(j) .eq. 0)) then
                  task = task + 1
                  proc = mod(task, numprocs)
                  if (proc .eq. rank) then
                     nmat = nmat + 1
                  end if
               end if
            end do
         end do
         nmat_tot = nmat
         if (store_surface_matrix) then
            call clear_stored_ref_mat(stored_ref)
            allocate (stored_ref(nmat))
         end if
      end if

      task = 0
      nmat = 0
      aout = 0.d0
      time1 = mstm_mpi_wtime()
      do i = 1, number_spheres
         i1 = sphere_offset(i) + 1
         i2 = sphere_offset(i) + sphere_block(i)
         if (one_side_only) then
            jstart = i
         else
            jstart = 1
         end if
         do j = jstart, number_spheres
            if ((host_sphere(i) .eq. 0) .and. (host_sphere(j) .eq. 0)) then
               task = task + 1
               proc = mod(task, numprocs)
               if (proc .eq. rank) then
                  j1 = sphere_offset(j) + 1
                  j2 = sphere_offset(j) + sphere_block(j)
                  nmat = nmat + 1
                  if (calcmat) then
!if(mstm_global_rank.eq.0) then
!write(*,'(2i6,l2)') i,j,calcmat
!flush(6)
!endif
                     rp = sphere_position(:, i) - sphere_position(:, j)
                     rmatsymm = (abs(rp(1)) .lt. 1.d-4 .and. abs(rp(2)) .lt. 1.d-4 &
                                 .and. (.not. periodic_lattice))
                     if (rmatsymm) then
                        rmatdim = 2 * atcdim(sphere_order(i), sphere_order(j))
                     else
                        rmatdim = sphere_block(i) * sphere_block(j)
                     end if
                     if (store_surface_matrix) then
                        stored_ref(nmat)%row_order = sphere_order(i)
                        stored_ref(nmat)%col_order = sphere_order(j)
                        stored_ref(nmat)%symmetrical = rmatsymm
                        allocate (stored_ref(nmat)%matrix(rmatdim))
                        loc_rmat => stored_ref(nmat)
                     else
                        rmat%row_order = sphere_order(i)
                        rmat%col_order = sphere_order(j)
                        rmat%symmetrical = rmatsymm
                        allocate (rmat%matrix(rmatdim))
                        loc_rmat => rmat
                     end if
                     call plane_interaction(sphere_order(i), sphere_order(j), &
                                            rp(1), rp(2), sphere_position(3, j), sphere_position(3, i), &
                                            loc_rmat%matrix, &
                                            index_model=2, lr_transformation=.true., &
                                            make_symmetric=rmatsymm)
                     if (store_surface_matrix .and. rank0 .eq. 0) then
                        time2 = mstm_mpi_wtime()
                        if (time2 - time1 .ge. 15.d0) then
                           write (run_print_unit, '('' assembling surf matrix '',i5,''/'',i5)') &
                              nmat, nmat_tot
                           flush (run_print_unit)
                           time1 = time2
                        end if
                     end if
                  else
                     loc_rmat => stored_ref(nmat)
                  end if
                  allocate (atempj(sphere_block(j)), atempi(sphere_block(i)))
                  if ((j .ne. i) .and. one_side_only) allocate (atempj2(sphere_block(j)), atempi2(sphere_block(i)))
                  do rhs = 1, nrhs
                     if (.not. rhslist(rhs)) cycle
                     if (.not. contran(rhs)) then
                        call surface_interaction_matrix_mult(sphere_order(j), sphere_order(i), ain(j1:j2, rhs), atempi, loc_rmat, 1)
                        aout(i1:i2, rhs) = aout(i1:i2, rhs) + atempi(1:sphere_block(i))
                     else
                        call surface_interaction_matrix_mult(sphere_order(i), sphere_order(j), ain(i1:i2, rhs), atempj, loc_rmat, 2)
                        aout(j1:j2, rhs) = aout(j1:j2, rhs) + atempj(1:sphere_block(j))
                     end if
                     if ((j .ne. i) .and. one_side_only) then
                        if (.not. contran(rhs)) then
                           call degree_transformation(sphere_order(i), &
                                                      ain(i1:i2, rhs), atempi)
                           call surface_interaction_matrix_mult(sphere_order(i), sphere_order(j), atempi, atempj, loc_rmat, 2)
                           call degree_transformation(sphere_order(j), &
                                                      atempj, atempj2)
                           aout(j1:j2, rhs) = aout(j1:j2, rhs) &
                                              + atempj2(1:sphere_block(j))
                        else
                           call degree_transformation(sphere_order(j), &
                                                      ain(j1:j2, rhs), atempj)
                           call surface_interaction_matrix_mult(sphere_order(j), sphere_order(i), atempj, atempi, loc_rmat, 1)
                           call degree_transformation(sphere_order(i), &
                                                      atempi, atempi2)
                           aout(i1:i2, rhs) = aout(i1:i2, rhs) + atempi2(1:sphere_block(i))
                        end if
                     end if
                  end do
                  deallocate (atempi, atempj)
                  if ((j .ne. i) .and. one_side_only) deallocate (atempi2, atempj2)
                  if (.not. store_surface_matrix) deallocate (rmat%matrix)
               end if
            end if
         end do
      end do
      if (store_surface_matrix) recalculate_surface_matrix = .false.
   end subroutine spheresurfaceinteraction

   subroutine surface_interaction_matrix_mult(nin, nout, ain, aout, rmat, dir)
      implicit none(type, external)
      integer, intent(in) :: nin, nout, dir
      integer :: bin, bout, m, m1, n, l, p, q, mnp, klq, i, moff
      complex(real64), intent(in) :: ain(2 * nin * (nin + 2))
      complex(real64), intent(out) :: aout(2 * nout * (nout + 2))
      type(surface_ref_data), intent(in) :: rmat
      bin = 2 * nin * (nin + 2)
      bout = 2 * nout * (nout + 2)
      aout = 0.d0
      if (rmat%symmetrical) then
         do m = -min(nin, nout), min(nin, nout)
            m1 = max(abs(m), 1)
            moff = 2 * moffset(m, nin, nout)
            do n = m1, nout
               do p = 1, 2
                  mnp = amnpaddress(m, n, p, nout, 2)
                  do l = m1, nin
                     do q = 1, 2
                        klq = amnpaddress(m, l, q, nin, 2)
                        if (dir .eq. 1) then
                           i = moff + p + 2 * (n - m1) + 2 * (nout - m1 + 1) * (q - 1 + 2 * (l - m1))
                        else
                           i = moff + q + 2 * (l - m1) + 2 * (nin - m1 + 1) * (p - 1 + 2 * (n - m1))
                        end if
                        aout(mnp) = aout(mnp) + rmat%matrix(i) * ain(klq)
                     end do
                  end do
               end do
            end do
         end do
      else
         if (dir .eq. 1) then
            do n = 1, bout
               do l = 1, bin
                  i = n + bout * (l - 1)
                  aout(n) = aout(n) + rmat%matrix(i) * ain(l)
               end do
            end do
         else
            do n = 1, bout
               do l = 1, bin
                  i = l + bin * (n - 1)
                  aout(n) = aout(n) + rmat%matrix(i) * ain(l)
               end do
            end do
         end if
      end if
   end subroutine surface_interaction_matrix_mult
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
      logical :: smopt, rhslist(nrhs), contran(nrhs), rot
      logical, save :: calcmat, firstrun = .true.
      logical, optional, intent(in) :: store_matrix_option, initial_run, &
                                       rhs_list(nrhs), con_tran(nrhs)
      integer, optional, intent(in) :: mpi_comm
      real(real64) :: rdist
      complex(real64), intent(in) :: ain(neqns, nrhs)
      complex(real64), intent(out) :: gout(neqns, nrhs)
      complex(real64) :: rimedium(2)
      type(translation_data), pointer :: loc_tranmat
      type(translation_data), target :: tranmat
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call mstm_mpi(mpi_command='size', mpi_size=numprocs, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      nsphere = number_spheres
      gout = 0.
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
               call clear_stored_trans_mat(stored_trans_mat)
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
                  call exteriorrefindex(i, rimedium)
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
                                             rhs_list, mpi_comm, con_tran)
      implicit none(type, external)
      integer, intent(in) :: neqns, nrhs
      integer :: rank, numprocs, nsphere, mpicomm, &
                 i, j, task, proc, extsurf, intsurf, ext1, ext2, &
                 int1, int2, noext, noint, rhs
      logical :: rhslist(nrhs), contran(nrhs)
      logical, optional, intent(in) :: rhs_list(nrhs), con_tran(nrhs)
      integer, optional, intent(in) :: mpi_comm
      complex(real64), intent(in) :: ain(neqns, nrhs)
      complex(real64), intent(out) :: bout(neqns, nrhs)
      type(translation_data) :: tranmat
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
      bout = 0.d0
      nsphere = number_spheres
      task = 0

      do i = 1, nsphere - 1
         do j = i + 1, nsphere
            if (host_sphere(j) .eq. i .or. host_sphere(i) .eq. j) then
               task = task + 1
               proc = mod(task, numprocs)
               if (proc .eq. rank) then
                  if (host_sphere(j) .eq. i) then
                     extsurf = j
                     intsurf = i
                  else
                     extsurf = i
                     intsurf = j
                  end if
                  noext = sphere_order(extsurf)
                  noint = sphere_order(intsurf)
                  call tranmat%configure(1, sphere_position(:, intsurf) - sphere_position(:, extsurf), &
                                         sphere_ref_index(:, intsurf), .true.)
                  ext1 = sphere_offset(extsurf) + 1
                  ext2 = ext1 - 1 + sphere_block(extsurf)
                  int1 = sphere_offset(intsurf) + 1 + sphere_block(intsurf)
                  int2 = int1 - 1 + sphere_block(intsurf)
                  do rhs = 1, nrhs
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
   end subroutine external_to_internal_expansion

   subroutine coefficient_translation(nodra, nmodea, nodrg, nmodeg, &
                                      acoef, gcoef, tranmat, shift_op, tran_op)
      implicit none(type, external)
      logical :: sop, top, rot
      logical, optional, intent(in) :: shift_op, tran_op
      integer, intent(in) :: nodra, nodrg, nmodea, nmodeg
      integer :: nblka, nblkg, lengtha, lengthg, nmode, shiftvec(2), n, m, im, nn1, nn2, n1, &
                 p, nmin, m1, offset, blocksize, nmax, tdim, vtype, nodrs, nodrt
      real(real64) :: r, ct, rtran(3)
      complex(real64), intent(in) :: acoef(*)
      complex(real64), intent(inout) :: gcoef(*)
      complex(real64) :: &
         a_t(0:nodra + 1, nodra, nmodea), g_t(0:nodrg + 1, nodrg, nmodeg), &
         g_shift(0:nodrg + 1, nodrg, nmodeg), &
         a_tt(-nodra:nodra, nodra, 2), g_tt(-nodrg:nodrg, nodrg, 2), &
         atc(max(nodra, nodrg), max(nodra, nodrg), 2), &
         a_t2(nodra * (nodra + 2), 2), g_t2(nodrg * (nodrg + 2), 2), rimed(2), ephi
      class(translation_data), intent(inout) :: tranmat
      nblka = nodra * (nodra + 2)
      nblkg = nodrg * (nodrg + 2)
      nmin = min(nodra, nodrg)
      nmax = max(nodra, nodrg)
      nmode = max(nmodea, nmodeg)
      if (present(shift_op)) then
         sop = shift_op
      else
         sop = .false.
      end if
      if (present(tran_op)) then
         top = tran_op
      else
         top = .false.
      end if
      rot = tranmat%rot_op

      if (.not. tranmat%matrix_calculated) then
         vtype = tranmat%vswf_type
         rimed = tranmat%refractive_index
         rtran = tranmat%translation_vector
         r = dot_product(rtran, rtran)
         if (r .lt. 1.d-12) then
            tranmat%zero_translation = .true.
         else
            tranmat%zero_translation = .false.
         end if
         if (sop) then
            nodrs = nodrg
            nodrt = nodra
         else
            nodrs = nodra
            nodrt = nodrg
         end if
         if (.not. tranmat%zero_translation) then
            if (rot) then
               tdim = atcdim(nodrt, nodrs)
               allocate (tranmat%rot_mat(-nmin:nmin, 0:nmax * (nmax + 2)))
               allocate (tranmat%phi_mat(-nmax:nmax))
               allocate (tranmat%z_mat(1:tdim))
               call cartosphere(rtran, r, ct, ephi)
               call rotcoef(ct, nmin, nmax, tranmat%rot_mat)
               call axialtrancoefrecurrence(vtype, r, rimed, nodrt, nodrs, &
                                            tdim, tranmat%z_mat)
               call ephicoef(ephi, nmax, tranmat%phi_mat)
            else
               allocate (tranmat%gen_mat(nodrt * (nodrt + 2), nodrs * (nodrs + 2), 2))
               call gentranmatrix(nodrs, nodrt, translation_vector=rtran, &
                                  refractive_index=rimed, ac_matrix=tranmat%gen_mat, vswf_type=vtype, &
                                  mode_s=2, mode_t=2)
            end if
         end if
         tranmat%matrix_calculated = .true.
      end if
      lengtha = nblka * nmodea
      lengthg = nblkg * nmodeg
      shiftvec = (/1, 1/)
      if (sop .neqv. top) then
         shiftvec = -shiftvec
      end if
      if (sop) then
         im = -1
      else
         im = 1
      end if
      if (tranmat%zero_translation) then
         if (tranmat%vswf_type .eq. 1) then
            call mtransfer(nodra, nodrg, acoef(1:lengtha), g_t)
            gcoef(1:lengthg) = gcoef(1:lengthg) &
                               + reshape(g_t(0:nodrg + 1, 1:nodrg, 1:nmodeg), (/lengthg/))
         end if
      else
         if (rot) then
            call shiftcoefficient(nodra, nmodea, shiftvec(1), shiftvec(2), &
                                  acoef(1:lengtha), a_t(0:nodra + 1, 1:nodra, 1:2))
            a_tt(0, 1:nodra, 1:2) = a_t(0, 1:nodra, 1:2)
            do m = 1, nodra
               a_tt(m, m:nodra, 1:2) = a_t(m, m:nodra, 1:2) * tranmat%phi_mat(im * m)
               a_tt(-m, m:nodra, 1:2) = a_t(m + 1:nodra + 1, m, 1:2) * tranmat%phi_mat(-im * m)
            end do
            do n = 1, nodra
               nn1 = n * (n + 1) - n
               nn2 = nn1 + (2 * n + 1) - 1
               n1 = min(n, nodrg)
               a_tt(-n1:n1, n, 1:2) &
                  = matmul(tranmat%rot_mat(-n1:n1, nn1:nn2), a_tt(-n:n, n, 1:2))
            end do
            do m = -nmin, nmin
               m1 = max(1, abs(m))
               if (sop) then
                  offset = moffset(m, nodra, nodrg)
                  blocksize = (nodrg - m1 + 1) * (nodra - m1 + 1) * 2
                  atc(m1:nodra, m1:nodrg, 1:2) = &
                     reshape(tranmat%z_mat(offset + 1:offset + blocksize), &
                             (/nodra - m1 + 1, nodrg - m1 + 1, 2/))
                  do p = 1, 2
                     g_tt(m, m1:nodrg, p) &
                        = matmul(a_tt(m, m1:nodra, p), atc(m1:nodra, m1:nodrg, p))
                  end do
               else
                  offset = moffset(m, nodrg, nodra)
                  blocksize = (nodrg - m1 + 1) * (nodra - m1 + 1) * 2
                  atc(m1:nodrg, m1:nodra, 1:2) = &
                     reshape(tranmat%z_mat(offset + 1:offset + blocksize), &
                             (/nodrg - m1 + 1, nodra - m1 + 1, 2/))
                  do p = 1, 2
                     g_tt(m, m1:nodrg, p) &
                        = matmul(atc(m1:nodrg, m1:nodra, p), a_tt(m, m1:nodra, p))
                  end do
               end if
            end do
            do n = 1, nodrg
               nn1 = n * (n + 1) - n
               nn2 = nn1 + (2 * n + 1) - 1
               n1 = min(n, nodra)
               g_tt(-n:n, n, 1) = matmul(g_tt(-n1:n1, n, 1), tranmat%rot_mat(-n1:n1, nn1:nn2))
               g_tt(-n:n, n, 2) = matmul(g_tt(-n1:n1, n, 2), tranmat%rot_mat(-n1:n1, nn1:nn2))
            end do
            g_t(0, 1:nodrg, 1:2) = g_tt(0, 1:nodrg, 1:2)
            do m = 1, nodrg
               g_t(m, m:nodrg, 1:2) = g_tt(m, m:nodrg, 1:2) * tranmat%phi_mat(-im * m)
               g_t(m + 1:nodrg + 1, m, 1:2) = g_tt(-m, m:nodrg, 1:2) * tranmat%phi_mat(im * m)
            end do
            call shiftcoefficient(nodrg, nmodeg, shiftvec(1), shiftvec(2), &
                                  g_t, g_shift)
            g_t = g_shift
         else
            call shiftcoefficient(nodra, nmodea, shiftvec(1), shiftvec(2), &
                                  acoef(1:lengtha), &
                                  a_t2(1:nblka, 1:2))
            if (sop) then
               g_t2(:, 1) = matmul(a_t2(1:nblka, 1), tranmat%gen_mat(1:nblka, 1:nblkg, 1))
            else
               g_t2(:, 1) = matmul(tranmat%gen_mat(1:nblkg, 1:nblka, 1), a_t2(1:nblka, 1))
            end if
            if (nmodeg .eq. 2 .and. nmodea .eq. 1) then
               if (sop) then
                  g_t2(:, 2) = matmul(a_t2(1:nblka, 1), tranmat%gen_mat(:, :, 2))
                  g_t2 = 0.5d0 * g_t2
               else
                  g_t2(:, 2) = matmul(tranmat%gen_mat(:, :, 2), a_t2(1:nblka, 1))
                  g_t2 = 0.5d0 * g_t2
               end if
            elseif (nmodeg .eq. 1 .and. nmodea .eq. 2) then
               if (sop) then
                  g_t2(:, 1) = g_t2(:, 1) + matmul(a_t2(1:nblka, 2), tranmat%gen_mat(:, :, 2))
               else
                  g_t2(:, 1) = g_t2(:, 1) + matmul(tranmat%gen_mat(:, :, 2), a_t2(1:nblka, 2))
               end if
            elseif (nmodeg .eq. 2 .and. nmodea .eq. 2) then
               if (sop) then
                  g_t2(:, 2) = matmul(a_t2(1:nblka, 2), tranmat%gen_mat(1:nblka, 1:nblkg, 2))
               else
                  g_t2(:, 2) = matmul(tranmat%gen_mat(1:nblkg, 1:nblka, 2), a_t2(1:nblka, 2))
               end if
            end if
            call shiftcoefficient(nodrg, nmodeg, shiftvec(1), shiftvec(2), &
                                  g_t2, g_t)
         end if
         gcoef(1:lengthg) &
            = gcoef(1:lengthg) &
              + reshape(g_t(0:nodrg + 1, 1:nodrg, 1:nmodeg), (/lengthg/))
      end if
   end subroutine coefficient_translation

   subroutine shiftcoefficient(nodr, nmode, msign, mflip, &
                               ain, aout)
      implicit none(type, external)
      integer, intent(in) :: nodr, nmode, msign, mflip
      integer :: m, n, im
      complex(real64), intent(in) :: ain(0:nodr + 1, nodr, nmode)
      complex(real64), intent(out) :: aout(0:nodr + 1, nodr, nmode)
      complex(real64) :: at(nmode)
      if (msign .eq. 1 .and. mflip .eq. 1) then
         aout = ain
      else
         aout(0, 1:nodr, 1:nmode) = ain(0, 1:nodr, 1:nmode)
         if (mflip .eq. -1) then
            im = 1
            do m = 1, nodr
               im = im * msign
               do n = m, nodr
                  at = ain(n + 1, m, 1:nmode)
                  aout(n + 1, m, 1:nmode) = im * ain(m, n, 1:nmode)
                  aout(m, n, 1:nmode) = im * at
               end do
            end do
         else
            im = 1
            do m = 1, nodr
               im = im * msign
               do n = m, nodr
                  aout(m, n, 1:nmode) = im * ain(m, n, 1:nmode)
                  aout(n + 1, m, 1:nmode) = im * ain(n + 1, m, 1:nmode)
               end do
            end do
         end if
      end if
   end subroutine shiftcoefficient

end module translation
