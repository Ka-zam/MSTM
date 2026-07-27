module translation_surface_interactions
   use angular_functions, only: axial_translation_offset, axial_translation_size
   use coefficient_indexing, only: polarized_mode_index
   use iso_fortran_env, only: real64
   use mpidefs, only: mpi_comm_world, mstm_mpi, mstm_mpi_wtime
   use periodic_lattice_subroutines, only: periodic_lattice, plane_boundary_lattice_interaction
   use spheredata, only: host_sphere, number_spheres, one_side_only, recalculate_surface_matrix, &
                         run_print_unit, sphere_block, sphere_layer, sphere_offset, sphere_order, sphere_position, &
                         store_surface_matrix, store_translation_matrix
   use surface_subroutines, only: layer_ref_index, plane_interaction, plane_surface_present
   use wave_functions, only: reverse_azimuthal_modes

   implicit none(type, external)
   private

   public :: periodic_lattice_sphere_interaction
   public :: sphere_surface_interaction

   type surface_interaction_cache
      private
      logical :: symmetrical = .false.
      complex(real64), allocatable :: matrix(:)
   contains
      procedure :: configure => configure_surface_interaction_cache
      procedure :: clear => clear_surface_interaction_cache
      final :: finalize_surface_interaction_cache
   end type surface_interaction_cache

   type periodic_lattice_interaction_cache
      private
      complex(real64), allocatable :: matrix(:)
   contains
      procedure :: configure => configure_periodic_lattice_cache
      procedure :: clear => clear_periodic_lattice_interaction_cache
      final :: finalize_periodic_lattice_interaction_cache
   end type periodic_lattice_interaction_cache

   type(surface_interaction_cache), target, allocatable :: stored_surface_interactions(:)
   type(periodic_lattice_interaction_cache), target, allocatable :: stored_periodic_interactions(:)

contains

   subroutine configure_surface_interaction_cache(self, symmetrical, matrix_size)
      class(surface_interaction_cache), intent(inout) :: self
      integer, intent(in) :: matrix_size
      logical, intent(in) :: symmetrical

      call self%clear()
      self%symmetrical = symmetrical
      allocate (self%matrix(matrix_size))
   end subroutine configure_surface_interaction_cache

   subroutine clear_surface_interaction_cache(self)
      class(surface_interaction_cache), intent(inout) :: self

      if (allocated(self%matrix)) deallocate (self%matrix)
      self%symmetrical = .false.
   end subroutine clear_surface_interaction_cache

   subroutine finalize_surface_interaction_cache(self)
      type(surface_interaction_cache), intent(inout) :: self

      call self%clear()
   end subroutine finalize_surface_interaction_cache

   subroutine configure_periodic_lattice_cache(self, matrix_size)
      class(periodic_lattice_interaction_cache), intent(inout) :: self
      integer, intent(in) :: matrix_size

      call self%clear()
      allocate (self%matrix(matrix_size))
   end subroutine configure_periodic_lattice_cache

   subroutine clear_periodic_lattice_interaction_cache(self)
      class(periodic_lattice_interaction_cache), intent(inout) :: self

      if (allocated(self%matrix)) deallocate (self%matrix)
   end subroutine clear_periodic_lattice_interaction_cache

   subroutine finalize_periodic_lattice_interaction_cache(self)
      type(periodic_lattice_interaction_cache), intent(inout) :: self

      call self%clear()
   end subroutine finalize_periodic_lattice_interaction_cache

   subroutine clear_surface_cache(mat)
      implicit none(type, external)
      type(surface_interaction_cache), allocatable, intent(inout) :: mat(:)
      if (.not. allocated(mat)) return
      deallocate (mat)
   end subroutine clear_surface_cache

   subroutine clear_periodic_lattice_cache(mat)
      implicit none(type, external)
      type(periodic_lattice_interaction_cache), allocatable, intent(inout) :: mat(:)
      if (.not. allocated(mat)) return
      deallocate (mat)
   end subroutine clear_periodic_lattice_cache

   subroutine periodic_lattice_sphere_interaction(neqns, nrhs, ain, aout, &
                                                  initial_run, rhs_list, mpi_comm, con_tran, store_matrix_option)
      implicit none(type, external)
      integer, intent(in) :: neqns, nrhs
      integer :: rank, numprocs, nmat, mpicomm, proc, i, j, rhs, &
                 i1, i2, j1, j2, task, rank0, nbi, nbj, rmatdim, nmat_tot
      logical :: initrun, rhslist(nrhs), calcmat, contran(nrhs), smopt
      logical, optional, intent(in) :: initial_run, rhs_list(nrhs), con_tran(nrhs), store_matrix_option
      integer, optional, intent(in) :: mpi_comm
      real(real64) :: rp(3), time1, time2
      complex(real64), intent(inout) :: ain(neqns, nrhs)
      complex(real64), intent(out) :: aout(neqns, nrhs)
      complex(real64) :: ri
      type(periodic_lattice_interaction_cache), target :: rmat
      type(periodic_lattice_interaction_cache), pointer :: loc_rmat

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
            call clear_periodic_lattice_cache(stored_periodic_interactions)
            allocate (stored_periodic_interactions(nmat))
         end if
      end if

      task = 0
      nmat = 0
      aout = 0.0_real64
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
                        call stored_periodic_interactions(nmat)%configure(rmatdim)
                        loc_rmat => stored_periodic_interactions(nmat)
                     else
                        call rmat%configure(rmatdim)
                        loc_rmat => rmat
                     end if
                     call plane_boundary_lattice_interaction(sphere_order(i), sphere_order(j), &
                                                             rp(1), rp(2), sphere_position(3, i), sphere_position(3, j), &
                                                    loc_rmat%matrix, include_source=.true., lr_transformation=.true., index_model=2)
                     if (smopt .and. store_translation_matrix .and. rank0 .eq. 0) then
                        time2 = mstm_mpi_wtime()
                        if (time2 - time1 .ge. 15.0_real64) then
                           write (run_print_unit, '('' assembling pl matrix '',i5,''/'',i5)') &
                              nmat, nmat_tot
                           flush (run_print_unit)
                           time1 = time2
                        end if
                     end if
                  else
                     loc_rmat => stored_periodic_interactions(nmat)
                  end if
                  do rhs = 1, nrhs
                     if (.not. rhslist(rhs)) cycle
                     if (plane_surface_present) then
                        if (.not. contran(rhs)) then
                           call apply_periodic_lattice_matrix(sphere_order(i), sphere_order(j), &
                                                              ain(j1:j2, rhs), aout(i1:i2, rhs), &
                                                              .false., pb_mat=loc_rmat%matrix)
                        else
                           call apply_periodic_lattice_matrix(sphere_order(i), sphere_order(j), &
                                                              aout(j1:j2, rhs), ain(i1:i2, rhs), &
                                                              .true., pb_mat=loc_rmat%matrix)
                        end if
                     else
                        if (.not. contran(rhs)) then
                           call apply_periodic_lattice_matrix(sphere_order(i), sphere_order(j), &
                                                              ain(j1:j2, rhs), aout(i1:i2, rhs), &
                                                              .false., fs_mat=loc_rmat%matrix)
                        else
                           call apply_periodic_lattice_matrix(sphere_order(i), sphere_order(j), &
                                                              aout(j1:j2, rhs), ain(i1:i2, rhs), &
                                                              .true., fs_mat=loc_rmat%matrix)
                        end if
                     end if
                  end do
                  if (.not. (smopt .and. store_translation_matrix)) call rmat%clear()
               end if
            end if
         end do
      end do
      if (store_surface_matrix) recalculate_surface_matrix = .false.
   end subroutine periodic_lattice_sphere_interaction

   subroutine apply_periodic_lattice_matrix(nodrt, nodrs, as, at, tran, fs_mat, pb_mat)
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
   end subroutine apply_periodic_lattice_matrix

   subroutine sphere_surface_interaction(neqns, nrhs, ain, aout, &
                                         initial_run, rhs_list, mpi_comm, con_tran)
      implicit none(type, external)
      integer, intent(in) :: neqns, nrhs
      integer :: rank, numprocs, nmat, mpicomm, proc, i, j, rhs, &
                 i1, i2, j1, j2, task, jstart, rmatdim, rank0, nmat_tot
      logical :: initrun, rhslist(nrhs), calcmat, contran(nrhs), rmatsymm
      logical, optional, intent(in) :: initial_run, rhs_list(nrhs), con_tran(nrhs)
      integer, optional, intent(in) :: mpi_comm
      real(real64) :: rp(3), time1, time2
      complex(real64), intent(inout) :: ain(neqns, nrhs)
      complex(real64), intent(out) :: aout(neqns, nrhs)
      complex(real64), allocatable :: atempi(:), atempj(:), &
                                      atempi2(:), atempj2(:)
      type(surface_interaction_cache), target :: rmat
      type(surface_interaction_cache), pointer :: loc_rmat

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
            call clear_surface_cache(stored_surface_interactions)
            allocate (stored_surface_interactions(nmat))
         end if
      end if

      task = 0
      nmat = 0
      aout = 0.0_real64
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
                     rp = sphere_position(:, i) - sphere_position(:, j)
                     rmatsymm = (abs(rp(1)) .lt. 1.0e-4_real64 .and. abs(rp(2)) .lt. 1.0e-4_real64 &
                                 .and. (.not. periodic_lattice))
                     if (rmatsymm) then
                        rmatdim = 2 * axial_translation_size(sphere_order(i), sphere_order(j))
                     else
                        rmatdim = sphere_block(i) * sphere_block(j)
                     end if
                     if (store_surface_matrix) then
                        call stored_surface_interactions(nmat)%configure(rmatsymm, rmatdim)
                        loc_rmat => stored_surface_interactions(nmat)
                     else
                        call rmat%configure(rmatsymm, rmatdim)
                        loc_rmat => rmat
                     end if
                     call plane_interaction(sphere_order(i), sphere_order(j), &
                                            rp(1), rp(2), sphere_position(3, j), sphere_position(3, i), &
                                            loc_rmat%matrix, &
                                            index_model=2, lr_transformation=.true., &
                                            make_symmetric=rmatsymm)
                     if (store_surface_matrix .and. rank0 .eq. 0) then
                        time2 = mstm_mpi_wtime()
                        if (time2 - time1 .ge. 15.0_real64) then
                           write (run_print_unit, '('' assembling surf matrix '',i5,''/'',i5)') &
                              nmat, nmat_tot
                           flush (run_print_unit)
                           time1 = time2
                        end if
                     end if
                  else
                     loc_rmat => stored_surface_interactions(nmat)
                  end if
                  allocate (atempj(sphere_block(j)), atempi(sphere_block(i)))
                  if ((j .ne. i) .and. one_side_only) allocate (atempj2(sphere_block(j)), atempi2(sphere_block(i)))
                  do rhs = 1, nrhs
                     if (.not. rhslist(rhs)) cycle
                     if (.not. contran(rhs)) then
                        call apply_surface_interaction_matrix(sphere_order(j), sphere_order(i), &
                                                              ain(j1:j2, rhs), atempi, loc_rmat, 1)
                        aout(i1:i2, rhs) = aout(i1:i2, rhs) + atempi(1:sphere_block(i))
                     else
                        call apply_surface_interaction_matrix(sphere_order(i), sphere_order(j), &
                                                              ain(i1:i2, rhs), atempj, loc_rmat, 2)
                        aout(j1:j2, rhs) = aout(j1:j2, rhs) + atempj(1:sphere_block(j))
                     end if
                     if ((j .ne. i) .and. one_side_only) then
                        if (.not. contran(rhs)) then
                           call reverse_azimuthal_modes(sphere_order(i), &
                                                        ain(i1:i2, rhs), atempi)
                           call apply_surface_interaction_matrix(sphere_order(i), sphere_order(j), &
                                                                 atempi, atempj, loc_rmat, 2)
                           call reverse_azimuthal_modes(sphere_order(j), &
                                                        atempj, atempj2)
                           aout(j1:j2, rhs) = aout(j1:j2, rhs) &
                                              + atempj2(1:sphere_block(j))
                        else
                           call reverse_azimuthal_modes(sphere_order(j), &
                                                        ain(j1:j2, rhs), atempj)
                           call apply_surface_interaction_matrix(sphere_order(j), sphere_order(i), &
                                                                 atempj, atempi, loc_rmat, 1)
                           call reverse_azimuthal_modes(sphere_order(i), &
                                                        atempi, atempi2)
                           aout(i1:i2, rhs) = aout(i1:i2, rhs) + atempi2(1:sphere_block(i))
                        end if
                     end if
                  end do
                  deallocate (atempi, atempj)
                  if ((j .ne. i) .and. one_side_only) deallocate (atempi2, atempj2)
                  if (.not. store_surface_matrix) call rmat%clear()
               end if
            end if
         end do
      end do
      if (store_surface_matrix) recalculate_surface_matrix = .false.
   end subroutine sphere_surface_interaction

   subroutine apply_surface_interaction_matrix(nin, nout, ain, aout, rmat, dir)
      implicit none(type, external)
      integer, intent(in) :: nin, nout, dir
      integer :: bin, bout, m, m1, n, l, p, q, mnp, klq, i, moff
      complex(real64), intent(in) :: ain(2 * nin * (nin + 2))
      complex(real64), intent(out) :: aout(2 * nout * (nout + 2))
      type(surface_interaction_cache), intent(in) :: rmat
      bin = 2 * nin * (nin + 2)
      bout = 2 * nout * (nout + 2)
      aout = 0.0_real64
      if (rmat%symmetrical) then
         do m = -min(nin, nout), min(nin, nout)
            m1 = max(abs(m), 1)
            moff = 2 * axial_translation_offset(m, nin, nout)
            do n = m1, nout
               do p = 1, 2
                  mnp = polarized_mode_index(m, n, p, nout, 2)
                  do l = m1, nin
                     do q = 1, 2
                        klq = polarized_mode_index(m, l, q, nin, 2)
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
   end subroutine apply_surface_interaction_matrix

end module translation_surface_interactions
