!
!  module solver: subroutines for solving interaction equations for fixed orientation
!  and T matrix problems
!
!
!  last revised: 15 January 2011
!                february 2013
!
module solver
   use, intrinsic :: iso_fortran_env, only: real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use direct_lu_solver, only: direct_lu_solver_t, solver_breakdown, solver_converged, &
                               solver_iteration_limit, solver_non_finite, solver_singular, &
                               solver_status_message
   use interaction_operators, only: create_interaction_operator, interaction_operator_t
   use parallel_runtime, only: mpi_comm_null, mpi_comm_world, mstm_global_rank, parallel_allreduce_sum, &
                               parallel_barrier, parallel_broadcast, parallel_communicator_create, parallel_group, &
                               parallel_group_include, parallel_rank, parallel_reduce_sum, parallel_size, &
                               parallel_split, parallel_wall_time
   use numerical_tables
   use runtime_support, only: open_output_file, runtime_failed, set_runtime_error, synchronize_runtime_status
   use coefficient_indexing, only: polarized_mode_index
   use sphere_data
   use mie
   use translation_expansions, only: general_interaction_matrix
   use random_orientation_scattering, only: random_orientation_scattering_matrix
   use scattering_efficiencies, only: configuration_efficiency_factors
   use scattering_interactions, only: distribute_from_common_origin, merge_to_common_origin, &
                                      phase_shift, sphere_interaction, sphere_plane_wave_coefficients
   implicit none

   real(real64) :: direct_matrix_assembly_time = 0.0_real64
   real(real64) :: direct_factorization_time = 0.0_real64
   real(real64) :: direct_condition_estimation_time = 0.0_real64
   real(real64) :: direct_backsolve_time = 0.0_real64

contains

   pure real(real64) function relative_residual_norm(residual, right_hand_side)
      complex(real64), intent(in) :: residual(:), right_hand_side(:)
      real(real64) :: denominator

      denominator = norm2(abs(right_hand_side))
      if (denominator <= tiny(1.0_real64)) then
         relative_residual_norm = norm2(abs(residual))
      else
         relative_residual_norm = norm2(abs(residual)) / denominator
      end if
   end function relative_residual_norm

   pure logical function complex_vector_is_finite(values)
      complex(real64), intent(in) :: values(:)

      complex_vector_is_finite = all(ieee_is_finite(real(values, kind=real64))) .and. &
         all(ieee_is_finite(aimag(values)))
   end function complex_vector_is_finite

   subroutine solve_t_matrix(solution_method, solution_eps, convergence_eps, &
                             max_iterations, t_matrix_file, procs_per_soln, mpi_comm, &
                             sphere_qeff, solution_status, sphere_excitation_list)
      implicit none
      logical :: firstrun, initialize, itersoln, continueloop, exlist(sphere_cluster%number_spheres)
      logical, optional :: sphere_excitation_list(sphere_cluster%number_spheres)
      integer :: file_unit, ierr, io_status, iter, niter, istat, rank, maxiter, rank0, &
                 numprocs, mpicomm, i, nblkt, l, k, q, ka, la, nsolns, &
                 n, m, p, nssoln, kq, ns
      integer, save :: pcomm, prank, pgroup, ppsoln, pcomm0
      integer, optional :: mpi_comm, procs_per_soln, max_iterations, solution_status
      real(real64) :: r0(3), maxerr, qeffi(3, sphere_cluster%number_spheres), &
                      dqeffi(3, sphere_cluster%number_spheres), qeff(3), qeffold(3), time0, time1, timepersoln, timeleft, &
                      solneps, conveps, solnerr, converr(1), qteff(3), &
                      ttime(0:6), dqteff(3), converri, qeffiold(3)
      real(real64), allocatable :: sexp(:, :), scexp(:, :)
      real(real64), optional :: sphere_qeff(3, sphere_cluster%number_spheres), &
                                solution_eps, convergence_eps
      complex(real64) :: amnpkq(sphere_cluster%number_eqns), pmnpkq(sphere_cluster%number_eqns), pmnpan(sphere_cluster%number_eqns)
      complex(real64), allocatable :: pmnp0(:, :, :), amnp0(:)
      class(interaction_operator_t), allocatable :: external_operator
      type(direct_lu_solver_t) :: direct_solver
      character(len=4) :: timeunit
      character(len=128) :: tmatrixfile
      character(len=256) :: io_message
      character(len=1), optional :: solution_method
      character(len=128), optional :: t_matrix_file
      data firstrun/.true./
      istat = 0
      if (present(solution_method)) then
         itersoln = solution_method .eq. 'i'
      else
         itersoln = .true.
      end if
      if (present(solution_eps)) then
         solneps = solution_eps
      else
         solneps = 1.d-6
      end if
      if (present(convergence_eps)) then
         conveps = convergence_eps
      else
         conveps = 1.d-6
      end if
      if (present(max_iterations)) then
         niter = max_iterations
      else
         niter = 100
      end if
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call create_interaction_operator(sphere_cluster%fft_translation_option, external_operator)
      if (present(procs_per_soln)) then
         ppsoln = procs_per_soln
      else
         ppsoln = 1
      end if
      if (present(t_matrix_file)) then
         tmatrixfile = t_matrix_file
      else
         tmatrixfile = 'tmatrix_temp.dat'
      end if
      if (present(sphere_excitation_list)) then
         exlist = sphere_excitation_list
      else
         exlist = .true.
      end if

      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank0)
      if (.not. itersoln) ppsoln = 1
      ppsoln = min(ppsoln, numprocs)
      nssoln = max(1, numprocs / ppsoln)

      if (firstrun) then
         firstrun = .false.
         if (numprocs .gt. 1) then
            pgroup = floor(dble(rank) / dble(ppsoln))
            call parallel_split( &
               mpi_color=pgroup, mpi_key=rank, &
               mpi_new_comm=pcomm, &
               mpi_comm=mpicomm)
            call parallel_rank(mpi_rank=prank, mpi_comm=pcomm)
            call parallel_split( &
               mpi_color=prank, mpi_key=rank, &
               mpi_new_comm=pcomm0, &
               mpi_comm=mpicomm)
         else
            pgroup = 0
            prank = 0
            pcomm = mpicomm
            pcomm0 = mpicomm
         end if
      end if

      nblkt = 2 * sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2)
      r0 = sphere_cluster%cluster_origin(:)
      qeffi = 0.d0
      qeffold = 0.d0
      qteff = 0.d0
      nsolns = 0
      initialize = .true.
      istat = 0
      if (rank .eq. 0) then
         call open_output_file(tmatrixfile, file_unit)
         if (.not. runtime_failed()) then
            write (file_unit, '(2i4)') sphere_cluster%t_matrix_order, sphere_cluster%t_matrix_order
            time0 = parallel_wall_time()
            close (file_unit)
         end if
      end if
      call synchronize_runtime_status(mpicomm)
      if (runtime_failed()) then
         if (present(solution_status)) solution_status = 2
         return
      end if

      maxiter = 0
      maxerr = 0.d0

      if (rank0 .eq. 0) then
         write (sphere_cluster%run_print_unit, '(''  n   # its  qext         qabs'',&
                   &''           error     sec/soln est. time rem.'')')
         flush (sphere_cluster%run_print_unit)
      end if

      do l = 1, sphere_cluster%t_matrix_order
         if (rank .eq. 0) then
            ttime(0) = parallel_wall_time()
         end if
         if (rank .eq. 0) then
            time1 = parallel_wall_time()
         end if

         allocate (pmnp0(0:l + 1, l, 2), amnp0(2 * l * (l + 2)))
         kq = 0
         continueloop = .true.
         do while (continueloop)
            ns = 0
            dqeffi = 0.d0
            dqteff = 0.d0
            do i = 1, nssoln
               ns = ns + 1
               kq = kq + 1
               k = -l + (kq - 1) / 2
               q = mod(kq - 1, 2) + 1
               if ((i - 1) .eq. pgroup) then
!if(prank.eq.0) then
!write(*,'('' pgroup,l,k,q:'',4i5)') pgroup,l,k,q
!flush(6)
!endif
                  pmnp0 = 0.d0
                  if (k .le. -1) then
                     ka = l + 1
                     la = -k
                  else
                     ka = k
                     la = l
                  end if
!
!  the t matrix is te-tm based; hence the following two lines
!
                  pmnp0(ka, la, 1) = .5d0
                  pmnp0(ka, la, 2) = -.5d0 * (-1)**q
                  call distribute_from_common_origin(l, pmnp0, pmnpkq, &
                                                     number_rhs=1, &
                                                     origin_position=r0, &
                                                     origin_host=0, &
                                                     vswf_type=1, &
                                                     mpi_comm=pcomm)
!                        sphere_translation_list=exlist)

                  call apply_mie_coefficients(sphere_cluster%number_eqns, 1, 1, pmnpkq, pmnpan)
                  amnpkq = pmnpan
                  if (niter .ne. 0) then
                     if (itersoln) then
                        call solve_complex_biconjugate_gradient(niter, solneps, pmnpan, amnpkq, 0, &
                                                                iter, solnerr, initialize_solver=initialize, &
                                                                mpi_comm=pcomm, solution_status=ierr, &
                                                                external_operator=external_operator)
                        maxiter = max(iter, maxiter)
                        maxerr = max(solnerr, maxerr)
                        if (ierr /= solver_converged) then
                           istat = ierr
                           deallocate (pmnp0, amnp0)
                           if (present(solution_status)) solution_status = istat
                           return
                        end if
                     else
                        call solve_direct_system(direct_solver, pmnpan, amnpkq, &
                                                 initialize_solver=initialize, solution_error=solnerr, &
                                                 number_iterations=iter, solution_eps=solneps, &
                                                 solution_status=ierr, mpi_comm=pcomm)
                        maxiter = max(iter, maxiter)
                        maxerr = max(solnerr, maxerr)
                        if (ierr /= solver_converged) then
                           istat = ierr
                           deallocate (pmnp0, amnp0)
                           if (present(solution_status)) solution_status = istat
                           return
                        end if
                     end if
                  end if
                  initialize = .false.
                  call configuration_efficiency_factors(sphere_cluster%number_spheres, 1, amnpkq, &
                                                        pmnpkq, dqeffi, mpi_comm=pcomm)
!                     dqeffi(3,:)=dqeffi(1,:)-dqeffi(2,:)
! patch 10-22.  used new formula for qeff(3) in qefficiencyfactor SR.  Assumes ri_medium is real
!                     dqeffi(3,:)=0.d0
                  amnp0 = 0.d0
                  call merge_to_common_origin(l, amnpkq, &
                                              amnp0, &
                                              number_rhs=1, &
                                              origin_position=r0, &
                                              mpi_comm=pcomm, &
                                              sphere_translation_list=exlist)
                  ka = polarized_mode_index(k, l, q, l, 2)
                  dqteff(1) = -2.d0 / sphere_cluster%vol_radius**2 * dble(amnp0(ka))
                  dqteff(3) = 2.d0 / sphere_cluster%vol_radius**2 &
                              * dble(sum(amnp0(:) * conjg(amnp0(:))))
                  dqteff(2) = dqteff(1) - dqteff(3)

               end if

               if (kq .eq. 2 * (2 * l + 1)) then
                  continueloop = .false.
                  exit
               end if
            end do
            if (rank .eq. 0) then
               ttime(1) = parallel_wall_time()
            end if

            if (prank .eq. 0) then
               call parallel_allreduce_sum(receive_buffer=dqteff, &
                                           mpi_number=3, mpi_comm=pcomm0)
               call parallel_allreduce_sum(receive_buffer=dqeffi, &
                                           mpi_number=3 * sphere_cluster%number_spheres, mpi_comm=pcomm0)

               qeffi = qeffi + dqeffi
               qteff = qteff + dqteff

               do i = 1, ns
                  if (i - 1 .eq. pgroup) then
                     call open_output_file(tmatrixfile, file_unit, append=.true.)
                     if (.not. runtime_failed()) then
                        do n = 1, l
                           do m = -n, n
                              do p = 1, 2
                                 ka = polarized_mode_index(m, n, p, l, 2)
                                 write (file_unit, '(''('',e18.10,'','',e18.10,'')'')') amnp0(ka)
                              end do
                           end do
                        end do
                        close (file_unit)
                     end if
                  end if
                  call parallel_barrier(mpi_comm=pcomm0)
               end do
            end if
         end do

!if(rank.eq.0) then
!ttime(2)=parallel_wall_time()
!write(*,'('' timings:'',2es14.5)') ttime(1)-ttime(0),ttime(2)-ttime(1)
!endif

         if (rank .eq. 0) then
            qeff = 0.d0
            do i = 1, sphere_cluster%number_spheres
               if (.not. exlist(i)) cycle
               qeff(:) = qeff(:) + qeffi(:, i) * sphere_cluster%sphere_radius(i)**2 / sphere_cluster%vol_radius**2
            end do

            converr = qeff(1) - qeffold(1)
            qeffold = qeff
            nsolns = nsolns + (l + l + 1) * 2

            if (rank0 .eq. 0) then
               timepersoln = (parallel_wall_time() - time1) / dble((l + l + 1) * 2)
               timeleft = timepersoln * (nblkt - nsolns)
               if (timeleft .gt. 3600.d0) then
                  timeleft = timeleft / 3600.d0
                  timeunit = ' hrs'
               elseif (timeleft .gt. 60.d0) then
                  timeleft = timeleft / 60.d0
                  timeunit = ' min'
               else
                  timeunit = ' sec'
               end if
               write (sphere_cluster%run_print_unit, '(i4,i5,4e13.5,2f8.2,a4)') l, maxiter, qeff(1:2), &
                  converr, dqeffi(1, 1), timepersoln, timeleft, timeunit
               flush (sphere_cluster%run_print_unit)
            end if
         end if
         deallocate (amnp0, pmnp0)

!            allocate(sexp(16,0:2*l),scexp(16,0:2*l))
!            call random_orientation_scattering_matrix(tmatrixfile,sexp, &
!               scexp,override_order=l, &
!               keep_quiet=.true., &
!               mpi_comm=mpicomm)
!            if(rank.eq.0) then
!               if(l.eq.1) then
!                  open(20,file='tmatsmexp.dat')
!               else
!                  open(20,file='tmatsmexp.dat',position='append')
!               endif
!               write(20,'(i5)') l
!               do n=0,2*l
!                  write(20,'(i5,2es14.6)') n,sexp(1,n),scexp(1,n)
!               enddo
!               close(20)
!            endif
!            deallocate(sexp,scexp)

         call parallel_broadcast(mpi_rank=0, send_buffer=converr, mpi_number=1, mpi_comm=mpicomm)

         if (conveps .gt. 0.d0 .and. converr(1) .lt. conveps) exit
      end do

      if (l .lt. sphere_cluster%t_matrix_order) then
         sphere_cluster%t_matrix_order = l
         if (rank .eq. 0) then
            open (newunit=file_unit, file=tmatrixfile, status='old', action='write', &
                  form='formatted', access='direct', recl=8, iostat=io_status, iomsg=io_message)
            if (io_status /= 0) then
               call set_runtime_error("Cannot update T-matrix file '"//trim(tmatrixfile)//"': "//trim(io_message), &
                                      io_status)
            else
               write (file_unit, '(2i4)', rec=1) l, l
               close (file_unit)
            end if
         end if
      end if

      call synchronize_runtime_status(mpicomm)
      if (runtime_failed()) istat = 2

      if (present(sphere_qeff)) sphere_qeff = qeffi
      if (present(solution_status)) solution_status = istat

   end subroutine solve_t_matrix

!
!  solution of interaction equations for a fixed orientation
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: modification of efficiency calculation, to calculate
!           polarized components
!  30 March 2011: took out gbfocus argument: this is not needed since positions are defined
!  relative to the gb focus.
!  20 April 2011: used 2-group MPI formulation
!  October 2011: adapted to far field approximation.
!  December 2011: changed efficiency factor calculation, adapted to generalized sphere
!                 configuration.
! february 2013: number rhs and mpi comm options added, completely rewritten.
!
!
   subroutine solve_fixed_orientation(alpha, sinc, dir, eps, niter, amnp, qeff, &
                                      qeffdim, maxerr, maxiter, iterwrite, istat, &
                                      mpi_comm, excited_spheres, solution_method, initialize_solver, &
                                      reciprocal_condition, incident_coefficients, number_rhs)
      implicit none
      logical :: firstrun, exsphere(sphere_cluster%number_spheres), dirsoln, initialize
      logical, save :: inp1, inp2
      logical, optional :: excited_spheres(sphere_cluster%number_spheres), initialize_solver
      integer :: ierr, iter, niter, istat, rank, maxiter, iterwrite, nsend, &
                 numprocs, mpicomm, prank, oddnumproc, &
                 groupsize, pgroup, mpigroup, syncgroup, i, p, dir, qeffdim, nrhs
      integer, save :: pcomm, synccomm1, synccomm2, p1, p2
      integer, allocatable :: grouplist(:)
      integer, optional :: mpi_comm, number_rhs
      real(real64) :: alpha, sinc, eps, serr, qeff(3, qeffdim, sphere_cluster%number_spheres), maxerr, direct_condition
      real(real64), optional, intent(out) :: reciprocal_condition
      complex(real64) :: amnp(sphere_cluster%number_eqns, 2)
      complex(real64), optional, intent(in) :: incident_coefficients(sphere_cluster%number_eqns, 2)
      complex(real64), allocatable :: pmnpan(:), pmnp0(:, :)
      class(interaction_operator_t), allocatable :: external_operator
      type(direct_lu_solver_t) :: direct_solver
      character(len=1), optional :: solution_method
      data firstrun/.true./
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      if (present(number_rhs)) then
         nrhs = number_rhs
      else
         nrhs = 2
      end if
      call create_interaction_operator(sphere_cluster%fft_translation_option, external_operator)
      if (present(excited_spheres)) then
         exsphere = excited_spheres
      else
         exsphere = .true.
      end if
      if (present(solution_method)) then
         dirsoln = (solution_method .ne. 'i')
      else
         dirsoln = .false.
      end if
      if (present(initialize_solver)) then
         initialize = initialize_solver
      else
         initialize = firstrun
      end if

!         if(mpicomm.eq.mpi_comm_null) return
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)

      if (initialize) then
         if (numprocs .gt. 1 .and. (.not. dirsoln) .and. nrhs .eq. 2) then
            oddnumproc = mod(numprocs, 2)
            pgroup = floor(dble(2 * rank) / dble(numprocs)) + 1
            p1 = pgroup
            p2 = p1
            call parallel_split( &
               mpi_color=pgroup, mpi_key=rank, &
               mpi_new_comm=pcomm, &
               mpi_comm=mpicomm)
            call parallel_rank(mpi_rank=prank, mpi_comm=pcomm)
            call parallel_group(mpi_group=mpigroup, mpi_comm=mpicomm)
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
            call parallel_group_include( &
               mpi_group=mpigroup, &
               mpi_size=groupsize, &
               mpi_new_group_list=grouplist, &
               mpi_new_group=syncgroup)
            call parallel_communicator_create( &
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
            call parallel_group_include( &
               mpi_group=mpigroup, &
               mpi_size=groupsize, &
               mpi_new_group_list=grouplist, &
               mpi_new_group=syncgroup)
            call parallel_communicator_create( &
               mpi_group=syncgroup, &
               mpi_comm=mpicomm, &
               mpi_new_comm=synccomm2)
            deallocate (grouplist)
         else
            p1 = 1
            p2 = nrhs
            pcomm = mpicomm
         end if
      end if
      firstrun = .false.

!write(*,*) ' step 3'
!flush(6)
      if (light_up) then
         write (*, '('' s8.2.1 '',i3)') mstm_global_rank
         flush (6)
      end if
!call parallel_barrier()

      allocate (pmnpan(sphere_cluster%number_eqns), pmnp0(sphere_cluster%number_eqns, 2))
      pmnp0 = (0.0_real64, 0.0_real64)
      if (present(reciprocal_condition)) reciprocal_condition = 1.0_real64
      if (present(incident_coefficients)) then
         pmnp0(:, 1:nrhs) = incident_coefficients(:, 1:nrhs)
      else
         call sphere_plane_wave_coefficients(alpha, sinc, dir, pmnp0, &
                                             excited_spheres=exsphere, mpi_comm=mpicomm)
      end if
!if(phase_shift_form) call phase_shift(pmnp0,1)

      do p = p1, p2
         istat = 0
         maxiter = 0
         maxerr = 0.
!
!  calculate the two solutions
!
         if (light_up) then
            write (*, '('' s8.2.2 '',i3)') mstm_global_rank
            flush (6)
         end if
         call apply_mie_coefficients(sphere_cluster%number_eqns, 1, 1, pmnp0(:, p), pmnpan)
         amnp(:, p) = pmnpan
         if (niter .ne. 0) then
!write(*,*) 'step 2, p:',p
!flush(6)
            if (light_up) then
               write (*, '('' s8.2.3 '',i3)') mstm_global_rank
               flush (6)
            end if
            if (dirsoln) then
               call solve_direct_system(direct_solver, pmnpan, amnp(:, p), initialize_solver=initialize, &
                                        number_iterations=iter, solution_error=serr, solution_eps=eps, &
                                        solution_status=ierr, mpi_comm=mpicomm, &
                                        reciprocal_condition=direct_condition)
               if (present(reciprocal_condition)) &
                  reciprocal_condition = min(reciprocal_condition, direct_condition)
               if (ierr /= solver_converged) then
                  istat = ierr
                  deallocate (pmnp0, pmnpan)
                  return
               end if
            else
               call solve_complex_biconjugate_gradient(niter, eps, pmnpan, amnp(:, p), iterwrite, &
                                                       iter, serr, mpi_comm=pcomm, initialize_solver=initialize, &
                                                       solution_status=ierr, external_operator=external_operator)
               if (ierr /= solver_converged) then
                  istat = ierr
                  deallocate (pmnp0, pmnpan)
                  return
               end if
            end if
         else
            iter = 0
            serr = 0.d0
            if (sphere_cluster%number_host_spheres .gt. 0) then
               pmnpan = 0.d0
               call sphere_interaction(sphere_cluster%number_eqns, 1, amnp(:, p), pmnpan, &
                                       initial_run=initialize, &
                                       skip_external_translation=.true., &
                                       mpi_comm=pcomm, external_operator=external_operator)
               call sphere_interaction(sphere_cluster%number_eqns, 1, pmnpan, pmnpan, &
                                       skip_external_translation=.true., &
                                       mpi_comm=pcomm, external_operator=external_operator)
               amnp(:, p) = amnp(:, p) + pmnpan
            end if
         end if
         initialize = .false.
         maxiter = max(iter, maxiter)
         maxerr = max(serr, maxerr)
      end do

!         call parallel_barrier()

      if (numprocs .gt. 1 .and. (.not. dirsoln) .and. nrhs .eq. 2) then
         nsend = sphere_cluster%number_eqns
         if (inp1) then
            call parallel_broadcast( &
               send_buffer=amnp(1:nsend, 1), &
               mpi_number=nsend, &
               mpi_rank=0, &
               mpi_comm=synccomm1)
         end if
!            call parallel_barrier()
         if (inp2) then
            call parallel_broadcast( &
               send_buffer=amnp(1:nsend, 2), &
               mpi_number=nsend, &
               mpi_rank=0, &
               mpi_comm=synccomm2)
         end if
      end if
!         call parallel_barrier()
!
!  efficiency factor calculations
!
      if (light_up) then
         write (*, '('' s8.2.4 '',i3)') mstm_global_rank
         flush (6)
      end if
!
! what is this doing here?
!call parallel_barrier()
      if (qeffdim .eq. 1) then
         i = 1
      else
         i = nrhs
      end if
      call configuration_efficiency_factors(sphere_cluster%number_spheres, i, amnp(:, 1:i), pmnp0(:, 1:i), &
                                            qeff(:, 1:2 * i - 1, :), mpi_comm=mpicomm)

!if(phase_shift_form) call phase_shift(amnp,-1)

      deallocate (pmnp0, pmnpan)
   end subroutine solve_fixed_orientation

   subroutine solve_direct_system(direct_solver, pnp, anp, initialize_solver, solution_error, &
                                  number_iterations, solution_eps, solution_status, mpi_comm, &
                                  max_refinements, reciprocal_condition)
      implicit none(type, external)
      class(direct_lu_solver_t), intent(inout) :: direct_solver
      logical :: initialize
      logical, optional, intent(in) :: initialize_solver
      integer :: rank, refinement, refinement_limit, mpicomm, status
      integer :: integer_buffer(2)
      integer, optional, intent(out) :: number_iterations, solution_status
      integer, optional, intent(in) :: mpi_comm, max_refinements
      real(real64) :: residual, seps, phase_start
      real(real64) :: real_buffer(2)
      real(real64), optional, intent(out) :: solution_error, reciprocal_condition
      real(real64), optional, intent(in) :: solution_eps
      complex(real64), intent(in) :: pnp(sphere_cluster%number_eqns)
      complex(real64), intent(out) :: anp(sphere_cluster%number_eqns)
      complex(real64), allocatable :: interaction_matrix(:, :)

      mpicomm = mpi_comm_world
      if (present(mpi_comm)) mpicomm = mpi_comm
      initialize = .not. direct_solver%is_factorized()
      if (present(initialize_solver)) initialize = initialize_solver
      seps = 1.0e-12_real64
      if (present(solution_eps)) seps = max(solution_eps, 0.0_real64)
      refinement_limit = 3
      if (present(max_refinements)) refinement_limit = max(0, max_refinements)
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)

      status = solver_converged
      refinement = 0
      residual = huge(1.0_real64)
      anp = 0.0_real64

      if (initialize) then
         direct_matrix_assembly_time = 0.0_real64
         direct_factorization_time = 0.0_real64
         direct_condition_estimation_time = 0.0_real64
         direct_backsolve_time = 0.0_real64
         allocate (interaction_matrix(sphere_cluster%number_eqns, sphere_cluster%number_eqns))
         pl_error_codes = 0
         phase_start = parallel_wall_time()
         call general_interaction_matrix(interaction_matrix, mie_mult=.true., mpi_comm=mpicomm)
         if (rank == 0) direct_matrix_assembly_time = parallel_wall_time() - phase_start
         call parallel_reduce_sum(mpi_rank=0, &
                                  receive_buffer=pl_error_codes, mpi_number=6, mpi_comm=mpicomm)
         if (rank == 0) then
            if (any(pl_error_codes /= 0)) then
               write (sphere_cluster%run_print_unit, '('' direct interaction-matrix error codes:'',10i4)') &
                  pl_error_codes, pl_rs_imax
            end if
            call direct_solver%factor(interaction_matrix, status)
            direct_factorization_time = direct_solver%factorization_time()
            direct_condition_estimation_time = direct_solver%condition_estimation_time()
         end if
         deallocate (interaction_matrix)
         integer_buffer(1) = status
         call parallel_broadcast(send_buffer=integer_buffer(1:1), mpi_number=1, mpi_rank=0, mpi_comm=mpicomm)
         status = integer_buffer(1)
      end if

      if (status == solver_converged .and. rank == 0) then
         call direct_solver%solve(pnp, anp, seps, refinement_limit, refinement, residual, status)
         direct_backsolve_time = direct_solver%backsolve_time()
      end if

      integer_buffer = [status, refinement]
      real_buffer = [residual, direct_solver%condition_estimate()]
      call parallel_broadcast(send_buffer=integer_buffer, mpi_number=2, mpi_rank=0, mpi_comm=mpicomm)
      call parallel_broadcast(send_buffer=real_buffer, mpi_number=2, mpi_rank=0, mpi_comm=mpicomm)
      call parallel_broadcast(send_buffer=anp, mpi_number=sphere_cluster%number_eqns, mpi_rank=0, mpi_comm=mpicomm)
      status = integer_buffer(1)
      refinement = integer_buffer(2)
      residual = real_buffer(1)

      if (present(solution_status)) solution_status = status
      if (present(solution_error)) solution_error = residual
      if (present(number_iterations)) number_iterations = refinement
      if (present(reciprocal_condition)) reciprocal_condition = real_buffer(2)
   end subroutine solve_direct_system
!
! iteration solver
! generalized complex biconjugate gradient method
! original code: Piotr Flatau, although not much remains.
! specialized to the multiple sphere problem
!
!
!  last revised: 15 January 2011
!  october 2011: translation calls modified
!  february 2013: number rhs option added, completely rewritten.
!
   subroutine solve_complex_biconjugate_gradient(niter, eps, pnp, anp, iterwrite, iter, errmax, &
                                                 initialize_solver, mpi_comm, solution_status, external_operator)
      implicit none
      logical :: firstrun, initialize, contran2(2)
      logical, save :: inp1, inp2
      logical, optional, intent(in) :: initialize_solver
      integer :: neqns, niter, iter, writetime, &
                 rank, iunit, iterwrite, numprocs, i, oddnumproc, &
                 mpicomm, prank, mpigroup, syncgroup, groupsize, rank0, status
      integer, save :: pgroup, pcomm, synccomm1, synccomm2
      integer, allocatable :: grouplist(:)
      integer, optional, intent(in) :: mpi_comm
      integer, optional, intent(out) :: solution_status
      class(interaction_operator_t), intent(inout) :: external_operator
      real(real64) :: eps, time1, time2, eerr, enorm, errmax, errmin, time0, &
                      breakdown_scale, residual_squared
      complex(real64) :: pnp(sphere_cluster%number_eqns), anp(sphere_cluster%number_eqns), cak, csk, cbk, csk2
      complex(real64), allocatable :: cr(:), cp(:), cw(:), cq(:), cap(:), caw(:), &
                                      capt(:), cawt(:), ctin(:, :), ctout(:, :)
      data firstrun/.true./
      data writetime/0/

      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      if (present(initialize_solver)) then
         initialize = initialize_solver
      else
         initialize = .true.
      end if
      iunit = sphere_cluster%run_print_unit
      neqns = sphere_cluster%number_eqns
      iter = 0
      errmax = 0.0_real64
      status = solver_converged
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
!         call parallel_rank(mpi_rank=rank0)
      rank0 = mstm_global_rank

      if (light_up) then
         write (*, '('' s8.2.3.1 '',i3)') mstm_global_rank
         flush (6)
      end if
      residual_squared = real(dot_product(pnp, pnp), kind=real64)
      if (.not. ieee_is_finite(residual_squared) .or. .not. complex_vector_is_finite(pnp)) then
         status = solver_non_finite
         if (present(solution_status)) solution_status = status
         return
      end if
      if (residual_squared <= tiny(1.0_real64)) then
         if (present(solution_status)) solution_status = status
         return
      end if
      allocate (cr(neqns), cp(neqns), cw(neqns), cq(neqns), &
                cap(neqns), caw(neqns), capt(neqns), &
                cawt(neqns))

      if (firstrun) then
         firstrun = .false.
         if (numprocs .gt. 1 .and. niter .gt. 0) then
            oddnumproc = mod(numprocs, 2)
            pgroup = floor(dble(2 * rank) / dble(numprocs)) + 1
            call parallel_split( &
               mpi_color=pgroup, mpi_key=rank, &
               mpi_new_comm=pcomm, &
               mpi_comm=mpicomm)
            call parallel_rank(mpi_rank=prank, mpi_comm=pcomm)
            call parallel_group(mpi_group=mpigroup, mpi_comm=mpicomm)
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
            call parallel_group_include( &
               mpi_group=mpigroup, &
               mpi_size=groupsize, &
               mpi_new_group_list=grouplist, &
               mpi_new_group=syncgroup)
            call parallel_communicator_create( &
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
            call parallel_group_include( &
               mpi_group=mpigroup, &
               mpi_size=groupsize, &
               mpi_new_group_list=grouplist, &
               mpi_new_group=syncgroup)
            call parallel_communicator_create( &
               mpi_group=syncgroup, &
               mpi_comm=mpicomm, &
               mpi_new_comm=synccomm2)
            deallocate (grouplist)
         else
            pcomm = mpicomm
         end if
      end if

      if (sphere_cluster%normalize_solution_error) then
         enorm = sqrt(residual_squared)
      else
         enorm = 1.d0
      end if
      errmin = 1.d10

      if (niter .eq. 0) then
         cp = anp
         cr = 0.
         call sphere_interaction(neqns, 1, cp, cr, &
                                 initial_run=initialize, &
                                 skip_external_translation=.true., &
                                 mpi_comm=mpicomm, external_operator=external_operator)
         cp = anp
         call sphere_interaction(neqns, 1, cp, anp, &
                                 skip_external_translation=.true., &
                                 mpi_comm=mpicomm, external_operator=external_operator)
         deallocate (cr, cp, cw, cq, cap, caw, capt, cawt)
         if (present(solution_status)) solution_status = status
         return
      end if
!
!  setting niter < 0 runs the following simple order--of--scattering solution
!
      if (niter .lt. 0) then
         cp = anp
         if (rank0 .eq. 0) time1 = parallel_wall_time()
         do iter = 1, -niter
            if (rank0 .eq. 0) time0 = parallel_wall_time()
            cr = 0.
            if (iter .eq. 1) then
               call sphere_interaction(neqns, 1, cp, cr, &
                                       initial_run=initialize, &
                                       mpi_comm=mpicomm, external_operator=external_operator)
            else
               call sphere_interaction(neqns, 1, cp, cr, &
                                       mpi_comm=mpicomm, external_operator=external_operator)
            end if
            anp = anp + cr
            cp = cr
            eerr = norm2(abs(cr)) / enorm
            errmax = eerr
            if (.not. ieee_is_finite(eerr)) then
               status = solver_non_finite
               exit
            end if
            if (eerr .lt. eps) then
               exit
            end if
            if (rank0 .eq. 0) time2 = parallel_wall_time()
            if (rank0 .eq. 0 .and. iterwrite .eq. 1 .and. time2 - time1 .gt. 5.d0) then
               write (iunit, '('' iter,err,tpi:'',i5,2e13.5)') iter, eerr, time2 - time0
               flush (iunit)
               time1 = time2
            end if
!if(rank0.eq.0) then
!write(*,'(i5,6es20.12)') iter,sum(abs(anp)),anp(sphere_cluster%number_eqns)
!endif
         end do
         call parallel_barrier(mpi_comm=mpicomm)
         deallocate (cr, cp, cw, cq, cap, caw, capt, cawt)
         if (present(solution_status)) solution_status = status
         return
      end if
!
! the following is the implementation of the complex biconjugate gradient
! iteration scheme
!
      cr = 0.d0
      call sphere_interaction(neqns, 1, anp, cr, initial_run=initialize, &
                              mpi_comm=pcomm, external_operator=external_operator)
      cr = pnp - anp + cr
      cq = conjg(cr)
      cw = cq
      cp = cr
      eerr = norm2(abs(cr)) / enorm
      errmax = eerr
      if (.not. ieee_is_finite(eerr) .or. .not. complex_vector_is_finite(cr)) then
         status = solver_non_finite
         deallocate (cr, cp, cw, cq, cap, caw, capt, cawt)
         if (present(solution_status)) solution_status = status
         return
      elseif (eerr <= eps) then
         deallocate (cr, cp, cw, cq, cap, caw, capt, cawt)
         if (present(solution_status)) solution_status = status
         return
      end if
      csk = dot_product(conjg(cr), cr)
      breakdown_scale = sqrt(real(dot_product(cq, cq), kind=real64) * &
                             real(dot_product(cr, cr), kind=real64))
      if (abs(csk) <= max(tiny(1.0_real64), 100.0_real64 * epsilon(1.0_real64) * breakdown_scale)) then
         status = solver_breakdown
         deallocate (cr, cp, cw, cq, cap, caw, capt, cawt)
         if (present(solution_status)) solution_status = status
         return
      end if
!
!  here starts the main iteration loop
!
      if (rank0 .eq. 0) time1 = parallel_wall_time()
      if (light_up) then
         write (*, '('' s8.2.3.2 '',i3)') mstm_global_rank
         flush (6)
      end if
      status = solver_iteration_limit
      do iter = 1, niter
         call parallel_barrier(mpi_comm=mpicomm)
         if (rank0 .eq. 0) time0 = parallel_wall_time()
         cak = 0.d0
         cawt = 0.d0
         capt = 0.d0
         cap = 0.d0
         caw = 0.d0

         if (light_up) then
            write (*, '('' s8.2.3.3 '',3i3)') mstm_global_rank, iter
            flush (6)
         end if
         if (numprocs .eq. 1) then
            allocate (ctin(neqns, 2), ctout(neqns, 2))
            ctin(:, 1) = cp(:)
            ctin(:, 2) = cw(:)
            contran2 = (/.false., .true./)
            ctout = 0.d0
            call sphere_interaction(neqns, 2, ctin, ctout, &
                                    con_tran=contran2, mpi_comm=mpicomm, &
                                    external_operator=external_operator)
            cap(:) = ctout(:, 1)
            caw(:) = ctout(:, 2)
            deallocate (ctin, ctout)
         else
            if (pgroup .eq. 1) then
               call sphere_interaction(neqns, 1, cp, cap, &
                                       mpi_comm=pcomm, external_operator=external_operator)
            else
               call sphere_interaction(neqns, 1, cw, caw, &
                                       con_tran=(/.true./), mpi_comm=pcomm, &
                                       external_operator=external_operator)
            end if
            call parallel_barrier(mpi_comm=mpicomm)
            if (inp2) then
               call parallel_broadcast( &
                  send_buffer=caw, &
                  mpi_number=neqns, &
                  mpi_rank=0, &
                  mpi_comm=synccomm2)
            end if
            if (inp1) then
               call parallel_broadcast( &
                  send_buffer=cap, &
                  mpi_number=neqns, &
                  mpi_rank=0, &
                  mpi_comm=synccomm1)
            end if
         end if

         call parallel_barrier(mpi_comm=mpicomm)

         if (light_up) then
            write (*, '('' s8.2.3.4 '',3i3)') mstm_global_rank, iter
            flush (6)
         end if
         cap = cp - cap
         caw = cw - caw
         cak = dot_product(cw, cap)
         breakdown_scale = sqrt(real(dot_product(cw, cw), kind=real64) * &
                                real(dot_product(cap, cap), kind=real64))
         if (.not. ieee_is_finite(abs(cak)) .or. &
             abs(cak) <= max(tiny(1.0_real64), 100.0_real64 * epsilon(1.0_real64) * breakdown_scale)) then
            status = solver_breakdown
            exit
         end if
         cak = csk / cak

         anp = anp + cak * cp
!if(rank0.eq.0) then
!write(*,'(i5,6es20.12)') iter,sum(abs(anp)),anp(sphere_cluster%number_eqns),cak
!endif
         cr = cr - cak * cap
         cq = cq - conjg(cak) * caw
         csk2 = dot_product(cq, cr)
         eerr = norm2(abs(cr)) / enorm
         errmax = eerr
         errmin = min(errmin, errmax)

         if (.not. ieee_is_finite(eerr) .or. .not. complex_vector_is_finite(anp) .or. &
             .not. complex_vector_is_finite(cr)) then
            status = solver_non_finite
            exit
         elseif (eerr <= eps) then
            status = solver_converged
            exit
         end if

         breakdown_scale = sqrt(real(dot_product(cq, cq), kind=real64) * &
                                real(dot_product(cr, cr), kind=real64))
         if (abs(csk2) <= max(tiny(1.0_real64), 100.0_real64 * epsilon(1.0_real64) * breakdown_scale)) then
            status = solver_breakdown
            exit
         end if
         cbk = csk2 / csk
         csk = csk2

         cp = cr + cbk * cp
         cw = cq + conjg(cbk) * cw

         if (light_up) then
            write (*, '('' s8.2.3.5 '',3i3)') mstm_global_rank, iter
            flush (6)
         end if
         if (rank0 .eq. 0) time2 = parallel_wall_time()
         if (rank0 .eq. 0 .and. iterwrite .eq. 1 .and. time2 - time1 .gt. 5.d0) then
            write (sphere_cluster%run_print_unit, '('' iter,err,min err, tpi:'',i5,2e12.4,e12.4)') &
               iter, errmax, errmin, time2 - time0
            flush (iunit)
            time1 = time2
         end if
      end do
      if (status == solver_iteration_limit) iter = niter
      deallocate (cr, cp, cw, cq, cap, caw, capt, cawt)
      if (present(solution_status)) solution_status = status
   end subroutine solve_complex_biconjugate_gradient

end module solver
