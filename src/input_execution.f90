module input_execution
   use, intrinsic :: iso_fortran_env, only: real64
   use angular_functions, only: cartesian_vectors_to_spherical, estimate_translation_order, &
                                rotate_expansion_coefficients
   use configuration_data
   use constants
   use effective_medium_analysis
   use fft_translation, only: fft_plan
   use input_parser
   use input_reporting
   use parallel_runtime, only: mpi_comm_world, mstm_global_rank, parallel_barrier, parallel_broadcast, &
                               parallel_rank, parallel_reduce_sum, parallel_size, parallel_split, parallel_wall_time
   use runtime_support, only: open_output_file, runtime_failed, set_runtime_error
   use random_orientation_scattering, only: evaluate_random_orientation_scattering_matrix
   use scattering_amplitudes, only: fixed_orientation_scattering_matrix_expansion, periodic_lattice_scattering
   use scattering_efficiencies, only: boundary_extinction, common_origin_hemispherical_scattering, &
                                      hemispherical_scattering, total_efficiency_factors
   use scattering_matrix_driver
   use solver, only: solver_converged, solver_status_message
   use wave_functions, only: compose_group_filename
   implicit none
contains

   subroutine execute_simulation(print_output, set_t_matrix_order, dry_run, mpi_comm)
      implicit none
      logical :: stopit, singleorigin, iframe, sett, printout, dryrun, averagerun
      logical, optional :: print_output, set_t_matrix_order, dry_run
      integer :: file_unit, n, istat, niter, rank, numprocs, i, nodrw, celldim(3), itemp(6), sx, sy, maxt, &
                 mpicomm, lochost
      integer, optional :: mpi_comm
      real(real64) :: alpha, time1, r0(3), rtran, costheta, &
                      csca, zext, targetvol, timet, tmin(3), tmax(3), rannum
      complex(real64) :: rimedium(2)
      character(len=256) :: timatrixfile
      if (present(dry_run)) then
         dryrun = dry_run
      else
         dryrun = .false.
      end if
      if (present(print_output)) then
         printout = print_output
      else
         printout = .true.
      end if
      if (present(set_t_matrix_order)) then
         sett = set_t_matrix_order
      else
         sett = .true.
      end if
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      averagerun = configuration_average .or. incidence_average
      calculate_up_down_scattering = input_calculate_up_down_scattering
      if (reflection_model) then
         calculate_up_down_scattering = .true.
         incident_frame = .false.
      end if

      first_run = .false.
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      if (solution_method(1:1) /= 'i' .and. solution_method(1:1) /= 'd') then
         call set_runtime_error("Unknown solution method '"//trim(solution_method)//"'; use iteration or direct")
         return
      end if
      if (max_iterations > 0 .and. solution_epsilon <= 0.0_real64) then
         call set_runtime_error('solution_epsilon must be positive')
         return
      end if
!         if(rank.ne.0) light_up=.false.
      local_rank = rank
      global_rank = rank
      if ((.not. configuration_average) .and. (.not. incidence_average)) n_configuration_groups = 1
!         random_configuration=(trim(sphere_data_input_file).eq.'random_configuration')
      if (random_configuration) then
         if (auto_target_radius .and. target_shape .eq. 2) then
            sphere_cluster%number_spheres = input_number_spheres
            target_dimensions(1:3) = (dble(sphere_cluster%number_spheres) / sphere_volume_fraction)**(1.d0 / 3.d0)
         else
            if (target_width_specified) then
               if (target_shape .eq. 0) then
                  target_dimensions(1:2) = target_width
                  target_dimensions(3) = target_thickness
               elseif (target_shape .eq. 1) then
                  target_dimensions(1:2) = target_width
                  target_dimensions(3) = target_thickness
               else
                  target_dimensions(1:3) = target_width
               end if
            end if
            call calculate_target_volume(target_dimensions, targetvol)
            if (number_spheres_specified) then
               sphere_cluster%number_spheres = input_number_spheres
               sphere_volume_fraction = dble(input_number_spheres) * four_pi_over_three / targetvol
            else
               sphere_cluster%number_spheres = ceiling(targetvol * sphere_volume_fraction) / (four_pi_over_three)
            end if
         end if
         if (target_shape .eq. 2 .and. random_configuration_host) then
            sphere_cluster%number_spheres = sphere_cluster%number_spheres + 1
         end if
      else
         sphere_cluster%number_spheres = input_number_spheres
      end if
      if (medium_ref_index_specified) then
         if (medium_reim_ref_index_specified) then
            layer_ref_index(0) = cmplx(medium_re_ref_index, medium_im_ref_index, kind=real64)
         else
            layer_ref_index(0) = medium_ref_index
         end if
      end if
      if (allocated(sphere_cluster%sphere_radius)) then
         deallocate (sphere_cluster%sphere_radius, &
                     sphere_cluster%sphere_position, &
                     sphere_cluster%sphere_ref_index, &
                     sphere_cluster%host_sphere, &
                     sphere_cluster%number_field_expansions, &
                     sphere_excitation_switch, &
                     sphere_index)
      end if
      allocate (sphere_cluster%sphere_radius(sphere_cluster%number_spheres), &
                sphere_cluster%sphere_position(3, sphere_cluster%number_spheres), &
                sphere_cluster%sphere_ref_index(2, 0:sphere_cluster%number_spheres), &
                sphere_cluster%host_sphere(sphere_cluster%number_spheres), &
                sphere_cluster%number_field_expansions(sphere_cluster%number_spheres), &
                sphere_excitation_switch(sphere_cluster%number_spheres), &
                sphere_index(sphere_cluster%number_spheres))
      if (random_configuration) then
         if (print_timings .and. mstm_global_rank .eq. 0) then
            write (sphere_cluster%run_print_unit, '('' generating random configuration:'')', advance='no')
            timet = parallel_wall_time()
         end if
         call generate_random_configuration(mpi_comm=mpicomm, skip_diffusion=dryrun)
         if (runtime_failed()) return
!            call generate_random_configuration(mpi_comm=mpicomm)
         if (print_timings .and. mstm_global_rank .eq. 0) then
            write (sphere_cluster%run_print_unit, '('' completed, time:'',es12.5,'' s'')') parallel_wall_time() - timet
         end if
         if (rank .eq. 0) then
!               if(print_random_configuration.and.(.not.configuration_average)) then
            if (print_random_configuration .and. mstm_global_rank .eq. 0) then
               call open_output_file(trim(random_configuration_output_file), file_unit)
               if (runtime_failed()) return
               do i = 1, sphere_cluster%number_spheres
                  write (file_unit, '(4es13.5)') sphere_cluster%sphere_position(:, i) / length_scale_factor, &
                     sphere_cluster%sphere_radius(i) / length_scale_factor
               end do
               close (file_unit)
            end if
         end if
      else
         call read_sphere_data_input_file(mpi_comm=mpicomm)
         if (runtime_failed()) return
      end if

      position_shift = (/x_shift, y_shift, z_shift/) * length_scale_factor
      if (shifted_sphere .gt. sphere_cluster%number_spheres) shifted_sphere = 0
      if (any(position_shift .ne. 0.d0)) then
         if (shifted_sphere .eq. 0) then
            do i = 1, sphere_cluster%number_spheres
               sphere_cluster%sphere_position(:, i) = sphere_cluster%sphere_position(:, i) + position_shift(:)
            end do
         else
           sphere_cluster%sphere_position(:, shifted_sphere) = sphere_cluster%sphere_position(:, shifted_sphere) + position_shift(:)
         end if
      end if

      sphere_cluster%sphere_ref_index(:, 0) = layer_ref_index(0)
      if (periodic_lattice) then
         if (random_configuration .and. target_shape .eq. 0) then
            cell_width(1:2) = target_dimensions(1:2) * 2.d0 * length_scale_factor
         else
            if (square_cell) then
               cell_width = input_cell_width_x * length_scale_factor
            else
               cell_width = input_cell_width * length_scale_factor
            end if
         end if
      end if

      plane_surface_present = number_plane_boundaries .gt. 0
      layer_thickness = input_layer_thickness * length_scale_factor
      call initialize_plane_boundaries()

      if (move_to_front .and. plane_surface_present) then
         zext = maxval(sphere_cluster%sphere_position(3, :) + sphere_cluster%sphere_radius(:))
         if (zext .gt. 0.d0) sphere_cluster%sphere_position(3, :) = sphere_cluster%sphere_position(3, :) - zext
      end if
      if (move_to_back .and. plane_surface_present) then
         zext = minval(sphere_cluster%sphere_position(3, :) - sphere_cluster%sphere_radius(:))
         if (zext .lt. plane_boundary_position(number_plane_boundaries)) &
            sphere_cluster%sphere_position(3, :) = sphere_cluster%sphere_position(3, :) - zext + plane_boundary_position(number_plane_boundaries)
      end if

      stopit = .false.
      if (random_orientation) then
         if (number_plane_boundaries .gt. 0) then
            if (rank .eq. 0) write (sphere_cluster%run_print_unit, '('' random orientation requires number_plane_boundaries=0'')')
            stopit = .true.
         end if
         if (periodic_lattice) then
            if (rank .eq. 0) write (sphere_cluster%run_print_unit, '('' random orientation and periodic lattice incompatible'')')
            stopit = .true.
         end if
      end if

     sphere_cluster%fft_translation_option = (input_fft_translation_option .and. sphere_cluster%number_spheres .ge. min_fft_nsphere)
      if (sphere_cluster%fft_translation_option) then
         if (number_plane_boundaries .gt. 0) then
            if (rank .eq. 0) write (sphere_cluster%run_print_unit, '('' fft option requires number_plane_boundaries=0'')')
            stopit = .true.
         end if
         if (periodic_lattice) then
            if (rank .eq. 0) write (sphere_cluster%run_print_unit, '('' fft option and periodic lattice incompatible'')')
            stopit = .true.
         end if
      end if

!         if(sphere_cluster%fft_translation_option) single_origin_expansion=.true.

      if (stopit) return
      if (light_up) then
         write (*, '('' s2 '',i3)') mstm_global_rank
         flush (6)
      end if
      sphere_cluster%cluster_origin = 0.d0
      call sphere_cluster%find_hosts()

      if (configuration_average .and. (target_shape .eq. 2) &
          .and. random_configuration_host .and. auto_target_radius) then
         sphere_cluster%host_sphere(1:sphere_cluster%number_spheres - 1) = sphere_cluster%number_spheres
         sphere_cluster%number_field_expansions(1:sphere_cluster%number_spheres - 1) = 1
         sphere_cluster%number_host_spheres = sphere_cluster%number_spheres
         sphere_cluster%host_sphere(sphere_cluster%number_spheres) = 0
         sphere_cluster%number_field_expansions(sphere_cluster%number_spheres) = 2
      end if

      if (light_up) then
         write (*, '('' s3 '',i3)') mstm_global_rank
         flush (6)
      end if
      call sphere_cluster%initialize_layers()
      call calculate_mie_coefficients(mie_epsilon)
      call initialize_numerical_tables(sphere_cluster%max_mie_order)
      if (light_up) then
         write (*, '('' s4 '',i3)') mstm_global_rank
         flush (6)
      end if

      singleorigin = number_plane_boundaries .eq. 0 .and. single_origin_expansion
      iframe = singleorigin .and. incident_frame

      sphere_cluster%cluster_origin = 0.d0
      if (singleorigin .or. random_orientation .or. .true.) then
         if (allocated(sphere_cluster%translation_order)) deallocate (sphere_cluster%translation_order)
         allocate (sphere_cluster%translation_order(sphere_cluster%number_spheres))
    sphere_cluster%translation_order(1:sphere_cluster%number_spheres) = sphere_cluster%sphere_order(1:sphere_cluster%number_spheres)
         sphere_cluster%cluster_origin = 0.d0
         if ((.not. configuration_average) .and. (.not. incidence_average)) then
            if (sphere_cluster%gaussian_beam_constant .eq. 0.d0) then
               n = 0
               do i = 1, sphere_cluster%number_spheres
                  if (sphere_cluster%host_sphere(i) .eq. 0) then
                     n = n + 1
                     sphere_cluster%cluster_origin = sphere_cluster%cluster_origin + sphere_cluster%sphere_position(:, i)
                  end if
               end do
               sphere_cluster%cluster_origin = sphere_cluster%cluster_origin / dble(n)
            else
               sphere_cluster%cluster_origin = sphere_cluster%gaussian_beam_focal_point * length_scale_factor
            end if
         end if
         if (sett) then
            maxt = sphere_cluster%max_t_matrix_order
         else
            maxt = sphere_cluster%t_matrix_order
         end if
         sphere_cluster%t_matrix_order = min(sphere_cluster%max_mie_order, maxt)
         do i = 1, sphere_cluster%number_spheres
            if (sphere_cluster%host_sphere(i) .eq. 0) then
               r0 = sphere_cluster%cluster_origin
               call exterior_refractive_index(i, rimedium)
               rtran = sqrt(sum((sphere_cluster%sphere_position(:, i) - r0(:))**2))
!                  if(rtran.gt.scattered_field_sample_length) cycle
               call estimate_translation_order(rtran, rimedium(1), sphere_cluster%sphere_order(i), &
                                               translation_epsilon, sphere_cluster%translation_order(i))
               sphere_cluster%translation_order(i) = min(sphere_cluster%translation_order(i), maxt)
               sphere_cluster%t_matrix_order = max(sphere_cluster%t_matrix_order, sphere_cluster%translation_order(i))
            end if
         end do
         if (.not. sett) sphere_cluster%t_matrix_order = maxt
      end if

      if (light_up) then
         write (*, '('' s5 '',i3)') mstm_global_rank
         flush (6)
      end if
      sphere_cluster%one_side_only = .false.
      sphere_cluster%vol_radius = 0.
      do i = 1, sphere_cluster%number_spheres
         if (sphere_cluster%host_sphere(i) .eq. 0) then
            sphere_cluster%vol_radius = sphere_cluster%vol_radius + sphere_cluster%sphere_radius(i)**3
         end if
      end do
      sphere_cluster%vol_radius = sphere_cluster%vol_radius**.333333

      if (periodic_lattice) then
         sphere_cluster%cross_section_radius = sqrt(product(cell_width) / pi)
      else
         if (sphere_cluster%gaussian_beam_constant .ne. 0.d0) then
            sphere_cluster%cross_section_radius = 1.d0 / sphere_cluster%gaussian_beam_constant / sqrt(2.d0)
         elseif (reflection_model) then
            if (random_configuration) then
               if (target_shape .eq. 0) then
                  sphere_cluster%cross_section_radius = length_scale_factor * 2.d0 * sqrt(product(target_dimensions(1:2)) / pi)
               elseif (target_shape .ge. 1) then
                  sphere_cluster%cross_section_radius = length_scale_factor * target_dimensions(1)
               end if
            else
               sphere_cluster%cross_section_radius = sqrt(product(sphere_cluster%sphere_max_position(1:2) - sphere_cluster%sphere_min_position(1:2)) / pi)
            end if
            sphere_cluster%cross_section_radius = min(sphere_cluster%cross_section_radius, length_scale_factor * excitation_radius)
         else
            sphere_cluster%cross_section_radius = sphere_cluster%vol_radius
         end if
      end if

      if (auto_absorption_sample_radius .and. random_configuration) then
         absorption_sample_radius = absorption_sample_radius_fraction * target_dimensions(1)
      end if

      if (sphere_cluster%fft_translation_option) then
         lochost = 0
         if (random_configuration) then
            tmin = -target_dimensions * length_scale_factor
            tmax = target_dimensions * length_scale_factor
            if (dryrun) call fft_plan%clear(clear_h=.true.)
            if (target_shape .eq. 2 .and. random_configuration_host) then
               lochost = sphere_cluster%number_spheres
            end if
         else
            tmin = sphere_cluster%sphere_min_position
            tmax = sphere_cluster%sphere_max_position
            if (.not. averagerun) call fft_plan%clear(clear_h=.true.)
         end if
         call fft_plan%configure(input_cell_volume_fraction, target_min=tmin, target_max=tmax, &
                                 d_specified=d_cell_specified, local_host=lochost, &
                                 requested_cell_size=input_d_cell, requested_node_order=input_node_order, &
                                 requested_neighbor_model=input_neighbor_node_model)
      end if

      sphere_excitation_switch = .true.
      do i = 1, sphere_cluster%number_spheres
         if (sphere_cluster%host_sphere(i) .ne. 0) cycle
         if (random_configuration) then
            if (target_shape .le. 1) then
               rtran = sqrt(sum((sphere_cluster%sphere_position(1:2, i) - sphere_cluster%cluster_origin(1:2))**2))
            else
               rtran = sqrt(sum((sphere_cluster%sphere_position(1:3, i) - sphere_cluster%cluster_origin(1:3))**2))
            end if
         else
            rtran = sqrt(sum((sphere_cluster%sphere_position(1:3, i) - sphere_cluster%cluster_origin(1:3))**2))
         end if
         if (excitation_radius .gt. 0.d0) then
            sphere_excitation_switch(i) = rtran .le. excitation_radius * length_scale_factor
         else
            sphere_excitation_switch(i) = i .le. -int(excitation_radius)
         end if
      end do
      if (excitation_radius .eq. 0.d0) then
         sphere_excitation_switch = .false.
         if (rank .eq. 0) then
            call random_seed()
            do
               call random_number(rannum)
               itemp(1) = 1 + floor(sphere_cluster%number_spheres * rannum)
               if (sphere_cluster%host_sphere(itemp(1)) .eq. 0) exit
            end do
         end if
         call parallel_broadcast(mpi_rank=0, &
                                 send_buffer=itemp(1), mpi_number=1, mpi_comm=mpicomm)
         sphere_excitation_switch(itemp(1)) = .true.
      end if
      if (light_up) then
         write (*, '('' s6 '',i3)') mstm_global_rank
         flush (6)
      end if
      if (random_orientation) then
         qeff_dim = 1
         if (calculate_scattering_matrix) then
            scat_mat_udim = floor(180.00001d0 / scattering_map_increment)
            scat_mat_mdim = 16
            scat_mat_ldim = 0
            scat_mat_amin = 0.d0
            scat_mat_amax = 180.d0
         end if
      else
         if (incident_beta_specified) then
            incident_beta = incident_beta_deg * degrees_to_radians
            if (incident_beta_deg .le. 90.d0) then
               incident_direction = 1
               incident_sin_beta = sind(incident_beta_deg) / dble(layer_ref_index(0))
            else
               incident_direction = 2
               incident_sin_beta = sind(incident_beta_deg) &
                                   / dble(layer_ref_index(number_plane_boundaries))
            end if
         else
            incident_beta = 0.d0
         end if
         if (incidence_average) then
            qeff_dim = 1
         else
            qeff_dim = 3
         end if
         alpha = incident_alpha_deg * degrees_to_radians
         call initialize_incident_field(alpha, incident_sin_beta, incident_direction)
         if (calculate_scattering_matrix) then
            if (allocated(scat_mat)) deallocate (scat_mat)
            if (periodic_lattice) then
               call periodic_lattice_scattering(amnp_s, pl_sca, dry_run=.true., num_dirs=number_rl_dirs)
               max_number_rl_dirs = maxval(number_rl_dirs)
               if (allocated(rl_vec)) deallocate (rl_vec)
               allocate (rl_vec(2, max_number_rl_dirs))
               scat_mat_udim = max_number_rl_dirs
               scat_mat_ldim = 1
               scat_mat_mdim = 32
            else
               if (scattering_map_model .eq. 0) then
                  if (number_plane_boundaries .eq. 0) then
                     scat_mat_udim = floor(180.00001d0 / scattering_map_increment)
                     if (azimuthal_average) then
                        scat_mat_ldim = 0
                     else
                        scat_mat_ldim = -scat_mat_udim
                     end if
                     scat_mat_mdim = 16
                     scat_mat_amax = 180.d0
                  else
                     scat_mat_udim = floor(90.00001d0 / scattering_map_increment)
!                        scat_mat_ldim=-scat_mat_udim
! 10-22 azimuthal_average applies to multiple plane boundaries
                     if (azimuthal_average) then
                        scat_mat_ldim = 0
                     else
                        scat_mat_ldim = -scat_mat_udim
                     end if
                     scat_mat_mdim = 32
                     scat_mat_amax = 90.d0
                  end if
                  scat_mat_amin = scat_mat_amax * (scat_mat_ldim / scat_mat_udim)
               else
                  i = 0
                  do sy = -scattering_map_dimension, scattering_map_dimension
                     do sx = -scattering_map_dimension, scattering_map_dimension
                        if (sx * sx + sy * sy .gt. scattering_map_dimension**2) cycle
                        i = i + 1
                     end do
                  end do
                  scat_mat_udim = i
                  scat_mat_ldim = 1
                  scat_mat_mdim = 32
               end if
            end if
         end if
         if (periodic_lattice .or. reflection_model .and. (target_shape .le. 1)) then
            sphere_cluster%cross_section_radius = sphere_cluster%cross_section_radius * sqrt(cos(incident_beta))
         end if
      end if
      if (light_up) then
         write (*, '('' s7 '',i3)') mstm_global_rank
         flush (6)
      end if
      if (allocated(boundary_sca)) deallocate (boundary_sca, boundary_ext)
!         allocate(boundary_sca(2,0:number_plane_boundaries+1),boundary_ext(2,0:number_plane_boundaries+1))
      allocate (boundary_sca(2, 0:1), boundary_ext(2, 0:1))

      if (allocated(q_eff)) deallocate (q_eff, q_eff_tot, q_vabs)
      allocate (q_eff(3, qeff_dim, sphere_cluster%number_spheres), q_eff_tot(3, qeff_dim), q_vabs(qeff_dim, sphere_cluster%number_spheres))
      if (calculate_scattering_matrix) then
         if (allocated(scat_mat)) deallocate (scat_mat)
         allocate (scat_mat(scat_mat_mdim, scat_mat_ldim:scat_mat_udim))
      end if
      if (light_up) then
         write (*, '('' s8 '',i3)') mstm_global_rank
         flush (6)
      end if
!call parallel_barrier()
      if (rank .eq. 0) time1 = parallel_wall_time()
      if (rank .eq. 0 .and. printout) then
         if (check_positions) call check_sphere_positions()
         call print_run_variables(sphere_cluster%run_print_unit)
         call open_output_file(output_file, file_unit, append=.true.)
         if (runtime_failed()) return
         call print_run_variables(file_unit)
         close (file_unit)
      end if

      if (dryrun) return

      if (random_orientation) then
         niter = max_iterations
         timatrixfile = 'titemp.dat'
         if (allocated(mean_t)) deallocate (mean_t)
         allocate (mean_t(2, sphere_cluster%t_matrix_order))
         mean_t = 0.d0
         call solve_t_matrix(solution_method=solution_method(1:1), &
                             solution_eps=solution_epsilon, &
                             convergence_eps=t_matrix_convergence_epsilon, &
                             max_iterations=niter, &
                             t_matrix_file=t_matrix_output_file, &
                             procs_per_soln=t_matrix_procs_per_solution, &
                             sphere_qeff=q_eff, &
                             solution_status=istat, &
                             mpi_comm=mpicomm, &
                             sphere_excitation_list=sphere_excitation_switch)
         if (runtime_failed()) return
         if (istat /= solver_converged) then
            call set_runtime_error('T-matrix solver failed: '//trim(solver_status_message(istat)), istat)
            return
         end if
         if (sphere_cluster%fft_translation_option) call fft_plan%clear(clear_h=.true.)
         if (calculate_scattering_matrix) then
            if (allocated(scat_mat_exp_coef)) deallocate (scat_mat_exp_coef)
            if (allocated(coh_scat_mat_exp_coef)) deallocate (coh_scat_mat_exp_coef)
            allocate (scat_mat_exp_coef(4, 4, 0:2 * sphere_cluster%t_matrix_order), coh_scat_mat_exp_coef(4, 4, 0:2 * sphere_cluster%t_matrix_order))
            nodrw = 2 * sphere_cluster%t_matrix_order
            call random_orientation_scattering_matrix(t_matrix_output_file, scat_mat_exp_coef, &
                                                      coh_scat_mat_exp_coef, &
                                                      beam_width=sphere_cluster%gaussian_beam_constant, &
                                                      number_processors=t_matrix_procs_per_solution, &
                                                      mean_t_matrix=mean_t, mpi_comm=mpicomm)
            coherent_scattering_ratio = coh_scat_mat_exp_coef(1, 1, 0)
            do i = scat_mat_ldim, scat_mat_udim
               costheta = cos(dble(i - scat_mat_ldim) * pi / dble(scat_mat_udim - scat_mat_ldim))
               call evaluate_random_orientation_scattering_matrix(costheta, scat_mat_exp_coef, nodrw, scat_mat(:, i))
            end do
         end if
         call total_efficiency_factors(sphere_cluster%number_spheres, qeff_dim, sphere_cluster%cross_section_radius, &
                                       q_eff, q_vabs, q_eff_tot)
      else
         if (light_up) then
            write (*, '('' s8.1 '',i3)') mstm_global_rank
            flush (6)
         end if
!call parallel_barrier()

         if (allocated(amnp_s)) deallocate (amnp_s)
         allocate (amnp_s(sphere_cluster%number_eqns, 2))
         amnp_s = 0.d0
         niter = max_iterations
         error_codes = 0
         pl_error_codes = 0
!            sphere_cluster%recalculate_surface_matrix=.true.
         if (light_up) then
            write (*, '('' s8.2 '',i3)') mstm_global_rank
            flush (6)
         end if
!call parallel_barrier()
         if (print_timings .and. mstm_global_rank .eq. 0) then
            write (sphere_cluster%run_print_unit, '('' generating solution:'')', advance='no')
            timet = parallel_wall_time()
         end if
         call solve_fixed_orientation(alpha, incident_sin_beta, incident_direction, solution_epsilon, niter, &
                                      amnp_s, q_eff, &
                                      qeff_dim, solution_error, solution_iterations, 1, istat, &
                                      mpi_comm=mpicomm, &
                                      excited_spheres=sphere_excitation_switch, &
                                      solution_method=solution_method(1:1), &
                                      initialize_solver=.true., &
                                      reciprocal_condition=solution_reciprocal_condition)
         if (runtime_failed()) return
         if (istat /= solver_converged) then
            call set_runtime_error('Fixed-orientation solver failed: '//trim(solver_status_message(istat)), istat)
            return
         end if
         if (print_timings .and. mstm_global_rank .eq. 0) then
            write (sphere_cluster%run_print_unit, '('' completed, time:'',es12.5,'' s'')') parallel_wall_time() - timet
         end if
         if (sphere_cluster%fft_translation_option) then
            call fft_plan%clear()
         end if
         if (mstm_global_rank .eq. 0. .and. ((.not. configuration_average) .and. (.not. incidence_average))) then
            write (sphere_cluster%run_print_unit, '('' solution completed: number iterations='',i5)') solution_iterations
            flush (sphere_cluster%run_print_unit)
         end if

         if (light_up) then
            write (*, '('' s8.3 '',i3)') mstm_global_rank
            flush (6)
         end if
!call parallel_barrier()
         if (print_timings .and. mstm_global_rank .eq. 0) then
            write (sphere_cluster%run_print_unit, '('' post processing solution:'')', advance='no')
            timet = parallel_wall_time()
         end if

         call total_efficiency_factors(sphere_cluster%number_spheres, qeff_dim, sphere_cluster%cross_section_radius, &
                                       q_eff, q_vabs, q_eff_tot)
!            q_eff_tot(3,:)=q_eff_tot(1,:)-q_eff_tot(2,:)
         csca = q_eff_tot(3, 1) * pi * sphere_cluster%cross_section_radius**2
         if (singleorigin) then
            if (allocated(amnp_0)) deallocate (amnp_0)
            allocate (amnp_0(2 * sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2), 2))
            amnp_0 = 0.d0
            if (light_up) then
               write (*, '('' s8.3.1 '',i3)') mstm_global_rank
               flush (6)
            end if
!call parallel_barrier()
            do i = 1, 2
               call merge_to_common_origin(sphere_cluster%t_matrix_order, amnp_s(:, i), amnp_0(:, i), &
                                           origin_position=sphere_cluster%cluster_origin, merge_procs=.true., &
                                           mpi_comm=mpicomm)
               if (iframe) then
                  call rotate_expansion_coefficients(alpha, incident_beta, 0.d0, sphere_cluster%t_matrix_order, &
                                                     sphere_cluster%t_matrix_order, amnp_0(:, i), 1)
               end if
            end do
            if (light_up) then
               write (*, '('' s8.3.2 '',i3)') mstm_global_rank
               flush (6)
            end if
!call parallel_barrier()
            if (singleorigin .and. azimuthal_average .and. (.not. numerical_azimuthal_average)) then
               if (allocated(scat_mat_exp_coef)) deallocate (scat_mat_exp_coef)
               allocate (scat_mat_exp_coef(16, 0:2 * sphere_cluster%t_matrix_order, 4))
               call fixed_orientation_scattering_matrix_expansion( &
                  sphere_cluster%t_matrix_order, amnp_0, scat_mat_exp_coef(:, :, 1), scat_mat_exp_coef(:, :, 2), &
                  scat_mat_exp_coef(:, :, 3), scat_mat_exp_coef(:, :, 4), mpi_comm=mpicomm)
            end if
            if (light_up) then
               write (*, '('' s8.3.4 '',i3)') mstm_global_rank
               flush (6)
            end if
!call parallel_barrier()
            if (calculate_scattering_matrix) then
               call compute_scattering_matrix(amnp_0, scat_mat, mpi_comm=mpicomm)
            end if
            if (light_up) then
               write (*, '('' s8.3.5 '',i3)') mstm_global_rank
               flush (6)
            end if
!call parallel_barrier()
            if (sphere_cluster%gaussian_beam_constant .eq. 0.d0) then
               boundary_ext = 0.d0
               ! A finite cluster without interfaces does not use plane-boundary
               ! extinction.  Avoid the spectral Green-function normalization at
               ! exact grazing incidence, where its longitudinal wave number is zero.
               if (number_plane_boundaries .gt. 0 .or. periodic_lattice .or. reflection_model) then
                  call boundary_extinction(amnp_0, alpha, incident_sin_beta, incident_direction, boundary_ext, &
                                           common_origin=singleorigin)
               end if
            end if
         else
            if (light_up) then
               write (*, '('' s8.3.5 '',i3)') mstm_global_rank
               flush (6)
            end if
!call parallel_barrier()
            if (light_up) then
               write (*, '('' s8.3.5 '',i3)') mstm_global_rank
               flush (6)
            end if
!call parallel_barrier()
            if (calculate_scattering_matrix) then
               call compute_scattering_matrix(amnp_s, scat_mat, mpi_comm=mpicomm)
            end if
            if (sphere_cluster%gaussian_beam_constant .eq. 0.d0) then
               boundary_ext = 0.d0
               if (number_plane_boundaries .gt. 0 .or. periodic_lattice .or. reflection_model) then
                  call boundary_extinction(amnp_s, alpha, incident_sin_beta, incident_direction, boundary_ext)
               end if
            end if
         end if

         if (calculate_up_down_scattering) then
            if (singleorigin) then
               call hemispherical_scattering(amnp_0, .true., numerical_hemispherical_integration, &
                                             boundary_sca, mpi_comm=mpicomm)
            else
               call hemispherical_scattering(amnp_s, .false., numerical_hemispherical_integration, &
                                             boundary_sca, mpi_comm=mpicomm)
            end if
         end if

         if (sphere_cluster%gaussian_beam_constant .ne. 0.d0) then
            boundary_ext = 0
            boundary_ext(1:2, 1) = -q_eff_tot(1, 2:3)
         end if

         if (periodic_lattice) then
            call periodic_lattice_scattering(amnp_s, pl_sca)
         elseif (reflection_model) then
!               pl_sca(:,1)=boundary_sca(:,number_plane_boundaries+1)
            pl_sca(:, 1) = boundary_sca(:, 1)
            pl_sca(:, 2) = -boundary_sca(:, 0)
         end if

         if (periodic_lattice .or. reflection_model) then
            call calculate_surface_absorptance()
         end if

         if (number_plane_boundaries .gt. 0 .and. .not. periodic_lattice) then
!               evan_sca(1:2)=q_eff_tot(1,2:3)-q_eff_tot(2,2:3)+boundary_sca(1:2,0) &
!                  - boundary_sca(:,number_plane_boundaries+1)
            evan_sca(1:2) = q_eff_tot(1, 2:3) - q_eff_tot(2, 2:3) + boundary_sca(1:2, 0) &
                            - boundary_sca(:, 1)
         end if

         if (print_timings .and. mstm_global_rank .eq. 0) then
            write (sphere_cluster%run_print_unit, '('' completed, time:'',es12.5,'' s'')') parallel_wall_time() - timet
         end if

      end if

      if (rank .eq. 0) solution_time = parallel_wall_time() - time1

      call gather_error_codes(mpicomm)
      if (light_up) then
         write (*, '('' s12 '',i3)') mstm_global_rank
         flush (6)
      end if
      if (mstm_global_rank .eq. 0 .and. printout) then
         call print_calculation_results(output_file)
      end if
      if (light_up) then
         write (*, '('' s13 '',i3)') mstm_global_rank
         flush (6)
      end if
      error_codes = 0
      pl_error_codes = 0

      if (calculate_near_field .and. (.not. random_orientation)) then
         celldim = ceiling((near_field_plane_vertices(:, 2) - near_field_plane_vertices(:, 1)) / near_field_step_size)
         celldim = max(celldim, (/1, 1, 1/))
         if (configuration_average) then
            if (allocated(e_field)) deallocate (e_field, h_field)
            if (rank .eq. 0) then
               allocate (e_field(3, 2, celldim(1), celldim(2), celldim(3)), &
                         h_field(3, 2, celldim(1), celldim(2), celldim(3)))
               e_field = 0.d0
               h_field = 0.d0
            end if
            call compute_near_field(amnp_s, alpha, incident_sin_beta, incident_direction, &
                                    near_field_plane_vertices, celldim, &
                                    incident_model=near_field_calculation_model, output_unit=0, &
                                    e_field_array=e_field, h_field_array=h_field, mpi_comm=mpicomm)
         else
            if (append_near_field_output_file) then
               call open_output_file(near_field_output_file, file_unit, append=.true.)
            else
               call open_output_file(near_field_output_file, file_unit)
            end if
            if (runtime_failed()) return
            append_near_field_output_file = .true.
            call compute_near_field(amnp_s, alpha, incident_sin_beta, incident_direction, &
                                    near_field_plane_vertices, celldim, &
                                    incident_model=near_field_calculation_model, output_unit=file_unit, output_header=.true.)
            close (file_unit)
         end if
         call gather_error_codes(mpicomm)
         if (rank .eq. 0) call print_error_codes(sphere_cluster%run_print_unit)
      end if

   end subroutine execute_simulation

   subroutine gather_error_codes(mpicomm)
      implicit none
      integer :: mpicomm, itemp(6)
      if (number_plane_boundaries .gt. 0) then
         itemp(1:4) = error_codes
         call parallel_reduce_sum(mpi_rank=0, mpi_number=4, &
                                  receive_buffer=error_codes, send_buffer=itemp(1:4), mpi_comm=mpicomm)
      else
         error_codes = 0
      end if
      if (periodic_lattice) then
         itemp(1:6) = pl_error_codes
         call parallel_reduce_sum(mpi_rank=0, mpi_number=6, &
                                  receive_buffer=pl_error_codes, send_buffer=itemp(1:6), mpi_comm=mpicomm)
      else
         pl_error_codes = 0
      end if
   end subroutine gather_error_codes

   subroutine run_configuration_average()
      implicit none
      logical :: singleorigin, iframe
      integer :: file_unit, rank, numprocs, m, n, p, mnp, griddim(3), ipos(3), ix, iy, iz, &
                 numprocsperconfig, configcolor, configgroup, configcomm, configrank, config0comm, nconfigave, nsend
      real(real64) :: time1, timet, diffac, csca(1), xspfit, rpos(3), rtemp(1)
      real(real64), allocatable :: texpcoef(:, :, :), spherical_position(:, :)
      complex(real64) :: ritemp(2), aneff, ctemp(1), rieff, e0
      complex(real64), allocatable :: pmnp0(:, :), anp0(:, :), edat(:)
      character(len=256) :: tmatchar1, tmatchar2
      data tmatchar1, tmatchar2/'tmat-', '.tmp'/
      first_run = .false.
      call parallel_rank(mpi_rank=rank)
      call parallel_size(mpi_size=numprocs)
!         if(rank.ne.0) light_up=.false.
      local_rank = rank
      global_rank = rank

      if (max_iterations .le. 1) then
         numprocsperconfig = 2
      else
         numprocsperconfig = 4
      end if
!numprocsperconfig=2

      n_configuration_groups = numprocs / numprocsperconfig
      n_configuration_groups = max(n_configuration_groups, 1)
      configcolor = floor(dble(n_configuration_groups * rank) / dble(numprocs))
      configgroup = configcolor
      call parallel_split( &
         mpi_color=configcolor, mpi_key=rank, &
         mpi_new_comm=configcomm)
      call parallel_rank( &
         mpi_rank=configrank, &
         mpi_comm=configcomm)
      configcolor = configrank
      call parallel_split( &
         mpi_color=configcolor, mpi_key=rank, &
         mpi_new_comm=config0comm)
      random_configuration = .true.
      singleorigin = number_plane_boundaries .eq. 0 .and. single_origin_expansion
!singleorigin=.true.
      iframe = singleorigin .and. incident_frame

      call execute_simulation(print_output=.false., set_t_matrix_order=.true., dry_run=.true.)
      if (runtime_failed()) return

      if (allocated(q_eff_ave)) deallocate (q_eff_ave, q_eff_tot_ave, q_vabs_ave, sphere_position_ave, boundary_sca_ave, &
                                            boundary_ext_ave, dif_boundary_sca)
      allocate (q_eff_ave(3, qeff_dim, sphere_cluster%number_spheres), q_eff_tot_ave(3, qeff_dim), q_vabs_ave(qeff_dim, sphere_cluster%number_spheres), &
                sphere_position_ave(3, sphere_cluster%number_spheres), boundary_sca_ave(2, 0:1), &
                boundary_ext_ave(2, 0:1), dif_boundary_sca(2, 0:1))
      q_eff_ave = 0.d0
      q_eff_tot_ave = 0.d0
      q_vabs_ave = 0.d0
      sphere_position_ave = 0.d0
      if (target_shape .eq. 2) allocate (spherical_position(3, sphere_cluster%number_spheres))
      pl_sca_ave = 0.d0
      boundary_sca_ave = 0.d0
      boundary_ext_ave = 0.d0
      solution_time_ave = 0.d0
      surface_absorptance_ave = 0.
      sphere_cluster%effective_medium_simulation = .false.
      if (singleorigin) then
         if (allocated(amnp_0_ave)) deallocate (amnp_0_ave, scat_mat_exp_coef_ave, coh_scat_mat_exp_coef)
         allocate (amnp_0_ave(2 * sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2), 2), &
                   scat_mat_exp_coef_ave(16, 0:2 * sphere_cluster%t_matrix_order, 4), &
                   coh_scat_mat_exp_coef(16, 0:2 * sphere_cluster%t_matrix_order, 4))
         amnp_0_ave = 0.d0
         scat_mat_exp_coef_ave = 0.d0
         tot_csca_ave = 0.d0
         dif_csca_ratio = 0.d0
         if (input_effective_medium_simulation .and. target_shape .eq. 2 .and. &
             (.not. random_configuration_host)) then
            sphere_cluster%effective_medium_simulation = .true.
            sphere_cluster%effective_ref_index = layer_ref_index(0)
            if (random_configuration_host_model .eq. 1) then
               sphere_cluster%effective_cluster_radius = target_dimensions(1) * length_scale_factor
            elseif (random_configuration_host_model .eq. 2) then
               sphere_cluster%effective_cluster_radius = sphere_cluster%vol_radius / (sphere_volume_fraction)**0.33333
            end if
         end if
      end if
      if (calculate_scattering_matrix) then
         if (allocated(scat_mat_ave)) deallocate (scat_mat_ave, dif_scat_mat)
         allocate (scat_mat_ave(scat_mat_mdim, scat_mat_ldim:scat_mat_udim), &
                   dif_scat_mat(scat_mat_mdim, scat_mat_ldim:scat_mat_udim))
         scat_mat_ave = 0.d0
      end if
      if (calculate_near_field .and. configrank .eq. 0) then
         griddim = ceiling((near_field_plane_vertices(:, 2) - near_field_plane_vertices(:, 1)) / near_field_step_size)
         griddim = max(griddim, (/1, 1, 1/))
         if (allocated(e_field_ave)) deallocate (e_field_ave, h_field_ave, s_field, s_field_ave)
         allocate (e_field_ave(3, 2, griddim(1), griddim(2), griddim(3)), &
                   h_field_ave(3, 2, griddim(1), griddim(2), griddim(3)), &
                   s_field(3, 2, griddim(1), griddim(2), griddim(3)), &
                   s_field_ave(3, 2, griddim(1), griddim(2), griddim(3)))
         e_field_ave = 0.d0
         h_field_ave = 0.d0
         s_field_ave = 0.d0
      end if

      nconfigave = 0
      do random_configuration_number = 1, ceiling(dble(number_configurations) / dble(n_configuration_groups))

         if (rank .eq. 0) then
            if (random_configuration_number .eq. 1) then
               call print_run_variables(sphere_cluster%run_print_unit)
               call open_output_file(output_file, file_unit, append=.true.)
               if (runtime_failed()) return
               call print_run_variables(file_unit)
               close (file_unit)
            end if
            write (sphere_cluster%run_print_unit, '('' configuration averaging, samples:'',i5,''-'',i5)') &
               (random_configuration_number - 1) * n_configuration_groups + 1, &
               random_configuration_number * n_configuration_groups
         end if

         if (rank .eq. 0) time1 = parallel_wall_time()

         call execute_simulation(print_output=.false., set_t_matrix_order=.false., mpi_comm=configcomm)
         if (runtime_failed()) return

         if (singleorigin .and. configrank .eq. 0) then
            call common_origin_scattering_cross_section(sphere_cluster%t_matrix_order, amnp_0, csca)
         end if

         if (configrank .eq. 0) then
            if (rank .eq. 0) solution_time = parallel_wall_time() - time1
            q_eff_ave = q_eff_ave + q_eff
            q_eff_tot_ave = q_eff_tot_ave + q_eff_tot
            q_vabs_ave = q_vabs_ave + q_vabs
            if (target_shape .eq. 2) then
              call cartesian_vectors_to_spherical(sphere_cluster%number_spheres, sphere_cluster%sphere_position, spherical_position)
               sphere_position_ave = sphere_position_ave + spherical_position
            else
               sphere_position_ave = sphere_position_ave + sphere_cluster%sphere_position
            end if
            pl_sca_ave = pl_sca_ave + pl_sca
            surface_absorptance_ave = surface_absorptance_ave + surface_absorptance
            if (calculate_up_down_scattering) boundary_sca_ave = boundary_sca_ave + boundary_sca
            boundary_ext_ave = boundary_ext_ave + boundary_ext
            if (singleorigin .and. azimuthal_average .and. (.not. numerical_azimuthal_average)) &
               scat_mat_exp_coef_ave = scat_mat_exp_coef_ave + scat_mat_exp_coef
            if (calculate_scattering_matrix) then
               scat_mat_ave = scat_mat_ave + scat_mat
            end if
            if (calculate_near_field) then
               e_field_ave = e_field_ave + e_field
               h_field_ave = h_field_ave + h_field
               s_field(1, :, :, :, :) = 0.5 * (e_field(2, :, :, :, :) * conjg(h_field(3, :, :, :, :)) &
                                               - e_field(3, :, :, :, :) * conjg(h_field(2, :, :, :, :)))
               s_field(2, :, :, :, :) = 0.5 * (-e_field(1, :, :, :, :) * conjg(h_field(3, :, :, :, :)) &
                                               + e_field(3, :, :, :, :) * conjg(h_field(1, :, :, :, :)))
               s_field(3, :, :, :, :) = 0.5 * (e_field(1, :, :, :, :) * conjg(h_field(2, :, :, :, :)) &
                                               - e_field(2, :, :, :, :) * conjg(h_field(1, :, :, :, :)))
               s_field_ave = s_field_ave + s_field
            end if
            if (rank .eq. 0) solution_time_ave = solution_time_ave + solution_time
         end if

         if (singleorigin) then
!               if(sphere_1_fixed) call subtract_fixed_sphere_from_common_origin()
!               if(sphere_1_fixed) then
!                  erase_sphere_1=.true.
!                  use_previous_configuration=.true.
!                  sphere_1_fixed=.false.
!                  call execute_simulation(print_output=.false.,set_t_matrix_order=.false.,mpi_comm=configcomm)
!                  erase_sphere_1=.false.
!                  use_previous_configuration=.false.
!                  sphere_1_fixed=.true.
!               endif
            amnp_0_ave = amnp_0_ave + amnp_0
            tot_csca_ave = tot_csca_ave + csca
         end if

         nsend = 3 * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=sphere_position_ave, &
            receive_buffer=sphere_cluster%sphere_position, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 3 * qeff_dim * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=q_eff_ave, &
            receive_buffer=q_eff, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 3 * qeff_dim
         call parallel_reduce_sum( &
            send_buffer=q_eff_tot_ave, &
            receive_buffer=q_eff_tot, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = qeff_dim * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=q_vabs_ave, &
            receive_buffer=q_vabs, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 4
         call parallel_reduce_sum( &
            send_buffer=pl_sca_ave, &
            receive_buffer=pl_sca, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 2
         call parallel_reduce_sum( &
            send_buffer=surface_absorptance_ave, &
            receive_buffer=surface_absorptance, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         if (calculate_up_down_scattering) then
            nsend = 4
            call parallel_reduce_sum( &
               send_buffer=boundary_sca_ave, &
               receive_buffer=boundary_sca, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
         call parallel_reduce_sum( &
            send_buffer=boundary_ext_ave, &
            receive_buffer=boundary_ext, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         if (singleorigin) then
            nsend = 4 * sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2)
            call parallel_reduce_sum( &
               send_buffer=amnp_0_ave, &
               receive_buffer=amnp_0, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend = 1
            call parallel_reduce_sum( &
               send_buffer=tot_csca_ave, &
               receive_buffer=csca, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
         if (calculate_scattering_matrix) then
            nsend = scat_mat_mdim * (scat_mat_udim - scat_mat_ldim + 1)
            call parallel_reduce_sum( &
               send_buffer=scat_mat_ave, &
               receive_buffer=scat_mat, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
         if (singleorigin .and. azimuthal_average .and. (.not. numerical_azimuthal_average)) then
            nsend = 16 * 4 * (2 * sphere_cluster%t_matrix_order + 1)
            call parallel_reduce_sum( &
               send_buffer=scat_mat_exp_coef_ave, &
               receive_buffer=scat_mat_exp_coef, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
         if (calculate_near_field) then
            nsend = 6 * product(griddim)
            call parallel_reduce_sum( &
               send_buffer=e_field_ave, &
               receive_buffer=e_field, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            call parallel_reduce_sum( &
               send_buffer=h_field_ave, &
               receive_buffer=h_field, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            call parallel_reduce_sum( &
               send_buffer=s_field_ave, &
               receive_buffer=s_field, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if

         nconfigave = nconfigave + n_configuration_groups
         diffac = (dble(sphere_cluster%number_spheres) - 1.d0) / dble(sphere_cluster%number_spheres)
         diffac = 1.d0
!            diffac=(dble(sphere_cluster%number_spheres))/dble(sphere_cluster%number_spheres-1)
!            diffac=(1.d0-(1.d0/dble(sphere_cluster%number_spheres))**.5d0)

         if (singleorigin) then
            if (rank .eq. 0 .and. print_timings) then
               timet = parallel_wall_time()
               write (sphere_cluster%run_print_unit, '('' calculating diffuse field:'')', advance='no')
            end if
            amnp_0 = amnp_0 / dble(nconfigave)
!
!  zero out azimuth order .ne. pm 1 for sphere targets 2/23
!
!               if(target_shape.eq.2) then
!                  do n=1,sphere_cluster%t_matrix_order
!                     do m=-n,n
!                        do p=1,2
!                           if(abs(m).ne.1) then
!                              amnp_0(polarized_mode_index(m,n,p,sphere_cluster%t_matrix_order,2),:)=0.d0
!                           endif
!                        enddo
!                     enddo
!                  enddo
!               endif
            csca = csca / dble(nconfigave)
            tot_csca = csca(1)
            call common_origin_scattering_cross_section(sphere_cluster%t_matrix_order, amnp_0, dif_csca_ratio)
            dif_csca = tot_csca - dif_csca_ratio(1)
            dif_csca_ratio = 1.d0 - dif_csca_ratio / csca
            if (azimuthal_average .and. (.not. numerical_azimuthal_average)) then
               allocate (texpcoef(16, 0:2 * sphere_cluster%t_matrix_order, 4))
               texpcoef = scat_mat_exp_coef / dble(nconfigave)
!                  scat_mat_exp_coef=scat_mat_exp_coef/dble(nconfigave)
               call fixed_orientation_scattering_matrix_expansion( &
                  sphere_cluster%t_matrix_order, amnp_0, coh_scat_mat_exp_coef(:, :, 1), coh_scat_mat_exp_coef(:, :, 2), &
                  coh_scat_mat_exp_coef(:, :, 3), coh_scat_mat_exp_coef(:, :, 4), mpi_comm=configcomm)
               scat_mat_exp_coef = coh_scat_mat_exp_coef
            end if
            if (calculate_scattering_matrix) then
               call compute_scattering_matrix(amnp_0, dif_scat_mat, mpi_comm=configcomm)
            end if
            if (azimuthal_average .and. (.not. numerical_azimuthal_average)) then
!                  scat_mat_exp_coef=texpcoef-scat_mat_exp_coef*diffac
               scat_mat_exp_coef = texpcoef
               deallocate (texpcoef)
            end if
            call hemispherical_scattering(amnp_0, .true., numerical_hemispherical_integration, &
                                          dif_boundary_sca, mpi_comm=configcomm)
!               call common_origin_hemispherical_scattering(amnp_0,dif_boundary_sca)
      if (rank .eq. 0 .and. print_timings) write (sphere_cluster%run_print_unit, '('' completed, '',es12.4,'' sec'')') parallel_wall_time() - timet
            if (rank .eq. 0 .and. target_shape .eq. 2 .and. (.not. random_configuration_host)) then
allocate (pmnp0(2 * sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2), 2), anp0(2, sphere_cluster%t_matrix_order))
               call generate_plane_wave_coefficients(0.d0, (1.d0, 0.d0), sphere_cluster%t_matrix_order, pmnp0, lr_tran=.false.)
               call open_output_file('anpeff.dat', file_unit)
               if (runtime_failed()) return
               write (file_unit, '(i5)') sphere_cluster%t_matrix_order
               do n = 1, sphere_cluster%t_matrix_order
                  do p = 1, 2
                     aneff = 0.d0
                     do m = -1, 1, 2
                        mnp = polarized_mode_index(m, n, p, sphere_cluster%t_matrix_order, 2)
                        aneff = aneff + 0.5d0 * sum(amnp_0(mnp, :) / pmnp0(mnp, :))
                     end do
                     aneff = aneff / 2.d0
                     write (file_unit, '(2i4,2es13.5)') n, p, aneff
                     anp0(p, n) = aneff
                  end do
               end do
               close (file_unit)
               deallocate (pmnp0)
!                  sphere_cluster%effective_ref_index=(1.1d0,0.01d0)
!                  fit_radius=target_dimensions(1)*length_scale_factor
               call fit_effective_refractive_index(anp0, sphere_cluster%effective_ref_index, fit_radius, fit_stat)
               deallocate (anp0)
            end if
         end if
         if (allocated(amnp_0)) deallocate (amnp_0)

         if (rank .eq. 0) then
            sphere_cluster%sphere_position = sphere_cluster%sphere_position / dble(nconfigave)
            q_eff = q_eff / dble(nconfigave)
            q_vabs = q_vabs / dble(nconfigave)
            q_eff_tot = q_eff_tot / dble(nconfigave)
            pl_sca = pl_sca / dble(nconfigave)
            surface_absorptance = surface_absorptance / dble(nconfigave)
            if (calculate_up_down_scattering) boundary_sca = boundary_sca / dble(nconfigave)
            boundary_ext = boundary_ext / dble(nconfigave)
!               if(singleorigin) dif_boundary_sca=boundary_sca-dif_boundary_sca
            if (singleorigin) dif_boundary_sca = boundary_sca &
                                                 - dif_boundary_sca * diffac
            if (calculate_scattering_matrix) then
               scat_mat = scat_mat / dble(nconfigave)
! experiment

!                  if(singleorigin) dif_scat_mat=scat_mat-dif_scat_mat*(dble(sphere_cluster%number_spheres-1)/dble(sphere_cluster%number_spheres))
!                  if(singleorigin) dif_scat_mat=scat_mat-dif_scat_mat*diffac

!                  dif_scat_mat=scat_mat-dif_scat_mat
            end if
            solution_time = solution_time_ave / dble(random_configuration_number)
            call print_calculation_results(output_file)
            if (calculate_near_field) then
               e_field = e_field / dble(nconfigave)
               h_field = h_field / dble(nconfigave)
               s_field = s_field / dble(nconfigave)
               call open_output_file(near_field_output_file, file_unit)
               if (runtime_failed()) return
               call write_near_field_output_header(griddim, file_unit, print_intersecting_spheres=.false.)
               allocate (edat(griddim(3)))
               edat = 0.d0
               do iz = 1, griddim(3)
                  do iy = 1, griddim(2)
                     do ix = 1, griddim(1)
                        edat(iz) = edat(iz) + e_field(1, 1, ix, iy, iz) + e_field(2, 2, ix, iy, iz)
                        ipos(:) = (/ix, iy, iz/)
                        rpos(:) = (dble(ipos(:)) - (/0.5d0, 0.5d0, 0.5d0/)) * grid_spacing(:) + grid_region(:, 1)
                        write (file_unit, '(33es12.4)') rpos(:), &
                           e_field(:, 1, ix, iy, iz), h_field(:, 1, ix, iy, iz), &
                           e_field(:, 2, ix, iy, iz), h_field(:, 2, ix, iy, iz), &
                           s_field(:, 1, ix, iy, iz), s_field(:, 2, ix, iy, iz)
                     end do
                  end do
               end do
               edat = edat / dble(griddim(1) * griddim(2) * 2.d0)
               call effective_refractive_index(griddim(3), edat, grid_spacing(3), rieff, e0)
               nf_eff_ref_index = rieff
               close (file_unit)
!                  write(*,'('' field fit ri:'',2es12.5)') rieff
               deallocate (edat)
            end if
         end if
         if (sphere_cluster%effective_medium_simulation) then
            ctemp = sphere_cluster%effective_ref_index
            call parallel_broadcast(mpi_rank=0, mpi_number=1, &
                                    send_buffer=ctemp)
            sphere_cluster%effective_ref_index = ctemp(1)
         end if
      end do

   end subroutine run_configuration_average

   subroutine run_random_orientation_configuration_average()
      implicit none
      integer :: file_unit, rank, numprocs, itemp(1), n, &
                 numprocsperconfig, configcolor, configgroup, configcomm, configrank, config0comm, nconfigave, nsend
      real(real64) :: time1, timet
      real(real64), allocatable :: tpos(:, :)
      complex(real64) :: ctemp(1)
      character(len=256) :: tmatchar1, tmatchar2
      data tmatchar1, tmatchar2/'tmat-', '.tmp'/
      first_run = .false.
      call parallel_rank(mpi_rank=rank)
      call parallel_size(mpi_size=numprocs)
!         if(rank.ne.0) light_up=.false.
      local_rank = rank
      global_rank = rank

      if (max_iterations .le. 0) then
         numprocsperconfig = 2
      else
         numprocsperconfig = 4
      end if
      n_configuration_groups = numprocs / numprocsperconfig
      n_configuration_groups = max(n_configuration_groups, 1)
      configcolor = floor(dble(n_configuration_groups * rank) / dble(numprocs))
      configgroup = configcolor
      call parallel_split( &
         mpi_color=configcolor, mpi_key=rank, &
         mpi_new_comm=configcomm)
      call parallel_rank( &
         mpi_rank=configrank, &
         mpi_comm=configcomm)
      configcolor = configrank
      call parallel_split( &
         mpi_color=configcolor, mpi_key=rank, &
         mpi_new_comm=config0comm)
      random_configuration = .true.
!singleorigin=.true.

      call execute_simulation(print_output=.false., set_t_matrix_order=.true., dry_run=.true.)
      if (runtime_failed()) return
      call compose_group_filename(tmatchar1, configgroup, tmatchar2, t_matrix_output_file)

      if (allocated(q_eff_ave)) deallocate (q_eff_ave, q_eff_tot_ave, q_vabs_ave, sphere_position_ave)
      if (allocated(mean_t_ave)) deallocate (mean_t_ave)
      if (allocated(scat_mat_exp_coef_ave)) deallocate (scat_mat_exp_coef_ave)
      if (allocated(coh_scat_mat_exp_coef_ave)) deallocate (coh_scat_mat_exp_coef_ave)
      allocate (q_eff_ave(3, qeff_dim, sphere_cluster%number_spheres), q_eff_tot_ave(3, qeff_dim), q_vabs_ave(qeff_dim, sphere_cluster%number_spheres), &
                sphere_position_ave(3, sphere_cluster%number_spheres), mean_t_ave(2, sphere_cluster%t_matrix_order), &
                scat_mat_exp_coef_ave(4, 4, 0:2 * sphere_cluster%t_matrix_order), coh_scat_mat_exp_coef_ave(4, 4, 0:2 * sphere_cluster%t_matrix_order))
      q_eff_ave = 0.d0
      q_eff_tot_ave = 0.d0
      q_vabs_ave = 0.d0
      sphere_position_ave = 0.d0
      mean_t_ave = 0.d0
      scat_mat_exp_coef_ave = 0.d0
      coh_scat_mat_exp_coef_ave = 0.d0
      if (calculate_scattering_matrix) then
         if (allocated(scat_mat_ave)) deallocate (scat_mat_ave)
         allocate (scat_mat_ave(scat_mat_mdim, scat_mat_ldim:scat_mat_udim))
         scat_mat_ave = 0.d0
      end if

      nconfigave = 0
      do random_configuration_number = 1, ceiling(dble(number_configurations) / dble(n_configuration_groups))

         if (rank .eq. 0) then
            if (random_configuration_number .eq. 1) then
               call print_run_variables(sphere_cluster%run_print_unit)
               call open_output_file(output_file, file_unit, append=.true.)
               if (runtime_failed()) return
               call print_run_variables(file_unit)
               close (file_unit)
            end if
            write (sphere_cluster%run_print_unit, '('' configuration averaging, samples:'',i5,''-'',i5)') &
               (random_configuration_number - 1) * n_configuration_groups + 1, &
               random_configuration_number * n_configuration_groups
         end if

         if (rank .eq. 0) time1 = parallel_wall_time()

         call execute_simulation(print_output=.false., set_t_matrix_order=.false., mpi_comm=configcomm)
         if (runtime_failed()) return

         if (configrank .eq. 0) then
            if (rank .eq. 0) solution_time = parallel_wall_time() - time1
            allocate (tpos(3, sphere_cluster%number_spheres))
            call cartesian_vectors_to_spherical(sphere_cluster%number_spheres, sphere_cluster%sphere_position(:, 1:sphere_cluster%number_spheres), &
                                                tpos(:, 1:sphere_cluster%number_spheres))
            q_eff_ave = q_eff_ave + q_eff
            q_eff_tot_ave = q_eff_tot_ave + q_eff_tot
            q_vabs_ave = q_vabs_ave + q_vabs
            sphere_position_ave = sphere_position_ave + tpos
            mean_t_ave(:, 1:sphere_cluster%t_matrix_order) = mean_t_ave(:, 1:sphere_cluster%t_matrix_order) + mean_t(:, 1:sphere_cluster%t_matrix_order)
            if (calculate_scattering_matrix) then
               scat_mat_ave = scat_mat_ave + scat_mat
               scat_mat_exp_coef_ave = scat_mat_exp_coef_ave + scat_mat_exp_coef
               coh_scat_mat_exp_coef_ave = coh_scat_mat_exp_coef_ave + coh_scat_mat_exp_coef
            end if
            if (rank .eq. 0) solution_time_ave = solution_time_ave + solution_time
            deallocate (tpos)
         end if

         nsend = 3 * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=sphere_position_ave, &
            receive_buffer=sphere_cluster%sphere_position, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 3 * qeff_dim * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=q_eff_ave, &
            receive_buffer=q_eff, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 3 * qeff_dim
         call parallel_reduce_sum( &
            send_buffer=q_eff_tot_ave, &
            receive_buffer=q_eff_tot, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = qeff_dim * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=q_vabs_ave, &
            receive_buffer=q_vabs, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 2 * sphere_cluster%t_matrix_order
         call parallel_reduce_sum( &
            send_buffer=mean_t_ave, &
            receive_buffer=mean_t, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         if (calculate_scattering_matrix) then
            nsend = scat_mat_mdim * (scat_mat_udim - scat_mat_ldim + 1)
            call parallel_reduce_sum( &
               send_buffer=scat_mat_ave, &
               receive_buffer=scat_mat, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend = 16 * (2 * sphere_cluster%t_matrix_order + 1)
            call parallel_reduce_sum( &
               send_buffer=scat_mat_exp_coef_ave, &
               receive_buffer=scat_mat_exp_coef, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            call parallel_reduce_sum( &
               send_buffer=coh_scat_mat_exp_coef_ave, &
               receive_buffer=coh_scat_mat_exp_coef, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if

         nconfigave = nconfigave + n_configuration_groups

         if (rank .eq. 0) then
            sphere_cluster%sphere_position = sphere_cluster%sphere_position / dble(nconfigave)
            q_eff = q_eff / dble(nconfigave)
            q_vabs = q_vabs / dble(nconfigave)
            q_eff_tot = q_eff_tot / dble(nconfigave)
            mean_t(:, 1:sphere_cluster%t_matrix_order) = mean_t(:, 1:sphere_cluster%t_matrix_order) / dble(nconfigave)
            if (calculate_scattering_matrix) then
               scat_mat = scat_mat / dble(nconfigave)
               scat_mat_exp_coef = scat_mat_exp_coef / dble(nconfigave)
               coh_scat_mat_exp_coef = coh_scat_mat_exp_coef / dble(nconfigave)
            end if
            solution_time = solution_time_ave / dble(random_configuration_number)
            call fit_effective_refractive_index(mean_t, sphere_cluster%effective_ref_index, fit_radius, fit_stat)
            call print_calculation_results(output_file)
         end if
!            if(random_configuration_host) then
!               ctemp=sphere_cluster%effective_ref_index
!               if(rank.eq.0) write(*,'('' new ri:'',2es12.4)') sphere_cluster%effective_ref_index
!               call parallel_broadcast(mpi_rank=0,mpi_number=1, &
!                 send_buffer=ctemp)
!               layer_ref_index(0)=ctemp(1)
!            endif
!if(rank.eq.0) then
!write(*,'(4es12.4)') layer_ref_index(0),sphere_cluster%sphere_ref_index(1,0)
!endif
      end do
   end subroutine run_random_orientation_configuration_average

   subroutine run_incidence_average()
      implicit none
      logical :: singleorigin, prancon, aa, soe, iframe, cuds
      integer :: file_unit, rank, numprocs, &
                 numprocsperconfig, configcolor, configgroup, configcomm, configrank, config0comm, nconfigave, nsend
      real(real64) :: time1, timet
      real(real64), allocatable :: texpcoef(:, :, :)
      character(len=256) :: sdatfile
      first_run = .false.
      call parallel_rank(mpi_rank=rank)
      call parallel_size(mpi_size=numprocs)
!         if(rank.ne.0) light_up=.false.
      local_rank = rank
      global_rank = rank
      sdatfile = sphere_data_input_file
      prancon = print_random_configuration
      print_random_configuration = .true.
      aa = azimuthal_average
      soe = single_origin_expansion
      iframe = incident_frame
      cuds = calculate_up_down_scattering
      azimuthal_average = .true.
      single_origin_expansion = .true.
      incident_frame = .true.
      calculate_up_down_scattering = .false.

      if (max_iterations .lt. 0) then
         numprocsperconfig = 2
      else
         numprocsperconfig = 4
      end if
      n_configuration_groups = numprocs / numprocsperconfig
      n_configuration_groups = max(n_configuration_groups, 1)
      configcolor = floor(dble(n_configuration_groups * rank) / dble(numprocs))
      configgroup = configcolor
      call parallel_split( &
         mpi_color=configcolor, mpi_key=rank, &
         mpi_new_comm=configcomm)
      call parallel_rank( &
         mpi_rank=configrank, &
         mpi_comm=configcomm)
      configcolor = configrank
      call parallel_split( &
         mpi_color=configcolor, mpi_key=rank, &
         mpi_new_comm=config0comm)
      singleorigin = number_plane_boundaries .eq. 0 .and. single_origin_expansion
      incident_beta_specified = .true.

      call execute_simulation(print_output=.false., set_t_matrix_order=.true., dry_run=.true.)
      if (runtime_failed()) return

      if (trim(sphere_data_input_file) .eq. 'random_configuration') then
         sphere_data_input_file = 'random_configuration.pos'
      end if

      if (allocated(q_eff_ave)) deallocate (q_eff_ave, q_eff_tot_ave, q_vabs_ave, boundary_sca_ave, &
                                            boundary_ext_ave, dif_boundary_sca)
      allocate (q_eff_ave(3, qeff_dim, sphere_cluster%number_spheres), q_eff_tot_ave(3, qeff_dim), q_vabs_ave(qeff_dim, sphere_cluster%number_spheres), &
                boundary_sca_ave(2, 0:1), boundary_ext_ave(2, 0:1), dif_boundary_sca(2, 0:1))
      q_eff_ave = 0.d0
      q_eff_tot_ave = 0.d0
      q_vabs_ave = 0.d0
      pl_sca_ave = 0.d0
      boundary_sca_ave = 0.d0
      boundary_ext_ave = 0.d0
      solution_time_ave = 0.d0
      if (singleorigin) then
         if (allocated(amnp_0_ave)) deallocate (amnp_0_ave, scat_mat_exp_coef_ave)
         allocate (amnp_0_ave(2 * sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2), 2), &
                   scat_mat_exp_coef_ave(16, 0:2 * sphere_cluster%t_matrix_order, 4))
         amnp_0_ave = 0.d0
         scat_mat_exp_coef_ave = 0.d0
      end if
      if (calculate_scattering_matrix) then
         if (allocated(scat_mat_ave)) deallocate (scat_mat_ave, dif_scat_mat)
         allocate (scat_mat_ave(scat_mat_mdim, scat_mat_ldim:scat_mat_udim), &
                   dif_scat_mat(scat_mat_mdim, scat_mat_ldim:scat_mat_udim))
         scat_mat_ave = 0.d0
      end if

      nconfigave = 0
      do incident_direction_number = 1, ceiling(dble(number_incident_directions) / dble(n_configuration_groups))

         if (rank .eq. 0) then
            if (incident_direction_number .eq. 1) then
               call print_run_variables(sphere_cluster%run_print_unit)
               call open_output_file(output_file, file_unit, append=.true.)
               if (runtime_failed()) return
               call print_run_variables(file_unit)
               close (file_unit)
            end if
            write (sphere_cluster%run_print_unit, '('' incidence averaging, samples:'',i5,''-'',i5)') &
               (incident_direction_number - 1) * n_configuration_groups + 1, &
               incident_direction_number * n_configuration_groups
         end if

         call sample_incident_direction(mpi_comm=configcomm)

         if (rank .eq. 0) time1 = parallel_wall_time()

         call execute_simulation(print_output=.false., set_t_matrix_order=.false., mpi_comm=configcomm)
         if (runtime_failed()) return

         if (singleorigin) then
            amnp_0_ave = amnp_0_ave + amnp_0
         end if

         if (configrank .eq. 0) then
            if (rank .eq. 0) solution_time = parallel_wall_time() - time1
            q_eff_ave = q_eff_ave + q_eff
            q_eff_tot_ave = q_eff_tot_ave + q_eff_tot
            q_vabs_ave = q_vabs_ave + q_vabs
            pl_sca_ave = pl_sca_ave + pl_sca
            if (calculate_up_down_scattering) boundary_sca_ave = boundary_sca_ave + boundary_sca
            boundary_ext_ave = boundary_ext_ave + boundary_ext
            if (singleorigin .and. azimuthal_average .and. (.not. numerical_azimuthal_average)) &
               scat_mat_exp_coef_ave = scat_mat_exp_coef_ave + scat_mat_exp_coef
            if (calculate_scattering_matrix) then
               scat_mat_ave = scat_mat_ave + scat_mat
            end if
            if (rank .eq. 0) solution_time_ave = solution_time_ave + solution_time
         end if

         nsend = 3 * qeff_dim * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=q_eff_ave, &
            receive_buffer=q_eff, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 3 * qeff_dim
         call parallel_reduce_sum( &
            send_buffer=q_eff_tot_ave, &
            receive_buffer=q_eff_tot, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = qeff_dim * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=q_vabs_ave, &
            receive_buffer=q_vabs, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 4
         call parallel_reduce_sum( &
            send_buffer=pl_sca_ave, &
            receive_buffer=pl_sca, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         if (calculate_up_down_scattering) then
            nsend = 4
            call parallel_reduce_sum( &
               send_buffer=boundary_sca_ave, &
               receive_buffer=boundary_sca, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
         call parallel_reduce_sum( &
            send_buffer=boundary_ext_ave, &
            receive_buffer=boundary_ext, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         if (singleorigin) then
            nsend = 4 * sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2)
            call parallel_reduce_sum( &
               send_buffer=amnp_0_ave, &
               receive_buffer=amnp_0, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
         if (calculate_scattering_matrix) then
            nsend = scat_mat_mdim * (scat_mat_udim - scat_mat_ldim + 1)
            call parallel_reduce_sum( &
               send_buffer=scat_mat_ave, &
               receive_buffer=scat_mat, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
         if (singleorigin .and. azimuthal_average .and. (.not. numerical_azimuthal_average)) then
            nsend = 16 * 4 * (2 * sphere_cluster%t_matrix_order + 1)
            call parallel_reduce_sum( &
               send_buffer=scat_mat_exp_coef_ave, &
               receive_buffer=scat_mat_exp_coef, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if

         nconfigave = nconfigave + n_configuration_groups

         if (singleorigin) then
            if (rank .eq. 0 .and. print_timings) then
               timet = parallel_wall_time()
               write (sphere_cluster%run_print_unit, '('' calculating diffuse field:'')', advance='no')
            end if
            amnp_0 = amnp_0 / dble(nconfigave)
            if (azimuthal_average .and. (.not. numerical_azimuthal_average)) then
               allocate (texpcoef(16, 0:2 * sphere_cluster%t_matrix_order, 4))
               texpcoef = scat_mat_exp_coef / dble(nconfigave)
               call fixed_orientation_scattering_matrix_expansion( &
                  sphere_cluster%t_matrix_order, amnp_0, scat_mat_exp_coef(:, :, 1), scat_mat_exp_coef(:, :, 2), &
                  scat_mat_exp_coef(:, :, 3), scat_mat_exp_coef(:, :, 4), mpi_comm=configcomm)
            end if
            if (calculate_scattering_matrix) then
               call compute_scattering_matrix(amnp_0, dif_scat_mat, mpi_comm=configcomm)
            end if
            if (azimuthal_average .and. (.not. numerical_azimuthal_average)) then
               scat_mat_exp_coef = texpcoef - scat_mat_exp_coef * (dble(sphere_cluster%number_spheres - 1) / dble(sphere_cluster%number_spheres))**2
               deallocate (texpcoef)
            end if
            call common_origin_hemispherical_scattering(amnp_0, dif_boundary_sca)
      if (rank .eq. 0 .and. print_timings) write (sphere_cluster%run_print_unit, '('' completed, '',es12.4,'' sec'')') parallel_wall_time() - timet
         end if
         if (allocated(amnp_0)) deallocate (amnp_0)

         if (rank .eq. 0) then
            q_eff = q_eff / dble(nconfigave)
            q_vabs = q_vabs / dble(nconfigave)
            q_eff_tot = q_eff_tot / dble(nconfigave)
            pl_sca = pl_sca / dble(nconfigave)
            if (calculate_up_down_scattering) boundary_sca = boundary_sca / dble(nconfigave)
            boundary_ext = boundary_ext / dble(nconfigave)
!               if(singleorigin) dif_boundary_sca=boundary_sca-dif_boundary_sca
            if (singleorigin) dif_boundary_sca = boundary_sca &
                                  - dif_boundary_sca * dble(sphere_cluster%number_spheres - 1) / dble(sphere_cluster%number_spheres)
            if (calculate_scattering_matrix) then
               scat_mat = scat_mat / dble(nconfigave)
! experiment
               if (singleorigin) dif_scat_mat = scat_mat - dif_scat_mat * (dble(sphere_cluster%number_spheres - 1) / dble(sphere_cluster%number_spheres))**2

!                  dif_scat_mat=scat_mat-dif_scat_mat
            end if
            solution_time = solution_time_ave / dble(random_configuration_number)
            call print_calculation_results(output_file)
         end if
      end do
      sphere_data_input_file = sdatfile
      print_random_configuration = prancon
      azimuthal_average = aa
      single_origin_expansion = soe
      incident_frame = iframe
      calculate_up_down_scattering = cuds
   end subroutine run_incidence_average

   subroutine common_origin_scattering_cross_section(n, a, c)
      implicit none
      integer :: n
      real(real64) :: c(1)
      complex(real64) :: a(4 * n * (n + 2))
      c(1) = sum(a(:) * conjg(a(:)))
   end subroutine common_origin_scattering_cross_section

   subroutine subtract_fixed_sphere_from_common_origin()
      implicit none
      integer :: i, i1, mnp0, mnp1, n, m, p
      real(real64) :: fn
      complex(real64) :: a(2, 2), b(2, 2)

      fn = dble(sphere_cluster%number_spheres - 1) / dble(sphere_cluster%number_spheres)
      fn = 1.d0
      do i = 1, sphere_cluster%number_spheres
         if (sum(sphere_cluster%sphere_position(:, i)**2) .lt. 1.d-7) then
            i1 = i
            exit
         end if
      end do
      do n = 1, sphere_cluster%sphere_order(i1)
         do m = -n, n
            do p = 1, 2
               mnp1 = polarized_mode_index(m, n, p, sphere_cluster%sphere_order(i1), 2)
               a(p, :) = amnp_s(mnp1 + sphere_cluster%sphere_offset(i1), :)
            end do
            b(1, :) = a(1, :) + a(2, :)
            b(2, :) = a(1, :) - a(2, :)
            do p = 1, 2
               mnp0 = polarized_mode_index(m, n, p, sphere_cluster%t_matrix_order, 2)
               amnp_0(mnp0, :) = amnp_0(mnp0, :) - fn * b(p, :)
!                  amnp_0=amnp_0*dble(sphere_cluster%number_spheres)/dble(sphere_cluster%number_spheres-1)
            end do
         end do
      end do
   end subroutine subtract_fixed_sphere_from_common_origin

   subroutine calculate_surface_absorptance()
      implicit none
      integer :: i
      real(real64) :: rsamp, r, asamp
      if (periodic_lattice) then
         surface_absorptance(1:2) = q_eff_tot(2, 2:3)
      else
         surface_absorptance = 0.d0
         asamp = absorption_sample_radius * length_scale_factor
         rsamp = min(sphere_cluster%cross_section_radius, asamp)
         do i = 1, sphere_cluster%number_spheres
            if (sphere_cluster%host_sphere(i) .ne. 0) cycle
            r = sqrt(sum((sphere_cluster%sphere_position(1:2, i) - sphere_cluster%cluster_origin(1:2))**2))
            if (r .le. asamp) then
               surface_absorptance(1:2) = surface_absorptance(1:2) &
                                          + q_eff(2, 2:3, i) * (sphere_cluster%sphere_radius(i) / rsamp)**2
            end if
         end do
      end if
   end subroutine calculate_surface_absorptance

   subroutine sample_incident_direction(mpi_comm)
      implicit none
      integer :: mpicomm, rank, numprocs
      integer, optional :: mpi_comm
      real(real64) :: rnum(2), sbuf(2), cbeta
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      if (rank .eq. 0) then
         call random_number(rnum)
         cbeta = 2.d0 * rnum(1) - 1.d0
         sbuf(1) = acosd(cbeta)
         sbuf(2) = 360.d0 * rnum(2)
      end if
      if (numprocs .gt. 1) then
         call parallel_broadcast(mpi_number=2, &
                                 mpi_rank=0, send_buffer=sbuf, mpi_comm=mpicomm)
      end if
      incident_beta_deg = sbuf(1)
      incident_alpha_deg = sbuf(2)
   end subroutine sample_incident_direction

end module input_execution
