module simulation_execution
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
   private
   public :: execute_simulation
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
      averagerun = simulation_config%configuration_average .or. simulation_config%incidence_average
      simulation_config%calculate_up_down_scattering = simulation_config%input_calculate_up_down_scattering
      if (simulation_config%reflection_model) then
         simulation_config%calculate_up_down_scattering = .true.
         simulation_config%incident_frame = .false.
      end if

      simulation_config%first_run = .false.
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      if (simulation_config%solver%solution_method(1:1) /= 'i' .and. simulation_config%solver%solution_method(1:1) /= 'd') then
   call set_runtime_error("Unknown solution method '"//trim(simulation_config%solver%solution_method)//"'; use iteration or direct")
         return
      end if
      if (simulation_config%solver%max_iterations > 0 .and. simulation_config%solver%solution_epsilon <= 0.0_real64) then
         call set_runtime_error('solution_epsilon must be positive')
         return
      end if
!         if(rank.ne.0) light_up=.false.
      simulation_result%local_rank = rank
      global_rank = rank
      if ((.not. simulation_config%configuration_average) .and. (.not. simulation_config%incidence_average)) simulation_config%number_configuration_groups = 1
!         simulation_config%random_configuration=(trim(simulation_config%output%sphere_data_file).eq.'simulation_config%random_configuration')
      if (simulation_config%random_configuration) then
         if (simulation_config%auto_target_radius .and. target_shape .eq. 2) then
            sphere_cluster%number_spheres = simulation_config%input_number_spheres
            target_dimensions(1:3) = (dble(sphere_cluster%number_spheres) / simulation_config%sphere_volume_fraction)**(1.d0 / 3.d0)
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
            if (simulation_config%number_spheres_specified) then
               sphere_cluster%number_spheres = simulation_config%input_number_spheres
            simulation_config%sphere_volume_fraction = dble(simulation_config%input_number_spheres) * four_pi_over_three / targetvol
            else
               sphere_cluster%number_spheres = ceiling(targetvol * simulation_config%sphere_volume_fraction) / (four_pi_over_three)
            end if
         end if
         if (target_shape .eq. 2 .and. simulation_config%random_configuration_host) then
            sphere_cluster%number_spheres = sphere_cluster%number_spheres + 1
         end if
      else
         sphere_cluster%number_spheres = simulation_config%input_number_spheres
      end if
      if (simulation_config%medium_ref_index_specified) then
         if (simulation_config%medium_reim_ref_index_specified) then
            layer_ref_index(0) = cmplx(simulation_config%medium_re_ref_index, simulation_config%medium_im_ref_index, kind=real64)
         else
            layer_ref_index(0) = simulation_config%medium_ref_index
         end if
      end if
      if (allocated(sphere_cluster%sphere_radius)) then
         deallocate (sphere_cluster%sphere_radius, &
                     sphere_cluster%sphere_position, &
                     sphere_cluster%sphere_ref_index, &
                     sphere_cluster%material_model, &
                     sphere_cluster%host_sphere, &
                     sphere_cluster%number_field_expansions, &
                     simulation_config%sphere_excitation_switch, &
                     simulation_result%sphere_index)
      end if
      allocate (sphere_cluster%sphere_radius(sphere_cluster%number_spheres), &
                sphere_cluster%sphere_position(3, sphere_cluster%number_spheres), &
                sphere_cluster%sphere_ref_index(2, 0:sphere_cluster%number_spheres), &
                sphere_cluster%material_model(sphere_cluster%number_spheres), &
                sphere_cluster%host_sphere(sphere_cluster%number_spheres), &
                sphere_cluster%number_field_expansions(sphere_cluster%number_spheres), &
                simulation_config%sphere_excitation_switch(sphere_cluster%number_spheres), &
                simulation_result%sphere_index(sphere_cluster%number_spheres))
      sphere_cluster%material_model = material_dielectric
      if (simulation_config%random_configuration) then
         if (simulation_config%output%print_timings .and. mstm_global_rank .eq. 0) then
            write (sphere_cluster%run_print_unit, '('' generating random configuration:'')', advance='no')
            timet = parallel_wall_time()
         end if
         call generate_random_configuration(mpi_comm=mpicomm, skip_diffusion=dryrun)
         if (runtime_failed()) return
!            call generate_random_configuration(mpi_comm=mpicomm)
         if (simulation_config%output%print_timings .and. mstm_global_rank .eq. 0) then
            write (sphere_cluster%run_print_unit, '('' completed, time:'',es12.5,'' s'')') parallel_wall_time() - timet
         end if
         if (rank .eq. 0) then
!               if(simulation_config%output%print_random_configuration.and.(.not.simulation_config%configuration_average)) then
            if (simulation_config%output%print_random_configuration .and. mstm_global_rank .eq. 0) then
               call open_output_file(trim(simulation_config%output%random_configuration_file), file_unit)
               if (runtime_failed()) return
               do i = 1, sphere_cluster%number_spheres
                  write (file_unit, '(4es13.5)') sphere_cluster%sphere_position(:, i) / simulation_config%length_scale_factor, &
                     sphere_cluster%sphere_radius(i) / simulation_config%length_scale_factor
               end do
               close (file_unit)
            end if
         end if
      else
         call read_sphere_data_input_file(mpi_comm=mpicomm)
         if (runtime_failed()) return
      end if

      simulation_result%position_shift = (/simulation_config%x_shift, simulation_config%y_shift, simulation_config%z_shift/) * simulation_config%length_scale_factor
      if (simulation_config%shifted_sphere .gt. sphere_cluster%number_spheres) simulation_config%shifted_sphere = 0
      if (any(simulation_result%position_shift .ne. 0.d0)) then
         if (simulation_config%shifted_sphere .eq. 0) then
            do i = 1, sphere_cluster%number_spheres
               sphere_cluster%sphere_position(:, i) = sphere_cluster%sphere_position(:, i) + simulation_result%position_shift(:)
            end do
         else
           sphere_cluster%sphere_position(:, simulation_config%shifted_sphere) = sphere_cluster%sphere_position(:, simulation_config%shifted_sphere) + simulation_result%position_shift(:)
         end if
      end if

      sphere_cluster%sphere_ref_index(:, 0) = layer_ref_index(0)
      if (periodic_lattice) then
         if (simulation_config%random_configuration .and. target_shape .eq. 0) then
            cell_width(1:2) = target_dimensions(1:2) * 2.d0 * simulation_config%length_scale_factor
         else
            if (simulation_config%square_cell) then
               cell_width = simulation_config%input_cell_width_x * simulation_config%length_scale_factor
            else
               cell_width = simulation_config%input_cell_width * simulation_config%length_scale_factor
            end if
         end if
      end if

      plane_surface_present = number_plane_boundaries .gt. 0
      layer_thickness = simulation_config%input_layer_thickness * simulation_config%length_scale_factor
      call initialize_plane_boundaries()

      if (simulation_config%move_to_front .and. plane_surface_present) then
         zext = maxval(sphere_cluster%sphere_position(3, :) + sphere_cluster%sphere_radius(:))
         if (zext .gt. 0.d0) sphere_cluster%sphere_position(3, :) = sphere_cluster%sphere_position(3, :) - zext
      end if
      if (simulation_config%move_to_back .and. plane_surface_present) then
         zext = minval(sphere_cluster%sphere_position(3, :) - sphere_cluster%sphere_radius(:))
         if (zext .lt. plane_boundary_position(number_plane_boundaries)) &
            sphere_cluster%sphere_position(3, :) = sphere_cluster%sphere_position(3, :) - zext + plane_boundary_position(number_plane_boundaries)
      end if

      stopit = .false.
      if (simulation_config%random_orientation) then
         if (number_plane_boundaries .gt. 0) then
            if (rank .eq. 0) write (sphere_cluster%run_print_unit, '('' random orientation requires number_plane_boundaries=0'')')
            stopit = .true.
         end if
         if (periodic_lattice) then
            if (rank .eq. 0) write (sphere_cluster%run_print_unit, '('' random orientation and periodic lattice incompatible'')')
            stopit = .true.
         end if
      end if

     sphere_cluster%fft_translation_option = (simulation_config%fft%enabled .and. sphere_cluster%number_spheres .ge. simulation_config%fft%minimum_spheres)
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

!         if(sphere_cluster%fft_translation_option) simulation_config%single_origin_expansion=.true.

      if (stopit) return
      if (light_up) then
         write (*, '('' s2 '',i3)') mstm_global_rank
         flush (6)
      end if
      sphere_cluster%cluster_origin = 0.d0
      call sphere_cluster%find_hosts()

      if (simulation_config%configuration_average .and. (target_shape .eq. 2) &
          .and. simulation_config%random_configuration_host .and. simulation_config%auto_target_radius) then
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
      call validate_material_configuration()
      if (runtime_failed()) return
      call calculate_mie_coefficients(simulation_config%solver%mie_epsilon)
      if (runtime_failed()) return
      call initialize_numerical_tables(sphere_cluster%max_mie_order)
      if (light_up) then
         write (*, '('' s4 '',i3)') mstm_global_rank
         flush (6)
      end if

      singleorigin = number_plane_boundaries .eq. 0 .and. simulation_config%single_origin_expansion
      iframe = singleorigin .and. simulation_config%incident_frame

      sphere_cluster%cluster_origin = 0.d0
      if (singleorigin .or. simulation_config%random_orientation .or. .true.) then
         if (allocated(sphere_cluster%translation_order)) deallocate (sphere_cluster%translation_order)
         allocate (sphere_cluster%translation_order(sphere_cluster%number_spheres))
    sphere_cluster%translation_order(1:sphere_cluster%number_spheres) = sphere_cluster%sphere_order(1:sphere_cluster%number_spheres)
         sphere_cluster%cluster_origin = 0.d0
         if ((.not. simulation_config%configuration_average) .and. (.not. simulation_config%incidence_average)) then
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
               sphere_cluster%cluster_origin = sphere_cluster%gaussian_beam_focal_point * simulation_config%length_scale_factor
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
                                               simulation_config%solver%translation_epsilon, sphere_cluster%translation_order(i))
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
         elseif (simulation_config%reflection_model) then
            if (simulation_config%random_configuration) then
               if (target_shape .eq. 0) then
     sphere_cluster%cross_section_radius = simulation_config%length_scale_factor * 2.d0 * sqrt(product(target_dimensions(1:2)) / pi)
               elseif (target_shape .ge. 1) then
                  sphere_cluster%cross_section_radius = simulation_config%length_scale_factor * target_dimensions(1)
               end if
            else
               sphere_cluster%cross_section_radius = sqrt(product(sphere_cluster%sphere_max_position(1:2) - sphere_cluster%sphere_min_position(1:2)) / pi)
            end if
            sphere_cluster%cross_section_radius = min(sphere_cluster%cross_section_radius, simulation_config%length_scale_factor * simulation_config%excitation_radius)
         else
            sphere_cluster%cross_section_radius = sphere_cluster%vol_radius
         end if
      end if

      if (simulation_config%auto_absorption_sample_radius .and. simulation_config%random_configuration) then
         simulation_config%absorption_sample_radius = simulation_config%absorption_sample_radius_fraction * target_dimensions(1)
      end if

      if (sphere_cluster%fft_translation_option) then
         lochost = 0
         if (simulation_config%random_configuration) then
            tmin = -target_dimensions * simulation_config%length_scale_factor
            tmax = target_dimensions * simulation_config%length_scale_factor
            if (dryrun) call fft_plan%clear(clear_h=.true.)
            if (target_shape .eq. 2 .and. simulation_config%random_configuration_host) then
               lochost = sphere_cluster%number_spheres
            end if
         else
            tmin = sphere_cluster%sphere_min_position
            tmax = sphere_cluster%sphere_max_position
            if (.not. averagerun) call fft_plan%clear(clear_h=.true.)
         end if
         call fft_plan%configure(simulation_config%fft%cell_volume_fraction, target_min=tmin, target_max=tmax, &
                                 d_specified=simulation_config%fft%cell_size_specified, local_host=lochost, &
                       requested_cell_size=simulation_config%fft%cell_size, requested_node_order=simulation_config%fft%node_order, &
                                 requested_neighbor_model=simulation_config%fft%neighbor_node_model)
      end if

      simulation_config%sphere_excitation_switch = .true.
      do i = 1, sphere_cluster%number_spheres
         if (sphere_cluster%host_sphere(i) .ne. 0) cycle
         if (simulation_config%random_configuration) then
            if (target_shape .le. 1) then
               rtran = sqrt(sum((sphere_cluster%sphere_position(1:2, i) - sphere_cluster%cluster_origin(1:2))**2))
            else
               rtran = sqrt(sum((sphere_cluster%sphere_position(1:3, i) - sphere_cluster%cluster_origin(1:3))**2))
            end if
         else
            rtran = sqrt(sum((sphere_cluster%sphere_position(1:3, i) - sphere_cluster%cluster_origin(1:3))**2))
         end if
         if (simulation_config%excitation_radius .gt. 0.d0) then
            simulation_config%sphere_excitation_switch(i) = rtran .le. simulation_config%excitation_radius * simulation_config%length_scale_factor
         else
            simulation_config%sphere_excitation_switch(i) = i .le. -int(simulation_config%excitation_radius)
         end if
      end do
      if (simulation_config%excitation_radius .eq. 0.d0) then
         simulation_config%sphere_excitation_switch = .false.
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
         simulation_config%sphere_excitation_switch(itemp(1)) = .true.
      end if
      if (light_up) then
         write (*, '('' s6 '',i3)') mstm_global_rank
         flush (6)
      end if
      if (simulation_config%random_orientation) then
         simulation_result%efficiency_dimension = 1
         if (simulation_config%calculate_scattering_matrix) then
            simulation_result%scattering_matrix_upper_bound = floor(180.00001d0 / simulation_config%scattering_map_increment)
            simulation_result%scattering_matrix_dimension = 16
            simulation_result%scattering_matrix_lower_bound = 0
            simulation_config%scattering_matrix_angle_minimum = 0.d0
            simulation_config%scattering_matrix_angle_maximum = 180.d0
         end if
      else
         if (simulation_config%incident_beta_specified) then
            simulation_result%incident_beta = simulation_config%incident_beta_degrees * degrees_to_radians
            if (simulation_config%incident_beta_degrees .le. 90.d0) then
               simulation_config%incident_direction = 1
               simulation_result%incident_sin_beta = sind(simulation_config%incident_beta_degrees) / dble(layer_ref_index(0))
            else
               simulation_config%incident_direction = 2
               simulation_result%incident_sin_beta = sind(simulation_config%incident_beta_degrees) &
                                                     / dble(layer_ref_index(number_plane_boundaries))
            end if
         else
            simulation_result%incident_beta = 0.d0
         end if
         if (simulation_config%incidence_average) then
            simulation_result%efficiency_dimension = 1
         else
            simulation_result%efficiency_dimension = 3
         end if
         alpha = simulation_config%incident_alpha_degrees * degrees_to_radians
         call initialize_incident_field(alpha, simulation_result%incident_sin_beta, simulation_config%incident_direction)
         if (simulation_config%calculate_scattering_matrix) then
            if (allocated(simulation_result%scattering_matrix)) deallocate (simulation_result%scattering_matrix)
            if (periodic_lattice) then
               call periodic_lattice_scattering(simulation_result%solution_coefficients, simulation_result%plane_scattering, dry_run=.true., num_dirs=simulation_result%reflection_transmission_direction_counts)
   simulation_result%maximum_reflection_transmission_directions = maxval(simulation_result%reflection_transmission_direction_counts)
    if (allocated(simulation_result%reflection_transmission_vectors)) deallocate (simulation_result%reflection_transmission_vectors)
       allocate (simulation_result%reflection_transmission_vectors(2, simulation_result%maximum_reflection_transmission_directions))
               simulation_result%scattering_matrix_upper_bound = simulation_result%maximum_reflection_transmission_directions
               simulation_result%scattering_matrix_lower_bound = 1
               simulation_result%scattering_matrix_dimension = 32
            else
               if (simulation_config%scattering_map_model .eq. 0) then
                  if (number_plane_boundaries .eq. 0) then
                   simulation_result%scattering_matrix_upper_bound = floor(180.00001d0 / simulation_config%scattering_map_increment)
                     if (simulation_config%azimuthal_average) then
                        simulation_result%scattering_matrix_lower_bound = 0
                     else
                        simulation_result%scattering_matrix_lower_bound = -simulation_result%scattering_matrix_upper_bound
                     end if
                     simulation_result%scattering_matrix_dimension = 16
                     simulation_config%scattering_matrix_angle_maximum = 180.d0
                  else
                    simulation_result%scattering_matrix_upper_bound = floor(90.00001d0 / simulation_config%scattering_map_increment)
!                        simulation_result%scattering_matrix_lower_bound=-simulation_result%scattering_matrix_upper_bound
! 10-22 simulation_config%azimuthal_average applies to multiple plane boundaries
                     if (simulation_config%azimuthal_average) then
                        simulation_result%scattering_matrix_lower_bound = 0
                     else
                        simulation_result%scattering_matrix_lower_bound = -simulation_result%scattering_matrix_upper_bound
                     end if
                     simulation_result%scattering_matrix_dimension = 32
                     simulation_config%scattering_matrix_angle_maximum = 90.d0
                  end if
                  simulation_config%scattering_matrix_angle_minimum = simulation_config%scattering_matrix_angle_maximum * (simulation_result%scattering_matrix_lower_bound / simulation_result%scattering_matrix_upper_bound)
               else
                  i = 0
                  do sy = -simulation_config%scattering_map_dimension, simulation_config%scattering_map_dimension
                     do sx = -simulation_config%scattering_map_dimension, simulation_config%scattering_map_dimension
                        if (sx * sx + sy * sy .gt. simulation_config%scattering_map_dimension**2) cycle
                        i = i + 1
                     end do
                  end do
                  simulation_result%scattering_matrix_upper_bound = i
                  simulation_result%scattering_matrix_lower_bound = 1
                  simulation_result%scattering_matrix_dimension = 32
               end if
            end if
         end if
         if (periodic_lattice .or. simulation_config%reflection_model .and. (target_shape .le. 1)) then
            sphere_cluster%cross_section_radius = sphere_cluster%cross_section_radius * sqrt(cos(simulation_result%incident_beta))
         end if
      end if
      if (light_up) then
         write (*, '('' s7 '',i3)') mstm_global_rank
         flush (6)
      end if
      if (allocated(simulation_result%boundary_scattering)) deallocate (simulation_result%boundary_scattering, simulation_result%boundary_extinction)
!         allocate(simulation_result%boundary_scattering(2,0:number_plane_boundaries+1),simulation_result%boundary_extinction(2,0:number_plane_boundaries+1))
      allocate (simulation_result%boundary_scattering(2, 0:1), simulation_result%boundary_extinction(2, 0:1))

      if (allocated(simulation_result%efficiency)) deallocate (simulation_result%efficiency, simulation_result%total_efficiency, simulation_result%volume_absorption)
      allocate (simulation_result%efficiency(3, simulation_result%efficiency_dimension, sphere_cluster%number_spheres), simulation_result%total_efficiency(3, simulation_result%efficiency_dimension), simulation_result%volume_absorption(simulation_result%efficiency_dimension, sphere_cluster%number_spheres))
      if (simulation_config%calculate_scattering_matrix) then
         if (allocated(simulation_result%scattering_matrix)) deallocate (simulation_result%scattering_matrix)
         allocate (simulation_result%scattering_matrix(simulation_result%scattering_matrix_dimension, simulation_result%scattering_matrix_lower_bound:simulation_result%scattering_matrix_upper_bound))
      end if
      if (light_up) then
         write (*, '('' s8 '',i3)') mstm_global_rank
         flush (6)
      end if
!call parallel_barrier()
      if (rank .eq. 0) time1 = parallel_wall_time()
      if (rank .eq. 0 .and. printout) then
         if (simulation_config%check_positions) call check_sphere_positions()
         call print_run_variables(sphere_cluster%run_print_unit)
         call open_output_file(simulation_config%output%output_file, file_unit, append=.true.)
         if (runtime_failed()) return
         call print_run_variables(file_unit)
         close (file_unit)
      end if

      if (dryrun) return

      if (simulation_config%random_orientation) then
         niter = simulation_config%solver%max_iterations
         timatrixfile = 'titemp.dat'
         if (allocated(simulation_result%mean_t_matrix)) deallocate (simulation_result%mean_t_matrix)
         allocate (simulation_result%mean_t_matrix(2, sphere_cluster%t_matrix_order))
         simulation_result%mean_t_matrix = 0.d0
         call solve_t_matrix(solution_method=simulation_config%solver%solution_method(1:1), &
                             solution_eps=simulation_config%solver%solution_epsilon, &
                             convergence_eps=simulation_config%solver%t_matrix_convergence_epsilon, &
                             max_iterations=niter, &
                             t_matrix_file=simulation_config%output%t_matrix_file, &
                             procs_per_soln=simulation_config%solver%t_matrix_procs_per_solution, &
                             sphere_qeff=simulation_result%efficiency, &
                             solution_status=istat, &
                             mpi_comm=mpicomm, &
                             sphere_excitation_list=simulation_config%sphere_excitation_switch)
         if (runtime_failed()) return
         if (istat /= solver_converged) then
            call set_runtime_error('T-matrix solver failed: '//trim(solver_status_message(istat)), istat)
            return
         end if
         if (sphere_cluster%fft_translation_option) call fft_plan%clear(clear_h=.true.)
         if (simulation_config%calculate_scattering_matrix) then
            if (allocated(simulation_result%scattering_matrix_expansion)) deallocate (simulation_result%scattering_matrix_expansion)
        if (allocated(simulation_result%coherent_scattering_expansion)) deallocate (simulation_result%coherent_scattering_expansion)
            allocate (simulation_result%scattering_matrix_expansion(4, 4, 0:2 * sphere_cluster%t_matrix_order), simulation_result%coherent_scattering_expansion(4, 4, 0:2 * sphere_cluster%t_matrix_order))
            nodrw = 2 * sphere_cluster%t_matrix_order
  call random_orientation_scattering_matrix(simulation_config%output%t_matrix_file, simulation_result%scattering_matrix_expansion, &
                                                      simulation_result%coherent_scattering_expansion, &
                                                      beam_width=sphere_cluster%gaussian_beam_constant, &
                                                      number_processors=simulation_config%solver%t_matrix_procs_per_solution, &
                                                      mean_t_matrix=simulation_result%mean_t_matrix, mpi_comm=mpicomm)
            simulation_result%coherent_scattering_ratio = simulation_result%coherent_scattering_expansion(1, 1, 0)
            do i = simulation_result%scattering_matrix_lower_bound, simulation_result%scattering_matrix_upper_bound
               costheta = cos(dble(i - simulation_result%scattering_matrix_lower_bound) * pi / dble(simulation_result%scattering_matrix_upper_bound - simulation_result%scattering_matrix_lower_bound))
               call evaluate_random_orientation_scattering_matrix(costheta, simulation_result%scattering_matrix_expansion, nodrw, simulation_result%scattering_matrix(:, i))
            end do
         end if
         call total_efficiency_factors(sphere_cluster%number_spheres, simulation_result%efficiency_dimension, sphere_cluster%cross_section_radius, &
                              simulation_result%efficiency, simulation_result%volume_absorption, simulation_result%total_efficiency)
      else
         if (light_up) then
            write (*, '('' s8.1 '',i3)') mstm_global_rank
            flush (6)
         end if
!call parallel_barrier()

         if (allocated(simulation_result%solution_coefficients)) deallocate (simulation_result%solution_coefficients)
         allocate (simulation_result%solution_coefficients(sphere_cluster%number_eqns, 2))
         simulation_result%solution_coefficients = 0.d0
         niter = simulation_config%solver%max_iterations
         error_codes = 0
         pl_error_codes = 0
!            sphere_cluster%recalculate_surface_matrix=.true.
         if (light_up) then
            write (*, '('' s8.2 '',i3)') mstm_global_rank
            flush (6)
         end if
!call parallel_barrier()
         if (simulation_config%output%print_timings .and. mstm_global_rank .eq. 0) then
            write (sphere_cluster%run_print_unit, '('' generating solution:'')', advance='no')
            timet = parallel_wall_time()
         end if
         call solve_fixed_orientation(alpha, simulation_result%incident_sin_beta, simulation_config%incident_direction, simulation_config%solver%solution_epsilon, niter, &
                                      simulation_result%solution_coefficients, simulation_result%efficiency, &
        simulation_result%efficiency_dimension, simulation_result%solution_error, simulation_result%solution_iterations, 1, istat, &
                                      mpi_comm=mpicomm, &
                                      excited_spheres=simulation_config%sphere_excitation_switch, &
                                      solution_method=simulation_config%solver%solution_method(1:1), &
                                      initialize_solver=.true., &
                                      reciprocal_condition=simulation_result%reciprocal_condition)
         if (runtime_failed()) return
         if (istat /= solver_converged) then
            call set_runtime_error('Fixed-orientation solver failed: '//trim(solver_status_message(istat)), istat)
            return
         end if
         if (simulation_config%output%print_timings .and. mstm_global_rank .eq. 0) then
            write (sphere_cluster%run_print_unit, '('' completed, time:'',es12.5,'' s'')') parallel_wall_time() - timet
         end if
         if (sphere_cluster%fft_translation_option) then
            call fft_plan%clear()
         end if
         if (mstm_global_rank .eq. 0. .and. ((.not. simulation_config%configuration_average) .and. (.not. simulation_config%incidence_average))) then
     write (sphere_cluster%run_print_unit, '('' solution completed: number iterations='',i5)') simulation_result%solution_iterations
            flush (sphere_cluster%run_print_unit)
         end if

         if (light_up) then
            write (*, '('' s8.3 '',i3)') mstm_global_rank
            flush (6)
         end if
!call parallel_barrier()
         if (simulation_config%output%print_timings .and. mstm_global_rank .eq. 0) then
            write (sphere_cluster%run_print_unit, '('' post processing solution:'')', advance='no')
            timet = parallel_wall_time()
         end if

         call total_efficiency_factors(sphere_cluster%number_spheres, simulation_result%efficiency_dimension, sphere_cluster%cross_section_radius, &
                              simulation_result%efficiency, simulation_result%volume_absorption, simulation_result%total_efficiency)
!            simulation_result%total_efficiency(3,:)=simulation_result%total_efficiency(1,:)-simulation_result%total_efficiency(2,:)
         csca = simulation_result%total_efficiency(3, 1) * pi * sphere_cluster%cross_section_radius**2
         if (singleorigin) then
            if (allocated(simulation_result%incident_coefficients)) deallocate (simulation_result%incident_coefficients)
      allocate (simulation_result%incident_coefficients(2 * sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2), 2))
            simulation_result%incident_coefficients = 0.d0
            if (light_up) then
               write (*, '('' s8.3.1 '',i3)') mstm_global_rank
               flush (6)
            end if
!call parallel_barrier()
            do i = 1, 2
               call merge_to_common_origin(sphere_cluster%t_matrix_order, simulation_result%solution_coefficients(:, i), simulation_result%incident_coefficients(:, i), &
                                           origin_position=sphere_cluster%cluster_origin, merge_procs=.true., &
                                           mpi_comm=mpicomm)
               if (iframe) then
                  call rotate_expansion_coefficients(alpha, simulation_result%incident_beta, 0.d0, sphere_cluster%t_matrix_order, &
                                                    sphere_cluster%t_matrix_order, simulation_result%incident_coefficients(:, i), 1)
               end if
            end do
            if (light_up) then
               write (*, '('' s8.3.2 '',i3)') mstm_global_rank
               flush (6)
            end if
!call parallel_barrier()
        if (singleorigin .and. simulation_config%azimuthal_average .and. (.not. simulation_config%numerical_azimuthal_average)) then
            if (allocated(simulation_result%scattering_matrix_expansion)) deallocate (simulation_result%scattering_matrix_expansion)
               allocate (simulation_result%scattering_matrix_expansion(16, 0:2 * sphere_cluster%t_matrix_order, 4))
               call fixed_orientation_scattering_matrix_expansion( &
                  sphere_cluster%t_matrix_order, simulation_result%incident_coefficients, simulation_result%scattering_matrix_expansion(:, :, 1), simulation_result%scattering_matrix_expansion(:, :, 2), &
   simulation_result%scattering_matrix_expansion(:, :, 3), simulation_result%scattering_matrix_expansion(:, :, 4), mpi_comm=mpicomm)
            end if
            if (light_up) then
               write (*, '('' s8.3.4 '',i3)') mstm_global_rank
               flush (6)
            end if
!call parallel_barrier()
            if (simulation_config%calculate_scattering_matrix) then
      call compute_scattering_matrix(simulation_result%incident_coefficients, simulation_result%scattering_matrix, mpi_comm=mpicomm)
            end if
            if (light_up) then
               write (*, '('' s8.3.5 '',i3)') mstm_global_rank
               flush (6)
            end if
!call parallel_barrier()
            if (sphere_cluster%gaussian_beam_constant .eq. 0.d0) then
               simulation_result%boundary_extinction = 0.d0
               ! A finite cluster without interfaces does not use plane-boundary
               ! extinction.  Avoid the spectral Green-function normalization at
               ! exact grazing incidence, where its longitudinal wave number is zero.
               if (number_plane_boundaries .gt. 0 .or. periodic_lattice .or. simulation_config%reflection_model) then
                  call boundary_extinction(simulation_result%incident_coefficients, alpha, simulation_result%incident_sin_beta, simulation_config%incident_direction, simulation_result%boundary_extinction, &
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
            if (simulation_config%calculate_scattering_matrix) then
      call compute_scattering_matrix(simulation_result%solution_coefficients, simulation_result%scattering_matrix, mpi_comm=mpicomm)
            end if
            if (sphere_cluster%gaussian_beam_constant .eq. 0.d0) then
               simulation_result%boundary_extinction = 0.d0
               if (number_plane_boundaries .gt. 0 .or. periodic_lattice .or. simulation_config%reflection_model) then
                  call boundary_extinction(simulation_result%solution_coefficients, alpha, simulation_result%incident_sin_beta, simulation_config%incident_direction, simulation_result%boundary_extinction)
               end if
            end if
         end if

         if (simulation_config%calculate_up_down_scattering) then
            if (singleorigin) then
               call hemispherical_scattering(simulation_result%incident_coefficients, .true., simulation_config%numerical_hemispherical_integration, &
                                             simulation_result%boundary_scattering, mpi_comm=mpicomm)
            else
               call hemispherical_scattering(simulation_result%solution_coefficients, .false., simulation_config%numerical_hemispherical_integration, &
                                             simulation_result%boundary_scattering, mpi_comm=mpicomm)
            end if
         end if

         if (sphere_cluster%gaussian_beam_constant .ne. 0.d0) then
            simulation_result%boundary_extinction = 0
            simulation_result%boundary_extinction(1:2, 1) = -simulation_result%total_efficiency(1, 2:3)
         end if

         if (periodic_lattice) then
            call periodic_lattice_scattering(simulation_result%solution_coefficients, simulation_result%plane_scattering)
         elseif (simulation_config%reflection_model) then
!               simulation_result%plane_scattering(:,1)=simulation_result%boundary_scattering(:,number_plane_boundaries+1)
            simulation_result%plane_scattering(:, 1) = simulation_result%boundary_scattering(:, 1)
            simulation_result%plane_scattering(:, 2) = -simulation_result%boundary_scattering(:, 0)
         end if

         if (periodic_lattice .or. simulation_config%reflection_model) then
            call calculate_surface_absorptance()
         end if

         if (number_plane_boundaries .gt. 0 .and. .not. periodic_lattice) then
!               simulation_result%evanescent_scattering(1:2)=simulation_result%total_efficiency(1,2:3)-simulation_result%total_efficiency(2,2:3)+simulation_result%boundary_scattering(1:2,0) &
!                  - simulation_result%boundary_scattering(:,number_plane_boundaries+1)
            simulation_result%evanescent_scattering(1:2) = simulation_result%total_efficiency(1, 2:3) - simulation_result%total_efficiency(2, 2:3) + simulation_result%boundary_scattering(1:2, 0) &
                                                           - simulation_result%boundary_scattering(:, 1)
         end if

         if (simulation_config%output%print_timings .and. mstm_global_rank .eq. 0) then
            write (sphere_cluster%run_print_unit, '('' completed, time:'',es12.5,'' s'')') parallel_wall_time() - timet
         end if

      end if

      if (rank .eq. 0) simulation_result%solution_time = parallel_wall_time() - time1

      call gather_error_codes(mpicomm)
      if (light_up) then
         write (*, '('' s12 '',i3)') mstm_global_rank
         flush (6)
      end if
      if (mstm_global_rank .eq. 0 .and. printout) then
         call print_calculation_results(simulation_config%output%output_file)
      end if
      if (light_up) then
         write (*, '('' s13 '',i3)') mstm_global_rank
         flush (6)
      end if
      error_codes = 0
      pl_error_codes = 0

      if (simulation_config%calculate_near_field .and. (.not. simulation_config%random_orientation)) then
         celldim = ceiling((simulation_config%near_field_plane_vertices(:, 2) - simulation_config%near_field_plane_vertices(:, 1)) / simulation_config%near_field_step_size)
         celldim = max(celldim, (/1, 1, 1/))
         if (simulation_config%configuration_average) then
    if (allocated(simulation_result%electric_field)) deallocate (simulation_result%electric_field, simulation_result%magnetic_field)
            if (rank .eq. 0) then
               allocate (simulation_result%electric_field(3, 2, celldim(1), celldim(2), celldim(3)), &
                         simulation_result%magnetic_field(3, 2, celldim(1), celldim(2), celldim(3)))
               simulation_result%electric_field = 0.d0
               simulation_result%magnetic_field = 0.d0
            end if
            call compute_near_field(simulation_result%solution_coefficients, alpha, simulation_result%incident_sin_beta, simulation_config%incident_direction, &
                                    simulation_config%near_field_plane_vertices, celldim, &
                                    incident_model=simulation_config%near_field_calculation_model, output_unit=0, &
                   e_field_array=simulation_result%electric_field, h_field_array=simulation_result%magnetic_field, mpi_comm=mpicomm)
         else
            if (simulation_config%output%append_near_field) then
               call open_output_file(simulation_config%output%near_field_file, file_unit, append=.true.)
            else
               call open_output_file(simulation_config%output%near_field_file, file_unit)
            end if
            if (runtime_failed()) return
            simulation_config%output%append_near_field = .true.
            call compute_near_field(simulation_result%solution_coefficients, alpha, simulation_result%incident_sin_beta, simulation_config%incident_direction, &
                                    simulation_config%near_field_plane_vertices, celldim, &
                         incident_model=simulation_config%near_field_calculation_model, output_unit=file_unit, output_header=.true.)
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

   subroutine calculate_surface_absorptance()
      implicit none
      integer :: i
      real(real64) :: rsamp, r, asamp
      if (periodic_lattice) then
         simulation_result%surface_absorptance(1:2) = simulation_result%total_efficiency(2, 2:3)
      else
         simulation_result%surface_absorptance = 0.d0
         asamp = simulation_config%absorption_sample_radius * simulation_config%length_scale_factor
         rsamp = min(sphere_cluster%cross_section_radius, asamp)
         do i = 1, sphere_cluster%number_spheres
            if (sphere_cluster%host_sphere(i) .ne. 0) cycle
            r = sqrt(sum((sphere_cluster%sphere_position(1:2, i) - sphere_cluster%cluster_origin(1:2))**2))
            if (r .le. asamp) then
               simulation_result%surface_absorptance(1:2) = simulation_result%surface_absorptance(1:2) &
                                            + simulation_result%efficiency(2, 2:3, i) * (sphere_cluster%sphere_radius(i) / rsamp)**2
            end if
         end do
      end if
   end subroutine calculate_surface_absorptance

end module simulation_execution
