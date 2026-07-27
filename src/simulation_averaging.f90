module simulation_averaging
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
   use simulation_execution, only: execute_simulation
   use solver, only: solver_converged, solver_status_message
   use wave_functions, only: compose_group_filename
   implicit none
   private
   public :: run_configuration_average, run_incidence_average, &
             run_random_orientation_configuration_average
contains

   subroutine run_configuration_average()
      implicit none
      logical :: singleorigin, iframe
      integer :: file_unit, rank, numprocs, m, n, p, mnp, griddim(3), ipos(3), ix, iy, iz, configuration_number, &
                 numprocsperconfig, configcolor, configgroup, configcomm, configrank, config0comm, nconfigave, nsend
      real(real64) :: time1, timet, diffac, csca(1), xspfit, rpos(3), rtemp(1)
      real(real64), allocatable :: texpcoef(:, :, :), spherical_position(:, :)
      complex(real64) :: ritemp(2), aneff, ctemp(1), rieff, e0
      complex(real64), allocatable :: pmnp0(:, :), anp0(:, :), edat(:)
      character(len=256) :: tmatchar1, tmatchar2
      data tmatchar1, tmatchar2/'tmat-', '.tmp'/
      simulation_config%first_run = .false.
      call parallel_rank(mpi_rank=rank)
      call parallel_size(mpi_size=numprocs)
!         if(rank.ne.0) light_up=.false.
      simulation_result%local_rank = rank
      global_rank = rank

      if (simulation_config%solver%max_iterations .le. 1) then
         numprocsperconfig = 2
      else
         numprocsperconfig = 4
      end if
!numprocsperconfig=2

      simulation_config%number_configuration_groups = numprocs / numprocsperconfig
      simulation_config%number_configuration_groups = max(simulation_config%number_configuration_groups, 1)
      configcolor = floor(dble(simulation_config%number_configuration_groups * rank) / dble(numprocs))
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
      simulation_config%random_configuration = .true.
      singleorigin = number_plane_boundaries .eq. 0 .and. simulation_config%single_origin_expansion
!singleorigin=.true.
      iframe = singleorigin .and. simulation_config%incident_frame

      call execute_simulation(print_output=.false., set_t_matrix_order=.true., dry_run=.true.)
      if (runtime_failed()) return

      if (allocated(simulation_result%average_efficiency)) deallocate (simulation_result%average_efficiency, simulation_result%average_total_efficiency, simulation_result%average_volume_absorption, simulation_result%average_sphere_position, simulation_result%average_boundary_scattering, &
                                       simulation_result%average_boundary_extinction, simulation_result%diffuse_boundary_scattering)
      allocate (simulation_result%average_efficiency(3, simulation_result%efficiency_dimension, sphere_cluster%number_spheres), simulation_result%average_total_efficiency(3, simulation_result%efficiency_dimension), simulation_result%average_volume_absorption(simulation_result%efficiency_dimension, sphere_cluster%number_spheres), &
                simulation_result%average_sphere_position(3, sphere_cluster%number_spheres), simulation_result%average_boundary_scattering(2, 0:1), &
                simulation_result%average_boundary_extinction(2, 0:1), simulation_result%diffuse_boundary_scattering(2, 0:1))
      simulation_result%average_efficiency = 0.d0
      simulation_result%average_total_efficiency = 0.d0
      simulation_result%average_volume_absorption = 0.d0
      simulation_result%average_sphere_position = 0.d0
      if (target_shape .eq. 2) allocate (spherical_position(3, sphere_cluster%number_spheres))
      simulation_result%average_plane_scattering = 0.d0
      simulation_result%average_boundary_scattering = 0.d0
      simulation_result%average_boundary_extinction = 0.d0
      simulation_result%average_solution_time = 0.d0
      simulation_result%average_surface_absorptance = 0.
      sphere_cluster%effective_medium_simulation = .false.
      if (singleorigin) then
         if (allocated(simulation_result%average_incident_coefficients)) deallocate (simulation_result%average_incident_coefficients, simulation_result%average_scattering_matrix_expansion, simulation_result%coherent_scattering_expansion)
         allocate (simulation_result%average_incident_coefficients(2 * sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2), 2), &
                   simulation_result%average_scattering_matrix_expansion(16, 0:2 * sphere_cluster%t_matrix_order, 4), &
                   simulation_result%coherent_scattering_expansion(16, 0:2 * sphere_cluster%t_matrix_order, 4))
         simulation_result%average_incident_coefficients = 0.d0
         simulation_result%average_scattering_matrix_expansion = 0.d0
         simulation_result%average_total_cross_section = 0.d0
         simulation_result%diffuse_cross_section_ratio = 0.d0
         if (simulation_config%effective_medium_simulation .and. target_shape .eq. 2 .and. &
             (.not. simulation_config%random_configuration_host)) then
            sphere_cluster%effective_medium_simulation = .true.
            sphere_cluster%effective_ref_index = layer_ref_index(0)
            if (simulation_config%random_configuration_host_model .eq. 1) then
               sphere_cluster%effective_cluster_radius = target_dimensions(1) * simulation_config%length_scale_factor
            elseif (simulation_config%random_configuration_host_model .eq. 2) then
           sphere_cluster%effective_cluster_radius = sphere_cluster%vol_radius / (simulation_config%sphere_volume_fraction)**0.33333
            end if
         end if
      end if
      if (simulation_config%calculate_scattering_matrix) then
         if (allocated(simulation_result%average_scattering_matrix)) deallocate (simulation_result%average_scattering_matrix, simulation_result%diffuse_scattering_matrix)
         allocate (simulation_result%average_scattering_matrix(simulation_result%scattering_matrix_dimension, simulation_result%scattering_matrix_lower_bound:simulation_result%scattering_matrix_upper_bound), &
                   simulation_result%diffuse_scattering_matrix(simulation_result%scattering_matrix_dimension, simulation_result%scattering_matrix_lower_bound:simulation_result%scattering_matrix_upper_bound))
         simulation_result%average_scattering_matrix = 0.d0
      end if
      if (simulation_config%calculate_near_field .and. configrank .eq. 0) then
         griddim = ceiling((simulation_config%near_field_plane_vertices(:, 2) - simulation_config%near_field_plane_vertices(:, 1)) / simulation_config%near_field_step_size)
         griddim = max(griddim, (/1, 1, 1/))
         if (allocated(simulation_result%average_electric_field)) deallocate (simulation_result%average_electric_field, simulation_result%average_magnetic_field, simulation_result%scattering_field, simulation_result%average_scattering_field)
         allocate (simulation_result%average_electric_field(3, 2, griddim(1), griddim(2), griddim(3)), &
                   simulation_result%average_magnetic_field(3, 2, griddim(1), griddim(2), griddim(3)), &
                   simulation_result%scattering_field(3, 2, griddim(1), griddim(2), griddim(3)), &
                   simulation_result%average_scattering_field(3, 2, griddim(1), griddim(2), griddim(3)))
         simulation_result%average_electric_field = 0.d0
         simulation_result%average_magnetic_field = 0.d0
         simulation_result%average_scattering_field = 0.d0
      end if

      nconfigave = 0
      do configuration_number = 1, ceiling(dble(simulation_config%number_configurations) / &
                                           dble(simulation_config%number_configuration_groups))
         simulation_result%random_configuration_number = configuration_number

         if (rank .eq. 0) then
            if (simulation_result%random_configuration_number .eq. 1) then
               call print_run_variables(sphere_cluster%run_print_unit)
               call open_output_file(simulation_config%output%output_file, file_unit, append=.true.)
               if (runtime_failed()) return
               call print_run_variables(file_unit)
               close (file_unit)
            end if
            write (sphere_cluster%run_print_unit, '('' configuration averaging, samples:'',i5,''-'',i5)') &
               (simulation_result%random_configuration_number - 1) * simulation_config%number_configuration_groups + 1, &
               simulation_result%random_configuration_number * simulation_config%number_configuration_groups
         end if

         if (rank .eq. 0) time1 = parallel_wall_time()

         call execute_simulation(print_output=.false., set_t_matrix_order=.false., mpi_comm=configcomm)
         if (runtime_failed()) return

         if (singleorigin .and. configrank .eq. 0) then
           call common_origin_scattering_cross_section(sphere_cluster%t_matrix_order, simulation_result%incident_coefficients, csca)
         end if

         if (configrank .eq. 0) then
            if (rank .eq. 0) simulation_result%solution_time = parallel_wall_time() - time1
            simulation_result%average_efficiency = simulation_result%average_efficiency + simulation_result%efficiency
        simulation_result%average_total_efficiency = simulation_result%average_total_efficiency + simulation_result%total_efficiency
     simulation_result%average_volume_absorption = simulation_result%average_volume_absorption + simulation_result%volume_absorption
            if (target_shape .eq. 2) then
              call cartesian_vectors_to_spherical(sphere_cluster%number_spheres, sphere_cluster%sphere_position, spherical_position)
               simulation_result%average_sphere_position = simulation_result%average_sphere_position + spherical_position
            else
              simulation_result%average_sphere_position = simulation_result%average_sphere_position + sphere_cluster%sphere_position
            end if
        simulation_result%average_plane_scattering = simulation_result%average_plane_scattering + simulation_result%plane_scattering
            simulation_result%average_surface_absorptance = simulation_result%average_surface_absorptance + simulation_result%surface_absorptance
            if (simulation_config%calculate_up_down_scattering) simulation_result%average_boundary_scattering = simulation_result%average_boundary_scattering + simulation_result%boundary_scattering
            simulation_result%average_boundary_extinction = simulation_result%average_boundary_extinction + simulation_result%boundary_extinction
           if (singleorigin .and. simulation_config%azimuthal_average .and. (.not. simulation_config%numerical_azimuthal_average)) &
               simulation_result%average_scattering_matrix_expansion = simulation_result%average_scattering_matrix_expansion + simulation_result%scattering_matrix_expansion
            if (simulation_config%calculate_scattering_matrix) then
     simulation_result%average_scattering_matrix = simulation_result%average_scattering_matrix + simulation_result%scattering_matrix
            end if
            if (simulation_config%calculate_near_field) then
              simulation_result%average_electric_field = simulation_result%average_electric_field + simulation_result%electric_field
              simulation_result%average_magnetic_field = simulation_result%average_magnetic_field + simulation_result%magnetic_field
               simulation_result%scattering_field(1, :, :, :, :) = 0.5 * (simulation_result%electric_field(2, :, :, :, :) * conjg(simulation_result%magnetic_field(3, :, :, :, :)) &
                         - simulation_result%electric_field(3, :, :, :, :) * conjg(simulation_result%magnetic_field(2, :, :, :, :)))
               simulation_result%scattering_field(2, :, :, :, :) = 0.5 * (-simulation_result%electric_field(1, :, :, :, :) * conjg(simulation_result%magnetic_field(3, :, :, :, :)) &
                         + simulation_result%electric_field(3, :, :, :, :) * conjg(simulation_result%magnetic_field(1, :, :, :, :)))
               simulation_result%scattering_field(3, :, :, :, :) = 0.5 * (simulation_result%electric_field(1, :, :, :, :) * conjg(simulation_result%magnetic_field(2, :, :, :, :)) &
                         - simulation_result%electric_field(2, :, :, :, :) * conjg(simulation_result%magnetic_field(1, :, :, :, :)))
        simulation_result%average_scattering_field = simulation_result%average_scattering_field + simulation_result%scattering_field
            end if
if (rank .eq. 0) simulation_result%average_solution_time = simulation_result%average_solution_time + simulation_result%solution_time
         end if

         if (singleorigin) then
!               if(sphere_1_fixed) call subtract_fixed_sphere_from_common_origin()
!               if(sphere_1_fixed) then
!                  simulation_config%erase_sphere_1=.true.
!                  simulation_config%use_previous_configuration=.true.
!                  sphere_1_fixed=.false.
!                  call execute_simulation(print_output=.false.,set_t_matrix_order=.false.,mpi_comm=configcomm)
!                  simulation_config%erase_sphere_1=.false.
!                  simulation_config%use_previous_configuration=.false.
!                  sphere_1_fixed=.true.
!               endif
            simulation_result%average_incident_coefficients = simulation_result%average_incident_coefficients + simulation_result%incident_coefficients
            simulation_result%average_total_cross_section = simulation_result%average_total_cross_section + csca
         end if

         nsend = 3 * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_sphere_position, &
            receive_buffer=sphere_cluster%sphere_position, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 3 * simulation_result%efficiency_dimension * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_efficiency, &
            receive_buffer=simulation_result%efficiency, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 3 * simulation_result%efficiency_dimension
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_total_efficiency, &
            receive_buffer=simulation_result%total_efficiency, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = simulation_result%efficiency_dimension * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_volume_absorption, &
            receive_buffer=simulation_result%volume_absorption, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 4
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_plane_scattering, &
            receive_buffer=simulation_result%plane_scattering, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 2
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_surface_absorptance, &
            receive_buffer=simulation_result%surface_absorptance, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         if (simulation_config%calculate_up_down_scattering) then
            nsend = 4
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_boundary_scattering, &
               receive_buffer=simulation_result%boundary_scattering, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_boundary_extinction, &
            receive_buffer=simulation_result%boundary_extinction, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         if (singleorigin) then
            nsend = 4 * sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2)
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_incident_coefficients, &
               receive_buffer=simulation_result%incident_coefficients, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend = 1
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_total_cross_section, &
               receive_buffer=csca, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
         if (simulation_config%calculate_scattering_matrix) then
            nsend = simulation_result%scattering_matrix_dimension * (simulation_result%scattering_matrix_upper_bound - simulation_result%scattering_matrix_lower_bound + 1)
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_scattering_matrix, &
               receive_buffer=simulation_result%scattering_matrix, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
        if (singleorigin .and. simulation_config%azimuthal_average .and. (.not. simulation_config%numerical_azimuthal_average)) then
            nsend = 16 * 4 * (2 * sphere_cluster%t_matrix_order + 1)
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_scattering_matrix_expansion, &
               receive_buffer=simulation_result%scattering_matrix_expansion, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
         if (simulation_config%calculate_near_field) then
            nsend = 6 * product(griddim)
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_electric_field, &
               receive_buffer=simulation_result%electric_field, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_magnetic_field, &
               receive_buffer=simulation_result%magnetic_field, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_scattering_field, &
               receive_buffer=simulation_result%scattering_field, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if

         nconfigave = nconfigave + simulation_config%number_configuration_groups
         diffac = (dble(sphere_cluster%number_spheres) - 1.d0) / dble(sphere_cluster%number_spheres)
         diffac = 1.d0
!            diffac=(dble(sphere_cluster%number_spheres))/dble(sphere_cluster%number_spheres-1)
!            diffac=(1.d0-(1.d0/dble(sphere_cluster%number_spheres))**.5d0)

         if (singleorigin) then
            if (rank .eq. 0 .and. simulation_config%output%print_timings) then
               timet = parallel_wall_time()
               write (sphere_cluster%run_print_unit, '('' calculating diffuse field:'')', advance='no')
            end if
            simulation_result%incident_coefficients = simulation_result%incident_coefficients / dble(nconfigave)
!
!  zero out azimuth order .ne. pm 1 for sphere targets 2/23
!
!               if(target_shape.eq.2) then
!                  do n=1,sphere_cluster%t_matrix_order
!                     do m=-n,n
!                        do p=1,2
!                           if(abs(m).ne.1) then
!                              simulation_result%incident_coefficients(polarized_mode_index(m,n,p,sphere_cluster%t_matrix_order,2),:)=0.d0
!                           endif
!                        enddo
!                     enddo
!                  enddo
!               endif
            csca = csca / dble(nconfigave)
            simulation_result%total_cross_section = csca(1)
            call common_origin_scattering_cross_section(sphere_cluster%t_matrix_order, simulation_result%incident_coefficients, simulation_result%diffuse_cross_section_ratio)
  simulation_result%diffuse_cross_section = simulation_result%total_cross_section - simulation_result%diffuse_cross_section_ratio(1)
            simulation_result%diffuse_cross_section_ratio = 1.d0 - simulation_result%diffuse_cross_section_ratio / csca
            if (simulation_config%azimuthal_average .and. (.not. simulation_config%numerical_azimuthal_average)) then
               allocate (texpcoef(16, 0:2 * sphere_cluster%t_matrix_order, 4))
               texpcoef = simulation_result%scattering_matrix_expansion / dble(nconfigave)
!                  simulation_result%scattering_matrix_expansion=simulation_result%scattering_matrix_expansion/dble(nconfigave)
               call fixed_orientation_scattering_matrix_expansion( &
                  sphere_cluster%t_matrix_order, simulation_result%incident_coefficients, simulation_result%coherent_scattering_expansion(:, :, 1), simulation_result%coherent_scattering_expansion(:, :, 2), &
                  simulation_result%coherent_scattering_expansion(:, :, 3), simulation_result%coherent_scattering_expansion(:, :, 4), mpi_comm=configcomm)
               simulation_result%scattering_matrix_expansion = simulation_result%coherent_scattering_expansion
            end if
            if (simulation_config%calculate_scattering_matrix) then
               call compute_scattering_matrix(simulation_result%incident_coefficients, simulation_result%diffuse_scattering_matrix, mpi_comm=configcomm)
            end if
            if (simulation_config%azimuthal_average .and. (.not. simulation_config%numerical_azimuthal_average)) then
!                  simulation_result%scattering_matrix_expansion=texpcoef-simulation_result%scattering_matrix_expansion*diffac
               simulation_result%scattering_matrix_expansion = texpcoef
               deallocate (texpcoef)
            end if
            call hemispherical_scattering(simulation_result%incident_coefficients, .true., simulation_config%numerical_hemispherical_integration, &
                                          simulation_result%diffuse_boundary_scattering, mpi_comm=configcomm)
!               call common_origin_hemispherical_scattering(simulation_result%incident_coefficients,simulation_result%diffuse_boundary_scattering)
      if (rank .eq. 0 .and. simulation_config%output%print_timings) write (sphere_cluster%run_print_unit, '('' completed, '',es12.4,'' sec'')') parallel_wall_time() - timet
            if (rank .eq. 0 .and. target_shape .eq. 2 .and. (.not. simulation_config%random_configuration_host)) then
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
                        aneff = aneff + 0.5d0 * sum(simulation_result%incident_coefficients(mnp, :) / pmnp0(mnp, :))
                     end do
                     aneff = aneff / 2.d0
                     write (file_unit, '(2i4,2es13.5)') n, p, aneff
                     anp0(p, n) = aneff
                  end do
               end do
               close (file_unit)
               deallocate (pmnp0)
!                  sphere_cluster%effective_ref_index=(1.1d0,0.01d0)
!                  simulation_result%fit_radius=target_dimensions(1)*simulation_config%length_scale_factor
               call fit_effective_refractive_index(anp0, sphere_cluster%effective_ref_index, simulation_result%fit_radius, simulation_result%fit_status)
               deallocate (anp0)
            end if
         end if
         if (allocated(simulation_result%incident_coefficients)) deallocate (simulation_result%incident_coefficients)

         if (rank .eq. 0) then
            sphere_cluster%sphere_position = sphere_cluster%sphere_position / dble(nconfigave)
            simulation_result%efficiency = simulation_result%efficiency / dble(nconfigave)
            simulation_result%volume_absorption = simulation_result%volume_absorption / dble(nconfigave)
            simulation_result%total_efficiency = simulation_result%total_efficiency / dble(nconfigave)
            simulation_result%plane_scattering = simulation_result%plane_scattering / dble(nconfigave)
            simulation_result%surface_absorptance = simulation_result%surface_absorptance / dble(nconfigave)
            if (simulation_config%calculate_up_down_scattering) simulation_result%boundary_scattering = simulation_result%boundary_scattering / dble(nconfigave)
            simulation_result%boundary_extinction = simulation_result%boundary_extinction / dble(nconfigave)
!               if(singleorigin) simulation_result%diffuse_boundary_scattering=simulation_result%boundary_scattering-simulation_result%diffuse_boundary_scattering
            if (singleorigin) simulation_result%diffuse_boundary_scattering = simulation_result%boundary_scattering &
                                                                            - simulation_result%diffuse_boundary_scattering * diffac
            if (simulation_config%calculate_scattering_matrix) then
               simulation_result%scattering_matrix = simulation_result%scattering_matrix / dble(nconfigave)
! experiment

!                  if(singleorigin) simulation_result%diffuse_scattering_matrix=simulation_result%scattering_matrix-simulation_result%diffuse_scattering_matrix*(dble(sphere_cluster%number_spheres-1)/dble(sphere_cluster%number_spheres))
!                  if(singleorigin) simulation_result%diffuse_scattering_matrix=simulation_result%scattering_matrix-simulation_result%diffuse_scattering_matrix*diffac

!                  simulation_result%diffuse_scattering_matrix=simulation_result%scattering_matrix-simulation_result%diffuse_scattering_matrix
            end if
     simulation_result%solution_time = simulation_result%average_solution_time / dble(simulation_result%random_configuration_number)
            call print_calculation_results(simulation_config%output%output_file)
            if (simulation_config%calculate_near_field) then
               simulation_result%electric_field = simulation_result%electric_field / dble(nconfigave)
               simulation_result%magnetic_field = simulation_result%magnetic_field / dble(nconfigave)
               simulation_result%scattering_field = simulation_result%scattering_field / dble(nconfigave)
               call open_output_file(simulation_config%output%near_field_file, file_unit)
               if (runtime_failed()) return
               call write_near_field_output_header(griddim, file_unit, print_intersecting_spheres=.false.)
               allocate (edat(griddim(3)))
               edat = 0.d0
               do iz = 1, griddim(3)
                  do iy = 1, griddim(2)
                     do ix = 1, griddim(1)
       edat(iz) = edat(iz) + simulation_result%electric_field(1, 1, ix, iy, iz) + simulation_result%electric_field(2, 2, ix, iy, iz)
                        ipos(:) = (/ix, iy, iz/)
                        rpos(:) = (dble(ipos(:)) - (/0.5d0, 0.5d0, 0.5d0/)) * grid_spacing(:) + grid_region(:, 1)
                        write (file_unit, '(33es12.4)') rpos(:), &
                           simulation_result%electric_field(:, 1, ix, iy, iz), simulation_result%magnetic_field(:, 1, ix, iy, iz), &
                           simulation_result%electric_field(:, 2, ix, iy, iz), simulation_result%magnetic_field(:, 2, ix, iy, iz), &
                          simulation_result%scattering_field(:, 1, ix, iy, iz), simulation_result%scattering_field(:, 2, ix, iy, iz)
                     end do
                  end do
               end do
               edat = edat / dble(griddim(1) * griddim(2) * 2.d0)
               call effective_refractive_index(griddim(3), edat, grid_spacing(3), rieff, e0)
               simulation_result%near_field_effective_ref_index = rieff
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
      integer :: file_unit, rank, numprocs, itemp(1), n, configuration_number, &
                 numprocsperconfig, configcolor, configgroup, configcomm, configrank, config0comm, nconfigave, nsend
      real(real64) :: time1, timet
      real(real64), allocatable :: tpos(:, :)
      complex(real64) :: ctemp(1)
      character(len=256) :: tmatchar1, tmatchar2
      data tmatchar1, tmatchar2/'tmat-', '.tmp'/
      simulation_config%first_run = .false.
      call parallel_rank(mpi_rank=rank)
      call parallel_size(mpi_size=numprocs)
!         if(rank.ne.0) light_up=.false.
      simulation_result%local_rank = rank
      global_rank = rank

      if (simulation_config%solver%max_iterations .le. 0) then
         numprocsperconfig = 2
      else
         numprocsperconfig = 4
      end if
      simulation_config%number_configuration_groups = numprocs / numprocsperconfig
      simulation_config%number_configuration_groups = max(simulation_config%number_configuration_groups, 1)
      configcolor = floor(dble(simulation_config%number_configuration_groups * rank) / dble(numprocs))
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
      simulation_config%random_configuration = .true.
!singleorigin=.true.

      call execute_simulation(print_output=.false., set_t_matrix_order=.true., dry_run=.true.)
      if (runtime_failed()) return
      call compose_group_filename(tmatchar1, configgroup, tmatchar2, simulation_config%output%t_matrix_file)

      if (allocated(simulation_result%average_efficiency)) deallocate (simulation_result%average_efficiency, simulation_result%average_total_efficiency, simulation_result%average_volume_absorption, simulation_result%average_sphere_position)
      if (allocated(simulation_result%average_mean_t_matrix)) deallocate (simulation_result%average_mean_t_matrix)
      if (allocated(simulation_result%average_scattering_matrix_expansion)) deallocate (simulation_result%average_scattering_matrix_expansion)
      if (allocated(simulation_result%average_coherent_scattering_expansion)) deallocate (simulation_result%average_coherent_scattering_expansion)
      allocate (simulation_result%average_efficiency(3, simulation_result%efficiency_dimension, sphere_cluster%number_spheres), simulation_result%average_total_efficiency(3, simulation_result%efficiency_dimension), simulation_result%average_volume_absorption(simulation_result%efficiency_dimension, sphere_cluster%number_spheres), &
                simulation_result%average_sphere_position(3, sphere_cluster%number_spheres), simulation_result%average_mean_t_matrix(2, sphere_cluster%t_matrix_order), &
                simulation_result%average_scattering_matrix_expansion(4, 4, 0:2 * sphere_cluster%t_matrix_order), simulation_result%average_coherent_scattering_expansion(4, 4, 0:2 * sphere_cluster%t_matrix_order))
      simulation_result%average_efficiency = 0.d0
      simulation_result%average_total_efficiency = 0.d0
      simulation_result%average_volume_absorption = 0.d0
      simulation_result%average_sphere_position = 0.d0
      simulation_result%average_mean_t_matrix = 0.d0
      simulation_result%average_scattering_matrix_expansion = 0.d0
      simulation_result%average_coherent_scattering_expansion = 0.d0
      if (simulation_config%calculate_scattering_matrix) then
         if (allocated(simulation_result%average_scattering_matrix)) deallocate (simulation_result%average_scattering_matrix)
         allocate (simulation_result%average_scattering_matrix(simulation_result%scattering_matrix_dimension, simulation_result%scattering_matrix_lower_bound:simulation_result%scattering_matrix_upper_bound))
         simulation_result%average_scattering_matrix = 0.d0
      end if

      nconfigave = 0
      do configuration_number = 1, ceiling(dble(simulation_config%number_configurations) / &
                                           dble(simulation_config%number_configuration_groups))
         simulation_result%random_configuration_number = configuration_number

         if (rank .eq. 0) then
            if (simulation_result%random_configuration_number .eq. 1) then
               call print_run_variables(sphere_cluster%run_print_unit)
               call open_output_file(simulation_config%output%output_file, file_unit, append=.true.)
               if (runtime_failed()) return
               call print_run_variables(file_unit)
               close (file_unit)
            end if
            write (sphere_cluster%run_print_unit, '('' configuration averaging, samples:'',i5,''-'',i5)') &
               (simulation_result%random_configuration_number - 1) * simulation_config%number_configuration_groups + 1, &
               simulation_result%random_configuration_number * simulation_config%number_configuration_groups
         end if

         if (rank .eq. 0) time1 = parallel_wall_time()

         call execute_simulation(print_output=.false., set_t_matrix_order=.false., mpi_comm=configcomm)
         if (runtime_failed()) return

         if (configrank .eq. 0) then
            if (rank .eq. 0) simulation_result%solution_time = parallel_wall_time() - time1
            allocate (tpos(3, sphere_cluster%number_spheres))
            call cartesian_vectors_to_spherical(sphere_cluster%number_spheres, sphere_cluster%sphere_position(:, 1:sphere_cluster%number_spheres), &
                                                tpos(:, 1:sphere_cluster%number_spheres))
            simulation_result%average_efficiency = simulation_result%average_efficiency + simulation_result%efficiency
        simulation_result%average_total_efficiency = simulation_result%average_total_efficiency + simulation_result%total_efficiency
     simulation_result%average_volume_absorption = simulation_result%average_volume_absorption + simulation_result%volume_absorption
            simulation_result%average_sphere_position = simulation_result%average_sphere_position + tpos
            simulation_result%average_mean_t_matrix(:, 1:sphere_cluster%t_matrix_order) = simulation_result%average_mean_t_matrix(:, 1:sphere_cluster%t_matrix_order) + simulation_result%mean_t_matrix(:, 1:sphere_cluster%t_matrix_order)
            if (simulation_config%calculate_scattering_matrix) then
     simulation_result%average_scattering_matrix = simulation_result%average_scattering_matrix + simulation_result%scattering_matrix
               simulation_result%average_scattering_matrix_expansion = simulation_result%average_scattering_matrix_expansion + simulation_result%scattering_matrix_expansion
               simulation_result%average_coherent_scattering_expansion = simulation_result%average_coherent_scattering_expansion + simulation_result%coherent_scattering_expansion
            end if
if (rank .eq. 0) simulation_result%average_solution_time = simulation_result%average_solution_time + simulation_result%solution_time
            deallocate (tpos)
         end if

         nsend = 3 * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_sphere_position, &
            receive_buffer=sphere_cluster%sphere_position, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 3 * simulation_result%efficiency_dimension * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_efficiency, &
            receive_buffer=simulation_result%efficiency, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 3 * simulation_result%efficiency_dimension
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_total_efficiency, &
            receive_buffer=simulation_result%total_efficiency, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = simulation_result%efficiency_dimension * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_volume_absorption, &
            receive_buffer=simulation_result%volume_absorption, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 2 * sphere_cluster%t_matrix_order
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_mean_t_matrix, &
            receive_buffer=simulation_result%mean_t_matrix, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         if (simulation_config%calculate_scattering_matrix) then
            nsend = simulation_result%scattering_matrix_dimension * (simulation_result%scattering_matrix_upper_bound - simulation_result%scattering_matrix_lower_bound + 1)
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_scattering_matrix, &
               receive_buffer=simulation_result%scattering_matrix, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend = 16 * (2 * sphere_cluster%t_matrix_order + 1)
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_scattering_matrix_expansion, &
               receive_buffer=simulation_result%scattering_matrix_expansion, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_coherent_scattering_expansion, &
               receive_buffer=simulation_result%coherent_scattering_expansion, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if

         nconfigave = nconfigave + simulation_config%number_configuration_groups

         if (rank .eq. 0) then
            sphere_cluster%sphere_position = sphere_cluster%sphere_position / dble(nconfigave)
            simulation_result%efficiency = simulation_result%efficiency / dble(nconfigave)
            simulation_result%volume_absorption = simulation_result%volume_absorption / dble(nconfigave)
            simulation_result%total_efficiency = simulation_result%total_efficiency / dble(nconfigave)
            simulation_result%mean_t_matrix(:, 1:sphere_cluster%t_matrix_order) = simulation_result%mean_t_matrix(:, 1:sphere_cluster%t_matrix_order) / dble(nconfigave)
            if (simulation_config%calculate_scattering_matrix) then
               simulation_result%scattering_matrix = simulation_result%scattering_matrix / dble(nconfigave)
               simulation_result%scattering_matrix_expansion = simulation_result%scattering_matrix_expansion / dble(nconfigave)
               simulation_result%coherent_scattering_expansion = simulation_result%coherent_scattering_expansion / dble(nconfigave)
            end if
     simulation_result%solution_time = simulation_result%average_solution_time / dble(simulation_result%random_configuration_number)
            call fit_effective_refractive_index(simulation_result%mean_t_matrix, sphere_cluster%effective_ref_index, simulation_result%fit_radius, simulation_result%fit_status)
            call print_calculation_results(simulation_config%output%output_file)
         end if
!            if(simulation_config%random_configuration_host) then
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
      integer :: file_unit, rank, numprocs, direction_number, &
                 numprocsperconfig, configcolor, configgroup, configcomm, configrank, config0comm, nconfigave, nsend
      real(real64) :: time1, timet
      real(real64), allocatable :: texpcoef(:, :, :)
      character(len=256) :: sdatfile
      simulation_config%first_run = .false.
      call parallel_rank(mpi_rank=rank)
      call parallel_size(mpi_size=numprocs)
!         if(rank.ne.0) light_up=.false.
      simulation_result%local_rank = rank
      global_rank = rank
      sdatfile = simulation_config%output%sphere_data_file
      prancon = simulation_config%output%print_random_configuration
      simulation_config%output%print_random_configuration = .true.
      aa = simulation_config%azimuthal_average
      soe = simulation_config%single_origin_expansion
      iframe = simulation_config%incident_frame
      cuds = simulation_config%calculate_up_down_scattering
      simulation_config%azimuthal_average = .true.
      simulation_config%single_origin_expansion = .true.
      simulation_config%incident_frame = .true.
      simulation_config%calculate_up_down_scattering = .false.

      if (simulation_config%solver%max_iterations .lt. 0) then
         numprocsperconfig = 2
      else
         numprocsperconfig = 4
      end if
      simulation_config%number_configuration_groups = numprocs / numprocsperconfig
      simulation_config%number_configuration_groups = max(simulation_config%number_configuration_groups, 1)
      configcolor = floor(dble(simulation_config%number_configuration_groups * rank) / dble(numprocs))
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
      singleorigin = number_plane_boundaries .eq. 0 .and. simulation_config%single_origin_expansion
      simulation_config%incident_beta_specified = .true.

      call execute_simulation(print_output=.false., set_t_matrix_order=.true., dry_run=.true.)
      if (runtime_failed()) return

      if (trim(simulation_config%output%sphere_data_file) .eq. 'random_configuration') then
         simulation_config%output%sphere_data_file = 'random_configuration.pos'
      end if

      if (allocated(simulation_result%average_efficiency)) deallocate (simulation_result%average_efficiency, simulation_result%average_total_efficiency, simulation_result%average_volume_absorption, simulation_result%average_boundary_scattering, &
                                       simulation_result%average_boundary_extinction, simulation_result%diffuse_boundary_scattering)
      allocate (simulation_result%average_efficiency(3, simulation_result%efficiency_dimension, sphere_cluster%number_spheres), simulation_result%average_total_efficiency(3, simulation_result%efficiency_dimension), simulation_result%average_volume_absorption(simulation_result%efficiency_dimension, sphere_cluster%number_spheres), &
                simulation_result%average_boundary_scattering(2, 0:1), simulation_result%average_boundary_extinction(2, 0:1), simulation_result%diffuse_boundary_scattering(2, 0:1))
      simulation_result%average_efficiency = 0.d0
      simulation_result%average_total_efficiency = 0.d0
      simulation_result%average_volume_absorption = 0.d0
      simulation_result%average_plane_scattering = 0.d0
      simulation_result%average_boundary_scattering = 0.d0
      simulation_result%average_boundary_extinction = 0.d0
      simulation_result%average_solution_time = 0.d0
      if (singleorigin) then
         if (allocated(simulation_result%average_incident_coefficients)) deallocate (simulation_result%average_incident_coefficients, simulation_result%average_scattering_matrix_expansion)
         allocate (simulation_result%average_incident_coefficients(2 * sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2), 2), &
                   simulation_result%average_scattering_matrix_expansion(16, 0:2 * sphere_cluster%t_matrix_order, 4))
         simulation_result%average_incident_coefficients = 0.d0
         simulation_result%average_scattering_matrix_expansion = 0.d0
      end if
      if (simulation_config%calculate_scattering_matrix) then
         if (allocated(simulation_result%average_scattering_matrix)) deallocate (simulation_result%average_scattering_matrix, simulation_result%diffuse_scattering_matrix)
         allocate (simulation_result%average_scattering_matrix(simulation_result%scattering_matrix_dimension, simulation_result%scattering_matrix_lower_bound:simulation_result%scattering_matrix_upper_bound), &
                   simulation_result%diffuse_scattering_matrix(simulation_result%scattering_matrix_dimension, simulation_result%scattering_matrix_lower_bound:simulation_result%scattering_matrix_upper_bound))
         simulation_result%average_scattering_matrix = 0.d0
      end if

      nconfigave = 0
      do direction_number = 1, ceiling(dble(simulation_config%number_incident_directions) / &
                                       dble(simulation_config%number_configuration_groups))
         simulation_result%incident_direction_number = direction_number

         if (rank .eq. 0) then
            if (simulation_result%incident_direction_number .eq. 1) then
               call print_run_variables(sphere_cluster%run_print_unit)
               call open_output_file(simulation_config%output%output_file, file_unit, append=.true.)
               if (runtime_failed()) return
               call print_run_variables(file_unit)
               close (file_unit)
            end if
            write (sphere_cluster%run_print_unit, '('' incidence averaging, samples:'',i5,''-'',i5)') &
               (simulation_result%incident_direction_number - 1) * simulation_config%number_configuration_groups + 1, &
               simulation_result%incident_direction_number * simulation_config%number_configuration_groups
         end if

         call sample_incident_direction(mpi_comm=configcomm)

         if (rank .eq. 0) time1 = parallel_wall_time()

         call execute_simulation(print_output=.false., set_t_matrix_order=.false., mpi_comm=configcomm)
         if (runtime_failed()) return

         if (singleorigin) then
            simulation_result%average_incident_coefficients = simulation_result%average_incident_coefficients + simulation_result%incident_coefficients
         end if

         if (configrank .eq. 0) then
            if (rank .eq. 0) simulation_result%solution_time = parallel_wall_time() - time1
            simulation_result%average_efficiency = simulation_result%average_efficiency + simulation_result%efficiency
        simulation_result%average_total_efficiency = simulation_result%average_total_efficiency + simulation_result%total_efficiency
     simulation_result%average_volume_absorption = simulation_result%average_volume_absorption + simulation_result%volume_absorption
        simulation_result%average_plane_scattering = simulation_result%average_plane_scattering + simulation_result%plane_scattering
            if (simulation_config%calculate_up_down_scattering) simulation_result%average_boundary_scattering = simulation_result%average_boundary_scattering + simulation_result%boundary_scattering
            simulation_result%average_boundary_extinction = simulation_result%average_boundary_extinction + simulation_result%boundary_extinction
           if (singleorigin .and. simulation_config%azimuthal_average .and. (.not. simulation_config%numerical_azimuthal_average)) &
               simulation_result%average_scattering_matrix_expansion = simulation_result%average_scattering_matrix_expansion + simulation_result%scattering_matrix_expansion
            if (simulation_config%calculate_scattering_matrix) then
     simulation_result%average_scattering_matrix = simulation_result%average_scattering_matrix + simulation_result%scattering_matrix
            end if
if (rank .eq. 0) simulation_result%average_solution_time = simulation_result%average_solution_time + simulation_result%solution_time
         end if

         nsend = 3 * simulation_result%efficiency_dimension * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_efficiency, &
            receive_buffer=simulation_result%efficiency, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 3 * simulation_result%efficiency_dimension
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_total_efficiency, &
            receive_buffer=simulation_result%total_efficiency, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = simulation_result%efficiency_dimension * sphere_cluster%number_spheres
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_volume_absorption, &
            receive_buffer=simulation_result%volume_absorption, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         nsend = 4
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_plane_scattering, &
            receive_buffer=simulation_result%plane_scattering, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         if (simulation_config%calculate_up_down_scattering) then
            nsend = 4
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_boundary_scattering, &
               receive_buffer=simulation_result%boundary_scattering, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
         call parallel_reduce_sum( &
            send_buffer=simulation_result%average_boundary_extinction, &
            receive_buffer=simulation_result%boundary_extinction, &
            mpi_rank=0, &
            mpi_number=nsend, &
            mpi_comm=config0comm)
         if (singleorigin) then
            nsend = 4 * sphere_cluster%t_matrix_order * (sphere_cluster%t_matrix_order + 2)
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_incident_coefficients, &
               receive_buffer=simulation_result%incident_coefficients, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
         if (simulation_config%calculate_scattering_matrix) then
            nsend = simulation_result%scattering_matrix_dimension * (simulation_result%scattering_matrix_upper_bound - simulation_result%scattering_matrix_lower_bound + 1)
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_scattering_matrix, &
               receive_buffer=simulation_result%scattering_matrix, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if
        if (singleorigin .and. simulation_config%azimuthal_average .and. (.not. simulation_config%numerical_azimuthal_average)) then
            nsend = 16 * 4 * (2 * sphere_cluster%t_matrix_order + 1)
            call parallel_reduce_sum( &
               send_buffer=simulation_result%average_scattering_matrix_expansion, &
               receive_buffer=simulation_result%scattering_matrix_expansion, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
         end if

         nconfigave = nconfigave + simulation_config%number_configuration_groups

         if (singleorigin) then
            if (rank .eq. 0 .and. simulation_config%output%print_timings) then
               timet = parallel_wall_time()
               write (sphere_cluster%run_print_unit, '('' calculating diffuse field:'')', advance='no')
            end if
            simulation_result%incident_coefficients = simulation_result%incident_coefficients / dble(nconfigave)
            if (simulation_config%azimuthal_average .and. (.not. simulation_config%numerical_azimuthal_average)) then
               allocate (texpcoef(16, 0:2 * sphere_cluster%t_matrix_order, 4))
               texpcoef = simulation_result%scattering_matrix_expansion / dble(nconfigave)
               call fixed_orientation_scattering_matrix_expansion( &
                  sphere_cluster%t_matrix_order, simulation_result%incident_coefficients, simulation_result%scattering_matrix_expansion(:, :, 1), simulation_result%scattering_matrix_expansion(:, :, 2), &
simulation_result%scattering_matrix_expansion(:, :, 3), simulation_result%scattering_matrix_expansion(:, :, 4), mpi_comm=configcomm)
            end if
            if (simulation_config%calculate_scattering_matrix) then
               call compute_scattering_matrix(simulation_result%incident_coefficients, simulation_result%diffuse_scattering_matrix, mpi_comm=configcomm)
            end if
            if (simulation_config%azimuthal_average .and. (.not. simulation_config%numerical_azimuthal_average)) then
               simulation_result%scattering_matrix_expansion = texpcoef - simulation_result%scattering_matrix_expansion * (dble(sphere_cluster%number_spheres - 1) / dble(sphere_cluster%number_spheres))**2
               deallocate (texpcoef)
            end if
 call common_origin_hemispherical_scattering(simulation_result%incident_coefficients, simulation_result%diffuse_boundary_scattering)
      if (rank .eq. 0 .and. simulation_config%output%print_timings) write (sphere_cluster%run_print_unit, '('' completed, '',es12.4,'' sec'')') parallel_wall_time() - timet
         end if
         if (allocated(simulation_result%incident_coefficients)) deallocate (simulation_result%incident_coefficients)

         if (rank .eq. 0) then
            simulation_result%efficiency = simulation_result%efficiency / dble(nconfigave)
            simulation_result%volume_absorption = simulation_result%volume_absorption / dble(nconfigave)
            simulation_result%total_efficiency = simulation_result%total_efficiency / dble(nconfigave)
            simulation_result%plane_scattering = simulation_result%plane_scattering / dble(nconfigave)
            if (simulation_config%calculate_up_down_scattering) simulation_result%boundary_scattering = simulation_result%boundary_scattering / dble(nconfigave)
            simulation_result%boundary_extinction = simulation_result%boundary_extinction / dble(nconfigave)
!               if(singleorigin) simulation_result%diffuse_boundary_scattering=simulation_result%boundary_scattering-simulation_result%diffuse_boundary_scattering
            if (singleorigin) simulation_result%diffuse_boundary_scattering = simulation_result%boundary_scattering &
     - simulation_result%diffuse_boundary_scattering * dble(sphere_cluster%number_spheres - 1) / dble(sphere_cluster%number_spheres)
            if (simulation_config%calculate_scattering_matrix) then
               simulation_result%scattering_matrix = simulation_result%scattering_matrix / dble(nconfigave)
! experiment
               if (singleorigin) simulation_result%diffuse_scattering_matrix = simulation_result%scattering_matrix - simulation_result%diffuse_scattering_matrix * (dble(sphere_cluster%number_spheres - 1) / dble(sphere_cluster%number_spheres))**2

!                  simulation_result%diffuse_scattering_matrix=simulation_result%scattering_matrix-simulation_result%diffuse_scattering_matrix
            end if
     simulation_result%solution_time = simulation_result%average_solution_time / dble(simulation_result%random_configuration_number)
            call print_calculation_results(simulation_config%output%output_file)
         end if
      end do
      simulation_config%output%sphere_data_file = sdatfile
      simulation_config%output%print_random_configuration = prancon
      simulation_config%azimuthal_average = aa
      simulation_config%single_origin_expansion = soe
      simulation_config%incident_frame = iframe
      simulation_config%calculate_up_down_scattering = cuds
   end subroutine run_incidence_average

   subroutine common_origin_scattering_cross_section(n, a, c)
      implicit none
      integer :: n
      real(real64) :: c(1)
      complex(real64) :: a(4 * n * (n + 2))
      c(1) = sum(a(:) * conjg(a(:)))
   end subroutine common_origin_scattering_cross_section

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
      simulation_config%incident_beta_degrees = sbuf(1)
      simulation_config%incident_alpha_degrees = sbuf(2)
   end subroutine sample_incident_direction

end module simulation_averaging
