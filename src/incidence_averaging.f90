module incidence_averaging
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
   public :: run_incidence_average
contains

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

end module incidence_averaging
