module random_orientation_averaging
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
   public :: run_random_orientation_configuration_average
contains

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
end module random_orientation_averaging
