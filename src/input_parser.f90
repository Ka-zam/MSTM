module input_parser
   use, intrinsic :: iso_fortran_env, only: real64
   use input_state
   use input_value_parsing
   use parallel_runtime, only: parallel_rank
   use runtime_support, only: open_output_file, runtime_failed, set_runtime_error
   implicit none
contains

   subroutine process_input_variable(varlabel, &
                                     var_value, var_type, &
                                     var_position, var_operation, var_status, &
                                     i_var_pointer, r_var_pointer, c_var_pointer)
      implicit none
      logical :: operate
      logical, pointer :: lvarvalue, lavarvalue(:)
      integer :: varpos, varstatus, varlen
      integer, optional :: var_position, var_status
      integer, pointer :: ivarvalue
      integer, optional, pointer :: i_var_pointer
      real(real64), pointer :: rvarvalue, ravarvalue(:)
      real(real64), optional, pointer :: r_var_pointer
      complex(real64), pointer :: cvarvalue, cavarvalue(:)
      complex(real64), optional, pointer :: c_var_pointer
      character(len=1) :: vartype
      character(len=1), optional :: var_type
      character(len=*), optional :: var_value, var_operation
      character(len=256) :: varop, sentvarvalue, varlabel
      character(len=256), pointer :: avarvalue

      if (present(var_operation)) then
         varop = trim(var_operation)
      else
         varop = ' '
      end if
      if (present(var_value)) then
         sentvarvalue = trim(var_value)
         operate = .true.
      else
         sentvarvalue = ' '
         operate = .false.
      end if
      if (present(var_position)) then
         varpos = var_position
      else
         varpos = 1
      end if
      varstatus = 0
      vartype = 'n'
      varlen = 1

      if (varlabel .eq. 'output_file') then
         vartype = 'a'
         avarvalue => simulation_config%output%output_file

      elseif (varlabel .eq. 'append_output_file') then
         vartype = 'l'
         lvarvalue => simulation_config%output%append

      elseif (varlabel .eq. 'copy_input_file') then
         vartype = 'l'
         lvarvalue => simulation_config%output%copy_input

      elseif (varlabel .eq. 'run_file') then
         vartype = 'a'
         avarvalue => simulation_config%output%run_file

      elseif (varlabel .eq. 'sphere_data_input_file') then
         vartype = 'a'
         avarvalue => simulation_config%output%sphere_data_file
         if (operate) then
            simulation_config%embedded_sphere_data = .false.
            simulation_config%number_sphere_data_records = 0
            simulation_config%sphere_data_source = trim(sentvarvalue)
            if (allocated(simulation_config%sphere_data_records)) &
               deallocate (simulation_config%sphere_data_records, simulation_config%sphere_data_record_lines)
         end if
         sphere_cluster%recalculate_surface_matrix = .true.

      elseif (varlabel .eq. 'max_iterations') then
         vartype = 'i'
         ivarvalue => simulation_config%solver%max_iterations

      elseif (varlabel .eq. 'solution_epsilon') then
         vartype = 'r'
         rvarvalue => simulation_config%solver%solution_epsilon

      elseif (varlabel .eq. 'normalize_solution_error') then
         vartype = 'l'
         lvarvalue => sphere_cluster%normalize_solution_error

      elseif (varlabel .eq. 'mie_epsilon') then
         vartype = 'r'
         rvarvalue => simulation_config%solver%mie_epsilon

      elseif (varlabel .eq. 'translation_epsilon') then
         vartype = 'r'
         rvarvalue => simulation_config%solver%translation_epsilon

      elseif (varlabel .eq. 'random_orientation') then
         vartype = 'l'
         lvarvalue => simulation_config%random_orientation

      elseif (varlabel .eq. 't_matrix_centered_on_1') then
         vartype = 'l'
         lvarvalue => simulation_config%t_matrix_centered_on_1

      elseif (varlabel .eq. 't_matrix_convergence_epsilon') then
         vartype = 'r'
         rvarvalue => simulation_config%solver%t_matrix_convergence_epsilon

      elseif (varlabel .eq. 'solution_method') then
         vartype = 'a'
         avarvalue => simulation_config%solver%solution_method

      elseif (varlabel .eq. 't_matrix_procs_per_solution') then
         vartype = 'i'
         ivarvalue => simulation_config%solver%t_matrix_procs_per_solution

      elseif (varlabel .eq. 'max_t_matrix_order') then
         vartype = 'i'
         ivarvalue => sphere_cluster%max_t_matrix_order

      elseif (varlabel .eq. 'fft_translation_option') then
         vartype = 'l'
         lvarvalue => simulation_config%fft%enabled

      elseif (varlabel .eq. 'node_order') then
         vartype = 'i'
         ivarvalue => simulation_config%fft%node_order

      elseif (varlabel .eq. 'min_fft_nsphere') then
         vartype = 'i'
         ivarvalue => simulation_config%fft%minimum_spheres

      elseif (varlabel .eq. 'neighbor_node_model') then
         vartype = 'i'
         ivarvalue => simulation_config%fft%neighbor_node_model

      elseif (varlabel .eq. 'cell_volume_fraction') then
         vartype = 'r'
         rvarvalue => simulation_config%fft%cell_volume_fraction
         simulation_config%fft%cell_size_specified = .false.

      elseif (varlabel .eq. 'd_cell') then
         vartype = 'r'
         rvarvalue => simulation_config%fft%cell_size
         simulation_config%fft%cell_size_specified = .true.

      elseif (varlabel .eq. 'incident_beta_deg') then
         vartype = 'r'
         rvarvalue => simulation_config%incident_beta_degrees
         simulation_config%incident_beta_specified = .true.

      elseif (varlabel .eq. 'incident_sin_beta') then
         vartype = 'r'
         rvarvalue => simulation_result%incident_sin_beta
         simulation_config%incident_beta_specified = .false.

      elseif (varlabel .eq. 'incident_direction') then
         vartype = 'i'
         ivarvalue => simulation_config%incident_direction

      elseif (varlabel .eq. 'incident_alpha_deg') then
         vartype = 'r'
         rvarvalue => simulation_config%incident_alpha_degrees

      elseif (varlabel .eq. 'excitation_type') then
         vartype = 'a'
         avarvalue => simulation_config%excitation_type

      elseif (varlabel .eq. 'electric_dipole_position') then
         vartype = 'r'
         varlen = 3
         ravarvalue => simulation_config%electric_dipole_position

      elseif (varlabel .eq. 'electric_dipole_moment') then
         vartype = 'c'
         varlen = 3
         cavarvalue => simulation_config%electric_dipole_moment

      elseif (varlabel .eq. 'gaussian_beam_constant') then
         vartype = 'r'
         rvarvalue => sphere_cluster%gaussian_beam_constant

      elseif (varlabel .eq. 'excitation_radius') then
         vartype = 'r'
         rvarvalue => simulation_config%excitation_radius

      elseif (varlabel .eq. 'interaction_radius') then
         vartype = 'r'
         rvarvalue => interaction_radius

      elseif (varlabel .eq. 'incidence_average') then
         vartype = 'l'
         lvarvalue => simulation_config%incidence_average

      elseif (varlabel .eq. 'number_incident_directions') then
         vartype = 'i'
         ivarvalue => simulation_config%number_incident_directions

      elseif (varlabel .eq. 'gaussian_beam_focal_point') then
         vartype = 'r'
         varlen = 3
         ravarvalue => sphere_cluster%gaussian_beam_focal_point(1:3)

      elseif (varlabel .eq. 'calculate_scattering_matrix') then
         vartype = 'l'
         lvarvalue => simulation_config%calculate_scattering_matrix

      elseif (varlabel .eq. 'single_origin_expansion') then
         vartype = 'l'
         lvarvalue => simulation_config%single_origin_expansion

      elseif (varlabel .eq. 'scattering_map_model') then
         vartype = 'i'
         ivarvalue => simulation_config%scattering_map_model

      elseif (varlabel .eq. 'scattering_map_dimension') then
         vartype = 'i'
         ivarvalue => simulation_config%scattering_map_dimension

      elseif (varlabel .eq. 'scattering_map_increment') then
         vartype = 'r'
         rvarvalue => simulation_config%scattering_map_increment

      elseif (varlabel .eq. 'azimuthal_average') then
         vartype = 'l'
         lvarvalue => simulation_config%azimuthal_average

      elseif (varlabel .eq. 'numerical_azimuthal_average') then
         vartype = 'l'
         lvarvalue => simulation_config%numerical_azimuthal_average

      elseif (varlabel .eq. 'incident_frame') then
         vartype = 'l'
         lvarvalue => simulation_config%incident_frame

      elseif (varlabel .eq. 'print_sphere_data') then
         vartype = 'l'
         lvarvalue => simulation_config%output%print_sphere_data

      elseif (varlabel .eq. 'number_spheres') then
         vartype = 'i'
         ivarvalue => simulation_config%input_number_spheres
         sphere_cluster%recalculate_surface_matrix = .true.
         simulation_config%number_spheres_specified = .true.

      elseif (varlabel .eq. 'length_scale_factor') then
         vartype = 'r'
         rvarvalue => simulation_config%length_scale_factor
         sphere_cluster%recalculate_surface_matrix = .true.

      elseif (varlabel .eq. 'ref_index_scale_factor') then
         vartype = 'c'
         cvarvalue => simulation_config%ref_index_scale_factor
         sphere_cluster%recalculate_surface_matrix = .true.

      elseif (varlabel .eq. 'number_plane_boundaries') then
         vartype = 'i'
         ivarvalue => number_plane_boundaries
         sphere_cluster%recalculate_surface_matrix = .true.

      elseif (varlabel .eq. 'maximum_integration_subdivisions') then
         vartype = 'i'
         ivarvalue => maximum_integration_subdivisions
         sphere_cluster%recalculate_surface_matrix = .true.

      elseif (varlabel .eq. 'integration_error_epsilon') then
         vartype = 'r'
         rvarvalue => integration_error_epsilon
         sphere_cluster%recalculate_surface_matrix = .true.

      elseif (varlabel .eq. 'integration_limit_epsilon') then
         vartype = 'r'
         rvarvalue => integration_limit_epsilon
         sphere_cluster%recalculate_surface_matrix = .true.

      elseif (varlabel .eq. 'minimum_initial_segment_size') then
         vartype = 'r'
         rvarvalue => minimum_initial_segment_size
         sphere_cluster%recalculate_surface_matrix = .true.

      elseif (varlabel .eq. 'gf_switch_factor') then
         vartype = 'r'
         rvarvalue => gf_switch_factor
         sphere_cluster%recalculate_surface_matrix = .true.

      elseif (varlabel .eq. 's_scale_constant') then
         vartype = 'r'
         rvarvalue => s_scale_constant
         sphere_cluster%recalculate_surface_matrix = .true.

      elseif (varlabel .eq. 'real_axis_integration_limit') then
         vartype = 'r'
         rvarvalue => real_axis_integration_limit
         sphere_cluster%recalculate_surface_matrix = .true.

      elseif (varlabel .eq. 'minimum_integration_spacing') then
         vartype = 'r'
         rvarvalue => minimum_integration_spacing
         sphere_cluster%recalculate_surface_matrix = .true.

      elseif (varlabel .eq. 'move_to_front') then
         vartype = 'l'
         lvarvalue => simulation_config%move_to_front

      elseif (varlabel .eq. 'move_to_back') then
         vartype = 'l'
         lvarvalue => simulation_config%move_to_back

      elseif (varlabel .eq. 'store_translation_matrix') then
         vartype = 'l'
         lvarvalue => sphere_cluster%store_translation_matrix

      elseif (varlabel .eq. 'store_surface_matrix') then
         vartype = 'l'
         lvarvalue => sphere_cluster%store_surface_matrix

      elseif (varlabel .eq. 'calculate_near_field') then
         vartype = 'l'
         lvarvalue => simulation_config%calculate_near_field

      elseif (varlabel .eq. 'store_surface_vector') then
         vartype = 'l'
         lvarvalue => store_surface_vector

      elseif (varlabel .eq. 'fast_near_field') then
         vartype = 'l'
         lvarvalue => fast_near_field

      elseif (varlabel .eq. 'near_field_output_file') then
         vartype = 'a'
         avarvalue => simulation_config%output%near_field_file
         simulation_config%output%append_near_field = .false.

      elseif (varlabel .eq. 'near_field_calculation_model') then
         vartype = 'i'
         ivarvalue => simulation_config%near_field_calculation_model

      elseif (varlabel .eq. 'near_field_expansion_order') then
         vartype = 'i'
         ivarvalue => near_field_expansion_order

      elseif (varlabel .eq. 'near_field_expansion_spacing') then
         vartype = 'r'
         rvarvalue => near_field_expansion_spacing

      elseif (varlabel .eq. 'near_field_step_size') then
         vartype = 'r'
         rvarvalue => simulation_config%near_field_step_size

      elseif (varlabel .eq. 'near_field_minimum_border') then
         vartype = 'r'
         varlen = 3
         ravarvalue => simulation_config%near_field_plane_vertices(1:3, 1)

      elseif (varlabel .eq. 'near_field_maximum_border') then
         vartype = 'r'
         varlen = 3
         ravarvalue => simulation_config%near_field_plane_vertices(1:3, 2)

      elseif (varlabel .eq. 'normalize_s11') then
         vartype = 'l'
         lvarvalue => simulation_config%output%normalize_s11

      elseif (varlabel .eq. 'periodic_lattice') then
         vartype = 'l'
         lvarvalue => periodic_lattice

      elseif (varlabel .eq. 'phase_shift_form') then
         vartype = 'l'
         lvarvalue => phase_shift_form

      elseif (varlabel .eq. 'finite_lattice') then
         vartype = 'l'
         lvarvalue => finite_lattice

      elseif (varlabel .eq. 'cell_width') then
         vartype = 'r'
         varlen = 2
         ravarvalue => simulation_config%input_cell_width(1:2)
         simulation_config%square_cell = .false.

      elseif (varlabel .eq. 'cell_width_x') then
         vartype = 'r'
         rvarvalue => simulation_config%input_cell_width_x
         simulation_config%square_cell = .true.

      elseif (varlabel .eq. 'random_configuration') then
         vartype = 'l'
         lvarvalue => simulation_config%random_configuration

      elseif (varlabel .eq. 'random_configuration_output_file') then
         vartype = 'a'
         avarvalue => simulation_config%output%random_configuration_file

      elseif (varlabel .eq. 'target_dimensions') then
         vartype = 'r'
         varlen = 3
         ravarvalue => target_dimensions(1:3)
         target_width_specified = .false.

      elseif (varlabel .eq. 'target_shape') then
         vartype = 'i'
         ivarvalue => target_shape

      elseif (varlabel .eq. 'target_width') then
         vartype = 'r'
         rvarvalue => target_width
         target_width_specified = .true.

      elseif (varlabel .eq. 'target_thickness') then
         vartype = 'r'
         rvarvalue => target_thickness
         target_width_specified = .true.
!
!         elseif(varlabel.eq.'psd_sigma') then
!            vartype='r'
!            rvarvalue=>psd_sigma

      elseif (varlabel .eq. 'number_components') then
         vartype = 'i'
         ivarvalue => number_components

      elseif (varlabel .eq. 'random_seed') then
         vartype = 'i'
         ivarvalue => random_seed_value

      elseif (varlabel .eq. 'max_diffusion_simulation_time') then
         vartype = 'r'
         rvarvalue => max_diffusion_simulation_time

      elseif (varlabel .eq. 'max_diffusion_cpu_time') then
         vartype = 'r'
         rvarvalue => max_diffusion_cpu_time

      elseif (varlabel .eq. 'max_collisions_per_sphere') then
         vartype = 'r'
         rvarvalue => max_collisions_per_sphere

      elseif (varlabel .eq. 'number_configurations') then
         vartype = 'i'
         ivarvalue => simulation_config%number_configurations

      elseif (varlabel .eq. 'sphere_volume_fraction') then
         vartype = 'r'
         rvarvalue => simulation_config%sphere_volume_fraction
         simulation_config%number_spheres_specified = .false.

      elseif (varlabel .eq. 'periodic_bc') then
         vartype = 'l'
         varlen = 3
         lavarvalue => periodic_bc

      elseif (varlabel .eq. 'wall_boundary_model') then
         vartype = 'i'
         ivarvalue => wall_boundary_model

      elseif (varlabel .eq. 'auto_target_radius') then
         vartype = 'l'
         lvarvalue => simulation_config%auto_target_radius

      elseif (varlabel .eq. 'target_radius_padding') then
         vartype = 'r'
         rvarvalue => simulation_config%target_radius_padding

      elseif (varlabel .eq. 'sphere_1_fixed') then
         vartype = 'l'
         lvarvalue => sphere_1_fixed

      elseif (varlabel .eq. 'random_lattice_configuration') then
         vartype = 'l'
         lvarvalue => random_lattice_configuration

      elseif (varlabel .eq. 'erase_sphere_1') then
         vartype = 'l'
         lvarvalue => simulation_config%erase_sphere_1

      elseif (varlabel .eq. 'configuration_average') then
         vartype = 'l'
         lvarvalue => simulation_config%configuration_average

      elseif (varlabel .eq. 'frozen_configuration') then
         vartype = 'l'
         lvarvalue => simulation_config%frozen_configuration

      elseif (varlabel .eq. 'random_configuration_host') then
         vartype = 'l'
         lvarvalue => simulation_config%random_configuration_host

      elseif (varlabel .eq. 'host_sphere_ref_index') then
         vartype = 'c'
         cvarvalue => simulation_config%host_sphere_ref_index

      elseif (varlabel .eq. 'random_configuration_host_model') then
         vartype = 'i'
         ivarvalue => simulation_config%random_configuration_host_model

      elseif (varlabel .eq. 'fit_for_radius') then
         vartype = 'l'
         lvarvalue => simulation_config%fit_for_radius

      elseif (varlabel .eq. 'effective_medium_simulation') then
         vartype = 'l'
         lvarvalue => simulation_config%effective_medium_simulation

      elseif (varlabel .eq. 'reflection_model') then
         vartype = 'l'
         lvarvalue => simulation_config%reflection_model

      elseif (varlabel .eq. 'absorption_sample_radius') then
         vartype = 'r'
         rvarvalue => simulation_config%absorption_sample_radius

      elseif (varlabel .eq. 'absorption_sample_radius_fraction') then
         vartype = 'r'
         rvarvalue => simulation_config%absorption_sample_radius_fraction

      elseif (varlabel .eq. 'auto_absorption_sample_radius') then
         vartype = 'l'
         lvarvalue => simulation_config%auto_absorption_sample_radius

      elseif (varlabel .eq. 'print_random_configuration') then
         vartype = 'l'
         lvarvalue => simulation_config%output%print_random_configuration

      elseif (varlabel .eq. 'print_timings') then
         vartype = 'l'
         lvarvalue => simulation_config%output%print_timings

      elseif (varlabel .eq. 'calculate_up_down_scattering') then
         vartype = 'l'
         lvarvalue => simulation_config%input_calculate_up_down_scattering

      elseif (varlabel .eq. 'numerical_hemispherical_integration') then
         vartype = 'l'
         lvarvalue => simulation_config%numerical_hemispherical_integration

      elseif (varlabel .eq. 'x_shift') then
         vartype = 'r'
         rvarvalue => simulation_config%x_shift

      elseif (varlabel .eq. 'y_shift') then
         vartype = 'r'
         rvarvalue => simulation_config%y_shift

      elseif (varlabel .eq. 'z_shift') then
         vartype = 'r'
         rvarvalue => simulation_config%z_shift

      elseif (varlabel .eq. 'shifted_sphere') then
         vartype = 'i'
         ivarvalue => simulation_config%shifted_sphere

      elseif (varlabel .eq. 'check_positions') then
         vartype = 'l'
         lvarvalue => simulation_config%check_positions

      elseif (varlabel .eq. 'medium_ref_index') then
         vartype = 'c'
         cvarvalue => simulation_config%medium_ref_index
         simulation_config%medium_ref_index_specified = .true.
         simulation_config%medium_reim_ref_index_specified = .false.

      elseif (varlabel .eq. 'medium_re_ref_index') then
         vartype = 'r'
         rvarvalue => simulation_config%medium_re_ref_index
         simulation_config%medium_ref_index_specified = .true.
         simulation_config%medium_reim_ref_index_specified = .true.

      elseif (varlabel .eq. 'medium_im_ref_index') then
         vartype = 'r'
         rvarvalue => simulation_config%medium_im_ref_index
         simulation_config%medium_ref_index_specified = .true.
         simulation_config%medium_reim_ref_index_specified = .true.

      elseif (varlabel .eq. 'light_up') then
         vartype = 'l'
         lvarvalue => light_up

      end if

      if (vartype .eq. 'n') then
         varstatus = 1
         if (present(var_status)) var_status = varstatus
         return
      end if
      if (present(var_type)) var_type = vartype
      if (present(i_var_pointer)) i_var_pointer => ivarvalue
      if (present(r_var_pointer)) r_var_pointer => rvarvalue
      if (present(c_var_pointer)) c_var_pointer => cvarvalue

      if (operate) then
         if (vartype .eq. 'i') then
            call apply_integer_input_value(sentvarvalue, &
                                           ivarvalue, var_operation=varop)
         elseif (vartype .eq. 'r') then
            if (varlen .eq. 1) then
               call apply_real_input_value(sentvarvalue, &
                                           rvarvalue, var_operation=varop)
            else
               call apply_real_array_input_value(sentvarvalue, &
                                                 ravarvalue, var_operation=varop, var_len=varlen)
            end if
         elseif (vartype .eq. 'c') then
            if (varlen .eq. 1) then
               call apply_complex_input_value(sentvarvalue, &
                                              cvarvalue, var_operation=varop)
            else
               call apply_complex_array_input_value(sentvarvalue, &
                                                    cavarvalue, var_operation=varop, var_len=varlen)
            end if
         elseif (vartype .eq. 'l') then
            if (varlen .eq. 1) then
               call apply_logical_input_value(sentvarvalue, &
                                              lvarvalue, var_operation=varop)
            else
               call apply_logical_array_input_value(sentvarvalue, &
                                                    lavarvalue, var_operation=varop, var_len=varlen)
            end if
         elseif (vartype .eq. 'a') then
            avarvalue = sentvarvalue
         end if
      end if
   end subroutine process_input_variable

   subroutine parse_input_data(inputfiledata, read_status)
      implicit none
      integer :: readok, n, record_line, spherenum, varstat, rank, stopit, istat
      integer, save :: inputline
      integer, optional :: read_status
      character(len=256) :: parmid, parmval, varop, inputfiledata(*)
      data inputline/1/

      call parallel_rank(mpi_rank=rank)
      readok = 0
      stopit = 0
      varstat = 0
      do while (readok .eq. 0)
         parmid = inputfiledata(inputline)
         inputline = inputline + 1
         if (trim(parmid) .eq. 'run_file') then
            parmval = inputfiledata(inputline)
            inputline = inputline + 1
            if (trim(parmval) .ne. ' ' .and. .not. simulation_config%validation_only) then
               if (rank .eq. 0) then
                  call open_output_file(trim(parmval), sphere_cluster%run_print_unit)
                  if (runtime_failed()) then
                     stopit = 1
                     exit
                  end if
               end if
            end if
            cycle
         end if

         if (parmid(1:1) .eq. '!' .or. parmid(1:1) .eq. '%') then
            cycle
         end if

         if (trim(parmid) .eq. 'loop_variable') then
            simulation_config%loop_job = .true.
            simulation_config%number_nested_loops = simulation_config%number_nested_loops + 1
            n = simulation_config%number_nested_loops
            parmid = inputfiledata(inputline)
            inputline = inputline + 1
            if (trim(parmid) .eq. 'sphere_number') then
               read (inputfiledata(inputline), *) spherenum
               inputline = inputline + 1
               simulation_config%loop_sphere_number(n) = spherenum
               parmid = inputfiledata(inputline)
               inputline = inputline + 1
            else
               spherenum = 1
            end if
            simulation_config%loop_variable_label(n) = parmid
            call process_input_variable(simulation_config%loop_variable_label(n), &
                                        var_type=simulation_config%loop_variable_type(n), var_position=spherenum)
            if (simulation_config%loop_variable_type(n) .eq. 'i') then
               read (inputfiledata(inputline), *) simulation_config%integer_loop_start(n), simulation_config%integer_loop_stop(n), simulation_config%integer_loop_step(n)
            elseif (simulation_config%loop_variable_type(n) .eq. 'r') then
               read (inputfiledata(inputline), *) simulation_config%real_loop_start(n), simulation_config%real_loop_stop(n), simulation_config%real_loop_step(n)
            elseif (simulation_config%loop_variable_type(n) .eq. 'c') then
               read (inputfiledata(inputline), *) simulation_config%complex_loop_start(n), simulation_config%complex_loop_stop(n), simulation_config%complex_loop_step(n)
            end if
            inputline = inputline + 1
            cycle

         elseif (trim(parmid) .eq. 'sphere_data') then
            if (allocated(simulation_config%sphere_data_records)) &
               deallocate (simulation_config%sphere_data_records, simulation_config%sphere_data_record_lines)
            allocate (simulation_config%sphere_data_records(simulation_config%input_number_spheres), &
                      simulation_config%sphere_data_record_lines(simulation_config%input_number_spheres))
            simulation_config%sphere_data_records = ''
            simulation_config%sphere_data_record_lines = 0
            simulation_config%output%sphere_data_file = '<embedded sphere_data>'
            simulation_config%sphere_data_source = simulation_config%input_file
            simulation_config%embedded_sphere_data = .true.
            n = 0
            do
               parmval = inputfiledata(inputline)
               if (trim(parmval) .eq. 'end_of_options') exit
               record_line = inputline
               inputline = inputline + 1
               if (trim(parmval) .eq. 'end_of_sphere_data') exit
               parmval = adjustl(parmval)
               if (len_trim(parmval) == 0) cycle
               if (parmval(1:1) .eq. '!' .or. parmval(1:1) .eq. '%') cycle
               n = n + 1
               if (n <= simulation_config%input_number_spheres) then
                  simulation_config%sphere_data_records(n) = parmval
                  simulation_config%sphere_data_record_lines(n) = record_line
               end if
            end do
            simulation_config%number_sphere_data_records = min(n, simulation_config%input_number_spheres)
            simulation_config%data_scaled = .false.
            sphere_cluster%recalculate_surface_matrix = .true.
            cycle

         elseif (trim(parmid) .eq. 'new_run') then
            simulation_config%repeat_run = .true.
            exit

         elseif (trim(parmid) .eq. 'end_of_options') then
            simulation_config%repeat_run = .false.
            readok = -1
            exit

         elseif (trim(parmid) .eq. 'layer_ref_index') then
            if (number_plane_boundaries .gt. max_number_plane_boundaries) then
               if (rank .eq. 0) write (sphere_cluster%run_print_unit, '('' max # plane boundaries exceeded:'',i3,''>'',i3)') &
                  number_plane_boundaries, max_number_plane_boundaries
               call set_runtime_error('Maximum number of plane boundaries exceeded')
               stopit = 1
               exit
            end if
            parmval = inputfiledata(inputline)
            inputline = inputline + 1
            read (parmval, *, iostat=istat) layer_ref_index(0)
            layer_ref_index(1:max(1, number_plane_boundaries)) = layer_ref_index(0)
            read (parmval, *, iostat=istat) layer_ref_index(0:number_plane_boundaries)
            sphere_cluster%recalculate_surface_matrix = .true.
            simulation_config%medium_ref_index_specified = .false.

         elseif (trim(parmid) .eq. 'layer_thickness') then
            parmval = inputfiledata(inputline)
            inputline = inputline + 1
            simulation_config%input_layer_thickness(1:max(1, number_plane_boundaries)) = 0.d0
            read (parmval, *, iostat=istat) simulation_config%input_layer_thickness(1:max(1, number_plane_boundaries))
            sphere_cluster%recalculate_surface_matrix = .true.

         elseif (trim(parmid) .eq. 'component_radii') then
            call read_real_list(component_radii, number_components)

         elseif (trim(parmid) .eq. 'component_number_fraction') then
            call read_real_list(component_number_fraction, number_components)

         elseif (trim(parmid) .eq. 'psd_sigma') then
            call read_real_list(psd_sigma, number_components)

         elseif (trim(parmid) .eq. 'component_ref_index') then
            call read_complex_list(simulation_config%component_ref_index, number_components)

         else
            varstat = 0
            call process_input_variable(parmid, &
                                        var_status=varstat)
            if (varstat .ne. 0) then
               stopit = 1
               call set_runtime_error('Unknown input parameter: '//trim(parmid))
               if (rank .eq. 0) then
                  write (sphere_cluster%run_print_unit, '('' unknown input parameter:'',a)') trim(parmid)
                  flush (sphere_cluster%run_print_unit)
               end if
               exit
            else
               parmval = inputfiledata(inputline)
               inputline = inputline + 1
               if (readok .ne. 0) cycle
               parmval = trim(parmval)
               varop = 'assign'
               call process_input_variable(parmid, var_value=parmval, &
                                           var_position=1, var_operation='assign', &
                                           var_status=varstat)
            end if
         end if
      end do
      if (stopit .eq. 1) varstat = 1
      if (present(read_status)) read_status = varstat

   contains
      subroutine read_real_list(listvar, listnum)
         implicit none
         integer :: listnum
         real(real64) :: listvar(*)
         parmval = inputfiledata(inputline)
         inputline = inputline + 1
         read (parmval, *, iostat=istat) listvar(1:listnum)
      end subroutine read_real_list

      subroutine read_complex_list(listvar, listnum)
         implicit none
         integer :: listnum
         complex(real64) :: listvar(*)
         parmval = inputfiledata(inputline)
         inputline = inputline + 1
         read (parmval, *, iostat=istat) listvar(1:listnum)
      end subroutine read_complex_list
   end subroutine parse_input_data
end module input_parser
