!june 18 original
!27 july: core volume fraction added to output

module input_state
   use specialfuncs
   use mpidefs
   use solver
   use spheredata
   use translation, only: interaction_radius
   use mie
   use nearfield
   use scatprops
   use fft_translation
   use surface_subroutines
   use periodic_lattice_subroutines
   use random_sphere_configuration
   implicit none
   logical :: loop_job, repeat_run, first_run, data_scaled, temporary_pos_file, &
              append_near_field_output_file, incident_beta_specified, number_spheres_specified, &
              square_cell, use_previous_configuration, calculate_up_down_scattering, &
              d_cell_specified, medium_ref_index_specified, medium_reim_ref_index_specified
   logical, target :: append_output_file, &
                      print_scattering_matrix, &
                      copy_input_file, calculate_near_field, &
                      move_to_front, move_to_back, random_orientation, &
                      t_matrix_centered_on_1, calculate_scattering_matrix, &
                      normalize_s11, print_sphere_data, single_origin_expansion, &
                      azimuthal_average, incident_frame, configuration_average, &
                      frozen_configuration, reflection_model, input_fft_translation_option, &
                      print_random_configuration, print_timings, &
                      input_calculate_up_down_scattering, incidence_average, auto_absorption_sample_radius, &
                      random_configuration, check_positions, random_configuration_host, fit_for_radius, &
                      numerical_azimuthal_average, numerical_hemispherical_integration, input_effective_medium_simulation, &
                      auto_target_radius, erase_sphere_1
   logical, allocatable :: sphere_excitation_switch(:)
   integer :: n_nest_loops, i_var_start(5), i_var_stop(5), i_var_step(5), &
              run_number, loop_sphere_number(5), qeff_dim, &
              scattering_map_directions, local_rank, scat_mat_ldim, scat_mat_udim, scat_mat_mdim, &
              ran_config_stat, ran_config_time_steps, n_configuration_groups, random_configuration_number, &
              solution_iterations, incident_direction_number, number_rl_dirs(2), max_number_rl_dirs, fit_stat
   integer, target :: max_iterations, t_matrix_procs_per_solution, &
                      scattering_map_model, scattering_map_dimension, near_field_calculation_model, &
                      incident_direction, number_configurations, min_fft_nsphere, input_node_order, &
                      number_incident_directions, shifted_sphere, number_excited_spheres, input_number_spheres, &
                      random_configuration_host_model
   integer, allocatable :: sphere_index(:)
   real(8) :: r_var_start(5), r_var_stop(5), r_var_step(5), diffuse_scattering_ratio, &
              coherent_scattering_ratio, hemispherical_sca(2, 2), evan_sca(2), prop_sca(2), &
              input_layer_thickness(max_number_plane_boundaries), &
              pl_sca(2, 2), scat_mat_amin, scat_mat_amax, pl_sca_ave(2, 2), solution_time, solution_time_ave, &
              incident_beta, solution_error, surface_absorptance(2), surface_absorptance_ave(2), &
              position_shift(3), tot_csca_ave(1), dif_csca_ratio(1), fit_radius, tot_csca, dif_csca
   real(8), allocatable :: q_eff(:, :, :), q_vabs(:, :), q_eff_tot(:, :), scat_mat(:, :), &
                           dif_scat_mat(:, :), sm_coef(:, :, :), sm_cf_coef(:, :, :), boundary_sca(:, :), boundary_ext(:, :), &
                           q_eff_ave(:, :, :), q_vabs_ave(:, :), q_eff_tot_ave(:, :), scat_mat_ave(:, :), &
                           boundary_sca_ave(:, :), boundary_ext_ave(:, :), sphere_position_ave(:, :), dif_boundary_sca(:, :), &
                         scat_mat_exp_coef(:, :, :), scat_mat_exp_coef_ave(:, :, :), rl_vec(:, :), coh_scat_mat_exp_coef(:, :, :), &
                           coh_scat_mat_exp_coef_ave(:, :, :), s_field(:, :, :, :, :), s_field_ave(:, :, :, :, :)
   real(8), target :: incident_beta_deg, incident_alpha_deg, solution_epsilon, &
                      mie_epsilon, length_scale_factor, near_field_plane_position, &
                      near_field_plane_vertices(3, 2), near_field_step_size, &
                      translation_epsilon, t_matrix_convergence_epsilon, &
                      scattering_map_increment, incident_sin_beta, input_cell_width(2), &
                      sphere_volume_fraction, input_cell_width_x, &
                      input_cell_volume_fraction, medium_re_ref_index, medium_im_ref_index, &
                      excitation_radius, absorption_sample_radius, absorption_sample_radius_fraction, &
                      x_shift, y_shift, z_shift, input_d_cell, target_radius_padding
   complex(8) :: c_var_start(5), c_var_stop(5), c_var_step(5), nf_eff_ref_index
   complex(8), target :: ref_index_scale_factor, host_sphere_ref_index, medium_ref_index, &
                         component_ref_index(4)
   complex(8), allocatable :: amnp_s(:, :), amnp_0_ave(:, :), amnp_0(:, :), e_field(:, :, :, :, :), &
                              e_field_ave(:, :, :, :, :), h_field(:, :, :, :, :), h_field_ave(:, :, :, :, :), &
                              mean_t(:, :), mean_t_ave(:, :)
   character(len=1) :: loop_var_type(5)
   character(len=20) :: run_date_and_time
   character(len=256) :: loop_var_label(5), input_file
   character(len=256), target :: output_file, run_file, t_matrix_output_file, &
                                 sphere_data_input_file, near_field_output_file, solution_method, &
                                 random_configuration_output_file
   integer :: effective_fit_order
   real(8) :: effective_fit_radius
   complex(8), allocatable :: effective_fit_coefficients(:, :)
   data loop_job, repeat_run/.false., .false./
   data append_output_file/.false./
   data copy_input_file/.false./
   data n_nest_loops/0/
   data run_number/0/
   data i_var_start, i_var_stop, i_var_step/5*0, 5*0, 5*0/
   data r_var_start, r_var_stop, r_var_step/5*0.d0, 5*0.d0, 5*0.d0/
   data c_var_start, c_var_stop, c_var_step/5*(0.d0, 0.d0), 5*(0.d0, 0.d0), 5*(0.d0, 0.d0)/
   data max_iterations/10000/
   data incident_beta_deg/0.d0/
   data incident_alpha_deg/0.d0/
   data incident_sin_beta, incident_direction, incident_frame/0.d0, 1, .false./
   data solution_epsilon/1.d-6/
   data mie_epsilon/1.d-6/
   data translation_epsilon/1.d-5/
   data t_matrix_convergence_epsilon/1.d-6/
   data output_file/'mstmtest.dat'/
   data run_file/'run1.dat'/
   data length_scale_factor/1.d0/
   data ref_index_scale_factor/(1.d0, 0.d0)/
   data move_to_front, move_to_back/.false., .false./
   data calculate_near_field/.false./
   data near_field_output_file/'nftest.dat'/
   data append_near_field_output_file/.false./
   data near_field_plane_vertices/-.5d0, 0.d0, -.5d0, .5d0, 0.d0, .5d0/
   data near_field_step_size/0.2d0/
   data data_scaled/.false./
   data temporary_pos_file/.false./
   data random_orientation/.false./
   data t_matrix_output_file/'tmattemp.dat'/
   data t_matrix_procs_per_solution/4/
   data t_matrix_centered_on_1/.false./
   data calculate_scattering_matrix/.true./
   data solution_method/'iteration'/
   data scattering_map_model/0/
   data scattering_map_dimension/15/
   data scattering_map_increment/1.d0/
   data normalize_s11, print_sphere_data, single_origin_expansion, azimuthal_average/.true., .true., .true., .false./
   data numerical_azimuthal_average/.false./
   data number_spheres_specified, configuration_average/.true., .false./
   data frozen_configuration, reflection_model, random_configuration_number/.false., .false., 1/
   data random_configuration, random_configuration_output_file/.false., 'random_configuration.pos'/
   data min_fft_nsphere, input_fft_translation_option, input_node_order, input_cell_volume_fraction/200, .false., -1, 0.d0/
   data d_cell_specified/.false./
   data print_random_configuration, print_timings/.false., .true./
   data input_calculate_up_down_scattering/.true./
   data incidence_average, number_incident_directions/.false., 16/
   data use_previous_configuration/.false./
   data absorption_sample_radius, excitation_radius/1.d10, 1.d10/
   data auto_absorption_sample_radius, absorption_sample_radius_fraction/.true., 0.8d0/
   data x_shift, y_shift, z_shift, shifted_sphere/0.d0, 0.d0, 0.d0, 0/
   data erase_sphere_1/.false./
   data check_positions/.true./
   data number_excited_spheres/1000000/
   data host_sphere_ref_index, random_configuration_host, random_configuration_host_model/(1.d0, 0.d0), .false., 1/
   data fit_for_radius/.true./
   data medium_ref_index, medium_ref_index_specified/(1.d0, 0.d0), .false./
   data medium_re_ref_index, medium_im_ref_index, medium_reim_ref_index_specified/1.d0, 0.d0, .false./
   data numerical_hemispherical_integration/.false./
   data input_effective_medium_simulation/.false./
   data auto_target_radius, target_radius_padding/.false., 5.d0/
   data component_ref_index/4*(1.d0, 0.d0)/

end module input_state
