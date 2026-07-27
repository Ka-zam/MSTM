module input_state
   use, intrinsic :: iso_fortran_env, only: real64
   use solver
   use sphere_data
   use translation_expansions, only: interaction_radius
   use mie
   use near_field
   use surface
   use periodic_lattice_operations
   use random_sphere_configuration
   implicit none(type, external)

   type, public :: solver_options_t
      integer :: max_iterations = 10000
      integer :: t_matrix_procs_per_solution = 4
      real(real64) :: solution_epsilon = 1.0e-6_real64
      real(real64) :: mie_epsilon = 1.0e-6_real64
      real(real64) :: translation_epsilon = 1.0e-5_real64
      real(real64) :: t_matrix_convergence_epsilon = 1.0e-6_real64
      character(len=256) :: solution_method = 'iteration'
   end type solver_options_t

   type, public :: fft_options_t
      logical :: enabled = .false.
      logical :: cell_size_specified = .false.
      integer :: minimum_spheres = 200
      integer :: node_order = -1
      integer :: neighbor_node_model = 2
      real(real64) :: cell_volume_fraction = 0.0_real64
      real(real64) :: cell_size = 0.0_real64
   end type fft_options_t

   type, public :: output_options_t
      logical :: append = .false.
      logical :: append_near_field = .false.
      logical :: print_scattering_matrix = .false.
      logical :: copy_input = .false.
      logical :: normalize_s11 = .true.
      logical :: print_sphere_data = .true.
      logical :: print_random_configuration = .false.
      logical :: print_timings = .true.
      character(len=256) :: output_file = 'mstmtest.dat'
      character(len=256) :: run_file = 'run1.dat'
      character(len=256) :: t_matrix_file = 'tmattemp.dat'
      character(len=256) :: sphere_data_file = ''
      character(len=256) :: near_field_file = 'nftest.dat'
      character(len=256) :: random_configuration_file = 'random_configuration.pos'
   end type output_options_t

   type, public :: simulation_config_t
      type(solver_options_t) :: solver
      type(fft_options_t) :: fft
      type(output_options_t) :: output
      logical :: loop_job = .false.
      logical :: repeat_run = .false.
      logical :: first_run = .false.
      logical :: data_scaled = .false.
      logical :: temporary_position_file = .false.
      logical :: incident_beta_specified = .false.
      logical :: number_spheres_specified = .true.
      logical :: square_cell = .false.
      logical :: use_previous_configuration = .false.
      logical :: calculate_up_down_scattering = .false.
      logical :: medium_ref_index_specified = .false.
      logical :: medium_reim_ref_index_specified = .false.
      logical :: calculate_near_field = .false.
      logical :: move_to_front = .false.
      logical :: move_to_back = .false.
      logical :: random_orientation = .false.
      logical :: t_matrix_centered_on_1 = .false.
      logical :: calculate_scattering_matrix = .true.
      logical :: single_origin_expansion = .true.
      logical :: azimuthal_average = .false.
      logical :: incident_frame = .false.
      logical :: configuration_average = .false.
      logical :: frozen_configuration = .false.
      logical :: reflection_model = .false.
      logical :: input_calculate_up_down_scattering = .true.
      logical :: incidence_average = .false.
      logical :: auto_absorption_sample_radius = .true.
      logical :: random_configuration = .false.
      logical :: check_positions = .true.
      logical :: random_configuration_host = .false.
      logical :: fit_for_radius = .true.
      logical :: numerical_azimuthal_average = .false.
      logical :: numerical_hemispherical_integration = .false.
      logical :: effective_medium_simulation = .false.
      logical :: auto_target_radius = .false.
      logical :: erase_sphere_1 = .false.
      integer :: number_nested_loops = 0
      integer :: integer_loop_start(5) = 0
      integer :: integer_loop_stop(5) = 0
      integer :: integer_loop_step(5) = 0
      integer :: loop_sphere_number(5) = 0
      integer :: scattering_map_model = 0
      integer :: scattering_map_dimension = 15
      integer :: near_field_calculation_model = 0
      integer :: incident_direction = 1
      integer :: number_configurations = 0
      integer :: number_incident_directions = 16
      integer :: shifted_sphere = 0
      integer :: number_excited_spheres = 1000000
      integer :: input_number_spheres = 0
      integer :: random_configuration_host_model = 1
      integer :: random_configuration_time_steps = 0
      integer :: number_configuration_groups = 0
      real(real64) :: real_loop_start(5) = 0.0_real64
      real(real64) :: real_loop_stop(5) = 0.0_real64
      real(real64) :: real_loop_step(5) = 0.0_real64
      real(real64) :: input_layer_thickness(max_number_plane_boundaries) = 0.0_real64
      real(real64) :: scattering_matrix_angle_minimum = 0.0_real64
      real(real64) :: scattering_matrix_angle_maximum = 0.0_real64
      real(real64) :: incident_beta_degrees = 0.0_real64
      real(real64) :: incident_alpha_degrees = 0.0_real64
      real(real64) :: length_scale_factor = 1.0_real64
      real(real64) :: near_field_plane_position = 0.0_real64
      real(real64) :: near_field_plane_vertices(3, 2) = reshape( &
                      [-0.5_real64, 0.0_real64, -0.5_real64, 0.5_real64, 0.0_real64, 0.5_real64], [3, 2])
      real(real64) :: near_field_step_size = 0.2_real64
      real(real64) :: scattering_map_increment = 1.0_real64
      real(real64) :: input_cell_width(2) = 0.0_real64
      real(real64) :: sphere_volume_fraction = 0.0_real64
      real(real64) :: input_cell_width_x = 0.0_real64
      real(real64) :: medium_re_ref_index = 1.0_real64
      real(real64) :: medium_im_ref_index = 0.0_real64
      real(real64) :: excitation_radius = 1.0e10_real64
      real(real64) :: absorption_sample_radius = 1.0e10_real64
      real(real64) :: absorption_sample_radius_fraction = 0.8_real64
      real(real64) :: x_shift = 0.0_real64
      real(real64) :: y_shift = 0.0_real64
      real(real64) :: z_shift = 0.0_real64
      real(real64) :: target_radius_padding = 5.0_real64
      complex(real64) :: complex_loop_start(5) = (0.0_real64, 0.0_real64)
      complex(real64) :: complex_loop_stop(5) = (0.0_real64, 0.0_real64)
      complex(real64) :: complex_loop_step(5) = (0.0_real64, 0.0_real64)
      complex(real64) :: ref_index_scale_factor = (1.0_real64, 0.0_real64)
      complex(real64) :: host_sphere_ref_index = (1.0_real64, 0.0_real64)
      complex(real64) :: medium_ref_index = (1.0_real64, 0.0_real64)
      complex(real64) :: component_ref_index(4) = (1.0_real64, 0.0_real64)
      character(len=1) :: loop_variable_type(5) = ''
      character(len=256) :: loop_variable_label(5) = ''
      character(len=256) :: input_file = ''
      logical, allocatable :: sphere_excitation_switch(:)
   end type simulation_config_t

   type, public :: simulation_result_t
      integer :: run_number = 0
      integer :: efficiency_dimension = 0
      integer :: scattering_map_directions = 0
      integer :: local_rank = 0
      integer :: scattering_matrix_lower_bound = 0
      integer :: scattering_matrix_upper_bound = 0
      integer :: scattering_matrix_dimension = 0
      integer :: random_configuration_status = 0
      integer :: random_configuration_number = 1
      integer :: solution_iterations = 0
      integer :: incident_direction_number = 0
      integer :: reflection_transmission_direction_counts(2) = 0
      integer :: maximum_reflection_transmission_directions = 0
      integer :: fit_status = 0
      integer :: effective_fit_order = 0
      integer, allocatable :: sphere_index(:)
      real(real64) :: diffuse_scattering_ratio = 0.0_real64
      real(real64) :: coherent_scattering_ratio = 0.0_real64
      real(real64) :: hemispherical_scattering(2, 2) = 0.0_real64
      real(real64) :: evanescent_scattering(2) = 0.0_real64
      real(real64) :: propagating_scattering(2) = 0.0_real64
      real(real64) :: plane_scattering(2, 2) = 0.0_real64
      real(real64) :: average_plane_scattering(2, 2) = 0.0_real64
      real(real64) :: solution_time = 0.0_real64
      real(real64) :: average_solution_time = 0.0_real64
      real(real64) :: incident_beta = 0.0_real64
      real(real64) :: incident_sin_beta = 0.0_real64
      real(real64) :: solution_error = 0.0_real64
      real(real64) :: reciprocal_condition = 1.0_real64
      real(real64) :: surface_absorptance(2) = 0.0_real64
      real(real64) :: average_surface_absorptance(2) = 0.0_real64
      real(real64) :: position_shift(3) = 0.0_real64
      real(real64) :: average_total_cross_section(1) = 0.0_real64
      real(real64) :: diffuse_cross_section_ratio(1) = 0.0_real64
      real(real64) :: fit_radius = 0.0_real64
      real(real64) :: total_cross_section = 0.0_real64
      real(real64) :: diffuse_cross_section = 0.0_real64
      real(real64) :: effective_fit_radius = 0.0_real64
      real(real64), allocatable :: efficiency(:, :, :), volume_absorption(:, :), total_efficiency(:, :)
      real(real64), allocatable :: scattering_matrix(:, :), diffuse_scattering_matrix(:, :)
      real(real64), allocatable :: scattering_coefficients(:, :, :), configuration_scattering_coefficients(:, :, :)
      real(real64), allocatable :: boundary_scattering(:, :), boundary_extinction(:, :)
      real(real64), allocatable :: average_efficiency(:, :, :), average_volume_absorption(:, :)
      real(real64), allocatable :: average_total_efficiency(:, :), average_scattering_matrix(:, :)
      real(real64), allocatable :: average_boundary_scattering(:, :), average_boundary_extinction(:, :)
      real(real64), allocatable :: average_sphere_position(:, :), diffuse_boundary_scattering(:, :)
      real(real64), allocatable :: scattering_matrix_expansion(:, :, :), average_scattering_matrix_expansion(:, :, :)
      real(real64), allocatable :: reflection_transmission_vectors(:, :)
      real(real64), allocatable :: coherent_scattering_expansion(:, :, :), average_coherent_scattering_expansion(:, :, :)
      real(real64), allocatable :: scattering_field(:, :, :, :, :), average_scattering_field(:, :, :, :, :)
      complex(real64) :: near_field_effective_ref_index = (0.0_real64, 0.0_real64)
      complex(real64), allocatable :: solution_coefficients(:, :), average_incident_coefficients(:, :)
      complex(real64), allocatable :: incident_coefficients(:, :), electric_field(:, :, :, :, :)
      complex(real64), allocatable :: average_electric_field(:, :, :, :, :), magnetic_field(:, :, :, :, :)
      complex(real64), allocatable :: average_magnetic_field(:, :, :, :, :), mean_t_matrix(:, :)
      complex(real64), allocatable :: average_mean_t_matrix(:, :), effective_fit_coefficients(:, :)
      character(len=20) :: run_date_and_time = ''
   end type simulation_result_t

   type(simulation_config_t), target, public :: simulation_config
   type(simulation_result_t), target, public :: simulation_result
end module input_state
