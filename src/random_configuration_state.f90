module random_configuration_state
   implicit none

   type coll_list
      logical :: wallcoll
      integer :: wall, sphere
      real(8) :: time, collpos(3)
   end type coll_list
   type l_list
      integer :: index
      type(l_list), pointer :: next
   end type l_list
   type c_list
      integer :: number_elements
      type(l_list), pointer :: members
   end type c_list
   logical :: target_width_specified
   logical, target :: sphere_1_fixed, periodic_bc(3), random_lattice_configuration
   integer :: cell_dim(3)
   integer, allocatable :: sphere_cell(:, :)
   integer, target :: target_shape, wall_boundary_model, max_number_time_steps, number_components
   real(8) :: fv_crit, time_step
   real(8) :: minimum_gap, d_cell, target_boundaries(3, 2)
   real(8), target :: target_dimensions(3), psd_sigma(4), target_width, target_thickness, max_collisions_per_sphere, &
                      max_diffusion_cpu_time, max_diffusion_simulation_time, component_radii(4), &
                      component_number_fraction(4)
   real(8) :: sim_timings(10), time_0
   type(c_list), allocatable :: cell_list(:, :, :)
   type(coll_list), allocatable :: coll_data(:)
   character(len=1) :: c_temp
   data fv_crit, time_step/0.25d0, .1d0/
   data minimum_gap, sphere_1_fixed, target_shape/1.0d-3, .false., 0/
   data periodic_bc/.true., .true., .true./
   data wall_boundary_model/1/
   data max_number_time_steps, max_collisions_per_sphere, max_diffusion_simulation_time/100, 3.d0, 5.d0/
   data max_diffusion_cpu_time/100.d0/
   data number_components, component_radii, component_number_fraction, psd_sigma/1, 4*1.d0, &
      1.d0, 0.d0, 0.d0, 0.d0, 4*0.d0/
end module random_configuration_state
