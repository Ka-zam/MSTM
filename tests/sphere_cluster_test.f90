program sphere_cluster_test
   use, intrinsic :: iso_fortran_env, only: real64
   use sphere_data, only: sphere_cluster
   use surface, only: number_plane_boundaries
   implicit none(type, external)

   sphere_cluster%number_spheres = 3
   sphere_cluster%cluster_origin = 0.0_real64
   number_plane_boundaries = 0
   allocate (sphere_cluster%sphere_radius(3), sphere_cluster%sphere_position(3, 3), &
             sphere_cluster%host_sphere(3), sphere_cluster%number_field_expansions(3))
   sphere_cluster%sphere_radius = [0.2_real64, 0.2_real64, 1.0_real64]
   sphere_cluster%sphere_position(:, 1) = [0.1_real64, 0.0_real64, 0.0_real64]
   sphere_cluster%sphere_position(:, 2) = [2.0_real64, 0.0_real64, 0.0_real64]
   sphere_cluster%sphere_position(:, 3) = 0.0_real64

   call sphere_cluster%find_hosts()
   call require(all(sphere_cluster%host_sphere == [3, 0, 0]), 'host classification is incorrect')
   call require(sphere_cluster%number_host_spheres == 1, 'nested sphere count is incorrect')

   call sphere_cluster%initialize_layers()
   call require(all(sphere_cluster%sphere_layer == 0), 'layer classification is incorrect')
   call require(all(sphere_cluster%sphere_links(3, 0)%indices == [1]), 'nested host index list is incorrect')
   call require(all(sphere_cluster%sphere_links(0, 0)%indices == [2, 3]), 'external sphere index list is incorrect')

contains

   subroutine require(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) error stop message
   end subroutine require
end program sphere_cluster_test
