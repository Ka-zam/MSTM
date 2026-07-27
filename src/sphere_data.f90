module sphere_data
   use, intrinsic :: iso_fortran_env, only: real64
   use numerical_tables
   use surface
   use periodic_lattice_operations
   type linked_sphere_list
      integer :: sphere
      type(linked_sphere_list), pointer :: next
   end type linked_sphere_list
   type host_list
      integer :: number
      type(linked_sphere_list), pointer :: sphere_list
   end type host_list
   type(host_list), allocatable :: sphere_links(:, :)
   logical :: one_side_only, recalculate_surface_matrix, any_optically_active
   logical, target :: store_translation_matrix, effective_medium_simulation, &
                      store_surface_matrix, fft_translation_option, normalize_solution_error
   logical, allocatable :: optically_active(:)
   integer, target :: number_spheres, run_print_unit, translation_switch_order, &
                      max_t_matrix_order
   integer :: max_mie_order, number_host_spheres, number_eqns, t_matrix_order, max_sphere_depth
   integer, allocatable, target :: host_sphere(:), sphere_order(:), sphere_block(:), sphere_offset(:)
   integer, allocatable :: translation_order(:), number_field_expansions(:), mie_offset(:), &
                           mie_block_offset(:), sphere_layer(:), sphere_depth(:)
   real(real64) :: cluster_origin(3), vol_radius, sphere_mean_position(3), area_mean_radius, &
                   sphere_min_position(3), sphere_max_position(3), mean_qext_mie, mean_qabs_mie, &
                   circumscribing_radius, cross_section_radius, effective_cluster_radius
   real(real64), target :: gaussian_beam_constant, gaussian_beam_focal_point(3)
   real(real64), allocatable :: qext_mie(:), qabs_mie(:)
   real(real64), allocatable, target :: sphere_radius(:), sphere_position(:, :)
   complex(real64) :: effective_ref_index
   complex(real64), allocatable :: an_mie(:), cn_mie(:), un_mie(:), vn_mie(:), dn_mie(:), an_inv_mie(:)
   complex(real64), allocatable, target :: sphere_ref_index(:, :)

   data run_print_unit/6/
   data store_translation_matrix, store_surface_matrix/.false., .true./
   data recalculate_surface_matrix/.true./
   data translation_switch_order/3/
   data fft_translation_option/.false./
   data max_t_matrix_order/100/
   data normalize_solution_error/.true./
   data gaussian_beam_constant, gaussian_beam_focal_point/0.d0, 0.d0, 0.d0, 0.d0/
   data effective_medium_simulation, effective_ref_index/.false., (1.d0, 0.d0)/

contains

   subroutine initialize_sphere_layers()
      implicit none
      type(linked_sphere_list), pointer :: slist
      integer :: i, j, l

      if (allocated(sphere_layer)) deallocate (sphere_layer)
      allocate (sphere_layer(number_spheres))
      sphere_layer = 0
      if (number_plane_boundaries .gt. 0) then
         do i = 1, number_spheres
            do j = 1, number_plane_boundaries
               if (sphere_position(3, i) .gt. plane_boundary_position(j)) then
                  sphere_layer(i) = j
               else
                  exit
               end if
            end do
         end do
      end if

      if (number_plane_boundaries .gt. 0) then
         top_boundary = max(plane_boundary_position(max(1, number_plane_boundaries)), sphere_max_position(3)) + 1.d-5
         bot_boundary = min(0.d0, sphere_min_position(3)) - 1.d-5
      else
         top_boundary = sphere_max_position(3) + 1.d-5
         bot_boundary = sphere_min_position(3) - 1.d-5
      end if
      top_boundary = max(plane_boundary_position(max(1, number_plane_boundaries)), sphere_max_position(3)) + 1.d-5
      bot_boundary = min(0.d0, sphere_min_position(3)) - 1.d-5

      if (allocated(sphere_links)) then
         do l = 0, ubound(sphere_links, 2) - 1
            do i = 0, ubound(sphere_links, 1) - 1
               call clear_host_list(sphere_links(i, l))
            end do
         end do
         deallocate (sphere_links)
      end if
      allocate (sphere_links(0:number_spheres, 0:number_plane_boundaries))
      do l = 0, number_plane_boundaries
         do i = 0, number_spheres
            sphere_links(i, l)%number = 0
            allocate (sphere_links(i, l)%sphere_list)
         end do
      end do
      do l = 0, number_plane_boundaries
         do j = 0, number_spheres
            slist => sphere_links(j, l)%sphere_list
            do i = 1, number_spheres
               if (host_sphere(i) .eq. j .and. sphere_layer(i) .eq. l) then
                  sphere_links(j, l)%number = sphere_links(j, l)%number + 1
                  slist%sphere = i
                  allocate (slist%next)
                  slist => slist%next
               end if
            end do
         end do
      end do

      if (allocated(sphere_depth)) deallocate (sphere_depth)
      allocate (sphere_depth(number_spheres))
      max_sphere_depth = 0
      do i = 1, number_spheres
         sphere_depth(i) = 0
         j = host_sphere(i)
         do while (j .ne. 0)
            sphere_depth(i) = sphere_depth(i) + 1
            j = host_sphere(j)
         end do
         max_sphere_depth = max(max_sphere_depth, sphere_depth(i))
      end do

   end subroutine initialize_sphere_layers

   subroutine clear_host_list(hlist)
      implicit none
      type(host_list) :: hlist
      type(linked_sphere_list), pointer :: slist, slist2
      integer :: n, i
      if (.not. associated(hlist%sphere_list)) return
      n = hlist%number
      slist => hlist%sphere_list
      do i = 1, n
         slist2 => slist%next
         deallocate (slist)
         slist => slist2
      end do
   end subroutine clear_host_list
!
!  find_host_spheres finds the host sphere of each sphere in the set.   host=0 for
!  external sphere.
!
!  december 2011
!  march 2013: something changed
!
   subroutine find_host_spheres()
      implicit none
      integer :: i, j, m
      real(real64) :: xij(3), rij, xspmin, zmax, zmin

      host_sphere = 0
      sphere_mean_position = 0.d0
      sphere_min_position = 1.d10
      sphere_max_position = -1.d10
      do i = 1, number_spheres
         sphere_mean_position = sphere_mean_position + sphere_position(:, i)
         do m = 1, 3
            sphere_min_position(m) = min(sphere_min_position(m), sphere_position(m, i))
            sphere_max_position(m) = max(sphere_max_position(m), sphere_position(m, i))
         end do
         xspmin = 1.d6
         do j = 1, number_spheres
            if (sphere_radius(j) .gt. sphere_radius(i)) then
               xij(:) = sphere_position(:, j) - sphere_position(:, i)
               rij = sqrt(dot_product(xij, xij))
!                  if(rij.le.sphere_radius(j)-sphere_radius(i).and.sphere_radius(j).lt.xspmin) then
               if (rij .le. sphere_radius(j) .and. sphere_radius(j) .lt. xspmin) then
                  host_sphere(i) = j
                  xspmin = sphere_radius(j)
               end if
            end if
         end do
      end do
      sphere_mean_position = sphere_mean_position / dble(number_spheres)
      number_field_expansions = 1
      number_host_spheres = 0
      circumscribing_radius = 0.d0
      do i = 1, number_spheres
         rij = sqrt(sum((sphere_position(:, i) - cluster_origin(:))**2)) + sphere_radius(i)
!            rij=sqrt(sum((sphere_position(:,i)-sphere_mean_position(:))**2))+sphere_radius(i)
         circumscribing_radius = max(circumscribing_radius, rij)
         j = host_sphere(i)
         if (j .ne. 0) then
            number_field_expansions(j) = 2
            number_host_spheres = number_host_spheres + 1
         end if
      end do
      zmax = maxval(sphere_position(3, 1:number_spheres))
      zmin = minval(sphere_position(3, 1:number_spheres))
      one_side_only = (zmax * zmin .gt. 0.d0)

   end subroutine find_host_spheres

end module sphere_data
