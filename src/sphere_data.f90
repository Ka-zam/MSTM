module sphere_data
   use, intrinsic :: iso_fortran_env, only: real64
   use numerical_tables
   use surface
   use periodic_lattice_operations
   implicit none(type, external)

   integer, parameter, public :: material_dielectric = 0
   integer, parameter, public :: material_pec = 1
   public :: sphere_record_is_pec, sphere_record_mentions_pec
   private :: lowercase

   type sphere_index_list
      integer, allocatable :: indices(:)
   end type sphere_index_list

   type, public :: sphere_cluster_t
      logical :: one_side_only = .false.
      logical :: recalculate_surface_matrix = .true.
      logical :: any_optically_active = .false.
      logical :: store_translation_matrix = .false.
      logical :: effective_medium_simulation = .false.
      logical :: store_surface_matrix = .true.
      logical :: fft_translation_option = .false.
      logical :: normalize_solution_error = .true.
      logical, allocatable :: optically_active(:)
      integer :: number_spheres = 0
      integer :: run_print_unit = 6
      integer :: translation_switch_order = 3
      integer :: max_t_matrix_order = 100
      integer :: max_mie_order = 0
      integer :: number_host_spheres = 0
      integer :: number_eqns = 0
      integer :: t_matrix_order = 0
      integer :: max_sphere_depth = 0
      integer, allocatable :: host_sphere(:), sphere_order(:), sphere_block(:), sphere_offset(:)
      integer, allocatable :: material_model(:)
      integer, allocatable :: translation_order(:), number_field_expansions(:), mie_offset(:), &
                              mie_block_offset(:), sphere_layer(:), sphere_depth(:)
      type(sphere_index_list), allocatable :: sphere_links(:, :)
      real(real64) :: cluster_origin(3) = 0.0_real64
      real(real64) :: vol_radius = 0.0_real64
      real(real64) :: sphere_mean_position(3) = 0.0_real64
      real(real64) :: area_mean_radius = 0.0_real64
      real(real64) :: sphere_min_position(3) = 0.0_real64
      real(real64) :: sphere_max_position(3) = 0.0_real64
      real(real64) :: mean_qext_mie = 0.0_real64
      real(real64) :: mean_qabs_mie = 0.0_real64
      real(real64) :: circumscribing_radius = 0.0_real64
      real(real64) :: cross_section_radius = 0.0_real64
      real(real64) :: effective_cluster_radius = 0.0_real64
      real(real64) :: gaussian_beam_constant = 0.0_real64
      real(real64) :: gaussian_beam_focal_point(3) = 0.0_real64
      real(real64), allocatable :: qext_mie(:), qabs_mie(:)
      real(real64), allocatable :: sphere_radius(:), sphere_position(:, :)
      complex(real64) :: effective_ref_index = (1.0_real64, 0.0_real64)
      complex(real64), allocatable :: an_mie(:), cn_mie(:), un_mie(:), vn_mie(:), dn_mie(:), an_inv_mie(:)
      complex(real64), allocatable :: sphere_ref_index(:, :)
   contains
      procedure, public :: initialize_layers => initialize_sphere_layers
      procedure, public :: find_hosts => find_host_spheres
      procedure, public :: is_pec => sphere_is_pec
   end type sphere_cluster_t

   type(sphere_cluster_t), target, public :: sphere_cluster

contains

   logical function sphere_record_is_pec(record)
      character(len=*), intent(in) :: record
      character(len=32) :: material_name
      integer :: io_status
      real(real64) :: geometry(4)

      read (record, *, iostat=io_status) geometry, material_name
      sphere_record_is_pec = io_status == 0 .and. lowercase(trim(material_name)) == 'pec'
   end function sphere_record_is_pec

   pure logical function sphere_record_mentions_pec(record)
      character(len=*), intent(in) :: record
      character(len=len(record)) :: lower_record
      integer :: i, record_length
      logical :: left_delimiter, right_delimiter

      lower_record = lowercase(record)
      record_length = len_trim(lower_record)
      sphere_record_mentions_pec = .false.
      do i = 1, record_length - 2
         if (lower_record(i:i + 2) /= 'pec') cycle
         left_delimiter = i == 1
         if (.not. left_delimiter) left_delimiter = index(' ,'//achar(9), lower_record(i - 1:i - 1)) > 0
         right_delimiter = i + 2 == record_length
         if (.not. right_delimiter) right_delimiter = index(' ,!'//achar(9), lower_record(i + 3:i + 3)) > 0
         if (left_delimiter .and. right_delimiter) then
            sphere_record_mentions_pec = .true.
            return
         end if
      end do
   end function sphere_record_mentions_pec

   pure function lowercase(value) result(lower_value)
      character(len=*), intent(in) :: value
      character(len=len(value)) :: lower_value
      integer :: character_code, i

      lower_value = value
      do i = 1, len(value)
         character_code = iachar(value(i:i))
         if (character_code >= iachar('A') .and. character_code <= iachar('Z')) &
            lower_value(i:i) = achar(character_code + iachar('a') - iachar('A'))
      end do
   end function lowercase

   pure logical function sphere_is_pec(self, sphere_index)
      class(sphere_cluster_t), intent(in) :: self
      integer, intent(in) :: sphere_index

      sphere_is_pec = allocated(self%material_model) .and. sphere_index >= 1
      if (sphere_is_pec) sphere_is_pec = sphere_index <= size(self%material_model)
      if (sphere_is_pec) sphere_is_pec = self%material_model(sphere_index) == material_pec
   end function sphere_is_pec

   subroutine initialize_sphere_layers(self)
      implicit none(type, external)
      class(sphere_cluster_t), intent(inout) :: self
      integer :: i, j, l, member, member_count

      if (allocated(self%sphere_layer)) deallocate (self%sphere_layer)
      allocate (self%sphere_layer(self%number_spheres))
      self%sphere_layer = 0
      if (number_plane_boundaries .gt. 0) then
         do i = 1, self%number_spheres
            do j = 1, number_plane_boundaries
               if (self%sphere_position(3, i) .gt. plane_boundary_position(j)) then
                  self%sphere_layer(i) = j
               else
                  exit
               end if
            end do
         end do
      end if

      if (number_plane_boundaries .gt. 0) then
         top_boundary = max(plane_boundary_position(max(1, number_plane_boundaries)), self%sphere_max_position(3)) + 1.d-5
         bot_boundary = min(0.d0, self%sphere_min_position(3)) - 1.d-5
      else
         top_boundary = self%sphere_max_position(3) + 1.d-5
         bot_boundary = self%sphere_min_position(3) - 1.d-5
      end if
      if (allocated(self%sphere_links)) deallocate (self%sphere_links)
      allocate (self%sphere_links(0:self%number_spheres, 0:number_plane_boundaries))
      do l = 0, number_plane_boundaries
         do j = 0, self%number_spheres
            member_count = count(self%host_sphere == j .and. self%sphere_layer == l)
            allocate (self%sphere_links(j, l)%indices(member_count))
            member = 0
            do i = 1, self%number_spheres
               if (self%host_sphere(i) .eq. j .and. self%sphere_layer(i) .eq. l) then
                  member = member + 1
                  self%sphere_links(j, l)%indices(member) = i
               end if
            end do
         end do
      end do

      if (allocated(self%sphere_depth)) deallocate (self%sphere_depth)
      allocate (self%sphere_depth(self%number_spheres))
      self%max_sphere_depth = 0
      do i = 1, self%number_spheres
         self%sphere_depth(i) = 0
         j = self%host_sphere(i)
         do while (j .ne. 0)
            self%sphere_depth(i) = self%sphere_depth(i) + 1
            j = self%host_sphere(j)
         end do
         self%max_sphere_depth = max(self%max_sphere_depth, self%sphere_depth(i))
      end do

   end subroutine initialize_sphere_layers

!  find_host_spheres finds the host sphere of each sphere in the set.   host=0 for
!  external sphere.
!
!  december 2011
!  march 2013: something changed
!
   subroutine find_host_spheres(self)
      implicit none(type, external)
      class(sphere_cluster_t), intent(inout) :: self
      integer :: i, j, m
      real(real64) :: xij(3), rij, xspmin, zmax, zmin

      self%host_sphere = 0
      self%sphere_mean_position = 0.d0
      self%sphere_min_position = 1.d10
      self%sphere_max_position = -1.d10
      do i = 1, self%number_spheres
         self%sphere_mean_position = self%sphere_mean_position + self%sphere_position(:, i)
         do m = 1, 3
            self%sphere_min_position(m) = min(self%sphere_min_position(m), self%sphere_position(m, i))
            self%sphere_max_position(m) = max(self%sphere_max_position(m), self%sphere_position(m, i))
         end do
         xspmin = 1.d6
         do j = 1, self%number_spheres
            if (self%sphere_radius(j) .gt. self%sphere_radius(i)) then
               xij(:) = self%sphere_position(:, j) - self%sphere_position(:, i)
               rij = sqrt(dot_product(xij, xij))
!                  if(rij.le.self%sphere_radius(j)-self%sphere_radius(i).and.self%sphere_radius(j).lt.xspmin) then
               if (rij .le. self%sphere_radius(j) .and. self%sphere_radius(j) .lt. xspmin) then
                  self%host_sphere(i) = j
                  xspmin = self%sphere_radius(j)
               end if
            end if
         end do
      end do
      self%sphere_mean_position = self%sphere_mean_position / dble(self%number_spheres)
      self%number_field_expansions = 1
      self%number_host_spheres = 0
      self%circumscribing_radius = 0.d0
      do i = 1, self%number_spheres
         rij = sqrt(sum((self%sphere_position(:, i) - self%cluster_origin(:))**2)) + self%sphere_radius(i)
!            rij=sqrt(sum((self%sphere_position(:,i)-self%sphere_mean_position(:))**2))+self%sphere_radius(i)
         self%circumscribing_radius = max(self%circumscribing_radius, rij)
         j = self%host_sphere(i)
         if (j .ne. 0) then
            self%number_field_expansions(j) = 2
            self%number_host_spheres = self%number_host_spheres + 1
         end if
      end do
      zmax = maxval(self%sphere_position(3, 1:self%number_spheres))
      zmin = minval(self%sphere_position(3, 1:self%number_spheres))
      self%one_side_only = (zmax * zmin .gt. 0.d0)

   end subroutine find_host_spheres

end module sphere_data
