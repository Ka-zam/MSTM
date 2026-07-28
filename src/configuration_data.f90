module configuration_data
   use, intrinsic :: iso_fortran_env, only: real64
   use input_state
   use parallel_runtime, only: mpi_comm_world, mstm_global_rank, parallel_barrier, parallel_broadcast, parallel_rank
   use runtime_support, only: open_input_file, runtime_failed, set_runtime_error, synchronize_runtime_status
   use sphere_data, only: material_dielectric, material_pec, sphere_cluster, sphere_record_is_pec, &
                          sphere_record_mentions_pec
   implicit none
   private
   public :: check_sphere_positions, generate_random_configuration, read_sphere_data_input_file, &
             validate_material_configuration
contains

   subroutine read_sphere_data_input_file(mpi_comm)
      implicit none
      integer :: n
      integer, optional :: mpi_comm

      if (.not. simulation_config%embedded_sphere_data) call load_external_sphere_data_records()
      if (runtime_failed()) return
      if (simulation_config%number_sphere_data_records < sphere_cluster%number_spheres) then
         call set_runtime_error("Sphere data source '"//trim(simulation_config%sphere_data_source)//"' contains "// &
                                trim(integer_string(simulation_config%number_sphere_data_records))// &
                                ' records; expected '//trim(integer_string(sphere_cluster%number_spheres)))
         return
      end if

      do n = 1, sphere_cluster%number_spheres
         call parse_sphere_data_record(n)
         if (runtime_failed()) return
      end do
      do n = 1, sphere_cluster%number_spheres
         sphere_cluster%sphere_radius(n) = sphere_cluster%sphere_radius(n) * simulation_config%length_scale_factor
         sphere_cluster%sphere_position(:, n) = sphere_cluster%sphere_position(:, n) * simulation_config%length_scale_factor
         if (.not. sphere_cluster%is_pec(n)) &
            sphere_cluster%sphere_ref_index(:, n) = sphere_cluster%sphere_ref_index(:, n) &
                                                    * simulation_config%ref_index_scale_factor
      end do
   end subroutine read_sphere_data_input_file

   subroutine load_external_sphere_data_records()
      integer :: input_unit, io_status, line_number, record_count
      character(len=256) :: record

      simulation_config%sphere_data_source = simulation_config%output%sphere_data_file
      simulation_config%number_sphere_data_records = 0
      if (allocated(simulation_config%sphere_data_records)) &
         deallocate (simulation_config%sphere_data_records, simulation_config%sphere_data_record_lines)
      allocate (simulation_config%sphere_data_records(sphere_cluster%number_spheres), &
                simulation_config%sphere_data_record_lines(sphere_cluster%number_spheres))
      simulation_config%sphere_data_records = ''
      simulation_config%sphere_data_record_lines = 0

      call open_input_file(simulation_config%output%sphere_data_file, input_unit)
      if (runtime_failed()) return
      line_number = 0
      record_count = 0
      do while (record_count < sphere_cluster%number_spheres)
         read (input_unit, '(a)', iostat=io_status) record
         if (io_status /= 0) exit
         line_number = line_number + 1
         record = adjustl(record)
         if (len_trim(record) == 0) cycle
         if (record(1:1) == '!' .or. record(1:1) == '%') cycle
         record_count = record_count + 1
         simulation_config%sphere_data_records(record_count) = record
         simulation_config%sphere_data_record_lines(record_count) = line_number
      end do
      close (input_unit)
      simulation_config%number_sphere_data_records = record_count
   end subroutine load_external_sphere_data_records

   subroutine parse_sphere_data_record(record_index)
      integer, intent(in) :: record_index
      integer :: io_status
      real(real64) :: geometry(4)
      complex(real64) :: refractive_index(2)
      character(len=32) :: material_name
      character(len=256) :: record

      record = simulation_config%sphere_data_records(record_index)
      sphere_cluster%sphere_radius(record_index) = 1.0_real64
      sphere_cluster%sphere_ref_index(:, record_index) = (1.0_real64, 0.0_real64)
      sphere_cluster%material_model(record_index) = material_dielectric

      read (record, *, iostat=io_status) sphere_cluster%sphere_position(:, record_index)
      if (io_status /= 0) then
         call set_sphere_data_record_error(record_index, 'expected numeric x, y, and z coordinates')
         return
      end if

      read (record, *, iostat=io_status) geometry
      if (io_status == 0) sphere_cluster%sphere_radius(record_index) = geometry(4)
      if (sphere_record_mentions_pec(record) .and. .not. sphere_record_is_pec(record)) then
         call set_sphere_data_record_error(record_index, 'PEC sphere records require x, y, z, radius, and PEC')
         return
      end if
      if (sphere_record_is_pec(record)) then
         read (record, *, iostat=io_status) geometry, material_name
         if (io_status /= 0) then
            call set_sphere_data_record_error(record_index, 'PEC sphere records require x, y, z, radius, and PEC')
            return
         end if
         sphere_cluster%sphere_radius(record_index) = geometry(4)
         sphere_cluster%material_model(record_index) = material_pec
         return
      end if

      read (record, *, iostat=io_status) geometry, refractive_index(1)
      if (io_status == 0) sphere_cluster%sphere_ref_index(1, record_index) = refractive_index(1)
      read (record, *, iostat=io_status) geometry, refractive_index
      if (io_status == 0) then
         sphere_cluster%sphere_ref_index(2, record_index) = refractive_index(2)
      else
         sphere_cluster%sphere_ref_index(2, record_index) = sphere_cluster%sphere_ref_index(1, record_index)
      end if
   end subroutine parse_sphere_data_record

   subroutine set_sphere_data_record_error(record_index, detail)
      integer, intent(in) :: record_index
      character(len=*), intent(in) :: detail
      character(len=512) :: message

      message = 'Sphere record '//trim(integer_string(record_index))//" in '"// &
                trim(simulation_config%sphere_data_source)//"'"
      if (allocated(simulation_config%sphere_data_record_lines)) then
         if (record_index <= size(simulation_config%sphere_data_record_lines)) &
            message = trim(message)//' at line '// &
                      trim(integer_string(simulation_config%sphere_data_record_lines(record_index)))
      end if
      message = trim(message)//': '//trim(detail)
      if (allocated(simulation_config%sphere_data_records)) then
         if (record_index <= size(simulation_config%sphere_data_records)) &
            message = trim(message)//'. Offending record: '//trim(simulation_config%sphere_data_records(record_index))
      end if
      call set_runtime_error(trim(message))
   end subroutine set_sphere_data_record_error

   subroutine generate_random_configuration(mpi_comm, skip_diffusion)
      implicit none
      logical :: skipdif, firstrun, frozen
      logical, optional :: skip_diffusion
      integer :: mpicomm, rank, nsend, nsphere, nspheresamp, i, generation_status(1)
      integer, optional :: mpi_comm
      integer, allocatable, save :: sphereindex(:)
      real(real64) :: targetdimensions(3), crad
      real(real64), allocatable, save :: sphereradius(:), sphereposition(:, :)
      data firstrun/.true./
      if (present(skip_diffusion)) then
         skipdif = skip_diffusion
      else
         skipdif = .false.
      end if
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      nsphere = sphere_cluster%number_spheres
      if (simulation_config%random_configuration_host) nsphere = nsphere - 1
      if (rank .eq. 0) then
        frozen = ((.not. firstrun) .and. (simulation_config%frozen_configuration .or. simulation_config%use_previous_configuration))
         if (frozen) then
            sphere_cluster%sphere_position(:, 1:nsphere) = sphereposition(:, 1:nsphere)
            sphere_cluster%sphere_radius(1:nsphere) = sphereradius(1:nsphere)
         else
            if (simulation_config%auto_target_radius .and. target_shape .eq. 2) then
               targetdimensions = target_dimensions + simulation_config%target_radius_padding
!                  nspheresamp=ceiling(simulation_config%sphere_volume_fraction*(targetdimensions(1)-1.d0)**3.d0)
               nspheresamp = ceiling(simulation_config%sphere_volume_fraction * (targetdimensions(1))**3.d0)
               nspheresamp = max(nspheresamp, nsphere)
               if (mstm_global_rank .eq. 0) then
                  write (sphere_cluster%run_print_unit, '('' set, sampled number spheres:'',2i6)') nsphere, nspheresamp
               end if
            else
               targetdimensions = target_dimensions
               nspheresamp = nsphere
            end if
            if (allocated(sphereradius)) deallocate (sphereradius, sphereposition, sphereindex)
            allocate (sphereradius(nspheresamp), sphereposition(3, nspheresamp), sphereindex(nspheresamp))
            call generate_random_sphere_cluster(nspheresamp, targetdimensions, sphereposition, sphereradius, &
                                               sphereindex, sphere_cluster%run_print_unit, simulation_result%random_configuration_status, simulation_config%random_configuration_time_steps, &
                                                skip_diffusion=skipdif, print_progress=.true.)
            firstrun = .false.
            sphere_cluster%sphere_position(:, 1:nsphere) = sphereposition(:, 1:nsphere)
            sphere_cluster%sphere_radius(1:nsphere) = sphereradius(1:nsphere)
            simulation_result%sphere_index(1:nsphere) = sphereindex(1:nsphere)
         end if
      end if
      call synchronize_runtime_status(mpicomm)
      if (runtime_failed()) return
      generation_status(1) = simulation_result%random_configuration_status
      call parallel_broadcast(send_buffer=generation_status, mpi_number=1, &
                              mpi_rank=0, mpi_comm=mpicomm)
      simulation_result%random_configuration_status = generation_status(1)
      if (generation_status(1) .ge. 3) then
         call set_runtime_error('Unable to generate random sphere configuration')
         return
      end if
!         call parallel_barrier()
      nsend = nsphere
      call parallel_broadcast(send_buffer=sphere_cluster%sphere_radius, &
                              mpi_number=nsend, mpi_rank=0, mpi_comm=mpicomm)
      call parallel_broadcast(send_buffer=simulation_result%sphere_index, &
                              mpi_number=nsend, mpi_rank=0, mpi_comm=mpicomm)
      nsend = nsphere * 3
      call parallel_broadcast(send_buffer=sphere_cluster%sphere_position, &
                              mpi_number=nsend, mpi_rank=0, mpi_comm=mpicomm)
      sphere_cluster%sphere_radius(1:nsphere) = sphere_cluster%sphere_radius(1:nsphere) * simulation_config%length_scale_factor
 sphere_cluster%sphere_position(:, 1:nsphere) = sphere_cluster%sphere_position(:, 1:nsphere) * simulation_config%length_scale_factor
      do i = 1, nsphere
         sphere_cluster%material_model(i) = material_dielectric
         sphere_cluster%sphere_ref_index(:, i) = simulation_config%ref_index_scale_factor &
                                                 * simulation_config%component_ref_index(simulation_result%sphere_index(i))
      end do
      if (simulation_config%erase_sphere_1) sphere_cluster%sphere_ref_index(:, 1) = (1.000001d0, 0.d0)
      if (simulation_config%random_configuration_host) then
         crad = maxval(sqrt(sum(sphere_cluster%sphere_position(:, 1:nsphere)**2, 1))) + 0.01d0
         if (simulation_config%auto_target_radius .and. target_shape .eq. 2) then
            sphere_cluster%sphere_radius(sphere_cluster%number_spheres) &
  = (sum(sphere_cluster%sphere_radius(1:sphere_cluster%number_spheres - 1)**3.) / simulation_config%sphere_volume_fraction)**0.33333
         else
            if (simulation_config%random_configuration_host_model .eq. 1) then
          sphere_cluster%sphere_radius(sphere_cluster%number_spheres) = target_dimensions(1) * simulation_config%length_scale_factor
            elseif (simulation_config%random_configuration_host_model .eq. 2) then
               sphere_cluster%sphere_radius(sphere_cluster%number_spheres) = (sum(sphere_cluster%sphere_radius(1:sphere_cluster%number_spheres - 1)**3.) / simulation_config%sphere_volume_fraction)**0.33333
            end if
         end if
         sphere_cluster%sphere_position(:, sphere_cluster%number_spheres) = 0.d0
         sphere_cluster%material_model(sphere_cluster%number_spheres) = material_dielectric
         sphere_cluster%sphere_ref_index(:, sphere_cluster%number_spheres) = simulation_config%host_sphere_ref_index
      end if
   end subroutine generate_random_configuration

   subroutine validate_material_configuration()
      integer :: i

      do i = 1, sphere_cluster%number_spheres
         if (.not. sphere_cluster%is_pec(i)) cycle
         if (sphere_cluster%host_sphere(i) /= 0) then
            call set_sphere_data_record_error(i, 'nested PEC sphere; only top-level solid PEC spheres are supported')
            return
         end if
         if (any(sphere_cluster%host_sphere == i)) then
            call set_sphere_data_record_error(i, 'PEC sphere contains another sphere; PEC hosts and cavities are unsupported')
            return
         end if
      end do
   end subroutine validate_material_configuration

   function integer_string(value) result(text)
      integer, intent(in) :: value
      character(len=32) :: text

      write (text, '(i0)') value
   end function integer_string

   subroutine check_sphere_positions()
      implicit none
      logical :: check
      integer :: i, j, imin, jmin
      real(real64) :: r, amax, amin, rmingap

      check = .true.
      rmingap = -1.d10
      do i = 1, sphere_cluster%number_spheres - 1
         do j = i + 1, sphere_cluster%number_spheres
            amax = max(sphere_cluster%sphere_radius(i), sphere_cluster%sphere_radius(j))
            amin = min(sphere_cluster%sphere_radius(i), sphere_cluster%sphere_radius(j))
            r = sqrt(sum((sphere_cluster%sphere_position(:, i) - sphere_cluster%sphere_position(:, j))**2))
            if (r .ge. amax + amin) then
               cycle
            else
               if (amin + amax - r .gt. rmingap) then
                  rmingap = amin + amax - r
                  imin = i
                  jmin = j
               end if
            end if
            if (r .le. amax - amin) cycle
            check = .false.
!               write(sphere_cluster%run_print_unit,'('' spheres '',i4,'' and '',i4,'' intersect'')') i,j
!               write(2,'('' spheres '',i4,'' and '',i4,'' intersect'')') i,j
         end do
      end do
      if (.not. check) then
         write (sphere_cluster%run_print_unit, '('' warning: sphere-sphere intersections detected, max overlap:'',es12.4, &
            &''  Results might be garbage!'')') rmingap
   write (sphere_cluster%run_print_unit, '('' positions:'',i5,3es12.4,i5,3es12.4)') imin, sphere_cluster%sphere_position(:, imin), &
            jmin, sphere_cluster%sphere_position(:, jmin)
         flush (sphere_cluster%run_print_unit)
      end if
      check = .true.
      rmingap = -1.d10
      do i = 1, sphere_cluster%number_spheres
         do j = 1, number_plane_boundaries
            if (abs(sphere_cluster%sphere_position(3, i) - plane_boundary_position(j)) .ge. sphere_cluster%sphere_radius(i)) cycle
            check = .false.
    rmingap = max(rmingap, sphere_cluster%sphere_radius(i) - abs(sphere_cluster%sphere_position(3, i) - plane_boundary_position(j)))
!               write(sphere_cluster%run_print_unit,'('' sphere '',i4,'' and plane boundary'',i4,'' intersect'')') i,j
!               write(2,'('' sphere '',i4,'' and plane boundary'',i4,'' intersect'')') i,j
         end do
      end do
      if (.not. check) then
         write (sphere_cluster%run_print_unit, '('' warning: sphere-plane boundary intersections detected, max overlap:'',es12.4, &
            &''  Results might be garbage!'')') rmingap
!            write(sphere_cluster%run_print_unit,'('' positions:'',i5,3es12.4,i5,3es12.4)') imin,sphere_cluster%sphere_position(:,imin), &
!               jmin,sphere_cluster%sphere_position(:,jmin)
         flush (sphere_cluster%run_print_unit)
      end if
   end subroutine check_sphere_positions
end module configuration_data
