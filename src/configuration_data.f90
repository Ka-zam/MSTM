module configuration_data
   use, intrinsic :: iso_fortran_env, only: real64
   use input_state
   use runtime_support, only: open_input_file, runtime_failed, set_runtime_error, synchronize_runtime_status
   implicit none
   private
   public :: check_sphere_positions, generate_random_configuration, read_sphere_data_input_file
contains

   subroutine read_sphere_data_input_file(mpi_comm)
      implicit none
      integer :: input_unit, mpicomm, rank, istat, n
      integer, optional :: mpi_comm
      real(real64) :: rtemp(4)
      complex(real64) :: ctemp(2)
      character(len=256) :: parmval

      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)

      call open_input_file(sphere_data_input_file, input_unit)
      if (runtime_failed()) return
      do n = 1, number_spheres
         sphere_radius(n) = 1.d0
         sphere_ref_index(1, n) = (1.d0, 0.d0)
         sphere_ref_index(2, n) = (0.d0, 0.d0)
         read (input_unit, '(a)', iostat=istat) parmval
         if (istat .ne. 0) then
            if (rank .eq. 0) then
               write (run_print_unit, '('' insufficient data in input file: '', i4,'' lines, need '',i4)') &
                  n, number_spheres
            end if
            call set_runtime_error('Insufficient data in sphere input file: '//trim(sphere_data_input_file))
            close (input_unit)
            return
         end if
         read (parmval, *, iostat=istat) sphere_position(:, n)
         if (istat .ne. 0) then
            if (rank .eq. 0) then
               write (run_print_unit, '('' read error in sphere data input file'')')
            end if
            call set_runtime_error('Invalid data in sphere input file: '//trim(sphere_data_input_file))
            close (input_unit)
            return
         end if
         read (parmval, *, iostat=istat) rtemp(1:4)
         if (istat .eq. 0) sphere_radius(n) = rtemp(4)
         read (parmval, *, iostat=istat) rtemp(1:4), ctemp(1)
         if (istat .eq. 0) sphere_ref_index(1, n) = ctemp(1)
         read (parmval, *, iostat=istat) rtemp(1:4), ctemp(1:2)
         if (istat .eq. 0) then
            sphere_ref_index(2, n) = ctemp(2)
         else
            sphere_ref_index(2, n) = sphere_ref_index(1, n)
         end if
      end do
      number_spheres = min(n, number_spheres)
      close (input_unit)
      do n = 1, number_spheres
         sphere_radius(n) = sphere_radius(n) * length_scale_factor
         sphere_position(:, n) = sphere_position(:, n) * length_scale_factor
         sphere_ref_index(:, n) = sphere_ref_index(:, n) * ref_index_scale_factor
      end do
   end subroutine read_sphere_data_input_file

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
      call mstm_mpi(mpi_command='rank', mpi_rank=rank, mpi_comm=mpicomm)
      nsphere = number_spheres
      if (random_configuration_host) nsphere = nsphere - 1
      if (rank .eq. 0) then
         frozen = ((.not. firstrun) .and. (frozen_configuration .or. use_previous_configuration))
         if (frozen) then
            sphere_position(:, 1:nsphere) = sphereposition(:, 1:nsphere)
            sphere_radius(1:nsphere) = sphereradius(1:nsphere)
         else
            if (auto_target_radius .and. target_shape .eq. 2) then
               targetdimensions = target_dimensions + target_radius_padding
!                  nspheresamp=ceiling(sphere_volume_fraction*(targetdimensions(1)-1.d0)**3.d0)
               nspheresamp = ceiling(sphere_volume_fraction * (targetdimensions(1))**3.d0)
               nspheresamp = max(nspheresamp, nsphere)
               if (mstm_global_rank .eq. 0) then
                  write (run_print_unit, '('' set, sampled number spheres:'',2i6)') nsphere, nspheresamp
               end if
            else
               targetdimensions = target_dimensions
               nspheresamp = nsphere
            end if
            if (allocated(sphereradius)) deallocate (sphereradius, sphereposition, sphereindex)
            allocate (sphereradius(nspheresamp), sphereposition(3, nspheresamp), sphereindex(nspheresamp))
            call generate_random_sphere_cluster(nspheresamp, targetdimensions, sphereposition, sphereradius, &
                                                sphereindex, run_print_unit, ran_config_stat, ran_config_time_steps, &
                                                skip_diffusion=skipdif, print_progress=.true.)
            firstrun = .false.
            sphere_position(:, 1:nsphere) = sphereposition(:, 1:nsphere)
            sphere_radius(1:nsphere) = sphereradius(1:nsphere)
            sphere_index(1:nsphere) = sphereindex(1:nsphere)
         end if
      end if
      call synchronize_runtime_status(mpicomm)
      if (runtime_failed()) return
      generation_status(1) = ran_config_stat
      call mstm_mpi(mpi_command='bcast', mpi_send_buf_i=generation_status, mpi_number=1, &
                    mpi_rank=0, mpi_comm=mpicomm)
      ran_config_stat = generation_status(1)
      if (generation_status(1) .ge. 3) then
         call set_runtime_error('Unable to generate random sphere configuration')
         return
      end if
!         call mstm_mpi(mpi_command='barrier')
      nsend = nsphere
      call mstm_mpi(mpi_command='bcast', mpi_send_buf_dp=sphere_radius, &
                    mpi_number=nsend, mpi_rank=0, mpi_comm=mpicomm)
      call mstm_mpi(mpi_command='bcast', mpi_send_buf_i=sphere_index, &
                    mpi_number=nsend, mpi_rank=0, mpi_comm=mpicomm)
      nsend = nsphere * 3
      call mstm_mpi(mpi_command='bcast', mpi_send_buf_dp=sphere_position, &
                    mpi_number=nsend, mpi_rank=0, mpi_comm=mpicomm)
      sphere_radius(1:nsphere) = sphere_radius(1:nsphere) * length_scale_factor
      sphere_position(:, 1:nsphere) = sphere_position(:, 1:nsphere) * length_scale_factor
      do i = 1, nsphere
         sphere_ref_index(:, i) = ref_index_scale_factor &
                                  * component_ref_index(sphere_index(i))
      end do
      if (erase_sphere_1) sphere_ref_index(:, 1) = (1.000001d0, 0.d0)
      if (random_configuration_host) then
         crad = maxval(sqrt(sum(sphere_position(:, 1:nsphere)**2, 1))) + 0.01d0
         if (auto_target_radius .and. target_shape .eq. 2) then
            sphere_radius(number_spheres) &
               = (sum(sphere_radius(1:number_spheres - 1)**3.) / sphere_volume_fraction)**0.33333
         else
            if (random_configuration_host_model .eq. 1) then
               sphere_radius(number_spheres) = target_dimensions(1) * length_scale_factor
            elseif (random_configuration_host_model .eq. 2) then
               sphere_radius(number_spheres) = (sum(sphere_radius(1:number_spheres - 1)**3.) / sphere_volume_fraction)**0.33333
            end if
         end if
         sphere_position(:, number_spheres) = 0.d0
         sphere_ref_index(:, number_spheres) = host_sphere_ref_index
      end if
   end subroutine generate_random_configuration

   subroutine check_sphere_positions()
      implicit none
      logical :: check
      integer :: i, j, imin, jmin
      real(real64) :: r, amax, amin, rmingap

      check = .true.
      rmingap = -1.d10
      do i = 1, number_spheres - 1
         do j = i + 1, number_spheres
            amax = max(sphere_radius(i), sphere_radius(j))
            amin = min(sphere_radius(i), sphere_radius(j))
            r = sqrt(sum((sphere_position(:, i) - sphere_position(:, j))**2))
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
!               write(run_print_unit,'('' spheres '',i4,'' and '',i4,'' intersect'')') i,j
!               write(2,'('' spheres '',i4,'' and '',i4,'' intersect'')') i,j
         end do
      end do
      if (.not. check) then
         write (run_print_unit, '('' warning: sphere-sphere intersections detected, max overlap:'',es12.4, &
            &''  Results might be garbage!'')') rmingap
         write (run_print_unit, '('' positions:'',i5,3es12.4,i5,3es12.4)') imin, sphere_position(:, imin), &
            jmin, sphere_position(:, jmin)
         flush (run_print_unit)
      end if
      check = .true.
      rmingap = -1.d10
      do i = 1, number_spheres
         do j = 1, number_plane_boundaries
            if (abs(sphere_position(3, i) - plane_boundary_position(j)) .ge. sphere_radius(i)) cycle
            check = .false.
            rmingap = max(rmingap, sphere_radius(i) - abs(sphere_position(3, i) - plane_boundary_position(j)))
!               write(run_print_unit,'('' sphere '',i4,'' and plane boundary'',i4,'' intersect'')') i,j
!               write(2,'('' sphere '',i4,'' and plane boundary'',i4,'' intersect'')') i,j
         end do
      end do
      if (.not. check) then
         write (run_print_unit, '('' warning: sphere-plane boundary intersections detected, max overlap:'',es12.4, &
            &''  Results might be garbage!'')') rmingap
!            write(run_print_unit,'('' positions:'',i5,3es12.4,i5,3es12.4)') imin,sphere_position(:,imin), &
!               jmin,sphere_position(:,jmin)
         flush (run_print_unit)
      end if
   end subroutine check_sphere_positions
end module configuration_data
