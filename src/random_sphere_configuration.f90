module random_sphere_configuration
   use, intrinsic :: iso_fortran_env, only: real64
   use constants
   use parallel_runtime, only: mstm_global_rank, parallel_rank, parallel_wall_time
   use random_configuration_dynamics
   use random_configuration_geometry
   use random_configuration_sorting
   use random_configuration_state
   use runtime_support, only: open_output_file, runtime_failed
   use angular_functions, only: euler_rotate_cartesian_vectors
   implicit none
   private
   public :: add_sphere_to_cluster, calculate_target_distribution_statistics, calculate_target_volume, &
             check_position_in_target, circumscribing_sphere, clear_cells, direct_overlap_test, &
             find_next_sphere_collision, find_next_wall_collision, generate_hexagonal_positions, &
             heap_sort_with_tolerance, initialize_cells, modify_cells, move_spheres, position_to_cell_index, &
             generate_random_sphere_cluster, resolve_collision_velocities, sample_layered_configuration, &
             sample_particle_radius, sample_random_velocities, sample_target_position, sort_sphere_positions, &
             sort_sphere_radii, swap_cell_contents, test_pair_collision
   public :: c_list, c_temp, cell_list, coll_data, coll_list, component_number_fraction, component_radii, &
             l_list, max_collisions_per_sphere, max_diffusion_cpu_time, max_diffusion_simulation_time, &
             max_number_time_steps, number_components, periodic_bc, psd_sigma, random_lattice_configuration, &
             random_seed_value, &
             sim_timings, sphere_1_fixed, sphere_cell, target_dimensions, target_shape, target_thickness, &
             target_width, target_width_specified, time_0, wall_boundary_model
contains

   subroutine generate_random_sphere_cluster(numberspheres, targetdimensions, sphereposition, sphereradius, &
                                             sphereindex, iunit, istatus, &
                                             ntsteps, skip_diffusion, use_saved_values, print_progress, simulation_file)
      implicit none
      logical :: fitok, allin, initial0, initial1, trystage1, skipdif, pprog, printsim, multicomp
      logical, optional :: skip_diffusion, use_saved_values, print_progress
      logical, save :: firstrun
      integer :: i, j, maxsamp0, maxsamp1, numberspheres, ncolls, ncollstot, maxns, ntsteps, istatus, iunit, &
                 rank, opair(2), simulation_unit, &
                 nscompi(4), sphereindex(numberspheres), n
      real(real64) :: samppos(3), sphereposition(3, numberspheres), sphereradius(numberspheres), &
                      spherevol, targetfv, u(3, numberspheres), wallboundaries(3, 2), targetvol, targetdimensions(3), &
                      targetstretch, collspersphere, mfp, time0, time1, sum1, sum2, sdev, mean, time2, rmin, rnum(3), &
                      radscale
      real(real64), allocatable, save :: saved_sphereradius(:), saved_sphereposition(:, :)
      character(len=255), optional :: simulation_file
      data maxsamp0, maxsamp1, firstrun/10000, 100, .true./
      if (present(use_saved_values)) then
         if (use_saved_values) then
            sphereposition = saved_sphereposition
            sphereradius = saved_sphereradius
            return
         end if
      end if
      if (present(skip_diffusion)) then
         skipdif = skip_diffusion
      else
         skipdif = .false.
      end if
      if (present(print_progress)) then
         pprog = print_progress
      else
         pprog = .false.
      end if
      printsim = (present(simulation_file) .and. mstm_global_rank .eq. 0)
      if (firstrun) then
         call initialize_random_seed(random_seed_value)
         firstrun = .false.
      end if
      allocate (coll_data(numberspheres))
      call parallel_rank(mpi_rank=rank)
      trystage1 = .false.

      nscompi(1:number_components) = numberspheres * component_number_fraction(1:number_components)
      if (number_components .gt. 1) nscompi(number_components) = nscompi(number_components) &
                                                                 + numberspheres - sum(nscompi(1:number_components))
      radscale = (sum(component_number_fraction(1:number_components) &
                      * component_radii(1:number_components)**3))**(1.d0 / 3.d0)

      spherevol = 0.d0
      n = 0
      do j = 1, number_components
         do i = 1, nscompi(j)
            n = n + 1
            sphereindex(n) = j
            call sample_particle_radius(psd_sigma(j), 2.5d0, sphereradius(n))
            sphereradius(n) = sphereradius(n) * component_radii(j) / radscale
            spherevol = spherevol + four_pi_over_three * sphereradius(n)**3
         end do
         !            if(psd_sigma.gt.0.1d0) then
         !               call sort_sphere_radii(numberspheres,sphereradius)
         !               trystage1=.true.
         !            endif
      end do
      call calculate_target_volume(targetdimensions, targetvol)
      targetfv = spherevol / targetvol
      targetstretch = (1.d0 / targetfv)**(1.d0 / 3.d0)
      targetstretch = max(targetstretch, 1.02d0)
      mfp = targetvol / dble(numberspheres) / 4.d0
      target_boundaries(:, 1) = -targetdimensions
      target_boundaries(:, 2) = targetdimensions
      wallboundaries = target_boundaries
      d_cell = 2.5d0 * maxval(sphereradius(1:numberspheres))
      allin = .false.
      istatus = 3
      sim_timings = 0.d0
      if ((targetfv .le. 0.25d0 .or. trystage1) .and. (.not. allin) .and. (.not. random_lattice_configuration)) then
         allin = .true.
         call initialize_cells(numberspheres)
         do i = 1, numberspheres
            do j = 1, maxsamp0
               call sample_target_position(samppos, sphereradius(i))
               if (sphere_1_fixed .and. i .eq. 1) samppos = 0.d0
               call add_sphere_to_cluster(sphereradius(i), samppos, i - 1, sphereradius, sphereposition, fitok)
               if (fitok) then
                  sphereposition(:, i) = samppos(:)
                  exit
               end if
            end do
            if (j .ge. maxsamp0) then
               allin = .false.
               exit
            end if
         end do
         if (allin) then
            ntsteps = min(ceiling(2.d0 / time_step), max_number_time_steps)
            if (mstm_global_rank .eq. 0 .and. pprog) then
               write (iunit, '('' target configuration computed using random sampling'')')
               flush (iunit)
            end if
            istatus = 0
         else
            call clear_cells()
         end if
      end if
      if (targetfv .lt. 0.6d0 .and. (.not. allin) .and. (.not. random_lattice_configuration)) then
         allin = .true.
         sum1 = 0.
         sum2 = 0.
         do j = 1, maxsamp1
            call initialize_cells(numberspheres)
            call sample_layered_configuration(numberspheres, sphereradius, sphereposition, wallboundaries, maxns)
            if (maxns .ge. numberspheres) exit
            sum1 = sum1 + maxns
            sum2 = sum2 + maxns * maxns
            mean = sum1 / dble(j)
            sdev = sqrt(dble(j) * sum2 - sum1 * sum1) / dble(j)
!if(rank.eq.0) then
!write(*,'(3i10,es12.4)') j,maxns,numberspheres,2.d0*sdev+mean
!flush(6)
!endif
            if (j .gt. 20 .and. 2.d0 * sdev + mean .lt. numberspheres) exit
            call clear_cells()
         end do
         if (maxns .lt. numberspheres) then
            allin = .false.
         else
            ntsteps = min(ceiling(mfp / time_step), max_number_time_steps)
            if (mstm_global_rank .eq. 0 .and. pprog) then
               write (iunit, '('' target configuration computed using layered sampling + diffusion, time steps:'',i5)') ntsteps
               flush (iunit)
            end if
            istatus = 1
         end if
      end if
      if ((.not. allin) .or. random_lattice_configuration) then
         do
            call initialize_cells(numberspheres)
            call generate_hexagonal_positions(numberspheres, sphereradius, sphereposition, &
                                              wallboundaries, targetstretch, allin, maxns)
!if(rank.eq.0) then
!write(*,'(i10,es12.4)') maxns,targetstretch
!flush(6)
!endif

            if (allin) exit
            call clear_cells()
            if (targetstretch .le. 1.02d0) then
               write (iunit, '('' MC configuration sampler failed'')')
               istatus = 3
               return
            end if
            targetstretch = targetstretch - 0.001
            targetstretch = max(targetstretch, 1.02d0)
         end do
         istatus = 2
         ntsteps = max_number_time_steps
         if (mstm_global_rank .eq. 0 .and. pprog) then
            write (iunit, '('' target configuration computed initial HCP + diffusion, time steps:'',i5)') ntsteps
            flush (iunit)
         end if
      end if
      do i = 1, numberspheres
         call check_position_in_target(sphereradius(i), sphereposition(:, i), wallboundaries, allin)
         if (.not. allin) write (iunit, '('' initially outside:'',i5,3es12.4)') i, sphereposition(:, i)
      end do
!call direct_overlap_test(numberspheres,sphereradius,sphereposition,allin,distance=rmin,pair=opair)
!if(allin) write(iunit,'('' initial overlap:'',2i4,f8.3)') opair,rmin
      ntsteps = max_number_time_steps
      if (printsim) then
         call open_output_file(trim(simulation_file), simulation_unit)
         printsim = .not. runtime_failed()
         if (printsim) then
            write (simulation_unit, '(i8)') numberspheres
            write (simulation_unit, '(es13.5)') 0.d0
            do i = 1, numberspheres
               write (simulation_unit, '(4es13.5)') sphereposition(:, i), sphereradius(i)
            end do
         end if
      end if
      if ((.not. skipdif) .and. max_number_time_steps .gt. 0 .and. max_diffusion_simulation_time .gt. 0.) then
         call sample_random_velocities(numberspheres, u)
         ncollstot = 0
         time1 = parallel_wall_time()
         time0 = time1
         do j = 1, max_number_time_steps
            call move_spheres(numberspheres, sphereradius, sphereposition, u, time_step, wallboundaries, &
                              number_wall_hits=ncolls)
            ncollstot = ncollstot + ncolls
            collspersphere = dble(ncollstot) / dble(numberspheres)
            time2 = parallel_wall_time()
            if (mstm_global_rank .eq. 0 .and. pprog .and. time2 - time1 .gt. 15.d0) then
               write (iunit, '('' diffusion step, collision per sphere:'',i8,es12.4)') j, collspersphere
               flush (iunit)
               time1 = parallel_wall_time()
            end if
            if (time2 - time0 .gt. max_diffusion_cpu_time) exit
            if (j * time_step .gt. max_diffusion_simulation_time) exit
            if (collspersphere .gt. max_collisions_per_sphere) exit
            if (printsim) then
               write (simulation_unit, '(es13.5)') j * time_step
               do i = 1, numberspheres
                  write (simulation_unit, '(4es13.5)') sphereposition(:, i), sphereradius(i)
               end do
            end if
         end do
         ntsteps = min(ntsteps, j)
         do i = 1, numberspheres
            call check_position_in_target(sphereradius(i), sphereposition(:, i), wallboundaries, allin)
            if (.not. allin) write (iunit, '('' outside:'',i5,3es12.4)') i, sphereposition(:, i)
         end do
!call direct_overlap_test(numberspheres,sphereradius,sphereposition,allin,distance=rmin,pair=opair)
!if(allin) write(iunit,'('' overlap:'',2i4,f8.3)') opair,rmin
      end if
      if (printsim) close (simulation_unit)
      if (target_shape .eq. 0 .or. target_shape .eq. 1) then
         call sort_sphere_positions(numberspheres, sphereradius, sphereposition, sphereindex, 3)
      else
         call sort_sphere_positions(numberspheres, sphereradius, sphereposition, sphereindex, 0)
      end if
      if (random_lattice_configuration) then
         call random_number(rnum(1:3))
         rnum(1) = two_pi * rnum(1)
         rnum(3) = two_pi * rnum(3)
         rnum(2) = dacos(-1.d0 + 2.d0 * rnum(2))
         call euler_rotate_cartesian_vectors(sphereposition(:, 1:numberspheres), rnum, 1, &
                                             sphereposition(:, 1:numberspheres), numberspheres)
      end if

      sphereposition(:, 1:numberspheres) = sphereposition(:, 1:numberspheres) * radscale
      sphereradius(1:numberspheres) = sphereradius(1:numberspheres) * radscale

      if (allocated(saved_sphereposition)) then
         deallocate (saved_sphereposition, saved_sphereradius)
      end if
      allocate (saved_sphereposition(3, numberspheres), saved_sphereradius(numberspheres))
      saved_sphereposition = sphereposition
      saved_sphereradius = sphereradius
      call clear_cells()
      deallocate (coll_data)
   end subroutine generate_random_sphere_cluster

   subroutine initialize_random_seed(seed_value)
      integer, intent(in) :: seed_value
      integer :: i, seed_size
      integer, allocatable :: seed(:)

      if (seed_value < 0) then
         call random_seed()
         return
      end if

      call random_seed(size=seed_size)
      allocate (seed(seed_size))
      seed = seed_value + 104729*[(i - 1, i=1, seed_size)]
      call random_seed(put=seed)
   end subroutine initialize_random_seed
end module random_sphere_configuration
