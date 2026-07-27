program output_regression_test
   use, intrinsic :: iso_fortran_env, only: error_unit, real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   character(len=64) :: case_name
   character(len=1024) :: comparison_directory_1, comparison_directory_2, work_directory
   integer :: failures

   call get_command_argument(1, case_name)
   call get_command_argument(2, work_directory)
   failures = 0

   select case (trim(case_name))
   case ('figure1')
      call check_figure1(trim(work_directory), failures)
   case ('effective_medium')
      call check_effective_medium(trim(work_directory), failures)
   case ('two_sphere_broadside')
      call check_two_sphere_broadside(trim(work_directory), failures)
   case ('rigid_transform_invariants')
      call check_rigid_transform_invariants(trim(work_directory), failures)
   case ('solver_validation')
      call check_solver_validation(trim(work_directory), failures)
   case ('solver_reciprocity')
      call check_solver_reciprocity(trim(work_directory), failures)
   case ('fft_translation_validation')
      call check_fft_translation_validation(trim(work_directory), failures)
   case ('random_orientation_average')
      call check_random_orientation_average(trim(work_directory), failures)
   case ('incidence_average')
      call check_incidence_average(trim(work_directory), failures)
   case ('periodic_surface')
      call check_periodic_surface(trim(work_directory), failures)
   case ('nested_sphere')
      call check_nested_sphere(trim(work_directory), failures)
   case ('near_field_modes')
      call check_near_field_modes(trim(work_directory), failures)
   case ('state_lifecycle')
      call check_state_lifecycle(trim(work_directory), failures)
   case ('parallel_solver_equivalence')
      call get_command_argument(3, comparison_directory_1)
      call get_command_argument(4, comparison_directory_2)
      call check_parallel_equivalence(trim(work_directory), trim(comparison_directory_1), &
                                      trim(comparison_directory_2), 'solver-validation.dat', 3, failures)
   case ('parallel_fft_equivalence')
      call get_command_argument(3, comparison_directory_1)
      call get_command_argument(4, comparison_directory_2)
      call check_parallel_equivalence(trim(work_directory), trim(comparison_directory_1), &
                                      trim(comparison_directory_2), 'fft-translation-validation.dat', 2, failures)
   case default
      write (error_unit, '(a)') 'Unknown regression case: '//trim(case_name)
      error stop 2
   end select

   if (failures /= 0) then
      write (error_unit, '(i0,a)') failures, ' numerical regression checks failed'
      error stop 1
   end if

contains

   subroutine check_figure1(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: efficiency(9), near_field_row(29)
      integer :: unit, i

      call open_regression_file(directory//'/mstm-2022b-fig1.dat', unit)
      call find_line(unit, 'total extinction, absorption, scattering efficiencies', .true.)
      read (unit, *) efficiency
      close (unit)

      call assert_close('Figure 1 total extinction', efficiency(1), 2.7566_real64, failure_count)
      call assert_close('Figure 1 total scattering', efficiency(3), 2.7566_real64, failure_count)
      call assert_close('Figure 1 parallel extinction', efficiency(4), 2.7566_real64, failure_count)

      call open_regression_file(directory//'/nf-fig1.dat', unit)
      do i = 1, 11
         read (unit, *)
      end do
      read (unit, *) near_field_row
      close (unit)

      call assert_close('Figure 1 near-field x', near_field_row(1), -9.95_real64, failure_count)
      call assert_close('Figure 1 first electric-field component', near_field_row(4), 0.48809_real64, failure_count)
      call assert_close('Figure 1 second electric-field component', near_field_row(5), 0.76756_real64, failure_count)
   end subroutine check_figure1

   subroutine check_effective_medium(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: efficiency(9), fit(4), coefficient_real, coefficient_imaginary
      integer :: unit, fit_status, order, coefficient_order, polarization

      call open_regression_file(directory//'/effective-medium-smoke.dat', unit)
      call find_line(unit, 'total extinction, absorption, scattering efficiencies', .true.)
      read (unit, *) efficiency
      call find_line(unit, 'mie fit effective ri, radius, RT ratio, status', .true.)
      read (unit, *) fit, fit_status
      close (unit)

      call assert_close('Effective-medium total extinction', efficiency(1), 6.0387e-3_real64, failure_count)
      call assert_close('Effective-medium total absorption', efficiency(2), 5.1562e-3_real64, failure_count)
      call assert_close('Effective-medium refractive-index real part', fit(1), 1.0680_real64, failure_count)
      call assert_close('Effective-medium refractive-index imaginary part', fit(2), -1.4990e-2_real64, failure_count)
      call assert_close('Effective-medium fitted radius', fit(3), 4.6107e-1_real64, failure_count)
      if (fit_status /= 1) then
         write (error_unit, '(a,i0,a)') 'Effective-medium fit status was ', fit_status, ', expected 1'
         failure_count = failure_count + 1
      end if

      call open_regression_file(directory//'/anpeff.dat', unit)
      read (unit, *) order
      read (unit, *) coefficient_order, polarization, coefficient_real, coefficient_imaginary
      close (unit)
      if (order /= 5 .or. coefficient_order /= 1 .or. polarization /= 1) then
         write (error_unit, '(a)') 'Effective-medium coefficient metadata changed'
         failure_count = failure_count + 1
      end if
      call assert_close('Effective-medium first coefficient real part', coefficient_real, &
                        6.76874e-4_real64, failure_count)
      call assert_close('Effective-medium first coefficient imaginary part', coefficient_imaginary, &
                        2.92138e-3_real64, failure_count)
   end subroutine check_effective_medium

   subroutine check_two_sphere_broadside(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: efficiency(9)
      integer :: unit

      call open_regression_file(directory//'/two-sphere-broadside.dat', unit)
      call find_line(unit, 'total extinction, absorption, scattering efficiencies', .true.)
      read (unit, *) efficiency
      close (unit)

      call assert_close('Broadside total extinction', efficiency(1), 3.6038e-2_real64, failure_count)
      call assert_close('Broadside parallel extinction', efficiency(4), 4.1080e-2_real64, failure_count)
      call assert_close('Broadside perpendicular extinction', efficiency(7), 3.0996e-2_real64, failure_count)
   end subroutine check_two_sphere_broadside

   subroutine check_rigid_transform_invariants(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: efficiency(9, 4)
      integer :: component, run, unit
      character(len=64) :: label

      call open_regression_file(directory//'/rigid-transform-invariants.dat', unit)
      do run = 1, 4
         call find_line(unit, 'total extinction, absorption, scattering efficiencies', .true.)
         read (unit, *) efficiency(:, run)
      end do
      close (unit)

      do component = 1, 9
         write (label, '(a,i0)') 'Translated efficiency component ', component
         call assert_close(trim(label), efficiency(component, 2), &
                           efficiency(component, 1), failure_count)
         write (label, '(a,i0)') 'Z-rotated efficiency component ', component
         call assert_close(trim(label), efficiency(component, 3), &
                           efficiency(component, 1), failure_count)
      end do
      do component = 1, 3
         write (label, '(a,i0)') 'Y-rotated unpolarized efficiency component ', component
         call assert_close(trim(label), efficiency(component, 4), &
                           efficiency(component, 1), failure_count)
      end do
   end subroutine check_rigid_transform_invariants

   subroutine check_solver_validation(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: efficiency(9, 3), residual(3), elapsed, reciprocal_condition
      integer :: iterations(3), run, unit

      call open_regression_file(directory//'/solver-validation.dat', unit)
      do run = 1, 3
         call find_line(unit, 'number iterations, error, solution time', .true.)
         read (unit, *) iterations(run), residual(run), elapsed
         if (run == 3) then
            call find_line(unit, 'direct reciprocal condition estimate', .true.)
            read (unit, *) reciprocal_condition
         end if
         call find_line(unit, 'total extinction, absorption, scattering efficiencies', .true.)
         read (unit, *) efficiency(:, run)
      end do
      close (unit)

      if (residual(1) > 1.0e-4_real64) then
         write (error_unit, '(a,es14.6)') 'Loose iterative residual exceeded tolerance: ', residual(1)
         failure_count = failure_count + 1
      end if
      if (residual(2) > 1.0e-10_real64) then
         write (error_unit, '(a,es14.6)') 'Tight iterative residual exceeded tolerance: ', residual(2)
         failure_count = failure_count + 1
      end if
      if (residual(3) > 1.0e-10_real64) then
         write (error_unit, '(a,es14.6)') 'Direct residual exceeded tolerance: ', residual(3)
         failure_count = failure_count + 1
      end if
      if (.not. ieee_is_finite(reciprocal_condition) .or. reciprocal_condition <= 0.0_real64 .or. &
          reciprocal_condition > 1.0_real64) then
         write (error_unit, '(a,es14.6)') 'Invalid direct reciprocal condition estimate: ', reciprocal_condition
         failure_count = failure_count + 1
      end if
      if (iterations(2) < iterations(1) .or. residual(2) > residual(1)) then
         write (error_unit, '(a)') 'Tighter tolerance did not improve iterative convergence'
         failure_count = failure_count + 1
      end if
      do run = 1, 9
         call assert_close('Direct/iterative efficiency', efficiency(run, 3), efficiency(run, 2), failure_count)
      end do
      do run = 1, 3
         call assert_close('Lossless unpolarized energy balance', efficiency(3, run), &
                           efficiency(1, run), failure_count)
         call assert_close('Lossless parallel energy balance', efficiency(6, run), &
                           efficiency(4, run), failure_count)
         call assert_close('Lossless perpendicular energy balance', efficiency(9, run), &
                           efficiency(7, run), failure_count)
      end do
   end subroutine check_solver_validation

   subroutine check_solver_reciprocity(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: direct_matrix(16), reciprocal_matrix_left(16), reciprocal_matrix_right(16)
      integer :: unit

      call open_regression_file(directory//'/solver-reciprocity.dat', unit)
      call find_line(unit, 'scattering matrix in incident plane', .true.)
      call find_scattering_row(unit, 60.0_real64, direct_matrix)
      call find_line(unit, 'scattering matrix in incident plane', .true.)
      call find_scattering_row(unit, -180.0_real64, reciprocal_matrix_left)
      call find_scattering_row(unit, 180.0_real64, reciprocal_matrix_right)
      close (unit)

      call assert_close('Reciprocal unpolarized differential scattering', &
                        0.5_real64 * (reciprocal_matrix_left(1) + reciprocal_matrix_right(1)), &
                        direct_matrix(1), failure_count)
   end subroutine check_solver_reciprocity

   subroutine check_fft_translation_validation(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: efficiency(9, 2)
      integer :: component, interaction_count, run, transform_count, unit

      call open_regression_file(directory//'/fft-translation-validation.dat', unit)
      do run = 1, 2
         call find_line(unit, 'total extinction, absorption, scattering efficiencies', .true.)
         read (unit, *) efficiency(:, run)
      end do
      rewind (unit)
      call find_line(unit, 'FFT interactions, 3-D transforms', .true.)
      read (unit, *) interaction_count, transform_count
      close (unit)

      if (interaction_count <= 0 .or. transform_count <= 0) then
         write (error_unit, '(a,2i8)') 'FFT path did not report work: ', interaction_count, transform_count
         failure_count = failure_count + 1
      end if
      do component = 1, 9
         call assert_close('FFT/pairwise efficiency', efficiency(component, 2), &
                           efficiency(component, 1), failure_count)
      end do
   end subroutine check_fft_translation_validation

   subroutine check_random_orientation_average(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: efficiency(3)
      integer :: unit

      call open_regression_file(directory//'/random-orientation-average.dat', unit)
      call find_line(unit, 'total extinction, absorption, scattering efficiencies', .true.)
      read (unit, *) efficiency
      close (unit)

      call require_finite('Random-orientation efficiencies', efficiency, failure_count)
      call assert_close('Random-orientation extinction', efficiency(1), 8.8930e-4_real64, failure_count)
      call assert_close('Random-orientation lossless energy balance', efficiency(1), &
                        efficiency(2) + efficiency(3), failure_count)
   end subroutine check_random_orientation_average

   subroutine check_incidence_average(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: efficiency(3)
      integer :: unit

      call open_regression_file(directory//'/incidence-average.dat', unit)
      call find_line(unit, 'total extinction, absorption, scattering efficiencies', .true.)
      read (unit, *) efficiency
      close (unit)

      call require_finite('Incidence-average efficiencies', efficiency, failure_count)
      call assert_close('Incidence-average analytical sphere extinction', efficiency(1), &
                        2.1510e-1_real64, failure_count)
      call assert_close('Incidence-average lossless energy balance', efficiency(1), &
                        efficiency(2) + efficiency(3), failure_count)
   end subroutine check_incidence_average

   subroutine check_periodic_surface(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: flux(9)
      integer :: polarization, unit

      call open_regression_file(directory//'/periodic-surface.dat', unit)
      call find_line(unit, 'unit cell reflectance, absorptance, transmittance', .true.)
      read (unit, *) flux
      close (unit)

      call require_finite('Periodic-surface fluxes', flux, failure_count)
      call assert_close('Periodic-surface unpolarized reflectance', flux(1), 7.8075e-3_real64, failure_count)
      do polarization = 0, 2
         call assert_close('Periodic-surface energy balance', &
                           sum(flux(3 * polarization + 1:3 * polarization + 3)), &
                           1.0_real64, failure_count)
      end do
   end subroutine check_periodic_surface

   subroutine check_nested_sphere(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: efficiency(9)
      integer :: polarization, unit

      call open_regression_file(directory//'/nested-sphere.dat', unit)
      call find_line(unit, 'total extinction, absorption, scattering efficiencies', .true.)
      read (unit, *) efficiency
      close (unit)

      call require_finite('Nested-sphere efficiencies', efficiency, failure_count)
      call assert_close('Nested-sphere unpolarized extinction', efficiency(1), &
                        1.8049e-1_real64, failure_count)
      do polarization = 0, 2
         call assert_close('Nested-sphere lossless energy balance', efficiency(3 * polarization + 1), &
                           efficiency(3 * polarization + 2) + efficiency(3 * polarization + 3), &
                           failure_count)
      end do
   end subroutine check_nested_sphere

   subroutine check_near_field_modes(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: scattered_field(29), total_field(29)
      integer :: unit

      call open_regression_file(directory//'/near-field-total.dat', unit)
      call skip_lines(unit, 7)
      read (unit, *) total_field
      close (unit)
      call open_regression_file(directory//'/near-field-scattered.dat', unit)
      call skip_lines(unit, 7)
      read (unit, *) scattered_field
      close (unit)

      call require_finite('Total near field', total_field, failure_count)
      call require_finite('Scattered near field', scattered_field, failure_count)
      call assert_close('Near-field sample x coordinate', total_field(1), 2.05_real64, failure_count)
      call assert_close('Near-field electric incident contribution', &
                        total_field(4) - scattered_field(4), 1.0_real64, failure_count)
      call assert_close('Near-field magnetic incident contribution', &
                        total_field(22) - scattered_field(22), -1.0_real64, failure_count)
   end subroutine check_near_field_modes

   subroutine check_state_lifecycle(directory, failure_count)
      character(len=*), intent(in) :: directory
      integer, intent(inout) :: failure_count
      real(real64) :: efficiency(9, 5)
      integer :: component, run, unit

      call open_regression_file(directory//'/state-lifecycle.dat', unit)
      do run = 1, 5
         call find_line(unit, 'total extinction, absorption, scattering efficiencies', .true.)
         read (unit, *) efficiency(:, run)
      end do
      close (unit)

      call require_finite('Repeated-run efficiencies', reshape(efficiency, [size(efficiency)]), failure_count)
      call assert_close('Resized single-sphere analytical result', efficiency(1, 4), &
                        2.1510e-1_real64, failure_count)
      do component = 1, 9
         call assert_close('FFT state transition', efficiency(component, 2), &
                           efficiency(component, 1), failure_count)
         call assert_close('Direct-solver state transition', efficiency(component, 3), &
                           efficiency(component, 1), failure_count)
         call assert_close('Restored-problem state transition', efficiency(component, 5), &
                           efficiency(component, 1), failure_count)
      end do
   end subroutine check_state_lifecycle

   subroutine check_parallel_equivalence(serial_directory, two_rank_directory, four_rank_directory, &
                                         output_name, number_runs, failure_count)
      character(len=*), intent(in) :: serial_directory, two_rank_directory, four_rank_directory, output_name
      integer, intent(in) :: number_runs
      integer, intent(inout) :: failure_count
      real(real64), allocatable :: four_rank_efficiency(:, :), serial_efficiency(:, :), two_rank_efficiency(:, :)
      integer :: component, run

      allocate (serial_efficiency(9, number_runs), two_rank_efficiency(9, number_runs), &
                four_rank_efficiency(9, number_runs))
      call read_efficiency_runs(serial_directory//'/'//output_name, serial_efficiency)
      call read_efficiency_runs(two_rank_directory//'/'//output_name, two_rank_efficiency)
      call read_efficiency_runs(four_rank_directory//'/'//output_name, four_rank_efficiency)

      do run = 1, number_runs
         do component = 1, 9
            call assert_parallel_close('Serial/two-rank efficiency', two_rank_efficiency(component, run), &
                                       serial_efficiency(component, run), failure_count)
            call assert_parallel_close('Serial/four-rank efficiency', four_rank_efficiency(component, run), &
                                       serial_efficiency(component, run), failure_count)
         end do
      end do
   end subroutine check_parallel_equivalence

   subroutine read_efficiency_runs(path, efficiency)
      character(len=*), intent(in) :: path
      real(real64), intent(out) :: efficiency(:, :)
      integer :: run, unit

      call open_regression_file(path, unit)
      do run = 1, size(efficiency, 2)
         call find_line(unit, 'total extinction, absorption, scattering efficiencies', .true.)
         read (unit, *) efficiency(:, run)
      end do
      close (unit)
   end subroutine read_efficiency_runs

   subroutine open_regression_file(path, unit)
      character(len=*), intent(in) :: path
      integer, intent(out) :: unit
      integer :: io_status
      character(len=512) :: io_message

      open (newunit=unit, file=path, status='old', action='read', iostat=io_status, iomsg=io_message)
      if (io_status /= 0) then
         write (error_unit, '(3a)') 'Cannot open ', trim(path), ': '//trim(io_message)
         error stop 2
      end if
   end subroutine open_regression_file

   subroutine find_line(unit, text, required)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: text
      logical, intent(in) :: required
      character(len=2048) :: line
      integer :: io_status

      do
         read (unit, '(a)', iostat=io_status) line
         if (io_status /= 0) exit
         if (index(line, text) /= 0) return
      end do
      if (required) then
         write (error_unit, '(3a)') 'Could not find "', trim(text), '" in regression output'
         error stop 2
      end if
   end subroutine find_line

   subroutine find_scattering_row(unit, requested_angle, matrix)
      integer, intent(in) :: unit
      real(real64), intent(in) :: requested_angle
      real(real64), intent(out) :: matrix(16)
      character(len=2048) :: line
      integer :: io_status
      real(real64) :: angle

      do
         read (unit, '(a)', iostat=io_status) line
         if (io_status /= 0) exit
         read (line, *, iostat=io_status) angle, matrix
         if (io_status == 0) then
            if (abs(angle - requested_angle) < 1.0e-6_real64) return
         end if
      end do
      write (error_unit, '(a,f8.2)') 'Could not find scattering angle ', requested_angle
      error stop 2
   end subroutine find_scattering_row

   subroutine skip_lines(unit, number_lines)
      integer, intent(in) :: unit, number_lines
      integer :: line

      do line = 1, number_lines
         read (unit, *)
      end do
   end subroutine skip_lines

   subroutine require_finite(label, values, failure_count)
      character(len=*), intent(in) :: label
      real(real64), intent(in) :: values(:)
      integer, intent(inout) :: failure_count

      if (.not. all(ieee_is_finite(values))) then
         write (error_unit, '(a)') trim(label)//' contain a non-finite value'
         failure_count = failure_count + 1
      end if
   end subroutine require_finite

   subroutine assert_close(label, actual, expected, failure_count)
      character(len=*), intent(in) :: label
      real(real64), intent(in) :: actual, expected
      integer, intent(inout) :: failure_count
      real(real64), parameter :: absolute_tolerance = 2.0e-7_real64
      real(real64), parameter :: relative_tolerance = 2.0e-4_real64

      if (.not. ieee_is_finite(actual)) then
         write (error_unit, '(a,es14.6)') trim(label)//' is not finite: actual=', actual
         failure_count = failure_count + 1
      elseif (abs(actual - expected) > absolute_tolerance + relative_tolerance * abs(expected)) then
         write (error_unit, '(a,2(a,es14.6))') trim(label)//' changed:', ' actual=', actual, ' expected=', expected
         failure_count = failure_count + 1
      end if
   end subroutine assert_close

   subroutine assert_parallel_close(label, actual, expected, failure_count)
      character(len=*), intent(in) :: label
      real(real64), intent(in) :: actual, expected
      integer, intent(inout) :: failure_count
      real(real64), parameter :: absolute_tolerance = 2.0e-10_real64
      real(real64), parameter :: relative_tolerance = 2.0e-8_real64

      if (.not. ieee_is_finite(actual)) then
         write (error_unit, '(a,es14.6)') trim(label)//' is not finite: actual=', actual
         failure_count = failure_count + 1
      elseif (abs(actual - expected) > absolute_tolerance + relative_tolerance * abs(expected)) then
         write (error_unit, '(a,2(a,es14.6))') trim(label)//' differs:', &
            ' actual=', actual, ' expected=', expected
         failure_count = failure_count + 1
      end if
   end subroutine assert_parallel_close

end program output_regression_test
