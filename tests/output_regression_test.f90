program output_regression_test
   use, intrinsic :: iso_fortran_env, only: error_unit, real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   character(len=64) :: case_name
   character(len=1024) :: work_directory
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

end program output_regression_test
