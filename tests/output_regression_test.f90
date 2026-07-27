program output_regression_test
   use, intrinsic :: iso_fortran_env, only: error_unit, real64
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

      call assert_close('Effective-medium total extinction', efficiency(1), 6.0097e-3_real64, failure_count)
      call assert_close('Effective-medium total absorption', efficiency(2), 5.1521e-3_real64, failure_count)
      call assert_close('Effective-medium refractive-index real part', fit(1), 1.0680_real64, failure_count)
      call assert_close('Effective-medium refractive-index imaginary part', fit(2), -1.4983e-2_real64, failure_count)
      call assert_close('Effective-medium fitted radius', fit(3), 4.6115e-1_real64, failure_count)
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
                        6.76879e-4_real64, failure_count)
      call assert_close('Effective-medium first coefficient imaginary part', coefficient_imaginary, &
                        2.92016e-3_real64, failure_count)
   end subroutine check_effective_medium

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

   subroutine assert_close(label, actual, expected, failure_count)
      character(len=*), intent(in) :: label
      real(real64), intent(in) :: actual, expected
      integer, intent(inout) :: failure_count
      real(real64), parameter :: absolute_tolerance = 2.0e-7_real64
      real(real64), parameter :: relative_tolerance = 2.0e-4_real64

      if (abs(actual - expected) > absolute_tolerance + relative_tolerance * abs(expected)) then
         write (error_unit, '(a,2(a,es14.6))') trim(label)//' changed:', ' actual=', actual, ' expected=', expected
         failure_count = failure_count + 1
      end if
   end subroutine assert_close

end program output_regression_test
