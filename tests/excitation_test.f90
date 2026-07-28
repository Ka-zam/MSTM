program excitation_test
   use, intrinsic :: iso_fortran_env, only: error_unit, real64
   use constants, only: imaginary_unit
   use excitation, only: electric_dipole_field, electric_dipole_mode_coefficients, &
                         electric_dipole_source_power
   use translation_operator, only: translation_operator_state
   use wave_functions, only: vector_spherical_wave_functions
   implicit none(type, external)

   integer, parameter :: target_order = 8
   integer, parameter :: target_block = target_order * (target_order + 2)
   integer :: failures
   real(real64) :: source_position(3), target_position(3), local_position(3), global_position(3)
   complex(real64) :: moment(3), scaled_moment(3), coefficients(0:2, 1, 2)
   complex(real64) :: regular_origin(3, 6), source_coefficients(6), target_coefficients(2 * target_block)
   complex(real64) :: regular_target(3, 2 * target_block)
   complex(real64) :: electric_direct(3), magnetic_direct(3), electric_translated(3), magnetic_translated(3)
   complex(real64) :: electric_scaled(3), magnetic_scaled(3), electric_rotated(3), magnetic_rotated(3)
   type(translation_operator_state) :: translation

   failures = 0
   moment = [(1.2_real64, -0.3_real64), (-0.4_real64, 0.7_real64), (0.8_real64, 0.2_real64)]
   call electric_dipole_mode_coefficients(moment, coefficients)
   call vector_spherical_wave_functions([0.0_real64, 0.0_real64, 0.0_real64], &
                                        [(1.0_real64, 0.0_real64), (1.0_real64, 0.0_real64)], &
                                        1, 1, regular_origin, index_model=2, lr_to_mode=.true.)
   call assert_close('Cartesian-to-mode normalization', matmul(regular_origin(:, 1:3), coefficients(:, 1, 1)), &
                     moment, 2.0e-14_real64, failures)
   call assert_close_real('source modal power', electric_dipole_source_power(moment), &
                          sum(abs(reshape(coefficients, [6]))**2), 2.0e-14_real64, failures)

   source_position = [0.0_real64, 0.0_real64, 0.0_real64]
   target_position = [0.0_real64, 0.0_real64, 3.0_real64]
   local_position = [0.07_real64, -0.04_real64, 0.09_real64]
   global_position = target_position + local_position
   source_coefficients(1:3) = 0.5_real64 * coefficients(:, 1, 1)
   source_coefficients(4:6) = source_coefficients(1:3)
   target_coefficients = (0.0_real64, 0.0_real64)
   call translation%configure(3, target_position - source_position, &
                              [(1.0_real64, 0.0_real64), (1.0_real64, 0.0_real64)], .false.)
   call translation%apply(1, 2, target_order, 2, source_coefficients, target_coefficients)
   call vector_spherical_wave_functions(local_position, &
                                        [(1.0_real64, 0.0_real64), (1.0_real64, 0.0_real64)], &
                                        target_order, 1, regular_target, index_model=2)
   electric_translated = matmul(regular_target, target_coefficients)
   magnetic_translated = (matmul(regular_target(:, 1:target_block), target_coefficients(1:target_block)) &
                          - matmul(regular_target(:, target_block + 1:), target_coefficients(target_block + 1:))) &
                         / imaginary_unit
   call electric_dipole_field(global_position - source_position, (1.0_real64, 0.0_real64), moment, &
                              electric_direct, magnetic_direct)
   call assert_close('outgoing-to-regular electric translation', electric_translated, electric_direct, &
                     2.0e-10_real64, failures)
   call assert_close('outgoing-to-regular magnetic translation', magnetic_translated, magnetic_direct, &
                     2.0e-10_real64, failures)

   scaled_moment = (2.0_real64, -0.5_real64) * moment
   call electric_dipole_field(global_position, (1.0_real64, 0.0_real64), scaled_moment, &
                              electric_scaled, magnetic_scaled)
   call assert_close('complex-amplitude electric scaling', electric_scaled, &
                     (2.0_real64, -0.5_real64) * electric_direct, 2.0e-14_real64, failures)
   call assert_close('complex-amplitude magnetic scaling', magnetic_scaled, &
                     (2.0_real64, -0.5_real64) * magnetic_direct, 2.0e-14_real64, failures)

   call electric_dipole_field([3.0_real64, 0.0_real64, 0.0_real64], (1.0_real64, 0.0_real64), &
                              [(0.0_real64, 0.0_real64), (0.0_real64, 0.0_real64), (1.0_real64, 0.0_real64)], &
                              electric_rotated, magnetic_rotated)
   call electric_dipole_field([0.0_real64, 0.0_real64, 3.0_real64], (1.0_real64, 0.0_real64), &
                              [(1.0_real64, 0.0_real64), (0.0_real64, 0.0_real64), (0.0_real64, 0.0_real64)], &
                              electric_direct, magnetic_direct)
   call assert_close_real('rotated electric-field norm', sum(abs(electric_rotated)**2), &
                          sum(abs(electric_direct)**2), 2.0e-14_real64, failures)
   call assert_close_real('rotated magnetic-field norm', sum(abs(magnetic_rotated)**2), &
                          sum(abs(magnetic_direct)**2), 2.0e-14_real64, failures)

   if (failures > 0) then
      write (error_unit, '(i0,a)') failures, ' excitation checks failed'
      error stop 1
   end if

contains

   subroutine assert_close(label, actual, expected, tolerance, failure_count)
      character(len=*), intent(in) :: label
      complex(real64), intent(in) :: actual(:), expected(:)
      real(real64), intent(in) :: tolerance
      integer, intent(inout) :: failure_count

      if (maxval(abs(actual - expected)) > tolerance * max(1.0_real64, maxval(abs(expected)))) then
         write (error_unit, '(a,es12.4)') trim(label)//' error: ', maxval(abs(actual - expected))
         failure_count = failure_count + 1
      end if
   end subroutine assert_close

   subroutine assert_close_real(label, actual, expected, tolerance, failure_count)
      character(len=*), intent(in) :: label
      real(real64), intent(in) :: actual, expected, tolerance
      integer, intent(inout) :: failure_count

      if (abs(actual - expected) > tolerance * max(1.0_real64, abs(expected))) then
         write (error_unit, '(a,2es12.4)') trim(label)//' actual/expected: ', actual, expected
         failure_count = failure_count + 1
      end if
   end subroutine assert_close_real

end program excitation_test
