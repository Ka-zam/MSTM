program mie_test
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: real64
   use mie, only: optically_active_mie_coefficients
   implicit none(type, external)

   integer, parameter :: maximum_order = 32
   real(real64), parameter :: size_parameter = 1.0_real64
   real(real64), parameter :: reference_efficiency = 0.2150975960428854_real64
   complex(real64) :: coefficients(2, 2, maximum_order)
   complex(real64) :: medium_index(2), sphere_index(2)
   integer :: order
   real(real64) :: absorption, extinction, scattering
   real(real64) :: coefficient_extinction, coefficient_scattering

   medium_index = cmplx(1.0_real64, 0.0_real64, kind=real64)
   sphere_index = cmplx(1.5_real64, 0.0_real64, kind=real64)
   coefficients = cmplx(0.0_real64, 0.0_real64, kind=real64)
   call optically_active_mie_coefficients(size_parameter, sphere_index, order, 1.0e-13_real64, &
                                          extinction, scattering, absorption, &
                                          anp_mie=coefficients, ri_medium=medium_index)

   call require_close(extinction, reference_efficiency, 2.0e-12_real64, &
                      'single-sphere extinction differs from the analytical Mie series')
   call require_close(scattering, reference_efficiency, 2.0e-12_real64, &
                      'single-sphere scattering differs from the analytical Mie series')
   call require_close(absorption, 0.0_real64, 2.0e-13_real64, &
                      'lossless Mie sphere does not conserve energy')

   call efficiencies_from_coefficients(coefficients(:, :, :order), coefficient_extinction, &
                                       coefficient_scattering)
   call require_close(extinction, coefficient_extinction, 2.0e-13_real64, &
                      'Mie coefficients violate the forward optical theorem')
   call require_close(scattering, coefficient_scattering, 2.0e-13_real64, &
                      'Mie coefficient norm disagrees with scattering efficiency')

   sphere_index = medium_index
   call optically_active_mie_coefficients(size_parameter, sphere_index, order, 1.0e-13_real64, &
                                          extinction, scattering, absorption, ri_medium=medium_index)
   call require_close(extinction, 0.0_real64, 2.0e-13_real64, &
                      'zero-contrast sphere has nonzero extinction')
   call require_close(scattering, 0.0_real64, 2.0e-13_real64, &
                      'zero-contrast sphere has nonzero scattering')
   call require_close(absorption, 0.0_real64, 2.0e-13_real64, &
                      'zero-contrast sphere has nonzero absorption')

   sphere_index = cmplx(1.5_real64, 0.05_real64, kind=real64)
   call optically_active_mie_coefficients(size_parameter, sphere_index, order, 1.0e-13_real64, &
                                          extinction, scattering, absorption, ri_medium=medium_index)
   call require(ieee_is_finite(extinction) .and. ieee_is_finite(scattering) .and. &
                ieee_is_finite(absorption), 'absorbing-sphere efficiencies are not finite')
   call require(absorption > 0.0_real64, 'passive absorbing sphere has nonpositive absorption')
   call require_close(extinction, scattering + absorption, 2.0e-13_real64, &
                      'absorbing Mie sphere does not conserve energy')

contains

   subroutine efficiencies_from_coefficients(mie_coefficients, qext, qsca)
      complex(real64), intent(in) :: mie_coefficients(:, :, :)
      real(real64), intent(out) :: qext, qsca
      integer :: degree

      qext = 0.0_real64
      qsca = 0.0_real64
      do degree = 1, size(mie_coefficients, 3)
         qext = qext - real(2 * degree + 1, real64) * &
                real(mie_coefficients(1, 1, degree) + mie_coefficients(2, 2, degree), real64)
         qsca = qsca + real(2 * degree + 1, real64) * &
                sum(abs(mie_coefficients(:, :, degree))**2)
      end do
      qext = 2.0_real64 * qext / size_parameter**2
      qsca = 2.0_real64 * qsca / size_parameter**2
   end subroutine efficiencies_from_coefficients

   subroutine require_close(actual, expected, tolerance, message)
      real(real64), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message

      if (.not. ieee_is_finite(actual) .or. abs(actual - expected) > tolerance) then
         write (*, '(a,2(a,es24.16))') trim(message), ': actual=', actual, ', expected=', expected
         error stop 1
      end if
   end subroutine require_close

   subroutine require(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) error stop message
   end subroutine require
end program mie_test
