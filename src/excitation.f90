module excitation
   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: imaginary_unit
   use wave_functions, only: vector_spherical_wave_functions
   implicit none(type, external)
   private

   character(len=*), parameter, public :: plane_wave_excitation = 'plane_wave'
   character(len=*), parameter, public :: electric_dipole_excitation = 'electric_dipole'

   public :: electric_dipole_field, electric_dipole_mode_coefficients, electric_dipole_source_power

contains

   pure subroutine electric_dipole_mode_coefficients(moment, coefficients)
      complex(real64), intent(in) :: moment(3)
      complex(real64), intent(out) :: coefficients(0:2, 1, 2)

      ! The first mode is electric/TM.  The normalization is chosen so that
      ! its regular n=1 field at the origin equals the Cartesian moment.
      coefficients = (0.0_real64, 0.0_real64)
      coefficients(0, 1, 1) = sqrt(3.0_real64) * moment(3)
      coefficients(1, 1, 1) = sqrt(1.5_real64) * (-moment(1) + imaginary_unit * moment(2))
      coefficients(2, 1, 1) = sqrt(1.5_real64) * (moment(1) + imaginary_unit * moment(2))
   end subroutine electric_dipole_mode_coefficients

   pure real(real64) function electric_dipole_source_power(moment)
      complex(real64), intent(in) :: moment(3)

      electric_dipole_source_power = 3.0_real64 * sum(abs(moment)**2)
   end function electric_dipole_source_power

   subroutine electric_dipole_field(displacement, refractive_index, moment, electric_field, magnetic_field, singular)
      real(real64), intent(in) :: displacement(3)
      complex(real64), intent(in) :: refractive_index, moment(3)
      complex(real64), intent(out) :: electric_field(3), magnetic_field(3)
      logical, intent(out), optional :: singular
      complex(real64) :: coefficients(0:2, 1, 2), wave_functions(3, 6)
      logical :: at_source

      at_source = sum(displacement**2) <= 1.0e-24_real64
      if (present(singular)) singular = at_source
      if (at_source) then
         electric_field = (0.0_real64, 0.0_real64)
         magnetic_field = (0.0_real64, 0.0_real64)
         return
      end if

      call electric_dipole_mode_coefficients(moment, coefficients)
      call vector_spherical_wave_functions(displacement, [refractive_index, refractive_index], 1, 3, &
                                           wave_functions, index_model=2, lr_to_mode=.true.)
      electric_field = matmul(wave_functions(:, 1:3), coefficients(0:2, 1, 1))
      magnetic_field = matmul(wave_functions(:, 4:6), coefficients(0:2, 1, 1)) &
                       * refractive_index / imaginary_unit
   end subroutine electric_dipole_field

end module excitation
