module excitation
   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: imaginary_unit, pi, two_pi
   use quadrature, only: gauss_legendre_rule
   use wave_functions, only: vector_spherical_wave_functions
   implicit none(type, external)
   private

   character(len=*), parameter, public :: plane_wave_excitation = 'plane_wave'
   character(len=*), parameter, public :: electric_dipole_excitation = 'electric_dipole'
   character(len=*), parameter, public :: magnetic_current_excitation = 'magnetic_current_segments'

   type, public :: magnetic_current_segment_t
      real(real64) :: start_point(3) = 0.0_real64
      real(real64) :: end_point(3) = 0.0_real64
      complex(real64) :: amplitude = (1.0_real64, 0.0_real64)
   end type magnetic_current_segment_t

   public :: electric_dipole_field, electric_dipole_mode_coefficients, electric_dipole_source_power, &
             magnetic_current_field, magnetic_current_mode_coefficients, magnetic_current_source_power

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

   pure subroutine magnetic_current_mode_coefficients(moment, coefficients)
      complex(real64), intent(in) :: moment(3)
      complex(real64), intent(out) :: coefficients(0:2, 1, 2)

      ! A differential magnetic current excites the magnetic/TE n=1 mode.
      ! It uses the same Cartesian normalization as the electric/TM mode.
      coefficients = (0.0_real64, 0.0_real64)
      coefficients(0, 1, 2) = sqrt(3.0_real64) * moment(3)
      coefficients(1, 1, 2) = sqrt(1.5_real64) * (-moment(1) + imaginary_unit * moment(2))
      coefficients(2, 1, 2) = sqrt(1.5_real64) * (moment(1) + imaginary_unit * moment(2))
   end subroutine magnetic_current_mode_coefficients

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

   subroutine magnetic_current_field(position, refractive_index, segments, quadrature_order, &
                                     electric_field, magnetic_field, singular)
      real(real64), intent(in) :: position(3)
      complex(real64), intent(in) :: refractive_index
      type(magnetic_current_segment_t), intent(in) :: segments(:)
      integer, intent(in) :: quadrature_order
      complex(real64), intent(out) :: electric_field(3), magnetic_field(3)
      logical, intent(out), optional :: singular
      integer :: segment_index, node_index
      real(real64) :: delta(3), source_position(3), parameter_nodes(quadrature_order), &
                      parameter_weights(quadrature_order)
      complex(real64) :: differential_moment(3), electric_contribution(3), magnetic_contribution(3)
      logical :: at_source

      electric_field = (0.0_real64, 0.0_real64)
      magnetic_field = (0.0_real64, 0.0_real64)
      at_source = .false.
      call gauss_legendre_rule(0.0_real64, 1.0_real64, parameter_nodes, parameter_weights, quadrature_order)

      do segment_index = 1, size(segments)
         delta = segments(segment_index)%end_point - segments(segment_index)%start_point
         if (point_segment_distance(position, segments(segment_index)%start_point, &
                                    segments(segment_index)%end_point) <= 1.0e-12_real64) then
            at_source = .true.
            cycle
         end if
         do node_index = 1, quadrature_order
            source_position = segments(segment_index)%start_point + parameter_nodes(node_index) * delta
            differential_moment = segments(segment_index)%amplitude * delta * parameter_weights(node_index)
            call magnetic_dipole_field(position - source_position, refractive_index, differential_moment, &
                                       electric_contribution, magnetic_contribution)
            electric_field = electric_field + electric_contribution
            magnetic_field = magnetic_field + magnetic_contribution
         end do
      end do

      if (at_source) then
         electric_field = (0.0_real64, 0.0_real64)
         magnetic_field = (0.0_real64, 0.0_real64)
      end if
      if (present(singular)) singular = at_source
   end subroutine magnetic_current_field

   real(real64) function magnetic_current_source_power(segments, refractive_index, angular_order)
      type(magnetic_current_segment_t), intent(in) :: segments(:)
      real(real64), intent(in) :: refractive_index
      integer, intent(in), optional :: angular_order
      integer :: azimuth_index, azimuth_order, polar_index, polar_order, segment_index
      real(real64) :: azimuth, delta(3), diameter, direction(3), integration_sum, midpoint(3), &
                      mu, reference_position(3), transverse_power
      real(real64), allocatable :: mu_nodes(:), mu_weights(:)
      complex(real64) :: far_field_vector(3), phase

      if (size(segments) == 0) then
         magnetic_current_source_power = 0.0_real64
         return
      end if
      if (present(angular_order)) then
         polar_order = angular_order
      else
         diameter = source_diameter(segments)
         polar_order = max(48, 16 + ceiling(4.0_real64 * refractive_index * diameter))
      end if
      azimuth_order = 2 * polar_order
      allocate (mu_nodes(polar_order), mu_weights(polar_order))
      call gauss_legendre_rule(-1.0_real64, 1.0_real64, mu_nodes, mu_weights, polar_order)
      reference_position = 0.5_real64 * (segments(1)%start_point + segments(1)%end_point)
      integration_sum = 0.0_real64

      do polar_index = 1, polar_order
         mu = mu_nodes(polar_index)
         do azimuth_index = 1, azimuth_order
            azimuth = two_pi * real(azimuth_index - 1, real64) / real(azimuth_order, real64)
            direction = [sqrt(max(0.0_real64, 1.0_real64 - mu * mu)) * cos(azimuth), &
                         sqrt(max(0.0_real64, 1.0_real64 - mu * mu)) * sin(azimuth), mu]
            far_field_vector = (0.0_real64, 0.0_real64)
            do segment_index = 1, size(segments)
               delta = segments(segment_index)%end_point - segments(segment_index)%start_point
               midpoint = 0.5_real64 * (segments(segment_index)%start_point + &
                                        segments(segment_index)%end_point)
               phase = exp(-imaginary_unit * refractive_index * dot_product(direction, midpoint - reference_position))
               far_field_vector = far_field_vector + segments(segment_index)%amplitude * delta * &
                                  normalized_sinc(0.5_real64 * refractive_index * dot_product(direction, delta)) * phase
            end do
            transverse_power = sum(abs(far_field_vector)**2) - &
                               abs(dot_product(direction, far_field_vector))**2
            integration_sum = integration_sum + mu_weights(polar_index) * transverse_power
         end do
      end do

      ! Divide the full angular radiation integral by 2*pi to match the
      ! modal-power convention used by electric_dipole_source_power().
      magnetic_current_source_power = 9.0_real64 * integration_sum / &
                                      (4.0_real64 * real(azimuth_order, real64))
   end function magnetic_current_source_power

   subroutine magnetic_dipole_field(displacement, refractive_index, moment, electric_field, magnetic_field)
      real(real64), intent(in) :: displacement(3)
      complex(real64), intent(in) :: refractive_index, moment(3)
      complex(real64), intent(out) :: electric_field(3), magnetic_field(3)
      complex(real64) :: coefficients(0:2, 1, 2), wave_functions(3, 6)

      call magnetic_current_mode_coefficients(moment, coefficients)
      call vector_spherical_wave_functions(displacement, [refractive_index, refractive_index], 1, 3, &
                                           wave_functions, index_model=2, lr_to_mode=.true.)
      electric_field = matmul(wave_functions(:, 4:6), coefficients(0:2, 1, 2))
      magnetic_field = matmul(wave_functions(:, 1:3), coefficients(0:2, 1, 2)) &
                       * refractive_index / imaginary_unit
   end subroutine magnetic_dipole_field

   pure real(real64) function normalized_sinc(argument)
      real(real64), intent(in) :: argument

      if (abs(argument) <= 1.0e-8_real64) then
         normalized_sinc = 1.0_real64 - argument**2 / 6.0_real64
      else
         normalized_sinc = sin(argument) / argument
      end if
   end function normalized_sinc

   pure real(real64) function point_segment_distance(point, start_point, end_point)
      real(real64), intent(in) :: point(3), start_point(3), end_point(3)
      real(real64) :: delta(3), projection

      delta = end_point - start_point
      projection = dot_product(point - start_point, delta) / max(tiny(1.0_real64), sum(delta**2))
      projection = max(0.0_real64, min(1.0_real64, projection))
      point_segment_distance = sqrt(sum((point - start_point - projection * delta)**2))
   end function point_segment_distance

   pure real(real64) function source_diameter(segments)
      type(magnetic_current_segment_t), intent(in) :: segments(:)
      integer :: first_segment, second_segment

      source_diameter = 0.0_real64
      do first_segment = 1, size(segments)
         do second_segment = 1, size(segments)
            source_diameter = max(source_diameter, &
                                  sqrt(sum((segments(first_segment)%start_point - &
                                            segments(second_segment)%start_point)**2)), &
                                  sqrt(sum((segments(first_segment)%start_point - &
                                            segments(second_segment)%end_point)**2)), &
                                  sqrt(sum((segments(first_segment)%end_point - &
                                            segments(second_segment)%start_point)**2)), &
                                  sqrt(sum((segments(first_segment)%end_point - &
                                            segments(second_segment)%end_point)**2)))
         end do
      end do
   end function source_diameter

end module excitation
