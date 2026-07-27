program concurrent_kernels_test
   use angular_functions, only: cartesian_vectors_to_spherical
   use iso_fortran_env, only: real64
   use random_configuration_geometry, only: circumscribing_sphere
   use scattering_amplitudes, only: amplitude_to_scattering_matrix
   use wave_functions, only: invert_two_by_two_matrix, reverse_azimuthal_modes

   implicit none(type, external)

   integer, parameter :: order = 2
   integer, parameter :: coefficient_count = 2 * order * (order + 2)
   real(real64), parameter :: tolerance = 1.0e-12_real64
   complex(real64) :: coefficients(coefficient_count), restored(coefficient_count), transformed(coefficient_count)
   complex(real64) :: inverse(2, 2), matrix(2, 2), identity(2, 2)
   complex(real64) :: amplitudes(4)
   real(real64) :: expected_scattering(4, 4), positions(3, 3), radii(3), spherical(3, 3), scattering(4, 4)
   integer :: coefficient

   positions = reshape([3.0_real64, 4.0_real64, 0.0_real64, &
                        0.0_real64, 0.0_real64, -2.0_real64, &
                        0.0_real64, 0.0_real64, 0.0_real64], shape(positions))
   call cartesian_vectors_to_spherical(size(positions, dim=2), positions, spherical)
   if (maxval(abs(spherical(:, 1) - [0.0_real64, atan2(4.0_real64, 3.0_real64), 5.0_real64])) > tolerance) &
      error stop 'Cartesian-to-spherical concurrent loop failed'
   if (maxval(abs(spherical(:, 2) - [-1.0_real64, 0.0_real64, 2.0_real64])) > tolerance) &
      error stop 'Cartesian-to-spherical axial case failed'
   if (maxval(abs(spherical(:, 3) - [1.0_real64, 0.0_real64, 0.0_real64])) > tolerance) &
      error stop 'Cartesian-to-spherical origin case failed'

   positions = reshape([3.0_real64, 4.0_real64, 0.0_real64, &
                        0.0_real64, 0.0_real64, 2.0_real64, &
                        -1.0_real64, -2.0_real64, -2.0_real64], shape(positions))
   radii = [1.0_real64, 0.5_real64, 2.0_real64]
   call circumscribing_sphere(size(radii), radii, positions, spherical(1, 1))
   if (abs(spherical(1, 1) - 6.0_real64) > tolerance) error stop 'Concurrent maximum reduction failed'

   matrix = reshape([cmplx(2.0_real64, 0.0_real64, kind=real64), &
                     cmplx(1.0_real64, -1.0_real64, kind=real64), &
                     cmplx(0.5_real64, 0.25_real64, kind=real64), &
                     cmplx(3.0_real64, 0.0_real64, kind=real64)], shape(matrix))
   call invert_two_by_two_matrix(matrix, inverse)
   identity = cmplx(0.0_real64, 0.0_real64, kind=real64)
   identity(1, 1) = cmplx(1.0_real64, 0.0_real64, kind=real64)
   identity(2, 2) = cmplx(1.0_real64, 0.0_real64, kind=real64)
   if (maxval(abs(matmul(matrix, inverse) - identity)) > tolerance) error stop 'Concurrent matrix inverse failed'

   do coefficient = 1, coefficient_count
      coefficients(coefficient) = cmplx(0.25_real64 * coefficient, -0.125_real64 * coefficient, kind=real64)
   end do
   call reverse_azimuthal_modes(order, coefficients, transformed)
   call reverse_azimuthal_modes(order, transformed, restored)
   if (maxval(abs(restored - coefficients)) > tolerance) error stop 'Concurrent degree transformation failed'

   amplitudes = cmplx(0.0_real64, 0.0_real64, kind=real64)
   amplitudes(1) = cmplx(1.0_real64, 0.0_real64, kind=real64)
   call amplitude_to_scattering_matrix(amplitudes, scattering)
   expected_scattering = 0.0_real64
   expected_scattering(1, 1) = 1.0_real64
   expected_scattering(1, 2) = -1.0_real64
   expected_scattering(2, 1) = -1.0_real64
   expected_scattering(2, 2) = 1.0_real64
   if (maxval(abs(scattering - expected_scattering)) > tolerance) &
      error stop 'Concurrent scattering-matrix construction failed'
end program concurrent_kernels_test
