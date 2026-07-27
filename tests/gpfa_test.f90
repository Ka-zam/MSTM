program gpfa_test
   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: two_pi
   use gpfa_controller, only: cgpfa
   use gpfa_setup, only: setgpfa
   implicit none

   integer, parameter :: transform_size = 30
   integer, parameter :: trig_size = 2 * (2 + 3 + 5)
   integer :: input_index, output_index
   ! The inherited GPFA kernels use legacy trigonometric recurrences whose
   ! observed double-precision error is approximately 1.7e-7 for this case.
   real(real64), parameter :: tolerance = 5.d-7
   real(real64) :: angle, max_error
   real(real64) :: real_data(transform_size), imag_data(transform_size)
   real(real64) :: real_input(transform_size), imag_input(transform_size)
   real(real64) :: trigs(trig_size)
   complex(real64) :: expected, actual

   do input_index = 1, transform_size
      real_input(input_index) = sin(0.2d0 * input_index) + 0.03d0 * input_index
      imag_input(input_index) = cos(0.3d0 * input_index) - 0.02d0 * input_index
   end do
   real_data = real_input
   imag_data = imag_input

   call setgpfa(trigs, transform_size)
   call cgpfa(real_data, imag_data, trigs, 1, transform_size, 1)

   max_error = 0.d0
   do output_index = 1, transform_size
      expected = (0.d0, 0.d0)
      do input_index = 1, transform_size
         angle = two_pi * dble((output_index - 1) * (input_index - 1)) / dble(transform_size)
         expected = expected + cmplx(real_input(input_index), imag_input(input_index), kind=real64) &
                    * cmplx(cos(angle), sin(angle), kind=real64)
      end do
      actual = cmplx(real_data(output_index), imag_data(output_index), kind=real64)
      max_error = max(max_error, abs(actual - expected))
   end do

   if (max_error > tolerance) then
      write (*, '(a,es12.4)') 'GPFA maximum error: ', max_error
      error stop 'GPFA result differs from direct DFT'
   end if
end program gpfa_test
