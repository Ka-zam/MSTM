program bessel_test
   use bessel_functions, only: bessel_integer_complex
   implicit none

   integer, parameter :: maximum_order = 8
   integer :: highest_order, order
   real(8), parameter :: tolerance = 1.d-12
   real(8) :: max_error, real_argument
   complex(8) :: actual(0:maximum_order), complex_argument, expected

   real_argument = 3.25d0
   call bessel_integer_complex(maximum_order, cmplx(real_argument, 0.d0, kind=8), highest_order, actual)
   max_error = maxval(abs(actual - cmplx(bessel_jn(0, maximum_order, real_argument), 0.d0, kind=8)))
   if (highest_order /= maximum_order .or. max_error > tolerance) then
      write (*, '(a,es12.4)') 'Real Bessel maximum error: ', max_error
      error stop 'Real Bessel branch differs from the Fortran intrinsic'
   end if

   complex_argument = cmplx(0.8d0, 0.3d0, kind=8)
   call bessel_integer_complex(maximum_order, complex_argument, highest_order, actual)
   max_error = 0.d0
   do order = 0, maximum_order
      expected = bessel_series(order, complex_argument)
      max_error = max(max_error, abs(actual(order) - expected))
   end do
   if (highest_order /= maximum_order .or. max_error > tolerance) then
      write (*, '(a,es12.4)') 'Complex Bessel maximum error: ', max_error
      error stop 'Complex Bessel branch differs from its converged series'
   end if

contains

   complex(8) function bessel_series(order, argument) result(value)
      integer, intent(in) :: order
      complex(8), intent(in) :: argument
      integer :: index
      complex(8) :: term

      term = (argument / 2.d0)**order
      do index = 2, order
         term = term / dble(index)
      end do
      value = term
      do index = 1, 80
         term = -term * argument * argument / (4.d0 * dble(index) * dble(order + index))
         value = value + term
         if (abs(term) <= epsilon(1.d0) * max(1.d0, abs(value))) exit
      end do
   end function bessel_series
end program bessel_test
