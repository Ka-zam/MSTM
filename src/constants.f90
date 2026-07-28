module constants
   use iso_fortran_env, only: real64
   implicit none
   private

   real(real64), parameter, public :: pi = acos(-1.0_real64)
   real(real64), parameter, public :: quarter_pi = pi / 4.0_real64
   real(real64), parameter, public :: three_quarters_pi = 3.0_real64 * quarter_pi
   real(real64), parameter, public :: two_pi = 2.0_real64 * pi
   real(real64), parameter, public :: four_pi = 4.0_real64 * pi
   real(real64), parameter, public :: four_pi_over_three = four_pi / 3.0_real64
   real(real64), parameter, public :: degrees_to_radians = pi / 180.0_real64
   real(real64), parameter, public :: sqrt_two_pi = sqrt(two_pi)
   complex(real64), parameter, public :: imaginary_unit = (0.0_real64, 1.0_real64)
end module constants
