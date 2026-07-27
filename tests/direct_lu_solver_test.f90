program direct_lu_solver_test
   use, intrinsic :: iso_fortran_env, only: real64
   use direct_lu_solver, only: direct_lu_solver_t, solver_converged, solver_singular
   implicit none(type, external)

   type(direct_lu_solver_t) :: solver
   complex(real64) :: matrix(3, 3), right_hand_side(3), solution(3), expected(3)
   integer :: refinements, status
   real(real64) :: residual

   matrix = reshape([cmplx(4.0_real64, 1.0_real64, real64), cmplx(1.0_real64, -1.0_real64, real64), &
                     cmplx(0.5_real64, 0.0_real64, real64), cmplx(1.0_real64, 2.0_real64, real64), &
                     cmplx(5.0_real64, 0.0_real64, real64), cmplx(-1.0_real64, 0.5_real64, real64), &
                     cmplx(0.0_real64, -0.5_real64, real64), cmplx(2.0_real64, 0.0_real64, real64), &
                     cmplx(6.0_real64, -1.0_real64, real64)], shape(matrix))
   expected = [cmplx(1.0_real64, -0.5_real64, real64), cmplx(-2.0_real64, 1.0_real64, real64), &
               cmplx(0.25_real64, 0.75_real64, real64)]
   right_hand_side = matmul(matrix, expected)

   call solver%factor(matrix, status)
   call require(status == solver_converged, 'factorization failed')
   call require(solver%is_factorized(), 'factorized state was not retained')
   call require(solver%condition_estimate() > 0.0_real64, 'condition estimate was not calculated')
   call solver%solve(right_hand_side, solution, 1.0e-13_real64, 3, refinements, residual, status)
   call require(status == solver_converged, 'first solve failed')
   call require(maxval(abs(solution - expected)) < 1.0e-12_real64, 'first solution is inaccurate')
   call require(residual < 1.0e-13_real64, 'first residual is too large')

   expected = 2.0_real64 * expected
   right_hand_side = matmul(matrix, expected)
   call solver%solve(right_hand_side, solution, 1.0e-13_real64, 3, refinements, residual, status)
   call require(status == solver_converged, 'reused solve failed')
   call require(maxval(abs(solution - expected)) < 1.0e-12_real64, 'reused solution is inaccurate')

   call solver%clear()
   call require(.not. solver%is_factorized(), 'clear did not release factorization state')
   matrix = 0.0_real64
   call solver%factor(matrix, status)
   call require(status == solver_singular, 'singular matrix was not rejected')

contains

   subroutine require(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) error stop message
   end subroutine require
end program direct_lu_solver_test
