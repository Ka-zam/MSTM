module direct_lu_solver
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: real64
   use parallel_runtime, only: parallel_wall_time
   implicit none(type, external)
   private

   integer, parameter, public :: solver_converged = 0
   integer, parameter, public :: solver_iteration_limit = 1
   integer, parameter, public :: solver_breakdown = 2
   integer, parameter, public :: solver_singular = 3
   integer, parameter, public :: solver_non_finite = 4

   type, public :: direct_lu_solver_t
      private
      complex(real64), allocatable :: original_matrix(:, :)
      complex(real64), allocatable :: factored_matrix(:, :)
      integer, allocatable :: pivots(:)
      real(real64) :: reciprocal_condition = 0.0_real64
      real(real64) :: factorization_seconds = 0.0_real64
      real(real64) :: condition_estimation_seconds = 0.0_real64
      real(real64) :: backsolve_seconds = 0.0_real64
      logical :: factorized = .false.
   contains
      procedure, public :: factor => factor_direct_matrix
      procedure, public :: solve => solve_direct_rhs
      procedure, public :: clear => clear_direct_solver
      procedure, public :: is_factorized => direct_solver_is_factorized
      procedure, public :: condition_estimate => direct_solver_condition_estimate
      procedure, public :: factorization_time => direct_solver_factorization_time
      procedure, public :: condition_estimation_time => direct_solver_condition_estimation_time
      procedure, public :: backsolve_time => direct_solver_backsolve_time
      final :: finalize_direct_solver
   end type direct_lu_solver_t

   public :: solver_status_message

   interface
      subroutine zgetrf(m, n, a, lda, ipiv, info)
         import real64
         integer, intent(in) :: m, n, lda
         complex(real64), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*), info
      end subroutine zgetrf

      subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import real64
         character(len=1), intent(in) :: trans
         integer, intent(in) :: n, nrhs, lda, ldb, ipiv(*)
         complex(real64), intent(in) :: a(lda, *)
         complex(real64), intent(inout) :: b(ldb, *)
         integer, intent(out) :: info
      end subroutine zgetrs

      function zlange(norm, m, n, a, lda, work) result(matrix_norm)
         import real64
         character(len=1), intent(in) :: norm
         integer, intent(in) :: m, n, lda
         complex(real64), intent(in) :: a(lda, *)
         real(real64), intent(out) :: work(*)
         real(real64) :: matrix_norm
      end function zlange

      subroutine zgecon(norm, n, a, lda, matrix_norm, reciprocal_condition, work, real_work, info)
         import real64
         character(len=1), intent(in) :: norm
         integer, intent(in) :: n, lda
         complex(real64), intent(in) :: a(lda, *)
         real(real64), intent(in) :: matrix_norm
         real(real64), intent(out) :: reciprocal_condition
         complex(real64), intent(out) :: work(*)
         real(real64), intent(out) :: real_work(*)
         integer, intent(out) :: info
      end subroutine zgecon
   end interface

contains

   pure function solver_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=48) :: message

      select case (status)
      case (solver_converged)
         message = 'converged'
      case (solver_iteration_limit)
         message = 'iteration limit reached'
      case (solver_breakdown)
         message = 'iterative solver breakdown'
      case (solver_singular)
         message = 'singular or numerically singular matrix'
      case (solver_non_finite)
         message = 'non-finite solver result'
      case default
         message = 'unknown solver status'
      end select
   end function solver_status_message

   subroutine factor_direct_matrix(self, matrix, status)
      class(direct_lu_solver_t), intent(inout) :: self
      complex(real64), intent(in) :: matrix(:, :)
      integer, intent(out) :: status
      integer :: matrix_size, lapack_status
      real(real64) :: matrix_norm, phase_start
      real(real64), allocatable :: norm_work(:), condition_real_work(:)
      complex(real64), allocatable :: condition_work(:)

      call self%clear()
      status = solver_converged
      if (size(matrix, 1) /= size(matrix, 2) .or. size(matrix, 1) == 0) then
         status = solver_breakdown
         return
      end if

      matrix_size = size(matrix, 1)
      allocate (self%original_matrix(matrix_size, matrix_size), &
                self%factored_matrix(matrix_size, matrix_size), self%pivots(matrix_size))
      allocate (norm_work(matrix_size), condition_work(2 * matrix_size), &
                condition_real_work(2 * matrix_size))
      self%original_matrix = matrix
      self%factored_matrix = matrix

      phase_start = parallel_wall_time()
      matrix_norm = zlange('1', matrix_size, matrix_size, self%original_matrix, matrix_size, norm_work)
      self%condition_estimation_seconds = parallel_wall_time() - phase_start

      phase_start = parallel_wall_time()
      call zgetrf(matrix_size, matrix_size, self%factored_matrix, matrix_size, self%pivots, lapack_status)
      self%factorization_seconds = parallel_wall_time() - phase_start
      if (lapack_status > 0) then
         status = solver_singular
      elseif (lapack_status < 0) then
         status = solver_breakdown
      else
         phase_start = parallel_wall_time()
         call zgecon('1', matrix_size, self%factored_matrix, matrix_size, matrix_norm, &
                     self%reciprocal_condition, condition_work, condition_real_work, lapack_status)
         self%condition_estimation_seconds = self%condition_estimation_seconds + &
                                             parallel_wall_time() - phase_start
         if (lapack_status /= 0) then
            status = solver_breakdown
         elseif (.not. ieee_is_finite(self%reciprocal_condition) .or. &
                 self%reciprocal_condition <= tiny(1.0_real64)) then
            status = solver_singular
         else
            self%factorized = .true.
         end if
      end if
      deallocate (norm_work, condition_work, condition_real_work)
   end subroutine factor_direct_matrix

   subroutine solve_direct_rhs(self, right_hand_side, solution, tolerance, maximum_refinements, &
                               refinements, residual, status)
      class(direct_lu_solver_t), intent(inout) :: self
      complex(real64), intent(in) :: right_hand_side(:)
      complex(real64), intent(out) :: solution(size(right_hand_side))
      real(real64), intent(in) :: tolerance
      integer, intent(in) :: maximum_refinements
      integer, intent(out) :: refinements, status
      real(real64), intent(out) :: residual
      integer :: lapack_status, matrix_size
      real(real64) :: phase_start
      complex(real64) :: correction(size(right_hand_side)), residual_vector(size(right_hand_side))

      matrix_size = size(right_hand_side)
      refinements = 0
      residual = huge(1.0_real64)
      solution = 0.0_real64
      if (.not. self%factorized .or. size(self%original_matrix, 1) /= matrix_size) then
         status = solver_breakdown
         return
      end if

      status = solver_converged
      solution = right_hand_side
      phase_start = parallel_wall_time()
      call zgetrs('N', matrix_size, 1, self%factored_matrix, matrix_size, self%pivots, &
                  solution, matrix_size, lapack_status)
      self%backsolve_seconds = self%backsolve_seconds + parallel_wall_time() - phase_start
      if (lapack_status /= 0) status = solver_breakdown

      if (status == solver_converged) then
         residual_vector = right_hand_side - matmul(self%original_matrix, solution)
         residual = relative_residual_norm(residual_vector, right_hand_side)
      end if
      do while (status == solver_converged .and. residual > tolerance .and. &
                refinements < maximum_refinements)
         correction = residual_vector
         phase_start = parallel_wall_time()
         call zgetrs('N', matrix_size, 1, self%factored_matrix, matrix_size, self%pivots, &
                     correction, matrix_size, lapack_status)
         self%backsolve_seconds = self%backsolve_seconds + parallel_wall_time() - phase_start
         if (lapack_status /= 0) then
            status = solver_breakdown
            exit
         end if
         solution = solution + correction
         refinements = refinements + 1
         residual_vector = right_hand_side - matmul(self%original_matrix, solution)
         residual = relative_residual_norm(residual_vector, right_hand_side)
      end do

      if (status == solver_converged) then
         if (.not. ieee_is_finite(residual) .or. .not. complex_vector_is_finite(solution)) then
            status = solver_non_finite
         elseif (residual > tolerance) then
            status = solver_iteration_limit
         end if
      end if
   end subroutine solve_direct_rhs

   subroutine clear_direct_solver(self)
      class(direct_lu_solver_t), intent(inout) :: self

      if (allocated(self%original_matrix)) deallocate (self%original_matrix)
      if (allocated(self%factored_matrix)) deallocate (self%factored_matrix)
      if (allocated(self%pivots)) deallocate (self%pivots)
      self%reciprocal_condition = 0.0_real64
      self%factorization_seconds = 0.0_real64
      self%condition_estimation_seconds = 0.0_real64
      self%backsolve_seconds = 0.0_real64
      self%factorized = .false.
   end subroutine clear_direct_solver

   subroutine finalize_direct_solver(self)
      type(direct_lu_solver_t), intent(inout) :: self

      call self%clear()
   end subroutine finalize_direct_solver

   pure logical function direct_solver_is_factorized(self)
      class(direct_lu_solver_t), intent(in) :: self

      direct_solver_is_factorized = self%factorized
   end function direct_solver_is_factorized

   pure real(real64) function direct_solver_condition_estimate(self)
      class(direct_lu_solver_t), intent(in) :: self

      direct_solver_condition_estimate = self%reciprocal_condition
   end function direct_solver_condition_estimate

   pure real(real64) function direct_solver_factorization_time(self)
      class(direct_lu_solver_t), intent(in) :: self

      direct_solver_factorization_time = self%factorization_seconds
   end function direct_solver_factorization_time

   pure real(real64) function direct_solver_condition_estimation_time(self)
      class(direct_lu_solver_t), intent(in) :: self

      direct_solver_condition_estimation_time = self%condition_estimation_seconds
   end function direct_solver_condition_estimation_time

   pure real(real64) function direct_solver_backsolve_time(self)
      class(direct_lu_solver_t), intent(in) :: self

      direct_solver_backsolve_time = self%backsolve_seconds
   end function direct_solver_backsolve_time

   pure real(real64) function relative_residual_norm(residual_vector, right_hand_side)
      complex(real64), intent(in) :: residual_vector(:), right_hand_side(:)
      real(real64) :: denominator

      denominator = sqrt(sum(abs(right_hand_side)**2))
      if (denominator <= tiny(1.0_real64)) then
         relative_residual_norm = sqrt(sum(abs(residual_vector)**2))
      else
         relative_residual_norm = sqrt(sum(abs(residual_vector)**2)) / denominator
      end if
   end function relative_residual_norm

   pure logical function complex_vector_is_finite(vector)
      complex(real64), intent(in) :: vector(:)

      complex_vector_is_finite = all(ieee_is_finite(real(vector, real64))) .and. &
         all(ieee_is_finite(aimag(vector)))
   end function complex_vector_is_finite
end module direct_lu_solver
