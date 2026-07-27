module interaction_operators
   use, intrinsic :: iso_fortran_env, only: real64
   use fft_translation, only: fft_plan, fft_translation_plan_t
   use translation_expansions, only: external_to_external_expansion
   implicit none(type, external)
   private

   type, abstract, public :: interaction_operator_t
   contains
      procedure(apply_interaction), deferred, public :: apply
   end type interaction_operator_t

   type, extends(interaction_operator_t), public :: pairwise_interaction_operator_t
   contains
      procedure, public :: apply => apply_pairwise_interaction
   end type pairwise_interaction_operator_t

   type, extends(interaction_operator_t), public :: fft_interaction_operator_t
      private
      type(fft_translation_plan_t), pointer :: plan => null()
   contains
      procedure, public :: apply => apply_fft_interaction
      procedure, public :: attach => attach_fft_plan
   end type fft_interaction_operator_t

   public :: create_interaction_operator

   abstract interface
      subroutine apply_interaction(self, equation_count, right_hand_side_count, input_coefficients, &
                                   output_coefficients, store_matrix, initial_run, right_hand_side_list, &
                                   communicator, conjugate_transpose)
         import interaction_operator_t, real64
         class(interaction_operator_t), intent(inout) :: self
         integer, intent(in) :: equation_count, right_hand_side_count, communicator
         logical, intent(in) :: store_matrix, initial_run
         logical, intent(in) :: right_hand_side_list(right_hand_side_count)
         logical, intent(in) :: conjugate_transpose(right_hand_side_count)
         complex(real64), intent(in) :: input_coefficients(equation_count, right_hand_side_count)
         complex(real64), intent(out) :: output_coefficients(equation_count, right_hand_side_count)
      end subroutine apply_interaction
   end interface

contains

   subroutine create_interaction_operator(use_fft, interaction_operator)
      logical, intent(in) :: use_fft
      class(interaction_operator_t), allocatable, intent(out) :: interaction_operator

      if (use_fft) then
         allocate (fft_interaction_operator_t :: interaction_operator)
         select type (interaction_operator)
         type is (fft_interaction_operator_t)
            call interaction_operator%attach(fft_plan)
         end select
      else
         allocate (pairwise_interaction_operator_t :: interaction_operator)
      end if
   end subroutine create_interaction_operator

   subroutine attach_fft_plan(self, plan)
      class(fft_interaction_operator_t), intent(inout) :: self
      type(fft_translation_plan_t), target, intent(inout) :: plan

      self%plan => plan
   end subroutine attach_fft_plan

   subroutine apply_pairwise_interaction(self, equation_count, right_hand_side_count, input_coefficients, &
                                         output_coefficients, store_matrix, initial_run, right_hand_side_list, &
                                         communicator, conjugate_transpose)
      class(pairwise_interaction_operator_t), intent(inout) :: self
      integer, intent(in) :: equation_count, right_hand_side_count, communicator
      logical, intent(in) :: store_matrix, initial_run
      logical, intent(in) :: right_hand_side_list(right_hand_side_count)
      logical, intent(in) :: conjugate_transpose(right_hand_side_count)
      complex(real64), intent(in) :: input_coefficients(equation_count, right_hand_side_count)
      complex(real64), intent(out) :: output_coefficients(equation_count, right_hand_side_count)

      call external_to_external_expansion(equation_count, right_hand_side_count, input_coefficients, &
                                          output_coefficients, store_matrix_option=store_matrix, &
                                          initial_run=initial_run, rhs_list=right_hand_side_list, &
                                          mpi_comm=communicator, con_tran=conjugate_transpose)
   end subroutine apply_pairwise_interaction

   subroutine apply_fft_interaction(self, equation_count, right_hand_side_count, input_coefficients, &
                                    output_coefficients, store_matrix, initial_run, right_hand_side_list, &
                                    communicator, conjugate_transpose)
      class(fft_interaction_operator_t), intent(inout) :: self
      integer, intent(in) :: equation_count, right_hand_side_count, communicator
      logical, intent(in) :: store_matrix, initial_run
      logical, intent(in) :: right_hand_side_list(right_hand_side_count)
      logical, intent(in) :: conjugate_transpose(right_hand_side_count)
      complex(real64), intent(in) :: input_coefficients(equation_count, right_hand_side_count)
      complex(real64), intent(out) :: output_coefficients(equation_count, right_hand_side_count)

      if (.not. associated(self%plan)) error stop 'FFT interaction operator has no translation plan'
      call self%plan%apply(equation_count, right_hand_side_count, input_coefficients, output_coefficients, &
                           store_matrix_option=store_matrix, initial_run=initial_run, &
                           rhs_list=right_hand_side_list, mpi_comm=communicator, &
                           con_tran=conjugate_transpose)
   end subroutine apply_fft_interaction
end module interaction_operators
