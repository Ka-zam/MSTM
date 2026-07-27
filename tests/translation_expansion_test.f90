program translation_expansion_test
   use iso_fortran_env, only: real64
   use mpidefs, only: mstm_mpi
   use numerical_tables, only: initialize_numerical_tables
   use translation, only: external_to_internal_expansion, nested_sphere_geometry_view

   implicit none(type, external)

   integer, parameter :: test_order = 2
   integer, parameter :: block_size = 2 * test_order * (test_order + 2)
   integer, parameter :: equation_count = 3 * block_size
   integer, parameter :: rhs_count = 2
   integer :: coefficient
   real(real64) :: tolerance
   complex(real64) :: input(equation_count, rhs_count)
   complex(real64) :: all_output(equation_count, rhs_count)
   complex(real64) :: first_output(equation_count, rhs_count)
   complex(real64) :: second_output(equation_count, rhs_count)
   type(nested_sphere_geometry_view) :: geometry

   call mstm_mpi(mpi_command='init')
   call configure_nested_spheres()

   do coefficient = 1, equation_count
      input(coefficient, 1) = cmplx(0.02_real64 * coefficient, -0.01_real64 * coefficient, kind=real64)
      input(coefficient, 2) = cmplx(-0.03_real64 * coefficient, 0.015_real64 * coefficient, kind=real64)
   end do

   call external_to_internal_expansion(equation_count, rhs_count, input, all_output, &
                                       rhs_list=[.true., .true.], geometry=geometry)
   call external_to_internal_expansion(equation_count, rhs_count, input, first_output, &
                                       rhs_list=[.true., .false.], geometry=geometry)
   call external_to_internal_expansion(equation_count, rhs_count, input, second_output, &
                                       rhs_list=[.false., .true.], geometry=geometry)

   tolerance = 100.0_real64 * epsilon(1.0_real64) * max(1.0_real64, maxval(abs(all_output)))
   if (maxval(abs(all_output(:, 1))) <= tolerance .or. maxval(abs(all_output(:, 2))) <= tolerance) then
      error stop 'Nested-sphere translation produced no output'
   end if
   if (maxval(abs(first_output(:, 1) - all_output(:, 1))) > tolerance &
       .or. maxval(abs(first_output(:, 2))) > tolerance) then
      error stop 'First selective RHS result is incorrect'
   end if
   if (maxval(abs(second_output(:, 1))) > tolerance &
       .or. maxval(abs(second_output(:, 2) - all_output(:, 2))) > tolerance) then
      error stop 'Second selective RHS result is incorrect'
   end if

   call mstm_mpi(mpi_command='finalize')

contains

   subroutine configure_nested_spheres()
      integer, target, save :: hosts(2), blocks(2), offsets(2), orders(2)
      real(real64), target, save :: positions(3, 2)
      complex(real64), target, save :: refractive_indices(2, 2)

      hosts = [2, 0]
      orders = test_order
      blocks = block_size
      offsets = [0, block_size]
      positions(:, 1) = [0.25_real64, -0.15_real64, 0.10_real64]
      positions(:, 2) = 0.0_real64
      refractive_indices = cmplx(1.0_real64, 0.0_real64, kind=real64)
      refractive_indices(:, 2) = cmplx(1.4_real64, 0.0_real64, kind=real64)

      call initialize_numerical_tables(2 * test_order)
      call geometry%configure(hosts, blocks, offsets, orders, positions, refractive_indices)
   end subroutine configure_nested_spheres

end program translation_expansion_test
