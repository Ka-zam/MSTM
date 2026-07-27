program mstm
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit, real64
   use input_execution, only: execute_simulation, run_configuration_average, run_incidence_average, &
                              run_random_orientation_configuration_average
   use input_parser, only: parse_input_data, process_input_variable
   use input_reporting, only: output_header
   use input_state, only: simulation_config, simulation_result
   use mstm_version_info, only: mstm_version
   use parallel_runtime, only: parallel_barrier, parallel_finalize, parallel_initialize, parallel_rank, parallel_size
   use runtime_support, only: clear_runtime_status, open_input_file, open_output_file, &
                              report_runtime_error, runtime_failed, set_runtime_error, &
                              synchronize_runtime_status
   implicit none
   logical :: input_exists
   integer :: input_unit, looplevel, output_unit_number, rank, numprocs, &
              readstat(1), i, istat, numberinputlines, n
   character(len=256) :: inputfile, inputline, oldoutputfile
   character(len=256), allocatable :: inputfiledata(:)
   data oldoutputfile/' '/

   call parallel_initialize()
   call clear_runtime_status()
   call parallel_rank(mpi_rank=rank)
   call parallel_size(mpi_size=numprocs)
   call parse_command_line(rank, inputfile)

   inquire (file=trim(inputfile), exist=input_exists)
   if (.not. input_exists) then
      if (rank .eq. 0) then
         write (error_unit, '(a)') "mstm: cannot open input file '"//trim(inputfile)//"'"
         write (error_unit, '(a)') "Try 'mstm --help' for usage information."
      end if
      call parallel_finalize()
      stop 2, quiet = .true.
   end if

   do i = 0, numprocs - 1
      if (rank .eq. i) then
         simulation_config%input_file = inputfile
         call open_input_file(trim(inputfile), input_unit)
         if (.not. runtime_failed()) then
            numberinputlines = 0
            istat = 0
            do while (istat .eq. 0)
               read (input_unit, '(a)', iostat=istat) inputline
               numberinputlines = numberinputlines + 1
               if (trim(inputline) .eq. 'end_of_options') exit
            end do
            if (istat .ne. 0) numberinputlines = numberinputlines + 1
            allocate (inputfiledata(numberinputlines))
            rewind (input_unit)
            istat = 0
            n = 0
            do while (istat .eq. 0)
               read (input_unit, '(a)', iostat=istat) inputline
               n = n + 1
               inputfiledata(n) = inputline
               if (trim(inputline) .eq. 'end_of_options') exit
            end do
            if (istat .ne. 0) inputfiledata(numberinputlines) = 'end_of_options'
            close (input_unit)
         end if
      end if
      call parallel_barrier()
   end do
   call synchronize_runtime_status()
   if (runtime_failed()) then
      if (rank .eq. 0) call report_runtime_error(error_unit)
      call parallel_finalize()
      stop 2, quiet = .true.
   end if

   simulation_config%repeat_run = .true.
   simulation_config%first_run = .true.

   simulation_loop: do while (simulation_config%repeat_run)
      readstat = 0
      do i = 0, numprocs - 1
         if (i .eq. rank) then
            call parse_input_data(inputfiledata, read_status=readstat(1))
            if (readstat(1) /= 0 .and. .not. runtime_failed()) &
               call set_runtime_error('Invalid simulation input')
         end if
         call parallel_barrier()
      end do
      call synchronize_runtime_status()
      if (runtime_failed()) exit simulation_loop
      if (oldoutputfile .ne. simulation_config%output%output_file) then
         simulation_config%first_run = .true.
         simulation_result%run_number = 0
         oldoutputfile = simulation_config%output%output_file
      end if
      if (rank .eq. 0) then
         if (simulation_config%first_run) then
            call open_output_file(simulation_config%output%output_file, output_unit_number, append=simulation_config%output%append)
            if (.not. runtime_failed()) then
               call output_header(output_unit_number, inputfile)
               close (output_unit_number)
            end if
         end if
      end if
      call synchronize_runtime_status()
      if (runtime_failed()) exit simulation_loop

      if (simulation_config%number_nested_loops .eq. 0) then
         simulation_result%run_number = simulation_result%run_number + 1
         if (simulation_config%configuration_average) then
            if (simulation_config%random_orientation) then
               call run_random_orientation_configuration_average()
            else
               call run_configuration_average()
            end if
         elseif (simulation_config%incidence_average) then
            call run_incidence_average()
         else
            call execute_simulation()
         end if
      else
         looplevel = 1
         call execute_nested_loop(looplevel, rank)
      end if
      call synchronize_runtime_status()
      if (runtime_failed()) exit simulation_loop
      simulation_config%number_nested_loops = 0
   end do simulation_loop
!      if(simulation_config%temporary_position_file.and.rank.eq.0) then
!         open(20,file='temp_pos.dat')
!         close(20,status='delete')
!      endif
   call parallel_barrier()
   call parallel_finalize()
   if (runtime_failed()) then
      if (rank .eq. 0) call report_runtime_error(error_unit)
      stop 2, quiet = .true.
   end if

contains

   subroutine parse_command_line(rank, input_file_name)
      implicit none
      integer, intent(in) :: rank
      character(len=*), intent(out) :: input_file_name
      integer :: number_arguments
      character(len=256) :: argument, extra_argument

      number_arguments = command_argument_count()
      if (number_arguments .eq. 0) then
         if (rank .eq. 0) call print_help()
         call parallel_finalize()
         stop 0, quiet = .true.
      end if

      call get_command_argument(1, argument)
      select case (trim(argument))
      case ('help', '--help', '-h')
         if (rank .eq. 0) call print_help()
         call parallel_finalize()
         stop 0, quiet = .true.
      case ('version', '--version', '-V')
         if (rank .eq. 0) write (output_unit, '(a)') 'MSTM '//mstm_version
         call parallel_finalize()
         stop 0, quiet = .true.
      case default
         if (argument(1:1) .eq. '-') then
            if (rank .eq. 0) then
               write (error_unit, '(a)') "mstm: unknown option '"//trim(argument)//"'"
               write (error_unit, '(a)') "Try 'mstm --help' for usage information."
            end if
            call parallel_finalize()
            stop 2, quiet = .true.
         end if
      end select

      if (number_arguments .gt. 1) then
         call get_command_argument(2, extra_argument)
         if (rank .eq. 0) then
            write (error_unit, '(a)') "mstm: unexpected argument '"//trim(extra_argument)//"'"
            write (error_unit, '(a)') "Usage: mstm [INPUT_FILE]"
         end if
         call parallel_finalize()
         stop 2, quiet = .true.
      end if

      input_file_name = trim(argument)
   end subroutine parse_command_line

   subroutine print_help()
      implicit none

      write (output_unit, '(a)') 'MSTM '//mstm_version
      write (output_unit, '(a)') 'Multiple Sphere T-Matrix electromagnetic scattering solver'
      write (output_unit, '(a)') ''
      write (output_unit, '(a)') 'Usage:'
      write (output_unit, '(a)') '  mstm [INPUT_FILE]'
      write (output_unit, '(a)') '  mstm help | --help | -h'
      write (output_unit, '(a)') '  mstm version | --version | -V'
      write (output_unit, '(a)') ''
      write (output_unit, '(a)') 'Arguments:'
      write (output_unit, '(a)') '  INPUT_FILE    Simulation input file.'
      write (output_unit, '(a)') ''
      write (output_unit, '(a)') 'Options:'
      write (output_unit, '(a)') '  -h, --help    Show this help and exit.'
      write (output_unit, '(a)') '  -V, --version Show version information and exit.'
      write (output_unit, '(a)') ''
      write (output_unit, '(a)') 'See docs/manual.md for input parameters and examples.'
   end subroutine print_help

   recursive subroutine execute_nested_loop(looplevel, rank)
      implicit none
      logical :: continueloop
      integer :: looplevel, varposition, loopindex, rank
      integer, pointer :: i_loop_var_pointer
      real(real64) :: maxdif, loopdif
      real(real64), pointer :: r_loop_var_pointer
      complex(real64), pointer :: c_loop_var_pointer
      character(len=256) :: varlabel
      character(len=1) :: vartype

      varlabel = simulation_config%loop_variable_label(looplevel)
      vartype = simulation_config%loop_variable_type(looplevel)
      varposition = simulation_config%loop_sphere_number(looplevel)
      if (vartype .eq. 'i') then
         maxdif = abs(simulation_config%integer_loop_stop(looplevel) - simulation_config%integer_loop_start(looplevel))
         call process_input_variable(varlabel, &
                                     var_position=varposition, &
                                     i_var_pointer=i_loop_var_pointer)
         i_loop_var_pointer = simulation_config%integer_loop_start(looplevel)
      elseif (vartype .eq. 'r') then
         maxdif = abs(simulation_config%real_loop_stop(looplevel) - simulation_config%real_loop_start(looplevel))
         call process_input_variable(varlabel, &
                                     var_position=varposition, &
                                     r_var_pointer=r_loop_var_pointer)
         r_loop_var_pointer = simulation_config%real_loop_start(looplevel)
      elseif (vartype .eq. 'c') then
         maxdif = abs(simulation_config%complex_loop_stop(looplevel) - simulation_config%complex_loop_start(looplevel))
         call process_input_variable(varlabel, &
                                     var_position=varposition, &
                                     c_var_pointer=c_loop_var_pointer)
         c_loop_var_pointer = simulation_config%complex_loop_start(looplevel)
      end if

      loopindex = 0
      continueloop = .true.
      do while (continueloop)
         loopindex = loopindex + 1
         if (looplevel .eq. simulation_config%number_nested_loops) then
            simulation_result%run_number = simulation_result%run_number + 1
!               if(loopindex.eq.1) then
!                  simulation_result%run_number=simulation_result%run_number+1
!               else
!                  if(.not.simulation_config%configuration_average) simulation_result%run_number=simulation_result%run_number+1
!               endif
            if (simulation_config%configuration_average) then
               if (simulation_config%random_orientation) then
                  call run_random_orientation_configuration_average()
               else
                  call run_configuration_average()
               end if
            elseif (simulation_config%incidence_average) then
               call run_incidence_average()
            else
               call execute_simulation()
            end if
            if (runtime_failed()) return
         else
            call execute_nested_loop(looplevel + 1, rank)
            if (runtime_failed()) return
         end if
         if (vartype .eq. 'i') then
            i_loop_var_pointer = i_loop_var_pointer + simulation_config%integer_loop_step(looplevel)
            loopdif = abs(i_loop_var_pointer - simulation_config%integer_loop_start(looplevel))
         elseif (vartype .eq. 'r') then
            r_loop_var_pointer = r_loop_var_pointer + simulation_config%real_loop_step(looplevel)
            loopdif = abs(r_loop_var_pointer - simulation_config%real_loop_start(looplevel))
         elseif (vartype .eq. 'c') then
            c_loop_var_pointer = c_loop_var_pointer + simulation_config%complex_loop_step(looplevel)
            loopdif = abs(c_loop_var_pointer - simulation_config%complex_loop_start(looplevel))
         end if
!            if(loopindex.gt.1000) exit
         if (loopdif - maxdif .gt. 1.d-6) continueloop = .false.
      end do
   end subroutine execute_nested_loop

end program mstm
