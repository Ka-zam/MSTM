program mstm
   use inputinterface
   use solver
   use mpidefs
   use specialfuncs
   use spheredata
   implicit none
   integer :: looplevel, rank, numprocs, &
              readstat(1), i, numargs, istat, numberinputlines, n
   character(len=256) :: inputfile, inputline, oldoutputfile
   character(len=256), allocatable :: inputfiledata(:)
   data oldoutputfile/' '/

   call mstm_mpi(mpi_command='init')
   call mstm_mpi(mpi_command='rank', mpi_rank=rank)
   call mstm_mpi(mpi_command='size', mpi_size=numprocs)

   do i = 0, numprocs - 1
      if (rank .eq. i) then
         numargs = command_argument_count()
         if (numargs .eq. 0) then
            inputfile = 'mstm.inp'
         else
            call get_command_argument(1, inputfile)
         end if
         input_file = inputfile
         open (2, file=inputfile)
         numberinputlines = 0
         istat = 0
         do while (istat .eq. 0)
            read (2, '(a)', iostat=istat) inputline
            numberinputlines = numberinputlines + 1
            if (trim(inputline) .eq. 'end_of_options') exit
         end do
         if (istat .ne. 0) numberinputlines = numberinputlines + 1
         allocate (inputfiledata(numberinputlines))
         rewind (2)
         istat = 0
         n = 0
         do while (istat .eq. 0)
            read (2, '(a)', iostat=istat) inputline
            n = n + 1
            inputfiledata(n) = inputline
            if (trim(inputline) .eq. 'end_of_options') exit
         end do
         if (istat .ne. 0) inputfiledata(numberinputlines) = 'end_of_options'
         close (2)
      end if
      call mstm_mpi(mpi_command='barrier')
   end do

   repeat_run = .true.
   first_run = .true.

   do while (repeat_run)
      readstat = 0
      do i = 0, numprocs - 1
         if (i .eq. rank) then
            call parse_input_data(inputfiledata, read_status=readstat(1))
         end if
         call mstm_mpi(mpi_command='barrier')
      end do
      if (oldoutputfile .ne. output_file) then
         first_run = .true.
         run_number = 0
         oldoutputfile = output_file
      end if
      if (rank .eq. 0) then
         if (first_run) then
            if (append_output_file) then
               open (2, file=output_file, position='append')
            else
               open (2, file=output_file)
               close (2, status='delete')
               open (2, file=output_file)
            end if
            call output_header(2, inputfile)
            close (2)
         end if
      end if

      if (n_nest_loops .eq. 0) then
         run_number = run_number + 1
         if (configuration_average) then
            if (random_orientation) then
               call run_random_orientation_configuration_average()
            else
               call run_configuration_average()
            end if
         elseif (incidence_average) then
            call run_incidence_average()
         else
            call execute_simulation()
         end if
      else
         looplevel = 1
         call execute_nested_loop(looplevel, rank)
      end if
      n_nest_loops = 0
   end do
!      if(temporary_pos_file.and.rank.eq.0) then
!         open(20,file='temp_pos.dat')
!         close(20,status='delete')
!      endif
   call mstm_mpi(mpi_command='barrier')
   call mstm_mpi(mpi_command='finalize')

contains

   recursive subroutine execute_nested_loop(looplevel, rank)
      implicit none
      logical :: continueloop
      integer :: looplevel, varposition, loopindex, rank
      integer, pointer :: i_loop_var_pointer
      real(8) :: maxdif, loopdif
      real(8), pointer :: r_loop_var_pointer
      complex(8), pointer :: c_loop_var_pointer
      character(len=256) :: varlabel
      character(len=1) :: vartype

      varlabel = loop_var_label(looplevel)
      vartype = loop_var_type(looplevel)
      varposition = loop_sphere_number(looplevel)
      if (vartype .eq. 'i') then
         maxdif = abs(i_var_stop(looplevel) - i_var_start(looplevel))
         call process_input_variable(varlabel, &
                                     var_position=varposition, &
                                     i_var_pointer=i_loop_var_pointer)
         i_loop_var_pointer = i_var_start(looplevel)
      elseif (vartype .eq. 'r') then
         maxdif = abs(r_var_stop(looplevel) - r_var_start(looplevel))
         call process_input_variable(varlabel, &
                                     var_position=varposition, &
                                     r_var_pointer=r_loop_var_pointer)
         r_loop_var_pointer = r_var_start(looplevel)
      elseif (vartype .eq. 'c') then
         maxdif = abs(c_var_stop(looplevel) - c_var_start(looplevel))
         call process_input_variable(varlabel, &
                                     var_position=varposition, &
                                     c_var_pointer=c_loop_var_pointer)
         c_loop_var_pointer = c_var_start(looplevel)
      end if

      loopindex = 0
      continueloop = .true.
      do while (continueloop)
         loopindex = loopindex + 1
         if (looplevel .eq. n_nest_loops) then
            run_number = run_number + 1
!               if(loopindex.eq.1) then
!                  run_number=run_number+1
!               else
!                  if(.not.configuration_average) run_number=run_number+1
!               endif
            if (configuration_average) then
               if (random_orientation) then
                  call run_random_orientation_configuration_average()
               else
                  call run_configuration_average()
               end if
            elseif (incidence_average) then
               call run_incidence_average()
            else
               call execute_simulation()
            end if
         else
            call execute_nested_loop(looplevel + 1, rank)
         end if
         if (vartype .eq. 'i') then
            i_loop_var_pointer = i_loop_var_pointer + i_var_step(looplevel)
            loopdif = abs(i_loop_var_pointer - i_var_start(looplevel))
         elseif (vartype .eq. 'r') then
            r_loop_var_pointer = r_loop_var_pointer + r_var_step(looplevel)
            loopdif = abs(r_loop_var_pointer - r_var_start(looplevel))
         elseif (vartype .eq. 'c') then
            c_loop_var_pointer = c_loop_var_pointer + c_var_step(looplevel)
            loopdif = abs(c_loop_var_pointer - c_var_start(looplevel))
         end if
!            if(loopindex.gt.1000) exit
         if (loopdif - maxdif .gt. 1.d-6) continueloop = .false.
      end do
   end subroutine execute_nested_loop

end program mstm
