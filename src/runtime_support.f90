module runtime_support
   use, intrinsic :: iso_fortran_env, only: real64
   use parallel_runtime, only: mpi_comm_world, mstm_mpi, mstm_mpi_max
   implicit none
   private
   public :: clear_runtime_status, open_input_file, open_output_file, open_update_file, &
             report_runtime_error, runtime_failed, set_runtime_error, synchronize_runtime_status, &
             write_elapsed_time

   integer, save :: runtime_status_code = 0
   character(len=512), save :: runtime_status_message = ''
contains

   subroutine clear_runtime_status()
      runtime_status_code = 0
      runtime_status_message = ''
   end subroutine clear_runtime_status

   logical function runtime_failed()
      runtime_failed = runtime_status_code /= 0
   end function runtime_failed

   subroutine set_runtime_error(message, status_code)
      character(len=*), intent(in) :: message
      integer, intent(in), optional :: status_code

      if (runtime_status_code /= 0) return
      runtime_status_code = 1
      if (present(status_code)) runtime_status_code = max(1, abs(status_code))
      runtime_status_message = message
   end subroutine set_runtime_error

   subroutine report_runtime_error(unit)
      integer, intent(in) :: unit

      if (runtime_status_code == 0) return
      if (len_trim(runtime_status_message) == 0) then
         write (unit, '(a,i0)') 'MSTM failed with status ', runtime_status_code
      else
         write (unit, '(a)') trim(runtime_status_message)
      end if
   end subroutine report_runtime_error

   subroutine synchronize_runtime_status(mpi_comm)
      integer, intent(in), optional :: mpi_comm
      integer :: communicator, global_status(1), local_status(1)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      local_status(1) = runtime_status_code
      global_status = local_status
      call mstm_mpi(mpi_command='allreduce', mpi_send_buf_i=local_status, &
                    mpi_recv_buf_i=global_status, mpi_number=1, mpi_operation=mstm_mpi_max, &
                    mpi_comm=communicator)
      runtime_status_code = global_status(1)
      if (global_status(1) /= 0 .and. len_trim(runtime_status_message) == 0) &
         runtime_status_message = 'MSTM failed on another parallel rank'
   end subroutine synchronize_runtime_status

   subroutine open_input_file(path, unit)
      character(len=*), intent(in) :: path
      integer, intent(out) :: unit
      integer :: io_status
      character(len=256) :: io_message

      open (newunit=unit, file=trim(path), status='old', action='read', iostat=io_status, iomsg=io_message)
      if (io_status /= 0) call set_runtime_error( &
         "Cannot open input file '"//trim(path)//"': "//trim(io_message), io_status)
   end subroutine open_input_file

   subroutine open_output_file(path, unit, append)
      character(len=*), intent(in) :: path
      integer, intent(out) :: unit
      logical, intent(in), optional :: append
      integer :: io_status
      character(len=256) :: io_message

      if (present(append)) then
         if (append) then
            open (newunit=unit, file=trim(path), status='unknown', action='write', position='append', &
                  iostat=io_status, iomsg=io_message)
         else
            open (newunit=unit, file=trim(path), status='replace', action='write', &
                  iostat=io_status, iomsg=io_message)
         end if
      else
         open (newunit=unit, file=trim(path), status='replace', action='write', &
               iostat=io_status, iomsg=io_message)
      end if
      if (io_status /= 0) call set_runtime_error( &
         "Cannot open output file '"//trim(path)//"': "//trim(io_message), io_status)
   end subroutine open_output_file

   subroutine open_update_file(path, unit)
      character(len=*), intent(in) :: path
      integer, intent(out) :: unit
      integer :: io_status
      character(len=256) :: io_message

      open (newunit=unit, file=trim(path), status='old', action='readwrite', iostat=io_status, iomsg=io_message)
      if (io_status /= 0) call set_runtime_error( &
         "Cannot update output file '"//trim(path)//"': "//trim(io_message), io_status)
   end subroutine open_update_file

   subroutine write_elapsed_time(iunit, char1, time, line_break)
      implicit none
      integer :: iunit
      real(real64) :: time, time2
      logical :: linebreak
      logical, optional :: line_break
      character(*) :: char1
      if (present(line_break)) then
         linebreak = line_break
      else
         linebreak = .true.
      end if
      if (time .gt. 3600.d0) then
         time2 = time / 3600.d0
         if (linebreak) then
            write (iunit, '(a,f9.3,'' hours'')') char1, time2
         else
            write (iunit, '(a,f9.3,'' hours'')', advance='no') char1, time2
         end if
      elseif (time .gt. 60.d0) then
         time2 = time / 60.d0
         if (linebreak) then
            write (iunit, '(a,f9.2,'' min'')') char1, time2
         else
            write (iunit, '(a,f9.2,'' min'')', advance='no') char1, time2
         end if
      else
         if (linebreak) then
            write (iunit, '(a,f9.2,'' sec'')') char1, time
         else
            write (iunit, '(a,f9.2,'' sec'')', advance='no') char1, time
         end if
      end if
      if (linebreak) flush (iunit)
   end subroutine write_elapsed_time
end module runtime_support
