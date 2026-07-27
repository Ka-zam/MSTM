module intrinsics
!
! Portable wrappers around standard Fortran intrinsic procedures.
!
!
!  last revised: 15 January 2011
!
   implicit none
contains
!
!   system clock
!
   real function mytime()
      implicit none
      call cpu_time(mytime)
   end function mytime
!
!   number of command-line arguments.
!
   integer function mstm_nargs()
      implicit none
      mstm_nargs = command_argument_count()
   end function mstm_nargs
!
!   command line argument retrieval
!
   subroutine mstm_getarg(char)
      implicit none
      character(*) :: char
      call get_command_argument(1, char)
   end subroutine mstm_getarg
end module intrinsics
