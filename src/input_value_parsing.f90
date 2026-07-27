module input_value_parsing
   implicit none
   private
   public :: set_string_to_cmplx_variable, set_string_to_int_variable, &
             set_string_to_logical_array_variable, set_string_to_logical_variable, &
             set_string_to_real_array_variable, set_string_to_real_variable
contains

   subroutine set_string_to_int_variable(sentvarvalue, &
                                         ivarvalue, var_operation)
      implicit none
      integer :: itemp
      integer, pointer :: ivarvalue
      character(len=256) :: sentvarvalue, varop, intfile
      character(len=*), optional :: var_operation
      if (present(var_operation)) then
         varop = var_operation(:index(var_operation, ' '))
      else
         varop = 'assign'
      end if
      write (intfile, '(a)') sentvarvalue
      read (intfile, *) itemp
      if (varop(1:6) .eq. 'assign') then
         ivarvalue = itemp
      elseif (varop(1:3) .eq. 'add') then
         ivarvalue = ivarvalue + itemp
      end if
   end subroutine set_string_to_int_variable

   subroutine set_string_to_real_variable(sentvarvalue, &
                                          rvarvalue, var_operation)
      implicit none
      real(8) :: rtemp
      real(8), pointer :: rvarvalue
      character(len=256) :: sentvarvalue, varop, intfile
      character(len=256), optional :: var_operation
      if (present(var_operation)) then
         varop = var_operation(:index(var_operation, ' '))
      else
         varop = 'assign'
      end if
      write (intfile, '(a)') sentvarvalue
      read (intfile, *) rtemp
      if (varop(1:6) .eq. 'assign') then
         rvarvalue = rtemp
      elseif (varop(1:3) .eq. 'add') then
         rvarvalue = rvarvalue + rtemp
      end if
   end subroutine set_string_to_real_variable

   subroutine set_string_to_real_array_variable(sentvarvalue, &
                                                rvarvalue, var_operation, var_len)
      implicit none
      integer :: varlen, i, ierr
      integer, optional :: var_len
      real(8) :: rtemp(4)
      real(8), pointer :: rvarvalue(:)
      character(len=256) :: sentvarvalue, varop, intfile
      character(len=256), optional :: var_operation
      if (present(var_operation)) then
         varop = var_operation(:index(var_operation, ' '))
      else
         varop = 'assign'
      end if
      if (present(var_len)) then
         varlen = var_len
      else
         varlen = 1
      end if
      write (intfile, '(a)') sentvarvalue
      do i = 1, varlen
         read (intfile, *, iostat=ierr) rtemp(1:i)
         if (ierr .ne. 0) then
            rtemp(i:varlen) = rtemp(i - 1)
            exit
         end if
      end do
!         read(intfile,*) rtemp(1:varlen)
      if (varop(1:6) .eq. 'assign') then
         rvarvalue = rtemp(1:varlen)
      elseif (varop(1:3) .eq. 'add') then
         rvarvalue = rvarvalue + rtemp
      end if
   end subroutine set_string_to_real_array_variable

   subroutine set_string_to_cmplx_variable(sentvarvalue, &
                                           cvarvalue, var_operation)
      implicit none
      complex(8) :: ctemp
      complex(8), pointer :: cvarvalue
      character(len=256) :: sentvarvalue, varop, intfile
      character(len=256), optional :: var_operation
      if (present(var_operation)) then
         varop = var_operation(:index(var_operation, ' '))
      else
         varop = 'assign'
      end if
      write (intfile, '(a)') sentvarvalue
      read (intfile, *) ctemp
      if (varop(1:6) .eq. 'assign') then
         cvarvalue = ctemp
      elseif (varop(1:3) .eq. 'add') then
         cvarvalue = cvarvalue + ctemp
      end if
   end subroutine set_string_to_cmplx_variable

   subroutine set_string_to_logical_variable(sentvarvalue, &
                                             lvarvalue, var_operation)
      implicit none
      logical :: ltemp
      logical, pointer :: lvarvalue
      character(len=256) :: sentvarvalue, varop, intfile
      character(len=256), optional :: var_operation
      if (present(var_operation)) then
         varop = var_operation(:index(var_operation, ' '))
      else
         varop = 'assign'
      end if
      write (intfile, '(a)') sentvarvalue
      read (intfile, *) ltemp
      if (varop(1:6) .eq. 'assign') then
         lvarvalue = ltemp
      end if
   end subroutine set_string_to_logical_variable

   subroutine set_string_to_logical_array_variable(sentvarvalue, &
                                                   lvarvalue, var_operation, var_len)
      implicit none
      logical :: ltemp(5)
      logical, pointer :: lvarvalue(:)
      integer :: i, varlen, ierr
      integer, optional :: var_len
      character(len=256) :: sentvarvalue, varop, intfile
      character(len=256), optional :: var_operation
      if (present(var_operation)) then
         varop = var_operation(:index(var_operation, ' '))
      else
         varop = 'assign'
      end if
      if (present(var_len)) then
         varlen = var_len
      else
         varlen = 1
      end if
      write (intfile, '(a)') sentvarvalue
      do i = 1, varlen
         read (intfile, *, iostat=ierr) ltemp(1:i)
         if (ierr .ne. 0) then
            ltemp(i:varlen) = ltemp(i - 1)
            exit
         end if
      end do
      if (varop(1:6) .eq. 'assign') then
         lvarvalue = ltemp(1:varlen)
      end if
   end subroutine set_string_to_logical_array_variable
end module input_value_parsing
