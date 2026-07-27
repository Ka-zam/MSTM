module runtime_support
   implicit none
contains

   subroutine timewrite(iunit, char1, time, line_break)
      use intrinsics
      implicit none
      integer :: iunit
      real(8) :: time, time2
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
   end subroutine timewrite
end module runtime_support
