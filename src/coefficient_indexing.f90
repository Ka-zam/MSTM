module coefficient_indexing
   implicit none
contains

   pure integer function mode_index(m, n, l, model)
      implicit none
      integer, intent(in) :: m, n, l, model
      if (model .eq. 1) then
         mode_index = n * (n + 1) + m
      elseif (m .ge. 0) then
         mode_index = (n - 1) * (l + 2) + m + 1
      else
         mode_index = -(m + 1) * (l + 2) + n + 2
      end if
   end function mode_index

   pure integer function polarized_mode_index(m, n, p, l, model)
      implicit none
      integer, intent(in) :: m, n, p, l, model
      if (model .eq. 1) then
         polarized_mode_index = 2 * (n * (n + 1) + m - 1) + p
      elseif (m .ge. 0) then
         polarized_mode_index = (n - 1) * (l + 2) + m + 1 + (p - 1) * l * (l + 2)
      else
         polarized_mode_index = -(m + 1) * (l + 2) + n + 2 + (p - 1) * l * (l + 2)
      end if
   end function polarized_mode_index
end module coefficient_indexing
