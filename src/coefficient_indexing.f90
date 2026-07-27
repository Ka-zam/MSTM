module coefficient_indexing
   implicit none
contains

   pure integer function amnaddress(m, n, l, model)
      implicit none
      integer, intent(in) :: m, n, l, model
      if (model .eq. 1) then
         amnaddress = n * (n + 1) + m
      elseif (m .ge. 0) then
         amnaddress = (n - 1) * (l + 2) + m + 1
      else
         amnaddress = -(m + 1) * (l + 2) + n + 2
      end if
   end function amnaddress

   pure integer function amnpaddress(m, n, p, l, model)
      implicit none
      integer, intent(in) :: m, n, p, l, model
      if (model .eq. 1) then
         amnpaddress = 2 * (n * (n + 1) + m - 1) + p
      elseif (m .ge. 0) then
         amnpaddress = (n - 1) * (l + 2) + m + 1 + (p - 1) * l * (l + 2)
      else
         amnpaddress = -(m + 1) * (l + 2) + n + 2 + (p - 1) * l * (l + 2)
      end if
   end function amnpaddress
end module coefficient_indexing
