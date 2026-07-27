module random_configuration_geometry
   use constants
   use random_configuration_sorting, only: heap_sort_with_tolerance
   use random_configuration_state
   implicit none
   private
   public :: add_sphere_to_cluster, calculate_target_distribution_statistics, calculate_target_volume, &
             check_position_in_target, circumscribing_sphere, clear_cells, direct_overlap_test, &
             generate_hexagonal_positions, initialize_cells, modify_cells, position_to_cell_index, &
             sample_layered_configuration, sample_target_position, sort_sphere_positions, sort_sphere_radii, &
             swap_cell_contents
contains

   pure subroutine direct_overlap_test(nsphere, radius, position, overlap, distance, pair)
      implicit none
      integer, intent(in) :: nsphere
      real(8), intent(in) :: radius(nsphere), position(3, nsphere)
      logical, intent(out) :: overlap
      real(8), optional, intent(out) :: distance
      integer, optional, intent(out) :: pair(2)
      integer :: i, j
      real(8) :: rij
      overlap = .false.
      do i = 1, nsphere - 1
         do j = i + 1, nsphere
            rij = sqrt(sum((position(:, i) - position(:, j))**2))
            if (rij .lt. radius(i) + radius(j)) then
               overlap = .true.
               if (present(distance)) distance = rij
               if (present(pair)) pair = (/i, j/)
               return
            end if
         end do
      end do
   end subroutine direct_overlap_test

   subroutine calculate_target_volume(targetdimensions, targetvol)
      implicit none
      integer :: i, ipbc(3)
      real(8) :: targetvol, targetdimensions(3)
      ipbc = wall_boundary_model
      do i = 1, 3
         if (periodic_bc(i)) ipbc(i) = 0
      end do
      if (target_shape .eq. 0) then
         targetvol = 8.d0 * product(targetdimensions(1:3) - dble(ipbc))
      elseif (target_shape .eq. 1) then
         targetvol = two_pi * ((targetdimensions(1) - wall_boundary_model)**2) * (targetdimensions(3) - ipbc(3))
      else
         targetvol = four_pi_over_three * (targetdimensions(1) - wall_boundary_model)**3
      end if
   end subroutine calculate_target_volume

   subroutine position_to_cell_index(pos, cell)
      implicit none
      integer :: cell(3)
      real(8) :: pos(3)
      cell = floor((pos(:) - target_boundaries(:, 1)) / (target_boundaries(:, 2) - target_boundaries(:, 1)) * dble(cell_dim(:))) + 1
      cell = max(cell, (/1, 1, 1/))
      cell = min(cell, cell_dim)
   end subroutine position_to_cell_index

   subroutine sample_target_position(pos, rad)
      implicit none
      integer :: i
      real(8) :: pos(3), rannum(3), r, phi, ct, st, rad, wshift(3)

      call random_number(rannum)
      if (target_shape .eq. 0) then
         do i = 1, 3
            if (periodic_bc(i)) then
               wshift(i) = 0.d0
            else
               wshift(i) = rad * wall_boundary_model + minimum_gap
            end if
         end do
      pos = target_boundaries(:, 1) + wshift(:) + (target_boundaries(:, 2) - target_boundaries(:, 1) - 2.d0 * wshift(:)) * rannum(:)
      elseif (target_shape .eq. 1) then
         wshift(1) = rad * wall_boundary_model + minimum_gap
         if (periodic_bc(3)) then
            wshift(3) = 0.d0
         else
            wshift(3) = rad * wall_boundary_model + minimum_gap
         end if
         r = (target_boundaries(1, 2) - wshift(1)) * rannum(1)**0.5d0
         phi = two_pi * rannum(2)
         pos(1) = r * cos(phi)
         pos(2) = r * sin(phi)
   pos(3) = target_boundaries(3, 1) + wshift(3) + (target_boundaries(3, 2) - target_boundaries(3, 1) - 2.d0 * wshift(3)) * rannum(3)
      else
         wshift(1) = rad * wall_boundary_model + minimum_gap
         r = (target_boundaries(1, 2) - wshift(1)) * rannum(1)**0.333333d0
         phi = two_pi * rannum(2)
         ct = -1.d0 + 2.d0 * rannum(3)
         st = sqrt(1.d0 - ct * ct)
         pos(1) = r * st * cos(phi)
         pos(2) = r * st * sin(phi)
         pos(3) = r * ct
      end if
   end subroutine sample_target_position

   subroutine clear_cells()
      implicit none
      integer :: n, ix, iy, iz, i
      type(l_list), pointer :: llist, llist2
      if (allocated(cell_list)) then
         do iz = 1, cell_dim(3)
            do iy = 1, cell_dim(2)
               do ix = 1, cell_dim(1)
                  n = cell_list(ix, iy, iz)%number_elements
                  if (.not. associated(cell_list(ix, iy, iz)%members)) cycle
                  llist => cell_list(ix, iy, iz)%members
                  do i = 1, n
                     llist2 => llist%next
                     deallocate (llist)
                     nullify (llist)
                     if (.not. associated(llist2)) exit
                     llist => llist2
                  end do
               end do
            end do
         end do
         deallocate (cell_list)
      end if
      if (allocated(sphere_cell)) deallocate (sphere_cell)
   end subroutine clear_cells

   subroutine initialize_cells(nsphere)
      implicit none
      integer :: nsphere
      if (allocated(sphere_cell)) deallocate (sphere_cell)
      allocate (sphere_cell(3, nsphere))
      sphere_cell(:, :) = 0
      cell_dim(:) = floor((target_boundaries(:, 2) - target_boundaries(:, 1) - 1.d-6) / d_cell) + 1
      if (allocated(cell_list)) deallocate (cell_list)
      allocate (cell_list(cell_dim(1), cell_dim(2), cell_dim(3)))
      cell_list(:, :, :)%number_elements = 0
   end subroutine initialize_cells

   subroutine swap_cell_contents(i, newcell)
      implicit none
      integer :: i, newcell(3), cell(3), n, l
      type(l_list), pointer :: llist, llist2, llistnew
      cell(:) = sphere_cell(:, i)
      n = cell_list(cell(1), cell(2), cell(3))%number_elements
      if (cell_list(cell(1), cell(2), cell(3))%members%index .eq. i) then
         llist2 => cell_list(cell(1), cell(2), cell(3))%members%next
         llistnew => cell_list(cell(1), cell(2), cell(3))%members
         cell_list(cell(1), cell(2), cell(3))%members => llist2
      else
         llist => cell_list(cell(1), cell(2), cell(3))%members
         do l = 1, n - 1
            if (llist%next%index .eq. i) then
               llist2 => llist%next%next
               llistnew => llist%next
               llist%next => llist2
               exit
            end if
            llist => llist%next
         end do
      end if
      cell_list(cell(1), cell(2), cell(3))%number_elements = n - 1
      cell = newcell
      n = cell_list(cell(1), cell(2), cell(3))%number_elements
      llist2 => cell_list(cell(1), cell(2), cell(3))%members
      cell_list(cell(1), cell(2), cell(3))%members => llistnew
      cell_list(cell(1), cell(2), cell(3))%members%next => llist2
      cell_list(cell(1), cell(2), cell(3))%number_elements = n + 1
      sphere_cell(:, i) = cell
   end subroutine swap_cell_contents

   subroutine modify_cells(nsphere, position, start_sphere, end_sphere)
      implicit none
      integer :: isphere, nsphere, cell(3), istart, iend
      integer, optional :: start_sphere, end_sphere
      real(8) :: position(3, nsphere)
      if (present(start_sphere)) then
         istart = start_sphere
      else
         istart = 1
      end if
      if (present(end_sphere)) then
         iend = end_sphere
      else
         iend = nsphere
      end if
      do isphere = istart, iend
         call position_to_cell_index(position(:, isphere), cell)
         if (any(sphere_cell(:, isphere) .ne. cell)) then
            call swap_cell_contents(isphere, cell)
         end if
      end do
   end subroutine modify_cells

   subroutine calculate_target_distribution_statistics(nsphere, sdev)
      implicit none
      integer :: nsphere, n, iz, iy, ix, nt, ncell
      real(8) :: sdev, nmean
      sdev = 0.d0
      ncell = product(cell_dim)
      nmean = dble(nsphere) / dble(ncell)
      nt = 0
      do iz = 1, cell_dim(3)
         do iy = 1, cell_dim(2)
            do ix = 1, cell_dim(1)
               n = cell_list(ix, iy, iz)%number_elements
               nt = nt + n
               sdev = sdev + (dble(n) / nmean - 1.d0) * (dble(n) / nmean - 1.d0)
            end do
         end do
      end do
      sdev = sqrt(sdev)
   end subroutine calculate_target_distribution_statistics

   subroutine add_sphere_to_cluster(newrad, newpos, nsphere, radius, position, fitok)
      implicit none
      logical :: fitok, bndok, pbc(3)
      integer :: nsphere, i, cell(3), n, m, ccell(3), scell(3), j
      real(8) :: radius(*), position(3, *), newrad, newpos(3), rij, tpos(3)
      type(l_list), pointer :: llist
      pbc = .false.
      if (target_shape .eq. 0) then
         pbc = periodic_bc
      elseif (target_shape .eq. 1) then
         pbc(3) = periodic_bc(3)
      end if
      call position_to_cell_index(newpos, ccell)
      fitok = .true.
      do m = 0, 26
         scell(1) = mod(m, 3) - 1
         scell(2) = mod(m / 3, 3) - 1
         scell(3) = mod(m / 9, 3) - 1
         cell = ccell + scell
         bndok = .true.
         tpos = newpos
         do i = 1, 3
            if (cell(i) .lt. 1 .or. cell(i) .gt. cell_dim(i)) then
               if (pbc(i)) then
                  if (cell(i) .lt. 1) then
                     cell(i) = cell_dim(i)
                     tpos(i) = tpos(i) + target_boundaries(i, 2) - target_boundaries(i, 1)
                  elseif (cell(i) .gt. cell_dim(i)) then
                     cell(i) = 1
                     tpos(i) = tpos(i) - target_boundaries(i, 2) + target_boundaries(i, 1)
                  end if
               else
                  bndok = .false.
                  exit
               end if
            end if
         end do
         if (.not. bndok) cycle
         n = cell_list(cell(1), cell(2), cell(3))%number_elements
         if (n .eq. 0) cycle
         llist => cell_list(cell(1), cell(2), cell(3))%members
         do j = 1, n
            i = llist%index
            rij = sqrt(sum((tpos(:) - position(:, i))**2))
            if (rij .lt. newrad + radius(i) + minimum_gap) then
               fitok = .false.
               return
            end if
            if (j .lt. n) llist => llist%next
         end do
      end do
      sphere_cell(:, nsphere + 1) = ccell(:)
      n = cell_list(ccell(1), ccell(2), ccell(3))%number_elements
      if (n .eq. 0) allocate (cell_list(ccell(1), ccell(2), ccell(3))%members)
      llist => cell_list(ccell(1), ccell(2), ccell(3))%members
      do i = 1, n
         if (i .eq. n) allocate (llist%next)
         llist => llist%next
      end do
      llist%index = nsphere + 1
      cell_list(ccell(1), ccell(2), ccell(3))%number_elements = n + 1
   end subroutine add_sphere_to_cluster

   subroutine sort_sphere_positions(nsphere, radius, position, cindex, sort_elem, make_positive)
      implicit none
      logical :: makepos
      logical, optional :: make_positive
      integer :: nsphere, i, ind(nsphere), selem, cindex(nsphere), tindex(nsphere)
      integer, optional :: sort_elem
      real(8) :: radius(nsphere), position(3, nsphere), r(nsphere), &
                 tpos(3, nsphere)
      if (present(sort_elem)) then
         selem = sort_elem
      else
         selem = 0
      end if
      if (present(make_positive)) then
         makepos = make_positive
      else
         makepos = .false.
      end if
      if (selem .eq. 0) then
         r(:) = sqrt(sum(position(:, :)**2, 1))
      else
         if (makepos) then
            r(:) = abs(position(selem, :))
         else
            r(:) = position(selem, :)
         end if
      end if
      ind(1) = 0
      call heap_sort_with_tolerance(nsphere, r, ind, 1.d-15)
      r = radius
      tpos = position
      tindex = cindex
      do i = 1, nsphere
         radius(i) = r(ind(i))
         position(:, i) = tpos(:, ind(i))
         cindex(i) = tindex(ind(i))
      end do
   end subroutine sort_sphere_positions

   subroutine sort_sphere_radii(nsphere, radius)
      implicit none
      integer :: nsphere, ind(nsphere)
      real(8) :: radius(nsphere)
      radius = -radius
      ind(1) = 0
      call heap_sort_with_tolerance(nsphere, radius, ind, 1.d-15)
      radius = -radius
   end subroutine sort_sphere_radii

   pure subroutine circumscribing_sphere(nsphere, radius, position, rcell)
      implicit none
      integer, intent(in) :: nsphere
      real(8), intent(in) :: radius(nsphere), position(3, nsphere)
      real(8), intent(out) :: rcell
      integer :: i
      rcell = 0.d0
      do concurrent(i=1:nsphere) reduce(max:rcell)
         rcell = max(rcell, sqrt(sum(position(:, i)**2)) + radius(i))
      end do
   end subroutine circumscribing_sphere

   subroutine check_position_in_target(rad, pos, wallbound, intarget)
      implicit none
      logical :: intarget
      integer :: i
      real(8) :: rad, pos(3), wallbound(3, 2), rho, wrad
      intarget = .true.
      wrad = rad * wall_boundary_model
      if (target_shape .eq. 0) then
         do i = 1, 3
            if (periodic_bc(i)) then
               intarget = (pos(i) .ge. wallbound(i, 1) .and. pos(i) .le. wallbound(i, 2))
            else
               intarget = (pos(i) - wrad .ge. wallbound(i, 1) .and. pos(i) + wrad .le. wallbound(i, 2))
            end if
            if (.not. intarget) return
         end do
      elseif (target_shape .eq. 1) then
         rho = sqrt(sum(pos(1:2)**2))
         if (rho + wrad .ge. wallbound(1, 2)) then
            intarget = .false.
            return
         end if
         if (periodic_bc(3)) then
            intarget = (pos(3) .ge. wallbound(3, 1) .and. pos(3) .le. wallbound(3, 2))
         else
            intarget = (pos(3) - wrad .ge. wallbound(3, 1) .and. pos(3) + wrad .le. wallbound(3, 2))
         end if
         if (.not. intarget) return
      else
         rho = sqrt(sum(pos(1:3)**2))
         if (rho + wrad .gt. wallbound(1, 2)) then
            intarget = .false.
            return
         end if
      end if
   end subroutine check_position_in_target

   subroutine sample_layered_configuration(nsphere, rad, pos, wallbound, nin)
      implicit none
      logical :: fitok, pbc(3)
      integer :: nsphere, nin, i, m, maxsamp
      real(8) :: rad(nsphere), pos(3, nsphere), wallbound(3, 2), r2, vtot, delv, wdist(3), &
                 dz, z1, z2, samp(3), rho, phi, r1, r, st, ct, rannum(3)
      data maxsamp/5000/
      if (target_shape .eq. 0) then
         pbc = periodic_bc
      elseif (target_shape .eq. 1) then
         pbc(1:2) = .false.
         pbc(3) = periodic_bc(3)
      else
         pbc = .false.
      end if
      if (target_shape .eq. 2) then
         r2 = 0.d0
         vtot = 0.d0
         delv = four_pi_over_three * (wallbound(1, 2) - dble(wall_boundary_model) - minimum_gap)**3 / dble(nsphere)
      end if
      nin = 0
      do i = 1, nsphere
         wdist = 0.d0
         do m = 1, 3
            if (.not. pbc(m)) wdist(m) = rad(i) * wall_boundary_model + minimum_gap
         end do
         if (target_shape .eq. 0) then
            dz = (wallbound(3, 2) - wallbound(3, 1) - 2.d0 * wdist(3)) / dble(nsphere)
            z1 = wallbound(3, 1) + wdist(3) + dble(i - 1) * dz
            z2 = z1 + dz
            z2 = min(z2, wallbound(3, 2) - wdist(3))
            do m = 1, maxsamp
               call random_number(rannum)
               samp = (/wallbound(1, 1) + wdist(1), wallbound(2, 1) + wdist(2), z1/) &
                      + ((/wallbound(1, 2) - wdist(1), wallbound(2, 2) - wdist(2), z2/) &
                         - (/wallbound(1, 1) + wdist(1), wallbound(2, 1) + wdist(2), z1/)) * rannum
               call check_position_in_target(rad(i) * wall_boundary_model, samp, wallbound, fitok)
               if (.not. fitok) cycle
               call add_sphere_to_cluster(rad(i), samp, i - 1, rad, pos, fitok)
               if (fitok) then
                  pos(:, i) = samp(:)
                  exit
               end if
            end do
            if (.not. fitok) return
         elseif (target_shape .eq. 1) then
            dz = (wallbound(3, 2) - wallbound(3, 1) - 2.d0 * wdist(3)) / dble(nsphere)
            z1 = wallbound(3, 1) + wdist(3) + dble(i - 1) * dz
            z2 = z1 + dz
            do m = 1, maxsamp
               call random_number(rannum)
               rho = (wallbound(1, 2) - wdist(1)) * sqrt(rannum(1))
               phi = two_pi * rannum(2)
               samp(1) = rho * cos(phi)
               samp(2) = rho * sin(phi)
               samp(3) = z1 + dz * rannum(3)
               call check_position_in_target(rad(i) * wall_boundary_model, samp, wallbound, fitok)
               if (.not. fitok) cycle
               call add_sphere_to_cluster(rad(i), samp, i - 1, rad, pos, fitok)
               if (fitok) then
                  pos(:, i) = samp(:)
                  exit
               end if
            end do
            if (.not. fitok) return
         else
            r1 = r2
            vtot = vtot + delv
            r2 = (vtot / four_pi_over_three)**(1.d0 / 3.d0)
            do m = 1, maxsamp
               if (sphere_1_fixed .and. i .eq. 1) then
                  pos = 0.d0
                  fitok = .true.
               else
                  call random_number(rannum)
                  if ((wallbound(1, 2) - wdist(1)) .le. r2) then
                     r2 = wallbound(1, 2) - wdist(1)
                     r2 = max(r2, 0.d0)
                  end if
                  r = (delv / four_pi_over_three * rannum(1) + r1**3)**(1.d0 / 3.d0)
                  ct = -1.d0 + 2.d0 * rannum(2)
                  st = sqrt(1.d0 - ct * ct)
                  phi = two_pi * rannum(3)
                  samp(1) = r * st * cos(phi)
                  samp(2) = r * st * sin(phi)
                  samp(3) = r * ct
                  call check_position_in_target(rad(i) * wall_boundary_model, samp, wallbound, fitok)
               end if
               if (.not. fitok) cycle
               call add_sphere_to_cluster(rad(i), samp, i - 1, rad, pos, fitok)
               if (fitok) then
                  pos(:, i) = samp(:)
                  exit
               end if
            end do
            if (.not. fitok) return
         end if
         nin = i
      end do
   end subroutine sample_layered_configuration

   subroutine generate_hexagonal_positions(nsphere, rad, pos, wallbound, s, allin, ns)
      implicit none
      logical :: intarget, fitok, allin
      integer :: nsphere, i, l, m, n, ns, l0, imax, m0, i2, n2, m2, i21, l1, n0, ns0
      real(8) :: rad(nsphere), pos(3, nsphere), wallbound(3, 2), s, cscale(3), tpos(3), &
                 tpos1(3), tpos2(3), trad, cscale2(3)
      data imax/200/
      cscale = (/2.d0, sqrt(3.d0), sqrt(8.d0 / 3.d0)/)
      cscale2 = cscale * cscale
      ns = 0
      allin = .true.
      do i = 0, imax
         if (mod(i, 2) .eq. 0) ns0 = ns
         i2 = i * i
         i21 = (i + 1) * (i + 1)
         n0 = ceiling(dble(i) / cscale(3))
         do n = -n0 - 1, n0 + 1
            n2 = n * n
            m0 = ceiling(sqrt(max(dble(i2) - cscale2(3) * dble(n2), 0.d0)) / cscale(2))
            do m = -m0 - 1, m0 + 1
               m2 = m * m
               l0 = floor(sqrt(max(dble(i2) - cscale2(3) * dble(n2) - cscale2(2) * dble(m2), 0.d0)) / cscale(1))
               l0 = max(l0, 1)
               l1 = ceiling(sqrt(max(dble(i21) - cscale2(3) * dble(n2) - cscale2(2) * dble(m2), 0.d0)) / cscale(1))
               tpos1 = s * (/dble(mod(abs(m + n), 2)), dble(mod(abs(n), 2)) / cscale(2), 0.d0/)
               do l = l0 - 1, l1 + 1
                  tpos2 = s * cscale * (/dble(l), dble(m), dble(n)/)
                  tpos = tpos1 + tpos2
                  trad = sqrt(sum(tpos**2)) / s
                  if (trad .ge. dble(i) .and. trad .lt. dble(i + 1)) then
                     trad = rad(ns + 1) * wall_boundary_model
                     intarget = .true.
                     call check_position_in_target(trad, tpos, wallbound, intarget)
                     if (intarget) then
                        call add_sphere_to_cluster(rad(ns + 1), tpos, ns, rad, pos, fitok)
                        if (fitok) then
                           ns = ns + 1
                           pos(:, ns) = tpos(:)
                           if (ns .eq. nsphere) return
                        end if
                     end if
                  end if
                  if (l .eq. 0) cycle
                  tpos2 = s * cscale * (/-dble(l), dble(m), dble(n)/)
                  tpos = tpos1 + tpos2
                  trad = sqrt(sum(tpos**2)) / s
                  if (trad .ge. dble(i) .and. trad .lt. dble(i + 1)) then
                     trad = rad(ns + 1) * wall_boundary_model
                     intarget = .true.
                     call check_position_in_target(trad, tpos, wallbound, intarget)
                     if (intarget) then
                        call add_sphere_to_cluster(rad(ns + 1), tpos, ns, rad, pos, fitok)
                        if (fitok) then
                           ns = ns + 1
                           pos(:, ns) = tpos(:)
                           if (ns .eq. nsphere) return
                        end if
                     end if
                  end if
               end do
            end do
         end do
         if (ns .eq. ns0 .and. mod(i, 2) .ne. 0) then
            if (ns .lt. nsphere) allin = .false.
            return
         end if
      end do
   end subroutine generate_hexagonal_positions
end module random_configuration_geometry
