module random_configuration_dynamics
   use constants
   use mpidefs, only: mstm_mpi_wtime
   use random_configuration_geometry, only: check_position_in_target, modify_cells, position_to_cell_index
   use random_configuration_state
   implicit none
   private
   public :: find_next_sphere_collision, find_next_wall_collision, move_spheres, &
             resolve_collision_velocities, sample_particle_radius, sample_random_velocities, test_pair_collision
contains

   subroutine move_spheres(nsphere, radius, pos, u, maxtime, wallboundaries, number_wall_hits)
      implicit none
      logical :: collision, wallcollision, pbc(3), intarget
      integer :: nsphere, i, is, js, collisionpair(2), iwall, iswall, m, nwhits, it
      integer, optional :: number_wall_hits
      real(8) :: pos(3, nsphere), radius(nsphere), maxtime, &
                 tcmin, tmove, u1new(3), u2new(3), u(1:3, nsphere), &
                 twallmin, rho, cp, sp, urho, uphi, tpos(3), &
                 u1pn(3), ct, st, r, tcoll, wallboundaries(3, 2), collpos(3)
      pbc = .false.
      if (target_shape .eq. 0) then
         pbc = periodic_bc
      elseif (target_shape .eq. 1) then
         pbc(3) = periodic_bc(3)
      end if
      tmove = maxtime
      i = 1
      nwhits = 0
      time_0 = mstm_mpi_wtime()
      call find_next_sphere_collision(nsphere, radius, pos, u, tmove, wallboundaries, &
                                      collision, tcoll, collisionpair)
      sim_timings(1) = sim_timings(1) + mstm_mpi_wtime() - time_0
      do while (tmove .gt. 0.d0)
         time_0 = mstm_mpi_wtime()
         call modify_cells(nsphere, pos)
         sim_timings(2) = sim_timings(2) + mstm_mpi_wtime() - time_0
         time_0 = mstm_mpi_wtime()
         tcoll = tmove
         collision = .false.
         do is = 1, nsphere
            if (coll_data(is)%time .lt. tcoll) then
               tcoll = coll_data(is)%time
               collision = .true.
               collisionpair(1) = is
               collisionpair(2) = coll_data(is)%sphere
               collpos(:) = coll_data(is)%collpos(:)
            end if
         end do
         sim_timings(3) = sim_timings(3) + mstm_mpi_wtime() - time_0
         time_0 = mstm_mpi_wtime()
         tcmin = tcoll
         call find_next_wall_collision(nsphere, radius, pos, u, tmove, wallboundaries, &
                                       twallmin, iswall, iwall)
         sim_timings(4) = sim_timings(4) + mstm_mpi_wtime() - time_0
         time_0 = mstm_mpi_wtime()
         wallcollision = (twallmin .lt. tcmin)
         tcmin = min(tcmin, twallmin)
         do is = 1, nsphere
            if (is .eq. 1 .and. sphere_1_fixed) cycle
            intarget = .false.
            do while (.not. intarget)
               tpos(1:3) = pos(1:3, is) + u(1:3, is) * tcmin
               do m = 1, 3
                  if (pbc(m)) then
                     if (tpos(m) .ge. wallboundaries(m, 2)) then
                        tpos(m) = tpos(m) - (wallboundaries(m, 2) - wallboundaries(m, 1))
                     elseif (tpos(m) .lt. wallboundaries(m, 1)) then
                        tpos(m) = tpos(m) + (wallboundaries(m, 2) - wallboundaries(m, 1))
                     end if
                  end if
               end do
               call check_position_in_target(radius(is), tpos, wallboundaries, intarget)
               if (.not. intarget) tcmin = .95 * tcmin
            end do
            if (intarget) then
               pos(:, is) = tpos
            else
               write (*, '('' out of target'')')
               write (*, '(8es12.4)') tpos, sqrt(sum(tpos**2)), tcmin, tcoll, twallmin
            end if
         end do
         sim_timings(5) = sim_timings(5) + mstm_mpi_wtime() - time_0
         time_0 = mstm_mpi_wtime()
         if (tcmin .lt. tmove) then
            nwhits = nwhits + 1
            if (wallcollision) then
               if (target_shape .eq. 0) then
                  u(iwall, iswall) = -u(iwall, iswall)
               elseif (target_shape .eq. 1) then
                  if (iwall .le. 2) then
                     rho = sqrt(pos(1, iswall) * pos(1, iswall) + pos(2, iswall) * pos(2, iswall))
                     cp = pos(1, iswall) / rho
                     sp = pos(2, iswall) / rho
                     urho = cp * u(1, iswall) + sp * u(2, iswall)
                     uphi = -sp * u(1, iswall) + cp * u(2, iswall)
                     u(1, iswall) = -cp * urho - sp * uphi
                     u(2, iswall) = -sp * urho + cp * uphi
                  else
                     u(iwall, iswall) = -u(iwall, iswall)
                  end if
               elseif (target_shape .eq. 2) then
                  rho = sqrt(pos(1, iswall) * pos(1, iswall) + pos(2, iswall) * pos(2, iswall))
                  if (rho .eq. 0.d0) then
                     cp = 1.d0
                     sp = 0.d0
                  else
                     cp = pos(1, iswall) / rho
                     sp = pos(2, iswall) / rho
                  end if
                  r = sqrt(rho * rho + pos(3, iswall) * pos(3, iswall))
                  if (r .eq. 0.d0) then
                     ct = 1.d0
                     st = 0.d0
                  else
                     ct = pos(3, iswall) / r
                     st = rho / r
                  end if
                  u1pn(1) = (u(1, iswall) * cp + u(2, iswall) * sp) * st + u(3, iswall) * ct
                  u1pn(2) = (u(1, iswall) * cp + u(2, iswall) * sp) * ct - u(3, iswall) * st
                  u1pn(3) = u(1, iswall) * sp - u(2, iswall) * cp
                  u1pn(1) = -u1pn(1)
                  u(1, iswall) = (u1pn(1) * st + u1pn(2) * ct) * cp + u1pn(3) * sp
                  u(2, iswall) = (u1pn(1) * st + u1pn(2) * ct) * sp - u1pn(3) * cp
                  u(3, iswall) = u1pn(1) * ct - u1pn(2) * st
               end if
            elseif (collision) then
               is = collisionpair(1)
               js = collisionpair(2)
               if (is .eq. 1 .and. sphere_1_fixed) then
                  call resolve_collision_velocities(1.d20, pos(1:3, is), u(1:3, is), 1.d0, &
                                                    pos(1:3, js), u(1:3, js), u1new, u2new)
               else
                  call resolve_collision_velocities(1.d0, collpos(1:3), u(1:3, is), 1.d0, &
                                                    pos(1:3, js), u(1:3, js), u1new, u2new)
               end if
               u(1:3, is) = u1new(1:3)
               u(1:3, js) = u2new(1:3)
            end if
         end if
         tmove = tmove - abs(tcmin)
         coll_data(1:nsphere)%time = coll_data(1:nsphere)%time - abs(tcmin)
         if (wallcollision) then
            call find_next_sphere_collision(nsphere, radius, pos, u, tmove, wallboundaries, &
                                            collision, tcoll, collisionpair, &
                                            start_sphere=iswall, end_sphere=iswall)
         elseif (collision) then
            is = collisionpair(1)
            js = collisionpair(2)
            call find_next_sphere_collision(nsphere, radius, pos, u, tmove, wallboundaries, &
                                            collision, tcoll, collisionpair, start_sphere=is, end_sphere=is)
            call find_next_sphere_collision(nsphere, radius, pos, u, tmove, wallboundaries, &
                                            collision, tcoll, collisionpair, start_sphere=js, end_sphere=js)
         end if
         i = i + 1
         if (sphere_1_fixed) u(:, 1) = 0.d0
         sim_timings(6) = sim_timings(6) + mstm_mpi_wtime() - time_0
      end do
      if (present(number_wall_hits)) number_wall_hits = nwhits
   end subroutine move_spheres

   subroutine find_next_wall_collision(nsphere, radius, pos, u, tmove, wallboundaries, twallmin, is, iswall, &
                                       start_sphere, end_sphere)
      implicit none
      integer :: nsphere, is, iwall, i, iswall, i1, i2
      integer, optional :: start_sphere, end_sphere
      real(8) :: pos(3, nsphere), radius(nsphere), tmove, wallboundaries(3, 2), &
                 twall, u(3, nsphere), twallmin, vel, rho, dist, cp, sp, ct, st, r, urho
      if (present(start_sphere)) then
         i1 = start_sphere
      else
         i1 = 1
      end if
      if (present(end_sphere)) then
         i2 = end_sphere
      else
         i2 = nsphere
      end if
      twallmin = tmove
      if (target_shape .eq. 0) then
         do iwall = 1, 3
            if (periodic_bc(iwall)) cycle
            do i = i1, i2
               vel = u(iwall, i)
               if (vel .lt. 0.d0) then
                  dist = -pos(iwall, i) + wallboundaries(iwall, 1) + radius(i) * wall_boundary_model + minimum_gap
                  twall = dist / vel
               elseif (vel .gt. 0.d0) then
                  dist = wallboundaries(iwall, 2) - pos(iwall, i) - radius(i) * wall_boundary_model - minimum_gap
                  twall = dist / vel
               else
                  twall = 1.d6
               end if
               if (twall .lt. twallmin) then
                  twallmin = twall
                  is = i
                  iswall = iwall
               end if
            end do
         end do
      elseif (target_shape .eq. 1) then
         do iwall = 2, 3
            if (iwall .eq. 3 .and. periodic_bc(iwall)) cycle
            do i = i1, i2
               if (iwall .lt. 3) then
                  rho = sqrt(pos(1, i) * pos(1, i) + pos(2, i) * pos(2, i))
                  if (rho .ne. 0.d0) then
                     vel = (pos(1, i) * u(1, i) + pos(2, i) * u(2, i)) / rho
                  else
                     vel = sqrt(u(1, i) * u(1, i) + u(2, i) * u(2, i))
                  end if
               else
                  vel = u(iwall, i)
               end if
               if (vel .lt. 0.d0) then
                  if (iwall .lt. 3) then
                     dist = -rho - wallboundaries(1, 2) + radius(i) * wall_boundary_model + minimum_gap
                  else
                     dist = -pos(iwall, i) + wallboundaries(iwall, 1) + radius(i) * wall_boundary_model + minimum_gap
                  end if
                  twall = dist / vel
               elseif (vel .gt. 0.d0) then
                  if (iwall .lt. 3) then
                     dist = wallboundaries(1, 2) - rho - radius(i) * wall_boundary_model - minimum_gap
                  else
                     dist = wallboundaries(iwall, 2) - pos(iwall, i) - radius(i) * wall_boundary_model - minimum_gap
                  end if
                  twall = dist / vel
               else
                  twall = 1.d6
               end if
               if (twall .lt. twallmin) then
                  twallmin = twall
                  is = i
                  iswall = iwall
               end if
            end do
         end do
      elseif (target_shape .eq. 2) then
         do i = i1, i2
            rho = sqrt(pos(1, i) * pos(1, i) + pos(2, i) * pos(2, i))
            if (rho .eq. 0.d0) then
               cp = 1.d0
               sp = 0.d0
            else
               cp = pos(1, i) / rho
               sp = pos(2, i) / rho
            end if
            r = sqrt(rho * rho + pos(3, i) * pos(3, i))
            if (r .eq. 0.d0) then
               vel = sqrt(dot_product(u(:, i), u(:, i)))
            else
               ct = pos(3, i) / r
               st = rho / r
               urho = cp * u(1, i) + sp * u(2, i)
               vel = urho * st + u(3, i) * ct
            end if
            if (vel .lt. 0.d0) then
               dist = -r - wallboundaries(1, 2) + radius(i) * wall_boundary_model + minimum_gap
               twall = dist / vel
            elseif (vel .gt. 0.d0) then
               dist = wallboundaries(1, 2) - r - radius(i) * wall_boundary_model - minimum_gap
               twall = dist / vel
            else
               twall = 1.d6
            end if
            if (twall .lt. twallmin) then
               twallmin = twall
               is = i
               iswall = 3
            end if
         end do
      end if
   end subroutine find_next_wall_collision

   subroutine find_next_sphere_collision(nsphere, radius, pos, u, maxtime, wallboundaries, collision, &
                                         tcmin, collisionpair, start_sphere, end_sphere, minimum_distance, collision_pos)
      implicit none
      logical :: collision, bndok, loccoll, pbc(3)
      integer :: nsphere, is, j, js, cell(3), ccell(3), scell(3), collisionpair(2), m, n, &
                 istart, iend, i
      integer, optional :: start_sphere, end_sphere
      real(8) :: pos(3, nsphere), radius(nsphere), maxtime, wallboundaries(3, 2), &
                 tcmin, rcol, tcollision, u(1:3, nsphere), mindist, tpos(3), collisionpos(3)
      real(8), optional :: minimum_distance, collision_pos(3)
      type(l_list), pointer :: llist
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
      if (present(minimum_distance)) then
         mindist = minimum_distance
      else
         mindist = minimum_gap
      end if
      if (target_shape .eq. 0) then
         pbc = periodic_bc
      elseif (target_shape .eq. 1) then
         pbc(3) = periodic_bc(3)
      end if
      tcmin = maxtime
      collision = .false.
      do is = istart, iend
         coll_data(is)%wallcoll = .false.
         coll_data(is)%time = maxtime
         coll_data(is)%collpos(:) = pos(:, is)
         call position_to_cell_index(pos(:, is), ccell)
         do m = 0, 26
            scell(1) = mod(m, 3) - 1
            scell(2) = mod(m / 3, 3) - 1
            scell(3) = mod(m / 9, 3) - 1
            cell = ccell + scell
            bndok = .true.
            tpos = pos(:, is)
            do i = 1, 3
               if (cell(i) .lt. 1 .or. cell(i) .gt. cell_dim(i)) then
                  if (pbc(i)) then
                     if (cell(i) .lt. 1) then
                        cell(i) = cell_dim(i)
                        tpos(i) = tpos(i) + wallboundaries(i, 2) - wallboundaries(i, 1)
                     elseif (cell(i) .gt. cell_dim(i)) then
                        cell(i) = 1
                        tpos(i) = tpos(i) - wallboundaries(i, 2) + wallboundaries(i, 1)
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
               js = llist%index
               if (js .ne. is) then
                  rcol = radius(is) + radius(js) + mindist
                  call test_pair_collision(tpos(1:3), u(1:3, is), pos(1:3, js), u(1:3, js), &
                                           rcol, loccoll, tcollision)
                  if (loccoll) then
                     if (coll_data(is)%time .gt. tcollision) then
                        coll_data(is)%time = tcollision
                        coll_data(is)%sphere = js
                        coll_data(is)%collpos(:) = tpos(:)
                     end if
                     if (tcollision .lt. tcmin) then
                        collision = .true.
                        tcmin = tcollision
                        collisionpair(1:2) = (/is, js/)
                        collisionpos = tpos
                     end if
                  end if
               end if
               if (j .lt. n) llist => llist%next
            end do
         end do
      end do
      if (present(collision_pos)) collision_pos = collisionpos
   end subroutine find_next_sphere_collision

   pure subroutine test_pair_collision(pos1, u1, pos2, u2, rcol, collision, tcollision)
      implicit none
      real(8), intent(in) :: pos1(3), u1(3), pos2(3), u2(3), rcol
      logical, intent(out) :: collision
      real(8), intent(out) :: tcollision
      real(8) :: urel(3), posrel(3), a, b, c, d
      urel = u2 - u1
      posrel = pos2 - pos1
      b = 2.d0 * dot_product(urel, posrel)
      if (b .ge. 0.d0) then
         collision = .false.
         return
      end if
      a = dot_product(urel, urel)
      c = max(dot_product(posrel, posrel) - rcol * rcol, 0.d0)
      if (c .eq. 0.d0) then
         collision = .true.
         tcollision = 0.d0
         return
      end if
      d = b * b - 4.d0 * a * c
      if (d .lt. 0.d0) then
         collision = .false.
         return
      end if
!         tc1=-(b+sqrt(d))/2.d0/a
!         tc2=-(b-sqrt(d))/2.d0/a
      tcollision = -(b + sqrt(b * b - 4.d0 * a * c)) / 2.d0 / a
      collision = .true.
   end subroutine test_pair_collision

   pure subroutine resolve_collision_velocities(mass1, pos1, u1, mass2, pos2, u2, u1new, u2new)
      implicit none
      real(8), intent(in) :: mass1, pos1(3), u1(3), mass2, pos2(3), u2(3)
      real(8), intent(out) :: u1new(3), u2new(3)
      real(8) :: posrel(3), rc, cosb, sinb, alpha, cosa, sina, rotmat(3, 3), u1p(3), u2p(3), &
                 u1pn(3), u2pn(3)
      posrel = pos2 - pos1
      rc = sqrt(dot_product(posrel, posrel))
      cosb = posrel(3) / rc
      sinb = sqrt((1.d0 - cosb) * (1.d0 + cosb))
      if (posrel(1) .eq. 0.d0 .and. posrel(2) .eq. 0.d0) then
         alpha = 0.d0
      else
         alpha = datan2(posrel(2), posrel(1))
      end if
      cosa = cos(alpha)
      sina = sin(alpha)
      rotmat = reshape((/cosa * cosb, -sina, cosa * sinb, sina * cosb, cosa, sina * sinb, -sinb, 0.d0, cosb/), (/3, 3/))
      u1p = matmul(rotmat, u1)
      u2p = matmul(rotmat, u2)
      u1pn(1:2) = u1p(1:2)
      u1pn(3) = ((mass1 - mass2) * u1p(3) + 2.d0 * mass2 * u2p(3)) / (mass1 + mass2)
      u2pn(1:2) = u2p(1:2)
      u2pn(3) = ((mass2 - mass1) * u2p(3) + 2.d0 * mass1 * u1p(3)) / (mass1 + mass2)
      u1new = matmul(transpose(rotmat), u1pn)
      u2new = matmul(transpose(rotmat), u2pn)
   end subroutine resolve_collision_velocities

   subroutine sample_random_velocities(nsphere, u)
      implicit none
      integer :: i, nsphere
      real(8) :: u(3, nsphere), cb, sb, alpha, ca, sa, rannum(2)
      do i = 1, nsphere
         call random_number(rannum)
         cb = -1.d0 + 2.d0 * rannum(1)
         sb = sqrt((1.d0 - cb) * (1.d0 + cb))
         alpha = two_pi * rannum(2)
         ca = cos(alpha)
         sa = sin(alpha)
         u(1:3, i) = (/ca * sb, sa * sb, cb/)
      end do
   end subroutine sample_random_velocities

   subroutine sample_particle_radius(sigma, maxradius, x)
      implicit none
      integer :: i
      real(8) :: sigma, maxradius, f1, fd, x, fmax, s2, xmax, &
                 t1, rannum(2)
      if (sigma .eq. 0.d0) then
         x = 1.d0
         return
      end if
      s2 = sigma * sigma
      f1 = 1.d0
      fd = 0.d0
      xmax = exp(-2.5d0 * s2)
      t1 = (log(xmax) + 1.5d0 * s2)
      fmax = exp(-t1 * t1 / (2.d0 * s2)) / sqrt_two_pi / xmax / sigma
      i = 0
      do while (f1 .gt. fd)
         i = i + 1
         call random_number(rannum)
         x = maxradius * rannum(1)
         f1 = fmax * rannum(2)
         t1 = (log(x) + 1.5d0 * s2)
         fd = exp(-t1 * t1 / (2.d0 * s2)) / sqrt_two_pi / x / sigma
      end do
   end subroutine sample_particle_radius
end module random_configuration_dynamics
