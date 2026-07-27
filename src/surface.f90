module surface
   use, intrinsic :: iso_fortran_env, only: real64
   use constants
   use numerical_tables
   use angular_functions, only: axial_translation_offset, axial_translation_size, &
                                complex_vector_spherical_harmonics, generate_plane_wave_coefficients
   use bessel_functions, only: bessel_integer_complex
   use coefficient_indexing, only: polarized_mode_index
   use quadrature, only: integrate_gauss_kronrod_adaptive, sort_unique_real_values
   use parallel_runtime, only: mstm_global_rank
   implicit none
   logical :: source_sum, include_direct_source, &
              pole_integration, plane_surface_present
   integer, parameter :: max_number_plane_boundaries = 10, max_singular_points = 100
   integer :: source_order, target_order, max_azimuth_mode, number_layers, source_layer, target_layer, &
              source_layer2, number_limits, max_gf_iterations, number_gf_iterations, source_order2, &
              error_codes(4), energy_kernel_region, number_singular_points, singular_point_polarization(max_singular_points)
   integer, target :: number_plane_boundaries, maximum_integration_subdivisions
   real(real64) :: source_z, target_z, azimuth_angle, incident_field_boundary, &
                   radial_distance, max_s, gf_error_epsilon, &
                   s_sc1, s_sc2, max_gf, max_bf, max_pi, max_picon, &
                   top_boundary, bot_boundary, source_z2, integration_error, singular_points(max_singular_points), &
                   incident_field_scale(2), pole_integration_radius, singular_gf_value(max_singular_points), &
                   incident_lateral_vector(2), g_cut, g_sing_mag
   real(real64), target :: layer_thickness(max_number_plane_boundaries), integration_limit_epsilon, &
                           integration_error_epsilon, real_axis_integration_limit, gf_switch_factor, s_scale_constant, &
                           minimum_integration_spacing, minimum_initial_segment_size
   real(real64), allocatable :: plane_boundary_position(:), real_axis_limits(:)
   complex(real64) :: source_ri, target_ri, pole_integration_s
   complex(real64), target :: layer_ref_index(0:max_number_plane_boundaries)
   complex(real64), allocatable :: source_coefficient(:, :), source_coefficient_1(:, :, :), source_coefficient_2(:, :, :)
   data real_axis_integration_limit, minimum_integration_spacing, minimum_initial_segment_size/100.d0, 1.d-5, 0.d0/
   data include_direct_source/.false./
   data max_gf_iterations, gf_switch_factor, number_gf_iterations, gf_error_epsilon/1000, .5d0, 0, 1.d-7/
   data s_scale_constant/0.01d0/
   data layer_thickness/max_number_plane_boundaries*1.d0/
   data layer_ref_index(0)/(1.d0, 0.d0)/
   data number_plane_boundaries/0/
   data integration_limit_epsilon, integration_error_epsilon, maximum_integration_subdivisions/1.d-15, 1.d-6, 25/
   data layer_ref_index(1:max_number_plane_boundaries)/max_number_plane_boundaries*(1.d0, 0.d0)/
   data pole_integration_radius/1.d-6/
   data g_cut, g_sing_mag/10.d0, 1.d5/

contains
   subroutine initialize_plane_boundaries()
      implicit none
      integer :: i
      real(real64) :: smax
      if (allocated(plane_boundary_position)) deallocate (plane_boundary_position)
      allocate (plane_boundary_position(1:max(number_plane_boundaries, 1)))
      plane_boundary_position(1) = 0.d0
      do i = 1, number_plane_boundaries - 1
         plane_boundary_position(i + 1) = plane_boundary_position(i) + layer_thickness(i)
      end do
      smax = maxval(dble(layer_ref_index(0:number_plane_boundaries)))
      top_boundary = plane_boundary_position(max(1, number_plane_boundaries)) + 1.d-8
      bot_boundary = -1.d-8
      if (number_plane_boundaries .gt. 1) then
         call find_green_function_singular_points(bot_boundary, top_boundary, g_cut, smax, &
                                                  number_singular_points, singular_points, singular_point_polarization, &
                                                  singular_gf_value)
      else
         number_singular_points = 0
         singular_gf_value = 1.d0
      end if
   end subroutine initialize_plane_boundaries

   subroutine find_green_function_singular_points(sourcez, targetz, gcut, smax, nsingpoints, singpoints, singpol, singval)
      implicit none
      integer :: p, i, nbrack, nsingpoints, singpol(*)
    real(real64) :: sourcez, targetz, gcut, smax, singpoints(*), sbrack(2, max_singular_points), s0, s1, s2, fmax, sfmax, singval(*)
      nsingpoints = 0
      do p = 1, 2
         call bracket_green_function_singular_points(sourcez, targetz, p, gcut, smax, nbrack, sbrack)
         do i = 1, nbrack
            s0 = sbrack(1, i)
            s2 = sbrack(2, i)
            s1 = 0.5d0 * (s0 + s2)
            call maximize_green_function(sourcez, targetz, p, s0, s1, s2, 1.d-9, 100, fmax, sfmax)
            if (fmax .lt. g_sing_mag) cycle
            nsingpoints = nsingpoints + 1
            singpoints(nsingpoints) = sfmax
            singpol(nsingpoints) = p
            singval(nsingpoints) = fmax
         end do
      end do
   end subroutine find_green_function_singular_points

   subroutine bracket_green_function_singular_points(sz, tz, p, gcut, smax, nbrack, sbrack)
      implicit none
      logical :: inbrack
      integer :: nbrack, p
      real(real64) :: sz, tz, smax, sbrack(2, *), dels, fm, gcut
      complex(real64) :: gf(2, 2, 2), skz, tkz, s
      data dels/1.d-3/
      inbrack = .false.
      nbrack = 0
      s = dels * 0.5d0
      do while (dble(s) .lt. smax)
         call layer_green_function(s, sz, tz, gf, skz, tkz)
         fm = dble(sum(abs(gf(:, :, p))))
         if (fm .gt. gcut) then
            if (.not. inbrack) then
               nbrack = nbrack + 1
               sbrack(1, nbrack) = dble(s)
               inbrack = .true.
            end if
         else
            if (inbrack) then
               inbrack = .false.
               sbrack(2, nbrack) = dble(s)
               if (nbrack .eq. max_singular_points) then
                  write (*, '('' max number GF singular points exceeded'')')
                  exit
               end if
            end if
         end if
         s = s + dels
      end do
   end subroutine bracket_green_function_singular_points

   subroutine maximize_green_function(sz, tz, p, ax, bx, cx, tol, maxsteps, gmax, xmax)
      implicit none
      integer :: n, maxsteps, p
      real(real64) :: ax, bx, cx, tol, xmax, r, c, x0, x3, x1, x2, f1, f2, f3, f0, gmax, sz, tz
      complex(real64) :: gf(2, 2, 2), skz, tkz, s
      data r, c/.61803399d0, .38196602d0/
      x0 = ax
      x3 = cx
      if (abs(cx - bx) .gt. abs(bx - ax)) then
         x1 = bx
         x2 = bx + c * (cx - bx)
      else
         x2 = bx
         x1 = bx - c * (bx - ax)
      end if
      s = x1
      call layer_green_function(s, sz, tz, gf, skz, tkz)
      f1 = dble(sum(abs(gf(:, :, p))))
      s = x2
      call layer_green_function(s, sz, tz, gf, skz, tkz)
      f2 = dble(sum(abs(gf(:, :, p))))
      n = 1
      do while (abs(x3 - x0) .gt. tol * (abs(x1) + abs(x2)) .and. n .le. maxsteps)
         n = n + 1
         if (f2 .gt. f1) then
            x0 = x1
            x1 = x2
            x2 = r * x1 + c * x3
            f0 = f1
            f1 = f2
            s = x2
            call layer_green_function(s, sz, tz, gf, skz, tkz)
            f2 = dble(sum(abs(gf(:, :, p))))
         else
            x3 = x2
            x2 = x1
            x1 = r * x2 + c * x0
            f3 = f2
            f2 = f1
            s = x1
            call layer_green_function(s, sz, tz, gf, skz, tkz)
            f1 = dble(sum(abs(gf(:, :, p))))
         end if
      end do
      if (f1 .gt. f2) then
         gmax = f1
         xmax = x1
      else
         gmax = f2
         xmax = x2
      end if
   end subroutine maximize_green_function

   subroutine layer_green_function_series(sourcelayer, targetlayer, sourcez, targetz, kz, omega, tm, gfs)
      implicit none
      integer :: n, sourcelayer, targetlayer, p, sourcedir, iter
      real(real64) :: scale, scale0, sourcez, targetz
      complex(real64) :: gf(2, number_plane_boundaries, 2), kz(0:number_plane_boundaries), omega(0:number_plane_boundaries), &
tm(2, 2, 0:number_plane_boundaries + 1, 2), gfs(2, 2, 2), gp(2, number_plane_boundaries, 2), gftot(2, number_plane_boundaries, 2), &
                         tfup, tfdn
      if (targetlayer .gt. 0) then
         tfup = exp((0.d0, 1.d0) * layer_ref_index(targetlayer) * kz(targetlayer) &
                    * (targetz - plane_boundary_position(targetlayer)))
      else
         tfup = 0.d0
      end if
      if (targetlayer .lt. number_plane_boundaries) then
         tfdn = exp((0.d0, 1.d0) * layer_ref_index(targetlayer) * kz(targetlayer) &
                    * (plane_boundary_position(targetlayer + 1) - targetz))
      else
         tfdn = 0.d0
      end if
      do sourcedir = 1, 2
         gfs(:, sourcedir, :) = 0.d0
         gf = 0.d0
         if (sourcedir .eq. 1) then
            if (sourcelayer .lt. number_plane_boundaries) then
               gf(:, sourcelayer + 1, :) = tm(:, 1, sourcelayer + 1, :) &
                                           * exp((0.d0, 1.d0) * layer_ref_index(sourcelayer) * kz(sourcelayer) &
                                                 * (plane_boundary_position(sourcelayer + 1) - sourcez))
            else
               cycle
            end if
         end if
         if (sourcedir .eq. 2) then
            if (sourcelayer .gt. 0) then
               gf(:, sourcelayer, :) = tm(:, 2, sourcelayer, :) &
                                       * exp((0.d0, 1.d0) * layer_ref_index(sourcelayer) * kz(sourcelayer) &
                                             * (sourcez - plane_boundary_position(sourcelayer)))
            else
               cycle
            end if
         end if
         gftot = gf
         if (targetlayer .gt. 0) then
            gfs(1, sourcedir, :) = gftot(1, targetlayer, :) * tfup
         end if
         if (targetlayer .lt. number_plane_boundaries) then
            gfs(2, sourcedir, :) = gftot(2, targetlayer + 1, :) * tfdn
         end if
         scale0 = sqrt(sum(abs(gftot(:, :, :))**2))
         do iter = 1, max_gf_iterations
            gp = 0
            do p = 1, 2
               do n = 2, number_plane_boundaries
                  gp(:, n, p) = omega(n - 1) * tm(:, 1, n, p) * gf(1, n - 1, p)
               end do
               do n = 1, number_plane_boundaries - 1
                  gp(:, n, p) = gp(:, n, p) + omega(n) * tm(:, 2, n, p) * gf(2, n + 1, p)
               end do
            end do
            gftot(:, :, :) = gftot(:, :, :) + gp(:, :, :)
            gf(:, :, :) = gp(:, :, :)
            if (targetlayer .gt. 0) then
               gfs(1, sourcedir, :) = gftot(1, targetlayer, :) * tfup
            end if
            if (targetlayer .lt. number_plane_boundaries) then
               gfs(2, sourcedir, :) = gftot(2, targetlayer + 1, :) * tfdn
            end if
            scale = sqrt(sum(abs(gftot(:, :, :))**2))
            if (abs(scale - scale0) / max(scale0, 1.d-12) .lt. gf_error_epsilon) then
               number_gf_iterations = max(number_gf_iterations, iter)
               exit
            end if
            scale0 = scale
         end do
         number_gf_iterations = max(number_gf_iterations, iter)
      end do
      if (number_gf_iterations .gt. max_gf_iterations) then
         error_codes(1) = 1
      end if
   end subroutine layer_green_function_series

   subroutine layer_green_function_recurrence(sourcelayer, targetlayer, sourcez, targetz, kz, omega, tm, gfs)
      implicit none
      integer :: n, sourcelayer, targetlayer, p, sourcedir
      real(real64) :: sourcez, targetz
      complex(real64) :: gf(2, 0:number_plane_boundaries + 1), kz(0:number_plane_boundaries), omega(0:number_plane_boundaries), &
                         tm(2, 2, 0:number_plane_boundaries + 1, 2), gfs(2, 2, 2), tfup, tfdn
!         complex(16) :: um(2,2,0:number_plane_boundaries),sm(2,2,0:number_plane_boundaries), &
!            sv(2,0:number_plane_boundaries+1),tempv(2),umt(2,2),bv(2,0:number_plane_boundaries+1),umt0(2,2)
      complex(real64) :: um(2, 2, 0:number_plane_boundaries), sm(2, 2, 0:number_plane_boundaries), &
                         sv(2, 0:number_plane_boundaries + 1), tempv(2), umt(2, 2), bv(2, 0:number_plane_boundaries + 1), umt0(2, 2)
      if (targetlayer .gt. 0) then
         tfup = exp((0.d0, 1.d0) * layer_ref_index(targetlayer) * kz(targetlayer) &
                    * (targetz - plane_boundary_position(targetlayer)))
      else
         tfup = 0.d0
      end if
      if (targetlayer .lt. number_plane_boundaries) then
         tfdn = exp((0.d0, 1.d0) * layer_ref_index(targetlayer) * kz(targetlayer) &
                    * (plane_boundary_position(targetlayer + 1) - targetz))
      else
         tfdn = 0.d0
      end if
      do p = 1, 2
         do n = 0, number_plane_boundaries
            um(1, 1, n) = omega(n) * tm(1, 1, n, p)
            um(1, 2, n) = omega(n) * tm(1, 2, n, p)
            um(2, 1, n) = -omega(n) * tm(2, 1, n + 1, p) * tm(1, 1, n, p) / tm(2, 2, n + 1, p)
            um(2, 2, n) = (1.d0 - omega(n) * omega(n) * tm(2, 1, n + 1, p) * tm(1, 2, n, p)) / (omega(n) * tm(2, 2, n + 1, p))
            sm(1, 1, n) = 1.d0
            sm(1, 2, n) = 0.d0
            sm(2, 1, n) = -tm(2, 1, n + 1, p) / tm(2, 2, n + 1, p)
            sm(2, 2, n) = -1.d0 / (omega(n) * tm(2, 2, n + 1, p))
         end do
         do sourcedir = 1, 2
            sv = 0.d0
            gf(:, :) = 0.d0
            gfs(:, sourcedir, p) = 0.d0
            if (sourcedir .eq. 1) then
               if (sourcelayer .eq. number_plane_boundaries) then
                  cycle
               else
                  sv(1, sourcelayer + 1) = exp((0.d0, 1.d0) * layer_ref_index(sourcelayer) * kz(sourcelayer) &
                                               * (plane_boundary_position(sourcelayer + 1) - sourcez))
               end if
            elseif (sourcedir .eq. 2) then
               if (sourcelayer .eq. 0) then
                  cycle
               else
                  sv(2, sourcelayer) = exp((0.d0, 1.d0) * layer_ref_index(sourcelayer) * kz(sourcelayer) &
                                           * (sourcez - plane_boundary_position(sourcelayer)))
               end if
            end if
            tempv(1) = sm(1, 1, sourcelayer) * sv(1, sourcelayer + 1)
            tempv(2) = sm(2, 1, sourcelayer) * sv(1, sourcelayer + 1) + sm(2, 2, sourcelayer) * sv(2, sourcelayer)
            bv(:, sourcelayer + 1) = tempv(:)
            do n = sourcelayer + 1, number_plane_boundaries
               bv(1, n + 1) = um(1, 1, n) * bv(1, n) + um(1, 2, n) * bv(2, n)
               bv(2, n + 1) = um(2, 1, n) * bv(1, n) + um(2, 2, n) * bv(2, n)
            end do
            umt(:, :) = um(:, :, 0)
            do n = 1, number_plane_boundaries
               umt0 = umt
               umt(1, 1) = um(1, 1, n) * umt0(1, 1) + um(1, 2, n) * umt0(2, 1)
               umt(2, 1) = um(2, 1, n) * umt0(1, 1) + um(2, 2, n) * umt0(2, 1)
               umt(1, 2) = um(1, 1, n) * umt0(1, 2) + um(1, 2, n) * umt0(2, 2)
               umt(2, 2) = um(2, 1, n) * umt0(1, 2) + um(2, 2, n) * umt0(2, 2)
            end do
            gf(2, 0) = -bv(2, number_plane_boundaries + 1) / umt(2, 2)
            do n = 0, number_plane_boundaries
               gf(1, n + 1) = um(1, 1, n) * gf(1, n) + um(1, 2, n) * gf(2, n) + sm(1, 1, n) * sv(1, n + 1)
               gf(2, n + 1) = um(2, 1, n) * gf(1, n) + um(2, 2, n) * gf(2, n) + sm(2, 1, n) * sv(1, n + 1) + sm(2, 2, n) * sv(2, n)
            end do
            if (targetlayer .gt. 0) then
               gfs(1, sourcedir, p) = tfup * sum(tm(1, :, targetlayer, p) * gf(:, targetlayer))
            end if
            if (targetlayer .lt. number_plane_boundaries) then
               gfs(2, sourcedir, p) = tfdn * sum(tm(2, :, targetlayer + 1, p) * gf(:, targetlayer + 1))
            end if
         end do
      end do
   end subroutine layer_green_function_recurrence

   subroutine layer_green_function(s, sourcez, targetz, gfs, sourcekz, targetkz, include_direct)
      implicit none
      logical :: incdir, prop
      logical, optional :: include_direct
      integer :: n, sourcelayer, targetlayer, i
      real(real64) :: sourcez, targetz, c
      complex(real64) :: kz(0:number_plane_boundaries), omega(0:number_plane_boundaries), &
                         den, tm(2, 2, 0:number_plane_boundaries + 1, 2), &
                         gfs(2, 2, 2), &
                         sourcekz, targetkz, s
      if (present(include_direct)) then
         incdir = include_direct
      else
         incdir = .false.
      end if
      c = 1.d0
      sourcelayer = find_layer_index(sourcez)
      targetlayer = find_layer_index(targetz)
      do n = 0, number_plane_boundaries
         kz(n) = sqrt((layer_ref_index(n) - s) * (layer_ref_index(n) + s)) / layer_ref_index(n)
      end do
      if (incdir) then
         if (sourcez .eq. targetz) then
            c = 0.5d0
         else
            c = 1.d0
         end if
      end if
      prop = dble(s) .le. dble(layer_ref_index(targetlayer))
      if (number_plane_boundaries .eq. 0) then
         sourcekz = kz(0)
         targetkz = kz(0)
         gfs = 0.d0
         if (incdir) then
            if (targetz .le. sourcez) then
               gfs(2, 2, :) = c * exp((0.d0, 1.d0) * layer_ref_index(targetlayer) * sourcekz * (abs(sourcez - targetz)))
            end if
            if (targetz .ge. sourcez) then
               gfs(1, 1, :) = c * exp((0.d0, 1.d0) * layer_ref_index(targetlayer) * sourcekz * (abs(sourcez - targetz)))
            end if
         end if
         return
      end if
      omega(0) = 1.d0
      omega(number_plane_boundaries) = 1.d0
      do n = 1, number_plane_boundaries - 1
         omega(n) = exp((0.d0, 1.d0) * layer_ref_index(n) * kz(n) * layer_thickness(n))
      end do
      do n = 1, number_plane_boundaries
         den = kz(n) * layer_ref_index(n - 1) + kz(n - 1) * layer_ref_index(n)
         tm(1, 1, n, 1) = 2.d0 * kz(n - 1) * layer_ref_index(n - 1) / den
         tm(1, 2, n, 1) = (kz(n) * layer_ref_index(n - 1) - kz(n - 1) * layer_ref_index(n)) / den
         tm(2, 1, n, 1) = -tm(1, 2, n, 1)
         tm(2, 2, n, 1) = 2.d0 * kz(n) * layer_ref_index(n) / den
         den = kz(n - 1) * layer_ref_index(n - 1) + kz(n) * layer_ref_index(n)
         tm(1, 1, n, 2) = 2.d0 * kz(n - 1) * layer_ref_index(n - 1) / den
         tm(1, 2, n, 2) = -(kz(n - 1) * layer_ref_index(n - 1) - kz(n) * layer_ref_index(n)) / den
         tm(2, 1, n, 2) = -tm(1, 2, n, 2)
         tm(2, 2, n, 2) = 2.d0 * kz(n) * layer_ref_index(n) / den
      end do
      sourcekz = kz(sourcelayer)
      targetkz = kz(targetlayer)
      gfs = 0.d0
      if (number_plane_boundaries .eq. 1) then
         gfs = 0.d0
         if (sourcelayer .eq. 0 .and. targetlayer .eq. 0) then
            gfs(2, 1, :) = exp((0.d0, 1.d0) * layer_ref_index(0) * sourcekz * (abs(sourcez) + abs(targetz))) &
                           * tm(2, 1, 1, :)
         elseif (sourcelayer .eq. 0 .and. targetlayer .eq. 1) then
            gfs(1, 1, :) = exp((0.d0, 1.d0) * (layer_ref_index(0) * sourcekz * abs(sourcez) &
                                               + layer_ref_index(1) * targetkz * abs(targetz))) * tm(1, 1, 1, :)
         elseif (sourcelayer .eq. 1 .and. targetlayer .eq. 1) then
            gfs(1, 2, :) = exp((0.d0, 1.d0) * layer_ref_index(1) * sourcekz * (abs(sourcez) + abs(targetz))) &
                           * tm(1, 2, 1, :)
         elseif (sourcelayer .eq. 1 .and. targetlayer .eq. 0) then
            gfs(2, 2, :) = exp((0.d0, 1.d0) * (layer_ref_index(1) * sourcekz * abs(sourcez) &
                                               + layer_ref_index(0) * targetkz * abs(targetz))) * tm(2, 2, 1, :)
         end if
      else
         do i = 1, 2
            tm(i, i, 0, :) = 1.d0
            tm(i, 3 - i, 0, :) = 0.d0
            tm(i, i, number_plane_boundaries + 1, :) = 1.d0
            tm(i, 3 - i, number_plane_boundaries + 1, :) = 0.d0
         end do
         if (maxval(abs(omega(1:number_plane_boundaries - 1))) .lt. gf_switch_factor) then
            number_gf_iterations = 0
            call layer_green_function_series(sourcelayer, targetlayer, sourcez, targetz, kz, omega, tm, gfs)
         else
            call layer_green_function_recurrence(sourcelayer, targetlayer, sourcez, targetz, kz, omega, tm, gfs)
         end if
      end if
      if (incdir .and. sourcelayer .eq. targetlayer) then
         if (targetz .le. sourcez) then
            gfs(2, 2, :) = gfs(2, 2, :) + c * exp((0.d0, 1.d0) * layer_ref_index(targetlayer) * sourcekz * (abs(sourcez - targetz)))
         end if
         if (targetz .ge. sourcez) then
            gfs(1, 1, :) = gfs(1, 1, :) &
                           + c * exp((0.d0, 1.d0) * layer_ref_index(targetlayer) * sourcekz * (abs(sourcez - targetz)))
         end if
      end if
   end subroutine layer_green_function

   subroutine real_axis_kernel(ntot, t, kernmat)
      implicit none
      integer :: m, m1, k, k1, mk, n, l, p, q, mn, kl, i, pol, dir, tsign, ssign, mmax, ntot
      integer, save :: count
      real(real64) :: t
      complex(real64) :: kernmat(ntot), temp, bfunc(0:source_order + target_order), &
                         s, const, const2, sr, dsdt, &
                         pivec(2, source_order * (source_order + 2), 2), picvec(2, target_order * (target_order + 2), 2), &
                         bvec(2, 2), sourcekz, targetkz, gfunc(2, 2, 2), sourceri, targetri, sources, targets
      data count/0/
      if (pole_integration) then
         dsdt = (0.d0, 1.d0) * pole_integration_radius * exp((0.d0, 1.d0) * t)
         s = pole_integration_s + pole_integration_radius * exp((0.d0, 1.d0) * t)
      else
         s = t * cmplx(1.d0, -s_sc1, kind=real64) &
             + cmplx(0.d0, -1.d0, kind=real64) * s_sc2
         dsdt = cmplx(1.d0, -s_sc1, kind=real64)
      end if
      call layer_green_function(s, source_z, target_z, gfunc, sourcekz, targetkz, include_direct_source)
      sourceri = layer_ref_index(source_layer)
      targetri = layer_ref_index(target_layer)
      sr = s * radial_distance
      sources = s / sourceri
      targets = s / targetri
      call complex_vector_spherical_harmonics(sourcekz, source_order, pivec, 1)
      call complex_vector_spherical_harmonics(targetkz, target_order, picvec, -1)
      bfunc = 0.d0
      mmax = target_order + source_order
      if (radial_distance .eq. 0.d0) then
         bfunc(0) = 1.d0
      else
         call bessel_integer_complex(target_order + source_order, sr, mmax, bfunc)
      end if
!         max_gf=max(max_gf,maxval(abs(gfunc)))
!         max_pi=max(max_pi,maxval(abs(pivec)))
!         max_picon=max(max_picon,maxval(abs(picvec)))
!         max_bf=max(max_bf,maxval(abs(bfunc)))
      const = four_pi * s / sourcekz / sourceri / sourceri * dsdt
      i = 0
      if (source_sum) then
         do m = -target_order, target_order
            bvec = 0.d0
            m1 = max(1, abs(m))
            do k = -source_order, source_order
               k1 = max(1, abs(k))
               mk = abs(k - m)
               if ((radial_distance .eq. 0.d0) .and. mk .ne. 0) cycle
               if (mk .gt. mmax) cycle
               const2 = exp((0.d0, 1.d0) * (k - m) * azimuth_angle) * ((0.d0, 1.d0)**mk) * bfunc(mk) * const
               do l = k1, source_order
                  kl = l * (l + 1) + k
                  do q = 1, 2
                     ssign = (-1)**(k + l + q - 1)
                     do pol = 1, 2
                        do dir = 1, 2
                           bvec(dir, pol) = bvec(dir, pol) + const2 * pivec(q, kl, pol) * (gfunc(dir, 1, pol) &
                                                               + ssign * (-1)**pol * gfunc(dir, 2, pol)) * source_coefficient(q, kl)
                        end do
                     end do
                  end do
               end do
            end do
            do n = m1, target_order
               mn = n * (n + 1) + m
               do p = 1, 2
                  tsign = (-1)**(m + n + p - 1)
                  i = i + 1
                  kernmat(i) = picvec(p, mn, 1) * (bvec(1, 1) - tsign * bvec(2, 1)) &
                               + picvec(p, mn, 2) * (bvec(1, 2) + tsign * bvec(2, 2))
               end do
            end do
         end do
      else
         do m = 0, target_order
            m1 = max(1, abs(m))
            do k = -source_order, source_order
               k1 = max(1, abs(k))
               mk = abs(k - m)
               if ((radial_distance .eq. 0.d0) .and. mk .ne. 0) then
                  cycle
               end if
               const2 = ((0.d0, 1.d0)**mk) * bfunc(mk) * const
               do n = m1, target_order
                  mn = m + n * (n + 1)
                  do l = k1, source_order
                     kl = k + l * (l + 1)
                     do p = 1, 2
                        tsign = (-1)**(m + n + p - 1)
                        do q = 1, 2
                           ssign = (-1)**(k + l + q - 1)
                           i = i + 1
                           temp = 0.d0
                           do pol = 1, 2
                              temp = temp + const2 * picvec(p, mn, pol) * pivec(q, kl, pol) &
                                     * (gfunc(1, 1, pol) + (-1)**pol * (tsign * gfunc(2, 1, pol) + ssign * gfunc(1, 2, pol)) &
                                        + tsign * ssign * gfunc(2, 2, pol))
                           end do
                           kernmat(i) = temp
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end if
   end subroutine real_axis_kernel

   subroutine energy_kernel(ntot, t, kernmat)
      implicit none
      integer :: m, m1, k, k1, mk, n, l, p, q, mn, kl, pol, tsign, ssign, mmax, sourceorder(2), &
                 nmax, spol, tlay, ntot
      integer, save :: count
      real(real64) :: t, sourcez(2), targetz, const3(3, 2)
      complex(real64) :: kernmat(ntot), bfunc(0:source_order + source_order2), &
                         s, const, const2, sr, cscale, &
                         pivec1(2, source_order * (source_order + 2), 2), &
                         pivec2(2, source_order2 * (source_order2 + 2), 2), bvec(2, 2), &
                         bvec1(2), bvec2(2), sourcekz(2), targetkz, gfunc1(2, 2, 2), gfunc2(2, 2, 2), sourceri(2), targetri
      data count/0/
      sourcez(1) = source_z
      sourcez(2) = source_z2
      targetz = target_z
      sourceri(1) = layer_ref_index(source_layer)
      sourceri(2) = layer_ref_index(source_layer2)
      sourceorder(1) = source_order
      sourceorder(2) = source_order2
      tlay = target_layer
      targetri = layer_ref_index(tlay)
      nmax = max(source_order, source_order2)
      if (energy_kernel_region .eq. 0) then
         s = sqrt(((1.d0, 0.d0) - t) * ((1.d0, 0.d0) + t)) * dble(targetri)
      else
         s = sqrt((1.d0, 0.d0) + t * t) * dble(targetri)
      end if
      call layer_green_function(s, sourcez(1), targetz, gfunc1, sourcekz(1), targetkz, include_direct=.true.)
      call layer_green_function(s, sourcez(2), targetz, gfunc2, sourcekz(2), targetkz, include_direct=.true.)
      cscale = t * dble(targetri)**2
      call complex_vector_spherical_harmonics(sourcekz(1), sourceorder(1), pivec1, 1)
      call complex_vector_spherical_harmonics(sourcekz(2), sourceorder(2), pivec2, 1)
      bfunc = 0.d0
      sr = s * radial_distance
      mmax = sum(sourceorder(:))
      if (radial_distance .eq. 0.d0) then
         bfunc(0) = 1.d0
      else
         call bessel_integer_complex(sum(sourceorder(:)), sr, mmax, bfunc)
      end if
      kernmat = 0.d0
      const = four_pi * cscale / sourcekz(1) / sourceri(1)**2 / conjg(sourcekz(2) * sourceri(2)**2)
      const3(1, 1) = (abs(targetkz)**2 + s * s / abs(targetri)**2) * dble(targetri * targetkz)
      const3(2, 1) = -const3(1, 1)
      const3(3, 1) = (abs(targetkz)**2 - s * s / abs(targetri)**2) * aimag(targetri * targetkz)
      const3(1, 2) = dble(targetri * targetkz)
      const3(2, 2) = -const3(1, 2)
      const3(3, 2) = -aimag(targetri * targetkz)
      do spol = 1, 2
         do pol = 1, 2
            bvec = 0.d0
            do m = -sourceorder(2), sourceorder(2)
               bvec1 = 0.d0
               bvec2 = 0.d0
               m1 = max(1, abs(m))
               do k = -sourceorder(1), sourceorder(1)
                  k1 = max(1, abs(k))
                  mk = abs(k - m)
                  if ((radial_distance .eq. 0.d0) .and. mk .ne. 0) cycle
                  if (mk .gt. mmax) cycle
                  const2 = exp((0.d0, 1.d0) * (k - m) * azimuth_angle) * ((0.d0, 1.d0)**mk) * bfunc(mk) * const
                  do l = k1, sourceorder(1)
                     kl = l * (l + 1) + k
                     do q = 1, 2
                        ssign = (-1)**(k + l + q - 1 + pol)
                        bvec1(:) = bvec1(:) + const2 * pivec1(q, kl, pol) * (gfunc1(:, 1, pol) &
                                                                    + ssign * gfunc1(:, 2, pol)) * source_coefficient_1(q, kl, spol)
                     end do
                  end do
               end do
               do n = m1, sourceorder(2)
                  mn = n * (n + 1) + m
                  do p = 1, 2
                     tsign = (-1)**(m + n + p - 1 + pol)
                     bvec2(:) = bvec2(:) + conjg(pivec2(p, mn, pol) * (gfunc2(:, 1, pol) &
                                                                       + tsign * gfunc2(:, 2, pol)) &
                                                 * source_coefficient_2(p, mn, spol))
                  end do
               end do
               do p = 1, 2
                  do q = 1, 2
                     bvec(p, q) = bvec(p, q) + bvec1(p) * bvec2(q)
                  end do
               end do
            end do
            kernmat(spol) = kernmat(spol) + const3(1, pol) * dble(bvec(1, 1)) + const3(2, pol) * dble(bvec(2, 2)) &
                            + const3(3, pol) * aimag(bvec(1, 2) - bvec(2, 1))
         end do
      end do
   end subroutine energy_kernel

   pure real(real64) function mode_vector_norm(n, m)
      implicit none
      integer, intent(in) :: n
      complex(real64), intent(in) :: m(n)
      mode_vector_norm = sqrt(dble(dot_product(m, m)))
   end function mode_vector_norm

   pure integer function find_layer_index(z)
      implicit none
      integer :: n
      real(real64), intent(in) :: z
      find_layer_index = 0
      do n = 1, number_plane_boundaries
         if (z .ge. plane_boundary_position(n)) then
            find_layer_index = n
         else
            return
         end if
      end do
   end function find_layer_index

   subroutine plane_boundary_interaction(ntot, ltot, x, y, sourcez, targetz, &
                                         interactionmatrix, source_vector, &
                                         index_model, lr_transformation, make_symmetric, propagating_directions_only)
      implicit none
      logical :: lrtran, makesymmetric, propdir
      logical, optional :: lr_transformation, make_symmetric, propagating_directions_only
      integer :: ntot, ltot, l, k, n, m, p, q, mn, kl, tranmat(2, 2), indexmodel, &
                 mnp, qtot, m1, k1, i, rmatdim, rmataddress
      integer, optional :: index_model
      real(real64) :: x, y, rho, sourcez, targetz, ssc, rail
      complex(real64) :: ephi, ephim, interactionmatrix(*), &
                         ctemp(2), c2temp(2, 2), c2tempm(2, 2)
      complex(real64), allocatable :: rmat(:)
      complex(real64), optional :: source_vector(2 * ltot * (ltot + 2))
      tranmat = reshape((/1, 1, 1, -1/), (/2, 2/))
      rmatdim = 0
      target_order = ntot
      source_order = ltot
      if (present(index_model)) then
         indexmodel = index_model
      else
         indexmodel = 1
      end if
      if (present(lr_transformation)) then
         lrtran = lr_transformation
      else
         lrtran = .true.
      end if
      if (present(make_symmetric)) then
         makesymmetric = make_symmetric
      else
         makesymmetric = .false.
      end if
      if (present(propagating_directions_only)) then
         propdir = propagating_directions_only
      else
         propdir = .false.
      end if
      if (present(source_vector)) then
         source_sum = .true.
         allocate (source_coefficient(2, ltot * (ltot + 2)))
         do n = 1, ltot
            do m = -n, n
               mn = n * (n + 1) + m
               do p = 1, 2
                  mnp = polarized_mode_index(m, n, p, ltot, indexmodel)
                  ctemp(p) = source_vector(mnp)
               end do
               if (lrtran) ctemp = matmul(tranmat, ctemp)
               source_coefficient(:, mn) = ctemp(:)
            end do
         end do
      else
         source_sum = .false.
      end if
      if (x .eq. 0.d0 .and. y .eq. 0.d0) then
         rho = 0.d0
         ephi = 1.d0
         azimuth_angle = 0.d0
      else
         rho = sqrt(x * x + y * y)
         ephi = cmplx(x, y, kind=real64) / rho
         azimuth_angle = atan2(y, x)
      end if
      if (rho .le. 0.0001d0) then
         max_azimuth_mode = 0
         rho = 0.d0
      else
         max_azimuth_mode = source_order + target_order
      end if
      source_z = sourcez
      target_z = targetz
      radial_distance = rho
      if (.not. source_sum) then
         if (makesymmetric) then
            max_azimuth_mode = 0
            rmatdim = 2 * axial_translation_size(ntot, ltot)
         else
            rmatdim = 4 * ntot * (ntot + 2) * ltot * (ltot + 2)
         end if
      end if
      source_layer = find_layer_index(sourcez)
      target_layer = find_layer_index(targetz)
      allocate (real_axis_limits(number_plane_boundaries + 1))
      real_axis_limits(1:number_plane_boundaries + 1) = abs(layer_ref_index(0:number_plane_boundaries))
      call sort_unique_real_values(number_plane_boundaries + 1, real_axis_limits, 1.d-10, number_limits)
      if (source_sum) then
         qtot = 2 * ntot * (ntot + 2)
      else
         qtot = 0
         do m = 0, ntot
            m1 = max(1, m)
            do k = -ltot, ltot
               if (abs(m - k) .gt. max_azimuth_mode) cycle
               k1 = max(abs(k), 1)
               do n = m1, ntot
                  do l = k1, ltot
                     qtot = qtot + 4
                  end do
               end do
            end do
         end do
      end if
      allocate (rmat(qtot))
      max_gf = 0.d0
      max_bf = 0.d0
      max_pi = 0.d0
      max_picon = 0.d0

      if (propdir) then
         ssc = s_scale_constant
         rail = real_axis_integration_limit
         s_scale_constant = 0.d0
         real_axis_integration_limit = dble(layer_ref_index(source_layer)) - 1.d-6
      end if

      call integrate_reflection_matrix_real_axis(qtot, rmat)

      if (propdir) then
         s_scale_constant = ssc
         real_axis_integration_limit = rail
      end if

      if (source_sum) then
         interactionmatrix(1:2 * ntot * (ntot + 2)) = 0.d0
         i = 0
         do m = -ntot, ntot
            m1 = max(1, abs(m))
            do n = m1, ntot
               mn = n * (n + 1) + m
               do p = 1, 2
                  i = i + 1
                  ctemp(p) = rmat(i)
               end do
               if (lrtran) ctemp = matmul(tranmat, ctemp) / 2.d0
               do p = 1, 2
                  mnp = polarized_mode_index(m, n, p, ntot, indexmodel)
                  interactionmatrix(mnp) = ctemp(p)
               end do
            end do
         end do
         deallocate (rmat, source_coefficient)
      else
         interactionmatrix(1:rmatdim) = 0.d0
         i = 0
         do m = 0, ntot
            m1 = max(1, m)
            do k = -ltot, ltot
               if (abs(m - k) .gt. max_azimuth_mode) cycle
               k1 = max(abs(k), 1)
               ephim = ephi**(k - m)
               do n = m1, ntot
                  do l = k1, ltot
                     do p = 1, 2
                        do q = 1, 2
                           i = i + 1
                           c2temp(p, q) = rmat(i) * ephim
                           if (m .ne. 0) c2tempm(p, q) = rmat(i) / ephim * ((-1)**(abs(m - k) + p + q))
                        end do
                     end do
                     if (lrtran) then
                        c2temp = matmul(tranmat, matmul(c2temp, tranmat)) / 2.d0
                        if (m .ne. 0) c2tempm = matmul(tranmat, matmul(c2tempm, tranmat)) / 2.d0
                     end if
                     do q = 1, 2
                        do p = 1, 2
                           if (makesymmetric) then
                              rmataddress = 2 * axial_translation_offset(m, ntot, ltot) + p &
                                            + 2 * (n - m1) + 2 * (ntot - m1 + 1) * (q - 1 + 2 * (l - m1))
                           else
                              mn = polarized_mode_index(m, n, p, ntot, indexmodel)
                              kl = polarized_mode_index(k, l, q, ltot, indexmodel)
                              rmataddress = mn + (kl - 1) * 2 * ntot * (ntot + 2)
                           end if
                           interactionmatrix(rmataddress) = c2temp(p, q)
                           if (m .ne. 0) then
                              if (makesymmetric) then
                                 rmataddress = 2 * axial_translation_offset(-m, ntot, ltot) + p &
                                               + 2 * (n - m1) + 2 * (ntot - m1 + 1) * (q - 1 + 2 * (l - m1))
                              else
                                 mn = polarized_mode_index(-m, n, p, ntot, indexmodel)
                                 kl = polarized_mode_index(-k, l, q, ltot, indexmodel)
                                 rmataddress = mn + (kl - 1) * 2 * ntot * (ntot + 2)
                              end if
                              interactionmatrix(rmataddress) = c2tempm(p, q)
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate (rmat)
      end if
      deallocate (real_axis_limits)
   end subroutine plane_boundary_interaction

   subroutine sphere_boundary_scattering(ntot1, rpos1, scoef1, ntot2, rpos2, scoef2, &
                                         targetz, qsca, lr_to_mode)
      implicit none
      logical :: lr2mode
      logical, optional :: lr_to_mode
      integer :: ntot1, ntot2, n, m, p, mnp, mn, tranmat(2, 2), nlimits, subdiv, nlimits0, ec
      real(real64) :: rpos1(3), rpos2(3), qsca(2), &
                      xyvec(2), limits(1:number_plane_boundaries + max_singular_points), s0, s1, errstep, targetz, riscale
      complex(real64) :: scoef1(2 * ntot1 * (ntot1 + 2), 2), scoef2(2 * ntot2 * (ntot2 + 2), 2), ctemp(2, 2), qmat(2)
      tranmat = reshape((/1, 1, 1, -1/), (/2, 2/))

      if (present(lr_to_mode)) then
         lr2mode = lr_to_mode
      else
         lr2mode = .true.
      end if
      allocate (source_coefficient_1(2, ntot1 * (ntot1 + 2), 2))
      allocate (source_coefficient_2(2, ntot2 * (ntot2 + 2), 2))
      do n = 1, ntot1
         do m = -n, n
            mn = n * (n + 1) + m
            do p = 1, 2
               mnp = polarized_mode_index(m, n, p, ntot1, 2)
               ctemp(p, :) = scoef1(mnp, :)
            end do
            if (lr2mode) ctemp = matmul(tranmat, ctemp)
            source_coefficient_1(:, mn, :) = ctemp(:, :)
         end do
      end do
      do n = 1, ntot2
         do m = -n, n
            mn = n * (n + 1) + m
            do p = 1, 2
               mnp = polarized_mode_index(m, n, p, ntot2, 2)
               ctemp(p, :) = scoef2(mnp, :)
            end do
            if (lr2mode) ctemp = matmul(tranmat, ctemp)
            source_coefficient_2(:, mn, :) = ctemp(:, :)
         end do
      end do
      source_layer = find_layer_index(rpos1(3))
      source_layer2 = find_layer_index(rpos2(3))
      target_layer = find_layer_index(targetz)
      source_order = ntot1
      source_order2 = ntot2
      source_z = rpos1(3)
      source_z2 = rpos2(3)
      target_z = targetz
      xyvec = rpos2(1:2) - rpos1(1:2)
      radial_distance = sqrt(sum(xyvec * xyvec))
      if (radial_distance .eq. 0.d0) then
         azimuth_angle = 0.d0
      else
         azimuth_angle = atan2(xyvec(2), xyvec(1))
      end if
      riscale = dble(layer_ref_index(target_layer))
      nlimits0 = 1
      limits(1) = 1.d0
      do n = 1, number_plane_boundaries + 1
         if (n - 1 .ne. target_layer .and. dble(layer_ref_index(n - 1)) .lt. riscale) then
            nlimits0 = nlimits0 + 1
            limits(nlimits0) = sqrt((1.d0 - dble(layer_ref_index(n - 1)) / riscale) &
                                    * (1.d0 + dble(layer_ref_index(n - 1)) / riscale))
         end if
      end do
      do n = 1, number_singular_points
         if (singular_points(n) .lt. riscale) then
            nlimits0 = nlimits0 + 1
            limits(nlimits0) = sqrt((1.d0 - singular_points(n) / riscale) &
                                    * (1.d0 + singular_points(n) / riscale))
         end if
      end do
      call sort_unique_real_values(nlimits0, limits, 1.d-10, nlimits)
      qsca = 0.d0
      integration_error = 0.d0
      s1 = 0.d0
      energy_kernel_region = 0
      do n = 1, nlimits
         qmat = 0.d0
         s0 = s1
         s1 = limits(n)
         subdiv = 0
         ec = 0
         call integrate_gauss_kronrod_adaptive(2, s0, s1, energy_kernel, qmat, subdiv, ec, &
                                           integration_error_epsilon, minimum_integration_spacing, maximum_integration_subdivisions)
         if (ec .eq. 1) error_codes(4) = 1
         if (ec .eq. 2) error_codes(3) = 1
         qsca = qsca + qmat
      end do
      nlimits0 = 0
      do n = 1, number_plane_boundaries + 1
         if (n - 1 .ne. target_layer .and. dble(layer_ref_index(n - 1)) .gt. riscale) then
            nlimits0 = nlimits0 + 1
            limits(nlimits0) = sqrt(dble(layer_ref_index(n - 1))**2 / riscale**2 + 1.d0)
         end if
      end do
      do n = 1, number_singular_points
         if (singular_points(n) .gt. riscale) then
            nlimits0 = nlimits0 + 1
            limits(nlimits0) = sqrt(1.d0 + singular_points(n)**2 / riscale**2)
         end if
      end do
      s1 = 0.d0
      energy_kernel_region = 1
      if (nlimits0 .gt. 0) then
         call sort_unique_real_values(nlimits0, limits, 1.d-10, nlimits)
         do n = 1, nlimits
            qmat = 0.d0
            s0 = s1
            s1 = limits(n)
            subdiv = 0
            ec = 0
            call integrate_gauss_kronrod_adaptive(2, s0, s1, energy_kernel, qmat, subdiv, ec, &
                                           integration_error_epsilon, minimum_integration_spacing, maximum_integration_subdivisions)
            if (ec .eq. 1) error_codes(4) = 1
            if (ec .eq. 2) error_codes(3) = 1
            qsca = qsca + qmat
         end do
      end if
      errstep = 1.d0
      do while (errstep .gt. integration_limit_epsilon .and. s1 .lt. real_axis_integration_limit)
         s0 = s1
         s1 = s1 + 0.5d0
         qmat = 0.d0
         subdiv = 0
         ec = 0
         call integrate_gauss_kronrod_adaptive(2, s0, s1, energy_kernel, qmat, subdiv, ec, &
                                           integration_error_epsilon, minimum_integration_spacing, maximum_integration_subdivisions)
         if (ec .eq. 1) error_codes(4) = 1
         if (ec .eq. 2) error_codes(3) = 1
         qsca = qsca + qmat
         errstep = abs(sum(qmat)) / abs(sum(qsca))
      end do
      deallocate (source_coefficient_1, source_coefficient_2)
   end subroutine sphere_boundary_scattering

   subroutine integrate_reflection_matrix_real_axis(qtot, rmat)
      implicit none
      integer :: qtot, n, limit, nseg, seg, n0, subdiv, ec
      real(real64) :: t1, t2, delt, dnorm, norm0, r, &
                      deltseg, dt1, dt2, errlim, t1t, t2t
      complex(real64) :: rmat(qtot), drmat(qtot)
      if (pole_integration) then
         t1 = 0.d0
         t2 = two_pi
         subdiv = 0.d0
         rmat = 0.d0
         ec = 0
         call integrate_gauss_kronrod_adaptive(qtot, t1, t2, real_axis_kernel, rmat, subdiv, ec, &
                                           integration_error_epsilon, minimum_integration_spacing, maximum_integration_subdivisions)
         return
      end if
      r = sqrt(sqrt(radial_distance**2 + (abs(source_z) + abs(target_z))**2))
      delt = .5d0 / r
      n0 = 1
      do n = 1, number_limits
         if (abs(abs(layer_ref_index(source_layer)) - real_axis_limits(n)) .le. 1.d-8) then
            n0 = n
            exit
         end if
      end do
      delt = min(delt, .5d0)
      delt = max(delt, minimum_initial_segment_size)
      s_sc1 = s_scale_constant / real_axis_limits(n0)
      s_sc2 = 0.d0
      rmat = 0.d0
      t2 = 0.d0
      do limit = 1, n0
         t1 = t2
         t2 = real_axis_limits(limit)
         nseg = ceiling((t2 - t1) / delt)
         deltseg = (t2 - t1) / dble(nseg)
         dt2 = t1
         do seg = 1, nseg
            dt1 = dt2
            dt2 = min(dt2 + deltseg, real_axis_integration_limit)
            drmat = 0.d0
            subdiv = 0
            ec = 0
            call integrate_gauss_kronrod_adaptive(qtot, dt1, dt2, real_axis_kernel, drmat, subdiv, ec, &
                                           integration_error_epsilon, minimum_integration_spacing, maximum_integration_subdivisions)
!if(mstm_global_rank.eq.0) write(*,'('' p1 '',2es12.4,2i6)') dt1,dt2,ec,subdiv
            if (ec .eq. 1) error_codes(4) = 1
            if (ec .eq. 2) error_codes(3) = 1
            rmat = rmat + drmat
            if (dt2 .ge. real_axis_integration_limit) return
         end do
      end do
      s_sc1 = 0.d0
      s_sc2 = s_scale_constant
      do limit = n0, number_limits - 1
         t1 = t2
         t2 = real_axis_limits(limit + 1)
         nseg = ceiling((t2 - t1) / delt)
         deltseg = (t2 - t1) / dble(nseg)
         dt2 = t1
         do seg = 1, nseg
            dt1 = dt2
            dt2 = min(dt2 + deltseg, real_axis_integration_limit)
            drmat = 0.d0
            subdiv = 0
            ec = 0
            call integrate_gauss_kronrod_adaptive(qtot, dt1, dt2, real_axis_kernel, drmat, subdiv, ec, &
                                           integration_error_epsilon, minimum_integration_spacing, maximum_integration_subdivisions)
!if(mstm_global_rank.eq.0) write(*,'('' p2 '',2es12.4,2i6)') dt1,dt2,ec,subdiv
            if (ec .eq. 1) error_codes(4) = 1
            if (ec .eq. 2) error_codes(3) = 1
            rmat = rmat + drmat
            if (dt2 .ge. real_axis_integration_limit) return
         end do
      end do
      delt = max(1.d0, minimum_initial_segment_size)
      do while (t2 .lt. real_axis_integration_limit)
         t1 = t2
         t2 = t2 + delt
         subdiv = 0
         t1t = t1
         t2t = t2
         ec = 0
         call integrate_gauss_kronrod_adaptive(qtot, t1t, t2t, real_axis_kernel, drmat, subdiv, ec, &
                                           integration_error_epsilon, minimum_integration_spacing, maximum_integration_subdivisions)
!if(mstm_global_rank.eq.0) write(*,'('' p3 '',2es12.4,2i6)') t1t,t2t,ec,subdiv
         if (ec .eq. 1) error_codes(4) = 1
         if (ec .eq. 2) error_codes(3) = 1
         rmat = rmat + drmat
         dnorm = mode_vector_norm(qtot, drmat)
         norm0 = max(mode_vector_norm(qtot, rmat), 0.01 * integration_limit_epsilon)
         errlim = dnorm / norm0
         if (errlim .lt. integration_limit_epsilon) return
      end do
      error_codes(2) = 1
   end subroutine integrate_reflection_matrix_real_axis

   subroutine boundary_energy_transfer(sinc, sdir, r, t, a, fres_r, fres_t)
      implicit none
      integer :: sdir, rdir, tdir
      real(real64) :: sinc, zref, ztra, r(2), t(2), a(2)
      complex(real64) :: riref, s, rkz, tkz, gfr(2, 2, 2), gft(2, 2, 2), ritra
      complex(real64), optional :: fres_r(2), fres_t(2)
      if (sdir .eq. 1) then
         riref = layer_ref_index(0)
         ritra = layer_ref_index(number_plane_boundaries)
         zref = -1.d-8
         ztra = plane_boundary_position(max(1, number_plane_boundaries)) + 1.d-8
         rdir = 2
         tdir = 1
      else
         riref = layer_ref_index(number_plane_boundaries)
         ritra = layer_ref_index(0)
         ztra = -1.d-8
         zref = plane_boundary_position(max(1, number_plane_boundaries)) + 1.d-8
         rdir = 1
         tdir = 2
      end if
      s = sinc
      if (number_plane_boundaries .eq. 0) then
         gfr = 0.d0
         gft = 1.d0
         rkz = 1.d0
         tkz = 1.d0
      else
         call layer_green_function(s, zref, zref, gfr, rkz, tkz)
         call layer_green_function(s, zref, ztra, gft, rkz, tkz)
      end if
      r(:) = abs(gfr(rdir, sdir, :))**2
      t(1) = abs(gft(tdir, sdir, 1))**2 * dble(tkz * conjg(ritra / riref) / rkz)
      t(2) = abs(gft(tdir, sdir, 2))**2 * dble(conjg(tkz * ritra / rkz / riref))
      a(:) = 1.d0 - r(:) - t(:)
      if (present(fres_r)) then
         fres_r(:) = gfr(rdir, sdir, :)
      end if
      if (present(fres_t)) then
         fres_t(:) = gft(tdir, sdir, :)
      end if
   end subroutine boundary_energy_transfer

   subroutine initialize_incident_field(alpha, sinc, sdir)
      implicit none
      integer :: sdir, ssign, p, k, q, klq
      real(real64) :: alpha, sinc, targetz, sourcez
      complex(real64) s, pmnp(6, 2), riinc, skz, tkz, gfs(2, 2, 2)
      incident_lateral_vector = (/sinc * cos(alpha), sinc * sin(alpha)/)
      if (number_plane_boundaries .eq. 0) then
         incident_field_scale = 1.d0
         incident_field_boundary = 0.d0
         return
      end if
      if (sdir .eq. 1) then
         riinc = layer_ref_index(0)
      else
         riinc = layer_ref_index(number_plane_boundaries)
      end if
      if (sinc .gt. dble(riinc)) then
         s = sinc
         if (sdir .eq. 1) then
            sourcez = -1.d-8
         else
            sourcez = plane_boundary_position(number_plane_boundaries) + 1.d-8
         end if
         targetz = 0.5d0 * plane_boundary_position(number_plane_boundaries)
         call layer_green_function(s, sourcez, targetz, gfs, skz, tkz)
         call generate_plane_wave_coefficients(alpha, tkz, 1, pmnp)
         do p = 1, 2
            do k = -1, 1
               do q = 1, 2
                  klq = polarized_mode_index(k, 1, q, 1, 2)
                  ssign = (-1)**(k + q + p)
                  pmnp(klq, p) = pmnp(klq, p) * (gfs(1, sdir, p) + ssign * gfs(2, sdir, p))
               end do
            end do
            incident_field_scale(p) = sqrt(sum(abs(pmnp(1:3, p)**2)))
         end do
      else
         if (sdir .eq. 1) then
            sourcez = bot_boundary
         else
            sourcez = top_boundary
         end if
         incident_field_scale(:) = 1.d0
      end if
      incident_field_scale(:) = maxval(incident_field_scale(1:2))
      incident_field_boundary = sourcez
   end subroutine initialize_incident_field

   subroutine layer_plane_wave_coefficients(alpha, sinc, sdir, rpos, nodr, pmnp, include_direct)
      implicit none
      logical, optional :: include_direct
      logical :: incdir, evanescent
      integer :: p, incregion, sdir, nodr, nblk, layer
      real(real64) :: alpha, ca, sa, rpos(3), sourcez, targetz, sinc
      complex(real64) :: pmnp(2 * nodr * (nodr + 2), 2), riinc, cbinc, &
                         s, phasefaclat, skz, tkz, gfs(2, 2, 2)
      complex(real64), allocatable :: pmnpinc(:, :), pmnpup(:, :), &
                                      pmnpdn(:, :), pmnptot(:, :)
      if (present(include_direct)) then
         incdir = include_direct
      else
         incdir = .true.
      end if
      nblk = 2 * nodr * (nodr + 2)
      ca = cos(alpha)
      sa = sin(alpha)
      layer = find_layer_index(rpos(3))
      if (sdir .eq. 1) then
         riinc = layer_ref_index(0)
         incregion = 0
      else
         riinc = layer_ref_index(number_plane_boundaries)
         incregion = number_plane_boundaries
      end if
      if (sinc .gt. dble(riinc)) then
         evanescent = .true.
         incdir = .false.
      else
         evanescent = .false.
      end if
      sourcez = incident_field_boundary
      s = sinc
      cbinc = sqrt((1.d0 - sinc / riinc) * (1.d0 + sinc / riinc)) * (3 - 2 * sdir)
      pmnp = 0.d0
      targetz = rpos(3)
      phasefaclat = exp((0.d0, 1.d0) * s * (ca * rpos(1) + sa * rpos(2)))
      allocate (pmnptot(nblk, 2))
      pmnptot = 0.d0
      if (layer .eq. incregion .and. incdir) then
         allocate (pmnpinc(nblk, 2))
         call generate_plane_wave_coefficients(alpha, cbinc, nodr, pmnpinc)
         pmnptot = pmnpinc * exp((0.d0, 1.d0) * riinc * cbinc * (rpos(3) - sourcez))
         deallocate (pmnpinc)
      end if
      if (number_plane_boundaries .gt. 0) then
         call layer_green_function(s, sourcez, targetz, gfs, skz, tkz)
         allocate (pmnpup(nblk, 2), pmnpdn(nblk, 2))
         call generate_plane_wave_coefficients(alpha, tkz, nodr, pmnpup)
         call generate_plane_wave_coefficients(alpha, -tkz, nodr, pmnpdn)
         do p = 1, 2
            pmnptot(:, p) = pmnptot(:, p) + pmnpup(:, p) * gfs(1, sdir, p) &
                            + pmnpdn(:, p) * gfs(2, sdir, p)
         end do
         deallocate (pmnpup, pmnpdn)
      else
         tkz = abs(cbinc)
         skz = abs(cbinc)
      end if
      do p = 1, 2
         pmnp(:, p) = pmnptot(:, p) * phasefaclat / incident_field_scale(p)
      end do
      deallocate (pmnptot)
   end subroutine layer_plane_wave_coefficients

   subroutine layer_vector_spherical_harmonics(s, alpha, targetz, tdir, rpos, nodr, pmnp)
      implicit none
      integer :: p, tdir, nodr, nblk, slayer, tlayer, k, l, q, klq, ssign
      real(real64) :: alpha, ca, sa, rpos(3), sourcez, targetz
      complex(real64) :: pmnp(2 * nodr * (nodr + 2), 2), targetri, sourceri, &
                         s, phasefaclat, skz, tkz, gfs(2, 2, 2)

      nblk = 2 * nodr * (nodr + 2)
      ca = cos(alpha)
      sa = sin(alpha)
      slayer = find_layer_index(rpos(3))
      tlayer = find_layer_index(targetz)
      targetri = layer_ref_index(tlayer)
      sourceri = layer_ref_index(slayer)
      sourcez = rpos(3)
      phasefaclat = exp(-(0.d0, 1.d0) * s * (ca * rpos(1) + sa * rpos(2)))
      call layer_green_function(s, sourcez, targetz, gfs, skz, tkz, include_direct=.true.)
      call generate_plane_wave_coefficients(alpha, conjg(skz), nodr, pmnp, lr_tran=.false.)
      do p = 1, 2
         do l = 1, nodr
            do k = -l, l
               do q = 1, 2
                  klq = polarized_mode_index(k, l, q, nodr, 2)
                  ssign = (-1)**(k + l + q - 1 + p)
                  pmnp(klq, p) = conjg(pmnp(klq, p)) * (gfs(tdir, 1, p) + ssign * gfs(tdir, 2, p))
               end do
            end do
         end do
      end do
      pmnp = pmnp * phasefaclat / 4.d0 / sourceri / sourceri / skz
   end subroutine layer_vector_spherical_harmonics

end module surface
