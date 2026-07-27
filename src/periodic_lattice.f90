module periodic_lattice_subroutines
   use numconstants
   use specialfuncs
   use surface_subroutines
   use mpidefs
   implicit none
   logical, target :: periodic_lattice, time_it, phase_shift_form, finite_lattice
   integer :: pl_max_subdivs, pl_rs_nmax, pl_error_codes(6), pl_fs_method, pl_rs_imax, &
              q1d_number_segments, q2d_number_segments, s_max_q2
   integer, target :: pl_integration_method
   real(8) :: lattice_integration_segment, pl_rs_eps, time_count(4), time_0
   real(8), target :: cell_width(2), rs_dz_min, pl_integration_error_epsilon, &
                      pl_integration_limit_epsilon
   integer :: qkernel_nodr, qkernel_integration_model
   real(8) :: qkernel_x, qkernel_y, qkernel_z, qkernel_width, &
              qkernel_k0y, qkernel_kz
   complex(8) :: qkernel_ref_index
   data lattice_integration_segment/1.d0/
   data pl_rs_nmax, pl_rs_eps, pl_rs_imax/200, 1.d-7, 0/
   data time_it, rs_dz_min, s_max_q2/.true., 100.d0, 100/
   data pl_integration_error_epsilon, pl_integration_limit_epsilon, pl_integration_method/1.d-7, 1.d-9, 1/
   data time_count/0.d0, 0.d0, 0.d0, 0.d0/
   data phase_shift_form, finite_lattice/.false., .false./

contains

   subroutine plane_boundary_lattice_interaction(nodrt, nodrs, x0, y0, zt, zs, &
                                                 matrix, source_vector, include_source, lr_transformation, index_model)
      implicit none
      logical :: incsrc, lrtran, rsincsrc, fsincsrc
      logical, optional :: include_source, lr_transformation
      integer :: nodrt, nodrs, i, ix, iy, nmax, nterms, n, p, q, imodl, m, l, k, mnp, klq, mn, kl, tranmat(2, 2), &
                 slay, tlay, np(2)
      integer, optional :: index_model
      real(8) :: x0, y0, zt, zs, w(2), wx, wy, k0x, k0y, kconst, kx, ky, eps, cerr, asum, asum0, x, y, wcrit
      complex(8) :: matrix(*), csum(2, 2), ri
      complex(8), optional :: source_vector(2 * nodrs * (nodrs + 2))
      complex(8), allocatable :: kernel(:, :, :, :), dkernel(:, :, :, :), fsmat(:, :, :), rsmat(:, :)

      tranmat = reshape((/1, 1, 1, -1/), (/2, 2/))
      if (present(include_source)) then
         incsrc = include_source
      else
         incsrc = .false.
      end if
      if (present(lr_transformation)) then
         lrtran = lr_transformation
      else
         lrtran = .true.
      end if
      if (present(index_model)) then
         imodl = index_model
      else
         imodl = 2
      end if
      slay = layer_id(zs)
      tlay = layer_id(zt)
      if (.not. plane_surface_present) then
         if (present(source_vector)) then
            call free_space_lattice_translation_matrix(nodrt, nodrs, (/x0, y0, zt - zs/), &
                                                       layer_ref_index(0), matrix, source_vector=source_vector, &
                                                       include_source=incsrc, lr_transformation=lrtran, index_model=imodl)
         else
            call free_space_lattice_translation_matrix(nodrt, nodrs, (/x0, y0, zt - zs/), &
                                                       layer_ref_index(0), matrix, &
                                                       include_source=incsrc, lr_transformation=lrtran, index_model=imodl)
         end if
         return
      elseif (slay .eq. tlay) then
         if (present(source_vector)) then
            call common_layer_lattice_translation_matrix(nodrt, nodrs, x0, y0, zt, zs, &
                                                         matrix, source_vector=source_vector, &
                                                         include_source=incsrc, lr_transformation=lrtran, index_model=imodl)
         else
            call common_layer_lattice_translation_matrix(nodrt, nodrs, x0, y0, zt, zs, &
                                                         matrix, include_source=incsrc, lr_transformation=lrtran, index_model=imodl)
         end if
         return
      end if
      w = cell_width
      wcrit = min(rs_dz_min, minval(cell_width) / 2.d0)
      k0x = incident_lateral_vector(1)
      k0y = incident_lateral_vector(2)
      np(1) = floor((x0 + cell_width(1) / 2.d0) / cell_width(1))
      np(2) = floor((y0 + cell_width(2) / 2.d0) / cell_width(2))
      x = x0 - cell_width(1) * dble(np(1))
      y = y0 - cell_width(2) * dble(np(2))
      nmax = pl_rs_nmax
      eps = pl_rs_eps
      if (incsrc) then
         ri = layer_ref_index(slay)
         if (slay .eq. tlay) then
            pl_fs_method = 0
            rsincsrc = .false.
            fsincsrc = .true.
         else
            pl_fs_method = 1
            rsincsrc = .true.
            fsincsrc = .false.
         end if
      else
         rsincsrc = .false.
         fsincsrc = .false.
      end if
      allocate (kernel(2, nodrt * (nodrt + 2), 2, nodrs * (nodrs + 2)), dkernel(2, nodrt * (nodrt + 2), 2, nodrs * (nodrs + 2)))
      wx = w(1)
      wy = w(2)
      kconst = 8.d0 * pi * pi / cell_width(1) / cell_width(2)
      kernel = 0.d0
      kx = incident_lateral_vector(1)
      ky = incident_lateral_vector(2)
      call plane_boundary_lattice_kernel(nodrt, nodrs, kx, ky, x, y, zt, zs, kernel, include_source=rsincsrc)
      do n = 1, pl_rs_nmax
         dkernel = 0.d0
         do i = 0, 8 * n - 1
            q = i / (2 * n)
            p = i - 2 * q * n
            if (q .eq. 0) then
               ix = n
               iy = -n + p
            elseif (q .eq. 1) then
               ix = n - p
               iy = n
            elseif (q .eq. 2) then
               ix = -n
               iy = n - p
            else
               ix = -n + p
               iy = -n
            end if
            kx = 2.d0 * pi * dble(ix) / cell_width(1) + incident_lateral_vector(1)
            ky = 2.d0 * pi * dble(iy) / cell_width(2) + incident_lateral_vector(2)
            call plane_boundary_lattice_kernel(nodrt, nodrs, kx, ky, x, y, zt, zs, dkernel, include_source=rsincsrc)
         end do
         kernel = kernel + dkernel
         asum0 = sum(abs(kernel))
         asum = sum(abs(dkernel))
         cerr = asum / asum0
         if (cerr .lt. pl_rs_eps) exit
      end do
      nterms = n
      pl_rs_imax = nterms
      if (nterms .ge. pl_rs_nmax) pl_error_codes(3) = 1
      deallocate (dkernel)
      allocate (rsmat(2 * nodrt * (nodrt + 2), 2 * nodrs * (nodrs + 2)))
      do l = 1, nodrs
         do k = -l, l
            kl = l * (l + 1) + k
            do n = 1, nodrt
               do m = -n, n
                  mn = n * (n + 1) + m
                  do q = 1, 2
                     do p = 1, 2
                        csum(p, q) = kernel(p, mn, q, kl)
                     end do
                  end do
                  if (lrtran) then
                     csum = matmul(tranmat, matmul(csum, tranmat)) / 2.d0
                  end if
                  do q = 1, 2
                     klq = amnpaddress(k, l, q, nodrs, imodl)
                     do p = 1, 2
                        mnp = amnpaddress(m, n, p, nodrt, imodl)
                        rsmat(mnp, klq) = csum(p, q) * kconst
                     end do
                  end do
               end do
            end do
         end do
      end do
      deallocate (kernel)
      if (slay .eq. tlay) then
         allocate (fsmat(nodrt * (nodrt + 2), nodrs * (nodrs + 2), 2))
         ri = layer_ref_index(slay)
         call free_space_lattice_translation_matrix(nodrt, nodrs, (/x, y, zt - zs/), &
                                                    ri, fsmat, include_source=fsincsrc, lr_transformation=lrtran, index_model=imodl)
         do l = 1, nodrs
            do k = -l, l
               kl = amnaddress(k, l, nodrs, imodl)
               do n = 1, nodrt
                  do m = -n, n
                     mn = amnaddress(m, n, nodrt, imodl)
                     do p = 1, 2
                        klq = amnpaddress(k, l, p, nodrs, imodl)
                        mnp = amnpaddress(m, n, p, nodrt, imodl)
                        rsmat(mnp, klq) = rsmat(mnp, klq) + fsmat(mn, kl, p)
                     end do
                  end do
               end do
            end do
         end do
         deallocate (fsmat)
      end if
      if (present(source_vector)) then
         matrix(1:2 * nodrt * (nodrt + 2)) &
            = matmul(rsmat(:, :), source_vector(:)) * exp((0.d0, 1.d0) * sum(incident_lateral_vector * dble(np) * cell_width))
      else
         matrix(1:4 * nodrt * (nodrt + 2) * nodrs * (nodrs + 2)) &
            = reshape(rsmat, (/4 * nodrt * (nodrt + 2) * nodrs * (nodrs + 2)/)) * exp((0.d0, 1.d0) &
                                                                             * sum(incident_lateral_vector * dble(np) * cell_width))
      end if
      deallocate (rsmat)
   end subroutine plane_boundary_lattice_interaction

   subroutine free_space_lattice_translation_matrix(nodrt, nodrs, rpos0, &
                                                    ri, matrix, source_vector, include_source, lr_transformation, index_model)
      implicit none
      logical :: incsrc, lrtran
      logical, optional :: include_source, lr_transformation
      integer :: nodrt, nodrs, wmax, l, m, mn, n, np(2), &
                 m1m, k, kl, nn1, ll1, v, w, wmin, vw, imodl, nterms, klp, p
      integer, optional :: index_model
      real(8) :: rpos(3), vc1(0:nodrs + nodrt), vc2(0:nodrs + nodrt), rpos0(3), wcrit
      complex(8) :: ci, c, a, b, ri, pshift, &
                    ysum(0:(nodrt + nodrs) * (nodrt + nodrs + 2)), &
                    swf(0:(nodrt + nodrs) * (nodrt + nodrs + 2)), matrix(*)
      complex(8), allocatable :: fsmat(:, :, :)
      complex(8), optional :: source_vector(nodrs * (nodrs + 2) * 2)
      data ci/(0.d0, 1.d0)/
      if (present(include_source)) then
         incsrc = include_source
      else
         incsrc = .true.
      end if
      if (present(lr_transformation)) then
         lrtran = lr_transformation
      else
         lrtran = .true.
      end if
      if (present(index_model)) then
         imodl = index_model
      else
         imodl = 2
      end if
      wcrit = min(rs_dz_min, minval(cell_width) / 2.d0)
      rpos(3) = rpos0(3)
      np(1) = floor((rpos0(1) + cell_width(1) / 2.d0) / cell_width(1))
      np(2) = floor((rpos0(2) + cell_width(2) / 2.d0) / cell_width(2))
      rpos(1) = rpos0(1) - dble(np(1)) * cell_width(1)
      rpos(2) = rpos0(2) - dble(np(2)) * cell_width(2)
      pshift = exp((0.d0, 1.d0) * sum(incident_lateral_vector * dble(np) * cell_width))
!
! phase shift form=.true. scales scattering coefficients with lateral phase shift.   Used only for testing purposes.
!
      if (phase_shift_form) pshift = exp((0.d0, -1.d0) * sum(incident_lateral_vector * rpos(1:2)))
      pl_max_subdivs = 0
      allocate (fsmat(nodrt * (nodrt + 2), nodrs * (nodrs + 2), 2))
!
! finite lattice uses cell-centered translation matrix H^i-j instead of L.   Included only for testing/expeimentation purposes
!
      if (finite_lattice) then
         call gentranmatrix(nodrs, nodrt, translation_vector=rpos, &
                            refractive_index=(/ri, ri/), ac_matrix=fsmat, vswf_type=3, &
                            index_model=2)
      else
         wmax = nodrt + nodrs
         ysum = 0.d0
         if (abs(rpos(3)) .lt. wcrit) then
            pl_fs_method = 0
            if (time_it) time_0 = mstm_mpi_wtime()
            call swf_lattice_sum(wmax, rpos(1), rpos(2), rpos(3), cell_width, incident_lateral_vector, &
                                 ri, ysum, include_source=incsrc)
            if (time_it) time_count(1) = mstm_mpi_wtime() - time_0 + time_count(1)
            if (pl_max_subdivs .ge. maximum_integration_subdivisions) pl_error_codes(1) = 1
         else
            pl_fs_method = 1
            if (time_it) time_0 = mstm_mpi_wtime()
            call reciprocal_space_swf_lattice_sum(wmax, rpos(1), rpos(2), rpos(3), cell_width, &
                                                  incident_lateral_vector, ri, pl_rs_nmax, pl_rs_eps, nterms, ysum)
            if (time_it) time_count(2) = mstm_mpi_wtime() - time_0 + time_count(2)
            if (.not. incsrc) then
               call scalar_wave_function(wmax, 3, rpos(1), rpos(2), rpos(3), ri, swf)
               ysum = ysum - swf
            end if
            if (nterms .ge. pl_rs_nmax) pl_error_codes(2) = 1
         end if
         do l = 1, nodrs
            ll1 = l * (l + 1)
            do n = 1, nodrt
               nn1 = n * (n + 1)
               wmax = n + l
               call vcfunc(-1, n, 1, l, vc2)
               c = -sqrt(4.d0 * pi) * ci**(n - l) * fnr(n + n + 1) * fnr(l + l + 1)
               do k = -l, l
                  kl = amnaddress(k, l, nodrs, imodl)
                  do m = -n, n
                     m1m = (-1)**m
                     mn = amnaddress(m, n, nodrt, imodl)
                     v = k - m
                     call vcfunc(-m, n, k, l, vc1)
                     wmin = max(abs(v), abs(n - l))
                     a = 0.
                     b = 0.
                     do w = wmax, wmin, -1
                        vw = w * (w + 1) + v
                        if (mod(wmax - w, 2) .eq. 0) then
                           a = a + (ci**w) * vc1(w) * vc2(w) * ysum(vw) / sqrt(dble(w + w + 1))
                        else
                           b = b + (ci**w) * vc1(w) * vc2(w) * ysum(vw) / sqrt(dble(w + w + 1))
                        end if
                     end do
                     if (lrtran) then
                        fsmat(mn, kl, 1) = m1m * c * (a + b)
                        fsmat(mn, kl, 2) = m1m * c * (a - b)
                     else
                        fsmat(mn, kl, 1) = m1m * c * a
                        fsmat(mn, kl, 2) = m1m * c * b
                     end if
                  end do
               end do
            end do
         end do
      end if
      if (present(source_vector)) then
         do p = 1, 2
            do n = 1, nodrt
               do m = -n, n
                  mn = amnaddress(m, n, nodrt, imodl)
                  c = 0.d0
                  do l = 1, nodrs
                     do k = -l, l
                        kl = amnaddress(k, l, nodrs, imodl)
                        klp = amnpaddress(k, l, p, nodrs, imodl)
                        c = c + fsmat(mn, kl, p) * source_vector(klp)
                     end do
                  end do
                  mn = amnpaddress(m, n, p, nodrt, imodl)
                  matrix(mn) = c * pshift
               end do
            end do
         end do
      else
matrix(1:2 * nodrt * (nodrt + 2) * nodrs * (nodrs + 2)) = reshape(fsmat, (/2 * nodrt * (nodrt + 2) * nodrs * (nodrs + 2)/)) * pshift
      end if
      deallocate (fsmat)
   end subroutine free_space_lattice_translation_matrix

   subroutine common_layer_lattice_translation_matrix(nodrt, nodrs, x0, y0, zt, zs, &
                                                      matrix, source_vector, include_source, lr_transformation, index_model)
      implicit none
      logical :: incsrc, lrtran
      logical, optional :: include_source, lr_transformation
      integer :: nodrt, nodrs, slay, nodrw, nterms, nblk, sdirs(2), tdirs(2), n, i, q, p, ix, iy, &
                 ll1, nn1, wmax, l, k, m, mn, kl, mnp, klq, tsign, ssign, iw, wmin, v, np(2), imodl, sdir, tdir, &
                 tranmat(2, 2), vw, w, nelem
      integer, optional :: index_model
      real(8) :: x, y, zt, zs, kx, ky, vc1(0:nodrt + nodrs), &
                 vcp1m1(0:nodrt + nodrs), vcm1m1(0:nodrt + nodrs), x0, y0, a0mag
      real(8), allocatable :: asum(:), asum0(:), cerr(:)
      complex(8) :: ri, c, tmat(-1:2), csum(2, 2), matrix(*)
      complex(8), allocatable :: qsum(:, :, :, :), dqsum(:, :, :, :), mat(:, :), fsmat(:, :, :)
      complex(8), optional :: source_vector(2 * nodrs * (nodrs + 2))
      tranmat = reshape((/1, 1, 1, -1/), (/2, 2/))
      if (present(include_source)) then
         incsrc = include_source
      else
         incsrc = .false.
      end if
      if (present(lr_transformation)) then
         lrtran = lr_transformation
      else
         lrtran = .true.
      end if
      if (present(index_model)) then
         imodl = index_model
      else
         imodl = 2
      end if
      np(1) = floor((x0 + cell_width(1) / 2.d0) / cell_width(1))
      np(2) = floor((y0 + cell_width(2) / 2.d0) / cell_width(2))
      x = x0 - cell_width(1) * dble(np(1))
      y = y0 - cell_width(2) * dble(np(2))
      slay = layer_id(zs)
      ri = layer_ref_index(slay)
      nodrw = nodrs + nodrt
      nblk = nodrw * (nodrw + 2)
      if (slay .eq. 0) then
         sdirs = (/1, 1/)
         tdirs = (/2, 2/)
      elseif (slay .eq. number_plane_boundaries) then
         sdirs = (/2, 2/)
         tdirs = (/1, 1/)
      else
         sdirs = (/1, 2/)
         tdirs = (/1, 2/)
      end if
      nelem = 3 * (nblk + 1) * (tdirs(2) - tdirs(1) + 1) * (sdirs(2) - sdirs(1) + 1)
      allocate (qsum(-1:1, 0:nblk, tdirs(1):tdirs(2), sdirs(1):sdirs(2)), &
                dqsum(-1:1, 0:nblk, tdirs(1):tdirs(2), sdirs(1):sdirs(2)), asum(nelem), asum0(nelem), cerr(nelem))
      kx = incident_lateral_vector(1)
      ky = incident_lateral_vector(2)
      qsum = 0.d0
      call common_layer_lattice_kernel(nodrw, kx, ky, x, y, zt, zs, tdirs, sdirs, qsum)
      do n = 1, pl_rs_nmax
         dqsum = 0.d0
         do i = 0, 8 * n - 1
            q = i / (2 * n)
            p = i - 2 * q * n
            if (q .eq. 0) then
               ix = n
               iy = -n + p
            elseif (q .eq. 1) then
               ix = n - p
               iy = n
            elseif (q .eq. 2) then
               ix = -n
               iy = n - p
            else
               ix = -n + p
               iy = -n
            end if
            kx = 2.d0 * pi * dble(ix) / cell_width(1) + incident_lateral_vector(1)
            ky = 2.d0 * pi * dble(iy) / cell_width(2) + incident_lateral_vector(2)
            call common_layer_lattice_kernel(nodrw, kx, ky, x, y, zt, zs, tdirs, sdirs, dqsum)
         end do
         qsum = qsum + dqsum
         asum0 = reshape(abs(qsum), (/nelem/))
         asum = reshape(abs(dqsum), (/nelem/))
         a0mag = maxval(asum0)
         cerr = 0.d0
         do i = 1, nelem
            if (asum0(i) .gt. 1.d-15 * a0mag) cerr(i) = asum(i) / asum0(i)
         end do
         if (maxval(cerr) .lt. pl_rs_eps) exit
      end do
      nterms = n
      if (nterms .ge. pl_rs_nmax) pl_error_codes(3) = 1
      deallocate (dqsum)
      qsum = qsum / cell_width(1) / cell_width(2) * pi**1.5d0
      allocate (mat(2 * nodrt * (nodrt + 2), 2 * nodrs * (nodrs + 2)))
      do l = 1, nodrs
         ll1 = l * (l + 1)
         do n = 1, nodrt
            nn1 = n * (n + 1)
            wmax = n + l
            call vcfunc(1, n, -1, l, vcp1m1)
            call vcfunc(-1, n, -1, l, vcm1m1)
            c = (0.d0, 1.d0)**(n - l) * fnr(n + n + 1) * fnr(l + l + 1)
            do k = -l, l
               kl = ll1 + k
               do m = -n, n
                  mn = nn1 + m
                  v = k - m
                  call vcfunc(-m, n, k, l, vc1)
                  wmin = max(abs(v), abs(n - l))
                  csum = 0.d0
                  do q = 1, 2
                     klq = amnpaddress(k, l, q, nodrs, imodl)
                     do p = 1, 2
                        mnp = amnpaddress(m, n, p, nodrt, imodl)
                        tmat = 0.d0
                        do sdir = sdirs(1), sdirs(2)
                           if (sdir .eq. 1) then
                              ssign = 1
                           else
                              ssign = (-1)**(k + l + q - 1)
                           end if
                           do tdir = tdirs(1), tdirs(2)
                              if (tdir .eq. 1) then
                                 tsign = ssign
                              else
                                 tsign = (-1)**(m + n + p - 1) * ssign
                              end if
                              do w = wmax, wmin, -1
                                 iw = (-1)**w
                                 vw = w * (w + 1) + v
                                 if (w .ge. 2) then
                                    tmat(-1) = tmat(-1) + vc1(w) * vcm1m1(w) * qsum(-1, vw, tdir, sdir) * tsign
                                    tmat(1) = tmat(1) + iw * vc1(w) * vcm1m1(w) * qsum(1, vw, tdir, sdir) * tsign
                                 end if
                                 tmat(0) = tmat(0) + vc1(w) * vcp1m1(w) * qsum(0, vw, tdir, sdir) * tsign
                                 tmat(2) = tmat(2) + iw * vc1(w) * vcp1m1(w) * qsum(0, vw, tdir, sdir) * tsign
                              end do
                           end do
                        end do
                        csum(p, q) = (-1)**m * c * ((-1)**q * tmat(-1) + (-1)**(p + q) * tmat(0) &
                                                    + (-1)**(n + l) * tmat(2) + (-1)**(n + l + p) * tmat(1))
                     end do
                  end do
                  if (lrtran) then
                     csum = matmul(tranmat, matmul(csum, tranmat)) / 2.d0
                  end if
                  do q = 1, 2
                     klq = amnpaddress(k, l, q, nodrs, imodl)
                     do p = 1, 2
                        mnp = amnpaddress(m, n, p, nodrt, imodl)
                        mat(mnp, klq) = csum(p, q)
                     end do
                  end do
               end do
            end do
         end do
      end do
      deallocate (qsum)
      allocate (fsmat(nodrt * (nodrt + 2), nodrs * (nodrs + 2), 2))
      call free_space_lattice_translation_matrix(nodrt, nodrs, (/x, y, zt - zs/), &
                                                 ri, fsmat, include_source=incsrc, lr_transformation=lrtran, index_model=imodl)
      do l = 1, nodrs
         do k = -l, l
            kl = amnaddress(k, l, nodrs, imodl)
            do n = 1, nodrt
               do m = -n, n
                  mn = amnaddress(m, n, nodrt, imodl)
                  do p = 1, 2
                     klq = amnpaddress(k, l, p, nodrs, imodl)
                     mnp = amnpaddress(m, n, p, nodrt, imodl)
                     mat(mnp, klq) = mat(mnp, klq) + fsmat(mn, kl, p)
                  end do
               end do
            end do
         end do
      end do
      deallocate (fsmat)
      if (present(source_vector)) then
         matrix(1:2 * nodrt * (nodrt + 2)) &
            = matmul(mat(:, :), source_vector(:)) * exp((0.d0, 1.d0) * sum(incident_lateral_vector * dble(np) * cell_width))
      else
         matrix(1:4 * nodrt * (nodrt + 2) * nodrs * (nodrs + 2)) &
            = reshape(mat, (/4 * nodrt * (nodrt + 2) * nodrs * (nodrs + 2)/)) * exp((0.d0, 1.d0) &
                                                                             * sum(incident_lateral_vector * dble(np) * cell_width))
      end if
      deallocate (mat)
   end subroutine common_layer_lattice_translation_matrix

   subroutine swf_lattice_sum(nodr, x, y, z, w, k0, ri, swfsum, include_source)
      implicit none
      logical :: addsource
      logical, optional :: include_source
      integer :: nodr, n, nn1, k, m, mn
      real(8) :: x, y, z, w(2), k0(2), drot(-nodr:nodr, 0:nodr * (nodr + 2))
      complex(8) :: ri, swfsum(0:nodr * (nodr + 2)), swf(0:nodr * (nodr + 2)), act(-nodr:nodr), csum

      if (present(include_source)) then
         if (x .ne. 0.d0 .or. y .ne. 0.d0 .or. z .ne. 0.d0) then
            addsource = include_source
         else
            addsource = .false.
         end if
      else
         addsource = .false.
      end if
      call swfyzlatticesum(nodr, z, -y, x, (/w(2), w(1)/), -k0(2), k0(1), ri, swfsum)
      call rotcoef(0.d0, nodr, nodr, drot)
      do n = 0, nodr
         nn1 = n * (n + 1)
         do k = -n, n
            act(k) = swfsum(nn1 + k)
         end do
         do m = -n, n
            mn = nn1 + m
            csum = 0.
            do k = -n, n
               csum = csum + drot(m, nn1 - k) * act(k)
            end do
            swfsum(nn1 + m) = ((-1)**n) * csum
         end do
      end do
      if (addsource) then
         call scalar_wave_function(nodr, 3, x, y, z, ri, swf)
         swfsum = swfsum + swf
      end if
   end subroutine swf_lattice_sum

   subroutine swfyzlatticesum(nodr, x, y, z, w, k0y, k0z, ri, swfyzsum)
      implicit none
      logical :: convrg
      integer :: nodr, m, n, s, ntz, l, smaxp, smaxn
      real(8) :: x, k0y, k0z, w(2), y, z, kz, mag, terr(0:nodr * (nodr + 2)), derr(0:nodr * (nodr + 2))
      complex(8) :: q1d(0:nodr * (nodr + 2)), q2d(-nodr:nodr), &
                    swfyzsum(0:nodr * (nodr + 2)), ri, ci, &
                    c, ct, ymn(0:nodr * (nodr + 2))
      complex(8), allocatable :: sum1(:, :)
      data ci/(0.d0, 1.d0)/
      ntz = max(s_max_q2, ceiling(w(2) / 2.d0 / pi))
      allocate (sum1(-nodr:nodr, 0:ntz))
      smaxp = ntz
      smaxn = -ntz
      swfyzsum = 0.d0

      call q1dbnosource(nodr, x, y, z, w(2), k0z, ri, q1d)
      do n = 0, nodr
         do m = -n, n
            swfyzsum(m + n * (n + 1)) = q1d(m + n * (n + 1))
         end do
      end do

      sum1 = 0.d0
      convrg = .false.
      do s = 0, ntz
         kz = k0z + 2.d0 * pi * dble(s) / w(2)
         call q2db(nodr, x, y, w(1), k0y, kz, ri, q2d)
         sum1(:, s) = q2d
         ct = kz / ri
         call crotcoef(ct, 0, nodr, ymn)
         terr = 0.d0
         derr = abs(swfyzsum(:))
         do n = 0, nodr
            c = -((0.d0, -1.d0)**n) / w(2) / ri / sqrt(4.d0 * pi / dble(n + n + 1))
            do m = -n, n
               l = m + n * (n + 1)
               swfyzsum(l) = swfyzsum(l) &
                             + c * ymn(l) * sum1(m, s) * exp((0.d0, 1.d0) * kz * z)
            end do
         end do
         terr = abs(swfyzsum(:))
         if (abs(kz) .gt. 1.d0) then
            do n = 0, nodr * (nodr + 2)
               if (terr(n) .ne. 0.d0) then
                  terr(n) = abs(terr(n) - derr(n)) / terr(n)
               end if
            end do
            if (maxval(terr) .lt. pl_rs_eps) then
               smaxp = s
               convrg = .true.
               exit
            end if
         end if
      end do
      if (.not. convrg) pl_error_codes(6) = 1
      if (k0z .eq. 0.d0) then
         do s = 1, smaxp
            kz = -2.d0 * pi * dble(s) / w(2)
            ct = kz / ri
            call crotcoef(ct, 0, nodr, ymn)
            do n = 0, nodr
               c = -((0.d0, -1.d0)**n) / w(2) / ri / sqrt(4.d0 * pi / dble(n + n + 1))
               do m = -n, n
                  l = m + n * (n + 1)
                  swfyzsum(l) = swfyzsum(l) &
                                + c * ymn(l) * sum1(m, s) * exp((0.d0, 1.d0) * kz * z)
               end do
            end do
         end do
      else
         convrg = .false.
         do s = 1, ntz
            kz = k0z - 2.d0 * pi * dble(s) / w(2)
            call q2db(nodr, x, y, w(1), k0y, kz, ri, q2d)
            sum1(:, s) = q2d
            ct = kz / ri
            call crotcoef(ct, 0, nodr, ymn)
            terr = 0.d0
            derr = abs(swfyzsum(:))
            do n = 0, nodr
               c = -((0.d0, -1.d0)**n) / w(2) / ri / sqrt(4.d0 * pi / dble(n + n + 1))
               do m = -n, n
                  l = m + n * (n + 1)
                  swfyzsum(l) = swfyzsum(l) &
                                + c * ymn(l) * sum1(m, s) * exp((0.d0, 1.d0) * kz * z)
               end do
            end do
            terr = abs(swfyzsum(:))
            if (abs(kz) .gt. 1.d0) then
               do n = 0, nodr * (nodr + 2)
                  if (terr(n) .ne. 0.d0) then
                     terr(n) = abs(terr(n) - derr(n)) / terr(n)
                  end if
               end do
               if (maxval(terr) .lt. pl_rs_eps) then
                  convrg = .true.
                  exit
               end if
            end if
         end do
         if (.not. convrg) pl_error_codes(6) = 1
      end if
      deallocate (sum1)
   end subroutine swfyzlatticesum

   subroutine swfyzlatticesum0(nodr, x, y, z, w, k0y, k0z, ri, swfyzsum)
      implicit none
      logical :: convrg
      integer :: nodr, m, n, s, ntz, l, smaxp, smaxn
      real(8) :: x, k0y, k0z, w(2), y, z, kz, mag, tmag
      complex(8) :: q1d(0:nodr * (nodr + 2)), q2d(-nodr:nodr), &
                    swfyzsum(0:nodr * (nodr + 2)), ri, ci, &
                    c, ct, ymn(0:nodr * (nodr + 2))
      complex(8), allocatable :: sum1(:, :)
      data ci/(0.d0, 1.d0)/
      ntz = max(20, ceiling(w(2) / 2.d0 / pi))
      allocate (sum1(-nodr:nodr, -ntz:ntz))
      smaxp = ntz
      smaxn = -ntz
      swfyzsum = 0.d0
      sum1 = 0.d0
      convrg = .false.
      do s = 0, ntz
         kz = k0z + 2.d0 * pi * dble(s) / w(2)
         call q2db(nodr, x, y, w(1), k0y, kz, ri, q2d)
         sum1(:, s) = q2d
         mag = sum(abs(q2d))
         if ((abs(kz) .gt. 1.d0) .and. (mag .lt. 1.d-8)) then
            smaxp = s
            convrg = .true.
            exit
         end if
      end do
      if (.not. convrg) pl_error_codes(6) = 1
      if (k0z .eq. 0.) then
         smaxn = -smaxp
         do s = 1, smaxp
            sum1(:, -s) = sum1(:, s)
         end do
      else
         convrg = .false.
         do s = 1, ntz
            kz = k0z - 2.d0 * pi * dble(s) / w(2)
            call q2db(nodr, x, y, w(1), k0y, kz, ri, q2d)
            sum1(:, -s) = q2d
            mag = sum(abs(q2d))
            if ((abs(kz) .gt. 1.d0) .and. (mag .lt. 1.d-8)) then
               smaxn = -s
               convrg = .true.
               exit
            end if
         end do
         if (.not. convrg) pl_error_codes(6) = 1
      end if
      call q1dbnosource(nodr, x, y, z, w(2), k0z, ri, q1d)
      do n = 0, nodr
         do m = -n, n
            swfyzsum(m + n * (n + 1)) = q1d(m + n * (n + 1))
         end do
      end do
      do s = smaxn, smaxp
         kz = k0z + 2.d0 * pi * dble(s) / w(2)
         ct = kz / ri
         call crotcoef(ct, 0, nodr, ymn)
         do n = 0, nodr
            c = -((0.d0, -1.d0)**n) / w(2) / ri / sqrt(4.d0 * pi / dble(n + n + 1))
            do m = -n, n
               l = m + n * (n + 1)
               swfyzsum(l) = swfyzsum(l) &
                             + c * ymn(l) * sum1(m, s) * exp((0.d0, 1.d0) * kz * z)
            end do
         end do
      end do
   end subroutine swfyzlatticesum0

   subroutine qkernel2d(ntot, t, qfunc)
      implicit none
      integer :: ntot, m
      real(8) :: t, tt, dt
      complex(8) :: qfunc(ntot), ci, u, efunc1, efunc2, pfunc1, pfunc2, c, du, v, rkz
      data ci/(0.d0, 1.d0)/
      if (qkernel_integration_model .eq. 0) then
         tt = t
         dt = 1.d0
      else
         tt = 1.d0 / t
         dt = tt / t
      end if
      rkz = sqrt((1.d0, 0.d0) - qkernel_kz * qkernel_kz / qkernel_ref_index**2)
      c = tt * tt - 2.d0 * ci * rkz
      u = tt * sqrt(c)
      du = (c + tt * tt) / sqrt(c)
      v = sqrt(1.d0 - qkernel_kz * qkernel_kz / qkernel_ref_index**2 - u * u)
      pfunc1 = (u - ci * v) / rkz
      pfunc2 = (u + ci * v) / rkz
      efunc1 = exp(ci * (qkernel_k0y * qkernel_width + qkernel_ref_index * v * &
                         (qkernel_width - qkernel_y))) / (exp(ci * (qkernel_k0y / qkernel_ref_index + v) * &
                                                              qkernel_width * qkernel_ref_index) - 1.d0)
      efunc2 = -exp(ci * v * qkernel_ref_index * (qkernel_width + qkernel_y)) / &
               (exp(ci * qkernel_k0y * qkernel_width) - exp(ci * v * qkernel_ref_index * qkernel_width))
      qfunc = 0.d0
      do m = -qkernel_nodr, qkernel_nodr
         qfunc(m + 1 + qkernel_nodr) = efunc1 * pfunc1**m + efunc2 * pfunc2**m
      end do
      qfunc = qfunc * exp(ci * u * qkernel_x * qkernel_ref_index) * du / v * dt
   end subroutine qkernel2d

   subroutine q2db(nodr, x, y, w, k0y, kz, ri, qint)
      implicit none
      integer :: nodr, ntot, subdiv, nseg, ec
      real(8) :: x, y, w, k0y, kz, t0, t1, cerr, dseg
      complex(8) :: qint(-nodr:nodr), qintt(-nodr:nodr), qintp(-nodr:nodr), ri
      qkernel_nodr = nodr
      qkernel_x = x
      qkernel_y = y
      qkernel_width = w
      qkernel_k0y = k0y
      qkernel_kz = kz
      qkernel_ref_index = ri
      qint = 0.
      ntot = 2 * nodr + 1
      qkernel_integration_model = 0
      cerr = 1.d0
      t0 = 0.d0
      nseg = 0
      dseg = lattice_integration_segment / qkernel_width
      do while (cerr .gt. pl_integration_limit_epsilon)
         t1 = t0
         t0 = t0 - dseg
         nseg = nseg + 1
         qintt = 0.
         subdiv = 0
         if (nseg .eq. 2 .and. pl_integration_method .eq. 0) then
            qkernel_integration_model = 1
            t0 = -1.d0 / dseg
            t1 = 0.d0
         end if
         ec = 0
         call gkintegrate(ntot, t0, t1, qkernel2d, qintt, subdiv, ec, &
                          pl_integration_error_epsilon, minimum_integration_spacing, maximum_integration_subdivisions)
         if (ec .ne. 0) pl_error_codes(4) = 1
         cerr = sum(abs(qintt))
         qint = qint + qintt
         cerr = cerr / sum(abs(qint))
         pl_max_subdivs = max(pl_max_subdivs, subdiv)
         if (nseg .eq. 2 .and. pl_integration_method .eq. 0) exit
      end do
      q2d_number_segments = nseg
      qkernel_integration_model = 0
      qintp = 0.d0
      cerr = 1.d0
      t1 = 0.d0
      nseg = 0
      do while (cerr .gt. pl_integration_limit_epsilon)
         t0 = t1
         t1 = t1 + dseg
         nseg = nseg + 1
         qintt = 0.
         subdiv = 0
         if (nseg .eq. 2 .and. pl_integration_method .eq. 0) then
            qkernel_integration_model = 1
            t1 = 1.d0 / dseg
            t0 = 0.d0
         end if
         ec = 0
         call gkintegrate(ntot, t0, t1, qkernel2d, qintt, subdiv, ec, &
                          pl_integration_error_epsilon, minimum_integration_spacing, maximum_integration_subdivisions)
         if (ec .ne. 0) pl_error_codes(4) = 1
         cerr = sum(abs(qintt))
         qintp = qintp + qintt
         cerr = cerr / sum(abs(qintp))
         pl_max_subdivs = max(pl_max_subdivs, subdiv)
         if (nseg .eq. 2 .and. pl_integration_method .eq. 0) exit
      end do
      q2d_number_segments = max(q2d_number_segments, nseg)
      qint = qint + qintp
   end subroutine q2db

   subroutine qkernel1d(ntot, t, qfunc)
      implicit none
      integer :: ntot, i, n, m, nodrtemp, mlim
      real(8) :: t, tt, dt, rho
      complex(8) :: qfunc(ntot), ymn1(0:qkernel_nodr * (qkernel_nodr + 2)), ci, u, st1, rhos, &
                    ephi, bfunc(-qkernel_nodr:qkernel_nodr), expzfunc1, expzfunc2, c
      data ci/(0.d0, 1.d0)/
      if (qkernel_integration_model .eq. 0) then
         tt = t
         dt = 1.d0
      else
         tt = 1.d0 / t
         dt = tt / t
      end if
      u = cmplx(1.d0, tt, kind=kind(0.0d0))
      rho = sqrt(qkernel_x * qkernel_x + qkernel_y * qkernel_y)
      st1 = sqrt((1.d0 - u) * (1.d0 + u))
      rhos = rho * st1 * qkernel_ref_index
      bfunc = 0.d0
      if (rho .eq. 0.) then
         ephi = 1.d0
         bfunc(0) = 1.d0
         nodrtemp = 0
      else
         ephi = cmplx(qkernel_x, qkernel_y, kind=kind(0.0d0)) / rho
         nodrtemp = qkernel_nodr
         call bessel_integer_complex(qkernel_nodr, rhos, nodrtemp, bfunc(0:qkernel_nodr))
         do m = 1, nodrtemp
            bfunc(-m) = (-1)**m * bfunc(m)
         end do
      end if
      expzfunc1 = exp(ci * u * qkernel_ref_index * (qkernel_width + qkernel_z)) / &
                  (exp(ci * qkernel_kz * qkernel_width) - exp(ci * qkernel_ref_index * u * qkernel_width))
      expzfunc2 = exp(ci * u * qkernel_ref_index * (qkernel_width - qkernel_z)) / &
                  (exp(-ci * qkernel_kz * qkernel_width) - exp(ci * qkernel_ref_index * u * qkernel_width))
      call crotcoef(u, 0, qkernel_nodr, ymn1)
      qfunc = 0.
      do n = 0, qkernel_nodr
         mlim = min(n, nodrtemp)
         do m = -mlim, mlim
            i = n * (n + 1) + m + 1
            c = -ci * ((-ci)**(n - m)) * (ephi**m) / 2.d0 / sqrt(pi / dble(n + n + 1))
            qfunc(i) = dt * c * ymn1(n * (n + 1) + m) * bfunc(m) * (expzfunc1 + ((-1)**(n + m)) * expzfunc2)
         end do
      end do
   end subroutine qkernel1d

   subroutine q1dbnosource(nodr, x, y, z, w, kz, ri, qint)
      implicit none
      integer :: nodr, ntot, subdiv, nseg, ec
      real(8) :: x, y, z, w, kz, t0, t1, cerr, dseg
      complex(8) :: qint(0:nodr * (nodr + 2)), qintt(0:nodr * (nodr + 2)), ri
      qkernel_nodr = nodr
      qkernel_x = x
      qkernel_y = y
      qkernel_z = z
      qkernel_width = w
      qkernel_kz = kz
      qkernel_ref_index = ri
      ntot = 1 + nodr * (nodr + 2)
      qkernel_integration_model = 0
      qint = 0.
      t1 = 0.d0
      cerr = 1.d0
      nseg = 0
      dseg = lattice_integration_segment / qkernel_width
      do while (cerr .gt. pl_integration_limit_epsilon)
         t0 = t1
         t1 = t1 + dseg
         nseg = nseg + 1
         subdiv = 0
         qintt = 0.
         if (nseg .eq. 2 .and. pl_integration_method .eq. 0) then
            qkernel_integration_model = 1
            t1 = 1.d0 / dseg
            t0 = 0.d0
         end if
         ec = 0
         call gkintegrate(ntot, t0, t1, qkernel1d, qintt, subdiv, ec, &
                          pl_integration_error_epsilon, minimum_integration_spacing, maximum_integration_subdivisions)
         if (ec .ne. 0) pl_error_codes(5) = 1
         cerr = sum(abs(qintt))
         qint = qint + qintt
         cerr = cerr / sum(abs(qint))
         pl_max_subdivs = max(pl_max_subdivs, subdiv)
         if (nseg .eq. 2 .and. pl_integration_method .eq. 0) exit
      end do
      q1d_number_segments = nseg
   end subroutine q1dbnosource

   subroutine reciprocal_space_swf_lattice_sum(nodr, x0, y0, z0, w, k0, ri, nmax, eps, nterms, wf)
      implicit none
      integer :: nodr, i, ix, iy, nmax, nterms, n, p, q
      real(8) :: x0, y0, z0, w(2), wx, wy, k0(2), kconst, kx, ky, eps, cerr(0:nodr * (nodr + 2)), &
                 asum(0:nodr * (nodr + 2)), asum0(0:nodr * (nodr + 2)), a0mag
      complex(8) :: ri, wf(0:nodr * (nodr + 2)), kswf(0:nodr * (nodr + 2)), dwf(0:nodr * (nodr + 2))
      wx = w(1)
      wy = w(2)
      kconst = 2.d0 * pi / wx / wy
      kx = k0(1)
      ky = k0(2)
      call reciprocal_scalar_wave_function(nodr, kx, ky, x0, y0, z0, ri, kswf)
      wf = kswf * kconst
      do n = 1, nmax
         dwf = 0.d0
         do i = 0, 8 * n - 1
            q = i / (2 * n)
            p = i - 2 * q * n
            if (q .eq. 0) then
               ix = n
               iy = -n + p
            elseif (q .eq. 1) then
               ix = n - p
               iy = n
            elseif (q .eq. 2) then
               ix = -n
               iy = n - p
            else
               ix = -n + p
               iy = -n
            end if
            kx = 2.d0 * pi * dble(ix) / wx + k0(1)
            ky = 2.d0 * pi * dble(iy) / wy + k0(2)
            call reciprocal_scalar_wave_function(nodr, kx, ky, x0, y0, z0, ri, kswf)
            dwf = dwf + kswf * kconst
         end do
         wf = wf + dwf
         asum0 = abs(wf)
         asum = abs(dwf)
         a0mag = maxval(asum0)
         cerr = 0.d0
         do i = 0, nodr * (nodr + 2)
            if (asum0(i) .gt. 1.d-15 * a0mag) cerr(i) = asum(i) / asum0(i)
         end do
         if (maxval(cerr) .lt. eps) exit
      end do
      nterms = n
   end subroutine reciprocal_space_swf_lattice_sum

   subroutine common_layer_lattice_kernel(nodr, kx, ky, x, y, zt, zs, tdirs, sdirs, kernel)
      implicit none
      integer :: nodr, n, m, k, mn, sdirs(2), tdirs(2), sdir, tdir, slay, i, is, id, it
      real(8) :: kx, ky, x, y, zs, zt, kr
      complex(8) :: s, gfunc(2, 2, 2), skz, tkz, &
                    drot(-2:2, 0:nodr * (nodr + 2)), ealpha, ri, c, c2, &
                    kernel(-1:1, 0:nodr * (nodr + 2), tdirs(1):tdirs(2), sdirs(1):sdirs(2))
      if (time_it) time_0 = mstm_mpi_wtime()
      slay = layer_id(zs)
      ri = layer_ref_index(slay)
      kr = sqrt(kx * kx + ky * ky)
      if (kr .eq. 0.d0) then
         ealpha = 1.d0
      else
         ealpha = cmplx(kx, ky, kind=kind(0.0d0)) / kr
      end if
      s = kr
      call layer_gf(s, zs, zt, gfunc, skz, tkz)
      call crotcoef(skz, 2, nodr, drot)
      c = exp((0.d0, 1.d0) * (kx * x + ky * y)) / ri / ri / skz / sqrt(4.d0 * pi)
      do k = -2, 2, 2
         i = k / 2
         is = (-1)**i
         do n = abs(k), nodr
            do m = -n, n
               mn = n * (n + 1) + m
               c2 = c * ealpha**m * drot(k, mn)
               do sdir = sdirs(1), sdirs(2)
                  do tdir = tdirs(1), tdirs(2)
                     if (sdir .eq. tdir) then
                        id = 1
                        it = -1
                     else
                        id = -1
                        it = 1
                     end if
           kernel(i, mn, tdir, sdir) = kernel(i, mn, tdir, sdir) + c2 * it * (gfunc(tdir, sdir, 1) + is * id * gfunc(tdir, sdir, 2))
                  end do
               end do
            end do
         end do
      end do
      if (time_it) time_count(4) = mstm_mpi_wtime() - time_0 + time_count(4)
   end subroutine common_layer_lattice_kernel

   subroutine plane_boundary_lattice_kernel(nodrt, nodrs, kx, ky, x, y, zt, zs, kernel, &
                                            include_source)
      implicit none
      logical :: incsrc
      logical, optional :: include_source
      integer :: nodrs, nodrt, n, m, p, k, l, q, mn, kl, ssign, tsign, pol
      real(8) :: kx, ky, x, y, zs, zt, kr
      complex(8) :: kernel(2, nodrt * (nodrt + 2), 2, nodrs * (nodrs + 2)), s, gfunc(2, 2, 2), skz, tkz, &
                    pivec(2, nodrs * (nodrs + 2), 2), picvec(2, nodrt * (nodrt + 2), 2), ealpha, ri, ekm, csum(2, 2), c
      if (present(include_source)) then
         incsrc = include_source
      else
         incsrc = .false.
      end if
      if (zt .eq. zs) incsrc = .false.
      if (time_it) time_0 = mstm_mpi_wtime()
      ri = layer_ref_index(layer_id(zs))
      kr = sqrt(kx * kx + ky * ky)
      if (kr .eq. 0.d0) then
         ealpha = 1.d0
      else
         ealpha = cmplx(kx, ky, kind=kind(0.0d0)) / kr
      end if
      s = kr
      call layer_gf(s, zs, zt, gfunc, skz, tkz, incsrc)
      call complexpivec(skz, nodrs, pivec, 1)
      call complexpivec(tkz, nodrt, picvec, -1)
      c = exp((0.d0, 1.d0) * (kx * x + ky * y)) / ri / ri / skz
      do n = 1, nodrt
         do m = -n, n
            mn = n * (n + 1) + m
            do l = 1, nodrs
               do k = -l, l
                  kl = l * (l + 1) + k
                  ekm = ealpha**(k - m)
                  csum = 0.d0
                  do p = 1, 2
                     tsign = (-1)**(m + n + p - 1)
                     do q = 1, 2
                        ssign = (-1)**(k + l + q - 1)
                        do pol = 1, 2
                           csum(p, q) = csum(p, q) + picvec(p, mn, pol) * pivec(q, kl, pol) &
                                        * (gfunc(1, 1, pol) + (-1)**pol * (tsign * gfunc(2, 1, pol) + ssign * gfunc(1, 2, pol)) &
                                           + tsign * ssign * gfunc(2, 2, pol))
                        end do
                     end do
                  end do
                  kernel(:, mn, :, kl) = kernel(:, mn, :, kl) + csum(:, :) * c * ekm
               end do
            end do
         end do
      end do
      if (time_it) time_count(3) = mstm_mpi_wtime() - time_0 + time_count(3)
   end subroutine plane_boundary_lattice_kernel

end module periodic_lattice_subroutines
