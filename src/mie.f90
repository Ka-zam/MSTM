module mie

   use, intrinsic :: iso_fortran_env, only: real64
contains
!
!  calculation of the max order of sphere expansions and storage of mie coefficients
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!  april 2012: all spheres are assumed OA: l/r formulation
!  february 2013: tmatrix file option.
!
   subroutine calculate_mie_coefficients(qeps)
      use sphere_data
      implicit none
      integer :: i, nodrn, nsphere, ntermstot, nblktot, nterms, &
                 n1, n2, j, n
      real(real64) :: qext, qsca, qeps, qabs
      complex(real64) :: rihost(2)
      complex(real64), allocatable :: anp(:, :, :), cnp(:, :, :), unp(:, :, :), &
                                      vnp(:, :, :), dnp(:, :, :), anpinv(:, :, :)
      nsphere = number_spheres
      if (allocated(sphere_order)) deallocate (sphere_order, mie_offset, sphere_block, &
                                               qext_mie, qabs_mie, optically_active, sphere_offset)
      allocate (sphere_order(nsphere), mie_offset(nsphere + 1), sphere_block(nsphere), &
                qext_mie(nsphere), qabs_mie(nsphere), optically_active(nsphere), &
                sphere_offset(nsphere + 1))
      ntermstot = 0
      nblktot = 0
      max_mie_order = 0
      mean_qext_mie = 0.d0
      mean_qabs_mie = 0.d0
      area_mean_radius = 0.d0
!
!
! march 2013: calculates orders, and
!             forces host to have at least same order as constituents.
!
      do i = 1, number_spheres
         call exterior_refractive_index(i, rihost)
         call optically_active_mie_coefficients(sphere_radius(i), sphere_ref_index(:, i), &
                                                nodrn, qeps, qext, qsca, qabs, &
                                                ri_medium=rihost)
         sphere_order(i) = nodrn
      end do
      do i = 1, number_spheres
         j = host_sphere(i)
         do while (j .ne. 0)
            sphere_order(j) = max(sphere_order(j), sphere_order(i))
            j = host_sphere(j)
         end do
      end do
!
!  calculate the order limits and efficiencies
!
      n = 0
      any_optically_active = .false.
      do i = 1, number_spheres
         call exterior_refractive_index(i, rihost)
         if (sphere_ref_index(1, i) .eq. sphere_ref_index(2, i)) then
            optically_active(i) = .false.
         else
            optically_active(i) = .true.
            any_optically_active = .true.
         end if
         nodrn = sphere_order(i)
         sphere_block(i) = 2 * nodrn * (nodrn + 2)
         call optically_active_mie_coefficients(sphere_radius(i), sphere_ref_index(:, i), &
                                                nodrn, 0.d0, qext, qsca, qabs, &
                                                ri_medium=rihost)
         nterms = 4 * nodrn
         mie_offset(i) = ntermstot
         ntermstot = ntermstot + nterms
         qext_mie(i) = qext
         qabs_mie(i) = qabs
         if (host_sphere(i) .eq. 0) then
            mean_qext_mie = mean_qext_mie + qext * sphere_radius(i)**2
            mean_qabs_mie = mean_qabs_mie + qabs * sphere_radius(i)**2
            area_mean_radius = area_mean_radius + sphere_radius(i)**2
            n = n + 1
         end if
         sphere_offset(i) = nblktot
         nblktot = nblktot + sphere_block(i) * number_field_expansions(i)
         max_mie_order = max(max_mie_order, sphere_order(i))
      end do
      n = max(n, 1)
      mie_offset(nsphere + 1) = ntermstot
      sphere_offset(nsphere + 1) = nblktot
      number_eqns = nblktot
      area_mean_radius = sqrt(area_mean_radius / dble(n))
      mean_qext_mie = mean_qext_mie / area_mean_radius**2 / dble(n)
      mean_qabs_mie = mean_qabs_mie / area_mean_radius**2 / dble(n)
!
! calculate the mie coefficients, and store in memory
!
      if (allocated(an_mie)) deallocate (an_mie, cn_mie, un_mie, vn_mie, dn_mie, &
                                         an_inv_mie)
      allocate (an_mie(ntermstot), cn_mie(ntermstot), un_mie(ntermstot), &
                vn_mie(ntermstot), dn_mie(ntermstot), an_inv_mie(ntermstot))
      do i = 1, number_spheres
         nodrn = sphere_order(i)
         call exterior_refractive_index(i, rihost)
         allocate (anp(2, 2, nodrn), cnp(2, 2, nodrn), unp(2, 2, nodrn), &
                   vnp(2, 2, nodrn), dnp(2, 2, nodrn), anpinv(2, 2, nodrn))
         call optically_active_mie_coefficients(sphere_radius(i), sphere_ref_index(:, i), &
                                                nodrn, 0.d0, qext, qsca, qabs, &
                                                anp_mie=anp, cnp_mie=cnp, dnp_mie=dnp, &
                                                unp_mie=unp, vnp_mie=vnp, anp_inv_mie=anpinv, &
                                                ri_medium=rihost)
         nterms = 4 * nodrn
         n1 = mie_offset(i) + 1
         n2 = mie_offset(i) + nterms
         an_mie(n1:n2) = reshape(anp(1:2, 1:2, 1:nodrn), (/nterms/))
         cn_mie(n1:n2) = reshape(cnp(1:2, 1:2, 1:nodrn), (/nterms/))
         dn_mie(n1:n2) = reshape(dnp(1:2, 1:2, 1:nodrn), (/nterms/))
         un_mie(n1:n2) = reshape(unp(1:2, 1:2, 1:nodrn), (/nterms/))
         vn_mie(n1:n2) = reshape(vnp(1:2, 1:2, 1:nodrn), (/nterms/))
         an_inv_mie(n1:n2) = reshape(anpinv(1:2, 1:2, 1:nodrn), (/nterms/))
         deallocate (anp, cnp, unp, vnp, dnp, anpinv)
      end do
   end subroutine calculate_mie_coefficients

   subroutine exterior_refractive_index(i, rihost)
      use sphere_data
      use surface
      implicit none
      integer :: i
      complex(real64) :: rihost(2)
      if (plane_surface_present .and. (host_sphere(i) .eq. 0)) then
         rihost = layer_ref_index(sphere_layer(i))
      else
         rihost = sphere_ref_index(:, host_sphere(i))
      end if
   end subroutine exterior_refractive_index
!
! transformation between lr and te tm basis
!
   subroutine left_right_to_mode_matrix(at, am)
      implicit none
      complex(real64) :: a(2, 2), am(2, 2), at(2, 2)
      a = at
      am(1, 1) = (a(1, 1) + a(1, 2) + a(2, 1) + a(2, 2)) / 2.
      am(1, 2) = (a(1, 1) - a(1, 2) + a(2, 1) - a(2, 2)) / 2.
      am(2, 1) = (a(1, 1) + a(1, 2) - a(2, 1) - a(2, 2)) / 2.
      am(2, 2) = (a(1, 1) - a(1, 2) - a(2, 1) + a(2, 2)) / 2.
   end subroutine left_right_to_mode_matrix
!
! optically active lorenz/mie coefficients
! original 30 March 2011
! April 2012: generalized LR formulation, generalized mie coefficients
!
   subroutine optically_active_mie_coefficients(x, ri, nodr0, qeps, qext, qsca, qabs, anp_mie, dnp_mie, &
                                                unp_mie, vnp_mie, cnp_mie, ri_medium, anp_inv_mie, dnp_eff_mie, anp_eff_mie)
      use bessel_functions, only: riccati_bessel, riccati_hankel
      use wave_functions, only: invert_two_by_two_matrix
      implicit none
      integer :: nstop, n, i, p, q, nodr0, s, t, ss, st
      real(real64) :: x, qeps, qext, qsca, fn1, err, qextt, qabs
      real(real64), allocatable :: qext1(:)
      complex(real64) :: ci, ri(2), xri(2, 2), rii(2, 2), ribulk(2), psip(2, 2), &
                         xip(2, 2), gmatinv(2, 2), bmatinv(2, 2), gmat(2, 2), bmat(2, 2), &
                         amat(2, 2), dmat(2, 2), umat(2, 2), vmat(2, 2), amatinv(2, 2), &
                         psipn(2, 2), xipn(2, 2)
      complex(real64), optional :: anp_mie(2, 2, *), dnp_mie(2, 2, *), ri_medium(2), &
                                   unp_mie(2, 2, *), vnp_mie(2, 2, *), cnp_mie(2, 2, *), &
                                   anp_inv_mie(2, 2, *), dnp_eff_mie(2, 2, *), anp_eff_mie(2, 2, *)
      complex(real64), allocatable :: psi(:, :, :), xi(:, :, :)
      data ci/(0.d0, 1.d0)/

      if (present(ri_medium)) then
         rii(:, 1) = ri_medium
      else
         rii(:, 1) = (/(1.d0, 0.d0), (1.d0, 0.d0)/)
      end if
      rii(:, 2) = ri
      ribulk(:) = 2.d0 / (1.d0 / rii(1, :) + 1.d0 / rii(2, :))
      xri = x * rii
      if (qeps .gt. 0.) then
         nstop = nint(x + 4.*x**(1./3.)) + 5.
      elseif (qeps .lt. 0.) then
         nstop = ceiling(-qeps)
         nodr0 = nstop
      else
         nstop = nodr0
      end if
      allocate (psi(0:nstop + 1, 2, 2), xi(0:nstop + 1, 2, 2), qext1(nstop))

      do i = 1, 2
         do p = 1, 2
            call riccati_bessel(nstop + 1, xri(p, i), psi(0, p, i))
            call riccati_hankel(nstop + 1, xri(p, i), xi(0, p, i))
         end do
      end do
      qabs = 0.d0
      qsca = 0.d0
      qext = 0.d0
! 2/6/2019: scaled psip,xip
      do n = 1, nstop
         psip(:, :) = psi(n - 1, :, :) - dble(n) * psi(n, :, :) / xri(:, :)
         xip(:, :) = xi(n - 1, :, :) - dble(n) * xi(n, :, :) / xri(:, :)
         psipn(:, :) = psi(n - 1, :, :) / psi(n, :, :) - dble(n) / xri(:, :)
         xipn(:, :) = xi(n - 1, :, :) / xi(n, :, :) - dble(n) / xri(:, :)
         do s = 1, 2
            ss = (-1)**s
            do t = 1, 2
               st = (-1)**t
               gmat(s, t) = (ss * ribulk(1) + st * ribulk(2)) / (rii(s, 2) * rii(t, 1)) &
                            * (psipn(s, 2) * xi(n, t, 1) - ss * st * xip(t, 1))
               amat(s, t) = (ss * ribulk(1) + st * ribulk(2)) / (rii(s, 2) * rii(t, 1)) &
                            * (psipn(s, 2) * psi(n, t, 1) - ss * st * psip(t, 1))
!                  umat(s,t)=(ss*ribulk(2)+st*ribulk(2))/(rii(s,2)*rii(t,2)) &
!                    *(xip(s,2)*psi(n,t,2)-ss*st*xi(n,s,2)*psip(t,2))
               umat(s, t) = (ss * ribulk(2) + st * ribulk(2)) / (rii(s, 2) * rii(t, 2)) * ci / psi(n, s, 2)
               bmat(s, t) = (ss * ribulk(2) + st * ribulk(1)) / (rii(s, 1) * rii(t, 2)) &
                            * (xipn(s, 1) * psi(n, t, 2) - ss * st * psip(t, 2))
!                  dmat(s,t)=(ss*ribulk(1)+st*ribulk(1))/(rii(s,1)*rii(t,1)) &
!                    *(psip(s,1)*xi(n,t,1)-ss*st*psi(n,s,1)*xip(t,1))
               dmat(s, t) = -(ss * ribulk(1) + st * ribulk(1)) / (rii(s, 1) * rii(t, 1)) * ci / xi(n, s, 1)
               vmat(s, t) = (ss * ribulk(2) + st * ribulk(1)) / (rii(s, 1) * rii(t, 2)) &
                            * (xipn(s, 1) * xi(n, t, 2) - ss * st * xip(t, 2))
            end do
         end do

!write(*,'(8e12.4)') gmat
!flush(6)

         call invert_two_by_two_matrix(gmat, gmatinv)
         call invert_two_by_two_matrix(bmat, bmatinv)

         amat = -matmul(gmatinv, amat)
         umat = -matmul(gmatinv, umat)
         dmat = -matmul(bmatinv, dmat)
         vmat = -matmul(bmatinv, vmat)

         if (present(anp_mie)) then
            anp_mie(:, :, n) = amat(:, :)
         end if
         if (present(dnp_mie)) then
            dnp_mie(:, :, n) = dmat(:, :)
         end if
         if (present(unp_mie)) then
            unp_mie(:, :, n) = umat(:, :)
         end if
         if (present(vnp_mie)) then
            vnp_mie(:, :, n) = vmat(:, :)
         end if
         if (present(cnp_mie)) then
            call invert_two_by_two_matrix(amat, amatinv)
            amatinv = matmul(dmat, amatinv)
            cnp_mie(:, :, n) = amatinv(:, :)
         end if
         if (present(anp_inv_mie)) then
            call invert_two_by_two_matrix(amat, amatinv)
            anp_inv_mie(:, :, n) = amatinv(:, :)
         end if
         if (present(dnp_eff_mie)) then
            call invert_two_by_two_matrix(umat, amatinv)
            amatinv = matmul(amatinv, amat)
            amatinv = matmul(vmat, amatinv)
            dnp_eff_mie(:, :, n) = dmat - amatinv
         end if
         if (present(anp_eff_mie)) then
            call invert_two_by_two_matrix(umat, amatinv)
            amatinv = matmul(amatinv, amat)
            anp_eff_mie(:, :, n) = -amatinv
         end if
         qext1(n) = 0.d0
         fn1 = n + n + 1
         do p = 1, 2
            do q = 1, 2
               qsca = qsca + fn1 * abs(amat(p, q)) * abs(amat(p, q))
            end do
            qext1(n) = qext1(n) - fn1 * dble(amat(p, p))
         end do
         qext = qext + qext1(n)
      end do

      if (qeps .gt. 0.d0) then
         qextt = qext
         qext = 0.
         do n = 1, nstop
            qext = qext + qext1(n)
            err = abs(1.d0 - qext / qextt)
            if (err .lt. qeps .or. n .eq. nstop) exit
         end do
         nodr0 = n
      end if
      qsca = 2./x / x * qsca
      qext = 2./x / x * qext
      qabs = qext - qsca
      nstop = min(n, nstop)
      return
   end subroutine optically_active_mie_coefficients

!
!  multiplies coefficients for sphere i by appropriate lm coefficient.
!  lr, oa model
!  april 2012
!
   subroutine apply_single_sphere_mie_coefficients(i, nodr, cx, cy, mie_coefficient)
      use sphere_data
      implicit none
      integer :: i, n, p, nodr, n1, n2, nterms
      complex(real64) :: cx(0:nodr + 1, nodr, 2), cy(0:nodr + 1, nodr, 2)
      complex(real64), allocatable :: an1(:, :, :)
      character(len=1), optional :: mie_coefficient
      character(len=1) :: miecoefficient

      if (present(mie_coefficient)) then
         miecoefficient = mie_coefficient
      else
         miecoefficient = 'a'
      end if
      nterms = 4 * nodr
      n1 = mie_offset(i) + 1
      n2 = mie_offset(i) + nterms
      allocate (an1(2, 2, nodr))
      if (miecoefficient .eq. 'a') then
         an1 = reshape(an_mie(n1:n2), (/2, 2, nodr/))
      elseif (miecoefficient .eq. 'c') then
         an1 = reshape(cn_mie(n1:n2), (/2, 2, nodr/))
      elseif (miecoefficient .eq. 'd') then
         an1 = reshape(dn_mie(n1:n2), (/2, 2, nodr/))
      elseif (miecoefficient .eq. 'u') then
         an1 = reshape(un_mie(n1:n2), (/2, 2, nodr/))
      elseif (miecoefficient .eq. 'v') then
         an1 = reshape(vn_mie(n1:n2), (/2, 2, nodr/))
      elseif (miecoefficient .eq. 'i') then
         an1 = reshape(an_inv_mie(n1:n2), (/2, 2, nodr/))
      end if
      do n = 1, nodr
         do p = 1, 2
            cy(n + 1, n:1:-1, p) = an1(p, 1, n) * cx(n + 1, n:1:-1, 1) &
                                   + an1(p, 2, n) * cx(n + 1, n:1:-1, 2)
            cy(0:n, n, p) = an1(p, 1, n) * cx(0:n, n, 1) &
                            + an1(p, 2, n) * cx(0:n, n, 2)
         end do
      end do
      deallocate (an1)
   end subroutine apply_single_sphere_mie_coefficients
!
! generalized mie coefficient mult:
!  (a,f) = (generalized mie matrix)*(g,b)
! idir not = 1 does the transpose.
! aout is written over in this one.
! february 2013: tmatrix file option
!
   subroutine apply_mie_coefficients(neqns, nrhs, idir, ain, aout, rhs_list)
      use sphere_data
      implicit none
      integer :: neqns, i, n, p, q, nodri, nblki, n1, n2, b11, b12, b21, b22, idir, &
                 nterms, j, nrhs
      complex(real64) :: ain(neqns, nrhs), aout(neqns, nrhs)
      complex(real64), allocatable :: gin_t(:, :, :), aout_t(:, :, :), &
                                      an1(:, :, :), dn1(:, :, :), un1(:, :, :), vn1(:, :, :), &
                                      bin_t(:, :, :), fout_t(:, :, :)
      logical, optional :: rhs_list(nrhs)
      logical :: rhslist(nrhs)
      if (present(rhs_list)) then
         rhslist = rhs_list
      else
         rhslist = .true.
      end if

      do j = 1, nrhs
         do i = 1, number_spheres
            nodri = sphere_order(i)
            nblki = 2 * nodri * (nodri + 2)
            allocate (gin_t(0:nodri + 1, nodri, 2), aout_t(0:nodri + 1, nodri, 2), &
                      an1(2, 2, nodri))
            b11 = sphere_offset(i) + 1
            b12 = sphere_offset(i) + nblki
            gin_t = reshape(ain(b11:b12, j), (/nodri + 2, nodri, 2/))
            nterms = 4 * nodri
            n1 = mie_offset(i) + 1
            n2 = mie_offset(i) + nterms
            an1 = reshape(an_mie(n1:n2), (/2, 2, nodri/))
            if (number_field_expansions(i) .eq. 1) then
               aout_t = 0.d0
               if (rhslist(j) .and. idir .eq. 1) then
                  do n = 1, nodri
                     do p = 1, 2
                        do q = 1, 2
                           aout_t(n + 1, n:1:-1, p) = aout_t(n + 1, n:1:-1, p) &
                                                      + an1(p, q, n) * gin_t(n + 1, n:1:-1, q)
                           aout_t(0:n, n, p) = aout_t(0:n, n, p) &
                                               + an1(p, q, n) * gin_t(0:n, n, q)
                        end do
                     end do
                  end do
               elseif (rhslist(j) .and. idir .ne. 1) then
                  do n = 1, nodri
                     do p = 1, 2
                        do q = 1, 2
                           aout_t(n + 1, n:1:-1, p) = aout_t(n + 1, n:1:-1, p) &
                                                      + an1(q, p, n) * gin_t(n + 1, n:1:-1, q)
                           aout_t(0:n, n, p) = aout_t(0:n, n, p) &
                                               + an1(q, p, n) * gin_t(0:n, n, q)
                        end do
                     end do
                  end do
               end if
               aout(b11:b12, j) &
                  = reshape(aout_t(0:nodri + 1, 1:nodri, 1:2), (/nblki/))
            else
               allocate (bin_t(0:nodri + 1, nodri, 2), &
                         fout_t(0:nodri + 1, nodri, 2), &
                         dn1(2, 2, nodri), un1(2, 2, nodri), vn1(2, 2, nodri))
               b21 = sphere_offset(i) + nblki + 1
               b22 = sphere_offset(i) + 2 * nblki
               bin_t = reshape(ain(b21:b22, j), (/nodri + 2, nodri, 2/))
               dn1 = reshape(dn_mie(n1:n2), (/2, 2, nodri/))
               un1 = reshape(un_mie(n1:n2), (/2, 2, nodri/))
               vn1 = reshape(vn_mie(n1:n2), (/2, 2, nodri/))
               aout_t = 0.d0
               fout_t = 0.d0
               if (rhslist(j) .and. idir .eq. 1) then
                  do n = 1, nodri
                     do p = 1, 2
                        do q = 1, 2
                           aout_t(n + 1, n:1:-1, p) = aout_t(n + 1, n:1:-1, p) &
                                                      + an1(p, q, n) * gin_t(n + 1, n:1:-1, q) &
                                                      + un1(p, q, n) * bin_t(n + 1, n:1:-1, q)
                           aout_t(0:n, n, p) = aout_t(0:n, n, p) &
                                               + an1(p, q, n) * gin_t(0:n, n, q) &
                                               + un1(p, q, n) * bin_t(0:n, n, q)
                           fout_t(n + 1, n:1:-1, p) = fout_t(n + 1, n:1:-1, p) &
                                                      + dn1(p, q, n) * gin_t(n + 1, n:1:-1, q) &
                                                      + vn1(p, q, n) * bin_t(n + 1, n:1:-1, q)
                           fout_t(0:n, n, p) = fout_t(0:n, n, p) &
                                               + dn1(p, q, n) * gin_t(0:n, n, q) &
                                               + vn1(p, q, n) * bin_t(0:n, n, q)
                        end do
                     end do
                  end do
               elseif (rhslist(j) .and. idir .ne. 1) then
                  do n = 1, nodri
                     do p = 1, 2
                        do q = 1, 2
                           aout_t(n + 1, n:1:-1, p) = aout_t(n + 1, n:1:-1, p) &
                                                      + an1(q, p, n) * gin_t(n + 1, n:1:-1, q) &
                                                      + dn1(q, p, n) * bin_t(n + 1, n:1:-1, q)
                           aout_t(0:n, n, p) = aout_t(0:n, n, p) &
                                               + an1(q, p, n) * gin_t(0:n, n, q) &
                                               + dn1(q, p, n) * bin_t(0:n, n, q)
                           fout_t(n + 1, n:1:-1, p) = fout_t(n + 1, n:1:-1, p) &
                                                      + un1(q, p, n) * gin_t(n + 1, n:1:-1, q) &
                                                      + vn1(q, p, n) * bin_t(n + 1, n:1:-1, q)
                           fout_t(0:n, n, p) = fout_t(0:n, n, p) &
                                               + un1(q, p, n) * gin_t(0:n, n, q) &
                                               + vn1(q, p, n) * bin_t(0:n, n, q)

                        end do
                     end do
                  end do
               end if
               aout(b11:b12, j) &
                  = reshape(aout_t(0:nodri + 1, 1:nodri, 1:2), (/nblki/))
               aout(b21:b22, j) &
                  = reshape(fout_t(0:nodri + 1, 1:nodri, 1:2), (/nblki/))
               deallocate (bin_t, fout_t, un1, vn1, dn1)
            end if
            deallocate (gin_t, aout_t, an1)
         end do
      end do
   end subroutine apply_mie_coefficients

end module mie
