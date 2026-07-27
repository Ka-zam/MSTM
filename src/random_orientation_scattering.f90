module random_orientation_scattering
   use, intrinsic :: iso_fortran_env, only: real32, real64
   use constants
   use fft_translation
   use mie
   use parallel_runtime, only: mpi_comm_world, parallel_allreduce_sum, parallel_barrier, parallel_broadcast, &
                               parallel_communicator_create, parallel_group, parallel_group_include, parallel_rank, &
                               parallel_size, parallel_wall_time
   use numerical_tables
   use periodic_lattice_operations
   use runtime_support, only: open_input_file, runtime_failed, synchronize_runtime_status, write_elapsed_time
   use angular_functions, only: rotation_coefficients, vector_coupling_coefficients
   use sphere_data
   use surface
   implicit none
   private
   public :: evaluate_random_orientation_scattering_matrix, random_orientation_scattering_matrix
contains

!
!  compute the coefficients for the GSF expansion of the random orientation
!  scattering matrix.
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!  January 2012: added computation of coherent field average
!  February 2013:  added number processors option.
!
   subroutine random_orientation_scattering_matrix(tmatrixfile, sm, smcf, beam_width, number_processors, mpi_comm, &
                                                   mean_t_matrix, override_order, keep_quiet)
      implicit none
      logical :: symmetrical, nkq
      logical, optional :: keep_quiet
      integer :: file_unit, nodr, nodrw, nodr2, m, n, p, k, l, q, t, v, u, w, nblk, kl, mn, nn1, tn, &
                 lmax, ll1, tvl, ku, k1, ns, ik, ik1, m1, nu, n1s, n1e, nu1, p1, n1max, &
                 in, n1, i, kt, nodrt, nodrrhs, mnm, klm, ikm, &
                 rank, numprocs, &
                 nblkrhs, numprocscalc, orig_group, new_group, &
                 new_comm, new_rank, nblkw, wv, sizedm, sizetm, nread, mpicomm, rank0
      integer, allocatable :: windex(:), vindex(:), wvindex(:), wvnum(:), group_list(:)
      integer, optional :: number_processors, mpi_comm, override_order
      real(real64) :: sm(4, 4, 0:*), fl, xv, fl2, cbeam, gbn, wvperproc, sum, &
                      time1, time2, smcf(4, 4, 0:*), qextt, qscatt
      real(real64), allocatable :: vc(:)
      real(real64), optional :: beam_width
      complex(real64) :: ci, cin, a, tct, tcp(2)
      complex(real64), allocatable :: aw(:, :, :), bw(:, :, :), cw(:), &
                                      dw(:), pp(:, :, :), bm(:, :, :), &
                                      am(:, :, :), fm(:, :, :, :, :), bmcf(:, :, :), &
                                      amcf(:, :, :), fmcf(:, :, :, :, :), awcf(:, :, :), &
                                      bwcf(:, :, :), cwcf(:), dwcf(:)
      complex(real64), allocatable :: dm(:, :, :, :, :, :), dmcf(:, :, :, :, :, :)
      complex(real64), optional :: mean_t_matrix(*)
      complex(real32), allocatable :: tc(:, :, :, :)
      character(len=30) :: tmatrixfile
      data ci/(0.d0, 1.d0)/
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      if (present(keep_quiet)) then
         nkq = .not. keep_quiet
      else
         nkq = .true.
      end if
      call parallel_rank(mpi_rank=rank, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=rank0)
      call parallel_size(mpi_size=numprocs, mpi_comm=mpicomm)
      if (present(beam_width)) then
         cbeam = beam_width
      else
         cbeam = 0.d0
      end if
      if (present(number_processors)) then
         numprocscalc = min(number_processors, numprocs)
      else
         numprocscalc = numprocs
      end if
      allocate (group_list(numprocscalc))
      group_list = (/(i, i=0, numprocscalc - 1)/)
      call parallel_group(mpi_group=orig_group, mpi_comm=mpicomm)
      call parallel_group_include(mpi_group=orig_group, &
                                  mpi_size=numprocscalc, mpi_new_group_list=group_list, &
                                  mpi_new_group=new_group, mpi_comm=mpicomm)
      call parallel_communicator_create(mpi_group=new_group, &
                                        mpi_new_comm=new_comm, mpi_comm=mpicomm)

      xv = sphere_cluster%cross_section_radius
      if (rank .le. numprocscalc - 1) then

         call parallel_rank(mpi_comm=new_comm, &
                            mpi_rank=new_rank)
         if (rank .eq. 0) time1 = parallel_wall_time()
!
!  read the T matrix from the file
!
         do i = 0, numprocscalc - 1
            if (i .eq. new_rank) then
               call open_input_file(tmatrixfile, file_unit)
               if (.not. runtime_failed()) then
                  read (file_unit, *) nodr, nodrrhs
                  if (present(override_order)) then
                     nodr = override_order
                     nodrrhs = nodr
                  end if
                  symmetrical = nodr .eq. nodrrhs
                  nodrt = nodr
                  nodr2 = nodr + nodr
                  nodrw = nodr2
                  nblk = nodr * (nodr + 2)
!                  nodrrhs=nodr
                  nblkrhs = nodrrhs * (nodrrhs + 2)
                  sizetm = 4 * nblk * nblkrhs
                  allocate (tc(2, nblk, 2, nblkrhs))
                  tc = (0., 0.)
                  do l = 1, nodrrhs
                     if (present(mean_t_matrix)) then
                        mean_t_matrix(2 * l - 1:2 * l) = 0.d0
                     end if
                     do k = -l, l
                        kl = l * (l + 1) + k
                        klm = l * (l + 1) - k
                        do q = 1, 2
                           if (symmetrical) then
                              nread = l
                           else
                              nread = nodr
                           end if
                           do n = 1, nread
                              do m = -n, n
                                 mn = n * (n + 1) + m
                                 mnm = n * (n + 1) - m
                                 do p = 1, 2
                                    read (file_unit, *) tcp(p)
                                    tc(p, mn, q, kl) = tcp(p)
                                 end do
                                 if (present(mean_t_matrix) .and. n .eq. l .and. m .eq. k) then
                                    mean_t_matrix(2 * (l - 1) + q) = mean_t_matrix(2 * (l - 1) + q) + tcp(q)
                                 end if
                                 if (n .lt. l .and. symmetrical) then
                                    ikm = (-1)**(m + k)
                                    do p = 1, 2
                                       tc(q, klm, p, mnm) = tc(p, mn, q, kl) * ikm
                                    end do
                                 end if
                              end do
                           end do
                        end do
                     end do
                     if (present(mean_t_matrix)) then
                        mean_t_matrix(2 * l - 1:2 * l) = mean_t_matrix(2 * l - 1:2 * l) / dble(2 * l + 1)
                     end if
                  end do
!if(rank.eq.0) then
!open(30,file='meantmatrix.dat')
!do n=1,nodrt
!do p=1,2
!tct=0.
!do m=-n,n
!mn=n*(n+1)+m
!tct=tct+tc(p,mn,p,mn)
!enddo
!tct=tct/dble(n+n+1)
!write(30,'(2i5,2es13.5)') n,p,tct
!enddo
!enddo
!close(30)
!endif
                  qextt = 0.d0
                  qscatt = 0.d0
                  do l = 1, nodrrhs
                     do k = -l, l
                        kl = l * (l + 1) + k
                        do q = 1, 2
                           do n = 1, nodr
                              do m = -n, n
                                 mn = n * (n + 1) + m
                                 do p = 1, 2
                                    qscatt = qscatt + cabs(tc(p, mn, q, kl))**2
                                 end do
                              end do
                           end do
                           if (l .le. nodr) qextt = qextt - real(tc(q, kl, q, kl))
                        end do
                     end do
                  end do
                  qextt = qextt * 2./sphere_cluster%cross_section_radius**2
                  qscatt = qscatt * 2./sphere_cluster%cross_section_radius**2
                  close (file_unit)
               end if
            end if
            call parallel_barrier(mpi_comm=new_comm)
            call synchronize_runtime_status(new_comm)
            if (runtime_failed()) return
         end do
         if (rank0 .eq. 0 .and. nkq) then
            write (sphere_cluster%run_print_unit, '('' t matrix ext, sca:'',2e13.5)') qextt, qscatt
         end if
!!
!!  send to the other processors
!!
!            if(numprocscalc.gt.1) then
!               call parallel_broadcast(send_buffer=tc, &
!                    mpi_number=sizetm,mpi_rank=0,mpi_comm=new_comm)
!            endif

         allocate (vc(0:4 * nodr + 2))
         allocate (aw(0:2, -1:1, 0:nodrw), bw(0:2, -1:1, 0:nodrw), cw(0:nodrw), &
                   dw(0:nodrw), pp(nodr, 2, 2), bm(2, nodr * (nodr + 2), 2), &
                   am(2, nodrrhs + 1, 2), fm(3, nodr, 2, nodr, 2), bmcf(2, nodr * (nodr + 2), 2), &
                   amcf(2, nodrrhs + 1, 2), fmcf(3, nodr, 2, nodr, 2), awcf(0:2, -1:1, 0:nodrw), &
                   bwcf(0:2, -1:1, 0:nodrw), cwcf(0:nodrw), dwcf(0:nodrw))
         allocate (dm(-nodr - 1:nodr + 1, 3, nodr, 2, nodr, 2), dmcf(-nodr - 1:nodr + 1, 3, nodr, 2, nodr, 2))
         if (rank0 .eq. 0 .and. nkq) then
            time2 = parallel_wall_time() - time1
            call write_elapsed_time(sphere_cluster%run_print_unit, ' t matrix read time:', time2)
            time1 = parallel_wall_time()
         end if
         dm = (0.d0, 0.d0)
         dmcf = (0.d0, 0.d0)
         sizedm = size(dm)
         call initialize_numerical_tables(nodr2)
!
!  compute the GB modified T matrix
!
         do n = 1, nodrrhs
            gbn = exp(-((dble(n) + .5d0) * cbeam)**2.)
            cin = ci**(n + 1)
            pp(n, 1, 1) = -.5d0 * cin * fnr(n + n + 1) * gbn
            pp(n, 2, 1) = -pp(n, 1, 1)
            pp(n, 1, 2) = -pp(n, 1, 1)
            pp(n, 2, 2) = pp(n, 2, 1)
         end do
         do n = 1, nodr
            nn1 = n * (n + 1)
            do m = -n, n
               mn = nn1 + m
               do p = 1, 2
                  do l = 1, nodrrhs
                     do k = -l, l
                        kl = l * (l + 1) + k
                        a = tc(p, mn, 1, kl)
                        tc(p, mn, 1, kl) = tc(p, mn, 1, kl) * pp(l, 1, 1) &
                                           + tc(p, mn, 2, kl) * pp(l, 1, 2)
                        tc(p, mn, 2, kl) = a * pp(l, 2, 1) + tc(p, mn, 2, kl) * pp(l, 2, 2)
                     end do
                  end do
               end do
            end do
         end do
!
!  determine the distribution of work load among the processors
!
         nblkw = nodr2 * (nodr2 + 2) + 1
         allocate (windex(nblkw), vindex(nblkw), wvindex(0:numprocscalc - 1), &
                   wvnum(0:numprocscalc - 1))
         w = 0
         do n = 0, nodr2
            do m = -n, n
               w = w + 1
               windex(w) = n
               vindex(w) = m
            end do
         end do
         wvperproc = dble(nblkw) / dble(numprocscalc)
         sum = 0.
         do i = 0, numprocscalc - 1
            wvindex(i) = floor(sum)
            sum = sum + wvperproc
         end do
         do i = 0, numprocscalc - 2
            wvnum(i) = wvindex(i + 1) - wvindex(i)
         end do
         wvnum(numprocscalc - 1) = nblkw - wvindex(numprocscalc - 1)
         if (rank0 .eq. 0 .and. nkq) then
            write (sphere_cluster%run_print_unit, '('' d matrix calculation, order+degree per proc.:'',f9.2)') &
               wvperproc
            flush (sphere_cluster%run_print_unit)
         end if
!
!  the big loop
!
         do wv = wvindex(rank) + 1, wvindex(rank) + wvnum(rank)
            w = windex(wv)
            v = vindex(wv)
            bm = (0.d0, 0.d0)
            bmcf = (0.d0, 0.d0)
            do n = 1, nodr
               nn1 = n * (n + 1)
               do l = max(1, abs(w - n)), min(nodrrhs, w + n)
                  am(1, l, 1) = 0.d0
                  am(1, l, 2) = 0.d0
                  am(2, l, 1) = 0.d0
                  am(2, l, 2) = 0.d0
                  amcf(1, l, 1) = 0.d0
                  amcf(1, l, 2) = 0.d0
                  amcf(2, l, 1) = 0.d0
                  amcf(2, l, 2) = 0.d0
               end do
               do t = -n, n
                  tn = nn1 + t
                  lmax = min(nodrrhs, w + n)
                  call vector_coupling_coefficients(v, w, -t, n, vc)
                  do l = max(1, abs(v - t), abs(n - w)), lmax
                     ll1 = l * (l + 1)
                     tvl = ll1 + t - v
                     do k = 1, 2
                        do p = 1, 2
                           am(k, l, p) = am(k, l, p) + vc(l) * tc(p, tn, k, tvl)
                           if (l .eq. n .and. v .eq. 0) then
! may 2019: orientation averaged T matrix is azimuthal independent
                              tct = 0.d0
                              do kt = -n, n
                                 tct = tct + tc(p, nn1 + kt, k, nn1 + kt)
                              end do
                              amcf(k, l, p) = amcf(k, l, p) + vc(l) * tct / dble(n + n + 1)
!                                 amcf(k,l,p)=amcf(k,l,p)+vc(l)*tcm(p,n)
!                                 amcf(k,l,p)=amcf(k,l,p)+vc(l)*tc(p,tn,k,tvl)
                           end if
                        end do
                     end do
                  end do
               end do
               do m = -n, n
                  mn = nn1 + m
                  do k = 1, 2
                     u = m - (-3 + 2 * k)
                     if (abs(u) .le. w) then
                        lmax = min(nodrrhs, w + n)
                        call vector_coupling_coefficients(-u, w, m, n, vc)
                        do l = max(1, abs(w - n)), lmax
                           fl = -(-1)**l * vc(l) / dble(l + l + 1)
                           do p = 1, 2
                              bm(k, mn, p) = bm(k, mn, p) + am(k, l, p) * fl
                              if (v .eq. 0) then
                                 bmcf(k, mn, p) = bmcf(k, mn, p) + amcf(k, l, p) * fl
                              end if
                           end do
                        end do
                     end if
                  end do
               end do
            end do

!               do k=1,2
!                  do n=1,nodr*(nodr+2)
!                     do p=1,2
!                        bmcf(k,n,p)=bm(k,n,p)-bmcf(k,n,p)
!                     enddo
!                  enddo
!               enddo

            do u = -min(w, nodr + 1), min(w, nodr + 1)
               do ku = 1, 3
                  if (ku .eq. 1) then
                     k = -1
                     k1 = -1
                  elseif (ku .eq. 2) then
                     k = 1
                     k1 = 1
                  else
                     k = 1
                     k1 = -1
                  end if
                  m = u + k
                  ns = max(1, abs(m))
                  ik = (k + 1) / 2 + 1
                  ik1 = (k1 + 1) / 2 + 1
                  m1 = u + k1
                  do n = ns, nodr
                     nu = n * (n + 1) + m
                     n1s = max(1, abs(m1), n - nodrw)
                     n1e = min(nodr, n + nodrw)
                     do n1 = n1s, n1e
                        cin = ci**(n - n1)
                        nu1 = n1 * (n1 + 1) + m1
                        fl = -fnr(n + n + 1) * fnr(n1 + n1 + 1) * dble(w + w + 1)
                        do p = 1, 2
                           do p1 = 1, 2
                              a = bm(ik, nu, p) * cin * fl * conjg(bm(ik1, nu1, p1))
                              dm(u, ku, n, p, n1, p1) = dm(u, ku, n, p, n1, p1) + a
                              if (v .eq. 0) then
                                 a = bmcf(ik, nu, p) * cin * fl * conjg(bmcf(ik1, nu1, p1))
                                 dmcf(u, ku, n, p, n1, p1) = dmcf(u, ku, n, p, n1, p1) + a
                              end if
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
         deallocate (tc)

         call parallel_allreduce_sum(receive_buffer=dm, &
                                     mpi_number=sizedm, mpi_comm=new_comm)
         call parallel_allreduce_sum(receive_buffer=dmcf, &
                                     mpi_number=sizedm, mpi_comm=new_comm)
         if (rank0 .eq. 0 .and. nkq) then
            time2 = parallel_wall_time() - time1
            call write_elapsed_time(sphere_cluster%run_print_unit, ' d matrix time:', time2)
            time1 = parallel_wall_time()
         end if
!
!  compute the expansion coefficients
!
         aw = 0.d0
         bw = 0.d0
         cw = 0.d0
         dw = 0.d0
         awcf = 0.d0
         bwcf = 0.d0
         cwcf = 0.d0
         dwcf = 0.d0
         do w = 0, nodrw
            do n = 1, nodr
               n1s = max(1, abs(n - w))
               n1e = min(nodr, n + w)
               do n1 = n1s, n1e
                  do k = 1, 3
                     do p = 1, 2
                        do p1 = 1, 2
                           fm(k, n, p, n1, p1) = 0.
                           fmcf(k, n, p, n1, p1) = 0.
                        end do
                     end do
                  end do
               end do
            end do
            do u = -nodr - 1, nodr + 1
               do k = -1, 1, 2
                  m = u + k
                  ik = (k + 1) / 2 + 1
                  ns = max(1, abs(m))
                  do n = ns, nodr
                     n1max = min(w + n, nodr)
                     call vector_coupling_coefficients(m, n, 0, w, vc)
                     do n1 = ns, nodr
                        if ((n + n1 .lt. w) .or. (abs(n - n1) .gt. w)) cycle
                        fl = -(-1)**n * vc(n1) * fnr(w + w + 1) / fnr(n1 + n1 + 1)
                        do p = 1, 2
                           do p1 = 1, 2
                              fm(ik, n, p, n1, p1) = fm(ik, n, p, n1, p1) &
                                                     + dm(u, ik, n, p, n1, p1) * fl
                              fmcf(ik, n, p, n1, p1) = fmcf(ik, n, p, n1, p1) &
                                                       + dmcf(u, ik, n, p, n1, p1) * fl
                           end do
                        end do
                     end do
                  end do
               end do
               if (w .lt. 2) cycle
               m = u + 1
               m1 = u - 1
               ns = max(1, abs(m))
               n1s = max(1, abs(m1))
               do n = ns, nodr
                  n1max = min(w + n, nodr)
                  call vector_coupling_coefficients(m, n, -2, w, vc)
                  do n1 = n1s, nodr
                     if ((n + n1 .lt. w) .or. (abs(n - n1) .gt. w)) cycle
                     fl = -(-1)**n * vc(n1) * fnr(w + w + 1) / fnr(n1 + n1 + 1)
                     do p = 1, 2
                        do p1 = 1, 2
                           fm(3, n, p, n1, p1) = fm(3, n, p, n1, p1) &
                                                 + dm(u, 3, n, p, n1, p1) * fl
                           fmcf(3, n, p, n1, p1) = fmcf(3, n, p, n1, p1) &
                                                   + dmcf(u, 3, n, p, n1, p1) * fl
                        end do
                     end do
                  end do
               end do
            end do
            do n = 1, nodr
               n1s = max(1, abs(n - w))
               n1e = min(nodr, n + w)
               in = (-1)**n
               n1max = min(w + n, nodr)
               call vector_coupling_coefficients(1, n, 0, w, vc)
               do n1 = n1s, n1e
                  fl = 2.d0 * in * vc(n1) * fnr(w + w + 1) / fnr(n1 + n1 + 1)
                  i = mod(n + n1 - w, 2) + 1
                  do p = 1, 2
                     p1 = (2 - i) * p + (i - 1) * (3 - p)
                     do k = -1, 1, 2
                        ik = (k + 1) / 2 + 1
                        aw(0, k, w) = aw(0, k, w) + fm(ik, n, p, n1, p1) * fl
                        bw(0, k, w) = bw(0, k, w) + fm(ik, n, p, n1, 3 - p1) * fl
                        awcf(0, k, w) = awcf(0, k, w) + fmcf(ik, n, p, n1, p1) * fl
                        bwcf(0, k, w) = bwcf(0, k, w) + fmcf(ik, n, p, n1, 3 - p1) * fl
                     end do
                     bw(2, 0, w) = bw(2, 0, w) + fm(3, n, p, n1, 3 - p1) * fl
                     aw(2, 0, w) = aw(2, 0, w) + fm(3, n, p, n1, p1) * fl
                     bwcf(2, 0, w) = bwcf(2, 0, w) + fmcf(3, n, p, n1, 3 - p1) * fl
                     awcf(2, 0, w) = awcf(2, 0, w) + fmcf(3, n, p, n1, p1) * fl
                  end do
               end do
               if (w .lt. 2) cycle
               call vector_coupling_coefficients(1, n, -2, w, vc)
               do n1 = n1s, n1e
                  fl = 2.d0 * in * vc(n1) * fnr(w + w + 1) / fnr(n1 + n1 + 1)
                  i = mod(n + n1 - w, 2) + 1
                  do p = 1, 2
                     p1 = (2 - i) * p + (i - 1) * (3 - p)
                     do k = -1, 1, 2
                        ik = (k + 1) / 2 + 1
                        aw(2, k, w) = aw(2, k, w) + fm(ik, n, p, n1, p1) * fl * (-1)**p1
                        bw(2, k, w) = bw(2, k, w) + fm(ik, n, p, n1, 3 - p1) * fl * (-1)**(3 - p1)
                        awcf(2, k, w) = awcf(2, k, w) + fmcf(ik, n, p, n1, p1) * fl * (-1)**p1
                        bwcf(2, k, w) = bwcf(2, k, w) &
                                        + fmcf(ik, n, p, n1, 3 - p1) * fl * (-1)**(3 - p1)
                     end do
                  end do
                  fl2 = 2.*(-1)**(n1 + w) * vc(n1) * fnr(w + w + 1) / fnr(n1 + n1 + 1)
                  do p = 1, 2
                     do p1 = 1, 2
                        cw(w) = cw(w) + fm(3, n, p, n1, p1) * fl * (-1)**p1
                        dw(w) = dw(w) + fm(3, n, p, n1, p1) * fl2 * (-1)**p
                        cwcf(w) = cwcf(w) + fmcf(3, n, p, n1, p1) * fl * (-1)**p1
                        dwcf(w) = dwcf(w) + fmcf(3, n, p, n1, p1) * fl2 * (-1)**p
                     end do
                  end do
               end do
            end do
         end do
!            do w=0,nodrw
!               do k=-1,1
!                  do i=0,2
!                     aw(i,k,w)=aw(i,k,w)*2./xv/xv
!                     bw(i,k,w)=bw(i,k,w)*2./xv/xv
!                     awcf(i,k,w)=awcf(i,k,w)*2./xv/xv
!                     bwcf(i,k,w)=bwcf(i,k,w)*2./xv/xv
!                  enddo
!               enddo
!               cw(w)=cw(w)*2./xv/xv
!               dw(w)=dw(w)*2./xv/xv
!               cwcf(w)=cwcf(w)*2./xv/xv
!               dwcf(w)=dwcf(w)*2./xv/xv
!            enddo
         do n = 0, nodrw
            sm(1, 1, n) = aw(0, -1, n) + aw(0, 1, n)
            sm(1, 2, n) = aw(2, -1, n) + aw(2, 1, n)
            sm(1, 3, n) = 2.d0 * aimag(aw(2, 0, n))
            sm(1, 4, n) = aw(0, 1, n) - aw(0, -1, n)
            sm(2, 2, n) = dw(n)
            sm(2, 3, n) = aimag(dw(n))
            sm(2, 4, n) = aw(2, 1, n) - aw(2, -1, n)
            sm(3, 1, n) = -aimag(bw(2, -1, n) + bw(2, 1, n))
            sm(3, 2, n) = aimag(cw(n))
            sm(3, 3, n) = cw(n)
            sm(3, 4, n) = aimag(bw(2, -1, n) - bw(2, 1, n))
            sm(4, 1, n) = bw(0, -1, n) + bw(0, 1, n)
            sm(4, 2, n) = 2.*bw(2, 0, n)
            sm(4, 4, n) = bw(0, 1, n) - bw(0, -1, n)

            smcf(1, 1, n) = awcf(0, -1, n) + awcf(0, 1, n)
            smcf(1, 2, n) = awcf(2, -1, n) + awcf(2, 1, n)
            smcf(1, 3, n) = 2.d0 * aimag(awcf(2, 0, n))
            smcf(1, 4, n) = awcf(0, 1, n) - awcf(0, -1, n)
            smcf(2, 2, n) = dwcf(n)
            smcf(2, 3, n) = aimag(dwcf(n))
            smcf(2, 4, n) = awcf(2, 1, n) - awcf(2, -1, n)
            smcf(3, 1, n) = -aimag(bwcf(2, -1, n) + bwcf(2, 1, n))
            smcf(3, 2, n) = aimag(cwcf(n))
            smcf(3, 3, n) = cwcf(n)
            smcf(3, 4, n) = aimag(bwcf(2, -1, n) - bwcf(2, 1, n))
            smcf(4, 1, n) = bwcf(0, -1, n) + bwcf(0, 1, n)
            smcf(4, 2, n) = 2.*bwcf(2, 0, n)
            smcf(4, 4, n) = bwcf(0, 1, n) - bwcf(0, -1, n)

         end do
! patch 10-22
         sm(:, :, 0:nodrw) = sm(:, :, 0:nodrw) / dble(layer_ref_index(0)) / 2.d0
         smcf(:, :, 0:nodrw) = smcf(:, :, 0:nodrw) / dble(layer_ref_index(0)) / 2.d0
!            sm(:,:,0:nodrw)=sm(:,:,0:nodrw)*four_pi
!            smcf(:,:,0:nodrw)=smcf(:,:,0:nodrw)*four_pi
!
!  normalization
!
!            qsca0=sm(1,1,0)
!            do n=0,nodrw
!               do i=1,4
!                  do j=1,4
!                     sm(i,j,n)=sm(i,j,n)/qsca0
!                     smcf(i,j,n)=smcf(i,j,n)/qsca0
!                  enddo
!               enddo
!            enddo
         if (rank0 .eq. 0 .and. nkq) then
            time2 = parallel_wall_time() - time1
            call write_elapsed_time(sphere_cluster%run_print_unit, ' scat matrix coef time:', time2)
         end if
         deallocate (windex, vindex, wvindex, wvnum, dm)
         deallocate (vc)
         deallocate (aw, bw, cw, dw, pp, bm, am, fm, bmcf, amcf, fmcf, awcf, &
                     bwcf, cwcf, dwcf)
      end if
      call parallel_barrier(mpi_comm=mpicomm)
   end subroutine random_orientation_scattering_matrix
!
!  calculation of the RO scattering matrix from the GSF expansion
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!
   subroutine evaluate_random_orientation_scattering_matrix(ct, smc, nodrexp, sm)
      implicit none
      integer :: nodrexp, n, nn0, nnp2, nnm2
      real(real64) :: smc(4, 4, 0:nodrexp), sm(4, 4), dc(-2:2, 0:nodrexp * (nodrexp + 2)), &
                      ct
!
!     dc is the normalized generalized spherical function
!     dc(k,n*(n+1)+m) = ((n-k)!(n+m)!/(n+k)!/(n-m)!)^(1/2) D^k_{mn},
!     where D^k_{mn} is defined in M&M JOSA 96
!
      call rotation_coefficients(ct, 2, nodrexp, dc)
      sm = 0.d0
      do n = 0, nodrexp
         nn0 = n * (n + 1)
         nnp2 = nn0 + 2
         nnm2 = nn0 - 2
         sm(1, 1) = sm(1, 1) + dc(0, nn0) * smc(1, 1, n)
         sm(1, 4) = sm(1, 4) + dc(0, nn0) * smc(1, 4, n)
         sm(4, 1) = sm(4, 1) + dc(0, nn0) * smc(4, 1, n)
         sm(4, 4) = sm(4, 4) + dc(0, nn0) * smc(4, 4, n)
         if (n .ge. 2) then
            sm(1, 2) = sm(1, 2) + dc(2, nn0) * smc(1, 2, n)
            sm(2, 4) = sm(2, 4) + dc(2, nn0) * smc(2, 4, n)
            sm(3, 4) = sm(3, 4) + dc(2, nn0) * smc(3, 4, n)
            sm(1, 3) = sm(1, 3) + dc(2, nn0) * smc(1, 3, n)
            sm(3, 1) = sm(1, 3) + dc(2, nn0) * smc(3, 1, n)
            sm(4, 2) = sm(4, 2) + dc(2, nn0) * smc(4, 2, n)
            sm(2, 2) = sm(2, 2) + dc(2, nnm2) * smc(2, 2, n) + dc(2, nnp2) * smc(3, 3, n)
            sm(2, 3) = sm(2, 3) + dc(2, nnp2) * smc(2, 3, n) + dc(2, nnp2) * smc(3, 2, n)
            sm(3, 3) = sm(3, 3) - dc(2, nnm2) * smc(2, 2, n) + dc(2, nnp2) * smc(3, 3, n)
            sm(3, 2) = sm(3, 2) + dc(2, nnp2) * smc(2, 3, n) - dc(2, nnp2) * smc(3, 2, n)
         end if
      end do
      sm(2, 1) = sm(1, 2)
      sm(4, 3) = -sm(3, 4)
!
!  discontiued scaling option: now done in main program
!
!            if(iscale.eq.1) then
!               do j=1,4
!                  do k=j,4
!                     if(j.ne.1.or.k.ne.1) then
!                        sm(j,k,i)=sm(j,k,i)/sm(1,1,i)
!                     endif
!                  enddo
!               enddo
!            endif
!
!    here are the VV and HH differential cross sections
!
!            gvv=.25*(sm(1,1)+sm(2,2)-2.*sm(1,2))
!            ghh=.25*(sm(1,1)+sm(2,2)+2.*sm(1,2))
!
      return
   end subroutine evaluate_random_orientation_scattering_matrix
end module random_orientation_scattering
