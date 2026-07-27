module gpfa_radix3
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
contains

   subroutine gpfa3f(a, b, trigs, inc, jump, n, mm, lot, isign)
      implicit none
      integer :: inc, jump, n, mm, lot, isign, lvr, inq, jstepx, ninc, ink, &
                 m, mh, nblox, left, nb, nvex, la, mu, ipass, jstep, jstepl, &
                 jjj, ja, nu, jb, jc, jd, j, l, kk, k, je, jf, jg, jh, laincl, ji, &
                 n3, istart, ll
      real(real64) :: a(*), b(*), trigs(*), s, t2, t1, t3, u2, u1, u3, co1, si1, &
                      co2, si2, c1, &
                      sin60
      data sin60/0.866025403784437/
      data lvr/64/
!
!     ***************************************************************
!     *                                                             *
!     *  n.b. lvr = length of vector registers, set to 128 for c90. *
!     *  reset to 64 for other cray machines, or to any large value *
!     *  (greater than or equal to lot) for a scalar computer.      *
!     *                                                             *
!     ***************************************************************
!
      n3 = 3**mm
      inq = n / n3
      jstepx = (n3 - n) * inc
      ninc = n * inc
      ink = inc * inq
      mu = mod(inq, 3)
      if (isign .eq. -1) mu = 3 - mu
      m = mm
      mh = (m + 1) / 2
      s = real(isign, kind=real64)
      c1 = sin60
      if (mu .eq. 2) c1 = -c1
!
      nblox = 1 + (lot - 1) / lvr
      left = lot
      s = real(isign, kind=real64)
      istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
      do nb = 1, nblox
!
         if (left .le. lvr) then
            nvex = left
         else if (left .lt. (2 * lvr)) then
            nvex = left / 2
            nvex = nvex + mod(nvex, 2)
         else
            nvex = lvr
         end if
         left = left - nvex
!
         la = 1
!
!  loop on type i radix-3 passes
!  -----------------------------
         do ipass = 1, mh
            jstep = (n * inc) / (3 * la)
            jstepl = jstep - ninc
!
!  k = 0 loop (no twiddle factors)
!  -------------------------------
            do jjj = 0, (n - 1) * inc, 3 * jstep
               ja = istart + jjj
!
!  "transverse" loop
!  -----------------
               do nu = 1, inq
                  jb = ja + jstepl
                  if (jb .lt. istart) jb = jb + ninc
                  jc = jb + jstepl
                  if (jc .lt. istart) jc = jc + ninc
                  j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
                  do l = 1, nvex
                     t1 = a(jb + j) + a(jc + j)
                     t2 = a(ja + j) - 0.5 * t1
                     t3 = c1 * (a(jb + j) - a(jc + j))
                     u1 = b(jb + j) + b(jc + j)
                     u2 = b(ja + j) - 0.5 * u1
                     u3 = c1 * (b(jb + j) - b(jc + j))
                     a(ja + j) = a(ja + j) + t1
                     b(ja + j) = b(ja + j) + u1
                     a(jb + j) = t2 - u3
                     b(jb + j) = u2 + t3
                     a(jc + j) = t2 + u3
                     b(jc + j) = u2 - t3
                     j = j + jump
                  end do
                  ja = ja + jstepx
                  if (ja .lt. istart) ja = ja + ninc
               end do
            end do
!
!  finished if n3 = 3
!  ------------------
            if (n3 .eq. 3) go to 490
            kk = 2 * la
!
!  loop on nonzero k
!  -----------------
            do k = ink, jstep - ink, ink
               co1 = trigs(kk + 1)
               si1 = s * trigs(kk + 2)
               co2 = trigs(2 * kk + 1)
               si2 = s * trigs(2 * kk + 2)
!
!  loop along transform
!  --------------------
               do jjj = k, (n - 1) * inc, 3 * jstep
                  ja = istart + jjj
!
!  "transverse" loop
!  -----------------
                  do nu = 1, inq
                     jb = ja + jstepl
                     if (jb .lt. istart) jb = jb + ninc
                     jc = jb + jstepl
                     if (jc .lt. istart) jc = jc + ninc
                     j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep,shortloop
                     do l = 1, nvex
                        t1 = a(jb + j) + a(jc + j)
                        t2 = a(ja + j) - 0.5 * t1
                        t3 = c1 * (a(jb + j) - a(jc + j))
                        u1 = b(jb + j) + b(jc + j)
                        u2 = b(ja + j) - 0.5 * u1
                        u3 = c1 * (b(jb + j) - b(jc + j))
                        a(ja + j) = a(ja + j) + t1
                        b(ja + j) = b(ja + j) + u1
                        a(jb + j) = co1 * (t2 - u3) - si1 * (u2 + t3)
                        b(jb + j) = si1 * (t2 - u3) + co1 * (u2 + t3)
                        a(jc + j) = co2 * (t2 + u3) - si2 * (u2 - t3)
                        b(jc + j) = si2 * (t2 + u3) + co2 * (u2 - t3)
                        j = j + jump
                     end do
!-----( end of loop across transforms )
                     ja = ja + jstepx
                     if (ja .lt. istart) ja = ja + ninc
                  end do
               end do
!-----( end of loop along transforms )
               kk = kk + 2 * la
            end do
!-----( end of loop on nonzero k )
            la = 3 * la
         end do
!-----( end of loop on type i radix-3 passes)
!
!  loop on type ii radix-3 passes
!  ------------------------------
!
         do ipass = mh + 1, m
            jstep = (n * inc) / (3 * la)
            jstepl = jstep - ninc
            laincl = la * ink - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
            do ll = 0, (la - 1) * ink, 3 * jstep
!
               do jjj = ll, (n - 1) * inc, 3 * la * ink
                  ja = istart + jjj
!
!  "transverse" loop
!  -----------------
                  do nu = 1, inq
                     jb = ja + jstepl
                     if (jb .lt. istart) jb = jb + ninc
                     jc = jb + jstepl
                     if (jc .lt. istart) jc = jc + ninc
                     jd = ja + laincl
                     if (jd .lt. istart) jd = jd + ninc
                     je = jd + jstepl
                     if (je .lt. istart) je = je + ninc
                     jf = je + jstepl
                     if (jf .lt. istart) jf = jf + ninc
                     jg = jd + laincl
                     if (jg .lt. istart) jg = jg + ninc
                     jh = jg + jstepl
                     if (jh .lt. istart) jh = jh + ninc
                     ji = jh + jstepl
                     if (ji .lt. istart) ji = ji + ninc
                     j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
                     do l = 1, nvex
                        t1 = a(jb + j) + a(jc + j)
                        t2 = a(ja + j) - 0.5 * t1
                        t3 = c1 * (a(jb + j) - a(jc + j))
                        a(jb + j) = a(jd + j)
                        u1 = b(jb + j) + b(jc + j)
                        u2 = b(ja + j) - 0.5 * u1
                        u3 = c1 * (b(jb + j) - b(jc + j))
                        b(jb + j) = b(jd + j)
                        a(ja + j) = a(ja + j) + t1
                        b(ja + j) = b(ja + j) + u1
                        a(jd + j) = t2 - u3
                        b(jd + j) = u2 + t3
                        a(jc + j) = t2 + u3
                        b(jc + j) = u2 - t3
!----------------------
                        t1 = a(je + j) + a(jf + j)
                        t2 = a(jb + j) - 0.5 * t1
                        t3 = c1 * (a(je + j) - a(jf + j))
                        a(jf + j) = a(jh + j)
                        u1 = b(je + j) + b(jf + j)
                        u2 = b(jb + j) - 0.5 * u1
                        u3 = c1 * (b(je + j) - b(jf + j))
                        b(jf + j) = b(jh + j)
                        a(jb + j) = a(jb + j) + t1
                        b(jb + j) = b(jb + j) + u1
                        a(je + j) = t2 - u3
                        b(je + j) = u2 + t3
                        a(jh + j) = t2 + u3
                        b(jh + j) = u2 - t3
!----------------------
                        t1 = a(jf + j) + a(ji + j)
                        t2 = a(jg + j) - 0.5 * t1
                        t3 = c1 * (a(jf + j) - a(ji + j))
                        t1 = a(jg + j) + t1
                        a(jg + j) = a(jc + j)
                        u1 = b(jf + j) + b(ji + j)
                        u2 = b(jg + j) - 0.5 * u1
                        u3 = c1 * (b(jf + j) - b(ji + j))
                        u1 = b(jg + j) + u1
                        b(jg + j) = b(jc + j)
                        a(jc + j) = t1
                        b(jc + j) = u1
                        a(jf + j) = t2 - u3
                        b(jf + j) = u2 + t3
                        a(ji + j) = t2 + u3
                        b(ji + j) = u2 - t3
                        j = j + jump
                     end do
!-----( end of loop across transforms )
                     ja = ja + jstepx
                     if (ja .lt. istart) ja = ja + ninc
                  end do
               end do
            end do
!-----( end of double loop for k=0 )
!
!  finished if last pass
!  ---------------------
            if (ipass .eq. m) go to 490
!
            kk = 2 * la
!
!     loop on nonzero k
!     -----------------
            do k = ink, jstep - ink, ink
               co1 = trigs(kk + 1)
               si1 = s * trigs(kk + 2)
               co2 = trigs(2 * kk + 1)
               si2 = s * trigs(2 * kk + 2)
!
!  double loop along first transform in block
!  ------------------------------------------
               do ll = k, (la - 1) * ink, 3 * jstep
!
                  do jjj = ll, (n - 1) * inc, 3 * la * ink
                     ja = istart + jjj
!
!  "transverse" loop
!  -----------------
                     do nu = 1, inq
                        jb = ja + jstepl
                        if (jb .lt. istart) jb = jb + ninc
                        jc = jb + jstepl
                        if (jc .lt. istart) jc = jc + ninc
                        jd = ja + laincl
                        if (jd .lt. istart) jd = jd + ninc
                        je = jd + jstepl
                        if (je .lt. istart) je = je + ninc
                        jf = je + jstepl
                        if (jf .lt. istart) jf = jf + ninc
                        jg = jd + laincl
                        if (jg .lt. istart) jg = jg + ninc
                        jh = jg + jstepl
                        if (jh .lt. istart) jh = jh + ninc
                        ji = jh + jstepl
                        if (ji .lt. istart) ji = ji + ninc
                        j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
                        do l = 1, nvex
                           t1 = a(jb + j) + a(jc + j)
                           t2 = a(ja + j) - 0.5 * t1
                           t3 = c1 * (a(jb + j) - a(jc + j))
                           a(jb + j) = a(jd + j)
                           u1 = b(jb + j) + b(jc + j)
                           u2 = b(ja + j) - 0.5 * u1
                           u3 = c1 * (b(jb + j) - b(jc + j))
                           b(jb + j) = b(jd + j)
                           a(ja + j) = a(ja + j) + t1
                           b(ja + j) = b(ja + j) + u1
                           a(jd + j) = co1 * (t2 - u3) - si1 * (u2 + t3)
                           b(jd + j) = si1 * (t2 - u3) + co1 * (u2 + t3)
                           a(jc + j) = co2 * (t2 + u3) - si2 * (u2 - t3)
                           b(jc + j) = si2 * (t2 + u3) + co2 * (u2 - t3)
!----------------------
                           t1 = a(je + j) + a(jf + j)
                           t2 = a(jb + j) - 0.5 * t1
                           t3 = c1 * (a(je + j) - a(jf + j))
                           a(jf + j) = a(jh + j)
                           u1 = b(je + j) + b(jf + j)
                           u2 = b(jb + j) - 0.5 * u1
                           u3 = c1 * (b(je + j) - b(jf + j))
                           b(jf + j) = b(jh + j)
                           a(jb + j) = a(jb + j) + t1
                           b(jb + j) = b(jb + j) + u1
                           a(je + j) = co1 * (t2 - u3) - si1 * (u2 + t3)
                           b(je + j) = si1 * (t2 - u3) + co1 * (u2 + t3)
                           a(jh + j) = co2 * (t2 + u3) - si2 * (u2 - t3)
                           b(jh + j) = si2 * (t2 + u3) + co2 * (u2 - t3)
!----------------------
                           t1 = a(jf + j) + a(ji + j)
                           t2 = a(jg + j) - 0.5 * t1
                           t3 = c1 * (a(jf + j) - a(ji + j))
                           t1 = a(jg + j) + t1
                           a(jg + j) = a(jc + j)
                           u1 = b(jf + j) + b(ji + j)
                           u2 = b(jg + j) - 0.5 * u1
                           u3 = c1 * (b(jf + j) - b(ji + j))
                           u1 = b(jg + j) + u1
                           b(jg + j) = b(jc + j)
                           a(jc + j) = t1
                           b(jc + j) = u1
                           a(jf + j) = co1 * (t2 - u3) - si1 * (u2 + t3)
                           b(jf + j) = si1 * (t2 - u3) + co1 * (u2 + t3)
                           a(ji + j) = co2 * (t2 + u3) - si2 * (u2 - t3)
                           b(ji + j) = si2 * (t2 + u3) + co2 * (u2 - t3)
                           j = j + jump
                        end do
!-----(end of loop across transforms)
                        ja = ja + jstepx
                        if (ja .lt. istart) ja = ja + ninc
                     end do
                  end do
               end do
!-----( end of double loop for this k )
               kk = kk + 2 * la
            end do
!-----( end of loop over values of k )
            la = 3 * la
         end do
!-----( end of loop on type ii radix-3 passes )
!-----( nvex transforms completed)
490      continue
         istart = istart + nvex * jump
      end do
!-----( end of loop on blocks of transforms )
!
      return
   end subroutine gpfa3f

end module gpfa_radix3
