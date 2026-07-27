module gpfa_radix2
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
contains

   subroutine gpfa2f(a, b, trigs, inc, jump, n, mm, lot, isign)
      implicit none
      integer :: inc, jump, n, mm, lot, isign, lvr, n2, inq, jstepx, ninc, ink, &
                 m2, m8, m, mh, nblox, left, nb, nvex, la, mu, ipass, jstep, jstepl, &
                 jjj, ja, nu, jb, jc, jd, j, l, kk, k, je, jf, jg, jh, laincl, ji, jj, jk, &
                 jl, jm, jn, jo, jp, istart, ll
      real(real64) :: a(*), b(*), trigs(*), s, ss, t0, t2, t1, t3, u0, u2, u1, u3, co1, si1, &
                      co2, si2, co3, si3, c1, c2, c3, co4, si4, co5, si5, co6, si6, co7, si7
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
      n2 = 2**mm
      inq = n / n2
      jstepx = (n2 - n) * inc
      ninc = n * inc
      ink = inc * inq
!
      m2 = 0
      m8 = 0
      m = 0
      if (mod(mm, 2) .eq. 0) then
         m = mm / 2
      else if (mod(mm, 4) .eq. 1) then
         m = (mm - 1) / 2
         m2 = 1
      else if (mod(mm, 4) .eq. 3) then
         m = (mm - 3) / 2
         m8 = 1
      end if
      mh = (m + 1) / 2
!
      nblox = 1 + (lot - 1) / lvr
      left = lot
      s = dble(isign)
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
!  loop on type i radix-4 passes
!  -----------------------------
         mu = mod(inq, 4)
         if (isign .eq. -1) mu = 4 - mu
         ss = 1.0
         if (mu .eq. 3) ss = -1.0
!
         if (mh .eq. 0) go to 200
!
         do ipass = 1, mh
            jstep = (n * inc) / (4 * la)
            jstepl = jstep - ninc
!
!  k = 0 loop (no twiddle factors)
!  -------------------------------
            do jjj = 0, (n - 1) * inc, 4 * jstep
               ja = istart + jjj
!
!     "transverse" loop
!     -----------------
               do nu = 1, inq
                  jb = ja + jstepl
                  if (jb .lt. istart) jb = jb + ninc
                  jc = jb + jstepl
                  if (jc .lt. istart) jc = jc + ninc
                  jd = jc + jstepl
                  if (jd .lt. istart) jd = jd + ninc
                  j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
                  do l = 1, nvex
                     t0 = a(ja + j) + a(jc + j)
                     t2 = a(ja + j) - a(jc + j)
                     t1 = a(jb + j) + a(jd + j)
                     t3 = ss * (a(jb + j) - a(jd + j))
                     u0 = b(ja + j) + b(jc + j)
                     u2 = b(ja + j) - b(jc + j)
                     u1 = b(jb + j) + b(jd + j)
                     u3 = ss * (b(jb + j) - b(jd + j))
                     a(ja + j) = t0 + t1
                     a(jc + j) = t0 - t1
                     b(ja + j) = u0 + u1
                     b(jc + j) = u0 - u1
                     a(jb + j) = t2 - u3
                     a(jd + j) = t2 + u3
                     b(jb + j) = u2 + t3
                     b(jd + j) = u2 - t3
                     j = j + jump
                  end do
                  ja = ja + jstepx
                  if (ja .lt. istart) ja = ja + ninc
               end do
            end do
!
!  finished if n2 = 4
!  ------------------

            if (n2 .eq. 4) go to 490
            kk = 2 * la
!
!  loop on nonzero k
!  -----------------
            do k = ink, jstep - ink, ink
               co1 = trigs(kk + 1)
               si1 = s * trigs(kk + 2)
               co2 = trigs(2 * kk + 1)
               si2 = s * trigs(2 * kk + 2)
               co3 = trigs(3 * kk + 1)
               si3 = s * trigs(3 * kk + 2)
!
!  loop along transform
!  --------------------
               do jjj = k, (n - 1) * inc, 4 * jstep
                  ja = istart + jjj
!
!     "transverse" loop
!     -----------------
                  do nu = 1, inq
                     jb = ja + jstepl
                     if (jb .lt. istart) jb = jb + ninc
                     jc = jb + jstepl
                     if (jc .lt. istart) jc = jc + ninc
                     jd = jc + jstepl
                     if (jd .lt. istart) jd = jd + ninc
                     j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep,shortloop
                     do l = 1, nvex
                        t0 = a(ja + j) + a(jc + j)
                        t2 = a(ja + j) - a(jc + j)
                        t1 = a(jb + j) + a(jd + j)
                        t3 = ss * (a(jb + j) - a(jd + j))
                        u0 = b(ja + j) + b(jc + j)
                        u2 = b(ja + j) - b(jc + j)
                        u1 = b(jb + j) + b(jd + j)
                        u3 = ss * (b(jb + j) - b(jd + j))
                        a(ja + j) = t0 + t1
                        b(ja + j) = u0 + u1
                        a(jb + j) = co1 * (t2 - u3) - si1 * (u2 + t3)
                        b(jb + j) = si1 * (t2 - u3) + co1 * (u2 + t3)
                        a(jc + j) = co2 * (t0 - t1) - si2 * (u0 - u1)
                        b(jc + j) = si2 * (t0 - t1) + co2 * (u0 - u1)
                        a(jd + j) = co3 * (t2 + u3) - si3 * (u2 - t3)
                        b(jd + j) = si3 * (t2 + u3) + co3 * (u2 - t3)
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
            la = 4 * la
         end do
!-----( end of loop on type i radix-4 passes)
!
!  central radix-2 pass
!  --------------------
200      continue
         if (m2 .eq. 0) go to 300
!
         jstep = (n * inc) / (2 * la)
         jstepl = jstep - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
         do jjj = 0, (n - 1) * inc, 2 * jstep
            ja = istart + jjj
!
!     "transverse" loop
!     -----------------
            do nu = 1, inq
               jb = ja + jstepl
               if (jb .lt. istart) jb = jb + ninc
               j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
               do l = 1, nvex
                  t0 = a(ja + j) - a(jb + j)
                  a(ja + j) = a(ja + j) + a(jb + j)
                  a(jb + j) = t0
                  u0 = b(ja + j) - b(jb + j)
                  b(ja + j) = b(ja + j) + b(jb + j)
                  b(jb + j) = u0
                  j = j + jump
               end do
!-----(end of loop across transforms)
               ja = ja + jstepx
               if (ja .lt. istart) ja = ja + ninc
            end do
         end do
!
!  finished if n2=2
!  ----------------
         if (n2 .eq. 2) go to 490
!
         kk = 2 * la
!
!  loop on nonzero k
!  -----------------
         do k = ink, jstep - ink, ink
            co1 = trigs(kk + 1)
            si1 = s * trigs(kk + 2)
!
!  loop along transforms
!  ---------------------
            do jjj = k, (n - 1) * inc, 2 * jstep
               ja = istart + jjj
!
!     "transverse" loop
!     -----------------
               do nu = 1, inq
                  jb = ja + jstepl
                  if (jb .lt. istart) jb = jb + ninc
                  j = 0
!
!  loop across transforms
!  ----------------------
                  if (kk .eq. n2 / 2) then
!cdir$ ivdep, shortloop
                     do l = 1, nvex
                        t0 = ss * (a(ja + j) - a(jb + j))
                        a(ja + j) = a(ja + j) + a(jb + j)
                        a(jb + j) = ss * (b(jb + j) - b(ja + j))
                        b(ja + j) = b(ja + j) + b(jb + j)
                        b(jb + j) = t0
                        j = j + jump
                     end do
!
                  else
!
!cdir$ ivdep, shortloop
                     do l = 1, nvex
                        t0 = a(ja + j) - a(jb + j)
                        a(ja + j) = a(ja + j) + a(jb + j)
                        u0 = b(ja + j) - b(jb + j)
                        b(ja + j) = b(ja + j) + b(jb + j)
                        a(jb + j) = co1 * t0 - si1 * u0
                        b(jb + j) = si1 * t0 + co1 * u0
                        j = j + jump
                     end do
!
                  end if
!
!-----(end of loop across transforms)
                  ja = ja + jstepx
                  if (ja .lt. istart) ja = ja + ninc
               end do
            end do
!-----(end of loop along transforms)
            kk = kk + 2 * la
         end do
!-----(end of loop on nonzero k)
!-----(end of radix-2 pass)
!
         la = 2 * la
         go to 400
!
!  central radix-8 pass
!  --------------------

300      continue
         if (m8 .eq. 0) go to 400
         jstep = (n * inc) / (8 * la)
         jstepl = jstep - ninc
         mu = mod(inq, 8)
         if (isign .eq. -1) mu = 8 - mu
         c1 = 1.0
         if (mu .eq. 3 .or. mu .eq. 7) c1 = -1.0
         c2 = sqrt(0.5)
         if (mu .eq. 3 .or. mu .eq. 5) c2 = -c2
         c3 = c1 * c2
!
!  stage 1
!  -------
         do k = 0, jstep - ink, ink
            do jjj = k, (n - 1) * inc, 8 * jstep
               ja = istart + jjj
!
!     "transverse" loop
!     -----------------
               do nu = 1, inq
                  jb = ja + jstepl
                  if (jb .lt. istart) jb = jb + ninc
                  jc = jb + jstepl
                  if (jc .lt. istart) jc = jc + ninc
                  jd = jc + jstepl
                  if (jd .lt. istart) jd = jd + ninc
                  je = jd + jstepl
                  if (je .lt. istart) je = je + ninc
                  jf = je + jstepl
                  if (jf .lt. istart) jf = jf + ninc
                  jg = jf + jstepl
                  if (jg .lt. istart) jg = jg + ninc
                  jh = jg + jstepl
                  if (jh .lt. istart) jh = jh + ninc
                  j = 0
!cdir$ ivdep, shortloop
                  do l = 1, nvex
                     t0 = a(ja + j) - a(je + j)
                     a(ja + j) = a(ja + j) + a(je + j)
                     t1 = c1 * (a(jc + j) - a(jg + j))
                     a(je + j) = a(jc + j) + a(jg + j)
                     t2 = a(jb + j) - a(jf + j)
                     a(jc + j) = a(jb + j) + a(jf + j)
                     t3 = a(jd + j) - a(jh + j)
                     a(jg + j) = a(jd + j) + a(jh + j)
                     a(jb + j) = t0
                     a(jf + j) = t1
                     a(jd + j) = c2 * (t2 - t3)
                     a(jh + j) = c3 * (t2 + t3)
                     u0 = b(ja + j) - b(je + j)
                     b(ja + j) = b(ja + j) + b(je + j)
                     u1 = c1 * (b(jc + j) - b(jg + j))
                     b(je + j) = b(jc + j) + b(jg + j)
                     u2 = b(jb + j) - b(jf + j)
                     b(jc + j) = b(jb + j) + b(jf + j)
                     u3 = b(jd + j) - b(jh + j)
                     b(jg + j) = b(jd + j) + b(jh + j)
                     b(jb + j) = u0
                     b(jf + j) = u1
                     b(jd + j) = c2 * (u2 - u3)
                     b(jh + j) = c3 * (u2 + u3)
                     j = j + jump
                  end do
                  ja = ja + jstepx
                  if (ja .lt. istart) ja = ja + ninc
               end do
            end do
         end do
!
!  stage 2
!  -------
!
!  k=0 (no twiddle factors)
!  ------------------------
         do jjj = 0, (n - 1) * inc, 8 * jstep
            ja = istart + jjj
!
!     "transverse" loop
!     -----------------
            do nu = 1, inq
               jb = ja + jstepl
               if (jb .lt. istart) jb = jb + ninc
               jc = jb + jstepl
               if (jc .lt. istart) jc = jc + ninc
               jd = jc + jstepl
               if (jd .lt. istart) jd = jd + ninc
               je = jd + jstepl
               if (je .lt. istart) je = je + ninc
               jf = je + jstepl
               if (jf .lt. istart) jf = jf + ninc
               jg = jf + jstepl
               if (jg .lt. istart) jg = jg + ninc
               jh = jg + jstepl
               if (jh .lt. istart) jh = jh + ninc
               j = 0
!cdir$ ivdep, shortloop
               do l = 1, nvex
                  t0 = a(ja + j) + a(je + j)
                  t2 = a(ja + j) - a(je + j)
                  t1 = a(jc + j) + a(jg + j)
                  t3 = c1 * (a(jc + j) - a(jg + j))
                  u0 = b(ja + j) + b(je + j)
                  u2 = b(ja + j) - b(je + j)
                  u1 = b(jc + j) + b(jg + j)
                  u3 = c1 * (b(jc + j) - b(jg + j))
                  a(ja + j) = t0 + t1
                  a(je + j) = t0 - t1
                  b(ja + j) = u0 + u1
                  b(je + j) = u0 - u1
                  a(jc + j) = t2 - u3
                  a(jg + j) = t2 + u3
                  b(jc + j) = u2 + t3
                  b(jg + j) = u2 - t3
                  t0 = a(jb + j) + a(jd + j)
                  t2 = a(jb + j) - a(jd + j)
                  t1 = a(jf + j) - a(jh + j)
                  t3 = a(jf + j) + a(jh + j)
                  u0 = b(jb + j) + b(jd + j)
                  u2 = b(jb + j) - b(jd + j)
                  u1 = b(jf + j) - b(jh + j)
                  u3 = b(jf + j) + b(jh + j)
                  a(jb + j) = t0 - u3
                  a(jh + j) = t0 + u3
                  b(jb + j) = u0 + t3
                  b(jh + j) = u0 - t3
                  a(jd + j) = t2 + u1
                  a(jf + j) = t2 - u1
                  b(jd + j) = u2 - t1
                  b(jf + j) = u2 + t1
                  j = j + jump
               end do
               ja = ja + jstepx
               if (ja .lt. istart) ja = ja + ninc
            end do
         end do
!
         if (n2 .eq. 8) go to 490
!
!  loop on nonzero k
!  -----------------
         kk = 2 * la
!
         do k = ink, jstep - ink, ink
!
            co1 = trigs(kk + 1)
            si1 = s * trigs(kk + 2)
            co2 = trigs(2 * kk + 1)
            si2 = s * trigs(2 * kk + 2)
            co3 = trigs(3 * kk + 1)
            si3 = s * trigs(3 * kk + 2)
            co4 = trigs(4 * kk + 1)
            si4 = s * trigs(4 * kk + 2)
            co5 = trigs(5 * kk + 1)
            si5 = s * trigs(5 * kk + 2)
            co6 = trigs(6 * kk + 1)
            si6 = s * trigs(6 * kk + 2)
            co7 = trigs(7 * kk + 1)
            si7 = s * trigs(7 * kk + 2)
!
            do jjj = k, (n - 1) * inc, 8 * jstep
               ja = istart + jjj
!
!     "transverse" loop
!     -----------------
               do nu = 1, inq
                  jb = ja + jstepl
                  if (jb .lt. istart) jb = jb + ninc
                  jc = jb + jstepl
                  if (jc .lt. istart) jc = jc + ninc
                  jd = jc + jstepl
                  if (jd .lt. istart) jd = jd + ninc
                  je = jd + jstepl
                  if (je .lt. istart) je = je + ninc
                  jf = je + jstepl
                  if (jf .lt. istart) jf = jf + ninc
                  jg = jf + jstepl
                  if (jg .lt. istart) jg = jg + ninc
                  jh = jg + jstepl
                  if (jh .lt. istart) jh = jh + ninc
                  j = 0
!cdir$ ivdep, shortloop
                  do l = 1, nvex
                     t0 = a(ja + j) + a(je + j)
                     t2 = a(ja + j) - a(je + j)
                     t1 = a(jc + j) + a(jg + j)
                     t3 = c1 * (a(jc + j) - a(jg + j))
                     u0 = b(ja + j) + b(je + j)
                     u2 = b(ja + j) - b(je + j)
                     u1 = b(jc + j) + b(jg + j)
                     u3 = c1 * (b(jc + j) - b(jg + j))
                     a(ja + j) = t0 + t1
                     b(ja + j) = u0 + u1
                     a(je + j) = co4 * (t0 - t1) - si4 * (u0 - u1)
                     b(je + j) = si4 * (t0 - t1) + co4 * (u0 - u1)
                     a(jc + j) = co2 * (t2 - u3) - si2 * (u2 + t3)
                     b(jc + j) = si2 * (t2 - u3) + co2 * (u2 + t3)
                     a(jg + j) = co6 * (t2 + u3) - si6 * (u2 - t3)
                     b(jg + j) = si6 * (t2 + u3) + co6 * (u2 - t3)
                     t0 = a(jb + j) + a(jd + j)
                     t2 = a(jb + j) - a(jd + j)
                     t1 = a(jf + j) - a(jh + j)
                     t3 = a(jf + j) + a(jh + j)
                     u0 = b(jb + j) + b(jd + j)
                     u2 = b(jb + j) - b(jd + j)
                     u1 = b(jf + j) - b(jh + j)
                     u3 = b(jf + j) + b(jh + j)
                     a(jb + j) = co1 * (t0 - u3) - si1 * (u0 + t3)
                     b(jb + j) = si1 * (t0 - u3) + co1 * (u0 + t3)
                     a(jh + j) = co7 * (t0 + u3) - si7 * (u0 - t3)
                     b(jh + j) = si7 * (t0 + u3) + co7 * (u0 - t3)
                     a(jd + j) = co3 * (t2 + u1) - si3 * (u2 - t1)
                     b(jd + j) = si3 * (t2 + u1) + co3 * (u2 - t1)
                     a(jf + j) = co5 * (t2 - u1) - si5 * (u2 + t1)
                     b(jf + j) = si5 * (t2 - u1) + co5 * (u2 + t1)
                     j = j + jump
                  end do
                  ja = ja + jstepx
                  if (ja .lt. istart) ja = ja + ninc
               end do
            end do
            kk = kk + 2 * la
         end do
!
         la = 8 * la
!
!  loop on type ii radix-4 passes
!  ------------------------------
400      continue
         mu = mod(inq, 4)
         if (isign .eq. -1) mu = 4 - mu
         ss = 1.0
         if (mu .eq. 3) ss = -1.0
!
         do ipass = mh + 1, m
            jstep = (n * inc) / (4 * la)
            jstepl = jstep - ninc
            laincl = la * ink - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
            do ll = 0, (la - 1) * ink, 4 * jstep
!
               do jjj = ll, (n - 1) * inc, 4 * la * ink
                  ja = istart + jjj
!
!     "transverse" loop
!     -----------------
                  do nu = 1, inq
                     jb = ja + jstepl
                     if (jb .lt. istart) jb = jb + ninc
                     jc = jb + jstepl
                     if (jc .lt. istart) jc = jc + ninc
                     jd = jc + jstepl
                     if (jd .lt. istart) jd = jd + ninc
                     je = ja + laincl
                     if (je .lt. istart) je = je + ninc
                     jf = je + jstepl
                     if (jf .lt. istart) jf = jf + ninc
                     jg = jf + jstepl
                     if (jg .lt. istart) jg = jg + ninc
                     jh = jg + jstepl
                     if (jh .lt. istart) jh = jh + ninc
                     ji = je + laincl
                     if (ji .lt. istart) ji = ji + ninc
                     jj = ji + jstepl
                     if (jj .lt. istart) jj = jj + ninc
                     jk = jj + jstepl
                     if (jk .lt. istart) jk = jk + ninc
                     jl = jk + jstepl
                     if (jl .lt. istart) jl = jl + ninc
                     jm = ji + laincl
                     if (jm .lt. istart) jm = jm + ninc
                     jn = jm + jstepl
                     if (jn .lt. istart) jn = jn + ninc
                     jo = jn + jstepl
                     if (jo .lt. istart) jo = jo + ninc
                     jp = jo + jstepl
                     if (jp .lt. istart) jp = jp + ninc
                     j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
                     do l = 1, nvex
                        t0 = a(ja + j) + a(jc + j)
                        t2 = a(ja + j) - a(jc + j)
                        t1 = a(jb + j) + a(jd + j)
                        t3 = ss * (a(jb + j) - a(jd + j))
                        a(jc + j) = a(ji + j)
                        u0 = b(ja + j) + b(jc + j)
                        u2 = b(ja + j) - b(jc + j)
                        u1 = b(jb + j) + b(jd + j)
                        u3 = ss * (b(jb + j) - b(jd + j))
                        a(jb + j) = a(je + j)
                        a(ja + j) = t0 + t1
                        a(ji + j) = t0 - t1
                        b(ja + j) = u0 + u1
                        b(jc + j) = u0 - u1
                        b(jd + j) = b(jm + j)
                        a(je + j) = t2 - u3
                        a(jd + j) = t2 + u3
                        b(jb + j) = u2 + t3
                        b(jm + j) = u2 - t3
!----------------------
                        t0 = a(jb + j) + a(jg + j)
                        t2 = a(jb + j) - a(jg + j)
                        t1 = a(jf + j) + a(jh + j)
                        t3 = ss * (a(jf + j) - a(jh + j))
                        a(jg + j) = a(jj + j)
                        u0 = b(je + j) + b(jg + j)
                        u2 = b(je + j) - b(jg + j)
                        u1 = b(jf + j) + b(jh + j)
                        u3 = ss * (b(jf + j) - b(jh + j))
                        b(je + j) = b(jb + j)
                        a(jb + j) = t0 + t1
                        a(jj + j) = t0 - t1
                        b(jg + j) = b(jj + j)
                        b(jb + j) = u0 + u1
                        b(jj + j) = u0 - u1
                        a(jf + j) = t2 - u3
                        a(jh + j) = t2 + u3
                        b(jf + j) = u2 + t3
                        b(jh + j) = u2 - t3
!----------------------
                        t0 = a(jc + j) + a(jk + j)
                        t2 = a(jc + j) - a(jk + j)
                        t1 = a(jg + j) + a(jl + j)
                        t3 = ss * (a(jg + j) - a(jl + j))
                        u0 = b(ji + j) + b(jk + j)
                        u2 = b(ji + j) - b(jk + j)
                        a(jl + j) = a(jo + j)
                        u1 = b(jg + j) + b(jl + j)
                        u3 = ss * (b(jg + j) - b(jl + j))
                        b(ji + j) = b(jc + j)
                        a(jc + j) = t0 + t1
                        a(jk + j) = t0 - t1
                        b(jl + j) = b(jo + j)
                        b(jc + j) = u0 + u1
                        b(jk + j) = u0 - u1
                        a(jg + j) = t2 - u3
                        a(jo + j) = t2 + u3
                        b(jg + j) = u2 + t3
                        b(jo + j) = u2 - t3
!----------------------
                        t0 = a(jm + j) + a(jl + j)
                        t2 = a(jm + j) - a(jl + j)
                        t1 = a(jn + j) + a(jp + j)
                        t3 = ss * (a(jn + j) - a(jp + j))
                        a(jm + j) = a(jd + j)
                        u0 = b(jd + j) + b(jl + j)
                        u2 = b(jd + j) - b(jl + j)
                        u1 = b(jn + j) + b(jp + j)
                        u3 = ss * (b(jn + j) - b(jp + j))
                        a(jn + j) = a(jh + j)
                        a(jd + j) = t0 + t1
                        a(jl + j) = t0 - t1
                        b(jd + j) = u0 + u1
                        b(jl + j) = u0 - u1
                        b(jn + j) = b(jh + j)
                        a(jh + j) = t2 - u3
                        a(jp + j) = t2 + u3
                        b(jh + j) = u2 + t3
                        b(jp + j) = u2 - t3
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
               co3 = trigs(3 * kk + 1)
               si3 = s * trigs(3 * kk + 2)
!
!  double loop along first transform in block
!  ------------------------------------------
               do ll = k, (la - 1) * ink, 4 * jstep
!
                  do jjj = ll, (n - 1) * inc, 4 * la * ink
                     ja = istart + jjj
!
!     "transverse" loop
!     -----------------
                     do nu = 1, inq
                        jb = ja + jstepl
                        if (jb .lt. istart) jb = jb + ninc
                        jc = jb + jstepl
                        if (jc .lt. istart) jc = jc + ninc
                        jd = jc + jstepl
                        if (jd .lt. istart) jd = jd + ninc
                        je = ja + laincl
                        if (je .lt. istart) je = je + ninc
                        jf = je + jstepl
                        if (jf .lt. istart) jf = jf + ninc
                        jg = jf + jstepl
                        if (jg .lt. istart) jg = jg + ninc
                        jh = jg + jstepl
                        if (jh .lt. istart) jh = jh + ninc
                        ji = je + laincl
                        if (ji .lt. istart) ji = ji + ninc
                        jj = ji + jstepl
                        if (jj .lt. istart) jj = jj + ninc
                        jk = jj + jstepl
                        if (jk .lt. istart) jk = jk + ninc
                        jl = jk + jstepl
                        if (jl .lt. istart) jl = jl + ninc
                        jm = ji + laincl
                        if (jm .lt. istart) jm = jm + ninc
                        jn = jm + jstepl
                        if (jn .lt. istart) jn = jn + ninc
                        jo = jn + jstepl
                        if (jo .lt. istart) jo = jo + ninc
                        jp = jo + jstepl
                        if (jp .lt. istart) jp = jp + ninc
                        j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
                        do l = 1, nvex
                           t0 = a(ja + j) + a(jc + j)
                           t2 = a(ja + j) - a(jc + j)
                           t1 = a(jb + j) + a(jd + j)
                           t3 = ss * (a(jb + j) - a(jd + j))
                           a(jc + j) = a(ji + j)
                           u0 = b(ja + j) + b(jc + j)
                           u2 = b(ja + j) - b(jc + j)
                           u1 = b(jb + j) + b(jd + j)
                           u3 = ss * (b(jb + j) - b(jd + j))
                           a(jb + j) = a(je + j)
                           a(ja + j) = t0 + t1
                           b(ja + j) = u0 + u1
                           a(je + j) = co1 * (t2 - u3) - si1 * (u2 + t3)
                           b(jb + j) = si1 * (t2 - u3) + co1 * (u2 + t3)
                           b(jd + j) = b(jm + j)
                           a(ji + j) = co2 * (t0 - t1) - si2 * (u0 - u1)
                           b(jc + j) = si2 * (t0 - t1) + co2 * (u0 - u1)
                           a(jd + j) = co3 * (t2 + u3) - si3 * (u2 - t3)
                           b(jm + j) = si3 * (t2 + u3) + co3 * (u2 - t3)
!----------------------------------------
                           t0 = a(jb + j) + a(jg + j)
                           t2 = a(jb + j) - a(jg + j)
                           t1 = a(jf + j) + a(jh + j)
                           t3 = ss * (a(jf + j) - a(jh + j))
                           a(jg + j) = a(jj + j)
                           u0 = b(je + j) + b(jg + j)
                           u2 = b(je + j) - b(jg + j)
                           u1 = b(jf + j) + b(jh + j)
                           u3 = ss * (b(jf + j) - b(jh + j))
                           b(je + j) = b(jb + j)
                           a(jb + j) = t0 + t1
                           b(jb + j) = u0 + u1
                           b(jg + j) = b(jj + j)
                           a(jf + j) = co1 * (t2 - u3) - si1 * (u2 + t3)
                           b(jf + j) = si1 * (t2 - u3) + co1 * (u2 + t3)
                           a(jj + j) = co2 * (t0 - t1) - si2 * (u0 - u1)
                           b(jj + j) = si2 * (t0 - t1) + co2 * (u0 - u1)
                           a(jh + j) = co3 * (t2 + u3) - si3 * (u2 - t3)
                           b(jh + j) = si3 * (t2 + u3) + co3 * (u2 - t3)
!----------------------------------------
                           t0 = a(jc + j) + a(jk + j)
                           t2 = a(jc + j) - a(jk + j)
                           t1 = a(jg + j) + a(jl + j)
                           t3 = ss * (a(jg + j) - a(jl + j))
                           u0 = b(ji + j) + b(jk + j)
                           u2 = b(ji + j) - b(jk + j)
                           a(jl + j) = a(jo + j)
                           u1 = b(jg + j) + b(jl + j)
                           u3 = ss * (b(jg + j) - b(jl + j))
                           b(ji + j) = b(jc + j)
                           a(jc + j) = t0 + t1
                           b(jc + j) = u0 + u1
                           b(jl + j) = b(jo + j)
                           a(jg + j) = co1 * (t2 - u3) - si1 * (u2 + t3)
                           b(jg + j) = si1 * (t2 - u3) + co1 * (u2 + t3)
                           a(jk + j) = co2 * (t0 - t1) - si2 * (u0 - u1)
                           b(jk + j) = si2 * (t0 - t1) + co2 * (u0 - u1)
                           a(jo + j) = co3 * (t2 + u3) - si3 * (u2 - t3)
                           b(jo + j) = si3 * (t2 + u3) + co3 * (u2 - t3)
!----------------------------------------
                           t0 = a(jm + j) + a(jl + j)
                           t2 = a(jm + j) - a(jl + j)
                           t1 = a(jn + j) + a(jp + j)
                           t3 = ss * (a(jn + j) - a(jp + j))
                           a(jm + j) = a(jd + j)
                           u0 = b(jd + j) + b(jl + j)
                           u2 = b(jd + j) - b(jl + j)
                           a(jn + j) = a(jh + j)
                           u1 = b(jn + j) + b(jp + j)
                           u3 = ss * (b(jn + j) - b(jp + j))
                           b(jn + j) = b(jh + j)
                           a(jd + j) = t0 + t1
                           b(jd + j) = u0 + u1
                           a(jh + j) = co1 * (t2 - u3) - si1 * (u2 + t3)
                           b(jh + j) = si1 * (t2 - u3) + co1 * (u2 + t3)
                           a(jl + j) = co2 * (t0 - t1) - si2 * (u0 - u1)
                           b(jl + j) = si2 * (t0 - t1) + co2 * (u0 - u1)
                           a(jp + j) = co3 * (t2 + u3) - si3 * (u2 - t3)
                           b(jp + j) = si3 * (t2 + u3) + co3 * (u2 - t3)
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
            la = 4 * la
         end do
!-----( end of loop on type ii radix-4 passes )
!-----( nvex transforms completed)
490      continue
         istart = istart + nvex * jump
      end do
!-----( end of loop on blocks of transforms )
!
      return
   end subroutine gpfa2f

end module gpfa_radix2
