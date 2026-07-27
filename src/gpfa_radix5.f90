module gpfa_radix5
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
contains

   subroutine gpfa5f(a, b, trigs, inc, jump, n, mm, lot, isign)
      implicit none
      integer :: inc, jump, n, mm, lot, isign, lvr, inq, jstepx, ninc, ink, &
                 m, mh, nblox, left, nb, nvex, la, mu, ipass, jstep, jstepl, &
                 jjj, ja, nu, jb, jc, jd, j, l, kk, k, je, jf, jg, jh, laincl, ji, jj, jk, &
                 jl, jm, jn, jo, jp, n5, jq, jr, js, jt, ju, jv, jw, jx, jy, istart, ll
      real(real64) :: a(*), b(*), trigs(*), s, t2, t1, t3, u2, u1, u3, co1, si1, &
                      co2, si2, co3, si3, c1, c2, c3, co4, si4, &
                      sin36, sin72, qrt5, t4, t5, t6, t7, t8, t9, t10, t11, &
                      ax, bx, u4, u5, u6, u7, u8, u9, u10, u11
      data sin36/0.587785252292473/, sin72/0.951056516295154/, &
         qrt5/0.559016994374947/
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
      n5 = 5**mm
      inq = n / n5
      jstepx = (n5 - n) * inc
      ninc = n * inc
      ink = inc * inq
      mu = mod(inq, 5)
      if (isign .eq. -1) mu = 5 - mu
!
      m = mm
      mh = (m + 1) / 2
      s = real(isign, kind=real64)
      c1 = qrt5
      c2 = sin72
      c3 = sin36
      if (mu .eq. 2 .or. mu .eq. 3) then
         c1 = -c1
         c2 = sin36
         c3 = sin72
      end if
      if (mu .eq. 3 .or. mu .eq. 4) c2 = -c2
      if (mu .eq. 2 .or. mu .eq. 4) c3 = -c3
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
!  loop on type i radix-5 passes
!  -----------------------------
         do ipass = 1, mh
            jstep = (n * inc) / (5 * la)
            jstepl = jstep - ninc
            kk = 0
!
!  loop on k
!  ---------
            do k = 0, jstep - ink, ink
!
               if (k .gt. 0) then
                  co1 = trigs(kk + 1)
                  si1 = s * trigs(kk + 2)
                  co2 = trigs(2 * kk + 1)
                  si2 = s * trigs(2 * kk + 2)
                  co3 = trigs(3 * kk + 1)
                  si3 = s * trigs(3 * kk + 2)
                  co4 = trigs(4 * kk + 1)
                  si4 = s * trigs(4 * kk + 2)
               end if
!
!  loop along transform
!  --------------------
               do jjj = k, (n - 1) * inc, 5 * jstep
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
                     j = 0
!
!  loop across transforms
!  ----------------------
                     if (k .eq. 0) then
!
!cdir$ ivdep, shortloop
                        do l = 1, nvex
                           t1 = a(jb + j) + a(je + j)
                           t2 = a(jc + j) + a(jd + j)
                           t3 = a(jb + j) - a(je + j)
                           t4 = a(jc + j) - a(jd + j)
                           t5 = t1 + t2
                           t6 = c1 * (t1 - t2)
                           t7 = a(ja + j) - 0.25 * t5
                           a(ja + j) = a(ja + j) + t5
                           t8 = t7 + t6
                           t9 = t7 - t6
                           t10 = c3 * t3 - c2 * t4
                           t11 = c2 * t3 + c3 * t4
                           u1 = b(jb + j) + b(je + j)
                           u2 = b(jc + j) + b(jd + j)
                           u3 = b(jb + j) - b(je + j)
                           u4 = b(jc + j) - b(jd + j)
                           u5 = u1 + u2
                           u6 = c1 * (u1 - u2)
                           u7 = b(ja + j) - 0.25 * u5
                           b(ja + j) = b(ja + j) + u5
                           u8 = u7 + u6
                           u9 = u7 - u6
                           u10 = c3 * u3 - c2 * u4
                           u11 = c2 * u3 + c3 * u4
                           a(jb + j) = t8 - u11
                           b(jb + j) = u8 + t11
                           a(je + j) = t8 + u11
                           b(je + j) = u8 - t11
                           a(jc + j) = t9 - u10
                           b(jc + j) = u9 + t10
                           a(jd + j) = t9 + u10
                           b(jd + j) = u9 - t10
                           j = j + jump
                        end do
!
                     else
!
!cdir$ ivdep,shortloop
                        do l = 1, nvex
                           t1 = a(jb + j) + a(je + j)
                           t2 = a(jc + j) + a(jd + j)
                           t3 = a(jb + j) - a(je + j)
                           t4 = a(jc + j) - a(jd + j)
                           t5 = t1 + t2
                           t6 = c1 * (t1 - t2)
                           t7 = a(ja + j) - 0.25 * t5
                           a(ja + j) = a(ja + j) + t5
                           t8 = t7 + t6
                           t9 = t7 - t6
                           t10 = c3 * t3 - c2 * t4
                           t11 = c2 * t3 + c3 * t4
                           u1 = b(jb + j) + b(je + j)
                           u2 = b(jc + j) + b(jd + j)
                           u3 = b(jb + j) - b(je + j)
                           u4 = b(jc + j) - b(jd + j)
                           u5 = u1 + u2
                           u6 = c1 * (u1 - u2)
                           u7 = b(ja + j) - 0.25 * u5
                           b(ja + j) = b(ja + j) + u5
                           u8 = u7 + u6
                           u9 = u7 - u6
                           u10 = c3 * u3 - c2 * u4
                           u11 = c2 * u3 + c3 * u4
                           a(jb + j) = co1 * (t8 - u11) - si1 * (u8 + t11)
                           b(jb + j) = si1 * (t8 - u11) + co1 * (u8 + t11)
                           a(je + j) = co4 * (t8 + u11) - si4 * (u8 - t11)
                           b(je + j) = si4 * (t8 + u11) + co4 * (u8 - t11)
                           a(jc + j) = co2 * (t9 - u10) - si2 * (u9 + t10)
                           b(jc + j) = si2 * (t9 - u10) + co2 * (u9 + t10)
                           a(jd + j) = co3 * (t9 + u10) - si3 * (u9 - t10)
                           b(jd + j) = si3 * (t9 + u10) + co3 * (u9 - t10)
                           j = j + jump
                        end do
!
                     end if
!
!-----( end of loop across transforms )
!
                     ja = ja + jstepx
                     if (ja .lt. istart) ja = ja + ninc
                  end do
               end do
!-----( end of loop along transforms )
               kk = kk + 2 * la
            end do
!-----( end of loop on nonzero k )
            la = 5 * la
         end do
!-----( end of loop on type i radix-5 passes)
!
         if (n .eq. 5) go to 490
!
!  loop on type ii radix-5 passes
!  ------------------------------
!
         do ipass = mh + 1, m
            jstep = (n * inc) / (5 * la)
            jstepl = jstep - ninc
            laincl = la * ink - ninc
            kk = 0
!
!     loop on k
!     ---------
            do k = 0, jstep - ink, ink
!
               if (k .gt. 0) then
                  co1 = trigs(kk + 1)
                  si1 = s * trigs(kk + 2)
                  co2 = trigs(2 * kk + 1)
                  si2 = s * trigs(2 * kk + 2)
                  co3 = trigs(3 * kk + 1)
                  si3 = s * trigs(3 * kk + 2)
                  co4 = trigs(4 * kk + 1)
                  si4 = s * trigs(4 * kk + 2)
               end if
!
!  double loop along first transform in block
!  ------------------------------------------
               do ll = k, (la - 1) * ink, 5 * jstep
!
                  do jjj = ll, (n - 1) * inc, 5 * la * ink
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
                        jf = ja + laincl
                        if (jf .lt. istart) jf = jf + ninc
                        jg = jf + jstepl
                        if (jg .lt. istart) jg = jg + ninc
                        jh = jg + jstepl
                        if (jh .lt. istart) jh = jh + ninc
                        ji = jh + jstepl
                        if (ji .lt. istart) ji = ji + ninc
                        jj = ji + jstepl
                        if (jj .lt. istart) jj = jj + ninc
                        jk = jf + laincl
                        if (jk .lt. istart) jk = jk + ninc
                        jl = jk + jstepl
                        if (jl .lt. istart) jl = jl + ninc
                        jm = jl + jstepl
                        if (jm .lt. istart) jm = jm + ninc
                        jn = jm + jstepl
                        if (jn .lt. istart) jn = jn + ninc
                        jo = jn + jstepl
                        if (jo .lt. istart) jo = jo + ninc
                        jp = jk + laincl
                        if (jp .lt. istart) jp = jp + ninc
                        jq = jp + jstepl
                        if (jq .lt. istart) jq = jq + ninc
                        jr = jq + jstepl
                        if (jr .lt. istart) jr = jr + ninc
                        js = jr + jstepl
                        if (js .lt. istart) js = js + ninc
                        jt = js + jstepl
                        if (jt .lt. istart) jt = jt + ninc
                        ju = jp + laincl
                        if (ju .lt. istart) ju = ju + ninc
                        jv = ju + jstepl
                        if (jv .lt. istart) jv = jv + ninc
                        jw = jv + jstepl
                        if (jw .lt. istart) jw = jw + ninc
                        jx = jw + jstepl
                        if (jx .lt. istart) jx = jx + ninc
                        jy = jx + jstepl
                        if (jy .lt. istart) jy = jy + ninc
                        j = 0
!
!  loop across transforms
!  ----------------------
                        if (k .eq. 0) then
!
!cdir$ ivdep, shortloop
                           do l = 1, nvex
                              t1 = a(jb + j) + a(je + j)
                              t2 = a(jc + j) + a(jd + j)
                              t3 = a(jb + j) - a(je + j)
                              t4 = a(jc + j) - a(jd + j)
                              a(jb + j) = a(jf + j)
                              t5 = t1 + t2
                              t6 = c1 * (t1 - t2)
                              t7 = a(ja + j) - 0.25 * t5
                              a(ja + j) = a(ja + j) + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(jc + j) = a(jk + j)
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              u1 = b(jb + j) + b(je + j)
                              u2 = b(jc + j) + b(jd + j)
                              u3 = b(jb + j) - b(je + j)
                              u4 = b(jc + j) - b(jd + j)
                              b(jb + j) = b(jf + j)
                              u5 = u1 + u2
                              u6 = c1 * (u1 - u2)
                              u7 = b(ja + j) - 0.25 * u5
                              b(ja + j) = b(ja + j) + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(jc + j) = b(jk + j)
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              a(jf + j) = t8 - u11
                              b(jf + j) = u8 + t11
                              a(je + j) = t8 + u11
                              b(je + j) = u8 - t11
                              a(jk + j) = t9 - u10
                              b(jk + j) = u9 + t10
                              a(jd + j) = t9 + u10
                              b(jd + j) = u9 - t10
!----------------------
                              t1 = a(jg + j) + a(jj + j)
                              t2 = a(jh + j) + a(ji + j)
                              t3 = a(jg + j) - a(jj + j)
                              t4 = a(jh + j) - a(ji + j)
                              a(jh + j) = a(jl + j)
                              t5 = t1 + t2
                              t6 = c1 * (t1 - t2)
                              t7 = a(jb + j) - 0.25 * t5
                              a(jb + j) = a(jb + j) + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(ji + j) = a(jq + j)
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              u1 = b(jg + j) + b(jj + j)
                              u2 = b(jh + j) + b(ji + j)
                              u3 = b(jg + j) - b(jj + j)
                              u4 = b(jh + j) - b(ji + j)
                              b(jh + j) = b(jl + j)
                              u5 = u1 + u2
                              u6 = c1 * (u1 - u2)
                              u7 = b(jb + j) - 0.25 * u5
                              b(jb + j) = b(jb + j) + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(ji + j) = b(jq + j)
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              a(jg + j) = t8 - u11
                              b(jg + j) = u8 + t11
                              a(jj + j) = t8 + u11
                              b(jj + j) = u8 - t11
                              a(jl + j) = t9 - u10
                              b(jl + j) = u9 + t10
                              a(jq + j) = t9 + u10
                              b(jq + j) = u9 - t10
!----------------------
                              t1 = a(jh + j) + a(jo + j)
                              t2 = a(jm + j) + a(jn + j)
                              t3 = a(jh + j) - a(jo + j)
                              t4 = a(jm + j) - a(jn + j)
                              a(jn + j) = a(jr + j)
                              t5 = t1 + t2
                              t6 = c1 * (t1 - t2)
                              t7 = a(jc + j) - 0.25 * t5
                              a(jc + j) = a(jc + j) + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(jo + j) = a(jw + j)
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              u1 = b(jh + j) + b(jo + j)
                              u2 = b(jm + j) + b(jn + j)
                              u3 = b(jh + j) - b(jo + j)
                              u4 = b(jm + j) - b(jn + j)
                              b(jn + j) = b(jr + j)
                              u5 = u1 + u2
                              u6 = c1 * (u1 - u2)
                              u7 = b(jc + j) - 0.25 * u5
                              b(jc + j) = b(jc + j) + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(jo + j) = b(jw + j)
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              a(jh + j) = t8 - u11
                              b(jh + j) = u8 + t11
                              a(jw + j) = t8 + u11
                              b(jw + j) = u8 - t11
                              a(jm + j) = t9 - u10
                              b(jm + j) = u9 + t10
                              a(jr + j) = t9 + u10
                              b(jr + j) = u9 - t10
!----------------------
                              t1 = a(ji + j) + a(jt + j)
                              t2 = a(jn + j) + a(js + j)
                              t3 = a(ji + j) - a(jt + j)
                              t4 = a(jn + j) - a(js + j)
                              a(jt + j) = a(jx + j)
                              t5 = t1 + t2
                              t6 = c1 * (t1 - t2)
                              t7 = a(jp + j) - 0.25 * t5
                              ax = a(jp + j) + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(jp + j) = a(jd + j)
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              a(jd + j) = ax
                              u1 = b(ji + j) + b(jt + j)
                              u2 = b(jn + j) + b(js + j)
                              u3 = b(ji + j) - b(jt + j)
                              u4 = b(jn + j) - b(js + j)
                              b(jt + j) = b(jx + j)
                              u5 = u1 + u2
                              u6 = c1 * (u1 - u2)
                              u7 = b(jp + j) - 0.25 * u5
                              bx = b(jp + j) + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(jp + j) = b(jd + j)
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              b(jd + j) = bx
                              a(ji + j) = t8 - u11
                              b(ji + j) = u8 + t11
                              a(jx + j) = t8 + u11
                              b(jx + j) = u8 - t11
                              a(jn + j) = t9 - u10
                              b(jn + j) = u9 + t10
                              a(js + j) = t9 + u10
                              b(js + j) = u9 - t10
!----------------------
                              t1 = a(jv + j) + a(jy + j)
                              t2 = a(jo + j) + a(jt + j)
                              t3 = a(jv + j) - a(jy + j)
                              t4 = a(jo + j) - a(jt + j)
                              a(jv + j) = a(jj + j)
                              t5 = t1 + t2
                              t6 = c1 * (t1 - t2)
                              t7 = a(ju + j) - 0.25 * t5
                              ax = a(ju + j) + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(ju + j) = a(je + j)
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              a(je + j) = ax
                              u1 = b(jv + j) + b(jy + j)
                              u2 = b(jo + j) + b(jt + j)
                              u3 = b(jv + j) - b(jy + j)
                              u4 = b(jo + j) - b(jt + j)
                              b(jv + j) = b(jj + j)
                              u5 = u1 + u2
                              u6 = c1 * (u1 - u2)
                              u7 = b(ju + j) - 0.25 * u5
                              bx = b(ju + j) + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(ju + j) = b(je + j)
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              b(je + j) = bx
                              a(jj + j) = t8 - u11
                              b(jj + j) = u8 + t11
                              a(jy + j) = t8 + u11
                              b(jy + j) = u8 - t11
                              a(jo + j) = t9 - u10
                              b(jo + j) = u9 + t10
                              a(jt + j) = t9 + u10
                              b(jt + j) = u9 - t10
                              j = j + jump
                           end do
!
                        else
!
!cdir$ ivdep, shortloop
                           do l = 1, nvex
                              t1 = a(jb + j) + a(je + j)
                              t2 = a(jc + j) + a(jd + j)
                              t3 = a(jb + j) - a(je + j)
                              t4 = a(jc + j) - a(jd + j)
                              a(jb + j) = a(jf + j)
                              t5 = t1 + t2
                              t6 = c1 * (t1 - t2)
                              t7 = a(ja + j) - 0.25 * t5
                              a(ja + j) = a(ja + j) + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(jc + j) = a(jk + j)
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              u1 = b(jb + j) + b(je + j)
                              u2 = b(jc + j) + b(jd + j)
                              u3 = b(jb + j) - b(je + j)
                              u4 = b(jc + j) - b(jd + j)
                              b(jb + j) = b(jf + j)
                              u5 = u1 + u2
                              u6 = c1 * (u1 - u2)
                              u7 = b(ja + j) - 0.25 * u5
                              b(ja + j) = b(ja + j) + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(jc + j) = b(jk + j)
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              a(jf + j) = co1 * (t8 - u11) - si1 * (u8 + t11)
                              b(jf + j) = si1 * (t8 - u11) + co1 * (u8 + t11)
                              a(je + j) = co4 * (t8 + u11) - si4 * (u8 - t11)
                              b(je + j) = si4 * (t8 + u11) + co4 * (u8 - t11)
                              a(jk + j) = co2 * (t9 - u10) - si2 * (u9 + t10)
                              b(jk + j) = si2 * (t9 - u10) + co2 * (u9 + t10)
                              a(jd + j) = co3 * (t9 + u10) - si3 * (u9 - t10)
                              b(jd + j) = si3 * (t9 + u10) + co3 * (u9 - t10)
!----------------------
                              t1 = a(jg + j) + a(jj + j)
                              t2 = a(jh + j) + a(ji + j)
                              t3 = a(jg + j) - a(jj + j)
                              t4 = a(jh + j) - a(ji + j)
                              a(jh + j) = a(jl + j)
                              t5 = t1 + t2
                              t6 = c1 * (t1 - t2)
                              t7 = a(jb + j) - 0.25 * t5
                              a(jb + j) = a(jb + j) + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(ji + j) = a(jq + j)
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              u1 = b(jg + j) + b(jj + j)
                              u2 = b(jh + j) + b(ji + j)
                              u3 = b(jg + j) - b(jj + j)
                              u4 = b(jh + j) - b(ji + j)
                              b(jh + j) = b(jl + j)
                              u5 = u1 + u2
                              u6 = c1 * (u1 - u2)
                              u7 = b(jb + j) - 0.25 * u5
                              b(jb + j) = b(jb + j) + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(ji + j) = b(jq + j)
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              a(jg + j) = co1 * (t8 - u11) - si1 * (u8 + t11)
                              b(jg + j) = si1 * (t8 - u11) + co1 * (u8 + t11)
                              a(jj + j) = co4 * (t8 + u11) - si4 * (u8 - t11)
                              b(jj + j) = si4 * (t8 + u11) + co4 * (u8 - t11)
                              a(jl + j) = co2 * (t9 - u10) - si2 * (u9 + t10)
                              b(jl + j) = si2 * (t9 - u10) + co2 * (u9 + t10)
                              a(jq + j) = co3 * (t9 + u10) - si3 * (u9 - t10)
                              b(jq + j) = si3 * (t9 + u10) + co3 * (u9 - t10)
!----------------------
                              t1 = a(jh + j) + a(jo + j)
                              t2 = a(jm + j) + a(jn + j)
                              t3 = a(jh + j) - a(jo + j)
                              t4 = a(jm + j) - a(jn + j)
                              a(jn + j) = a(jr + j)
                              t5 = t1 + t2
                              t6 = c1 * (t1 - t2)
                              t7 = a(jc + j) - 0.25 * t5
                              a(jc + j) = a(jc + j) + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(jo + j) = a(jw + j)
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              u1 = b(jh + j) + b(jo + j)
                              u2 = b(jm + j) + b(jn + j)
                              u3 = b(jh + j) - b(jo + j)
                              u4 = b(jm + j) - b(jn + j)
                              b(jn + j) = b(jr + j)
                              u5 = u1 + u2
                              u6 = c1 * (u1 - u2)
                              u7 = b(jc + j) - 0.25 * u5
                              b(jc + j) = b(jc + j) + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(jo + j) = b(jw + j)
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              a(jh + j) = co1 * (t8 - u11) - si1 * (u8 + t11)
                              b(jh + j) = si1 * (t8 - u11) + co1 * (u8 + t11)
                              a(jw + j) = co4 * (t8 + u11) - si4 * (u8 - t11)
                              b(jw + j) = si4 * (t8 + u11) + co4 * (u8 - t11)
                              a(jm + j) = co2 * (t9 - u10) - si2 * (u9 + t10)
                              b(jm + j) = si2 * (t9 - u10) + co2 * (u9 + t10)
                              a(jr + j) = co3 * (t9 + u10) - si3 * (u9 - t10)
                              b(jr + j) = si3 * (t9 + u10) + co3 * (u9 - t10)
!----------------------
                              t1 = a(ji + j) + a(jt + j)
                              t2 = a(jn + j) + a(js + j)
                              t3 = a(ji + j) - a(jt + j)
                              t4 = a(jn + j) - a(js + j)
                              a(jt + j) = a(jx + j)
                              t5 = t1 + t2
                              t6 = c1 * (t1 - t2)
                              t7 = a(jp + j) - 0.25 * t5
                              ax = a(jp + j) + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(jp + j) = a(jd + j)
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              a(jd + j) = ax
                              u1 = b(ji + j) + b(jt + j)
                              u2 = b(jn + j) + b(js + j)
                              u3 = b(ji + j) - b(jt + j)
                              u4 = b(jn + j) - b(js + j)
                              b(jt + j) = b(jx + j)
                              u5 = u1 + u2
                              u6 = c1 * (u1 - u2)
                              u7 = b(jp + j) - 0.25 * u5
                              bx = b(jp + j) + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(jp + j) = b(jd + j)
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              b(jd + j) = bx
                              a(ji + j) = co1 * (t8 - u11) - si1 * (u8 + t11)
                              b(ji + j) = si1 * (t8 - u11) + co1 * (u8 + t11)
                              a(jx + j) = co4 * (t8 + u11) - si4 * (u8 - t11)
                              b(jx + j) = si4 * (t8 + u11) + co4 * (u8 - t11)
                              a(jn + j) = co2 * (t9 - u10) - si2 * (u9 + t10)
                              b(jn + j) = si2 * (t9 - u10) + co2 * (u9 + t10)
                              a(js + j) = co3 * (t9 + u10) - si3 * (u9 - t10)
                              b(js + j) = si3 * (t9 + u10) + co3 * (u9 - t10)
!----------------------
                              t1 = a(jv + j) + a(jy + j)
                              t2 = a(jo + j) + a(jt + j)
                              t3 = a(jv + j) - a(jy + j)
                              t4 = a(jo + j) - a(jt + j)
                              a(jv + j) = a(jj + j)
                              t5 = t1 + t2
                              t6 = c1 * (t1 - t2)
                              t7 = a(ju + j) - 0.25 * t5
                              ax = a(ju + j) + t5
                              t8 = t7 + t6
                              t9 = t7 - t6
                              a(ju + j) = a(je + j)
                              t10 = c3 * t3 - c2 * t4
                              t11 = c2 * t3 + c3 * t4
                              a(je + j) = ax
                              u1 = b(jv + j) + b(jy + j)
                              u2 = b(jo + j) + b(jt + j)
                              u3 = b(jv + j) - b(jy + j)
                              u4 = b(jo + j) - b(jt + j)
                              b(jv + j) = b(jj + j)
                              u5 = u1 + u2
                              u6 = c1 * (u1 - u2)
                              u7 = b(ju + j) - 0.25 * u5
                              bx = b(ju + j) + u5
                              u8 = u7 + u6
                              u9 = u7 - u6
                              b(ju + j) = b(je + j)
                              u10 = c3 * u3 - c2 * u4
                              u11 = c2 * u3 + c3 * u4
                              b(je + j) = bx
                              a(jj + j) = co1 * (t8 - u11) - si1 * (u8 + t11)
                              b(jj + j) = si1 * (t8 - u11) + co1 * (u8 + t11)
                              a(jy + j) = co4 * (t8 + u11) - si4 * (u8 - t11)
                              b(jy + j) = si4 * (t8 + u11) + co4 * (u8 - t11)
                              a(jo + j) = co2 * (t9 - u10) - si2 * (u9 + t10)
                              b(jo + j) = si2 * (t9 - u10) + co2 * (u9 + t10)
                              a(jt + j) = co3 * (t9 + u10) - si3 * (u9 - t10)
                              b(jt + j) = si3 * (t9 + u10) + co3 * (u9 - t10)
                              j = j + jump
                           end do
!
                        end if
!
!-----(end of loop across transforms)
!
                        ja = ja + jstepx
                        if (ja .lt. istart) ja = ja + ninc
                     end do
                  end do
               end do
!-----( end of double loop for this k )
               kk = kk + 2 * la
            end do
!-----( end of loop over values of k )
            la = 5 * la
         end do
!-----( end of loop on type ii radix-5 passes )
!-----( nvex transforms completed)
490      continue
         istart = istart + nvex * jump
      end do
!-----( end of loop on blocks of transforms )
!
      return
   end subroutine gpfa5f

end module gpfa_radix5
