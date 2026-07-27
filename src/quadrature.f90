module quadrature_functions
   implicit none
contains

   subroutine qng(n, a, b, epsabs, epsrel, f, resultf, abserr, neval, ier)
      implicit none
      integer :: n, ier, k, l, neval, ipx
      real(8) a, absc, abserr, b, centr, dhlgth, epsabs, epsrel, hlgth, w10(5), w21a(5), w21b(6), &
         w43a(10), w43b(12), w87a(21), w87b(23), x1(5), x2(5), x3(11), x4(22)
      complex(8) :: fcentr(n), fval(n), fval1(n), fval2(n), fv1(5, n), fv2(5, n), fv3(5, n), fv4(5, n), &
                    resultf(n), res10(n), res21(n), res43(n), res87(n), savfun(21, n)
      external :: f
      data x1(1), x1(2), x1(3), x1(4), x1(5)/ &
         9.739065285171717E-01, 8.650633666889845E-01, &
         6.794095682990244E-01, 4.333953941292472E-01, &
         1.488743389816312E-01/
      data x2(1), x2(2), x2(3), x2(4), x2(5)/ &
         9.956571630258081E-01, 9.301574913557082E-01, &
         7.808177265864169E-01, 5.627571346686047E-01, &
         2.943928627014602E-01/
      data x3(1), x3(2), x3(3), x3(4), x3(5), x3(6), x3(7), x3(8), x3(9), x3(10), &
         x3(11)/ &
         9.993333609019321E-01, 9.874334029080889E-01, &
         9.548079348142663E-01, 9.001486957483283E-01, &
         8.251983149831142E-01, 7.321483889893050E-01, &
         6.228479705377252E-01, 4.994795740710565E-01, &
         3.649016613465808E-01, 2.222549197766013E-01, &
         7.465061746138332E-02/
      data x4(1), x4(2), x4(3), x4(4), x4(5), x4(6), x4(7), x4(8), x4(9), x4(10), &
         x4(11), x4(12), x4(13), x4(14), x4(15), x4(16), x4(17), x4(18), x4(19), &
         x4(20), x4(21), x4(22)/9.999029772627292E-01, &
         9.979898959866787E-01, 9.921754978606872E-01, &
         9.813581635727128E-01, 9.650576238583846E-01, &
         9.431676131336706E-01, 9.158064146855072E-01, &
         8.832216577713165E-01, 8.457107484624157E-01, &
         8.035576580352310E-01, 7.570057306854956E-01, &
         7.062732097873218E-01, 6.515894665011779E-01, &
         5.932233740579611E-01, 5.314936059708319E-01, &
         4.667636230420228E-01, 3.994248478592188E-01, &
         3.298748771061883E-01, 2.585035592021616E-01, &
         1.856953965683467E-01, 1.118422131799075E-01, &
         3.735212339461987E-02/
      data w10(1), w10(2), w10(3), w10(4), w10(5)/ &
         6.667134430868814E-02, 1.494513491505806E-01, &
         2.190863625159820E-01, 2.692667193099964E-01, &
         2.955242247147529E-01/
      data w21a(1), w21a(2), w21a(3), w21a(4), w21a(5)/ &
         3.255816230796473E-02, 7.503967481091995E-02, &
         1.093871588022976E-01, 1.347092173114733E-01, &
         1.477391049013385E-01/
      data w21b(1), w21b(2), w21b(3), w21b(4), w21b(5), w21b(6)/ &
         1.169463886737187E-02, 5.475589657435200E-02, &
         9.312545458369761E-02, 1.234919762620659E-01, &
         1.427759385770601E-01, 1.494455540029169E-01/
      data w43a(1), w43a(2), w43a(3), w43a(4), w43a(5), w43a(6), w43a(7), &
         w43a(8), w43a(9), w43a(10)/1.629673428966656E-02, &
         3.752287612086950E-02, 5.469490205825544E-02, &
         6.735541460947809E-02, 7.387019963239395E-02, &
         5.768556059769796E-03, 2.737189059324884E-02, &
         4.656082691042883E-02, 6.174499520144256E-02, &
         7.138726726869340E-02/
      data w43b(1), w43b(2), w43b(3), w43b(4), w43b(5), w43b(6), w43b(7), &
         w43b(8), w43b(9), w43b(10), w43b(11), w43b(12)/ &
         1.844477640212414E-03, 1.079868958589165E-02, &
         2.189536386779543E-02, 3.259746397534569E-02, &
         4.216313793519181E-02, 5.074193960018458E-02, &
         5.837939554261925E-02, 6.474640495144589E-02, &
         6.956619791235648E-02, 7.282444147183321E-02, &
         7.450775101417512E-02, 7.472214751740301E-02/
      data w87a(1), w87a(2), w87a(3), w87a(4), w87a(5), w87a(6), w87a(7), &
         w87a(8), w87a(9), w87a(10), w87a(11), w87a(12), w87a(13), w87a(14), &
         w87a(15), w87a(16), w87a(17), w87a(18), w87a(19), w87a(20), w87a(21)/ &
         8.148377384149173E-03, 1.876143820156282E-02, &
         2.734745105005229E-02, 3.367770731163793E-02, &
         3.693509982042791E-02, 2.884872430211531E-03, &
         1.368594602271270E-02, 2.328041350288831E-02, &
         3.087249761171336E-02, 3.569363363941877E-02, &
         9.152833452022414E-04, 5.399280219300471E-03, &
         1.094767960111893E-02, 1.629873169678734E-02, &
         2.108156888920384E-02, 2.537096976925383E-02, &
         2.918969775647575E-02, 3.237320246720279E-02, &
         3.478309895036514E-02, 3.641222073135179E-02, &
         3.725387550304771E-02/
      data w87b(1), w87b(2), w87b(3), w87b(4), w87b(5), w87b(6), w87b(7), &
         w87b(8), w87b(9), w87b(10), w87b(11), w87b(12), w87b(13), w87b(14), &
         w87b(15), w87b(16), w87b(17), w87b(18), w87b(19), w87b(20), w87b(21), &
         w87b(22), w87b(23)/2.741455637620724E-04, &
         1.807124155057943E-03, 4.096869282759165E-03, &
         6.758290051847379E-03, 9.549957672201647E-03, &
         1.232944765224485E-02, 1.501044734638895E-02, &
         1.754896798624319E-02, 1.993803778644089E-02, &
         2.219493596101229E-02, 2.433914712600081E-02, &
         2.637450541483921E-02, 2.828691078877120E-02, &
         3.005258112809270E-02, 3.164675137143993E-02, &
         3.305041341997850E-02, 3.425509970422606E-02, &
         3.526241266015668E-02, 3.607698962288870E-02, &
         3.669860449845609E-02, 3.712054926983258E-02, &
         3.733422875193504E-02, 3.736107376267902E-02/
!
!  Test on validity of parameters.
!
      resultf = 0.0D+00
      abserr = 0.0D+00
      neval = 0

      hlgth = 5.0D-01 * (b - a)
      dhlgth = abs(hlgth)
      centr = 5.0D-01 * (b + a)
      call f(n, centr, fcentr)
      neval = 21
      ier = 1

      do l = 1, 3
         if (l == 1) then
            res10 = 0.0D+00
            res21 = w21b(6) * fcentr
            do k = 1, 5
               absc = hlgth * x1(k)
               call f(n, centr + absc, fval1)
               call f(n, centr - absc, fval2)
               fval = fval1 + fval2
               res10 = res10 + w10(k) * fval
               res21 = res21 + w21a(k) * fval
               savfun(k, :) = fval(:)
               fv1(k, :) = fval1(:)
               fv2(k, :) = fval2(:)
            end do
            ipx = 5
            do k = 1, 5
               ipx = ipx + 1
               absc = hlgth * x2(k)
               call f(n, centr + absc, fval1)
               call f(n, centr - absc, fval2)
               fval = fval1 + fval2
               res21 = res21 + w21b(k) * fval
               savfun(ipx, :) = fval(:)
               fv3(k, :) = fval1(:)
               fv4(k, :) = fval2(:)
            end do
            resultf = res21 * hlgth
            abserr = maxval(abs((res21 - res10) * hlgth))
         elseif (l == 2) then
            res43 = w43b(12) * fcentr
            neval = 43
            do k = 1, 10
               res43 = res43 + savfun(k, :) * w43a(k)
            end do
            do k = 1, 11
               ipx = ipx + 1
               absc = hlgth * x3(k)
               call f(n, centr + absc, fval1)
               call f(n, centr - absc, fval2)
               fval = fval1 + fval2
               res43 = res43 + fval * w43b(k)
               savfun(ipx, :) = fval(:)
            end do
            resultf = res43 * hlgth
            abserr = maxval(abs((res43 - res21) * hlgth))
         elseif (l == 3) then
            res87 = w87b(23) * fcentr
            neval = 87
            do k = 1, 21
               res87 = res87 + savfun(k, :) * w87a(k)
            end do
            do k = 1, 22
               absc = hlgth * x4(k)
               call f(n, centr + absc, fval1)
               call f(n, centr - absc, fval2)
               res87 = res87 + w87b(k) * (fval1 + fval2)
            end do
            resultf = res87 * hlgth
            abserr = maxval(abs((res87 - res43) * hlgth))
         end if
         if (abserr <= max(epsabs, epsrel * maxval(abs(resultf)))) then
            ier = 0
         end if
         if (ier == 0) then
            exit
         end if
      end do
   end subroutine qng

   recursive subroutine gkintegrate(ntot, t0, t1, qsub, qint, subdiv, &
                                    errorcodes, inteps, mindiv, maxnumdiv)
      implicit none
      integer, intent(in) :: ntot
      integer, intent(inout) :: subdiv, errorcodes
      integer :: nsteps, ier, subdiv1, subdiv2, maxnumdiv, ec1, ec2
      real(8), intent(in) :: t0, t1
      real(8) :: t00, tmid, t11, errstep, inteps, mindiv
      complex(8), intent(out) :: qint(ntot)
      complex(8) :: qint1(ntot), qint2(ntot)
      external :: qsub

      errorcodes = 0
      call qng(ntot, t0, t1, inteps, inteps, qsub, qint, errstep, nsteps, ier)
      if (abs(t1 - t0) .lt. mindiv) then
         errorcodes = 2
!            write(*,'('' min delta: subdiv,t0,t1,err,steps:'',i5,5es12.4,i4)') &
!                subdiv,t0,t1,t1-t0,maxval(abs(qint)),errstep,nsteps
         return
      end if
      if (ier .ne. 0) then
         if (subdiv .ge. maxnumdiv) then
            errorcodes = 1
!               write(*,'('' max sub: subdiv,t0,t1,err,steps:'',i5,5es12.4,i4)') &
!                  subdiv,t0,t1,t1-t0,maxval(abs(qint)),errstep,nsteps
         else
            subdiv = subdiv + 1
            subdiv1 = subdiv
            subdiv2 = subdiv
            t00 = t0
            tmid = (t0 + t1) * 0.5d0
            t11 = t1
            call gkintegrate(ntot, t00, tmid, qsub, qint1, subdiv1, ec1, inteps, mindiv, maxnumdiv)
            call gkintegrate(ntot, tmid, t11, qsub, qint2, subdiv2, ec2, inteps, mindiv, maxnumdiv)
            subdiv = max(subdiv1, subdiv2)
            errorcodes = max(ec1, ec2)
            qint = qint1 + qint2
         end if
      end if
      return
   end subroutine gkintegrate

   subroutine gauleg(x1, x2, x, w, n)
      implicit none
      integer :: n, m, j, i
      real(8) :: x1, x2, x(n), w(n), xm, xl, z, p1, p2, p3, pp, z1, dj
      real(8), parameter :: eps = 3.d-14
      m = (n + 1) / 2
      xm = 0.5d0 * (x2 + x1)
      xl = 0.5d0 * (x2 - x1)
      do i = 1, m
         z = cos(3.141592654d0 * (i - .25d0) / (n + .5d0))
         z1 = z - 1.d0
         do while (abs(z - z1) .gt. eps)
            p1 = 1.d0
            p2 = 0.d0
            do j = 1, n
               dj = dble(j)
               p3 = p2
               p2 = p1
               p1 = ((2.d0 * dj - 1.d0) * z * p2 - (dj - 1.d0) * p3) / dj
            end do
            pp = n * (z * p1 - p2) / (z * z - 1.d0)
            z1 = z
            z = z1 - p1 / pp
         end do
         x(i) = xm - xl * z
         x(n + 1 - i) = xm + xl * z
         w(i) = 2.d0 * xl / ((1.d0 - z * z) * pp * pp)
         w(n + 1 - i) = w(i)
      end do
      return
   end subroutine gauleg

   subroutine realsort(nlimits0, limits, eps, nlimits)
      implicit none
      integer :: nlimits0, nlimits, imin(1), n
      real(8) :: limits(1:nlimits0), rtemp(1:nlimits0), eps
      rtemp(1:nlimits0) = limits(1:nlimits0)
      imin = minloc(rtemp(1:nlimits0))
      nlimits = 1
      limits(nlimits) = rtemp(imin(1))
      rtemp(imin(1)) = 1.d10
      do n = 2, nlimits0
         imin = minloc(rtemp(1:nlimits0))
         if (rtemp(imin(1)) - limits(nlimits) .gt. eps) then
            nlimits = nlimits + 1
            limits(nlimits) = rtemp(imin(1))
         end if
         rtemp(imin(1)) = 1.d10
      end do
   end subroutine realsort

end module quadrature_functions
