module translation_operator
   use angular_functions, only: atcdim, axialtrancoefrecurrence, cartosphere, ephicoef, gentranmatrix, moffset, rotcoef
   use iso_fortran_env, only: real64
   use wave_functions, only: mtransfer

   implicit none(type, external)
   private

   public :: shiftcoefficient, translation_data

   type translation_data
      private
      logical :: matrix_calculated = .false.
      logical :: rot_op = .false.
      logical :: zero_translation = .false.
      integer :: vswf_type = 0
      real(real64) :: translation_vector(3) = 0.0_real64
      real(real64), allocatable :: rot_mat(:, :)
      complex(real64) :: refractive_index(2) = (1.0_real64, 0.0_real64)
      complex(real64), allocatable :: phi_mat(:), z_mat(:), gen_mat(:, :, :)
   contains
      procedure, public :: configure => configure_translation
      procedure, public, pass(tranmat) :: apply => coefficient_translation
      procedure, public :: clear => clear_translation
      final :: finalize_translation
   end type translation_data

contains

   subroutine configure_translation(self, vswf_type, translation_vector, refractive_index, use_rotation)
      class(translation_data), intent(inout) :: self
      integer, intent(in) :: vswf_type
      real(real64), intent(in) :: translation_vector(3)
      complex(real64), intent(in) :: refractive_index(2)
      logical, intent(in) :: use_rotation

      call self%clear()
      self%vswf_type = vswf_type
      self%translation_vector = translation_vector
      self%refractive_index = refractive_index
      self%rot_op = use_rotation
   end subroutine configure_translation

   subroutine clear_translation(self)
      class(translation_data), intent(inout) :: self

      if (allocated(self%rot_mat)) deallocate (self%rot_mat)
      if (allocated(self%phi_mat)) deallocate (self%phi_mat)
      if (allocated(self%z_mat)) deallocate (self%z_mat)
      if (allocated(self%gen_mat)) deallocate (self%gen_mat)
      self%matrix_calculated = .false.
      self%zero_translation = .false.
      self%vswf_type = 0
      self%translation_vector = 0.0_real64
      self%refractive_index = (1.0_real64, 0.0_real64)
      self%rot_op = .false.
   end subroutine clear_translation

   subroutine finalize_translation(self)
      type(translation_data), intent(inout) :: self

      call self%clear()
   end subroutine finalize_translation

   subroutine coefficient_translation(nodra, nmodea, nodrg, nmodeg, &
                                      acoef, gcoef, tranmat, shift_op, tran_op)
      implicit none(type, external)
      logical :: sop, top, rot
      logical, optional, intent(in) :: shift_op, tran_op
      integer, intent(in) :: nodra, nodrg, nmodea, nmodeg
      integer :: nblka, nblkg, lengtha, lengthg, nmode, shiftvec(2), n, m, im, nn1, nn2, n1, &
                 p, nmin, m1, offset, blocksize, nmax, tdim, vtype, nodrs, nodrt
      real(real64) :: r, ct, rtran(3)
      complex(real64), intent(in) :: acoef(*)
      complex(real64), intent(inout) :: gcoef(*)
      complex(real64) :: &
         a_t(0:nodra + 1, nodra, nmodea), g_t(0:nodrg + 1, nodrg, nmodeg), &
         g_shift(0:nodrg + 1, nodrg, nmodeg), &
         a_tt(-nodra:nodra, nodra, 2), g_tt(-nodrg:nodrg, nodrg, 2), &
         atc(max(nodra, nodrg), max(nodra, nodrg), 2), &
         a_t2(nodra * (nodra + 2), 2), g_t2(nodrg * (nodrg + 2), 2), rimed(2), ephi
      class(translation_data), intent(inout) :: tranmat

      nblka = nodra * (nodra + 2)
      nblkg = nodrg * (nodrg + 2)
      nmin = min(nodra, nodrg)
      nmax = max(nodra, nodrg)
      nmode = max(nmodea, nmodeg)
      if (present(shift_op)) then
         sop = shift_op
      else
         sop = .false.
      end if
      if (present(tran_op)) then
         top = tran_op
      else
         top = .false.
      end if
      rot = tranmat%rot_op

      if (.not. tranmat%matrix_calculated) then
         vtype = tranmat%vswf_type
         rimed = tranmat%refractive_index
         rtran = tranmat%translation_vector
         r = dot_product(rtran, rtran)
         if (r .lt. 1.0e-12_real64) then
            tranmat%zero_translation = .true.
         else
            tranmat%zero_translation = .false.
         end if
         if (sop) then
            nodrs = nodrg
            nodrt = nodra
         else
            nodrs = nodra
            nodrt = nodrg
         end if
         if (.not. tranmat%zero_translation) then
            if (rot) then
               tdim = atcdim(nodrt, nodrs)
               allocate (tranmat%rot_mat(-nmin:nmin, 0:nmax * (nmax + 2)))
               allocate (tranmat%phi_mat(-nmax:nmax))
               allocate (tranmat%z_mat(1:tdim))
               call cartosphere(rtran, r, ct, ephi)
               call rotcoef(ct, nmin, nmax, tranmat%rot_mat)
               call axialtrancoefrecurrence(vtype, r, rimed, nodrt, nodrs, &
                                            tdim, tranmat%z_mat)
               call ephicoef(ephi, nmax, tranmat%phi_mat)
            else
               allocate (tranmat%gen_mat(nodrt * (nodrt + 2), nodrs * (nodrs + 2), 2))
               call gentranmatrix(nodrs, nodrt, translation_vector=rtran, &
                                  refractive_index=rimed, ac_matrix=tranmat%gen_mat, vswf_type=vtype, &
                                  mode_s=2, mode_t=2)
            end if
         end if
         tranmat%matrix_calculated = .true.
      end if
      lengtha = nblka * nmodea
      lengthg = nblkg * nmodeg
      shiftvec = (/1, 1/)
      if (sop .neqv. top) then
         shiftvec = -shiftvec
      end if
      if (sop) then
         im = -1
      else
         im = 1
      end if
      if (tranmat%zero_translation) then
         if (tranmat%vswf_type .eq. 1) then
            call mtransfer(nodra, nodrg, acoef(1:lengtha), g_t)
            gcoef(1:lengthg) = gcoef(1:lengthg) &
                               + reshape(g_t(0:nodrg + 1, 1:nodrg, 1:nmodeg), (/lengthg/))
         end if
      else
         if (rot) then
            call shiftcoefficient(nodra, nmodea, shiftvec(1), shiftvec(2), &
                                  acoef(1:lengtha), a_t(0:nodra + 1, 1:nodra, 1:2))
            a_tt(0, 1:nodra, 1:2) = a_t(0, 1:nodra, 1:2)
            do m = 1, nodra
               a_tt(m, m:nodra, 1:2) = a_t(m, m:nodra, 1:2) * tranmat%phi_mat(im * m)
               a_tt(-m, m:nodra, 1:2) = a_t(m + 1:nodra + 1, m, 1:2) * tranmat%phi_mat(-im * m)
            end do
            do n = 1, nodra
               nn1 = n * (n + 1) - n
               nn2 = nn1 + (2 * n + 1) - 1
               n1 = min(n, nodrg)
               a_tt(-n1:n1, n, 1:2) &
                  = matmul(tranmat%rot_mat(-n1:n1, nn1:nn2), a_tt(-n:n, n, 1:2))
            end do
            do m = -nmin, nmin
               m1 = max(1, abs(m))
               if (sop) then
                  offset = moffset(m, nodra, nodrg)
                  blocksize = (nodrg - m1 + 1) * (nodra - m1 + 1) * 2
                  atc(m1:nodra, m1:nodrg, 1:2) = &
                     reshape(tranmat%z_mat(offset + 1:offset + blocksize), &
                             (/nodra - m1 + 1, nodrg - m1 + 1, 2/))
                  do p = 1, 2
                     g_tt(m, m1:nodrg, p) &
                        = matmul(a_tt(m, m1:nodra, p), atc(m1:nodra, m1:nodrg, p))
                  end do
               else
                  offset = moffset(m, nodrg, nodra)
                  blocksize = (nodrg - m1 + 1) * (nodra - m1 + 1) * 2
                  atc(m1:nodrg, m1:nodra, 1:2) = &
                     reshape(tranmat%z_mat(offset + 1:offset + blocksize), &
                             (/nodrg - m1 + 1, nodra - m1 + 1, 2/))
                  do p = 1, 2
                     g_tt(m, m1:nodrg, p) &
                        = matmul(atc(m1:nodrg, m1:nodra, p), a_tt(m, m1:nodra, p))
                  end do
               end if
            end do
            do n = 1, nodrg
               nn1 = n * (n + 1) - n
               nn2 = nn1 + (2 * n + 1) - 1
               n1 = min(n, nodra)
               g_tt(-n:n, n, 1) = matmul(g_tt(-n1:n1, n, 1), tranmat%rot_mat(-n1:n1, nn1:nn2))
               g_tt(-n:n, n, 2) = matmul(g_tt(-n1:n1, n, 2), tranmat%rot_mat(-n1:n1, nn1:nn2))
            end do
            g_t(0, 1:nodrg, 1:2) = g_tt(0, 1:nodrg, 1:2)
            do m = 1, nodrg
               g_t(m, m:nodrg, 1:2) = g_tt(m, m:nodrg, 1:2) * tranmat%phi_mat(-im * m)
               g_t(m + 1:nodrg + 1, m, 1:2) = g_tt(-m, m:nodrg, 1:2) * tranmat%phi_mat(im * m)
            end do
            call shiftcoefficient(nodrg, nmodeg, shiftvec(1), shiftvec(2), &
                                  g_t, g_shift)
            g_t = g_shift
         else
            call shiftcoefficient(nodra, nmodea, shiftvec(1), shiftvec(2), &
                                  acoef(1:lengtha), &
                                  a_t2(1:nblka, 1:2))
            if (sop) then
               g_t2(:, 1) = matmul(a_t2(1:nblka, 1), tranmat%gen_mat(1:nblka, 1:nblkg, 1))
            else
               g_t2(:, 1) = matmul(tranmat%gen_mat(1:nblkg, 1:nblka, 1), a_t2(1:nblka, 1))
            end if
            if (nmodeg .eq. 2 .and. nmodea .eq. 1) then
               if (sop) then
                  g_t2(:, 2) = matmul(a_t2(1:nblka, 1), tranmat%gen_mat(:, :, 2))
                  g_t2 = 0.5_real64 * g_t2
               else
                  g_t2(:, 2) = matmul(tranmat%gen_mat(:, :, 2), a_t2(1:nblka, 1))
                  g_t2 = 0.5_real64 * g_t2
               end if
            elseif (nmodeg .eq. 1 .and. nmodea .eq. 2) then
               if (sop) then
                  g_t2(:, 1) = g_t2(:, 1) + matmul(a_t2(1:nblka, 2), tranmat%gen_mat(:, :, 2))
               else
                  g_t2(:, 1) = g_t2(:, 1) + matmul(tranmat%gen_mat(:, :, 2), a_t2(1:nblka, 2))
               end if
            elseif (nmodeg .eq. 2 .and. nmodea .eq. 2) then
               if (sop) then
                  g_t2(:, 2) = matmul(a_t2(1:nblka, 2), tranmat%gen_mat(1:nblka, 1:nblkg, 2))
               else
                  g_t2(:, 2) = matmul(tranmat%gen_mat(1:nblkg, 1:nblka, 2), a_t2(1:nblka, 2))
               end if
            end if
            call shiftcoefficient(nodrg, nmodeg, shiftvec(1), shiftvec(2), &
                                  g_t2, g_t)
         end if
         gcoef(1:lengthg) &
            = gcoef(1:lengthg) &
              + reshape(g_t(0:nodrg + 1, 1:nodrg, 1:nmodeg), (/lengthg/))
      end if
   end subroutine coefficient_translation

   subroutine shiftcoefficient(nodr, nmode, msign, mflip, &
                               ain, aout)
      implicit none(type, external)
      integer, intent(in) :: nodr, nmode, msign, mflip
      integer :: m, n, im
      complex(real64), intent(in) :: ain(0:nodr + 1, nodr, nmode)
      complex(real64), intent(out) :: aout(0:nodr + 1, nodr, nmode)
      complex(real64) :: at(nmode)

      if (msign .eq. 1 .and. mflip .eq. 1) then
         aout = ain
      else
         aout(0, 1:nodr, 1:nmode) = ain(0, 1:nodr, 1:nmode)
         if (mflip .eq. -1) then
            im = 1
            do m = 1, nodr
               im = im * msign
               do n = m, nodr
                  at = ain(n + 1, m, 1:nmode)
                  aout(n + 1, m, 1:nmode) = im * ain(m, n, 1:nmode)
                  aout(m, n, 1:nmode) = im * at
               end do
            end do
         else
            im = 1
            do m = 1, nodr
               im = im * msign
               do n = m, nodr
                  aout(m, n, 1:nmode) = im * ain(m, n, 1:nmode)
                  aout(n + 1, m, 1:nmode) = im * ain(n + 1, m, 1:nmode)
               end do
            end do
         end if
      end if
   end subroutine shiftcoefficient

end module translation_operator
