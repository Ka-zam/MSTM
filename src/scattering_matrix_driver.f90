module scattering_matrix_driver
   use input_state
   use scattering, only: common_origin_scattering_matrix, evaluate_fixed_orientation_scattering_matrix, &
                         multiple_origin_scattering_matrix, &
                         numerical_scattering_matrix_azimuthal_average_multiple_origin, &
                         numerical_scattering_matrix_azimuthal_average_single_origin, periodic_lattice_scattering
   implicit none
   private
   public :: compute_scattering_matrix
contains

   subroutine compute_scattering_matrix(amnp, scatmat, mpi_comm)
      implicit none
      logical :: singleorigin, iframe
      integer :: i, sy, sx, mpicomm
      integer, optional :: mpi_comm
      real(8) :: scatmat(scat_mat_mdim, scat_mat_ldim:scat_mat_udim), costheta, phi, csca, &
                 ky, kx, sintheta, ctm
      complex(8) :: amnp(*), ampmat(2, 2)
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      singleorigin = number_plane_boundaries .eq. 0 .and. single_origin_expansion
      iframe = singleorigin .and. incident_frame
      csca = (q_eff_tot(1, 1) - q_eff_tot(2, 1)) * pi * cross_section_radius**2
      csca = two_pi

      if (periodic_lattice) then
         call periodic_lattice_scattering(amnp, pl_sca, scat_mat=scatmat, krho_vec=rl_vec)
         return
      end if

      if (scattering_map_model .eq. 0) then
         do i = scat_mat_ldim, scat_mat_udim
            costheta = cosd(scat_mat_amin + (scat_mat_amax - scat_mat_amin) &
                            * dble(i - scat_mat_ldim) / dble(scat_mat_udim - scat_mat_ldim))
            if (costheta .eq. 1.d0) costheta = 0.9999999d0
            if (costheta .eq. -1.d0) costheta = -0.9999999d0
            if (i .lt. 0) then
               phi = incident_alpha_deg * degrees_to_radians + pi
            else
               phi = incident_alpha_deg * degrees_to_radians
            end if
            if (number_plane_boundaries .eq. 0) then
               if (singleorigin) then
                  if (azimuthal_average) then
                     if (.not. numerical_azimuthal_average) then
!                        call evaluate_fixed_orientation_scattering_matrix(12,s00,s02,sp22,sm22,costheta,scatmat(:,i),normalize_s11=.false.)
                        call evaluate_fixed_orientation_scattering_matrix( &
                           t_matrix_order, scat_mat_exp_coef(:, :, 1), scat_mat_exp_coef(:, :, 2), &
                           scat_mat_exp_coef(:, :, 3), scat_mat_exp_coef(:, :, 4), &
                           costheta, scatmat(:, i), normalize_s11=.false.)
                     else
                        call numerical_scattering_matrix_azimuthal_average_single_origin( &
                           amnp, t_matrix_order, costheta, scatmat(:, i), &
                           rotate_plane=.true., normalize_s11=.false.)
                     end if
                  else
                     call common_origin_scattering_matrix(amnp, t_matrix_order, costheta, phi, ampmat, scatmat(:, i), &
                                                          rotate_plane=iframe, normalize_s11=.false.)
                  end if
               else
                  if (azimuthal_average) then
                     call numerical_scattering_matrix_azimuthal_average_multiple_origin( &
                        amnp, costheta, scatmat(:, i), rotate_plane=.true.)
                  else
                     call multiple_origin_scattering_matrix(amnp, costheta, phi, csca, ampmat, scatmat(:, i), &
                                                            rotate_plane=.true.)
                  end if
               end if
            else
               ctm = -costheta
               if (azimuthal_average) then
                  call numerical_scattering_matrix_azimuthal_average_multiple_origin(amnp, ctm, scatmat(1:16, i))
                  call numerical_scattering_matrix_azimuthal_average_multiple_origin( &
                     amnp, costheta, scatmat(17:32, i))
               else
                  call multiple_origin_scattering_matrix(amnp, ctm, phi, csca, ampmat, scatmat(1:16, i))
                  call multiple_origin_scattering_matrix(amnp, costheta, phi, csca, ampmat, scatmat(17:32, i))
               end if
            end if
         end do
      else
         i = 0
         do sy = -scattering_map_dimension, scattering_map_dimension
            ky = dble(sy) / dble(scattering_map_dimension)
            do sx = -scattering_map_dimension, scattering_map_dimension
               if (sx * sx + sy * sy .gt. scattering_map_dimension**2) cycle
               kx = dble(sx) / dble(scattering_map_dimension)
               sintheta = kx * kx + ky * ky
               sintheta = min(sintheta, .99999d0)
               i = i + 1
               if (sx .eq. 0 .and. sy .eq. 0) then
                  phi = 0.d0
               else
                  phi = datan2(ky, kx)
               end if
               costheta = -sqrt(1.d0 - sintheta)
               if (singleorigin) then
                  call common_origin_scattering_matrix(amnp, t_matrix_order, costheta, phi, ampmat, &
                                                       scatmat(1:16, i), &
                                                       rotate_plane=iframe, normalize_s11=.false.)
               else
                  call multiple_origin_scattering_matrix(amnp, costheta, phi, csca, ampmat, scatmat(1:16, i), &
                                                         rotate_plane=incident_frame)
               end if
               costheta = -costheta
               if (singleorigin) then
                  call common_origin_scattering_matrix(amnp, t_matrix_order, costheta, phi, ampmat, &
                                                       scatmat(17:32, i), &
                                                       rotate_plane=iframe, normalize_s11=.false.)
               else
                  call multiple_origin_scattering_matrix(amnp, costheta, phi, csca, ampmat, &
                                                         scatmat(17:32, i), rotate_plane=incident_frame)
               end if
            end do
         end do
      end if
   end subroutine compute_scattering_matrix
end module scattering_matrix_driver
