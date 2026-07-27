module scattering_matrix_driver
   use, intrinsic :: iso_fortran_env, only: real64
   use input_state
   use parallel_runtime, only: mpi_comm_world
   use scattering_amplitudes, only: common_origin_scattering_matrix, &
                                    evaluate_fixed_orientation_scattering_matrix, multiple_origin_scattering_matrix, &
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
      real(real64) :: scatmat(simulation_result%scattering_matrix_dimension, simulation_result%scattering_matrix_lower_bound:simulation_result%scattering_matrix_upper_bound), costheta, phi, csca, &
                      ky, kx, sintheta, ctm
      complex(real64) :: amnp(*), ampmat(2, 2)
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      singleorigin = number_plane_boundaries .eq. 0 .and. simulation_config%single_origin_expansion
      iframe = singleorigin .and. simulation_config%incident_frame
      csca = (simulation_result%total_efficiency(1, 1) - simulation_result%total_efficiency(2, 1)) * pi * sphere_cluster%cross_section_radius**2
      csca = two_pi

      if (periodic_lattice) then
         call periodic_lattice_scattering(amnp, simulation_result%plane_scattering, scat_mat=scatmat, &
                                          krho_vec=simulation_result%reflection_transmission_vectors)
         return
      end if

      if (simulation_config%scattering_map_model .eq. 0) then
         do i = simulation_result%scattering_matrix_lower_bound, simulation_result%scattering_matrix_upper_bound
            costheta = cosd(simulation_config%scattering_matrix_angle_minimum + (simulation_config%scattering_matrix_angle_maximum - simulation_config%scattering_matrix_angle_minimum) &
                            * dble(i - simulation_result%scattering_matrix_lower_bound) / dble(simulation_result%scattering_matrix_upper_bound - simulation_result%scattering_matrix_lower_bound))
            if (costheta .eq. 1.d0) costheta = 0.9999999d0
            if (costheta .eq. -1.d0) costheta = -0.9999999d0
            if (i .lt. 0) then
               phi = simulation_config%incident_alpha_degrees * degrees_to_radians + pi
            else
               phi = simulation_config%incident_alpha_degrees * degrees_to_radians
            end if
            if (number_plane_boundaries .eq. 0) then
               if (singleorigin) then
                  if (simulation_config%azimuthal_average) then
                     if (.not. simulation_config%numerical_azimuthal_average) then
!                        call evaluate_fixed_orientation_scattering_matrix(12,s00,s02,sp22,sm22,costheta,scatmat(:,i),simulation_config%output%normalize_s11=.false.)
                        call evaluate_fixed_orientation_scattering_matrix( &
                           sphere_cluster%t_matrix_order, simulation_result%scattering_matrix_expansion(:, :, 1), simulation_result%scattering_matrix_expansion(:, :, 2), &
                   simulation_result%scattering_matrix_expansion(:, :, 3), simulation_result%scattering_matrix_expansion(:, :, 4), &
                           costheta, scatmat(:, i), normalize_s11=.false.)
                     else
                        call numerical_scattering_matrix_azimuthal_average_single_origin( &
                           amnp, sphere_cluster%t_matrix_order, costheta, scatmat(:, i), &
                           rotate_plane=.true., normalize_s11=.false.)
                     end if
                  else
                   call common_origin_scattering_matrix(amnp, sphere_cluster%t_matrix_order, costheta, phi, ampmat, scatmat(:, i), &
                                                          rotate_plane=iframe, normalize_s11=.false.)
                  end if
               else
                  if (simulation_config%azimuthal_average) then
                     call numerical_scattering_matrix_azimuthal_average_multiple_origin( &
                        amnp, costheta, scatmat(:, i), rotate_plane=.true.)
                  else
                     call multiple_origin_scattering_matrix(amnp, costheta, phi, csca, ampmat, scatmat(:, i), &
                                                            rotate_plane=.true.)
                  end if
               end if
            else
               ctm = -costheta
               if (simulation_config%azimuthal_average) then
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
         do sy = -simulation_config%scattering_map_dimension, simulation_config%scattering_map_dimension
            ky = dble(sy) / dble(simulation_config%scattering_map_dimension)
            do sx = -simulation_config%scattering_map_dimension, simulation_config%scattering_map_dimension
               if (sx * sx + sy * sy .gt. simulation_config%scattering_map_dimension**2) cycle
               kx = dble(sx) / dble(simulation_config%scattering_map_dimension)
               sintheta = kx * kx + ky * ky
               sintheta = min(sintheta, .99999d0)
               i = i + 1
               if (sx .eq. 0 .and. sy .eq. 0) then
                  phi = 0.d0
               else
                  phi = atan2(ky, kx)
               end if
               costheta = -sqrt(1.d0 - sintheta)
               if (singleorigin) then
                  call common_origin_scattering_matrix(amnp, sphere_cluster%t_matrix_order, costheta, phi, ampmat, &
                                                       scatmat(1:16, i), &
                                                       rotate_plane=iframe, normalize_s11=.false.)
               else
                  call multiple_origin_scattering_matrix(amnp, costheta, phi, csca, ampmat, scatmat(1:16, i), &
                                                         rotate_plane=simulation_config%incident_frame)
               end if
               costheta = -costheta
               if (singleorigin) then
                  call common_origin_scattering_matrix(amnp, sphere_cluster%t_matrix_order, costheta, phi, ampmat, &
                                                       scatmat(17:32, i), &
                                                       rotate_plane=iframe, normalize_s11=.false.)
               else
                  call multiple_origin_scattering_matrix(amnp, costheta, phi, csca, ampmat, &
                                                         scatmat(17:32, i), rotate_plane=simulation_config%incident_frame)
               end if
            end do
         end do
      end if
   end subroutine compute_scattering_matrix
end module scattering_matrix_driver
