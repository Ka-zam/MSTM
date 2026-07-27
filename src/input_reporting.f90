module input_reporting
   use, intrinsic :: iso_fortran_env, only: real64
   use runtime_support, only: open_update_file, runtime_failed, set_runtime_error
   use constants
   use effective_medium_analysis, only: diffuse_scattering_effective_refractive_index, &
                                        effective_extinction_coefficient_ratio
   use fft_translation, only: fft_plan, fft_translation_metrics_t
   use input_state
   implicit none
contains

   subroutine output_header(iunit, inputfile)
      implicit none
      integer :: iunit
      character(len=8) :: rundate
      character(len=10) :: runtime
      character(len=256) :: inputfile
      call date_and_time(date=rundate, time=runtime)
      run_date_and_time = trim(rundate)//' '//trim(runtime)
      write (iunit, '(''****************************************************'')')
      write (iunit, '(''****************************************************'')')
      write (iunit, '('' mstm calculation results'')')
      write (iunit, '('' date, time:'')')
      write (iunit, '(a)') run_date_and_time
      write (iunit, '('' input file:'')')
      write (iunit, '(a)') trim(inputfile)
   end subroutine output_header

   subroutine print_run_variables(iunit)
      implicit none
      integer :: iunit, i, n, fft_cell_dimensions(3), fft_neighbor_count, fft_order
      real(real64) :: cb, r(2), t(2), a(2), tvol, svol, fft_volume_fraction, fft_cell_size

      write (iunit, '(''****************************************************'')')
      write (iunit, '('' input variables for run '',i5)') run_number
      if (random_configuration) then
         write (iunit, '('' sphere positions randomly generated'')')
         if (target_shape .eq. 0) then
            write (iunit, '('' rectangular target, half-widths in x,y,z'')')
            write (iunit, '(3es12.4)') target_dimensions * length_scale_factor
         elseif (target_shape .eq. 1) then
            write (iunit, '('' cylindrical target, radius, half-thickness'')')
            write (iunit, '(3es12.4)') target_dimensions(1) * length_scale_factor, &
               target_dimensions(3) * length_scale_factor
         else
            write (iunit, '('' spherical target, radius'')')
            write (iunit, '(3es12.4)') target_dimensions(1) * length_scale_factor
            if (random_configuration_host) then
               write (iunit, '('' target enclosed in host sphere w/ radius, ref index'')')
               write (iunit, '(3es12.4)') sphere_radius(number_spheres), &
                  sphere_ref_index(1, number_spheres)
            end if
         end if
         write (iunit, '('' number of components:'')')
         write (iunit, '(i3)') number_components

         write (iunit, '('' sphere log-normal PSD sigma:'')')
         write (iunit, '(4es12.4)') psd_sigma(1:number_components)

         if (number_components .gt. 1) then
            write (iunit, '('' sphere component radii:'')')
            write (iunit, '(4es12.4)') component_radii(1:number_components)
            write (iunit, '('' sphere component number_fraction:'')')
            write (iunit, '(4es12.4)') component_number_fraction(1:number_components)
            write (iunit, '('' sphere component refractive index:'')')
            write (iunit, '(8es12.4)') component_ref_index(1:number_components)
         end if

         write (iunit, '('' sphere volume fraction'')')
         if ((.not. number_spheres_specified) .or. auto_target_radius) then
            write (iunit, '(es12.4)') sphere_volume_fraction
         else
            call calculate_target_volume(target_dimensions, tvol)
            tvol = tvol * length_scale_factor**3
            svol = four_pi_over_three * sum(sphere_radius(:)**3)
            write (iunit, '(es12.4)') svol / tvol
         end if
         if (ran_config_stat .eq. 0) then
            write (iunit, '('' target configuration computed using random sampling + diffusion'')')
         elseif (ran_config_stat .eq. 1) then
            write (iunit, '('' target configuration computed using layered sampling + diffusion'')')
         elseif (ran_config_stat .eq. 2) then
            write (iunit, '('' target configuration computed initial HCP + diffusion'')')
         end if
         write (iunit, '('' number diffusion time steps:'',i5)') ran_config_time_steps
      else
         write (iunit, '('' sphere data input file:'')')
         write (iunit, '(a)') trim(sphere_data_input_file)
      end if
      write (iunit, '('' number spheres'')')
      write (iunit, '(i7)') number_spheres
      write (iunit, '('' length, ref index scale factors'')')
      write (iunit, '(3es15.7)') length_scale_factor, ref_index_scale_factor
      write (iunit, '('' volume cluster radius, area mean sphere radius, circumscribing radius, cross section radius'')')
      write (iunit, '(4es15.7)') vol_radius, area_mean_radius, circumscribing_radius, cross_section_radius
      if (print_sphere_data) then
         write (iunit, '('' sphere properties and associations'')')
         if (any_optically_active) then
            write (iunit, '(''   sphere    host   layer radius     x       y       z    '', &
              &''     ref indx(L)             ref_indx(R)'')')
         else
            write (iunit, '(''   sphere    host   layer radius     x       y       z           ref indx'')')
         end if
         do n = 1, number_spheres
            if (any_optically_active) then
               write (iunit, '(3i8,4f8.3,4es12.4)') n, host_sphere(n), sphere_layer(n), sphere_radius(n), &
                  sphere_position(:, n), sphere_ref_index(1:2, n)
            else
               write (iunit, '(3i8,4f8.3,4es12.4)') n, host_sphere(n), sphere_layer(n), sphere_radius(n), &
                  sphere_position(:, n), sphere_ref_index(1, n)
            end if
         end do
      end if

      if (random_orientation) then
         write (iunit, '('' random orientation, estimated t matrix order:'')')
         write (iunit, '(i6)') t_matrix_order
      else
         if (gaussian_beam_constant .ne. 0.d0) then
            write (iunit, '('' incident Gaussian beam: 1/beam width, focal point'')')
            write (iunit, '(4es12.4)') gaussian_beam_constant, gaussian_beam_focal_point
         else
            write (iunit, '('' incident plane wave'')')
         end if
         if (incidence_average) then
            write (iunit, '('' Monte Carlo average over incident directions'')')
         else
            if (incident_beta_specified) then
               write (iunit, '('' incident alpha, beta(deg)'')')
               write (iunit, '(2es12.4)') incident_alpha_deg, incident_beta_deg
            else
               write (iunit, '('' incident alpha(deg), incident sin(beta), incident direction'')')
               write (iunit, '(2es12.4,i3)') incident_alpha_deg, incident_sin_beta, 3 - 2 * incident_direction
            end if
         end if
         if (single_origin_expansion) then
            write (iunit, '('' t matrix order:'')')
            write (iunit, '(i6)') t_matrix_order
         end if
      end if
      if (reflection_model) then
         if (auto_absorption_sample_radius) then
            write (iunit, '('' particle layer reflectance/absorptance model, absorption sample radius fraction'')')
            write (iunit, '(es12.4)') absorption_sample_radius_fraction
         else
            write (iunit, '('' particle layer reflectance/absorptance model, absorption sample radius'')')
            write (iunit, '(es12.4)') absorption_sample_radius * length_scale_factor
         end if
      end if
      write (iunit, '('' layer 0 refractive index'')')
      write (iunit, '(4es12.4)') layer_ref_index(0)
      write (iunit, '('' number of plane boundaries '')')
      write (iunit, '(i3)') number_plane_boundaries
      if (number_plane_boundaries .gt. 0) then
         write (iunit, '('' boundary, position, refractive index'')')
         do i = 1, number_plane_boundaries
            write (iunit, '(i3,3es12.4)') i, plane_boundary_position(i), layer_ref_index(i)
         end do
         if (.not. incidence_average) then
            cb = cosd(incident_beta_deg)
            call boundary_energy_transfer(incident_sin_beta, incident_direction, r, t, a)
            write (iunit, '('' Fresnel boundary reflectance, transmittance, absorptance (par, perp)'')')
            write (iunit, '(3es12.4)') r(1), t(1), a(1)
            write (iunit, '(3es12.4)') r(2), t(2), a(2)
         end if
         if (number_singular_points .gt. 0) then
            write (iunit, '('' GF singular points (in s)'')')
            do i = 1, number_singular_points
               write (iunit, '(2i5,es12.4,es20.10)') i, singular_point_polarization(i), &
                  singular_gf_value(i), singular_points(i)
            end do
         end if
      end if
      if (periodic_lattice) then
         write (iunit, '('' periodic lattice cell width, incident lateral vector '')')
         write (iunit, '(4es15.7)') cell_width, incident_lateral_vector
      end if
      write (iunit, '('' max_iterations,solution_epsilon, mie_epsilon, interaction radius'')')
      write (iunit, '(i10,3es12.4)') max_iterations, solution_epsilon, mie_epsilon, interaction_radius
      write (iunit, '('' maximum Mie order, number of equations:'')')
      write (iunit, '(i4,i10)') max_mie_order, number_eqns
      write (iunit, '('' mean sphere Mie extinction, absorption efficiencies, albedo'')')
      write (iunit, '(3es12.4)') mean_qext_mie, mean_qabs_mie, 1.d0 - mean_qabs_mie / mean_qext_mie
      if (fft_translation_option) then
         call fft_plan%configuration(fft_cell_dimensions, fft_volume_fraction, fft_cell_size, &
                                     fft_neighbor_count, fft_order)
         write (iunit, '('' fft translation option implemented'')')
         write (iunit, '('' cell width, cell volume fraction, cell dimension:'')')
         write (iunit, '(3i8,2es12.4)') fft_cell_dimensions, fft_volume_fraction, fft_cell_size
         write (iunit, '('' number of neighbor nodes, node order:'')')
         write (iunit, '(2i5)') fft_neighbor_count, fft_order
      end if

      write (iunit, '(''****************************************************'')')
      write (iunit, '('' calculation results for run '')')
      write (iunit, *) run_number
      flush (iunit)
   end subroutine print_run_variables

   subroutine print_error_codes(outunit)
      implicit none
      integer :: outunit
      if (number_plane_boundaries .eq. 0) then
       if (maxval(error_codes) .ne. 0) write (outunit, '('' warning: problems encountered with surface interaction calculations'')')
         if (error_codes(1) .ne. 0) write (outunit, '('' iterative surface GF algorithm did not converge'')')
         if (error_codes(2) .ne. 0) then
            write (outunit, '('' integration along real s axis did not converge'')')
            write (outunit, '('' increase real_axis_integration_limit from current value of '',es10.2)') &
               real_axis_integration_limit
            write (outunit, '('' or increase integration_limit_epsilon from current value of '',es10.2, &
            &  '' and see if results change'')') integration_limit_epsilon
         end if
         if (error_codes(3) .ne. 0) write (outunit, '('' subdivided integration interval below 1d-12'')')
         if (error_codes(4) .ne. 0) then
            write (outunit, '('' maximum subdivision reached in GK integration algorithm'')')
            write (outunit, '('' increase maximum_integration_subdivisions from current value of '', i5, &
            &'' and see if results change'')') maximum_integration_subdivisions
         end if
      end if
      if (periodic_lattice) then
       if (maxval(pl_error_codes) .ne. 0) write (outunit, '('' warning: problems encountered with periodic lattice calculations'')')
         if (pl_error_codes(1) .ne. 0) write (outunit, '('' integration formulas for FS PL DMGF did not converge''&
         &'' in subroutine scalar_wave_function_lattice_sum'')')
         if (pl_error_codes(2) .ne. 0) write (outunit, '('' RS series for PL DMGF did not converge in subroutine''&
         &'' reciprocal_space_scalar_wave_function_lattice_sum'')')
         if (pl_error_codes(3) .ne. 0) write (outunit, '('' reciprocal space series for PL DMGF did not converge''&
             &'' in subroutine plane_boundary_lattice_interaction'')')
         if (pl_error_codes(4) .ne. 0) &
            write (outunit, '('' integration did not converge in subroutine integrate_lattice_term_2d'')')
         if (pl_error_codes(5) .ne. 0) &
            write (outunit, '('' integration did not converge in subroutine integrate_lattice_term_1d_without_source'')')
         if (pl_error_codes(6) .ne. 0) &
            write (outunit, '('' series in s did not converge in subroutine scalar_wave_function_yz_lattice_sum'')')
      end if
   end subroutine print_error_codes

   subroutine print_calculation_results(fout)
      implicit none
      integer :: io_status, outunit, n, i, j, smvec(6), sx, sy, s, smvec0(16), nsmat, smvecp(16)
      real(real64) :: smt(16), kx, ky, s11scale, r(2), t(2), a(2), scacoef, scarat, &
                      abscoef, absrat, rl(2), al(2), tl(2), tvol, imrieff, qeeff
      character(len=2) :: smlabel(16)
      character(len=256) :: fout, chartemp
      type(fft_translation_metrics_t) :: fft_metrics
      smvec = (/1, 5, 6, 11, 15, 16/)
      smvec0 = (/1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16/)
      smlabel = (/'11', '21', '31', '41', '12', '22', '32', '42', '13', '23', '33', '43', '14', '24', '34', '44'/)
      if (fout(1:7) .eq. 'console') then
         outunit = 6
      else
         call open_update_file(trim(fout), outunit)
         if (runtime_failed()) return
         do
            read (outunit, '(a)', iostat=io_status) chartemp
            if (io_status /= 0) then
               call set_runtime_error('Cannot find current run in output file: '//trim(fout))
               close (outunit)
               return
            end if
            if (trim(chartemp) .eq. run_date_and_time) exit
         end do
         do
            read (outunit, '(a)', iostat=io_status) chartemp
            if (io_status /= 0) then
               call set_runtime_error('Cannot find calculation section in output file: '//trim(fout))
               close (outunit)
               return
            end if
            if (chartemp(1:28) .eq. ' calculation results for run') then
               read (outunit, *) i
               if (i .eq. run_number) exit
            end if
         end do
      end if
      if (configuration_average) then
         write (outunit, '('' averages collected: sample:'')')
         write (outunit, '(i5)') random_configuration_number * n_configuration_groups
      end if
      if (incidence_average) then
         write (outunit, '('' incidence averages collected: sample:'')')
         write (outunit, '(i5)') incident_direction_number * n_configuration_groups
      end if

      write (outunit, '('' number iterations, error, solution time '')')
      write (outunit, '(i6,3es12.4)') solution_iterations, solution_error, solution_time
      if (fft_translation_option .and. print_timings) then
         fft_metrics = fft_plan%performance_metrics()
         write (outunit, '('' FFT interactions, 3-D transforms'')')
         write (outunit, '(2i12)') fft_metrics%interaction_calls, fft_metrics%transform_calls
         write (outunit, '('' FFT initialize, sphere-node, node-node, node-sphere, local times'')')
         write (outunit, '(5es12.4)') fft_metrics%initialization_time, fft_metrics%sphere_to_node_time, &
            fft_metrics%node_to_node_time, fft_metrics%node_to_sphere_time, fft_metrics%local_interaction_time
         write (outunit, '('' FFT 3-D transform time'')')
         write (outunit, '(es12.4)') fft_metrics%transform_time
      end if
      if (solution_method(1:1) /= 'i') then
         write (outunit, '('' direct reciprocal condition estimate'')')
         write (outunit, '(es12.4)') solution_reciprocal_condition
         if (print_timings) then
            write (outunit, '('' direct matrix, factorization, condition, backsolve times'')')
            write (outunit, '(4es12.4)') direct_matrix_assembly_time, direct_factorization_time, &
               direct_condition_estimation_time, direct_backsolve_time
         end if
      end if
      call print_error_codes(outunit)
      if (number_plane_boundaries .gt. 0) then
         call boundary_energy_transfer(incident_sin_beta, incident_direction, r, t, a)
      else
         r = 0.d0
         a = 0.d0
         t = 1.d0
      end if
      if (random_orientation) then
         write (outunit, '('' calculated t matrix order:'')')
         write (outunit, '(i6)') t_matrix_order
      end if

      if (print_sphere_data) then
         write (outunit, '('' sphere extinction, absorption, volume absorption efficiencies (unpolarized incidence)'')')
         if (configuration_average) then
            if (target_shape .le. 1) then
               write (outunit, '(''   sphere z_ave       Qext        Qabs        Qvabs'')')
            else
               write (outunit, '(''   sphere r_ave       Qext        Qabs        Qvabs'')')
            end if
         else
            write (outunit, '(''   sphere   Qext        Qabs        Qvabs'')')
         end if
         do n = 1, number_spheres
            if (configuration_average) then
               write (outunit, '(i8,4es12.4)') n, sphere_position(3, n), q_eff(1:2, 1, n), q_vabs(1, n)
            else
               write (outunit, '(i8,3es12.4)') n, q_eff(1:2, 1, n), q_vabs(1, n)
            end if
         end do
      end if

      if (periodic_lattice .or. reflection_model) then
         rl = pl_sca(:, 2) - boundary_ext(:, 0) + r(:)
!            tl=pl_sca(:,1)+boundary_ext(:,number_plane_boundaries+1)+t(:)
!            al=1.d0-rl-tl
         al = surface_absorptance(:)
         tl = 1.d0 - rl - al
         write (outunit, '('' unit cell reflectance, absorptance, transmittance (unpol, par, perp)'')')
         write (outunit, '(9es12.4)') 0.5d0 * sum(rl), 0.5d0 * sum(al), 0.5d0 * sum(tl), &
            (rl(i), al(i), tl(i), i=1, 2)
!            write(outunit,'('' unit cell reflectance, absorptance, transmittance (unpol, par, perp)'')')
!               write(outunit,'(9es12.4)') 0.5d0*sum(pl_sca(:,2)-boundary_ext(:,0)+r(:)),q_eff_tot(2,1), &
!                 .5d0*sum(pl_sca(:,1)+boundary_ext(:,number_plane_boundaries+1)+t(:)), &
!                 (pl_sca(i,2)-boundary_ext(i,0)+r(i),q_eff_tot(2,i+1), &
!                  pl_sca(i,1)+boundary_ext(i,number_plane_boundaries+1)+t(i),i=1,2)

!            write(outunit,'('' down, up scattering fraction, total scattering cross section (unpol)'')')
!               write(outunit,'(9es12.4)') 0.5d0*sum(pl_sca(:,2)),.5d0*sum(pl_sca(:,1)), &
!                 pi*(0.5d0*sum(pl_sca(:,2))+.5d0*sum(pl_sca(:,1)))*cross_section_radius**2
         if (configuration_average .and. single_origin_expansion) then
            scacoef = (-0.5d0 * sum(dif_boundary_sca(:, 0)) + 0.5d0 * sum(dif_boundary_sca(:, 1))) &
                      * pi * cross_section_radius**2
            scarat = scacoef / (q_eff_tot(3, 1) * pi * cross_section_radius**2)
            write (outunit, '('' down, up diffuse scattering fraction, optically thin dependent/independent ratio'')')
            write (outunit, '(9es12.4)') - 0.5d0 * sum(dif_boundary_sca(:, 0)), &
               0.5d0 * sum(dif_boundary_sca(:, number_plane_boundaries + 1)), &
               scarat
            call effective_extinction_coefficient_ratio(scacoef, abscoef, scarat, absrat)
            write (outunit, '('' dimensionless extinction, absorption coefficients, dependent/independent ratios'')')
            write (outunit, '(4es12.4)') scacoef, abscoef, scarat, absrat
         end if
      else
         if (qeff_dim .eq. 1) then
            write (outunit, '('' total extinction, absorption, scattering efficiencies (unpolarized incidence)'')')
            write (outunit, '(3es12.4)') q_eff_tot(1:3, 1)
         else
            write (outunit, '('' total extinction, absorption, scattering efficiencies (unpol, par, perp incidence)'')')
            write (outunit, '(9es12.4)') q_eff_tot(1:3, 1:3)
         end if
         if (.not. random_orientation) then
            if (number_plane_boundaries .gt. 0) then
               write (outunit, '(''  down and up extinction efficiencies (unpol, par, perp)'')')
               write (outunit, '(16es12.4)') &
                  0.5d0 * sum(boundary_ext(:, 0)), 0.5d0 * sum(-boundary_ext(:, 1)), &
                  (boundary_ext(i, 0), -boundary_ext(i, 1), i=1, 2)
            end if
            if (calculate_up_down_scattering) then
               write (outunit, '(''  down and up hemispherical scattering efficiencies (unpol, par, perp) '')')
               write (outunit, '(16es12.4)') &
                  0.5d0 * sum(-boundary_sca(:, 0)), 0.5d0 * sum(boundary_sca(:, 1)), &
                  (-boundary_sca(i, 0), boundary_sca(i, 1), i=1, 2)
            end if
            if (number_plane_boundaries .gt. 0 .and. number_singular_points .gt. 0) then
               write (outunit, '(''  waveguide scattering efficiencies (unpol, par, perp)  '')')
               write (outunit, '(16es12.4)') &
                  0.5d0 * sum(evan_sca(:)), (evan_sca(i), i=1, 2)
            end if
         end if
      end if
      if (configuration_average .and. single_origin_expansion .and. .not. random_orientation) then
         call diffuse_scattering_effective_refractive_index(imrieff, qeeff, scarat)
         write (outunit, '('' diffuse/total scattering ratio, extinction coefficient, extinction prob., albedo, RT ratio'')')
         write (outunit, '(6es12.4)') dif_csca_ratio, 2.d0 * imrieff, qeeff, &
            dif_csca_ratio * q_eff_tot(3, 1) / (dif_csca_ratio * q_eff_tot(3, 1) + q_eff_tot(2, 1)), scarat
         write (outunit, '('' NI formulation ratio'')')
         write (outunit, '(5es12.4)') tot_csca, dif_csca, &
            pi * (dif_csca) / (pi * cross_section_radius**2 * q_eff_tot(1, 1))
      end if
      if (configuration_average .and. target_shape .eq. 2 .and. (.not. random_configuration_host)) then
         scarat = mean_qext_mie * area_mean_radius**2 * 3.d0 * dble(number_spheres) / (4.d0 * fit_radius**3)
         scarat = 2 * aimag(effective_ref_index) / scarat
         write (outunit, '('' mie fit effective ri, radius, RT ratio, status'')')
         write (outunit, '(4es12.4,i5)') effective_ref_index, fit_radius, scarat, fit_stat
      end if
      if (configuration_average .and. target_shape .eq. 2 .and. calculate_near_field) then
         write (outunit, '('' field fit effective ri'')')
         write (outunit, '(2es12.4)') nf_eff_ref_index
      end if

      if (calculate_scattering_matrix .and. .not. periodic_lattice) then
         if (normalize_s11) then
            s11scale = 1.d0 / ((cross_section_radius**2) * pi * q_eff_tot(3, 1))
         else
!               s11scale=(cross_section_radius**2)*q_eff_tot(3,1)
            s11scale = two_pi
!               if(.not.random_orientation) s11scale=pi*s11scale
         end if
         if (((.not. any_optically_active) .and. random_orientation) .or. azimuthal_average) then
            nsmat = 6
            smvecp(1:6) = smvec(1:6)
         else
            nsmat = 16
            smvecp(1:16) = smvec0(1:16)
         end if
         if (random_orientation) then
            write (outunit, '('' total scattering'')')
            call print_scattering_matrix_header()
            do i = scat_mat_ldim, scat_mat_udim
               smt = scaled_scattering_matrix(scat_mat(1:16, i))
               smt(1) = smt(1) * s11scale
               call print_scattering_matrix_row(i, smt)
            end do
         else
            if (scattering_map_model .eq. 0) then
               if (number_plane_boundaries .eq. 0) then
                  if (azimuthal_average) then
                     write (outunit, '('' azimuthal averaged scattering matrix'')')
                  else
                     if (incident_frame) then
                        write (outunit, '('' scattering matrix in incident plane: 0 deg = incident direction'')')
                     else
                        write (outunit, '('' scattering matrix in incident plane: 0 deg= z axis'')')
                     end if
                  end if
                  call print_scattering_matrix_header()
                  do i = scat_mat_ldim, scat_mat_udim
                     smt = scaled_scattering_matrix(scat_mat(1:16, i))
                     smt(1) = smt(1) * s11scale
                     call print_scattering_matrix_row(i, smt)
                  end do
                  if ((configuration_average .or. incidence_average) .and. single_origin_expansion) then
!write(*,'(es12.4)') s11scale
                     write (outunit, '('' diffuse scattering matrix '')')
                     call print_scattering_matrix_header(no_numbers=.true.)
                     do i = scat_mat_ldim, scat_mat_udim
                        smt = scaled_scattering_matrix(dif_scat_mat(1:16, i))
                        smt(1) = smt(1) * s11scale
                        call print_scattering_matrix_row(i, smt)
                     end do
                  end if
               else
                  write (outunit, '('' scattering matrix in incident plane'')')
                  write (outunit, '('' reflection'')')
                  call print_scattering_matrix_header()
                  do i = scat_mat_ldim, scat_mat_udim
                     smt = scaled_scattering_matrix(scat_mat(1:16, i))
                     smt(1) = smt(1) * s11scale
                     call print_scattering_matrix_row(i, smt)
                  end do
                  write (outunit, '('' transmission'')')
                  call print_scattering_matrix_header(no_numbers=.true.)
                  do i = scat_mat_ldim, scat_mat_udim
                     smt = scaled_scattering_matrix(scat_mat(17:32, i))
                     smt(1) = smt(1) * s11scale
                     call print_scattering_matrix_row(i, smt)
                  end do
               end if
            else
               write (outunit, '('' 2-D scattering in backward and forward hemispheres.   Number points:'')')
               write (outunit, '(i10)') scat_mat_udim
               write (outunit, '('' backward hemisphere scattering'')')
               write (outunit, '(''    kx      ky '')', advance='no')
               do j = 1, 4
                  do i = 1, 4
                     write (outunit, '(''     '',2i1,''     '')', advance='no') i, j
                  end do
               end do
               write (outunit, *)
               s = 0
               do sy = -scattering_map_dimension, scattering_map_dimension
                  ky = dble(sy) / dble(scattering_map_dimension)
                  do sx = -scattering_map_dimension, scattering_map_dimension
                     kx = dble(sx) / dble(scattering_map_dimension)
                     if (sx * sx + sy * sy .gt. scattering_map_dimension**2) cycle
                     s = s + 1
                     smt = scaled_scattering_matrix(scat_mat(1:16, s))
                     smt(1) = smt(1) * s11scale
                     write (outunit, '(2f9.5)', advance='no') kx, ky
                     do j = 1, 16
                        write (outunit, '(es12.4)', advance='no') smt(smvec0(j))
                     end do
                     write (outunit, *)
                  end do
               end do
               write (outunit, '('' forward hemisphere scattering'')')
               write (outunit, '(''    kx      ky '')', advance='no')
               do j = 1, 4
                  do i = 1, 4
                     write (outunit, '(''     '',2i1,''     '')', advance='no') i, j
                  end do
               end do
               write (outunit, *)
               s = 0
               do sy = -scattering_map_dimension, scattering_map_dimension
                  ky = dble(sy) / dble(scattering_map_dimension)
                  do sx = -scattering_map_dimension, scattering_map_dimension
                     if (sx * sx + sy * sy .gt. scattering_map_dimension**2) cycle
                     kx = dble(sx) / dble(scattering_map_dimension)
                     s = s + 1
                     smt = scaled_scattering_matrix(scat_mat(17:32, s))
                     smt(1) = smt(1) * s11scale
                     write (outunit, '(2f9.5)', advance='no') kx, ky
                     do j = 1, 16
                        write (outunit, '(es12.4)', advance='no') smt(smvec0(j))
                     end do
                     write (outunit, *)
                  end do
               end do
            end if
         end if
         if (azimuthal_average .and. (.not. random_orientation) .and. (.not. numerical_azimuthal_average)) then
!               s11scale=1.d0/scat_mat_exp_coef(1,0,1)
!               scat_mat_exp_coef=scat_mat_exp_coef*s11scale
            write (outunit, '('' azimuthal averaged scattering matrix expansion coefficients, total field'')')
            write (outunit, '(''    n         11            44            12            34           22p           22m'')')
            do n = 0, 2 * t_matrix_order
               write (outunit, '(i5,6es14.6)') n, scat_mat_exp_coef(1, n, 1), scat_mat_exp_coef(16, n, 1), &
                  0.5d0 * (scat_mat_exp_coef(2, n, 2) + scat_mat_exp_coef(5, n, 2)), &
                  0.5d0 * (scat_mat_exp_coef(12, n, 2) + scat_mat_exp_coef(15, n, 2)), &
                  scat_mat_exp_coef(6, n, 3), scat_mat_exp_coef(6, n, 4)
            end do
            if (configuration_average) then
               write (outunit, '('' azimuthal averaged scattering matrix expansion coefficients, coherent field'')')
               write (outunit, '(''    n         11            44            12            34           22p           22m'')')
               do n = 0, 2 * t_matrix_order
                  write (outunit, '(i5,6es14.6)') n, coh_scat_mat_exp_coef(1, n, 1), coh_scat_mat_exp_coef(16, n, 1), &
                     0.5d0 * (coh_scat_mat_exp_coef(2, n, 2) + coh_scat_mat_exp_coef(5, n, 2)), &
                     0.5d0 * (coh_scat_mat_exp_coef(12, n, 2) + coh_scat_mat_exp_coef(15, n, 2)), &
                     coh_scat_mat_exp_coef(6, n, 3), coh_scat_mat_exp_coef(6, n, 4)
               end do
            end if
         end if
         if (random_orientation) then
            write (outunit, '('' orientation averaged scattering matrix expansion coefficients'')')
!               write(outunit,'(''    w  a11           a22           a33           '',&
!                &''a23           a32           a44           a12         '',&
!                &''a34           a13           a24           a14'')')
!               do n=0,2*t_matrix_order
!                  write(outunit,'(i5,11es14.6)') w,scat_mat_exp_coef(1,1,n),scat_mat_exp_coef(2,2,n),&
!                    scat_mat_exp_coef(3,3,n),scat_mat_exp_coef(2,3,n),scat_mat_exp_coef(3,2,n),scat_mat_exp_coef(4,4,n),&
!                    scat_mat_exp_coef(1,2,n),scat_mat_exp_coef(3,4,n),scat_mat_exp_coef(1,3,n),scat_mat_exp_coef(2,4,n),&
!                    scat_mat_exp_coef(1,4,n)
!               enddo
            write (outunit, '(''    n         11            44            12            34           22            33'')')
            do n = 0, 2 * t_matrix_order
               write (outunit, '(i5,11es14.6)') n, scat_mat_exp_coef(1, 1, n), scat_mat_exp_coef(4, 4, n), &
                  scat_mat_exp_coef(1, 2, n), scat_mat_exp_coef(3, 4, n), &
                  scat_mat_exp_coef(2, 2, n), scat_mat_exp_coef(3, 3, n)
            end do
            write (outunit, '('' orientation averaged coherent scattering matrix expansion coefficients'')')
            write (outunit, '(''    n         11            44            12            34           22            33'')')
            do n = 0, 2 * t_matrix_order
               write (outunit, '(i5,11es14.6)') n, coh_scat_mat_exp_coef(1, 1, n), coh_scat_mat_exp_coef(4, 4, n), &
                  coh_scat_mat_exp_coef(1, 2, n), coh_scat_mat_exp_coef(3, 4, n), &
                  coh_scat_mat_exp_coef(2, 2, n), coh_scat_mat_exp_coef(3, 3, n)
            end do
         end if
      end if

      if (calculate_scattering_matrix .and. periodic_lattice) then
         write (outunit, '('' scattering by periodic lattice at reciprocal lattice directions'')')
         write (outunit, '('' backward hemisphere scattering'')')
         write (outunit, '('' number directions, number SM elements'')')
         write (outunit, '(2i6)') number_rl_dirs(1), 16
         write (outunit, '(''    kx      ky '')', advance='no')
         do j = 1, 4
            do i = 1, 4
               write (outunit, '(''     '',2i1,''     '')', advance='no') i, j
            end do
         end do
         write (outunit, *)
         do i = 1, number_rl_dirs(1)
            smt = scaled_scattering_matrix(scat_mat(1:16, i))
            write (outunit, '(2f9.5)', advance='no') rl_vec(1:2, i) / dble(layer_ref_index(0))
            do j = 1, 16
               write (outunit, '(es12.4)', advance='no') smt(smvec0(j))
            end do
            write (outunit, *)
         end do
         write (outunit, '('' forward hemisphere scattering'')')
         write (outunit, '('' number directions, number SM elements'')')
         write (outunit, '(2i6)') number_rl_dirs(2), 16
         write (outunit, '(''    kx      ky '')', advance='no')
         do j = 1, 4
            do i = 1, 4
               write (outunit, '(''     '',2i1,''     '')', advance='no') i, j
            end do
         end do
         write (outunit, *)
         do i = 1, number_rl_dirs(2)
            smt = scaled_scattering_matrix(scat_mat(17:32, i))
            write (outunit, '(2f9.5)', advance='no') &
               rl_vec(1:2, i) / dble(layer_ref_index(number_plane_boundaries))
            do j = 1, 16
               write (outunit, '(es12.4)', advance='no') smt(smvec0(j))
            end do
            write (outunit, *)
         end do
      end if

      if (random_orientation .and. configuration_average) then
         write (outunit, '('' mean t matrix elements (p=1,2)'')')
         do i = 1, t_matrix_order
            write (outunit, '(i5,4es12.4)') i, mean_t(1, i), mean_t(2, i)
         end do
      end if

      if (outunit .ne. 6) close (outunit)

   contains
      subroutine print_scattering_matrix_header(no_numbers)
         implicit none
         logical, optional :: no_numbers
         if (present(no_numbers)) then
            if (.not. no_numbers) then
               write (outunit, '('' number directions, number SM elements:'')')
               write (outunit, '(2i6)') scat_mat_udim - scat_mat_ldim + 1, nsmat
            end if
         else
            write (outunit, '('' number directions, number SM elements:'')')
            write (outunit, '(2i6)') scat_mat_udim - scat_mat_ldim + 1, nsmat
         end if
         write (outunit, '(''   theta'')', advance='no')
         do i = 1, nsmat
            write (outunit, '(''     '',a2,''     '')', advance='no') smlabel(smvecp(i))
         end do
         write (outunit, *)
      end subroutine print_scattering_matrix_header

      subroutine print_scattering_matrix_row(i, smt)
         implicit none
         integer :: i
         real(real64) :: smt(16)
!            write(outunit,'(f8.2)',advance='no') dble(i)*180.d0/dble(scat_mat_udim-scat_mat_ldim)
         write (outunit, '(f8.2)', advance='no') scat_mat_amin &
            + dble(i - scat_mat_ldim) / dble(scat_mat_udim - scat_mat_ldim) * (scat_mat_amax - scat_mat_amin)
         do j = 1, nsmat
            write (outunit, '(es12.4)', advance='no') smt(smvecp(j))
         end do
         write (outunit, *)
      end subroutine print_scattering_matrix_row

   end subroutine print_calculation_results

   pure function scaled_scattering_matrix(s)
      real(real64), intent(in) :: s(16)
      real(real64) :: scaled_scattering_matrix(16)
      if (s(1) .eq. 0.d0) then
         scaled_scattering_matrix = 0.d0
      else
         scaled_scattering_matrix(1) = s(1)
         scaled_scattering_matrix(2:16) = s(2:16) / s(1)
      end if
   end function scaled_scattering_matrix

   pure subroutine scattering_matrix_to_phase_matrix(smat, u, up, phi, smatrot)
      implicit none
      real(real64), intent(in) :: smat(4, 4), u, up, phi
      real(real64), intent(out) :: smatrot(4, 4)
      real(real64) :: s, sp, us, ss, csig, ssig, c2sig, s2sig, mat1(4, 4), mat2(4, 4)

      s = sqrt(1.d0 - u * u)
      sp = sqrt(1.d0 - up * up)
      us = u * up + s * sp * cos(phi)
      ss = sqrt(1.d0 - us * us)
      if (up .eq. 1.d0) then
         csig = cos(phi)
         ssig = -sin(phi)
      elseif (up .eq. -1.d0) then
         csig = -cos(phi)
         ssig = -sin(phi)
      else
         csig = (-u + up * us) / (sp * ss)
         ssig = -s * sin(phi) / ss
      end if
      c2sig = 2.d0 * csig * csig - 1.d0
      s2sig = 2.d0 * ssig * csig
      mat1(:, 1) = (/1.d0, 0.d0, 0.d0, 0.d0/)
      mat1(:, 2) = (/0.d0, c2sig, -s2sig, 0.d0/)
      mat1(:, 3) = (/0.d0, s2sig, c2sig, 0.d0/)
      mat1(:, 4) = (/0.d0, 0.d0, 0.d0, 1.d0/)
      if (u .eq. 1.d0) then
         csig = cos(phi)
         ssig = -sin(phi)
      elseif (u .eq. -1.d0) then
         csig = -cos(phi)
         ssig = -sin(phi)
      else
         csig = (-up + u * us) / (s * ss)
         ssig = -sp * sin(phi) / ss
      end if
      c2sig = 2.d0 * csig * csig - 1.d0
      s2sig = 2.d0 * ssig * csig
      mat2(:, 1) = (/1.d0, 0.d0, 0.d0, 0.d0/)
      mat2(:, 2) = (/0.d0, c2sig, -s2sig, 0.d0/)
      mat2(:, 3) = (/0.d0, s2sig, c2sig, 0.d0/)
      mat2(:, 4) = (/0.d0, 0.d0, 0.d0, 1.d0/)
      mat1 = matmul(smat, mat1)
      smatrot = matmul(mat2, mat1)
   end subroutine scattering_matrix_to_phase_matrix
end module input_reporting
