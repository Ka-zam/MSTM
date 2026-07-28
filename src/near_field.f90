module near_field
   use, intrinsic :: iso_fortran_env, only: real32, real64
   use excitation, only: electric_dipole_field, magnetic_current_field, magnetic_current_segment_t
   use parallel_runtime, only: mpi_comm_world, mstm_global_rank, parallel_rank, parallel_reduce_sum, &
                               parallel_size, parallel_wall_time
   use wave_functions, only: vector_spherical_wave_functions
   use sphere_data
   use mie
   use numerical_tables
   use surface
   use periodic_lattice_operations
   use scattering_interactions, only: layered_gaussian_beam_coefficients
   use translation_operator, only: translation_operator_state
   use runtime_support, only: set_runtime_error
   implicit none
   type grid_info
      logical :: initialized
      integer :: cellnum, host, layer
      type(cell_info), pointer :: cellinfo => null()
   end type grid_info
   type cell_info
      logical :: outside_spheres
      integer :: ncell(3), layer, host, order, nispheres
      real(real64) :: rcell(3)
      integer, allocatable :: sphere_indices(:)
      complex(real64), pointer :: vector(:, :, :) => null()
      complex(real64), pointer :: reg_source_vector(:, :, :) => null()
      complex(real64), pointer :: gb_vector(:, :, :) => null()
   end type cell_info
   type linked_cell_list
      type(cell_info) :: cellinfo
      type(linked_cell_list), pointer :: next => null()
   end type linked_cell_list
   type vector_storage
      complex(real64), pointer :: vector(:, :, :) => null()
   end type vector_storage
   type linked_sphere_data
      integer :: sphere, host
      real(real64) :: position(3), radius
      type(linked_sphere_data), pointer :: next
   end type linked_sphere_data

   type(linked_cell_list), pointer, private :: cell_info_list
   type(vector_storage), allocatable, private :: internal_field_vector(:)
   logical :: incident_gb
   logical, target :: store_surface_vector, fast_near_field
   integer, private :: local_rank, local_numprocs, local_run_number, total_cells, number_intersecting_spheres
   integer, target :: near_field_expansion_order
   real(real64) :: grid_region(3, 2), grid_spacing(3)
   real(real64), target :: near_field_expansion_spacing
   type(linked_sphere_data), pointer, private :: intersecting_spheres
   complex(real64), private :: vwf_0(3, 3, 2)
   data near_field_expansion_order, near_field_expansion_spacing, store_surface_vector/10, 5.d0, .true./
   data local_run_number, fast_near_field/1, .true./

contains

   subroutine compute_near_field(amnp, alpha, sinc, dir, gridregion, griddim, incident_model, output_unit, &
                                 e_field_array, h_field_array, e_field_ave_array, output_header, mpi_comm, &
                                 dipole_position, dipole_moment, magnetic_segments, magnetic_quadrature_order)
      implicit none
      logical, optional :: output_header
      logical :: dipole_incident, localized_source_incident, magnetic_current_incident
      integer :: incmodel, i, p, outputunit, nblk, l, griddim(3), ix, iy, iz, &
                 layer, host, ipos(3), cellnum, nsend, mpicomm, totpoints, point, dir
      integer, optional :: incident_model, output_unit, mpi_comm, magnetic_quadrature_order
      real(real64) :: gridregion(3, 2), rpos(3), alpha, time1, time0, rtemp, sinc
      real(real64), intent(in), optional :: dipole_position(3)
      complex(real64) :: amnp(sphere_cluster%number_eqns, 2), evec(3, 2), hvec(3, 2), evec1(3, 2), hvec1(3, 2)
      complex(real64), intent(in), optional :: dipole_moment(3)
      type(magnetic_current_segment_t), intent(in), optional :: magnetic_segments(:)
      complex(real64), allocatable :: vector(:, :, :)
      complex(real64), optional :: e_field_ave_array(3, 2, griddim(3))
      complex(real64), target, optional :: e_field_array(3, 2, griddim(1), griddim(2), griddim(3)), &
                                           h_field_array(3, 2, griddim(1), griddim(2), griddim(3))
      complex(real32) :: earray(3, 2, griddim(1), griddim(2)), &
                         harray(3, 2, griddim(1), griddim(2))
      type(grid_info), allocatable :: gridinfo(:, :, :)
      type(cell_info), pointer :: cellinfo
      type(linked_cell_list), pointer :: clist1, clist2
!         complex(real32), pointer :: earray(:,:,:,:,:),harray(:,:,:,:,:)

      if (present(incident_model)) then
         incmodel = incident_model
      else
         incmodel = 1
      end if
      if (present(output_unit)) then
         outputunit = output_unit
      else
         outputunit = 0
      end if
      if (present(mpi_comm)) then
         mpicomm = mpi_comm
      else
         mpicomm = mpi_comm_world
      end if
      incident_gb = sphere_cluster%gaussian_beam_constant .ne. 0.d0
      dipole_incident = present(dipole_position) .and. present(dipole_moment)
      magnetic_current_incident = present(magnetic_segments) .and. present(magnetic_quadrature_order)
      localized_source_incident = dipole_incident .or. magnetic_current_incident

      call parallel_size(mpi_size=local_numprocs, mpi_comm=mpicomm)
      call parallel_rank(mpi_rank=local_rank, mpi_comm=mpicomm)
      do i = 1, 3
         if (griddim(i) .eq. 1) then
            grid_spacing(i) = 0.d0
            rtemp = 0.5d0 * (gridregion(i, 1) + gridregion(i, 2))
            grid_region(i, :) = rtemp
         else
            grid_spacing(i) = (gridregion(i, 2) - gridregion(i, 1)) / dble(griddim(i))
            grid_region(i, :) = gridregion(i, :)
         end if
      end do
      totpoints = product(griddim)

      if (dipole_incident) then
         do iz = 1, griddim(3)
            do iy = 1, griddim(2)
               do ix = 1, griddim(1)
                  ipos = [ix, iy, iz]
                  rpos = (dble(ipos) - 0.5_real64) * grid_spacing + grid_region(:, 1)
                  if (sum((rpos - dipole_position)**2) <= 1.0e-24_real64) then
                     call set_runtime_error('Near-field grid contains the singular electric-dipole source point')
                     return
                  end if
               end do
            end do
         end do
      elseif (magnetic_current_incident) then
         do iz = 1, griddim(3)
            do iy = 1, griddim(2)
               do ix = 1, griddim(1)
                  ipos = [ix, iy, iz]
                  rpos = (dble(ipos) - 0.5_real64) * grid_spacing + grid_region(:, 1)
                  do i = 1, size(magnetic_segments)
                     if (distance_to_segment(rpos, magnetic_segments(i)%start_point, &
                                             magnetic_segments(i)%end_point) <= 1.0e-12_real64) then
                        call set_runtime_error('Near-field grid contains a singular magnetic-current segment point')
                        return
                     end if
                  end do
               end do
            end do
         end do
      end if

      call vector_spherical_wave_functions((/0.d0, 0.d0, 0.d0/), &
                                           (/(1.d0, 0.d0), (1.d0, 0.d0)/), 1, 1, vwf_0, index_model=2)

      allocate (gridinfo(griddim(1), griddim(2), griddim(3)))
      call initialize_grid_points(griddim, gridinfo)

      if (allocated(internal_field_vector)) then
         l = ubound(internal_field_vector, 1)
         do i = 1, l
            if (associated(internal_field_vector(i)%vector)) deallocate (internal_field_vector(i)%vector)
         end do
         deallocate (internal_field_vector)
      end if
      allocate (internal_field_vector(sphere_cluster%number_spheres))
      do i = 1, sphere_cluster%number_spheres
         nblk = sphere_cluster%sphere_order(i) * (sphere_cluster%sphere_order(i) + 2)
         allocate (vector(nblk, 2, 2))
         allocate (internal_field_vector(i)%vector(nblk, 2, 2))
         do p = 1, 2
            if (sphere_cluster%is_pec(i)) then
               vector(1:nblk, 1:2, p) = 0.0_real64
            elseif (sphere_cluster%number_field_expansions(i) .eq. 1) then
               call apply_single_sphere_mie_coefficients(i, sphere_cluster%sphere_order(i), &
                    amnp(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i), p), &
                                                         vector(1:nblk, 1:2, p), 'c')
            else
               vector(1:nblk, 1:2, p) &
                  = reshape(amnp(sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i) + 1:sphere_cluster%sphere_offset(i) + 2 * sphere_cluster%sphere_block(i), p), &
                            (/nblk, 2/))
            end if
         end do
         internal_field_vector(i)%vector = vector
         deallocate (vector)
      end do

      if (present(output_header)) then
         if (output_header .and. outputunit .ne. 0) call write_near_field_output_header(griddim, outputunit)
      end if
      if (local_rank .eq. 0 .and. outputunit .ne. 0) then
         write (sphere_cluster%run_print_unit, '('' near field calculation'')')
         write (sphere_cluster%run_print_unit, '('' grid minimum x,y,z:'',3es12.4)') grid_region(:, 1)
         write (sphere_cluster%run_print_unit, '('' grid maximum x,y,z:'',3es12.4)') grid_region(:, 2)
         write (sphere_cluster%run_print_unit, '('' grid x,y,z dimensions, total points:'',3i5,i12)') griddim(:), product(griddim)
         write (sphere_cluster%run_print_unit, '('' total reexpansion cells:'',i5)') total_cells
         flush (sphere_cluster%run_print_unit)
      end if

!clist=>cell_info_list
!open(20,file='ctest.dat')
!do i=1,total_cells
!write(20,'('' cell:'',i5)') i
!write(20,'('' ncell,host,layer,order:'',6i5)') clist%ncell,clist%host,clist%layer,clist%order
!write(20,'('' rcell:'',3es12.4)') clist%rcell
!write(20,'('' vector stored:'',l2)') associated(clist%vector)
!write(20,*)
!clist=>clist%next
!if(.not.associated(clist)) exit
!enddo

      time0 = parallel_wall_time()
      point = 0
      if (local_rank .eq. 0 .and. present(e_field_ave_array)) e_field_ave_array = 0.d0
      do iz = 1, griddim(3)
         earray = 0.
         harray = 0.
         do iy = 1, griddim(2)
            do ix = 1, griddim(1)
               cellnum = gridinfo(ix, iy, iz)%cellnum
               point = point + 1
               time1 = parallel_wall_time()
!                  if(time1-time0.gt.15.d0.and.local_rank.eq.0.and.(.not.present(e_field_ave_array))) then
               if (time1 - time0 .gt. 15.d0 .and. mstm_global_rank .eq. 0) then
                  write (sphere_cluster%run_print_unit, '('' completed field point calculation '',i8,''/'',i8)') point, totpoints
                  flush (sphere_cluster%run_print_unit)
                  time0 = time1
               end if
               if (mod(cellnum, local_numprocs) .ne. local_rank) cycle
!
               cellinfo => gridinfo(ix, iy, iz)%cellinfo
               layer = cellinfo%layer
               host = cellinfo%host
               ipos(:) = (/ix, iy, iz/)
               rpos(:) = (dble(ipos(:)) - (/0.5d0, 0.5d0, 0.5d0/)) * grid_spacing(:) + grid_region(:, 1)
!write(*,*) 'nf 1', layer,host
!flush(6)
!if(local_rank.eq.0) then
!write(*,'(3i5,3es12.4)') ix,iz,host,rpos
!flush(6)
!endif
               if (host /= 0 .and. sphere_cluster%is_pec(host)) then
                  evec = 0.0_real64
                  hvec = 0.0_real64
               else
                  if ((.not. periodic_lattice) .and. fast_near_field) then
!                     if(host.eq.0.and.(.not.periodic_lattice).and.(number_plane_boundaries.eq.0)) then
                     call calculate_source_field_fast(rpos, amnp, cellinfo, evec, hvec)
                  else
                     call calculate_source_field(rpos, amnp, host, layer, evec, hvec)
                  end if
!write(*,*) 'nf 2'
!flush(6)
                  if (host .eq. 0 .and. (number_plane_boundaries .gt. 0 .or. periodic_lattice)) then
                     call calculate_surface_field(rpos, amnp, cellinfo, evec1, hvec1)
                     evec(:, :) = evec(:, :) + evec1(:, :)
                     hvec(:, :) = hvec(:, :) + hvec1(:, :)
                  end if
!write(*,*) 'nf 3'
!flush(6)
                  if (incmodel .ne. 2 .and. host .eq. 0) then
                     call calculate_incident_field(rpos, layer, alpha, sinc, dir, cellinfo, evec1, hvec1, &
                                                   dipole_position, dipole_moment, magnetic_segments, &
                                                   magnetic_quadrature_order)
                     evec(:, :) = evec(:, :) + evec1(:, :)
                     hvec(:, :) = hvec(:, :) + hvec1(:, :)
                  elseif (incmodel .eq. 2 .and. host .ne. 0) then
                     call calculate_incident_field(rpos, layer, alpha, sinc, dir, cellinfo, evec1, hvec1, &
                                                   dipole_position, dipole_moment, magnetic_segments, &
                                                   magnetic_quadrature_order)
                     evec(:, :) = evec(:, :) - evec1(:, :)
                     hvec(:, :) = hvec(:, :) - hvec1(:, :)
                  end if
               end if
               earray(:, :, ix, iy) = evec
               harray(:, :, ix, iy) = hvec
            end do
         end do

         if (local_numprocs .gt. 1) then
            nsend = 3 * 2 * product(griddim(1:2))
            call parallel_reduce_sum(receive_buffer=earray, &
                                     mpi_number=nsend, mpi_rank=0, mpi_comm=mpicomm)
            call parallel_reduce_sum(receive_buffer=harray, &
                                     mpi_number=nsend, mpi_rank=0, mpi_comm=mpicomm)
         end if
         if (local_rank .eq. 0 .and. present(e_field_ave_array)) then
            do iy = 1, griddim(2)
               do ix = 1, griddim(1)
                  e_field_ave_array(:, :, iz) = e_field_ave_array(:, :, iz) &
                                                + earray(:, :, ix, iy)
               end do
            end do
            e_field_ave_array(:, :, iz) = e_field_ave_array(:, :, iz) / dble(griddim(1) * griddim(2))
         end if

         if (outputunit .ne. 0 .and. local_rank .eq. 0) then
            do iy = 1, griddim(2)
               do ix = 1, griddim(1)
                  ipos(:) = (/ix, iy, iz/)
                  rpos(:) = (dble(ipos(:)) - (/0.5d0, 0.5d0, 0.5d0/)) * grid_spacing(:) + grid_region(:, 1)
                  if (localized_source_incident) then
                     write (outputunit, '(15es12.4)') rpos(:), earray(:, 1, ix, iy), harray(:, 1, ix, iy)
                  else
                     write (outputunit, '(27es12.4)') rpos(:), earray(:, 1, ix, iy), harray(:, 1, ix, iy), &
                        earray(:, 2, ix, iy), harray(:, 2, ix, iy)
                  end if
               end do
            end do
         end if
         if (local_rank .eq. 0 .and. present(e_field_array)) then
            e_field_array(:, :, :, :, iz) = earray
         end if
         if (local_rank .eq. 0 .and. present(h_field_array)) then
            h_field_array(:, :, :, :, iz) = harray
         end if
      end do

      clist1 => cell_info_list
      do while (associated(clist1))
         clist2 => clist1%next
         deallocate (clist1)
         nullify (clist1)
         clist1 => clist2
      end do
      if (associated(cell_info_list)) then
         nullify (cell_info_list)
      end if
      deallocate (gridinfo)

   end subroutine compute_near_field

   subroutine calculate_source_field_fast(rpos, sourcevec, cellinfo, evec, hvec)
      implicit none
      integer :: nodr, nblk, i, p, j, cellhost, celllayer
      real(real64) :: rpos(3), rtran(3)
      complex(real64) :: sourcevec(sphere_cluster%number_eqns, 2), evec(3, 2), hvec(3, 2), ri2(2), ri
      complex(real64), allocatable :: vwf(:, :, :), svec(:, :)
      type(cell_info), pointer :: cellinfo

      evec = 0.d0
      hvec = 0.d0
      if (.not. associated(cellinfo%reg_source_vector)) then
         call calculate_stored_source_vector(sourcevec, cellinfo)
      end if
      cellhost = cellinfo%host
      celllayer = cellinfo%layer
      if (cellhost .eq. 0) then
         ri = layer_ref_index(celllayer)
         ri2 = ri
      else
         ri2 = sphere_cluster%sphere_ref_index(:, cellhost)
         ri = 2.d0 / (1.d0 / ri2(1) + 1.d0 / ri2(2))
      end if
      if (cellinfo%nispheres .gt. 0) then
         do i = 1, cellinfo%nispheres
            j = cellinfo%sphere_indices(i)
            nodr = sphere_cluster%sphere_order(j)
            nblk = nodr * (nodr + 2)
            allocate (vwf(3, nblk, 2), svec(nblk, 2))
            rtran(:) = rpos(:) - sphere_cluster%sphere_position(:, j)
            call vector_spherical_wave_functions(rtran, ri2, nodr, 3, vwf, index_model=2)
            do p = 1, 2
svec = reshape(sourcevec(sphere_cluster%sphere_offset(j) + 1:sphere_cluster%sphere_offset(j) + sphere_cluster%sphere_block(j), p), &
                              (/nblk, 2/))
               evec(:, p) = evec(:, p) + (matmul(vwf(:, :, 1), svec(:, 1)) + matmul(vwf(:, :, 2), svec(:, 2)))
               hvec(:, p) = hvec(:, p) + (matmul(vwf(:, :, 1), svec(:, 1)) - matmul(vwf(:, :, 2), svec(:, 2))) * ri / (0.d0, 1.d0)
            end do
            deallocate (vwf, svec)
         end do
      end if
      if (cellinfo%outside_spheres) then
         nodr = cellinfo%order
         nblk = nodr * (nodr + 2)
         allocate (vwf(3, nblk, 2))
         rtran(:) = rpos(:) - cellinfo%rcell(:)
         call vector_spherical_wave_functions(rtran, ri2, nodr, 1, vwf, index_model=2)
         do p = 1, 2
            evec(:, p) = evec(:, p) + (matmul(vwf(:, :, 1), cellinfo%reg_source_vector(:, 1, p)) &
                                       + matmul(vwf(:, :, 2), cellinfo%reg_source_vector(:, 2, p)))
            hvec(:, p) = hvec(:, p) + (matmul(vwf(:, :, 1), cellinfo%reg_source_vector(:, 1, p)) &
                                       - matmul(vwf(:, :, 2), cellinfo%reg_source_vector(:, 2, p))) * ri / (0.d0, 1.d0)
         end do
         deallocate (vwf)
      end if
   end subroutine calculate_source_field_fast

   subroutine calculate_source_field(rpos, sourcevec, host, layer, evec, hvec)
      implicit none
      integer :: layer, host, nodr, nblk, i, p, num, j, np(2), np0(2)
      real(real64) :: rpos(3), rtran(3)
      complex(real64) :: sourcevec(sphere_cluster%number_eqns, 2), evec(3, 2), hvec(3, 2), ri, ri2(2), pshift, pshift0
      complex(real64), allocatable :: vwf(:, :, :), svec(:, :)

      evec = 0.d0
      hvec = 0.d0
!return
      num = size(sphere_cluster%sphere_links(host, layer)%indices)
      if (host .eq. 0) then
         ri = layer_ref_index(layer)
         ri2 = ri
      else
         ri2 = sphere_cluster%sphere_ref_index(:, host)
         ri = 2.d0 / (1.d0 / ri2(1) + 1.d0 / ri2(2))
      end if

      np0 = 0
      pshift0 = 1.d0
      if (host .ne. 0 .and. periodic_lattice) then
         i = host
         do while (sphere_cluster%host_sphere(i) .ne. 0)
            i = sphere_cluster%host_sphere(i)
         end do
         rtran = rpos - sphere_cluster%sphere_position(:, i)
         np0(1:2) = floor((rtran(1:2) + cell_width(1:2) / 2.d0) / cell_width(1:2))
         pshift0 = exp((0.d0, 1.d0) * sum(incident_lateral_vector * dble(np0) * cell_width))
      end if

      if (num .gt. 0) then
         do i = 1, num
            j = sphere_cluster%sphere_links(host, layer)%indices(i)
            nodr = sphere_cluster%sphere_order(j)
            nblk = nodr * (nodr + 2)
            allocate (vwf(3, nblk, 2), svec(nblk, 2))
            rtran = rpos - sphere_cluster%sphere_position(:, j)
            if (host .eq. 0 .and. periodic_lattice) then
               np(1:2) = floor((rtran(1:2) + cell_width(1:2) / 2.d0) / cell_width(1:2))
               rtran(1:2) = rtran(1:2) - cell_width(1:2) * dble(np(1:2))
               pshift = exp((0.d0, 1.d0) * sum(incident_lateral_vector * dble(np) * cell_width)) * pshift0
            else
               pshift = pshift0
            end if
            call vector_spherical_wave_functions(rtran, ri2, nodr, 3, vwf, index_model=2)
            do p = 1, 2
svec = reshape(sourcevec(sphere_cluster%sphere_offset(j) + 1:sphere_cluster%sphere_offset(j) + sphere_cluster%sphere_block(j), p), &
                              (/nblk, 2/))
               evec(:, p) = evec(:, p) + (matmul(vwf(:, :, 1), svec(:, 1)) + matmul(vwf(:, :, 2), svec(:, 2))) * pshift
        hvec(:, p) = hvec(:, p) + (matmul(vwf(:, :, 1), svec(:, 1)) - matmul(vwf(:, :, 2), svec(:, 2))) * ri / (0.d0, 1.d0) * pshift
            end do
            deallocate (vwf, svec)
         end do
      end if
      if (host .ne. 0) then
         nodr = sphere_cluster%sphere_order(host)
         nblk = nodr * (nodr + 2)
         allocate (vwf(3, nblk, 2), svec(nblk, 2))
         rtran = rpos - sphere_cluster%sphere_position(:, host)
         rtran(1:2) = rtran(1:2) - cell_width(1:2) * dble(np0(1:2))
         call vector_spherical_wave_functions(rtran, ri2, nodr, 1, vwf, index_model=2)
         do p = 1, 2
            svec(:, :) = internal_field_vector(host)%vector(:, :, p)
            evec(:, p) = evec(:, p) + (matmul(vwf(:, :, 1), svec(:, 1)) + matmul(vwf(:, :, 2), svec(:, 2))) * pshift0
       hvec(:, p) = hvec(:, p) + (matmul(vwf(:, :, 1), svec(:, 1)) - matmul(vwf(:, :, 2), svec(:, 2))) * ri / (0.d0, 1.d0) * pshift0
         end do
      end if
   end subroutine calculate_source_field

   subroutine calculate_surface_field(rpos, sourcevec, cellinfo, evec, hvec)
      implicit none
      logical :: storecalc
      integer :: layer, nodr, nblk, i, p
      real(real64) :: rpos(3), rtran(3), rcell(3)
      complex(real64) :: sourcevec(sphere_cluster%number_eqns, 2), evec(3, 2), hvec(3, 2), ri, rvec(3, 2)
      complex(real64), allocatable :: vwf(:, :, :)
      type(cell_info), pointer :: cellinfo

      layer = cellinfo%layer
      ri = layer_ref_index(layer)
      if (store_surface_vector) then
         storecalc = .true.
      else
         storecalc = .false.
      end if
!evec=0.
!hvec=0.
!return
      if (store_surface_vector) then
         nodr = cellinfo%order
         rcell = cellinfo%rcell
         nblk = nodr * (nodr + 2)
         allocate (vwf(3, nblk, 2))
         if (.not. associated(cellinfo%vector)) then
            call calculate_stored_surface_vector(nodr, rcell, sourcevec, cellinfo%vector)
         end if
         rtran(:) = rpos(:) - rcell(:)
         call vector_spherical_wave_functions(rtran, (/ri, ri/), nodr, 1, vwf, index_model=2)
         do p = 1, 2
            evec(:, p) = matmul(vwf(:, :, 1), cellinfo%vector(1:nblk, 1, p)) &
                         + matmul(vwf(:, :, 2), cellinfo%vector(1:nblk, 2, p))
            hvec(:, p) = (matmul(vwf(:, :, 1), cellinfo%vector(1:nblk, 1, p)) &
                          - matmul(vwf(:, :, 2), cellinfo%vector(1:nblk, 2, p))) * ri / (0.d0, 1.d0)
         end do
         deallocate (vwf)
      else
         evec = 0.d0
         hvec = 0.d0
         do i = 1, sphere_cluster%number_spheres
            if (sphere_cluster%host_sphere(i) .ne. 0) cycle
            rtran(:) = rpos(:) - sphere_cluster%sphere_position(:, i)
!write(*,*) rtran
!flush(6)
            do p = 1, 2
               if (periodic_lattice) then
                  call plane_boundary_lattice_interaction(1, sphere_cluster%sphere_order(i), rtran(1), rtran(2), rpos(3), &
                 sphere_cluster%sphere_position(3, i), rvec, source_vector=sourcevec(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i), p), &
                                                          include_source=.false., lr_transformation=.true., index_model=2)
               else
                  call plane_boundary_interaction(1, sphere_cluster%sphere_order(i), &
                                                  rtran(1), rtran(2), sphere_cluster%sphere_position(3, i), rpos(3), &
                                        rvec, source_vector=sourcevec(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i), p), &
                                                  index_model=2, lr_transformation=.true., make_symmetric=.false.)
               end if
               evec(:, p) = evec(:, p) + matmul(vwf_0(:, :, 1), rvec(:, 1)) + matmul(vwf_0(:, :, 2), rvec(:, 2))
             hvec(:, p) = hvec(:, p) + (matmul(vwf_0(:, :, 1), rvec(:, 1)) - matmul(vwf_0(:, :, 2), rvec(:, 2))) * ri / (0.d0, 1.d0)
            end do
         end do
      end if
   end subroutine calculate_surface_field

   subroutine calculate_incident_field(rpos, layer, alpha, sinc, dir, cellinfo, evec, hvec, &
                                       dipole_position, dipole_moment, magnetic_segments, magnetic_quadrature_order)
      implicit none
      integer :: p, layer, dir, nodr, nblk
      integer, intent(in), optional :: magnetic_quadrature_order
      real(real64) :: alpha, sinc, rpos(3), rcell(3), rtran(3)
      real(real64), intent(in), optional :: dipole_position(3)
      complex(real64) :: riinc, pmnp(3, 2, 2), evec(3, 2), hvec(3, 2)
      complex(real64), intent(in), optional :: dipole_moment(3)
      type(magnetic_current_segment_t), intent(in), optional :: magnetic_segments(:)
      complex(real64), allocatable :: vwf(:, :, :)
      type(cell_info), pointer :: cellinfo
      riinc = layer_ref_index(layer)
      if (present(dipole_position) .and. present(dipole_moment)) then
         evec = (0.0_real64, 0.0_real64)
         hvec = (0.0_real64, 0.0_real64)
         call electric_dipole_field(rpos - dipole_position, riinc, dipole_moment, evec(:, 1), hvec(:, 1))
      elseif (present(magnetic_segments) .and. present(magnetic_quadrature_order)) then
         evec = (0.0_real64, 0.0_real64)
         hvec = (0.0_real64, 0.0_real64)
         call magnetic_current_field(rpos, riinc, magnetic_segments, magnetic_quadrature_order, &
                                     evec(:, 1), hvec(:, 1))
      elseif (incident_gb) then
         if (store_surface_vector) then
            nodr = cellinfo%order
            rcell = cellinfo%rcell
            nblk = nodr * (nodr + 2)
            allocate (vwf(3, nblk, 2))
!write(*,'(6f10.4)') rpos,rcell
            if (.not. associated(cellinfo%gb_vector)) then
!write(*,'('' allocating gb vec'' )')
               allocate (cellinfo%gb_vector(nblk, 2, 2))
               call layered_gaussian_beam_coefficients(alpha, sinc, dir, rcell, nodr, &
                                                       cellinfo%gb_vector, include_direct=.true., &
                                                       include_indirect=.true.)
            end if
            rtran(:) = rpos(:) - rcell(:)
            call vector_spherical_wave_functions(rtran, (/riinc, riinc/), nodr, 1, vwf, index_model=2)
            do p = 1, 2
               evec(:, p) = matmul(vwf(:, :, 1), cellinfo%gb_vector(1:nblk, 1, p)) &
                            + matmul(vwf(:, :, 2), cellinfo%gb_vector(1:nblk, 2, p))
               hvec(:, p) = (matmul(vwf(:, :, 1), cellinfo%gb_vector(1:nblk, 1, p)) &
                             - matmul(vwf(:, :, 2), cellinfo%gb_vector(1:nblk, 2, p))) * riinc / (0.d0, 1.d0)
            end do
            deallocate (vwf)
         else
            call layered_gaussian_beam_coefficients(alpha, sinc, dir, rcell, 1, pmnp, include_direct=.true., &
                                                    include_indirect=.true.)
            do p = 1, 2
               evec(:, p) = matmul(vwf_0(:, :, 1), pmnp(:, 1, p)) + matmul(vwf_0(:, :, 2), pmnp(:, 2, p))
               hvec(:, p) = (matmul(vwf_0(:, :, 1), pmnp(:, 1, p)) - matmul(vwf_0(:, :, 2), pmnp(:, 2, p))) * riinc / (0.d0, 1.d0)
            end do
         end if
      else
         call layer_plane_wave_coefficients(alpha, sinc, dir, rpos, 1, pmnp)
         do p = 1, 2
            evec(:, p) = matmul(vwf_0(:, :, 1), pmnp(:, 1, p)) + matmul(vwf_0(:, :, 2), pmnp(:, 2, p))
            hvec(:, p) = (matmul(vwf_0(:, :, 1), pmnp(:, 1, p)) - matmul(vwf_0(:, :, 2), pmnp(:, 2, p))) * riinc / (0.d0, 1.d0)
         end do
      end if
   end subroutine calculate_incident_field

   pure real(real64) function distance_to_segment(point, start_point, end_point)
      real(real64), intent(in) :: point(3), start_point(3), end_point(3)
      real(real64) :: delta(3), projection

      delta = end_point - start_point
      projection = dot_product(point - start_point, delta) / max(tiny(1.0_real64), sum(delta**2))
      projection = max(0.0_real64, min(1.0_real64, projection))
      distance_to_segment = sqrt(sum((point - start_point - projection * delta)**2))
   end function distance_to_segment

   subroutine calculate_stored_surface_vector(nodr, rc, sourcevec, storedvector)
      implicit none
      integer :: nodr, nblk, i, p
      real(real64) :: rc(3), rhovec(2)
      complex(real64) :: sourcevec(sphere_cluster%number_eqns, 2)
      complex(real64), allocatable :: vector(:, :, :)
      complex(real64), pointer :: storedvector(:, :, :)

      nblk = nodr * (nodr + 2)
      allocate (vector(nblk, 2, 2), storedvector(nblk, 2, 2))
      storedvector(1:nblk, 1:2, 1:2) = 0.d0
      do i = 1, sphere_cluster%number_spheres
         if (sphere_cluster%host_sphere(i) .ne. 0) cycle
         rhovec(1:2) = rc(1:2) - sphere_cluster%sphere_position(1:2, i)
         vector = 0.d0
         do p = 1, 2
            if (periodic_lattice) then
               call plane_boundary_lattice_interaction(nodr, sphere_cluster%sphere_order(i), &
                                                       rhovec(1), rhovec(2), rc(3), sphere_cluster%sphere_position(3, i), &
                                                       vector(:, :, p), &
 source_vector=sourcevec(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i), p), &
                                                       include_source=.false., lr_transformation=.true., index_model=2)
            else
               call plane_boundary_interaction(nodr, sphere_cluster%sphere_order(i), &
                                               rhovec(1), rhovec(2), sphere_cluster%sphere_position(3, i), rc(3), &
                                               vector(:, :, p), &
 source_vector=sourcevec(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i), p), &
                                               index_model=2, lr_transformation=.true., &
                                               make_symmetric=.false.)
            end if
         end do
         storedvector(1:nblk, 1:2, 1:2) = storedvector(1:nblk, 1:2, 1:2) + vector(1:nblk, 1:2, 1:2)
      end do
      deallocate (vector)
   end subroutine calculate_stored_surface_vector

   subroutine calculate_stored_source_vector(sourcevec, cellinfo)
      implicit none
      integer :: nodr, nblk, i, p, cellhost, vswf_type
      real(real64) :: rc(3), r
      complex(real64) :: sourcevec(sphere_cluster%number_eqns, 2), ri2(2)
      type(cell_info), pointer :: cellinfo
      type(translation_operator_state) :: tranmat

      cellhost = cellinfo%host
      if (cellhost .eq. 0) then
         ri2(1:2) = layer_ref_index(cellinfo%layer)
      else
         ri2 = sphere_cluster%sphere_ref_index(:, cellhost)
      end if
      nodr = cellinfo%order
      nblk = nodr * (nodr + 2)
      allocate (cellinfo%reg_source_vector(nblk, 2, 2))
      cellinfo%reg_source_vector(1:nblk, 1:2, 1:2) = 0.d0
      cellinfo%nispheres = 0
      cellinfo%outside_spheres = .false.
      do i = 1, sphere_cluster%number_spheres
         if (sphere_cluster%host_sphere(i) /= cellhost) cycle
         rc = cellinfo%rcell - sphere_cluster%sphere_position(:, i)
         if (sqrt(sum(rc**2)) <= 2.0_real64 * near_field_expansion_spacing) then
            cellinfo%nispheres = cellinfo%nispheres + 1
         end if
      end do
      if (allocated(cellinfo%sphere_indices)) deallocate (cellinfo%sphere_indices)
      allocate (cellinfo%sphere_indices(cellinfo%nispheres))
      cellinfo%nispheres = 0
      do i = 1, sphere_cluster%number_spheres
         if (sphere_cluster%host_sphere(i) .eq. cellhost .or. i .eq. cellhost) then
            rc = cellinfo%rcell(:) - sphere_cluster%sphere_position(:, i)
            r = sqrt(sum(rc**2))
            if (r .le. 2.d0 * near_field_expansion_spacing .and. sphere_cluster%host_sphere(i) .eq. cellhost) then
               cellinfo%nispheres = cellinfo%nispheres + 1
               cellinfo%sphere_indices(cellinfo%nispheres) = i
            else
               cellinfo%outside_spheres = .true.
               if (sphere_cluster%host_sphere(i) .eq. cellhost) then
                  vswf_type = 3
               else
                  vswf_type = 1
               end if
               call tranmat%configure(vswf_type, rc, ri2, &
                                      max(sphere_cluster%sphere_order(i), nodr) .ge. sphere_cluster%translation_switch_order)
               do p = 1, 2
                  if (sphere_cluster%host_sphere(i) .eq. cellhost) then
                     call tranmat%apply(sphere_cluster%sphere_order(i), 2, nodr, 2, &
               sourcevec(sphere_cluster%sphere_offset(i) + 1:sphere_cluster%sphere_offset(i) + sphere_cluster%sphere_block(i), p), &
                                        cellinfo%reg_source_vector(:, :, p))
                  else
                     call tranmat%apply(sphere_cluster%sphere_order(i), 2, nodr, 2, &
                                        internal_field_vector(cellhost)%vector(:, :, p), &
                                        cellinfo%reg_source_vector(:, :, p))
                  end if
               end do
               call tranmat%clear()
            end if
         end if
      end do
   end subroutine calculate_stored_source_vector

   subroutine initialize_grid_points(griddim, gridinfo)
      implicit none
      logical :: ingrid
      integer :: depth, i, l, layer, griddim(3), ix, iy, iz, nodr, nbound, zbsign, ncell(3), jy, jx, jlim(2)
      type(grid_info) :: gridinfo(griddim(1), griddim(2), griddim(3))
      real(real64) :: x, y, z, pbcellsize(1:max(1, number_plane_boundaries)), zbound, zbcell, rcell(3), cellsize, spos(3)
      type(cell_info) :: cellinfo
      type(linked_sphere_data), pointer :: slist

!         call clear_cell_info_list()
      total_cells = 0
      do iz = 1, griddim(3)
         do iy = 1, griddim(2)
            do ix = 1, griddim(1)
               gridinfo(ix, iy, iz)%initialized = .false.
               gridinfo(ix, iy, iz)%cellnum = 0
            end do
         end do
      end do

      number_intersecting_spheres = 0
      allocate (intersecting_spheres)
      slist => intersecting_spheres
      if (periodic_lattice) then
         jlim = (/-1, 1/)
      else
         jlim = 0
      end if
      do depth = sphere_cluster%max_sphere_depth, 0, -1
         do i = 1, sphere_cluster%number_spheres
            if (sphere_cluster%sphere_depth(i) .ne. depth) cycle
            do jy = jlim(1), jlim(2)
               do jx = jlim(1), jlim(2)
                  spos = sphere_cluster%sphere_position(:, i)
                  if (periodic_lattice) &
                     spos(1:2) = spos(1:2) + dble((/jx, jy/)) * cell_width
                  call sphere_to_grid_points(i, spos, griddim, gridinfo, ingrid)
                  if (ingrid) then
                     slist%sphere = i
                     slist%host = sphere_cluster%host_sphere(i)
                     slist%position = spos
                     slist%radius = sphere_cluster%sphere_radius(i)
                     number_intersecting_spheres = number_intersecting_spheres + 1
                     allocate (slist%next)
                     slist => slist%next
                  end if
               end do
            end do
         end do
      end do
      do l = 1, max(1, number_plane_boundaries)
         pbcellsize(l) = near_field_expansion_spacing
         if (plane_surface_present) then
            do i = 1, sphere_cluster%number_spheres
               if ((sphere_cluster%sphere_layer(i) .eq. l - 1 .or. sphere_cluster%sphere_layer(i) .eq. l) .and. sphere_cluster%host_sphere(i) .eq. 0) then
                  pbcellsize(l) = min(pbcellsize(l), 1.d0 * abs(plane_boundary_position(l) - sphere_cluster%sphere_position(3, i)))
               end if
            end do
         end if
      end do

      do iz = 1, griddim(3)
         z = (dble(iz) - 0.5d0) * grid_spacing(3) + grid_region(3, 1)
         layer = find_layer_index(z)
         if (plane_surface_present) then
            if (layer .eq. 0) then
               nbound = 1
               zbound = plane_boundary_position(1) - z
               zbcell = pbcellsize(1)
            elseif (layer .eq. number_plane_boundaries) then
               nbound = number_plane_boundaries
               zbound = z - plane_boundary_position(number_plane_boundaries)
               zbcell = pbcellsize(number_plane_boundaries)
            else
               nbound = layer
               zbound = z - plane_boundary_position(layer)
               if (zbound .gt. plane_boundary_position(layer + 1) - z) then
                  zbound = plane_boundary_position(layer + 1) - z
                  nbound = layer + 1
               end if
               zbcell = min(pbcellsize(nbound), &
                            plane_boundary_position(layer + 1) - plane_boundary_position(layer))
            end if
            if (nbound .eq. layer) then
               zbsign = 1
            else
               zbsign = -1
            end if
            if (zbound .le. pbcellsize(nbound)) then
               cellsize = pbcellsize(nbound)
               rcell(3) = plane_boundary_position(nbound) + 0.5d0 * zbsign * zbcell
            else
               cellsize = near_field_expansion_spacing
               zbound = plane_boundary_position(nbound) + zbsign * pbcellsize(nbound)
               rcell(3) = zbound &
                          + (floor(abs(z - zbound) / near_field_expansion_spacing) + 0.5d0) &
                          * zbsign * near_field_expansion_spacing
            end if
         else
            cellsize = near_field_expansion_spacing
            rcell(3) = (dble(floor((z - grid_region(3, 1)) / cellsize)) + 0.5d0) * cellsize + grid_region(3, 1)
         end if
         if (griddim(3) .eq. 1) then
            ncell(3) = 0
            rcell(3) = grid_region(3, 1)
         else
            ncell(3) = floor(rcell(3) / grid_spacing(3))
         end if
         nodr = max(ceiling(near_field_expansion_order * 0.9999d0 * (cellsize / near_field_expansion_spacing)), 1)
!            nodr=near_field_expansion_order

         do iy = 1, griddim(2)
            y = (dble(iy) - 0.5d0) * grid_spacing(2) + grid_region(2, 1)
            if (griddim(2) .eq. 1) then
               rcell(2) = grid_region(2, 1)
               ncell(2) = 0
            else
               rcell(2) = (dble(floor((y - grid_region(2, 1)) / cellsize)) + 0.5d0) * cellsize + grid_region(2, 1)
               ncell(2) = floor(rcell(2) / grid_spacing(2))
            end if
            do ix = 1, griddim(1)
               x = (dble(ix) - 0.5d0) * grid_spacing(1) + grid_region(1, 1)
               if (griddim(1) .eq. 1) then
                  rcell(1) = grid_region(1, 1)
                  ncell(1) = 0
               else
                  rcell(1) = (dble(floor((x - grid_region(1, 1)) / cellsize)) + 0.5d0) * cellsize + grid_region(1, 1)
                  ncell(1) = floor(rcell(1) / grid_spacing(1))
               end if
               if (.not. gridinfo(ix, iy, iz)%initialized) then
                  gridinfo(ix, iy, iz)%initialized = .true.
                  gridinfo(ix, iy, iz)%host = 0
                  gridinfo(ix, iy, iz)%layer = layer
               end if
               cellinfo%ncell = ncell
               cellinfo%rcell(:) = rcell
               cellinfo%order = nodr
               cellinfo%host = gridinfo(ix, iy, iz)%host
               cellinfo%layer = gridinfo(ix, iy, iz)%layer
               call find_or_add_cell_info(cellinfo, gridinfo(ix, iy, iz)%cellinfo, cell_info_list)
               gridinfo(ix, iy, iz)%cellnum = total_cells
!write(*,'(5i3,10es12.4)') ix,iz,ncell(1),ncell(3),total_cells,rcell(1),rcell(3)
!flush(6)
            end do
         end do
      end do

   end subroutine initialize_grid_points

   subroutine sphere_to_grid_points(sphere, spos, griddim, gridinfo, ingrid)
      implicit none
      logical :: ingrid, skip
      integer :: sphere, ic(3), i, limits(3, 2), ix, iy, iz, griddim(3), n1, n2
      real(real64) :: rpos(3), rc(3), r, r1, r2, spos(3)
      type(grid_info) :: gridinfo(griddim(1), griddim(2), griddim(3))
      type(cell_info) :: cellinfo

      ingrid = .false.
      skip = .false.
      do i = 1, 3
         if (grid_spacing(i) .eq. 0.d0) then
            if (abs(grid_region(i, 1) - spos(i)) .gt. sphere_cluster%sphere_radius(sphere)) then
               skip = .true.
            else
               limits(i, 1) = 1
               limits(i, 2) = 1
            end if
         else
            r1 = spos(i) - sphere_cluster%sphere_radius(sphere)
            if (r1 .gt. grid_region(i, 2)) then
               skip = .true.
            end if
            r2 = spos(i) + sphere_cluster%sphere_radius(sphere)
            if (r2 .lt. grid_region(i, 1)) then
               skip = .true.
            end if
            n1 = ceiling((r1 - grid_region(i, 1)) / grid_spacing(i))
            n2 = ceiling((r2 - grid_region(i, 1)) / grid_spacing(i))
            limits(i, 1) = max(1, n1)
            limits(i, 2) = min(griddim(i), n2)
         end if
      end do
      if (skip) return

      do iz = limits(3, 1), limits(3, 2)
         do iy = limits(2, 1), limits(2, 2)
            do ix = limits(1, 1), limits(1, 2)
               ic = (/ix, iy, iz/)
               rpos(:) = (dble(ic(:)) - (/0.5d0, 0.5d0, 0.5d0/)) * grid_spacing(:) + grid_region(:, 1)
               rc(:) = rpos(:) - spos(:)
               r = sqrt(dot_product(rc(:), rc(:)))
               if (r .gt. sphere_cluster%sphere_radius(sphere)) cycle
               if (gridinfo(ic(1), ic(2), ic(3))%initialized) cycle
               ingrid = .true.
               gridinfo(ic(1), ic(2), ic(3))%host = sphere
               gridinfo(ic(1), ic(2), ic(3))%layer = sphere_cluster%sphere_layer(sphere)
               gridinfo(ic(1), ic(2), ic(3))%initialized = .true.

!!                  cellinfo%ncell=(/0,0,0/)
!!                  cellinfo%rcell(:)=sphere_cluster%sphere_position(:,sphere)
!!                  cellinfo%order=sphere_cluster%sphere_order(sphere)
!                  cellinfo%host=sphere
!                  cellinfo%order=near_field_expansion_order
!                  cellinfo%layer=sphere_cluster%sphere_layer(sphere)
!                  call find_or_add_cell_info(cellinfo,gridinfo(ic(1),ic(2),ic(3))%cellinfo,cell_info_list)
!                  gridinfo(ic(1),ic(2),ic(3))%cellnum=total_cells
            end do
         end do
      end do
   end subroutine sphere_to_grid_points

   subroutine find_or_add_cell_info(delem, elem, list)
      implicit none
      logical :: inlist
      type(linked_cell_list), pointer :: list, tlist
      type(cell_info), pointer :: elem
      type(cell_info) :: delem

      if (associated(list)) then
         tlist => list
         inlist = .true.
         do while (delem%ncell(1) .ne. tlist%cellinfo%ncell(1) &
                   .or. delem%ncell(2) .ne. tlist%cellinfo%ncell(2) &
                   .or. delem%ncell(3) .ne. tlist%cellinfo%ncell(3) &
                   .or. delem%host .ne. tlist%cellinfo%host &
                   .or. delem%layer .ne. tlist%cellinfo%layer &
                   .or. delem%order .ne. tlist%cellinfo%order)
            if (.not. associated(tlist%next)) then
               inlist = .false.
               exit
            else
               tlist => tlist%next
            end if
         end do
         if (inlist) then
            elem => tlist%cellinfo
         else
            total_cells = total_cells + 1
            allocate (tlist%next)
            tlist => tlist%next
            tlist%cellinfo%ncell(:) = delem%ncell(:)
            tlist%cellinfo%rcell(:) = delem%rcell(:)
            tlist%cellinfo%host = delem%host
            tlist%cellinfo%layer = delem%layer
            tlist%cellinfo%order = delem%order
            elem => tlist%cellinfo
         end if
      else
         total_cells = 1
         allocate (list)
         list%cellinfo%ncell(:) = delem%ncell(:)
         list%cellinfo%rcell(:) = delem%rcell(:)
         list%cellinfo%host = delem%host
         list%cellinfo%layer = delem%layer
         list%cellinfo%order = delem%order
         elem => list%cellinfo
      end if
   end subroutine find_or_add_cell_info

   subroutine write_near_field_output_header(griddim, outputunit, print_intersecting_spheres)
      implicit none
      logical :: pis
      logical, optional :: print_intersecting_spheres
      integer :: griddim(3), outputunit, n, j, l1, l2
      type(linked_sphere_data), pointer :: slist
      if (present(print_intersecting_spheres)) then
         pis = print_intersecting_spheres
      else
         pis = .true.
      end if
      write (outputunit, '('' run number:'')')
      write (outputunit, '(i5)') local_run_number
      local_run_number = local_run_number + 1
      if (pis) then
         n = number_intersecting_spheres
      else
         n = 0
      end if
      write (outputunit, '(i5)') n
      slist => intersecting_spheres
      do j = 1, n
         write (outputunit, '(4es12.4)') slist%position, slist%radius
         if (j .lt. n) slist => slist%next
      end do
      l1 = find_layer_index(grid_region(3, 1))
      l2 = find_layer_index(grid_region(3, 2))
      write (outputunit, '(i5)') l2 - l1
      do j = l1 + 1, l2
         write (outputunit, '(es12.4)') plane_boundary_position(j)
      end do
      write (outputunit, '(3es12.4)') grid_region(:, 1)
      write (outputunit, '(3es12.4)') grid_region(:, 2)
      write (outputunit, '(3i5)') griddim(:)
   end subroutine write_near_field_output_header

end module near_field
