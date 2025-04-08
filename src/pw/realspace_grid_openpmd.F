!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Generate Gaussian cube files
! **************************************************************************************************
MODULE realspace_grid_openpmd
   USE cp_files, ONLY: close_file, &
                       open_file
   USE cp_log_handling, ONLY: cp_logger_get_default_io_unit
   USE cp_output_handling_openpmd, ONLY: cp_openpmd_get_value_unit_nr, &
                                         cp_openpmd_per_call_value_type
   USE kinds, ONLY: default_string_length, &
                    dp
   USE message_passing, ONLY: &
      file_amode_rdonly, file_offset, mp_comm_type, mp_file_descriptor_type, mp_file_type, &
      mp_file_type_free, mp_file_type_hindexed_make_chv, mp_file_type_set_view_chv, &
      mpi_character_size
! TODO deal properly with MPI
#if defined(__parallel)
#if defined(__MPI_F08)
   USE mpi_f08, ONLY: mpi_allreduce, mpi_integer, mpi_bxor, mpi_comm, mpi_allgather
#else
   USE mpi, ONLY: mpi_allreduce, mpi_integer, mpi_bxor, mpi_comm, mpi_allgather
#endif
#endif

   USE openpmd_api, ONLY: &
      openpmd_dynamic_memory_view_resolve_1d_double, &
      openpmd_dynamic_memory_view_resolve_3d_double, openpmd_dynamic_memory_view_type_1d, &
      openpmd_dynamic_memory_view_type_3d, openpmd_iteration_get_mesh, &
      openpmd_iteration_get_particle_species, openpmd_mesh_as_record_component, &
      openpmd_mesh_set_axis_labels, openpmd_mesh_set_grid_global_offset, &
      openpmd_mesh_set_grid_spacing, openpmd_mesh_set_position, openpmd_mesh_type, &
      openpmd_particle_species_get_record, openpmd_particle_species_type, &
      openpmd_record_component_make_constant_zero, openpmd_record_component_reset_dataset, &
      openpmd_record_component_store_chunk_span_1d_double, &
      openpmd_record_component_store_chunk_span_3d_double, openpmd_record_component_type, &
      openpmd_record_get_component, openpmd_record_type, openpmd_type_double, openpmd_type_int, openpmd_iteration_as_attributable, openpmd_attributable_series_flush
   USE pw_grid_types, ONLY: PW_MODE_LOCAL
   USE pw_types, ONLY: pw_r3d_rs_type
   USE util, ONLY: sort_unique
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pw_to_cube, pw_to_openpmd, cube_to_pw, pw_to_simple_volumetric

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'realspace_grid_openpmd'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   LOGICAL, PRIVATE                     :: parses_linebreaks = .FALSE., &
                                           parse_test = .TRUE.

   TYPE cp_openpmd_write_buffer_1d
      REAL(KIND=dp), POINTER :: buffer(:)
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param particles_z ...
!> \param res_atom_types ...
!> \param res_atom_counts ...
!> \param res_len ...
! **************************************************************************************************
   SUBROUTINE pw_get_atom_types(particles_z, res_atom_types, res_atom_counts, res_len)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: particles_z
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: res_atom_types, res_atom_counts
      INTEGER, INTENT(OUT)                               :: res_len

      INTEGER                                            :: current_atom_number, i
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: particles_z_sorted
      LOGICAL                                            :: unique

      ALLOCATE (particles_z_sorted(SIZE(particles_z)))
      particles_z_sorted(:) = particles_z(:)
      CALL sort_unique(particles_z_sorted, unique)

      ALLOCATE (res_atom_types(MIN(118, SIZE(particles_z))))
      ALLOCATE (res_atom_counts(MIN(118, SIZE(particles_z))))
      current_atom_number = -1
      res_len = 0
      DO i = 1, SIZE(particles_z_sorted)
         IF (particles_z_sorted(i) /= current_atom_number) THEN
            res_len = res_len + 1
            current_atom_number = particles_z_sorted(i)
            res_atom_types(res_len) = current_atom_number
            res_atom_counts(res_len) = 1
         ELSE
            res_atom_counts(res_len) = res_atom_counts(res_len) + 1
         END IF
      END DO

   END SUBROUTINE pw_get_atom_types

! **************************************************************************************************
!> \brief ...
!> \param particles_z ...
!> \param particles_r ...
!> \param atom_type ...
!> \param atom_count ...
!> \param openpmd_data ...
!> \param do_write_data ...
! **************************************************************************************************
   SUBROUTINE pw_write_particle_species(particles_z, particles_r, atom_type, atom_count, openpmd_data, do_write_data)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: particles_z
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: particles_r
      INTEGER, INTENT(IN)                                :: atom_type, atom_count
      TYPE(cp_openpmd_per_call_value_type)               :: openpmd_data
      LOGICAL                                            :: do_write_data

      CHARACTER(len=1), DIMENSION(3), PARAMETER          :: dims = (/"x", "y", "z"/)

      CHARACTER(len=3)                                   :: atom_type_as_string
      CHARACTER(len=default_string_length)               :: species_name
      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(1)                              :: global_extent, global_offset, &
                                                            local_extent
      TYPE(cp_openpmd_write_buffer_1d), DIMENSION(3)     :: write_buffers
      TYPE(openpmd_dynamic_memory_view_type_1d), &
         DIMENSION(3)                                    :: unresolved_write_buffers
      TYPE(openpmd_particle_species_type)                :: species
      TYPE(openpmd_record_component_type)                :: position_component, &
                                                            position_offset_component
      TYPE(openpmd_record_type)                          :: position, position_offset

      global_extent(1) = atom_count
      IF (do_write_data) THEN
         global_offset(1) = 0
         local_extent(1) = atom_count
      ELSE
         global_offset(1) = 0
         local_extent(1) = 0
      END IF

      WRITE (atom_type_as_string, '(I3)') atom_type
      species_name = TRIM(openpmd_data%name_prefix)//"-"//ADJUSTL(atom_type_as_string)

      species = openpmd_iteration_get_particle_species(openpmd_data%iteration, TRIM(species_name))

      position_offset = openpmd_particle_species_get_record(species, "positionOffset")
      position = openpmd_particle_species_get_record(species, "position")
      DO k = 1, SIZE(dims)
         position_offset_component = openpmd_record_get_component(position_offset, dims(k))
         CALL openpmd_record_component_make_constant_zero(position_offset_component, openpmd_type_int, global_extent)
         position_component = openpmd_record_get_component(position, dims(k))
         CALL openpmd_record_component_reset_dataset(position_component, openpmd_type_double, global_extent)
         unresolved_write_buffers(k) = &
            openpmd_record_component_store_chunk_span_1d_double(position_component, global_offset, local_extent)
         IF (do_write_data) THEN
            write_buffers(k)%buffer => &
               openpmd_dynamic_memory_view_resolve_1d_double(unresolved_write_buffers(k), DEALLOCATE=.FALSE.)
         END IF
      END DO

      IF (do_write_data) THEN
         ! Resolve Spans for a second time to allow for internal reallocations in BP4 engine of ADIOS2
         DO k = 1, SIZE(dims)
            write_buffers(k)%buffer = &
               openpmd_dynamic_memory_view_resolve_1d_double(unresolved_write_buffers(k), DEALLOCATE=.TRUE.)
         END DO
         j = 1
         DO i = 1, SIZE(particles_z)
            IF (particles_z(i) == atom_type) THEN
               DO k = 1, 3
                  write_buffers(k)%buffer(j) = particles_r(k, i)
               END DO
               j = j + 1
            END IF
         END DO
      END IF
   END SUBROUTINE pw_write_particle_species

! **************************************************************************************************
!> \brief ...
!> \param particles_z ...
!> \param particles_r ...
!> \param atom_types ...
!> \param atom_counts ...
!> \param num_atom_types ...
!> \param openpmd_data ...
!> \param gid ...
! **************************************************************************************************
   SUBROUTINE pw_write_particles(particles_z, particles_r, atom_types, atom_counts, num_atom_types, openpmd_data, gid)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: particles_z
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: particles_r
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_types, atom_counts
      INTEGER, INTENT(IN), TARGET                        :: num_atom_types
      TYPE(cp_openpmd_per_call_value_type)               :: openpmd_data
      TYPE(mp_comm_type), OPTIONAL                       :: gid

      INTEGER                                            :: i, mpi_rank
      LOGICAL                                            :: do_write_data

      IF (PRESENT(gid)) THEN
         CALL gid%get_rank(mpi_rank)
         do_write_data = mpi_rank == 0
      ELSE
         do_write_data = .TRUE.
      END IF
      DO i = 1, num_atom_types
         CALL pw_write_particle_species(particles_z, particles_r, atom_types(i), atom_counts(i), openpmd_data, do_write_data)
      END DO
   END SUBROUTINE pw_write_particles

! TODO move this to own module
! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param unit_nr ...
!> \param title ...
!> \param particles_r ...
!> \param particles_z ...
!> \param stride ...
!> \param zero_tails ...
!> \param silent ...
!> \param mpi_io ...
! **************************************************************************************************
   SUBROUTINE pw_to_openpmd(pw, unit_nr, title, particles_r, particles_z, stride, zero_tails, &
                            silent, mpi_io)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw
      INTEGER                                            :: unit_nr
      CHARACTER(*), INTENT(IN), OPTIONAL                 :: title
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: particles_r
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: particles_z
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: stride
      LOGICAL, INTENT(IN), OPTIONAL                      :: zero_tails, silent, mpi_io

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_to_openpmd'
      INTEGER, PARAMETER                                 :: entry_len = 13, num_entries_line = 6

      INTEGER :: count1, count2, count3, handle, i, I1, I2, I3, iat, j, L1, L2, L3, msglen, &
         my_rank, my_stride(3), np, num_atom_types, num_linebreak, num_pe, size_of_z, U1, U2, U3
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_counts, atom_types
      INTEGER, DIMENSION(3)                              :: global_extent, local_extent, offset
      LOGICAL                                            :: be_silent, my_zero_tails, parallel_write
      REAL(KIND=dp), DIMENSION(3)                        :: grid_spacing
      REAL(KIND=dp), POINTER                             :: write_buffer(:, :, :)
      TYPE(cp_openpmd_per_call_value_type)               :: openpmd_data
      TYPE(mp_comm_type)                                 :: gid
      TYPE(openpmd_dynamic_memory_view_type_3d)          :: unresolved_write_buffer
      TYPE(openpmd_mesh_type)                            :: mesh
      TYPE(openpmd_record_component_type)                :: scalar_mesh

      CALL timeset(routineN, handle)

      my_zero_tails = .FALSE.
      be_silent = .FALSE.
      parallel_write = .FALSE.
      gid = pw%pw_grid%para%group
      IF (PRESENT(zero_tails)) my_zero_tails = zero_tails
      IF (PRESENT(silent)) be_silent = silent
      IF (PRESENT(mpi_io)) parallel_write = mpi_io
      my_stride = 1
      IF (PRESENT(stride)) THEN
         IF (SIZE(stride) /= 1 .AND. SIZE(stride) /= 3) &
            CALL cp_abort(__LOCATION__, "STRIDE keyword can accept only 1 "// &
                          "(the same for X,Y,Z) or 3 values. Correct your input file.")
         IF (SIZE(stride) == 1) THEN
            DO i = 1, 3
               my_stride(i) = stride(1)
            END DO
         ELSE
            my_stride = stride(1:3)
         END IF
         CPASSERT(my_stride(1) > 0)
         CPASSERT(my_stride(2) > 0)
         CPASSERT(my_stride(3) > 0)
      END IF

      openpmd_data = cp_openpmd_get_value_unit_nr(unit_nr)

      IF (.NOT. parallel_write .OR. parallel_write) THEN

         CPASSERT(PRESENT(particles_z) .EQV. PRESENT(particles_r))
         np = 0
         IF (PRESENT(particles_z)) THEN
            CALL pw_get_atom_types(particles_z, atom_types, atom_counts, num_atom_types)
            CPASSERT(SIZE(particles_z) == SIZE(particles_r, dim=2))
            np = SIZE(particles_z)
         END IF

         DO i = 1, 3
            DO j = 1, 3
               IF (i == j) THEN
                  grid_spacing(i) = pw%pw_grid%dh(i, i)*REAL(my_stride(i), dp)
               ELSE
                  CPASSERT(pw%pw_grid%dh(i, j) == 0)
               END IF
            END DO
         END DO

         IF (PRESENT(particles_z)) THEN
            IF (parallel_write) THEN
               CALL pw_write_particles(particles_z, particles_r, atom_types, atom_counts, num_atom_types, openpmd_data, gid)
            ELSE
               CALL pw_write_particles(particles_z, particles_r, atom_types, atom_counts, num_atom_types, openpmd_data)
            END IF
         END IF

         DO iat = 1, 3
            global_extent(iat) = (pw%pw_grid%npts(iat) + my_stride(iat) - 1)/my_stride(iat)
            ! '+ 1' because upper end is inclusive, '- 1' for upper gaussian bracket
            offset(iat) = ((pw%pw_grid%bounds_local(1, iat) - pw%pw_grid%bounds(1, iat) + my_stride(iat) - 1)/my_stride(iat))
            ! refer local_extent to the global offset first in order to have consistent rounding
         local_extent(iat) = ((pw%pw_grid%bounds_local(2, iat) + 1 - pw%pw_grid%bounds(1, iat) + my_stride(iat) - 1)/my_stride(iat))
         END DO
         local_extent = local_extent - offset

         mesh = openpmd_iteration_get_mesh(openpmd_data%iteration, TRIM(openpmd_data%name_prefix))
         CALL openpmd_mesh_set_axis_labels(mesh, (/"x", "y", "z"/))
         CALL openpmd_mesh_set_position(mesh, (/0.5_dp, 0.5_dp, 0.5_dp/))
         CALL openpmd_mesh_set_grid_global_offset(mesh, pw%pw_grid%bounds(1, :)*grid_spacing/my_stride)
         CALL openpmd_mesh_set_grid_spacing(mesh, grid_spacing)
         scalar_mesh = openpmd_mesh_as_record_component(mesh)
         CALL openpmd_record_component_reset_dataset(scalar_mesh, openpmd_type_double, global_extent)

         ! shortcut
         L1 = pw%pw_grid%bounds_local(1, 1)
         L2 = pw%pw_grid%bounds_local(1, 2)
         L3 = pw%pw_grid%bounds_local(1, 3)
         U1 = pw%pw_grid%bounds_local(2, 1)
         U2 = pw%pw_grid%bounds_local(2, 2)
         U3 = pw%pw_grid%bounds_local(2, 3)

         my_rank = pw%pw_grid%para%group%mepos
         num_pe = pw%pw_grid%para%group%num_pe

         count3 = 0
         DO I3 = L3, U3, my_stride(3)
            ! maybe add an overload to provide `buf` here for HDF5, might have better performance
            ! for intermittent flushing
            ! or just call the buffer in the outer function if memory is no problem...
            unresolved_write_buffer = openpmd_record_component_store_chunk_span_3d_double( &
                                      scalar_mesh, &
                                      (/offset(1), offset(2), offset(3) + count3/), &
                                      (/local_extent(1), local_extent(2), 1/))
            write_buffer => openpmd_dynamic_memory_view_resolve_3d_double(unresolved_write_buffer, DEALLOCATE=.TRUE.)

            count2 = 0
            DO I2 = L2, U2, my_stride(2)
               ! This loop deals with ray (:, count2, count3) of the local subspace
               ! The write buffer itself has been allocated for slice (:, :, count3)
               count1 = 0
               DO I1 = L1, U1, my_stride(1)
                  write_buffer(count1 + 1, count2 + 1, 1) = pw%array(I1, I2, I3)
                  count1 = count1 + 1
               END DO
               count2 = count2 + 1
            END DO
            count3 = count3 + 1
         END DO
      ELSE
         size_of_z = CEILING(REAL(pw%pw_grid%bounds(2, 3) - pw%pw_grid%bounds(1, 3) + 1, dp)/REAL(my_stride(3), dp))
         num_linebreak = size_of_z/num_entries_line
         IF (MODULO(size_of_z, num_entries_line) /= 0) &
            num_linebreak = num_linebreak + 1
         msglen = (size_of_z*entry_len + num_linebreak)*mpi_character_size
         ! CALL mp_unit%set_handle(unit_nr)
         ! CALL pw_to_cube_parallel(pw, mp_unit, title, particles_r, particles_z, my_stride, my_zero_tails, msglen)
      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_to_openpmd

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param unit_nr ...
!> \param title ...
!> \param particles_r ...
!> \param particles_z ...
!> \param stride ...
!> \param zero_tails ...
!> \param silent ...
!> \param mpi_io ...
! **************************************************************************************************
   SUBROUTINE pw_to_cube(pw, unit_nr, title, particles_r, particles_z, stride, zero_tails, &
                         silent, mpi_io)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(*), INTENT(IN), OPTIONAL                 :: title
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: particles_r
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: particles_z
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: stride
      LOGICAL, INTENT(IN), OPTIONAL                      :: zero_tails, silent, mpi_io

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_to_cube'
      INTEGER, PARAMETER                                 :: entry_len = 13, num_entries_line = 6

      INTEGER :: checksum, dest, handle, i, I1, I2, I3, iat, ip, L1, L2, L3, msglen, my_rank, &
         my_stride(3), np, num_linebreak, num_pe, rank(2), size_of_z, source, tag, U1, U2, U3
      LOGICAL                                            :: be_silent, my_zero_tails, parallel_write
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: buf
      TYPE(mp_comm_type)                                 :: gid
      TYPE(mp_file_type)                                 :: mp_unit

      CALL timeset(routineN, handle)

      my_zero_tails = .FALSE.
      be_silent = .FALSE.
      parallel_write = .FALSE.
      IF (PRESENT(zero_tails)) my_zero_tails = zero_tails
      IF (PRESENT(silent)) be_silent = silent
      IF (PRESENT(mpi_io)) parallel_write = mpi_io
      my_stride = 1
      IF (PRESENT(stride)) THEN
         IF (SIZE(stride) /= 1 .AND. SIZE(stride) /= 3) &
            CALL cp_abort(__LOCATION__, "STRIDE keyword can accept only 1 "// &
                          "(the same for X,Y,Z) or 3 values. Correct your input file.")
         IF (SIZE(stride) == 1) THEN
            DO i = 1, 3
               my_stride(i) = stride(1)
            END DO
         ELSE
            my_stride = stride(1:3)
         END IF
         CPASSERT(my_stride(1) > 0)
         CPASSERT(my_stride(2) > 0)
         CPASSERT(my_stride(3) > 0)
      END IF

      IF (.NOT. parallel_write) THEN
         IF (unit_nr > 0) THEN
            ! this format seems to work for e.g. molekel and gOpenmol
            ! latest version of VMD can read non orthorhombic cells
            WRITE (unit_nr, '(a11)') "-Quickstep-"
            IF (PRESENT(title)) THEN
               WRITE (unit_nr, *) TRIM(title)
            ELSE
               WRITE (unit_nr, *) "No Title"
            END IF

            CPASSERT(PRESENT(particles_z) .EQV. PRESENT(particles_r))
            np = 0
            IF (PRESENT(particles_z)) THEN
               CPASSERT(SIZE(particles_z) == SIZE(particles_r, dim=2))
               ! cube files can only be written for 99999 particles due to a format limitation (I5)
               ! so we limit the number of particles written.
               np = MIN(99999, SIZE(particles_z))
            END IF

            WRITE (unit_nr, '(I5,3f12.6)') np, 0.0_dp, 0._dp, 0._dp !start of cube

            WRITE (unit_nr, '(I5,3f12.6)') (pw%pw_grid%npts(1) + my_stride(1) - 1)/my_stride(1), &
               pw%pw_grid%dh(1, 1)*REAL(my_stride(1), dp), pw%pw_grid%dh(2, 1)*REAL(my_stride(1), dp), &
               pw%pw_grid%dh(3, 1)*REAL(my_stride(1), dp)
            WRITE (unit_nr, '(I5,3f12.6)') (pw%pw_grid%npts(2) + my_stride(2) - 1)/my_stride(2), &
               pw%pw_grid%dh(1, 2)*REAL(my_stride(2), dp), pw%pw_grid%dh(2, 2)*REAL(my_stride(2), dp), &
               pw%pw_grid%dh(3, 2)*REAL(my_stride(2), dp)
            WRITE (unit_nr, '(I5,3f12.6)') (pw%pw_grid%npts(3) + my_stride(3) - 1)/my_stride(3), &
               pw%pw_grid%dh(1, 3)*REAL(my_stride(3), dp), pw%pw_grid%dh(2, 3)*REAL(my_stride(3), dp), &
               pw%pw_grid%dh(3, 3)*REAL(my_stride(3), dp)

            IF (PRESENT(particles_z)) THEN
               DO iat = 1, np
                  WRITE (unit_nr, '(I5,4f12.6)') particles_z(iat), 0._dp, particles_r(:, iat)
               END DO
            END IF
         END IF

         ! shortcut
         L1 = pw%pw_grid%bounds(1, 1)
         L2 = pw%pw_grid%bounds(1, 2)
         L3 = pw%pw_grid%bounds(1, 3)
         U1 = pw%pw_grid%bounds(2, 1)
         U2 = pw%pw_grid%bounds(2, 2)
         U3 = pw%pw_grid%bounds(2, 3)

         ALLOCATE (buf(L3:U3))

         my_rank = pw%pw_grid%para%group%mepos
         gid = pw%pw_grid%para%group
         num_pe = pw%pw_grid%para%group%num_pe
         tag = 1

         rank (1) = unit_nr
         rank (2) = my_rank
         checksum = 0
         IF (unit_nr > 0) checksum = 1

         CALL gid%sum(checksum)
         CPASSERT(checksum == 1)

         CALL gid%maxloc(rank)
         CPASSERT(rank(1) > 0)

         dest = rank(2)
         DO I1 = L1, U1, my_stride(1)
            DO I2 = L2, U2, my_stride(2)

               ! cycling through the CPUs, check if the current ray (I1,I2) is local to that CPU
               IF (pw%pw_grid%para%mode .NE. PW_MODE_LOCAL) THEN
                  DO ip = 0, num_pe - 1
                     IF (pw%pw_grid%para%bo(1, 1, ip, 1) <= I1 - L1 + 1 .AND. pw%pw_grid%para%bo(2, 1, ip, 1) >= I1 - L1 + 1 .AND. &
                         pw%pw_grid%para%bo(1, 2, ip, 1) <= I2 - L2 + 1 .AND. pw%pw_grid%para%bo(2, 2, ip, 1) >= I2 - L2 + 1) THEN
                        source = ip
                     END IF
                  END DO
               ELSE
                  source = dest
               END IF

               IF (source == dest) THEN
                  IF (my_rank == source) THEN
                     buf(:) = pw%array(I1, I2, :)
                  END IF
               ELSE
                  IF (my_rank == source) THEN
                     buf(:) = pw%array(I1, I2, :)
                     CALL gid%send(buf, dest, tag)
                  END IF
                  IF (my_rank == dest) THEN
                     CALL gid%recv(buf, source, tag)
                  END IF
               END IF

               IF (unit_nr > 0) THEN
                  IF (my_zero_tails) THEN
                     DO I3 = L3, U3
                        IF (buf(I3) < 1.E-7_dp) buf(I3) = 0.0_dp
                     END DO
                  END IF
                  WRITE (unit_nr, '(6E13.5)') (buf(I3), I3=L3, U3, my_stride(3))
               END IF

               ! this double loop generates so many messages that it can overload
               ! the message passing system, e.g. on XT3
               ! we therefore put a barrier here that limits the amount of message
               ! that flies around at any given time.
               ! if ever this routine becomes a bottleneck, we should go for a
               ! more complicated rewrite
               CALL gid%sync()

            END DO
         END DO

         DEALLOCATE (buf)
      ELSE
         size_of_z = CEILING(REAL(pw%pw_grid%bounds(2, 3) - pw%pw_grid%bounds(1, 3) + 1, dp)/REAL(my_stride(3), dp))
         num_linebreak = size_of_z/num_entries_line
         IF (MODULO(size_of_z, num_entries_line) /= 0) &
            num_linebreak = num_linebreak + 1
         msglen = (size_of_z*entry_len + num_linebreak)*mpi_character_size
         CALL mp_unit%set_handle(unit_nr)
         CALL pw_to_cube_parallel(pw, mp_unit, title, particles_r, particles_z, my_stride, my_zero_tails, msglen)
      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_to_cube

! **************************************************************************************************
!> \brief  Computes the external density on the grid
!>         hacked from external_read_density
!> \param grid     pw to read from cube file
!> \param filename name of cube file
!> \param scaling  scale values before storing
!> \param parallel_read ...
!> \param silent ...
!> \par History
!>      Created [M.Watkins] (01.2014)
!>      Use blocking, collective MPI read for parallel simulations [Nico Holmberg] (05.2017)
! **************************************************************************************************
   SUBROUTINE cube_to_pw(grid, filename, scaling, parallel_read, silent)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: grid
      CHARACTER(len=*), INTENT(in)                       :: filename
      REAL(kind=dp), INTENT(in)                          :: scaling
      LOGICAL, INTENT(in)                                :: parallel_read
      LOGICAL, INTENT(in), OPTIONAL                      :: silent

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cube_to_pw'
      INTEGER, PARAMETER                                 :: entry_len = 13, num_entries_line = 6

      INTEGER                                            :: extunit, handle, i, j, k, msglen, &
                                                            my_rank, nat, ndum, num_linebreak, &
                                                            num_pe, output_unit, size_of_z, tag
      INTEGER, DIMENSION(3)                              :: lbounds, lbounds_local, npoints, &
                                                            npoints_local, ubounds, ubounds_local
      LOGICAL                                            :: be_silent
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: buffer
      REAL(kind=dp), DIMENSION(3)                        :: dr, rdum
      TYPE(mp_comm_type)                                 :: gid

      output_unit = cp_logger_get_default_io_unit()

      CALL timeset(routineN, handle)

      be_silent = .FALSE.
      IF (PRESENT(silent)) THEN
         be_silent = silent
      END IF
      !get rs grids and parallel environment
      gid = grid%pw_grid%para%group
      my_rank = grid%pw_grid%para%group%mepos
      num_pe = grid%pw_grid%para%group%num_pe
      tag = 1

      lbounds_local = grid%pw_grid%bounds_local(1, :)
      ubounds_local = grid%pw_grid%bounds_local(2, :)
      size_of_z = ubounds_local(3) - lbounds_local(3) + 1

      IF (.NOT. parallel_read) THEN
         npoints = grid%pw_grid%npts
         lbounds = grid%pw_grid%bounds(1, :)
         ubounds = grid%pw_grid%bounds(2, :)

         DO i = 1, 3
            dr(i) = grid%pw_grid%dh(i, i)
         END DO

         npoints_local = grid%pw_grid%npts_local
         !pw grids at most pencils - all processors have a full set of z data for x,y
         ALLOCATE (buffer(lbounds(3):ubounds(3)))

         IF (my_rank == 0) THEN
            IF (output_unit > 0 .AND. .NOT. be_silent) THEN
               WRITE (output_unit, FMT="(/,T2,A,/,/,T2,A,/)") "Reading the cube file:     ", TRIM(filename)
            END IF

            CALL open_file(file_name=filename, &
                           file_status="OLD", &
                           file_form="FORMATTED", &
                           file_action="READ", &
                           unit_number=extunit)

            !skip header comments
            DO i = 1, 2
               READ (extunit, *)
            END DO
            READ (extunit, *) nat, rdum
            DO i = 1, 3
               READ (extunit, *) ndum, rdum
               IF ((ndum /= npoints(i) .OR. (ABS(rdum(i) - dr(i)) > 1e-4)) .AND. &
                   output_unit > 0) THEN
                  WRITE (output_unit, *) "Restart from density | ERROR! | CUBE FILE NOT COINCIDENT WITH INTERNAL GRID ", i
                  WRITE (output_unit, *) "Restart from density | ", ndum, " DIFFERS FROM ", npoints(i)
                  WRITE (output_unit, *) "Restart from density | ", rdum, " DIFFERS FROM ", dr(i)
               END IF
            END DO
            !ignore atomic position data - read from coord or topology instead
            DO i = 1, nat
               READ (extunit, *)
            END DO
         END IF

         !master sends all data to everyone
         DO i = lbounds(1), ubounds(1)
            DO j = lbounds(2), ubounds(2)
               IF (my_rank .EQ. 0) THEN
                  READ (extunit, *) (buffer(k), k=lbounds(3), ubounds(3))
               END IF
               CALL gid%bcast(buffer(lbounds(3):ubounds(3)), 0)

               !only use data that is local to me - i.e. in slice of pencil I own
               IF ((lbounds_local(1) .LE. i) .AND. (i .LE. ubounds_local(1)) .AND. (lbounds_local(2) .LE. j) &
                   .AND. (j .LE. ubounds_local(2))) THEN
                  !allow scaling of external potential values by factor 'scaling' (SCALING_FACTOR in input file)
                  grid%array(i, j, lbounds(3):ubounds(3)) = buffer(lbounds(3):ubounds(3))*scaling
               END IF

            END DO
         END DO

         IF (my_rank == 0) CALL close_file(unit_number=extunit)

         CALL gid%sync()
      ELSE
         ! Parallel routine needs as input the byte size of each grid z-slice
         ! This is a hack to prevent compilation errors with gcc -Wall (up to versions 6.3)
         ! related to allocatable-length string declaration CHARACTER(LEN=:), ALLOCATABLE, DIMENSION(:) :: string
         ! Each data line of a Gaussian cube contains max 6 entries with last line potentially containing less if nz % 6 /= 0
         ! Thus, this size is simply the number of entries multiplied by the entry size + the number of line breaks
         num_linebreak = size_of_z/num_entries_line
         IF (MODULO(size_of_z, num_entries_line) /= 0) &
            num_linebreak = num_linebreak + 1
         msglen = (size_of_z*entry_len + num_linebreak)*mpi_character_size
         CALL cube_to_pw_parallel(grid, filename, scaling, msglen, silent=silent)
      END IF

      CALL timestop(handle)

   END SUBROUTINE cube_to_pw

! **************************************************************************************************
!> \brief Reads a realspace potential/density from a cube file using collective MPI I/O and
!>        stores it in grid.
!> \param grid     pw to read from cube file
!> \param filename name of cube file
!> \param scaling  scale values before storing
!> \param msglen   the size of each grid slice along z-axis in bytes
!> \param silent ...
!> \par History
!>      Created [Nico Holmberg] (05.2017)
! **************************************************************************************************
   SUBROUTINE cube_to_pw_parallel(grid, filename, scaling, msglen, silent)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: grid
      CHARACTER(len=*), INTENT(in)                       :: filename
      REAL(kind=dp), INTENT(in)                          :: scaling
      INTEGER, INTENT(in)                                :: msglen
      LOGICAL, INTENT(in), OPTIONAL                      :: silent

      INTEGER, PARAMETER                                 :: entry_len = 13, num_entries_line = 6

      INTEGER, DIMENSION(3)                              :: lbounds, lbounds_local, npoints, &
                                                            npoints_local, ubounds, ubounds_local
      INTEGER, ALLOCATABLE, DIMENSION(:), TARGET         :: blocklengths
      INTEGER(kind=file_offset), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: displacements
      INTEGER(kind=file_offset)                          :: BOF
      INTEGER :: extunit_handle, i, ientry, islice, j, k, my_rank, nat, ndum, nslices, num_pe, &
         offset_global, output_unit, pos, readstat, size_of_z, tag
      CHARACTER(LEN=msglen), ALLOCATABLE, DIMENSION(:)   :: readbuffer
      CHARACTER(LEN=msglen)                              :: tmp
      LOGICAL                                            :: be_silent, should_read(2)
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: buffer
      REAL(kind=dp), DIMENSION(3)                        :: dr, rdum
      TYPE(mp_comm_type)                                 :: gid
      TYPE(mp_file_descriptor_type)                      :: mp_file_desc
      TYPE(mp_file_type)                                 :: extunit

      output_unit = cp_logger_get_default_io_unit()

      be_silent = .FALSE.
      IF (PRESENT(silent)) THEN
         be_silent = silent
      END IF

      !get rs grids and parallel envnment
      gid = grid%pw_grid%para%group
      my_rank = grid%pw_grid%para%group%mepos
      num_pe = grid%pw_grid%para%group%num_pe
      tag = 1

      DO i = 1, 3
         dr(i) = grid%pw_grid%dh(i, i)
      END DO

      npoints = grid%pw_grid%npts
      lbounds = grid%pw_grid%bounds(1, :)
      ubounds = grid%pw_grid%bounds(2, :)

      npoints_local = grid%pw_grid%npts_local
      lbounds_local = grid%pw_grid%bounds_local(1, :)
      ubounds_local = grid%pw_grid%bounds_local(2, :)
      size_of_z = ubounds_local(3) - lbounds_local(3) + 1
      nslices = (ubounds_local(1) - lbounds_local(1) + 1)*(ubounds_local(2) - lbounds_local(2) + 1)
      islice = 1

      ! Read header information and determine byte offset of cube data on master process
      IF (my_rank == 0) THEN
         IF (output_unit > 0 .AND. .NOT. be_silent) THEN
            WRITE (output_unit, FMT="(/,T2,A,/,/,T2,A,/)") "Reading the cube file:     ", TRIM(filename)
         END IF

         CALL open_file(file_name=filename, &
                        file_status="OLD", &
                        file_form="FORMATTED", &
                        file_action="READ", &
                        file_access="STREAM", &
                        unit_number=extunit_handle)

         !skip header comments
         DO i = 1, 2
            READ (extunit_handle, *)
         END DO
         READ (extunit_handle, *) nat, rdum
         DO i = 1, 3
            READ (extunit_handle, *) ndum, rdum
            IF ((ndum /= npoints(i) .OR. (ABS(rdum(i) - dr(i)) > 1e-4)) .AND. &
                output_unit > 0) THEN
               WRITE (output_unit, *) "Restart from density | ERROR! | CUBE FILE NOT COINCIDENT WITH INTERNAL GRID ", i
               WRITE (output_unit, *) "Restart from density | ", ndum, " DIFFERS FROM ", npoints(i)
               WRITE (output_unit, *) "Restart from density | ", rdum, " DIFFERS FROM ", dr(i)
            END IF
         END DO
         !ignore atomic position data - read from coord or topology instead
         DO i = 1, nat
            READ (extunit_handle, *)
         END DO
         ! Get byte offset
         INQUIRE (extunit_handle, POS=offset_global)
         CALL close_file(unit_number=extunit_handle)
      END IF
      ! Sync offset and start parallel read
      CALL gid%bcast(offset_global, grid%pw_grid%para%group%source)
      BOF = offset_global
      CALL extunit%open(groupid=gid, filepath=filename, amode_status=file_amode_rdonly)
      ! Determine byte offsets for each grid z-slice which are local to a process
      ALLOCATE (displacements(nslices))
      displacements = 0
      DO i = lbounds(1), ubounds(1)
         should_read(:) = .TRUE.
         IF (i .LT. lbounds_local(1)) THEN
            should_read(1) = .FALSE.
         ELSE IF (i .GT. ubounds_local(1)) THEN
            EXIT
         END IF
         DO j = lbounds(2), ubounds(2)
            should_read(2) = .TRUE.
            IF (j .LT. lbounds_local(2) .OR. j .GT. ubounds_local(2)) THEN
               should_read(2) = .FALSE.
            END IF
            IF (ALL(should_read .EQV. .TRUE.)) THEN
               IF (islice > nslices) CPABORT("Index out of bounds.")
               displacements(islice) = BOF
               islice = islice + 1
            END IF
            ! Update global byte offset
            BOF = BOF + msglen
         END DO
      END DO
      ! Size of each z-slice is msglen
      ALLOCATE (blocklengths(nslices))
      blocklengths(:) = msglen
      ! Create indexed MPI type using calculated byte offsets as displacements and use it as a file view
      mp_file_desc = mp_file_type_hindexed_make_chv(nslices, blocklengths, displacements)
      BOF = 0
      CALL mp_file_type_set_view_chv(extunit, BOF, mp_file_desc)
      ! Collective read of cube
      ALLOCATE (readbuffer(nslices))
      readbuffer(:) = ''
      CALL extunit%read_all(msglen, nslices, readbuffer, mp_file_desc)
      CALL mp_file_type_free(mp_file_desc)
      CALL extunit%close()
      ! Convert cube values string -> real
      i = lbounds_local(1)
      j = lbounds_local(2)
      ALLOCATE (buffer(lbounds(3):ubounds(3)))
      buffer = 0.0_dp
      ! Test if the compiler supports parsing lines with line breaks in them
      IF (parse_test) THEN
         READ (readbuffer(1), *, IOSTAT=readstat) (buffer(k), k=lbounds(3), ubounds(3))
         IF (readstat == 0) THEN
            parses_linebreaks = .TRUE.
         ELSE
            parses_linebreaks = .FALSE.
         END IF
         parse_test = .FALSE.
         buffer = 0.0_dp
      END IF
      DO islice = 1, nslices
         IF (parses_linebreaks) THEN
            ! Use list-directed conversion if supported
            ! Produces faster code, but maybe the latter method could be optimized
            READ (readbuffer(islice), *) (buffer(k), k=lbounds(3), ubounds(3))
         ELSE
            ! Convert values by going through the z-slice one value at a time
            tmp = readbuffer(islice)
            pos = 1
            ientry = 1
            DO k = lbounds_local(3), ubounds_local(3)
               IF (MODULO(ientry, num_entries_line) == 0 .OR. k == ubounds_local(3)) THEN
                  ! Last value on line, dont read line break
                  READ (tmp(pos:pos + (entry_len - 2)), '(E12.5)') buffer(k)
                  pos = pos + (entry_len + 1)
               ELSE
                  READ (tmp(pos:pos + (entry_len - 1)), '(E13.5)') buffer(k)
                  pos = pos + entry_len
               END IF
               ientry = ientry + 1
            END DO
         END IF
         ! Optionally scale cube file values
         grid%array(i, j, lbounds(3):ubounds(3)) = scaling*buffer(lbounds(3):ubounds(3))
         j = j + 1
         IF (j > ubounds_local(2)) THEN
            j = lbounds_local(2)
            i = i + 1
         END IF
      END DO
      DEALLOCATE (readbuffer)
      DEALLOCATE (blocklengths, displacements)
      IF (debug_this_module) THEN
         ! Check that cube was correctly read using intrinsic read on master who sends data to everyone
         buffer = 0.0_dp
         IF (my_rank == 0) THEN
            IF (output_unit > 0 .AND. .NOT. be_silent) THEN
               WRITE (output_unit, FMT="(/,T2,A,/,/,T2,A)") "Reading the cube file:     ", filename
            END IF

            CALL open_file(file_name=filename, &
                           file_status="OLD", &
                           file_form="FORMATTED", &
                           file_action="READ", &
                           unit_number=extunit_handle)

            !skip header comments
            DO i = 1, 2
               READ (extunit_handle, *)
            END DO
            READ (extunit_handle, *) nat, rdum
            DO i = 1, 3
               READ (extunit_handle, *) ndum, rdum
               IF ((ndum /= npoints(i) .OR. (ABS(rdum(i) - dr(i)) > 1e-4)) .AND. &
                   output_unit > 0) THEN
                  WRITE (output_unit, *) "Restart from density | ERROR! | CUBE FILE NOT COINCIDENT WITH INTERNAL GRID ", i
                  WRITE (output_unit, *) "Restart from density | ", ndum, " DIFFERS FROM ", npoints(i)
                  WRITE (output_unit, *) "Restart from density | ", rdum, " DIFFERS FROM ", dr(i)
               END IF
            END DO
            !ignore atomic position data - read from coord or topology instead
            DO i = 1, nat
               READ (extunit_handle, *)
            END DO
         END IF

         !master sends all data to everyone
         DO i = lbounds(1), ubounds(1)
            DO j = lbounds(2), ubounds(2)
               IF (my_rank .EQ. 0) THEN
                  READ (extunit_handle, *) (buffer(k), k=lbounds(3), ubounds(3))
               END IF
               CALL gid%bcast(buffer(lbounds(3):ubounds(3)), 0)

               !only use data that is local to me - i.e. in slice of pencil I own
               IF ((lbounds_local(1) .LE. i) .AND. (i .LE. ubounds_local(1)) .AND. (lbounds_local(2) .LE. j) &
                   .AND. (j .LE. ubounds_local(2))) THEN
                  !allow scaling of external potential values by factor 'scaling' (SCALING_FACTOR in input file)
                  IF (ANY(grid%array(i, j, lbounds(3):ubounds(3)) /= buffer(lbounds(3):ubounds(3))*scaling)) &
                     CALL cp_abort(__LOCATION__, &
                                   "Error in parallel read of input cube file.")
               END IF

            END DO
         END DO

         IF (my_rank == 0) CALL close_file(unit_number=extunit_handle)

         CALL gid%sync()
      END IF
      DEALLOCATE (buffer)

   END SUBROUTINE cube_to_pw_parallel

! **************************************************************************************************
!> \brief Writes a realspace potential to a cube file using collective MPI I/O.
!> \param grid        the pw to output to the cube file
!> \param unit_nr     the handle associated with the cube file
!> \param title       title of the cube file
!> \param particles_r Cartersian coordinates of the system
!> \param particles_z atomic masses of atoms in the system
!> \param stride      every stride(i)th value of the potential is outputted (i=x,y,z)
!> \param zero_tails  flag that determines if small values of the potential should be zeroed
!> \param msglen      the size of each grid slice along z-axis in bytes
!> \par History
!>      Created [Nico Holmberg] (11.2017)
! **************************************************************************************************
   SUBROUTINE pw_to_cube_parallel(grid, unit_nr, title, particles_r, particles_z, stride, zero_tails, &
                                  msglen)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: grid
      TYPE(mp_file_type), INTENT(IN)                     :: unit_nr
      CHARACTER(*), INTENT(IN), OPTIONAL                 :: title
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: particles_r
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: particles_z
      INTEGER, INTENT(IN)                                :: stride(3)
      LOGICAL, INTENT(IN)                                :: zero_tails
      INTEGER, INTENT(IN)                                :: msglen

      INTEGER, PARAMETER                                 :: entry_len = 13, header_len = 41, &
                                                            header_len_z = 53, num_entries_line = 6

      CHARACTER(LEN=entry_len)                           :: value
      CHARACTER(LEN=header_len)                          :: header
      CHARACTER(LEN=header_len_z)                        :: header_z
      INTEGER, DIMENSION(3)                              :: lbounds, lbounds_local, ubounds, &
                                                            ubounds_local
      INTEGER, ALLOCATABLE, DIMENSION(:), TARGET         :: blocklengths
      INTEGER(kind=file_offset), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: displacements
      INTEGER(kind=file_offset)                          :: BOF
      INTEGER                                            :: counter, i, islice, j, k, last_z, &
                                                            my_rank, np, nslices, size_of_z
      CHARACTER(LEN=msglen), ALLOCATABLE, DIMENSION(:)   :: writebuffer
      CHARACTER(LEN=msglen)                              :: tmp
      LOGICAL                                            :: should_write(2)
      TYPE(mp_comm_type)                                 :: gid
      TYPE(mp_file_descriptor_type)                      :: mp_desc

      !get rs grids and parallel envnment
      gid = grid%pw_grid%para%group
      my_rank = grid%pw_grid%para%group%mepos

      ! Shortcut
      lbounds = grid%pw_grid%bounds(1, :)
      ubounds = grid%pw_grid%bounds(2, :)
      lbounds_local = grid%pw_grid%bounds_local(1, :)
      ubounds_local = grid%pw_grid%bounds_local(2, :)
      ! Determine the total number of z-slices and the number of values per slice
      size_of_z = CEILING(REAL(ubounds_local(3) - lbounds_local(3) + 1, dp)/REAL(stride(3), dp))
      islice = 1
      DO i = lbounds(1), ubounds(1), stride(1)
         should_write(:) = .TRUE.
         IF (i .LT. lbounds_local(1)) THEN
            should_write(1) = .FALSE.
         ELSE IF (i .GT. ubounds_local(1)) THEN
            EXIT
         END IF
         DO j = lbounds(2), ubounds(2), stride(2)
            should_write(2) = .TRUE.
            IF (j .LT. lbounds_local(2) .OR. j .GT. ubounds_local(2)) THEN
               should_write(2) = .FALSE.
            END IF
            IF (ALL(should_write .EQV. .TRUE.)) THEN
               islice = islice + 1
            END IF
         END DO
      END DO
      nslices = islice - 1
      DO k = lbounds(3), ubounds(3), stride(3)
         IF (k + stride(3) > ubounds(3)) last_z = k
      END DO
      islice = 1
      ! Determine initial byte offset (0 or EOF if data is appended)
      CALL unit_nr%get_position(BOF)
      ! Write header information on master process and update byte offset accordingly
      IF (my_rank == 0) THEN
         ! this format seems to work for e.g. molekel and gOpenmol
         ! latest version of VMD can read non orthorhombic cells
         CALL unit_nr%write_at(BOF, "-Quickstep-"//NEW_LINE("C"))
         BOF = BOF + LEN("-Quickstep-"//NEW_LINE("C"))*mpi_character_size
         IF (PRESENT(title)) THEN
            CALL unit_nr%write_at(BOF, TRIM(title)//NEW_LINE("C"))
            BOF = BOF + LEN(TRIM(title)//NEW_LINE("C"))*mpi_character_size
         ELSE
            CALL unit_nr%write_at(BOF, "No Title"//NEW_LINE("C"))
            BOF = BOF + LEN("No Title"//NEW_LINE("C"))*mpi_character_size
         END IF

         CPASSERT(PRESENT(particles_z) .EQV. PRESENT(particles_r))
         np = 0
         IF (PRESENT(particles_z)) THEN
            CPASSERT(SIZE(particles_z) == SIZE(particles_r, dim=2))
            ! cube files can only be written for 99999 particles due to a format limitation (I5)
            ! so we limit the number of particles written.
            np = MIN(99999, SIZE(particles_z))
         END IF

         WRITE (header, '(I5,3f12.6)') np, 0.0_dp, 0._dp, 0._dp !start of cube
         CALL unit_nr%write_at(BOF, header//NEW_LINE("C"))
         BOF = BOF + LEN(header//NEW_LINE("C"))*mpi_character_size

         WRITE (header, '(I5,3f12.6)') (grid%pw_grid%npts(1) + stride(1) - 1)/stride(1), &
            grid%pw_grid%dh(1, 1)*REAL(stride(1), dp), grid%pw_grid%dh(2, 1)*REAL(stride(1), dp), &
            grid%pw_grid%dh(3, 1)*REAL(stride(1), dp)
         CALL unit_nr%write_at(BOF, header//NEW_LINE("C"))
         BOF = BOF + LEN(header//NEW_LINE("C"))*mpi_character_size

         WRITE (header, '(I5,3f12.6)') (grid%pw_grid%npts(2) + stride(2) - 1)/stride(2), &
            grid%pw_grid%dh(1, 2)*REAL(stride(2), dp), grid%pw_grid%dh(2, 2)*REAL(stride(2), dp), &
            grid%pw_grid%dh(3, 2)*REAL(stride(2), dp)
         CALL unit_nr%write_at(BOF, header//NEW_LINE("C"))
         BOF = BOF + LEN(header//NEW_LINE("C"))*mpi_character_size

         WRITE (header, '(I5,3f12.6)') (grid%pw_grid%npts(3) + stride(3) - 1)/stride(3), &
            grid%pw_grid%dh(1, 3)*REAL(stride(3), dp), grid%pw_grid%dh(2, 3)*REAL(stride(3), dp), &
            grid%pw_grid%dh(3, 3)*REAL(stride(3), dp)
         CALL unit_nr%write_at(BOF, header//NEW_LINE("C"))
         BOF = BOF + LEN(header//NEW_LINE("C"))*mpi_character_size

         IF (PRESENT(particles_z)) THEN
            DO i = 1, np
               WRITE (header_z, '(I5,4f12.6)') particles_z(i), 0._dp, particles_r(:, i)
               CALL unit_nr%write_at(BOF, header_z//NEW_LINE("C"))
               BOF = BOF + LEN(header_z//NEW_LINE("C"))*mpi_character_size
            END DO
         END IF
      END IF
      ! Sync offset
      CALL gid%bcast(BOF, grid%pw_grid%para%group%source)
      ! Determine byte offsets for each grid z-slice which are local to a process
      ! and convert z-slices to cube format compatible strings
      ALLOCATE (displacements(nslices))
      displacements = 0
      ALLOCATE (writebuffer(nslices))
      writebuffer(:) = ''
      DO i = lbounds(1), ubounds(1), stride(1)
         should_write(:) = .TRUE.
         IF (i .LT. lbounds_local(1)) THEN
            should_write(1) = .FALSE.
         ELSE IF (i .GT. ubounds_local(1)) THEN
            EXIT
         END IF
         DO j = lbounds(2), ubounds(2), stride(2)
            should_write(2) = .TRUE.
            IF (j .LT. lbounds_local(2) .OR. j .GT. ubounds_local(2)) THEN
               should_write(2) = .FALSE.
            END IF
            IF (ALL(should_write .EQV. .TRUE.)) THEN
               IF (islice > nslices) CPABORT("Index out of bounds.")
               displacements(islice) = BOF
               tmp = ''
               counter = 0
               DO k = lbounds(3), ubounds(3), stride(3)
                  IF (zero_tails .AND. grid%array(i, j, k) < 1.E-7_dp) THEN
                     WRITE (value, '(E13.5)') 0.0_dp
                  ELSE
                     WRITE (value, '(E13.5)') grid%array(i, j, k)
                  END IF
                  tmp = TRIM(tmp)//TRIM(value)
                  counter = counter + 1
                  IF (MODULO(counter, num_entries_line) == 0 .OR. k == last_z) &
                     tmp = TRIM(tmp)//NEW_LINE('C')
               END DO
               writebuffer(islice) = tmp
               islice = islice + 1
            END IF
            ! Update global byte offset
            BOF = BOF + msglen
         END DO
      END DO
      ! Create indexed MPI type using calculated byte offsets as displacements
      ! Size of each z-slice is msglen
      ALLOCATE (blocklengths(nslices))
      blocklengths(:) = msglen
      mp_desc = mp_file_type_hindexed_make_chv(nslices, blocklengths, displacements)
      ! Use the created type as a file view
      ! NB. The vector 'displacements' contains the absolute offsets of each z-slice i.e.
      ! they are given relative to the beginning of the file. The global offset to
      ! set_view must therefore be set to 0
      BOF = 0
      CALL mp_file_type_set_view_chv(unit_nr, BOF, mp_desc)
      ! Collective write of cube
      CALL unit_nr%write_all(msglen, nslices, writebuffer, mp_desc)
      ! Clean up
      CALL mp_file_type_free(mp_desc)
      DEALLOCATE (writebuffer)
      DEALLOCATE (blocklengths, displacements)

   END SUBROUTINE pw_to_cube_parallel

! **************************************************************************************************
!> \brief Prints a simple grid file: X Y Z value
!> \param pw ...
!> \param unit_nr ...
!> \param stride ...
!> \param pw2 ...
!> \par History
!>      Created [Vladimir Rybkin] (08.2018)
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE pw_to_simple_volumetric(pw, unit_nr, stride, pw2)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: stride
      TYPE(pw_r3d_rs_type), INTENT(IN), OPTIONAL         :: pw2

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_to_simple_volumetric'

      INTEGER                                            :: checksum, dest, handle, i, I1, I2, I3, &
                                                            ip, L1, L2, L3, my_rank, my_stride(3), &
                                                            ngrids, npoints, num_pe, rank(2), &
                                                            source, tag, U1, U2, U3
      LOGICAL                                            :: DOUBLE
      REAL(KIND=dp)                                      :: x, y, z
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: buf, buf2
      TYPE(mp_comm_type)                                 :: gid

      CALL timeset(routineN, handle)

      ! Check if we write two grids
      DOUBLE = .FALSE.
      IF (PRESENT(pw2)) DOUBLE = .TRUE.

      my_stride = 1
      IF (PRESENT(stride)) THEN
         IF (SIZE(stride) /= 1 .AND. SIZE(stride) /= 3) &
            CALL cp_abort(__LOCATION__, "STRIDE keyword can accept only 1 "// &
                          "(the same for X,Y,Z) or 3 values. Correct your input file.")
         IF (SIZE(stride) == 1) THEN
            DO i = 1, 3
               my_stride(i) = stride(1)
            END DO
         ELSE
            my_stride = stride(1:3)
         END IF
         CPASSERT(my_stride(1) > 0)
         CPASSERT(my_stride(2) > 0)
         CPASSERT(my_stride(3) > 0)
      END IF

      ! shortcut
      L1 = pw%pw_grid%bounds(1, 1)
      L2 = pw%pw_grid%bounds(1, 2)
      L3 = pw%pw_grid%bounds(1, 3)
      U1 = pw%pw_grid%bounds(2, 1)
      U2 = pw%pw_grid%bounds(2, 2)
      U3 = pw%pw_grid%bounds(2, 3)

      ! Write the header: number of points and number of spins
      ngrids = 1
      IF (DOUBLE) ngrids = 2
      npoints = ((pw%pw_grid%npts(1) + my_stride(1) - 1)/my_stride(1))* &
                ((pw%pw_grid%npts(2) + my_stride(2) - 1)/my_stride(1))* &
                ((pw%pw_grid%npts(3) + my_stride(3) - 1)/my_stride(1))
      IF (unit_nr > 1) WRITE (unit_nr, '(I7,I5)') npoints, ngrids

      ALLOCATE (buf(L3:U3))
      IF (DOUBLE) ALLOCATE (buf2(L3:U3))

      my_rank = pw%pw_grid%para%group%mepos
      gid = pw%pw_grid%para%group
      num_pe = pw%pw_grid%para%group%num_pe
      tag = 1

      rank (1) = unit_nr
      rank (2) = my_rank
      checksum = 0
      IF (unit_nr > 0) checksum = 1

      CALL gid%sum(checksum)
      CPASSERT(checksum == 1)

      CALL gid%maxloc(rank)
      CPASSERT(rank(1) > 0)

      dest = rank(2)
      DO I1 = L1, U1, my_stride(1)
         DO I2 = L2, U2, my_stride(2)

            ! cycling through the CPUs, check if the current ray (I1,I2) is local to that CPU
            IF (pw%pw_grid%para%mode .NE. PW_MODE_LOCAL) THEN
               DO ip = 0, num_pe - 1
                  IF (pw%pw_grid%para%bo(1, 1, ip, 1) <= I1 - L1 + 1 .AND. pw%pw_grid%para%bo(2, 1, ip, 1) >= I1 - L1 + 1 .AND. &
                      pw%pw_grid%para%bo(1, 2, ip, 1) <= I2 - L2 + 1 .AND. pw%pw_grid%para%bo(2, 2, ip, 1) >= I2 - L2 + 1) THEN
                     source = ip
                  END IF
               END DO
            ELSE
               source = dest
            END IF

            IF (source == dest) THEN
               IF (my_rank == source) THEN
                  buf(:) = pw%array(I1, I2, :)
                  IF (DOUBLE) buf2(:) = pw2%array(I1, I2, :)
               END IF
            ELSE
               IF (my_rank == source) THEN
                  buf(:) = pw%array(I1, I2, :)
                  CALL gid%send(buf, dest, tag)
                  IF (DOUBLE) THEN
                     buf2(:) = pw2%array(I1, I2, :)
                     CALL gid%send(buf2, dest, tag)
                  END IF
               END IF
               IF (my_rank == dest) THEN
                  CALL gid%recv(buf, source, tag)
                  IF (DOUBLE) CALL gid%recv(buf2, source, tag)
               END IF
            END IF

            IF (.NOT. DOUBLE) THEN
               DO I3 = L3, U3, my_stride(3)
                  x = pw%pw_grid%dh(1, 1)*I1 + &
                      pw%pw_grid%dh(2, 1)*I2 + &
                      pw%pw_grid%dh(3, 1)*I3

                  y = pw%pw_grid%dh(1, 2)*I1 + &
                      pw%pw_grid%dh(2, 2)*I2 + &
                      pw%pw_grid%dh(3, 2)*I3

                  z = pw%pw_grid%dh(1, 3)*I1 + &
                      pw%pw_grid%dh(2, 3)*I2 + &
                      pw%pw_grid%dh(3, 3)*I3

                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, '(6E13.5, 6E13.5, 6E13.5, 6E13.5)') x, y, z, buf(I3)
                  END IF
               END DO

            ELSE

               DO I3 = L3, U3, my_stride(3)
                  x = pw%pw_grid%dh(1, 1)*I1 + &
                      pw%pw_grid%dh(2, 1)*I2 + &
                      pw%pw_grid%dh(3, 1)*I3

                  y = pw%pw_grid%dh(1, 2)*I1 + &
                      pw%pw_grid%dh(2, 2)*I2 + &
                      pw%pw_grid%dh(3, 2)*I3

                  z = pw%pw_grid%dh(1, 3)*I1 + &
                      pw%pw_grid%dh(2, 3)*I2 + &
                      pw%pw_grid%dh(3, 3)*I3

                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, '(6E13.5, 6E13.5, 6E13.5, 6E13.5)') x, y, z, buf(I3), buf2(I3)
                  END IF
               END DO

            END IF ! Double

            ! this double loop generates so many messages that it can overload
            ! the message passing system, e.g. on XT3
            ! we therefore put a barrier here that limits the amount of message
            ! that flies around at any given time.
            ! if ever this routine becomes a bottleneck, we should go for a
            ! more complicated rewrite
            CALL gid%sync()

         END DO
      END DO

      DEALLOCATE (buf)
      IF (DOUBLE) DEALLOCATE (buf2)

      CALL timestop(handle)

   END SUBROUTINE pw_to_simple_volumetric

END MODULE realspace_grid_openpmd
