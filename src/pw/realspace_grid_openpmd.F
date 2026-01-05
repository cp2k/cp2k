!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Generate Gaussian cube files
! **************************************************************************************************
MODULE realspace_grid_openpmd

#ifdef __OPENPMD
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
#if defined(__parallel)
#if defined(__MPI_F08)
   USE mpi_f08, ONLY: mpi_allreduce, mpi_integer, mpi_bxor, mpi_allgather
#else
   USE mpi, ONLY: mpi_allreduce, mpi_integer, mpi_bxor, mpi_allgather
#endif
#endif

   USE openpmd_api, ONLY: &
      openpmd_attributable_type, &
      openpmd_dynamic_memory_view_type_1d, &
      openpmd_dynamic_memory_view_type_3d, &
      openpmd_mesh_type, &
      openpmd_particle_species_type, &
      openpmd_record_component_type, &
      openpmd_record_type, openpmd_type_double, openpmd_type_int
   USE pw_grid_types, ONLY: PW_MODE_LOCAL
   USE pw_types, ONLY: pw_r3d_rs_type
   USE util, ONLY: sort_unique

#else

   USE pw_types, ONLY: pw_r3d_rs_type
   USE kinds, ONLY: dp

#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC ::pw_to_openpmd

#ifdef __OPENPMD
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'realspace_grid_openpmd'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   LOGICAL, PRIVATE                     :: parses_linebreaks = .FALSE., &
                                           parse_test = .TRUE.

   TYPE cp_openpmd_write_buffer_1d
      REAL(KIND=dp), POINTER :: buffer(:)
   END TYPE cp_openpmd_write_buffer_1d
#endif

CONTAINS

#ifdef __OPENPMD
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
!> \param particles_zeff ...
!> \param atom_type ...
!> \param atom_count ...
!> \param openpmd_data ...
!> \param do_write_data ...
! **************************************************************************************************
   SUBROUTINE pw_write_particle_species( &
      particles_z, &
      particles_r, &
      particles_zeff, &
      atom_type, &
      atom_count, &
      openpmd_data, &
      do_write_data &
      )
      INTEGER, DIMENSION(:), INTENT(IN)                  :: particles_z
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: particles_r
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: particles_zeff
      INTEGER, INTENT(IN)                                :: atom_type, atom_count
      TYPE(cp_openpmd_per_call_value_type)               :: openpmd_data
      LOGICAL                                            :: do_write_data

      CHARACTER(len=1), DIMENSION(3), PARAMETER          :: dims = ["x", "y", "z"]

      CHARACTER(len=3)                                   :: atom_type_as_string
      CHARACTER(len=default_string_length)               :: species_name
      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(1)                              :: global_extent, global_offset, &
                                                            local_extent
      TYPE(cp_openpmd_write_buffer_1d)                   :: charge_write_buffer
      TYPE(cp_openpmd_write_buffer_1d), DIMENSION(3)     :: write_buffers
      TYPE(openpmd_attributable_type)                    :: attr
      TYPE(openpmd_dynamic_memory_view_type_1d)          :: unresolved_charge_write_buffer
      TYPE(openpmd_dynamic_memory_view_type_1d), &
         DIMENSION(3)                                    :: unresolved_write_buffers
      TYPE(openpmd_particle_species_type)                :: species
      TYPE(openpmd_record_component_type)                :: charge_component, position_component, &
                                                            position_offset_component
      TYPE(openpmd_record_type)                          :: charge, position, position_offset

! TODO: The charge is probably constant per species?
! If yes, we could use a constant component and save storage space

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

      species = openpmd_data%iteration%get_particle_species(TRIM(species_name))

      position_offset = species%get_record("positionOffset")
      position = species%get_record("position")
      DO k = 1, SIZE(dims)
         position_offset_component = position_offset%get_component(dims(k))
         CALL position_offset_component%make_constant_zero(openpmd_type_int, global_extent)
         position_component = position%get_component(dims(k))
         CALL position_component%reset_dataset(openpmd_type_double, global_extent)
         IF (do_write_data) THEN
            unresolved_write_buffers(k) = &
               position_component%store_chunk_span_1d_double(global_offset, local_extent)
            write_buffers(k)%buffer => unresolved_write_buffers(k)%resolve_double(DEALLOCATE=.FALSE.)
         END IF
      END DO

      IF (PRESENT(particles_zeff)) THEN
         charge = species%get_record("charge")
         charge_component = charge%as_record_component()
         CALL charge_component%reset_dataset(openpmd_type_double, global_extent)
         IF (do_write_data) THEN
            unresolved_charge_write_buffer = charge_component%store_chunk_span_1d_double(global_offset, local_extent)
            charge_write_buffer%buffer => unresolved_charge_write_buffer%resolve_double(DEALLOCATE=.FALSE.)
         END IF
      END IF

      IF (do_write_data) THEN
         ! Resolve Spans for a second time to allow for internal reallocations in BP4 engine of ADIOS2
         DO k = 1, SIZE(dims)
            write_buffers(k)%buffer = unresolved_write_buffers(k)%resolve_double(DEALLOCATE=.TRUE.)
         END DO
         IF (PRESENT(particles_zeff)) THEN
            charge_write_buffer%buffer = unresolved_charge_write_buffer%resolve_double(DEALLOCATE=.TRUE.)
         END IF
         j = 1
         DO i = 1, SIZE(particles_z)
            IF (particles_z(i) == atom_type) THEN
               DO k = 1, 3
                  write_buffers(k)%buffer(j) = particles_r(k, i)
               END DO
               IF (PRESENT(particles_zeff)) THEN
                  charge_write_buffer%buffer(j) = particles_zeff(i)
               END IF
               j = j + 1
            END IF
         END DO
      END IF
      attr = openpmd_data%iteration%as_attributable()
      CALL attr%series_flush("hdf5.independent_stores = true")
   END SUBROUTINE pw_write_particle_species

! **************************************************************************************************
!> \brief ...
!> \param particles_z ...
!> \param particles_r ...
!> \param particles_zeff ...
!> \param atom_types ...
!> \param atom_counts ...
!> \param num_atom_types ...
!> \param openpmd_data ...
!> \param gid ...
! **************************************************************************************************
   SUBROUTINE pw_write_particles( &
      particles_z, &
      particles_r, &
      particles_zeff, &
      atom_types, &
      atom_counts, &
      num_atom_types, &
      openpmd_data, &
      gid &
      )
      INTEGER, DIMENSION(:), INTENT(IN)                  :: particles_z
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: particles_r
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: particles_zeff
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
         CALL pw_write_particle_species( &
            particles_z, &
            particles_r, &
            particles_zeff, &
            atom_types(i), &
            atom_counts(i), &
            openpmd_data, &
            do_write_data &
            )
      END DO
   END SUBROUTINE pw_write_particles

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param unit_nr ...
!> \param title ...
!> \param particles_r ...
!> \param particles_z ...
!> \param particles_zeff ...
!> \param stride ...
!> \param zero_tails ...
!> \param silent ...
!> \param mpi_io ...
! **************************************************************************************************
   SUBROUTINE pw_to_openpmd( &
      pw, &
      unit_nr, &
      title, &
      particles_r, &
      particles_z, &
      particles_zeff, &
      stride, &
      zero_tails, &
      silent, &
      mpi_io &
      )
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw
      INTEGER                                            :: unit_nr
      CHARACTER(*), INTENT(IN), OPTIONAL                 :: title
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: particles_r
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: particles_z
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: particles_zeff
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: stride
      LOGICAL, INTENT(IN), OPTIONAL                      :: zero_tails, silent, mpi_io

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_to_openpmd'
      INTEGER, PARAMETER                                 :: entry_len = 13, num_entries_line = 6

      INTEGER                                            :: count1, count2, count3, handle, i, I1, &
                                                            I2, I3, iat, L1, L2, L3, my_rank, &
                                                            my_stride(3), np, num_atom_types, &
                                                            num_pe, U1, U2, U3
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_counts, atom_types
      INTEGER, DIMENSION(3)                              :: global_extent, local_extent, offset
      LOGICAL                                            :: be_silent, my_zero_tails, parallel_write
      REAL(KIND=dp), DIMENSION(3)                        :: grid_spacing
      REAL(KIND=dp), POINTER                             :: write_buffer(:, :, :)
      TYPE(cp_openpmd_per_call_value_type)               :: openpmd_data
      TYPE(mp_comm_type)                                 :: gid
      TYPE(openpmd_attributable_type)                    :: attr
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

      CPASSERT(PRESENT(particles_z) .EQV. PRESENT(particles_r))
      np = 0
      IF (PRESENT(particles_z)) THEN
         CALL pw_get_atom_types(particles_z, atom_types, atom_counts, num_atom_types)
         CPASSERT(SIZE(particles_z) == SIZE(particles_r, dim=2))
         np = SIZE(particles_z)
      END IF

      DO i = 1, 3
         ! Notes:
         ! 1. This loses information on the rotation of the mesh, the mesh is stored
         !    without reference to a global coordinate system
         ! 2. This assumes that the coordinate system is not sheared
         grid_spacing(i) = SQRT(SUM(pw%pw_grid%dh(:, i)**2))*REAL(my_stride(i), dp)
      END DO

      IF (PRESENT(particles_z)) THEN
         IF (parallel_write) THEN
            CALL pw_write_particles( &
               particles_z, &
               particles_r, &
               particles_zeff, &
               atom_types, &
               atom_counts, &
               num_atom_types, &
               openpmd_data, &
               gid &
               )
         ELSE
            CALL pw_write_particles( &
               particles_z, &
               particles_r, &
               particles_zeff, &
               atom_types, &
               atom_counts, &
               num_atom_types, &
               openpmd_data &
               )
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

      mesh = openpmd_data%iteration%get_mesh(TRIM(openpmd_data%name_prefix))
      CALL mesh%set_axis_labels(["x", "y", "z"])
      CALL mesh%set_position([0.5_dp, 0.5_dp, 0.5_dp])
      CALL mesh%set_grid_global_offset([0._dp, 0._dp, 0._dp])
      CALL mesh%set_grid_spacing(grid_spacing)
      scalar_mesh = mesh%as_record_component()
      CALL scalar_mesh%reset_dataset(openpmd_type_double, global_extent)

      ! shortcut
      ! need to adjust L1/U1 for uneven distributions across MPI ranks
      ! (when working with a stride, we might have to skip the first n values)
      ! so keep this consistent with the offset and local_extent computed above
      ! L1 = pw%pw_grid%bounds_local(1, 1)
      L1 = offset(1)*my_stride(1)
      L2 = pw%pw_grid%bounds_local(1, 2)
      L3 = pw%pw_grid%bounds_local(1, 3)
      ! U1 = pw%pw_grid%bounds_local(2, 1)
      ! offset + local_extent is the start index for the next rank already
      ! since the indexes are inclusive, subtract 1 from the boundary index
      U1 = (offset(1) + local_extent(1) - 1)*my_stride(1)
      U2 = pw%pw_grid%bounds_local(2, 2)
      U3 = pw%pw_grid%bounds_local(2, 3)

      my_rank = pw%pw_grid%para%group%mepos
      num_pe = pw%pw_grid%para%group%num_pe

      IF (ALL(my_stride == 1)) THEN
         CALL scalar_mesh%store_chunk(pw%array(L1:U1, L2:U2, L3:U3), offset)
         ! Are there some conditions under which we can skip this flush?
         attr = openpmd_data%iteration%as_attributable()
         CALL attr%series_flush("hdf5.independent_stores = false")
      ELSE
         count3 = 0
         DO I3 = L3, U3, my_stride(3)
            ! maybe add an overload to provide `buf` here for HDF5, might have better performance
            ! for intermittent flushing
            ! or just call the buffer in the outer function if memory is no problem...
            unresolved_write_buffer = scalar_mesh%store_chunk_span_3d_double( &
                                      [offset(1), offset(2), offset(3) + count3], &
                                      [local_extent(1), local_extent(2), 1])
            write_buffer => unresolved_write_buffer%resolve_double(DEALLOCATE=.TRUE.)

            ! Sanity checks: ensure buffer is associated and matches expected shape
            CPASSERT(ASSOCIATED(write_buffer))
            CPASSERT(SIZE(write_buffer, 1) == local_extent(1))
            CPASSERT(SIZE(write_buffer, 2) == local_extent(2))
            CPASSERT(SIZE(write_buffer, 3) == 1)

            count2 = 0
            DO I2 = L2, U2, my_stride(2)
               ! This loop deals with ray (:, count2, count3) of the local subspace
               ! The write buffer itself has been allocated for slice (:, :, count3)
               count1 = 0
               DO I1 = L1, U1, my_stride(1)
                  write_buffer(count1 + 1, count2 + 1, 1) = pw%array(I1, I2, I3)
                  ! Debug: print the target indices in write_buffer to the command line
                  ! WRITE(*,*) 'write_buffer index:', count1 + 1, ',', count2 + 1, ',', 1
                  count1 = count1 + 1
               END DO
               count2 = count2 + 1
            END DO
            count3 = count3 + 1
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_to_openpmd

#else

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param unit_nr ...
!> \param title ...
!> \param particles_r ...
!> \param particles_z ...
!> \param particles_zeff ...
!> \param stride ...
!> \param zero_tails ...
!> \param silent ...
!> \param mpi_io ...
! **************************************************************************************************
   SUBROUTINE pw_to_openpmd( &
      pw, &
      unit_nr, &
      title, &
      particles_r, &
      particles_z, &
      particles_zeff, &
      stride, &
      zero_tails, &
      silent, &
      mpi_io &
      )
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw
      INTEGER                                            :: unit_nr
      CHARACTER(*), INTENT(IN), OPTIONAL                 :: title
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: particles_r
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: particles_z
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: particles_zeff
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: stride
      LOGICAL, INTENT(IN), OPTIONAL                      :: zero_tails, silent, mpi_io

      MARK_USED(pw)
      MARK_USED(unit_nr)
      MARK_USED(title)
      MARK_USED(particles_r)
      MARK_USED(particles_z)
      MARK_USED(particles_zeff)
      MARK_USED(stride)
      MARK_USED(zero_tails)
      MARK_USED(silent)
      MARK_USED(mpi_io)
      CPABORT("CP2K compiled without the openPMD-api")

   END SUBROUTINE pw_to_openpmd

#endif

END MODULE realspace_grid_openpmd
