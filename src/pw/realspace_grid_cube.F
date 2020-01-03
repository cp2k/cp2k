!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Generate Gaussian cube files
! **************************************************************************************************
MODULE realspace_grid_cube
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: &
        file_amode_rdonly, file_offset, mp_bcast, mp_file_close, mp_file_descriptor_type, &
        mp_file_get_position, mp_file_open, mp_file_read_all_chv, mp_file_type_free, &
        mp_file_type_hindexed_make_chv, mp_file_type_set_view_chv, mp_file_write_all_chv, &
        mp_file_write_at, mp_maxloc, mp_recv, mp_send, mp_sum, mp_sync, mpi_character_size
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
   USE pw_types,                        ONLY: pw_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pw_to_cube, cube_to_pw, pw_to_simple_volumetric

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'realspace_grid_cube'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   LOGICAL, PRIVATE                     :: parses_linebreaks = .FALSE., &
                                           parse_test = .TRUE.

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param unit_nr ...
!> \param title ...
!> \param particles_r ...
!> \param particles_z ...
!> \param stride ...
!> \param zero_tails ...
!> \param mpi_io ...
! **************************************************************************************************
   SUBROUTINE pw_to_cube(pw, unit_nr, title, particles_r, particles_z, stride, zero_tails, mpi_io)
      TYPE(pw_type), POINTER                             :: pw
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(*), INTENT(IN), OPTIONAL                 :: title
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: particles_r
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: particles_z
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: stride
      LOGICAL, INTENT(IN), OPTIONAL                      :: zero_tails, mpi_io

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_to_cube', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: entry_len = 13, num_entries_line = 6

      INTEGER :: checksum, dest, gid, handle, i, I1, I2, I3, iat, ip, L1, L2, L3, msglen, my_rank, &
         my_stride(3), np, num_linebreak, num_pe, rank(2), size_of_z, source, tag, U1, U2, U3
      LOGICAL                                            :: my_zero_tails, parallel_write
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: buf

      CALL timeset(routineN, handle)

      my_zero_tails = .FALSE.
      parallel_write = .FALSE.
      IF (PRESENT(zero_tails)) my_zero_tails = zero_tails
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
            ENDIF

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
         ENDIF

         ! shortcut
         L1 = pw%pw_grid%bounds(1, 1)
         L2 = pw%pw_grid%bounds(1, 2)
         L3 = pw%pw_grid%bounds(1, 3)
         U1 = pw%pw_grid%bounds(2, 1)
         U2 = pw%pw_grid%bounds(2, 2)
         U3 = pw%pw_grid%bounds(2, 3)

         ALLOCATE (buf(L3:U3))

         my_rank = pw%pw_grid%para%my_pos
         gid = pw%pw_grid%para%group
         num_pe = pw%pw_grid%para%group_size
         tag = 1

         rank(1) = unit_nr
         rank(2) = my_rank
         checksum = 0
         IF (unit_nr > 0) checksum = 1

         CALL mp_sum(checksum, gid)
         CPASSERT(checksum == 1)

         CALL mp_maxloc(rank, gid)
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
                     ENDIF
                  ENDDO
               ELSE
                  source = dest
               ENDIF

               IF (source == dest) THEN
                  IF (my_rank == source) THEN
                     buf(:) = pw%cr3d(I1, I2, :)
                  ENDIF
               ELSE
                  IF (my_rank == source) THEN
                     buf(:) = pw%cr3d(I1, I2, :)
                     CALL mp_send(buf, dest, tag, gid)
                  ENDIF
                  IF (my_rank == dest) THEN
                     CALL mp_recv(buf, source, tag, gid)
                  ENDIF
               ENDIF

               IF (unit_nr > 0) THEN
                  IF (my_zero_tails) THEN
                     DO I3 = L3, U3
                        IF (buf(I3) < 1.E-7_dp) buf(I3) = 0.0_dp
                     END DO
                  END IF
                  WRITE (unit_nr, '(6E13.5)') (buf(I3), I3=L3, U3, my_stride(3))
               ENDIF

               ! this double loop generates so many messages that it can overload
               ! the message passing system, e.g. on XT3
               ! we therefore put a barrier here that limits the amount of message
               ! that flies around at any given time.
               ! if ever this routine becomes a bottleneck, we should go for a
               ! more complicated rewrite
               CALL mp_sync(gid)

            ENDDO
         ENDDO

         DEALLOCATE (buf)
      ELSE
         size_of_z = CEILING(REAL(pw%pw_grid%bounds(2, 3) - pw%pw_grid%bounds(1, 3) + 1, dp)/REAL(my_stride(3), dp))
         num_linebreak = size_of_z/num_entries_line
         IF (MODULO(size_of_z, num_entries_line) /= 0) &
            num_linebreak = num_linebreak + 1
         msglen = (size_of_z*entry_len + num_linebreak)*mpi_character_size
         CALL pw_to_cube_parallel(pw, unit_nr, title, particles_r, particles_z, my_stride, my_zero_tails, msglen)
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
!> \par History
!>      Created [M.Watkins] (01.2014)
!>      Use blocking, collective MPI read for parallel simulations [Nico Holmberg] (05.2017)
! **************************************************************************************************
   SUBROUTINE cube_to_pw(grid, filename, scaling, parallel_read)

      TYPE(pw_type), POINTER                             :: grid
      CHARACTER(len=*), INTENT(in)                       :: filename
      REAL(kind=dp), INTENT(in)                          :: scaling
      LOGICAL, INTENT(in)                                :: parallel_read

      CHARACTER(len=*), PARAMETER :: routineN = 'cube_to_pw', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: entry_len = 13, num_entries_line = 6

      INTEGER                                            :: extunit, gid, handle, i, ip, j, k, &
                                                            master, msglen, my_rank, nat, ndum, &
                                                            num_linebreak, num_pe, output_unit, &
                                                            size_of_z, tag
      INTEGER, DIMENSION(3)                              :: lbounds, lbounds_local, npoints, &
                                                            npoints_local, ubounds, ubounds_local
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: buffer
      REAL(kind=dp), DIMENSION(3)                        :: dr, rdum

      output_unit = cp_logger_get_default_io_unit()

      CALL timeset(routineN, handle)

      !get rs grids and parallel envnment
      gid = grid%pw_grid%para%group
      my_rank = grid%pw_grid%para%my_pos
      num_pe = grid%pw_grid%para%group_size
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
            IF (output_unit > 0) THEN
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
            ENDDO
            READ (extunit, *) nat, rdum
            DO i = 1, 3
               READ (extunit, *) ndum, rdum
               IF ((ndum /= npoints(i) .OR. (ABS(rdum(i) - dr(i)) > 1e-4)) .AND. &
                   output_unit > 0) THEN
                  WRITE (output_unit, *) "Restart from density | ERROR! | CUBE FILE NOT COINCIDENT WITH INTERNAL GRID ", i
                  WRITE (output_unit, *) "Restart from density | ", ndum, " DIFFERS FROM ", npoints(i)
                  WRITE (output_unit, *) "Restart from density | ", rdum, " DIFFERS FROM ", dr(i)
               ENDIF
            ENDDO
            !ignore atomic postion data - read from coord or topology instead
            DO i = 1, nat
               READ (extunit, *)
            ENDDO
         ENDIF

         !master sends all data to everyone
         DO i = lbounds(1), ubounds(1)
            DO j = lbounds(2), ubounds(2)
               IF (my_rank .EQ. 0) THEN
                  READ (extunit, *) (buffer(k), k=lbounds(3), ubounds(3))
                  IF (num_pe .GT. 1) THEN
                     DO ip = 1, num_pe - 1
                        CALL mp_send(buffer(lbounds(3):ubounds(3)), ip, tag, gid)
                     ENDDO
                  ENDIF
               ELSE
                  master = 0
                  CALL mp_recv(buffer(lbounds(3):ubounds(3)), master, tag, gid)
               ENDIF
               CALL mp_sync(gid)

               !only use data that is local to me - i.e. in slice of pencil I own
               IF ((lbounds_local(1) .LE. i) .AND. (i .LE. ubounds_local(1)) .AND. (lbounds_local(2) .LE. j) &
                   .AND. (j .LE. ubounds_local(2))) THEN
                  !allow scaling of external potential values by factor 'scaling' (SCALING_FACTOR in input file)
                  grid%cr3d(i, j, lbounds(3):ubounds(3)) = buffer(lbounds(3):ubounds(3))*scaling
               ENDIF

            ENDDO
         ENDDO

         IF (my_rank == 0) CALL close_file(unit_number=extunit)

         CALL mp_sync(gid)
      ELSE
         ! Parallel routine needs as input the byte size of each grid z-slice
         ! This is a hack to prevent compilation errors with gcc -Wall (upto versions 6.3)
         ! related to allocatable-length string declaration CHARACTER(LEN=:), ALLOCATABLE, DIMENSION(:) :: string
         ! Each data line of a Gaussian cube contains max 6 entries with last line potentially containing less if nz % 6 /= 0
         ! Thus, this size is simply the number of entries multiplied by the entry size + the number of line breaks
         num_linebreak = size_of_z/num_entries_line
         IF (MODULO(size_of_z, num_entries_line) /= 0) &
            num_linebreak = num_linebreak + 1
         msglen = (size_of_z*entry_len + num_linebreak)*mpi_character_size
         CALL cube_to_pw_parallel(grid, filename, scaling, msglen)
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
!> \par History
!>      Created [Nico Holmberg] (05.2017)
! **************************************************************************************************
   SUBROUTINE cube_to_pw_parallel(grid, filename, scaling, msglen)

      TYPE(pw_type), POINTER                             :: grid
      CHARACTER(len=*), INTENT(in)                       :: filename
      REAL(kind=dp), INTENT(in)                          :: scaling
      INTEGER, INTENT(in)                                :: msglen

      CHARACTER(len=*), PARAMETER :: routineN = 'cube_to_pw_parallel', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: entry_len = 13, num_entries_line = 6

      INTEGER, DIMENSION(3)                              :: lbounds, lbounds_local, npoints, &
                                                            npoints_local, ubounds, ubounds_local
      INTEGER, ALLOCATABLE, DIMENSION(:), TARGET         :: blocklengths
      INTEGER(kind=file_offset), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: displacements
      INTEGER(kind=file_offset)                          :: BOF
      INTEGER :: extunit, gid, i, ientry, ip, islice, j, k, master, my_rank, nat, ndum, nslices, &
         num_pe, offset_global, output_unit, pos, readstat, size_of_z, tag
      CHARACTER(LEN=msglen), ALLOCATABLE, DIMENSION(:)   :: readbuffer
      CHARACTER(LEN=msglen)                              :: tmp
      LOGICAL                                            :: should_read(2)
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: buffer
      REAL(kind=dp), DIMENSION(3)                        :: dr, rdum
      TYPE(mp_file_descriptor_type)                      :: mp_file_type

      output_unit = cp_logger_get_default_io_unit()

      !get rs grids and parallel envnment
      gid = grid%pw_grid%para%group
      my_rank = grid%pw_grid%para%my_pos
      num_pe = grid%pw_grid%para%group_size
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
         IF (output_unit > 0) THEN
            WRITE (output_unit, FMT="(/,T2,A,/,/,T2,A,/)") "Reading the cube file:     ", TRIM(filename)
         END IF

         CALL open_file(file_name=filename, &
                        file_status="OLD", &
                        file_form="FORMATTED", &
                        file_action="READ", &
                        file_access="STREAM", &
                        unit_number=extunit)

         !skip header comments
         DO i = 1, 2
            READ (extunit, *)
         ENDDO
         READ (extunit, *) nat, rdum
         DO i = 1, 3
            READ (extunit, *) ndum, rdum
            IF ((ndum /= npoints(i) .OR. (ABS(rdum(i) - dr(i)) > 1e-4)) .AND. &
                output_unit > 0) THEN
               WRITE (output_unit, *) "Restart from density | ERROR! | CUBE FILE NOT COINCIDENT WITH INTERNAL GRID ", i
               WRITE (output_unit, *) "Restart from density | ", ndum, " DIFFERS FROM ", npoints(i)
               WRITE (output_unit, *) "Restart from density | ", rdum, " DIFFERS FROM ", dr(i)
            ENDIF
         ENDDO
         !ignore atomic postion data - read from coord or topology instead
         DO i = 1, nat
            READ (extunit, *)
         ENDDO
         ! Get byte offset
         INQUIRE (extunit, POS=offset_global)
         CALL close_file(unit_number=extunit)
      ENDIF
      ! Sync offset and start parallel read
      CALL mp_bcast(offset_global, grid%pw_grid%para%group_head_id, gid)
      BOF = offset_global
      CALL mp_file_open(groupid=gid, fh=extunit, filepath=filename, amode_status=file_amode_rdonly)
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
      mp_file_type = mp_file_type_hindexed_make_chv(nslices, blocklengths, displacements)
      BOF = 0
      CALL mp_file_type_set_view_chv(extunit, BOF, mp_file_type)
      ! Collective read of cube
      ALLOCATE (readbuffer(nslices))
      readbuffer(:) = ''
      CALL mp_file_read_all_chv(extunit, msglen, nslices, readbuffer, mp_file_type)
      CALL mp_file_type_free(mp_file_type)
      CALL mp_file_close(extunit)
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
         grid%cr3d(i, j, lbounds(3):ubounds(3)) = scaling*buffer(lbounds(3):ubounds(3))
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
            IF (output_unit > 0) THEN
               WRITE (output_unit, FMT="(/,T2,A,/,/,T2,A)") "Reading the cube file:     ", filename
            END IF

            CALL open_file(file_name=filename, &
                           file_status="OLD", &
                           file_form="FORMATTED", &
                           file_action="READ", &
                           unit_number=extunit)

            !skip header comments
            DO i = 1, 2
               READ (extunit, *)
            ENDDO
            READ (extunit, *) nat, rdum
            DO i = 1, 3
               READ (extunit, *) ndum, rdum
               IF ((ndum /= npoints(i) .OR. (ABS(rdum(i) - dr(i)) > 1e-4)) .AND. &
                   output_unit > 0) THEN
                  WRITE (output_unit, *) "Restart from density | ERROR! | CUBE FILE NOT COINCIDENT WITH INTERNAL GRID ", i
                  WRITE (output_unit, *) "Restart from density | ", ndum, " DIFFERS FROM ", npoints(i)
                  WRITE (output_unit, *) "Restart from density | ", rdum, " DIFFERS FROM ", dr(i)
               ENDIF
            ENDDO
            !ignore atomic postion data - read from coord or topology instead
            DO i = 1, nat
               READ (extunit, *)
            ENDDO
         ENDIF

         !master sends all data to everyone
         DO i = lbounds(1), ubounds(1)
            DO j = lbounds(2), ubounds(2)
               IF (my_rank .EQ. 0) THEN
                  READ (extunit, *) (buffer(k), k=lbounds(3), ubounds(3))
                  IF (num_pe .GT. 1) THEN
                     DO ip = 1, num_pe - 1
                        CALL mp_send(buffer(lbounds(3):ubounds(3)), ip, tag, gid)
                     ENDDO
                  ENDIF
               ELSE
                  master = 0
                  CALL mp_recv(buffer(lbounds(3):ubounds(3)), master, tag, gid)
               ENDIF
               CALL mp_sync(gid)

               !only use data that is local to me - i.e. in slice of pencil I own
               IF ((lbounds_local(1) .LE. i) .AND. (i .LE. ubounds_local(1)) .AND. (lbounds_local(2) .LE. j) &
                   .AND. (j .LE. ubounds_local(2))) THEN
                  !allow scaling of external potential values by factor 'scaling' (SCALING_FACTOR in input file)
                  IF (ANY(grid%cr3d(i, j, lbounds(3):ubounds(3)) /= buffer(lbounds(3):ubounds(3))*scaling)) &
                     CALL cp_abort(__LOCATION__, &
                                   "Error in parallel read of input cube file.")
               ENDIF

            ENDDO
         ENDDO

         IF (my_rank == 0) CALL close_file(unit_number=extunit)

         CALL mp_sync(gid)
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

      TYPE(pw_type), POINTER                             :: grid
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(*), INTENT(IN), OPTIONAL                 :: title
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: particles_r
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: particles_z
      INTEGER, INTENT(IN)                                :: stride(3)
      LOGICAL, INTENT(IN)                                :: zero_tails
      INTEGER, INTENT(IN)                                :: msglen

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_to_cube_parallel', &
         routineP = moduleN//':'//routineN
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
      INTEGER                                            :: counter, gid, i, islice, j, k, last_z, &
                                                            my_rank, np, nslices, size_of_z
      CHARACTER(LEN=msglen), ALLOCATABLE, DIMENSION(:)   :: writebuffer
      CHARACTER(LEN=msglen)                              :: tmp
      LOGICAL                                            :: should_write(2)
      TYPE(mp_file_descriptor_type)                      :: mp_file_type

      !get rs grids and parallel envnment
      gid = grid%pw_grid%para%group
      my_rank = grid%pw_grid%para%my_pos

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
      CALL mp_file_get_position(unit_nr, BOF)
      ! Write header information on master process and update byte offset accordingly
      IF (my_rank == 0) THEN
         ! this format seems to work for e.g. molekel and gOpenmol
         ! latest version of VMD can read non orthorhombic cells
         CALL mp_file_write_at(unit_nr, BOF, "-Quickstep-"//NEW_LINE("C"))
         BOF = BOF + LEN("-Quickstep-"//NEW_LINE("C"))*mpi_character_size
         IF (PRESENT(title)) THEN
            CALL mp_file_write_at(unit_nr, BOF, TRIM(title)//NEW_LINE("C"))
            BOF = BOF + LEN(TRIM(title)//NEW_LINE("C"))*mpi_character_size
         ELSE
            CALL mp_file_write_at(unit_nr, BOF, "No Title"//NEW_LINE("C"))
            BOF = BOF + LEN("No Title"//NEW_LINE("C"))*mpi_character_size
         ENDIF

         CPASSERT(PRESENT(particles_z) .EQV. PRESENT(particles_r))
         np = 0
         IF (PRESENT(particles_z)) THEN
            CPASSERT(SIZE(particles_z) == SIZE(particles_r, dim=2))
            ! cube files can only be written for 99999 particles due to a format limitation (I5)
            ! so we limit the number of particles written.
            np = MIN(99999, SIZE(particles_z))
         END IF

         WRITE (header, '(I5,3f12.6)') np, 0.0_dp, 0._dp, 0._dp !start of cube
         CALL mp_file_write_at(unit_nr, BOF, header//NEW_LINE("C"))
         BOF = BOF + LEN(header//NEW_LINE("C"))*mpi_character_size

         WRITE (header, '(I5,3f12.6)') (grid%pw_grid%npts(1) + stride(1) - 1)/stride(1), &
            grid%pw_grid%dh(1, 1)*REAL(stride(1), dp), grid%pw_grid%dh(2, 1)*REAL(stride(1), dp), &
            grid%pw_grid%dh(3, 1)*REAL(stride(1), dp)
         CALL mp_file_write_at(unit_nr, BOF, header//NEW_LINE("C"))
         BOF = BOF + LEN(header//NEW_LINE("C"))*mpi_character_size

         WRITE (header, '(I5,3f12.6)') (grid%pw_grid%npts(2) + stride(2) - 1)/stride(2), &
            grid%pw_grid%dh(1, 2)*REAL(stride(2), dp), grid%pw_grid%dh(2, 2)*REAL(stride(2), dp), &
            grid%pw_grid%dh(3, 2)*REAL(stride(2), dp)
         CALL mp_file_write_at(unit_nr, BOF, header//NEW_LINE("C"))
         BOF = BOF + LEN(header//NEW_LINE("C"))*mpi_character_size

         WRITE (header, '(I5,3f12.6)') (grid%pw_grid%npts(3) + stride(3) - 1)/stride(3), &
            grid%pw_grid%dh(1, 3)*REAL(stride(3), dp), grid%pw_grid%dh(2, 3)*REAL(stride(3), dp), &
            grid%pw_grid%dh(3, 3)*REAL(stride(3), dp)
         CALL mp_file_write_at(unit_nr, BOF, header//NEW_LINE("C"))
         BOF = BOF + LEN(header//NEW_LINE("C"))*mpi_character_size

         IF (PRESENT(particles_z)) THEN
            DO i = 1, np
               WRITE (header_z, '(I5,4f12.6)') particles_z(i), 0._dp, particles_r(:, i)
               CALL mp_file_write_at(unit_nr, BOF, header_z//NEW_LINE("C"))
               BOF = BOF + LEN(header_z//NEW_LINE("C"))*mpi_character_size
            END DO
         END IF
      ENDIF
      ! Sync offset
      CALL mp_bcast(BOF, grid%pw_grid%para%group_head_id, gid)
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
                  IF (zero_tails .AND. grid%cr3d(i, j, k) < 1.E-7_dp) THEN
                     WRITE (value, '(E13.5)') 0.0_dp
                  ELSE
                     WRITE (value, '(E13.5)') grid%cr3d(i, j, k)
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
      mp_file_type = mp_file_type_hindexed_make_chv(nslices, blocklengths, displacements)
      ! Use the created type as a file view
      ! NB. The vector 'displacements' contains the absolute offsets of each z-slice i.e.
      ! they are given relative to the beginning of the file. The global offset to
      ! set_view must therefore be set to 0
      BOF = 0
      CALL mp_file_type_set_view_chv(unit_nr, BOF, mp_file_type)
      ! Collective write of cube
      CALL mp_file_write_all_chv(unit_nr, msglen, nslices, writebuffer, mp_file_type)
      ! Clean up
      CALL mp_file_type_free(mp_file_type)
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
      TYPE(pw_type), POINTER                             :: pw
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: stride
      TYPE(pw_type), OPTIONAL, POINTER                   :: pw2

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_to_simple_volumetric', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: checksum, dest, gid, handle, i, I1, I2, &
                                                            I3, ip, L1, L2, L3, my_rank, &
                                                            my_stride(3), ngrids, npoints, num_pe, &
                                                            rank(2), source, tag, U1, U2, U3
      LOGICAL                                            :: DOUBLE
      REAL(KIND=dp)                                      :: x, y, z
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: buf, buf2

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

      my_rank = pw%pw_grid%para%my_pos
      gid = pw%pw_grid%para%group
      num_pe = pw%pw_grid%para%group_size
      tag = 1

      rank(1) = unit_nr
      rank(2) = my_rank
      checksum = 0
      IF (unit_nr > 0) checksum = 1

      CALL mp_sum(checksum, gid)
      CPASSERT(checksum == 1)

      CALL mp_maxloc(rank, gid)
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
                  ENDIF
               ENDDO
            ELSE
               source = dest
            ENDIF

            IF (source == dest) THEN
               IF (my_rank == source) THEN
                  buf(:) = pw%cr3d(I1, I2, :)
                  IF (DOUBLE) buf2(:) = pw2%cr3d(I1, I2, :)
               ENDIF
            ELSE
               IF (my_rank == source) THEN
                  buf(:) = pw%cr3d(I1, I2, :)
                  CALL mp_send(buf, dest, tag, gid)
                  IF (DOUBLE) THEN
                     buf2(:) = pw2%cr3d(I1, I2, :)
                     CALL mp_send(buf2, dest, tag, gid)
                  ENDIF
               ENDIF
               IF (my_rank == dest) THEN
                  CALL mp_recv(buf, source, tag, gid)
                  IF (DOUBLE) CALL mp_recv(buf2, source, tag, gid)
               ENDIF
            ENDIF

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
                  ENDIF
               ENDDO

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
                  ENDIF
               ENDDO

            ENDIF ! Double

            ! this double loop generates so many messages that it can overload
            ! the message passing system, e.g. on XT3
            ! we therefore put a barrier here that limits the amount of message
            ! that flies around at any given time.
            ! if ever this routine becomes a bottleneck, we should go for a
            ! more complicated rewrite
            CALL mp_sync(gid)

         ENDDO
      ENDDO

      DEALLOCATE (buf)
      IF (DOUBLE) DEALLOCATE (buf2)

      CALL timestop(handle)

   END SUBROUTINE pw_to_simple_volumetric

END MODULE realspace_grid_cube
