!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
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
        file_amode_rdonly, file_offset, mp_comm_type, mp_file_descriptor_type, mp_file_type, &
        mp_file_type_free, mp_file_type_hindexed_make_chv, mp_file_type_set_view_chv, &
        mpi_character_size
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
   USE pw_types,                        ONLY: pw_r3d_rs_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: cube_reader_type, cube_value_format, cube_values_format, cube_to_pw, pw_to_cube, &
             pw_to_simple_volumetric

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'realspace_grid_cube'
   INTEGER, PARAMETER, PRIVATE          :: cube_entry_len = 14, &
                                           cube_num_entries_line = 6
   CHARACTER(len=*), PARAMETER          :: cube_value_format = '(1X,ES13.5E3)', &
                                           cube_values_format = '(6(1X,ES13.5E3))'
   CHARACTER(len=*), PARAMETER, PRIVATE :: cube_whitespace = ' '//ACHAR(9)//ACHAR(10)//ACHAR(11)//ACHAR(12)//ACHAR(13)
   CHARACTER(len=*), PARAMETER, PRIVATE :: missing_zeff_warning = &
                                           "Effective nuclear charges were not supplied; "// &
                                           "zeros will be written to the Cube atom records."

   TYPE cube_reader_type
      PRIVATE
      CHARACTER(LEN=:), ALLOCATABLE :: line
      INTEGER                       :: pos = 1, length = 0, complete = 0
   CONTAINS
      PROCEDURE :: read => cube_read_values
   END TYPE cube_reader_type

CONTAINS

! **************************************************************************************************
!> \brief Read whitespace-separated cube values, retaining any unused part of the current record.
!> \param reader input record and position
!> \param unit_nr formatted file positioned at the cube data
!> \param buffer parsed values
! **************************************************************************************************
   SUBROUTINE cube_read_values(reader, unit_nr, buffer)
      CLASS(cube_reader_type), INTENT(INOUT)             :: reader
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: buffer

      CHARACTER(LEN=:), ALLOCATABLE                      :: longer_line
      INTEGER                                            :: count, first, nchar, nvalues, readstat

      IF (.NOT. ALLOCATED(reader%line)) ALLOCATE (CHARACTER(LEN=512) :: reader%line)
      nvalues = 0
      DO WHILE (nvalues < SIZE(buffer))
         IF (reader%pos > reader%complete) THEN
            ! Carry just the unfinished token into the next read, even for a single very long line.
            reader%line(:reader%length - reader%complete) = reader%line(reader%complete + 1:reader%length)
            reader%length = reader%length - reader%complete
            DO
               IF (reader%length == LEN(reader%line)) THEN
                  ALLOCATE (CHARACTER(LEN=2*reader%length) :: longer_line)
                  longer_line(:reader%length) = reader%line
                  CALL MOVE_ALLOC(longer_line, reader%line)
               END IF
               READ (unit_nr, '(A)', ADVANCE='NO', SIZE=nchar, IOSTAT=readstat) reader%line(reader%length + 1:)
               IF (readstat > 0) CPABORT("Error reading cube data.")
               reader%length = reader%length + nchar
               IF (readstat < 0) THEN
                  IF (IS_IOSTAT_END(readstat) .AND. reader%length == 0) THEN
                     CPABORT("Unexpected end of cube data.")
                  END IF
                  reader%complete = reader%length
                  EXIT
               END IF
               reader%complete = SCAN(reader%line(:reader%length), cube_whitespace, BACK=.TRUE.)
               IF (reader%complete > 0) EXIT
            END DO
            reader%pos = 1
         END IF

         first = reader%pos
         count = cube_value_span(reader%line(:reader%complete), reader%pos, SIZE(buffer) - nvalues)
         CALL cube_parse_values(reader%line(first:reader%pos - 1), buffer(nvalues + 1:nvalues + count))
         nvalues = nvalues + count
      END DO

   END SUBROUTINE cube_read_values

! **************************************************************************************************
!> \brief Count up to limit values and normalize whitespace for list-directed internal reads.
!> \param text cube data
!> \param pos current character, advanced to the next unread value or end
!> \param limit maximum number of values to consume
!> \return number of values in the span
! **************************************************************************************************
   INTEGER FUNCTION cube_value_span(text, pos, limit) RESULT(count)
      CHARACTER(LEN=*), INTENT(INOUT)                    :: text
      INTEGER, INTENT(INOUT)                             :: pos
      INTEGER, INTENT(IN)                                :: limit

      LOGICAL                                            :: in_value

      count = 0
      in_value = .FALSE.
      DO WHILE (pos <= LEN(text))
         SELECT CASE (text(pos:pos))
         CASE (' ', ACHAR(9):ACHAR(13))
            text(pos:pos) = ' '
            in_value = .FALSE.
         CASE DEFAULT
            IF (.NOT. in_value) THEN
               IF (count == limit) EXIT
               count = count + 1
               in_value = .TRUE.
            END IF
         END SELECT
         pos = pos + 1
      END DO
   END FUNCTION cube_value_span

! **************************************************************************************************
!> \brief Convert a span of whitespace-separated values.
!> \param text cube data
!> \param values parsed values
! **************************************************************************************************
   SUBROUTINE cube_parse_values(text, values)
      CHARACTER(LEN=*), INTENT(IN)                       :: text
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: values

      INTEGER                                            :: stat

      IF (SIZE(values) == 0) RETURN
      ! Null values, repeat counts and list terminators would silently alter the grid.
      IF (SCAN(text, ',/*;!') /= 0) CPABORT("Bad value while reading cube data.")
      READ (text, *, IOSTAT=stat) values
      IF (stat /= 0) CPABORT("Bad value while reading cube data.")
   END SUBROUTINE cube_parse_values

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param unit_nr ...
!> \param title ...
!> \param particles_r ...
!> \param particles_z ...
!> \param particles_zeff ...
!> \param stride ...
!> \param max_file_size_mb ...
!> \param zero_tails ...
!> \param silent ...
!> \param mpi_io ...
! **************************************************************************************************
   SUBROUTINE pw_to_cube(pw, unit_nr, title, particles_r, particles_z, particles_zeff, &
                         stride, max_file_size_mb, zero_tails, silent, mpi_io)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: pw
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(*), INTENT(IN), OPTIONAL                 :: title
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: particles_r
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: particles_z
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: particles_zeff
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: stride
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: max_file_size_mb
      LOGICAL, INTENT(IN), OPTIONAL                      :: zero_tails, silent, mpi_io

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_to_cube'

      INTEGER :: checksum, dest, handle, i, I1, I2, I3, iat, ip, L1, L2, L3, msglen, my_rank, &
         my_stride(3), np, num_linebreak, num_pe, rank(2), size_of_z, source, tag, U1, U2, U3
      LOGICAL                                            :: be_silent, my_zero_tails, parallel_write
      REAL(KIND=dp)                                      :: compression_factor, my_max_file_size_mb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: buf
      TYPE(mp_comm_type)                                 :: gid
      TYPE(mp_file_type)                                 :: mp_unit

      CALL timeset(routineN, handle)

      my_zero_tails = .FALSE.
      be_silent = .FALSE.
      parallel_write = .FALSE.
      my_max_file_size_mb = 0.0_dp
      IF (PRESENT(zero_tails)) my_zero_tails = zero_tails
      IF (PRESENT(silent)) be_silent = silent
      IF (PRESENT(mpi_io)) parallel_write = mpi_io
      IF (PRESENT(max_file_size_mb)) my_max_file_size_mb = max_file_size_mb
      CPASSERT(my_max_file_size_mb >= 0)
      IF (.NOT. parallel_write .AND. PRESENT(particles_z) .AND. &
          .NOT. PRESENT(particles_zeff)) THEN
         CPWARN(missing_zeff_warning)
      END IF

      my_stride = 1
      IF (PRESENT(stride)) THEN
         IF (SIZE(stride) /= 1 .AND. SIZE(stride) /= 3) THEN
            CALL cp_abort(__LOCATION__, "STRIDE keyword can accept only 1 "// &
                          "(the same for X,Y,Z) or 3 values. Correct your input file.")
         END IF
         IF (SIZE(stride) == 1) THEN
            DO i = 1, 3
               my_stride(i) = stride(1)
            END DO
         ELSE
            my_stride = stride(1:3)
         END IF
      END IF

      IF (my_max_file_size_mb > 0) THEN
         ! Estimate the grid data size from the width of a value in bytes.
         compression_factor = 1.0e-06_dp*cube_entry_len*PRODUCT(REAL(pw%pw_grid%npts, dp))/my_max_file_size_mb
         my_stride(:) = INT(compression_factor**(1.0_dp/3.0_dp)) + 1
      END IF

      CPASSERT(my_stride(1) > 0)
      CPASSERT(my_stride(2) > 0)
      CPASSERT(my_stride(3) > 0)

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
               IF (PRESENT(particles_zeff)) THEN
                  DO iat = 1, np
                     WRITE (unit_nr, '(I5,4f12.6)') particles_z(iat), particles_zeff(iat), particles_r(:, iat)
                  END DO
               ELSE
                  DO iat = 1, np
                     WRITE (unit_nr, '(I5,4f12.6)') particles_z(iat), 0._dp, particles_r(:, iat)
                  END DO
               END IF
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
               IF (pw%pw_grid%para%mode /= PW_MODE_LOCAL) THEN
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
                  WRITE (unit_nr, cube_values_format) (buf(I3), I3=L3, U3, my_stride(3))
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
         num_linebreak = (size_of_z + cube_num_entries_line - 1)/cube_num_entries_line
         msglen = (size_of_z*cube_entry_len + num_linebreak)*mpi_character_size
         CALL mp_unit%set_handle(unit_nr)
         CALL pw_to_cube_parallel(pw, mp_unit, title, particles_r, particles_z, particles_zeff, &
                                  my_stride, my_zero_tails, msglen)
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
      CHARACTER(len=*), INTENT(IN)                       :: filename
      REAL(KIND=dp), INTENT(IN)                          :: scaling
      LOGICAL, INTENT(IN)                                :: parallel_read
      LOGICAL, INTENT(IN), OPTIONAL                      :: silent

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cube_to_pw'

      INTEGER                                            :: extunit, handle, i, j, my_rank, nat, &
                                                            ndum, output_unit
      INTEGER(KIND=file_offset)                          :: data_position
      INTEGER, DIMENSION(3)                              :: lbounds, lbounds_local, npoints, &
                                                            ubounds, ubounds_local
      LOGICAL                                            :: be_silent
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: buffer
      REAL(KIND=dp), DIMENSION(3)                        :: rdum
      TYPE(cube_reader_type)                             :: reader
      TYPE(mp_comm_type)                                 :: gid

      CALL timeset(routineN, handle)

      output_unit = cp_logger_get_default_io_unit()
      be_silent = .FALSE.
      IF (PRESENT(silent)) be_silent = silent
      gid = grid%pw_grid%para%group
      my_rank = gid%mepos
      npoints = grid%pw_grid%npts
      lbounds = grid%pw_grid%bounds(1, :)
      ubounds = grid%pw_grid%bounds(2, :)
      lbounds_local = grid%pw_grid%bounds_local(1, :)
      ubounds_local = grid%pw_grid%bounds_local(2, :)

      IF (my_rank == 0) THEN
         IF (output_unit > 0 .AND. .NOT. be_silent) THEN
            WRITE (output_unit, FMT="(/,T2,A,/,/,T2,A,/)") "Reading the cube file:     ", TRIM(filename)
         END IF
         CALL open_file(file_name=filename, file_status="OLD", file_form="FORMATTED", &
                        file_action="READ", file_access="STREAM", unit_number=extunit)

         ! Skip comment lines and read the grid dimensions.
         READ (extunit, *)
         READ (extunit, *)
         READ (extunit, *) nat, rdum
         DO i = 1, 3
            READ (extunit, *) ndum, rdum
            IF ((ndum /= npoints(i) .OR. ABS(rdum(i) - grid%pw_grid%dh(i, i)) > 1e-4_dp) .AND. &
                output_unit > 0) THEN
               WRITE (output_unit, *) "Restart from density | ERROR! | CUBE FILE NOT COINCIDENT WITH INTERNAL GRID ", i
               WRITE (output_unit, *) "Restart from density | ", ndum, " DIFFERS FROM ", npoints(i)
               WRITE (output_unit, *) "Restart from density | ", rdum, " DIFFERS FROM ", grid%pw_grid%dh(i, i)
            END IF
         END DO
         DO i = 1, nat
            READ (extunit, *)
         END DO
      END IF

      IF (parallel_read) THEN
         IF (my_rank == 0) THEN
            INQUIRE (UNIT=extunit, POS=data_position)
            CALL close_file(unit_number=extunit)
         END IF
         CALL gid%bcast(data_position, 0)
         ! Fortran stream positions start at 1 and MPI file offsets start at 0.
         CALL cube_to_pw_parallel(grid, filename, scaling, data_position - 1)
      ELSE
         ALLOCATE (buffer(lbounds(3):ubounds(3)))
         DO i = lbounds(1), ubounds(1)
            DO j = lbounds(2), ubounds(2)
               IF (my_rank == 0) CALL reader%read(extunit, buffer)
               CALL gid%bcast(buffer, 0)
               IF (lbounds_local(1) <= i .AND. i <= ubounds_local(1) .AND. &
                   lbounds_local(2) <= j .AND. j <= ubounds_local(2)) THEN
                  grid%array(i, j, lbounds(3):ubounds(3)) = scaling*buffer
               END IF
            END DO
         END DO
         IF (my_rank == 0) CALL close_file(unit_number=extunit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE cube_to_pw

! **************************************************************************************************
!> \brief Read consecutive file chunks and distribute their values to the owning grid pencils.
!> \param grid target grid
!> \param filename name of cube file
!> \param scaling scale values before storing
!> \param data_offset byte offset after the cube header
! **************************************************************************************************
   SUBROUTINE cube_to_pw_parallel(grid, filename, scaling, data_offset)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: grid
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      REAL(KIND=dp), INTENT(IN)                          :: scaling
      INTEGER(KIND=file_offset), INTENT(IN)              :: data_offset

      CHARACTER(LEN=:), ALLOCATABLE                      :: text
      INTEGER                                            :: count, dest, i, j, n, np, nz, packed, &
                                                            pos, rank, received, xhi, xlo
      INTEGER(KIND=file_offset)                          :: base, block_size, chunk_end, consumed, &
                                                            file_size, first, hi, lo, plane, &
                                                            plane_size, total
      INTEGER(KIND=file_offset), ALLOCATABLE, &
         DIMENSION(:)                                    :: counts, offsets
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rcount, rdispl, scount, sdispl
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: local, sendbuf, values
      TYPE(mp_comm_type)                                 :: gid
      TYPE(mp_file_type)                                 :: file

      gid = grid%pw_grid%para%group
      rank = gid%mepos
      np = gid%num_pe
      nz = grid%pw_grid%npts(3)
      plane_size = INT(grid%pw_grid%npts(2), file_offset)*nz
      total = grid%pw_grid%ngpts
      ALLOCATE (offsets(0:np), counts(0:np - 1), local(SIZE(grid%array)))
      ALLOCATE (rcount(0:np - 1), rdispl(0:np - 1), scount(0:np - 1), sdispl(0:np - 1))

      CALL file%open(groupid=gid, filepath=filename, amode_status=file_amode_rdonly)
      CALL file%get_size(file_size)
      ! Read up to 8 MiB per rank per round, extending each chunk to a token boundary.
      block_size = MIN(8_file_offset*1024*1024, (file_size - data_offset + np - 1)/np)
      base = data_offset
      consumed = 0
      received = 0
      DO WHILE (consumed < total)
         IF (base >= file_size) CPABORT("Unexpected end of cube data.")
         chunk_end = MIN(file_size, base + (rank + 1)*block_size)
         offsets(0) = base
         CALL gid%allgather(cube_chunk_end(file, chunk_end, file_size), offsets(1:))
         ALLOCATE (CHARACTER(LEN=offsets(rank + 1) - offsets(rank)) :: text)
         CALL file%read_at_all(offsets(rank), text)
         pos = 1
         count = cube_value_span(text, pos, LEN(text))
         CALL gid%allgather(INT(count, file_offset), counts)
         ! An appended frame may start in this chunk. Read only the requested grid.
         first = MIN(consumed + SUM(counts(:rank - 1)), total)
         n = INT(MIN(INT(count, file_offset), total - first))
         IF (n < count) THEN
            pos = 1
            n = cube_value_span(text, pos, n)
         END IF
         ALLOCATE (values(n), sendbuf(n))
         CALL cube_parse_values(text(:pos - 1), values)
         DEALLOCATE (text)

         ! At fixed x, each destination owns a contiguous run of y and z in file order.
         xlo = INT(first/plane_size) + 1
         xhi = INT((first + n - 1)/plane_size) + 1
         packed = 0
         DO dest = 0, np - 1
            sdispl(dest) = packed
            ASSOCIATE (bo => grid%pw_grid%para%bo(:, :, dest, 1))
               DO i = MAX(bo(1, 1), xlo), MIN(bo(2, 1), xhi)
                  plane = (i - 1)*plane_size
                  lo = MAX(first, plane + INT(bo(1, 2) - 1, file_offset)*nz)
                  hi = MIN(first + n, plane + INT(bo(2, 2), file_offset)*nz)
                  IF (hi <= lo) CYCLE
                  sendbuf(packed + 1:packed + hi - lo) = values(lo - first + 1:hi - first)
                  packed = packed + INT(hi - lo)
               END DO
            END ASSOCIATE
            scount(dest) = packed - sdispl(dest)
         END DO
         CALL gid%alltoall(scount, rcount, 1)
         rdispl(0) = 0
         DO i = 1, np - 1
            rdispl(i) = rdispl(i - 1) + rcount(i - 1)
         END DO
         CALL gid%alltoall(sendbuf, scount, sdispl, local(received + 1:), rcount, rdispl)
         received = received + SUM(rcount)
         consumed = consumed + SUM(counts)
         base = offsets(np)
         DEALLOCATE (sendbuf, values)
      END DO
      CALL file%close()
      ! Received chunks are already in file order, with z varying fastest.
      pos = 0
      ASSOCIATE (bounds => grid%pw_grid%bounds_local)
         DO i = bounds(1, 1), bounds(2, 1)
            DO j = bounds(1, 2), bounds(2, 2)
               grid%array(i, j, :) = scaling*local(pos + 1:pos + nz)
               pos = pos + nz
            END DO
         END DO
      END ASSOCIATE
   END SUBROUTINE cube_to_pw_parallel

! **************************************************************************************************
!> \brief Advance a byte offset to a value separator or EOF.
!> \param file open cube file
!> \param offset nominal end of a chunk
!> \param file_size size in bytes
!> \return exclusive end offset
! **************************************************************************************************
   FUNCTION cube_chunk_end(file, offset, file_size) RESULT(last)
      TYPE(mp_file_type), INTENT(IN)                     :: file
      INTEGER(KIND=file_offset), INTENT(IN)              :: offset, file_size
      INTEGER(KIND=file_offset)                          :: last

      CHARACTER(LEN=128)                                 :: probe
      INTEGER                                            :: n, pos

      last = offset
      DO WHILE (last < file_size)
         n = INT(MIN(INT(LEN(probe), file_offset), file_size - last))
         CALL file%read_at(last, probe(:n))
         pos = SCAN(probe(:n), cube_whitespace)
         IF (pos > 0) THEN
            last = last + pos - 1
            EXIT
         END IF
         last = last + n
      END DO
   END FUNCTION cube_chunk_end

! **************************************************************************************************
!> \brief Writes a realspace potential to a cube file using collective MPI I/O.
!> \param grid        the pw to output to the cube file
!> \param unit_nr     the handle associated with the cube file
!> \param title       title of the cube file
!> \param particles_r Cartersian coordinates of the system
!> \param particles_z atomic masses of atoms in the system
!> \param particles_zeff effective atomic charges of atoms in the system
!> \param stride      every stride(i)th value of the potential is outputted (i=x,y,z)
!> \param zero_tails  flag that determines if small values of the potential should be zeroed
!> \param msglen      the size of each grid slice along z-axis in bytes
!> \par History
!>      Created [Nico Holmberg] (11.2017)
! **************************************************************************************************
   SUBROUTINE pw_to_cube_parallel(grid, unit_nr, title, particles_r, particles_z, particles_zeff, &
                                  stride, zero_tails, msglen)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: grid
      TYPE(mp_file_type), INTENT(IN)                     :: unit_nr
      CHARACTER(*), INTENT(IN), OPTIONAL                 :: title
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: particles_r
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: particles_z
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: particles_zeff
      INTEGER, INTENT(IN)                                :: stride(3)
      LOGICAL, INTENT(IN)                                :: zero_tails
      INTEGER, INTENT(IN)                                :: msglen

      INTEGER, PARAMETER                                 :: header_len = 41, header_len_z = 53

      CHARACTER(LEN=cube_entry_len)                      :: value
      CHARACTER(LEN=header_len)                          :: header
      CHARACTER(LEN=header_len_z)                        :: header_z
      INTEGER, DIMENSION(3)                              :: lbounds, lbounds_local, ubounds, &
                                                            ubounds_local
      INTEGER, ALLOCATABLE, DIMENSION(:), TARGET         :: blocklengths
      INTEGER(kind=file_offset), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: displacements
      INTEGER(kind=file_offset)                          :: BOF
      INTEGER                                            :: counter, i, islice, j, k, my_rank, np, &
                                                            nslices
      CHARACTER(LEN=msglen), ALLOCATABLE, DIMENSION(:)   :: writebuffer
      CHARACTER(LEN=msglen)                              :: tmp
      LOGICAL                                            :: should_write(2)
      TYPE(mp_comm_type)                                 :: gid
      TYPE(mp_file_descriptor_type)                      :: mp_desc

      !get rs grids and parallel envnment
      gid = grid%pw_grid%para%group
      my_rank = grid%pw_grid%para%group%mepos
      IF (PRESENT(particles_z) .AND. .NOT. PRESENT(particles_zeff)) THEN
         CPWARN(missing_zeff_warning)
      END IF

      ! Shortcut
      lbounds = grid%pw_grid%bounds(1, :)
      ubounds = grid%pw_grid%bounds(2, :)
      lbounds_local = grid%pw_grid%bounds_local(1, :)
      ubounds_local = grid%pw_grid%bounds_local(2, :)
      ! Determine the number of local z-slices.
      islice = 1
      DO i = lbounds(1), ubounds(1), stride(1)
         should_write(:) = .TRUE.
         IF (i < lbounds_local(1)) THEN
            should_write(1) = .FALSE.
         ELSE IF (i > ubounds_local(1)) THEN
            EXIT
         END IF
         DO j = lbounds(2), ubounds(2), stride(2)
            should_write(2) = .TRUE.
            IF (j < lbounds_local(2) .OR. j > ubounds_local(2)) THEN
               should_write(2) = .FALSE.
            END IF
            IF (ALL(should_write .EQV. .TRUE.)) THEN
               islice = islice + 1
            END IF
         END DO
      END DO
      nslices = islice - 1
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
            IF (PRESENT(particles_zeff)) THEN
               DO i = 1, np
                  WRITE (header_z, '(I5,4f12.6)') particles_z(i), particles_zeff(i), particles_r(:, i)
                  CALL unit_nr%write_at(BOF, header_z//NEW_LINE("C"))
                  BOF = BOF + LEN(header_z//NEW_LINE("C"))*mpi_character_size
               END DO
            ELSE
               DO i = 1, np
                  WRITE (header_z, '(I5,4f12.6)') particles_z(i), 0._dp, particles_r(:, i)
                  CALL unit_nr%write_at(BOF, header_z//NEW_LINE("C"))
                  BOF = BOF + LEN(header_z//NEW_LINE("C"))*mpi_character_size
               END DO
            END IF
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
         IF (i < lbounds_local(1)) THEN
            should_write(1) = .FALSE.
         ELSE IF (i > ubounds_local(1)) THEN
            EXIT
         END IF
         DO j = lbounds(2), ubounds(2), stride(2)
            should_write(2) = .TRUE.
            IF (j < lbounds_local(2) .OR. j > ubounds_local(2)) THEN
               should_write(2) = .FALSE.
            END IF
            IF (ALL(should_write .EQV. .TRUE.)) THEN
               IF (islice > nslices) CPABORT("Index out of bounds.")
               displacements(islice) = BOF
               tmp = ''
               counter = 0
               DO k = lbounds(3), ubounds(3), stride(3)
                  IF (zero_tails .AND. grid%array(i, j, k) < 1.E-7_dp) THEN
                     WRITE (value, cube_value_format) 0.0_dp
                  ELSE
                     WRITE (value, cube_value_format) grid%array(i, j, k)
                  END IF
                  tmp = TRIM(tmp)//TRIM(value)
                  counter = counter + 1
                  IF (MODULO(counter, cube_num_entries_line) == 0 .OR. ubounds(3) - k < stride(3)) THEN
                     tmp = TRIM(tmp)//NEW_LINE('C')
                  END IF
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
         IF (SIZE(stride) /= 1 .AND. SIZE(stride) /= 3) THEN
            CALL cp_abort(__LOCATION__, "STRIDE keyword can accept only 1 "// &
                          "(the same for X,Y,Z) or 3 values. Correct your input file.")
         END IF
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
            IF (pw%pw_grid%para%mode /= PW_MODE_LOCAL) THEN
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
                     WRITE (unit_nr, '(6(1X,ES12.4E3), 6(1X,ES12.4E3), 6(1X,ES12.4E3), 6(1X,ES12.4E3))') x, y, z, buf(I3)
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
                     WRITE (unit_nr, '(6(1X,ES12.4E3), 6(1X,ES12.4E3), 6(1X,ES12.4E3), 6(1X,ES12.4E3))') x, y, z, buf(I3), buf2(I3)
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

END MODULE realspace_grid_cube
