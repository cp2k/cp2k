!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Read explicit Wannier90 overlap connections, including nonuniform closed loops.
! **************************************************************************************************
MODULE wannier90_nnkpts
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wannier90_nnkpts'
   PUBLIC :: read_wannier90_nnkpts
CONTAINS

! **************************************************************************************************
!> \brief Read and validate the cell, points, and directed connections of an nnkp file.
!> \param filename input file
!> \param real_lattice direct lattice vectors as columns, in Angstrom
!> \param recip_lattice reciprocal lattice vectors as columns, in inverse Angstrom (including 2*pi)
!> \param kpt fractional k-points, without coordinate wrapping
!> \param nnlist target point for each directed connection
!> \param nncell reciprocal lattice translation for each connection
! **************************************************************************************************
   SUBROUTINE read_wannier90_nnkpts(filename, real_lattice, recip_lattice, kpt, nnlist, nncell)
      CHARACTER(len=*), INTENT(IN)                       :: filename
      REAL(KIND=dp), INTENT(IN)                          :: real_lattice(3, 3), recip_lattice(3, 3)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: kpt(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: nnlist(:, :), nncell(:, :, :)

      CHARACTER(len=1024)                                :: line
      CHARACTER(len=64)                                  :: block, marker
      INTEGER                                            :: i, ios, j, n, nexclude, nneigh, origin, &
                                                            shift(3), TARGET, unit
      INTEGER, ALLOCATABLE                               :: counts(:)
      LOGICAL                                            :: have_real, have_recip
      REAL(KIND=dp)                                      :: lattice(3, 3)

      have_real = .FALSE.
      have_recip = .FALSE.
      CALL open_file(filename, unit_number=unit, file_status="OLD", file_action="READ")
      DO
         CALL next_line(unit, line, ios)
         IF (ios < 0) EXIT
         IF (ios /= 0) CPABORT("NNKP_FILE: read error.")
         marker = ""
         block = ""
         READ (line, *, IOSTAT=ios) marker, block
         IF (ios /= 0) CYCLE
         CALL lowercase(marker)
         CALL lowercase(block)
         IF (marker /= "begin") CYCLE
         SELECT CASE (TRIM(block))
         CASE ("real_lattice", "recip_lattice")
            DO i = 1, 3
               CALL required_line(unit, line)
               READ (line, *, IOSTAT=ios) lattice(:, i)
               IF (ios /= 0) CPABORT("NNKP_FILE: invalid lattice vector.")
            END DO
            IF (.NOT. ALL(ieee_is_finite(lattice))) CPABORT("NNKP_FILE: nonfinite lattice.")
            IF (block == "real_lattice") THEN
               IF (have_real) CPABORT("NNKP_FILE: duplicate real_lattice block.")
               IF (MAXVAL(ABS(lattice - real_lattice)) > 1.e-6_dp) THEN
                  CPABORT("NNKP_FILE: real lattice differs from the CP2K cell (Angstrom).")
               END IF
               have_real = .TRUE.
            ELSE
               IF (have_recip) CPABORT("NNKP_FILE: duplicate recip_lattice block.")
               IF (MAXVAL(ABS(lattice - recip_lattice)) > 1.e-6_dp) THEN
                  CPABORT("NNKP_FILE: reciprocal lattice differs from the CP2K cell.")
               END IF
               have_recip = .TRUE.
            END IF
            CALL end_block(unit, block)
         CASE ("kpoints")
            IF (ALLOCATED(kpt)) CPABORT("NNKP_FILE: duplicate kpoints block.")
            CALL required_line(unit, line)
            READ (line, *, IOSTAT=ios) n
            IF (ios /= 0) CPABORT("NNKP_FILE: invalid point count.")
            IF (n < 1) CPABORT("NNKP_FILE: point count must be positive.")
            ALLOCATE (kpt(3, n))
            DO i = 1, n
               CALL required_line(unit, line)
               READ (line, *, IOSTAT=ios) kpt(:, i)
               IF (ios /= 0) CPABORT("NNKP_FILE: invalid k-point.")
            END DO
            IF (.NOT. ALL(ieee_is_finite(kpt))) CPABORT("NNKP_FILE: nonfinite k-point.")
            CALL end_block(unit, block)
         CASE ("nnkpts")
            IF (.NOT. ALLOCATED(kpt)) CPABORT("NNKP_FILE: kpoints must precede nnkpts.")
            IF (ALLOCATED(nnlist)) CPABORT("NNKP_FILE: duplicate nnkpts block.")
            CALL required_line(unit, line)
            READ (line, *, IOSTAT=ios) nneigh
            IF (ios /= 0) CPABORT("NNKP_FILE: invalid neighbour count.")
            IF (nneigh < 1) CPABORT("NNKP_FILE: neighbour count must be positive.")
            n = SIZE(kpt, 2)
            ALLOCATE (nnlist(n, nneigh), nncell(3, n, nneigh), counts(n))
            counts = 0
            DO i = 1, n*nneigh
               CALL required_line(unit, line)
               READ (line, *, IOSTAT=ios) origin, TARGET, shift
               IF (ios /= 0) CPABORT("NNKP_FILE: invalid connection.")
               IF (MIN(origin, TARGET) < 1 .OR. MAX(origin, TARGET) > n) THEN
                  CPABORT("NNKP_FILE: connection point index out of range.")
               END IF
               IF (counts(origin) == nneigh) CPABORT("NNKP_FILE: too many neighbours for a point.")
               DO j = 1, counts(origin)
                  IF (nnlist(origin, j) == TARGET .AND. ALL(nncell(:, origin, j) == shift)) THEN
                     CPABORT("NNKP_FILE: duplicate connection.")
                  END IF
               END DO
               counts(origin) = counts(origin) + 1
               nnlist(origin, counts(origin)) = TARGET
               nncell(:, origin, counts(origin)) = shift
            END DO
            IF (ANY(counts /= nneigh)) CPABORT("NNKP_FILE: missing connections.")
            DEALLOCATE (counts)
            CALL end_block(unit, block)
         CASE ("exclude_bands")
            CALL required_line(unit, line)
            READ (line, *, IOSTAT=ios) nexclude
            IF (ios /= 0) CPABORT("NNKP_FILE: invalid exclude_bands count.")
            IF (nexclude /= 0) CPABORT("NNKP_FILE: specify band exclusions using CP2K EXCLUDE_BANDS.")
            CALL end_block(unit, block)
         CASE DEFAULT
            ! Skip projections and other Wannier90 blocks not required for overlaps.
            DO
               CALL required_line(unit, line)
               marker = ""
               READ (line, *, IOSTAT=ios) marker
               IF (ios /= 0) CPABORT("NNKP_FILE: invalid line in unused block.")
               CALL lowercase(marker)
               IF (marker == "end") EXIT
            END DO
         END SELECT
      END DO
      CALL close_file(unit)
      IF (.NOT. have_real .OR. .NOT. have_recip) CPABORT("NNKP_FILE: missing lattice blocks.")
      IF (.NOT. ALLOCATED(kpt) .OR. .NOT. ALLOCATED(nnlist)) THEN
         CPABORT("NNKP_FILE: missing kpoints or nnkpts block.")
      END IF
   END SUBROUTINE read_wannier90_nnkpts

! **************************************************************************************************
!> \brief Read a nonempty line, removing comments.
!> \param unit input unit
!> \param line resulting line
!> \param ios IO status
! **************************************************************************************************
   SUBROUTINE next_line(unit, line, ios)
      INTEGER, INTENT(IN)                                :: unit
      CHARACTER(len=*), INTENT(OUT)                      :: line
      INTEGER, INTENT(OUT)                               :: ios

      INTEGER                                            :: i

      DO
         READ (unit, '(A)', IOSTAT=ios) line
         IF (ios /= 0) RETURN
         i = SCAN(line, "!#")
         IF (i > 0) line(i:) = ""
         line = ADJUSTL(line)
         IF (LEN_TRIM(line) > 0) RETURN
      END DO
   END SUBROUTINE next_line

! **************************************************************************************************
!> \brief Read a required nonempty line.
!> \param unit input unit
!> \param line resulting line
! **************************************************************************************************
   SUBROUTINE required_line(unit, line)
      INTEGER, INTENT(IN)                                :: unit
      CHARACTER(len=*), INTENT(OUT)                      :: line

      INTEGER                                            :: ios

      CALL next_line(unit, line, ios)
      IF (ios /= 0) CPABORT("NNKP_FILE: unexpected end of file or read error.")
   END SUBROUTINE required_line

! **************************************************************************************************
!> \brief Check a block terminator.
!> \param unit input unit
!> \param block expected block name
! **************************************************************************************************
   SUBROUTINE end_block(unit, block)
      INTEGER, INTENT(IN)                                :: unit
      CHARACTER(len=*), INTENT(IN)                       :: block

      CHARACTER(len=1024)                                :: line
      CHARACTER(len=64)                                  :: marker, name
      INTEGER                                            :: ios

      CALL required_line(unit, line)
      READ (line, *, IOSTAT=ios) marker, name
      IF (ios /= 0) CPABORT("NNKP_FILE: missing block terminator.")
      CALL lowercase(marker)
      CALL lowercase(name)
      IF (marker /= "end" .OR. name /= block) CPABORT("NNKP_FILE: incorrect block terminator.")
   END SUBROUTINE end_block

! **************************************************************************************************
!> \brief Convert ASCII keyword to lowercase.
!> \param string keyword
! **************************************************************************************************
   SUBROUTINE lowercase(string)
      CHARACTER(len=*), INTENT(INOUT)                    :: string

      INTEGER                                            :: code, i

      DO i = 1, LEN_TRIM(string)
         code = IACHAR(string(i:i))
         IF (code >= IACHAR('A') .AND. code <= IACHAR('Z')) string(i:i) = ACHAR(code + 32)
      END DO
   END SUBROUTINE lowercase
END MODULE wannier90_nnkpts
