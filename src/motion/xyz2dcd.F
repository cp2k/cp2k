!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM xyz2dcd

! Version: 1.0
! Author:  Matthias Krack (MK)
! History: - Creation (30.03.2015,MK)
!
! Note: The input coordinates and the cell vectors should be in Angstrom.

! Uncomment the following line if this module is available (e.g. with gfortran) and comment the corresponding variable declarations below
! USE ISO_FORTRAN_ENV, ONLY: error_unit,input_unit,output_unit

   IMPLICIT NONE

   ! Comment the following lines if the ISO_FORTRAN_ENV is used (see above)
   INTEGER, PARAMETER                                 :: default_error_unit = 0, &
                                                         default_input_unit = 5, &
                                                         default_output_unit = 6
   INTEGER                                            :: error_unit = default_error_unit, &
                                                         output_unit = default_output_unit
   ! End Comment

   ! Parameters
   CHARACTER(LEN=*), PARAMETER :: routineN = "xyz2dcd", &
                                  version_info = routineN//" v1.0 (30.06.2020, Matthias Krack)"

   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 200), &
                         sp = SELECTED_REAL_KIND(6, 30)
   INTEGER, PARAMETER :: default_string_length = 240, &
                         cell_file_unit = 10, &
                         dcd_file_unit = 11, &
                         xyz_file_unit = default_input_unit

   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp
!  REAL(KIND=dp), PARAMETER :: angstrom = 0.52917720859_dp            ! [a.u.] -> [Angstrom]
   REAL(KIND=dp), PARAMETER :: degree = 180.0_dp/pi                   ! [rad]  -> [degree]

   ! Variables
  CHARACTER(LEN=default_string_length)               :: arg, cell_file_name, dcd_file_name, message, remark1, remark2, remark_xyz, &
                                                         string, xyz_file_name
   CHARACTER(LEN=5), DIMENSION(:), ALLOCATABLE        :: atomic_label
   CHARACTER(LEN=80), DIMENSION(2)                    :: remark_dcd
   INTEGER                                            :: first_frame, i, iarg, iatom, iskip, istat, istep, j, &
                                                         last_frame, narg, natom, nframe, nframe_read, nremark, stride
   LOGICAL                                            :: apply_pbc, debug, dump_frame, exists, have_cell_file, have_cell_info, &
                                                         info, pbc0, print_atomic_displacements, print_scaled_coordinates, &
                                                         print_scaled_pbc_coordinates, trace_atoms
   REAL(KIND=dp)                                      :: alpha, beta, dt, eps_out_of_box, gamma, tstep
   REAL(KIND=dp), DIMENSION(3)                        :: a, abc, b, c
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: atomic_displacement
   REAL(KIND=dp), DIMENSION(3, 3)                      :: h, hinv
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE         :: r, r_pbc, r0, s, s_pbc

   apply_pbc = .FALSE.
   debug = .FALSE.
   dump_frame = .TRUE.
   have_cell_file = .FALSE.
   have_cell_info = .FALSE.
   info = .FALSE.
   trace_atoms = .FALSE.
   pbc0 = .FALSE.
   print_scaled_coordinates = .FALSE.
   print_atomic_displacements = .FALSE.
   print_scaled_pbc_coordinates = .FALSE.
   first_frame = 1
   last_frame = 1000000 ! Hard limit of 1 Mio frames in total
   stride = 1
   nframe = 0
   nframe_read = 0
   nremark = 0
   cell_file_name = ""
   dcd_file_name = ""
   xyz_file_name = ""
   remark_dcd(:) = ""
   remark_xyz = ""
   eps_out_of_box = -HUGE(0.0_dp)
   h(:, :) = 0.0_dp
   hinv(:, :) = 0.0_dp

   ! Scan argument list and digest it

   narg = command_argument_count()

   IF (narg == 0) THEN
      CALL print_help()
      CALL abort_program(routineN, "No input file(s) specified")
   END IF

   iarg = 0
   arg_loop: DO
      iarg = iarg + 1
      CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)
      SELECT CASE (arg)
      CASE ("-abc")
         DO i = 1, 3
            iarg = iarg + 1
            CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)
            READ (UNIT=arg, FMT=*, ERR=100) h(i, i)
         END DO
         have_cell_info = .TRUE.
         CYCLE arg_loop
100      CALL abort_program(routineN, "Reading -abc arguments (3 reals are expected)")
      CASE ("-cell")
         DO i = 1, 3
            DO j = 1, 3
               iarg = iarg + 1
               CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)
               READ (UNIT=arg, FMT=*, ERR=101) h(j, i)
            END DO
         END DO
         have_cell_info = .TRUE.
         CYCLE arg_loop
101      CALL abort_program(routineN, "Reading -cell arguments (9 reals are expected)")
      CASE ("-cell_file", "-cf")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=cell_file_name, STATUS=istat)
         have_cell_file = .TRUE.
         have_cell_info = .TRUE.
         CYCLE arg_loop
      CASE ("-df", "-dcd_file")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=dcd_file_name, STATUS=istat)
         CYCLE arg_loop
      CASE ("-debug", "-d")
         debug = .TRUE.
         info = .TRUE.
         CYCLE arg_loop
      CASE ("-displacements", "-disp")
         print_atomic_displacements = .TRUE.
         CYCLE arg_loop
      CASE ("-eo")
         error_unit = output_unit
         CYCLE arg_loop
      CASE ("-first_frame", "-first", "-ff")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)
         READ (UNIT=arg, FMT=*, ERR=102) first_frame
         IF (first_frame <= 0) THEN
            CALL abort_program(routineN, "Invalid number for first frame specified: "// &
                               "first_frame must be greater than zero")
         END IF
         CYCLE arg_loop
102      CALL abort_program(routineN, "Invalid number for first frame specified "// &
                            "(an integer number greater than zero is expected)")
      CASE ("-help", "-h")
         CALL print_help()
         STOP
      CASE ("-info", "-i")
         info = .TRUE.
         CYCLE arg_loop
      CASE ("-last_frame", "-last", "-lf")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)
         READ (UNIT=arg, FMT=*, ERR=103) last_frame
         IF (last_frame <= 0) THEN
            CALL abort_program(routineN, "Invalid number for last frame specified: "// &
                               "last_frame must be greater than zero")
         END IF
         CYCLE arg_loop
103      CALL abort_program(routineN, "Invalid number for last frame specified "// &
                            "(an integer number greater than zero is expected)")
      CASE ("-pbc")
         apply_pbc = .TRUE.
         pbc0 = .FALSE.
         CYCLE arg_loop
      CASE ("-pbc0")
         apply_pbc = .TRUE.
         pbc0 = .TRUE.
         CYCLE arg_loop
      CASE ("-scaled_coordinates", "-sc")
         print_scaled_coordinates = .TRUE.
         CYCLE arg_loop
      CASE ("-scaled_pbc_coordinates", "-spc")
         print_scaled_pbc_coordinates = .TRUE.
         CYCLE arg_loop
      CASE ("-stride")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)
         READ (UNIT=arg, FMT=*, ERR=104) stride
         IF (stride < 1) THEN
            CALL abort_program(routineN, "Invalid stride for frame dump specified: stride must be greater than zero")
         END IF
         CYCLE arg_loop
104      CALL abort_program(routineN, "Invalid stride for frame dump specified "// &
                            "(an integer number greater than 0 is expected)")
      CASE ("-trace_atoms")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)
         READ (UNIT=arg, FMT=*, ERR=105) eps_out_of_box
         IF (eps_out_of_box <= 0.0_dp) THEN
            CALL abort_program(routineN, "Invalid threshold value for -trace_atoms flag specified")
         END IF
         trace_atoms = .TRUE.
         CYCLE arg_loop
105      CALL abort_program(routineN, "Invalid threshold value for -trace_atoms flag specified")
      CASE DEFAULT
         IF (arg(1:1) == "-") THEN
            CALL print_help()
            CALL abort_program(routineN, "Unknown command line flag """//TRIM(arg)//""" found")
         END IF
         xyz_file_name = arg
         EXIT arg_loop
      END SELECT
   END DO arg_loop

   ! Check flag compatibility
   IF (.NOT. have_cell_info) THEN
      CALL abort_program(routineN, "No cell information available. Neither -abc, -cell, nor -cell_file flag found")
   END IF
   IF (first_frame > last_frame) THEN
      CALL abort_program(routineN, "Number of first frame greater than number of last frame")
   END IF
   IF (.NOT. apply_pbc .AND. trace_atoms) THEN
      CALL abort_program(routineN, "The -trace_atoms flag requires the specification of a -pbc flag")
   END IF
   IF (print_scaled_coordinates .AND. print_scaled_pbc_coordinates) THEN
      CALL abort_program(routineN, "The -sc flag and the -spc flag are incompatible")
   END IF
   IF (.NOT. apply_pbc .AND. print_scaled_coordinates) THEN
      CALL abort_program(routineN, "The -sc flag requires the specification of a -pbc flag")
   END IF
   IF (.NOT. apply_pbc .AND. print_scaled_pbc_coordinates) THEN
      CALL abort_program(routineN, "The -spc flag requires the specification of a -pbc flag")
   END IF

   ! Open cell input file (if specified)
   IF (have_cell_file) THEN
      INQUIRE (FILE=cell_file_name, EXIST=exists)
      IF (.NOT. exists) CALL abort_program(routineN, "The specified cell file <"// &
                                           TRIM(cell_file_name)//"> does not exist")
      OPEN (UNIT=cell_file_unit, &
            FILE=cell_file_name, &
            STATUS="OLD", &
            ACCESS="SEQUENTIAL", &
            FORM="FORMATTED", &
            POSITION="REWIND", &
            ACTION="READ", &
            IOSTAT=istat)
      IF (istat /= 0) CALL abort_program(routineN, "The cell file <"// &
                                         TRIM(cell_file_name)//"> could not be opened")
      IF (info) WRITE (UNIT=error_unit, FMT="(A)") "# Reading cell file : "//TRIM(cell_file_name)
   END IF

   ! Open XYZ input file
   INQUIRE (FILE=xyz_file_name, EXIST=exists)
   IF (.NOT. exists) CALL abort_program(routineN, "The specified XYZ file <"// &
                                        TRIM(xyz_file_name)//"> does not exist")
   OPEN (UNIT=xyz_file_unit, &
         FILE=xyz_file_name, &
         STATUS="OLD", &
         ACCESS="SEQUENTIAL", &
         FORM="FORMATTED", &
         POSITION="REWIND", &
         ACTION="READ", &
         IOSTAT=istat)
   IF (istat /= 0) CALL abort_program(routineN, "The XYZ file <"// &
                                      TRIM(xyz_file_name)//"> could not be opened")
   IF (info) WRITE (UNIT=error_unit, FMT="(A)") "# Reading XYZ file : "//TRIM(xyz_file_name)

   ! Read the first two lines of the XYZ file
   READ (UNIT=xyz_file_unit, FMT="(A)", IOSTAT=istat) arg
   IF (istat /= 0) THEN
      CALL abort_program(routineN, "Reading the first line of the current frame from the XYZ file <"// &
                         TRIM(xyz_file_name)//"> failed")
   END IF
   IF (arg(1:1) == "#") THEN
      READ (UNIT=arg, FMT=*) string, natom
   ELSE
      READ (UNIT=arg, FMT=*) natom
   END IF
   IF (istat /= 0) THEN
      CALL abort_program(routineN, "Reading the number of atoms from the XYZ file failed")
   END IF
   READ (UNIT=xyz_file_unit, FMT="(A)", IOSTAT=istat) remark_xyz
   IF (istat /= 0) CALL abort_program(routineN, "Reading the second line from the XYZ file <"// &
                                      TRIM(xyz_file_name)//"> failed")
   REWIND (UNIT=xyz_file_unit)

   ! Open DCD output file
   IF (LEN_TRIM(dcd_file_name) == 0) THEN
      i = LEN_TRIM(xyz_file_name)
      IF (xyz_file_name(i - 2:i) == "xyz") THEN
         dcd_file_name = xyz_file_name(1:i - 3)//"dcd"
      ELSE
         dcd_file_name = xyz_file_name(1:i)//".dcd"
      END IF
   END IF
   INQUIRE (FILE=dcd_file_name, EXIST=exists)
   IF (exists) CALL abort_program(routineN, "The DCD file "// &
                                  TRIM(dcd_file_name)//" exists already")
   OPEN (UNIT=dcd_file_unit, &
         FILE=dcd_file_name, &
         STATUS="UNKNOWN", &
         ACCESS="SEQUENTIAL", &
         FORM="UNFORMATTED", &
         ACTION="WRITE", &
         IOSTAT=istat)
   IF (istat /= 0) CALL abort_program(routineN, "The unformatted DCD output file  "// &
                                      TRIM(dcd_file_name)//" could not be opened")
   IF (info) WRITE (UNIT=error_unit, FMT="(A)") "# Writing DCD file : "//TRIM(dcd_file_name)

   istep = 0
   iskip = 0
   dt = 0.0_dp

   ! Write DCD file header
   WRITE (UNIT=dcd_file_unit) "CORD", 0, istep, iskip, 0, 0, 0, 0, 0, 0, REAL(dt, KIND=sp), &
      1, 0, 0, 0, 0, 0, 0, 0, 0, 24
   remark1 = "REMARK CORD"//" DCD file created by "//TRIM(version_info)
   remark2 = "REMARK "//TRIM(ADJUSTL(remark_xyz))
   WRITE (UNIT=dcd_file_unit) 2, remark1(1:80), remark2(1:80)
   WRITE (UNIT=dcd_file_unit) natom

   ! Allocate work arrays
   ALLOCATE (r(natom, 3), STAT=istat)
   IF (istat /= 0) CALL abort_program(routineN, "Allocation of the array r failed")
   r(:, :) = 0.0_dp

   ALLOCATE (atomic_label(natom), STAT=istat)
   IF (istat /= 0) CALL abort_program(routineN, "Allocation of the vector atomic_label failed")
   atomic_label(:) = ""

   ! Loop over all frames in the XYZ file
   frame_loop: DO

      nframe = nframe + 1

      IF (nframe < first_frame) THEN
         dump_frame = .FALSE.
      ELSE
         IF (MODULO(nframe - first_frame, stride) == 0) THEN
            dump_frame = .TRUE.
         ELSE
            dump_frame = .FALSE.
         END IF
      END IF

      ! Read unit cell information, if available
      IF (have_cell_file) THEN
         DO
            READ (UNIT=cell_file_unit, FMT=*, IOSTAT=istat) arg
            IF (istat < 0) EXIT frame_loop
            IF (istat /= 0) THEN
               CALL abort_program(routineN, "Reading line from cell file")
            END IF
            IF (arg(1:1) == "#") THEN
               CYCLE
            ELSE
               BACKSPACE (UNIT=cell_file_unit)
               EXIT
            END IF
         END DO
         READ (UNIT=cell_file_unit, FMT=*, IOSTAT=istat) istep, tstep, ((h(j, i), j=1, 3), i=1, 3)
         IF (istat /= 0) THEN
            CALL abort_program(routineN, "Reading information from cell file")
         END IF
      END IF

      ! Initialise or update cell info
      IF (have_cell_file .OR. (nframe == 1)) THEN
         a(1:3) = h(1:3, 1)
         b(1:3) = h(1:3, 2)
         c(1:3) = h(1:3, 3)
         abc(1) = SQRT(DOT_PRODUCT(a(1:3), a(1:3)))
         abc(2) = SQRT(DOT_PRODUCT(b(1:3), b(1:3)))
         abc(3) = SQRT(DOT_PRODUCT(c(1:3), c(1:3)))
         alpha = angle(b(1:3), c(1:3))*degree
         beta = angle(a(1:3), c(1:3))*degree
         gamma = angle(a(1:3), b(1:3))*degree
      END IF

      ! Read first line of the current frame in the XYZ file
      READ (UNIT=xyz_file_unit, FMT="(A)", IOSTAT=istat) arg
      IF (istat < 0) EXIT frame_loop
      IF (istat /= 0) THEN
         CALL abort_program(routineN, "Reading the first line of the current frame from the XYZ file failed")
      END IF
      IF (arg(1:1) == "#") THEN
         READ (UNIT=arg, FMT=*) string, natom
      ELSE
         READ (UNIT=arg, FMT=*) natom
      END IF
      IF (istat /= 0) THEN
         CALL abort_program(routineN, "Reading the number of atoms from the XYZ file failed")
      END IF
      IF (natom /= SIZE(r, 1)) THEN
         CALL abort_program(routineN, "Number of atoms changed for the current frame")
      END IF

      ! Read second line of the current frame in the XYZ file
      READ (UNIT=xyz_file_unit, FMT="(A)", IOSTAT=istat) remark_xyz
      IF (istat /= 0) CALL abort_program(routineN, "Reading the second line from the XYZ file failed")

      ! Optionally print some info
      IF (info .AND. dump_frame) THEN
         WRITE (UNIT=error_unit, FMT="(A,/,A,I0)") &
            "#", "# Frame number     : ", nframe
         WRITE (UNIT=error_unit, FMT="(A,/,(A,F12.6))") &
            "#", &
            "# a [Angstrom]     : ", abc(1), &
            "# b [Angstrom]     : ", abc(2), &
            "# c [Angstrom]     : ", abc(3), &
            "# alpha [degree]   : ", alpha, &
            "# beta  [degree]   : ", beta, &
            "# gamma [degree]   : ", gamma
      END IF

      IF (info) THEN
         WRITE (UNIT=output_unit, FMT="(T2,I0)") natom
         WRITE (UNIT=output_unit, FMT="(A)") TRIM(ADJUSTL(remark_xyz))
      END IF

      ! Read in the atomic positions of the current frame from the XYZ file
      iatom = 0
      DO
         READ (UNIT=xyz_file_unit, FMT=*) arg
         IF (arg(1:1) == "#") THEN
            CYCLE
         ELSE
            BACKSPACE (UNIT=xyz_file_unit)
         END IF
         iatom = iatom + 1
         READ (UNIT=xyz_file_unit, FMT=*, IOSTAT=istat) atomic_label(iatom), r(iatom, 1:3)
         IF (istat /= 0) THEN
            message = ""
            WRITE (UNIT=message, FMT="(A,I0,A,I0,A)") &
               "Reading line ", iatom + 2, " of the current frame from XYZ file (atom ", iatom, ") failed"
            CALL abort_program(routineN, TRIM(message))
         END IF
         CALL uppercase(atomic_label(iatom) (1:1))
         IF (LEN_TRIM(atomic_label(iatom)) > 1) CALL lowercase(atomic_label(iatom) (2:2))
         atomic_label(iatom) = TRIM(ADJUSTL(atomic_label(iatom)))
         IF (iatom == natom) EXIT
      END DO

      ! Store the atomic positions of the first frame in the current XYZ input file for
      ! the output of the atomic displacements
      IF ((nframe == 1) .AND. print_atomic_displacements) THEN
         IF (.NOT. ALLOCATED(r0)) THEN
            ALLOCATE (r0(natom, 3), STAT=istat)
            IF (istat /= 0) CALL abort_program(routineN, "Allocation of the array r0 failed")
         END IF
         r0(:, :) = r(:, :)
         IF (.NOT. ALLOCATED(atomic_displacement)) THEN
            ALLOCATE (atomic_displacement(natom), STAT=istat)
            IF (istat /= 0) THEN
               CALL abort_program(routineN, "Allocation of the vector atomic_displacement failed")
            END IF
         END IF
         atomic_displacement(:) = 0.0_dp
      END IF

      IF (dump_frame) THEN
         ! Apply periodic boundary conditions before dumping the coordinates if requested
         IF (apply_pbc) THEN
            IF (.NOT. ALLOCATED(r_pbc)) THEN
               ALLOCATE (r_pbc(natom, 3), STAT=istat)
               IF (istat /= 0) CALL abort_program(routineN, "Allocation of the array r_pbc failed")
            END IF
            IF (.NOT. ALLOCATED(s)) THEN
               ALLOCATE (s(natom, 3), STAT=istat)
               IF (istat /= 0) CALL abort_program(routineN, "Allocation of the array s failed")
            END IF
            IF (.NOT. ALLOCATED(s_pbc)) THEN
               ALLOCATE (s_pbc(natom, 3), STAT=istat)
               IF (istat /= 0) CALL abort_program(routineN, "Allocation of the array s_pbc failed")
            END IF
            CALL pbc(r, r_pbc, s, s_pbc, h, hinv, debug, info, pbc0)
            CALL write_out_of_box_atoms(atomic_label, r, s, eps_out_of_box, h)
            ! Overwrite input coordinate with the PBCed coordinates for printing
            r(:, :) = r_pbc(:, :)
         END IF ! apply_pbc
         ! Calculate the atomic displacements with respect to the first frame,
         ! i.e. set of atomic positions
         IF (print_atomic_displacements) THEN
            DO iatom = 1, natom
               atomic_displacement(iatom) = SQRT((r(iatom, 1) - r0(iatom, 1))**2 + &
                                                 (r(iatom, 2) - r0(iatom, 2))**2 + &
                                                 (r(iatom, 3) - r0(iatom, 3))**2)
            END DO
         END IF
         IF (info) THEN
            IF (print_scaled_coordinates) THEN
               DO iatom = 1, natom
                  WRITE (UNIT=output_unit, FMT="(A5,3(1X,F14.6))") ADJUSTL(atomic_label(iatom)), s(iatom, 1:3)
               END DO
            ELSE IF (print_scaled_pbc_coordinates) THEN
               DO iatom = 1, natom
                  WRITE (UNIT=output_unit, FMT="(A5,3(1X,F14.6))") ADJUSTL(atomic_label(iatom)), s_pbc(iatom, 1:3)
               END DO
            ELSE
               DO iatom = 1, natom
                  WRITE (UNIT=output_unit, FMT="(A5,3(1X,F14.6))") ADJUSTL(atomic_label(iatom)), r(iatom, 1:3)
               END DO
            END IF
         END IF
         ! Dump the cell information in DCD format
         WRITE (UNIT=dcd_file_unit) abc(1), gamma, abc(2), beta, alpha, abc(3)
         ! Dump the atomic coordinates in DCD format
         DO i = 1, 3
            WRITE (UNIT=dcd_file_unit) REAL(r(1:natom, i), KIND=sp)
         END DO
         nframe_read = nframe_read + 1
      END IF ! dump_frame

      ! Exit loop and stop processing, if the last (requested) frame was encountered
      IF (nframe >= last_frame) EXIT

   END DO frame_loop

   nframe = nframe - 1

   ! Close files
   IF (have_cell_file) CLOSE (UNIT=cell_file_unit)
   CLOSE (UNIT=dcd_file_unit)
   CLOSE (UNIT=xyz_file_unit)

   IF (info) THEN
      WRITE (UNIT=error_unit, FMT="(A,/,A,I0)") &
         "#", &
         "# Frames processed : ", nframe_read
      WRITE (UNIT=error_unit, FMT="(A)") &
         "#", &
         "# Normal termination of "//TRIM(version_info)
   END IF

   ! Cleanup
   IF (ALLOCATED(atomic_label)) DEALLOCATE (atomic_label)
   IF (ALLOCATED(atomic_displacement)) DEALLOCATE (atomic_displacement)
   IF (ALLOCATED(r)) DEALLOCATE (r)
   IF (ALLOCATED(r0)) DEALLOCATE (r0)
   IF (ALLOCATED(r_pbc)) DEALLOCATE (r_pbc)
   IF (ALLOCATED(s)) DEALLOCATE (s)
   IF (ALLOCATED(s_pbc)) DEALLOCATE (s_pbc)

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param routine ...
!> \param message ...
! **************************************************************************************************
   SUBROUTINE abort_program(routine, message)
      ! Abort the program after printing an error message to stderr

      CHARACTER(LEN=*), INTENT(IN)                       :: routine, message

      CHARACTER(LEN=2*default_string_length)             :: error_message

      error_message = "*** ERROR in "//TRIM(routine)//": "//TRIM(message)//" ***"
      WRITE (UNIT=default_error_unit, FMT="(/,A,/)") TRIM(error_message)
      STOP "*** ABNORMAL PROGRAM TERMINATION of xyz2dcd v1.0 ***"

   END SUBROUTINE abort_program

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \return ...
! **************************************************************************************************
   FUNCTION angle(a, b) RESULT(angle_ab)
      ! Calculation of the angle between the vectors a and b. The angle is returned in radians.

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: a, b
      REAL(KIND=dp)                                      :: angle_ab

      REAL(KIND=dp), PARAMETER                           :: eps_geo = 1.0E-6_dp

      REAL(KIND=dp)                                      :: length_of_a, length_of_b
      REAL(KIND=dp), DIMENSION(SIZE(a, 1))               :: a_norm, b_norm

      length_of_a = SQRT(DOT_PRODUCT(a, a))
      length_of_b = SQRT(DOT_PRODUCT(b, b))

      IF ((length_of_a > eps_geo) .AND. (length_of_b > eps_geo)) THEN
         a_norm(:) = a(:)/length_of_a
         b_norm(:) = b(:)/length_of_b
         angle_ab = ACOS(MIN(MAX(DOT_PRODUCT(a_norm, b_norm), -1.0_dp), 1.0_dp))
      ELSE
         angle_ab = 0.0_dp
      END IF

   END FUNCTION angle

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param c ...
!> \param alpha ...
!> \param beta ...
!> \param gamma ...
!> \param h ...
! **************************************************************************************************
   SUBROUTINE build_h_matrix(a, b, c, alpha, beta, gamma, h)
      ! Calculate the h matrix of a simulation cell given the cell edge lengths a, b,
      ! and c in Angstrom and the angles alpha (<b,c)), beta (<(a,c)), and gamma (<(a,b))
      ! in degree.

      REAL(KIND=dp), INTENT(IN)                          :: a, b, c, alpha, beta, gamma
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: h

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_h_matrix'

      REAL(KIND=dp)                                      :: cosa, cosb, cosg, sing

      cosa = COS(alpha/degree)
      IF (ABS(cosa) < EPSILON(0.0_dp)) cosa = 0.0_dp

      cosb = COS(beta/degree)
      IF (ABS(cosb) < EPSILON(0.0_dp)) cosb = 0.0_dp

      cosg = COS(gamma/degree)
      IF (ABS(cosg) < EPSILON(0.0_dp)) cosg = 0.0_dp

      sing = SIN(gamma/degree)
      IF (ABS(sing) < EPSILON(0.0_dp)) sing = 0.0_dp

      h(1, 1) = 1.0_dp
      h(2, 1) = 0.0_dp
      h(3, 1) = 0.0_dp

      h(1, 2) = cosg
      h(2, 2) = sing
      h(3, 2) = 0.0_dp

      h(1, 3) = cosb
      h(2, 3) = (cosa - cosg*cosb)/sing
      IF ((1.0_dp - h(1, 3)**2 - h(2, 3)**2) < 0.0_dp) THEN
         CALL abort_program(routineN, "Build of the h matrix failed, check cell information")
      END IF
      h(3, 3) = SQRT(1.0_dp - h(1, 3)**2 - h(2, 3)**2)

      h(:, 1) = a*h(:, 1)
      h(:, 2) = b*h(:, 2)
      h(:, 3) = c*h(:, 3)

   END SUBROUTINE build_h_matrix

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \return ...
! **************************************************************************************************
   FUNCTION det_3x3(a) RESULT(det_a)
      ! Returns the determinante of the 3x3 matrix a.

      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: a
      REAL(KIND=dp)                                      :: det_a

      det_a = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) + &
              a(1, 2)*(a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3)) + &
              a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))

   END FUNCTION det_3x3

! **************************************************************************************************
!> \brief ...
!> \param h ...
!> \param hinv ...
!> \param deth ...
! **************************************************************************************************
   SUBROUTINE invert_matrix_3x3(h, hinv, deth)
      ! Calculate the inverse hinv and the determinant deth of the 3x3 matrix h.

      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: h
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: hinv
      REAL(KIND=dp), INTENT(OUT)                         :: deth

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'invert_matrix_3x3'

      deth = det_3x3(h)

      ! Numerics
      deth = ABS(deth)
      IF (deth < 1.0E-10_dp) THEN
         CALL abort_program(routineN, "Invalid h matrix for cell found; det(h) < 1.0E-10")
      END IF

      hinv(1, 1) = (h(2, 2)*h(3, 3) - h(3, 2)*h(2, 3))/deth
      hinv(2, 1) = (h(2, 3)*h(3, 1) - h(3, 3)*h(2, 1))/deth
      hinv(3, 1) = (h(2, 1)*h(3, 2) - h(3, 1)*h(2, 2))/deth

      hinv(1, 2) = (h(1, 3)*h(3, 2) - h(3, 3)*h(1, 2))/deth
      hinv(2, 2) = (h(1, 1)*h(3, 3) - h(3, 1)*h(1, 3))/deth
      hinv(3, 2) = (h(1, 2)*h(3, 1) - h(3, 2)*h(1, 1))/deth

      hinv(1, 3) = (h(1, 2)*h(2, 3) - h(2, 2)*h(1, 3))/deth
      hinv(2, 3) = (h(1, 3)*h(2, 1) - h(2, 3)*h(1, 1))/deth
      hinv(3, 3) = (h(1, 1)*h(2, 2) - h(2, 1)*h(1, 2))/deth

   END SUBROUTINE invert_matrix_3x3

! **************************************************************************************************
!> \brief ...
!> \param string ...
! **************************************************************************************************
   SUBROUTINE lowercase(string)
      ! Convert all letters in a string to lowercase
      CHARACTER(LEN=*), INTENT(INOUT)                    :: string

      INTEGER                                            :: i, iascii

      DO i = 1, LEN_TRIM(string)
         iascii = ICHAR(string(i:i))
         IF ((iascii >= 65) .AND. (iascii <= 90)) THEN
            string(i:i) = CHAR(iascii + 32)
         END IF
      END DO

   END SUBROUTINE lowercase

! **************************************************************************************************
!> \brief ...
!> \param r ...
!> \param r_pbc ...
!> \param s ...
!> \param s_pbc ...
!> \param h ...
!> \param hinv ...
!> \param debug ...
!> \param info ...
!> \param pbc0 ...
! **************************************************************************************************
   SUBROUTINE pbc(r, r_pbc, s, s_pbc, h, hinv, debug, info, pbc0)
      ! Apply the periodic boundary conditions (PBC) to the Cartesian coordinate array
      ! r given the cell edge lengths a, b, and c in Angstrom and the angles alpha (<b,c)),
      ! beta (<(a,c)), and gamma (<(a,b)) in degree.
      ! On output r_pbc is updated with the "PBCed" input coordinates, s with the scaled
      ! input coordinates, and s_pbc with scaled "PBCed" coordinates.
      ! If pbc0 is true then fold to the range [-l/2,+l/2[ (origin at box centre) else fold
      ! to the range [0,l[ (origin at lower left box corner).

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: r
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: r_pbc, s, s_pbc
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: h
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: hinv
      LOGICAL, INTENT(IN)                                :: debug, info, pbc0

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'pbc'

      INTEGER                                            :: i, natom
      LOGICAL                                            :: orthorhombic
      REAL(KIND=dp)                                      :: deth

      natom = SIZE(r, 1)
      IF (SIZE(r, 2) /= 3) CALL abort_program(routineN, "Array dimension for r must be 3")

      orthorhombic = ((h(2, 1) == 0.0_dp) .AND. &
                      (h(3, 1) == 0.0_dp) .AND. &
                      (h(1, 2) == 0.0_dp) .AND. &
                      (h(3, 2) == 0.0_dp) .AND. &
                      (h(1, 3) == 0.0_dp) .AND. &
                      (h(2, 3) == 0.0_dp))

      ! Build inverse of h
      hinv(:, :) = 0.0_dp
      CALL invert_matrix_3x3(h, hinv, deth)

      IF (info) THEN
         WRITE (UNIT=error_unit, FMT="(A)") "#"
         IF (orthorhombic) THEN
            WRITE (UNIT=error_unit, FMT="(A)") "# Cell symmetry   : orthorhombic"
         ELSE
            WRITE (UNIT=error_unit, FMT="(A)") "# Cell symmetry   : non-orthorhombic"
         END IF
         IF (debug) THEN
            WRITE (UNIT=error_unit, FMT="(A)") "#"
            WRITE (UNIT=error_unit, FMT="(A,3F12.6,A)") "#          / ", h(1, :), " \"
            WRITE (UNIT=error_unit, FMT="(A,3F12.6,A)") "# h      = | ", h(2, :), " |"
            WRITE (UNIT=error_unit, FMT="(A,3F12.6,A)") "#          \ ", h(3, :), " /"
            WRITE (UNIT=error_unit, FMT="(A)") "#"
            WRITE (UNIT=error_unit, FMT="(A,3F12.6,A)") "#          / ", hinv(1, :), " \"
            WRITE (UNIT=error_unit, FMT="(A,3F12.6,A)") "# Inv(h) = | ", hinv(2, :), " |"
            WRITE (UNIT=error_unit, FMT="(A,3F12.6,A)") "#          \ ", hinv(3, :), " /"
            WRITE (UNIT=error_unit, FMT="(A)") "#"
            WRITE (UNIT=error_unit, FMT="(A,F0.6)") "# det(h) = ", deth
         END IF
      END IF

      ! Calculate scaled coordinates and wrap back all atoms, i.e. apply the PBC
      IF (orthorhombic) THEN
         ! Try to save some flops in the case of an orthorhombic box
         DO i = 1, 3
            s(:, i) = r(:, i)*hinv(i, i)
         END DO
      ELSE
         s(:, :) = MATMUL(r(:, :), TRANSPOSE(hinv(:, :)))
      END IF

      IF (pbc0) THEN
         s_pbc(:, :) = s(:, :) - ANINT(s(:, :))
      ELSE
         s_pbc(:, :) = s(:, :) - FLOOR(s(:, :))
      END IF

      IF (orthorhombic) THEN
         DO i = 1, 3
            r_pbc(:, i) = s_pbc(:, i)*h(i, i)
         END DO
      ELSE
         r_pbc(:, :) = MATMUL(s_pbc(:, :), TRANSPOSE(h(:, :)))
      END IF

   END SUBROUTINE pbc

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE print_help()
      ! Print the program flags for help

      WRITE (UNIT=*, FMT="(T2,A)") &
         "", &
         "Program flags for "//TRIM(version_info)//":", &
         "", &
         " -abc <3 reals>                 : Cell vector lengths in Angstrom", &
         " -cell <9 reals>                : Cell vectors: a(x) a(y) a(z) b(x) b(y) b(z) c(x) c(y) c(z) in Angstrom", &
         " -cell_file, -cf <file>         : Name of the cell file in CP2K format", &
         " -debug, -d                     : Print debug information", &
         " -eo                            : Write standard output and standard error to the same logical unit", &
         " -first_frame, -ff <int>        : Number of the first frame which is dumped", &
         " -help, -h                      : Print this information", &
         " -info, -i                      : Print additional information for each frame (see also -debug flag)", &
         " -last_frame, -lf <int>         : Number of the last frame which is dumped", &
         " -pbc                           : Apply the periodic boundary conditions (PBC) to each frame before it is dumped", &
         "                                  (origin at lower left)", &
         " -pbc0                          : Apply the periodic boundary conditions (PBC) to each frame before it is dumped", &
         "                                  (origin at box centre)", &
         " -scaled_coordinates, -sc       : Print the scaled coordinates", &
      " -scaled_pbc_coordinates, -spc  : Print the scaled coordinates after periodic boundary conditions (PBC) have been applied", &
         " -stride <int>                  : Stride for frame dump (allows to skip frames, e.g. by dumping each 10th frame)", &
        " -trace_atoms <real>            : Print the atoms which left the simulation box given a threshold value in scaled units", &
         ""

      WRITE (UNIT=*, FMT="(T2,A)") &
         "Usage examples:", &
         "", &
         " xyz2dcd <optional flags> <cell information: -abc <3 reals>, -cell <9 reals>, or -cell_file <file>> <XYZ file>", &
         "", &
         "Specific usage examples:", &
         "", &
         " xyz2dcd -abc 27.341 27.341 27.341 project-pos-1.xyz", &
         " xyz2dcd -cell 27.341 0.0 0.0 0.0 27.341 0.0 0.0 0.0 27.341 project-pos-1.xyz", &
         " xyz2dcd -cell_file project-1.cell project-pos-1.xyz", &
         ""

      WRITE (UNIT=*, FMT="(T2,A)") &
         "Notes:", &
         "", &
         " - The -info and the -debug flags provide a more detailed output which is especially handy for tracing problems", &
         " - The input coordinates and cell vectors should be in Angstrom", &
         ""

   END SUBROUTINE print_help

! **************************************************************************************************
!> \brief ...
!> \param string ...
! **************************************************************************************************
   SUBROUTINE uppercase(string)
      ! Convert all letters in a string to uppercase

      CHARACTER(LEN=*), INTENT(INOUT)                    :: string

      INTEGER                                            :: i, iascii

      DO i = 1, LEN_TRIM(string)
         iascii = ICHAR(string(i:i))
         IF ((iascii >= 97) .AND. (iascii <= 122)) THEN
            string(i:i) = CHAR(iascii - 32)
         END IF
      END DO

   END SUBROUTINE uppercase

! **************************************************************************************************
!> \brief ...
!> \param atomic_label ...
!> \param r ...
!> \param s ...
!> \param eps_out_of_box ...
!> \param h ...
! **************************************************************************************************
   SUBROUTINE write_out_of_box_atoms(atomic_label, r, s, eps_out_of_box, h)
      ! Print a list of all atoms which have left the simulation box

      CHARACTER(LEN=5), DIMENSION(:), INTENT(IN)         :: atomic_label
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: r, s
      REAL(KIND=dp), INTENT(IN)                          :: eps_out_of_box
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: h

      INTEGER                                            :: i, iatom, natom, ncount
      REAL(KIND=dp)                                      :: rl, s_max, s_min, sl
      REAL(KIND=dp), DIMENSION(3)                        :: dr, ds

      ! Quick return, if no action is requested
      IF (eps_out_of_box <= 0.0_dp) RETURN

      s_max = 1.0_dp + eps_out_of_box
      s_min = -eps_out_of_box
      natom = SIZE(s, 1)
      ncount = 0
      DO iatom = 1, natom
         IF (ANY(s(iatom, :) < s_min) .OR. &
             ANY(s(iatom, :) > s_max)) THEN
            ncount = ncount + 1
            IF (ncount == 1) THEN
               WRITE (UNIT=error_unit, FMT="(A)") &
                  "#", &
                  "# Atoms out of box:", &
                  "# Atom index label              x              y              z           |dr|           |ds|"
            END IF
            ds(:) = s(iatom, :)
            DO i = 1, 3
               IF (s(iatom, i) < 0.0_dp) ds(i) = 0.0_dp
               IF (s(iatom, i) >= 1.0_dp) ds(i) = 1.0_dp
            END DO
            ds(:) = s(iatom, :) - ds(:)
            sl = SQRT(ds(1)**2 + ds(2)**2 + ds(3)**2)
            dr(:) = MATMUL(h(:, :), ds(:))
            rl = SQRT(dr(1)**2 + dr(2)**2 + dr(3)**2)
            WRITE (UNIT=error_unit, FMT="(A,I10,1X,A5,5(1X,F14.6))") &
               "# ", iatom, ADJUSTR(atomic_label(iatom)), r(iatom, :), rl, sl
         END IF
      END DO
      WRITE (UNIT=error_unit, FMT="(A,I0,A)") "# ", ncount, " atom(s) out of box"

   END SUBROUTINE write_out_of_box_atoms

END PROGRAM xyz2dcd
