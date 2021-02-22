!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM dumpdcd

! Version: 3.1
! Author:  Matthias Krack (MK)
! History: - Creation (13.02.2012,MK)
!          - XYZ file option added (13.02.2012,MK)
!          - Flags added for first and last frame and -o for output (24.05.2012,MK)
!          - XMOL flag added (25.05.2012,MK)
!          - PBC flag added (04.06.2012,MK)
!          - Stride flag added (05.06.2012,MK)
!          - Tracing of atoms added to detect the atoms which left the box (06.06.2012,MK)
!          - Keep input coordinates for further processing steps (15.06.2012,MK)
!          - VEL to CORD (-vel2cord flag) hack added (25.06.2012,MK)
!          - Added -displacement (-disp) flag (26.06.2012,MK)
!          - Dump the atomic displacements (CORD file) or temperatures (VEL file as x-coordinates of a DCD file (28.06.2012,MK)
!          - FRC to CORD (-frc2cord flag) hack added (28.01.2016,MK)
!          - Option -eformat added (15.02.2016,MK)
!          - Fix box unit string for velocity dump (07.06.2020,MK)
!          - Dump cell information to comment line of XMOL file if requested and available (15.02.2021,MK)
!
! Note: For -ekin a XYZ file is required to obtain the atomic labels.
!       The -info and the -debug flags provide a more detailed output which is especially handy for tracing problems.
!       The output in DCD format is in binary format.
!       The input coordinates should be in Angstrom. Velocities and forces are expected to be in atomic units.

! Uncomment the following line if this module is available (e.g. with gfortran) and comment the corresponding variable declarations below
! USE ISO_FORTRAN_ENV, ONLY: error_unit,input_unit,output_unit

   IMPLICIT NONE

   ! Comment the following lines if the ISO_FORTRAN_ENV is used (see above)
   INTEGER, PARAMETER                                 :: default_error_unit = 0, &
                                                         default_input_unit = 5, &
                                                         default_output_unit = 6
   INTEGER                                            :: error_unit = default_error_unit, &
                                                         input_unit = default_input_unit, &
                                                         output_unit = default_output_unit
   ! End Comment

   ! Parameters
   CHARACTER(LEN=*), PARAMETER :: routineN = "dumpdcd", &
                                  version_info = routineN//" v3.1 (15.02.2021)"

   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 200), &
                         sp = SELECTED_REAL_KIND(6, 30)
   INTEGER, PARAMETER :: default_string_length = 240, &
                         cell_input_unit = 10, &
                         xyz_input_unit = 11

   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp
   REAL(KIND=dp), PARAMETER :: angstrom = 0.52917720859_dp, &            ! [a.u.] -> [Angstrom]
                               degree = 180.0_dp/pi, &                   ! [rad]  -> [degree]
                               kelvin = 315774.647902944_dp, &           ! [a.u.] -> [K]
                               massunit = 1822.88484264550_dp            ! [u]    -> [a.u.]

   ! Variables
   CHARACTER(LEN=4)                                   :: id_dcd
   CHARACTER(LEN=10)                                  :: unit_string
   CHARACTER(LEN=17)                                  :: fmt_string
   CHARACTER(LEN=default_string_length)               :: arg, cell_file_name, dcd_file_name, message, out_file_name, &
                                                         output_format, remark_xyz, string, xyz_file_name
   CHARACTER(LEN=5), DIMENSION(:), ALLOCATABLE        :: atomic_label
   CHARACTER(LEN=80), DIMENSION(2)                    :: remark_dcd
   INTEGER                                            :: first_frame, have_unit_cell, i, iarg, iatom, &
                                                         iframe, iskip, istat, istep, last_frame, narg, &
                                                         natom_dcd, natom_xyz, ndcd_file, nframe, &
                                                         nframe_read, nremark, stride
   LOGICAL                                            :: apply_pbc, debug, dump_frame, eformat, ekin, eo, &
                                                         have_atomic_labels, have_cell_file, ignore_warnings, info, &
                                                         opened, output_format_dcd, output_format_xmol, pbc0, &
                                                         print_atomic_displacements, print_scaled_coordinates, &
                                                         print_scaled_pbc_coordinates, trace_atoms, vel2cord
   REAL(KIND=sp)                                      :: dt_dcd
   REAL(KIND=dp)                                      :: a, a_dcd, alpha, alpha_dcd, b, b_dcd, beta, beta_dcd, c, c_dcd, &
                                                         cell_volume, dt, eps_angle, eps_geo, eps_out_of_box, &
                                                         first_step_time, gamma, gamma_dcd, md_time_step, ndt, step_time, &
                                                         tavg, tavg_frame
   INTEGER, DIMENSION(16)                             :: idum
   REAL(KIND=dp), DIMENSION(3)                        :: rdum
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: atomic_displacement, atomic_mass, atomic_temperature
   REAL(KIND=dp), DIMENSION(3, 3)                     :: h, hinv, hmat
   REAL(KIND=sp), DIMENSION(:, :), ALLOCATABLE        :: r
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: r_pbc, r0, s, s_pbc

   apply_pbc = .FALSE.
   pbc0 = .FALSE.
   debug = .FALSE.
   dump_frame = .TRUE.
   eformat = .FALSE.
   ekin = .FALSE.
   eo = .FALSE.
   ignore_warnings = .FALSE.
   info = .FALSE.
   trace_atoms = .FALSE.
   vel2cord = .FALSE.
   print_atomic_displacements = .FALSE.
   print_scaled_coordinates = .FALSE.
   print_scaled_pbc_coordinates = .FALSE.
   first_frame = 1
   last_frame = 1000000 ! Hard limit of 1 Mio frames in total
   stride = 1
   ndcd_file = 0
   nframe = 0
   nframe_read = 0
   nremark = 0
   have_unit_cell = 0
   have_atomic_labels = .FALSE.
   have_cell_file = .FALSE.
   idum(:) = 0
   rdum(:) = 0.0_dp
   id_dcd = ""
   cell_file_name = ""
   dcd_file_name = ""
   out_file_name = ""
   xyz_file_name = ""
   output_format = "default"
   output_format_dcd = .FALSE.
   output_format_xmol = .FALSE.
   remark_dcd(:) = ""
   remark_xyz = ""
   fmt_string = ""
   dt = 0.0_dp
   dt_dcd = 0.0_sp
   a = 0.0_dp
   b = 0.0_dp
   c = 0.0_dp
   alpha = 0.0_dp
   beta = 0.0_dp
   gamma = 0.0_dp
   eps_out_of_box = -HUGE(0.0_dp)
   first_step_time = 0.0_dp
   md_time_step = 0.0_dp
   step_time = 0.0_dp
   tavg = 0.0_dp
   tavg_frame = 0.0_dp

   narg = command_argument_count()

   IF (narg == 0) THEN
      CALL print_help()
      CALL abort_program(routineN, "No input file(s) specified")
   END IF

   iarg = 0

   dcd_file_loop: DO

      iarg = iarg + 1

      CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)

      SELECT CASE (arg)
      CASE ("-cell_file")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=cell_file_name, STATUS=istat)
         have_cell_file = .TRUE.
         CYCLE dcd_file_loop
      CASE ("-debug", "-d")
         debug = .TRUE.
         info = .TRUE.
         CYCLE dcd_file_loop
      CASE ("-displacements", "-disp")
         print_atomic_displacements = .TRUE.
         CYCLE dcd_file_loop
      CASE ("-eformat")
         eformat = .TRUE.
         CYCLE dcd_file_loop
      CASE ("-ekin")
         ekin = .TRUE.
         CYCLE dcd_file_loop
      CASE ("-eo")
         eo = .TRUE.
         CYCLE dcd_file_loop
      CASE ("-first_frame", "-first", "-ff")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)
         READ (UNIT=arg, FMT=*, ERR=100) first_frame
         IF (first_frame <= 0) THEN
            CALL abort_program(routineN, "Invalid number for first frame specified: "// &
                               "first_frame must be greater than zero")
         END IF
         CYCLE dcd_file_loop
100      CALL abort_program(routineN, "Invalid number for first frame specified "// &
                            "(an integer number greater than zero is expected)")
      CASE ("-help", "-h")
         CALL print_help()
         STOP
      CASE ("-ignore_warnings")
         ignore_warnings = .TRUE.
         CYCLE dcd_file_loop
      CASE ("-info", "-i")
         info = .TRUE.
         CYCLE dcd_file_loop
      CASE ("-last_frame", "-last", "-lf")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)
         READ (UNIT=arg, FMT=*, ERR=101) last_frame
         IF (last_frame <= 0) THEN
            CALL abort_program(routineN, "Invalid number for last frame specified: "// &
                               "last_frame must be greater than zero")
         END IF
         CYCLE dcd_file_loop
101      CALL abort_program(routineN, "Invalid number for last frame specified "// &
                            "(an integer number greater than zero is expected)")
      CASE ("-md_time_step")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)
         READ (UNIT=arg, FMT=*, ERR=102) md_time_step
         IF (md_time_step <= 0.0_dp) THEN
            CALL abort_program(routineN, "Invalid (negative) MD time step specified")
         END IF
         CYCLE dcd_file_loop
102      CALL abort_program(routineN, "Invalid MD time step specified")
      CASE ("-o", "-output")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=out_file_name, STATUS=istat)
         CYCLE dcd_file_loop
      CASE ("-output_format", "-of")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=output_format, STATUS=istat)
         CALL uppercase(output_format)
         SELECT CASE (output_format)
         CASE ("DCD")
            output_format_dcd = .TRUE.
            output_format_xmol = .FALSE.
         CASE ("XMOL")
            output_format_dcd = .FALSE.
            output_format_xmol = .TRUE.
         CASE DEFAULT
            CALL abort_program(routineN, "Invalid output format type specified")
         END SELECT
         CYCLE dcd_file_loop
      CASE ("-pbc")
         apply_pbc = .TRUE.
         pbc0 = .FALSE.
         CYCLE dcd_file_loop
      CASE ("-pbc0")
         apply_pbc = .TRUE.
         pbc0 = .TRUE.
         CYCLE dcd_file_loop
      CASE ("-scaled_coordinates", "-sc")
         print_scaled_coordinates = .TRUE.
         CYCLE dcd_file_loop
      CASE ("-scaled_pbc_coordinates", "-spc")
         print_scaled_pbc_coordinates = .TRUE.
         CYCLE dcd_file_loop
      CASE ("-stride")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)
         READ (UNIT=arg, FMT=*, ERR=104) stride
         IF (stride < 1) THEN
            CALL abort_program(routineN, "Invalid stride for frame dump specified: stride must be greater than zero")
         END IF
         CYCLE dcd_file_loop
104      CALL abort_program(routineN, "Invalid stride for frame dump specified "// &
                            "(an integer number greater than 0 is expected)")
      CASE ("-trace_atoms")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=arg, STATUS=istat)
         READ (UNIT=arg, FMT=*, ERR=108) eps_out_of_box
         IF (eps_out_of_box <= 0.0_dp) THEN
            CALL abort_program(routineN, "Invalid threshold value for -trace_atoms flag specified")
         END IF
         trace_atoms = .TRUE.
         CYCLE dcd_file_loop
108      CALL abort_program(routineN, "Invalid threshold value for -trace_atoms flag specified")
      CASE ("-vel2cord", "-v2c", "-frc2cord", "-f2c")
         vel2cord = .TRUE.
         CYCLE dcd_file_loop
      CASE ("-xyz", "xyz_file")
         iarg = iarg + 1
         CALL get_command_argument(NUMBER=iarg, VALUE=xyz_file_name, STATUS=istat)
         have_atomic_labels = .TRUE.
         CYCLE dcd_file_loop
      CASE DEFAULT
         IF (arg(1:1) == "-") THEN
            CALL print_help()
            CALL abort_program(routineN, "Unknown command line flag """//TRIM(arg)//""" found")
         END IF
         dcd_file_name = arg
      END SELECT

      ! Check flag compatibility
      IF (first_frame > last_frame) THEN
         CALL abort_program(routineN, "Number of first frame greater than number of last frame")
      END IF
      IF ((.NOT. have_atomic_labels) .AND. output_format_xmol) THEN
         CALL abort_program(routineN, "The output format XMOL requires a valid xyz file (-xyz flag)")
      END IF
      IF (output_format_xmol .AND. ekin) THEN
         CALL abort_program(routineN, "Output format XMOL and the -ekin flag are incompatible")
      END IF
      IF (output_format_xmol .AND. print_atomic_displacements) THEN
         CALL abort_program(routineN, "Output format XMOL and the -displacements flag are incompatible")
      END IF
      IF (output_format_xmol .AND. eo) THEN
         CALL abort_program(routineN, "The -eo flag is incompatible with the output format XMOL")
      END IF
      IF (.NOT. have_atomic_labels .AND. ekin) THEN
         CALL abort_program(routineN, "ekin flag requires also the specification of a valid xyz file (-xyz flag)")
      END IF
      IF (.NOT. apply_pbc .AND. trace_atoms) THEN
         CALL abort_program(routineN, "The -trace_atoms flag requires the specification of a -pbc flag")
      END IF
      IF (ekin .AND. print_atomic_displacements) THEN
         CALL abort_program(routineN, "The -ekin flag and the -displacements flag are incompatible")
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

      ! Use optionally an ES format
      IF (eformat) THEN
         fmt_string = "(A5,3(1X,ES14.6))"
      ELSE
         fmt_string = "(A5,3(1X, F14.6))"
      END IF

      ! Open output units as requested
      IF (output_format_dcd) THEN
         ! Set default output file name if no file name was specified
         IF (LEN_TRIM(out_file_name) == 0) out_file_name = "output.dcd"
         ! Check if a new output file name was specified, if yes then close the old unit
         INQUIRE (UNIT=output_unit, NAME=string)
         IF (TRIM(string) /= TRIM(out_file_name)) CLOSE (UNIT=output_unit)
         INQUIRE (UNIT=output_unit, OPENED=opened)
         IF (.NOT. opened) THEN
            OPEN (UNIT=output_unit, &
                  FILE=out_file_name, &
                  STATUS="UNKNOWN", &
                  ACCESS="SEQUENTIAL", &
                  FORM="UNFORMATTED", &
                  ACTION="WRITE", &
                  IOSTAT=istat)
            IF (istat /= 0) CALL abort_program(routineN, "The unformatted output file could not be opened")
         END IF
      ELSE
         IF (eo) error_unit = output_unit
         IF (LEN_TRIM(out_file_name) > 0) THEN
            ! Check if a new output file name was specified, if yes then close the old unit
            INQUIRE (UNIT=output_unit, NAME=string)
            IF (TRIM(string) /= TRIM(out_file_name)) CLOSE (UNIT=output_unit)
            INQUIRE (UNIT=output_unit, OPENED=opened)
            IF (.NOT. opened) THEN
               OPEN (UNIT=output_unit, &
                     FILE=out_file_name, &
                     STATUS="UNKNOWN", &
                     ACCESS="SEQUENTIAL", &
                     FORM="FORMATTED", &
                     POSITION="REWIND", &
                     ACTION="WRITE", &
                     IOSTAT=istat)
               IF (istat /= 0) CALL abort_program(routineN, "The formatted output file could not be opened")
            END IF
         END IF
      END IF

      ! Avoid reading from and writing to the same file
      IF (TRIM(dcd_file_name) == TRIM(out_file_name)) THEN
         CALL abort_program(routineN, "Input and output file name cannot be the same")
      END IF

      ! Read reference xyz input file if requested
      IF (have_atomic_labels) THEN
         string = ""
         ! Check if a new XYZ input file name was specified, if yes then close the old unit
         INQUIRE (UNIT=xyz_input_unit, NAME=string)
         IF (TRIM(string) /= TRIM(xyz_file_name)) CLOSE (UNIT=xyz_input_unit)
         INQUIRE (UNIT=xyz_input_unit, OPENED=opened)
         ! Read new XYZ input file if needed and update reference information
         IF (.NOT. opened) THEN
            OPEN (UNIT=xyz_input_unit, &
                  FILE=xyz_file_name, &
                  STATUS="OLD", &
                  ACCESS="SEQUENTIAL", &
                  FORM="FORMATTED", &
                  POSITION="REWIND", &
                  ACTION="READ", &
                  IOSTAT=istat)
            IF (istat /= 0) CALL abort_program(routineN, "The XYZ file could not be opened")
            IF (info) WRITE (UNIT=error_unit, FMT="(A)") "#", "# Reading XYZ file: "//TRIM(xyz_file_name)
            READ (UNIT=xyz_input_unit, FMT="(A)", IOSTAT=istat) arg
            IF (istat /= 0) THEN
               CALL abort_program(routineN, "Reading line 1 of the XYZ file (number of atoms) failed")
            END IF
            IF (arg(1:1) == "#") THEN
               READ (UNIT=arg, FMT=*) string, natom_xyz
            ELSE
               READ (UNIT=arg, FMT=*) natom_xyz
            END IF
            IF (istat /= 0) THEN
               CALL abort_program(routineN, "Reading line 1 of the XYZ file (number of atoms) failed")
            END IF
            IF (ALLOCATED(atomic_label)) DEALLOCATE (atomic_label)
            ALLOCATE (atomic_label(natom_xyz), STAT=istat)
            IF (istat /= 0) CALL abort_program(routineN, "Allocation of the vector atomic_label failed")
            atomic_label(:) = ""
            IF (ekin) THEN
               IF (ALLOCATED(atomic_mass)) DEALLOCATE (atomic_mass)
               ALLOCATE (atomic_mass(natom_xyz), STAT=istat)
               IF (istat /= 0) CALL abort_program(routineN, "Allocation of the vector atomic_mass failed")
            END IF
            READ (UNIT=xyz_input_unit, FMT="(A)", IOSTAT=istat) remark_xyz
            IF (istat /= 0) CALL abort_program(routineN, "Reading line 2 of the XYZ file (remark line) failed")
            iatom = 0
            DO
               READ (UNIT=xyz_input_unit, FMT=*) arg
               IF (arg(1:1) == "#") THEN
                  CYCLE
               ELSE
                  BACKSPACE (UNIT=xyz_input_unit)
               END IF
               iatom = iatom + 1
               READ (UNIT=xyz_input_unit, FMT=*, IOSTAT=istat) atomic_label(iatom), rdum(1:3)
               IF (istat /= 0) THEN
                  message = ""
                  WRITE (UNIT=message, FMT="(A,I0,A,I0,A)") &
                     "Reading line ", iatom + 2, " of the XYZ file (atom ", iatom, ") failed"
                  CALL abort_program(routineN, TRIM(message))
               END IF
               CALL uppercase(atomic_label(iatom))
               IF (LEN_TRIM(atomic_label(iatom)) > 1) CALL lowercase(atomic_label(iatom) (2:2))
               atomic_label(iatom) = TRIM(ADJUSTL(atomic_label(iatom)))
               IF (ekin) atomic_mass(iatom) = get_atomic_mass(atomic_label(iatom))
               IF (iatom == natom_xyz) EXIT
            END DO
            IF (info) THEN
               WRITE (UNIT=error_unit, FMT="(A,I0)") &
                  "# Number of atoms : ", natom_xyz, &
                  "# Remark          : "//TRIM(ADJUSTL(remark_xyz))
            END IF
         END IF
      END IF

      ! Read cell information file if requested
      IF (have_cell_file) THEN
         string = ""
         ! Check if a new cell file name was specified, if yes then close the old unit
         INQUIRE (UNIT=cell_input_unit, NAME=string)
         IF (TRIM(string) /= TRIM(cell_file_name)) CLOSE (UNIT=cell_input_unit)
         INQUIRE (UNIT=cell_input_unit, OPENED=opened)
         ! Read new cell file if needed and update reference information
         IF (.NOT. opened) THEN
            OPEN (UNIT=cell_input_unit, &
                  FILE=cell_file_name, &
                  STATUS="OLD", &
                  ACCESS="SEQUENTIAL", &
                  FORM="FORMATTED", &
                  POSITION="REWIND", &
                  ACTION="READ", &
                  IOSTAT=istat)
            IF (istat /= 0) CALL abort_program(routineN, "The cell file could not be opened")
            IF (info) WRITE (UNIT=error_unit, FMT="(A)") "#", "# Reading cell file: "//TRIM(cell_file_name)
         END IF
      END IF

      ! Increment DCD file counter
      ndcd_file = ndcd_file + 1

      ! Open input file in DCD format
      OPEN (UNIT=input_unit, &
            FILE=dcd_file_name, &
            STATUS="OLD", &
            FORM="UNFORMATTED", &
            POSITION="REWIND", &
            ACTION="READ", &
            IOSTAT=istat)
      IF (istat /= 0) CALL abort_program(routineN, "The DCD file could not be opened")
      IF (info) THEN
         WRITE (UNIT=error_unit, FMT="(A,/,(A,I0))") &
            "#", &
            "# DCD file number : ", ndcd_file, &
            "# Reading DCD file: "//TRIM(dcd_file_name)
      END IF

      ! Read 1st record of DCD file
      READ (UNIT=input_unit) id_dcd, idum(1), istep, iskip, idum(2:7), dt_dcd, have_unit_cell, idum(8:16)
      IF (info) THEN
         WRITE (UNIT=error_unit, FMT="(A,2(/,A,I0),/,A,F9.3,/,A,I0)") &
            "# DCD id string   : "//id_dcd, &
            "# Step            : ", istep, &
            "# Print frequency : ", iskip, &
            "# Time step [fs]  : ", dt_dcd, &
            "# Unit cell       : ", have_unit_cell
      END IF

      IF (ekin .AND. TRIM(ADJUSTL(id_dcd)) /= "VEL") THEN
         CALL abort_program(routineN, "ekin flag requires a DCD file with VELocities")
      END IF
      IF (apply_pbc .AND. (have_unit_cell /= 1)) THEN
         CALL abort_program(routineN, "pbc flags require that unit cell information is available")
      END IF

      IF ((TRIM(ADJUSTL(id_dcd)) == "FRC") .OR. (TRIM(ADJUSTL(id_dcd)) == "VEL")) THEN
         unit_string = "[a.u.]"
         IF (apply_pbc) CALL abort_program(routineN, "pbc flags require a DCD file with COoRDinates")
      ELSE IF (TRIM(ADJUSTL(id_dcd)) == "CORD") THEN
         unit_string = "[Angstrom]"
      ELSE
         CALL abort_program(routineN, "Unknown DCD id found (use -debug or -info flag for details)")
      END IF

      ! Read 2nd record of DCD file
      READ (UNIT=input_unit) nremark, remark_dcd(1), remark_dcd(2)
      IF (info) THEN
         DO i = 1, nremark
            WRITE (UNIT=error_unit, FMT="(A,I1,A)") &
               "# Remark ", i, "        : "//TRIM(remark_dcd(i))
         END DO
      END IF

      ! Read 3rd record of DCD file
      READ (UNIT=input_unit) natom_dcd
      IF (info) THEN
         WRITE (UNIT=error_unit, FMT="(A,I0)") &
            "# Number of atoms : ", natom_dcd
      END IF

      IF (have_atomic_labels) THEN
         IF (natom_dcd /= natom_xyz) CALL abort_program(routineN, "Number of atoms in XYZ and DCD file differ")
      ELSE
         IF (.NOT. ALLOCATED(atomic_label)) THEN
            ALLOCATE (atomic_label(natom_dcd), STAT=istat)
            IF (istat /= 0) CALL abort_program(routineN, "Allocation of the vector atomic_label failed")
            atomic_label(:) = ""
         END IF
      END IF

      ! Dump output in DCD format => dcd2dcd functionality
      IF (output_format_dcd) THEN
         IF (vel2cord) id_dcd = "CORD"
         ! Overwrite MD time step, if supplied
         IF (md_time_step > 0.0_dp) dt_dcd = REAL(md_time_step, KIND=sp)
         ! Note, dt_dcd is REAL*4 and the rest is INTEGER*4
         WRITE (UNIT=output_unit) id_dcd, idum(1), istep, iskip, idum(2:7), dt_dcd, have_unit_cell, idum(8:16)
         WRITE (UNIT=output_unit) nremark, remark_dcd(1), remark_dcd(2)
         WRITE (UNIT=output_unit) natom_dcd
      END IF

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
         IF (have_unit_cell == 1) THEN
            READ (UNIT=input_unit, IOSTAT=istat) a, gamma, b, beta, alpha, c
            IF (istat < 0) EXIT frame_loop
         END IF
         IF (have_cell_file) THEN
            ! Read cell information from cell file if available
            DO
               READ (UNIT=cell_input_unit, FMT=*, IOSTAT=istat) string
               IF (istat /= 0) THEN
                  message = ""
                  WRITE (UNIT=message, FMT="(A)") &
                     "Reading frame (step) from cell file "//TRIM(cell_file_name)//" failed"
                  CALL abort_program(routineN, TRIM(message))
               END IF
               IF (string(1:1) == "#") THEN
                  CYCLE
               ELSE
                  BACKSPACE (UNIT=cell_input_unit)
                  EXIT
               END IF
            END DO
            READ (UNIT=cell_input_unit, FMT=*, IOSTAT=istat) iframe, step_time, hmat(1:3, 1:3), cell_volume
            IF (istat /= 0) THEN
               WRITE (UNIT=error_unit, FMT="(/,T2,A,I0)") &
                  "IOSTAT =   ", istat, &
                  "Line read: """//TRIM(string)//""""
               message = ""
               WRITE (UNIT=message, FMT="(A)") &
                  "Invalid cell information read from cell file "//TRIM(cell_file_name)
               CALL abort_program(routineN, TRIM(message))
            END IF
            ! Store step time of the first (selected) MD step
            IF (nframe == first_frame) first_step_time = step_time
            ! Perform checks only for the selected frames
            IF ((nframe >= first_frame) .AND. (nframe <= last_frame)) THEN
               ! Save cell information from DCD header, if available
               IF (have_unit_cell == 1) THEN
                  a_dcd = a
                  b_dcd = b
                  c_dcd = c
                  alpha_dcd = alpha
                  beta_dcd = beta
                  gamma_dcd = gamma
               END IF
               ! Overwrite cell information from DCD header
               a = SQRT(DOT_PRODUCT(hmat(1:3, 1), hmat(1:3, 1)))
               b = SQRT(DOT_PRODUCT(hmat(1:3, 2), hmat(1:3, 2)))
               c = SQRT(DOT_PRODUCT(hmat(1:3, 3), hmat(1:3, 3)))
               ! Lattice vectors from DCD headers of VEL files are atomic units, but the CP2K cell file is in Angstrom
               IF (unit_string == "[a.u.]") THEN
                  a = a/angstrom
                  b = b/angstrom
                  c = c/angstrom
               END IF
               alpha = angle(hmat(1:3, 2), hmat(1:3, 3))*degree
               beta = angle(hmat(1:3, 3), hmat(1:3, 1))*degree
               gamma = angle(hmat(1:3, 1), hmat(1:3, 2))*degree
               IF (have_unit_cell == 1) THEN
                  ! Check consistency of DCD and cell file information
                  IF (.NOT. ignore_warnings) THEN
                     IF (md_time_step > 0.0_dp) THEN
                        dt = md_time_step
                        string = "Time step (requested)   ="
                     ELSE
                        dt = REAL(dt_dcd, KIND=dp)
                        string = "Time step (DCD header)  ="
                     END IF
                     ndt = (step_time - first_step_time)/dt
                     IF (ABS(ndt - ANINT(ndt)) > 1.0E-8_dp) THEN
                        WRITE (UNIT=error_unit, FMT="(/,T2,A,I8,/,(T2,A,F15.6,A))") &
                           "Step number (CELL file) = ", iframe, &
                           "Step time (CELL file)   = ", step_time, " fs", &
                           TRIM(string)//" ", dt, " fs"
                        WRITE (UNIT=error_unit, FMT="(/,T2,A)") &
                           "*** WARNING: MD step time in cell file is not a multiple of the MD time step in the DCD file header ***"
                     END IF
                  END IF
                  eps_geo = 1.0E-7_dp
                  IF (ABS(a - a_dcd) > eps_geo) THEN
                     WRITE (UNIT=error_unit, FMT="(/,(T2,A,F14.6))") &
                        "a (CELL file)  =   ", a, &
                        "a (DCD header) =   ", a_dcd
                     CALL abort_program(routineN, "Cell file and DCD file information for lattice constant ""a"" differ")
                  END IF
                  IF (ABS(b - b_dcd) > eps_geo) THEN
                     WRITE (UNIT=error_unit, FMT="(/,(T2,A,F14.6))") &
                        "b (CELL file)  =   ", b, &
                        "b (DCD header) =   ", b_dcd
                     CALL abort_program(routineN, "Cell file and DCD file information for lattice constant ""b"" differ")
                  END IF
                  IF (ABS(c - c_dcd) > eps_geo) THEN
                     WRITE (UNIT=error_unit, FMT="(/,(T2,A,F14.6))") &
                        "c (CELL file)  =   ", c, &
                        "c (DCD header) =   ", c_dcd
                     CALL abort_program(routineN, "Cell file and DCD file information for lattice constant ""c"" differ")
                  END IF
                  eps_angle = 1.0E-4_dp
                  IF (ABS(alpha - alpha_dcd) > eps_angle) THEN
                     WRITE (UNIT=error_unit, FMT="(/,(T2,A,F14.6))") &
                        "alpha (CELL file)  =   ", alpha, &
                        "alpha (DCD header) =   ", alpha_dcd
                     CALL abort_program(routineN, "Cell file and DCD file information for cell angle ""alpha"" differ")
                  END IF
                  IF (ABS(beta - beta_dcd) > eps_angle) THEN
                     WRITE (UNIT=error_unit, FMT="(/,(T2,A,F14.6))") &
                        "beta (CELL file)  =   ", beta, &
                        "beta (DCD header) =   ", beta_dcd
                     CALL abort_program(routineN, "Cell file and DCD file information for cell angle ""beta"" differ")
                  END IF
                  IF (ABS(gamma - gamma_dcd) > eps_angle) THEN
                     WRITE (UNIT=error_unit, FMT="(/,(T2,A,F14.6))") &
                        "gamma (CELL file)  =   ", gamma, &
                        "gamma (DCD header) =   ", gamma_dcd
                     CALL abort_program(routineN, "Cell file and DCD file information for cell angle ""beta"" differ")
                  END IF
               END IF
            END IF
         END IF

         IF ((info .OR. trace_atoms) .AND. dump_frame) THEN
            WRITE (UNIT=error_unit, FMT="(A,/,A,I0)") &
               "#", "# Frame number    : ", nframe
         END IF

         ! Print unit cell information, if available
         IF (info .AND. dump_frame) THEN
            IF (have_unit_cell == 1) THEN
               WRITE (UNIT=error_unit, FMT="(A,/,(A,T19,A,F12.6))") &
                  "#", &
                  "# a "//TRIM(unit_string), ": ", a, &
                  "# b "//TRIM(unit_string), ": ", b, &
                  "# c "//TRIM(unit_string), ": ", c
               WRITE (UNIT=error_unit, FMT="(A,F12.6)") &
                  "# alpha [degree]  : ", alpha, &
                  "# beta  [degree]  : ", beta, &
                  "# gamma [degree]  : ", gamma
            END IF
         END IF

         ! Allocate the array for the current atomic positions if needed
         IF (.NOT. ALLOCATED(r)) THEN
            ALLOCATE (r(natom_dcd, 3), STAT=istat)
            IF (istat /= 0) CALL abort_program(routineN, "Allocation of the array r failed")
         END IF

         ! Read in the atomic positions of the current frame
         READ (UNIT=input_unit, IOSTAT=istat) r(1:natom_dcd, 1)
         IF (istat < 0) EXIT frame_loop
         READ (UNIT=input_unit) r(1:natom_dcd, 2)
         READ (UNIT=input_unit) r(1:natom_dcd, 3)

         ! Store the atomic positions of the first frame in the current DCD input file for
         ! the output of the atomic displacements
         IF ((nframe == 1) .AND. print_atomic_displacements) THEN
            IF (.NOT. ALLOCATED(r0)) THEN
               ALLOCATE (r0(natom_dcd, 3), STAT=istat)
               IF (istat /= 0) CALL abort_program(routineN, "Allocation of the array r0 failed")
            END IF
            r0(:, :) = r(:, :)
            IF (.NOT. ALLOCATED(atomic_displacement)) THEN
               ALLOCATE (atomic_displacement(natom_dcd), STAT=istat)
               IF (istat /= 0) THEN
                  CALL abort_program(routineN, "Allocation of the vector atomic_displacement failed")
               END IF
            END IF
            atomic_displacement(:) = 0.0_dp
         END IF

         IF (ekin .AND. TRIM(ADJUSTL(id_dcd)) == "VEL") THEN

            ! Dump the kinetic energy of each as an "atomic temperature"
            IF (dump_frame) THEN

               IF (.NOT. ALLOCATED(atomic_temperature)) THEN
                  ALLOCATE (atomic_temperature(natom_dcd), STAT=istat)
                  IF (istat /= 0) THEN
                     CALL abort_program(routineN, "Allocation of the vector atomic_temperature failed")
                  END IF
                  atomic_temperature(:) = 0.0_dp
               END IF

               IF (info) THEN
                  WRITE (UNIT=error_unit, FMT="(A)") &
                     "#", &
                     "# Temperature [K]"
               END IF

               tavg_frame = 0.0_dp
               DO iatom = 1, natom_dcd
                  atomic_temperature(iatom) = atomic_mass(iatom)*(r(iatom, 1)*r(iatom, 1) + &
                                                                  r(iatom, 2)*r(iatom, 2) + &
                                                                  r(iatom, 3)*r(iatom, 3))*kelvin/3.0_dp
                  tavg_frame = tavg_frame + atomic_temperature(iatom)
               END DO
               tavg_frame = tavg_frame/REAL(natom_dcd, KIND=dp)

               IF (output_format_dcd) THEN
                  IF (have_unit_cell == 1) THEN
                     WRITE (UNIT=output_unit) a, gamma, b, beta, alpha, c
                  END IF
                  ! This is a hack for VMD: dump the atomic temperatures as the x-coordinates
                  ! of a DCD VEL file
                  WRITE (UNIT=output_unit) REAL(atomic_temperature(:), KIND=sp)
                  ! The y and z coordinates are filled with zeros
                  atomic_temperature(:) = 0.0_dp
                  WRITE (UNIT=output_unit) REAL(atomic_temperature(:), KIND=sp)
                  WRITE (UNIT=output_unit) REAL(atomic_temperature(:), KIND=sp)
               ELSE
                  DO iatom = 1, natom_dcd
                     WRITE (UNIT=output_unit, FMT="(A5,5X,F25.3)") ADJUSTL(atomic_label(iatom)), atomic_temperature
                  END DO
               END IF

               IF (info) THEN
                  WRITE (UNIT=error_unit, FMT="(A,F12.3)") &
                     "# T [K] this frame: ", tavg_frame
               END IF

               tavg = tavg + tavg_frame

            END IF ! dump_frame

         ELSE

            IF (dump_frame) THEN

               ! Apply periodic boundary conditions before dumping the coordinates if requested
               IF (apply_pbc) THEN
                  IF (.NOT. ALLOCATED(r_pbc)) THEN
                     ALLOCATE (r_pbc(natom_dcd, 3), STAT=istat)
                     IF (istat /= 0) CALL abort_program(routineN, "Allocation of the array r_pbc failed")
                  END IF
                  IF (.NOT. ALLOCATED(s)) THEN
                     ALLOCATE (s(natom_dcd, 3), STAT=istat)
                     IF (istat /= 0) CALL abort_program(routineN, "Allocation of the array s failed")
                  END IF
                  IF (.NOT. ALLOCATED(s_pbc)) THEN
                     ALLOCATE (s_pbc(natom_dcd, 3), STAT=istat)
                     IF (istat /= 0) CALL abort_program(routineN, "Allocation of the array s_pbc failed")
                  END IF
                  CALL pbc(r, r_pbc, s, s_pbc, a, b, c, alpha, beta, gamma, debug, info, pbc0, h, hinv)
                  CALL write_out_of_box_atoms(atomic_label, r, s, eps_out_of_box, h)
                  ! Overwrite input coordinate with the PBCed coordinates for printing
                  r(:, :) = REAL(r_pbc(:, :), KIND=sp)
               END IF ! apply_pbc

               ! Calculate the atomic displacements with respect to the first frame,
               ! i.e. set of atomic positions
               IF (print_atomic_displacements) THEN
                  DO iatom = 1, natom_dcd
                     atomic_displacement(iatom) = SQRT((r(iatom, 1) - r0(iatom, 1))**2 + &
                                                       (r(iatom, 2) - r0(iatom, 2))**2 + &
                                                       (r(iatom, 3) - r0(iatom, 3))**2)
                  END DO
               END IF

               ! Dump the atomic coordinates in DCD or XMOL/MOLDEN format
               IF (output_format_dcd) THEN
                  IF (have_unit_cell == 1) THEN
                     WRITE (UNIT=output_unit) a, gamma, b, beta, alpha, c
                  END IF
                  IF (print_atomic_displacements) THEN
                     ! This is a hack for VMD: dump the atomic displacements as the x-coordinates
                     ! of a DCD CORD file
                     WRITE (UNIT=output_unit) REAL(atomic_displacement(:), KIND=sp)
                     ! The y and z coordinates are filled with zeros
                     atomic_displacement(:) = 0.0_dp
                     WRITE (UNIT=output_unit) REAL(atomic_displacement(:), KIND=sp)
                     WRITE (UNIT=output_unit) REAL(atomic_displacement(:), KIND=sp)
                  ELSE
                     ! Dump DCD file information
                     DO i = 1, 3
                        WRITE (UNIT=output_unit) r(:, i)
                     END DO
                  END IF
               ELSE
                  IF (print_atomic_displacements) THEN
                     IF (info) THEN
                        WRITE (UNIT=error_unit, FMT="(A,2(/,A),T15,A)") &
                           "#", &
                           "# Displacements   :", &
                           "# "//TRIM(unit_string), "|r|"
                     END IF
                     IF (have_atomic_labels) THEN
                        DO iatom = 1, natom_dcd
                           WRITE (UNIT=output_unit, FMT="(A5,5X,F25.6)") ADJUSTL(atomic_label(iatom)), atomic_displacement(iatom)
                        END DO
                     ELSE
                        DO iatom = 1, natom_dcd
                           WRITE (UNIT=output_unit, FMT="(10X,F25.6)") atomic_displacement(iatom)
                        END DO
                     END IF
                  ELSE
                     IF (output_format_xmol) THEN
                        IF (have_cell_file .OR. (have_unit_cell == 1)) THEN
                           WRITE (UNIT=output_unit, FMT="(T2,I0,/,3(1X,F14.6),3(1X,F9.3))") &
                              natom_dcd, a, b, c, alpha, beta, gamma
                        ELSE
                           WRITE (UNIT=output_unit, FMT="(T2,I0,/,A,2(I0,A))") &
                              natom_dcd, "Frame: ", nframe_read + 1, ", Step: ", nframe, ", "//TRIM(remark_xyz)
                        END IF
                        DO iatom = 1, natom_dcd
                           WRITE (UNIT=output_unit, FMT=fmt_string) ADJUSTL(atomic_label(iatom)), r(iatom, 1:3)
                        END DO
                     ELSE
                        IF (info) THEN
                           WRITE (UNIT=error_unit, FMT="(A,T20,A)") &
                              "# "//TRIM(unit_string), "x              y              z"
                        END IF
                        IF (print_scaled_coordinates) THEN
                           DO iatom = 1, natom_dcd
                              WRITE (UNIT=output_unit, FMT=fmt_string) ADJUSTL(atomic_label(iatom)), s(iatom, 1:3)
                           END DO
                        ELSE IF (print_scaled_pbc_coordinates) THEN
                           DO iatom = 1, natom_dcd
                              WRITE (UNIT=output_unit, FMT=fmt_string) ADJUSTL(atomic_label(iatom)), s_pbc(iatom, 1:3)
                           END DO
                        ELSE
                           DO iatom = 1, natom_dcd
                              WRITE (UNIT=output_unit, FMT=fmt_string) ADJUSTL(atomic_label(iatom)), r(iatom, 1:3)
                           END DO
                        END IF
                     END IF
                  END IF
               END IF ! output_format_dcd

            END IF ! dump_frame

         END IF

         IF (dump_frame) nframe_read = nframe_read + 1

         ! Exit loop and stop processing, if the last (requested) frame was encountered
         IF (nframe >= last_frame) THEN
            nframe = nframe - 1
            CLOSE (UNIT=input_unit)
            EXIT dcd_file_loop
         END IF

      END DO frame_loop

      nframe = nframe - 1

      CLOSE (UNIT=input_unit)

      IF (iarg >= narg) EXIT dcd_file_loop

   END DO dcd_file_loop

   IF (info) THEN
      WRITE (UNIT=error_unit, FMT="(A,/,A,I0)") &
         "#", &
         "# Frames processed: ", nframe_read
   END IF

   IF (ekin .AND. TRIM(ADJUSTL(id_dcd)) == "VEL") THEN
      IF (info) THEN
         WRITE (UNIT=error_unit, FMT="(A,/,A,F12.3)") &
            "#", &
            "# T [K] all frames: ", tavg/REAL(nframe_read, KIND=dp)
      END IF
   END IF

   ! Print version information
   IF (info) THEN
      WRITE (UNIT=error_unit, FMT="(A)") &
         "#", &
         "# Normal termination of "//TRIM(version_info)
   END IF

   ! Close files
   IF (have_atomic_labels) CLOSE (UNIT=xyz_input_unit)
   IF (LEN_TRIM(out_file_name) > 0) CLOSE (UNIT=output_unit)

   ! Cleanup
   IF (ALLOCATED(atomic_label)) DEALLOCATE (atomic_label)
   IF (ALLOCATED(atomic_displacement)) DEALLOCATE (atomic_displacement)
   IF (ALLOCATED(atomic_mass)) DEALLOCATE (atomic_mass)
   IF (ALLOCATED(atomic_temperature)) DEALLOCATE (atomic_temperature)
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
      STOP "*** ABNORMAL PROGRAM TERMINATION of dumpdcd v3.1 ***"

   END SUBROUTINE abort_program

! **************************************************************************************************
!> \brief  Calculation of the angle between the vectors a and b.
!>         The angle is returned in radians.
!> \param a ...
!> \param b ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION angle(a, b) RESULT(angle_ab)

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
!> \param element_symbol ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_atomic_mass(element_symbol) RESULT(amass)
      ! Get the atomic mass amass for an element

      CHARACTER(LEN=*), INTENT(IN)                       :: element_symbol
      REAL(KIND=dp)                                      :: amass

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_atomic_mass'

      SELECT CASE (TRIM(element_symbol))
      CASE ("O ")
         amass = 15.9994_dp
      CASE ("U ")
         amass = 238.02891_dp
      CASE DEFAULT
         CALL abort_program(routineN, "Unknown element symbol found")
      END SELECT

      amass = amass*massunit

   END FUNCTION get_atomic_mass

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
!> \param a ...
!> \param b ...
!> \param c ...
!> \param alpha ...
!> \param beta ...
!> \param gamma ...
!> \param debug ...
!> \param info ...
!> \param pbc0 ...
!> \param h ...
!> \param hinv ...
! **************************************************************************************************
   SUBROUTINE pbc(r, r_pbc, s, s_pbc, a, b, c, alpha, beta, gamma, debug, info, pbc0, h, hinv)
      ! Apply the periodic boundary conditions (PBC) to the Cartesian coordinate array
      ! r given the cell edge lengths a, b, and c in Angstrom and the angles alpha (<b,c)),
      ! beta (<(a,c)), and gamma (<(a,b)) in degree.
      ! On output r_pbc is updated with the "PBCed" input coordinates, s with the scaled
      ! input coordinates, and s_pbc with scaled "PBCed" coordinates.
      ! If pbc0 is true then fold to the range [-l/2,+l/2[ (origin at box centre) else fold
      ! to the range [0,l[ (origin at lower left box corner).

      REAL(KIND=sp), DIMENSION(:, :), INTENT(IN)         :: r
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: r_pbc, s, s_pbc
      REAL(KIND=dp), INTENT(IN)                          :: a, b, c, alpha, beta, gamma
      LOGICAL, INTENT(IN)                                :: debug, info, pbc0
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: h, hinv

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'pbc'

      INTEGER                                            :: i, natom
      LOGICAL                                            :: orthorhombic
      REAL(KIND=dp)                                      :: deth

      natom = SIZE(r, 1)
      IF (SIZE(r, 2) /= 3) CALL abort_program(routineN, "Array dimension for r must be 3")

      ! Build h matrix

      h(:, :) = 0.0_dp

      IF ((ABS(alpha - 90.0_dp) < EPSILON(0.0_dp)) .AND. &
          (ABS(beta - 90.0_dp) < EPSILON(0.0_dp)) .AND. &
          (ABS(gamma - 90.0_dp) < EPSILON(0.0_dp))) THEN
         orthorhombic = .TRUE.
         h(1, 1) = a
         h(2, 2) = b
         h(3, 3) = c
      ELSE
         orthorhombic = .FALSE.
         CALL build_h_matrix(a, b, c, alpha, beta, gamma, h)
      END IF

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
         " -cell_file                     : Input file with cell information in CP2K format (.cell)", &
         " -debug, -d                     : Print debug information", &
         " -ekin                          : Dump just the ""temperature"" of each atom", &
         " -eformat                       : Print coordinates in scientific format", &
         " -eo                            : Write standard output and standard error to the same logical unit", &
         " -first_frame, -ff <int>        : Number of the first frame which is dumped", &
         " -help, -h                      : Print this information", &
         " -ignore_warnings               : Do not print warning messages, e.g. about inconsistencies", &
         " -info, -i                      : Print additional information for each frame (see also -debug flag)", &
         " -last_frame, -lf <int>         : Number of the last frame which is dumped", &
         " -md_time_step <real>           : Use this MD time step instead of the one from the DCD header", &
         " -output, -o <file_name>        : Name of the output file (default is stdout)", &
         " -output_format, -of <DCD|XMOL> : Output format for dump", &
         " -pbc                           : Apply the periodic boundary conditions (PBC) to each frame before it is dumped", &
         "                                  (origin at lower left)", &
         " -pbc0                          : Apply the periodic boundary conditions (PBC) to each frame before it is dumped", &
         "                                  (origin at box centre)", &
         " -stride <int>                  : Stride for frame dump (allows to skip frames, e.g. by dumping each 10th frame)", &
         " -trace_atoms <real>            : Print atoms which left the simulation box given a threshold value in scaled units", &
         " -vel2cord, -v2c                : Dump a VELocity DCD file as COoRDinate DCD file (hack which allows a digest by VMD)", &
         " -xyz_file, -xyz <file_name>    : Name of a reference XYZ file in XMOL format that provides the atomic labels", &
         ""

      WRITE (UNIT=*, FMT="(T2,A)") &
         "Usage examples:", &
         "", &
         " dumpdcd <optional flags> <DCD file(s)>", &
         "", &
         "Specific usage examples:", &
         "", &
         " dumpdcd project-pos-1.dcd (without atomic labels from XYZ file)", &
         " dumpdcd -xyz project.xyz project-pos-1.dcd (single DCD file)", &
         " dumpdcd -xyz project.xyz project-pos-1.dcd project-pos-2.dcd ... (multiple DCD files are dumped consecutively)", &
         " dumpdcd -xyz project.xyz -cell_file project-1.cell -of xmol project-pos-1.dcd", &
         " dumpdcd -info -xyz project.xyz project-pos-1.dcd project-pos-2.dcd (print additional information)", &
         " dumpdcd -debug -xyz project.xyz project-pos-1.dcd project-pos-2.dcd (print debug information)", &
         " dumpdcd -ekin -d -xyz project.xyz project-vel-1.dcd (print the ""temperature"" of each atom)", &
         " dumpdcd -ekin -xyz project.xyz project-vel-1.dcd (print just the temperature of each atom)", &
         " dumpdcd -first_frame 5 -last_frame 10 project-pos-1.dcd (just dump frame 5 to 10, ie. 6 frames in total)", &
         " dumpdcd -o outfile.xyz project-pos-1.dcd (write output to the file ""outfile.xyz"" instead of stdout)", &
         " dumpdcd -o test.xyz -output_format xmol -xyz ref.xyz -first 10 -last 10 test.dcd (dump 10th frame in XMOL format)", &
         " dumpdcd -of dcd -ff 10 -lf 20 test.dcd (dump the frames 10 to 20 in DCD format to the default output file output.dcd)", &
         " dumpdcd -o part.dcd -of dcd -ff 1 -lf 3 test.dcd (dump the frames 1 to 3 in DCD format to the output file part.dcd)", &
         " dumpdcd -o part.dcd -of dcd -first 10 -lf 100 -stride 10 test.dcd (dump the frames 10,..., 100 to the file part.dcd)", &
         " dumpdcd -output new.dcd -output_format dcd -pbc old.dcd (dump all frames applying PBC to the output file new.dcd)", &
         " dumpdcd -o new.dcd -of dcd -pbc -trace_atoms 0.02 old.dcd (all atoms more than 2% out of the box are listed)", &
         " dumpdcd -o new.dcd -e out_of_box.log -of dcd -pbc -trace_atoms 0.1 old.dcd (atoms >10% out of the box are listed)", &
         " dumpdcd -o new.dcd -of dcd -vel2cord old.dcd (dump old.dcd as new.dcd and change only the DCD id from VEL to CORD)", &
         " dumpdcd -o new.dcd -of dcd -frc2cord old.dcd (dump old.dcd as new.dcd and change only the DCD id from FRC to CORD)", &
         " dumpdcd -i -disp project-pos-1.dcd (dump the displacements of all atoms w.r.t. their positions in the first frame)", &
         " dumpdcd -i -of dcd -disp project-pos-1.dcd (dump the atomic displacements as x-coordinates of a DCD CORD file)", &
         " dumpdcd -i -of dcd -ekin -v2c -xyz project.xyz project-vel-1.dcd (dump the atomic temperatures as x-coordinates of a ", &
         "                                                                   DCD CORD file -> hack for VMD)", &
         ""

      WRITE (UNIT=*, FMT="(T2,A)") &
         "Notes:", &
         "", &
         " - For -ekin a XYZ file is required to obtain the atomic labels", &
         " - The -info and the -debug flags provide a more detailed output which is especially handy for tracing problems", &
         " - The output in DCD format is in binary format", &
         " - The input coordinates should be in Angstrom. Velocities and forces are expected to be in atomic units", &
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
      REAL(KIND=sp), DIMENSION(:, :), INTENT(IN)         :: r
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: s
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

END PROGRAM dumpdcd
