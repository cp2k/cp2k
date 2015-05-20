PROGRAM dumpdcd

! Copyright (C) 2012 - 2015  Matthias Krack (MK)
!
! Version: 3.0
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
!
! Note: For -ekin a XYZ file is required to obtain the atomic labels.
!       The -info and the -debug flags provide a more detailed output which is especially handy for tracing problems.
!       The output in DCD format is in binary format.
!       The input coordinates should be in Angstrom. Velocities and forces are expected to be in atomic units.

! Uncomment the following line if this module is available (e.g. with gfortran) and comment the corresponding variable declarations below
! USE ISO_FORTRAN_ENV, ONLY: error_unit,input_unit,output_unit

  IMPLICIT NONE

  ! Comment the following lines if the ISO_FORTRAN_ENV is used (see above)
  INTEGER, PARAMETER                                 :: default_error_unit  = 0,&
                                                        default_input_unit  = 5,&
                                                        default_output_unit = 6
  INTEGER                                            :: error_unit  = default_error_unit,&
                                                        input_unit  = default_input_unit,&
                                                        output_unit = default_output_unit
  ! End Comment

  ! Parameters
  CHARACTER(LEN=*), PARAMETER :: routine_name = "dumpdcd",&
                                 version_info = routine_name//" v3.0 (31.01.2014, Matthias Krack)"

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,200),&
                        sp = SELECTED_REAL_KIND(6,30)
  INTEGER, PARAMETER :: default_string_length        = 240,&
                        xyz_input_unit               = 10

  REAL(KIND=dp), PARAMETER :: pi       = 3.14159265358979323846264338_dp
  REAL(KIND=dp), PARAMETER :: angstrom = 0.52917720859_dp,&              ! [a.u.] -> [Angstrom]
                              degree   = 180.0_dp/pi,&                   ! [rad]  -> [degree]
                              kelvin   = 315774.647902944_dp,&           ! [a.u.] -> [K]
                              massunit = 1822.88484264550_dp             ! [u]    -> [a.u.]

  ! Variables
  CHARACTER(LEN=4)                                   :: id_dcd
  CHARACTER(LEN=10)                                  :: unit_string
  CHARACTER(LEN=default_string_length)               :: arg,dcd_file_name,message,out_file_name,output_format,&
                                                        remark_xyz,string,xyz_file_name
  CHARACTER(LEN=5), DIMENSION(:), ALLOCATABLE        :: atomic_label
  CHARACTER(LEN=80), DIMENSION(2)                    :: remark_dcd
  INTEGER                                            :: first_frame,have_unit_cell,i,iarg,&
                                                        iatom,iskip,istat,istep,last_frame,narg,&
                                                        natom_dcd,natom_xyz,ndcd_file,nframe,&
                                                        nframe_read,nremark,stride
  LOGICAL                                            :: apply_pbc,debug,dump_frame,ekin,eo,have_atomic_labels,&
                                                        info,opened,output_format_dcd,output_format_xmol,&
                                                        pbc0,print_atomic_displacements,print_scaled_coordinates,&
                                                        print_scaled_pbc_coordinates,trace_atoms,vel2cord
  REAL(KIND=sp)                                      :: dt
  REAL(KIND=dp)                                      :: a,alpha,b,beta,c,eps_out_of_box,gamma,tavg,tavg_frame,x
  INTEGER, DIMENSION(16)                             :: idum
  REAL(KIND=dp), DIMENSION(3)                        :: rdum
  REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: atomic_displacement,atomic_mass,atomic_temperature
  REAL(KIND=dp), DIMENSION(3,3)                      :: h,hinv
  REAL(KIND=sp), DIMENSION(:,:), ALLOCATABLE         :: r
  REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE         :: r_pbc,r0,s,s_pbc

  apply_pbc = .FALSE.
  pbc0 = .FALSE.
  debug = .FALSE.
  dump_frame = .TRUE.
  ekin = .FALSE.
  eo = .FALSE.
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
  idum(:) = 0
  rdum(:) = 0.0_dp
  id_dcd = ""
  dcd_file_name = ""
  out_file_name = ""
  xyz_file_name = ""
  output_format = "default"
  output_format_dcd = .FALSE.
  output_format_xmol = .FALSE.
  remark_dcd(:) = ""
  remark_xyz = ""
  dt = 0.0_sp
  a = 0.0_dp
  b = 0.0_dp
  c = 0.0_dp
  alpha = 0.0_dp
  beta  = 0.0_dp
  gamma = 0.0_dp
  eps_out_of_box = -HUGE(0.0_dp)
  tavg = 0.0_dp
  tavg_frame = 0.0_dp

  narg = command_argument_count()

  IF (narg == 0) THEN
    CALL print_help()
    CALL abort_program(routine_name,"No input file(s) specified")
  END IF

  iarg = 0

  dcd_file_loop: DO

    iarg = iarg + 1

    CALL get_command_argument(NUMBER=iarg,VALUE=arg,STATUS=istat)

    SELECT CASE (arg)
    CASE ("-debug","-d")
      debug = .TRUE.
      info = .TRUE.
      CYCLE dcd_file_loop
    CASE ("-displacements","-disp")
      print_atomic_displacements = .TRUE.
      CYCLE dcd_file_loop
    CASE ("-ekin")
      ekin = .TRUE.
      CYCLE dcd_file_loop
    CASE ("-eo")
      eo = .TRUE.
      CYCLE dcd_file_loop
    CASE ("-first_frame","-first","-ff")
      iarg = iarg + 1
      CALL get_command_argument(NUMBER=iarg,VALUE=arg,STATUS=istat)
      READ (UNIT=arg,FMT=*,ERR=100) first_frame
      IF (first_frame <= 0) THEN
        CALL abort_program(routine_name,"Invalid number for first frame specified: "//&
                                        "first_frame must be greater than zero")
      END IF
      CYCLE dcd_file_loop
100   CALL abort_program(routine_name,"Invalid number for first frame specified "//&
                                      "(an integer number greater than zero is expected)")
    CASE ("-help","-h")
      CALL print_help()
      STOP 
    CASE ("-info","-i")
      info = .TRUE.
      CYCLE dcd_file_loop
    CASE ("-last_frame","-last","-lf")
      iarg = iarg + 1
      CALL get_command_argument(NUMBER=iarg,VALUE=arg,STATUS=istat)
      READ (UNIT=arg,FMT=*,ERR=101) last_frame
      IF (last_frame <= 0) THEN
        CALL abort_program(routine_name,"Invalid number for last frame specified: "//&
                                        "last_frame must be greater than zero")
      END IF
      CYCLE dcd_file_loop
101   CALL abort_program(routine_name,"Invalid number for last frame specified "//&
                                      "(an integer number greater than zero is expected)")
    CASE ("-o","-output")
      iarg = iarg + 1
      CALL get_command_argument(NUMBER=iarg,VALUE=out_file_name,STATUS=istat)
      CYCLE dcd_file_loop
    CASE ("-output_format","-of")
      iarg = iarg + 1
      CALL get_command_argument(NUMBER=iarg,VALUE=output_format,STATUS=istat)
      CALL uppercase(output_format)
      SELECT CASE (output_format)
      CASE ("DCD")
        output_format_dcd = .TRUE.
        output_format_xmol = .FALSE.
      CASE ("XMOL")
        output_format_dcd = .FALSE.
        output_format_xmol = .TRUE.
      CASE DEFAULT
        CALL abort_program(routine_name,"Invalid output format type specified")
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
    CASE ("-scaled_coordinates","-sc")
      print_scaled_coordinates = .TRUE.
      CYCLE dcd_file_loop
    CASE ("-scaled_pbc_coordinates","-spc")
      print_scaled_pbc_coordinates = .TRUE.
      CYCLE dcd_file_loop
    CASE ("-stride")
      iarg = iarg + 1
      CALL get_command_argument(NUMBER=iarg,VALUE=arg,STATUS=istat)
      READ (UNIT=arg,FMT=*,ERR=104) stride
      IF (stride < 1) THEN
        CALL abort_program(routine_name,"Invalid stride for frame dump specified: stride must be greater than zero")
      END IF
      CYCLE dcd_file_loop
104   CALL abort_program(routine_name,"Invalid stride for frame dump specified "//&
                                      "(an integer number greater than 0 is expected)")
    CASE ("-trace_atoms")
      iarg = iarg + 1
      CALL get_command_argument(NUMBER=iarg,VALUE=arg,STATUS=istat)
      READ (UNIT=arg,FMT=*,ERR=108) eps_out_of_box
      IF (eps_out_of_box <= 0.0_dp) THEN
        CALL abort_program(routine_name,"Invalid threshold value for -trace_atoms flag specified")
      END IF
      trace_atoms = .TRUE.
      CYCLE dcd_file_loop
108   CALL abort_program(routine_name,"Invalid threshold value for -trace_atoms flag specified")
    CASE ("-vel2cord","-v2c")
      vel2cord = .TRUE.
      CYCLE dcd_file_loop
    CASE ("-xyz","xyz_file")
      iarg = iarg + 1
      CALL get_command_argument(NUMBER=iarg,VALUE=xyz_file_name,STATUS=istat)
      have_atomic_labels = .TRUE.
      CYCLE dcd_file_loop
    CASE DEFAULT
      IF (arg(1:1) == "-") THEN
        CALL print_help()
        CALL abort_program(routine_name,"Unknown command line flag """//TRIM(arg)//""" found")
      END IF
      dcd_file_name = arg
    END SELECT

    ! Check flag compatibility
    IF (first_frame > last_frame) THEN
      CALL abort_program(routine_name,"Number of first frame greater than number of last frame")
    END IF
    IF ((.NOT.have_atomic_labels).AND.output_format_xmol) THEN
      CALL abort_program(routine_name,"The output format XMOL requires a valid xyz file (-xyz flag)")
    END IF
    IF (output_format_xmol.AND.ekin) THEN
      CALL abort_program(routine_name,"Output format XMOL and the -ekin flag are incompatible")
    END IF
    IF (output_format_xmol.AND.print_atomic_displacements) THEN
      CALL abort_program(routine_name,"Output format XMOL and the -displacements flag are incompatible")
    END IF
    IF (output_format_xmol.AND.eo) THEN
      CALL abort_program(routine_name,"The -eo flag is incompatible with the output format XMOL")
    END IF
    IF (.NOT.have_atomic_labels.AND.ekin) THEN
      CALL abort_program(routine_name,"ekin flag requires also the specification of a valid xyz file (-xyz flag)")
    END IF
    IF (.NOT.apply_pbc.AND.trace_atoms) THEN
      CALL abort_program(routine_name,"The -trace_atoms flag requires the specification of a -pbc flag")
    END IF
    IF (ekin.AND.print_atomic_displacements) THEN
      CALL abort_program(routine_name,"The -ekin flag and the -displacements flag are incompatible")
    END IF
    IF (print_scaled_coordinates.AND.print_scaled_pbc_coordinates) THEN
      CALL abort_program(routine_name,"The -sc flag and the -spc flag are incompatible")
    END IF
    IF (.NOT.apply_pbc.AND.print_scaled_coordinates) THEN
      CALL abort_program(routine_name,"The -sc flag requires the specification of a -pbc flag")
    END IF
    IF (.NOT.apply_pbc.AND.print_scaled_pbc_coordinates) THEN
      CALL abort_program(routine_name,"The -spc flag requires the specification of a -pbc flag")
    END IF

    ! Open output units as requested
    IF (output_format_dcd) THEN
      ! Set default output file name if no file name was specified
      IF (LEN_TRIM(out_file_name) == 0) out_file_name = "output.dcd"
      ! Check if a new output file name was specified, if yes then close the old unit
      INQUIRE (UNIT=output_unit,NAME=string)
      IF (TRIM(string) /= TRIM(out_file_name)) CLOSE (UNIT=output_unit)
      INQUIRE (UNIT=output_unit,OPENED=opened)
      IF (.NOT.opened) THEN
        OPEN (UNIT=output_unit,&
              FILE=out_file_name,&
              STATUS="UNKNOWN",&
              ACCESS="SEQUENTIAL",&
              FORM="UNFORMATTED",&
              ACTION="WRITE",&
              IOSTAT=istat)
        IF (istat /= 0) CALL abort_program(routine_name,"The unformatted output file could not be opened")
      END IF
    ELSE
      IF (eo) error_unit = output_unit
      IF (LEN_TRIM(out_file_name) > 0) THEN
        ! Check if a new output file name was specified, if yes then close the old unit
        INQUIRE (UNIT=output_unit,NAME=string)
        IF (TRIM(string) /= TRIM(out_file_name)) CLOSE (UNIT=output_unit)
        INQUIRE (UNIT=output_unit,OPENED=opened)
        IF (.NOT.opened) THEN
          OPEN (UNIT=output_unit,&
                FILE=out_file_name,&
                STATUS="UNKNOWN",&
                ACCESS="SEQUENTIAL",&
                FORM="FORMATTED",&
                POSITION="REWIND",&
                ACTION="WRITE",&
                IOSTAT=istat)
          IF (istat /= 0) CALL abort_program(routine_name,"The formatted output file could not be opened")
        END IF
      END IF
    END IF

    ! Avoid reading from and writing to the same file
    IF (TRIM(dcd_file_name) == TRIM(out_file_name)) THEN
      CALL abort_program(routine_name,"Input and output file name cannot be the same")
    END IF

    ! Read reference xyz input file if requested
    IF (have_atomic_labels) THEN
      string = ""
      ! Check if a new XYZ input file name was specified, if yes then close the old unit
      INQUIRE (UNIT=xyz_input_unit,NAME=string)
      IF (TRIM(string) /= TRIM(xyz_file_name)) CLOSE (UNIT=xyz_input_unit)
      INQUIRE (UNIT=xyz_input_unit,OPENED=opened)
      ! Read new XYZ input file if needed and update reference information
      IF (.NOT.opened) THEN
        OPEN (UNIT=xyz_input_unit,&
              FILE=xyz_file_name,&
              STATUS="OLD",&
              ACCESS="SEQUENTIAL",&
              FORM="FORMATTED",&
              POSITION="REWIND",&
              ACTION="READ",&
              IOSTAT=istat)
        IF (istat /= 0) CALL abort_program(routine_name,"The XYZ file could not be opened")
        IF (info) WRITE (UNIT=error_unit,FMT="(A)") "#","# Reading XYZ file: "//TRIM(xyz_file_name)
        READ (UNIT=xyz_input_unit,FMT="(A)",IOSTAT=istat) arg
        IF (istat /= 0) THEN
          CALL abort_program(routine_name,"Reading line 1 of the XYZ file (number of atoms) failed")
        END IF
        IF (arg(1:1) == "#") THEN
          READ (UNIT=arg,FMT=*,IOSTAT=istat) string,natom_xyz
        ELSE
          READ (UNIT=arg,FMT=*,IOSTAT=istat) natom_xyz
        END IF
        IF (istat /= 0) THEN
          CALL abort_program(routine_name,"Reading line 1 of the XYZ file (number of atoms) failed")
        END IF
        IF (ALLOCATED(atomic_label)) THEN
          DEALLOCATE (atomic_label,STAT=istat)
          IF (istat /= 0) CALL abort_program(routine_name,"Deallocation of the vector atomic_label failed")
        END IF
        ALLOCATE (atomic_label(natom_xyz),STAT=istat)
        IF (istat /= 0) CALL abort_program(routine_name,"Allocation of the vector atomic_label failed")
        atomic_label(:) = ""
        IF (ekin) THEN
          IF (ALLOCATED(atomic_mass)) THEN
            DEALLOCATE (atomic_mass,STAT=istat)
            IF (istat /= 0) CALL abort_program(routine_name,"Deallocation of the vector atomic_mass failed")
          END IF
          ALLOCATE (atomic_mass(natom_xyz),STAT=istat)
          IF (istat /= 0) CALL abort_program(routine_name,"Allocation of the vector atomic_mass failed")
        END IF
        READ (UNIT=xyz_input_unit,FMT="(A)",IOSTAT=istat) remark_xyz
        IF (istat /= 0) CALL abort_program(routine_name,"Reading line 2 of the XYZ file (remark line) failed")
        iatom = 0
        DO
          READ (UNIT=xyz_input_unit,FMT=*,IOSTAT=istat) arg
          IF (arg(1:1) == "#") THEN
            CYCLE
          ELSE
            BACKSPACE (UNIT=xyz_input_unit)
          END IF
          iatom = iatom + 1
          READ (UNIT=xyz_input_unit,FMT=*,IOSTAT=istat) atomic_label(iatom),rdum(1:3)
          IF (istat /= 0) THEN
            message = ""
            WRITE (UNIT=message,FMT="(A,I0,A,I0,A)")&
             "Reading line ",iatom+2," of the XYZ file (atom ",iatom,") failed"
            CALL abort_program(routine_name,TRIM(message))
          END IF
          CALL uppercase(atomic_label(iatom))
          IF (LEN_TRIM(atomic_label(iatom)) > 1) CALL lowercase(atomic_label(iatom)(2:2))
          atomic_label(iatom) = TRIM(ADJUSTL(atomic_label(iatom)))
          IF (ekin) atomic_mass(iatom) = get_atomic_mass(atomic_label(iatom))
          IF (iatom == natom_xyz) EXIT
        END DO
        IF (info) THEN
          WRITE (UNIT=error_unit,FMT="(A,I0)")&
           "# Number of atoms : ",natom_xyz,&
           "# Remark          : "//TRIM(ADJUSTL(remark_xyz))
        END IF
      END IF
    END IF

    ! Increment DCD file counter
    ndcd_file = ndcd_file + 1

    ! Open input file in DCD format
    OPEN (UNIT=input_unit,&
          FILE=dcd_file_name,&
          STATUS="OLD",&
          FORM="UNFORMATTED",&
          POSITION="REWIND",&
          ACTION="READ",&
          IOSTAT=istat)
    IF (istat /= 0) CALL abort_program(routine_name,"The DCD file could not be opened")
    IF (info) THEN
      WRITE (UNIT=error_unit,FMT="(A,/,(A,I0))")&
       "#",&
       "# DCD file number : ",ndcd_file,&
       "# Reading DCD file: "//TRIM(dcd_file_name)
    END IF

    ! Read 1st record of DCD file
    READ (UNIT=input_unit) id_dcd,idum(1),istep,iskip,idum(2:7),dt,have_unit_cell,idum(8:16)
    IF (info) THEN
      WRITE (UNIT=error_unit,FMT="(A,2(/,A,I0),/,A,F9.3,/,A,I0)")&
       "# DCD id string   : "//id_dcd,&
       "# Step            : ",istep,&
       "# Print frequency : ",iskip,&
       "# Time step [fs]  : ",dt,&
       "# Unit cell       : ",have_unit_cell
    END IF

    IF (ekin.AND.TRIM(ADJUSTL(id_dcd)) /= "VEL") THEN
      CALL abort_program(routine_name,"ekin flag requires a DCD file with VELocities")
    END IF
    IF (apply_pbc.AND.(have_unit_cell /= 1)) THEN
      CALL abort_program(routine_name,"pbc flags require that unit cell information is available")
    END IF

    IF (TRIM(ADJUSTL(id_dcd)) == "VEL") THEN
      unit_string = "[a.u.]"
      IF (apply_pbc) CALL abort_program(routine_name,"pbc flags require a DCD file with COoRDinates")
    ELSE IF (TRIM(ADJUSTL(id_dcd)) == "CORD") THEN
      unit_string = "[Angstrom]"
    ELSE
      CALL abort_program(routine_name,"Unknown DCD id found (use -debug or -info flag for details)")
    END IF

    ! Read 2nd record of DCD file
    READ (UNIT=input_unit) nremark,remark_dcd(1),remark_dcd(2)
    IF (info) THEN
      DO i=1,nremark
        WRITE (UNIT=error_unit,FMT="(A,I1,A)")&
         "# Remark ",i,"        : "//TRIM(remark_dcd(i))
      END DO
    END IF

    ! Read 3rd record of DCD file
    READ (UNIT=input_unit) natom_dcd
    IF (info) THEN
      WRITE (UNIT=error_unit,FMT="(A,I0)")&
       "# Number of atoms : ",natom_dcd
    END IF

    IF (have_atomic_labels) THEN
      IF (natom_dcd /= natom_xyz) CALL abort_program(routine_name,"Number of atoms in XYZ and DCD file differ")
    ELSE
      IF (.NOT.ALLOCATED(atomic_label)) THEN
        ALLOCATE (atomic_label(natom_dcd),STAT=istat)
        IF (istat /= 0) CALL abort_program(routine_name,"Allocation of the vector atomic_label failed")
        atomic_label(:) = ""
      END IF
    END IF

    ! Dump output in DCD format => dcd2dcd functionality
    IF (output_format_dcd) THEN
      IF (vel2cord) id_dcd = "CORD"
      ! Note, dt is REAL*4 and the rest is INTEGER*4
      WRITE (UNIT=output_unit) id_dcd,idum(1),istep,iskip,idum(2:7),dt,have_unit_cell,idum(8:16)
      WRITE (UNIT=output_unit) nremark,remark_dcd(1),remark_dcd(2)
      WRITE (UNIT=output_unit) natom_dcd
    END IF

    frame_loop: DO

      nframe = nframe + 1

      IF (nframe < first_frame) THEN
        dump_frame = .FALSE.
      ELSE
        IF (MODULO(nframe-first_frame,stride) == 0) THEN
          dump_frame = .TRUE.
        ELSE
          dump_frame = .FALSE.
        END IF
      END IF

      ! Read unit cell information, if available
      IF (have_unit_cell == 1) THEN
        READ (UNIT=input_unit,IOSTAT=istat) a,gamma,b,beta,alpha,c
        IF (istat < 0) EXIT frame_loop
      END IF

      IF ((info.OR.trace_atoms).AND.dump_frame) THEN
        WRITE (UNIT=error_unit,FMT="(A,/,A,I0)")&
         "#","# Frame number    : ",nframe
      END IF

      ! Print unit cell information, if available
      IF (info.AND.dump_frame) THEN
        IF (have_unit_cell == 1) THEN
          WRITE (UNIT=error_unit,FMT="(A,/,(A,F12.6))")&
           "#",&
           "# a [Angstrom]    : ",a,&
           "# b [Angstrom]    : ",b,&
           "# c [Angstrom]    : ",c,&
           "# alpha [degree]  : ",alpha,&
           "# beta  [degree]  : ",beta,&
           "# gamma [degree]  : ",gamma
        END IF
      END IF

      ! Allocate the array for the current atomic positions if needed
      IF (.NOT.ALLOCATED(r)) THEN
        ALLOCATE (r(natom_dcd,3),STAT=istat)
        IF (istat /= 0) CALL abort_program(routine_name,"Allocation of the array r failed")
      END IF

      ! Read in the atomic positions of the current frame
      READ (UNIT=input_unit,IOSTAT=istat) r(1:natom_dcd,1)
      IF (istat < 0) EXIT frame_loop
      READ (UNIT=input_unit) r(1:natom_dcd,2)
      READ (UNIT=input_unit) r(1:natom_dcd,3)

      ! Store the atomic positions of the first frame in the current DCD input file for
      ! the output of the atomic displacements
      IF ((nframe == 1).AND.print_atomic_displacements) THEN
        IF (.NOT.ALLOCATED(r0)) THEN
          ALLOCATE (r0(natom_dcd,3),STAT=istat)
          IF (istat /= 0) CALL abort_program(routine_name,"Allocation of the array r0 failed")
        END IF
        r0(:,:) = r(:,:)
        IF (.NOT.ALLOCATED(atomic_displacement)) THEN
          ALLOCATE (atomic_displacement(natom_dcd),STAT=istat)
          IF (istat /= 0) THEN
            CALL abort_program(routine_name,"Allocation of the vector atomic_displacement failed")
          END IF
        END IF
        atomic_displacement(:) = 0.0_dp
      END IF

      IF (ekin.AND.TRIM(ADJUSTL(id_dcd)) == "VEL") THEN

        ! Dump the kinetic energy of each as an "atomic temperature"
        IF (dump_frame) THEN

          IF (.NOT.ALLOCATED(atomic_temperature)) THEN
            ALLOCATE (atomic_temperature(natom_dcd),STAT=istat)
            IF (istat /= 0) THEN
              CALL abort_program(routine_name,"Allocation of the vector atomic_temperature failed")
            END IF
            atomic_temperature(:) = 0.0_dp
          END IF

          IF (info) THEN
            WRITE (UNIT=error_unit,FMT="(A)")&
             "#",&
             "# Temperature [K]"
          END IF

          tavg_frame = 0.0_dp
          DO iatom=1,natom_dcd
             atomic_temperature(iatom) = atomic_mass(iatom)*(r(iatom,1)*r(iatom,1) +&
                                                             r(iatom,2)*r(iatom,2) +&
                                                             r(iatom,3)*r(iatom,3))*kelvin/3.0_dp
             tavg_frame = tavg_frame + atomic_temperature(iatom)
          END DO
          tavg_frame = tavg_frame/REAL(natom_dcd,KIND=dp)

          IF (output_format_dcd) THEN
            IF (have_unit_cell == 1) THEN
              WRITE (UNIT=output_unit) a,gamma,b,beta,alpha,c
            END IF
            ! This is a hack for VMD: dump the atomic temperatures as the x-coordinates
            ! of a DCD VEL file
            WRITE (UNIT=output_unit) REAL(atomic_temperature(:),KIND=sp)
            ! The y and z coordinates are filled with zeros
            atomic_temperature(:) = 0.0_dp
            WRITE (UNIT=output_unit) REAL(atomic_temperature(:),KIND=sp)
            WRITE (UNIT=output_unit) REAL(atomic_temperature(:),KIND=sp)
          ELSE
            DO iatom=1,natom_dcd
              WRITE (UNIT=output_unit,FMT="(A5,5X,F25.3)") ADJUSTL(atomic_label(iatom)),atomic_temperature
            END DO
          END IF

          IF (info) THEN
            WRITE (UNIT=error_unit,FMT="(A,F12.3)")&
             "# T [K] this frame: ",tavg_frame
          END IF

          tavg = tavg + tavg_frame

        END IF ! dump_frame

      ELSE

        IF (dump_frame) THEN

          ! Apply periodic boundary conditions before dumping the coordinates if requested
          IF (apply_pbc) THEN
            IF (.NOT.ALLOCATED(r_pbc)) THEN
              ALLOCATE (r_pbc(natom_dcd,3),STAT=istat)
              IF (istat /= 0) CALL abort_program(routine_name,"Allocation of the array r_pbc failed")
            END IF
            IF (.NOT.ALLOCATED(s)) THEN
              ALLOCATE (s(natom_dcd,3),STAT=istat)
              IF (istat /= 0) CALL abort_program(routine_name,"Allocation of the array s failed")
            END IF
            IF (.NOT.ALLOCATED(s_pbc)) THEN
              ALLOCATE (s_pbc(natom_dcd,3),STAT=istat)
              IF (istat /= 0) CALL abort_program(routine_name,"Allocation of the array s_pbc failed")
            END IF
            CALL pbc(r,r_pbc,s,s_pbc,a,b,c,alpha,beta,gamma,debug,info,pbc0,h,hinv)
            CALL write_out_of_box_atoms(atomic_label,r,s,eps_out_of_box,h)
            ! Overwrite input coordinate with the PBCed coordinates for printing
            r(:,:) = REAL(r_pbc(:,:),KIND=sp)
          END IF ! apply_pbc

          ! Calculate the atomic displacements with respect to the first frame,
          ! i.e. set of atomic positions
          IF (print_atomic_displacements) THEN
            DO iatom=1,natom_dcd
              atomic_displacement(iatom) = SQRT((r(iatom,1) - r0(iatom,1))**2 +&
                                                (r(iatom,2) - r0(iatom,2))**2 +&
                                                (r(iatom,3) - r0(iatom,3))**2)
            END DO
          END IF

          ! Dump the atomic coordinates in DCD or XMOL/MOLDEN format
          IF (output_format_dcd) THEN
            IF (have_unit_cell == 1) THEN
              WRITE (UNIT=output_unit) a,gamma,b,beta,alpha,c
            END IF
            IF (print_atomic_displacements) THEN
              ! This is a hack for VMD: dump the atomic displacements as the x-coordinates
              ! of a DCD CORD file
              WRITE (UNIT=output_unit) REAL(atomic_displacement(:),KIND=sp)
              ! The y and z coordinates are filled with zeros
              atomic_displacement(:) = 0.0_dp
              WRITE (UNIT=output_unit) REAL(atomic_displacement(:),KIND=sp)
              WRITE (UNIT=output_unit) REAL(atomic_displacement(:),KIND=sp)
            ELSE
              ! Dump DCD file information
              DO i=1,3
                WRITE (UNIT=output_unit) r(:,i)
              END DO
            END IF
          ELSE
            IF (print_atomic_displacements) THEN
              IF (info) THEN
                WRITE (UNIT=error_unit,FMT="(A,2(/,A),T15,A)")&
                 "#",&
                 "# Displacements   :",&
                 "# "//TRIM(unit_string),"|r|"
              END IF
              IF (have_atomic_labels) THEN
                DO iatom=1,natom_dcd
                  WRITE (UNIT=output_unit,FMT="(A5,5X,F25.6)") ADJUSTL(atomic_label(iatom)),atomic_displacement(iatom)
                END DO
              ELSE
                DO iatom=1,natom_dcd
                  WRITE (UNIT=output_unit,FMT="(10X,F25.6)") atomic_displacement(iatom)
                END DO
              END IF
            ELSE
              IF (output_format_xmol) THEN
                WRITE (UNIT=output_unit,FMT="(T2,I0,/,A)") natom_dcd,TRIM(remark_xyz)
                DO iatom=1,natom_dcd
                  WRITE (UNIT=output_unit,FMT="(A5,3(1X,F14.6))") ADJUSTL(atomic_label(iatom)),r(iatom,1:3)
                END DO
              ELSE
                IF (info) THEN
                  WRITE (UNIT=error_unit,FMT="(A,T20,A)")&
                   "# "//TRIM(unit_string),"x              y              z"
                END IF
                IF (print_scaled_coordinates) THEN
                  DO iatom=1,natom_dcd
                    WRITE (UNIT=output_unit,FMT="(A5,3(1X,F14.6))") ADJUSTL(atomic_label(iatom)),s(iatom,1:3)
                  END DO
                ELSE IF (print_scaled_pbc_coordinates) THEN
                  DO iatom=1,natom_dcd
                    WRITE (UNIT=output_unit,FMT="(A5,3(1X,F14.6))") ADJUSTL(atomic_label(iatom)),s_pbc(iatom,1:3)
                  END DO
                ELSE
                  DO iatom=1,natom_dcd
                    WRITE (UNIT=output_unit,FMT="(A5,3(1X,F14.6))") ADJUSTL(atomic_label(iatom)),r(iatom,1:3)
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
    WRITE (UNIT=error_unit,FMT="(A,/,A,I0)")&
     "#",&
     "# Frames processed: ",nframe_read
  END IF

  IF (ekin.AND.TRIM(ADJUSTL(id_dcd)) == "VEL") THEN
    IF (info) THEN
      WRITE (UNIT=error_unit,FMT="(A,/,A,F12.3)")&
       "#",&
       "# T [K] all frames: ",tavg/REAL(nframe_read,KIND=dp)
    END IF
  END IF

  ! Print version information
  IF (info) THEN
    WRITE (UNIT=error_unit,FMT="(A)")&
     "#",&
     "# Normal termination of "//TRIM(version_info)
  END IF

  ! Close files
  IF (have_atomic_labels) CLOSE (UNIT=xyz_input_unit)
  IF (LEN_TRIM(out_file_name) > 0) CLOSE (UNIT=output_unit)

  ! Cleanup
  IF (ALLOCATED(atomic_label)) THEN
    DEALLOCATE (atomic_label,STAT=istat)
    IF (istat /= 0) THEN
      CALL abort_program(routine_name,"Deallocation of the vector atomic_label failed")
    END IF
  END IF

  IF (ALLOCATED(atomic_displacement)) THEN
    DEALLOCATE (atomic_displacement,STAT=istat)
    IF (istat /= 0) THEN
      CALL abort_program(routine_name,"Deallocation of the vector atomic_displacement failed")
    END IF
  END IF

  IF (ALLOCATED(atomic_mass)) THEN
    DEALLOCATE (atomic_mass,STAT=istat)
    IF (istat /= 0) THEN
      CALL abort_program(routine_name,"Deallocation of the vector atomic_mass failed")
    END IF
  END IF

  IF (ALLOCATED(atomic_temperature)) THEN
    DEALLOCATE (atomic_temperature,STAT=istat)
    IF (istat /= 0) THEN
      CALL abort_program(routine_name,"Deallocation of the vector atomic_temperature failed")
    END IF
  END IF

  IF (ALLOCATED(r)) THEN
    DEALLOCATE (r,STAT=istat)
    IF (istat /= 0) THEN
      CALL abort_program(routine_name,"Deallocation of the array r failed")
    END IF
  END IF

  IF (ALLOCATED(r0)) THEN
    DEALLOCATE (r0,STAT=istat)
    IF (istat /= 0) THEN
      CALL abort_program(routine_name,"Deallocation of the array r0 failed")
    END IF
  END IF

  IF (ALLOCATED(r_pbc)) THEN
    DEALLOCATE (r_pbc,STAT=istat)
    IF (istat /= 0) THEN
      CALL abort_program(routine_name,"Deallocation of the array r_pbc failed")
    END IF
  END IF

  IF (ALLOCATED(s)) THEN
    DEALLOCATE (s,STAT=istat)
    IF (istat /= 0) THEN
      CALL abort_program(routine_name,"Deallocation of the array s failed")
    END IF
  END IF

  IF (ALLOCATED(s_pbc)) THEN
    DEALLOCATE (s_pbc,STAT=istat)
    IF (istat /= 0) THEN
      CALL abort_program(routine_name,"Deallocation of the array s_pbc failed")
    END IF
  END IF

CONTAINS

  SUBROUTINE abort_program(routine,message)
    ! Abort the program after printing an error message to stderr

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)           :: routine,message

    CHARACTER(LEN=2*default_string_length) :: error_message

    error_message = "*** ERROR in "//TRIM(routine)//": "//TRIM(message)//" ***"
    WRITE (UNIT=default_error_unit,FMT="(/,A,/)") TRIM(error_message)
    STOP "*** ABNORMAL PROGRAM TERMINATION of dumpdcd v3.0 ***"

  END SUBROUTINE abort_program

  SUBROUTINE build_h_matrix(a,b,c,alpha,beta,gamma,h)
    ! Calculate the h matrix of a simulation cell given the cell edge lengths a, b,
    ! and c in Angstrom and the angles alpha (<b,c)), beta (<(a,c)), and gamma (<(a,b))
    ! in degree.

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER                :: routine_name = "build_h_matrix"

    REAL(KIND=dp), INTENT(IN)                  :: a,alpha,b,beta,c,gamma
    REAL(KIND=dp), DIMENSION(3,3), INTENT(OUT) :: h

    REAL(KIND=dp)                              :: cosa,cosb,cosg,sing

    cosa = COS(alpha/degree)
    IF (ABS(cosa) < EPSILON(0.0_dp)) cosa = 0.0_dp

    cosb = COS(beta/degree)
    IF (ABS(cosb) < EPSILON(0.0_dp)) cosb = 0.0_dp

    cosg = COS(gamma/degree)
    IF (ABS(cosg) < EPSILON(0.0_dp)) cosg = 0.0_dp

    sing = SIN(gamma/degree)
    IF (ABS(sing) < EPSILON(0.0_dp)) sing = 0.0_dp

    h(1,1) = 1.0_dp
    h(2,1) = 0.0_dp
    h(3,1) = 0.0_dp

    h(1,2) = cosg
    h(2,2) = sing
    h(3,2) = 0.0_dp

    h(1,3) = cosb
    h(2,3) = (cosa - cosg*cosb)/sing
    IF ((1.0_dp - h(1,3)**2 - h(2,3)**2) < 0.0_dp) THEN
      CALL abort_program(routine_name,"Build of the h matrix failed, check cell information")
    END IF
    h(3,3) = SQRT(1.0_dp - h(1,3)**2 - h(2,3)**2)

    h(:,1) = a*h(:,1)
    h(:,2) = b*h(:,2)
    h(:,3) = c*h(:,3)

  END SUBROUTINE build_h_matrix

  FUNCTION det_3x3(a) RESULT(det_a)
    ! Returns the determinante of the 3x3 matrix a.

    IMPLICIT NONE

    REAL(KIND=dp), DIMENSION(3,3), INTENT(IN) :: a

    REAL(KIND=dp) :: det_a

    det_a = a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) +&
            a(1,2)*(a(2,3)*a(3,1) - a(2,1)*a(3,3)) +&
            a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))

  END FUNCTION det_3x3

  FUNCTION get_atomic_mass(element_symbol) RESULT(amass)
    ! Get the atomic mass amass for an element

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: element_symbol

    CHARACTER(LEN=*), PARAMETER  :: routine_name = "get_atomic_mass"

    REAL(KIND=dp)                :: amass

    SELECT CASE (TRIM(element_symbol))
    CASE ("O ")
      amass = 15.9994_dp
    CASE ("U ")
      amass = 238.02891_dp
    CASE DEFAULT
      CALL abort_program(routine_name,"Unknown element symbol found")
    END SELECT

    amass = amass*massunit

  END FUNCTION get_atomic_mass

  SUBROUTINE invert_matrix_3x3(h,hinv,deth)
    ! Calculate the inverse hinv and the determinant deth of the 3x3 matrix h.

    IMPLICIT NONE

    REAL(KIND=dp), DIMENSION(3,3), INTENT(IN)  :: h
    REAL(KIND=dp), DIMENSION(3,3), INTENT(OUT) :: hinv
    REAL(KIND=dp), INTENT(OUT)                 :: deth

    CHARACTER(LEN=*), PARAMETER                :: routine_name = "invert_matrix_3x3"

    deth = det_3x3(h)

    ! Numerics
    deth = ABS(deth)
    IF (deth < 1.0E-10_dp) THEN
      CALL abort_program(routine_name,"Invalid h matrix for cell found; det(h) < 1.0E-10")
    END IF

    hinv(1,1) = (h(2,2)*h(3,3) - h(3,2)*h(2,3))/deth
    hinv(2,1) = (h(2,3)*h(3,1) - h(3,3)*h(2,1))/deth
    hinv(3,1) = (h(2,1)*h(3,2) - h(3,1)*h(2,2))/deth

    hinv(1,2) = (h(1,3)*h(3,2) - h(3,3)*h(1,2))/deth
    hinv(2,2) = (h(1,1)*h(3,3) - h(3,1)*h(1,3))/deth
    hinv(3,2) = (h(1,2)*h(3,1) - h(3,2)*h(1,1))/deth

    hinv(1,3) = (h(1,2)*h(2,3) - h(2,2)*h(1,3))/deth
    hinv(2,3) = (h(1,3)*h(2,1) - h(2,3)*h(1,1))/deth
    hinv(3,3) = (h(1,1)*h(2,2) - h(2,1)*h(1,2))/deth

  END SUBROUTINE invert_matrix_3x3

  SUBROUTINE lowercase(string)
    ! Convert all letters in a string to lowercase
    CHARACTER(LEN=*), INTENT(INOUT)          :: string

    INTEGER                                  :: i, iascii

    DO i=1,LEN_TRIM(string)
       iascii = ICHAR(string(i:i))
       IF ((iascii >= 65).AND.(iascii <= 90)) THEN
          string(i:i) = CHAR(iascii + 32)
       END IF
    END DO

  END SUBROUTINE lowercase

  SUBROUTINE pbc(r,r_pbc,s,s_pbc,a,b,c,alpha,beta,gamma,debug,info,pbc0,h,hinv)
    ! Apply the periodic boundary conditions (PBC) to the Cartesian coordinate array
    ! r given the cell edge lengths a, b, and c in Angstrom and the angles alpha (<b,c)),
    ! beta (<(a,c)), and gamma (<(a,b)) in degree.
    ! On output r_pbc is updated with the "PBCed" input coordinates, s with the scaled
    ! input coordinates, and s_pbc with scaled "PBCed" coordinates.
    ! If pbc0 is true then fold to the range [-l/2,+l/2[ (origin at box centre) else fold
    ! to the range [0,l[ (origin at lower left box corner).

    IMPLICIT NONE

    REAL(KIND=sp), DIMENSION(:,:), INTENT(IN)    :: r
    REAL(KIND=dp), DIMENSION(:,:), INTENT(OUT)   :: r_pbc,s,s_pbc
    REAL(KIND=dp), INTENT(IN)                    :: a,alpha,b,beta,c,gamma
    LOGICAL, INTENT(IN)                          :: debug,info,pbc0
    REAL(KIND=dp), DIMENSION(3,3), INTENT(OUT)   :: h,hinv

    CHARACTER(LEN=*), PARAMETER                  :: routine_name = "pbc"

    REAL(KIND=dp)                                :: deth
    INTEGER                                      :: i,natom
    LOGICAL                                      :: orthorhombic

    natom= SIZE(r,1)
    IF (SIZE(r,2) /= 3) CALL abort_program(routine_name,"Array dimension for r must be 3")

    ! Build h matrix

    h(:,:) = 0.0_dp

    IF ((alpha == 90.0_dp).AND.(beta == 90.0_dp).AND.(gamma == 90.0_dp)) THEN
      orthorhombic = .TRUE.
      h(1,1) = a
      h(2,2) = b
      h(3,3) = c
    ELSE
      orthorhombic = .FALSE.
      CALL build_h_matrix(a,b,c,alpha,beta,gamma,h)
    END IF

    ! Build inverse of h
    hinv(:,:) = 0.0_dp
    CALL invert_matrix_3x3(h,hinv,deth)

    IF (info) THEN
      WRITE (UNIT=error_unit,FMT="(A)")          "#"
      IF (orthorhombic) THEN
        WRITE (UNIT=error_unit,FMT="(A)")        "# Cell symmetry   : orthorhombic"
      ELSE
        WRITE (UNIT=error_unit,FMT="(A)")        "# Cell symmetry   : non-orthorhombic"
      END IF
      IF (debug) THEN
        WRITE (UNIT=error_unit,FMT="(A)")          "#"
        WRITE (UNIT=error_unit,FMT="(A,3F12.6,A)") "#          / ",h(1,:)," \"
        WRITE (UNIT=error_unit,FMT="(A,3F12.6,A)") "# h      = | ",h(2,:)," |"
        WRITE (UNIT=error_unit,FMT="(A,3F12.6,A)") "#          \ ",h(3,:)," /"
        WRITE (UNIT=error_unit,FMT="(A)")          "#"
        WRITE (UNIT=error_unit,FMT="(A,3F12.6,A)") "#          / ",hinv(1,:)," \"
        WRITE (UNIT=error_unit,FMT="(A,3F12.6,A)") "# Inv(h) = | ",hinv(2,:)," |"
        WRITE (UNIT=error_unit,FMT="(A,3F12.6,A)") "#          \ ",hinv(3,:)," /"
        WRITE (UNIT=error_unit,FMT="(A)")          "#"
        WRITE (UNIT=error_unit,FMT="(A,F0.6)")     "# det(h) = ",deth
      END IF
    END IF

    ! Calculate scaled coordinates and wrap back all atoms, i.e. apply the PBC
    IF (orthorhombic) THEN
      ! Try to save some flops in the case of an orthorhombic box
      DO i=1,3
        s(:,i) = r(:,i)*hinv(i,i)
      END DO
    ELSE
      s(:,:) = MATMUL(r(:,:),TRANSPOSE(hinv(:,:)))
    END IF

    IF (pbc0) THEN
      s_pbc(:,:) = s(:,:) - ANINT(s(:,:))
    ELSE
      s_pbc(:,:) = s(:,:) - FLOOR(s(:,:))
    END IF

    IF (orthorhombic) THEN
      DO i=1,3
        r_pbc(:,i) = s_pbc(:,i)*h(i,i)
      END DO
    ELSE
      r_pbc(:,:) = MATMUL(s_pbc(:,:),TRANSPOSE(h(:,:)))
    END IF

  END SUBROUTINE pbc

  SUBROUTINE print_help()
    ! Print the program flags for help

    IMPLICIT NONE

    WRITE (UNIT=*,FMT="(T2,A)")&
     "",&
     "Program flags for "//TRIM(version_info)//":",&
     "",&
     " -debug, -d                     : Print debug information",&
     " -ekin                          : Dump just the ""temperature"" of each atom",&
     " -eo                            : Write standard output and standard error to the same logical unit",&
     " -first_frame, -ff <int>        : Number of the first frame which is dumped",&
     " -help, -h                      : Print this information",&
     " -info, -i                      : Print additional information for each frame (see also -debug flag)",&
     " -last_frame, -lf <int>         : Number of the last frame which is dumped",&
     " -output, -o <file_name>        : Name of the output file (default is stdout)",&
     " -output_format, -of <DCD|XMOL> : Output format for dump",&
     " -pbc                           : Apply the periodic boundary conditions (PBC) to each frame before it is dumped",&
     "                                  (origin at lower left)",&
     " -pbc0                          : Apply the periodic boundary conditions (PBC) to each frame before it is dumped",&
     "                                  (origin at box centre)",&
     " -stride <int>                  : Stride for frame dump (allows to skip frames, e.g. by dumping each 10th frame)",&
     " -trace_atoms <real>            : Print the atoms which left the simulation box given a threshold value in scaled units",&
     " -vel2cord, -v2c                : Dump a VELocity DCD file as COoRDinate DCD file (hack which allows a digest by VMD)",&
     " -xyz_file, -xyz <file_name>    : Name of a reference XYZ file in XMOL format that provides the atomic labels",&
     ""

    WRITE (UNIT=*,FMT="(T2,A)")&
     "Usage examples:",&
     "",&
     " dumpdcd <optional flags> <DCD file(s)>",&
     "",&
     "Specific usage examples:",&
     "",&
     " dumpdcd UO2-2x2x2-pos-1.dcd (without atomic labels from XYZ file)",&
     " dumpdcd -xyz UO2-2x2x2.xyz UO2-2x2x2-pos-1.dcd (single DCD file)",&
     " dumpdcd -xyz UO2-2x2x2.xyz UO2-2x2x2-pos-1.dcd UO2-2x2x2-pos-2.dcd ... (multiple DCD files are dumped consecutively)",&
     " dumpdcd -info -xyz UO2-2x2x2.xyz UO2-2x2x2-pos-1.dcd UO2-2x2x2-pos-2.dcd (print additional information)",&
     " dumpdcd -debug -xyz UO2-2x2x2.xyz UO2-2x2x2-pos-1.dcd UO2-2x2x2-pos-2.dcd (print debug information)",&
     " dumpdcd -ekin -d -xyz UO2-2x2x2.xyz UO2-2x2x2-vel-1.dcd (print the ""temperature"" of each atom)",&
     " dumpdcd -ekin -xyz UO2-2x2x2.xyz UO2-2x2x2-vel-1.dcd (print just the temperature of each atom)",&
     " dumpdcd -first_frame 5 -last_frame 10 UO2-2x2x2-pos-1.dcd (just dump frame 5 to 10, ie. 6 frames in total)",&
     " dumpdcd -o outfile.xyz UO2-2x2x2-pos-1.dcd (write output to the file ""outfile.xyz"" instead of stdout)",&
     " dumpdcd -o test.xyz -output_format xmol -xyz ref.xyz -first 10 -last 10 test.dcd (dump 10th frame in XMOL format)",&
     " dumpdcd -of dcd -ff 10 -lf 20 test.dcd (dump the frames 10 to 20 in DCD format to the default output file output.dcd)",&
     " dumpdcd -o part.dcd -of dcd -ff 1 -lf 3 test.dcd (dump the frames 1 to 3 in DCD format to the output file part.dcd)",&
     " dumpdcd -o part.dcd -of dcd -first 10 -lf 100 -stride 10 test.dcd (dump the frames 10, 20, ..., 100 to the file part.dcd)",&
     " dumpdcd -output new.dcd -output_format dcd -pbc old.dcd (dump all frames applying PBC to the output file new.dcd)",&
     " dumpdcd -o new.dcd -of dcd -pbc -trace_atoms 0.02 old.dcd (all atoms more than 2% out of the box are listed)",&
     " dumpdcd -o new.dcd -e out_of_box.log -of dcd -pbc -trace_atoms 0.1 old.dcd (atoms more than 10% out of the box are listed)",&
     " dumpdcd -o new.dcd -of dcd -vel2cord old.dcd (dump old.dcd as new.dcd and change only the DCD id string from VEL to CORD)",&
     " dumpdcd -i -disp UO2-2x2x2-pos-1.dcd (dump the displacements of all atoms w.r.t. their positions in the first frame)",&
     " dumpdcd -i -of dcd -disp UO2-2x2x2-pos-1.dcd (dump the atomic displacements as x-coordinates of a DCD CORD file)",&
     " dumpdcd -i -of dcd -ekin -v2c -xyz UO2-2x2x2.xyz UO2-2x2x2-vel-1.dcd (dump the atomic temperatures as x-coordinates of a ",&
     "                                                                       DCD CORD file -> hack for VMD)",&
     ""

    WRITE (UNIT=*,FMT="(T2,A)")&
     "Notes:",&
     "",&
     " - For -ekin a XYZ file is required to obtain the atomic labels",&
     " - The -info and the -debug flags provide a more detailed output which is especially handy for tracing problems",&
     " - The output in DCD format is in binary format",&
     " - The input coordinates should be in Angstrom. Velocities and forces are expected to be in atomic units",&
     ""

  END SUBROUTINE print_help

  SUBROUTINE uppercase(string)
    ! Convert all letters in a string to uppercase

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(INOUT) :: string

    INTEGER                         :: i,iascii

    DO i=1,LEN_TRIM(string)
       iascii = ICHAR(string(i:i))
       IF ((iascii >= 97).AND.(iascii <= 122)) THEN
         string(i:i) = CHAR(iascii - 32)
       END IF
    END DO

  END SUBROUTINE uppercase

  SUBROUTINE write_out_of_box_atoms(atomic_label,r,s,eps_out_of_box,h)
    ! Print a list of all atoms which have left the simulation box

    IMPLICIT NONE

    CHARACTER(LEN=5), DIMENSION(:), INTENT(IN) :: atomic_label
    REAL(KIND=sp), DIMENSION(:,:), INTENT(IN)  :: r
    REAL(KIND=dp), DIMENSION(:,:), INTENT(IN)  :: s
    REAL(KIND=dp), INTENT(IN)                  :: eps_out_of_box
    REAL(KIND=dp), DIMENSION(3,3), INTENT(IN)  :: h

    REAL(KIND=dp), DIMENSION(3)                :: ds,dr
    REAL(KIND=dp)                              :: rl,s_max,s_min,sl
    INTEGER                                    :: i,iatom,natom,ncount

    ! Quick return, if no action is requested
    IF (eps_out_of_box <= 0.0_dp) RETURN

    s_max = 1.0_dp + eps_out_of_box
    s_min = -eps_out_of_box
    natom = SIZE(s,1)
    ncount = 0
    DO iatom=1,natom
      IF (ANY(s(iatom,:) < s_min).OR.&
          ANY(s(iatom,:) > s_max)) THEN
        ncount = ncount + 1
        IF (ncount == 1) THEN
          WRITE (UNIT=error_unit,FMT="(A)")&
           "#",&
           "# Atoms out of box:",&
           "# Atom index label              x              y              z           |dr|           |ds|"
        END IF
        ds(:) = s(iatom,:)
        DO i=1,3
          IF (s(iatom,i)  < 0.0_dp) ds(i) = 0.0_dp
          IF (s(iatom,i) >= 1.0_dp) ds(i) = 1.0_dp
        END DO
        ds(:) = s(iatom,:) - ds(:)
        sl = SQRT(ds(1)**2 + ds(2)**2 + ds(3)**2)
        dr(:) = MATMUL(h(:,:),ds(:))
        rl = SQRT(dr(1)**2 + dr(2)**2 + dr(3)**2)
        WRITE (UNIT=error_unit,FMT="(A,I10,1X,A5,5(1X,F14.6))")&
         "# ",iatom,ADJUSTR(atomic_label(iatom)),r(iatom,:),rl,sl
      END IF
    END DO
    WRITE (UNIT=error_unit,FMT="(A,I0,A)") "# ",ncount," atom(s) out of box"

  END SUBROUTINE write_out_of_box_atoms

END PROGRAM dumpdcd
