MODULE atom_type

  USE kinds, ONLY : dbl
  USE radial_grids
  USE global_types, ONLY : global_environment_type
  USE string_utilities, ONLY : uppercase
  USE termination, ONLY : stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: atom_info, atom_init
  
  TYPE atom_info
     ! the properties
     INTEGER :: z                          ! atomic number
     INTEGER, DIMENSION(0:9) :: nel        ! number of electrons
     CHARACTER(8) :: method                ! name of the method to be used
     REAL(DBL) :: precision                ! precision to be achieved
     INTEGER :: maxiter                    ! maximum number of iterations
     TYPE (radial_grid_type) :: rg         ! the grid to be used for computations
     ! the results
     REAL(DBL) :: gse                               ! ground state energy
     REAL(DBL), DIMENSION(:), POINTER :: oes        ! energy of the orbitals
     REAL(DBL), DIMENSION(:,:), POINTER :: orbitals ! the electronic orbitals
  END TYPE atom_info

  CONTAINS

  !-----------------------------------------------------------------------------!
  !!****** atom_type/atom_init [1.0] *
  !!
  !!   FUNCTION
  !!     Fetches information about the atom to be calcultated and fills an 
  !!     atom record with it. Reads from the file input in the working directory.
  !!
  !!   HISTORY
  !!     TCH, 10-OCT-2001, began work
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE atom_init(this, globenv, print)

    USE parser

    IMPLICIT NONE

    ! arguments
    TYPE(atom_info), INTENT(OUT) :: this                 ! atom to be initialized
    TYPE(global_environment_type), INTENT(IN) :: globenv ! containing program environment options
    LOGICAL, INTENT(IN), OPTIONAL :: print

    ! variables
    INTEGER :: err                                       ! error code
    INTEGER :: l                                         ! counter
    INTEGER :: iw                                        ! environment output device
    INTEGER :: len                                       ! length of the fields in file
    CHARACTER(255) :: string                             ! 
    CHARACTER(9) :: label                                ! label of the input file section to be read

    CHARACTER(8) :: rgtype                               ! variables for radial grid
    INTEGER :: n                                         !
    REAL(DBL) :: rmax                                    !
    REAL(DBL), DIMENSION(0:5) :: b                       ! 
    REAL(DBL) :: z                                       !

    LOGICAL :: z_set, nel_set                            ! flags to check if parameters
    LOGICAL, DIMENSION(0:5) :: b_set                       ! were supplied in the file
    LOGICAL :: pflag
    
    pflag = print
    label = '&ATOMINFO'
    iw = globenv%scr

    z_set = .false.
    nel_set = .false.
    b_set = .false.

    ! some default values
    this%nel = 0
    this%maxiter = 10
    this%precision = 1.E-9_dbl
    b = 1.0_dbl; b(1) = 0.04_dbl
    n = 100
    rmax = 10._dbl

    ! open file for reading
    CALL parser_init(globenv%input_file_name, label, err, globenv)
    IF (err /= 0 ) THEN
       IF (globenv%ionode) WRITE (iw,'(a)') 'No input section &ATOMINFO found '
    ELSE
       CALL read_line()
       DO WHILE (test_next()/='X')

          len = 8; CALL cfield(string, len)
          CALL uppercase(string)

          SELECT CASE (string)
          CASE DEFAULT
             CALL p_error()
             CALL stop_parser("atom_type/atom_init", "unkown option")
          CASE ('Z')
             IF (test_next() == 'N') THEN 
                this%z = get_int(); z_set = .true.
                IF (pflag) WRITE (6,'(1X,(A),T47,I3)') "Option Z read with value ", this%z
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option Z needs numeric argument")
             END IF

          CASE ('NOELEC')
             l = 0
             DO
                IF (test_next() /= 'N') EXIT
                this%nel(l) = get_int(); nel_set = .true.
                IF (pflag) &
                     WRITE (6,'(1X,(A),I1,(A),T47,I3)') &
                     "Option NOELEC read with value (l=", l, ") ", this%nel(l)
                l = l+1
             END DO
             IF (.NOT. nel_set) CALL stop_parser("atom_type/atom_init", &
                  "option NOELEC needs numeric argument")

          CASE ('METHOD')
             IF (test_next() == 'C') THEN
                len = 8; CALL cfield(string, len)
                CALL uppercase(string)
                this%method = string
                IF (pflag) WRITE (6,'(1X,(A),T42,A8)') "Option MEHTOD read with value ", TRIM(this%method)
             ELSE
                CALL stop_parser("atom_type/atom_init", &
                     "option METHOD needs string argument")
             END IF

          CASE ('PRECI')
             IF (test_next() == 'N') THEN
                this%precision =  DBLE(get_real())
                IF (pflag) WRITE (6,'(1X,(A),T40,ES10.3)') "Option PRECI read with value ", this%precision
             ELSE
                CALL stop_parser("atom_type/atom_init", &
                     "option PRECI needs numeric argument")
             END IF
             
          CASE ('MAXITER')
             IF (test_next() == 'N') THEN 
                this%maxiter = get_int()
                IF (pflag) WRITE (6,'(1X,(A),T47,I3)') "Option MAXITER read with value ", this%maxiter
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option MAXITER needs numeric argument")
             END IF

          CASE ('GRIDTYPE')
             IF (test_next() == 'C') THEN
                len = 8; CALL cfield(rgtype, len)
                CALL uppercase(rgtype)
                IF (pflag) WRITE (6,'(1X,(A),T40,A10)') "Option GRIDTYPE read with value ", TRIM(rgtype)
             ELSE
                CALL stop_parser("atom_type/atom_init", &
                     "option GRIDTYPE needs string argument")
             END IF
             
          CASE ('GRIDN')
             IF (test_next() == 'N') THEN 
                n = get_int()
                IF (pflag) WRITE (6,'(1X,(A),T46,I4)') "Option GRIDN read with value ", n
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option GRIDN needs numeric argument")
             END IF

          CASE ('GRIDRMAX')
             IF (test_next() == 'N') THEN 
                rmax = get_real()
                IF (pflag) WRITE (6,'(1X,(A),T40,ES10.3)') "Option GRIDRMAX read with value ", rmax
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option GRIDRMAX needs numeric argument")
             END IF

          CASE ('GRIDB0')
             IF (test_next() == 'N') THEN 
                b(0) = get_real(); b_set(0) = .TRUE.
                IF (pflag) WRITE (6,'(1X,(A),T40,ES10.3)') "Option GRIDB0 read with value ", b(0)
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option GRIDB0 needs numeric argument")
             END IF

          CASE ('GRIDB1')
             IF (test_next() == 'N') THEN 
                b(1) = get_real(); b_set(1) = .TRUE.
                IF (pflag) WRITE (6,'(1X,(A),T40,ES10.3)') "Option GRIDB1 read with value ", b(1)
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option GRIDB1 needs numeric argument")
             END IF

          CASE ('GRIDB2')
             IF (test_next() == 'N') THEN 
                b(2) = get_real(); b_set(2) = .TRUE.
                IF (pflag) WRITE (6,'(1X,(A),T40,ES10.3)') "Option GRIDB2 read with value ", b(2)
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option GRIDB2 needs numeric argument")
             END IF

          CASE ('GRIDB3')
             IF (test_next() == 'N') THEN 
                b(3) = get_real(); b_set(3) = .TRUE.
                IF (pflag) WRITE (6,'(1X,(A),T40,ES10.3)') "Option GRIDB3 read with value ", b(3)
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option GRIDB3 needs numeric argument")
             END IF

          CASE ('GRIDB4')
             IF (test_next() == 'N') THEN 
                b(4) = get_real(); b_set(4) = .TRUE.
                IF (pflag) WRITE (6,'(1X,(A),T40,ES10.3)') "Option GRIDB4 read with value ", b(4)
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option GRIDB4 needs numeric argument")
             END IF

          CASE ('GRIDB5')
             IF (test_next() == 'N') THEN 
                b(5) = get_real(); b_set(5) = .TRUE.
                IF (pflag) WRITE (6,'(1X,(A),T40,ES10.3)') "Option GRIDB5 read with value ", b(5)
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option GRIDB5 needs numeric argument")
             END IF

          END SELECT
          CALL read_line()
       END DO
    END IF

    ! close the file
    CALL parser_end()

    ! check the supplied values
    IF (.NOT. z_set .OR. (this%z < 1)) THEN
       CALL stop_program("atom_type/atom_init", "no or bad value provided for z")
    END IF
    IF (.NOT. nel_set .OR. SUM(this%nel) < 1) THEN
       CALL stop_program("atom_type/atom_init", &
            "no or bad value provided for the number of electrons")
    END IF
    IF (n < 10) THEN
       CALL stop_program("atom_type/atom_init", &
            "bad GRIDN value specified, must be at least 10")
    END IF
    IF (this%precision < 0._dbl) THEN
       CALL stop_program("atom_type/atom_init", &
            "bad precision value prodivded, must be greater than 0")
    END IF
    IF (rmax < 0._dbl) THEN
       CALL stop_program("atom_type/atom_init", &
            "bad rmax value prodivded, must be greater than 0")
    END IF

    ! now create the radial grid with the complete informations
    SELECT CASE (rgtype)
    CASE DEFAULT
       CALL stop_parser("atomic_code/atom_init", &
            "unknown gridtype specified")
    CASE ('LINEAR')
       CALL init_radial_grid ( this%rg, "LINEAR", n, b(0), b(1:5), rmax )
    CASE ('LOG')
       CALL init_radial_grid ( this%rg, "LOG", n, b(0), b(1:5), rmax )
    CASE ('POWER')
       CALL init_radial_grid ( this%rg, "POWER", n, b(0), b(1:5), rmax )
    CASE ('RATIONAL')
       CALL init_radial_grid ( this%rg, "RATIONAL", n, b(0), b(1:5), rmax )
    CASE ('ALG')
       CALL init_radial_grid ( this%rg, "ALGEBRAIC SQUARE ROOT", n, b(0), b(1:5), rmax )
    END SELECT

    IF (pflag) WRITE (6,*)

    IF (this%rg%n /= n) CALL stop_program("atom_type/atom_init","internal error")

  END SUBROUTINE atom_init

END MODULE atom_type
