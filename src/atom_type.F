MODULE atom_type

  USE kinds, ONLY : dbl
  USE radial_grids, ONLY : radial_grid_type, init_radial_grid
  USE global_types, ONLY : global_environment_type
  USE termination, ONLY : stop_program
  USE string_utilities, ONLY : uppercase

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: atom_info, atom_init, atom_destroy
  
  TYPE atom_info

     ! the ATOMic properties
     INTEGER :: z                          ! atomic number
     INTEGER, DIMENSION(1:55) :: elec_conf ! configuration of electrons in orbits

     ! the ATOMic results
     REAL(DBL) :: ground_state_energy               
     REAL(DBL), DIMENSION(55) :: eigenvalues        ! energy eigenvalues of the orbitals
     REAL(DBL), DIMENSION(:,:), POINTER :: orbitals ! the electronic orbitals

     ! program parameters
     TYPE (radial_grid_type) :: rg         ! the grid to be used for computations

     CHARACTER(8) :: method                ! name of the method to be used
     REAL(DBL), DIMENSION(5) :: mpara      ! parameters for the method (specific)

     REAL(DBL) :: precision                ! precision to be achieved
     INTEGER :: maxiter                    ! maximum number of iterations
     INTEGER :: maxshoot                   ! maximum number of shooting tries
     REAL(DBL) :: dma                      ! density mixing angle
     REAL(DBL) :: ncv                      ! node correction value
     REAL(DBL) :: lddrsv                   ! ldd rescale factor

     CHARACTER(255) :: wf_logfile
     LOGICAL :: log_wf

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

    USE functionals, ONLY : init_functionals
    USE mathconstants, ONLY : pi
    USE parser, ONLY : p_error, parser_end, parser_init, get_next, &
         test_next, read_line, stop_parser, search_label
    USE termination, ONLY : stop_memory

    IMPLICIT NONE

    ! arguments
    TYPE(atom_info), INTENT(OUT) :: this                 ! atom to be initialized
    TYPE(global_environment_type), INTENT(IN) :: globenv ! containing program environment options
    LOGICAL, INTENT(IN), OPTIONAL :: print

    ! variables
    INTEGER :: err                     ! error code
    INTEGER :: l                       ! angular quantum number
    INTEGER :: n                       ! main quantum number
    INTEGER :: iw                      ! environment output device
    INTEGER :: len                     ! length of the fields in file
    INTEGER :: index                   ! a simple index
    CHARACTER(255) :: string           ! 
    CHARACTER(20) :: label             ! label of the input file section to be read
    INTEGER, DIMENSION(0:9) :: nel     ! number of electrons specified in file
                                       
    CHARACTER(8) :: rgtype             ! variables for radial grid
    INTEGER :: ngp                     !
    INTEGER :: orbit
                                       
    REAL(DBL) :: rmax                  !
    REAL(DBL), DIMENSION(0:5) :: b     ! 
    REAL(DBL) :: z                     !
                                       
    LOGICAL, DIMENSION(0:5) :: b_set   ! were supplied in the file
    LOGICAL :: nel_set, mpara_set(5), z_set 
    LOGICAL :: pflag
    
    pflag = .FALSE.
    IF (PRESENT(print)) pflag = print
    label = "&ATOMINFO"
    iw = globenv%scr

    b_set = .FALSE.
    nel_set = .FALSE.
    mpara_set = .FALSE.
    z_set = .FALSE.
    
    ! some default values
    nel = 0
    this%elec_conf = 0
    this%mpara = 0.0_dbl
    this%maxiter = 40
    this%maxshoot = 1500
    this%precision = 1.E-9_dbl
    b = 1.0_dbl; b(1) = 0.08_dbl
    ngp = 100
    rmax = 10._dbl

    this%dma = 0.5_dbl
    this%ncv = 0.1_dbl
    this%lddrsv = 2.0_dbl

    this%log_wf = .FALSE.

    ! open file for reading
    CALL parser_init( globenv%input_file_name, globenv )
    CALL search_label( label, err )
    IF (err /= 0 ) THEN
       IF (globenv%ionode) WRITE (iw,'(A)') "No input section &ATOMINFO found"
    ELSE
       CALL read_line()
       DO WHILE (test_next() /= 'X')

          len = 8; CALL get_next(string, len)
          CALL uppercase(string)

          SELECT CASE (string)

          CASE DEFAULT
             CALL p_error()
             CALL stop_parser("atom_type/atom_init", "unkown option")

          CASE ('Z')
             IF (test_next() == 'N') THEN 
                CALL get_next ( this%z ); z_set = .true.
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option Z needs numeric argument")
             END IF

          CASE ('DMA')
             IF (test_next() == 'N') THEN 
                CALL get_next ( this%dma );
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option DMA needs numeric argument")
             END IF

          CASE ('NCV')
             IF (test_next() == 'N') THEN 
                CALL get_next ( this%ncv );
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option NCV needs numeric argument")
             END IF

          CASE ('LDDRSV')
             IF (test_next() == 'N') THEN 
                CALL get_next ( this%lddrsv );
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option LDDRSV needs numeric argument")
             END IF

          CASE ('NOELEC')
             l = 0
             DO
                IF (test_next() /= 'N') EXIT
                CALL get_next ( nel(l) ); nel_set = .true.
                l = l+1
             END DO
             IF (.NOT. nel_set ) CALL stop_parser("atom_type/atom_init", &
                  "option NOELEC needs numeric argument")

          CASE ('PRECI')
             IF (test_next() == 'N') THEN
                CALL get_next ( this%precision )
             ELSE
                CALL stop_parser("atom_type/atom_init", &
                     "option PRECI needs numeric argument")
             END IF
             
          CASE ('MAXITER')
             IF (test_next() == 'N') THEN 
                CALL get_next ( this%maxiter )
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option MAXITER needs numeric argument")
             END IF

          CASE ('MAXSHOOT')
             IF (test_next() == 'N') THEN 
                CALL get_next ( this%maxshoot )
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option MAXSHOOT needs numeric argument")
             END IF

          CASE ('GRIDTYPE')
             IF (test_next() == 'C') THEN
                len = 8; CALL get_next(rgtype, len)
                CALL uppercase(rgtype)
             ELSE
                CALL stop_parser("atom_type/atom_init", &
                     "option GRIDTYPE needs string argument")
             END IF
             
          CASE ('GRIDN')
             IF (test_next() == 'N') THEN 
                CALL get_next ( ngp )
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option GRIDN needs numeric argument")
             END IF

          CASE ('GRIDRMAX')
             IF (test_next() == 'N') THEN 
                CALL get_next ( rmax )
             ELSE 
                CALL stop_parser("atom_type/atom_init", &
                     "option GRIDRMAX needs numeric argument")
             END IF

          CASE ('GRIDB')
             index = 0
             DO
                IF (test_next() /= 'N') EXIT
                IF (index > 5) &
                     CALL stop_parser("atom_type/atom_init", &
                     "option GRIDB can only have up to 6 arguments")
                CALL get_next ( b(index) )
                index = index + 1
             ENDDO

          CASE ('WFLOGF')
             IF (test_next() == 'C') THEN
                len = 255; CALL get_next ( this%wf_logfile, len )
                this%log_wf = .TRUE.
             ELSE
                CALL stop_parser("atom_type/atom_init", &
                     "option WFLOGF needs a filename as an argument")
             END IF

          CASE ('METHOD')
             IF (test_next() == 'C') THEN
                len = 8; CALL get_next(this%method, len)
                CALL uppercase(this%method)
             ELSE
                CALL stop_parser("atom_type/atom_init", &
                     "option METHOD needs string argument")
             END IF

          CASE ('MPARA')
             index = 1
             DO
                IF (test_next() /= 'N') EXIT
                IF (index > 5) &
                     CALL stop_parser("atom_type/atom_init", &
                     "option MPARA can only have up to 5 arguments")
                CALL get_next ( this%mpara(index) ); mpara_set(index) = .TRUE.
                index = index + 1
             END DO


          END SELECT
          CALL read_line()
       END DO
    END IF

    ! close the file
    CALL parser_end()

    IF (pflag) THEN
       WRITE (6,'(1X,A,T47,I3)') &
            "Option Z read with value ", this%z
       WRITE (6,'(1X,A,10I3)') &
            "Option NOELEC read with values ", nel

       WRITE (6,'(1X,A,T42,A8)') &
            "Option MEHTOD read with value ", TRIM(this%method)
       WRITE (6,'(1X,A,5ES10.2)') &
            "Option MPARA read with values  ", this%mpara

       WRITE (6,'(1X,A,T40,ES10.3)') &
            "Option PRECI read with value ", this%precision

       WRITE (6,'(1X,A,T47,I3)') &
            "Option MAXITER read with value ", this%maxiter
       WRITE (6,'(1X,A,T45,I5)') &
            "Option MAXSHOOT read with value ", this%maxshoot
       WRITE (6,'(1X,A,T46,F4.2)') &
            "Option DMA read with value ", this%dma
       WRITE (6,'(1X,A,T45,F5.1)') &
            "Option NCV read with value ", this%ncv
       WRITE (6,'(1X,A,T45,F5.1)') &
            "Option LDDRSV read with value ", this%lddrsv

       WRITE (6,'(1X,A,T40,A10)') &
            "Option GRIDTYPE read with value ", TRIM(rgtype)
       WRITE (6,'(1X,A,T46,I4)') &
            "Option GRIDN read with value ", ngp
       WRITE (6,'(1X,A,T45,F5.2)') &
            "Option GRIDRMAX read with value ", rmax
       WRITE (6,'(1X,A,3ES10.3)') &
            "Option GRIDB read with values ", b(0:2)
       
       WRITE (6,*)
    END IF


    ! check the supplied values
    IF (.NOT. z_set .OR. (this%z < 1)) THEN
       CALL stop_program("atom_type/atom_init", "no or bad value provided for z")
    END IF
    IF (.NOT. nel_set .OR. SUM(nel) < 1) THEN
       CALL stop_program("atom_type/atom_init", &
            "no or bad value provided for the number of electrons")
    END IF
    IF (ngp < 10) THEN
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
       CALL init_radial_grid ( this%rg, "LINEAR", ngp, b(0), b(1:5), rmax )
    CASE ('LOG')
       CALL init_radial_grid ( this%rg, "LOG", ngp, b(0), b(1:5), rmax )
    CASE ('POWER')
       CALL init_radial_grid ( this%rg, "POWER", ngp, b(0), b(1:5), rmax )
    CASE ('RATIONAL')
       CALL init_radial_grid ( this%rg, "RATIONAL", ngp, b(0), b(1:5), rmax )
    CASE ('ALG')
       CALL init_radial_grid ( this%rg, "ALGEBRAIC SQUARE ROOT", ngp, b(0), b(1:5), rmax )
    END SELECT

    ! allocate memory for the electron orbitals
    IF (ASSOCIATED(this%orbitals)) DEALLOCATE (this%orbitals, STAT=err)
    ALLOCATE (this%orbitals(this%rg%n,0:55), STAT=err)
    IF (err /=0) CALL &
         stop_memory("atom_type/atom_init","atom_info%orbitals",55*this%rg%n)
    
    ! further necessary initialisations

    ! method specific
    SELECT CASE (this%method)
       
    CASE DEFAULT
       CALL stop_program("atom_type/atom_init", "Method not implemented")
       
    CASE ('SIMPLE')
       !! nothing to do
       
    CASE ('XA')
       IF (.NOT.mpara_set(1)) this%mpara(1) = 0.6666666666667_dbl
       
       !! the constant factor 3/2 alpha (3/pi)^1/3 
       !! is stored in parameter 2 
       this%mpara(2) = 1.5_dbl * this%mpara(1) * (3.0_dbl/pi)**(1.0_dbl/3.0_dbl)

    END SELECT

    ! from the read electron numbers compute the electron configuration of the atom
    DO n=1, 10
       DO l=0, n-1
          orbit = (n*(n-1))/2 + l + 1
          this%elec_conf(orbit) = MIN(2*(2*l+1), nel(l) - (n-(l+1))*2*(2*l+1))
          IF (this%elec_conf(orbit) < 1) this%elec_conf(orbit) = 0
       END DO
    END DO

    CALL init_functionals( xalpha = this%mpara(1) )

  END SUBROUTINE atom_init

  !! Destructor
  SUBROUTINE atom_destroy(this)
    IMPLICIT NONE
    TYPE(atom_info), INTENT(INOUT) :: this

    DEALLOCATE (this%orbitals)
  END SUBROUTINE atom_destroy

END MODULE atom_type
