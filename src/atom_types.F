!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE atom_types

!TODEL_TL  USE parser,                          ONLY: get_next,&
!TODEL_TL                                             parser_end,&
!TODEL_TL                                             parser_init,&
!TODEL_TL                                             read_line,&
!TODEL_TL                                             search_label,&
!TODEL_TL                                             stop_parser,&
!TODEL_TL                                             test_next
  USE global_types,                    ONLY: GLOBAL_ENVIRONMENT_TYPE
  USE kinds,                           ONLY: dp
  USE radial_grids,                    ONLY: RADIAL_GRID_TYPE,&
                                             init_radial_grid

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ATOM_TYPE, &
       ELEC_ORBIT_TYPE, &
       ELEC_SHELL_TYPE, &
       ATOM_CALCULATION_TYPE, &
       XA_PARA_TYPE, &
       init, &
       print_atom

  ! **************************** TYPE DECLARATIONS *****************************

  TYPE ELEC_ORBIT_TYPE
     INTEGER :: ne                              ! number of electrons in this orbit
     INTEGER :: n                               ! main quantum number of this orbit
     INTEGER :: l                               ! angular quantum number of this orbit
     REAL(KIND=dp) :: ev                            ! energy eigenvalue of this orbit
     REAL(KIND=dp), DIMENSION(:), POINTER :: wf     ! the wavefunction of this orbit
  END TYPE ELEC_ORBIT_TYPE
  
  TYPE ELEC_SHELL_TYPE
     INTEGER :: n
     INTEGER :: no
     TYPE(ELEC_ORBIT_TYPE), DIMENSION(:), POINTER :: orbits
  END TYPE ELEC_SHELL_TYPE
  
  TYPE ATOM_TYPE
     INTEGER :: z                                ! number of protons
     INTEGER :: ne                               ! number of electrons
     TYPE(ELEC_SHELL_TYPE), DIMENSION(:), POINTER :: shells
     INTEGER :: ns                               ! number of electronic shells
     REAL(KIND=dp) :: gse                            ! ground state energy
  END TYPE ATOM_TYPE
  
  TYPE XA_PARA_TYPE
     REAL(KIND=dp) :: alpha
     REAL(KIND=dp) :: c                              ! v_xa = - c * rho**(1/3)
  END TYPE XA_PARA_TYPE
  
  TYPE ATOM_CALCULATION_TYPE
     CHARACTER(8)            :: method
     TYPE(XA_PARA_TYPE)      :: xapara
     REAL(KIND=dp)               :: PRECISION
     TYPE(RADIAL_GRID_TYPE)  :: rg
     INTEGER                 :: maxiter
     INTEGER                 :: maxshoot
  END TYPE ATOM_CALCULATION_TYPE
  
  ! ****************************************************************************
  
CONTAINS

  ! ************************** SUBROUTINES AND FUNCTIONS ***********************
  
  FUNCTION is_valid ( atom )
    
    
    TYPE(ATOM_TYPE), INTENT(IN)              :: atom
    LOGICAL                                  :: is_valid

    INTEGER                                  :: l, n, ne

    is_valid = .FALSE.
    ne = 0
    
    DO n=1, atom%ns
       IF (atom%shells(n)%n /= n) RETURN
       IF (atom%shells(n)%no > n) RETURN
       DO l=0, atom%shells(n)%no-1
          ne = ne + atom%shells(n)%orbits(l)%ne
          IF (atom%shells(n)%orbits(l)%n /= n) RETURN
          IF (atom%shells(n)%orbits(l)%l /= l) RETURN
       END DO
    END DO
    IF (ne /= atom%ne) RETURN
    
    ! if we reach this then the atom is valid
    is_valid = .TRUE.
    
  END FUNCTION is_valid
  
  ! ****************************************************************************
  
  SUBROUTINE init ( atom, ac, globenv )
    !
    ! Reads the init file and initializes the atom and atom_calculation
    ! objects with the information therein.
    !
    TYPE(ATOM_TYPE), INTENT(OUT)             :: atom
    TYPE(ATOM_CALCULATION_TYPE), INTENT(OUT) :: ac
    TYPE(GLOBAL_ENVIRONMENT_TYPE), POINTER   :: globenv

! variables
!TODEL_TL    CALL parser_init(globenv%input_file_name, globenv)

    CALL handle_ac(ac, globenv)
    CALL handle_grid(ac, globenv)
    CALL handle_atom(atom, ac, globenv)
    CALL handle_method_specific(ac, globenv)
    
    CALL print_atom(atom)
    CALL print_parameters(ac)

!TODEL_TL    CALL parser_end()

  END SUBROUTINE init
  
  ! ****************************************************************************
  
  SUBROUTINE handle_atom ( atom, ac, globenv )
    TYPE(ATOM_TYPE), INTENT(OUT)             :: atom
    TYPE(ATOM_CALCULATION_TYPE), INTENT(IN)  :: ac
    TYPE(GLOBAL_ENVIRONMENT_TYPE), POINTER   :: globenv

    CHARACTER(8)                             :: label
    INTEGER                                  :: len, n, np, out

! variables

    np = ac%rg%n
    len = 8
    out = globenv%scr
    label = "&ATOM"
!TODEL_TL    CALL search_label(label, err)
!TODEL_TL    IF ( err /= 0 ) THEN
!TODEL_TL       IF ( globenv%ionode ) WRITE ( out,'(A)' ) "No &ATOM section found"
!TODEL_TL    ELSE
!TODEL_TL       CALL read_line ()
!TODEL_TL       DO WHILE ( test_next () /= 'X' )
!TODEL_TL          
!TODEL_TL          CALL get_next(string, len)
!TODEL_TL          CALL uppercase(string)
!TODEL_TL          
!TODEL_TL          SELECT CASE ( string )
!TODEL_TL             
!TODEL_TL          CASE DEFAULT
!TODEL_TL             
!TODEL_TL          CASE ('Z')
!TODEL_TL             IF (test_next() == 'N') THEN 
!TODEL_TL                CALL get_next(atom%z)
!TODEL_TL             ELSE 
!TODEL_TL                CALL stop_parser("atom_types/init", "Z needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL             
!TODEL_TL          CASE ('N')
!TODEL_TL             IF (test_next () == 'N') THEN
!TODEL_TL                CALL get_next(atom%ne)
!TODEL_TL             ELSE
!TODEL_TL                CALL stop_parser("atom_types/init", "N needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL             
!TODEL_TL          END SELECT
!TODEL_TL          CALL read_line ()
!TODEL_TL       END DO
!TODEL_TL    END IF
!TODEL_TL    
!TODEL_TL!!!
!TODEL_TL!!! Fill the electron shells of the atom
!TODEL_TL!!!
!TODEL_TL    ! find the number of shells
!TODEL_TL    n = 1
!TODEL_TL    DO 
!TODEL_TL       IF ((2*(n*(n+1)*(2*n+1))/6) >= atom%ne) EXIT
!TODEL_TL       n = n + 1
!TODEL_TL    END DO
!TODEL_TL    atom%ns = n
!TODEL_TL    ALLOCATE(atom%shells(atom%ns))
!TODEL_TL
!TODEL_TL    avail = atom%ne
!TODEL_TL    ! loop over shells
!TODEL_TL    DO n=1, atom%ns
!TODEL_TL       atom%shells(n)%n = n
!TODEL_TL
!TODEL_TL       ! find number of different angular momenta in this shell
!TODEL_TL       l=0
!TODEL_TL       DO  
!TODEL_TL          IF ((avail <= 2*(l+1)**2) .OR. (l == n-1)) EXIT
!TODEL_TL          l = l + 1
!TODEL_TL       END DO
!TODEL_TL       atom%shells(n)%no = l+1
!TODEL_TL       ALLOCATE(atom%shells(n)%orbits(0:atom%shells(n)%no-1))
!TODEL_TL
!TODEL_TL       ! loop over these momenta
!TODEL_TL       DO l=0, atom%shells(n)%no-1
!TODEL_TL          IF (avail > 0) THEN
!TODEL_TL
!TODEL_TL             set = MIN(avail, 2*(2*l+1))
!TODEL_TL
!TODEL_TL             atom%shells(n)%no = l+1
!TODEL_TL             atom%shells(n)%orbits(l)%ne = set
!TODEL_TL             atom%shells(n)%orbits(l)%n = n
!TODEL_TL             atom%shells(n)%orbits(l)%l = l
!TODEL_TL
!TODEL_TL             SELECT CASE (ac%method)
!TODEL_TL             CASE DEFAULT
!TODEL_TL                atom%shells(n)%orbits(l)%ev = &
!TODEL_TL                     -0.5_dp * (REAL(atom%z, KIND=dp)/REAL(n, KIND=dp))**2
!TODEL_TL             CASE ('SIMPLE')
!TODEL_TL                atom%shells(n)%orbits(l)%ev = &
!TODEL_TL                     -0.7_dp * (REAL(atom%z, KIND=dp)/REAL(n, KIND=dp))**2
!TODEL_TL             CASE ('XA')
!TODEL_TL                atom%shells(n)%orbits(l)%ev = &
!TODEL_TL                     -0.5_dp * (REAL(avail, KIND=dp)/REAL(n, KIND=dp))**2
!TODEL_TL             END SELECT
!TODEL_TL
!TODEL_TL!!! fill the orbits with (sensible) estimations of eigenvalues and eigenfunctions
!TODEL_TL             ALLOCATE(atom%shells(n)%orbits(l)%wf(1:np))
!TODEL_TL             DO i=1, np
!TODEL_TL                x = 2.25_dp * ac%rg%r(i) / REAL((n*(n-1))/2+l+1,KIND=dp)**2
!TODEL_TL                atom%shells(n)%orbits(l)%wf(i) = exp(-x)*x
!TODEL_TL             END DO
!TODEL_TL             CALL radial_normalize(atom%shells(n)%orbits(l)%wf(1:np), ac%rg)
!TODEL_TL
!TODEL_TL             avail = avail - set
!TODEL_TL
!TODEL_TL          ELSE
!TODEL_TL             EXIT
!TODEL_TL          END IF
!TODEL_TL       END DO
!TODEL_TL    END DO
!TODEL_TL    
!TODEL_TL    IF (.NOT.is_valid(atom)) THEN
!TODEL_TL       WRITE (6,*) "Atom entries invalid"
!TODEL_TL    END IF
    
  END SUBROUTINE handle_atom
  
  ! ****************************************************************************
  
  SUBROUTINE handle_ac ( ac, globenv )
    TYPE(ATOM_CALCULATION_TYPE), INTENT(OUT) :: ac
    TYPE(GLOBAL_ENVIRONMENT_TYPE), POINTER   :: globenv

    CHARACTER(8)                             :: label
    INTEGER                                  :: len, out

! variables

    len = 8
    out = globenv%scr
    label = "&CALC"
!TODEL_TL    CALL search_label(label, err)
!TODEL_TL    IF ( err /= 0 ) THEN
!TODEL_TL       IF ( globenv%ionode ) WRITE ( out,'(A)' ) "No &CALC section found"
!TODEL_TL    ELSE
!TODEL_TL       CALL read_line ()
!TODEL_TL       DO WHILE ( test_next () /= 'X' )
!TODEL_TL          
!TODEL_TL          CALL get_next ( string, len )
!TODEL_TL          CALL uppercase ( string )
!TODEL_TL          
!TODEL_TL          SELECT CASE ( string )
!TODEL_TL             
!TODEL_TL          CASE DEFAULT
!TODEL_TL
!TODEL_TL          CASE ('METHOD')
!TODEL_TL             IF (test_next() == 'C') THEN
!TODEL_TL                CALL get_next(ac%method, len)
!TODEL_TL                CALL uppercase(ac%method)
!TODEL_TL             ELSE
!TODEL_TL                CALL stop_parser("atom_types/init", "METHOD needs string argument")
!TODEL_TL             END IF
!TODEL_TL             
!TODEL_TL          CASE ('MAXITER')
!TODEL_TL             IF (test_next() == 'N') THEN 
!TODEL_TL                CALL get_next(ac%maxiter)
!TODEL_TL             ELSE 
!TODEL_TL                CALL stop_parser("atom_types/init", "MAXITER needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL
!TODEL_TL          CASE ('MAXSHOOT')
!TODEL_TL             IF (test_next() == 'N') THEN 
!TODEL_TL                CALL get_next(ac%maxshoot)
!TODEL_TL             ELSE 
!TODEL_TL                CALL stop_parser("atom_types/init", "MAXSHOOT needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL
!TODEL_TL          CASE ('PRECIS')
!TODEL_TL             IF (test_next() == 'N') THEN 
!TODEL_TL                CALL get_next(ac%precision)
!TODEL_TL             ELSE 
!TODEL_TL                CALL stop_parser("atom_types/init", "PRECIS needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL
!TODEL_TL          END SELECT
!TODEL_TL          CALL read_line ()
!TODEL_TL       END DO
!TODEL_TL    END IF

  END SUBROUTINE handle_ac
  
  ! ****************************************************************************
  
  SUBROUTINE handle_grid ( ac, globenv )
    TYPE(ATOM_CALCULATION_TYPE), INTENT(OUT) :: ac
    TYPE(GLOBAL_ENVIRONMENT_TYPE), POINTER   :: globenv

    CHARACTER(8)                             :: gridtype, label
    INTEGER                                  :: len, np, out
    REAL(KIND=dp)                            :: rmax
    REAL(KIND=dp), DIMENSION(0:5)            :: b

! variables

    len = 8
    out = globenv%scr
    label = "&GRID"
!TODEL_TL    CALL search_label(label, err)
!TODEL_TL    IF ( err /= 0 ) THEN
!TODEL_TL       IF ( globenv%ionode ) WRITE ( out,'(A)' ) "No &GRID section found"
!TODEL_TL    ELSE
!TODEL_TL       CALL read_line ()
!TODEL_TL       DO WHILE ( test_next () /= 'X' )
!TODEL_TL          
!TODEL_TL          CALL get_next ( string, len )
!TODEL_TL          CALL uppercase ( string )
!TODEL_TL          
!TODEL_TL          SELECT CASE ( string )
!TODEL_TL             
!TODEL_TL          CASE DEFAULT
!TODEL_TL
!TODEL_TL          CASE ('GRIDTYPE')
!TODEL_TL             IF (test_next() == 'C') THEN
!TODEL_TL                CALL get_next(gridtype, len)
!TODEL_TL                CALL uppercase (gridtype)
!TODEL_TL             ELSE
!TODEL_TL                CALL stop_parser("atom_types/init", "GRIDTYPE needs string argument" )
!TODEL_TL             END IF
!TODEL_TL             
!TODEL_TL          CASE ('NP')
!TODEL_TL             IF (test_next () == 'N') THEN 
!TODEL_TL                CALL get_next(np)
!TODEL_TL             ELSE 
!TODEL_TL                CALL stop_parser("atom_types/init", "NP needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL             
!TODEL_TL          CASE ('RMAX')
!TODEL_TL             IF (test_next() == 'N') THEN 
!TODEL_TL                CALL get_next(rmax)
!TODEL_TL             ELSE 
!TODEL_TL                CALL stop_parser("atom_types/init", "GRIDRMAX needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL             
!TODEL_TL          CASE ('B0')
!TODEL_TL             IF (test_next() == 'N') THEN 
!TODEL_TL                CALL get_next(b(0))
!TODEL_TL             ELSE 
!TODEL_TL                CALL stop_parser("atom_types/init", "B0 needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL
!TODEL_TL          CASE ('B1')
!TODEL_TL             IF (test_next() == 'N') THEN 
!TODEL_TL                CALL get_next(b(1))
!TODEL_TL             ELSE 
!TODEL_TL                CALL stop_parser("atom_types/init", "B1 needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL
!TODEL_TL          CASE ('B2')
!TODEL_TL             IF (test_next() == 'N') THEN 
!TODEL_TL                CALL get_next(b(2))
!TODEL_TL             ELSE 
!TODEL_TL                CALL stop_parser("atom_types/init", "B2 needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL
!TODEL_TL          CASE ('B3')
!TODEL_TL             IF (test_next() == 'N') THEN 
!TODEL_TL                CALL get_next(b(3))
!TODEL_TL             ELSE 
!TODEL_TL                CALL stop_parser("atom_types/init", "B3 needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL
!TODEL_TL          CASE ('B4')
!TODEL_TL             IF (test_next() == 'N') THEN 
!TODEL_TL                CALL get_next(b(4))
!TODEL_TL             ELSE 
!TODEL_TL                CALL stop_parser("atom_types/init", "B4 needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL
!TODEL_TL          CASE ('B5')
!TODEL_TL             IF (test_next() == 'N') THEN 
!TODEL_TL                CALL get_next(b(5))
!TODEL_TL             ELSE 
!TODEL_TL                CALL stop_parser("atom_types/init", "B5 needs numeric argument")
!TODEL_TL             END IF
!TODEL_TL
!TODEL_TL          END SELECT
!TODEL_TL          CALL read_line ()
!TODEL_TL       END DO
!TODEL_TL    END IF
    
!!! now create the radial grid with the complete informations
    SELECT CASE (gridtype)
    CASE DEFAULT
!TODEL_TL       CALL stop_parser("atom_types/init", "unknown gridtype specified")
    CASE ('LINEAR')
       CALL init_radial_grid ( ac%rg, "LINEAR", np, b(0), b(1:5), rmax )
    CASE ('LOG')
       CALL init_radial_grid ( ac%rg, "LOG", np, b(0), b(1:5), rmax )
    CASE ('POWER')
       CALL init_radial_grid ( ac%rg, "POWER", np, b(0), b(1:5), rmax )
    CASE ('RATIONAL')
       CALL init_radial_grid ( ac%rg, "RATIONAL", np, b(0), b(1:5), rmax )
    CASE ('ALG')
       CALL init_radial_grid ( ac%rg, "ALGEBRAIC SQUARE ROOT", np, b(0), b(1:5), rmax )
    END SELECT

  END SUBROUTINE handle_grid
  
!!! ****************************************************************************
  
  SUBROUTINE handle_method_specific ( ac, globenv )
    TYPE(ATOM_CALCULATION_TYPE), &
      INTENT(INOUT)                          :: ac
    TYPE(GLOBAL_ENVIRONMENT_TYPE), POINTER   :: globenv

    CHARACTER(8)                             :: label
    INTEGER                                  :: len, out

! variables

    len = 8
    out = globenv%scr
    SELECT CASE (ac%method)

    CASE ('SIMPLE')

       ! nothing to do

    CASE ('XA')

       ! defaults
       ac%xapara%alpha = 0.667_dp

       label = "&XA"
!TODEL_TL       CALL search_label(label, err)
!TODEL_TL       IF ( err /= 0 ) THEN
!TODEL_TL          ! not necessary to specify so can be left out
!TODEL_TL       ELSE
!TODEL_TL          CALL read_line ()
!TODEL_TL          DO WHILE ( test_next () /= 'X' )
!TODEL_TL
!TODEL_TL             CALL get_next ( string, len )
!TODEL_TL             CALL uppercase ( string )
!TODEL_TL
!TODEL_TL             SELECT CASE ( string )
!TODEL_TL
!TODEL_TL             CASE DEFAULT
!TODEL_TL                CALL stop_parser("atom_types/init", "unrecognized option: ")
!TODEL_TL
!TODEL_TL             CASE ('ALPHA')
!TODEL_TL                IF (test_next() == 'N') THEN
!TODEL_TL                   CALL get_next(ac%xapara%alpha)
!TODEL_TL                ELSE
!TODEL_TL                   CALL stop_parser("atom_types/init", "ALPHA needs numeric argument")
!TODEL_TL                END IF
!TODEL_TL                CALL xalpha_init(1.e-10_dp,ac%xapara%alpha)
!TODEL_TL                ac%xapara%c = 1.5_dp * ac%xapara%alpha * (3.0_dp/pi)**(1.0_dp/3.0_dp)
!TODEL_TL             END SELECT
!TODEL_TL
!TODEL_TL             CALL read_line ()
!TODEL_TL          END DO
!TODEL_TL       END IF

    CASE DEFAULT
!TODEL_TL       CALL stop_parser("atom_types/init", "unrecognized method: ")
    END SELECT

  END SUBROUTINE handle_method_specific

!!! ****************************************************************************

  SUBROUTINE print_atom(atom)
    TYPE(ATOM_TYPE), INTENT(IN)              :: atom

    CHARACTER(11)                            :: e(0:7), t(0:7)
    INTEGER                                  :: l, n

    WRITE (6,*)
    WRITE (6,'(1X,A)') "============================================================="
    WRITE (6,'(1X,A)') "STATE of the atom"
    WRITE (6,'(1X,A)') "                                                             "
    WRITE (6,'(1X,A,I3)') "Z = ", atom%z
    WRITE (6,'(1X,A,I3)') "N = ", atom%ne
    WRITE (6,*)
    DO n=1, atom%ns
       t = "          "
       e = "          "
       DO l=0, MIN(atom%shells(n)%no-1,6)
          WRITE (t(l),'(A,I2,A,I2)') "      ", l+1, ":", atom%shells(n)%orbits(l)%ne
          WRITE (e(l),'(A,F10.5)') " ", atom%shells(n)%orbits(l)%ev
       END DO
       WRITE (6,'(7A11)') t(0),t(1),t(2),t(3),t(4),t(5),t(6)
       WRITE (6,'(7A11)') e(0),e(1),e(2),e(3),e(4),e(5),e(6)
    END DO
    WRITE (6,'(1X,A)') "============================================================="
    WRITE (6,*)
  END SUBROUTINE print_atom

!!! ****************************************************************************

  SUBROUTINE print_parameters(ac)

    
    TYPE(ATOM_CALCULATION_TYPE), INTENT(IN)  :: ac

    WRITE (6,*)
    WRITE (6,'(1X,A)') "============================================================="
    WRITE (6,'(1X,A)') "ATOM CALCULATION PARAMETERS"
    WRITE (6,*)
    WRITE (6,'(1X,A,A)') "method : ", ac%method
    SELECT CASE (ac%method)
    CASE DEFAULT
!TODEL_TL       CALL stop_parser("atom_types/print_parameters","Unkown calculation method "//ac%method)
    CASE ('SIMPLE')
       ! nothing to do
    CASE ('XA')
       WRITE (6,'(1X,A,F12.7)')  "-> alpha : ", ac%xapara%alpha
    END SELECT
    WRITE (6,'(1X,A)') "============================================================="
    WRITE (6,*)

  END SUBROUTINE print_parameters

!!! ****************************************************************************

END MODULE atom_types
