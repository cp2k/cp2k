!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE atom_types

  USE global_types,                    ONLY: GLOBAL_ENVIRONMENT_TYPE
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE parser,                          ONLY: get_next,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE radial_grids,                    ONLY: RADIAL_GRID_TYPE,&
                                             init_radial_grid
  USE radial_util,                     ONLY: radial_normalize
  USE string_utilities,                ONLY: uppercase
  USE xc_xalpha_new,                   only: xalpha_init
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
     REAL(KIND=dp)               :: precision
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
    TYPE(GLOBAL_ENVIRONMENT_TYPE), &
      INTENT(IN)                             :: globenv

! variables

    CALL parser_init(globenv%input_file_name, globenv)
    
    CALL handle_ac(ac, globenv)
    CALL handle_grid(ac, globenv)
    CALL handle_atom(atom, ac, globenv)
    CALL handle_method_specific(ac, globenv)
    
    CALL print_atom(atom)
    CALL print_parameters(ac)

    CALL parser_end()

  END SUBROUTINE init
  
  ! ****************************************************************************
  
  SUBROUTINE handle_atom ( atom, ac, globenv )
    TYPE(ATOM_TYPE), INTENT(OUT)             :: atom
    TYPE(ATOM_CALCULATION_TYPE), INTENT(IN)  :: ac
    TYPE(GLOBAL_ENVIRONMENT_TYPE), &
      INTENT(IN)                             :: globenv

    CHARACTER(8)                             :: label, string
    INTEGER                                  :: avail, err, i, l, len, n, np, &
                                                out, set
    REAL(KIND=dp)                                :: x

! variables

    np = ac%rg%n
    len = 8
    out = globenv%scr
    label = "&ATOM"
    CALL search_label(label, err)
    IF ( err /= 0 ) THEN
       IF ( globenv%ionode ) WRITE ( out,'(A)' ) "No &ATOM section found"
    ELSE
       CALL read_line ()
       DO WHILE ( test_next () /= 'X' )
          
          CALL get_next(string, len)
          CALL uppercase(string)
          
          SELECT CASE ( string )
             
          CASE DEFAULT
             
          CASE ('Z')
             IF (test_next() == 'N') THEN 
                CALL get_next(atom%z)
             ELSE 
                CALL stop_parser("atom_types/init", "Z needs numeric argument")
             END IF
             
          CASE ('N')
             IF (test_next () == 'N') THEN
                CALL get_next(atom%ne)
             ELSE
                CALL stop_parser("atom_types/init", "N needs numeric argument")
             END IF
             
          END SELECT
          CALL read_line ()
       END DO
    END IF
    
!!!
!!! Fill the electron shells of the atom
!!!
    ! find the number of shells
    n = 1
    DO 
       IF ((2*(n*(n+1)*(2*n+1))/6) >= atom%ne) EXIT
       n = n + 1
    END DO
    atom%ns = n
    ALLOCATE(atom%shells(atom%ns))

    avail = atom%ne
    ! loop over shells
    DO n=1, atom%ns
       atom%shells(n)%n = n

       ! find number of different angular momenta in this shell
       l=0
       DO  
          IF ((avail <= 2*(l+1)**2) .OR. (l == n-1)) EXIT
          l = l + 1
       END DO
       atom%shells(n)%no = l+1
       ALLOCATE(atom%shells(n)%orbits(0:atom%shells(n)%no-1))

       ! loop over these momenta
       DO l=0, atom%shells(n)%no-1
          IF (avail > 0) THEN

             set = MIN(avail, 2*(2*l+1))

             atom%shells(n)%no = l+1
             atom%shells(n)%orbits(l)%ne = set
             atom%shells(n)%orbits(l)%n = n
             atom%shells(n)%orbits(l)%l = l

             SELECT CASE (ac%method)
             CASE DEFAULT
                atom%shells(n)%orbits(l)%ev = &
                     -0.5_dp * (REAL(atom%z, KIND=dp)/REAL(n, KIND=dp))**2
             CASE ('SIMPLE')
                atom%shells(n)%orbits(l)%ev = &
                     -0.7_dp * (REAL(atom%z, KIND=dp)/REAL(n, KIND=dp))**2
             CASE ('XA')
                atom%shells(n)%orbits(l)%ev = &
                     -0.5_dp * (REAL(avail, KIND=dp)/REAL(n, KIND=dp))**2
             END SELECT

!!! fill the orbits with (sensible) estimations of eigenvalues and eigenfunctions
             ALLOCATE(atom%shells(n)%orbits(l)%wf(1:np))
             DO i=1, np
                x = 2.25_dp * ac%rg%r(i) / REAL((n*(n-1))/2+l+1,KIND=dp)**2
                atom%shells(n)%orbits(l)%wf(i) = exp(-x)*x
             END DO
             CALL radial_normalize(atom%shells(n)%orbits(l)%wf(1:np), ac%rg)

             avail = avail - set

          ELSE
             EXIT
          END IF
       END DO
    END DO
    
    IF (.NOT.is_valid(atom)) THEN
       WRITE (6,*) "Atom entries invalid"
    END IF
    
  END SUBROUTINE handle_atom
  
  ! ****************************************************************************
  
  SUBROUTINE handle_ac ( ac, globenv )
    TYPE(ATOM_CALCULATION_TYPE), INTENT(OUT) :: ac
    TYPE(GLOBAL_ENVIRONMENT_TYPE), &
      INTENT(IN)                             :: globenv

    CHARACTER(8)                             :: label, string
    INTEGER                                  :: err, len, out

! variables

    len = 8
    out = globenv%scr
    label = "&CALC"
    CALL search_label(label, err)
    IF ( err /= 0 ) THEN
       IF ( globenv%ionode ) WRITE ( out,'(A)' ) "No &CALC section found"
    ELSE
       CALL read_line ()
       DO WHILE ( test_next () /= 'X' )
          
          CALL get_next ( string, len )
          CALL uppercase ( string )
          
          SELECT CASE ( string )
             
          CASE DEFAULT

          CASE ('METHOD')
             IF (test_next() == 'C') THEN
                CALL get_next(ac%method, len)
                CALL uppercase(ac%method)
             ELSE
                CALL stop_parser("atom_types/init", "METHOD needs string argument")
             END IF
             
          CASE ('MAXITER')
             IF (test_next() == 'N') THEN 
                CALL get_next(ac%maxiter)
             ELSE 
                CALL stop_parser("atom_types/init", "MAXITER needs numeric argument")
             END IF

          CASE ('MAXSHOOT')
             IF (test_next() == 'N') THEN 
                CALL get_next(ac%maxshoot)
             ELSE 
                CALL stop_parser("atom_types/init", "MAXSHOOT needs numeric argument")
             END IF

          CASE ('PRECIS')
             IF (test_next() == 'N') THEN 
                CALL get_next(ac%precision)
             ELSE 
                CALL stop_parser("atom_types/init", "PRECIS needs numeric argument")
             END IF

          END SELECT
          CALL read_line ()
       END DO
    END IF

  END SUBROUTINE handle_ac
  
  ! ****************************************************************************
  
  SUBROUTINE handle_grid ( ac, globenv )
    TYPE(ATOM_CALCULATION_TYPE), INTENT(OUT) :: ac
    TYPE(GLOBAL_ENVIRONMENT_TYPE), &
      INTENT(IN)                             :: globenv

    CHARACTER(8)                             :: gridtype, label, string
    INTEGER                                  :: err, len, np, out
    REAL(KIND=dp)                                :: rmax
    REAL(KIND=dp), DIMENSION(0:5)                :: b

! variables

    len = 8
    out = globenv%scr
    label = "&GRID"
    CALL search_label(label, err)
    IF ( err /= 0 ) THEN
       IF ( globenv%ionode ) WRITE ( out,'(A)' ) "No &GRID section found"
    ELSE
       CALL read_line ()
       DO WHILE ( test_next () /= 'X' )
          
          CALL get_next ( string, len )
          CALL uppercase ( string )
          
          SELECT CASE ( string )
             
          CASE DEFAULT

          CASE ('GRIDTYPE')
             IF (test_next() == 'C') THEN
                CALL get_next(gridtype, len)
                CALL uppercase (gridtype)
             ELSE
                CALL stop_parser("atom_types/init", "GRIDTYPE needs string argument" )
             END IF
             
          CASE ('NP')
             IF (test_next () == 'N') THEN 
                CALL get_next(np)
             ELSE 
                CALL stop_parser("atom_types/init", "NP needs numeric argument")
             END IF
             
          CASE ('RMAX')
             IF (test_next() == 'N') THEN 
                CALL get_next(rmax)
             ELSE 
                CALL stop_parser("atom_types/init", "GRIDRMAX needs numeric argument")
             END IF
             
          CASE ('B0')
             IF (test_next() == 'N') THEN 
                CALL get_next(b(0))
             ELSE 
                CALL stop_parser("atom_types/init", "B0 needs numeric argument")
             END IF

          CASE ('B1')
             IF (test_next() == 'N') THEN 
                CALL get_next(b(1))
             ELSE 
                CALL stop_parser("atom_types/init", "B1 needs numeric argument")
             END IF

          CASE ('B2')
             IF (test_next() == 'N') THEN 
                CALL get_next(b(2))
             ELSE 
                CALL stop_parser("atom_types/init", "B2 needs numeric argument")
             END IF

          CASE ('B3')
             IF (test_next() == 'N') THEN 
                CALL get_next(b(3))
             ELSE 
                CALL stop_parser("atom_types/init", "B3 needs numeric argument")
             END IF

          CASE ('B4')
             IF (test_next() == 'N') THEN 
                CALL get_next(b(4))
             ELSE 
                CALL stop_parser("atom_types/init", "B4 needs numeric argument")
             END IF

          CASE ('B5')
             IF (test_next() == 'N') THEN 
                CALL get_next(b(5))
             ELSE 
                CALL stop_parser("atom_types/init", "B5 needs numeric argument")
             END IF

          END SELECT
          CALL read_line ()
       END DO
    END IF
    
!!! now create the radial grid with the complete informations
    SELECT CASE (gridtype)
    CASE DEFAULT
       CALL stop_parser("atom_types/init", "unknown gridtype specified")
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
    TYPE(ATOM_CALCULATION_TYPE), INTENT(INOUT) :: ac
    TYPE(GLOBAL_ENVIRONMENT_TYPE), &
      INTENT(IN)                             :: globenv

    CHARACTER(8)                             :: label, string
    INTEGER                                  :: err, len, out

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
       CALL search_label(label, err)
       IF ( err /= 0 ) THEN
          ! not necessary to specify so can be left out
       ELSE
          CALL read_line ()
          DO WHILE ( test_next () /= 'X' )

             CALL get_next ( string, len )
             CALL uppercase ( string )

             SELECT CASE ( string )

             CASE DEFAULT
                CALL stop_parser("atom_types/init", "unrecognized option: ")

             CASE ('ALPHA')
                IF (test_next() == 'N') THEN
                   CALL get_next(ac%xapara%alpha)
                ELSE
                   CALL stop_parser("atom_types/init", "ALPHA needs numeric argument")
                END IF
                CALL xalpha_init(1.e-10_dp,ac%xapara%alpha)
                ac%xapara%c = 1.5_dp * ac%xapara%alpha * (3.0_dp/pi)**(1.0_dp/3.0_dp)
             END SELECT

             CALL read_line ()
          END DO
       END IF

    CASE DEFAULT
       CALL stop_parser("atom_types/init", "unrecognized method: ")
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
       CALL stop_parser("atom_types/print_parameters","Unkown calculation method "//ac%method)
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
