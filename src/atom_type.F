!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/atom_type *
!!
!!   NAME
!!     atom_type
!!
!!   FUNCTION
!!     This module contains the >> AM_ATOM_INFO << structure, which stores
!!     information about an atom. The structure comes with accessor methods and
!!     follows Fawzis PAO method (at least I think so :). The structures content
!!     is public because other modules in cp2k do not yet support data hiding
!!     (radial_grids) mainly.
!!
!!   SYNOPSIS
!!     am_dealloc_ref(atom)
!!     am_init(atom, globenv)
!!     am_get_elec_conf(atom)
!!     am_get_elec_conf(atom, index)
!!
!!   MODIFICATION HISTORY
!!     11-FEB-2002, TCH, working version
!!
!!   SEE ALSO
!!     atom_module
!!
!!****
MODULE atom_type

  USE kinds, ONLY : dbl
  USE radial_grids, ONLY : radial_grid_type, init_radial_grid
  USE global_types, ONLY : global_environment_type
  USE mathconstants, ONLY : pi
  USE parser, ONLY : p_error, parser_end, parser_init, get_next, &
       test_next, read_line, stop_parser, search_label
  USE radial_util, ONLY : radial_integrate
  USE string_utilities, ONLY : uppercase
  USE termination, ONLY : stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: AM_ATOM_INFO, am_init, am_dealloc_ref, &
       am_get_elec_conf, am_get_z

  INTERFACE am_init
     MODULE PROCEDURE atom_init
  END INTERFACE

  INTERFACE am_dealloc_ref
     MODULE PROCEDURE am_dealloc_ref_atom
  END INTERFACE

  INTERFACE am_get_elec_conf
     MODULE PROCEDURE am_get_elec_conf_m, am_get_elec_conf_s
  END INTERFACE


!!****s* cp2k/atom_type/AM_ATOM_INFO *
!!
!!   NAME
!!     AM_ATOM_INFO
!!
!!   FUNCTION
!!     A structure to hold information (physical) about the
!!     atom to be computed
!!
!!   MODIFICATION HISTORY
!!     12-FEB-2002, TCH, working version
!!
!!   SOURCE
!!     
  TYPE AM_ATOM_INFO

     ! Cannot make this private because the attribute rg must be public.
     TYPE (radial_grid_type) :: rg

     INTEGER :: z
     INTEGER, DIMENSION(1:55) :: elec_conf
     REAL(DBL) :: ground_state_energy               
     REAL(DBL), DIMENSION(55) :: eigenvalues
     REAL(DBL), DIMENSION(:,:), POINTER :: orbitals

     CHARACTER(255) :: wf_logfile
     LOGICAL :: log_wf

  END TYPE AM_ATOM_INFO
!!****

CONTAINS

  SUBROUTINE am_init_orbitals ( obj )
    
    IMPLICIT NONE

    TYPE ( AM_ATOM_INFO ), INTENT ( OUT ) :: obj

    INTEGER :: errno, np, orbit
    REAL ( DBL ) :: res

    np = obj%rg%n

    obj%orbitals = 0.0_dbl
    DO orbit = 1, SIZE ( obj%elec_conf )
       IF ( obj%elec_conf( orbit ) == 0 ) CYCLE
       obj%orbitals( :, orbit ) = &
            EXP ( -2.25_dbl*( obj%rg%r( 1:np ) / DBLE ( orbit ) )**2 ) &
            * ( 2.25_dbl / SQRT ( 4.0_dbl*pi ) ) * obj%rg%r( 1:np )**2 &
            / DBLE ( orbit )**2
       CALL radial_integrate ( np, 4.0_dbl*pi*obj%orbitals(:,orbit)**2, obj%rg, res, errno )
       obj%orbitals( :, orbit ) = obj%orbitals( :, orbit) / SQRT ( res )
    ENDDO

  END SUBROUTINE am_init_orbitals

  SUBROUTINE atom_init(obj, globenv )

    IMPLICIT NONE

    ! arguments
    TYPE ( AM_ATOM_INFO ), INTENT ( OUT ) :: obj
    TYPE ( GLOBAL_ENVIRONMENT_TYPE ), INTENT ( IN ) :: globenv

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

    REAL ( DBL ) :: rmax                  !
    REAL ( DBL ), DIMENSION(0:5) :: b     ! 
    REAL ( DBL ) :: z                     !

    LOGICAL, DIMENSION(0:5) :: b_set   ! were supplied in the file
    LOGICAL :: nel_set, z_set 

    label = "&ATOMINFO"
    iw = globenv%scr

    b_set = .FALSE.
    nel_set = .FALSE.
    z_set = .FALSE.

    ! some default values
    nel = 0
    obj%elec_conf = 0
    obj%eigenvalues = 0.0_dbl
    b = 1.0_dbl; b(1) = 0.08_dbl
    ngp = 100
    rmax = 10._dbl

    obj%log_wf = .FALSE.

    ! open file for reading
    CALL parser_init ( globenv%input_file_name, globenv )
    CALL search_label ( label, err )
    IF ( err /= 0 ) THEN
       IF ( globenv%ionode ) WRITE ( iw,'(A)' ) "No input section &ATOMINFO found"
    ELSE
       CALL read_line ()
       DO WHILE ( test_next () /= 'X' )

          len = 8; CALL get_next ( string, len )
          CALL uppercase ( string )

          SELECT CASE ( string )

          CASE DEFAULT

          CASE ( 'Z' )
             IF ( test_next () == 'N' ) THEN 
                CALL get_next ( obj%z ); z_set = .true.
             ELSE 
                CALL stop_parser( "atom_type/atom_init", &
                     "option Z needs numeric argument" )
             END IF

          CASE ( 'NOELEC' )
             l = 0
             DO
                IF ( test_next () /= 'N' ) EXIT
                CALL get_next ( nel( l ) ); nel_set = .true.
                l = l+1
             END DO
             IF (.NOT. nel_set ) CALL stop_parser( "atom_type/atom_init", &
                  "option NOELEC needs numeric argument" )

          CASE ( 'GRIDTYPE' )
             IF ( test_next () == 'C' ) THEN
                len = 8; CALL get_next( rgtype, len )
                CALL uppercase ( rgtype )
             ELSE
                CALL stop_parser ( "atom_type/atom_init", &
                     "option GRIDTYPE needs string argument" )
             END IF

          CASE ( 'GRIDN' )
             IF ( test_next () == 'N' ) THEN 
                CALL get_next ( ngp )
             ELSE 
                CALL stop_parser ( "atom_type/atom_init", &
                     "option GRIDN needs numeric argument" )
             END IF

          CASE ( 'GRIDRMAX' )
             IF ( test_next () == 'N' ) THEN 
                CALL get_next ( rmax )
             ELSE 
                CALL stop_parser ( "atom_type/atom_init", &
                     "option GRIDRMAX needs numeric argument" )
             END IF

          CASE ('GRIDB')
             index = 0
             DO
                IF ( test_next () /= 'N' ) EXIT
                IF ( index > 5 ) &
                     CALL stop_parser ( "atom_type/atom_init", &
                     "option GRIDB can only have up to 6 arguments" )
                CALL get_next ( b( index ) )
                index = index + 1
             ENDDO

          CASE ( 'WFLOGF' )
             IF ( test_next () == 'C' ) THEN
                len = 255; CALL get_next ( obj%wf_logfile, len )
                obj%log_wf = .TRUE.
             ELSE
                CALL stop_parser ( "atom_type/atom_init", &
                     "option WFLOGF needs a filename as an argument" )
             END IF

          END SELECT
          CALL read_line ()
       END DO
    END IF

    ! close the file
    CALL parser_end ()

    ! check the supplied values
    IF (.NOT. z_set .OR. (obj%z < 1)) THEN
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
       CALL init_radial_grid ( obj%rg, "LINEAR", ngp, b(0), b(1:5), rmax )
    CASE ('LOG')
       CALL init_radial_grid ( obj%rg, "LOG", ngp, b(0), b(1:5), rmax )
    CASE ('POWER')
       CALL init_radial_grid ( obj%rg, "POWER", ngp, b(0), b(1:5), rmax )
    CASE ('RATIONAL')
       CALL init_radial_grid ( obj%rg, "RATIONAL", ngp, b(0), b(1:5), rmax )
    CASE ('ALG')
       CALL init_radial_grid ( obj%rg, "ALGEBRAIC SQUARE ROOT", ngp, b(0), b(1:5), rmax )
    END SELECT

    ! allocate memory for the electron orbitals
    IF ( ASSOCIATED ( obj%orbitals ) ) DEALLOCATE ( obj%orbitals, STAT=err )
    ALLOCATE (obj%orbitals(obj%rg%n,0:55), STAT=err)
    IF (err /=0) CALL &
         stop_memory("atom_type/atom_init","am_atom_info%orbitals",55*obj%rg%n)

    ! from the read electron numbers compute the electron configuration of the atom
    DO n=1, 10
       DO l=0, n-1
          orbit = (n*(n-1))/2 + l + 1
          obj%elec_conf( orbit ) = MIN(2*(2*l+1), nel(l) - (n-(l+1))*2*(2*l+1))
          IF ( obj%elec_conf( orbit ) < 0 ) obj%elec_conf( orbit ) = 0
       END DO
    END DO

    CALL am_init_orbitals ( obj )
    
  END SUBROUTINE atom_init

  ! ***********************************************************************

  SUBROUTINE am_dealloc_ref_atom ( obj )

    IMPLICIT NONE

    TYPE ( AM_ATOM_INFO ), INTENT ( INOUT ) :: obj

    DEALLOCATE ( obj%orbitals )

  END SUBROUTINE am_dealloc_ref_atom
  
  ! ***********************************************************************

  FUNCTION am_get_z ( obj ) RESULT ( res )
    
    IMPLICIT NONE

    TYPE ( AM_ATOM_INFO ), INTENT ( IN ) :: obj

    INTEGER :: res

    res = obj%z

  END FUNCTION am_get_z

  ! ***********************************************************************

  FUNCTION am_get_elec_conf_s ( obj, index ) RESULT ( res )
    
    IMPLICIT NONE

    TYPE ( AM_ATOM_INFO ), INTENT ( IN ) :: obj
    INTEGER, INTENT ( IN ) :: index

    INTEGER :: res

    res = obj%elec_conf( index )

  END FUNCTION am_get_elec_conf_s

  ! ***********************************************************************
  
  FUNCTION am_get_elec_conf_m ( obj ) RESULT ( res )
    
    IMPLICIT NONE

    TYPE ( AM_ATOM_INFO ), INTENT ( IN ) :: obj

    INTEGER, DIMENSION ( 1:SIZE ( obj%elec_conf ) ) :: res

    res = obj%elec_conf

  END FUNCTION am_get_elec_conf_m

  ! ***********************************************************************
  
END MODULE atom_type
