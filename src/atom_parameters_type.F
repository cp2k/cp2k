!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****m* cp2k/atom_parameters_type
!!
!!   NAME
!!     atom_parameters_type
!!
!!   FUNCTION
!!     Contains the AM_PARAMETERS_INFO type which is used to 
!!     keep track of all possible program parameters
!!
!!   SYNOPSIS
!!     am_init ( obj, globenv )
!!     am_dealloc_ref ( obj )
!!     am_get_dma ( obj )
!!     am_set_dma ( obj, dma )
!!     am_get_ldd_rescale_value ( obj )
!!     am_get_mothod_parameters ( obj, index )
!!     am_get_method ( obj )
!!     am_get_maxiter ( obj )
!!     am_get_maxshoot ( obj )
!!     am_get_ncv ( obj )
!!     am_get_presision ( obj )
!!
!!   MODIFICATION HISTORY
!!     TCH, 02-FEB-2002, working version
!!
!!   SEE ALSO
!!     
!!****
MODULE atom_parameters_type

  USE kinds, ONLY : dbl

  IMPLICIT NONE

  PRIVATE

  PUBLIC AM_PARAMETERS_INFO, am_init, am_dealloc_ref

  !! getters and setters
  PUBLIC &
       am_get_dma, am_set_dma, &
       am_get_method_parameters, &
       am_get_method, am_get_maxshoot, &
       am_get_ncv, am_get_precision, &
       am_get_maxiter, am_print_parameters, &
       am_get_ldd_rescale_value
  
!!****f* cp2k/atom_parameters_type/am_init
!!
!!   NAME
!!     am_init
!!
!!   FUNCTION
!!     Initializes the parameters object. This function is overloaded with
!!     an interface, so every object in the atom_module framework is treated
!!     the same way.
!!
!!   ARGUMENTS
!!     obj - TYPE ( AM_PARAMETERS_INFO ) - the object to be initialized itself
!!     globenv - TYPE ( GLOBAL_ENVIRONMENT_TYPE ) - global program parameters
!!
!!   MODIFICATION HISTORY
!!     12-FEB-2002, TCH, working version
!!
!!   SEE ALSO
!!     AM_PARAMETERS_INFO
!!
!!****
  INTERFACE am_init
     MODULE PROCEDURE am_init_parameters
  END INTERFACE

!!****f* cp2k/atom_parameters_type/am_dealloc_ref *
!!
!!   NAME
!!     am_dealloc_ref
!!
!!   FUNCTION
!!     Destroys the object and frees up memory.
!!  
!!   ARGUMENTS
!!     obj - TYPE ( AM_PARAMETERS_INFO ) - the object to be initialized itself
!!
!!   MODIFICATION HISTORY
!!     12-FEB-2002, TCH, working version
!!
!!   SEE ALSO
!!     AM_PARAMETERS_INFO
!!
!!****
  INTERFACE am_dealloc_ref
     MODULE PROCEDURE am_dealloc_ref_parameters
  END INTERFACE
  
!!****s* cp2k/atom_parameters_type/AM_PARAMETERS_INFO
!!
!!   NAME
!!     AM_PARAMETERS_INFO
!!
!!   FUNCTION
!!     This is the structure that holds the program parameters
!!
!!   MODIFICATION HISTORY
!!     02-FEB-2002, TCH, working version
!!
!!   SOURCE
!!
  TYPE AM_PARAMETERS_INFO

     PRIVATE

     CHARACTER ( 8 ) :: method
     REAL ( DBL ), DIMENSION ( 5 ) :: method_parameters 

     REAL ( DBL ) :: density_mixing_angle
     REAL ( DBL ) :: ldd_rescale_value
     INTEGER :: maxiter                       
     INTEGER :: maxshoot                      
     REAL ( DBL ) :: node_correction_value
     REAL ( DBL ) :: precision                

     CHARACTER(255) :: wf_logfile
     LOGICAL :: log_wf

  END TYPE AM_PARAMETERS_INFO
!!****

CONTAINS
  
  !! ********************************************************************


  SUBROUTINE am_init_parameters ( obj, globenv )

    USE global_types, ONLY : GLOBAL_ENVIRONMENT_TYPE
    USE mathconstants, ONLY : pi
    USE parser, ONLY : parser_init, get_next, read_line, &
         search_label, test_next, parser_end, stop_parser
    USE string_utilities, ONLY : uppercase
    USE termination, ONLY : stop_program
    USE functionals, ONLY : init_functionals

    IMPLICIT NONE
    
    TYPE ( AM_PARAMETERS_INFO ), INTENT ( INOUT ) :: obj
    TYPE ( GLOBAL_ENVIRONMENT_TYPE ), INTENT ( IN ) :: globenv

    CHARACTER ( len = 8 ) :: label
    CHARACTER ( len = 20 ) :: string
    INTEGER :: errno, index
    INTEGER :: len
    LOGICAL, DIMENSION ( 5 ) :: method_parameters_set
    LOGICAL :: method_set
    
    ! initial values if not set in config file
    obj%density_mixing_angle = 0.5_dbl
    obj%ldd_rescale_value = 1.0_dbl
    obj%maxiter = 50
    obj%maxshoot = 1500
    obj%method_parameters( 1 ) = 0.6666667_dbl
    obj%node_correction_value = 0.5_dbl
    obj%precision = 1.0E-6_dbl
    
    ! other initialisations
    method_parameters_set = .FALSE.
    method_set = .FALSE.
    
    CALL parser_init( globenv%input_file_name, globenv )
    label = "&ATOMINFO"
    CALL search_label( label, errno )
    IF ( errno /= 0 ) THEN
       IF ( globenv%ionode ) WRITE ( globenv%scr, '(A)' ) "No input section &ATOMINFO found"
    ELSE
       CALL read_line()
       DO WHILE ( test_next() /= 'X' )

          len = 8; CALL get_next ( string, len )
          CALL uppercase( string )

          SELECT CASE ( string )
             
          CASE DEFAULT

          CASE ('DMA')
             IF (test_next() == 'N') THEN 
                CALL get_next ( obj%density_mixing_angle )
             ELSE 
                CALL stop_parser("atom_parameters_type/am_init_parameters", &
                     "option DMA needs numeric argument")
             END IF

          CASE ('LDDRSV')
             IF (test_next() == 'N') THEN 
                CALL get_next ( obj%ldd_rescale_value )
             ELSE 
                CALL stop_parser("atom_parameters_type/am_init_parameters", &
                     "option LDDRSV needs numeric argument")
             END IF

          CASE ('MAXITER')
             IF (test_next() == 'N') THEN 
                CALL get_next ( obj%maxiter )
             ELSE 
                CALL stop_parser("atom_parameters_type/am_init_parameters", &
                     "option MAXITER needs numeric argument")
             END IF

          CASE ('MAXSHOOT')
             IF (test_next() == 'N') THEN 
                CALL get_next ( obj%maxshoot )
             ELSE 
                CALL stop_parser("atom_parameters_type/am_init_parameters", &
                     "option MAXSHOOT needs numeric argument")
             END IF

          CASE ('METHOD')
             IF (test_next() == 'C') THEN
                len = 8; CALL get_next( obj%method, len )
                CALL uppercase( obj%method )
                method_set = .TRUE.
             ELSE
                CALL stop_parser( "atom_parameters_type/am_init_parameters", &
                     "option METHOD needs string argument" )
             END IF

          CASE ('MPARA')
             IF ( test_next() /= 'N' ) THEN
                CALL stop_parser( "atom_parameters_type/am_init_parameters", &
                     "option MPARA needs at least a numeric argument" )
             ENDIF
             index = 1
             DO
                IF ( test_next() /= 'N' ) EXIT
                IF ( index > 5 ) THEN
                   CALL stop_parser( "atom_parameters_type/am_init_parameters", &
                        "option MPARA can only have up to 5 arguments" )
                ENDIF
                CALL get_next ( obj%method_parameters( index ) )
                method_parameters_set( index ) = .TRUE.
                index = index + 1
             END DO

          CASE ('NCV')
             IF (test_next() == 'N') THEN 
                CALL get_next ( obj%node_correction_value )
             ELSE 
                CALL stop_parser("atom_parameters_type/am_init_parameters", &
                     "option NCV needs numeric argument")
             END IF

          CASE ('PRECI')
             IF (test_next() == 'N') THEN 
                CALL get_next ( obj%precision )
             ELSE 
                CALL stop_parser("atom_parameters_type/am_init_parameters", &
                     "option PRECISION needs numeric argument")
             END IF


          END SELECT
          
          CALL read_line()
       ENDDO
    ENDIF

    CALL parser_end()

    IF ( .NOT. method_set ) THEN
       CALL stop_program ( "atom_parameters_type/am_init_parameters", &
            "the option METHOD must at least be specified" )
    ENDIF

    CALL init_functionals( xalpha = obj%method_parameters( 1 ) )
    IF ( obj%method == "XA" ) THEN
       obj%method_parameters( 2 ) = &
            1.5_dbl * obj%method_parameters( 1 ) * (3.0_dbl/pi)**(1.0_dbl/3.0_dbl)
    ENDIF
    
  END SUBROUTINE am_init_parameters

  !! ********************************************************************

  SUBROUTINE am_dealloc_ref_parameters ( obj )

    IMPLICIT NONE

    TYPE ( AM_PARAMETERS_INFO ), INTENT ( INOUT ) :: obj

  END SUBROUTINE am_dealloc_ref_parameters

  !! ********************************************************************

  FUNCTION am_get_dma ( obj )

    IMPLICIT NONE

    TYPE ( AM_PARAMETERS_INFO ), INTENT ( IN ) :: obj
    
    REAL ( DBL ) :: am_get_dma

    am_get_dma = obj%density_mixing_angle

  END FUNCTION am_get_dma
  
  ! ********************************************************************

  SUBROUTINE am_set_dma ( obj, dma )

    IMPLICIT NONE

    TYPE ( AM_PARAMETERS_INFO ), INTENT ( INOUT ) :: obj

    REAL ( DBL ) :: dma

    obj%density_mixing_angle = dma

  END SUBROUTINE am_set_dma

  ! ********************************************************************

  FUNCTION am_get_method_parameters ( obj, index ) RESULT ( res ) 
    
    USE termination, ONLY : stop_program
    
    IMPLICIT NONE
    
    TYPE ( AM_PARAMETERS_INFO ), INTENT ( IN ) :: obj
    INTEGER, INTENT ( IN ) :: index

    REAL ( DBL ) :: res

    IF ( index < 1 .OR. index > 5 ) THEN
       CALL stop_program("this be better done with an assert"," a d")
    ENDIF
    res = obj%method_parameters( index )

  END FUNCTION am_get_method_parameters

  ! ********************************************************************

  FUNCTION am_get_method ( obj ) RESULT ( res )
    
    IMPLICIT NONE

    TYPE ( AM_PARAMETERS_INFO ), INTENT ( IN ) :: obj

    CHARACTER ( len = 8 ) :: res

    res = obj%method

  END FUNCTION am_get_method

  ! ********************************************************************

  FUNCTION am_get_maxshoot ( obj ) RESULT ( res )

    IMPLICIT NONE

    TYPE ( AM_PARAMETERS_INFO ), INTENT ( IN ) :: obj

    INTEGER :: res

    res = obj%maxshoot

  END FUNCTION am_get_maxshoot

  ! ********************************************************************

  FUNCTION am_get_ldd_rescale_value ( obj ) RESULT ( res )

    IMPLICIT NONE

    TYPE ( AM_PARAMETERS_INFO ), INTENT ( IN ) :: obj

    REAL ( DBL ) :: res

    res = obj%ldd_rescale_value

  END FUNCTION am_get_ldd_rescale_value

  ! ********************************************************************

  FUNCTION am_get_ncv ( obj ) RESULT ( res )

    IMPLICIT NONE

    TYPE ( AM_PARAMETERS_INFO ), INTENT ( IN ) :: obj

    REAL ( DBL ) :: res

    res = obj%node_correction_value

  END FUNCTION am_get_ncv

  ! ********************************************************************

  FUNCTION am_get_precision ( obj ) RESULT ( res )

    IMPLICIT NONE

    TYPE ( AM_PARAMETERS_INFO ), INTENT ( IN ) :: obj

    REAL ( DBL ) :: res

    res = obj%precision

  END FUNCTION am_get_precision

  ! ********************************************************************

  FUNCTION am_get_maxiter ( obj ) RESULT ( res )

    IMPLICIT NONE

    TYPE ( AM_PARAMETERS_INFO ), INTENT ( IN ) :: obj

    INTEGER :: res

    res = obj%maxiter

  END FUNCTION am_get_maxiter

  ! ********************************************************************

  SUBROUTINE am_print_parameters ( obj )
    
    IMPLICIT NONE

    TYPE ( AM_PARAMETERS_INFO ), INTENT ( IN ) :: obj

    !!! nothing for now

  END SUBROUTINE am_print_parameters

  ! ********************************************************************

END MODULE atom_parameters_type
