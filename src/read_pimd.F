!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/read_pimd [1.0] *
!!
!!   NAME
!!     read_pimd
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     HAF
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE read_pimd
  
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: read_pimd_section, pimd_parameters_type

  
  TYPE pimd_parameters_type
     INTEGER :: beads
     INTEGER :: ndelta
     REAL ( dbl ) :: wp
     REAL ( dbl ) :: beta
     CHARACTER ( LEN = 20 ) :: transformation
  END TYPE pimd_parameters_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** md/read_pimd_section [1.0] *
!!
!!   NAME
!!     read_pimd_section
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     HAF
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &pimd ... &end                                                   I
!!I                                                                           I
!!I beads :   n            Number of beads                                    I
!!I ndelta:   n            Number of NHC+harmonic coupling integrations per   I
!!I                        full time step dt                                  I
!!I transformation: [none,normal,stage]  coordinate transformation            I
!!I                                                                           I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!   SOURCE
!******************************************************************************

SUBROUTINE read_pimd_section ( pimdpar, mdpar )
  
  
    TYPE(pimd_parameters_type), INTENT(OUT)  :: pimdpar
    TYPE(global_environment_type), &
      INTENT(IN)                             :: mdpar

    CHARACTER(LEN=20)                        :: string
    CHARACTER(LEN=5)                         :: label
    INTEGER                                  :: ierror, ilen, iw

!------------------------------------------------------------------------------
!..defaults

  pimdpar % beads = 1
  pimdpar % ndelta = 2
  pimdpar % transformation = 'NONE'

#if defined(__parallel)
  CALL stop_program('read_pimd_section', &
     'PIMD is not checked to work in parallel yet')
#endif
  
  iw = mdpar % scr
     
!..parse the input section
  label = '&PIMD'
  CALL parser_init(mdpar % input_file_name,mdpar)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF (ierror /= 0 ) THEN
     IF( mdpar % ionode ) &
        WRITE ( iw, '( a )' ) ' No input section &PIMD found '
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 6
        CALL get_next ( string, ilen )
        CALL uppercase ( string )
        
        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_parser ( 'read_pimd_section', 'unknown option' )
           
        CASE ( 'BEADS' )
           CALL get_next ( pimdpar % beads )

        CASE ( 'NDELTA' )
           CALL get_next ( pimdpar % ndelta )
           
        CASE ( 'TRANSF' )
           ilen = 0
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           pimdpar % transformation = string ( 1:ilen )
           
        END SELECT
        CALL read_line
     END DO
        
  END IF
  CALL parser_end
!..end of parsing the input section
     
!..write some information to output
  IF (mdpar % ionode .AND. mdpar % print_level>=0) THEN
     WRITE ( iw, '( A )' ) ' PIMD| Path Integral Molecular Dynamics Protocol'
     WRITE ( iw, '( A,T71,I10 )' ) ' PIMD| Number of Beads ', &
          pimdpar % beads
     WRITE ( iw, '( A,T71,I10 )' ) &
        ' PIMD| Number of kinetic+NHC integrations ', pimdpar % ndelta
     WRITE ( iw, '( A,T61,A )' ) ' PIMD| Transformation Type ', &
          ADJUSTR ( pimdpar % transformation )
     WRITE (iw,*)
  END IF
  
END SUBROUTINE read_pimd_section

!!*****
!******************************************************************************

END MODULE read_pimd

!******************************************************************************
