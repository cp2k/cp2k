!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/brillouin [1.0] *
!!
!!   NAME
!!     brillouin
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
!
    MODULE brillouin
!
!------------------------------------------------------------------------------!
!
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

      IMPLICIT NONE
!
      PRIVATE
!
      PUBLIC :: kpoint_type, brillouin_info, kpoint_input, kpoint_dealloc

      CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="brillouin"

!
      TYPE kpoint_type
        CHARACTER (len=20) :: scheme
        INTEGER :: nk(3)
        REAL(KIND=dp) :: shift(3)
        LOGICAL :: symmetry
        INTEGER :: wfn_type
        INTEGER :: nkpt
        REAL(KIND=dp), DIMENSION (:), POINTER :: weight
        REAL(KIND=dp), DIMENSION (:,:), POINTER :: xk
      END TYPE kpoint_type
!!*****
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!!****** brillouin/kpoint_input [1.0] *
!!
!!   NAME
!!     kpoint_input
!!
!!   SYNOPSIS
!!     Subroutine kpoint_input(kp, inpar)
!!       Implicit None
!!       Type(kpoint_type), Intent (OUT):: kp
!!       Type(global_environment_type), Intent (IN):: inpar
!!     End Subroutine kpoint_input
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   INPUTS
!!I----------------------------------------------------------------------------I
!!I SECTION: &kpoint... &end                                                   I
!!I                                                                            I
!!I scheme           [Gamma, Monkhorst-Pack, MacDonald, General]               I
!!I   { nx ny nz }                                                             I
!!I   { nx ny nz  sx sy sz }                                                   I
!!I   { nkpt  x1 y1 z1 w1 ... xn yn zn wn }                                    I
!!I symmetry         [on, off]                                                 I
!!I wavefunction     [real, complex]                                           I
!!I                                                                            I
!!I----------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************
      SUBROUTINE kpoint_input(kp,inpar)


    TYPE(kpoint_type), INTENT(OUT)           :: kp
    TYPE(global_environment_type), &
      INTENT(IN)                             :: inpar

    CHARACTER(len=20)                        :: str2, string
    CHARACTER(len=6)                         :: label
    INTEGER                                  :: i, ierror, ilen, isos, iw

!------------------------------------------------------------------------------!
!..defaults

        kp%scheme = 'NULL'
        kp%symmetry = .FALSE.
        kp%wfn_type = 0
        kp%nkpt = 0
        kp%nk = 0
        kp%shift = 0.0_dp
        NULLIFY(kp%weight, kp%xk)
        iw = inpar%scr

!..parse the input section
        label = '&KPOINT'
        CALL parser_init(inpar%input_file_name,inpar)
        CALL search_label(label,ierror,ignore_case=.TRUE.)
        IF (ierror/=0) THEN
          IF (inpar%ionode .AND. inpar%print_level >= 2) &
            WRITE (iw,'(a,a)') ' No input section &KPOINT found on file ', &
              TRIM(inpar%input_file_name)
        ELSE
          CALL read_line
          DO WHILE (test_next()/='X')
            ilen = 8
            CALL get_next(string,ilen)
            CALL uppercase(string)
            SELECT CASE (string)
            CASE DEFAULT
              CALL p_error()
              CALL stop_parser('kpoint_input','unknown option')
            CASE ('SCHEME')
              ilen = 20
              CALL get_next(str2,ilen)
              CALL uppercase(str2)
              SELECT CASE (str2)
              CASE DEFAULT
                CALL p_error()
                CALL stop_program("kpoint_input","Scheme: unknown option")
              CASE ('GAMMA')
                kp%scheme = str2
              CASE ('MONKHORST-PACK')
                kp%scheme = str2
                CALL get_next ( kp%nk(1) )
                CALL get_next ( kp%nk(2) )
                CALL get_next ( kp%nk(3) )
              CASE ('MACDONALD')
                kp%scheme = str2
                CALL get_next ( kp%nk(1) )
                CALL get_next ( kp%nk(2) )
                CALL get_next ( kp%nk(3) )
                CALL get_next ( kp%shift(1) )
                CALL get_next ( kp%shift(2) )
                CALL get_next ( kp%shift(3) )
              CASE ('GENERAL')
                kp%scheme = str2
                CALL get_next ( kp%nkpt )
                ALLOCATE (kp%xk(3,kp%nkpt),STAT=isos)
                IF (isos/=0) CALL stop_memory('kpoint_input', &
                  'kp%xk',3*kp%nkpt)
                ALLOCATE (kp%weight(kp%nkpt),STAT=isos)
                IF (isos/=0) CALL stop_memory('kpoint_input', &
                  'kp%weight',kp%nkpt)
                DO i = 1, kp%nkpt
                  CALL get_next ( kp%xk(1,i) )
                  CALL get_next ( kp%xk(2,i) )
                  CALL get_next ( kp%xk(3,i) )
                  CALL get_next ( kp%weight(i) )
                END DO
              END SELECT
            CASE ('SYMMETRY')
              ilen = 3
              CALL get_next(str2,ilen)
              CALL uppercase(str2)
              IF (str2=='OFF') kp%symmetry = .FALSE.
              IF (str2=='ON') kp%symmetry = .TRUE.
            CASE ('WAVEFUNC')
              ilen = 4
              CALL get_next(str2,ilen)
              CALL uppercase(str2)
              SELECT CASE (str2)
              CASE DEFAULT
                CALL p_error()
                CALL stop_parser('kpoint_input',&
                   'Wavefunctions: unknown option')
              CASE ('REAL')
                kp%wfn_type = 0
              CASE ('COMP')
                kp%wfn_type = 1
              END SELECT
            END SELECT
          END DO
        END IF
        CALL parser_end

!..update defaults
        IF (kp%scheme=='NULL' .OR. kp%scheme=='GAMMA') THEN
          kp%nkpt = 1
          ALLOCATE (kp%xk(3,1),STAT=isos)
          IF (isos/=0) CALL stop_memory('kpoint_input','kp%xk',3)
          ALLOCATE (kp%weight(1),STAT=isos)
          IF (isos/=0) CALL stop_memory('kpoint_input','kp%weight',1)
          kp%xk(:,:) = 0.0_dp
          kp%weight(1) = 1.0_dp
          kp%scheme = 'GAMMA'
        END IF

      END SUBROUTINE kpoint_input
!!*****
!------------------------------------------------------------------------------!
!!****** brillouin/brillouin_info [1.0] *
!!
!!   NAME
!!     brillouin_info
!!
!!   SYNOPSIS
!!     Subroutine brillouin_info(kp, inpar)
!!       Implicit None
!!       Type(kpoint_type), Intent (IN):: kp
!!       Type(global_environment_type), Intent (IN):: inpar
!!     End Subroutine brillouin_info
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (27.11.2001) import environment variable
!!
!!   SOURCE
!******************************************************************************
      SUBROUTINE brillouin_info(kp,inpar)
    TYPE(kpoint_type), INTENT(IN)            :: kp
    TYPE(global_environment_type), &
      INTENT(IN)                             :: inpar

    INTEGER                                  :: i, punit

        IF ( inpar%ionode ) THEN
           punit = inpar%scr
           IF ( inpar%print_level >= 0 ) THEN
              IF (kp%scheme=='GAMMA') THEN
                WRITE (punit,*)
                WRITE (punit,'(A,T57,A)') &
                      ' BRILLOUIN|', ' Gamma-point calculation'
                WRITE (punit,'(A,T76,A)') &
                      ' BRILLOUIN| Wavefunction type', ' REAL'
              ELSE
                WRITE (punit,*)
                WRITE (punit,'(A,T61,A)') ' BRILLOUIN| K-point scheme ', &
                  ADJUSTR(kp%scheme)
                IF (kp%scheme=='MONKHORST-PACK') THEN
                  WRITE (punit,'(A,T66,3I5)') &
                      ' BRILLOUIN| K-Point grid', kp%nk
                ELSE IF (kp%scheme=='MACDONALD') THEN
                  WRITE (punit,'(A,T66,3I5)') &
                      ' BRILLOUIN| K-Point grid', kp%nk
                  WRITE (punit,'(A,T51,3F10.4)') &
                      ' BRILLOUIN| K-Point shift', kp%shift
                END IF
                IF (kp%symmetry) THEN
                  WRITE (punit,'(A,T76,A)') &
                      ' BRILLOUIN| K-Point symmetry', '   ON'
                ELSE
                  WRITE (punit,'(A,T76,A)') &
                      ' BRILLOUIN| K-Point symmetry', '  OFF'
                END IF
                IF (kp%wfn_type==0) THEN
                  WRITE (punit,'(A,T76,A)') &
                      ' BRILLOUIN| Wavefunction type', ' REAL'
                ELSE
                  WRITE (punit,'(A,T73,A)') &
                      ' BRILLOUIN| Wavefunction type', ' COMPLEX'
                END IF
                WRITE (punit,'(A,T71,I10)') &
                    ' BRILLOUIN| Number of K-points ', kp%nkpt
                WRITE (punit,'(A,T30,A,T48,A,T63,A,T78,A)') &
                    ' BRILLOUIN| Number ', 'Weight', 'X', 'Y', 'Z'
                DO i = 1, kp%nkpt
                  WRITE (punit,'(A,I5,3X,4F15.5)') &
                      ' BRILLOUIN| ', i, kp%weight(i), &
                      kp%xk(1,i), kp%xk(2,i), kp%xk(3,i)
                END DO
              END IF
           END IF
        END IF

      END SUBROUTINE brillouin_info

!!*****
!------------------------------------------------------------------------------!
!!****f* brillouin/kpoint_dealloc [1.0] *
!!
!!   NAME
!!     kpoint_dealloc
!!
!!   SYNOPSIS
!!     Subroutine kpoint_dealloc(kp, error)
!!       Type(kpoint_type), Pointer:: kp
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine kpoint_dealloc
!!
!!   FUNCTION
!!     deallocates the given kpoint info
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - kp: the object to deallocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE kpoint_dealloc(kp,error)
    TYPE(kpoint_type), POINTER               :: kp
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kpoint_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(kp)) THEN
     IF (ASSOCIATED(kp%weight)) THEN
        DEALLOCATE(kp%weight,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
     IF (ASSOCIATED(kp%xk)) THEN
        DEALLOCATE(kp%xk,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
     DEALLOCATE(kp,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
END SUBROUTINE kpoint_dealloc
!***************************************************************************

    END MODULE brillouin

!------------------------------------------------------------------------------!
