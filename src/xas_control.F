!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!***** cp2k/xas_control [1.0] *
!!
!!   NAME
!!     cp_control_types
!!
!!   FUNCTION
!!     Defines control structures, which contain the parameters and the
!!     settings for the calculations.
!!
!!   SOURCE
!******************************************************************************

MODULE xas_control

 USE cp_error_handling,                ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: xas_none, xas_tp_hh, xas_tp_fh, xas_tddft,&
                                             xas_1s_type, xas_2s_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE qs_loc_control,                  ONLY: localized_wfn_control_create,&
                                             localized_wfn_control_release,&
                                             localized_wfn_control_type,&
                                             read_loc_control_new,&
                                             read_loc_control_old
  USE qs_parser,                       ONLY: read_object,&
                                             stop_parser,&
                                             test_object
  USE string_utilities,                ONLY: uppercase,&
                                             xstring
  USE termination,                     ONLY: stop_memory,&
                                             stop_program


  IMPLICIT NONE

  PRIVATE


!!****t* xas_control/xas_control_type
!!
!! NAME
!!   xas_control_type
!!
!! FUNCTION
!!   A type that holds controling information for a xas calculation
!!
!! ATTRIBUTES
!!
!! NOTES  
!!  
!! 
!!***

  TYPE xas_control_type
    INTEGER                             :: ref_count
    CHARACTER(LEN=50)                   :: xas_rootfname
    INTEGER                             :: each_step
    INTEGER                             :: max_scf
    INTEGER                             :: nexc_atoms
    INTEGER                             :: nexc_search
    INTEGER                             :: state_type
    INTEGER                             :: xas_method
    INTEGER                             :: lb_atom, ub_atom
    INTEGER, DIMENSION(:), POINTER      :: exc_atoms 
    LOGICAL                             :: localized_wfn
    REAL(dp)                            :: eps_scf
    TYPE(localized_wfn_control_type), &
      POINTER                           :: localized_wfn_control
  END TYPE xas_control_type


  CHARACTER(LEN=*), PARAMETER :: moduleN = "xas_control"

! *** Public data types ***

  PUBLIC :: xas_control_type

! *** Public subroutines ***

  PUBLIC :: read_xas_control_new, read_xas_control_old, xas_control_create, &
            xas_control_release, xas_control_retain

!!***
! *****************************************************************************

CONTAINS

!   ***************************************************************************

  SUBROUTINE read_xas_control_new(xas_control,globenv, error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(global_environment_type)            :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER ::  routineN = "read_xas_control_new", &
      routineP = moduleN//"/"//routineN

    INTEGER                                  :: i,istat
    INTEGER,DIMENSION(:) , POINTER           :: bounds,list
    LOGICAL                                  :: failure, was_present
    TYPE(section_vals_type), POINTER         :: xas_section, loc_section


    failure=.FALSE.
    was_present = .FALSE.

    nullify(xas_section, loc_section)
    xas_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL",error=error)
    xas_section => section_vals_get_subs_vals(xas_section,"DFT",error=error)
    xas_section => section_vals_get_subs_vals(xas_section,"XAS",error=error)


    NULLIFY(bounds,list)

    CALL section_vals_val_get(xas_section,"CONVERGENCE", &
         r_val=xas_control%eps_scf,error=error)

    CALL section_vals_val_get(xas_section,"METHOD",&
         i_val=xas_control%xas_method,error=error)

    CALL section_vals_val_get(xas_section,"MAXSTEP",&
         i_val=xas_control%max_scf,error=error)

    CALL section_vals_val_get(xas_section,"EACH_STEP",&
         i_val=xas_control%each_step,error=error)

    CALL section_vals_val_get(xas_section,"STATE_TYPE",&
         i_val=xas_control%state_type,error=error)

    CALL section_vals_val_get(xas_section,"STATE_SEARCH",&
         i_val=xas_control%nexc_search,error=error)

    CALL section_vals_val_get(xas_section,"FNAME",&
         c_val=xas_control%xas_rootfname,error=error)

    CALL section_vals_val_get(xas_section,"ATOMS_LU_BOUNDS",&
         i_vals=bounds,error=error) 
    xas_control%lb_atom = bounds(1) 
    xas_control%ub_atom = bounds(2) 
    IF(xas_control%lb_atom > 0 .AND. xas_control%ub_atom > 0) THEN
      xas_control%nexc_atoms = xas_control%ub_atom - xas_control%lb_atom + 1 
      
      ALLOCATE(xas_control%exc_atoms(xas_control%nexc_atoms),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO i = 1,xas_control%nexc_atoms
        xas_control%exc_atoms(i) = xas_control%lb_atom + (i-1)
      END DO
    END IF

    IF(.NOT. ASSOCIATED(xas_control%exc_atoms)) THEN
      CALL section_vals_val_get(xas_section,"ATOMS_LIST",&
           i_vals=list,error=error)
      IF(ASSOCIATED(list)) THEN
        xas_control%nexc_atoms = SIZE(list)
        ALLOCATE(xas_control%exc_atoms(xas_control%nexc_atoms),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DO i = 1,xas_control%nexc_atoms
          xas_control%exc_atoms(i) = list(i)
        END DO
      END IF
    END IF

    IF(.NOT. ASSOCIATED(xas_control%exc_atoms)) THEN
      xas_control%nexc_atoms = 1
      xas_control%lb_atom = 1
      xas_control%ub_atom = 1
      ALLOCATE(xas_control%exc_atoms(1),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO i = 1,xas_control%nexc_atoms
        xas_control%exc_atoms(i) = xas_control%lb_atom + (i-1)
      END DO
    END IF

    CALL section_vals_val_get(xas_section,"LOCALIZE_WFN",&
         l_val=xas_control%localized_wfn,error=error)

    loc_section =>section_vals_get_subs_vals(xas_section,"LOCALIZE",error=error)
    call section_vals_get(loc_section,explicit=was_present,error=error)
    if (was_present) then
      CALL read_loc_control_new(xas_control%localized_wfn_control,loc_section,error=error)
    end if


  END SUBROUTINE read_xas_control_new


!!****f* xas_control/read_xas_control_old
!!
!!   NAME
!!      read_xas_control_old
!!   FUNCTION
!!      read the xas_control parameters from input, in the old fashion
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     05.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************

  SUBROUTINE read_xas_control_old(xas_control,string,error)

    TYPE(xas_control_type), POINTER           :: xas_control
    CHARACTER(LEN=60)                         :: string
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error

    CHARACTER(LEN=*), PARAMETER ::  routineN = "read_xas_control_old", &
      routineP = moduleN//"/"//routineN

    CHARACTER(LEN=60)                         :: str1
    INTEGER                                   :: i, istat, nat_now
    LOGICAL :: failure  

    failure = .FALSE.

    SELECT CASE (TRIM(string))
    CASE ("NONE")
      xas_control%xas_method=xas_none
    CASE ("TP_HH")
       xas_control%xas_method=xas_tp_hh
    CASE ("TP_FH")
       xas_control%xas_method=xas_tp_fh
    CASE ("TDDFT")
       xas_control%xas_method=xas_tddft
    CASE ("LOCALIZED")
       xas_control%localized_wfn=.TRUE.
       CALL read_loc_control_old(xas_control%localized_wfn_control,test_object())
    CASE ("MAXSTEP")
       CALL read_object(xas_control%max_scf)
    CASE ("CONVERGENCE")
       CALL read_object(xas_control%eps_scf)
    CASE ("EACH_STEP")
       CALL read_object(xas_control%each_step)
    CASE ("FNAME")
       CALL read_object(xas_control%xas_rootfname)
    CASE ("STATE_TYPE")
       CALL read_object(str1,lower_to_upper=.TRUE.)
       SELECT CASE (str1)
       CASE("1S")
          xas_control%state_type = xas_1s_type
       CASE("2S")
          xas_control%state_type = xas_2s_type
       CASE DEFAULT
         CALL stop_parser(routineP,"unexpected data following XAS STATE_TYPE")
       END SELECT
    CASE ("STATE_SEARCH")
       CALL read_object(xas_control%nexc_search)
    CASE ("ATOMS")
       CALL read_object(str1,lower_to_upper=.TRUE.)
       SELECT CASE (str1)
       CASE("ITOJ")
         CALL read_object(xas_control%lb_atom)
         CALL read_object(xas_control%ub_atom)
         xas_control%nexc_atoms = &
             xas_control%ub_atom - xas_control%lb_atom + 1
         IF(ASSOCIATED(xas_control%exc_atoms))  THEN
            DEALLOCATE(xas_control%exc_atoms, STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         END IF
         ALLOCATE(xas_control%exc_atoms(xas_control%nexc_atoms),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DO i = 1,xas_control%nexc_atoms
           xas_control%exc_atoms(i) = xas_control%lb_atom + (i-1)
         END DO
       CASE("LIST")
         IF(.NOT.ASSOCIATED(xas_control%exc_atoms)) THEN
            ALLOCATE(xas_control%exc_atoms(50),STAT=istat)
            CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            xas_control%exc_atoms = 0
            nat_now = 0
         ELSE
            nat_now = SIZE(xas_control%exc_atoms,1)
            CALL reallocate(xas_control%exc_atoms,1,50+nat_now)
         ENDIF
         DO
           IF(test_object()=="EOL") THEN
               EXIT
           ELSE
             nat_now = nat_now + 1
             IF(nat_now > SIZE(xas_control%exc_atoms,1))THEN
               CALL reallocate(xas_control%exc_atoms,1,nat_now)
             END IF
             CALL read_object(xas_control%exc_atoms(nat_now))
           END IF
         END DO
         CALL reallocate(xas_control%exc_atoms,1,nat_now)
         xas_control%nexc_atoms = nat_now
       CASE DEFAULT
         CALL stop_parser(routineP,"unexpected data following XAS ATOMS")
       END SELECT
       
    CASE DEFAULT
      CALL stop_program(routineP,"UNKNOWN XAS KEYWORD "//&
           TRIM(string))
    END SELECT

  END SUBROUTINE read_xas_control_old


!****f* xas_control/xas_control_create/retain/release *
!!
!!   NAME
!!      xas_control_create/retain/release
!!   FUNCTION
!!      create retain release the xas_control_type
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE xas_control_create(xas_control,error)

    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
    TYPE(xas_control_type), POINTER           :: xas_control
 
    CHARACTER(len=*), PARAMETER :: routineN = 'xas_control_create', &
       routineP = moduleN//':'//routineN

    INTEGER ::  stat
    LOGICAL :: failure
    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(xas_control),cp_failure_level,routineP,error,failure)
    ALLOCATE(xas_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    xas_control%ref_count= 1
    xas_control%each_step= 1
    xas_control%xas_method = xas_tp_hh
    xas_control%xas_rootfname = ""
    xas_control%nexc_atoms = 1
    xas_control%nexc_atoms = -1
    xas_control%lb_atom = 1
    xas_control%ub_atom = 1
    xas_control%state_type = xas_1s_type    
    xas_control%max_scf = 150
    xas_control%eps_scf = 0.5d-6
    xas_control%localized_wfn = .FALSE. 
    NULLIFY(xas_control%localized_wfn_control)
    CALL localized_wfn_control_create(xas_control%localized_wfn_control,error=error)
    NULLIFY(xas_control%exc_atoms)

  END SUBROUTINE xas_control_create

  SUBROUTINE xas_control_release(xas_control,error)

    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
    TYPE(xas_control_type), POINTER           :: xas_control

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_control_release', &
       routineP = moduleN//':'//routineN
    INTEGER :: istat
    LOGICAL :: failure
    failure=.FALSE.
    IF(ASSOCIATED(xas_control)) THEN
      CPPrecondition(xas_control%ref_count>0,cp_failure_level,routineP,error,failure)
      xas_control%ref_count=xas_control%ref_count-1
      IF (xas_control%ref_count==0) THEN
        IF (ASSOCIATED(xas_control%localized_wfn_control)) THEN
          CALL localized_wfn_control_release(xas_control%localized_wfn_control,error=error)
        END IF
        IF (ASSOCIATED(xas_control%exc_atoms)) THEN
          DEALLOCATE(xas_control%exc_atoms,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        DEALLOCATE(xas_control,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
    END IF

  END SUBROUTINE xas_control_release

  SUBROUTINE xas_control_retain(xas_control,error)

    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
    TYPE(xas_control_type), POINTER           :: xas_control

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_control_retain', &
       routineP = moduleN//':'//routineN
    INTEGER :: istat
    LOGICAL :: failure
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(xas_control),cp_failure_level,routineP,error,failure)
    xas_control%ref_count=xas_control%ref_count+1

  END SUBROUTINE xas_control_retain

END MODULE xas_control
