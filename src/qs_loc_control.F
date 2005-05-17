!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!***** cp2k/qs_loc_control [1.0] *
!!
!!   NAME
!!     qs_loc_control
!!
!!   FUNCTION
!!     Defines the control structure, which contain the parameters and the
!!     settings for the calculations of the localization of wavefunctions.
!!
!!   SOURCE
!******************************************************************************

MODULE qs_loc_control

  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_keywords,                     ONLY: find_keyword,&
                                             keyword_type,&
                                             list_keywords
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_loc_min,&
                                             do_loc_max,&
                                             do_loc_jacobi,&
                                             do_loc_crazy,&
                                             do_loc_direct,&
                                             do_loc_none,&
                                             op_loc_berry,&
                                             op_loc_boys,&
                                             op_loc_pipek,&
                                             state_loc_all,&
                                             state_loc_upton,&
                                             state_loc_fromitoj,&
                                             state_loc_list,&
                                             state_loc_none
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE memory_utilities,                ONLY: reallocate
  USE qs_parser,                       ONLY: read_object,&
                                             stop_parser,&
                                             test_object
  USE string_utilities,                ONLY: uppercase,&
                                             xstring
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE


!!***** qs_loc_control/localized_wfn_control_type
!!
!! NAME
!!   localized_wfn_control_type
!!
!! FUNCTION
!!   A type that holds controling information for the
!!   calculation of the spread of wfn and theoptimization of
!!   the spread functional 
!!
!! ATTRIBUTES
!!   ref_count:
!!   each_step:  how often the calculation is done
!!   localization_method: which algorithm is used for the optimization
!!   operator_type: how the spread is defined
!!   nloc_states: number of states on which the spread is computed
!!   set_of_states: how to choose the states
!!   lu_bound_states: lower and upper bounds of the set of states
!!   print_cubes: 
!!   print_centers:
!!   print_spreads:
!!   loc_states: list of states on which the spread is computed
!!   centers_set:  arrais containing centers and spreads of the selected wfn
!!   centers_file_name: output file names
!!   spreads_file_name: output file names 
!!
!! NOTES
!!
!!
!!***
  TYPE localized_wfn_control_type
    INTEGER                         :: ref_count
    INTEGER                         :: iter_count
    INTEGER                         :: each_step
    INTEGER                         :: min_or_max
    INTEGER                         :: localization_method
    INTEGER                         :: operator_type
    INTEGER, DIMENSION(2)           :: nloc_states
    INTEGER                         :: set_of_states
    INTEGER, DIMENSION(2,2)         :: lu_bound_states
    LOGICAL                         :: print_cubes
    LOGICAL                         :: print_centers
    LOGICAL                         :: print_spreads
    INTEGER, POINTER, DIMENSION(:)  :: loc_states
    TYPE(cp_2d_r_p_type), DIMENSION(2)  :: centers_set
    CHARACTER ( LEN = 50 )          :: loc_rootfname
    CHARACTER ( LEN = 50 ),DIMENSION(2):: centers_file_name
    CHARACTER ( LEN = 50 ),DIMENSION(2):: spreads_file_name
  END TYPE localized_wfn_control_type

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_loc_control"

! *** Public data types ***

  PUBLIC :: localized_wfn_control_type

! *** Public subroutines ***

  PUBLIC :: localized_wfn_control_create, &
            localized_wfn_control_release,&
            localized_wfn_control_retain,&
            read_loc_control_old,&
            read_loc_control_new,&
            set_loc_wfn_lists,&
            set_loc_centers

!!***
! *****************************************************************************

CONTAINS

!!****f* qs_loc_control/localized_wfn_control_create/retain/release *
!!
!!   NAME
!!      localized_wfn_control_create/retain/release
!!   FUNCTION
!!      create retain release the localized_wfn_control_type
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE localized_wfn_control_create(localized_wfn_control,error)
   TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
   TYPE(localized_wfn_control_type), POINTER :: localized_wfn_control
   CHARACTER(len=*), PARAMETER :: routineN = 'localized_wfn_control_create', &
      routineP = module_name//':'//routineN

   INTEGER :: stat
   LOGICAL :: failure
   failure=.FALSE.

   CPPrecondition(.NOT.ASSOCIATED(localized_wfn_control),cp_failure_level,routineP,error,failure)
   ALLOCATE(localized_wfn_control,stat=stat)
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

   localized_wfn_control%ref_count= 1
   localized_wfn_control%iter_count= 0
   localized_wfn_control%each_step= 1
   localized_wfn_control%min_or_max= do_loc_min
   localized_wfn_control%localization_method =do_loc_jacobi
   localized_wfn_control%operator_type = op_loc_berry
   localized_wfn_control%set_of_states = state_loc_all
   localized_wfn_control%nloc_states=0
   localized_wfn_control%lu_bound_states=0
   localized_wfn_control%loc_rootfname=""
   localized_wfn_control%print_cubes = .FALSE.
   localized_wfn_control%print_centers = .FALSE.
   localized_wfn_control%print_spreads = .FALSE.
   NULLIFY(localized_wfn_control%loc_states)
   NULLIFY(localized_wfn_control%centers_set(1)%array)
   NULLIFY(localized_wfn_control%centers_set(2)%array)

  END SUBROUTINE localized_wfn_control_create

  SUBROUTINE localized_wfn_control_release(localized_wfn_control,error)

    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(localized_wfn_control_type), POINTER :: localized_wfn_control
    CHARACTER(len=*), PARAMETER :: routineN = 'localized_wfn_control_release', &
       routineP = module_name//':'//routineN
    INTEGER :: istat
    LOGICAL :: failure
    failure=.FALSE.
    IF(ASSOCIATED(localized_wfn_control)) THEN
      CPPrecondition(localized_wfn_control%ref_count>0,cp_failure_level,routineP,error,failure)
      localized_wfn_control%ref_count=localized_wfn_control%ref_count-1
      IF (localized_wfn_control%ref_count==0) THEN
         IF (ASSOCIATED(localized_wfn_control%loc_states)) &
             DEALLOCATE(localized_wfn_control%loc_states,STAT=istat)
         IF (ASSOCIATED(localized_wfn_control%centers_set(1)%array)) &
             DEALLOCATE(localized_wfn_control%centers_set(1)%array,STAT=istat)
         IF (ASSOCIATED(localized_wfn_control%centers_set(2)%array)) &
             DEALLOCATE(localized_wfn_control%centers_set(2)%array,STAT=istat)
         localized_wfn_control%ref_count=0
         DEALLOCATE(localized_wfn_control,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ENDIF
    END IF
  END SUBROUTINE localized_wfn_control_release

  SUBROUTINE localized_wfn_control_retain(localized_wfn_control,error)
   TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
   TYPE(localized_wfn_control_type), POINTER :: localized_wfn_control
   CHARACTER(len=*), PARAMETER :: routineN = 'localized_wfn_control_retain', &
      routineP = module_name//':'//routineN
   LOGICAL :: failure
   failure=.FALSE.
   CPPrecondition(ASSOCIATED(localized_wfn_control),cp_failure_level,routineP,error,failure)

   localized_wfn_control%ref_count=localized_wfn_control%ref_count+1
  END SUBROUTINE localized_wfn_control_retain

!!****f* qs_loc_control/read_loc_control_new *
!!
!!   NAME
!!      read_loc_control_new
!!   FUNCTION
!!      read the controlparameter from input, using the new input scheme
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     05.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************

  SUBROUTINE read_loc_control_new(localized_wfn_control,loc_section,error)

    TYPE(localized_wfn_control_type), POINTER :: localized_wfn_control
    TYPE(section_vals_type), POINTER          :: loc_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error

    CHARACTER(LEN=*), PARAMETER ::  routine_name = "read_loc_control_new", &
      routineP = module_name//"/"//routine_name

    INTEGER                                   :: nloc_now, i, istat
    INTEGER,DIMENSION(:) , POINTER            :: bounds,list
    LOGICAL                                   :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(loc_section),cp_failure_level,routineP,error,failure)

    IF(.NOT. failure)  THEN

      nullify (bounds, list ) 

      CALL section_vals_val_get(loc_section,"EACH_STEP",&
           i_val=localized_wfn_control%each_step,error=error)

      CALL section_vals_val_get(loc_section,"MIN_OR_MAX",&
           i_val=localized_wfn_control%min_or_max,error=error)

      CALL section_vals_val_get(loc_section,"METHOD",&
           i_val=localized_wfn_control%localization_method,error=error)

      CALL section_vals_val_get(loc_section,"OPERATOR",&
           i_val=localized_wfn_control%operator_type,error=error)

      CALL section_vals_val_get(loc_section,"SET_OF_STATES",&
           i_val=localized_wfn_control%set_of_states, error=error)

      IF(localized_wfn_control%set_of_states==state_loc_all) THEN

      ELSE IF(localized_wfn_control%set_of_states==state_loc_upton) THEN
         CALL section_vals_val_get(loc_section,"SET_OF_STATES",&
              i_val=localized_wfn_control%lu_bound_states(2,1), error=error)
         localized_wfn_control%lu_bound_states(1,1) = 1
         localized_wfn_control%nloc_states(1) = &
                     localized_wfn_control%lu_bound_states(2,1) 
      ELSE IF(localized_wfn_control%set_of_states==state_loc_fromitoj) THEN
         CALL section_vals_val_get(loc_section,"LU_BOUNDS",&
              i_vals=bounds,error=error)
              localized_wfn_control%lu_bound_states(1,1) = bounds(1)
              localized_wfn_control%lu_bound_states= bounds(2)

         localized_wfn_control%nloc_states(1) = &
                     localized_wfn_control%lu_bound_states(2,1)-&
                     localized_wfn_control%lu_bound_states(1,1)+1 
      ELSE IF(localized_wfn_control%set_of_states==state_loc_list) THEN
         CALL section_vals_val_get(loc_section,"LISTSTATES",&
              i_vals=list,error=error)
         CPPrecondition(ASSOCIATED(list),cp_failure_level,routineP,error,failure)
         localized_wfn_control%nloc_states(1) = SIZE(list,1)
         ALLOCATE(localized_wfn_control%loc_states(SIZE(list,1)),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DO i = 1,localized_wfn_control%nloc_states(1)
           localized_wfn_control%loc_states(i) = list(i) 
         END DO
         
      ELSE IF(localized_wfn_control%set_of_states==state_loc_none) THEN
         localized_wfn_control%nloc_states(:) = 0
      END IF

      CALL section_vals_val_get(loc_section,"FNAME",&
           c_val=localized_wfn_control%loc_rootfname,error=error)


    END IF



  END SUBROUTINE read_loc_control_new

!!****f* qs_loc_control/read_loc_control_old *
!!
!!   NAME
!!      read_loc_control_old
!!   FUNCTION
!!      read the controlparameter from input, in the old fashion
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************

  SUBROUTINE read_loc_control_old(localized_wfn_control,test_result, error)

    TYPE(localized_wfn_control_type), POINTER :: localized_wfn_control
    CHARACTER(LEN=3)                          :: test_result
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error

    CHARACTER(LEN=*), PARAMETER ::  routine_name = "read_loc_control_old", &
      routineP = module_name//"/"//routine_name


    CHARACTER(LEN=60)                         :: keyword, string, string2
    INTEGER                                   :: nloc_now, i, istat
    LOGICAL                                   :: failure

    failure = .FALSE.

    SELECT CASE (test_result)
    CASE ("STR")
      CALL read_object(string,lower_to_upper=.TRUE.)
      SELECT CASE (string)
      CASE ("METHOD")
        CALL read_object(string2,lower_to_upper=.TRUE.)
        SELECT CASE (string2)
        CASE ("JACOBI")
           localized_wfn_control%localization_method = do_loc_jacobi
        CASE ("CRAZY")
           localized_wfn_control%localization_method = do_loc_crazy
        CASE ("DIRECT")
           localized_wfn_control%localization_method = do_loc_direct
        CASE ("NONE")
           localized_wfn_control%localization_method = do_loc_none
        CASE DEFAULT
           CALL stop_parser(routineP,"unexpected data following LOCALIZED_WFN METHOD")
        END SELECT
      CASE ("MIN_OR_MAX")
        CALL read_object(string2,lower_to_upper=.TRUE.)
        SELECT CASE (string2)
        CASE ("SPREADMIN")
           localized_wfn_control%min_or_max = do_loc_min
        CASE ("SPREADMAX")
           localized_wfn_control%min_or_max = do_loc_max
        CASE DEFAULT
           CALL stop_parser(routineP,"unexpected data following LOCALIZED_WFN MIN_OR_MAX")
        END SELECT
      CASE ("OPERATOR")
        CALL read_object(string2,lower_to_upper=.TRUE.)
        SELECT CASE (string2)
        CASE("BERRY")
           localized_wfn_control%operator_type = op_loc_berry
        CASE("BOYS")
           localized_wfn_control%operator_type = op_loc_boys
        CASE("PIPEK")
           localized_wfn_control%operator_type = op_loc_pipek
        CASE DEFAULT
           CALL stop_parser(routineP,"unexpected data following LOCALIZED_WFN OPERATOR ")
        END SELECT 
      CASE ("SET_OF_STATES")
        CALL read_object(string2,lower_to_upper=.TRUE.)
        SELECT CASE (string2)
        CASE("ALL")
           localized_wfn_control%set_of_states = state_loc_all
        CASE("NONE")
           localized_wfn_control%set_of_states = state_loc_none
        CASE("UPTON")
           localized_wfn_control%set_of_states = state_loc_upton
           CALL read_object(localized_wfn_control%lu_bound_states(2,1))    
           localized_wfn_control%lu_bound_states(1,1) = 1
           localized_wfn_control%nloc_states(1) = &
                     localized_wfn_control%lu_bound_states(2,1) 
        CASE("FROMITOJ")
           localized_wfn_control%set_of_states = state_loc_fromitoj
           CALL read_object(localized_wfn_control%lu_bound_states(1,1))
           CALL read_object(localized_wfn_control%lu_bound_states(2,1))
           localized_wfn_control%nloc_states(1) = &
                     localized_wfn_control%lu_bound_states(2,1)-&
                     localized_wfn_control%lu_bound_states(1,1)+1 
        CASE("LIST")
           localized_wfn_control%set_of_states = state_loc_list 
           IF(.NOT.ASSOCIATED(localized_wfn_control%loc_states)) THEN
             ALLOCATE(localized_wfn_control%loc_states(50),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             localized_wfn_control%loc_states = 0
             nloc_now = 0
           ELSE
             nloc_now = SIZE(localized_wfn_control%loc_states,1)
             CALL reallocate(localized_wfn_control%loc_states,1,50+nloc_now)
           ENDIF  
           DO 
             IF(test_object()=="EOL") THEN
                 EXIT
             ELSE
               nloc_now = nloc_now + 1
               IF(nloc_now > SIZE(localized_wfn_control%loc_states,1))THEN
                 CALL reallocate(localized_wfn_control%loc_states,1,nloc_now)
               END IF  
               CALL read_object(localized_wfn_control%loc_states(nloc_now))
             END IF 
           END DO 
           CALL reallocate(localized_wfn_control%loc_states,1,nloc_now)
           localized_wfn_control%nloc_states(1) = nloc_now
        CASE DEFAULT
           CALL stop_parser(routineP,"unexpected data following LOCALIZED_WFN ")
        END SELECT
      CASE ("EACH_STEP")
        CALL read_object(localized_wfn_control%each_step)
      CASE ("FNAME")
        CALL read_object(localized_wfn_control%loc_rootfname)
      END SELECT
    CASE ("EOL")
      ! do nothing
    CASE DEFAULT
      CALL stop_parser(routineP,"unexpected data following LOCALIZED_WFN") 
    END SELECT

  END SUBROUTINE read_loc_control_old

!!****f* qs_loc_control/set_loc_centers *
!!
!!   NAME
!!      set_loc_centers
!!   FUNCTION
!!      create the center and spread array and the file names for the output
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************

  SUBROUTINE set_loc_centers(localized_wfn_control,nmoloc,nspins,globenv,error)

    TYPE(localized_wfn_control_type)          :: localized_wfn_control
    INTEGER, DIMENSION(2), INTENT(IN)         :: nmoloc
    INTEGER, INTENT(IN)                       :: nspins
    TYPE(global_environment_type), POINTER    :: globenv
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_loc_centers', &
       routineP = module_name//':'//routineN

    INTEGER :: ia, ie, ispin, istat
    LOGICAL :: failure

    failure = .FALSE.

    IF(localized_wfn_control%loc_rootfname==" ") THEN
      CALL xstring ( globenv % project_name, ia, ie )
      DO ispin = 1,nspins
         localized_wfn_control%spreads_file_name(ispin) = globenv%project_name(ia:ie)//&
                 '_spreads_s'//TRIM(ADJUSTL(cp_to_string(ispin)))
         localized_wfn_control%centers_file_name(ispin) = globenv%project_name(ia:ie)//&
                 '_centers_s'//TRIM(ADJUSTL(cp_to_string(ispin)))
      END DO
    ELSE 
      CALL xstring (localized_wfn_control%loc_rootfname, ia, ie )
         localized_wfn_control%spreads_file_name(ispin) = &
                 localized_wfn_control%loc_rootfname(ia:ie)//&
                 '_spreads_s'//TRIM(ADJUSTL(cp_to_string(ispin)))
         localized_wfn_control%centers_file_name(ispin) = &
                 localized_wfn_control%loc_rootfname(ia:ie)//&
                 '_centers_s'//TRIM(ADJUSTL(cp_to_string(ispin)))
    END IF
    DO ispin = 1,nspins
      ALLOCATE(localized_wfn_control%centers_set(ispin)%array(6,nmoloc(ispin)), STAT = istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      localized_wfn_control%spreads_file_name(ispin) = globenv%project_name(ia:ie)//&
                 '_spreads_s'//TRIM(ADJUSTL(cp_to_string(ispin)))
      localized_wfn_control%centers_file_name(ispin) = globenv%project_name(ia:ie)//&
                 '_centers_s'//TRIM(ADJUSTL(cp_to_string(ispin)))
    END DO

    localized_wfn_control%print_centers = globenv%print%localization
    localized_wfn_control%print_spreads = globenv%print%localization
    localized_wfn_control%print_cubes = globenv%print%cubes

  END SUBROUTINE set_loc_centers
 
! *****************************************************************************


!!****f* qs_loc_control/set_loc_wfn_lists *
!!
!!   NAME
!!      loc_wfn_lists
!!   FUNCTION
!!      create the lists of mos that are taken into account 
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************

  SUBROUTINE set_loc_wfn_lists(localized_wfn_control,nmoloc,nmo,nspins,error)

    TYPE(localized_wfn_control_type)          :: localized_wfn_control
    INTEGER, DIMENSION(2), INTENT(IN)         :: nmoloc,nmo
    INTEGER, INTENT(IN)                       :: nspins
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_loc_wfn_lists', &
       routineP = module_name//':'//routineN

    INTEGER :: i, iloc, ispin, istat, max_iloc, max_nmoloc
    LOGICAL :: failure

    failure = .FALSE.

    localized_wfn_control%nloc_states(1:2) = nmoloc(1:2)
    max_nmoloc = max(nmoloc(1),nmoloc(2))


    IF(localized_wfn_control%set_of_states == state_loc_list) THEN
      CPPrecondition(ASSOCIATED(localized_wfn_control%loc_states),cp_failure_level,routineP,error,failure)
      IF(max_nmoloc<SIZE(localized_wfn_control%loc_states,1)) THEN
         CALL reallocate(localized_wfn_control%loc_states,1,max_nmoloc) 
      END IF
      DO ispin = 1,nspins
         localized_wfn_control%lu_bound_states(1,ispin) = 1
         localized_wfn_control%lu_bound_states(2,ispin) = nmoloc(ispin)
         DO i = 1,nmoloc(ispin)  
           iloc = localized_wfn_control%loc_states(i)
           CPPostcondition(iloc<=nmo(ispin),cp_failure_level,routineP,error,failure)
         ENDDO
      END DO 
    ELSE
      IF (ASSOCIATED(localized_wfn_control%loc_states)) THEN
          DEALLOCATE(localized_wfn_control%loc_states,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      IF(localized_wfn_control%set_of_states == state_loc_fromitoj) THEN
        DO ispin = 1,nspins
           localized_wfn_control%lu_bound_states(1,ispin) = &
               localized_wfn_control%lu_bound_states(1,1) 
           localized_wfn_control%lu_bound_states(2,ispin) = &
              localized_wfn_control%lu_bound_states(1,1) +  nmoloc(ispin) -1
           max_iloc = localized_wfn_control%lu_bound_states(2,ispin)
           CPPostcondition(max_iloc<=nmo(ispin),cp_failure_level,routineP,error,failure)
        END DO
      ELSE
        DO ispin = 1,nspins
           localized_wfn_control%lu_bound_states(1,ispin) = 1
           localized_wfn_control%lu_bound_states(2,ispin) = nmoloc(ispin)
        END DO 
      END IF
    END IF 

  END SUBROUTINE set_loc_wfn_lists
 
! ******************************************************************************


END MODULE qs_loc_control
