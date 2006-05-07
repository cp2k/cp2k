!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!


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
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_should_output
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       do_loc_crazy, do_loc_direct, do_loc_jacobi, do_loc_max, do_loc_min, &
       do_loc_none, op_loc_berry, op_loc_boys, op_loc_pipek, state_loc_all, &
       state_loc_fromitoj, state_loc_list, state_loc_none, state_loc_upton
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE memory_utilities,                ONLY: reallocate
  USE string_utilities,                ONLY: xstring
  USE kinds,                           ONLY: default_string_length,&
                                             default_path_length
#include "cp_common_uses.h"

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
    CHARACTER ( LEN = default_path_length )               :: loc_rootfname
    CHARACTER ( LEN = default_path_length ),DIMENSION(2):: centers_file_name
    CHARACTER ( LEN = default_path_length ),DIMENSION(2):: spreads_file_name
  END TYPE localized_wfn_control_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_loc_control'

! *** Public data types ***

  PUBLIC :: localized_wfn_control_type

! *** Public subroutines ***

  PUBLIC :: localized_wfn_control_create, &
            localized_wfn_control_release,&
            localized_wfn_control_retain,&
            read_loc_section,&
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
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'localized_wfn_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

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
   localized_wfn_control%print_cubes = .FALSE.
   localized_wfn_control%print_centers = .FALSE.
   localized_wfn_control%print_spreads = .FALSE.
   NULLIFY(localized_wfn_control%loc_states)
   NULLIFY(localized_wfn_control%centers_set(1)%array)
   NULLIFY(localized_wfn_control%centers_set(2)%array)

  END SUBROUTINE localized_wfn_control_create

  SUBROUTINE localized_wfn_control_release(localized_wfn_control,error)

    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'localized_wfn_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

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
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'localized_wfn_control_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

   failure=.FALSE.
   CPPrecondition(ASSOCIATED(localized_wfn_control),cp_failure_level,routineP,error,failure)

   localized_wfn_control%ref_count=localized_wfn_control%ref_count+1
  END SUBROUTINE localized_wfn_control_retain

!!****f* qs_loc_control/read_loc_section *
!!
!!   NAME
!!      read_loc_section
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

  SUBROUTINE read_loc_section(localized_wfn_control,loc_section,error)

    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(section_vals_type), POINTER         :: loc_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_loc_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat, output_unit
    INTEGER, DIMENSION(:), POINTER           :: bounds, list
    LOGICAL                                  :: failure, ionode
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(loc_section),cp_failure_level,routineP,error,failure)
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
      output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
      output_unit = -1
    END IF

    IF(.NOT. failure)  THEN

       NULLIFY (bounds, list ) 

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
         CALL section_vals_val_get(loc_section,"UP_BOUND",&
              i_val=localized_wfn_control%lu_bound_states(2,1), error=error)
         localized_wfn_control%lu_bound_states(1,1) = 1
         localized_wfn_control%nloc_states(1) = &
                     localized_wfn_control%lu_bound_states(2,1) 
      ELSE IF(localized_wfn_control%set_of_states==state_loc_fromitoj) THEN
         CALL section_vals_val_get(loc_section,"LU_BOUNDS",&
              i_vals=bounds,error=error)
              localized_wfn_control%lu_bound_states(1,1) = bounds(1)
              localized_wfn_control%lu_bound_states(2,1) = bounds(2)

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

           
      IF(ionode) THEN

        WRITE (UNIT=output_unit,FMT="(/,T3,A)")&
                 "Localization of a subset of KS orbitals"
        IF(localized_wfn_control%min_or_max==do_loc_min) &
           WRITE (UNIT=output_unit,FMT="(T2,A,T54,A)")&
          "LOC|","Minimize the orbital spread"
        IF(localized_wfn_control%min_or_max==do_loc_max) &
            WRITE (UNIT=output_unit,FMT="(T2,A,T54,A)")&
          "LOC|","Maximize the orbital spread"
        IF(localized_wfn_control%operator_type==op_loc_berry)  &
           WRITE (UNIT=output_unit,FMT="(T2,A,T65,A)")&
          "LOC|  Definition of the spread functional by ","     Berry Phase"
        IF(localized_wfn_control%operator_type==op_loc_pipek)  &
           WRITE (UNIT=output_unit,FMT="(T2,A,T65,A)")&
          "LOC|  Definition of the spread functional by ","Pipek definition"
        IF(localized_wfn_control%operator_type==op_loc_boys)  &
           WRITE (UNIT=output_unit,FMT="(T2,A,T65,A)")&
          "LOC|  Definition of the spread functional by "," Boys definition"
        IF(localized_wfn_control%operator_type==do_loc_jacobi)  &
           WRITE (UNIT=output_unit,FMT="(T2,A,T63,A)")&
          "LOC|  Optimal Rotation found by ",             "     Jacobi method"
        IF(localized_wfn_control%operator_type==do_loc_crazy)  &
           WRITE (UNIT=output_unit,FMT="(T2,A,T63,A)")&
          "LOC|  Optimal Rotation found by ","Crazy angle method"
        IF(localized_wfn_control%operator_type==do_loc_direct)  &
           WRITE (UNIT=output_unit,FMT="(T2,A,T63,A)")&
          "LOC|  Optimal Rotation found by ","  Direct SD method"
        IF(localized_wfn_control%operator_type==do_loc_none)  &
           WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
          "LOC|   ","No stationary state is located"
        IF(localized_wfn_control%set_of_states==state_loc_all) &
           WRITE (UNIT=output_unit,FMT="(T2,A,T57,A)")&
          "LOC|  Set of states: " , "All the occupied orbitals"
        IF(localized_wfn_control%set_of_states==state_loc_upton) &
           WRITE (UNIT=output_unit,FMT="(T2,A,T53,A,I6)")&
          "LOC|  Set of states: " , "All the orbitals up to",&
           localized_wfn_control%nloc_states(1)
        IF(localized_wfn_control%set_of_states==state_loc_fromitoj) &
           WRITE (UNIT=output_unit,FMT="(T2,A,T46,A,I6,A,I6)")&
          "LOC|  Set of states: " , "The orbitals between",&
           localized_wfn_control%lu_bound_states(1,1), " and ",&
           localized_wfn_control%lu_bound_states(2,1)
        IF(localized_wfn_control%set_of_states==state_loc_list) &
           WRITE (UNIT=output_unit,FMT="(T2,A,T50,A,I6,A)")&
          "LOC|  Set of states: " , "List of ",&
          localized_wfn_control%nloc_states(1), " orbitals from input"
      END IF
    END IF

  END SUBROUTINE read_loc_section

!!****f* qs_loc_control/set_loc_centers *
!!
!!   NAME
!!      set_loc_centers
!!   FUNCTION
!!      create the center and spread array and the file names for the output
!!   NOTES
!!     would be better to get rid of globenv
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************

  SUBROUTINE set_loc_centers(localized_wfn_control,nmoloc,nspins,globenv,error)

    TYPE(localized_wfn_control_type)         :: localized_wfn_control
    INTEGER, DIMENSION(2), INTENT(IN)        :: nmoloc
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_loc_centers', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ie, ispin, istat
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
 
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
      DO ispin = 1,nspins
         localized_wfn_control%spreads_file_name(ispin) = &
                 localized_wfn_control%loc_rootfname(ia:ie)//&
                 '_spreads_s'//TRIM(ADJUSTL(cp_to_string(ispin)))
         localized_wfn_control%centers_file_name(ispin) = &
                 localized_wfn_control%loc_rootfname(ia:ie)//&
                 '_centers_s'//TRIM(ADJUSTL(cp_to_string(ispin)))
       END DO 
    END IF

    DO ispin = 1,nspins
      ALLOCATE(localized_wfn_control%centers_set(ispin)%array(6,nmoloc(ispin)), STAT = istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END DO

    localized_wfn_control%print_centers = BTEST(cp_print_key_should_output(logger%iter_info,globenv%input_file,&
         "FORCE_EVAL%DFT%PRINT%LOCALIZATION",error=error),cp_p_file)
    localized_wfn_control%print_spreads = BTEST(cp_print_key_should_output(logger%iter_info,globenv%input_file,&
         "FORCE_EVAL%DFT%PRINT%LOCALIZATION",error=error),cp_p_file)
    localized_wfn_control%print_cubes =   BTEST(cp_print_key_should_output(logger%iter_info,globenv%input_file,&
         "FORCE_EVAL%DFT%PRINT%LOCALIZATION%CUBES",error=error),cp_p_file)

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

    TYPE(localized_wfn_control_type)         :: localized_wfn_control
    INTEGER, DIMENSION(2), INTENT(IN)        :: nmoloc, nmo
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_loc_wfn_lists', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iloc, ispin, istat, &
                                                max_iloc, max_nmoloc
    LOGICAL                                  :: failure

    failure = .FALSE.

    localized_wfn_control%nloc_states(1:2) = nmoloc(1:2)
    max_nmoloc = MAX(nmoloc(1),nmoloc(2))


    IF(localized_wfn_control%set_of_states == state_loc_list) THEN
      CPPrecondition(ASSOCIATED(localized_wfn_control%loc_states),cp_failure_level,routineP,error,failure)
      IF(max_nmoloc<SIZE(localized_wfn_control%loc_states,1)) THEN
         CALL reallocate(localized_wfn_control%loc_states,1,max_nmoloc) 
      END IF
      DO ispin = 1,nspins
         localized_wfn_control%lu_bound_states(1,ispin) = 1
         localized_wfn_control%lu_bound_states(2,ispin) = nmoloc(ispin)
         IF(nmoloc(ispin)< 1) localized_wfn_control%lu_bound_states(1,ispin) = 0
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
           IF(nmoloc(ispin)< 1) localized_wfn_control%lu_bound_states(1,ispin) = 0
        END DO 
      END IF
    END IF 

  END SUBROUTINE set_loc_wfn_lists
 
! ******************************************************************************


END MODULE qs_loc_control
