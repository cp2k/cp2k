!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Defines the control structure, which contain the parameters and the
!>      settings for the calculations of the localization of wavefunctions.
!> \par History
!>      01.2008 Teodoro Laino [tlaino] - University of Zurich
!>        - Merging the two localization codes and updating to new structures
!> \author MI (04.2005)
! *****************************************************************************
MODULE qs_loc_control

  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_loc_crazy, do_loc_direct, do_loc_jacobi, do_loc_l1_norm_sd, &
       do_loc_none, op_loc_berry, op_loc_boys, state_loc_all, state_loc_list, &
       state_loc_none, state_loc_range
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
!> \brief A type that holds controling information for the
!>    calculation of the spread of wfn and theoptimization of
!>    the spread functional
!> \param ref_count:
!> \param localization_method which algorithm is used for the optimization
!> \param operator_type how the spread is defined
!> \param nloc_states number of states on which the spread is computed
!> \param set_of_states how to choose the states
!> \param lu_bound_states lower and upper bounds of the set of states
!>    print_cubes:
!>    print_centers:
!>    print_spreads:
!> \param loc_states list of states on which the spread is computed
!> \param centers_set arrais containing centers and spreads of the selected wfn
!> \param centers_file_name output file names
!> \param spreads_file_name output file names
! *****************************************************************************
  TYPE localized_wfn_control_type
    INTEGER                              :: ref_count
    INTEGER                              :: min_or_max
    INTEGER                              :: localization_method
    INTEGER                              :: operator_type
    INTEGER, DIMENSION(2)                :: nloc_states
    INTEGER                              :: set_of_states
    INTEGER, DIMENSION(2,2)              :: lu_bound_states
    INTEGER                              :: max_iter
    INTEGER                              :: out_each
    REAL(KIND=dp)                        :: eps_localization
    LOGICAL                              :: print_cubes, jacobi_fallback
    LOGICAL                              :: print_centers
    LOGICAL                              :: print_spreads
    INTEGER, POINTER, DIMENSION(:)       :: loc_states
    TYPE(cp_2d_r_p_type), DIMENSION(2)   :: centers_set
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

CONTAINS

! *****************************************************************************
!> \brief create the localized_wfn_control_type
!> \par History
!>      04.2005 created [MI]
! *****************************************************************************
  SUBROUTINE localized_wfn_control_create(localized_wfn_control,error)
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'localized_wfn_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(localized_wfn_control),cp_failure_level,routineP,error,failure)
    ALLOCATE(localized_wfn_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    localized_wfn_control%ref_count= 1
    localized_wfn_control%nloc_states=0
    localized_wfn_control%lu_bound_states=0
    localized_wfn_control%print_cubes = .FALSE.
    localized_wfn_control%print_centers = .FALSE.
    localized_wfn_control%print_spreads = .FALSE.
    NULLIFY(localized_wfn_control%loc_states)
    NULLIFY(localized_wfn_control%centers_set(1)%array)
    NULLIFY(localized_wfn_control%centers_set(2)%array)
  END SUBROUTINE localized_wfn_control_create

! *****************************************************************************
!> \brief release the localized_wfn_control_type
!> \par History
!>      04.2005 created [MI]
! *****************************************************************************
  SUBROUTINE localized_wfn_control_release(localized_wfn_control,error)

    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief retain the localized_wfn_control_type
!> \par History
!>      04.2005 created [MI]
! *****************************************************************************
  SUBROUTINE localized_wfn_control_retain(localized_wfn_control,error)
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'localized_wfn_control_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(localized_wfn_control),cp_failure_level,routineP,error,failure)

    localized_wfn_control%ref_count=localized_wfn_control%ref_count+1
  END SUBROUTINE localized_wfn_control_retain

! *****************************************************************************
!> \brief read the controlparameter from input, using the new input scheme
!> \par History
!>      05.2005 created [MI]
! *****************************************************************************
  SUBROUTINE read_loc_section(localized_wfn_control,loc_section,localize,&
       do_xas,nloc_xas,error)

    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(section_vals_type), POINTER         :: loc_section
    LOGICAL, INTENT(OUT)                     :: localize
    LOGICAL, INTENT(IN), OPTIONAL            :: do_xas
    INTEGER, INTENT(IN), OPTIONAL            :: nloc_xas
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_loc_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ir, n_list, n_rep, &
                                                n_state, output_unit
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: failure, my_do_xas
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    my_do_xas=.FALSE.
    IF(PRESENT(do_xas)) my_do_xas=do_xas
    CPPrecondition(ASSOCIATED(loc_section),cp_failure_level,routineP,error,failure)
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    IF(.NOT. failure)  THEN
       CALL section_vals_val_get(loc_section,"_SECTION_PARAMETERS_",l_val=localize,error=error)
       IF(localize) THEN
          NULLIFY (list)
          localized_wfn_control%lu_bound_states = 0
          localized_wfn_control%nloc_states     = 0
          n_state = 0

          CALL section_vals_val_get(loc_section,"MAX_ITER",&
               i_val=localized_wfn_control%max_iter,error=error)
          CALL section_vals_val_get(loc_section,"OUT_ITER_EACH",&
               i_val=localized_wfn_control%out_each,error=error)
          CALL section_vals_val_get(loc_section,"EPS_LOCALIZATION",&
               r_val=localized_wfn_control%eps_localization,error=error)
          CALL section_vals_val_get(loc_section,"MIN_OR_MAX",&
               i_val=localized_wfn_control%min_or_max,error=error)
          CALL section_vals_val_get(loc_section,"JACOBI_FALLBACK",&
               l_val=localized_wfn_control%jacobi_fallback,error=error)
          CALL section_vals_val_get(loc_section,"METHOD",&
               i_val=localized_wfn_control%localization_method,error=error)
          CALL section_vals_val_get(loc_section,"OPERATOR",&
               i_val=localized_wfn_control%operator_type,error=error)

          ! List of States
          CALL section_vals_val_get(loc_section,"LIST",  n_rep_val=n_rep,error=error)
          n_list = 0
          DO ir = 1,n_rep
             NULLIFY(list)
             CALL section_vals_val_get(loc_section,"LIST",i_rep_val=ir,i_vals=list,error=error)
             IF(ASSOCIATED(list)) THEN
                CALL reallocate(localized_wfn_control%loc_states,1,n_list+SIZE(list))
                DO i = 1, SIZE(list)
                   localized_wfn_control%loc_states(n_list+i) = list(i)
                END DO ! i
                localized_wfn_control%nloc_states(1) = n_list + SIZE(list)
                n_list = n_list +  SIZE(list)
             END IF
          END DO ! ir
          IF(n_list/=0) localized_wfn_control%set_of_states=state_loc_list

          ! All States or XAS specific states
          IF(localized_wfn_control%nloc_states(1) == 0)  THEN
             IF(my_do_xas) THEN
                localized_wfn_control%set_of_states = state_loc_range
                localized_wfn_control%nloc_states(1) =  nloc_xas
                localized_wfn_control%lu_bound_states(1,1) = 1
                localized_wfn_control%lu_bound_states(2,1) = nloc_xas
                localized_wfn_control%nloc_states(2) =  0
                localized_wfn_control%lu_bound_states(1,2) = 0
                localized_wfn_control%lu_bound_states(2,2) = 0
             ELSE
                localized_wfn_control%set_of_states = state_loc_all
             END IF
          END IF

          localized_wfn_control%print_centers = &
            BTEST(cp_print_key_should_output(logger%iter_info,loc_section,&
                  "PRINT%WANNIER_CENTERS",error=error),cp_p_file)
          localized_wfn_control%print_spreads = &
            BTEST(cp_print_key_should_output(logger%iter_info,loc_section,&
                  "PRINT%WANNIER_SPREADS",error=error),cp_p_file)
          localized_wfn_control%print_cubes = &
            BTEST(cp_print_key_should_output(logger%iter_info,loc_section,&
                  "PRINT%WANNIER_CUBES",error=error),cp_p_file)

          output_unit = cp_print_key_unit_nr(logger,loc_section,&
                                             "PRINT%PROGRAM_RUN_INFO",&
                                             extension=".Log",&
                                             error=error)

          IF (output_unit > 0) THEN

            WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
              "LOCALIZE| The spread relative to a set of orbitals is computed"

            SELECT CASE (localized_wfn_control%set_of_states)
            CASE (state_loc_all)
              WRITE (UNIT=output_unit,FMT="(T2,A)")&
                "LOCALIZE| Orbitals to be localized: All the occupied"
            CASE (state_loc_range)
              WRITE (UNIT=output_unit,FMT="(T2,A,T65,2I8)")&
                "LOCALIZE| Orbitals to be localized: Those with index between ",&
                localized_wfn_control%lu_bound_states(1,1),&
                localized_wfn_control%lu_bound_states(2,1)
            CASE (state_loc_list)
              WRITE (UNIT=output_unit,FMT="(T2,A)")&
                "LOCALIZE| Orbitals to be localized: Those with index in the list"
                WRITE (UNIT=output_unit,FMT=*) localized_wfn_control%loc_states
            CASE DEFAULT
              WRITE (UNIT=output_unit,FMT="(T2,A)")&
                "LOCALIZE| Orbitals to be localized: None "
            END SELECT

            SELECT CASE (localized_wfn_control%operator_type)
            CASE(op_loc_berry)
              WRITE (UNIT=output_unit,FMT="(T2,A)")&
                "LOCALIZE| Spread defined by the Berry phase operator "
            CASE(op_loc_boys)
              WRITE (UNIT=output_unit,FMT="(T2,A)")&
                "LOCALIZE| Spread defined by the Boys phase operator "
            CASE DEFAULT
              WRITE (UNIT=output_unit,FMT="(T2,A)")&
                "LOCALIZE| Spread defined by the Pipek phase operator "
            END SELECT

            SELECT CASE (localized_wfn_control%localization_method)
            CASE (do_loc_jacobi)
              WRITE (UNIT=output_unit,FMT="(T2,A)")&
                "LOCALIZE| The optimal unitary transformation is generated by the Jacobi algorithm"
            CASE (do_loc_crazy)
              WRITE (UNIT=output_unit,FMT="(T2,A)")&
                "LOCALIZE| The optimal unitary transformation is gnerated by the Crazy angle algorithm"
            CASE (do_loc_direct)
              WRITE (UNIT=output_unit,FMT="(T2,A)")&
                "LOCALIZE| The optimal unitary transformation is gnerated by the direct algorithm"
            CASE (do_loc_l1_norm_sd)
              WRITE (UNIT=output_unit,FMT="(T2,A)")&
                "LOCALIZE| The optimal unitary transformation is gnerated by the "//&
                " steepest descent algorithm applyed on an approximate l1 norm"
            CASE (do_loc_none)
              WRITE (UNIT=output_unit,FMT="(T2,A)")&
                "LOCALIZE| No unitary transformation is applied"
            END SELECT

          END IF ! process has output_unit

          CALL cp_print_key_finished_output(output_unit,logger,loc_section,&
                                            "PRINT%PROGRAM_RUN_INFO",error=error)

       ELSE

          localized_wfn_control%localization_method = do_loc_none
          localized_wfn_control%localization_method = state_loc_none
          localized_wfn_control%print_centers  = .FALSE.
          localized_wfn_control%print_spreads  = .FALSE.
          localized_wfn_control%print_cubes  = .FALSE.

       END IF

    END IF

  END SUBROUTINE read_loc_section

! *****************************************************************************
!> \brief create the center and spread array and the file names for the output
!> \note
!>      would be better to get rid of globenv
!> \par History
!>      04.2005 created [MI]
! *****************************************************************************
  SUBROUTINE set_loc_centers(localized_wfn_control,nmoloc,nspins,error)

    TYPE(localized_wfn_control_type)         :: localized_wfn_control
    INTEGER, DIMENSION(2), INTENT(IN)        :: nmoloc
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_loc_centers', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    DO ispin = 1,nspins
       ALLOCATE(localized_wfn_control%centers_set(ispin)%array(6,nmoloc(ispin)), STAT = istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       localized_wfn_control%centers_set(ispin)%array = 0.0_dp
    END DO

  END SUBROUTINE set_loc_centers

! *****************************************************************************
!> \brief create the lists of mos that are taken into account
!> \par History
!>      04.2005 created [MI]
! *****************************************************************************
  SUBROUTINE set_loc_wfn_lists(localized_wfn_control,nmoloc,nmo,nspins,error)

    TYPE(localized_wfn_control_type)         :: localized_wfn_control
    INTEGER, DIMENSION(2), INTENT(IN)        :: nmoloc, nmo
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_loc_wfn_lists', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iloc, ispin, istat, &
                                                max_iloc, max_nmoloc
    LOGICAL                                  :: failure

    failure = .FALSE.

    localized_wfn_control%nloc_states(1:2) = nmoloc(1:2)
    max_nmoloc = MAX(nmoloc(1),nmoloc(2))

    IF(localized_wfn_control%set_of_states == state_loc_list) THEN
       ! List
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
       IF(localized_wfn_control%set_of_states == state_loc_range) THEN
          ! Range
          DO ispin = 1,nspins
             localized_wfn_control%lu_bound_states(1,ispin) = &
                  localized_wfn_control%lu_bound_states(1,1)
             localized_wfn_control%lu_bound_states(2,ispin) = &
                  localized_wfn_control%lu_bound_states(1,1) +  nmoloc(ispin) -1
             max_iloc = localized_wfn_control%lu_bound_states(2,ispin)
             CPPostcondition(max_iloc<=nmo(ispin),cp_failure_level,routineP,error,failure)
          END DO
       ELSE
          ! All
          ALLOCATE(localized_wfn_control%loc_states(max_nmoloc),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO ispin = 1,nspins
             localized_wfn_control%lu_bound_states(1,ispin) = 1
             localized_wfn_control%lu_bound_states(2,ispin) = nmoloc(ispin)
             IF(nmoloc(ispin)< 1) localized_wfn_control%lu_bound_states(1,ispin) = 0
             DO i = 1, nmoloc(ispin)
                localized_wfn_control%loc_states(i) = i
             END DO
          END DO
       END IF
    END IF

  END SUBROUTINE set_loc_wfn_lists

END MODULE qs_loc_control
