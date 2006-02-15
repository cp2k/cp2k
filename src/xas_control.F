!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!


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

  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: xas_1s_type,&
                                             xas_dscf,&
                                             xas_tddft,&
                                             xas_tp_fh,&
                                             xas_tp_hh,&
                                             xas_tp_xfh,&
                                             xas_tp_xhh,&
                                             xes_tp_val
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE qs_loc_control,                  ONLY: localized_wfn_control_create,&
                                             localized_wfn_control_release,&
                                             localized_wfn_control_type,&
                                             read_loc_section
#include "cp_common_uses.h"

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
    INTEGER                             :: scf_method
    INTEGER                             :: nexc_atoms
    INTEGER                             :: nexc_search
    INTEGER                             :: state_type
    INTEGER                             :: xas_method
    INTEGER                             :: dipole_form
    INTEGER                             :: lb_atom, ub_atom
    INTEGER                             :: added_mos
    INTEGER                             :: max_iter_added
    INTEGER                             :: ngauss
    INTEGER                             :: stride
    INTEGER, DIMENSION(:), POINTER      :: exc_atoms 
    LOGICAL                             :: emission
    LOGICAL                             :: localized_wfn
    LOGICAL                             :: cubes, do_centers
    LOGICAL                             :: fermidist
    INTEGER, DIMENSION(:), POINTER      :: list_cubes
    REAL(dp)                            :: eps_scf
    REAL(dp)                            :: eps_diis
    REAL(dp)                            :: p_mix
    REAL(dp)                            :: eps_added
    REAL(dp)                            :: smear
    REAL(dp)                            :: xes_core_occupation
    REAL(dp)                            :: occ_estate, occ_homo 
    TYPE(localized_wfn_control_type), &
      POINTER                           :: localized_wfn_control
  END TYPE xas_control_type


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_control'

! *** Public data types ***

  PUBLIC :: xas_control_type

! *** Public subroutines ***

  PUBLIC :: read_xas_control, write_xas_control, xas_control_create, &
            xas_control_release, xas_control_retain

!!***
! *****************************************************************************

CONTAINS

!   ***************************************************************************
!!
!!   NAME
!!      read_xas_control
!!   FUNCTION
!!      read from input the instructions for a xes/xas calculation
!!   NOTES
!!
!!   INPUTS
!!      xas_control :  control variables
!!      globenv : global variables
!!      error
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE read_xas_control(xas_control,globenv, error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(global_environment_type)            :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_xas_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ir, istat, n_rep, ncubes, &
                                                nex_at
    INTEGER, DIMENSION(:), POINTER           :: bounds, list
    LOGICAL                                  :: failure, was_present
    TYPE(section_vals_type), POINTER         :: loc_section, xas_section

    failure=.FALSE.
    was_present = .FALSE.

    NULLIFY(xas_section, loc_section)
    xas_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL",error=error)
    xas_section => section_vals_get_subs_vals(xas_section,"DFT",error=error)
    xas_section => section_vals_get_subs_vals(xas_section,"XAS",error=error)

    NULLIFY(bounds,list)

    CALL section_vals_val_get(xas_section,"CONVERGENCE", &
         r_val=xas_control%eps_scf,error=error)

    CALL section_vals_val_get(xas_section,"EPS_DIIS", &
         r_val=xas_control%eps_diis,error=error)

    CALL section_vals_val_get(xas_section,"MIXING", &
         r_val=xas_control%p_mix,error=error)

    CALL section_vals_val_get(xas_section,"METHOD",&
         i_val=xas_control%xas_method,error=error)

    CALL section_vals_val_get(xas_section,"DIPOLE_FORM",&
         i_val=xas_control%dipole_form,error=error)

    CALL section_vals_val_get(xas_section,"MAXSTEP",&
         i_val=xas_control%max_scf,error=error)

    CALL section_vals_val_get(xas_section,"SCF_OPTIMIZER",&
         i_val=xas_control%scf_method,error=error)

    CALL section_vals_val_get(xas_section,"EACH_STEP",&
         i_val=xas_control%each_step,error=error)

    CALL section_vals_val_get(xas_section,"STATE_TYPE",&
         i_val=xas_control%state_type,error=error)

    CALL section_vals_val_get(xas_section,"STATE_SEARCH",&
         i_val=xas_control%nexc_search,error=error)

    CALL section_vals_val_get(xas_section,"FNAME",&
         c_val=xas_control%xas_rootfname,error=error)

    CALL section_vals_val_get(xas_section,"SMEAR", &
         r_val=xas_control%smear,error=error)

    CALL section_vals_val_get(xas_section,"FERMI", &
         l_val=xas_control%fermidist,error=error)

    CALL section_vals_val_get(xas_section,"XES_CORE", &
         r_val=xas_control%xes_core_occupation,error=error)

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

! It should be further generalized
    IF(.NOT. ASSOCIATED(xas_control%exc_atoms)) THEN
      CALL section_vals_val_get(xas_section,"ATOMS_LIST",&
          n_rep_val=n_rep, error=error)
      
      nex_at = 0
      DO ir = 1,n_rep
        NULLIFY(list)
        CALL section_vals_val_get(xas_section,"ATOMS_LIST",&
             i_rep_val=ir,i_vals=list,error=error)

        IF(ASSOCIATED(list)) THEN
          CALL reallocate(xas_control%exc_atoms,1,nex_at+SIZE(list))
          DO i = 1,SIZE(list)
            xas_control%exc_atoms(i+nex_at) = list(i)

          END DO
          xas_control%nexc_atoms =  nex_at + SIZE(list)
          nex_at = nex_at + SIZE(list)
        END IF
      END DO ! ir
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

    CALL section_vals_val_get(xas_section,"ADDED_MOS",&
         i_val=xas_control%added_mos,error=error)

    CALL section_vals_val_get(xas_section,"MAX_ITER_ADDED",&
         i_val=xas_control%max_iter_added,error=error)

    CALL section_vals_val_get(xas_section,"EPS_ADDED", &
         r_val=xas_control%eps_added,error=error)

    CALL section_vals_val_get(xas_section,"NGAUSS",&
         i_val=xas_control%ngauss,error=error)

    CALL section_vals_val_get(xas_section,"EMISSION",&
         l_val=xas_control%emission,error=error)
         
    CALL section_vals_val_get(xas_section,"LOCALIZE_WFN",&
         l_val=xas_control%localized_wfn,error=error)

    loc_section =>section_vals_get_subs_vals(xas_section,"LOCALIZE",error=error)
    CALL section_vals_get(loc_section,explicit=was_present,error=error)
    IF (was_present) THEN
      CALL read_loc_section(xas_control%localized_wfn_control,loc_section,error=error)
    END IF

    NULLIFY(bounds,list)
    CALL section_vals_val_get(xas_section,"CUBES",&
         l_val=xas_control%cubes,error=error)
    IF(xas_control%cubes) THEN
      CALL section_vals_val_get(xas_section,"CUBES_LU_BOUNDS",&
           i_vals=bounds,error=error) 
      ncubes = bounds(2) - bounds(1)  + 1
      IF(ncubes > 0 ) THEN
        ALLOCATE( xas_control%list_cubes(ncubes),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        DO i = 1,ncubes
          xas_control%list_cubes(i) = bounds(1) + (i-1)
        END DO
      END IF

! It should be further generalized
      IF(.NOT. ASSOCIATED(xas_control%list_cubes)) THEN
        CALL section_vals_val_get(xas_section,"CUBES_LIST",&
             n_rep_val=n_rep,error=error)
        ncubes = 0
        DO ir = 1,n_rep
          NULLIFY(list)
          CALL section_vals_val_get(xas_section,"CUBES_LIST",&
             i_rep_val=ir,i_vals=list,error=error)
          IF(ASSOCIATED(list)) THEN
            CALL reallocate(xas_control%list_cubes,1,ncubes+ SIZE(list))
            DO i = 1, SIZE(list)
              xas_control%list_cubes(i+ncubes) = list(i)
            END DO
            ncubes = ncubes + SIZE(list)
          END IF
        END DO  ! ir

      END IF

      CALL section_vals_val_get(xas_section,"CUBES_STRIDE",&
           i_val=xas_control%stride,error=error)
      ! If the states are not given from input, don't print any
      IF(.NOT. ASSOCIATED(xas_control%list_cubes)) xas_control%cubes = .FALSE.
    ELSE
      NULLIFY(xas_control%list_cubes)
    ENDIF

  END SUBROUTINE read_xas_control
!   ***************************************************************************
!!
!!   NAME
!!      write_xas_control
!!   FUNCTION
!!      write on the instructions for a xes/xas calculation
!!   NOTES
!!
!!   INPUTS
!!      xas_control :  control variables
!!      globenv : global variables
!!      error
!!
!!   MODIFICATION HISTORY
!!     12.2005 created [MI]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE write_xas_control(xas_control,dft_section, error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_xas_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: output_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

!   ---------------------------------------------------------------------------

    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,dft_section,&
         "PRINT%DFT_CONTROL_PARAMETERS",extension=".Log",error=error)
    IF (output_unit>0) THEN
      SELECT CASE(xas_control%xas_method)
        CASE(xas_tp_hh)
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T40,A)")&
                "XAS| Method:",&
                   "      Transition potential with half hole" 
        CASE(xas_tp_xhh)
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T40,A)")&
                "XAS| Method:",&
                   "      Transition potential with excited half hole" 
        CASE(xas_tp_fh)
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T40,A)")&
                "XAS| Method:",&
                   "      Transition potential with full hole" 
        CASE(xas_tp_xfh)
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T40,A)")&
                "XAS| Method:",&
                   "      Transition potential with excited full hole" 
        CASE(xes_tp_val)
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T40,A)")&
                "XAS| Method:",&
                   " Only XES with full core and hole in lumo" 
        CASE(xas_tddft)
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T40,A)")&
                "XAS| Method:",&
                   "                                    TDDFT" 
        CASE(xas_dscf)
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T40,A)")&
                "XAS| Method:",&
                   "         DSCF for the first excited state" 
        CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown xas method "//TRIM(ADJUSTL(cp_to_string(xas_control%xas_method))),&
             error,failure)
      END SELECT
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,dft_section,&
         "PRINT%DFT_CONTROL_PARAMETERS",error=error)
  END SUBROUTINE write_xas_control

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

    TYPE(xas_control_type), POINTER          :: xas_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(xas_control),cp_failure_level,routineP,error,failure)
    ALLOCATE(xas_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    xas_control%ref_count= 1
    xas_control%each_step= 1
    xas_control%xas_method = xas_tp_hh
    xas_control%occ_estate = 0.5_dp
    xas_control%occ_homo = 1.0_dp
    xas_control%xas_rootfname = ""
    xas_control%nexc_atoms = 1
    xas_control%nexc_search = -1
    xas_control%lb_atom = 1
    xas_control%ub_atom = 1
    xas_control%state_type = xas_1s_type    
    xas_control%max_scf = 150
    xas_control%eps_scf = 0.5d-6
    xas_control%localized_wfn = .FALSE. 
    xas_control%added_mos=0
    xas_control%xes_core_occupation=1.0_dp
    NULLIFY(xas_control%localized_wfn_control)
    CALL localized_wfn_control_create(xas_control%localized_wfn_control,error=error)
    NULLIFY(xas_control%exc_atoms)
    xas_control%cubes = .FALSE.
    xas_control%do_centers = .FALSE.
    NULLIFY(xas_control%list_cubes)

  END SUBROUTINE xas_control_create

  SUBROUTINE xas_control_release(xas_control,error)

    TYPE(xas_control_type), POINTER          :: xas_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

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
        IF (ASSOCIATED(xas_control%list_cubes)) THEN
          DEALLOCATE(xas_control%list_cubes,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        DEALLOCATE(xas_control,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
    END IF

  END SUBROUTINE xas_control_release

  SUBROUTINE xas_control_retain(xas_control,error)

    TYPE(xas_control_type), POINTER          :: xas_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_control_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(xas_control),cp_failure_level,routineP,error,failure)
    xas_control%ref_count=xas_control%ref_count+1

  END SUBROUTINE xas_control_retain

END MODULE xas_control
