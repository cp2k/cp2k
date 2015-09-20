!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Defines control structures, which contain the parameters and the
!>      settings for the calculations.
! *****************************************************************************
MODULE xas_control

  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_constants,                 ONLY: xas_1s_type,&
                                             xas_dscf,&
                                             xas_tp_fh,&
                                             xas_tp_flex,&
                                             xas_tp_hh,&
                                             xas_tp_xfh,&
                                             xas_tp_xhh,&
                                             xes_tp_val
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
!> \brief A type that holds controling information for a xas calculation
! *****************************************************************************
  TYPE xas_control_type
    INTEGER                             :: ref_count
    INTEGER                             :: nexc_atoms
    INTEGER                             :: nexc_search
    INTEGER                             :: state_type
    INTEGER                             :: xas_method
    INTEGER                             :: dipole_form
    INTEGER                             :: added_mos
    INTEGER                             :: max_iter_added
    INTEGER                             :: ngauss
    INTEGER                             :: stride
    INTEGER, DIMENSION(:), POINTER      :: exc_atoms
    LOGICAL                             :: cubes, do_centers
    LOGICAL                             :: xas_restart
    INTEGER, DIMENSION(:), POINTER      :: list_cubes
!
    REAL(dp)                            :: eps_added
    REAL(dp)                            :: xes_core_occupation
    REAL(dp)                            :: xes_homo_occupation
    REAL(dp)                            :: nel_tot, xas_core_occupation
  END TYPE xas_control_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_control'

! *** Public data types ***

  PUBLIC :: xas_control_type

! *** Public subroutines ***

  PUBLIC :: read_xas_control, write_xas_control, xas_control_create,&
            xas_control_release

CONTAINS

! *****************************************************************************
!> \brief read from input the instructions for a xes/xas calculation
!> \param xas_control control variables
!>       error
!> \param xas_section ...
!> \par History
!>      04.2005 created [MI]
! *****************************************************************************
  SUBROUTINE read_xas_control(xas_control,xas_section)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(section_vals_type), POINTER         :: xas_section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_xas_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ir, istat, n_rep, nex_at
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: failure, hempty, was_present

    failure=.FALSE.
    was_present = .FALSE.


    NULLIFY(list)

    CALL section_vals_val_get(xas_section,"METHOD",&
         i_val=xas_control%xas_method)

    CALL section_vals_val_get(xas_section,"DIPOLE_FORM",&
         i_val=xas_control%dipole_form)

    CALL section_vals_val_get(xas_section,"RESTART",&
         l_val=xas_control%xas_restart)

    CALL section_vals_val_get(xas_section,"STATE_TYPE",&
         i_val=xas_control%state_type)

    CALL section_vals_val_get(xas_section,"STATE_SEARCH",&
         i_val=xas_control%nexc_search)

    CALL section_vals_val_get(xas_section,"XAS_CORE", &
         r_val=xas_control%xas_core_occupation)

    CALL section_vals_val_get(xas_section,"XAS_TOT_EL", &
         r_val=xas_control%nel_tot)

    CALL section_vals_val_get(xas_section,"XES_CORE", &
         r_val=xas_control%xes_core_occupation)

    CALL section_vals_val_get(xas_section,"XES_EMPTY_HOMO", &
         l_val=hempty)
    IF(hempty) THEN
         xas_control%xes_homo_occupation = 0
    ELSE
         xas_control%xes_homo_occupation = 1
    END IF

! It should be further generalized
    IF(.NOT. ASSOCIATED(xas_control%exc_atoms)) THEN
      CALL section_vals_val_get(xas_section,"ATOMS_LIST",&
          n_rep_val=n_rep)

      IF(n_rep > 0) THEN
        nex_at = 0
        DO ir = 1,n_rep
        NULLIFY(list)
        CALL section_vals_val_get(xas_section,"ATOMS_LIST",&
             i_rep_val=ir,i_vals=list)

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
    END IF

    IF(.NOT. ASSOCIATED(xas_control%exc_atoms)) THEN
      xas_control%nexc_atoms = 1
      ALLOCATE(xas_control%exc_atoms(1),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,failure)
      xas_control%exc_atoms(1) = 1
    END IF

    CALL section_vals_val_get(xas_section,"ADDED_MOS",&
         i_val=xas_control%added_mos)

    CALL section_vals_val_get(xas_section,"MAX_ITER_ADDED",&
         i_val=xas_control%max_iter_added)

    CALL section_vals_val_get(xas_section,"EPS_ADDED", &
         r_val=xas_control%eps_added)

    CALL section_vals_val_get(xas_section,"NGAUSS",&
         i_val=xas_control%ngauss)


  END SUBROUTINE read_xas_control

! *****************************************************************************
!> \brief write on the instructions for a xes/xas calculation
!> \param xas_control control variables
!>       error
!> \param dft_section ...
!> \par History
!>      12.2005 created [MI]
! *****************************************************************************
  SUBROUTINE write_xas_control(xas_control,dft_section)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(section_vals_type), POINTER         :: dft_section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_xas_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: output_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_get_default_logger()
    output_unit = cp_print_key_unit_nr(logger,dft_section,&
         "PRINT%DFT_CONTROL_PARAMETERS",extension=".Log")
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
        CASE(xas_tp_flex)
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T25,A)")&
                "XAS| Method:",&
                   "      Transition potential with occupation of core state given from input"
        CASE(xas_dscf)
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T40,A)")&
                "XAS| Method:",&
                   "         DSCF for the first excited state"
        CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unknown xas method "//TRIM(ADJUSTL(cp_to_string(xas_control%xas_method))),&
             failure)
      END SELECT
      IF(xas_control%xas_restart) THEN
          WRITE (UNIT=output_unit,FMT="(/,T2,A,T30,A)")&
                "XAS|"," Orbitals read from atom-specific restart file when available"
      END IF
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,dft_section,&
         "PRINT%DFT_CONTROL_PARAMETERS")
  END SUBROUTINE write_xas_control

!****f* xas_control/xas_control_create/retain/release *

! *****************************************************************************
!> \brief create retain release the xas_control_type
!> \param xas_control ...
!> \par History
!>      04.2005 created [MI]
! *****************************************************************************
  SUBROUTINE xas_control_create(xas_control)

    TYPE(xas_control_type), POINTER          :: xas_control

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(xas_control),cp_failure_level,routineP,failure)
    ALLOCATE(xas_control,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    xas_control%ref_count= 1
    xas_control%xas_method = xas_tp_hh
    xas_control%nexc_atoms = 1
    xas_control%nexc_search = -1
    xas_control%state_type = xas_1s_type
    xas_control%xas_restart = .FALSE.
    xas_control%added_mos=0
    xas_control%xes_core_occupation=1.0_dp
    xas_control%xes_homo_occupation=1.0_dp
    NULLIFY(xas_control%exc_atoms)
    xas_control%cubes = .FALSE.
    xas_control%do_centers = .FALSE.
    NULLIFY(xas_control%list_cubes)

  END SUBROUTINE xas_control_create

! *****************************************************************************
!> \brief ...
!> \param xas_control ...
! *****************************************************************************
  SUBROUTINE xas_control_release(xas_control)

    TYPE(xas_control_type), POINTER          :: xas_control

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF(ASSOCIATED(xas_control)) THEN
      CPPrecondition(xas_control%ref_count>0,cp_failure_level,routineP,failure)
      xas_control%ref_count=xas_control%ref_count-1
      IF (xas_control%ref_count==0) THEN
        IF (ASSOCIATED(xas_control%exc_atoms)) THEN
          DEALLOCATE(xas_control%exc_atoms,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,failure)
        END IF
        IF (ASSOCIATED(xas_control%list_cubes)) THEN
          DEALLOCATE(xas_control%list_cubes,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,failure)
        END IF
        DEALLOCATE(xas_control,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,failure)
      END IF
    END IF

  END SUBROUTINE xas_control_release

! *****************************************************************************
!> \brief ...
!> \param xas_control ...
! *****************************************************************************
  SUBROUTINE xas_control_retain(xas_control)

    TYPE(xas_control_type), POINTER          :: xas_control

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_control_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(xas_control),cp_failure_level,routineP,failure)
    xas_control%ref_count=xas_control%ref_count+1

  END SUBROUTINE xas_control_retain

END MODULE xas_control
