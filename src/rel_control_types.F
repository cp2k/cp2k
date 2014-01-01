!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief parameters that control a relativistic calculation
!> \par History
!>      09.2002 created [fawzi] (as scf_control_types.F)
!>      10.2008 modifed for relativistic control types (Jens Thar)
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE rel_control_types

  USE f77_blas
  USE input_constants,                 ONLY: rel_none,&
                                             rel_pot_full,&
                                             rel_trans_full,&
                                             rel_zora_full
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rel_control_types'
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.

  ! Public data types

  PUBLIC :: rel_control_type

  ! Public subroutines

  PUBLIC :: rel_c_create,&
            rel_c_read_parameters,&
            rel_c_release,&
            rel_c_retain

! *****************************************************************************
!> \brief contains the parameters needed by a relativistic calculation
!> \param method used relativistic method (NONE, DKH)
!> \param DKH_order order of the DKH transformation (2,3)
!> \param transformation used blocks of the full matrix (FULL, MOLECULE, ATOM)
!> \param z_cutoff considered range of the Coulomb interaction
!> \param potential nuclear electron Coulomb potential (FULL, ERFC)
!> \par History
!>      09.2002 created [fawzi] for scf_control_type
!>      10.2008 copied to rel_control_type [JT]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE rel_control_type
    INTEGER                               :: rel_method
    INTEGER                               :: rel_DKH_order
    INTEGER                               :: rel_ZORA_type
    INTEGER                               :: rel_transformation
    INTEGER                               :: rel_z_cutoff
    INTEGER                               :: rel_potential
    INTEGER                               :: ref_count
  END TYPE rel_control_type

CONTAINS

! *****************************************************************************
!> \brief allocates and initializes an rel control object with the default values
!> \param rel_control the object to initialize
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2002 created [fawzi] for scf_control_type
!>      10.2008 copied to rel_control_type [JT]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE rel_c_create(rel_control,error)

    TYPE(rel_control_type), POINTER          :: rel_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rel_c_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    ALLOCATE (rel_control,STAT=stat)
    CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      ! Load the default values

      rel_control%rel_method = rel_none
      rel_control%rel_DKH_order = 2
      rel_control%rel_ZORA_type = rel_zora_full
      rel_control%rel_transformation = rel_trans_full
      rel_control%rel_z_cutoff = 1
      rel_control%rel_potential = rel_pot_full
      rel_control%ref_count = 1

    END IF

  END SUBROUTINE rel_c_create

! *****************************************************************************
!> \brief retains the given rel_control (see cp2k/doc/ReferenceCounting.html)
!> \param rel_control the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2002 created [fawzi] for scf_control_type
!>      10.2008 copied to rel_control_type [JT]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE rel_c_retain(rel_control,error)

    TYPE(rel_control_type), POINTER          :: rel_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rel_c_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(rel_control),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
      CPPrecondition(rel_control%ref_count > 0,cp_failure_level,routineP,error,failure)
      rel_control%ref_count = rel_control%ref_count + 1
    END IF

  END SUBROUTINE rel_c_retain

! *****************************************************************************
!> \brief releases the given rel_control (see cp2k/doc/ReferenceCounting.html)
!> \param rel_control the object to free
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      at the moment does nothing
!> \par History
!>      09.2002 created [fawzi] for scf_control_type
!>      10.2008 copied to rel_control_type [JT]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE rel_c_release(rel_control,error)

    TYPE(rel_control_type), POINTER          :: rel_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rel_c_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(rel_control)) THEN
      CPPrecondition(rel_control%ref_count>0,cp_failure_level,routineP,error,failure)
      rel_control%ref_count = rel_control%ref_count - 1
      IF (rel_control%ref_count < 1) THEN
        DEALLOCATE(rel_control,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
      END IF
    END IF

    NULLIFY (rel_control)

  END SUBROUTINE rel_c_release

! *****************************************************************************
!> \brief reads the parameters of the relativistic section into the given rel_control
!> \param rel_control the object that wil contain the values read
!> \param error controls log and error handling
!> \par History
!>      05.2001 created [Matthias] for scf_control_type
!>      09.2002 created separated scf_control type [fawzi]
!>      10.2008 copied to rel_control_type [JT]
!> \author Matthias Krack
! *****************************************************************************
  SUBROUTINE rel_c_read_parameters(rel_control,dft_section,error)

    TYPE(rel_control_type), POINTER          :: rel_control
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rel_c_read_parameters', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: rel_section

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(rel_control),cp_failure_level,routineP,error,failure)
    CPPrecondition((rel_control%ref_count > 0),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_section),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      rel_section => section_vals_get_subs_vals(dft_section,"RELATIVISTIC",&
                                                error=error)
      CALL section_vals_val_get(rel_section,"method",&
                                i_val=rel_control%rel_method,&
                                error=error)
      CALL section_vals_val_get(rel_section,"DKH_order",&
                                i_val=rel_control%rel_DKH_order,&
                                error=error)
      CALL section_vals_val_get(rel_section,"ZORA_TYPE",&
                                i_val=rel_control%rel_zora_type,&
                                error=error)
      CALL section_vals_val_get(rel_section,"transformation",&
                                i_val=rel_control%rel_transformation,&
                                error=error)
      CALL section_vals_val_get(rel_section,"z_cutoff",&
                                i_val=rel_control%rel_z_cutoff,&
                                error=error)
      CALL section_vals_val_get(rel_section,"potential",&
                                i_val=rel_control%rel_potential,&
                                error=error)

    END IF

  END SUBROUTINE rel_c_read_parameters

END MODULE rel_control_types
