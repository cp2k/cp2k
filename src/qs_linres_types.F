!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/qs_linres_types [1.0] *
!!
!!   NAME
!!     qs_linres_types
!!
!!   FUNCTION
!!     Type definitiona for linear response calculations
!!
!!   AUTHOR
!!     MI 
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!******************************************************************************
#include "cp_prep_globals.h"

MODULE qs_linres_types

  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_error_handling,               ONLY: cp_a_l,&
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
  USE kinds,                           ONLY: dp
  USE qs_loc_control,                  ONLY: localized_wfn_control_create ,&
                                             localized_wfn_control_release,&
                                             localized_wfn_control_retain ,&
                                             localized_wfn_control_type
  IMPLICIT NONE

  PRIVATE

!****s* qs_linres_types/linres_control_type
!!
!! NAME
!!   linres_control_type
!!
!! FUNCTION
!!   General settings for linear response calculations
!!
!! ATTRIBUTES
!!   - property : which quantity is to be calculated by LR
!!   - opt_method : method to optimize the psi1 by minimization of the second order term of the energy
!!   - preconditioner : which kind of preconditioner should be used, if any
!!   - localized_psi0 : don't use the canonical psi0, but the maximally localized wavefunctions
!!   - do_kernel  : the kernel is zero if the rho1 is zero as for the magnetic field perturbation 
!!   - tolerance : convergence criterium for the optimization of the psi1
!!
!! NOTES
!!   
!!
!!***
  TYPE linres_control_type
     INTEGER                                   :: ref_count
     INTEGER                                   :: property
     INTEGER                                   :: preconditioner
     INTEGER                                   :: opt_method
     LOGICAL                                   :: localized_psi0
     LOGICAL                                   :: do_kernel
     REAL(KIND=dp)                             :: tolerance
     TYPE(localized_wfn_control_type),POINTER  :: localized_wfn_control
  END TYPE linres_control_type

  TYPE nmr_env_type
     INTEGER                                   :: ref_count
     TYPE(cp_2d_r_p_type), DIMENSION(:),POINTER  :: centers_set

  END TYPE nmr_env_type

  CHARACTER(LEN=*), PARAMETER :: module_name ="qs_linres_types"

! *** Public data types ***

  PUBLIC ::  linres_control_type, nmr_env_type


! *** Public subroutines ***

  PUBLIC :: linres_control_create, linres_control_retain, linres_control_release

!!***
! *****************************************************************************

CONTAINS  

  SUBROUTINE linres_control_create(linres_control,error)
 
     TYPE(linres_control_type), POINTER        ::  linres_control
     TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error
     CHARACTER(len=*), PARAMETER :: routineN = 'linres_control_create', &
      routineP = module_name//':'//routineN
     INTEGER                                   :: istat
     LOGICAL                                   :: failure


     failure =.FALSE.

     CPPrecondition(.NOT.ASSOCIATED(linres_control),cp_failure_level,routineP,error,failure)
     IF (.NOT. failure) THEN
       ALLOCATE (linres_control,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       linres_control%ref_count=1
       CALL localized_wfn_control_create(linres_control%localized_wfn_control, error=error)

     END IF

  END SUBROUTINE linres_control_create

  SUBROUTINE linres_control_release(linres_control,error)
 
     TYPE(linres_control_type), POINTER        ::  linres_control
     TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error
     CHARACTER(len=*), PARAMETER :: routineN = 'linres_control_release', &
      routineP = module_name//':'//routineN
     INTEGER                                   :: istat
     LOGICAL                                   :: failure

     failure =.FALSE.

     CPPrecondition(ASSOCIATED(linres_control),cp_failure_level,routineP,error,failure)
     IF (.NOT. failure) THEN
       CPPostcondition(linres_control%ref_count>0,cp_failure_level,routineP,error,failure)
       linres_control%ref_count=linres_control%ref_count-1
       IF(linres_control%ref_count<1)THEN
         CALL localized_wfn_control_release(linres_control%localized_wfn_control, error=error)
         DEALLOCATE(linres_control,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF
     END IF
     NULLIFY(linres_control)
  END SUBROUTINE linres_control_release

  SUBROUTINE linres_control_retain(linres_control,error)
 
     TYPE(linres_control_type), POINTER        ::  linres_control
     TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error
     CHARACTER(len=*), PARAMETER :: routineN = 'linres_control_retain', &
      routineP = module_name//':'//routineN
     INTEGER                                   :: istat
     LOGICAL                                   :: failure

     failure =.FALSE.

     CPPrecondition(ASSOCIATED(linres_control),cp_failure_level,routineP,error,failure)
     IF (.NOT. failure) THEN
       CPPostcondition(linres_control%ref_count>0,cp_failure_level,routineP,error,failure)
       linres_control%ref_count=linres_control%ref_count+1
     END IF
     
  END SUBROUTINE linres_control_retain
       
END MODULE qs_linres_types

