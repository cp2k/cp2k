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
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type
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
  USE qs_rho_atom_types,               ONLY: rho_atom_p_type
  USE qs_rho_types,                    ONLY: qs_rho_p_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix_set,&
                                             real_matrix_p_type
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
!! AUTHOR
!!   MI 
!!
!!***
  TYPE linres_control_type
     INTEGER                                   :: ref_count
     INTEGER                                   :: property
     INTEGER                                   :: preconditioner
     INTEGER                                   :: opt_method
     INTEGER                                   :: max_scf
     LOGICAL                                   :: localized_psi0
     LOGICAL                                   :: do_kernel
     REAL(KIND=dp)                             :: eps_scf
     TYPE(localized_wfn_control_type),POINTER  :: localized_wfn_control
  END TYPE linres_control_type

!****s* qs_linres_types/nmr_env_type
!!
!! NAME
!!   nmr_env_type
!!
!! FUNCTION
!!
!! ATTRIBUTES
!!   - ref_count
!!   - centers_set : poiter to the  centers of the maximally localized psi0
!!   - spreads_set : pointer to the spreads of the maximally localized psi0 (to evaluate thei extension)
!!   - op_p_ao     : sparse matrices for the px, py and pz operator appliet to the atomic orbitals
!!   - op_rxp_ao   : the rxp operator should be psi0 dependent, therefore it is used a s a temporary matrix
!!   - jp1_ao      : current density matrices. If the current density has to be calculated,
!!                   we need to store at least two set of current density matrices, 
!!                   each set has the x y and z components.
!!   - jrho1_set   : current density on the global grid, if gapw this is only the soft part
!!   - jrho1_atom_set : current density on the local atomic grids (only if gapw)
!!
!! NOTES
!!   
!! AUTHOR
!!   MI 
!!***

  TYPE nmr_env_type
     INTEGER                                     :: ref_count
     LOGICAL                                     :: current_density, &
                                                    full_nmr
     REAL(dp),DIMENSION(3,3)                     :: chemical_shift
     TYPE(cp_2d_r_p_type), DIMENSION(:),POINTER  :: centers_set
     TYPE(cp_fm_pool_p_type), DIMENSION(:), &
       POINTER                                   :: ao_mo_fm_pools
     TYPE(real_matrix_p_type), DIMENSION(:), &
          POINTER                                :: op_p_ao, op_rmd_ao
     TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER :: psi1_p, psi1_rxp, psi1_D
     TYPE(real_matrix_p_type), DIMENSION(:,:), &
          POINTER                                :: jp1_ao
     TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER :: p_psi0, rxp_psi0
     TYPE(qs_rho_p_type),DIMENSION(:),   POINTER :: jrho1_set
     TYPE(rho_atom_p_type), DIMENSION(:), POINTER:: jrho1_atom_set
  END TYPE nmr_env_type

  CHARACTER(LEN=*), PARAMETER :: module_name ="qs_linres_types"

! *** Public data types ***

  PUBLIC ::  linres_control_type, nmr_env_type


! *** Public subroutines ***

  PUBLIC :: linres_control_create, linres_control_retain, linres_control_release,&
            nmr_env_create

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

  SUBROUTINE nmr_env_create(nmr_env,error)

    TYPE(nmr_env_type)                         :: nmr_env
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error
     CHARACTER(len=*), PARAMETER :: routineN = 'nmr_env_create', &
      routineP = module_name//':'//routineN
     INTEGER                                   :: istat
     LOGICAL                                   :: failure

     failure =.FALSE.

     CPPrecondition(nmr_env%ref_count==0, cp_failure_level,routineP,error,failure)
     IF(.NOT. failure) THEN
       nmr_env%ref_count = 1
       NULLIFY(nmr_env%centers_set)
       NULLIFY(nmr_env%ao_mo_fm_pools)
       NULLIFY(nmr_env%op_p_ao)
       NULLIFY(nmr_env%op_rmd_ao)
       NULLIFY(nmr_env%psi1_p)
       NULLIFY(nmr_env%psi1_rxp)
       NULLIFY(nmr_env%psi1_D)
       NULLIFY(nmr_env%jp1_ao)
       NULLIFY(nmr_env%p_psi0)
       NULLIFY(nmr_env%rxp_psi0)
       NULLIFY(nmr_env%jrho1_set)
       NULLIFY(nmr_env%jrho1_atom_set)
     END IF 


  END SUBROUTINE nmr_env_create
       
END MODULE qs_linres_types

