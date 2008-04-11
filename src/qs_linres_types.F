!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Type definitiona for linear response calculations
!> \author MI 
! *****************************************************************************
MODULE qs_linres_types

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type
  USE kinds,                           ONLY: dp
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_loc_control,                  ONLY: localized_wfn_control_release,&
                                             localized_wfn_control_type
  USE qs_rho_atom_types,               ONLY: rho_atom_coeff,&
                                             rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_p_type,&
                                             qs_rho_release
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

!****s* qs_linres_types/linres_control_type

! *****************************************************************************
!> \brief General settings for linear response calculations
!> \param property which quantity is to be calculated by LR
!> \param opt_method method to optimize the psi1 by minimization of the second order term of the energy
!> \param preconditioner which kind of preconditioner should be used, if any
!> \param localized_psi 0 : don't use the canonical psi0, but the maximally localized wavefunctions
!> \param do_kernel the kernel is zero if the rho1 is zero as for the magnetic field perturbation 
!> \param tolerance convergence criterium for the optimization of the psi1
!> \author MI 
! *****************************************************************************
  TYPE linres_control_type
     INTEGER                                   :: ref_count
     INTEGER                                   :: property
     INTEGER                                   :: preconditioner_type
     REAL(dp)                                  :: energy_gap
     INTEGER                                   :: opt_method
     INTEGER                                   :: ls_method
     REAL(dp)                                  :: ds_min
     INTEGER                                   :: max_scf
     LOGICAL                                   :: localized_psi0
     LOGICAL                                   :: do_kernel
     LOGICAL                                   :: converged
     LOGICAL                                   :: check_ortho
     REAL(KIND=dp)                             :: eps_scf
     TYPE(localized_wfn_control_type),POINTER  :: localized_wfn_control
     CHARACTER(LEN=8)                          :: flag
  END TYPE linres_control_type

!****s* qs_linres_types/current_env_type

! *****************************************************************************
!> \param ref_coun t
!> \param full_nmr true if the full correction is calculated
!> \param simplenmr_done , fullnmr_done : flags that indicate what has been
!>                    already calculated: used for restart
!> \param centers_set centers of the maximally localized psi0
!> \param spreads_set spreads of the maximally localized psi0 
!> \param op_p_ao sparse matrixes for the px, py and pz operator 
!>                    in the contracted basis set representation 
!> \param op_rmd_ao sparse matrixes used to construct the (r-d)xp operator
!>                    in the contracted basis set representation
!> \param p_psi 0      : full matrixes, operator p applied to psi0
!> \param rxp_psi 0    : full matrixes, operator (r-d)xp applied to psi0
!> \param psi 1_p      : response wavefunctions to the perturbation given by 
!>                    H1=p (xyz)  applied to psi0
!> \param psi 1_rxp    : response wavefunctions to the perturbation given by 
!>                    H1=(r-d_i)xp applied to psi0_i where d_i is the center
!> \param psi 1_D      : response wavefunctions to the perturbation given by 
!>                    H1=(d_j-d_i)xp applied to psi0_i where d_i is the center
!>                    and d_j is the center of psi0_j and psi1_D_j is the result
!>                    This operator has to be used in nstate scf calculations,
!>                    one for each psi1_D_j vector
!> \param chemical_shift the tensor for each atom
!> \param chi_tensor the susceptibility tensor
!> \param jp 1_ao :      current density matrices. 
!> \param jp 2_ao :      density matrices used to compute the r-dependent terms
!> \param jrho 1_set   : current density on the global grid, if gapw this is only the soft part
!> \param jrho 1_atom_set : current density on the local atomic grids (only if gapw)
!> \author MI 
! *****************************************************************************
  TYPE current_env_type
     LOGICAL                                     :: full,simple_done(6),simple_converged(6)
     LOGICAL, DIMENSION(:,:), POINTER            :: full_done
     INTEGER                                     :: ref_count,nao,nstates(2),gauge
     INTEGER, DIMENSION(:    ), POINTER          :: list_cubes
     INTEGER, DIMENSION(:,:,:), POINTER          :: statetrueindex
     CHARACTER(LEN=30)                           :: gauge_name
     REAL(dp), DIMENSION(:,:), POINTER           :: basisfun_center
     TYPE(cp_2d_r_p_type),    DIMENSION(:  ), POINTER :: centers_set
     TYPE(cp_fm_pool_p_type), DIMENSION(:  ), POINTER :: ao_mo_fm_pools
     TYPE(cp_fm_p_type),      DIMENSION(:,:), POINTER :: psi1_p,psi1_rxp,psi1_D,p_psi0,rxp_psi0
     TYPE(jrho_atom_type),    DIMENSION(:  ), POINTER :: jrho1_atom_set
     TYPE(qs_rho_p_type),     DIMENSION(:,:), POINTER :: jrho1_set
     TYPE(real_matrix_p_type),DIMENSION(:  ), POINTER :: op_p_ao,op_rmd_ao,jp0_ao,jp1_ao
     TYPE(real_matrix_p_type),DIMENSION(:,:), POINTER :: jp2_ao
  END TYPE current_env_type

  TYPE nmr_env_type
     INTEGER                             :: ref_count, n_nics
     INTEGER, DIMENSION(:), POINTER      :: cs_atom_list
     INTEGER, DIMENSION(:), POINTER      :: do_calc_cs_atom
     LOGICAL                             :: do_nics
     REAL(dp)                            :: shift_factor, shift_factor_gapw, chi_factor, &
          &                                 chi_SI2shiftppm, chi_SI2ppmcgs, &
          &                                 chi_tensor(3,3), chi_tensor_loc(3,3)
     REAL(dp), DIMENSION(:,:  ), POINTER :: r_nics
     REAL(dp), DIMENSION(:,:,:), POINTER :: chemical_shift, chemical_shift_loc, &
          &                                 chemical_shift_nics_loc, chemical_shift_nics
  END TYPE nmr_env_type

! *****************************************************************************
  TYPE epr_env_type
     !INTEGER                                     :: gauge
     !CHARACTER(LEN=30)                           :: gauge_name
     INTEGER                                     :: ref_count!, nao, nstates(2)
     !INTEGER, DIMENSION(:), POINTER              :: list_cubes
     !INTEGER, DIMENSION(:,:,:), POINTER          :: statetrueindex
     !LOGICAL                                     :: full_epr,simpleepr_done(6),&
     !                                               simpleepr_converged(6),&
     !                                               restart_epr
     !LOGICAL, DIMENSION(:,:), POINTER            :: fullepr_done
     !REAL(dp), DIMENSION(:,:), POINTER           :: basisfun_center
     REAL(dp)                                    :: g_free_factor, g_soo_chicorr_factor, g_soo_factor,&
          &                                         g_so_factor, g_so_factor_gapw, g_zke_factor, g_zke
     REAL(dp), DIMENSION(:,:), POINTER           :: g_total, g_so, g_soo, bind0
     !TYPE(cp_2d_r_p_type), DIMENSION(:),POINTER  :: centers_set
     !TYPE(cp_fm_pool_p_type), DIMENSION(:), &
     !  POINTER                                   :: ao_mo_fm_pools
     !TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER :: psi1_p, psi1_rxp, psi1_D
     !TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER :: p_psi0, rxp_psi0
     !TYPE(jrho_atom_type),DIMENSION(:),&
     !  POINTER                                   :: jrho1_atom_set
     !TYPE(qs_rho_p_type),DIMENSION(:,:), POINTER :: jrho1_set
     !TYPE(real_matrix_p_type), DIMENSION(:), &
     !     POINTER                                :: op_p_ao, op_rmd_ao
     !TYPE(real_matrix_p_type), DIMENSION(:,:), &
     !     POINTER                                :: jp2_ao
     !TYPE(real_matrix_p_type), DIMENSION(:), &
     !     POINTER                                :: jp0_ao, jp1_ao
     TYPE(qs_rho_p_type),      DIMENSION(:,:), POINTER :: nablavks_set
     TYPE(nablavks_atom_type), DIMENSION(:  ), POINTER :: nablavks_atom_set
     TYPE(qs_rho_p_type),      DIMENSION(:,:), POINTER :: bind_set
     TYPE(rho_atom_coeff),     DIMENSION(:,:), POINTER :: bind_atom_set
     TYPE(rho_atom_type),      DIMENSION(:  ), POINTER :: vks_atom_set
  END TYPE epr_env_type

! *****************************************************************************
  TYPE nablavks_atom_type
     TYPE(rho_atom_coeff), DIMENSION(:,:),&
          POINTER                                :: nablavks_vec_rad_h,&
                                                    nablavks_vec_rad_s
  END TYPE nablavks_atom_type

! *****************************************************************************
  TYPE jrho_atom_p_type
     TYPE(jrho_atom_type), POINTER       :: jrho_atom
  END TYPE jrho_atom_p_type 

! *****************************************************************************
  TYPE jrho_atom_type
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER      :: cjc_h, cjc_s, cjc0_h, cjc0_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER      :: cjc_ii_h, cjc_ii_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER      :: cjc_iii_h, cjc_iii_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER      :: jrho_rad_h,&
                                               jrho_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:,:,:),&
                               POINTER      :: jrho_vec_rad_h,&
                                               jrho_vec_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_h, jrho_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_a_h, jrho_a_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_b_h, jrho_b_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_a_h_ii, jrho_a_s_ii
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_b_h_ii, jrho_b_s_ii
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_a_h_iii, jrho_a_s_iii
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                               :: jrho_b_h_iii, jrho_b_s_iii
  END TYPE jrho_atom_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_types'

! *** Public data types ***

  PUBLIC ::  allocate_jrho_atom_rad, deallocate_jrho_atom_set, get_nmr_env, &
             get_current_env, &
             init_jrho_atom_set, linres_control_type, init_nablavks_atom_set, &
             nmr_env_type, jrho_atom_type, set_nmr_env, set2zero_jrho_atom_rad, &
             get_epr_env, epr_env_type, set_epr_env, nablavks_atom_type, &
             deallocate_nablavks_atom_set, current_env_type, set_current_env

! *** Public subroutines ***

  PUBLIC :: linres_control_create, linres_control_retain, linres_control_release,&
            nmr_env_create, epr_env_create, current_env_create

CONTAINS  

! *****************************************************************************
  SUBROUTINE linres_control_create(linres_control,error)
 
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'linres_control_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

     failure =.FALSE.

     CPPrecondition(.NOT.ASSOCIATED(linres_control),cp_failure_level,routineP,error,failure)
     IF (.NOT. failure) THEN
       ALLOCATE (linres_control,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       linres_control%ref_count=1
       NULLIFY(linres_control%localized_wfn_control)
!       CALL localized_wfn_control_create(linres_control%localized_wfn_control, error=error)

     END IF

  END SUBROUTINE linres_control_create

! *****************************************************************************
  SUBROUTINE linres_control_release(linres_control,error)
 
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'linres_control_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

     failure =.FALSE.
     
     IF (ASSOCIATED(linres_control)) THEN
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

! *****************************************************************************
  SUBROUTINE linres_control_retain(linres_control,error)
 
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'linres_control_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     failure =.FALSE.

     CPPrecondition(ASSOCIATED(linres_control),cp_failure_level,routineP,error,failure)
     IF (.NOT. failure) THEN
       CPPostcondition(linres_control%ref_count>0,cp_failure_level,routineP,error,failure)
       linres_control%ref_count=linres_control%ref_count+1
     END IF
     
  END SUBROUTINE linres_control_retain

! *****************************************************************************
  SUBROUTINE current_env_create(current_env,error)

    TYPE(current_env_type)                   :: current_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'current_env_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPPrecondition(current_env%ref_count==0, cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN
       current_env%ref_count = 1
       NULLIFY(current_env%full_done)
       NULLIFY(current_env%list_cubes)
       NULLIFY(current_env%statetrueindex)
       NULLIFY(current_env%basisfun_center)
       NULLIFY(current_env%centers_set)
       NULLIFY(current_env%ao_mo_fm_pools)
       NULLIFY(current_env%psi1_p)
       NULLIFY(current_env%psi1_rxp)
       NULLIFY(current_env%psi1_D)
       NULLIFY(current_env%p_psi0)
       NULLIFY(current_env%rxp_psi0)
       NULLIFY(current_env%jrho1_atom_set)
       NULLIFY(current_env%jrho1_set)
       NULLIFY(current_env%op_p_ao)
       NULLIFY(current_env%op_rmd_ao)
       NULLIFY(current_env%jp2_ao)
       NULLIFY(current_env%jp0_ao)
       NULLIFY(current_env%jp1_ao)
    END IF

  END SUBROUTINE current_env_create
! *****************************************************************************
  SUBROUTINE nmr_env_create(nmr_env,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nmr_env_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     failure =.FALSE.

     CPPrecondition(nmr_env%ref_count==0, cp_failure_level,routineP,error,failure)
     IF(.NOT. failure) THEN
       nmr_env%ref_count = 1
       NULLIFY(nmr_env%chemical_shift)
       NULLIFY(nmr_env%chemical_shift_loc)
       NULLIFY(nmr_env%chemical_shift_nics_loc)
       NULLIFY(nmr_env%chemical_shift_nics)
       NULLIFY(nmr_env%r_nics)
       NULLIFY(nmr_env%cs_atom_list)
       NULLIFY(nmr_env%do_calc_cs_atom)
     END IF 

  END SUBROUTINE nmr_env_create
       
! *****************************************************************************
  SUBROUTINE epr_env_create(epr_env,error)
    
    TYPE(epr_env_type)                       :: epr_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'epr_env_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPPrecondition(epr_env%ref_count==0, cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN
       epr_env%ref_count = 1
       !NULLIFY(epr_env%centers_set)
       !NULLIFY(epr_env%ao_mo_fm_pools)
       !NULLIFY(epr_env%op_p_ao)
       !NULLIFY(epr_env%op_rmd_ao)
       !NULLIFY(epr_env%psi1_p)
       !NULLIFY(epr_env%psi1_rxp)
       !NULLIFY(epr_env%psi1_D)
       !NULLIFY(epr_env%jp0_ao)
       !NULLIFY(epr_env%jp1_ao)
       !NULLIFY(epr_env%jp2_ao)
       !NULLIFY(epr_env%p_psi0)
       !NULLIFY(epr_env%rxp_psi0)
       !NULLIFY(epr_env%jrho1_set)
       !NULLIFY(epr_env%jrho1_atom_set)
       !NULLIFY(epr_env%fullepr_done)
       !NULLIFY(epr_env%basisfun_center)
       !NULLIFY(epr_env%statetrueindex)
       !NULLIFY(epr_env%list_cubes)
       NULLIFY(epr_env%nablavks_set)
       NULLIFY(epr_env%nablavks_atom_set)
       NULLIFY(epr_env%bind_set)
       NULLIFY(epr_env%bind_atom_set)
       NULLIFY(epr_env%g_total)
       NULLIFY(epr_env%g_so)
       NULLIFY(epr_env%g_soo)
       NULLIFY(epr_env%bind0)
       NULLIFY(epr_env%vks_atom_set)
     END IF 

  END SUBROUTINE epr_env_create
       
  SUBROUTINE get_current_env(current_env,simple_done,simple_converged,full_done,ref_count,nao,&
       &                 nstates,gauge,list_cubes,statetrueindex,gauge_name,basisfun_center,&
       &                 centers_set,ao_mo_fm_pools,psi1_p,psi1_rxp,psi1_D,p_psi0,rxp_psi0,&
       &                 jrho1_atom_set,jrho1_set,op_p_ao,op_rmd_ao,jp2_ao,jp0_ao,jp1_ao,&
       &                 error)
    
    TYPE(current_env_type), OPTIONAL         :: current_env
    LOGICAL, OPTIONAL                        :: simple_done(6), &
                                                simple_converged(6)
    LOGICAL, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: full_done
    INTEGER, OPTIONAL                        :: ref_count, nao, nstates(2), &
                                                gauge
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: list_cubes
    INTEGER, DIMENSION(:, :, :), OPTIONAL, &
      POINTER                                :: statetrueindex
    CHARACTER(LEN=30), OPTIONAL              :: gauge_name
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: basisfun_center
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: centers_set
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ao_mo_fm_pools
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: psi1_p, psi1_rxp, psi1_D, &
                                                p_psi0, rxp_psi0
    TYPE(jrho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: jrho1_atom_set
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: jrho1_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: op_p_ao, op_rmd_ao
    TYPE(real_matrix_p_type), &
      DIMENSION(:, :), OPTIONAL, POINTER     :: jp2_ao
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: jp0_ao, jp1_ao
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_current_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!

    failure =.FALSE.

    CPPrecondition(current_env%ref_count>0, cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN
       IF(PRESENT(simple_done     )) simple_done(1:6)      = current_env%simple_done(1:6)
       IF(PRESENT(simple_converged)) simple_converged(1:6) = current_env%simple_converged(1:6)
       IF(PRESENT(full_done       )) full_done             => current_env%full_done
       IF(PRESENT(ref_count       )) ref_count             =  current_env%ref_count
       IF(PRESENT(nao             )) nao                   =  current_env%nao
       IF(PRESENT(nstates         )) nstates(1:2)          =  current_env%nstates(1:2)
       IF(PRESENT(gauge           )) gauge                 =  current_env%gauge
       IF(PRESENT(list_cubes      )) list_cubes            => current_env%list_cubes
       IF(PRESENT(statetrueindex  )) statetrueindex        => current_env%statetrueindex
       IF(PRESENT(gauge_name      )) gauge_name            =  current_env%gauge_name
       IF(PRESENT(basisfun_center )) basisfun_center       => current_env%basisfun_center
       IF(PRESENT(centers_set     )) centers_set           => current_env%centers_set
       IF(PRESENT(ao_mo_fm_pools  )) ao_mo_fm_pools        => current_env%ao_mo_fm_pools
       IF(PRESENT(psi1_p          )) psi1_p                => current_env%psi1_p
       IF(PRESENT(psi1_rxp        )) psi1_rxp              => current_env%psi1_rxp
       IF(PRESENT(psi1_D          )) psi1_D                => current_env%psi1_D
       IF(PRESENT(p_psi0          )) p_psi0                => current_env%p_psi0
       IF(PRESENT(rxp_psi0        )) rxp_psi0              => current_env%rxp_psi0
       IF(PRESENT(jrho1_atom_set  )) jrho1_atom_set        => current_env%jrho1_atom_set
       IF(PRESENT(jrho1_set       )) jrho1_set             => current_env%jrho1_set
       IF(PRESENT(op_rmd_ao       )) op_rmd_ao             => current_env%op_rmd_ao
       IF(PRESENT(op_p_ao         )) op_p_ao               => current_env%op_p_ao
       IF(PRESENT(jp2_ao          )) jp2_ao                => current_env%jp2_ao
       IF(PRESENT(jp1_ao          )) jp1_ao                => current_env%jp1_ao
       IF(PRESENT(jp0_ao          )) jp0_ao                => current_env%jp0_ao
    ENDIF

  END SUBROUTINE get_current_env

! *****************************************************************************
  SUBROUTINE get_nmr_env(nmr_env, n_nics, chi_tensor, cs_atom_list, do_calc_cs_atom, &
             chi_tensor_loc, r_nics, chemical_shift,&
             chemical_shift_loc, chemical_shift_nics_loc, chemical_shift_nics,&
             error)

    TYPE(nmr_env_type)                       :: nmr_env
    INTEGER, INTENT(OUT), OPTIONAL           :: n_nics
    REAL(dp), INTENT(out), OPTIONAL          :: chi_tensor(3,3)
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: cs_atom_list, do_calc_cs_atom
    REAL(dp), INTENT(out), OPTIONAL          :: chi_tensor_loc(3,3)
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: r_nics
    REAL(dp), DIMENSION(:, :, :), OPTIONAL, &
      POINTER                                :: chemical_shift, &
                                                chemical_shift_loc, &
                                                chemical_shift_nics_loc, &
                                                chemical_shift_nics
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_nmr_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPPrecondition(nmr_env%ref_count>0, cp_failure_level,routineP,error,failure)

    IF(PRESENT(n_nics             )) n_nics              =  nmr_env%n_nics
    IF(PRESENT(cs_atom_list       )) cs_atom_list        => nmr_env%cs_atom_list
    IF(PRESENT(do_calc_cs_atom    )) do_calc_cs_atom     => nmr_env%do_calc_cs_atom
    IF(PRESENT(chi_tensor         )) chi_tensor          =  nmr_env%chi_tensor
    IF(PRESENT(chi_tensor_loc     )) chi_tensor_loc      =  nmr_env%chi_tensor_loc
    IF(PRESENT(chemical_shift     )) chemical_shift      => nmr_env%chemical_shift
    IF(PRESENT(chemical_shift_loc )) chemical_shift_loc  => nmr_env%chemical_shift_loc
    IF(PRESENT(chemical_shift_nics)) chemical_shift_nics => nmr_env%chemical_shift_nics
    IF(PRESENT(r_nics             )) r_nics              => nmr_env%r_nics
    IF(PRESENT(chemical_shift_nics_loc)) chemical_shift_nics_loc => nmr_env%chemical_shift_nics_loc

  END SUBROUTINE get_nmr_env

! *****************************************************************************
  SUBROUTINE set_current_env(current_env,jrho1_atom_set,jrho1_set,error)

    TYPE(current_env_type)                   :: current_env
    TYPE(jrho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: jrho1_atom_set
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: jrho1_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_current_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iB, idir
    LOGICAL                                  :: failure

    failure =.FALSE.

    CPPrecondition(current_env%ref_count>0, cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN

       IF(PRESENT(jrho1_atom_set)) THEN
          IF(ASSOCIATED(current_env%jrho1_atom_set)) THEN
             CALL deallocate_jrho_atom_set(current_env%jrho1_atom_set,error=error)
          ENDIF
          current_env%jrho1_atom_set => jrho1_atom_set
       END IF
       
       IF(PRESENT(jrho1_set)) THEN
          IF(ASSOCIATED(current_env%jrho1_set)) THEN
             DO iB = 1,3
                DO idir = 1,3
                   CALL qs_rho_release(current_env%jrho1_set(idir,iB)%rho,error=error)
                END DO
             END DO
          END IF
          current_env%jrho1_set => jrho1_set
       END IF
    END IF

  END SUBROUTINE set_current_env
! *****************************************************************************
  SUBROUTINE set_nmr_env(nmr_env,shift_factor,chi_factor,chi_SI2shiftppm,chi_SI2ppmcgs,&
             chi_tensor,chi_tensor_loc,error)

    TYPE(nmr_env_type)                       :: nmr_env
    REAL(dp), INTENT(IN), OPTIONAL :: shift_factor, chi_factor, &
      chi_SI2shiftppm, chi_SI2ppmcgs, chi_tensor(3,3), chi_tensor_loc(3,3)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_nmr_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iB, idir
    LOGICAL                                  :: failure

    failure =.FALSE.

    CPPrecondition(nmr_env%ref_count>0, cp_failure_level,routineP,error,failure)

    IF(PRESENT(shift_factor   )) nmr_env%chi_factor      = chi_factor
    IF(PRESENT(shift_factor   )) nmr_env%chi_factor      = chi_factor
    IF(PRESENT(chi_SI2shiftppm)) nmr_env%chi_SI2shiftppm = chi_SI2shiftppm
    IF(PRESENT(chi_SI2ppmcgs  )) nmr_env%chi_SI2ppmcgs   = chi_SI2ppmcgs
    IF(PRESENT(chi_tensor)) THEN
       DO iB=1,3
          DO idir=1,3
             nmr_env%chi_tensor(idir,iB) = chi_tensor(idir,iB)
          END DO
       END DO
     END IF
     IF(PRESENT(chi_tensor_loc))   THEN
        DO iB=1,3
           DO idir=1,3
              nmr_env%chi_tensor_loc(idir,iB) = chi_tensor_loc(idir,iB)
           END DO
        END DO
     END IF

  END SUBROUTINE set_nmr_env

! *****************************************************************************
  SUBROUTINE get_epr_env(epr_env, &
             g_total, g_so, g_soo, bind0, &
             nablavks_set, nablavks_atom_set, &
             bind_set, bind_atom_set, &
             vks_atom_set, error)

    TYPE(epr_env_type)                       :: epr_env
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: g_total, g_so, g_soo, bind0
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: nablavks_set
    TYPE(nablavks_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: nablavks_atom_set
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: bind_set
    TYPE(rho_atom_coeff), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: bind_atom_set
    TYPE(rho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: vks_atom_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_epr_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!INTEGER, INTENT(OUT), OPTIONAL           :: nao, nstates(2)
!INTEGER, DIMENSION(:), OPTIONAL, POINTER :: list_cubes
!INTEGER, DIMENSION(:, :, :), OPTIONAL, &
!  POINTER                                :: statetrueindex
!, &
!basisfun_center
!TYPE(cp_2d_r_p_type), DIMENSION(:), &
!  OPTIONAL, POINTER                      :: centers_set
!TYPE(cp_fm_p_type), DIMENSION(:, :), &
!  OPTIONAL, POINTER                      :: psi1_p, psi1_rxp, psi1_D, &
!                                            p_psi0, rxp_psi0
!TYPE(jrho_atom_type), DIMENSION(:), &
!  OPTIONAL, POINTER                      :: jrho1_atom_set
!TYPE(qs_rho_p_type), DIMENSION(:, :), &
!  OPTIONAL, POINTER                      :: jrho1_set
!TYPE(real_matrix_p_type), DIMENSION(:), &
!  OPTIONAL, POINTER                      :: op_p_ao, op_rmd_ao, jp0_ao, &
!                                            jp1_ao
!TYPE(real_matrix_p_type), &
!  DIMENSION(:, :), OPTIONAL, POINTER     :: jp2_ao

     failure =.FALSE.

     CPPrecondition(epr_env%ref_count>0, cp_failure_level,routineP,error,failure)

     !IF(PRESENT(nao)) nao=epr_env%nao
     !IF(PRESENT(nstates)) nstates(1:2)=epr_env%nstates(1:2)
     !IF(PRESENT(list_cubes)) list_cubes => epr_env%list_cubes
     !IF(PRESENT(statetrueindex)) statetrueindex => epr_env%statetrueindex
     !IF(PRESENT(centers_set)) centers_set => epr_env%centers_set
     !IF(PRESENT(basisfun_center)) basisfun_center => epr_env%basisfun_center
     !IF(PRESENT(centers_set)) centers_set => epr_env%centers_set
     !IF(PRESENT(psi1_p)) psi1_p => epr_env%psi1_p
     !IF(PRESENT(psi1_rxp)) psi1_rxp => epr_env%psi1_rxp
     !IF(PRESENT(psi1_D)) psi1_D => epr_env%psi1_D
     !IF(PRESENT(p_psi0)) p_psi0 => epr_env%p_psi0
     !IF(PRESENT(rxp_psi0)) rxp_psi0 => epr_env%rxp_psi0
     !IF(PRESENT(jrho1_atom_set)) jrho1_atom_set => epr_env%jrho1_atom_set
     !IF(PRESENT(jrho1_set)) jrho1_set => epr_env%jrho1_set
     !IF(PRESENT(op_rmd_ao)) op_rmd_ao => epr_env%op_rmd_ao
     !IF(PRESENT(op_p_ao)) op_p_ao => epr_env%op_p_ao
     !IF(PRESENT(jp2_ao)) jp2_ao => epr_env%jp2_ao
     !IF(PRESENT(jp1_ao)) jp1_ao => epr_env%jp1_ao
     !IF(PRESENT(jp0_ao)) jp0_ao => epr_env%jp0_ao
     IF(PRESENT(g_total)) g_total => epr_env%g_total
     IF(PRESENT(g_so)) g_so => epr_env%g_so
     IF(PRESENT(g_soo)) g_soo => epr_env%g_soo
     IF(PRESENT(bind0)) bind0 => epr_env%bind0
     IF(PRESENT(nablavks_set)) nablavks_set => epr_env%nablavks_set
     IF(PRESENT(nablavks_atom_set)) nablavks_atom_set => epr_env%nablavks_atom_set
     IF(PRESENT(bind_set)) bind_set => epr_env%bind_set
     IF(PRESENT(bind_atom_set)) bind_atom_set => epr_env%bind_atom_set

  END SUBROUTINE get_epr_env

! *****************************************************************************
  SUBROUTINE set_epr_env(epr_env,g_free_factor,g_soo_chicorr_factor,&
                         g_soo_factor,g_so_factor,g_so_factor_gapw,&
                         g_zke_factor,&!jrho1_atom_set,jrho1_set,&
                         nablavks_set,nablavks_atom_set,&
                         error)

    TYPE(epr_env_type)                       :: epr_env
    REAL(dp), INTENT(IN), OPTIONAL :: g_free_factor, g_soo_chicorr_factor, &
      g_soo_factor, g_so_factor, g_so_factor_gapw, g_zke_factor
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: nablavks_set
    TYPE(nablavks_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: nablavks_atom_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_epr_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idir, ispin
    LOGICAL                                  :: failure

!TYPE(jrho_atom_type), DIMENSION(:), &
!  OPTIONAL, POINTER                      :: jrho1_atom_set
!jrho1_set,

    failure =.FALSE.

    CPPrecondition(epr_env%ref_count>0, cp_failure_level,routineP,error,failure)

    IF(PRESENT(g_free_factor)) epr_env%g_free_factor=g_free_factor
    IF(PRESENT(g_zke_factor)) epr_env%g_zke_factor=g_zke_factor
    IF(PRESENT(g_so_factor)) epr_env%g_so_factor=g_so_factor
    IF(PRESENT(g_so_factor_gapw)) epr_env%g_so_factor_gapw=g_so_factor_gapw
    IF(PRESENT(g_soo_factor)) epr_env%g_soo_factor=g_soo_factor
    IF(PRESENT(g_soo_chicorr_factor)) epr_env%g_soo_chicorr_factor=g_soo_chicorr_factor

    !IF(PRESENT(jrho1_atom_set)) THEN
    !   IF(ASSOCIATED(epr_env%jrho1_atom_set)) THEN
    !      CALL deallocate_jrho_atom_set(epr_env%jrho1_atom_set,error=error)
    !   ENDIF
    !   epr_env%jrho1_atom_set => jrho1_atom_set
    !ENDIF

    !IF(PRESENT(jrho1_set)) THEN
    !   IF(ASSOCIATED(epr_env%jrho1_set)) THEN
    !       DO iB = 1,3
    !         DO idir = 1,3
    !           CALL qs_rho_release(epr_env%jrho1_set(idir,iB)%rho,error=error)
    !         END DO
    !       END DO
    !    END IF
    !    epr_env%jrho1_set => jrho1_set
    !ENDIF

    IF(PRESENT(nablavks_set)) THEN
        IF(ASSOCIATED(epr_env%nablavks_set)) THEN
           DO ispin = 1,2
             DO idir = 1,3
               CALL qs_rho_release(epr_env%nablavks_set(idir,ispin)%rho,error=error)
             END DO
           END DO
        END IF
        epr_env%nablavks_set => nablavks_set
    ENDIF

    IF(PRESENT(nablavks_atom_set)) THEN
       IF(ASSOCIATED(epr_env%nablavks_atom_set)) THEN
          CALL deallocate_nablavks_atom_set(epr_env%nablavks_atom_set,error=error)
       ENDIF
       epr_env%nablavks_atom_set => nablavks_atom_set
    ENDIF

  END SUBROUTINE set_epr_env

! *****************************************************************************
  SUBROUTINE allocate_nablavks_atom_set(nablavks_atom_set,natom,error)

    TYPE(nablavks_atom_type), DIMENSION(:), &
      POINTER                                :: nablavks_atom_set
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_nablavks_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    ALLOCATE(nablavks_atom_set(natom), STAT=istat)
    CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)

    DO iat = 1,natom
       NULLIFY(nablavks_atom_set(iat)%nablavks_vec_rad_h)
       NULLIFY(nablavks_atom_set(iat)%nablavks_vec_rad_s)
    ENDDO
  END SUBROUTINE allocate_nablavks_atom_set

! *****************************************************************************
  SUBROUTINE deallocate_nablavks_atom_set(nablavks_atom_set,error)

    TYPE(nablavks_atom_type), DIMENSION(:), &
      POINTER                                :: nablavks_atom_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_nablavks_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iat, idir, istat, n, natom
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(nablavks_atom_set),cp_failure_level,routineP,error,failure)
    natom = SIZE(nablavks_atom_set)

    DO iat = 1,natom
       IF(ASSOCIATED(nablavks_atom_set(iat)%nablavks_vec_rad_h)) THEN
          IF(ASSOCIATED(nablavks_atom_set(iat)%nablavks_vec_rad_h(1,1)%r_coef)) THEN
             n = SIZE(nablavks_atom_set(iat)%nablavks_vec_rad_h,2)
             DO i=1,n
                DO idir = 1,3
                   DEALLOCATE(nablavks_atom_set(iat)%nablavks_vec_rad_h(idir,i)%r_coef,STAT=istat)
                   CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                   DEALLOCATE(nablavks_atom_set(iat)%nablavks_vec_rad_s(idir,i)%r_coef,STAT=istat)
                   CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                ENDDO
             ENDDO
          ENDIF   
          DEALLOCATE(nablavks_atom_set(iat)%nablavks_vec_rad_h,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(nablavks_atom_set(iat)%nablavks_vec_rad_s,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
       ENDIF
    ENDDO
    DEALLOCATE(nablavks_atom_set, STAT=istat)
    CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
  END SUBROUTINE deallocate_nablavks_atom_set

! *****************************************************************************
  SUBROUTINE allocate_jrho_atom_set(jrho_atom_set,natom,error)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho_atom_set
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_jrho_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    ALLOCATE(jrho_atom_set(natom), STAT=istat)
    CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
     
    DO iat = 1,natom
       NULLIFY(jrho_atom_set(iat)%cjc0_h)
       NULLIFY(jrho_atom_set(iat)%cjc0_s)
       NULLIFY(jrho_atom_set(iat)%cjc_h)
       NULLIFY(jrho_atom_set(iat)%cjc_s)
       NULLIFY(jrho_atom_set(iat)%cjc_ii_h)
       NULLIFY(jrho_atom_set(iat)%cjc_ii_s)
       NULLIFY(jrho_atom_set(iat)%cjc_iii_h)
       NULLIFY(jrho_atom_set(iat)%cjc_iii_s)
       NULLIFY(jrho_atom_set(iat)%jrho_rad_h)
       NULLIFY(jrho_atom_set(iat)%jrho_rad_s)
       NULLIFY(jrho_atom_set(iat)%jrho_vec_rad_h)
       NULLIFY(jrho_atom_set(iat)%jrho_vec_rad_s)
       NULLIFY(jrho_atom_set(iat)%jrho_h)
       NULLIFY(jrho_atom_set(iat)%jrho_s)
       NULLIFY(jrho_atom_set(iat)%jrho_a_h)
       NULLIFY(jrho_atom_set(iat)%jrho_a_s)
       NULLIFY(jrho_atom_set(iat)%jrho_b_h)
       NULLIFY(jrho_atom_set(iat)%jrho_b_s)
       NULLIFY(jrho_atom_set(iat)%jrho_a_h_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_a_s_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_b_h_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_b_s_ii)
       NULLIFY(jrho_atom_set(iat)%jrho_a_h_iii)
       NULLIFY(jrho_atom_set(iat)%jrho_a_s_iii)
       NULLIFY(jrho_atom_set(iat)%jrho_b_h_iii)
       NULLIFY(jrho_atom_set(iat)%jrho_b_s_iii)
    ENDDO
  END SUBROUTINE allocate_jrho_atom_set

! *****************************************************************************
  SUBROUTINE deallocate_jrho_atom_set(jrho_atom_set,error)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho_atom_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_jrho_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iat, iB, idir, istat, n, &
                                                natom
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(jrho_atom_set),cp_failure_level,routineP,error,failure)
    natom = SIZE(jrho_atom_set)

    DO iat = 1,natom

       IF(ASSOCIATED(jrho_atom_set(iat)%cjc_h)) THEN
          IF(ASSOCIATED(jrho_atom_set(iat)%cjc_h(1)%r_coef)) THEN
             n = SIZE(jrho_atom_set(iat)%cjc_h)
             DO i = 1,n
                DEALLOCATE(jrho_atom_set(iat)%cjc0_h(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%cjc0_s(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%cjc_h(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%cjc_s(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%cjc_ii_h(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%cjc_ii_s(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%cjc_iii_h(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%cjc_iii_s(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
             END DO
          END IF
          DEALLOCATE(jrho_atom_set(iat)%cjc0_h,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%cjc0_s,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%cjc_h,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%cjc_s,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%cjc_ii_h,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%cjc_ii_s,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%cjc_iii_h,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%cjc_iii_s,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
       END IF
       
       IF(ASSOCIATED(jrho_atom_set(iat)%jrho_a_h)) THEN
          IF(ASSOCIATED(jrho_atom_set(iat)%jrho_a_h(1)%r_coef)) THEN
             n = SIZE(jrho_atom_set(iat)%jrho_a_h)
             DO i = 1,n
                DEALLOCATE(jrho_atom_set(iat)%jrho_h(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_s(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_a_h(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_a_s(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_b_h(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_b_s(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_a_h_ii(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_a_s_ii(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_b_h_ii(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_b_s_ii(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_a_h_iii(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_a_s_iii(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_b_h_iii(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                DEALLOCATE(jrho_atom_set(iat)%jrho_b_s_iii(i)%r_coef,STAT=istat)
                CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
             END DO
          END IF
          DEALLOCATE(jrho_atom_set(iat)%jrho_h,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_s,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_a_h,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_a_s,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_b_h,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_b_s,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_a_h_ii,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_a_s_ii,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_b_h_ii,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_b_s_ii,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_a_h_iii,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_a_s_iii,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_b_h_iii,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_b_s_iii,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
       END IF

       IF(ASSOCIATED(jrho_atom_set(iat)%jrho_vec_rad_h)) THEN
          NULLIFY(jrho_atom_set(iat)%jrho_rad_h,jrho_atom_set(iat)%jrho_rad_s)
          n = SIZE(jrho_atom_set(iat)%jrho_vec_rad_h,3)
          DO i = 1,n
             DO iB=1,3
                DO idir = 1,3
                   DEALLOCATE(jrho_atom_set(iat)%jrho_vec_rad_h(idir,iB,i)%r_coef,STAT=istat)
                   CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                   DEALLOCATE(jrho_atom_set(iat)%jrho_vec_rad_s(idir,iB,i)%r_coef,STAT=istat)
                   CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
                END DO
             END DO
          END DO
          DEALLOCATE(jrho_atom_set(iat)%jrho_vec_rad_h,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_vec_rad_s,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
       ELSE
          n = SIZE(jrho_atom_set(iat)%jrho_rad_h,1)
          DO i = 1,n
             DEALLOCATE(jrho_atom_set(iat)%jrho_rad_h(i)%r_coef,STAT=istat)
             CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
             DEALLOCATE(jrho_atom_set(iat)%jrho_rad_s(i)%r_coef,STAT=istat)
             CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          END DO
          DEALLOCATE(jrho_atom_set(iat)%jrho_rad_h,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
          DEALLOCATE(jrho_atom_set(iat)%jrho_rad_s,STAT=istat)
          CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)
       END IF
    END DO
    DEALLOCATE(jrho_atom_set, STAT=istat)
    CPPrecondition(istat==0, cp_failure_level,routineP,error,failure)

  END SUBROUTINE deallocate_jrho_atom_set

! *****************************************************************************
  SUBROUTINE  allocate_jrho_atom_rad(jrho1_atom,ispin,nr,max_iso_not0,error)
  
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    INTEGER, INTENT(IN)                      :: ispin, nr, max_iso_not0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_jrho_atom_rad', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat
    LOGICAL                                  :: failure

    CALL timeset(routineN,"I"," ",handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(jrho1_atom),cp_failure_level,routineP,error,failure)

    IF(.NOT.failure) THEN
       ALLOCATE (jrho1_atom%jrho_h(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_h(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_s(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_s(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_a_h(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_a_h(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_a_s(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_a_s(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_b_h(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_b_h(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_b_s(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_b_s(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_a_h_ii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_a_h_ii(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_a_s_ii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_a_s_ii(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_b_h_ii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_b_h_ii(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_b_s_ii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_b_s_ii(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_a_h_iii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_a_h_iii(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_a_s_iii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_a_s_iii(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_b_h_iii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_b_h_iii(ispin)%r_coef = 0.0_dp
       ALLOCATE (jrho1_atom%jrho_b_s_iii(ispin)%r_coef(nr,max_iso_not0),STAT=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       jrho1_atom%jrho_b_s_iii(ispin)%r_coef = 0.0_dp
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE allocate_jrho_atom_rad

! *****************************************************************************
  SUBROUTINE set2zero_jrho_atom_rad(jrho1_atom,ispin,error)
    !
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set2zero_jrho_atom_rad', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!
!

    failure = .FALSE.
    !
    CPPrecondition(ASSOCIATED(jrho1_atom),cp_failure_level,routineP,error,failure)
    !
    jrho1_atom%jrho_h(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_s(ispin)%r_coef = 0.0_dp
    !
    jrho1_atom%jrho_a_h(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_a_s(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_h(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_s(ispin)%r_coef = 0.0_dp
    !
    jrho1_atom%jrho_a_h_ii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_a_s_ii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_h_ii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_s_ii(ispin)%r_coef = 0.0_dp
    !
    jrho1_atom%jrho_a_h_iii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_a_s_iii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_h_iii(ispin)%r_coef = 0.0_dp
    jrho1_atom%jrho_b_s_iii(ispin)%r_coef = 0.0_dp
    !
  END SUBROUTINE set2zero_jrho_atom_rad

! *****************************************************************************
  SUBROUTINE get_jrho_atom(jrho1_atom_set,iatom,cjc_h,cjc_s,cjc_ii_h,cjc_ii_s,&
             cjc_iii_h,cjc_iii_s,jrho_rad_h,jrho_rad_s,jrho_vec_rad_h,jrho_vec_rad_s,error)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(rho_atom_coeff), DIMENSION(:), &
      OPTIONAL, POINTER                      :: cjc_h, cjc_s, cjc_ii_h, &
                                                cjc_ii_s, cjc_iii_h, &
                                                cjc_iii_s, jrho_rad_h, &
                                                jrho_rad_s
    TYPE(rho_atom_coeff), &
      DIMENSION(:, :, :), OPTIONAL, POINTER  :: jrho_vec_rad_h, jrho_vec_rad_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_jrho_atom', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(jrho1_atom_set), cp_failure_level,routineP,error,failure)

    IF(PRESENT(cjc_h         )) cjc_h          => jrho1_atom_set(iatom)%cjc_h
    IF(PRESENT(cjc_s         )) cjc_s          => jrho1_atom_set(iatom)%cjc_s
    IF(PRESENT(cjc_ii_h      )) cjc_ii_h       => jrho1_atom_set(iatom)%cjc_ii_h
    IF(PRESENT(cjc_ii_s      )) cjc_ii_s       => jrho1_atom_set(iatom)%cjc_ii_s
    IF(PRESENT(cjc_iii_h     )) cjc_iii_h      => jrho1_atom_set(iatom)%cjc_iii_h
    IF(PRESENT(cjc_iii_s     )) cjc_iii_s      => jrho1_atom_set(iatom)%cjc_iii_s
    IF(PRESENT(jrho_rad_h    )) jrho_rad_h     => jrho1_atom_set(iatom)%jrho_rad_h
    IF(PRESENT(jrho_rad_s    )) jrho_rad_s     => jrho1_atom_set(iatom)%jrho_rad_s
    IF(PRESENT(jrho_vec_rad_h)) jrho_vec_rad_h => jrho1_atom_set(iatom)%jrho_vec_rad_h
    IF(PRESENT(jrho_vec_rad_s)) jrho_vec_rad_s => jrho1_atom_set(iatom)%jrho_vec_rad_s

  END SUBROUTINE get_jrho_atom

! *****************************************************************************
  SUBROUTINE init_jrho_atom_set(jrho1_atom_set,atomic_kind_set,nspins,store_current,error)

    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, INTENT(IN)                      :: nspins
    LOGICAL, INTENT(IN)                      :: store_current
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_jrho_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iat, iatom, iB, idir, ikind, ispin, istat, &
      max_iso_not0, maxso, na, nat, natom, nkind, nr, nset, nsotot
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, paw_atom
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(harmonics_atom_type), POINTER       :: harmonics

    CALL timeset(routineN,"I"," ",handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)

    IF(ASSOCIATED(jrho1_atom_set)) THEN
       CALL deallocate_jrho_atom_set(jrho1_atom_set,error=error)
    END IF

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom = natom)

    CALL allocate_jrho_atom_set(jrho1_atom_set,natom,error=error)

    nkind = SIZE(atomic_kind_set)

    DO ikind = 1,nkind

       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            orb_basis_set=orb_basis_set, &
                            atom_list=atom_list,natom=nat, &
                            paw_atom=paw_atom,&
                            harmonics=harmonics,&
                            grid_atom=grid_atom)

       na = grid_atom%ng_sphere 
       nr = grid_atom%nr 

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                              maxso=maxso, nset=nset)
       nsotot = maxso * nset 
       max_iso_not0 = harmonics%max_iso_not0
       DO iat = 1,nat
          iatom = atom_list(iat) 
          !*** allocate the radial density for each LM,for each atom ***

          IF(store_current) THEN
             ALLOCATE (jrho1_atom_set(iatom)%jrho_vec_rad_h(3,3,nspins),STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE (jrho1_atom_set(iatom)%jrho_vec_rad_s(3,3,nspins),STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             DO ispin = 1,nspins
                DO iB = 1,3
                   DO idir = 1,3
                      NULLIFY(jrho1_atom_set(iatom)%jrho_vec_rad_h(idir,iB,ispin)%r_coef)
                      NULLIFY(jrho1_atom_set(iatom)%jrho_vec_rad_s(idir,iB,ispin)%r_coef)
                      ALLOCATE(jrho1_atom_set(iatom)%jrho_vec_rad_h(idir,iB,ispin)%r_coef(nr,na),STAT=istat)
                      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                      ALLOCATE(jrho1_atom_set(iatom)%jrho_vec_rad_s(idir,iB,ispin)%r_coef(nr,na),STAT=istat)
                      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                   END DO
                END DO
             END DO  ! ispin
          ELSE
             ALLOCATE (jrho1_atom_set(iatom)%jrho_rad_h(nspins),STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE (jrho1_atom_set(iatom)%jrho_rad_s(nspins),STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
 
             DO ispin = 1,nspins
                NULLIFY(jrho1_atom_set(iatom)%jrho_rad_h(ispin)%r_coef)
                NULLIFY(jrho1_atom_set(iatom)%jrho_rad_s(ispin)%r_coef)

                ALLOCATE(jrho1_atom_set(iatom)%jrho_rad_h(ispin)%r_coef(nr,na),STAT=istat)
                CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                jrho1_atom_set(iatom)%jrho_rad_h(ispin)%r_coef = 0.0_dp
                ALLOCATE(jrho1_atom_set(iatom)%jrho_rad_s(ispin)%r_coef(nr,na),STAT=istat)
                CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                jrho1_atom_set(iatom)%jrho_rad_s(ispin)%r_coef = 0.0_dp
             END DO  ! ispin 
          END IF
          
          ALLOCATE (jrho1_atom_set(iatom)%jrho_h(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_s(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

          ALLOCATE (jrho1_atom_set(iatom)%jrho_a_h(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_a_s(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_b_h(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_b_s(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_a_h_ii(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_a_s_ii(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_b_s_ii(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_b_h_ii(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_a_h_iii(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_a_s_iii(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_b_s_iii(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%jrho_b_h_iii(nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc0_h(nspins), STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc0_s(nspins), STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc_h(nspins), STAT=istat) 
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc_s(nspins), STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc_ii_h(nspins), STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc_ii_s(nspins), STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc_iii_h(nspins), STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (jrho1_atom_set(iatom)%cjc_iii_s(nspins), STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

          DO ispin = 1,nspins
             NULLIFY(jrho1_atom_set(iatom)%jrho_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_s(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_a_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_a_s(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_b_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_b_s(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_a_h_ii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_a_s_ii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_b_h_ii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_b_s_ii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_a_h_iii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_a_s_iii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_b_h_iii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%jrho_b_s_iii(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc0_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc0_s(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef)
             NULLIFY(jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef)

             IF(paw_atom) THEN
                ALLOCATE(jrho1_atom_set(iatom)%cjc0_h(ispin)%r_coef(1:nsotot,1:nsotot),&
                     &   STAT=istat)
                CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(jrho1_atom_set(iatom)%cjc0_s(ispin)%r_coef(1:nsotot,1:nsotot),&
                     &   STAT=istat)
                CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(jrho1_atom_set(iatom)%cjc_h(ispin)%r_coef(1:nsotot,1:nsotot),&
                     &   STAT=istat)
                CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(jrho1_atom_set(iatom)%cjc_s(ispin)%r_coef(1:nsotot,1:nsotot),&
                     &   STAT=istat)
                CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(jrho1_atom_set(iatom)%cjc_ii_h(ispin)%r_coef(1:nsotot,1:nsotot),&
                     &   STAT=istat)
                CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(jrho1_atom_set(iatom)%cjc_ii_s(ispin)%r_coef(1:nsotot,1:nsotot),&
                     &   STAT=istat)
                CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(jrho1_atom_set(iatom)%cjc_iii_h(ispin)%r_coef(1:nsotot,1:nsotot),&
                     &   STAT=istat)
                CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(jrho1_atom_set(iatom)%cjc_iii_s(ispin)%r_coef(1:nsotot,1:nsotot),&
                     &   STAT=istat)
                CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             ENDIF

          ENDDO ! ispin

       END DO  ! iat

    END DO  ! ikind

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE init_jrho_atom_set

! *****************************************************************************
  SUBROUTINE init_nablavks_atom_set(nablavks_atom_set,atomic_kind_set,nspins,error)

    TYPE(nablavks_atom_type), DIMENSION(:), &
      POINTER                                :: nablavks_atom_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_nablavks_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iat, iatom, idir, ikind, ispin, istat, max_iso_not0, &
      maxso, na, nat, natom, nkind, nr, nset, nsotot
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, paw_atom
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(harmonics_atom_type), POINTER       :: harmonics

    CALL timeset(routineN,"I"," ",handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)

    IF(ASSOCIATED(nablavks_atom_set)) THEN
      CALL deallocate_nablavks_atom_set(nablavks_atom_set,error=error)
    END IF

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom = natom)

    CALL allocate_nablavks_atom_set(nablavks_atom_set,natom,error=error)

    nkind = SIZE(atomic_kind_set)

    DO ikind = 1,nkind

       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            orb_basis_set=orb_basis_set, &
                            atom_list=atom_list,natom=nat, &
                            paw_atom=paw_atom,&
                            harmonics=harmonics,&
                            grid_atom=grid_atom)

       na = grid_atom%ng_sphere 
       nr = grid_atom%nr 

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                              maxso=maxso, nset=nset)
       nsotot = maxso * nset 
       max_iso_not0 = harmonics%max_iso_not0
       DO iat = 1,nat
          iatom = atom_list(iat) 
          !*** allocate the radial density for each LM,for each atom ***

            ALLOCATE (nablavks_atom_set(iatom)%nablavks_vec_rad_h(3,nspins),STAT=istat)
            CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
            ALLOCATE (nablavks_atom_set(iatom)%nablavks_vec_rad_s(3,nspins),STAT=istat)
            CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
            DO ispin = 1,nspins
              DO idir = 1,3
                NULLIFY(nablavks_atom_set(iatom)%nablavks_vec_rad_h(idir,ispin)%r_coef)
                NULLIFY(nablavks_atom_set(iatom)%nablavks_vec_rad_s(idir,ispin)%r_coef)
                ALLOCATE(nablavks_atom_set(iatom)%nablavks_vec_rad_h(idir,ispin)%r_coef(nr,na),STAT=istat)
                CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(nablavks_atom_set(iatom)%nablavks_vec_rad_s(idir,ispin)%r_coef(nr,na),STAT=istat)
                CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
              END DO
            END DO  ! ispin
       END DO  ! iat

    END DO  ! ikind

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE init_nablavks_atom_set

END MODULE qs_linres_types

