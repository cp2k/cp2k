!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types and set/get functions for auxiliary denisty matrix methods
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE admm_types
  USE bibliography,                    ONLY: Guidon2010,&
                                             cite_reference
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: do_admm_block_density_matrix
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_release,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  admm_env_create, admm_env_release, admm_type, eigvals_p_type, eigvals_type,&
          admm_create_block_list, admm_block


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_types'

  TYPE eigvals_type
    REAL(dp), DIMENSION(:), POINTER          :: DATA
  END TYPE

  TYPE eigvals_p_type
    TYPE(eigvals_type), POINTER              :: eigvals
  END TYPE

  TYPE admm_block
    INTEGER, DIMENSION(:), POINTER            :: list
  END TYPE


! *****************************************************************************
!> \brief stores some data used in wavefunction fitting
!> \param S overlap matrix for auxiliary fit basis set
!> \param P overlap matrix for mixed aux_fit/orb basis set
!> \param A contains inv(S)*P
!> \param B contains transpose(P)*inv(S)*P = transpose(P)*A
!> \param lambda contains transpose(mo_coeff_aux_fit)*B*mo_coeff_aux_fit
!> \param lambda_inv_sqrt contains inv(SQRT(lambda))
!> \param R contains eigenvectors of lambda
!> \param work_aux_aux temporary matrix
!> \param work_orb_nmo temporary matrix
!> \param work_nmo_nmo1 temporary matrix
!> \param work_nmo_nmo2 temporary matrix
!> \param work_aux_nmo temporary matrix
!> \param H contains KS_matrix * mo_coeff for auxiliary basis set
!> \param K contains KS matrix for auxiliary basis set
!> \param M contains matrix holding the 2nd order residues
!> \param nao_orb number of atomic orbitals in orb basis set
!> \param nao_aux_fit number of atomic orbitals in aux basis set
!> \param nmo number of molecular orbitals per spin
!> \param eigvals_lamda eigenvalues of lambda matrix
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************

  TYPE admm_type
    TYPE(cp_fm_type), POINTER                :: S_inv, S
    TYPE(cp_fm_type), POINTER                :: Q
    TYPE(cp_fm_type), POINTER                :: A
    TYPE(cp_fm_type), POINTER                :: B
    TYPE(cp_fm_type), POINTER                :: work_orb_orb, work_orb_orb2
    TYPE(cp_fm_type), POINTER                :: work_aux_orb, work_aux_orb2
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: lambda, lambda_inv
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: lambda_inv_sqrt
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: R, R_purify
    TYPE(cp_fm_type), POINTER                :: work_aux_aux, work_aux_aux2, work_aux_aux3
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: work_orb_nmo
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: work_nmo_nmo1
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: R_schur_R_t
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: work_nmo_nmo2
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: work_aux_nmo
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: work_aux_nmo2
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: H, H_corr
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs_tmp
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: K
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: M, M_purify, P_to_be_purified
    INTEGER                                  :: nao_orb, nao_aux_fit, nmo(2)
    TYPE(eigvals_p_type), DIMENSION(:), &
      POINTER                                :: eigvals_lambda, eigvals_P_to_be_purified
    TYPE(section_vals_type), POINTER         :: xc_section_primary
    TYPE(section_vals_type), POINTER         :: xc_section_aux
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: lambda_inv2, C_hat, P_tilde
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: ks_to_be_merged
    INTEGER                                  :: method_id
    INTEGER                                  :: purification_method
    INTEGER                                  :: block_purification_method
    INTEGER                                  :: block_projection_method
    TYPE(admm_block), DIMENSION(:), &
      POINTER                                :: blocks
    INTEGER, DIMENSION(:,:), POINTER         :: block_map
  END TYPE



  CONTAINS

! *****************************************************************************
!> \brief creates ADDE environment, initializes the basic types
!>
!> \param mos the MO's of the orbital basis set
!> \param mos_aux_fit the MO's of the auxiliary fitting basis set
!> \param para_env The parallel env
!> \param admm_env The ADMM env
!> \param para_env The parallel env
!> \param error
!>
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************

  SUBROUTINE admm_env_create(mos, mos_aux_fit, para_env, admm_env, error)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, istat, nao_aux_fit, &
                                                nao_orb, nmo, nspins
    LOGICAL                                  :: failure
    TYPE(cp_fm_struct_type), POINTER :: fm_struct_aux_aux, fm_struct_aux_nmo, &
      fm_struct_aux_orb, fm_struct_nmo_nmo, fm_struct_orb_aux, &
      fm_struct_orb_nmo, fm_struct_orb_orb
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL cite_reference(Guidon2010)

    ALLOCATE(admm_env, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    nspins = SIZE(mos)
    CALL get_mo_set(mos(1)%mo_set,mo_coeff=mo_coeff, nmo=nmo, nao=nao_orb)
    CALL get_mo_set(mos_aux_fit(1)%mo_set,nao=nao_aux_fit)
    admm_env%nmo = 0
    admm_env%nao_aux_fit = nao_aux_fit
    admm_env%nao_orb = nao_orb
    CALL cp_fm_struct_create(fm_struct_aux_aux,&
                             context=mo_coeff%matrix_struct%context,&
                             nrow_global=nao_aux_fit,&
                             ncol_global=nao_aux_fit,&
                             para_env=para_env,&
                             error=error)
    CALL cp_fm_struct_create(fm_struct_aux_orb,&
                             context=mo_coeff%matrix_struct%context,&
                             nrow_global=nao_aux_fit,&
                             ncol_global=nao_orb,&
                             para_env=para_env,&
                             error=error)
    CALL cp_fm_struct_create(fm_struct_orb_aux,&
                             context=mo_coeff%matrix_struct%context,&
                             nrow_global=nao_orb,&
                             ncol_global=nao_aux_fit,&
                             para_env=para_env,&
                             error=error)
    CALL cp_fm_struct_create(fm_struct_orb_orb,&
                             context=mo_coeff%matrix_struct%context,&
                             nrow_global=nao_orb,&
                             ncol_global=nao_orb,&
                             para_env=para_env,&
                             error=error)

    CALL cp_fm_create(admm_env%S,fm_struct_aux_aux,name="aux_fit_overlap",error=error)
    CALL cp_fm_create(admm_env%S_inv,fm_struct_aux_aux,name="aux_fit_overlap_inv",error=error)
    CALL cp_fm_create(admm_env%Q,fm_struct_aux_orb,name="mixed_overlap",error=error)
    CALL cp_fm_create(admm_env%A,fm_struct_aux_orb,name="work_A",error=error)
    CALL cp_fm_create(admm_env%B,fm_struct_orb_orb,name="work_B",error=error)
    CALL cp_fm_create(admm_env%work_orb_orb,fm_struct_orb_orb,name="work_orb_orb",error=error)
    CALL cp_fm_create(admm_env%work_orb_orb2,fm_struct_orb_orb,name="work_orb_orb",error=error)
    CALL cp_fm_create(admm_env%work_aux_orb,fm_struct_aux_orb,name="work_aux_orb",error=error)
    CALL cp_fm_create(admm_env%work_aux_orb2,fm_struct_aux_orb,name="work_aux_orb2",error=error)

    CALL cp_fm_create(admm_env%work_aux_aux,fm_struct_aux_aux,name="work_aux_aux",error=error)
    CALL cp_fm_create(admm_env%work_aux_aux2,fm_struct_aux_aux,name="work_aux_aux2",error=error)
    CALL cp_fm_create(admm_env%work_aux_aux3,fm_struct_aux_aux,name="work_aux_aux3",error=error)

    ALLOCATE(admm_env%lambda_inv(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%lambda(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%lambda_inv_sqrt(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%R(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%R_purify(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%work_orb_nmo(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%work_nmo_nmo1(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%R_schur_R_t(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%work_nmo_nmo2(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%eigvals_lambda(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%eigvals_P_to_be_purified(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%H(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%K(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%M(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%M_purify(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%P_to_be_purified(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%work_aux_nmo(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%work_aux_nmo2(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%mo_derivs_tmp(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%H_corr(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%ks_to_be_merged(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%lambda_inv2(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%C_hat(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(admm_env%P_tilde(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO ispin = 1,nspins
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff, nmo=nmo)
      admm_env%nmo(ispin)=nmo
      CALL cp_fm_struct_create(fm_struct_aux_nmo,&
                               context=mo_coeff%matrix_struct%context,&
                               nrow_global=nao_aux_fit,&
                               ncol_global=nmo,&
                               para_env=para_env,&
                               error=error)
      CALL cp_fm_struct_create(fm_struct_orb_nmo,&
                               context=mo_coeff%matrix_struct%context,&
                               nrow_global=nao_orb,&
                               ncol_global=nmo,&
                               para_env=para_env,&
                               error=error)
      CALL cp_fm_struct_create(fm_struct_nmo_nmo,&
                               context=mo_coeff%matrix_struct%context,&
                               nrow_global=nmo,&
                               ncol_global=nmo,&
                               para_env=para_env,&
                               error=error)

      CALL cp_fm_create(admm_env%work_orb_nmo(ispin)%matrix,fm_struct_orb_nmo,name="work_orb_nmo",error=error)
      CALL cp_fm_create(admm_env%work_nmo_nmo1(ispin)%matrix,fm_struct_nmo_nmo,name="work_nmo_nmo1",error=error)
      CALL cp_fm_create(admm_env%R_schur_R_t(ispin)%matrix,fm_struct_nmo_nmo,name="R_schur_R_t",error=error)
      CALL cp_fm_create(admm_env%work_nmo_nmo2(ispin)%matrix,fm_struct_nmo_nmo,name="work_nmo_nmo2",error=error)
      CALL cp_fm_create(admm_env%lambda(ispin)%matrix,fm_struct_nmo_nmo,name="lambda",error=error)
      CALL cp_fm_create(admm_env%lambda_inv(ispin)%matrix,fm_struct_nmo_nmo,name="lambda_inv",error=error)
      CALL cp_fm_create(admm_env%lambda_inv_sqrt(ispin)%matrix,fm_struct_nmo_nmo,name="lambda_inv_sqrt",error=error)
      CALL cp_fm_create(admm_env%R(ispin)%matrix,fm_struct_nmo_nmo,name="R",error=error)
      CALL cp_fm_create(admm_env%R_purify(ispin)%matrix,fm_struct_aux_aux,name="R_purify",error=error)
      CALL cp_fm_create(admm_env%K(ispin)%matrix,fm_struct_aux_aux,name="K",error=error)
      CALL cp_fm_create(admm_env%H(ispin)%matrix,fm_struct_aux_nmo,name="H",error=error)
      CALL cp_fm_create(admm_env%H_corr(ispin)%matrix,fm_struct_orb_orb,name="H_corr",error=error)
      CALL cp_fm_create(admm_env%M(ispin)%matrix,fm_struct_nmo_nmo,name="M",error=error)
      CALL cp_fm_create(admm_env%M_purify(ispin)%matrix,fm_struct_aux_aux,name="M aux",error=error)
      CALL cp_fm_create(admm_env%P_to_be_purified(ispin)%matrix,fm_struct_aux_aux,name="P_to_be_purified",error=error)
      CALL cp_fm_create(admm_env%work_aux_nmo(ispin)%matrix,fm_struct_aux_nmo,name="work_aux_nmo",error=error)
      CALL cp_fm_create(admm_env%work_aux_nmo2(ispin)%matrix,fm_struct_aux_nmo,name="work_aux_nmo2",error=error)
      CALL cp_fm_create(admm_env%mo_derivs_tmp(ispin)%matrix,fm_struct_orb_nmo,name="mo_derivs_tmp",error=error)

      CALL cp_fm_create(admm_env%lambda_inv2(ispin)%matrix,fm_struct_nmo_nmo,name="lambda_inv2",error=error)
      CALL cp_fm_create(admm_env%C_hat(ispin)%matrix,fm_struct_aux_nmo,name="C_hat",error=error)
      CALL cp_fm_create(admm_env%P_tilde(ispin)%matrix,fm_struct_aux_aux,name="P_tilde",error=error)

      CALL cp_fm_create(admm_env%ks_to_be_merged(ispin)%matrix,fm_struct_orb_orb,name="KS_to_be_merged ",error=error)

      ALLOCATE(admm_env%eigvals_lambda(ispin)%eigvals, STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(admm_env%eigvals_P_to_be_purified(ispin)%eigvals, STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(admm_env%eigvals_lambda(ispin)%eigvals%data(nmo), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(nao_aux_fit), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      admm_env%eigvals_lambda(ispin)%eigvals%data = 0.0_dp
      admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data = 0.0_dp
      CALL cp_fm_struct_release(fm_struct_aux_nmo,error=error)
      CALL cp_fm_struct_release(fm_struct_orb_nmo,error=error)
      CALL cp_fm_struct_release(fm_struct_nmo_nmo,error=error)
    END DO

    CALL cp_fm_struct_release(fm_struct_aux_aux,error=error)
    CALL cp_fm_struct_release(fm_struct_aux_orb,error=error)
    CALL cp_fm_struct_release(fm_struct_orb_aux,error=error)
    CALL cp_fm_struct_release(fm_struct_orb_orb,error=error)
  END SUBROUTINE admm_env_create


! *****************************************************************************
!> \brief releases the ADMM environment, cleans up all types
!>
!> \param admm_env The ADMM env
!> \param error
!>
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE admm_env_release(admm_env, error)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, istat
    LOGICAL                                  :: failure

    CALL cp_fm_release(admm_env%S, error=error)
    CALL cp_fm_release(admm_env%S_inv, error=error)
    CALL cp_fm_release(admm_env%Q, error=error)
    CALL cp_fm_release(admm_env%A, error=error)
    CALL cp_fm_release(admm_env%B, error=error)
    CALL cp_fm_release(admm_env%work_orb_orb, error=error)
    CALL cp_fm_release(admm_env%work_orb_orb2, error=error)
    CALL cp_fm_release(admm_env%work_aux_aux, error=error)
    CALL cp_fm_release(admm_env%work_aux_aux2, error=error)
    CALL cp_fm_release(admm_env%work_aux_aux3, error=error)
    CALL cp_fm_release(admm_env%work_aux_orb, error=error)
    CALL cp_fm_release(admm_env%work_aux_orb2, error=error)
    DO ispin = 1,SIZE(admm_env%lambda)
      CALL cp_fm_release(admm_env%lambda(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%lambda_inv(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%lambda_inv_sqrt(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%lambda_inv2(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%C_hat(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%P_tilde(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%R(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%R_purify(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%H(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%H_corr(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%K(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%M(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%M_purify(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%P_to_be_purified(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%work_orb_nmo(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%work_nmo_nmo1(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%R_schur_R_t(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%work_nmo_nmo2(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%work_aux_nmo(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%work_aux_nmo2(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%mo_derivs_tmp(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%ks_to_be_merged(ispin)%matrix, error=error)
      CALL cp_fm_release(admm_env%lambda_inv2(ispin)%matrix, error=error)
      DEALLOCATE(admm_env%eigvals_lambda(ispin)%eigvals%data,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(admm_env%eigvals_lambda(ispin)%eigvals,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(admm_env%eigvals_P_to_be_purified(ispin)%eigvals,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END DO
    DEALLOCATE(admm_env%eigvals_lambda,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%eigvals_P_to_be_purified,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%lambda,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%lambda_inv,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%lambda_inv_sqrt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%R,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%R_purify,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%M,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%M_purify,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%P_to_be_purified,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%H,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%H_corr,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%K,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%work_orb_nmo,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%work_nmo_nmo1,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%R_schur_R_t,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%work_nmo_nmo2,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%work_aux_nmo,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%work_aux_nmo2,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%mo_derivs_tmp,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%ks_to_be_merged,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%lambda_inv2,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%C_hat,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(admm_env%P_tilde,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF( admm_env%method_id == do_admm_block_density_matrix) THEN
      DEALLOCATE(admm_env%block_map, STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    
    IF(ASSOCIATED(admm_env%xc_section_primary))&
       CALL section_vals_release(admm_env%xc_section_primary,error)
    IF(ASSOCIATED(admm_env%xc_section_aux))&
       CALL section_vals_release(admm_env%xc_section_aux,error)

    DEALLOCATE(admm_env, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE admm_env_release


  SUBROUTINE admm_create_block_list(admm_block_section, admm_env, natom, error)
    TYPE(section_vals_type), POINTER         :: admm_block_section
    TYPE(admm_type), POINTER                 :: admm_env
    INTEGER                                  :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_create_block_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, irep, j, jatom, &
                                                list_size, n_rep, stat
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: list_section

    NULLIFY(list_section)
    list_section => section_vals_get_subs_vals(admm_block_section,"BLOCK",error=error)
    CALL section_vals_get(list_section,n_repetition=n_rep,error=error)

    ALLOCATE(admm_env%blocks(n_rep), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO irep = 1, n_rep
      CALL section_vals_val_get(list_section,"LIST", i_rep_section=irep, &
                   i_vals=tmplist, error=error)
      list_size = SIZE(tmplist)
      ALLOCATE(admm_env%blocks(irep)%list(list_size), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      admm_env%blocks(irep)%list = tmplist
    END DO

    ALLOCATE(admm_env%block_map(natom,natom), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    admm_env%block_map = 0
    DO irep = 1,n_rep
      DO i = 1,SIZE(admm_env%blocks(irep)%list)
        iatom = admm_env%blocks(irep)%list(i)
        DO j = 1,SIZE(admm_env%blocks(irep)%list)
          jatom = admm_env%blocks(irep)%list(j)
          admm_env%block_map(iatom,jatom) = 1
        END DO
      END DO
    END DO

    DO irep = 1, n_rep
      DEALLOCATE(admm_env%blocks(irep)%list, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    DEALLOCATE(admm_env%blocks, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


  END SUBROUTINE admm_create_block_list

END MODULE admm_types

