!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types and set/get functions for wfn_fitting 
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE wfn_fitting_types
  USE ai_overlap_new,                  ONLY: overlap
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  wfn_fitting_env_create, wfn_fitting_env_release, wfn_fitting_type, eigvals_p_type, eigvals_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wfn_fitting_types'

  TYPE eigvals_type
    REAL(dp), DIMENSION(:), POINTER          :: DATA
  END TYPE

  TYPE eigvals_p_type
    TYPE(eigvals_type), POINTER              :: eigvals
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

  TYPE wfn_fitting_type
    TYPE(cp_fm_type), POINTER                :: S
    TYPE(cp_fm_type), POINTER                :: P
    TYPE(cp_fm_type), POINTER                :: A
    TYPE(cp_fm_type), POINTER                :: B
    TYPE(cp_fm_type), POINTER                :: work_orb_orb
    TYPE(cp_fm_type), POINTER                :: work_aux_orb
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: lambda
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: lambda_inv_sqrt
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: R
    TYPE(cp_fm_type), POINTER                :: work_aux_aux
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
      POINTER                                :: H
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs_tmp
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: K
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: M
    INTEGER                                  :: nao_orb, nao_aux_fit, nmo(2)
    TYPE(eigvals_p_type), DIMENSION(:), &
      POINTER                                :: eigvals_lambda
    TYPE(section_vals_type), POINTER         :: xc_section_primary
    TYPE(section_vals_type), POINTER         :: xc_section_aux
  END TYPE

  
  CONTAINS

! *****************************************************************************
!> \brief creates wfn fitting environment, initializes the basic types
!>
!> \param mos the MO's of the orbital basis set
!> \param mos_aux_fit the MO's of the auxiliary fitting basis set
!> \param para_env The parallel env
!> \param wfn_fitting_env The wavefunction fitting env
!> \param para_env The parallel env
!> \param error
!>
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************

  SUBROUTINE wfn_fitting_env_create(mos, mos_aux_fit, para_env, wfn_fitting_env, error)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(wfn_fitting_type), POINTER          :: wfn_fitting_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'wfn_fitting_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, istat, nao_aux_fit, &
                                                nao_orb, nmo, nspins
    LOGICAL                                  :: failure
    TYPE(cp_fm_struct_type), POINTER :: fm_struct_aux_aux, fm_struct_aux_nmo, &
      fm_struct_aux_orb, fm_struct_nmo_nmo, fm_struct_orb_aux, &
      fm_struct_orb_nmo, fm_struct_orb_orb
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    ALLOCATE(wfn_fitting_env, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    nspins = SIZE(mos)
    CALL get_mo_set(mos(1)%mo_set,mo_coeff=mo_coeff, nmo=nmo, nao=nao_orb)
    CALL get_mo_set(mos_aux_fit(1)%mo_set,nao=nao_aux_fit)
    wfn_fitting_env%nmo = 0
    wfn_fitting_env%nao_aux_fit = nao_aux_fit
    wfn_fitting_env%nao_orb = nao_orb
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

    CALL cp_fm_create(wfn_fitting_env%S,fm_struct_aux_aux,name="aux_fit_overlap",error=error)
    CALL cp_fm_create(wfn_fitting_env%P,fm_struct_aux_orb,name="mixed_overlap",error=error)
    CALL cp_fm_create(wfn_fitting_env%A,fm_struct_aux_orb,name="work_A",error=error)
    CALL cp_fm_create(wfn_fitting_env%B,fm_struct_orb_orb,name="work_B",error=error)
    CALL cp_fm_create(wfn_fitting_env%work_orb_orb,fm_struct_orb_orb,name="work_orb_orb",error=error)
    CALL cp_fm_create(wfn_fitting_env%work_aux_orb,fm_struct_aux_orb,name="work_aux_orb",error=error)

    CALL cp_fm_create(wfn_fitting_env%work_aux_aux,fm_struct_aux_aux,name="work_aux_aux",error=error)

    ALLOCATE(wfn_fitting_env%lambda(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%lambda_inv_sqrt(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%R(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%work_orb_nmo(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%work_nmo_nmo1(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%R_schur_R_t(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%work_nmo_nmo2(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%eigvals_lambda(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%H(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%K(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%M(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%work_aux_nmo(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%work_aux_nmo2(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wfn_fitting_env%mo_derivs_tmp(nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)




    DO ispin = 1,nspins
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff, nmo=nmo)
      wfn_fitting_env%nmo(ispin)=nmo
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

      CALL cp_fm_create(wfn_fitting_env%work_orb_nmo(ispin)%matrix,fm_struct_orb_nmo,name="work_orb_nmo",error=error)
      CALL cp_fm_create(wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,fm_struct_nmo_nmo,name="work_nmo_nmo1",error=error)
      CALL cp_fm_create(wfn_fitting_env%R_schur_R_t(ispin)%matrix,fm_struct_nmo_nmo,name="R_schur_R_t",error=error)
      CALL cp_fm_create(wfn_fitting_env%work_nmo_nmo2(ispin)%matrix,fm_struct_nmo_nmo,name="work_nmo_nmo2",error=error)
      CALL cp_fm_create(wfn_fitting_env%lambda(ispin)%matrix,fm_struct_nmo_nmo,name="lambda",error=error)
      CALL cp_fm_create(wfn_fitting_env%lambda_inv_sqrt(ispin)%matrix,fm_struct_nmo_nmo,name="lambda_inv_sqrt",error=error)
      CALL cp_fm_create(wfn_fitting_env%R(ispin)%matrix,fm_struct_nmo_nmo,name="R",error=error)
      CALL cp_fm_create(wfn_fitting_env%K(ispin)%matrix,fm_struct_aux_aux,name="K",error=error)
      CALL cp_fm_create(wfn_fitting_env%H(ispin)%matrix,fm_struct_aux_nmo,name="H",error=error)
      CALL cp_fm_create(wfn_fitting_env%M(ispin)%matrix,fm_struct_nmo_nmo,name="M",error=error)
      CALL cp_fm_create(wfn_fitting_env%work_aux_nmo(ispin)%matrix,fm_struct_aux_nmo,name="work_aux_nmo",error=error)
      CALL cp_fm_create(wfn_fitting_env%work_aux_nmo2(ispin)%matrix,fm_struct_aux_nmo,name="work_aux_nmo2",error=error)
      CALL cp_fm_create(wfn_fitting_env%mo_derivs_tmp(ispin)%matrix,fm_struct_orb_nmo,name="mo_derivs_tmp",error=error)
      ALLOCATE(wfn_fitting_env%eigvals_lambda(ispin)%eigvals, STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(wfn_fitting_env%eigvals_lambda(ispin)%eigvals%data(nmo), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      wfn_fitting_env%eigvals_lambda(ispin)%eigvals%data = 0.0_dp
      CALL cp_fm_struct_release(fm_struct_aux_nmo,error=error)
      CALL cp_fm_struct_release(fm_struct_orb_nmo,error=error)
      CALL cp_fm_struct_release(fm_struct_nmo_nmo,error=error)
    END DO

    CALL cp_fm_struct_release(fm_struct_aux_aux,error=error)
    CALL cp_fm_struct_release(fm_struct_aux_orb,error=error)
    CALL cp_fm_struct_release(fm_struct_orb_aux,error=error)
    CALL cp_fm_struct_release(fm_struct_orb_orb,error=error) 
  END SUBROUTINE wfn_fitting_env_create


! *****************************************************************************
!> \brief releases the wfn fitting environment, cleans up all types
!>
!> \param wfn_fitting_env The wavefunction fitting env
!> \param error
!>
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE wfn_fitting_env_release(wfn_fitting_env, error)
    TYPE(wfn_fitting_type), POINTER          :: wfn_fitting_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'wfn_fitting_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, istat
    LOGICAL                                  :: failure

    CALL cp_fm_release(wfn_fitting_env%S, error=error)
    CALL cp_fm_release(wfn_fitting_env%P, error=error)
    CALL cp_fm_release(wfn_fitting_env%A, error=error)
    CALL cp_fm_release(wfn_fitting_env%B, error=error)
    CALL cp_fm_release(wfn_fitting_env%work_orb_orb, error=error)
    CALL cp_fm_release(wfn_fitting_env%work_aux_aux, error=error)
    CALL cp_fm_release(wfn_fitting_env%work_aux_orb, error=error)
    DO ispin = 1,SIZE(wfn_fitting_env%lambda)
      CALL cp_fm_release(wfn_fitting_env%lambda(ispin)%matrix, error=error)
      CALL cp_fm_release(wfn_fitting_env%lambda_inv_sqrt(ispin)%matrix, error=error)
      CALL cp_fm_release(wfn_fitting_env%R(ispin)%matrix, error=error)
      CALL cp_fm_release(wfn_fitting_env%H(ispin)%matrix, error=error)
      CALL cp_fm_release(wfn_fitting_env%K(ispin)%matrix, error=error)
      CALL cp_fm_release(wfn_fitting_env%M(ispin)%matrix, error=error)
      CALL cp_fm_release(wfn_fitting_env%work_orb_nmo(ispin)%matrix, error=error)
      CALL cp_fm_release(wfn_fitting_env%work_nmo_nmo1(ispin)%matrix, error=error)
      CALL cp_fm_release(wfn_fitting_env%R_schur_R_t(ispin)%matrix, error=error)
      CALL cp_fm_release(wfn_fitting_env%work_nmo_nmo2(ispin)%matrix, error=error)
      CALL cp_fm_release(wfn_fitting_env%work_aux_nmo(ispin)%matrix, error=error)
      CALL cp_fm_release(wfn_fitting_env%work_aux_nmo2(ispin)%matrix, error=error)
      CALL cp_fm_release(wfn_fitting_env%mo_derivs_tmp(ispin)%matrix, error=error)
      DEALLOCATE(wfn_fitting_env%eigvals_lambda(ispin)%eigvals%data,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(wfn_fitting_env%eigvals_lambda(ispin)%eigvals,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END DO
    DEALLOCATE(wfn_fitting_env%eigvals_lambda,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%lambda,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%lambda_inv_sqrt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%R,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%M,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%H,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%K,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%work_orb_nmo,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%work_nmo_nmo1,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%R_schur_R_t,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%work_nmo_nmo2,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%work_aux_nmo,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%work_aux_nmo2,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wfn_fitting_env%mo_derivs_tmp,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)


    CALL section_vals_release(wfn_fitting_env%xc_section_primary,error)
    CALL section_vals_release(wfn_fitting_env%xc_section_aux,error)

    DEALLOCATE(wfn_fitting_env, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE wfn_fitting_env_release

END MODULE wfn_fitting_types

