!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains methods used in the context of density fitting 
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE wfn_fitting_methods 
  
  USE ai_overlap_new,                  ONLY: overlap
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_schur_product,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm
  USE kinds,                           ONLY: dp
  USE mathconstants
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE wfn_fitting_types,               ONLY: wfn_fitting_env_create,&
                                             wfn_fitting_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC fit_mo_coeffs, merge_ks_times_mo_coeff

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wfn_fitting_methods'

!***
  
  CONTAINS

! *****************************************************************************
!> \brief Calculates the MO coefficients for the auxiliary fitting basis set  
!>        by minimizing int (psi_i - psi_aux_i)^2 using Lagrangian Multipliers
!>      
!> \param wfn_fitting_env The wavefunction fitting env
!> \param para_env The parallel env  
!> \param matrix_s_aux_fit the overlap matrix of the auxiliary fitting basis set
!> \param matrix_s_mixed the mixed overlap matrix of the auxiliary fitting basis
!>        set and the orbital basis set
!> \param mos the MO's of the orbital basis set
!> \param mos_aux_fit the MO's of the auxiliary fitting basis set
!> \param geometry_did_change flag to indicate if the geomtry changed
!> \param error 
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE fit_mo_coeffs(wfn_fitting_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                           mos, mos_aux_fit, geometry_did_change, error) 

    TYPE(wfn_fitting_type), POINTER          :: wfn_fitting_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    LOGICAL, INTENT(IN)                      :: geometry_did_change
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fit_mo_coeffs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ispin, istat, nao_aux_fit, &
                                                nao_orb, nmo, nspins
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eig_work
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit

    IF (.NOT.(ASSOCIATED(wfn_fitting_env) )) THEN
      CALL wfn_fitting_env_create(mos, mos_aux_fit, &
                                  para_env, wfn_fitting_env,&
                                  error)
    END IF

    nao_aux_fit = wfn_fitting_env%nao_aux_fit
    nao_orb = wfn_fitting_env%nao_orb
    nspins = SIZE(mos)


    ! *** This part only depends on overlap matrices ==> needs only to be calculated if the geometry changed

    IF( geometry_did_change ) THEN
      CALL copy_sm_to_fm(matrix_s_aux_fit(1)%matrix,wfn_fitting_env%S,error)
      CALL cp_fm_upper_to_full(wfn_fitting_env%S,wfn_fitting_env%work_aux_aux,error=error)
      CALL copy_sm_to_fm(matrix_s_mixed(1)%matrix,wfn_fitting_env%P,error)

      !! Calculate S'_inverse
      CALL cp_fm_cholesky_decompose(wfn_fitting_env%S,error=error)
      CALL cp_fm_cholesky_invert(wfn_fitting_env%S,error=error) 
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(wfn_fitting_env%S,wfn_fitting_env%work_aux_aux,error=error)
      !! Calculate A=S'^(-1)*P
      CALL cp_fm_gemm('N','N',nao_aux_fit,nao_orb,nao_aux_fit,&
                    1.0_dp,wfn_fitting_env%S,wfn_fitting_env%P,0.0_dp,&
                    wfn_fitting_env%A,error)

      !! B=P^(T)*A
      CALL cp_fm_gemm('T','N',nao_orb,nao_orb,nao_aux_fit,&
                      1.0_dp,wfn_fitting_env%P,wfn_fitting_env%A,0.0_dp,&
                      wfn_fitting_env%B,error)
   END IF

    ! *** Calculate the mo_coeffs for the fitting basis
    DO ispin=1,nspins
      nmo = wfn_fitting_env%nmo(ispin)
      !! Lambda = C^(T)*B*C
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
      CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)
       
      CALL cp_fm_gemm('N','N',nao_orb,nmo,nao_orb,&
                      1.0_dp,wfn_fitting_env%B,mo_coeff,0.0_dp,&
                      wfn_fitting_env%work_orb_nmo(ispin)%matrix,error)
      CALL cp_fm_gemm('T','N',nmo,nmo,nao_orb,&
                      1.0_dp,mo_coeff,wfn_fitting_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                      wfn_fitting_env%lambda(ispin)%matrix,error)
      CALL cp_fm_to_fm(wfn_fitting_env%lambda(ispin)%matrix, wfn_fitting_env%work_nmo_nmo1(ispin)%matrix, error=error) 
      CALL cp_fm_syevd(wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,wfn_fitting_env%R(ispin)%matrix,&
                       wfn_fitting_env%eigvals_lambda(ispin)%eigvals%data,error)
      ALLOCATE(eig_work(nmo), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO i=1,nmo
        eig_work(i) = 1.0_dp/SQRT(wfn_fitting_env%eigvals_lambda(ispin)%eigvals%data(i))
      END DO
      CALL cp_fm_to_fm(wfn_fitting_env%R(ispin)%matrix, wfn_fitting_env%work_nmo_nmo1(ispin)%matrix, error=error)
      CALL cp_fm_column_scale(wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,eig_work)       
      CALL cp_fm_gemm('N','T',nmo,nmo,nmo,&
                      1.0_dp,wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,wfn_fitting_env%R(ispin)%matrix,0.0_dp,&
                      wfn_fitting_env%lambda_inv_sqrt(ispin)%matrix,error)
      CALL cp_fm_gemm('N','N',nao_orb,nmo,nmo,&
                      1.0_dp,mo_coeff,wfn_fitting_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                      wfn_fitting_env%work_orb_nmo(ispin)%matrix,error)
      CALL cp_fm_gemm('N','N',nao_aux_fit,nmo,nao_orb,&
                      1.0_dp,wfn_fitting_env%A,wfn_fitting_env%work_orb_nmo(ispin)%matrix, 0.0_dp,&
                      mo_coeff_aux_fit,error)
      DEALLOCATE(eig_work, STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END DO
  END SUBROUTINE fit_mo_coeffs


! *****************************************************************************
!> \brief Calculates the product Kohn-Sham-Matrix x mo_coeff for the auxiliary 
!>        basis set and transforms it into the orbital basis. This is needed
!>        in order to use OT
!>      
!> \param ispin which spin to transform
!> \param wfn_fitting_env The wavefunction fitting env
!> \param mo_coeff the MO coefficients from the orbital basis set  
!> \param mo_coeff_aux_fit the MO coefficients from the auxiliary fitting basis set
!> \param mo_derivs KS x mo_coeff from the orbital basis set to which we add the
!>        auxiliary basis set part
!> \param matrix_ks_aux_fit the Kohn-Sham matrix from the auxiliary fitting basis set
!> \param error 
!>
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE merge_ks_times_mo_coeff(ispin, wfn_fitting_env, mo_coeff, mo_coeff_aux_fit, mo_derivs, &
                                     matrix_ks_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(wfn_fitting_type), POINTER          :: wfn_fitting_env
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: i, j, nao_aux_fit, nao_orb, &
                                                nmo
    REAL(dp)                                 :: eig_diff, pole
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eigvals_lambda

    nao_aux_fit = wfn_fitting_env%nao_aux_fit
    nao_orb = wfn_fitting_env%nao_orb
    nmo = wfn_fitting_env%nmo(ispin)

    CALL copy_sm_to_fm(matrix_ks_aux_fit(ispin)%matrix,wfn_fitting_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(wfn_fitting_env%K(ispin)%matrix,wfn_fitting_env%work_aux_aux,error=error)

    CALL cp_fm_gemm('N','N', nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,wfn_fitting_env%K(ispin)%matrix,mo_coeff_aux_fit,0.0_dp,&
                    wfn_fitting_env%H(ispin)%matrix,error)

    ! *** Add first term 
    CALL cp_fm_gemm('N','T', nao_aux_fit, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%H(ispin)%matrix,wfn_fitting_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_aux_nmo(ispin)%matrix,error)
    CALL cp_fm_gemm('T','N', nao_orb, nmo, nao_aux_fit,&
                    1.0_dp,wfn_fitting_env%A,wfn_fitting_env%work_aux_nmo(ispin)%matrix,1.0_dp,&
                    mo_derivs(ispin)%matrix,error)


    !! Construct Matrix M for Hadamard Product
    pole = 0.0_dp
    DO i=1,nmo
      DO j=i,nmo
        eig_diff = ABS( wfn_fitting_env%eigvals_lambda(ispin)%eigvals%data(i) -&
                        wfn_fitting_env%eigvals_lambda(ispin)%eigvals%data(j) )
        ! *** two eigenvalues could be the degenerated. In that case use 2nd order formula for the poles
        IF( i==j .OR. ( eig_diff == 0.0_dp ) ) THEN
          pole = -0.5_dp/SQRT(wfn_fitting_env%eigvals_lambda(ispin)%eigvals%data(i))**3
          CALL cp_fm_set_element(wfn_fitting_env%M(ispin)%matrix,i,j,pole,error)
        ELSE
          pole = 1.0_dp/SQRT(wfn_fitting_env%eigvals_lambda(ispin)%eigvals%data(i))
          pole = pole - 1.0_dp/SQRT(wfn_fitting_env%eigvals_lambda(ispin)%eigvals%data(j))
          pole = pole/(wfn_fitting_env%eigvals_lambda(ispin)%eigvals%data(i)-&
                       wfn_fitting_env%eigvals_lambda(ispin)%eigvals%data(j))
          CALL cp_fm_set_element(wfn_fitting_env%M(ispin)%matrix,i,j,pole,error) 
        END IF
      END DO
    END DO
    CALL cp_fm_upper_to_full(wfn_fitting_env%M(ispin)%matrix,wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,error=error)    

    !! 2nd term to be added to fm_H
  
    !! Part 1: B^(T)*C* R*[R^(T)*c^(T)*A^(T)*H_aux_fit*R x M]*R^(T)
    !! Part 2: B*C*(R*[R^(T)*c^(T)*A^(T)*H_aux_fit*R x M]*R^(T))^(T)

    !! H'*R
    CALL cp_fm_gemm('N','N', nao_aux_fit, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%H(ispin)%matrix,wfn_fitting_env%R(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_aux_nmo(ispin)%matrix,error)
    !! A^(T)*H'*R
    CALL cp_fm_gemm('T','N', nao_orb, nmo, nao_aux_fit,&
                    1.0_dp,wfn_fitting_env%A,wfn_fitting_env%work_aux_nmo(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_orb_nmo(ispin)%matrix,error)
    !! c^(T)*A^(T)*H'*R
    CALL cp_fm_gemm('T','N', nmo, nmo, nao_orb,&
                    1.0_dp,mo_coeff,wfn_fitting_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,error)
    !! R^(T)*c^(T)*A^(T)*H'*R
    CALL cp_fm_gemm('T','N', nmo, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%R(ispin)%matrix,wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_nmo_nmo2(ispin)%matrix,error)
    !! R^(T)*c^(T)*A^(T)*H'*R x M
    CALL cp_fm_schur_product(wfn_fitting_env%work_nmo_nmo2(ispin)%matrix,&
                             wfn_fitting_env%M(ispin)%matrix,wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,error) 
    !! R* (R^(T)*c^(T)*A^(T)*H'*R x M)
    CALL cp_fm_gemm('N','N', nmo, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%R(ispin)%matrix,wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_nmo_nmo2(ispin)%matrix,error)

    !! R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)
    CALL cp_fm_gemm('N','T', nmo, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%work_nmo_nmo2(ispin)%matrix,wfn_fitting_env%R(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,error)

    !! B^(T)*c
    CALL cp_fm_gemm('T','N', nao_orb, nmo, nao_orb,&
                    1.0_dp,wfn_fitting_env%B,mo_coeff,0.0_dp,&
                    wfn_fitting_env%work_orb_nmo(ispin)%matrix,error)

    !! Add first term to fm_H
    !! B^(T)*c* R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)
    CALL cp_fm_gemm('N','N', nao_orb, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%work_orb_nmo(ispin)%matrix,wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,1.0_dp,&
                    mo_derivs(ispin)%matrix,error)

    !! Add second term to fm_H
    !! B*C *[ R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)]^(T)
    CALL cp_fm_gemm('N','T', nao_orb, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%work_orb_nmo(ispin)%matrix,wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,1.0_dp,&
                    mo_derivs(ispin)%matrix,error)

  END SUBROUTINE merge_ks_times_mo_coeff 

END MODULE wfn_fitting_methods
