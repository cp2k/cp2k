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
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
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
  USE cp_sm_fm_interactions,           ONLY: copy_fm_to_sm,&
                                             copy_sm_to_fm
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE input_constants,                 ONLY: use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE kinds,                           ONLY: dp
  USE mathconstants
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_set_type, &
       neighbor_list_type, neighbor_node_type, next
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             deallocate_matrix,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE wfn_fitting_types,               ONLY: wfn_fitting_env_create,&
                                             wfn_fitting_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC fit_mo_coeffs, &
         merge_mo_derivs, &
         calc_mixed_overlap_force

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

    INTEGER                                  :: handle, i, ispin, istat, &
                                                nao_aux_fit, nao_orb, nmo, &
                                                nspins
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eig_work
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit

    CALL timeset(routineN,handle)

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
    
    CALL timestop(handle)

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
  SUBROUTINE merge_mo_derivs(ispin, wfn_fitting_env, mo_coeff, mo_coeff_aux_fit, mo_derivs, &
                             mo_derivs_aux_fit, matrix_ks_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(wfn_fitting_type), POINTER          :: wfn_fitting_env
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs, mo_derivs_aux_fit
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_mo_derivs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, nao_aux_fit, &
                                                nao_orb, nmo
    REAL(dp)                                 :: eig_diff, pole
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eigvals_lambda

    CALL timeset(routineN,handle)

    nao_aux_fit = wfn_fitting_env%nao_aux_fit
    nao_orb = wfn_fitting_env%nao_orb
    nmo = wfn_fitting_env%nmo(ispin)

    CALL copy_sm_to_fm(matrix_ks_aux_fit(ispin)%matrix,wfn_fitting_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(wfn_fitting_env%K(ispin)%matrix,wfn_fitting_env%work_aux_aux,error=error)

    CALL cp_fm_gemm('N','N', nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,wfn_fitting_env%K(ispin)%matrix,mo_coeff_aux_fit,0.0_dp,&
                    wfn_fitting_env%H(ispin)%matrix,error)

    CALL cp_fm_to_fm(wfn_fitting_env%H(ispin)%matrix, mo_derivs_aux_fit(ispin)%matrix, error=error) 

    ! *** Add first term 
    CALL cp_fm_gemm('N','T', nao_aux_fit, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%H(ispin)%matrix,wfn_fitting_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_aux_nmo(ispin)%matrix,error)
    CALL cp_fm_gemm('T','N', nao_orb, nmo, nao_aux_fit,&
                    1.0_dp,wfn_fitting_env%A,wfn_fitting_env%work_aux_nmo(ispin)%matrix,1.0_dp,&
                    mo_derivs(ispin)%matrix,error)


    ! *** Construct Matrix M for Hadamard Product
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

    ! *** 2nd term to be added to fm_H
  
    !! Part 1: B^(T)*C* R*[R^(T)*c^(T)*A^(T)*H_aux_fit*R x M]*R^(T)
    !! Part 2: B*C*(R*[R^(T)*c^(T)*A^(T)*H_aux_fit*R x M]*R^(T))^(T)

    ! *** H'*R
    CALL cp_fm_gemm('N','N', nao_aux_fit, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%H(ispin)%matrix,wfn_fitting_env%R(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_aux_nmo(ispin)%matrix,error)
    ! *** A^(T)*H'*R
    CALL cp_fm_gemm('T','N', nao_orb, nmo, nao_aux_fit,&
                    1.0_dp,wfn_fitting_env%A,wfn_fitting_env%work_aux_nmo(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_orb_nmo(ispin)%matrix,error)
    ! *** c^(T)*A^(T)*H'*R
    CALL cp_fm_gemm('T','N', nmo, nmo, nao_orb,&
                    1.0_dp,mo_coeff,wfn_fitting_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,error)
    ! *** R^(T)*c^(T)*A^(T)*H'*R
    CALL cp_fm_gemm('T','N', nmo, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%R(ispin)%matrix,wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_nmo_nmo2(ispin)%matrix,error)
    ! *** R^(T)*c^(T)*A^(T)*H'*R x M
    CALL cp_fm_schur_product(wfn_fitting_env%work_nmo_nmo2(ispin)%matrix,&
                             wfn_fitting_env%M(ispin)%matrix,wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,error) 
    ! *** R* (R^(T)*c^(T)*A^(T)*H'*R x M)
    CALL cp_fm_gemm('N','N', nmo, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%R(ispin)%matrix,wfn_fitting_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_nmo_nmo2(ispin)%matrix,error)

    ! *** R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)
    CALL cp_fm_gemm('N','T', nmo, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%work_nmo_nmo2(ispin)%matrix,wfn_fitting_env%R(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%R_schur_R_t(ispin)%matrix,error)

    ! *** B^(T)*c
    CALL cp_fm_gemm('T','N', nao_orb, nmo, nao_orb,&
                    1.0_dp,wfn_fitting_env%B,mo_coeff,0.0_dp,&
                    wfn_fitting_env%work_orb_nmo(ispin)%matrix,error)

    ! *** Add first term to fm_H
    ! *** B^(T)*c* R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)
    CALL cp_fm_gemm('N','N', nao_orb, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%work_orb_nmo(ispin)%matrix,wfn_fitting_env%R_schur_R_t(ispin)%matrix,1.0_dp,&
                    mo_derivs(ispin)%matrix,error)

    ! *** Add second term to fm_H
    ! *** B*C *[ R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)]^(T)
    CALL cp_fm_gemm('N','T', nao_orb, nmo, nmo,&
                    1.0_dp,wfn_fitting_env%work_orb_nmo(ispin)%matrix,wfn_fitting_env%R_schur_R_t(ispin)%matrix,1.0_dp,&
                    mo_derivs(ispin)%matrix,error)

    CALL timestop(handle)

  END SUBROUTINE merge_mo_derivs 

! *****************************************************************************
!> \brief Calculates contribution of forces due to basis transformation
!> 
!>        dE/dR = dE/dC'*dC'/dR
!>        dE/dC = Ks'*c'*occ = H'
!>
!>        dC'/dR = - tr(A*lambda^(-1/2)*H'^(T)*S^(-1) * dS'/dR)
!>                 - tr(A*C*Y^(T)*C^(T)*P^(T)*A^(T) * dS'/dR)
!>                 + tr(C*lambda^(-1/2)*H'^(T)*S^(-1) * dP/dR)
!>                 + tr(A*C*Y^(T)*c^(T) * dP/dR)
!>                 + tr(C*Y^(T)*C^(T)*A^(T) * dP/dR)
!>
!>        where
!>
!>        A = S'^(-1)*P
!>        lambda = C^(T)*B*C
!>        B = P^(T)*A
!>        Y = R*[ (R^(T)*C^(T)*A^(T)*H'*R) xx M ]*R^(T)
!>        lambda = R*D*R^(T)
!>        Mij = Poles-Matrix (see above)
!>        xx = schur product
!>        
!> \param qs_env the QS environment      
!> \param ispin which spin to transform
!> \param wfn_fitting_env The wavefunction fitting env
!> \param mo_coeff the MO coefficients from the orbital basis set  
!> \param matrix_s_aux_fit overlap matrix from auxiliary fitting basis
!> \param matrix_s_aux_fit_vs_orb mixed orbital/aux fit overlap matrix 
!> \param error 
!>
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE calc_mixed_overlap_force(qs_env, para_env, ispin, wfn_fitting_env, mo_coeff,  &
                                      matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(wfn_fitting_type), POINTER          :: wfn_fitting_env
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_mixed_overlap_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nao_orb, &
                                                nblock, neighbor_list_id, &
                                                nmo, nrow
    INTEGER, DIMENSION(:), POINTER           :: first_col, first_row, &
                                                last_col, last_row
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(real_matrix_type), POINTER          :: matrix_w_p, matrix_w_s

    CALL timeset(routineN,handle)

    NULLIFY(matrix_w_p, matrix_w_s, distribution_2d)

    nao_aux_fit = wfn_fitting_env%nao_aux_fit
    nao_orb = wfn_fitting_env%nao_orb
    nmo = wfn_fitting_env%nmo(ispin)


    ! *** Create sparse work matrices
    CALL get_qs_env(qs_env=qs_env,&
                    distribution_2d=distribution_2d,&
                    neighbor_list_id=neighbor_list_id,&
                    error=error) 

    CALL get_matrix_info(matrix=matrix_s_aux_fit(1)%matrix,&
                         nblock_row=nblock,&
                         nrow=nrow,&
                         first_row=first_row,&
                         last_row=last_row,&
                         first_col=first_col,&
                         last_col=last_col)

    CALL allocate_matrix(matrix=matrix_w_s,&
                         nblock_row=nblock,&
                         nblock_col=nblock,&
                         nrow=nrow,&
                         ncol=nrow,&
                         first_row=first_row,&
                         last_row=last_row,&
                         first_col=first_col,&
                         last_col=last_col,&
                         matrix_name="W MATRIX AUX S",&
                         matrix_symmetry="none",&
                         sparsity_id=neighbor_list_id, &
                         distribution_2d=distribution_2d,&
                         error=error)
    ! *** allocate all blocks
    CALL allocate_w_matrix_blocks(qs_env, &
                                  matrix_w_s, &
                                  qs_env%sab_aux_fit_asymm, &
                                  do_symmetric=.FALSE., &
                                  error=error)

    CALL replicate_matrix_structure(matrix_s_aux_fit_vs_orb(1)%matrix,&
                                    matrix_w_p,&
                                    "W MATRIX AUX P",error=error)

    ! *** S'^(-T)*H'
    CALL cp_fm_gemm('T', 'N', nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,wfn_fitting_env%S,qs_env%mo_derivs_aux_fit(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_aux_nmo(ispin)%matrix,error)
    ! *** S'^(-T)*H'*Lambda^(-T/2)
    CALL cp_fm_gemm('N', 'T',  nao_aux_fit,nmo, nmo,&
                    1.0_dp,wfn_fitting_env%work_aux_nmo(ispin)%matrix,wfn_fitting_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_aux_nmo2(ispin)%matrix,error)
    ! *** C*Lambda^(-1/2)*H'^(T)*S'^(-1) minus sign due to force = -dE/dR
    CALL cp_fm_gemm('N', 'T',  nao_aux_fit, nao_orb, nmo,&
                    -1.0_dp,wfn_fitting_env%work_aux_nmo2(ispin)%matrix,mo_coeff,0.0_dp,&
                    wfn_fitting_env%work_aux_orb,error)

    ! *** A*C*Lambda^(-1/2)*H'^(T)*S'^(-1), minus sign to recover from above
    CALL cp_fm_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nao_orb,&
                    -1.0_dp,wfn_fitting_env%work_aux_orb,wfn_fitting_env%A,0.0_dp,&
                    wfn_fitting_env%work_aux_aux,error)

    ! *** C*Y
    CALL cp_fm_gemm('N', 'N',  nao_orb, nmo, nmo,&
                    1.0_dp,mo_coeff,wfn_fitting_env%R_schur_R_t(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_orb_nmo(ispin)%matrix,error)
    ! *** C*Y^(T)*C^(T)
    CALL cp_fm_gemm('N', 'T',  nao_orb, nao_orb, nmo,&
                    1.0_dp,mo_coeff,wfn_fitting_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_orb_orb,error)
    ! *** A*C*Y^(T)*C^(T) Add to work aux_orb, minus sign due to force = -dE/dR
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                    -1.0_dp,wfn_fitting_env%A,wfn_fitting_env%work_orb_orb,1.0_dp,&
                    wfn_fitting_env%work_aux_orb,error)

    ! *** C*Y^(T)
    CALL cp_fm_gemm('N', 'T',  nao_orb, nmo, nmo,&
                    1.0_dp,mo_coeff,wfn_fitting_env%R_schur_R_t(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_orb_nmo(ispin)%matrix,error)
    ! *** C*Y*C^(T)
    CALL cp_fm_gemm('N', 'T',  nao_orb, nao_orb, nmo,&
                    1.0_dp,mo_coeff,wfn_fitting_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                    wfn_fitting_env%work_orb_orb,error)
    ! *** A*C*Y*C^(T) Add to work aux_orb, minus sign due to -dE/dR
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                    -1.0_dp,wfn_fitting_env%A,wfn_fitting_env%work_orb_orb,1.0_dp,&
                    wfn_fitting_env%work_aux_orb,error)

    ! *** copy to sparse matrix
    CALL copy_fm_to_sm(wfn_fitting_env%work_aux_orb, matrix_w_p)
 
    ! *** A*C*Y^(T)*C^(T)
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                    1.0_dp,wfn_fitting_env%A,wfn_fitting_env%work_orb_orb,0.0_dp,&
                    wfn_fitting_env%work_aux_orb,error)
    ! *** A*C*Y^(T)*C^(T)*A^(T) add to aux_aux, minus sign cancels
    CALL cp_fm_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nao_orb,&
                    1.0_dp,wfn_fitting_env%work_aux_orb,wfn_fitting_env%A,1.0_dp,&
                    wfn_fitting_env%work_aux_aux,error)

    ! *** copy to sparse matrix
    CALL copy_fm_to_sm(wfn_fitting_env%work_aux_aux, matrix_w_s)

    ! *** This can be done in one call w_total = w_alpha + w_beta

    CALL calc_overlap_force_general(qs_env,para_env,matrix_w_s,&
                                   use_aux_fit_basis_set, use_aux_fit_basis_set, qs_env%sab_aux_fit_asymm, &
                                   do_symmetric=.FALSE., error=error)

    CALL calc_overlap_force_general(qs_env,para_env,matrix_w_p,&
                                    use_aux_fit_basis_set, use_orb_basis_set, qs_env%sab_aux_fit_vs_orb, &
                                    do_symmetric=.FALSE., error=error)
   
    ! *** Deallocated weighted density matrices
    CALL deallocate_matrix(matrix_w_s,error)
    CALL deallocate_matrix(matrix_w_p,error)

    CALL timestop(handle)

  END SUBROUTINE calc_mixed_overlap_force

! *****************************************************************************
!> \brief Calculates dExc/dR from contribution of wfn_fitting and adds the
!>        forces to force(:)%overlap_wfn_fit(:,:). The routines basically
!>        calculates dE/dR= trace(W^(T)*d(S)/dR) using a decontraction of
!>        W and S in order to assign the corresponding contributions to the
!>        correct kinds
!>      
!> \param qs_env the qs environment 
!> \param para_env the parallel environment
!> \param matrix_w the energy weighted density matrix
!> \param basis_set_id_a basis set for bra decontraction
!> \param basis_set_id_b basis set for ket decontraction
!> \param sab_orb the corresponding neighbor list
!> \param do_symmetric flag to indicate if the involved matrices are symmetric
!>        or not. We only have symmetry if basis_set_id_a == basis_set_id_b AND
!>        W=W^(T). Furthermore, the neighbor list symmetry must be consistent
!>        with this flag!
!> \param error 
!>
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE calc_overlap_force_general(qs_env,para_env,matrix_w,&
                                        basis_set_id_a, basis_set_id_b, sab_orb, &
                                        do_symmetric, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(real_matrix_type), POINTER          :: matrix_w
    INTEGER, INTENT(IN)                      :: basis_set_id_a, basis_set_id_b
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    LOGICAL, INTENT(IN)                      :: do_symmetric
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_overlap_force_general', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iab, iatom, icol, ikind, ilist, inode, &
      irow, iset, jatom, jkind, jset, last_jatom, ldai, ldsab, maxblock, &
      maxco, maxcoa, maxcob, maxdco, maxdcoa, maxdcob, maxder, maxl, maxlgto, &
      maxlgtoa, maxlgtob, maxlppl, maxlppla, maxlpplb, maxsgf, maxsgfa, &
      maxsgfb, n, natom, ncoa, ncoa_sum, ncob, nder, neighbor_list_id, nkind, &
      nlist, nnode, nrow, nseta, nsetb, sgfa, sgfb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, first_sgf
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, new_atom_b, &
                                                return_s_derivatives
    REAL(KIND=dp)                            :: dab, rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: ai_work
    REAL(KIND=dp), DIMENSION(3)              :: force_a, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab, rpgfa, rpgfb, sab, &
                                                sphi_a, sphi_b, w_block, &
                                                work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: sdab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set_a, &
                                                orb_basis_set_b
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

    failure = .FALSE.
    CALL timeset(routineN,handle)


    NULLIFY (atomic_kind_set)
    NULLIFY (force)
    NULLIFY (pab)
    NULLIFY (particle_set)
    NULLIFY (sab)
    NULLIFY (sdab)
    NULLIFY (work)


    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    neighbor_list_id=neighbor_list_id,&
                    force=force,&
                    error=error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)


    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

    nder = 1
    maxder = ncoset(nder)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
    
    ! *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxcoa,&
                             maxlgto=maxlgtoa,&
                             maxlppl=maxlppla,&
                             maxsgf=maxsgfa,&
                             basis_set_id=basis_set_id_a)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxcob,&
                             maxlgto=maxlgtob,&
                             maxlppl=maxlpplb,&
                             maxsgf=maxsgfb,&
                             basis_set_id=basis_set_id_b)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxdcoa,&
                             maxder=MAX(1,nder),&
                             basis_set_id=basis_set_id_a)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxdcob,&
                             maxder=MAX(1,nder),&
                             basis_set_id=basis_set_id_b)


    maxdco = MAX(maxdcoa,maxdcob)
    maxco = MAX(maxcoa,maxcob)
    maxlgto = MAX(maxlgtoa,maxlgtob)
    maxlppl = MAX(maxlppla,maxlpplb)
    maxsgf = MAX(maxsgfa,maxsgfb)
    maxl = maxlgto
    maxder = nder
    CALL init_orbital_pointers(maxl+nder+1)

    ldsab = MAX(maxco,maxsgf)
    maxblock = MAX(4,maxder) ! 4=size(matrix_s_aux_fit)
    ldai = ncoset(maxl+nder+1)

    ALLOCATE(sab(ldsab,ldsab*maxblock))
    ALLOCATE(sdab(maxdco,maxco,4))
    ALLOCATE(work(ldsab,ldsab*maxder))
    ALLOCATE(pab(maxco,maxco))


    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      SELECT CASE (basis_set_id_a)
      CASE (use_orb_basis_set)
        CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set_a)
      CASE (use_aux_fit_basis_set)
        CALL get_atomic_kind(atomic_kind=atomic_kind,aux_fit_basis_set=orb_basis_set_a)
      END SELECT

      IF (.NOT.ASSOCIATED(orb_basis_set_a)) CYCLE
      CALL get_gto_basis_set(gto_basis_set=orb_basis_set_a,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             nco_sum=ncoa_sum,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             pgf_radius=rpgfa,&
                             set_radius=set_radius_a,&
                             sphi=sphi_a,&
                             zet=zeta)
      DO jkind=1,nkind
        atomic_kind => atomic_kind_set(jkind)
        SELECT CASE (basis_set_id_b)
        CASE (use_orb_basis_set)
          CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set_b)
        CASE (use_aux_fit_basis_set)
          CALL get_atomic_kind(atomic_kind=atomic_kind,aux_fit_basis_set=orb_basis_set_b)
        END SELECT


        IF (.NOT.ASSOCIATED(orb_basis_set_b)) CYCLE
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set_b,&
                               first_sgf=first_sgfb,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               nsgf_set=nsgfb,&
                               pgf_radius=rpgfb,&
                               set_radius=set_radius_b,&
                               sphi=sphi_b,&
                               zet=zetb)

        iab = ikind + nkind*(jkind - 1)

        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE
        CALL get_neighbor_list_set(neighbor_list_set=sab_orb(iab)%neighbor_list_set,nlist=nlist)

        NULLIFY ( sab_orb_neighbor_list )
        ALLOCATE (ai_work(ldai,ldai,MAX(1,ncoset(maxlppl)),ncoset(nder+1)),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ilist=1,nlist
          IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
             sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
          ELSE
             sab_orb_neighbor_list => next(sab_orb_neighbor_list)
          END IF
          sab_orb_neighbor_list_local => sab_orb_neighbor_list
          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,atom=iatom,nnode=nnode)

          atom_a = atom_of_kind(iatom)

          last_jatom = 0
          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode
            CALL get_neighbor_node(sab_orb_neighbor_node,neighbor=jatom,r=rab)

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
            END IF

            atom_b = atom_of_kind(jatom)

            IF (new_atom_b) THEN
              IF ( do_symmetric ) THEN
                IF (iatom <= jatom) THEN
                   irow = iatom
                   icol = jatom
                ELSE
                   irow = jatom
                   icol = iatom
                END IF
              ELSE
                irow = iatom
                icol = jatom
              END IF
              CALL get_block_node(matrix=matrix_w,&
                                  block_row=irow,&
                                  block_col=icol,&
                                  BLOCK=w_block)

              IF (.NOT.ASSOCIATED(w_block)) THEN
                 CALL stop_program(routineN,moduleN,__LINE__,&
                      "A weighted density matrix block is missing",para_env)
              END IF
            END IF

            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = SQRT(rab2)
            nrow = 1
            DO iset=1,nseta
              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)
              DO jset=1,nsetb
                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)
                IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN
                  IF ( iatom /= jatom ) THEN
                    ! *** Decontract W matrix block ***
                    IF( do_symmetric ) THEN
                      IF (iatom <= jatom) THEN
                        CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                    w_block(sgfa,sgfb),SIZE(w_block,1),&
                                    0.0_dp,work(1,1),SIZE(work,1))
                      ELSE
                        CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                                   1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                   w_block(sgfb,sgfa),SIZE(w_block,1),&
                                   0.0_dp,work(1,1),SIZE(work,1))
                      END IF
                    ELSE
                      CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 w_block(sgfa,sgfb),SIZE(w_block,1),&
                                 0.0_dp,work(1,1),SIZE(work,1))
                    END IF  
                    CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                               1.0_dp,work(1,1),SIZE(work,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_dp,pab(1,1),SIZE(pab,1))

                    ! *** Calculate the primitive overlap integrals ***
                    ! *** and the corresponding force contribution  ***
                    CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,sab,nder,return_s_derivatives,&
                                 ai_work,ldai,sdab,pab,force_a)
                    IF( do_symmetric ) THEN
                      force(ikind)%overlap_wfn_fit(:,atom_a)=force(ikind)%overlap_wfn_fit(:,atom_a) - 2.0_dp*force_a(:)
                      force(jkind)%overlap_wfn_fit(:,atom_b)=force(jkind)%overlap_wfn_fit(:,atom_b) + 2.0_dp*force_a(:)
                    ELSE
                      force(ikind)%overlap_wfn_fit(:,atom_a)=force(ikind)%overlap_wfn_fit(:,atom_a) - force_a(:)
                      force(jkind)%overlap_wfn_fit(:,atom_b)=force(jkind)%overlap_wfn_fit(:,atom_b) + force_a(:)
                    END IF
                  END IF
                END IF
              END DO
              nrow = nrow + ncoa
            END DO
            sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          END DO
        END DO
        DEALLOCATE (ai_work,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END DO
    END DO

    DEALLOCATE (sab,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (sdab,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (atom_of_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (pab,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE calc_overlap_force_general

! *****************************************************************************
!> \brief Allocates all blocks of a distributed weighted density matrix W wrt
!>        symmetry and neighbor list
!>      
!> \param qs_env the qs environment 
!> \param matrix_w the energy weighted density matrix
!> \param neighbor_list_sab the corresponding neighbor list
!> \param do_symmetric flag to indicate if W and the neighbor list are symmetric
!>        both have to be consistent!
!> \param error 
!>
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE allocate_w_matrix_blocks(qs_env, matrix_w, neighbor_list_sab, do_symmetric, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_type), POINTER          :: matrix_w
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: neighbor_list_sab
    LOGICAL, INTENT(IN)                      :: do_symmetric
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: iab, iatom, icol, ikind, &
                                                ilist, inode, irow, jatom, &
                                                jkind, last_jatom, nkind, &
                                                nlist, nnode
    LOGICAL                                  :: new_atom_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    TYPE(neighbor_node_type), POINTER        :: neighbor_node

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    error=error)

    nkind = SIZE(atomic_kind_set)

    
    DO ikind=1,nkind
      DO jkind=1,nkind
        iab = ikind + nkind*(jkind - 1)

        IF (.NOT.ASSOCIATED(neighbor_list_sab(iab)%neighbor_list_set)) CYCLE

        neighbor_list_set => neighbor_list_sab(iab)%neighbor_list_set

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   nlist=nlist)

        neighbor_list => first_list(neighbor_list_set)

        DO ilist=1,nlist

          CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                 atom=iatom,&
                                 nnode=nnode)

          last_jatom = 0

          neighbor_node => first_node(neighbor_list)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=neighbor_node,&
                                   neighbor=jatom)

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
            END IF

            IF (new_atom_b) THEN
              IF( do_symmetric ) THEN
                IF (iatom <= jatom) THEN
                  irow = iatom
                  icol = jatom
                ELSE
                  irow = jatom
                  icol = iatom
                END IF
              ELSE
                irow = iatom
                icol = jatom
              END IF
              NULLIFY (block)
              CALL add_block_node(matrix=matrix_w,&
                                  block_row=irow,&
                                  block_col=icol,&
                                  BLOCK=block,error=error)
            END IF
            neighbor_node => next(neighbor_node)

          END DO

          neighbor_list => next(neighbor_list)

        END DO

      END DO

    END DO

  END SUBROUTINE allocate_w_matrix_blocks

END MODULE wfn_fitting_methods
