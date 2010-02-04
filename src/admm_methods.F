!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains methods used in the context of density fitting 
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE admm_methods 
  USE admm_types,                      ONLY: admm_env_create,&
                                             admm_type
  USE ai_overlap_new,                  ONLY: overlap
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             dbcsr_alloc_block_from_nbl,&
                                             dbcsr_deallocate_matrix,&
                                             dbcsr_obj_p_type
  USE cp_dbcsr_output,                 ONLY: dbcsr_write_sparse_matrix
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_schur_product,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_methods,                   ONLY: dbcsr_col_block_sizes,&
                                             dbcsr_distribution,&
                                             dbcsr_init,&
                                             dbcsr_row_block_sizes,&
                                             dbcsr_uses_special_memory
  USE dbcsr_operations,                ONLY: dbcsr_add,&
                                             dbcsr_copy,&
                                             dbcsr_get_block_p,&
                                             dbcsr_set
  USE dbcsr_types,                     ONLY: dbcsr_obj,&
                                             dbcsr_type_normal,&
                                             dbcsr_type_symmetric
  USE dbcsr_util,                      ONLY: dbcsr_create,&
                                             dbcsr_finalize
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE input_constants,                 ONLY: do_hfx_potential_coulomb,&
                                             do_hfx_potential_short,&
                                             do_hfx_potential_truncated,&
                                             use_aux_fit_basis_set,&
                                             use_orb_basis_set,&
                                             xc_funct_no_shortcut
  USE input_section_types,             ONLY: section_vals_duplicate,&
                                             section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE mathconstants
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC fit_mo_coeffs, &
         merge_ks_matrix, &
         merge_mo_derivs, &
         remove_ks_matrix, &
         calc_mixed_overlap_force, &
         create_admm_xc_section,&
         print_matlab_matrix

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_methods'

!***
  
  CONTAINS

! *****************************************************************************
!> \brief Calculates the MO coefficients for the auxiliary fitting basis set  
!>        by minimizing int (psi_i - psi_aux_i)^2 using Lagrangian Multipliers
!>      
!> \param admm_env The ADMM env
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
  SUBROUTINE fit_mo_coeffs(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                           mos, mos_aux_fit, geometry_did_change, error) 

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
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
    TYPE(section_vals_type), POINTER         :: input, xc_section

    CALL timeset(routineN,handle)

    IF (.NOT.(ASSOCIATED(admm_env) )) THEN
      CALL admm_env_create(mos, mos_aux_fit, &
                                  para_env, admm_env,&
                                  error)
      CALL get_qs_env(qs_env, input=input, error=error)
      xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
      CALL create_admm_xc_section(qs_env, xc_section, admm_env, error)
    END IF

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)


    ! *** This part only depends on overlap matrices ==> needs only to be calculated if the geometry changed

    IF( geometry_did_change ) THEN
      CALL copy_dbcsr_to_fm(matrix_s_aux_fit(1)%matrix,admm_env%S_inv,error)
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      CALL cp_fm_to_fm(admm_env%S_inv, admm_env%S, error=error) 

      CALL copy_dbcsr_to_fm(matrix_s_mixed(1)%matrix,admm_env%Q,error)


      !! Calculate S'_inverse
      CALL cp_fm_cholesky_decompose(admm_env%S_inv,error=error)
      CALL cp_fm_cholesky_invert(admm_env%S_inv,error=error) 
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      !! Calculate A=S'^(-1)*P
      CALL cp_fm_gemm('N','N',nao_aux_fit,nao_orb,nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%Q,0.0_dp,&
                    admm_env%A,error)

      !! B=Q^(T)*A
      CALL cp_fm_gemm('T','N',nao_orb,nao_orb,nao_aux_fit,&
                      1.0_dp,admm_env%Q,admm_env%A,0.0_dp,&
                      admm_env%B,error)
 
   END IF

    ! *** Calculate the mo_coeffs for the fitting basis
    DO ispin=1,nspins
      nmo = admm_env%nmo(ispin)
      IF( nmo == 0 ) CYCLE
      !! Lambda = C^(T)*B*C
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
      CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)
       
      CALL cp_fm_gemm('N','N',nao_orb,nmo,nao_orb,&
                      1.0_dp,admm_env%B,mo_coeff,0.0_dp,&
                      admm_env%work_orb_nmo(ispin)%matrix,error)
      CALL cp_fm_gemm('T','N',nmo,nmo,nao_orb,&
                      1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                      admm_env%lambda(ispin)%matrix,error)
      CALL cp_fm_to_fm(admm_env%lambda(ispin)%matrix, admm_env%work_nmo_nmo1(ispin)%matrix, error=error) 
      CALL cp_fm_syevd(admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%R(ispin)%matrix,&
                       admm_env%eigvals_lambda(ispin)%eigvals%data,error)
      ALLOCATE(eig_work(nmo), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO i=1,nmo
        eig_work(i) = 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(i))
      END DO
      CALL cp_fm_to_fm(admm_env%R(ispin)%matrix, admm_env%work_nmo_nmo1(ispin)%matrix, error=error)
      CALL cp_fm_column_scale(admm_env%work_nmo_nmo1(ispin)%matrix,eig_work)       
      CALL cp_fm_gemm('N','T',nmo,nmo,nmo,&
                      1.0_dp,admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%R(ispin)%matrix,0.0_dp,&
                      admm_env%lambda_inv_sqrt(ispin)%matrix,error)
      CALL cp_fm_gemm('N','N',nao_orb,nmo,nmo,&
                      1.0_dp,mo_coeff,admm_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                      admm_env%work_orb_nmo(ispin)%matrix,error)
      CALL cp_fm_gemm('N','N',nao_aux_fit,nmo,nao_orb,&
                      1.0_dp,admm_env%A,admm_env%work_orb_nmo(ispin)%matrix, 0.0_dp,&
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
!> \param admm_env The ADMM env
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
  SUBROUTINE merge_mo_derivs(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, mo_derivs, &
                             mo_derivs_aux_fit, matrix_ks_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs, mo_derivs_aux_fit
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_mo_derivs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, nao_aux_fit, &
                                                nao_orb, nmo, stat
    INTEGER, SAVE                            :: counter = 0
    LOGICAL                                  :: failure
    REAL(dp)                                 :: eig_diff, pole, tmp32, tmp52, &
                                                tmp72, tmp92
    REAL(dp), DIMENSION(:), POINTER          :: occupation_numbers, &
                                                scaling_factor

    failure = .FALSE.

    CALL timeset(routineN,handle)

    counter = counter + 1

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)

    CALL cp_fm_gemm('N','N', nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,mo_coeff_aux_fit,0.0_dp,&
                    admm_env%H(ispin)%matrix,error)

    CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
    ALLOCATE(scaling_factor(SIZE(occupation_numbers)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    scaling_factor = 2.0_dp*occupation_numbers

    CALL cp_fm_column_scale(admm_env%H(ispin)%matrix,scaling_factor)

    CALL cp_fm_to_fm(admm_env%H(ispin)%matrix, mo_derivs_aux_fit(ispin)%matrix, error=error) 

    ! *** Add first term 
    CALL cp_fm_gemm('N','T', nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%H(ispin)%matrix,admm_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix,error)
    CALL cp_fm_gemm('T','N', nao_orb, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%A,admm_env%work_aux_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix,error)


    ! *** Construct Matrix M for Hadamard Product
    pole = 0.0_dp
    DO i=1,nmo
      DO j=i,nmo
        eig_diff = ( admm_env%eigvals_lambda(ispin)%eigvals%data(i) -&
                     admm_env%eigvals_lambda(ispin)%eigvals%data(j) )
        ! *** two eigenvalues could be the degenerated. In that case use 2nd order formula for the poles
        IF( ABS(eig_diff) < 0.0001_dp ) THEN
          tmp32 = 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(j))**3
          tmp52 = tmp32/admm_env%eigvals_lambda(ispin)%eigvals%data(j)*eig_diff
          tmp72 = tmp52/admm_env%eigvals_lambda(ispin)%eigvals%data(j)*eig_diff
          tmp92 = tmp72/admm_env%eigvals_lambda(ispin)%eigvals%data(j)*eig_diff

          pole = -0.5_dp*tmp32 + 3.0_dp/8.0_dp*tmp52 - 5.0_dp/16.0_dp*tmp72 + 35.0_dp/128.0_dp*tmp92
          CALL cp_fm_set_element(admm_env%M(ispin)%matrix,i,j,pole,error)
        ELSE
          pole = 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(i))
          pole = pole - 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(j))
          pole = pole/(admm_env%eigvals_lambda(ispin)%eigvals%data(i)-&
                       admm_env%eigvals_lambda(ispin)%eigvals%data(j))
          CALL cp_fm_set_element(admm_env%M(ispin)%matrix,i,j,pole,error) 
        END IF
      END DO
    END DO
    CALL cp_fm_upper_to_full(admm_env%M(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,error=error)    

    ! *** 2nd term to be added to fm_H
  
    !! Part 1: B^(T)*C* R*[R^(T)*c^(T)*A^(T)*H_aux_fit*R x M]*R^(T)
    !! Part 2: B*C*(R*[R^(T)*c^(T)*A^(T)*H_aux_fit*R x M]*R^(T))^(T)

    ! *** H'*R
    CALL cp_fm_gemm('N','N', nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%H(ispin)%matrix,admm_env%R(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix,error)
    ! *** A^(T)*H'*R
    CALL cp_fm_gemm('T','N', nao_orb, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%A,admm_env%work_aux_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_nmo(ispin)%matrix,error)
    ! *** c^(T)*A^(T)*H'*R
    CALL cp_fm_gemm('T','N', nmo, nmo, nao_orb,&
                    1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%work_nmo_nmo1(ispin)%matrix,error)
    ! *** R^(T)*c^(T)*A^(T)*H'*R
    CALL cp_fm_gemm('T','N', nmo, nmo, nmo,&
                    1.0_dp,admm_env%R(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                    admm_env%work_nmo_nmo2(ispin)%matrix,error)
    ! *** R^(T)*c^(T)*A^(T)*H'*R x M
    CALL cp_fm_schur_product(admm_env%work_nmo_nmo2(ispin)%matrix,&
                             admm_env%M(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,error) 
    ! *** R* (R^(T)*c^(T)*A^(T)*H'*R x M)
    CALL cp_fm_gemm('N','N', nmo, nmo, nmo,&
                    1.0_dp,admm_env%R(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                    admm_env%work_nmo_nmo2(ispin)%matrix,error)

    ! *** R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)
    CALL cp_fm_gemm('N','T', nmo, nmo, nmo,&
                    1.0_dp,admm_env%work_nmo_nmo2(ispin)%matrix,admm_env%R(ispin)%matrix,0.0_dp,&
                    admm_env%R_schur_R_t(ispin)%matrix,error)

    ! *** B^(T)*c
    CALL cp_fm_gemm('T','N', nao_orb, nmo, nao_orb,&
                    1.0_dp,admm_env%B,mo_coeff,0.0_dp,&
                    admm_env%work_orb_nmo(ispin)%matrix,error)

    ! *** Add first term to fm_H
    ! *** B^(T)*c* R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)
    CALL cp_fm_gemm('N','N', nao_orb, nmo, nmo,&
                    1.0_dp,admm_env%work_orb_nmo(ispin)%matrix,admm_env%R_schur_R_t(ispin)%matrix,1.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix,error)

    ! *** Add second term to fm_H
    ! *** B*C *[ R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)]^(T)
    CALL cp_fm_gemm('N','T', nao_orb, nmo, nmo,&
                    1.0_dp,admm_env%work_orb_nmo(ispin)%matrix,admm_env%R_schur_R_t(ispin)%matrix,1.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix,error)

    DO i = 1,SIZE(scaling_factor)
      scaling_factor(i) = 1.0_dp/scaling_factor(i)
    END DO

    CALL cp_fm_column_scale(admm_env%mo_derivs_tmp(ispin)%matrix,scaling_factor)

    CALL cp_fm_scale_and_add(1.0_dp,mo_derivs(ispin)%matrix,1.0_dp,admm_env%mo_derivs_tmp(ispin)%matrix,error)

    DEALLOCATE(scaling_factor, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE merge_mo_derivs 

! *****************************************************************************
!> \brief Calculates contribution of forces due to basis transformation
!> 
!>        dE/dR = dE/dC'*dC'/dR
!>        dE/dC = Ks'*c'*occ = H'
!>
!>        dC'/dR = - tr(A*lambda^(-1/2)*H'^(T)*S^(-1) * dS'/dR)
!>                 - tr(A*C*Y^(T)*C^(T)*Q^(T)*A^(T) * dS'/dR)
!>                 + tr(C*lambda^(-1/2)*H'^(T)*S^(-1) * dQ/dR)
!>                 + tr(A*C*Y^(T)*c^(T) * dQ/dR)
!>                 + tr(C*Y^(T)*C^(T)*A^(T) * dQ/dR)
!>
!>        where
!>
!>        A = S'^(-1)*Q
!>        lambda = C^(T)*B*C
!>        B = Q^(T)*A
!>        Y = R*[ (R^(T)*C^(T)*A^(T)*H'*R) xx M ]*R^(T)
!>        lambda = R*D*R^(T)
!>        Mij = Poles-Matrix (see above)
!>        xx = schur product
!>        
!> \param qs_env the QS environment      
!> \param ispin which spin to transform
!> \param admm_env The ADMM env
!> \param mo_coeff the MO coefficients from the orbital basis set  
!> \param matrix_s_aux_fit overlap matrix from auxiliary fitting basis
!> \param matrix_s_aux_fit_vs_orb mixed orbital/aux fit overlap matrix 
!> \param error 
!>
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE calc_mixed_overlap_force(qs_env, para_env, ispin, admm_env, mo_coeff,  &
                                      matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, logger,&
                                      iw, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER                                  :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_mixed_overlap_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nao_orb, &
                                                neighbor_list_id, nmo
    TYPE(dbcsr_obj), POINTER                 :: matrix_w_q, matrix_w_s
    TYPE(distribution_2d_type), POINTER      :: distribution_2d

    CALL timeset(routineN,handle)

    NULLIFY(matrix_w_q, matrix_w_s, distribution_2d)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)


    ! *** Create sparse work matrices
    CALL get_qs_env(qs_env=qs_env,&
                    distribution_2d=distribution_2d,&
                    neighbor_list_id=neighbor_list_id,&
                    error=error) 


    ALLOCATE(matrix_w_s)
    CALL dbcsr_init (matrix_w_s, error)
    CALL dbcsr_create(matrix_w_s, 'W MATRIX AUX S', &
         dbcsr_distribution(matrix_s_aux_fit(1)%matrix), dbcsr_type_normal, dbcsr_row_block_sizes(matrix_s_aux_fit(1)%matrix),&
         dbcsr_col_block_sizes(matrix_s_aux_fit(1)%matrix), matrix_s_aux_fit(1)%matrix%m%nblks, matrix_s_aux_fit(1)%matrix%m%nze,&
         matrix_s_aux_fit(1)%matrix%m%data_type, dbcsr_uses_special_memory(matrix_s_aux_fit(1)%matrix), error=error)
    CALL dbcsr_finalize(matrix_w_s, error=error)
    CALL dbcsr_alloc_block_from_nbl(matrix_w_s,qs_env%sab_aux_fit_asymm,error=error)

    ALLOCATE(matrix_w_q)
    CALL dbcsr_init(matrix_w_q)
    CALL dbcsr_copy(matrix_w_q,matrix_s_aux_fit_vs_orb(1)%matrix,&
                    "W MATRIX AUX Q",error=error)

    ! *** S'^(-T)*H'
    CALL cp_fm_gemm('T', 'N', nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,qs_env%mo_derivs_aux_fit(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix,error)
    ! *** S'^(-T)*H'*Lambda^(-T/2)
    CALL cp_fm_gemm('N', 'T',  nao_aux_fit,nmo, nmo,&
                    1.0_dp,admm_env%work_aux_nmo(ispin)%matrix,admm_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo2(ispin)%matrix,error)
    ! *** C*Lambda^(-1/2)*H'^(T)*S'^(-1) minus sign due to force = -dE/dR
    CALL cp_fm_gemm('N', 'T',  nao_aux_fit, nao_orb, nmo,&
                    -1.0_dp,admm_env%work_aux_nmo2(ispin)%matrix,mo_coeff,0.0_dp,&
                    admm_env%work_aux_orb,error)

    ! *** A*C*Lambda^(-1/2)*H'^(T)*S'^(-1), minus sign to recover from above
    CALL cp_fm_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nao_orb,&
                    -1.0_dp,admm_env%work_aux_orb,admm_env%A,0.0_dp,&
                    admm_env%work_aux_aux,error)

    ! *** C*Y
    CALL cp_fm_gemm('N', 'N',  nao_orb, nmo, nmo,&
                    1.0_dp,mo_coeff,admm_env%R_schur_R_t(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_nmo(ispin)%matrix,error)
    ! *** C*Y^(T)*C^(T)
    CALL cp_fm_gemm('N', 'T',  nao_orb, nao_orb, nmo,&
                    1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_orb,error)
    ! *** A*C*Y^(T)*C^(T) Add to work aux_orb, minus sign due to force = -dE/dR
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                    -1.0_dp,admm_env%A,admm_env%work_orb_orb,1.0_dp,&
                    admm_env%work_aux_orb,error)

    ! *** C*Y^(T)
    CALL cp_fm_gemm('N', 'T',  nao_orb, nmo, nmo,&
                    1.0_dp,mo_coeff,admm_env%R_schur_R_t(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_nmo(ispin)%matrix,error)
    ! *** C*Y*C^(T)
    CALL cp_fm_gemm('N', 'T',  nao_orb, nao_orb, nmo,&
                    1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_orb,error)
    ! *** A*C*Y*C^(T) Add to work aux_orb, minus sign due to -dE/dR
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                    -1.0_dp,admm_env%A,admm_env%work_orb_orb,1.0_dp,&
                    admm_env%work_aux_orb,error)

    ! *** copy to sparse matrix
    CALL copy_fm_to_dbcsr(admm_env%work_aux_orb, matrix_w_q,keep_sparsity=.TRUE.)
 
    ! *** A*C*Y^(T)*C^(T)
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                    1.0_dp,admm_env%A,admm_env%work_orb_orb,0.0_dp,&
                    admm_env%work_aux_orb,error)
    ! *** A*C*Y^(T)*C^(T)*A^(T) add to aux_aux, minus sign cancels
    CALL cp_fm_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nao_orb,&
                    1.0_dp,admm_env%work_aux_orb,admm_env%A,1.0_dp,&
                    admm_env%work_aux_aux,error)

    ! *** copy to sparse matrix
    CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, matrix_w_s,keep_sparsity=.TRUE.)

    ! *** This can be done in one call w_total = w_alpha + w_beta

    CALL calc_overlap_force_general(qs_env,para_env,matrix_w_s,&
                                   use_aux_fit_basis_set, use_aux_fit_basis_set, qs_env%sab_aux_fit_asymm, &
                                   do_symmetric=.FALSE., error=error)

    CALL calc_overlap_force_general(qs_env,para_env,matrix_w_q,&
                                    use_aux_fit_basis_set, use_orb_basis_set, qs_env%sab_aux_fit_vs_orb, &
                                    do_symmetric=.FALSE., error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT",&
            extension=".Log",error=error)
       CALL dbcsr_write_sparse_matrix(matrix_w_s,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT", error=error)
    END IF
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT",&
            extension=".Log",error=error)
       CALL dbcsr_write_sparse_matrix(matrix_w_q,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT", error=error)
    END IF

    ! *** Deallocated weighted density matrices
    CALL dbcsr_deallocate_matrix(matrix_w_s,error)
    CALL dbcsr_deallocate_matrix(matrix_w_q,error)

    CALL timestop(handle)

  END SUBROUTINE calc_mixed_overlap_force

! *****************************************************************************
!> \brief Calculates dExc/dR from contribution of ADMM and adds the
!>        forces to force(:)%overlap_admm(:,:). The routines basically
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
    TYPE(dbcsr_obj), POINTER                 :: matrix_w
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
      maxsgfb, natom, ncoa, ncob, nder, neighbor_list_id, nkind, nlist, &
      nnode, nrow, nseta, nsetb, sgfa, sgfb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found, new_atom_b, &
                                                return_s_derivatives, tr
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
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
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

    return_s_derivatives = .TRUE.

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    neighbor_list_id=neighbor_list_id,&
                    force=force,&
                    distribution_2d=distribution_2d,&
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
              CALL dbcsr_get_block_p(matrix=matrix_w,row=irow,col=icol,&
                   block=w_block,tr=tr,found=found)

              IF (.NOT.ASSOCIATED(w_block)) CYCLE
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
                      force(ikind)%overlap_admm(:,atom_a)=force(ikind)%overlap_admm(:,atom_a) - 2.0_dp*force_a(:)
                      force(jkind)%overlap_admm(:,atom_b)=force(jkind)%overlap_admm(:,atom_b) + 2.0_dp*force_a(:)
                    ELSE
                      force(ikind)%overlap_admm(:,atom_a)=force(ikind)%overlap_admm(:,atom_a) - force_a(:)
                      force(jkind)%overlap_admm(:,atom_b)=force(jkind)%overlap_admm(:,atom_b) + force_a(:)
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
!> \brief This routine modifies the xc section depending on the potential type
!>        used for the HF exchange and the resulting correction term. Currently 
!>        three types of corrections are implemented:
!>        
!>        coulomb:     Ex,hf = Ex,hf' + (PBEx-PBEx')
!>        shortrange:  Ex,hf = Ex,hf' + (XWPBEX-XWPBEX')
!>        truncated:   Ex,hf = Ex,hf' + ( (XWPBEX0-PBE_HOLE_TC_LR) -(XWPBEX0-PBE_HOLE_TC_LR)' )
!>        
!>        with ' denoting the auxiliary basis set and
!>      
!>          PBEx:           PBE exchange functional
!>          XWPBEX:         PBE exchange hole for short-range potential (erfc(omega*r)/r)
!>          XWPBEX0:        PBE exchange hole for standard coulomb potential
!>          PBE_HOLE_TC_LR: PBE exchange hole for longrange truncated coulomb potential
!>      
!>      
!> \param qs_env the qs environment 
!> \param xc_section the original xc_section
!> \param admm_env the ADMM environment
!> \param error 
!>
!> \par History
!>      12.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE create_admm_xc_section(qs_env, xc_section, admm_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: hfx_potential_type, ifun, nfun
    LOGICAL                                  :: funct_found
    REAL(dp)                                 :: cutoff_radius, hfx_fraction, &
                                                omega, scale_x
    TYPE(section_vals_type), POINTER         :: xc_fun, xc_fun_section

    NULLIFY(admm_env%xc_section_aux, admm_env%xc_section_primary)

    !! ** Duplicate existing xc-section

    CALL section_vals_duplicate(xc_section,admm_env%xc_section_aux,error=error)
    CALL section_vals_duplicate(xc_section,admm_env%xc_section_primary,error=error)

    !** Now modify the auxiliary basis
    !** First remove all functionals
    xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_aux,"XC_FUNCTIONAL",error=error)

    !* Overwrite possible shortcut
    CALL section_vals_val_set(xc_fun_section,"_SECTION_PARAMETERS_",&
                              i_val=xc_funct_no_shortcut,error=error)

    !** Get number of Functionals in the list
    ifun = 0
    nfun = 0
    DO
      ifun = ifun+1
      xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
      IF (.NOT.ASSOCIATED(xc_fun)) EXIT
      nfun = nfun + 1
    END DO

    ifun = 0
    DO ifun = 1,nfun
      xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=1,error=error)
      IF (.NOT.ASSOCIATED(xc_fun)) EXIT
      CALL section_vals_remove_values(xc_fun, error=error)
    END DO

    hfx_potential_type = qs_env%x_data(1,1)%potential_parameter%potential_type
    hfx_fraction = qs_env%x_data(1,1)%general_parameter%fraction

    !! ** Add functionals evaluated with auxiliary basis
    SELECT CASE (hfx_potential_type)
    CASE (do_hfx_potential_coulomb)
      CALL section_vals_val_set(xc_fun_section,"PBE%_SECTION_PARAMETERS_",&
                                l_val=.TRUE.,error=error)
      CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_X",&
                                r_val=-hfx_fraction,error=error)
      CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_C",&
                                r_val=0.0_dp,error=error)
    CASE (do_hfx_potential_short)
      omega =  qs_env%x_data(1,1)%potential_parameter%omega
      CALL section_vals_val_set(xc_fun_section,"XWPBE%_SECTION_PARAMETERS_",&
                                l_val=.TRUE.,error=error)
      CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                r_val=-hfx_fraction,error=error)
      CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                r_val=0.0_dp,error=error)
      CALL section_vals_val_set(xc_fun_section,"XWPBE%OMEGA",&
                                r_val=omega,error=error)
    CASE (do_hfx_potential_truncated)
      cutoff_radius = qs_env%x_data(1,1)%potential_parameter%cutoff_radius
      CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_",&
                                l_val=.TRUE.,error=error)
      CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%SCALE_X",&
                                r_val=hfx_fraction,error=error)
      CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%CUTOFF_RADIUS",&
                                r_val=cutoff_radius,error=error)
      CALL section_vals_val_set(xc_fun_section,"XWPBE%_SECTION_PARAMETERS_",&
                                l_val=.TRUE.,error=error)
      CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                r_val=0.0_dp,error=error)
      CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                r_val=-hfx_fraction,error=error)
    END SELECT


    !** Now modify the functionals for the primary basis
    xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary,"XC_FUNCTIONAL",error=error)
    !* Overwrite possible shortcut
    CALL section_vals_val_set(xc_fun_section,"_SECTION_PARAMETERS_",&
                              i_val=xc_funct_no_shortcut,error=error)


    SELECT CASE (hfx_potential_type)
    CASE (do_hfx_potential_coulomb)
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        IF( xc_fun%section%name == "PBE" ) THEN
          funct_found = .TRUE.
        END IF
      END DO
      IF( .NOT. funct_found ) THEN
        CALL section_vals_val_set(xc_fun_section,"PBE%_SECTION_PARAMETERS_",&
                                  l_val=.TRUE.,error=error)
        CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_X",&
                                  r_val=hfx_fraction,error=error)
        CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_C",&
                                  r_val=0.0_dp,error=error)
      ELSE
        CALL section_vals_val_get(xc_fun_section,"PBE%SCALE_X",&
                                  r_val=scale_x,error=error)
        scale_x = scale_x + hfx_fraction
        CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_X",&
                                  r_val=scale_x,error=error)
      END IF
    CASE (do_hfx_potential_short)
      omega =  qs_env%x_data(1,1)%potential_parameter%omega
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        IF( xc_fun%section%name == "XWPBE" ) THEN
          funct_found = .TRUE.
        END IF
      END DO
      IF( .NOT. funct_found ) THEN
        CALL section_vals_val_set(xc_fun_section,"XWPBE%_SECTION_PARAMETERS_",&
                                  l_val=.TRUE.,error=error)
        CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                  r_val=hfx_fraction,error=error)
        CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                  r_val=0.0_dp,error=error)
        CALL section_vals_val_set(xc_fun_section,"XWPBE%OMEGA",&
                                  r_val=omega,error=error)
      ELSE
        CALL section_vals_val_get(xc_fun_section,"XWPBE%SCALE_X",&
                                  r_val=scale_x,error=error)
        scale_x = scale_x + hfx_fraction
        CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                  r_val=scale_x,error=error)
      END IF
    CASE (do_hfx_potential_truncated)
      cutoff_radius =  qs_env%x_data(1,1)%potential_parameter%cutoff_radius
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        IF( xc_fun%section%name == "PBE_HOLE_T_C_LR" ) THEN
          funct_found = .TRUE.
        END IF
      END DO
      IF( .NOT. funct_found ) THEN
        CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_",&
                                  l_val=.TRUE.,error=error)
        CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%SCALE_X",&
                                  r_val=-hfx_fraction,error=error)
        CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%CUTOFF_RADIUS",&
                                  r_val=cutoff_radius,error=error)

      ELSE
        CALL section_vals_val_get(xc_fun_section,"PBE_HOLE_T_C_LR%SCALE_X",&
                                  r_val=scale_x,error=error)
        scale_x = scale_x - hfx_fraction
        CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%SCALE_X",&
                                  r_val=scale_x,error=error)
        CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%CUTOFF_RADIUS",&
                                  r_val=cutoff_radius,error=error)
      END IF
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        IF( xc_fun%section%name == "XWPBE" ) THEN
          funct_found = .TRUE.
        END IF
      END DO
      IF( .NOT. funct_found ) THEN
        CALL section_vals_val_set(xc_fun_section,"XWPBE%_SECTION_PARAMETERS_",&
                                  l_val=.TRUE.,error=error)
        CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                  r_val=hfx_fraction,error=error)
        CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                  r_val=0.0_dp,error=error)

      ELSE
        CALL section_vals_val_get(xc_fun_section,"XWPBE%SCALE_X0",&
                                  r_val=scale_x,error=error)
        scale_x = scale_x + hfx_fraction
        CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                  r_val=scale_x,error=error)
      END IF

    END SELECT


    IF( 0==1 ) THEN
      WRITE(*,*) "primary"
      xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary,"XC_FUNCTIONAL",error=error)
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT

        scale_x=-1000.0_dp
        IF(xc_fun%section%name /= "LYP" .AND. xc_fun%section%name /= "VWN") THEN
          CALL section_vals_val_get(xc_fun,"SCALE_X",&
                                      r_val=scale_x,error=error)
        END IF
        IF( xc_fun%section%name == "XWPBE" ) THEN
          CALL section_vals_val_get(xc_fun,"SCALE_X0",&
                                    r_val=hfx_fraction,error=error)

           WRITE(*,*) xc_fun%section%name, scale_x, hfx_fraction 
        ELSE
          WRITE(*,*) xc_fun%section%name, scale_x
        END IF
      END DO

      WRITE(*,*) "auxiliary"
      xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_aux,"XC_FUNCTIONAL",error=error)
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        scale_x=-1000.0_dp
        IF(xc_fun%section%name /= "LYP" .AND. xc_fun%section%name /= "VWN") THEN
          CALL section_vals_val_get(xc_fun,"SCALE_X",&
                                      r_val=scale_x,error=error)
        END IF
        IF( xc_fun%section%name == "XWPBE" ) THEN
          CALL section_vals_val_get(xc_fun,"SCALE_X0",&
                                    r_val=hfx_fraction,error=error)

           WRITE(*,*) xc_fun%section%name, scale_x, hfx_fraction
        ELSE
          WRITE(*,*) xc_fun%section%name, scale_x
        END IF
      END DO
    END IF

  END SUBROUTINE create_admm_xc_section

  SUBROUTINE merge_ks_matrix(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, &
                             matrix_ks, matrix_ks_aux_fit, matrix_s, matrix_p_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_s, matrix_p_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nao_orb, &
                                                nmo
    INTEGER, SAVE                            :: counter = 0
    LOGICAL                                  :: failure
    TYPE(dbcsr_obj), POINTER                 :: matrix_k_tilde

    failure = .FALSE.

    CALL timeset(routineN,handle)

    counter = counter + 1
    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    IF ( 1==0 ) THEN
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,admm_env%work_orb_orb,error)
      CALL cp_fm_upper_to_full(admm_env%work_orb_orb,admm_env%work_orb_orb2,error=error)
      CALL print_matlab_matrix(admm_env%work_orb_orb, "S_real.dat", counter)
      CALL print_matlab_matrix(mo_coeff,"C.dat",counter)
      CALL print_matlab_matrix(mo_coeff_aux_fit,"C_tilde.dat",counter)
      CALL print_matlab_matrix(admm_env%Q,"Q.dat",counter)
      CALL print_matlab_matrix(admm_env%S,"S.dat",counter)
    END IF

    !! Calculate Lambda^{-2}
    CALL cp_fm_to_fm(admm_env%lambda(ispin)%matrix, admm_env%work_nmo_nmo1(ispin)%matrix, error=error) 
    CALL cp_fm_cholesky_decompose(admm_env%work_nmo_nmo1(ispin)%matrix,error=error)
    CALL cp_fm_cholesky_invert(admm_env%work_nmo_nmo1(ispin)%matrix,error=error) 
    !! Symmetrize the guy
    CALL cp_fm_upper_to_full(admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%lambda_inv2(ispin)%matrix,error=error)
    !! Take square
    CALL cp_fm_gemm('N', 'T',  nmo, nmo, nmo,&
                    1.0_dp,admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                    admm_env%lambda_inv2(ispin)%matrix,error)

    !! ** C_hat = AC
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nmo, nao_orb,&
                    1.0_dp,admm_env%A,mo_coeff,0.0_dp,&
                    admm_env%C_hat(ispin)%matrix,error)

    !! ** P_tilde
    CALL cp_fm_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nmo,&
                    1.0_dp,mo_coeff_aux_fit,mo_coeff_aux_fit,0.0_dp,&
                    admm_env%P_tilde(ispin)%matrix,error)
    

    !! ** C_hat*Lambda^{-2}
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%C_hat(ispin)%matrix,admm_env%lambda_inv2(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix,error)

    !! ** C_hat*Lambda^{-2}*C_hat^T
    CALL cp_fm_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nmo,&
                    1.0_dp,admm_env%work_aux_nmo(ispin)%matrix,admm_env%C_hat(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux,error)


    !! ** S*C_hat*Lambda^{-2}*C_hat^T 
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%S,admm_env%work_aux_aux,0.0_dp,&
                    admm_env%work_aux_aux2,error)


    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)

    !! ** S*C_hat*Lambda^{-2}*C_hat^T*H_tilde 
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux2,admm_env%K(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux,error)

    !! ** P_tilde*S
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%P_tilde(ispin)%matrix,admm_env%S,0.0_dp,&
                    admm_env%work_aux_aux2,error)


    !! ** -S*C_hat*Lambda^{-2}*C_hat^T*H_tilde*P_tilde*S
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    -1.0_dp,admm_env%work_aux_aux,admm_env%work_aux_aux2,0.0_dp,&
                    admm_env%work_aux_aux3,error)

    
    !! ** -S*C_hat*Lambda^{-2}*C_hat^T*H_tilde*P_tilde*S+S*C_hat*Lambda^{-2}*C_hat^T*H_tilde
    CALL cp_fm_scale_and_add(1.0_dp,admm_env%work_aux_aux3,1.0_dp,admm_env%work_aux_aux,error)


    !! first_part*A
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux3,admm_env%A,0.0_dp,&
                    admm_env%work_aux_orb,error)

    !! + first_part^T*A
    CALL cp_fm_gemm('T', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux3,admm_env%A,1.0_dp,&
                    admm_env%work_aux_orb,error)
    


    !! A^T*(first+seccond)=H
    CALL cp_fm_gemm('T', 'N',  nao_orb, nao_orb, nao_aux_fit,&
                    1.0_dp,admm_env%A,admm_env%work_aux_orb,0.0_dp,&
                    admm_env%work_orb_orb,error)


    NULLIFY(matrix_k_tilde)
    ALLOCATE(matrix_k_tilde)
    CALL dbcsr_init (matrix_k_tilde, error)
    CALL dbcsr_create(matrix_k_tilde, 'MATRIX K_tilde', &
         dbcsr_distribution(matrix_ks(ispin)%matrix), dbcsr_type_symmetric, dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
         dbcsr_col_block_sizes(matrix_ks(ispin)%matrix), matrix_ks(ispin)%matrix%m%nblks, matrix_ks(ispin)%matrix%m%nze,&
         matrix_ks(ispin)%matrix%m%data_type, dbcsr_uses_special_memory(matrix_ks(ispin)%matrix), error=error)
    CALL dbcsr_finalize(matrix_k_tilde, error=error)


    CALL cp_fm_to_fm(admm_env%work_orb_orb, admm_env%ks_to_be_merged(ispin)%matrix, error=error) 

    CALL dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix, error=error)
    CALL dbcsr_set(matrix_k_tilde, 0.0_dp, error)
    CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE.)

    CALL cp_fm_gemm('N', 'N',  nao_orb, nmo, nao_orb,&
                    1.0_dp,admm_env%work_orb_orb,mo_coeff,0.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix,error)
    

    CALL dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp, error)   

    CALL copy_dbcsr_to_fm(matrix_ks(1)%matrix,admm_env%work_orb_orb,error)
    CALL cp_fm_upper_to_full(admm_env%work_orb_orb,admm_env%work_orb_orb2,error=error) 
    CALL print_matlab_matrix(admm_env%work_orb_orb, "KS_full.dat", counter)


    CALL dbcsr_deallocate_matrix(matrix_k_tilde,error) 

    CALL timestop(handle)

  END SUBROUTINE merge_ks_matrix

  SUBROUTINE print_matlab_matrix(matrix, filename, counter)
    TYPE(cp_fm_type), POINTER                :: matrix
    CHARACTER(LEN=*)                         :: filename
    INTEGER                                  :: counter

    CHARACTER(LEN=20)                        :: string
    CHARACTER(LEN=50)                        :: filename_count
    INTEGER                                  :: i, j

    WRITE(filename_count,FMT='(A,I0)') filename, counter
    OPEN (unit=120,file=TRIM(filename_count))
    j=SIZE(matrix%local_data,2)
    WRITE(string,FMT='(A1,I4,A7)') "(",j,"F20.16)"
    DO i=1,SIZE(matrix%local_data,1)
      WRITE(unit=120,FMT=string) matrix%local_data(i,:)
    END DO
    CLOSE (unit=120)
  END SUBROUTINE
  
  SUBROUTINE remove_ks_matrix(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, &
                             matrix_ks, matrix_ks_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'remove_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(dbcsr_obj), POINTER                 :: matrix_k_tilde

    failure = .FALSE.

    CALL timeset(routineN,handle)


    NULLIFY(matrix_k_tilde)
    ALLOCATE(matrix_k_tilde)
    CALL dbcsr_init (matrix_k_tilde, error)
    CALL dbcsr_create(matrix_k_tilde, 'MATRIX K_tilde', &
         dbcsr_distribution(matrix_ks(ispin)%matrix), dbcsr_type_symmetric, dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
         dbcsr_col_block_sizes(matrix_ks(ispin)%matrix), matrix_ks(ispin)%matrix%m%nblks, matrix_ks(ispin)%matrix%m%nze,&
         matrix_ks(ispin)%matrix%m%data_type, dbcsr_uses_special_memory(matrix_ks(ispin)%matrix), error=error)
    CALL dbcsr_finalize(matrix_k_tilde, error=error)


    CALL cp_fm_to_fm(admm_env%work_orb_orb, admm_env%ks_to_be_merged(ispin)%matrix, error=error) 

    CALL dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix, error=error)
    CALL dbcsr_set(matrix_k_tilde, 0.0_dp, error)
    CALL copy_fm_to_dbcsr(admm_env%ks_to_be_merged(ispin)%matrix, matrix_k_tilde, keep_sparsity=.TRUE.)

    CALL dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, -1.0_dp, error)   

    CALL dbcsr_deallocate_matrix(matrix_k_tilde,error) 

    CALL timestop(handle)

  END SUBROUTINE remove_ks_matrix

END MODULE admm_methods
