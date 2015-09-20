!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE cp_control_types,                ONLY: admm_control_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: do_admm_blocked_projection,&
                                             do_admm_blocking_purify_full,&
                                             do_admm_charge_constrained_projection
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: admm_env_create, admm_env_release, admm_type


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_types'

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
!> \param gsi contains ratio N_dens_m/N_aux_dens_m 
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  TYPE admm_type
    TYPE(cp_fm_type), POINTER                :: S_inv => Null(),&
                                                S => Null(),&
                                                Q => Null(),&
                                                A => Null(),&
                                                B => Null(),&
                                                work_orb_orb => Null(),&
                                                work_orb_orb2 => Null(),&
                                                work_orb_orb3 => Null(),&
                                                work_aux_orb => Null(),&
                                                work_aux_orb2 => Null(),&
                                                work_aux_orb3 => Null(),&
                                                work_aux_aux => Null(),&
                                                work_aux_aux2 => Null(),&
                                                work_aux_aux3 => Null(),&
                                                work_aux_aux4 => Null(),&
                                                work_aux_aux5 => Null()

    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: lambda => Null(),&
                                                lambda_inv => Null(),&
                                                lambda_inv_sqrt => Null(),&
                                                R => Null(),&
                                                R_purify => Null(),&
                                                work_orb_nmo => Null(),&
                                                work_nmo_nmo1 => Null(),&
                                                R_schur_R_t => Null(),&
                                                work_nmo_nmo2 => Null(),&
                                                work_aux_nmo => Null(),&
                                                work_aux_nmo2 => Null(),&
                                                H => Null(),&
                                                H_corr => Null(),&
                                                mo_derivs_tmp => Null(),&
                                                K => Null(),&
                                                M => Null(),&
                                                M_purify => Null(),&
                                                P_to_be_purified => Null(),&
                                                lambda_inv2 => Null(),&
                                                C_hat => Null(),&
                                                P_tilde => Null(),&
                                                ks_to_be_merged => Null()
    TYPE(eigvals_p_type), DIMENSION(:), &
      POINTER                                :: eigvals_lambda => Null(),&
                                                eigvals_P_to_be_purified => Null()
    TYPE(section_vals_type), POINTER         :: xc_section_primary => Null(),&
                                                xc_section_aux => Null()
    REAL(KIND=dp)                            :: gsi(3), &   
                                                lambda_merlot(2), &
                                                n_large_basis(3) 
    INTEGER                                  :: nao_orb, nao_aux_fit, nmo(2)
    INTEGER                                  :: purification_method
    LOGICAL                                  :: charge_constrain = .FALSE.  
    INTEGER                                  :: scaling_model, &
                                                aux_exch_func 
    LOGICAL                                  :: block_dm  = .FALSE.
    LOGICAL                                  :: block_fit = .FALSE.
    INTEGER, DIMENSION(:,:), POINTER         :: block_map => Null()
  END TYPE



  CONTAINS

! *****************************************************************************
!> \brief creates ADDE environment, initializes the basic types
!>
!> \param admm_env The ADMM env
!> \param admm_control ...
!> \param mos the MO's of the orbital basis set
!> \param mos_aux_fit the MO's of the auxiliary fitting basis set
!> \param para_env The parallel env
!> \param natoms ...
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE admm_env_create(admm_env, admm_control, mos, mos_aux_fit, para_env, natoms)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(admm_control_type), POINTER         :: admm_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: natoms

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, iblock, ispin, j, &
                                                jatom, nao_aux_fit, nao_orb, &
                                                nmo, nspins
    TYPE(cp_fm_struct_type), POINTER :: fm_struct_aux_aux, fm_struct_aux_nmo, &
      fm_struct_aux_orb, fm_struct_nmo_nmo, fm_struct_orb_aux, &
      fm_struct_orb_nmo, fm_struct_orb_orb
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL cite_reference(Guidon2010)

    ALLOCATE(admm_env)

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
                             para_env=para_env)
    CALL cp_fm_struct_create(fm_struct_aux_orb,&
                             context=mo_coeff%matrix_struct%context,&
                             nrow_global=nao_aux_fit,&
                             ncol_global=nao_orb,&
                             para_env=para_env)
    CALL cp_fm_struct_create(fm_struct_orb_aux,&
                             context=mo_coeff%matrix_struct%context,&
                             nrow_global=nao_orb,&
                             ncol_global=nao_aux_fit,&
                             para_env=para_env)
    CALL cp_fm_struct_create(fm_struct_orb_orb,&
                             context=mo_coeff%matrix_struct%context,&
                             nrow_global=nao_orb,&
                             ncol_global=nao_orb,&
                             para_env=para_env)

    CALL cp_fm_create(admm_env%S,fm_struct_aux_aux,name="aux_fit_overlap")
    CALL cp_fm_create(admm_env%S_inv,fm_struct_aux_aux,name="aux_fit_overlap_inv")
    CALL cp_fm_create(admm_env%Q,fm_struct_aux_orb,name="mixed_overlap")
    CALL cp_fm_create(admm_env%A,fm_struct_aux_orb,name="work_A")
    CALL cp_fm_create(admm_env%B,fm_struct_orb_orb,name="work_B")
    CALL cp_fm_create(admm_env%work_orb_orb,fm_struct_orb_orb,name="work_orb_orb")
    CALL cp_fm_create(admm_env%work_orb_orb2,fm_struct_orb_orb,name="work_orb_orb")
    CALL cp_fm_create(admm_env%work_orb_orb3,fm_struct_orb_orb,name="work_orb_orb3")
    CALL cp_fm_create(admm_env%work_aux_orb,fm_struct_aux_orb,name="work_aux_orb")
    CALL cp_fm_create(admm_env%work_aux_orb2,fm_struct_aux_orb,name="work_aux_orb2")
    CALL cp_fm_create(admm_env%work_aux_orb3,fm_struct_aux_orb,name="work_aux_orb3")
    CALL cp_fm_create(admm_env%work_aux_aux,fm_struct_aux_aux,name="work_aux_aux")
    CALL cp_fm_create(admm_env%work_aux_aux2,fm_struct_aux_aux,name="work_aux_aux2")
    CALL cp_fm_create(admm_env%work_aux_aux3,fm_struct_aux_aux,name="work_aux_aux3")
    CALL cp_fm_create(admm_env%work_aux_aux4,fm_struct_aux_aux,name="work_aux_aux4")
    CALL cp_fm_create(admm_env%work_aux_aux5,fm_struct_aux_aux,name="work_aux_aux5")

    ALLOCATE(admm_env%lambda_inv(nspins))
    ALLOCATE(admm_env%lambda(nspins))
    ALLOCATE(admm_env%lambda_inv_sqrt(nspins))
    ALLOCATE(admm_env%R(nspins))
    ALLOCATE(admm_env%R_purify(nspins))
    ALLOCATE(admm_env%work_orb_nmo(nspins))
    ALLOCATE(admm_env%work_nmo_nmo1(nspins))
    ALLOCATE(admm_env%R_schur_R_t(nspins))
    ALLOCATE(admm_env%work_nmo_nmo2(nspins))
    ALLOCATE(admm_env%eigvals_lambda(nspins))
    ALLOCATE(admm_env%eigvals_P_to_be_purified(nspins))
    ALLOCATE(admm_env%H(nspins))
    ALLOCATE(admm_env%K(nspins))
    ALLOCATE(admm_env%M(nspins))
    ALLOCATE(admm_env%M_purify(nspins))
    ALLOCATE(admm_env%P_to_be_purified(nspins))
    ALLOCATE(admm_env%work_aux_nmo(nspins))
    ALLOCATE(admm_env%work_aux_nmo2(nspins))
    ALLOCATE(admm_env%mo_derivs_tmp(nspins))
    ALLOCATE(admm_env%H_corr(nspins))
    ALLOCATE(admm_env%ks_to_be_merged(nspins))
    ALLOCATE(admm_env%lambda_inv2(nspins))
    ALLOCATE(admm_env%C_hat(nspins))
    ALLOCATE(admm_env%P_tilde(nspins))

    DO ispin = 1,nspins
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff, nmo=nmo)
      admm_env%nmo(ispin)=nmo
      CALL cp_fm_struct_create(fm_struct_aux_nmo,&
                               context=mo_coeff%matrix_struct%context,&
                               nrow_global=nao_aux_fit,&
                               ncol_global=nmo,&
                               para_env=para_env)
      CALL cp_fm_struct_create(fm_struct_orb_nmo,&
                               context=mo_coeff%matrix_struct%context,&
                               nrow_global=nao_orb,&
                               ncol_global=nmo,&
                               para_env=para_env)
      CALL cp_fm_struct_create(fm_struct_nmo_nmo,&
                               context=mo_coeff%matrix_struct%context,&
                               nrow_global=nmo,&
                               ncol_global=nmo,&
                               para_env=para_env)

      CALL cp_fm_create(admm_env%work_orb_nmo(ispin)%matrix,fm_struct_orb_nmo,name="work_orb_nmo")
      CALL cp_fm_create(admm_env%work_nmo_nmo1(ispin)%matrix,fm_struct_nmo_nmo,name="work_nmo_nmo1")
      CALL cp_fm_create(admm_env%R_schur_R_t(ispin)%matrix,fm_struct_nmo_nmo,name="R_schur_R_t")
      CALL cp_fm_create(admm_env%work_nmo_nmo2(ispin)%matrix,fm_struct_nmo_nmo,name="work_nmo_nmo2")
      CALL cp_fm_create(admm_env%lambda(ispin)%matrix,fm_struct_nmo_nmo,name="lambda")
      CALL cp_fm_create(admm_env%lambda_inv(ispin)%matrix,fm_struct_nmo_nmo,name="lambda_inv")
      CALL cp_fm_create(admm_env%lambda_inv_sqrt(ispin)%matrix,fm_struct_nmo_nmo,name="lambda_inv_sqrt")
      CALL cp_fm_create(admm_env%R(ispin)%matrix,fm_struct_nmo_nmo,name="R")
      CALL cp_fm_create(admm_env%R_purify(ispin)%matrix,fm_struct_aux_aux,name="R_purify")
      CALL cp_fm_create(admm_env%K(ispin)%matrix,fm_struct_aux_aux,name="K")
      CALL cp_fm_create(admm_env%H(ispin)%matrix,fm_struct_aux_nmo,name="H")
      CALL cp_fm_create(admm_env%H_corr(ispin)%matrix,fm_struct_orb_orb,name="H_corr")
      CALL cp_fm_create(admm_env%M(ispin)%matrix,fm_struct_nmo_nmo,name="M")
      CALL cp_fm_create(admm_env%M_purify(ispin)%matrix,fm_struct_aux_aux,name="M aux")
      CALL cp_fm_create(admm_env%P_to_be_purified(ispin)%matrix,fm_struct_aux_aux,name="P_to_be_purified")
      CALL cp_fm_create(admm_env%work_aux_nmo(ispin)%matrix,fm_struct_aux_nmo,name="work_aux_nmo")
      CALL cp_fm_create(admm_env%work_aux_nmo2(ispin)%matrix,fm_struct_aux_nmo,name="work_aux_nmo2")
      CALL cp_fm_create(admm_env%mo_derivs_tmp(ispin)%matrix,fm_struct_orb_nmo,name="mo_derivs_tmp")

      CALL cp_fm_create(admm_env%lambda_inv2(ispin)%matrix,fm_struct_nmo_nmo,name="lambda_inv2")
      CALL cp_fm_create(admm_env%C_hat(ispin)%matrix,fm_struct_aux_nmo,name="C_hat")
      CALL cp_fm_create(admm_env%P_tilde(ispin)%matrix,fm_struct_aux_aux,name="P_tilde")

      CALL cp_fm_create(admm_env%ks_to_be_merged(ispin)%matrix,fm_struct_orb_orb,name="KS_to_be_merged ")

      ALLOCATE(admm_env%eigvals_lambda(ispin)%eigvals)
      ALLOCATE(admm_env%eigvals_P_to_be_purified(ispin)%eigvals)
      ALLOCATE(admm_env%eigvals_lambda(ispin)%eigvals%data(nmo))
      ALLOCATE(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(nao_aux_fit))
      admm_env%eigvals_lambda(ispin)%eigvals%data = 0.0_dp
      admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data = 0.0_dp
      CALL cp_fm_struct_release(fm_struct_aux_nmo)
      CALL cp_fm_struct_release(fm_struct_orb_nmo)
      CALL cp_fm_struct_release(fm_struct_nmo_nmo)
    END DO

    CALL cp_fm_struct_release(fm_struct_aux_aux)
    CALL cp_fm_struct_release(fm_struct_aux_orb)
    CALL cp_fm_struct_release(fm_struct_orb_aux)
    CALL cp_fm_struct_release(fm_struct_orb_orb)

    ! copy settings from admm_control
    admm_env%purification_method = admm_control%purification_method
    admm_env%scaling_model = admm_control%scaling_model 
    admm_env%aux_exch_func = admm_control%aux_exch_func
    admm_env%charge_constrain = (admm_control%method == do_admm_charge_constrained_projection)     
    admm_env%block_dm  = ((admm_control%method == do_admm_blocking_purify_full) .OR. &
                          (admm_control%method == do_admm_blocked_projection ))
    admm_env%block_fit = admm_control%method == do_admm_blocked_projection


    IF((admm_control%method == do_admm_blocking_purify_full)  .OR. &
       (admm_control%method == do_admm_blocked_projection)) THEN
       ! create block map
       ALLOCATE(admm_env%block_map(natoms,natoms))
       admm_env%block_map(:,:) = 0
       DO iblock = 1, SIZE(admm_control%blocks)
         DO i = 1,SIZE(admm_control%blocks(iblock)%list)
           iatom = admm_control%blocks(iblock)%list(i)
           DO j = 1,SIZE(admm_control%blocks(iblock)%list)
             jatom = admm_control%blocks(iblock)%list(j)
             admm_env%block_map(iatom,jatom) = 1
           END DO
         END DO
       END DO
     ENDIF
  END SUBROUTINE admm_env_create


! *****************************************************************************
!> \brief releases the ADMM environment, cleans up all types
!>
!> \param admm_env The ADMM env
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE admm_env_release(admm_env)
    TYPE(admm_type), POINTER                 :: admm_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin

    CALL cp_fm_release(admm_env%S)
    CALL cp_fm_release(admm_env%S_inv)
    CALL cp_fm_release(admm_env%Q)
    CALL cp_fm_release(admm_env%A)
    CALL cp_fm_release(admm_env%B)
    CALL cp_fm_release(admm_env%work_orb_orb)
    CALL cp_fm_release(admm_env%work_orb_orb2)
    CALL cp_fm_release(admm_env%work_orb_orb3)
    CALL cp_fm_release(admm_env%work_aux_aux)
    CALL cp_fm_release(admm_env%work_aux_aux2)
    CALL cp_fm_release(admm_env%work_aux_aux3)
    CALL cp_fm_release(admm_env%work_aux_aux4)
    CALL cp_fm_release(admm_env%work_aux_aux5)
    CALL cp_fm_release(admm_env%work_aux_orb)
    CALL cp_fm_release(admm_env%work_aux_orb2)
    CALL cp_fm_release(admm_env%work_aux_orb3)
    DO ispin = 1,SIZE(admm_env%lambda)
      CALL cp_fm_release(admm_env%lambda(ispin)%matrix)
      CALL cp_fm_release(admm_env%lambda_inv(ispin)%matrix)
      CALL cp_fm_release(admm_env%lambda_inv_sqrt(ispin)%matrix)
      CALL cp_fm_release(admm_env%lambda_inv2(ispin)%matrix)
      CALL cp_fm_release(admm_env%C_hat(ispin)%matrix)
      CALL cp_fm_release(admm_env%P_tilde(ispin)%matrix)
      CALL cp_fm_release(admm_env%R(ispin)%matrix)
      CALL cp_fm_release(admm_env%R_purify(ispin)%matrix)
      CALL cp_fm_release(admm_env%H(ispin)%matrix)
      CALL cp_fm_release(admm_env%H_corr(ispin)%matrix)
      CALL cp_fm_release(admm_env%K(ispin)%matrix)
      CALL cp_fm_release(admm_env%M(ispin)%matrix)
      CALL cp_fm_release(admm_env%M_purify(ispin)%matrix)
      CALL cp_fm_release(admm_env%P_to_be_purified(ispin)%matrix)
      CALL cp_fm_release(admm_env%work_orb_nmo(ispin)%matrix)
      CALL cp_fm_release(admm_env%work_nmo_nmo1(ispin)%matrix)
      CALL cp_fm_release(admm_env%R_schur_R_t(ispin)%matrix)
      CALL cp_fm_release(admm_env%work_nmo_nmo2(ispin)%matrix)
      CALL cp_fm_release(admm_env%work_aux_nmo(ispin)%matrix)
      CALL cp_fm_release(admm_env%work_aux_nmo2(ispin)%matrix)
      CALL cp_fm_release(admm_env%mo_derivs_tmp(ispin)%matrix)
      CALL cp_fm_release(admm_env%ks_to_be_merged(ispin)%matrix)
      CALL cp_fm_release(admm_env%lambda_inv2(ispin)%matrix)
      DEALLOCATE(admm_env%eigvals_lambda(ispin)%eigvals%data)
      DEALLOCATE(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data)
      DEALLOCATE(admm_env%eigvals_lambda(ispin)%eigvals)
      DEALLOCATE(admm_env%eigvals_P_to_be_purified(ispin)%eigvals)
    END DO
    DEALLOCATE(admm_env%eigvals_lambda)
    DEALLOCATE(admm_env%eigvals_P_to_be_purified)
    DEALLOCATE(admm_env%lambda)
    DEALLOCATE(admm_env%lambda_inv)
    DEALLOCATE(admm_env%lambda_inv_sqrt)
    DEALLOCATE(admm_env%R)
    DEALLOCATE(admm_env%R_purify)
    DEALLOCATE(admm_env%M)
    DEALLOCATE(admm_env%M_purify)
    DEALLOCATE(admm_env%P_to_be_purified)
    DEALLOCATE(admm_env%H)
    DEALLOCATE(admm_env%H_corr)
    DEALLOCATE(admm_env%K)
    DEALLOCATE(admm_env%work_orb_nmo)
    DEALLOCATE(admm_env%work_nmo_nmo1)
    DEALLOCATE(admm_env%R_schur_R_t)
    DEALLOCATE(admm_env%work_nmo_nmo2)
    DEALLOCATE(admm_env%work_aux_nmo)
    DEALLOCATE(admm_env%work_aux_nmo2)
    DEALLOCATE(admm_env%mo_derivs_tmp)
    DEALLOCATE(admm_env%ks_to_be_merged)
    DEALLOCATE(admm_env%lambda_inv2)
    DEALLOCATE(admm_env%C_hat)
    DEALLOCATE(admm_env%P_tilde)

    IF(ASSOCIATED(admm_env%block_map)) &
      DEALLOCATE(admm_env%block_map)

    IF(ASSOCIATED(admm_env%xc_section_primary))&
       CALL section_vals_release(admm_env%xc_section_primary)
    IF(ASSOCIATED(admm_env%xc_section_aux))&
       CALL section_vals_release(admm_env%xc_section_aux)

    DEALLOCATE(admm_env)

  END SUBROUTINE admm_env_release



END MODULE admm_types

