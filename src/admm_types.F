!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types and set/get functions for auxiliary density matrix methods
!> \par History
!>      05.2008 created [Manuel Guidon]
!>      12.2019 Made GAPW compatiblae [Augustin Bussy]
!> \author Manuel Guidon
! **************************************************************************************************
MODULE admm_types
   USE admm_dm_types,                   ONLY: admm_dm_release,&
                                              admm_dm_type
   USE bibliography,                    ONLY: Guidon2010,&
                                              cite_reference
   USE cp_control_types,                ONLY: admm_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none,&
                                              do_admm_blocked_projection,&
                                              do_admm_blocking_purify_full,&
                                              do_admm_charge_constrained_projection,&
                                              do_admm_exch_scaling_merlot,&
                                              do_admm_exch_scaling_none,&
                                              do_admm_purify_none
   USE input_section_types,             ONLY: section_vals_release,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_kind_types,                   ONLY: deallocate_qs_kind_set,&
                                              qs_kind_type
   USE qs_local_rho_types,              ONLY: local_rho_set_release,&
                                              local_rho_type
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_oce_types,                    ONLY: deallocate_oce_set,&
                                              oce_matrix_type
   USE qs_rho_types,                    ONLY: qs_rho_release,&
                                              qs_rho_type
   USE task_list_types,                 ONLY: deallocate_task_list,&
                                              task_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: admm_env_create, admm_env_release, admm_type, admm_gapw_type, set_admm_env, get_admm_env

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_types'

   TYPE eigvals_type
      REAL(dp), DIMENSION(:), POINTER          :: DATA => NULL()
   END TYPE

   TYPE eigvals_p_type
      TYPE(eigvals_type), POINTER              :: eigvals => NULL()
   END TYPE

! **************************************************************************************************
!> \brief A subtype of the admm_env that contains the extra data needed for an ADMM GAPW calculation
!> \param admm_kind_set gets its own qs_kind set to store all relevant basis/grid/etc info
!> \param local_rho_set caontains soft and hard AUX_FIT atomoc densities
!> \param task_list the task list used for all soft density pw operations
!> \param oce stores the precomputed oce integrals
! **************************************************************************************************
   TYPE admm_gapw_type
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: admm_kind_set => Null()
      TYPE(local_rho_type), POINTER                      :: local_rho_set => Null()
      TYPE(task_list_type), POINTER                      :: task_list => Null()
      TYPE(oce_matrix_type), POINTER                     :: oce => Null()
   END TYPE admm_gapw_type

! **************************************************************************************************
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
!> \param admm_gapw_env the type containing ADMM GAPW specific data
!> \param do_gapw an internal logical switch for GAPW
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   TYPE admm_type
      TYPE(cp_fm_type), POINTER                :: S_inv => Null(), &
                                                  S => Null(), &
                                                  Q => Null(), &
                                                  A => Null(), &
                                                  B => Null(), &
                                                  work_orb_orb => Null(), &
                                                  work_orb_orb2 => Null(), &
                                                  work_orb_orb3 => Null(), &
                                                  work_aux_orb => Null(), &
                                                  work_aux_orb2 => Null(), &
                                                  work_aux_orb3 => Null(), &
                                                  work_aux_aux => Null(), &
                                                  work_aux_aux2 => Null(), &
                                                  work_aux_aux3 => Null(), &
                                                  work_aux_aux4 => Null(), &
                                                  work_aux_aux5 => Null()

      TYPE(cp_fm_type), DIMENSION(:), &
         POINTER                                :: lambda => Null(), &
                                                   lambda_inv => Null(), &
                                                   lambda_inv_sqrt => Null(), &
                                                   R => Null(), &
                                                   R_purify => Null(), &
                                                   work_orb_nmo => Null(), &
                                                   work_nmo_nmo1 => Null(), &
                                                   R_schur_R_t => Null(), &
                                                   work_nmo_nmo2 => Null(), &
                                                   work_aux_nmo => Null(), &
                                                   work_aux_nmo2 => Null(), &
                                                   H => Null(), &
                                                   H_corr => Null(), &
                                                   mo_derivs_tmp => Null(), &
                                                   K => Null(), &
                                                   M => Null(), &
                                                   M_purify => Null(), &
                                                   P_to_be_purified => Null(), &
                                                   lambda_inv2 => Null(), &
                                                   C_hat => Null(), &
                                                   P_tilde => Null(), &
                                                   ks_to_be_merged => Null()
      TYPE(eigvals_p_type), DIMENSION(:), &
         POINTER                                :: eigvals_lambda => Null(), &
                                                   eigvals_P_to_be_purified => Null()
      TYPE(section_vals_type), POINTER         :: xc_section_primary => Null(), &
                                                  xc_section_aux => Null()
      REAL(KIND=dp)                            :: gsi(3) = 0.0_dp, &
                                                  lambda_merlot(2) = 0.0_dp, &
                                                  n_large_basis(3) = 0.0_dp
      INTEGER                                  :: nao_orb = 0, nao_aux_fit = 0, nmo(2) = 0
      INTEGER                                  :: purification_method = do_admm_purify_none
      LOGICAL                                  :: charge_constrain = .FALSE., do_admmp = .FALSE., &
                                                  do_admmq = .FALSE., do_admms = .FALSE.
      INTEGER                                  :: scaling_model = do_admm_exch_scaling_none, &
                                                  aux_exch_func = do_admm_aux_exch_func_none
      LOGICAL                                  :: aux_exch_func_param = .FALSE.
      REAL(KIND=dp), DIMENSION(3)              :: aux_x_param = 0.0_dp
      LOGICAL                                  :: block_dm = .FALSE.
      LOGICAL                                  :: block_fit = .FALSE.
      INTEGER, DIMENSION(:, :), POINTER        :: block_map => Null()
      TYPE(admm_gapw_type), POINTER            :: admm_gapw_env => NULL()
      LOGICAL                                  :: do_gapw = .FALSE.
      TYPE(admm_dm_type), POINTER              :: admm_dm => Null()

      TYPE(mo_set_type), DIMENSION(:), &
         POINTER                               :: mos_aux_fit => NULL()
      TYPE(neighbor_list_set_p_type), &
         DIMENSION(:), POINTER                 :: sab_aux_fit => NULL(), sab_aux_fit_asymm => NULL(), sab_aux_fit_vs_orb => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), &
         POINTER                               :: matrix_ks_aux_fit => NULL(), matrix_ks_aux_fit_im => NULL(), &
                                                  matrix_ks_aux_fit_dft => NULL(), matrix_ks_aux_fit_hfx => NULL(), &
                                                  matrix_s_aux_fit => NULL(), matrix_s_aux_fit_vs_orb => NULL()
      TYPE(qs_rho_type), POINTER               :: rho_aux_fit => NULL(), rho_aux_fit_buffer => NULL()
      TYPE(task_list_type), POINTER            :: task_list_aux_fit => NULL()
      TYPE(cp_fm_type), DIMENSION(:), &
         POINTER                               :: mo_derivs_aux_fit => NULL()

   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief creates ADMM environment, initializes the basic types
!>
!> \param admm_env The ADMM env
!> \param admm_control ...
!> \param mos the MO's of the orbital basis set
!> \param para_env The parallel env
!> \param natoms ...
!> \param nao_aux_fit ...
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE admm_env_create(admm_env, admm_control, mos, para_env, natoms, nao_aux_fit)

      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(admm_control_type), POINTER                   :: admm_control
      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: natoms, nao_aux_fit

      INTEGER                                            :: i, iatom, iblock, ispin, j, jatom, &
                                                            nao_orb, nmo, nspins
      TYPE(cp_fm_struct_type), POINTER :: fm_struct_aux_aux, fm_struct_aux_nmo, fm_struct_aux_orb, &
         fm_struct_nmo_nmo, fm_struct_orb_aux, fm_struct_orb_nmo, fm_struct_orb_orb
      TYPE(cp_fm_type), POINTER                          :: mo_coeff

      CALL cite_reference(Guidon2010)

      ALLOCATE (admm_env)

      nspins = SIZE(mos)
      CALL get_mo_set(mos(1), mo_coeff=mo_coeff, nmo=nmo, nao=nao_orb)
      admm_env%nmo = 0
      admm_env%nao_aux_fit = nao_aux_fit
      admm_env%nao_orb = nao_orb
      CALL cp_fm_struct_create(fm_struct_aux_aux, &
                               context=mo_coeff%matrix_struct%context, &
                               nrow_global=nao_aux_fit, &
                               ncol_global=nao_aux_fit, &
                               para_env=para_env)
      CALL cp_fm_struct_create(fm_struct_aux_orb, &
                               context=mo_coeff%matrix_struct%context, &
                               nrow_global=nao_aux_fit, &
                               ncol_global=nao_orb, &
                               para_env=para_env)
      CALL cp_fm_struct_create(fm_struct_orb_aux, &
                               context=mo_coeff%matrix_struct%context, &
                               nrow_global=nao_orb, &
                               ncol_global=nao_aux_fit, &
                               para_env=para_env)
      CALL cp_fm_struct_create(fm_struct_orb_orb, &
                               context=mo_coeff%matrix_struct%context, &
                               nrow_global=nao_orb, &
                               ncol_global=nao_orb, &
                               para_env=para_env)

      NULLIFY (admm_env%S, admm_env%S_inv, admm_env%Q, admm_env%A, admm_env%B, &
               admm_env%work_orb_orb, admm_env%work_orb_orb2, admm_env%work_orb_orb3, &
               admm_env%work_aux_orb, admm_env%work_aux_orb2, admm_env%work_aux_orb3, &
               admm_env%work_aux_aux, admm_env%work_aux_aux2, admm_env%work_aux_aux3, &
               admm_env%work_aux_aux4, admm_env%work_aux_aux5)
      ALLOCATE (admm_env%S, admm_env%S_inv, admm_env%Q, admm_env%A, admm_env%B, &
                admm_env%work_orb_orb, admm_env%work_orb_orb2, admm_env%work_orb_orb3, &
                admm_env%work_aux_orb, admm_env%work_aux_orb2, admm_env%work_aux_orb3, &
                admm_env%work_aux_aux, admm_env%work_aux_aux2, admm_env%work_aux_aux3, &
                admm_env%work_aux_aux4, admm_env%work_aux_aux5)
      CALL cp_fm_create(admm_env%S, fm_struct_aux_aux, name="aux_fit_overlap")
      CALL cp_fm_create(admm_env%S_inv, fm_struct_aux_aux, name="aux_fit_overlap_inv")
      CALL cp_fm_create(admm_env%Q, fm_struct_aux_orb, name="mixed_overlap")
      CALL cp_fm_create(admm_env%A, fm_struct_aux_orb, name="work_A")
      CALL cp_fm_create(admm_env%B, fm_struct_orb_orb, name="work_B")
      CALL cp_fm_create(admm_env%work_orb_orb, fm_struct_orb_orb, name="work_orb_orb")
      CALL cp_fm_create(admm_env%work_orb_orb2, fm_struct_orb_orb, name="work_orb_orb")
      CALL cp_fm_create(admm_env%work_orb_orb3, fm_struct_orb_orb, name="work_orb_orb3")
      CALL cp_fm_create(admm_env%work_aux_orb, fm_struct_aux_orb, name="work_aux_orb")
      CALL cp_fm_create(admm_env%work_aux_orb2, fm_struct_aux_orb, name="work_aux_orb2")
      CALL cp_fm_create(admm_env%work_aux_orb3, fm_struct_aux_orb, name="work_aux_orb3")
      CALL cp_fm_create(admm_env%work_aux_aux, fm_struct_aux_aux, name="work_aux_aux")
      CALL cp_fm_create(admm_env%work_aux_aux2, fm_struct_aux_aux, name="work_aux_aux2")
      CALL cp_fm_create(admm_env%work_aux_aux3, fm_struct_aux_aux, name="work_aux_aux3")
      CALL cp_fm_create(admm_env%work_aux_aux4, fm_struct_aux_aux, name="work_aux_aux4")
      CALL cp_fm_create(admm_env%work_aux_aux5, fm_struct_aux_aux, name="work_aux_aux5")

      ALLOCATE (admm_env%lambda_inv(nspins))
      ALLOCATE (admm_env%lambda(nspins))
      ALLOCATE (admm_env%lambda_inv_sqrt(nspins))
      ALLOCATE (admm_env%R(nspins))
      ALLOCATE (admm_env%R_purify(nspins))
      ALLOCATE (admm_env%work_orb_nmo(nspins))
      ALLOCATE (admm_env%work_nmo_nmo1(nspins))
      ALLOCATE (admm_env%R_schur_R_t(nspins))
      ALLOCATE (admm_env%work_nmo_nmo2(nspins))
      ALLOCATE (admm_env%eigvals_lambda(nspins))
      ALLOCATE (admm_env%eigvals_P_to_be_purified(nspins))
      ALLOCATE (admm_env%H(nspins))
      ALLOCATE (admm_env%K(nspins))
      ALLOCATE (admm_env%M(nspins))
      ALLOCATE (admm_env%M_purify(nspins))
      ALLOCATE (admm_env%P_to_be_purified(nspins))
      ALLOCATE (admm_env%work_aux_nmo(nspins))
      ALLOCATE (admm_env%work_aux_nmo2(nspins))
      ALLOCATE (admm_env%mo_derivs_tmp(nspins))
      ALLOCATE (admm_env%H_corr(nspins))
      ALLOCATE (admm_env%ks_to_be_merged(nspins))
      ALLOCATE (admm_env%lambda_inv2(nspins))
      ALLOCATE (admm_env%C_hat(nspins))
      ALLOCATE (admm_env%P_tilde(nspins))

      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, nmo=nmo)
         admm_env%nmo(ispin) = nmo
         CALL cp_fm_struct_create(fm_struct_aux_nmo, &
                                  context=mo_coeff%matrix_struct%context, &
                                  nrow_global=nao_aux_fit, &
                                  ncol_global=nmo, &
                                  para_env=para_env)
         CALL cp_fm_struct_create(fm_struct_orb_nmo, &
                                  context=mo_coeff%matrix_struct%context, &
                                  nrow_global=nao_orb, &
                                  ncol_global=nmo, &
                                  para_env=para_env)
         CALL cp_fm_struct_create(fm_struct_nmo_nmo, &
                                  context=mo_coeff%matrix_struct%context, &
                                  nrow_global=nmo, &
                                  ncol_global=nmo, &
                                  para_env=para_env)

         CALL cp_fm_create(admm_env%work_orb_nmo(ispin), fm_struct_orb_nmo, name="work_orb_nmo")
         CALL cp_fm_create(admm_env%work_nmo_nmo1(ispin), fm_struct_nmo_nmo, name="work_nmo_nmo1")
         CALL cp_fm_create(admm_env%R_schur_R_t(ispin), fm_struct_nmo_nmo, name="R_schur_R_t")
         CALL cp_fm_create(admm_env%work_nmo_nmo2(ispin), fm_struct_nmo_nmo, name="work_nmo_nmo2")
         CALL cp_fm_create(admm_env%lambda(ispin), fm_struct_nmo_nmo, name="lambda")
         CALL cp_fm_create(admm_env%lambda_inv(ispin), fm_struct_nmo_nmo, name="lambda_inv")
         CALL cp_fm_create(admm_env%lambda_inv_sqrt(ispin), fm_struct_nmo_nmo, name="lambda_inv_sqrt")
         CALL cp_fm_create(admm_env%R(ispin), fm_struct_nmo_nmo, name="R")
         CALL cp_fm_create(admm_env%R_purify(ispin), fm_struct_aux_aux, name="R_purify")
         CALL cp_fm_create(admm_env%K(ispin), fm_struct_aux_aux, name="K")
         CALL cp_fm_create(admm_env%H(ispin), fm_struct_aux_nmo, name="H")
         CALL cp_fm_create(admm_env%H_corr(ispin), fm_struct_orb_orb, name="H_corr")
         CALL cp_fm_create(admm_env%M(ispin), fm_struct_nmo_nmo, name="M")
         CALL cp_fm_create(admm_env%M_purify(ispin), fm_struct_aux_aux, name="M aux")
         CALL cp_fm_create(admm_env%P_to_be_purified(ispin), fm_struct_aux_aux, name="P_to_be_purified")
         CALL cp_fm_create(admm_env%work_aux_nmo(ispin), fm_struct_aux_nmo, name="work_aux_nmo")
         CALL cp_fm_create(admm_env%work_aux_nmo2(ispin), fm_struct_aux_nmo, name="work_aux_nmo2")
         CALL cp_fm_create(admm_env%mo_derivs_tmp(ispin), fm_struct_orb_nmo, name="mo_derivs_tmp")
         CALL cp_fm_create(admm_env%lambda_inv2(ispin), fm_struct_nmo_nmo, name="lambda_inv2")
         CALL cp_fm_create(admm_env%C_hat(ispin), fm_struct_aux_nmo, name="C_hat")
         CALL cp_fm_create(admm_env%P_tilde(ispin), fm_struct_aux_aux, name="P_tilde")
         CALL cp_fm_create(admm_env%ks_to_be_merged(ispin), fm_struct_orb_orb, name="KS_to_be_merged ")

         ALLOCATE (admm_env%eigvals_lambda(ispin)%eigvals)
         ALLOCATE (admm_env%eigvals_P_to_be_purified(ispin)%eigvals)
         ALLOCATE (admm_env%eigvals_lambda(ispin)%eigvals%data(nmo))
         ALLOCATE (admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(nao_aux_fit))
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

      ! Copy settings from admm_control
      CPASSERT(ASSOCIATED(admm_control))
      admm_env%purification_method = admm_control%purification_method
      admm_env%scaling_model = admm_control%scaling_model
      admm_env%aux_exch_func = admm_control%aux_exch_func
      admm_env%charge_constrain = (admm_control%method == do_admm_charge_constrained_projection)
      admm_env%block_dm = ((admm_control%method == do_admm_blocking_purify_full) .OR. &
                           (admm_control%method == do_admm_blocked_projection))
      admm_env%block_fit = admm_control%method == do_admm_blocked_projection
      admm_env%aux_exch_func_param = admm_control%aux_exch_func_param
      admm_env%aux_x_param(:) = admm_control%aux_x_param(:)

      !ADMMP, ADMMQ, ADMMS
      IF ((.NOT. admm_env%charge_constrain) .AND. (admm_env%scaling_model == do_admm_exch_scaling_merlot)) &
         admm_env%do_admmp = .TRUE.
      IF (admm_env%charge_constrain .AND. (admm_env%scaling_model == do_admm_exch_scaling_none)) &
         admm_env%do_admmq = .TRUE.
      IF (admm_env%charge_constrain .AND. (admm_env%scaling_model == do_admm_exch_scaling_merlot)) &
         admm_env%do_admms = .TRUE.

      IF ((admm_control%method == do_admm_blocking_purify_full) .OR. &
          (admm_control%method == do_admm_blocked_projection)) THEN
         ! Create block map
         ALLOCATE (admm_env%block_map(natoms, natoms))
         admm_env%block_map(:, :) = 0
         DO iblock = 1, SIZE(admm_control%blocks)
            DO i = 1, SIZE(admm_control%blocks(iblock)%list)
               iatom = admm_control%blocks(iblock)%list(i)
               DO j = 1, SIZE(admm_control%blocks(iblock)%list)
                  jatom = admm_control%blocks(iblock)%list(j)
                  admm_env%block_map(iatom, jatom) = 1
               END DO
            END DO
         END DO
      END IF

      NULLIFY (admm_env%admm_gapw_env)
      admm_env%do_gapw = .FALSE.

      NULLIFY (admm_env%mos_aux_fit, admm_env%sab_aux_fit, admm_env%sab_aux_fit_asymm, admm_env%sab_aux_fit_vs_orb)
      NULLIFY (admm_env%matrix_ks_aux_fit, admm_env%matrix_ks_aux_fit_im, admm_env%matrix_ks_aux_fit_dft)
      NULLIFY (admm_env%matrix_ks_aux_fit_hfx, admm_env%matrix_s_aux_fit, admm_env%matrix_s_aux_fit_vs_orb)
      NULLIFY (admm_env%rho_aux_fit, admm_env%rho_aux_fit_buffer, admm_env%task_list_aux_fit, admm_env%mo_derivs_aux_fit)

   END SUBROUTINE admm_env_create

! **************************************************************************************************
!> \brief releases the ADMM environment, cleans up all types
!>
!> \param admm_env The ADMM env
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE admm_env_release(admm_env)

      TYPE(admm_type), POINTER                           :: admm_env

      INTEGER                                            :: ispin

      CPASSERT(ASSOCIATED(admm_env))

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
      DEALLOCATE (admm_env%S, admm_env%S_inv, admm_env%Q, admm_env%A, admm_env%B, &
                  admm_env%work_orb_orb, admm_env%work_orb_orb2, admm_env%work_orb_orb3, &
                  admm_env%work_aux_orb, admm_env%work_aux_orb2, admm_env%work_aux_orb3, &
                  admm_env%work_aux_aux, admm_env%work_aux_aux2, admm_env%work_aux_aux3, &
                  admm_env%work_aux_aux4, admm_env%work_aux_aux5)
      CALL cp_fm_release(admm_env%lambda)
      CALL cp_fm_release(admm_env%lambda_inv)
      CALL cp_fm_release(admm_env%lambda_inv_sqrt)
      CALL cp_fm_release(admm_env%lambda_inv2)
      CALL cp_fm_release(admm_env%C_hat)
      CALL cp_fm_release(admm_env%P_tilde)
      CALL cp_fm_release(admm_env%R)
      CALL cp_fm_release(admm_env%R_purify)
      CALL cp_fm_release(admm_env%H)
      CALL cp_fm_release(admm_env%H_corr)
      CALL cp_fm_release(admm_env%K)
      CALL cp_fm_release(admm_env%M)
      CALL cp_fm_release(admm_env%M_purify)
      CALL cp_fm_release(admm_env%P_to_be_purified)
      CALL cp_fm_release(admm_env%work_orb_nmo)
      CALL cp_fm_release(admm_env%work_nmo_nmo1)
      CALL cp_fm_release(admm_env%R_schur_R_t)
      CALL cp_fm_release(admm_env%work_nmo_nmo2)
      CALL cp_fm_release(admm_env%work_aux_nmo)
      CALL cp_fm_release(admm_env%work_aux_nmo2)
      CALL cp_fm_release(admm_env%mo_derivs_tmp)
      CALL cp_fm_release(admm_env%ks_to_be_merged)
      CALL cp_fm_release(admm_env%lambda_inv2)
      DO ispin = 1, SIZE(admm_env%eigvals_lambda)
         DEALLOCATE (admm_env%eigvals_lambda(ispin)%eigvals%data)
         DEALLOCATE (admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data)
         DEALLOCATE (admm_env%eigvals_lambda(ispin)%eigvals)
         DEALLOCATE (admm_env%eigvals_P_to_be_purified(ispin)%eigvals)
      END DO
      DEALLOCATE (admm_env%eigvals_lambda)
      DEALLOCATE (admm_env%eigvals_P_to_be_purified)

      IF (ASSOCIATED(admm_env%block_map)) &
         DEALLOCATE (admm_env%block_map)

      IF (ASSOCIATED(admm_env%xc_section_primary)) &
         CALL section_vals_release(admm_env%xc_section_primary)
      IF (ASSOCIATED(admm_env%xc_section_aux)) &
         CALL section_vals_release(admm_env%xc_section_aux)

      IF (ASSOCIATED(admm_env%admm_gapw_env)) CALL admm_gapw_env_release(admm_env%admm_gapw_env)
      IF (ASSOCIATED(admm_env%admm_dm)) CALL admm_dm_release(admm_env%admm_dm)

      IF (ASSOCIATED(admm_env%mos_aux_fit)) THEN
         DO ispin = 1, SIZE(admm_env%mos_aux_fit)
            CALL deallocate_mo_set(admm_env%mos_aux_fit(ispin))
         END DO
         DEALLOCATE (admm_env%mos_aux_fit)
      END IF
      CALL cp_fm_release(admm_env%mo_derivs_aux_fit)

      IF (ASSOCIATED(admm_env%sab_aux_fit)) CALL release_neighbor_list_sets(admm_env%sab_aux_fit)
      IF (ASSOCIATED(admm_env%sab_aux_fit_vs_orb)) CALL release_neighbor_list_sets(admm_env%sab_aux_fit_vs_orb)
      IF (ASSOCIATED(admm_env%sab_aux_fit_asymm)) CALL release_neighbor_list_sets(admm_env%sab_aux_fit_asymm)

      IF (ASSOCIATED(admm_env%matrix_ks_aux_fit)) CALL dbcsr_deallocate_matrix_set(admm_env%matrix_ks_aux_fit)
      IF (ASSOCIATED(admm_env%matrix_ks_aux_fit_im)) CALL dbcsr_deallocate_matrix_set(admm_env%matrix_ks_aux_fit_im)
      IF (ASSOCIATED(admm_env%matrix_ks_aux_fit_dft)) CALL dbcsr_deallocate_matrix_set(admm_env%matrix_ks_aux_fit_dft)
      IF (ASSOCIATED(admm_env%matrix_ks_aux_fit_hfx)) CALL dbcsr_deallocate_matrix_set(admm_env%matrix_ks_aux_fit_hfx)
      IF (ASSOCIATED(admm_env%matrix_s_aux_fit)) CALL dbcsr_deallocate_matrix_set(admm_env%matrix_s_aux_fit)
      IF (ASSOCIATED(admm_env%matrix_s_aux_fit_vs_orb)) CALL dbcsr_deallocate_matrix_set(admm_env%matrix_s_aux_fit_vs_orb)

      IF (ASSOCIATED(admm_env%rho_aux_fit)) THEN
         CALL qs_rho_release(admm_env%rho_aux_fit)
         DEALLOCATE (admm_env%rho_aux_fit)
      END IF
      IF (ASSOCIATED(admm_env%rho_aux_fit_buffer)) THEN
         CALL qs_rho_release(admm_env%rho_aux_fit_buffer)
         DEALLOCATE (admm_env%rho_aux_fit_buffer)
      END IF

      IF (ASSOCIATED(admm_env%task_list_aux_fit)) CALL deallocate_task_list(admm_env%task_list_aux_fit)

      DEALLOCATE (admm_env)

   END SUBROUTINE admm_env_release

! **************************************************************************************************
!> \brief Release the ADMM GAPW stuff
!> \param admm_gapw_env ...
! **************************************************************************************************
   SUBROUTINE admm_gapw_env_release(admm_gapw_env)

      TYPE(admm_gapw_type), POINTER                      :: admm_gapw_env

      IF (ASSOCIATED(admm_gapw_env%admm_kind_set)) THEN
         CALL deallocate_qs_kind_set(admm_gapw_env%admm_kind_set)
      END IF

      IF (ASSOCIATED(admm_gapw_env%local_rho_set)) THEN
         CALL local_rho_set_release(admm_gapw_env%local_rho_set)
      END IF

      IF (ASSOCIATED(admm_gapw_env%task_list)) THEN
         CALL deallocate_task_list(admm_gapw_env%task_list)
      END IF

      IF (ASSOCIATED(admm_gapw_env%oce)) THEN
         CALL deallocate_oce_set(admm_gapw_env%oce)
      END IF

      DEALLOCATE (admm_gapw_env)

   END SUBROUTINE admm_gapw_env_release

! **************************************************************************************************
!> \brief Get routine for the ADMM env
!> \param admm_env ...
!> \param mo_derivs_aux_fit ...
!> \param mos_aux_fit ...
!> \param sab_aux_fit ...
!> \param sab_aux_fit_asymm ...
!> \param sab_aux_fit_vs_orb ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_aux_fit_vs_orb ...
!> \param task_list_aux_fit ...
!> \param matrix_ks_aux_fit ...
!> \param matrix_ks_aux_fit_im ...
!> \param matrix_ks_aux_fit_dft ...
!> \param matrix_ks_aux_fit_hfx ...
!> \param rho_aux_fit ...
!> \param rho_aux_fit_buffer ...
!> \param admm_dm ...
! **************************************************************************************************
   SUBROUTINE get_admm_env(admm_env, mo_derivs_aux_fit, mos_aux_fit, sab_aux_fit, sab_aux_fit_asymm, &
                           sab_aux_fit_vs_orb, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, &
                           task_list_aux_fit, matrix_ks_aux_fit, matrix_ks_aux_fit_im, &
                           matrix_ks_aux_fit_dft, matrix_ks_aux_fit_hfx, rho_aux_fit, &
                           rho_aux_fit_buffer, admm_dm)

      TYPE(admm_type), INTENT(IN), POINTER               :: admm_env
      TYPE(cp_fm_type), DIMENSION(:), OPTIONAL, POINTER  :: mo_derivs_aux_fit
      TYPE(mo_set_type), DIMENSION(:), OPTIONAL, POINTER :: mos_aux_fit
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_aux_fit, sab_aux_fit_asymm, &
                                                            sab_aux_fit_vs_orb
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_s_aux_fit, matrix_s_aux_fit_vs_orb
      TYPE(task_list_type), OPTIONAL, POINTER            :: task_list_aux_fit
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_ks_aux_fit, matrix_ks_aux_fit_im, &
                                                            matrix_ks_aux_fit_dft, &
                                                            matrix_ks_aux_fit_hfx
      TYPE(qs_rho_type), OPTIONAL, POINTER               :: rho_aux_fit, rho_aux_fit_buffer
      TYPE(admm_dm_type), OPTIONAL, POINTER              :: admm_dm

      CPASSERT(ASSOCIATED(admm_env))

      IF (PRESENT(mo_derivs_aux_fit)) mo_derivs_aux_fit => admm_env%mo_derivs_aux_fit
      IF (PRESENT(mos_aux_fit)) mos_aux_fit => admm_env%mos_aux_fit
      IF (PRESENT(sab_aux_fit)) sab_aux_fit => admm_env%sab_aux_fit
      IF (PRESENT(sab_aux_fit_asymm)) sab_aux_fit_asymm => admm_env%sab_aux_fit_asymm
      IF (PRESENT(sab_aux_fit_vs_orb)) sab_aux_fit_vs_orb => admm_env%sab_aux_fit_vs_orb
      IF (PRESENT(matrix_s_aux_fit)) matrix_s_aux_fit => admm_env%matrix_s_aux_fit
      IF (PRESENT(matrix_s_aux_fit_vs_orb)) matrix_s_aux_fit_vs_orb => admm_env%matrix_s_aux_fit_vs_orb
      IF (PRESENT(task_list_aux_fit)) task_list_aux_fit => admm_env%task_list_aux_fit
      IF (PRESENT(matrix_ks_aux_fit)) matrix_ks_aux_fit => admm_env%matrix_ks_aux_fit
      IF (PRESENT(matrix_ks_aux_fit_im)) matrix_ks_aux_fit_im => admm_env%matrix_ks_aux_fit_im
      IF (PRESENT(matrix_ks_aux_fit_dft)) matrix_ks_aux_fit_dft => admm_env%matrix_ks_aux_fit_dft
      IF (PRESENT(matrix_ks_aux_fit_hfx)) matrix_ks_aux_fit_hfx => admm_env%matrix_ks_aux_fit_hfx
      IF (PRESENT(rho_aux_fit)) rho_aux_fit => admm_env%rho_aux_fit
      IF (PRESENT(rho_aux_fit_buffer)) rho_aux_fit_buffer => admm_env%rho_aux_fit_buffer
      IF (PRESENT(admm_dm)) admm_dm => admm_env%admm_dm

   END SUBROUTINE get_admm_env

! **************************************************************************************************
!> \brief Set routine for the ADMM env
!> \param admm_env ...
!> \param mo_derivs_aux_fit ...
!> \param mos_aux_fit ...
!> \param sab_aux_fit ...
!> \param sab_aux_fit_asymm ...
!> \param sab_aux_fit_vs_orb ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_aux_fit_vs_orb ...
!> \param task_list_aux_fit ...
!> \param matrix_ks_aux_fit ...
!> \param matrix_ks_aux_fit_im ...
!> \param matrix_ks_aux_fit_dft ...
!> \param matrix_ks_aux_fit_hfx ...
!> \param rho_aux_fit ...
!> \param rho_aux_fit_buffer ...
!> \param admm_dm ...
! **************************************************************************************************
   SUBROUTINE set_admm_env(admm_env, mo_derivs_aux_fit, mos_aux_fit, sab_aux_fit, sab_aux_fit_asymm, &
                           sab_aux_fit_vs_orb, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, &
                           task_list_aux_fit, matrix_ks_aux_fit, matrix_ks_aux_fit_im, &
                           matrix_ks_aux_fit_dft, matrix_ks_aux_fit_hfx, rho_aux_fit, &
                           rho_aux_fit_buffer, admm_dm)

      TYPE(admm_type), INTENT(INOUT), POINTER            :: admm_env
      TYPE(cp_fm_type), DIMENSION(:), OPTIONAL, POINTER  :: mo_derivs_aux_fit
      TYPE(mo_set_type), DIMENSION(:), OPTIONAL, POINTER :: mos_aux_fit
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_aux_fit, sab_aux_fit_asymm, &
                                                            sab_aux_fit_vs_orb
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_s_aux_fit, matrix_s_aux_fit_vs_orb
      TYPE(task_list_type), OPTIONAL, POINTER            :: task_list_aux_fit
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_ks_aux_fit, matrix_ks_aux_fit_im, &
                                                            matrix_ks_aux_fit_dft, &
                                                            matrix_ks_aux_fit_hfx
      TYPE(qs_rho_type), OPTIONAL, POINTER               :: rho_aux_fit, rho_aux_fit_buffer
      TYPE(admm_dm_type), OPTIONAL, POINTER              :: admm_dm

      CPASSERT(ASSOCIATED(admm_env))

      IF (PRESENT(mo_derivs_aux_fit)) admm_env%mo_derivs_aux_fit => mo_derivs_aux_fit
      IF (PRESENT(mos_aux_fit)) admm_env%mos_aux_fit => mos_aux_fit
      IF (PRESENT(sab_aux_fit)) admm_env%sab_aux_fit => sab_aux_fit
      IF (PRESENT(sab_aux_fit_asymm)) admm_env%sab_aux_fit_asymm => sab_aux_fit_asymm
      IF (PRESENT(sab_aux_fit_vs_orb)) admm_env%sab_aux_fit_vs_orb => sab_aux_fit_vs_orb
      IF (PRESENT(matrix_s_aux_fit)) admm_env%matrix_s_aux_fit => matrix_s_aux_fit
      IF (PRESENT(matrix_s_aux_fit_vs_orb)) admm_env%matrix_s_aux_fit_vs_orb => matrix_s_aux_fit_vs_orb
      IF (PRESENT(task_list_aux_fit)) admm_env%task_list_aux_fit => task_list_aux_fit
      IF (PRESENT(matrix_ks_aux_fit)) admm_env%matrix_ks_aux_fit => matrix_ks_aux_fit
      IF (PRESENT(matrix_ks_aux_fit_im)) admm_env%matrix_ks_aux_fit_im => matrix_ks_aux_fit_im
      IF (PRESENT(matrix_ks_aux_fit_dft)) admm_env%matrix_ks_aux_fit_dft => matrix_ks_aux_fit_dft
      IF (PRESENT(matrix_ks_aux_fit_hfx)) admm_env%matrix_ks_aux_fit_hfx => matrix_ks_aux_fit_hfx
      IF (PRESENT(rho_aux_fit)) admm_env%rho_aux_fit => rho_aux_fit
      IF (PRESENT(rho_aux_fit_buffer)) admm_env%rho_aux_fit_buffer => rho_aux_fit_buffer
      IF (PRESENT(admm_dm)) admm_env%admm_dm => admm_dm

   END SUBROUTINE set_admm_env

END MODULE admm_types

