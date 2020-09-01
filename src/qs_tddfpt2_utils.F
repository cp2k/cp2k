!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_utils
   USE cell_types,                      ONLY: cell_type
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type,&
                                              get_blacs_info
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_triangular_invert
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_type,&
                                              fm_pool_create,&
                                              fm_pool_release
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_p_type,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_get_info,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE input_constants,                 ONLY: &
        cholesky_dbcsr, cholesky_inverse, cholesky_off, cholesky_restore, oe_gllb, oe_lb, oe_none, &
        oe_saop, oe_shift
   USE input_section_types,             ONLY: section_vals_create,&
                                              section_vals_get_subs_vals,&
                                              section_vals_release,&
                                              section_vals_retain,&
                                              section_vals_set_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_min,&
                                              mp_sum
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_scf_methods,                  ONLY: eigensolver
   USE qs_scf_post_gpw,                 ONLY: make_lumo_gpw
   USE qs_scf_types,                    ONLY: ot_method_nr,&
                                              qs_scf_env_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE util,                            ONLY: sort
   USE xc_pot_saop,                     ONLY: add_saop_pot
   USE xtb_matrices,                    ONLY: build_xtb_ks_matrix
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_utils'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: tddfpt_init_ground_state_mos, tddfpt_release_ground_state_mos
   PUBLIC :: tddfpt_guess_vectors, tddfpt_init_mos, tddfpt_oecorr
   PUBLIC :: tddfpt_total_number_of_states

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Prepare MOs for TDDFPT Calculations
!> \param qs_env  Quickstep environment
!> \param gs_mos  ...
! **************************************************************************************************
   SUBROUTINE tddfpt_init_mos(qs_env, gs_mos)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt_init_mos'

      INTEGER                                            :: handle, ispin, nmo_avail, nmo_occ, &
                                                            nmo_virt, nspins
      REAL(kind=dp), DIMENSION(:), POINTER               :: evals_virt_spin
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: evals_virt
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_virt
      TYPE(cp_fm_type), POINTER                          :: mos_virt_spin
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, blacs_env=blacs_env, cell=cell, dft_control=dft_control, &
                      matrix_ks=matrix_ks, matrix_s=matrix_s, mos=mos, scf_env=scf_env)
      tddfpt_control => dft_control%tddfpt2_control

      CPASSERT(.NOT. ASSOCIATED(gs_mos))
      ! obtain occupied and virtual (unoccupied) ground-state Kohn-Sham orbitals
      nspins = dft_control%nspins
      ALLOCATE (gs_mos(nspins))

      ! when the number of unoccupied orbitals is limited and OT has been used
      ! for the ground-state DFT calculation,
      ! compute the missing unoccupied orbitals using OT as well.
      NULLIFY (evals_virt, evals_virt_spin, mos_virt, mos_virt_spin)
      IF (ASSOCIATED(scf_env)) THEN
         IF (scf_env%method == ot_method_nr .AND. tddfpt_control%nlumo > 0) THEN
            ! As OT with ADDED_MOS/=0 is currently not implemented, the following block is equivalent to:
            ! nmo_virt = tddfpt_control%nlumo
            ! number of already computed unoccupied orbitals (added_mos) .
            nmo_virt = HUGE(0)
            DO ispin = 1, nspins
               CALL get_mo_set(mos(ispin)%mo_set, nmo=nmo_avail, homo=nmo_occ)
               nmo_virt = MIN(nmo_virt, nmo_avail - nmo_occ)
            END DO
            ! number of unoccupied orbitals to compute
            nmo_virt = tddfpt_control%nlumo - nmo_virt
            IF (nmo_virt > 0) THEN
               ALLOCATE (evals_virt(nspins), mos_virt(nspins))
               ! the number of actually computed unoccupied orbitals will be stored as 'nmo_avail'
               CALL make_lumo_gpw(qs_env, scf_env, mos_virt, evals_virt, nmo_virt, nmo_avail)
            END IF
         END IF
      END IF

      DO ispin = 1, nspins
         IF (ASSOCIATED(evals_virt)) THEN
            evals_virt_spin => evals_virt(ispin)%array
         ELSE
            NULLIFY (evals_virt_spin)
         END IF
         IF (ASSOCIATED(mos_virt)) THEN
            mos_virt_spin => mos_virt(ispin)%matrix
         ELSE
            NULLIFY (mos_virt_spin)
         END IF
         CALL tddfpt_init_ground_state_mos(gs_mos=gs_mos(ispin), mo_set=mos(ispin)%mo_set, &
                                           nlumo=tddfpt_control%nlumo, &
                                           blacs_env=blacs_env, cholesky_method=cholesky_restore, &
                                           matrix_ks=matrix_ks(ispin)%matrix, matrix_s=matrix_s(1)%matrix, &
                                           mos_virt=mos_virt_spin, evals_virt=evals_virt_spin)
      END DO

      IF (ASSOCIATED(evals_virt)) THEN
         DO ispin = 1, SIZE(evals_virt)
            IF (ASSOCIATED(evals_virt(ispin)%array)) DEALLOCATE (evals_virt(ispin)%array)
         END DO
         DEALLOCATE (evals_virt)
      END IF

      IF (ASSOCIATED(mos_virt)) THEN
         DO ispin = 1, SIZE(mos_virt)
            IF (ASSOCIATED(mos_virt(ispin)%matrix)) CALL cp_fm_release(mos_virt(ispin)%matrix)
         END DO
         DEALLOCATE (mos_virt)
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_init_mos

! **************************************************************************************************
!> \brief Generate all virtual molecular orbitals for a given spin by diagonalising
!>        the corresponding Kohn-Sham matrix.
!> \param gs_mos           structure to store occupied and virtual molecular orbitals
!>                         (allocated and initialised on exit)
!> \param mo_set           ground state molecular orbitals for a given spin
!> \param nlumo            number of unoccupied states to consider (-1 means all states)
!> \param blacs_env        BLACS parallel environment
!> \param cholesky_method  Cholesky method to compute the inverse overlap matrix
!> \param matrix_ks        Kohn-Sham matrix for a given spin
!> \param matrix_s         overlap matrix
!> \param mos_virt         precomputed (OT) expansion coefficients of virtual molecular orbitals
!>                         (in addition to the ADDED_MOS, if present). NULL when no OT is in use.
!> \param evals_virt       orbital energies of precomputed (OT) virtual molecular orbitals.
!>                         NULL when no OT is in use.
!> \par History
!>    * 05.2016 created as tddfpt_lumos() [Sergey Chulkov]
!>    * 06.2016 renamed, altered prototype [Sergey Chulkov]
!>    * 04.2019 limit the number of unoccupied states, orbital energy correction [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_init_ground_state_mos(gs_mos, mo_set, nlumo, blacs_env, cholesky_method, matrix_ks, matrix_s, &
                                           mos_virt, evals_virt)
      TYPE(tddfpt_ground_state_mos)                      :: gs_mos
      TYPE(mo_set_type), POINTER                         :: mo_set
      INTEGER, INTENT(in)                                :: nlumo
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      INTEGER, INTENT(in)                                :: cholesky_method
      TYPE(dbcsr_type), POINTER                          :: matrix_ks, matrix_s
      TYPE(cp_fm_type), POINTER                          :: mos_virt
      REAL(kind=dp), DIMENSION(:), POINTER               :: evals_virt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init_ground_state_mos'
      REAL(kind=dp), PARAMETER                           :: eps_dp = EPSILON(0.0_dp)

      INTEGER :: cholesky_method_inout, handle, icol_global, icol_local, imo, irow_global, &
         irow_local, nao, ncol_local, nelectrons, nmo_occ, nmo_scf, nmo_virt, nrow_local, sign_int
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: minrow_neg_array, minrow_pos_array, &
                                                            sum_sign_array
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: do_eigen
      REAL(kind=dp)                                      :: element, maxocc
      REAL(kind=dp), DIMENSION(:), POINTER               :: mo_evals_extended, mo_occ_extended, &
                                                            mo_occ_scf
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: my_block
      TYPE(cp_fm_pool_type), POINTER                     :: wfn_fm_pool
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fm_struct, ao_mo_occ_fm_struct, &
                                                            ao_mo_virt_fm_struct, wfn_fm_struct
      TYPE(cp_fm_type), POINTER                          :: matrix_ks_fm, mo_coeff_extended, &
                                                            ortho_fm, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mo_set_type), POINTER                         :: mos_extended

      CALL timeset(routineN, handle)

      CALL get_blacs_info(blacs_env, para_env=para_env)

      CALL get_mo_set(mo_set, nao=nao, nmo=nmo_scf, homo=nmo_occ, maxocc=maxocc, &
                      nelectron=nelectrons, occupation_numbers=mo_occ_scf)

      nmo_virt = nao - nmo_occ
      IF (nlumo >= 0) &
         nmo_virt = MIN(nmo_virt, nlumo)

      IF (nmo_virt <= 0) &
         CALL cp_abort(__LOCATION__, &
                       'At least one unoccupied molecular orbital is required to calculate excited states.')

      do_eigen = .FALSE.
      ! diagonalise the Kohn-Sham matrix one more time if the number of available unoccupied states are too small
      IF (ASSOCIATED(evals_virt)) THEN
         CPASSERT(ASSOCIATED(mos_virt))
         IF (nmo_virt > nmo_scf - nmo_occ + SIZE(evals_virt)) do_eigen = .TRUE.
      ELSE
         IF (nmo_virt > nmo_scf - nmo_occ) do_eigen = .TRUE.
      END IF

      ! ++ allocate storage space for gs_mos
      NULLIFY (ao_mo_occ_fm_struct, ao_mo_virt_fm_struct)
      CALL cp_fm_struct_create(ao_mo_occ_fm_struct, nrow_global=nao, ncol_global=nmo_occ, context=blacs_env)
      CALL cp_fm_struct_create(ao_mo_virt_fm_struct, nrow_global=nao, ncol_global=nmo_virt, context=blacs_env)

      NULLIFY (gs_mos%mos_occ, gs_mos%mos_virt, gs_mos%evals_occ_matrix)
      CALL cp_fm_create(gs_mos%mos_occ, ao_mo_occ_fm_struct)
      CALL cp_fm_create(gs_mos%mos_virt, ao_mo_virt_fm_struct)

      ALLOCATE (gs_mos%evals_occ(nmo_occ))
      ALLOCATE (gs_mos%evals_virt(nmo_virt))
      ALLOCATE (gs_mos%phases_occ(nmo_occ))

      CALL cp_fm_struct_release(ao_mo_virt_fm_struct)
      !CALL cp_fm_struct_release(ao_mo_occ_fm_struct)

      ! ++ nullify pointers
      NULLIFY (ao_ao_fm_struct, wfn_fm_struct, wfn_fm_pool)
      NULLIFY (mos_extended, mo_coeff_extended, mo_evals_extended, mo_occ_extended)
      CALL cp_fm_struct_create(ao_ao_fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)

      IF (do_eigen) THEN
         ! ++ set of molecular orbitals
         CALL cp_fm_struct_create(wfn_fm_struct, nrow_global=nao, ncol_global=nmo_occ + nmo_virt, context=blacs_env)
         CALL fm_pool_create(wfn_fm_pool, wfn_fm_struct)

         CALL allocate_mo_set(mos_extended, nao, nmo_occ + nmo_virt, nelectrons, &
                              REAL(nelectrons, dp), maxocc, flexible_electron_count=0.0_dp)
         CALL init_mo_set(mos_extended, fm_pool=wfn_fm_pool, name="mos-extended")
         CALL fm_pool_release(wfn_fm_pool)
         CALL cp_fm_struct_release(wfn_fm_struct)
         CALL get_mo_set(mos_extended, mo_coeff=mo_coeff_extended, &
                         eigenvalues=mo_evals_extended, occupation_numbers=mo_occ_extended)

         ! use the explicit loop in order to avoid temporary arrays.
         !
         ! The assignment statement : mo_occ_extended(1:nmo_scf) = mo_occ_scf(1:nmo_scf)
         ! implies temporary arrays as a compiler does not know in advance that the pointers
         ! on both sides of the statement point to non-overlapped memory regions
         DO imo = 1, nmo_scf
            mo_occ_extended(imo) = mo_occ_scf(imo)
         END DO
         mo_occ_extended(nmo_scf + 1:) = 0.0_dp

         ! ++ allocate temporary matrices
         NULLIFY (matrix_ks_fm, ortho_fm, work_fm)
         CALL cp_fm_create(matrix_ks_fm, ao_ao_fm_struct)
         CALL cp_fm_create(ortho_fm, ao_ao_fm_struct)
         CALL cp_fm_create(work_fm, ao_ao_fm_struct)
         CALL cp_fm_struct_release(ao_ao_fm_struct)

         ! some stuff from the subroutine general_eigenproblem()
         CALL copy_dbcsr_to_fm(matrix_s, ortho_fm)
         CALL copy_dbcsr_to_fm(matrix_ks, matrix_ks_fm)

         IF (cholesky_method == cholesky_dbcsr) THEN
            CPABORT('CHOLESKY DBCSR_INVERSE is not implemented in TDDFT.')
         ELSE IF (cholesky_method == cholesky_off) THEN
            CPABORT('CHOLESKY OFF is not implemented in TDDFT.')
         ELSE
            CALL cp_fm_cholesky_decompose(ortho_fm)
            IF (cholesky_method == cholesky_inverse) THEN
               CALL cp_fm_triangular_invert(ortho_fm)
            END IF

            ! need to store 'cholesky_method' in a temporary variable, as the subroutine eigensolver()
            ! will update this variable
            cholesky_method_inout = cholesky_method
            CALL eigensolver(matrix_ks_fm=matrix_ks_fm, mo_set=mos_extended, ortho=ortho_fm, &
                             work=work_fm, cholesky_method=cholesky_method_inout, &
                             do_level_shift=.FALSE., level_shift=0.0_dp, matrix_u_fm=null(), use_jacobi=.FALSE.)
         END IF

         ! -- clean up needless matrices
         CALL cp_fm_release(work_fm)
         CALL cp_fm_release(ortho_fm)
         CALL cp_fm_release(matrix_ks_fm)
      ELSE
         CALL get_mo_set(mo_set, mo_coeff=mo_coeff_extended, &
                         eigenvalues=mo_evals_extended, occupation_numbers=mo_occ_extended)
      END IF

      ! compute the phase of molecular orbitals;
      ! matrix work_fm holds occupied molecular orbital coefficients distributed among all the processors
      !CALL cp_fm_struct_create(ao_mo_occ_fm_struct, nrow_global=nao, ncol_global=nmo_occ, context=blacs_env)
      CALL cp_fm_create(work_fm, ao_mo_occ_fm_struct)
      CALL cp_fm_struct_release(ao_mo_occ_fm_struct)

      CALL cp_fm_to_fm(mo_coeff_extended, work_fm, ncol=nmo_occ, source_start=1, target_start=1)
      CALL cp_fm_get_info(work_fm, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, local_data=my_block)

      ALLOCATE (minrow_neg_array(nmo_occ), minrow_pos_array(nmo_occ), sum_sign_array(nmo_occ))
      minrow_neg_array(:) = nao
      minrow_pos_array(:) = nao
      sum_sign_array(:) = 0
      DO icol_local = 1, ncol_local
         icol_global = col_indices(icol_local)

         DO irow_local = 1, nrow_local
            element = my_block(irow_local, icol_local)

            sign_int = 0
            IF (element >= eps_dp) THEN
               sign_int = 1
            ELSE IF (element <= -eps_dp) THEN
               sign_int = -1
            END IF

            sum_sign_array(icol_global) = sum_sign_array(icol_global) + sign_int

            irow_global = row_indices(irow_local)
            IF (sign_int > 0) THEN
               IF (minrow_pos_array(icol_global) > irow_global) &
                  minrow_pos_array(icol_global) = irow_global
            ELSE IF (sign_int < 0) THEN
               IF (minrow_neg_array(icol_global) > irow_global) &
                  minrow_neg_array(icol_global) = irow_global
            END IF
         END DO
      END DO

      CALL mp_sum(sum_sign_array, para_env%group)
      CALL mp_min(minrow_neg_array, para_env%group)
      CALL mp_min(minrow_pos_array, para_env%group)

      DO icol_local = 1, nmo_occ
         IF (sum_sign_array(icol_local) > 0) THEN
            ! most of the expansion coefficients are positive => MO's phase = +1
            gs_mos%phases_occ(icol_local) = 1.0_dp
         ELSE IF (sum_sign_array(icol_local) < 0) THEN
            ! most of the expansion coefficients are negative => MO's phase = -1
            gs_mos%phases_occ(icol_local) = -1.0_dp
         ELSE
            ! equal number of positive and negative expansion coefficients
            IF (minrow_pos_array(icol_local) <= minrow_neg_array(icol_local)) THEN
               ! the first positive expansion coefficient has a lower index then
               ! the first negative expansion coefficient; MO's phase = +1
               gs_mos%phases_occ(icol_local) = 1.0_dp
            ELSE
               ! MO's phase = -1
               gs_mos%phases_occ(icol_local) = -1.0_dp
            END IF
         END IF
      END DO

      DEALLOCATE (minrow_neg_array, minrow_pos_array, sum_sign_array)

      ! return the requested occupied and virtual molecular orbitals and corresponding orbital energies
      CALL cp_fm_to_fm(mo_coeff_extended, gs_mos%mos_occ, ncol=nmo_occ, source_start=1, target_start=1)
      gs_mos%evals_occ(1:nmo_occ) = mo_evals_extended(1:nmo_occ)

      IF (ASSOCIATED(evals_virt) .AND. (.NOT. do_eigen) .AND. nmo_virt > nmo_scf - nmo_occ) THEN
         CALL cp_fm_to_fm(mo_coeff_extended, gs_mos%mos_virt, ncol=nmo_scf - nmo_occ, &
                          source_start=nmo_occ + 1, target_start=1)
         CALL cp_fm_to_fm(mos_virt, gs_mos%mos_virt, ncol=nmo_virt - (nmo_scf - nmo_occ), &
                          source_start=1, target_start=nmo_scf - nmo_occ + 1)
         gs_mos%evals_virt(1:nmo_scf - nmo_occ) = evals_virt(nmo_occ + 1:nmo_occ + nmo_scf)
         gs_mos%evals_virt(nmo_scf - nmo_occ + 1:nmo_virt) = evals_virt(1:nmo_virt - (nmo_scf - nmo_occ))
      ELSE
         CALL cp_fm_to_fm(mo_coeff_extended, gs_mos%mos_virt, ncol=nmo_virt, source_start=nmo_occ + 1, target_start=1)
         gs_mos%evals_virt(1:nmo_virt) = mo_evals_extended(nmo_occ + 1:nmo_occ + nmo_virt)
      END IF

      CALL cp_fm_release(work_fm)

      IF (do_eigen) CALL deallocate_mo_set(mos_extended)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_init_ground_state_mos

! **************************************************************************************************
!> \brief Release molecular orbitals.
!> \param gs_mos  structure that holds occupied and virtual molecular orbitals
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_release_ground_state_mos(gs_mos)
      TYPE(tddfpt_ground_state_mos), INTENT(inout)       :: gs_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_release_ground_state_mos'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ALLOCATED(gs_mos%phases_occ)) &
         DEALLOCATE (gs_mos%phases_occ)

      IF (ALLOCATED(gs_mos%evals_virt)) &
         DEALLOCATE (gs_mos%evals_virt)

      IF (ALLOCATED(gs_mos%evals_occ)) &
         DEALLOCATE (gs_mos%evals_occ)

      IF (ASSOCIATED(gs_mos%evals_occ_matrix)) &
         CALL cp_fm_release(gs_mos%evals_occ_matrix)

      IF (ASSOCIATED(gs_mos%mos_virt)) &
         CALL cp_fm_release(gs_mos%mos_virt)

      IF (ASSOCIATED(gs_mos%mos_occ)) &
         CALL cp_fm_release(gs_mos%mos_occ)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_release_ground_state_mos

! **************************************************************************************************
!> \brief Callculate orbital corrected KS matrix for TDDFPT
!> \param qs_env  Quickstep environment
!> \param gs_mos ...
!> \param matrix_ks_oep ...
!> \param do_hfx ...
! **************************************************************************************************
   SUBROUTINE tddfpt_oecorr(qs_env, gs_mos, matrix_ks_oep, do_hfx)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks_oep
      LOGICAL, INTENT(IN)                                :: do_hfx

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt_oecorr'

      INTEGER                                            :: handle, iounit, ispin, nao, nmo_occ, &
                                                            nspins
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: ao_mo_occ_fm_struct, &
                                                            mo_occ_mo_occ_fm_struct
      TYPE(cp_fm_type), POINTER                          :: work_fm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: xc_fun_empty, xc_fun_original
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env, blacs_env=blacs_env, dft_control=dft_control, matrix_ks=matrix_ks)
      tddfpt_control => dft_control%tddfpt2_control

      ! obtain corrected KS-matrix
      ! We should 'save' the energy values?
      nspins = SIZE(gs_mos)
      NULLIFY (matrix_ks_oep)
      IF (tddfpt_control%oe_corr /= oe_none) THEN
         IF (iounit > 0) THEN
            WRITE (iounit, "(1X,A)") "", &
               "-------------------------------------------------------------------------------", &
               "-                    Orbital Eigenvalue Correction Started                    -", &
               "-------------------------------------------------------------------------------"
         END IF

         CALL cp_warn(__LOCATION__, &
                      "Orbital energy correction potential is an experimental feature. "// &
                      "Use it with extreme care")

         IF (do_hfx) THEN
            CALL cp_abort(__LOCATION__, &
                          "Implementation of orbital energy correction XC-potentials is "// &
                          "currently incompatible with exact-exchange functionals")
         END IF

         CALL dbcsr_allocate_matrix_set(matrix_ks_oep, nspins)
         DO ispin = 1, nspins
            CALL dbcsr_init_p(matrix_ks_oep(ispin)%matrix)
            CALL dbcsr_copy(matrix_ks_oep(ispin)%matrix, matrix_ks(ispin)%matrix)
         END DO

         ! KS-matrix without XC-terms
         xc_fun_original => section_vals_get_subs_vals(qs_env%input, "DFT%XC%XC_FUNCTIONAL")
         CALL section_vals_retain(xc_fun_original)
         NULLIFY (xc_fun_empty)
         CALL section_vals_create(xc_fun_empty, xc_fun_original%section)
         CALL section_vals_set_subs_vals(qs_env%input, "DFT%XC%XC_FUNCTIONAL", xc_fun_empty)
         CALL section_vals_release(xc_fun_empty)

         IF (dft_control%qs_control%semi_empirical) THEN
            CPABORT("TDDFPT with SE not possible")
         ELSEIF (dft_control%qs_control%dftb) THEN
            CPABORT("TDDFPT with DFTB not possible")
         ELSEIF (dft_control%qs_control%xtb) THEN
            CALL build_xtb_ks_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE., &
                                     ext_ks_matrix=matrix_ks_oep)
         ELSE
            CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE., &
                                              ext_ks_matrix=matrix_ks_oep)
         END IF

         IF (tddfpt_control%oe_corr == oe_saop .OR. &
             tddfpt_control%oe_corr == oe_lb .OR. &
             tddfpt_control%oe_corr == oe_gllb) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(T2,A)") " Orbital energy correction of SAOP type "
            END IF
            CALL add_saop_pot(matrix_ks_oep, qs_env, tddfpt_control%oe_corr)
         ELSE IF (tddfpt_control%oe_corr == oe_shift) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(T2,A,T71,F10.3)") &
                  " Virtual Orbital Eigenvalue Shift [eV] ", tddfpt_control%ev_shift*evolt
               WRITE (iounit, "(T2,A,T71,F10.3)") &
                  " Open Shell Orbital Eigenvalue Shift [eV] ", tddfpt_control%eos_shift*evolt
            END IF
            CALL ev_shift_operator(qs_env, gs_mos, matrix_ks_oep, &
                                   tddfpt_control%ev_shift, tddfpt_control%eos_shift)
         ELSE
            CALL cp_abort(__LOCATION__, &
                          "Unimplemented orbital energy correction potential")
         END IF
         CALL section_vals_set_subs_vals(qs_env%input, "DFT%XC%XC_FUNCTIONAL", xc_fun_original)
         CALL section_vals_release(xc_fun_original)

         ! compute 'evals_occ_matrix'
         CALL dbcsr_get_info(matrix_ks(1)%matrix, nfullrows_total=nao)
         NULLIFY (mo_occ_mo_occ_fm_struct)
         DO ispin = 1, nspins
            nmo_occ = SIZE(gs_mos(ispin)%evals_occ)
            CALL cp_fm_struct_create(mo_occ_mo_occ_fm_struct, nrow_global=nmo_occ, ncol_global=nmo_occ, &
                                     context=blacs_env)
            CALL cp_fm_create(gs_mos(ispin)%evals_occ_matrix, mo_occ_mo_occ_fm_struct)
            CALL cp_fm_struct_release(mo_occ_mo_occ_fm_struct)
            ! work_fm is a temporary [nao x nmo_occ] matrix
            CALL cp_fm_struct_create(ao_mo_occ_fm_struct, nrow_global=nao, ncol_global=nmo_occ, &
                                     context=blacs_env)
            CALL cp_fm_create(work_fm, ao_mo_occ_fm_struct)
            CALL cp_fm_struct_release(ao_mo_occ_fm_struct)
            CALL cp_dbcsr_sm_fm_multiply(matrix_ks_oep(ispin)%matrix, gs_mos(ispin)%mos_occ, &
                                         work_fm, ncol=nmo_occ, alpha=1.0_dp, beta=0.0_dp)
            CALL cp_gemm('T', 'N', nmo_occ, nmo_occ, nao, 1.0_dp, gs_mos(ispin)%mos_occ, work_fm, &
                         0.0_dp, gs_mos(ispin)%evals_occ_matrix)
            CALL cp_fm_release(work_fm)
         END DO
         IF (iounit > 0) THEN
            WRITE (iounit, "(1X,A)") &
               "-------------------------------------------------------------------------------"
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_oecorr

! **************************************************************************************************
!> \brief Compute the number of possible singly excited states (occ -> virt)
!> \param gs_mos          occupied and virtual molecular orbitals optimised for the ground state
!> \return the number of possible single excitations
!> \par History
!>    * 01.2017 created [Sergey Chulkov]
! **************************************************************************************************
   PURE FUNCTION tddfpt_total_number_of_states(gs_mos) RESULT(nstates_total)
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      INTEGER(kind=int_8)                                :: nstates_total

      INTEGER                                            :: ispin, nspins

      nstates_total = 0
      nspins = SIZE(gs_mos)

      DO ispin = 1, nspins
         nstates_total = nstates_total + &
                         SIZE(gs_mos(ispin)%evals_occ, kind=int_8)* &
                         SIZE(gs_mos(ispin)%evals_virt, kind=int_8)
      END DO
   END FUNCTION tddfpt_total_number_of_states

! **************************************************************************************************
!> \brief Create a shift operator on virtual/open shell space
!>        Shift operator = Edelta*Q  Q: projector on virtual space (1-PS)
!>                                      projector on open shell space PosS
!> \param qs_env the qs_env that is perturbed by this p_env
!> \param gs_mos  ...
!> \param matrix_ks ...
!> \param ev_shift ...
!> \param eos_shift ...
!> \par History
!>      02.04.2019 adapted for TDDFT use from p_env (JGH)
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ev_shift_operator(qs_env, gs_mos, matrix_ks, ev_shift, eos_shift)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      REAL(KIND=dp), INTENT(IN)                          :: ev_shift, eos_shift

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ev_shift_operator'

      INTEGER                                            :: handle, ispin, n_spins, na, nb, nhomo, &
                                                            nl, nos, nrow, nu, nvirt
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), POINTER                          :: cmos, coeff, cvec
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: smat
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      CALL timeset(routineN, handle)

      n_spins = SIZE(gs_mos)
      CPASSERT(n_spins == SIZE(matrix_ks))

      IF (eos_shift /= 0.0_dp .AND. n_spins > 1) THEN
         CPABORT("eos_shift not implemented")
         CALL get_qs_env(qs_env, mos=mos, matrix_s=matrix_s)
         smat => matrix_s(1)%matrix
         CALL cp_fm_get_info(gs_mos(1)%mos_occ, ncol_global=na)
         CALL cp_fm_get_info(gs_mos(2)%mos_occ, ncol_global=nb)
         nl = MIN(na, nb)
         nu = MAX(na, nb)
         ! open shell orbital shift
         DO ispin = 1, n_spins
            coeff => gs_mos(ispin)%mos_occ
            CALL cp_fm_get_info(coeff, matrix_struct=fmstruct, ncol_global=nhomo)
            IF (nhomo == nu) THEN
               ! downshift with -eos_shift using occupied orbitals
               nos = nu - nl
               CALL cp_fm_create(cmos, fmstruct)
               CALL cp_fm_get_info(coeff, nrow_global=nrow)
               CALL cp_fm_to_fm_submat(coeff, cmos, nrow, nos, 1, nl + 1, 1, 1)
               CALL cp_fm_create(cvec, fmstruct)
               CALL cp_dbcsr_sm_fm_multiply(smat, cmos, cvec, nos, 1.0_dp, 0.0_dp)
               CALL cp_dbcsr_plus_fm_fm_t(matrix_ks(ispin)%matrix, matrix_v=cvec, ncol=nos, &
                                          alpha=-eos_shift, keep_sparsity=.TRUE.)
               CALL cp_fm_release(cmos)
               CALL cp_fm_release(cvec)
            ELSE
               ! upshift with eos_shift using virtual orbitals
               coeff => gs_mos(ispin)%mos_virt
               CALL cp_fm_get_info(coeff, matrix_struct=fmstruct, ncol_global=nvirt)
               nos = nu - nhomo
               CPASSERT(nvirt >= nos)
               CALL cp_fm_create(cvec, fmstruct)
               CALL cp_dbcsr_sm_fm_multiply(smat, coeff, cvec, nos, 1.0_dp, 0.0_dp)
               CALL cp_dbcsr_plus_fm_fm_t(matrix_ks(ispin)%matrix, matrix_v=cvec, ncol=nos, &
                                          alpha=eos_shift, keep_sparsity=.TRUE.)
               CALL cp_fm_release(cvec)
            END IF
         END DO
         ! virtual shift
         IF (ev_shift /= 0.0_dp) THEN
            DO ispin = 1, n_spins
               CALL dbcsr_add(matrix_ks(ispin)%matrix, smat, &
                              alpha_scalar=1.0_dp, beta_scalar=ev_shift)
               coeff => gs_mos(ispin)%mos_occ
               CALL cp_fm_get_info(coeff, matrix_struct=fmstruct, ncol_global=nhomo)
               CALL cp_fm_create(cvec, fmstruct)
               CALL cp_dbcsr_sm_fm_multiply(smat, coeff, cvec, nhomo, 1.0_dp, 0.0_dp)
               CALL cp_dbcsr_plus_fm_fm_t(matrix_ks(ispin)%matrix, matrix_v=cvec, ncol=nhomo, &
                                          alpha=-ev_shift, keep_sparsity=.TRUE.)
               CALL cp_fm_release(cvec)
               IF (nhomo < nu) THEN
                  nos = nu - nhomo
                  coeff => gs_mos(ispin)%mos_virt
                  CALL cp_fm_get_info(coeff, matrix_struct=fmstruct, ncol_global=nvirt)
                  CPASSERT(nvirt >= nos)
                  CALL cp_fm_create(cvec, fmstruct)
                  CALL cp_dbcsr_sm_fm_multiply(smat, coeff, cvec, nos, 1.0_dp, 0.0_dp)
                  CALL cp_dbcsr_plus_fm_fm_t(matrix_ks(ispin)%matrix, matrix_v=cvec, ncol=nos, &
                                             alpha=-ev_shift, keep_sparsity=.TRUE.)
                  CALL cp_fm_release(cvec)
               END IF
            END DO
         END IF
      ELSE
         ! virtual shift
         IF (ev_shift /= 0.0_dp) THEN
            CALL get_qs_env(qs_env, mos=mos, matrix_s=matrix_s)
            smat => matrix_s(1)%matrix
            DO ispin = 1, n_spins
               CALL dbcsr_add(matrix_ks(ispin)%matrix, smat, &
                              alpha_scalar=1.0_dp, beta_scalar=ev_shift)
               coeff => gs_mos(ispin)%mos_occ
               CALL cp_fm_get_info(coeff, matrix_struct=fmstruct, ncol_global=nhomo)
               CALL cp_fm_create(cvec, fmstruct)
               CALL cp_dbcsr_sm_fm_multiply(smat, coeff, cvec, nhomo, 1.0_dp, 0.0_dp)
               CALL cp_dbcsr_plus_fm_fm_t(matrix_ks(ispin)%matrix, matrix_v=cvec, ncol=nhomo, &
                                          alpha=-ev_shift, keep_sparsity=.TRUE.)
               CALL cp_fm_release(cvec)
            END DO
         END IF
      END IF
      ! set eigenvalues
      IF (eos_shift == 0.0_dp .OR. n_spins == 1) THEN
         DO ispin = 1, n_spins
            IF (ALLOCATED(gs_mos(ispin)%evals_virt)) THEN
               gs_mos(ispin)%evals_virt(:) = gs_mos(ispin)%evals_virt(:) + ev_shift
            END IF
         END DO
      ELSE
         CALL cp_fm_get_info(gs_mos(1)%mos_occ, ncol_global=na)
         CALL cp_fm_get_info(gs_mos(2)%mos_occ, ncol_global=nb)
         nl = MIN(na, nb)
         nu = MAX(na, nb)
         nos = nu - nl
         IF (na == nu) THEN
            IF (ALLOCATED(gs_mos(1)%evals_occ)) THEN
               gs_mos(1)%evals_occ(nl + 1:nu) = gs_mos(1)%evals_occ(nl + 1:nu) - eos_shift
            END IF
            IF (ALLOCATED(gs_mos(1)%evals_virt)) THEN
               gs_mos(1)%evals_virt(:) = gs_mos(1)%evals_virt(:) + ev_shift
            END IF
            IF (ALLOCATED(gs_mos(2)%evals_virt)) THEN
               gs_mos(2)%evals_virt(1:nos) = gs_mos(2)%evals_virt(1:nos) + eos_shift
               gs_mos(2)%evals_virt(nos + 1:) = gs_mos(2)%evals_virt(nos + 1:) + ev_shift
            END IF
         ELSE
            IF (ALLOCATED(gs_mos(1)%evals_virt)) THEN
               gs_mos(1)%evals_virt(1:nos) = gs_mos(1)%evals_virt(1:nos) + eos_shift
               gs_mos(1)%evals_virt(nos + 1:) = gs_mos(1)%evals_virt(nos + 1:) + ev_shift
            END IF
            IF (ALLOCATED(gs_mos(2)%evals_occ)) THEN
               gs_mos(2)%evals_occ(nl + 1:nu) = gs_mos(2)%evals_occ(nl + 1:nu) - eos_shift
            END IF
            IF (ALLOCATED(gs_mos(2)%evals_virt)) THEN
               gs_mos(2)%evals_virt(:) = gs_mos(2)%evals_virt(:) + ev_shift
            END IF
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE ev_shift_operator

! **************************************************************************************************
!> \brief Generate missed guess vectors.
!> \param evects   guess vectors distributed across all processors (initialised on exit)
!> \param evals    guessed transition energies (initialised on exit)
!> \param gs_mos   occupied and virtual molecular orbitals optimised for the ground state
!> \param log_unit output unit
!> \par History
!>    * 05.2016 created as tddfpt_guess() [Sergey Chulkov]
!>    * 06.2016 renamed, altered prototype, supports spin-polarised density [Sergey Chulkov]
!>    * 01.2017 simplified prototype, do not compute all possible singly-excited states
!>              [Sergey Chulkov]
!> \note \parblock
!>       Based on the subroutine co_initial_guess() which was originally created by
!>       Thomas Chassaing on 06.2003.
!>
!>       Only not associated guess vectors 'evects(spin, state)%matrix' are allocated and
!>       initialised; associated vectors assumed to be initialised elsewhere (e.g. using
!>       a restart file).
!>       \endparblock
! **************************************************************************************************
   SUBROUTINE tddfpt_guess_vectors(evects, evals, gs_mos, log_unit)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(inout) :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(inout)         :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      INTEGER, INTENT(in)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_guess_vectors'

      CHARACTER(len=5)                                   :: spin_label
      INTEGER                                            :: handle, imo_occ, imo_virt, ind, ispin, &
                                                            istate, jspin, nspins, nstates, &
                                                            nstates_occ_virt_alpha, &
                                                            nstates_selected
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ_avail, nmo_occ_selected, &
                                                            nmo_virt_selected
      REAL(kind=dp)                                      :: e_occ
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_minus_occ
      TYPE(cp_fm_struct_p_type), DIMENSION(maxspins)     :: fm_struct_evects

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nstates = SIZE(evects, 2)

      IF (debug_this_module) THEN
         CPASSERT(nstates > 0)
         CPASSERT(nspins == 1 .OR. nspins == 2)
         CPASSERT(SIZE(gs_mos) == nspins)
      END IF

      DO ispin = 1, nspins
         ! number of occupied orbitals for each spin component
         nmo_occ_avail(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         ! number of occupied and virtual orbitals which can potentially
         ! contribute to the excited states in question.
         nmo_occ_selected(ispin) = MIN(nmo_occ_avail(ispin), nstates)
         nmo_virt_selected(ispin) = MIN(SIZE(gs_mos(ispin)%evals_virt), nstates)

         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, matrix_struct=fm_struct_evects(ispin)%struct)
      END DO

      ! TO DO: the variable 'nstates_selected' should probably be declared as INTEGER(kind=int_8),
      !        however we need a special version of the subroutine sort() in order to do so
      nstates_selected = DOT_PRODUCT(nmo_occ_selected(1:nspins), nmo_virt_selected(1:nspins))

      ALLOCATE (inds(nstates_selected))
      ALLOCATE (e_virt_minus_occ(nstates_selected))

      istate = 0
      DO ispin = 1, nspins
         DO imo_occ = 1, nmo_occ_selected(ispin)
            ! Here imo_occ enumerate Occupied orbitals in inverse order (from the last to the first element)
            e_occ = gs_mos(ispin)%evals_occ(nmo_occ_avail(ispin) - imo_occ + 1)

            DO imo_virt = 1, nmo_virt_selected(ispin)
               istate = istate + 1
               e_virt_minus_occ(istate) = gs_mos(ispin)%evals_virt(imo_virt) - e_occ
            END DO
         END DO
      END DO

      IF (debug_this_module) THEN
         CPASSERT(istate == nstates_selected)
      END IF

      CALL sort(e_virt_minus_occ, nstates_selected, inds)

      IF (nspins == 1) THEN
         ispin = 1
         spin_label = '     '
      END IF

      nstates_occ_virt_alpha = nmo_occ_selected(1)*nmo_virt_selected(1)
      IF (log_unit > 0) THEN
         WRITE (log_unit, "(1X,A)") "", &
            "-------------------------------------------------------------------------------", &
            "-                            TDDFPT Initial Guess                             -", &
            "-------------------------------------------------------------------------------"
         WRITE (log_unit, '(T11,A)') "State         Occupied      ->      Virtual          Excitation"
         WRITE (log_unit, '(T11,A)') "number         orbital              orbital          energy (eV)"
         WRITE (log_unit, '(1X,79("-"))')
      END IF

      DO istate = 1, nstates
         IF (ASSOCIATED(evects(1, istate)%matrix)) THEN
            IF (log_unit > 0) &
               WRITE (log_unit, '(T7,I8,T28,A19,T60,F14.5)') &
               istate, "***  restarted  ***", evals(istate)*evolt
         ELSE
            ind = inds(istate) - 1
            IF (nspins > 1) THEN
               IF (ind < nstates_occ_virt_alpha) THEN
                  ispin = 1
                  spin_label = '(alp)'
               ELSE
                  ispin = 2
                  ind = ind - nstates_occ_virt_alpha
                  spin_label = '(bet)'
               END IF
            END IF

            imo_occ = nmo_occ_avail(ispin) - ind/nmo_virt_selected(ispin)
            imo_virt = MOD(ind, nmo_virt_selected(ispin)) + 1
            evals(istate) = e_virt_minus_occ(istate)

            IF (log_unit > 0) &
               WRITE (log_unit, '(T7,I8,T24,I8,T37,A5,T45,I8,T54,A5,T60,F14.5)') &
               istate, imo_occ, spin_label, nmo_occ_avail(ispin) + imo_virt, spin_label, e_virt_minus_occ(istate)*evolt

            DO jspin = 1, nspins
               ! .NOT. ASSOCIATED(evects(jspin, istate)%matrix))
               CALL cp_fm_create(evects(jspin, istate)%matrix, fm_struct_evects(jspin)%struct)
               CALL cp_fm_set_all(evects(jspin, istate)%matrix, 0.0_dp)

               IF (jspin == ispin) &
                  CALL cp_fm_to_fm(gs_mos(ispin)%mos_virt, evects(ispin, istate)%matrix, &
                                   ncol=1, source_start=imo_virt, target_start=imo_occ)
            END DO
         END IF
      END DO

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T7,A,T50,I24)') 'Number of active states:', tddfpt_total_number_of_states(gs_mos)
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      DEALLOCATE (e_virt_minus_occ)
      DEALLOCATE (inds)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_guess_vectors

END MODULE qs_tddfpt2_utils
