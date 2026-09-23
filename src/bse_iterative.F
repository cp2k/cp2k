!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Iterative solution of the Bethe-Salpeter equation: the block Davidson solvers on top of
!>        the matrix-free application of A and B, and the post-processing of their eigenvectors.
!>        The solvers live in bse_davidson, apart from this driver, so that a later Lanczos-Haydock
!>        solver for the spectrum can share the matrix-free A and B of bse_matvec and this driver
!> \par History
!>      09.2026 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_iterative

   USE bse_davidson,                    ONLY: abba_indefinite,&
                                              bse_davidson_abba_mk,&
                                              bse_davidson_abba_os,&
                                              bse_davidson_refcheck,&
                                              bse_davidson_tda
   USE bse_full_diag,                   ONLY: create_A_and_B,&
                                              postprocess_bse
   USE bse_matvec,                      ONLY: bse_matvec_create,&
                                              bse_matvec_env_type,&
                                              bse_matvec_release,&
                                              bse_matvec_selfcheck
   USE bse_print,                       ONLY: print_output_header
   USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE input_constants,                 ONLY: bse_abba_mk_davidson,&
                                              bse_abba_olsen_stratmann,&
                                              bse_screening_alpha,&
                                              bse_screening_rpa,&
                                              bse_screening_tdhf,&
                                              bse_singlet,&
                                              bse_triplet
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_iterative'

   PUBLIC :: solve_bse_iteratively

CONTAINS

! **************************************************************************************************
!> \brief Lowest excitations of a closed-shell reference from the block Davidson solvers on top of
!>        the matrix-free application of A and B, post-processed as in the full diagonalization.
!>        TDA (bse_davidson_tda): sum_jb A_ia,jb X_jb^n = Ω^n X_ia^n. Full problem, by ABBA_SOLVER:
!>        MK_DAVIDSON, sum_jb [(A+B)(A-B)]_ia,jb x_jb^n = (Ω^n)^2 x_ia^n with x^n = X^n - Y^n;
!>        OLSEN_STRATMANN, sum_jb (A+B)_ia,jb y_jb^n = Ω^n x_ia^n and sum_jb (A-B)_ia,jb x_jb^n = Ω^n y_ia^n
!>        with y^n = X^n + Y^n. Neither A nor B is formed: bse_matvec_apply contracts the RI slabs,
!>        (A Z)_ia = (ε_a-ε_i) Z_ia + α sum_P B^P_ia sum_jb B^P_jb Z_jb - w sum_Pjb \bar{B}^P_ij B^P_ab Z_jb
!>        (B Z)_ia =                  α sum_P B^P_ia sum_jb B^P_jb Z_jb - w sum_Pjb \bar{B}^P_ib B^P_ja Z_jb
!>        with α, w and \bar{B} set by create_matvec_env from SPIN_CONFIG and SCREENING_IN_W.
!>        The five slabs are released once the solver holds its RI-sliced copies
!> \param fm_mat_S_ia_trunc_arr B^P_ia of the active window, N_RI x homo_red*virt_red; entry 1 is used
!> \param fm_mat_S_ij_trunc_arr B^P_ij, N_RI x homo_red^2
!> \param fm_mat_S_ab_trunc_arr B^P_ab, N_RI x virt_red^2
!> \param fm_mat_S_bar_ia_bse_arr \bar{B}^P_ia = sum_Q [1+Q(0)]^-1_PQ B^Q_ia, the W-screened slab of W_ib,aj
!> \param fm_mat_S_bar_ij_bse_arr \bar{B}^P_ij, the W-screened slab of W_ij,ab
!> \param Eigenval_reduced quasiparticle energies of the homo_red + virt_red active levels in Hartree
!> \param homo_red_arr occupied levels in the active window, per spin; entry 1 is used
!> \param virt_red_arr virtual levels in the active window, per spin; entry 1 is used
!> \param homo occupied levels of the full reference, per spin; postprocess_bse indexes the MOs with it
!> \param dimen_RI number of RI functions; only the debug construction of the explicit A and B uses it
!> \param do_tda solve the TDA problem
!> \param do_abba solve the full problem; after the TDA in a TDA+ABBA run, seeded with X^TDA
!> \param mp2_env mp2_env%bse carries the BSE_ITERAT settings and BSE_DEBUG_PRINT
!> \param qs_env for the TDDFPT control flags and postprocess_bse
!> \param mo_coeff MO coefficients of the reference, handed to postprocess_bse
!> \param unit_nr output unit, positive on the writing rank only
! **************************************************************************************************
   SUBROUTINE solve_bse_iteratively(fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, &
                                    fm_mat_S_ab_trunc_arr, fm_mat_S_bar_ia_bse_arr, &
                                    fm_mat_S_bar_ij_bse_arr, Eigenval_reduced, &
                                    homo_red_arr, virt_red_arr, homo, dimen_RI, &
                                    do_tda, do_abba, mp2_env, qs_env, mo_coeff, unit_nr)

      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT) :: fm_mat_S_ia_trunc_arr, &
         fm_mat_S_ij_trunc_arr, fm_mat_S_ab_trunc_arr, fm_mat_S_bar_ia_bse_arr, &
         fm_mat_S_bar_ij_bse_arr
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval_reduced
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo_red_arr, virt_red_arr, homo
      INTEGER, INTENT(IN)                                :: dimen_RI
      LOGICAL, INTENT(IN)                                :: do_tda, do_abba
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_bse_iteratively'
      INTEGER, PARAMETER                                 :: n_ov_max_debug = 20000

      CHARACTER(LEN=16)                                  :: margin_str, n_ov_str, solver_str
      INTEGER                                            :: abba_status, handle, n_ov
      LOGICAL                                            :: do_debug
      REAL(KIND=dp)                                      :: ab_margin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: exc_ens, exc_ens_abba
      TYPE(bse_matvec_env_type)                          :: mv_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_ia
      TYPE(cp_fm_type)                                   :: fm_A_BSE, fm_B_BSE, fm_X_abba, fm_X_ia, &
                                                            fm_Y_abba, fm_Y_ia
      TYPE(cp_fm_type), POINTER                          :: fm_tda_guess
      TYPE(cp_fm_type), TARGET                           :: fm_X
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      IF (dft_control%tddfpt2_control%do_bse .OR. dft_control%tddfpt2_control%do_bse_w_only .OR. &
          dft_control%tddfpt2_control%do_bse_gw_only) THEN
         ! the TDDFPT path reads the explicit A and W matrices that create_A_and_B (FULLDIAG) stashes in the
         ! excited-state environment; ITERDIAG never forms them
         CALL cp_abort(__LOCATION__, &
                       "BSE_DIAG_METHOD ITERDIAG cannot provide the BSE kernel of TDDFPT "// &
                       "(DO_BSE, DO_BSE_W_ONLY, DO_BSE_GW_ONLY). Use BSE_DIAG_METHOD FULLDIAG.")
      END IF

      para_env => fm_mat_S_ia_trunc_arr(1)%matrix_struct%para_env
      n_ov = homo_red_arr(1)*virt_red_arr(1)
      ! the debug checks hold the explicit N_ov x N_ov matrices
      do_debug = mp2_env%bse%bse_debug_print .AND. n_ov <= n_ov_max_debug
      IF (mp2_env%bse%bse_debug_print .AND. .NOT. do_debug .AND. unit_nr > 0) THEN
         WRITE (n_ov_str, '(I0)') n_ov_max_debug
         CALL cp_warn(__LOCATION__, &
                      "BSE_DEBUG_PRINT: the checks of the iterative solver against the explicit A "// &
                      "and B are skipped above N_ov = "//TRIM(n_ov_str)//".")
      END IF

      CALL create_matvec_env(mv_env, fm_mat_S_ia_trunc_arr(1), fm_mat_S_ij_trunc_arr(1), &
                             fm_mat_S_ab_trunc_arr(1), fm_mat_S_bar_ia_bse_arr(1), &
                             fm_mat_S_bar_ij_bse_arr(1), Eigenval_reduced, &
                             homo_red_arr(1), virt_red_arr(1), do_abba, mp2_env, unit_nr)

      IF (do_debug) THEN
         CALL create_A_and_B(fm_mat_S_ia_trunc_arr, fm_mat_S_ij_trunc_arr, fm_mat_S_ab_trunc_arr, &
                             fm_mat_S_bar_ia_bse_arr, fm_mat_S_bar_ij_bse_arr, fm_A_BSE, fm_B_BSE, &
                             do_abba, Eigenval_reduced, unit_nr, homo_red_arr, virt_red_arr, dimen_RI, &
                             mp2_env, para_env, qs_env)
         IF (do_abba) THEN
            CALL bse_matvec_selfcheck(mv_env, fm_A_BSE, unit_nr, fm_B_BSE)
         ELSE
            CALL bse_matvec_selfcheck(mv_env, fm_A_BSE, unit_nr)
         END IF
      END IF

      ! the solver works on its RI-sliced copies; only the grid of the slabs is kept, for the eigenvectors
      blacs_env_ia => fm_mat_S_ia_trunc_arr(1)%matrix_struct%context
      CALL blacs_env_ia%retain()
      CALL cp_fm_release(fm_mat_S_ia_trunc_arr(1))
      CALL cp_fm_release(fm_mat_S_ij_trunc_arr(1))
      CALL cp_fm_release(fm_mat_S_ab_trunc_arr(1))
      CALL cp_fm_release(fm_mat_S_bar_ia_bse_arr(1))
      CALL cp_fm_release(fm_mat_S_bar_ij_bse_arr(1))

      IF (do_tda) THEN
         CALL print_output_header(homo_red_arr(1), virt_red_arr(1), homo(1), .TRUE., mp2_env, unit_nr)
         CALL bse_davidson_tda(mv_env, mp2_env, unit_nr, exc_ens, fm_X)
         IF (do_debug) CALL bse_davidson_refcheck(fm_A_BSE, exc_ens, fm_X, mp2_env, unit_nr)

         CALL copy_to_ia_context(fm_X, blacs_env_ia, fm_X_ia)
         CALL postprocess_bse(exc_ens, fm_X_ia, mp2_env, qs_env, mo_coeff, &
                              homo_red_arr(1), virt_red_arr(1), homo(1), unit_nr, .TRUE.)
         CALL cp_fm_release(fm_X_ia)
         DEALLOCATE (exc_ens)
      END IF

      IF (do_abba) THEN
         ! the TDA vectors are the initial guess of a combined run; a null pointer passed to the
         ! optional argument counts as absent (Fortran 2008)
         NULLIFY (fm_tda_guess)
         IF (do_tda) fm_tda_guess => fm_X
         CALL print_output_header(homo_red_arr(1), virt_red_arr(1), homo(1), .FALSE., mp2_env, unit_nr)
         SELECT CASE (mp2_env%bse%abba_solver)
         CASE (bse_abba_mk_davidson)
            solver_str = 'MK_DAVIDSON'
            CALL bse_davidson_abba_mk(mv_env, mp2_env, unit_nr, exc_ens_abba, fm_X_abba, fm_Y_abba, &
                                      abba_status, ab_margin, fm_tda_guess)
         CASE (bse_abba_olsen_stratmann)
            solver_str = 'OLSEN_STRATMANN'
            CALL bse_davidson_abba_os(mv_env, mp2_env, unit_nr, exc_ens_abba, fm_X_abba, fm_Y_abba, &
                                      abba_status, ab_margin, fm_tda_guess)
         CASE DEFAULT
            CPABORT("Unknown BSE_ITERAT%ABBA_SOLVER")
         END SELECT
         IF (abba_status == abba_indefinite) THEN
            WRITE (margin_str, '(F16.6)') ab_margin*evolt
            IF (.NOT. do_tda) THEN
               CALL cp_abort(__LOCATION__, &
                             "Matrix (A-B) is not positive definite on the trial space of the "// &
                             TRIM(solver_str)//" solver (eigenvalue "//TRIM(ADJUSTL(margin_str))//" eV). "// &
                             "The full ABBA problem is ill-defined for this reference.")
            END IF
            ! the TDA results of a combined run stand on their own
            IF (unit_nr > 0) THEN
               CALL cp_warn(__LOCATION__, &
                            "BSE: A-B is not positive definite on the trial space (margin "// &
                            TRIM(ADJUSTL(margin_str))//" eV); the ABBA calculation is skipped, "// &
                            "TDA results above are complete.")
               WRITE (unit_nr, '(T2,A4)') 'BSE|'
               WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.6)') 'BSE|', &
                  'A-B is not positive definite on the trial space (eV)', ab_margin*evolt
               WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', &
                  'The ABBA calculation is skipped, the TDA results above are complete'
            END IF
         ELSE
            IF (do_debug) THEN
               CALL bse_davidson_refcheck(fm_A_BSE, exc_ens_abba, fm_X_abba, mp2_env, unit_nr, &
                                          fm_B_BSE, fm_Y_abba)
            END IF

            CALL copy_to_ia_context(fm_X_abba, blacs_env_ia, fm_X_ia)
            CALL copy_to_ia_context(fm_Y_abba, blacs_env_ia, fm_Y_ia)
            CALL cp_fm_release(fm_X_abba)
            CALL cp_fm_release(fm_Y_abba)
            CALL postprocess_bse(exc_ens_abba, fm_X_ia, mp2_env, qs_env, mo_coeff, &
                                 homo_red_arr(1), virt_red_arr(1), homo(1), unit_nr, .FALSE., fm_Y_ia)
            CALL cp_fm_release(fm_X_ia)
            CALL cp_fm_release(fm_Y_ia)
            DEALLOCATE (exc_ens_abba)
         END IF
      END IF

      IF (do_tda) CALL cp_fm_release(fm_X)
      IF (do_debug) THEN
         CALL cp_fm_release(fm_A_BSE)
         IF (do_abba) CALL cp_fm_release(fm_B_BSE)
      END IF
      CALL cp_blacs_env_release(blacs_env_ia)
      CALL bse_matvec_release(mv_env)

      CALL timestop(handle)

   END SUBROUTINE solve_bse_iteratively

! **************************************************************************************************
!> \brief Copies eigenvectors from the process grid of the solver to the one of the RI slabs, from
!>        which postprocess_bse takes its grid
!> \param fm_vec eigenvectors on the npe x 1 process grid of the solver
!> \param blacs_env_ia process grid of the RI slabs
!> \param fm_vec_ia the copy on blacs_env_ia, created here
! **************************************************************************************************
   SUBROUTINE copy_to_ia_context(fm_vec, blacs_env_ia, fm_vec_ia)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_vec
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_ia
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_vec_ia

      INTEGER                                            :: n_exc, n_ov
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL cp_fm_get_info(fm_vec, nrow_global=n_ov, ncol_global=n_exc)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=blacs_env_ia%para_env, context=blacs_env_ia, &
                               nrow_global=n_ov, ncol_global=n_exc)
      CALL cp_fm_create(fm_vec_ia, fm_struct, name="fm_davidson_ia")
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_to_fm_submat_general(fm_vec, fm_vec_ia, n_ov, n_exc, 1, 1, 1, 1, blacs_env_ia)

   END SUBROUTINE copy_to_ia_context

! **************************************************************************************************
!> \brief Sets up the matrix-free application of A and B with the prefactors and the screened
!>        slabs that create_A and create_B use for the chosen SPIN_CONFIG and SCREENING_IN_W:
!>        α = 2 (SINGLET) or 0 (TRIPLET); w = 1 with the W-screened \bar{B} (default screening),
!>        w = SCREENING_FACTOR with \bar{B} = B (ALPHA), w = 1 with \bar{B} = B (TDHF), w = 0 (RPA)
!> \param mv_env the environment created here
!> \param fm_S_ia B^P_ia, N_RI x homo_red*virt_red
!> \param fm_S_ij B^P_ij, N_RI x homo_red^2
!> \param fm_S_ab B^P_ab, N_RI x virt_red^2
!> \param fm_S_bar_ia W-screened \bar{B}^P_ia
!> \param fm_S_bar_ij W-screened \bar{B}^P_ij
!> \param Eigenval_reduced quasiparticle energies of the homo_red + virt_red active levels in Hartree
!> \param homo_red occupied levels in the active window
!> \param virt_red virtual levels in the active window
!> \param do_abba B is applied too, so the \bar{B}^P_ia slab is sliced as well
!> \param mp2_env SPIN_CONFIG, SCREENING_IN_W, SCREENING_FACTOR and MATVEC_BLOCK_COLUMNS
!> \param unit_nr output unit, positive on the writing rank only
! **************************************************************************************************
   SUBROUTINE create_matvec_env(mv_env, fm_S_ia, fm_S_ij, fm_S_ab, fm_S_bar_ia, fm_S_bar_ij, &
                                Eigenval_reduced, homo_red, virt_red, do_abba, mp2_env, unit_nr)

      TYPE(bse_matvec_env_type), INTENT(OUT)             :: mv_env
      TYPE(cp_fm_type), INTENT(IN), TARGET               :: fm_S_ia, fm_S_ij, fm_S_ab, fm_S_bar_ia, &
                                                            fm_S_bar_ij
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval_reduced
      INTEGER, INTENT(IN)                                :: homo_red, virt_red
      LOGICAL, INTENT(IN)                                :: do_abba
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr

      REAL(KIND=dp)                                      :: alpha, w_fac
      TYPE(cp_fm_type), POINTER                          :: fm_W_ia, fm_W_ij

      alpha = 0.0_dp
      SELECT CASE (mp2_env%bse%bse_spin_config)
      CASE (bse_singlet)
         alpha = 2.0_dp
      CASE (bse_triplet)
         alpha = 0.0_dp
      END SELECT

      SELECT CASE (mp2_env%bse%screening_method)
      CASE (bse_screening_alpha)
         w_fac = mp2_env%bse%screening_factor
      CASE (bse_screening_rpa)
         w_fac = 0.0_dp
      CASE DEFAULT
         w_fac = 1.0_dp
      END SELECT

      ! TDHF and ALPHA screening contract the bare slabs, all others the W-multiplied ones
      IF (mp2_env%bse%screening_method == bse_screening_tdhf .OR. &
          mp2_env%bse%screening_method == bse_screening_alpha) THEN
         fm_W_ij => fm_S_ij
         fm_W_ia => fm_S_ia
      ELSE
         fm_W_ij => fm_S_bar_ij
         fm_W_ia => fm_S_bar_ia
      END IF
      CALL bse_matvec_create(mv_env, fm_S_ia, fm_W_ij, fm_S_ab, fm_W_ia, Eigenval_reduced, &
                             homo_red, virt_red, alpha, w_fac, do_abba, unit_nr, &
                             block_cols=mp2_env%bse%matvec_block_columns)

   END SUBROUTINE create_matvec_env

END MODULE bse_iterative
