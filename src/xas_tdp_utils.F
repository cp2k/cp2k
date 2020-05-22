!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for X-ray absorption spectroscopy using TDDFPT
!> \author AB (01.2018)
! **************************************************************************************************

MODULE xas_tdp_utils
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_power
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale,&
                                              cp_fm_transpose,&
                                              cp_fm_upper_to_full
   USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                              cp_fm_geeig
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, &
        cp_fm_release, cp_fm_set_element, cp_fm_to_fm_submat, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_block_p, &
        dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_reserve_all_blocks, dbcsr_set, &
        dbcsr_type
   USE input_constants,                 ONLY: tddfpt_singlet,&
                                              tddfpt_spin_cons,&
                                              tddfpt_spin_flip,&
                                              tddfpt_triplet,&
                                              xas_dip_len
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: get_diag
   USE physcon,                         ONLY: a_fine
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE xas_tdp_kernel,                  ONLY: kernel_coulomb_xc,&
                                              kernel_exchange
   USE xas_tdp_types,                   ONLY: donor_state_type,&
                                              xas_tdp_control_type,&
                                              xas_tdp_env_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_utils'

   PUBLIC :: setup_xas_tdp_prob, solve_xas_tdp_prob, include_rcs_soc, include_os_soc

   !A helper type for SOC
   TYPE dbcsr_soc_package_type
      TYPE(dbcsr_type), POINTER     :: dbcsr_sg
      TYPE(dbcsr_type), POINTER     :: dbcsr_tp
      TYPE(dbcsr_type), POINTER     :: dbcsr_sc
      TYPE(dbcsr_type), POINTER     :: dbcsr_sf
      TYPE(dbcsr_type), POINTER     :: dbcsr_prod
      TYPE(dbcsr_type), POINTER     :: dbcsr_ovlp
      TYPE(dbcsr_type), POINTER     :: dbcsr_tmp
      TYPE(dbcsr_type), POINTER     :: dbcsr_work
   END TYPE dbcsr_soc_package_type

CONTAINS

! **************************************************************************************************
!> \brief Builds the matrix that defines the XAS TDDFPT generalized eigenvalue problem to be solved
!>        for excitation energies omega. The problem has the form omega*G*C = M*C, where C contains
!>        the response orbitals coefficients. The matrix M and the metric G are stored in the given
!>        donor_state
!> \param donor_state the donor_state for which the problem is restricted
!> \param qs_env ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \note  the matrix M is symmetric and has the form | M_d   M_o |
!>                                                   | M_o   M_d |,
!>       -In the SPIN-RESTRICTED case:
!>        depending on whether we consider singlet or triplet excitation, the diagonal (M_d) and
!>        off-diagonal (M_o) parts of M differ:
!>        - For singlet: M_d = A + 2B + C_aa + C_ab - D
!>                       M_o = 2B + C_aa + C_ab - E
!>        - For triplet: M_d = A + C_aa - C_ab - D
!>                       M_o = C_aa - C_ab - E
!>        where other subroutines computes the matrices A, B, E, D and G, which are:
!>        - A: the ground-state contribution: F_pq*delta_IJ - epsilon_IJ*S_pq
!>        - B: the Coulomb kernel ~(pI|Jq)
!>        - C: the xc kernel c_aa (double derivatibe wrt to n_alpha) and C_ab (wrt n_alpha and n_beta)
!>        - D: the on-digonal exact exchange kernel ~(pq|IJ)
!>        - E: the off-diagonal exact exchange kernel ~(pJ|Iq)
!>        - G: the metric  S_pq*delta_IJ
!>        For the xc functionals, C_aa + C_ab or C_aa - C_ab are stored in the same matrix
!>        In the above definitions, I,J label the donnor MOs and p,q the sgfs of the basis
!>
!>       -In the SPIN-UNRESTRICTED, spin-conserving case:
!>        the on- and off-diagonal elements of M are:
!>                     M_d = A + B + C -D
!>                     M_o = B + C - E
!>        where the submatrices A, B, C, D and E are:
!>        - A: the groun-state contribution: (F_pq*delta_IJ - epsilon_IJ*S_pq) * delta_ab
!>        - B: the Coulomb kernel: (pI_a|J_b q)
!>        - C: the xc kernel: (pI_a|fxc_ab|J_b q)
!>        - D: the on-diagonal exact-exchange kernel: (pq|I_a J_b) delta_ab
!>        - E: the off-diagonal exact-exchange kernel: (pJ_b|I_a q) delta_ab
!>        - G: the metric S_pq*delta_IJ*delta_ab
!>        p,q label the sgfs, I,J the donro MOs and a,b the spins
!>
!>       -In both above cases, the matrix M is always  projected onto the unperturbed unoccupied
!>        ground state: M <= Q * M * Q^T = (1 - SP) * M * (1 - PS)
!>
!>       -In the SPIN-FLIP case:
!>        Only the TDA is implemented, that is, there are only on-diagonal elements:
!>                    M_d = A + C - D
!>        where the submatrices A, C and D are:
!>        - A: the ground state-contribution: (F_pq*delta_IJ - epsilon_IJ*S_pq) * delta_ab, but here,
!>                                            the alph-alpha quadrant has the beta Fock matrix and
!>                                            the beta-beta quadrant has the alpha Fock matrix
!>        - C: the SF xc kernel: (pI_a|fxc|J_bq), fxc = 1/m * (vxc_a -vxc_b)
!>        - D: the on-diagonal exact-exchange kernel: (pq|I_a J_b) delta_ab
!>        To ensure that all excitation start from a given spin to the opposite, we then multiply
!>        by a Q projector where we swap the alpha-alpha and beta-beta spin-quadrants
!>
!>        All possibilities: TDA or full-TDDFT, singlet or triplet, xc or hybrid, etc are treated
!>        in the same routine to avoid recomputing stuff
!>        Under TDA, only the on-diagonal elements of M are computed
!>        In the case of non-TDA, one turns the problem Hermitian
! **************************************************************************************************
   SUBROUTINE setup_xas_tdp_prob(donor_state, qs_env, xas_tdp_env, xas_tdp_control)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_xas_tdp_prob', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: submat_blk_size
      LOGICAL                                            :: do_coul, do_hfx, do_os, do_sc, do_sf, &
                                                            do_sg, do_tda, do_tp, do_xc
      REAL(dp)                                           :: eps_filter, sx
      TYPE(dbcsr_distribution_type), POINTER             :: submat_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ex_ker, xc_ker
      TYPE(dbcsr_type)                                   :: matrix_a, matrix_a_sf, matrix_b, proj_Q, &
                                                            proj_Q_sf, work
      TYPE(dbcsr_type), POINTER :: matrix_c_sc, matrix_c_sf, matrix_c_sg, matrix_c_tp, matrix_d, &
         matrix_e_sc, sc_matrix_tdp, sf_matrix_tdp, sg_matrix_tdp, tp_matrix_tdp

      NULLIFY (sg_matrix_tdp, tp_matrix_tdp, submat_dist, submat_blk_size, matrix_c_sf)
      NULLIFY (matrix_c_sg, matrix_c_tp, matrix_c_sc, matrix_d, matrix_e_sc)
      NULLIFY (sc_matrix_tdp, sf_matrix_tdp, ex_ker, xc_ker)

      CALL timeset(routineN, handle)

!  Initialization
      do_os = xas_tdp_control%do_uks .OR. xas_tdp_control%do_roks
      do_sc = xas_tdp_control%do_spin_cons
      do_sf = xas_tdp_control%do_spin_flip
      do_sg = xas_tdp_control%do_singlet
      do_tp = xas_tdp_control%do_triplet
      do_xc = xas_tdp_control%do_xc
      do_hfx = xas_tdp_control%do_hfx
      do_coul = xas_tdp_control%do_coulomb
      do_tda = xas_tdp_control%tamm_dancoff
      sx = xas_tdp_control%sx
      eps_filter = xas_tdp_control%eps_filter
      IF (do_sc) THEN
         ALLOCATE (donor_state%sc_matrix_tdp)
         sc_matrix_tdp => donor_state%sc_matrix_tdp
      END IF
      IF (do_sf) THEN
         ALLOCATE (donor_state%sf_matrix_tdp)
         sf_matrix_tdp => donor_state%sf_matrix_tdp
      END IF
      IF (do_sg) THEN
         ALLOCATE (donor_state%sg_matrix_tdp)
         sg_matrix_tdp => donor_state%sg_matrix_tdp
      END IF
      IF (do_tp) THEN
         ALLOCATE (donor_state%tp_matrix_tdp)
         tp_matrix_tdp => donor_state%tp_matrix_tdp
      END IF

!  Get the dist and block size of all matrices A, B, C, etc
      CALL compute_submat_dist_and_blk_size(donor_state, do_os, qs_env)
      submat_dist => donor_state%dbcsr_dist
      submat_blk_size => donor_state%blk_size

!  Allocate and compute all the matrices A, B, C, etc we will need

      ! The projector(s) on the unoccupied unperturbed ground state 1-SP and associated work matrix
      IF (do_sg .OR. do_tp .OR. do_sc) THEN !spin-conserving
         CALL get_q_projector(proj_Q, donor_state, do_os, xas_tdp_env)
      END IF
      IF (do_sf) THEN !spin-flip
         CALL get_q_projector(proj_Q_sf, donor_state, do_os, xas_tdp_env, do_sf=.TRUE.)
      END IF
      CALL dbcsr_create(matrix=work, matrix_type="N", dist=submat_dist, name="WORK", &
                        row_blk_size=submat_blk_size, col_blk_size=submat_blk_size)

      ! The ground state contribution(s)
      IF (do_sg .OR. do_tp .OR. do_sc) THEN !spin-conserving
         CALL build_gs_contribution(matrix_a, donor_state, do_os, qs_env)
      END IF
      IF (do_sf) THEN !spin-flip
         CALL build_gs_contribution(matrix_a_sf, donor_state, do_os, qs_env, do_sf=.TRUE.)
      END IF

      ! The Coulomb and XC kernels. Internal analysis to know which matrix to compute
      CALL dbcsr_allocate_matrix_set(xc_ker, 4)
      ALLOCATE (xc_ker(1)%matrix, xc_ker(2)%matrix, xc_ker(3)%matrix, xc_ker(4)%matrix)
      CALL kernel_coulomb_xc(matrix_b, xc_ker, donor_state, xas_tdp_env, xas_tdp_control, qs_env)
      matrix_c_sg => xc_ker(1)%matrix; matrix_c_tp => xc_ker(2)%matrix
      matrix_c_sc => xc_ker(3)%matrix; matrix_c_sf => xc_ker(4)%matrix

      ! The exact exchange. Internal analysis to know which matrices to compute
      CALL dbcsr_allocate_matrix_set(ex_ker, 2)
      ALLOCATE (ex_ker(1)%matrix, ex_ker(2)%matrix)
      CALL kernel_exchange(ex_ker, donor_state, xas_tdp_env, xas_tdp_control, qs_env)
      matrix_d => ex_ker(1)%matrix; matrix_e_sc => ex_ker(2)%matrix

      ! Build the metric G, also need its inverse in case of full-TDDFT
      IF (do_tda) THEN
         ALLOCATE (donor_state%metric(1))
         CALL build_metric(donor_state%metric, donor_state, qs_env, do_os)
      ELSE
         ALLOCATE (donor_state%metric(2))
         CALL build_metric(donor_state%metric, donor_state, qs_env, do_os, do_inv=.TRUE.)
      END IF

!  Build the eigenvalue problem, depending on the case (TDA, singlet, triplet, hfx, etc ...)
      IF (do_tda) THEN

         IF (do_sc) THEN ! open-shell spin-conserving under TDA

            ! The final matrix is M = A + B + C - D
            CALL dbcsr_copy(sc_matrix_tdp, matrix_a, name="OS MATRIX TDP")
            IF (do_coul) CALL dbcsr_add(sc_matrix_tdp, matrix_b, 1.0_dp, 1.0_dp)

            IF (do_xc) CALL dbcsr_add(sc_matrix_tdp, matrix_c_sc, 1.0_dp, 1.0_dp) !xc kernel
            IF (do_hfx) CALL dbcsr_add(sc_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx) !scaled hfx

            ! The product with the Q projector
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, sc_matrix_tdp, 0.0_dp, work, filter_eps=eps_filter)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, sc_matrix_tdp, filter_eps=eps_filter)

         END IF !do_sc

         IF (do_sf) THEN ! open-shell spin-flip under TDA

            ! The final matrix is M = A + C - D
            CALL dbcsr_copy(sf_matrix_tdp, matrix_a_sf, name="OS MATRIX TDP")

            IF (do_xc) CALL dbcsr_add(sf_matrix_tdp, matrix_c_sf, 1.0_dp, 1.0_dp) !xc kernel
            IF (do_hfx) CALL dbcsr_add(sf_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx) !scaled hfx

            ! Take the product with the (spin-flip) Q projector
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q_sf, sf_matrix_tdp, 0.0_dp, work, filter_eps=eps_filter)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q_sf, 0.0_dp, sf_matrix_tdp, filter_eps=eps_filter)

         END IF !do_sf

         IF (do_sg) THEN ! singlets under TDA

            ! The final matrix is M = A + 2B + (C_aa + C_ab) - D
            CALL dbcsr_copy(sg_matrix_tdp, matrix_a, name="SINGLET MATRIX TDP")
            IF (do_coul) CALL dbcsr_add(sg_matrix_tdp, matrix_b, 1.0_dp, 2.0_dp)

            IF (do_xc) CALL dbcsr_add(sg_matrix_tdp, matrix_c_sg, 1.0_dp, 1.0_dp) ! xc kernel
            IF (do_hfx) CALL dbcsr_add(sg_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx) ! scaled hfx

            ! Take the product with the Q projector:
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, sg_matrix_tdp, 0.0_dp, work, filter_eps=eps_filter)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, sg_matrix_tdp, filter_eps=eps_filter)

         END IF !do_sg (TDA)

         IF (do_tp) THEN ! triplets under TDA

            ! The final matrix is M =  A + (C_aa - C_ab) - D
            CALL dbcsr_copy(tp_matrix_tdp, matrix_a, name="TRIPLET MATRIX TDP")

            IF (do_xc) CALL dbcsr_add(tp_matrix_tdp, matrix_c_tp, 1.0_dp, 1.0_dp) ! xc_kernel
            IF (do_hfx) CALL dbcsr_add(tp_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx) ! scaled hfx

            ! Take the product with the Q projector:
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, tp_matrix_tdp, 0.0_dp, work, filter_eps=eps_filter)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, tp_matrix_tdp, filter_eps=eps_filter)

         END IF !do_tp (TDA)

      ELSE ! not TDA

         ! In the case of full-TDDFT, the problem is turned Hermitian with the help of auxiliary
         ! matrices AUX = (A-D+E)^(+-0.5) that are stored in donor_state
         CALL build_aux_matrix(1.0E-8_dp, sx, matrix_a, matrix_d, matrix_e_sc, do_hfx, proj_Q, &
                               work, donor_state, eps_filter, qs_env)

         IF (do_sc) THEN !full-TDDFT open-shell spin-conserving

            ! The final matrix is the sum of the on- and off-diagonal elements as in the description
            ! M = A + 2B + 2C - D - E
            CALL dbcsr_copy(sc_matrix_tdp, matrix_a, name="OS MATRIX TDP")
            IF (do_coul) CALL dbcsr_add(sc_matrix_tdp, matrix_b, 1.0_dp, 2.0_dp)

            IF (do_hfx) THEN !scaled hfx
               CALL dbcsr_add(sc_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx)
               CALL dbcsr_add(sc_matrix_tdp, matrix_e_sc, 1.0_dp, -1.0_dp*sx)
            END IF
            IF (do_xc) THEN
               CALL dbcsr_add(sc_matrix_tdp, matrix_c_sc, 1.0_dp, 2.0_dp)
            END IF

            ! Take the product with the Q projector
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, sc_matrix_tdp, 0.0_dp, work, filter_eps=eps_filter)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, sc_matrix_tdp, filter_eps=eps_filter)

            ! Take the product with the inverse metric
            ! M <= G^-1 * M * G^-1
            CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, sc_matrix_tdp, &
                                0.0_dp, work, filter_eps=eps_filter)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, work, donor_state%metric(2)%matrix, 0.0_dp, &
                                sc_matrix_tdp, filter_eps=eps_filter)

         END IF

         IF (do_sg) THEN ! full-TDDFT singlets

            ! The final matrix is the sum of the on- and off-diagonal elements as in the description
            ! M = A + 4B + 2(C_aa + C_ab) - D - E
            CALL dbcsr_copy(sg_matrix_tdp, matrix_a, name="SINGLET MATRIX TDP")
            IF (do_coul) CALL dbcsr_add(sg_matrix_tdp, matrix_b, 1.0_dp, 4.0_dp)

            IF (do_hfx) THEN !scaled hfx
               CALL dbcsr_add(sg_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx)
               CALL dbcsr_add(sg_matrix_tdp, matrix_e_sc, 1.0_dp, -1.0_dp*sx)
            END IF
            IF (do_xc) THEN !xc kernel
               CALL dbcsr_add(sg_matrix_tdp, matrix_c_sg, 1.0_dp, 2.0_dp)
            END IF

            ! Take the product with the Q projector
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, sg_matrix_tdp, 0.0_dp, work, filter_eps=eps_filter)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, sg_matrix_tdp, filter_eps=eps_filter)

            ! Take the product with the inverse metric
            ! M <= G^-1 * M * G^-1
            CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, sg_matrix_tdp, &
                                0.0_dp, work, filter_eps=eps_filter)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, work, donor_state%metric(2)%matrix, 0.0_dp, &
                                sg_matrix_tdp, filter_eps=eps_filter)

         END IF ! singlets

         IF (do_tp) THEN ! full-TDDFT triplets

            ! The final matrix is the sum of the on- and off-diagonal elements as in the description
            ! M = A + 2(C_aa - C_ab) - D - E
            CALL dbcsr_copy(tp_matrix_tdp, matrix_a, name="TRIPLET MATRIX TDP")

            IF (do_hfx) THEN !scaled hfx
               CALL dbcsr_add(tp_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx)
               CALL dbcsr_add(tp_matrix_tdp, matrix_e_sc, 1.0_dp, -1.0_dp*sx)
            END IF
            IF (do_xc) THEN
               CALL dbcsr_add(tp_matrix_tdp, matrix_c_tp, 1.0_dp, 2.0_dp)
            END IF

            ! Take the product with the Q projector
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, tp_matrix_tdp, 0.0_dp, work, filter_eps=eps_filter)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, tp_matrix_tdp, filter_eps=eps_filter)

            ! Take the product with the inverse metric
            ! M <= G^-1 * M * G^-1
            CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, tp_matrix_tdp, &
                                0.0_dp, work, filter_eps=eps_filter)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, work, donor_state%metric(2)%matrix, 0.0_dp, &
                                tp_matrix_tdp, filter_eps=eps_filter)

         END IF ! triplets

      END IF ! test on TDA

!  Clean-up
      CALL dbcsr_release(matrix_a)
      CALL dbcsr_release(matrix_a_sf)
      CALL dbcsr_release(matrix_b)
      CALL dbcsr_release(proj_Q)
      CALL dbcsr_release(proj_Q_sf)
      CALL dbcsr_release(work)
      CALL dbcsr_deallocate_matrix_set(ex_ker)
      CALL dbcsr_deallocate_matrix_set(xc_ker)

      CALL timestop(handle)

   END SUBROUTINE setup_xas_tdp_prob

! **************************************************************************************************
!> \brief Solves the XAS TDP generalized eigenvalue problem omega*C = matrix_tdp*C using standard
!>        full diagonalization methods. The problem is Hermitian (made that way even if not TDA)
!> \param donor_state ...
!> \param xas_tdp_control ...
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \param ex_type whether we deal with singlets, triplets, spin-conserving open-shell or spin-flip
!> \note The computed eigenvalues and eigenvectors are stored in the donor_state
!>       The eigenvectors are the LR-coefficients. In case of TDA, c^- is stored. In the general
!>       case, the sum c^+ + c^- is stored.
!>      - Spin-restricted:
!>       In case both singlets and triplets are considered, this routine must be called twice. This
!>       is the choice that was made because the body of the routine is exactly the same in both cases
!>       Note that for singlet we solve for u = 1/sqrt(2)*(c_alpha + c_beta) = sqrt(2)*c
!>       and that for triplets we solve for v = 1/sqrt(2)*(c_alpha - c_beta) = sqrt(2)*c
!>      - Spin-unrestricted:
!>       The problem is solved for the LR coefficients c_pIa as they are (not linear combination)
!>       The routine might be called twice (once for spin-conservign, one for spin-flip)
! **************************************************************************************************
   SUBROUTINE solve_xas_tdp_prob(donor_state, xas_tdp_control, xas_tdp_env, qs_env, ex_type)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ex_type

      CHARACTER(len=*), PARAMETER :: routineN = 'solve_xas_tdp_prob', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: first_ex, handle, i, imo, ispin, nao, &
                                                            ndo_mo, nelectron, nevals, nocc, nrow, &
                                                            nspins
      LOGICAL                                            :: do_os, do_range, do_sf
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: scaling, tmp_evals
      REAL(dp), DIMENSION(:), POINTER                    :: lr_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: ex_struct, fm_struct
      TYPE(cp_fm_type), POINTER                          :: c_diff, c_sum, lhs_matrix, lr_coeffs, &
                                                            rhs_matrix, work
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type)                                   :: tmp_mat
      TYPE(dbcsr_type), POINTER                          :: matrix_tdp

      CALL timeset(routineN, handle)

      NULLIFY (para_env, blacs_env, fm_struct, rhs_matrix, matrix_tdp, lhs_matrix, work)
      NULLIFY (c_diff, c_sum, ex_struct, lr_evals, lr_coeffs)
      CPASSERT(ASSOCIATED(xas_tdp_env))
      MARK_USED(xas_tdp_env)

      do_os = .FALSE.
      do_sf = .FALSE.
      IF (ex_type == tddfpt_spin_cons) THEN
         matrix_tdp => donor_state%sc_matrix_tdp
         do_os = .TRUE.
      ELSE IF (ex_type == tddfpt_spin_flip) THEN
         matrix_tdp => donor_state%sf_matrix_tdp
         do_os = .TRUE.
         do_sf = .TRUE.
      ELSE IF (ex_type == tddfpt_singlet) THEN
         matrix_tdp => donor_state%sg_matrix_tdp
      ELSE IF (ex_type == tddfpt_triplet) THEN
         matrix_tdp => donor_state%tp_matrix_tdp
      END IF
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env, nelectron_total=nelectron)

!     Initialization
      nspins = 1; IF (do_os) nspins = 2
      CALL cp_fm_get_info(donor_state%gs_coeffs, nrow_global=nao)
      CALL dbcsr_get_info(matrix_tdp, nfullrows_total=nrow)
      ndo_mo = donor_state%ndo_mo
      nocc = nelectron/2; IF (do_os) nocc = nelectron
      nocc = ndo_mo*nocc
      first_ex = nocc + 1 !where to find the first proper eigenvalue

      !solve by energy_range or number of states ?
      do_range = .FALSE.
      IF (xas_tdp_control%e_range > 0.0_dp) do_range = .TRUE.

      ! create the fm infrastructure
      ALLOCATE (tmp_evals(nrow))
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nrow, &
                               para_env=para_env, ncol_global=nrow)
      CALL cp_fm_create(c_sum, fm_struct)
      CALL cp_fm_create(rhs_matrix, fm_struct)
      CALL cp_fm_create(work, fm_struct)

!     Test on TDA
      IF (xas_tdp_control%tamm_dancoff) THEN

!        Get the main matrix_tdp as an fm
         CALL copy_dbcsr_to_fm(matrix_tdp, rhs_matrix)

!        Get the metric as a fm
         CALL cp_fm_create(lhs_matrix, fm_struct)
         CALL copy_dbcsr_to_fm(donor_state%metric(1)%matrix, lhs_matrix)

         !Diagonalisation (Cholesky decomposition). In TDA, c_sum = c^-
         CALL cp_fm_geeig(rhs_matrix, lhs_matrix, c_sum, tmp_evals, work)

!        TDA specific clean-up
         CALL cp_fm_release(lhs_matrix)

      ELSE ! not TDA

!        Need to multiply the current matrix_tdp with the auxiliary matrix
!        tmp_mat =  (A-D+E)^0.5 * M * (A-D+E)^0.5
         CALL dbcsr_create(matrix=tmp_mat, template=matrix_tdp, matrix_type="N")
         CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%matrix_aux, matrix_tdp, &
                             0.0_dp, tmp_mat, filter_eps=xas_tdp_control%eps_filter)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, tmp_mat, donor_state%matrix_aux, &
                             0.0_dp, tmp_mat, filter_eps=xas_tdp_control%eps_filter)

!        Get the matrix as a fm
         CALL copy_dbcsr_to_fm(tmp_mat, rhs_matrix)

!        Solve the "turned-Hermitian" eigenvalue problem
         CALL choose_eigv_solver(rhs_matrix, work, tmp_evals)

!        Currently, work = (A-D+E)^0.5 (c^+ - c^-) and tmp_evals = omega^2
!        Put tiny almost zero eigenvalues to zero (corresponding to occupied MOs)
         WHERE (tmp_evals < 1.0E-4_dp) tmp_evals = 0.0_dp

!        Retrieve c_diff = (c^+ - c^-) for normalization
!        (c^+ - c^-) = 1/omega^2 * M * (A-D+E)^0.5 * work
         CALL cp_fm_create(c_diff, fm_struct)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_tdp, donor_state%matrix_aux, &
                             0.0_dp, tmp_mat, filter_eps=xas_tdp_control%eps_filter)
         CALL cp_dbcsr_sm_fm_multiply(tmp_mat, work, c_diff, ncol=nrow)

         ALLOCATE (scaling(nrow))
         scaling = 0.0_dp
         WHERE (ABS(tmp_evals) > 1.0E-8_dp) scaling = 1.0_dp/tmp_evals
         CALL cp_fm_column_scale(c_diff, scaling)

!        Normalize with the metric: c_diff * G * c_diff = +- 1
         scaling = 0.0_dp
         CALL get_normal_scaling(scaling, c_diff, donor_state)
         CALL cp_fm_column_scale(c_diff, scaling)

!        Get the actual eigenvalues
         tmp_evals = SQRT(tmp_evals)

!        Get c_sum = (c^+ + c^-), which appears in all transition density related expressions
!        c_sum = -1/omega G^-1 * (A-D+E) * (c^+ - c^-)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%matrix_aux, donor_state%matrix_aux, &
                             0.0_dp, tmp_mat, filter_eps=xas_tdp_control%eps_filter)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, tmp_mat, &
                             0.0_dp, tmp_mat, filter_eps=xas_tdp_control%eps_filter)
         CALL cp_dbcsr_sm_fm_multiply(tmp_mat, c_diff, c_sum, ncol=nrow)
         WHERE (tmp_evals .NE. 0) scaling = -1.0_dp/tmp_evals
         CALL cp_fm_column_scale(c_sum, scaling)

!        Full TDDFT specific clean-up
         CALL cp_fm_release(c_diff)
         CALL dbcsr_release(tmp_mat)
         DEALLOCATE (scaling)

      END IF ! TDA

!     Full matrix clean-up
      CALL cp_fm_release(rhs_matrix)
      CALL cp_fm_release(work)

!  Reorganize the eigenvalues, we want a lr_evals array with the proper dimension and where the
!  first element is the first eval. Need a case study on do_range
      IF (do_range) THEN

         WHERE (tmp_evals > tmp_evals(first_ex) + xas_tdp_control%e_range) tmp_evals = 0.0_dp
         nevals = MAXLOC(tmp_evals, 1) - nocc

         ALLOCATE (lr_evals(nevals))
         lr_evals(:) = tmp_evals(first_ex:nocc + nevals)

      ELSE

         !Determine the number of evals to keep base on N_EXCITED
         nevals = nspins*nao - nocc/ndo_mo
         IF (xas_tdp_control%n_excited > 0 .AND. xas_tdp_control%n_excited < nevals) THEN
            nevals = xas_tdp_control%n_excited
         END IF
         nevals = ndo_mo*nevals !as in input description, multiply by # of donor MOs

         ALLOCATE (lr_evals(nevals))
         lr_evals(:) = tmp_evals(first_ex:nocc + nevals)
      END IF

!  Reorganize the eigenvectors in array of cp_fm so that each ndo_mo columns corresponds to an
!  excited state. Makes later calls to those easier and more efficient
!  In case of open-shell, we store the coeffs in the same logic as the matrix => first block where
!  the columns are the c_Ialpha and second block with columns as c_Ibeta
      CALL cp_fm_struct_create(ex_struct, nrow_global=nao, ncol_global=ndo_mo*nspins*nevals, &
                               para_env=para_env, context=blacs_env)
      CALL cp_fm_create(lr_coeffs, ex_struct)

      DO i = 1, nevals
         DO ispin = 1, nspins
            DO imo = 1, ndo_mo

               CALL cp_fm_to_fm_submat(msource=c_sum, mtarget=lr_coeffs, &
                                       nrow=nao, ncol=1, s_firstrow=((ispin - 1)*ndo_mo + imo - 1)*nao + 1, &
                                       s_firstcol=first_ex + i - 1, t_firstrow=1, &
                                       t_firstcol=(i - 1)*ndo_mo*nspins + (ispin - 1)*ndo_mo + imo)
            END DO !imo
         END DO !ispin
      END DO !istate

      IF (ex_type == tddfpt_spin_cons) THEN
         donor_state%sc_coeffs => lr_coeffs
         donor_state%sc_evals => lr_evals
      ELSE IF (ex_type == tddfpt_spin_flip) THEN
         donor_state%sf_coeffs => lr_coeffs
         donor_state%sf_evals => lr_evals
      ELSE IF (ex_type == tddfpt_singlet) THEN
         donor_state%sg_coeffs => lr_coeffs
         donor_State%sg_evals => lr_evals
      ELSE IF (ex_type == tddfpt_triplet) THEN
         donor_state%tp_coeffs => lr_coeffs
         donor_state%tp_evals => lr_evals
      END IF

!  Clean-up
      CALL cp_fm_release(c_sum)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_struct_release(ex_struct)

!  Perform a partial clean-up of the donor_state
      CALL dbcsr_release(matrix_tdp)

      CALL timestop(handle)

   END SUBROUTINE solve_xas_tdp_prob

! **************************************************************************************************
!> \brief Returns the scaling to apply to normalize the LR eigenvectors.
!> \param scaling the scaling array to apply
!> \param lr_coeffs the linear response coefficients as a fm
!> \param donor_state ...
!> \note The LR coeffs are normalized when c^T G c = +- 1, G is the metric, c = c^- for TDA and
!>       c = c^+ - c^- for the full problem
! **************************************************************************************************
   SUBROUTINE get_normal_scaling(scaling, lr_coeffs, donor_state)

      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: scaling
      TYPE(cp_fm_type), POINTER                          :: lr_coeffs
      TYPE(donor_state_type), POINTER                    :: donor_state

      CHARACTER(len=*), PARAMETER :: routineN = 'get_normal_scaling', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: nrow, nscal, nvals
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: norm_struct, work_struct
      TYPE(cp_fm_type), POINTER                          :: fm_norm, work
      TYPE(cp_para_env_type), POINTER                    :: para_env

      NULLIFY (para_env, blacs_env, norm_struct, work, fm_norm, work_struct)

!  Creating the matrix structures and initializing the work matrices
      CALL cp_fm_get_info(lr_coeffs, context=blacs_env, para_env=para_env, &
                          matrix_struct=work_struct, ncol_global=nvals, nrow_global=nrow)
      CALL cp_fm_struct_create(norm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=nvals, ncol_global=nvals)

      CALL cp_fm_create(work, work_struct)
      CALL cp_fm_create(fm_norm, norm_struct)

!  Taking c^T * G * C
      CALL cp_dbcsr_sm_fm_multiply(donor_state%metric(1)%matrix, lr_coeffs, work, ncol=nvals)
      CALL cp_gemm('T', 'N', nvals, nvals, nrow, 1.0_dp, lr_coeffs, work, 0.0_dp, fm_norm)

!  Computing the needed scaling
      ALLOCATE (diag(nvals))
      CALL cp_fm_get_diag(fm_norm, diag)
      WHERE (ABS(diag) > 1.0E-8_dp) diag = 1.0_dp/SQRT(ABS(diag))

      nscal = SIZE(scaling)
      scaling(1:nscal) = diag(1:nscal)

!  Clean-up
      CALL cp_fm_release(work)
      CALL cp_fm_release(fm_norm)
      CALL cp_fm_struct_release(norm_struct)

   END SUBROUTINE get_normal_scaling

! **************************************************************************************************
!> \brief This subroutine computes the row/column block structure as well as the dbcsr ditrinution
!>        for the submatrices making up the generalized XAS TDP eigenvalue problem. They all share
!>        the same properties, which are based on the replication of the KS matrix. Stored in the
!>        donor_state_type
!> \param donor_state ...
!> \param do_os whether this is a open-shell calculation
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_submat_dist_and_blk_size(donor_state, do_os, qs_env)

      TYPE(donor_state_type), POINTER                    :: donor_state
      LOGICAL, INTENT(IN)                                :: do_os
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_submat_dist_and_blk_size', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: group, i, nao, nblk_row, ndo_mo, nspins, &
                                                            scol_dist, srow_dist
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, col_dist_sub, row_blk_size, &
                                                            row_dist, row_dist_sub, submat_blk_size
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist, submat_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks

      NULLIFY (matrix_ks, dbcsr_dist, row_blk_size, row_dist, col_dist, pgrid, col_dist_sub)
      NULLIFY (row_dist_sub, submat_dist, submat_blk_size)

!  The submatrices are indexed by M_{pi sig,qj tau}, where p,q label basis functions and i,j donor
!  MOs and sig,tau the spins. For each spin and donor MOs combination, one has a submatrix of the
!  size of the KS matrix (nao x nao) with the same dbcsr block structure

!  Initialization
      ndo_mo = donor_state%ndo_mo
      CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_get_info(matrix_ks(1)%matrix, row_blk_size=row_blk_size)
      CALL dbcsr_distribution_get(dbcsr_dist, row_dist=row_dist, col_dist=col_dist, group=group, &
                                  pgrid=pgrid)
      nao = SUM(row_blk_size)
      nblk_row = SIZE(row_blk_size)
      srow_dist = SIZE(row_dist)
      scol_dist = SIZE(col_dist)
      nspins = 1; IF (do_os) nspins = 2

!  Creation if submatrix block size and col/row distribution
      ALLOCATE (submat_blk_size(ndo_mo*nspins*nblk_row))
      ALLOCATE (row_dist_sub(ndo_mo*nspins*srow_dist))
      ALLOCATE (col_dist_sub(ndo_mo*nspins*scol_dist))

      DO i = 1, ndo_mo*nspins
         submat_blk_size((i - 1)*nblk_row + 1:i*nblk_row) = row_blk_size
         row_dist_sub((i - 1)*srow_dist + 1:i*srow_dist) = row_dist
         col_dist_sub((i - 1)*scol_dist + 1:i*scol_dist) = col_dist
      END DO

!  Create the submatrix dbcsr distribution
      ALLOCATE (submat_dist)
      CALL dbcsr_distribution_new(submat_dist, group=group, pgrid=pgrid, row_dist=row_dist_sub, &
                                  col_dist=col_dist_sub)

      donor_state%dbcsr_dist => submat_dist
      donor_state%blk_size => submat_blk_size

!  Clean-up
      DEALLOCATE (col_dist_sub, row_dist_sub)

   END SUBROUTINE compute_submat_dist_and_blk_size

! **************************************************************************************************
!> \brief Returns the projector on the unperturbed unoccupied ground state Q = 1 - SP on the block
!>        diagonal of a matrix with the standard size and distribution.
!> \param proj_Q the matrix with the projector
!> \param donor_state ...
!> \param do_os whether it is open-shell calculation
!> \param xas_tdp_env ...
!> \param do_sf whether the projector should be prepared for spin-flip excitations
!> \note In the spin-flip case, the alpha spins are sent to beta and vice-versa. The structure of
!>       the projector is changed by swapping the alpha-alpha with the beta-beta block, which
!>       naturally take the spin change into account. Only for open-shell.
! **************************************************************************************************
   SUBROUTINE get_q_projector(proj_Q, donor_state, do_os, xas_tdp_env, do_sf)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: proj_Q
      TYPE(donor_state_type), POINTER                    :: donor_state
      LOGICAL, INTENT(IN)                                :: do_os
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_sf

      CHARACTER(len=*), PARAMETER :: routineN = 'get_q_projector', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, handle, iblk, imo, ispin, jblk, &
                                                            nblk_row, ndo_mo, nspins
      INTEGER, DIMENSION(:), POINTER                     :: blk_size_q, row_blk_size
      LOGICAL                                            :: found_block, my_dosf
      REAL(dp), DIMENSION(:), POINTER                    :: work_block
      TYPE(dbcsr_distribution_type), POINTER             :: dist_q
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type), POINTER                          :: one_sp

      NULLIFY (work_block, one_sp, row_blk_size, dist_q, blk_size_q)

      CALL timeset(routineN, handle)

!  Initialization
      nspins = 1; IF (do_os) nspins = 2
      ndo_mo = donor_state%ndo_mo
      one_sp => xas_tdp_env%q_projector(1)%matrix
      CALL dbcsr_get_info(one_sp, row_blk_size=row_blk_size)
      nblk_row = SIZE(row_blk_size)
      my_dosf = .FALSE.
      IF (PRESENT(do_sf)) my_dosf = do_sf
      dist_q => donor_state%dbcsr_dist
      blk_size_q => donor_state%blk_size

      ! the projector is not symmetric
      CALL dbcsr_create(matrix=proj_Q, name="PROJ Q", matrix_type="N", dist=dist_q, &
                        row_blk_size=blk_size_q, col_blk_size=blk_size_q)

!  Fill the projector by looping over 1-SP and duplicating blocks. (all on the spin-MO block diagonal)
      DO ispin = 1, nspins
         one_sp => xas_tdp_env%q_projector(ispin)%matrix

         !if spin-flip, swap the alpha-alpha and beta-beta blocks
         IF (my_dosf) one_sp => xas_tdp_env%q_projector(3 - ispin)%matrix

         CALL dbcsr_iterator_start(iter, one_sp)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

            ! get the block
            CALL dbcsr_get_block_p(one_sp, iblk, jblk, work_block, found_block)

            IF (found_block) THEN

               DO imo = 1, ndo_mo
                  CALL dbcsr_put_block(proj_Q, ((ispin - 1)*ndo_mo + imo - 1)*nblk_row + iblk, &
                                       ((ispin - 1)*ndo_mo + imo - 1)*nblk_row + jblk, work_block)
               END DO

            END IF
            NULLIFY (work_block)

         END DO !iterator
         CALL dbcsr_iterator_stop(iter)
      END DO !ispin

      CALL dbcsr_finalize(proj_Q)

      CALL timestop(handle)

   END SUBROUTINE get_q_projector

! **************************************************************************************************
!> \brief Builds the matrix containing the ground state contribution to the matrix_tdp (aka matrix A)
!>         => A_{pis,qjt} = (F_pq*delta_ij - epsilon_ij*S_pq) delta_st, where:
!>         F is the KS matrix
!>         S is the overlap matrix
!>         epsilon_ij is the donor MO eigenvalue
!>         i,j labels the MOs, p,q the AOs and s,t the spins
!> \param matrix_a  pointer to a DBCSR matrix containing A
!> \param donor_state ...
!> \param do_os ...
!> \param qs_env ...
!> \param do_sf whether the ground state contribution should accommodate spin-flip
!> \note Even localized non-canonical MOs are diagonalized in their subsapce => eps_ij = eps_ii*delta_ij
! **************************************************************************************************
   SUBROUTINE build_gs_contribution(matrix_a, donor_state, do_os, qs_env, do_sf)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_a
      TYPE(donor_state_type), POINTER                    :: donor_state
      LOGICAL, INTENT(IN)                                :: do_os
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_sf

      CHARACTER(len=*), PARAMETER :: routineN = 'build_gs_contribution', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, handle, iblk, imo, ispin, jblk, &
                                                            nblk_row, ndo_mo, nspins
      INTEGER, DIMENSION(:), POINTER                     :: blk_size_a, row_blk_size
      LOGICAL                                            :: found_block, my_dosf
      REAL(dp), DIMENSION(:), POINTER                    :: energy_evals, work_block
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist, dist_a
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: m_ks, matrix_ks, matrix_s
      TYPE(dbcsr_type)                                   :: work_matrix

      NULLIFY (matrix_ks, dbcsr_dist, row_blk_size, work_block, energy_evals, matrix_s, m_ks)
      NULLIFY (dist_a, blk_size_a)

      !  Note: matrix A is symmetric and block diagonal. If ADMM, the ks matrix is the total one,
      !        and it is corrected for eigenvalues (done at xas_tdp_init)

      CALL timeset(routineN, handle)

!  Initialization
      nspins = 1; IF (do_os) nspins = 2
      ndo_mo = donor_state%ndo_mo
      CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size)
      nblk_row = SIZE(row_blk_size)
      dist_a => donor_state%dbcsr_dist
      blk_size_a => donor_state%blk_size

!  Prepare the KS matrix pointer
      ALLOCATE (m_ks(nspins))
      m_ks(1)%matrix => matrix_ks(1)%matrix
      IF (do_os) m_ks(2)%matrix => matrix_ks(2)%matrix

! If spin-flip, swap the KS alpha-alpha and beta-beta quadrants.
      my_dosf = .FALSE.
      IF (PRESENT(do_sf)) my_dosf = do_sf
      IF (my_dosf .AND. do_os) THEN
         m_ks(1)%matrix => matrix_ks(2)%matrix
         m_ks(2)%matrix => matrix_ks(1)%matrix
      END IF

!  Creating the symmetric matrix A (and work)
      CALL dbcsr_create(matrix=matrix_a, name="MATRIX A", matrix_type="S", dist=dist_a, &
                        row_blk_size=blk_size_a, col_blk_size=blk_size_a)
      CALL dbcsr_create(matrix=work_matrix, name="WORK MAT", matrix_type="S", dist=dist_a, &
                        row_blk_size=blk_size_a, col_blk_size=blk_size_a)

      DO ispin = 1, nspins

!     Loop over the blocks of KS and put them on the spin-MO block diagonal of matrix A
         CALL dbcsr_iterator_start(iter, m_ks(ispin)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

!           Get the block
            CALL dbcsr_get_block_p(m_ks(ispin)%matrix, iblk, jblk, work_block, found_block)

            IF (found_block) THEN

!              The KS matrix only appears on diagonal of matrix A => loop over II donor MOs
               DO imo = 1, ndo_mo

!                 Put the block as it is
                  CALL dbcsr_put_block(matrix_a, ((ispin - 1)*ndo_mo + imo - 1)*nblk_row + iblk, &
                                       ((ispin - 1)*ndo_mo + imo - 1)*nblk_row + jblk, work_block)

               END DO !imo
            END IF !found_block
            NULLIFY (work_block)
         END DO ! iteration on KS blocks
         CALL dbcsr_iterator_stop(iter)

!     Loop over the blocks of S and put them on the block diagonal of work
         energy_evals => donor_state%energy_evals(:, ispin)

         CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

!           Get the block
            CALL dbcsr_get_block_p(matrix_s(1)%matrix, iblk, jblk, work_block, found_block)

            IF (found_block) THEN

!              Add S matrix on block diagonal as epsilon_ii*S_pq
               DO imo = 1, ndo_mo

                  CALL dbcsr_put_block(work_matrix, ((ispin - 1)*ndo_mo + imo - 1)*nblk_row + iblk, &
                                       ((ispin - 1)*ndo_mo + imo - 1)*nblk_row + jblk, &
                                       energy_evals(imo)*work_block)
               END DO !imo
            END IF !found block
            NULLIFY (work_block)
         END DO ! iteration on S blocks
         CALL dbcsr_iterator_stop(iter)

      END DO !ispin
      CALL dbcsr_finalize(matrix_a)
      CALL dbcsr_finalize(work_matrix)

!  Take matrix_a = matrix_a - work
      CALL dbcsr_add(matrix_a, work_matrix, 1.0_dp, -1.0_dp)
      CALL dbcsr_finalize(matrix_a)

!  Clean-up
      CALL dbcsr_release(work_matrix)
      DEALLOCATE (m_ks)

      CALL timestop(handle)

   END SUBROUTINE build_gs_contribution

! **************************************************************************************************
!> \brief Creates the metric (aka  matrix G) needed for the generalized eigenvalue problem and inverse
!>         => G_{pis,qjt} = S_pq*delta_ij*delta_st
!> \param matrix_g dbcsr matrix containing G
!> \param donor_state ...
!> \param qs_env ...
!> \param do_os if open-shell calculation
!> \param do_inv if the inverse of G should be computed
! **************************************************************************************************
   SUBROUTINE build_metric(matrix_g, donor_state, qs_env, do_os, do_inv)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_g
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: do_os
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_inv

      CHARACTER(len=*), PARAMETER :: routineN = 'build_metric', routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, handle, i, iblk, jblk, nao, &
                                                            nblk_row, ndo_mo, nspins
      INTEGER, DIMENSION(:), POINTER                     :: blk_size_g, row_blk_size
      LOGICAL                                            :: found_block, my_do_inv
      REAL(dp), DIMENSION(:), POINTER                    :: work_block
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: dist_g
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: matrix_sinv

      NULLIFY (matrix_s, row_blk_size, work_block, para_env, blacs_env, dist_g, blk_size_g)

      CALL timeset(routineN, handle)

!  Initialization
      nspins = 1; IF (do_os) nspins = 2
      ndo_mo = donor_state%ndo_mo
      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size, nfullrows_total=nao)
      nblk_row = SIZE(row_blk_size)
      my_do_inv = .FALSE.
      IF (PRESENT(do_inv)) my_do_inv = do_inv
      dist_g => donor_state%dbcsr_dist
      blk_size_g => donor_state%blk_size

!  Creating the symmetric  matrices G and G^-1 with the right size and distribution
      ALLOCATE (matrix_g(1)%matrix)
      CALL dbcsr_create(matrix=matrix_g(1)%matrix, name="MATRIX G", matrix_type="S", dist=dist_g, &
                        row_blk_size=blk_size_g, col_blk_size=blk_size_g)

!  Fill the matrices G by looping over the block of S and putting them on the diagonal
      CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

!        Get the block
         CALL dbcsr_get_block_p(matrix_s(1)%matrix, iblk, jblk, work_block, found_block)

         IF (found_block) THEN

!           Go over the diagonal of G => donor MOs ii, spin ss
            DO i = 1, ndo_mo*nspins
               CALL dbcsr_put_block(matrix_g(1)%matrix, (i - 1)*nblk_row + iblk, (i - 1)*nblk_row + jblk, work_block)
            END DO

         END IF
         NULLIFY (work_block)

      END DO ! dbcsr_iterator
      CALL dbcsr_iterator_stop(iter)

!  Finalize
      CALL dbcsr_finalize(matrix_g(1)%matrix)

!  If the inverse of G is required, do the same as above with the inverse
      IF (my_do_inv) THEN

         CPASSERT(SIZE(matrix_g) == 2)

         ! Create the matrix
         ALLOCATE (matrix_g(2)%matrix)
         CALL dbcsr_create(matrix=matrix_g(2)%matrix, name="MATRIX GINV", matrix_type="S", &
                           dist=dist_g, row_blk_size=blk_size_g, col_blk_size=blk_size_g)

         ! Invert the overlap matrix
         CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
         CALL dbcsr_copy(matrix_sinv, matrix_s(1)%matrix)
         CALL cp_dbcsr_cholesky_decompose(matrix_sinv, para_env=para_env, blacs_env=blacs_env)
         CALL cp_dbcsr_cholesky_invert(matrix_sinv, para_env=para_env, blacs_env=blacs_env, upper_to_full=.TRUE.)

!     Fill the matrices G^-1 by looping over the block of S^-1 and putting them on the diagonal
         CALL dbcsr_iterator_start(iter, matrix_sinv)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

!           Get the block
            CALL dbcsr_get_block_p(matrix_sinv, iblk, jblk, work_block, found_block)

            IF (found_block) THEN

!              Go over the diagonal of G => donor MOs ii spin ss
               DO i = 1, ndo_mo*nspins
                  CALL dbcsr_put_block(matrix_g(2)%matrix, (i - 1)*nblk_row + iblk, (i - 1)*nblk_row + jblk, work_block)
               END DO

            END IF
            NULLIFY (work_block)

         END DO ! dbcsr_iterator
         CALL dbcsr_iterator_stop(iter)

         !  Finalize
         CALL dbcsr_finalize(matrix_g(2)%matrix)

         !  Clean-up
         CALL dbcsr_release(matrix_sinv)
      END IF !do_inv

      CALL timestop(handle)

   END SUBROUTINE build_metric

! **************************************************************************************************
!> \brief Builds the auxiliary matrix (A-D+E)^+0.5 needed for the transofrmation of the
!>        full-TDDFT problem into an Hermitian one
!> \param threshold a threshold for allowed negative eigenvalues
!> \param sx the amount of exact exchange
!> \param matrix_a the ground state contribution matrix A
!> \param matrix_d the on-diagonal exchange kernel matrix (ab|IJ)
!> \param matrix_e the off-diagonal exchange kernel matrix (aJ|Ib)
!> \param do_hfx if exact exchange is included
!> \param proj_Q ...
!> \param work ...
!> \param donor_state ...
!> \param eps_filter for the dbcsr multiplication
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE build_aux_matrix(threshold, sx, matrix_a, matrix_d, matrix_e, do_hfx, proj_Q, &
                               work, donor_state, eps_filter, qs_env)

      REAL(dp), INTENT(IN)                               :: threshold, sx
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_a, matrix_d, matrix_e
      LOGICAL, INTENT(IN)                                :: do_hfx
      TYPE(dbcsr_type), INTENT(INOUT)                    :: proj_Q, work
      TYPE(donor_state_type), POINTER                    :: donor_state
      REAL(dp), INTENT(IN)                               :: eps_filter
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_aux_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ndep
      REAL(dp)                                           :: evals(2)
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type)                                   :: tmp_mat

      NULLIFY (blacs_env, para_env)

      CALL timeset(routineN, handle)

      CALL dbcsr_copy(tmp_mat, matrix_a)
      IF (do_hfx) THEN
         CALL dbcsr_add(tmp_mat, matrix_d, 1.0_dp, -1.0_dp*sx) !scaled hfx
         CALL dbcsr_add(tmp_mat, matrix_e, 1.0_dp, 1.0_dp*sx)
      END IF

      ! Take the product with the Q projector:
      CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, tmp_mat, 0.0_dp, work, filter_eps=eps_filter)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, tmp_mat, filter_eps=eps_filter)

      ! Actually computing and storing the auxiliary matrix
      ALLOCATE (donor_state%matrix_aux)
      CALL dbcsr_create(matrix=donor_state%matrix_aux, template=matrix_a, name="MAT AUX")

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      ! good quality sqrt
      CALL cp_dbcsr_power(tmp_mat, 0.5_dp, threshold, ndep, para_env, blacs_env, eigenvalues=evals)

      CALL dbcsr_copy(donor_state%matrix_aux, tmp_mat)

      ! Warn the user if matrix not positive semi-definite
      IF (evals(1) < 0.0_dp .AND. ABS(evals(1)) > threshold) THEN
         CPWARN("The full TDDFT problem might not have been soundly turned Hermitian. Try TDA.")
      END IF

      ! clean-up
      CALL dbcsr_release(tmp_mat)

      CALL timestop(handle)

   END SUBROUTINE build_aux_matrix

! **************************************************************************************************
!> \brief Includes the SOC effects on the precomputed spin-conserving and spin-flip excitations
!>        from an open-shell calculation (UKS or ROKS). This is a perturbative treatment
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note Using AMEWs, build an hermitian matrix with all excited states SOC coupling + the
!>       excitation energies on the diagonal. Then diagonalize it to get the new excitation
!>       energies and corresponding linear combinations of lr_coeffs.
!>       The AMEWs are normalized
!>       Only for open-shell calculations
! **************************************************************************************************
   SUBROUTINE include_os_soc(donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'include_os_soc', routineP = moduleN//':'//routineN

      INTEGER                                            :: group, handle, homo, iex, isc, isf, nao, &
                                                            ndo_mo, ndo_so, nex, npcols, nprows, &
                                                            nsc, nsf, ntot, tas(2), tbs(2)
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_dist, row_blk_size, &
                                                            row_dist, row_dist_new
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      LOGICAL                                            :: do_roks, do_uks
      REAL(dp)                                           :: eps_filter, gs_sum, soc
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag, tmp_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: domo_soc_x, domo_soc_y, domo_soc_z, &
                                                            gsex_block
      REAL(dp), DIMENSION(:), POINTER                    :: sc_evals, sf_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type), POINTER                         :: evecs_cfm, pert_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: full_struct, gsex_struct, prod_struct, &
                                                            vec_struct, work_struct
      TYPE(cp_fm_type), POINTER :: gs_coeffs, gsex_fm, img_fm, mo_coeff, prod_work, real_fm, &
         sc_coeffs, sf_coeffs, vec_soc_x, vec_soc_y, vec_soc_z, vec_work, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: coeffs_dist, dbcsr_dist, prod_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_soc_package_type)                       :: dbcsr_soc_package
      TYPE(dbcsr_type), POINTER                          :: dbcsr_ovlp, dbcsr_prod, dbcsr_sc, &
                                                            dbcsr_sf, dbcsr_tmp, dbcsr_work, &
                                                            orb_soc_x, orb_soc_y, orb_soc_z
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      NULLIFY (gs_coeffs, sc_coeffs, sf_coeffs, matrix_s, orb_soc_x, orb_soc_y, orb_soc_z, mos)
      NULLIFY (real_fm, img_fm, full_struct, para_env, blacs_env, mo_coeff, vec_soc_x, work_fm)
      NULLIFY (sc_evals, sf_evals, vec_work, prod_work, vec_struct, prod_struct, vec_soc_y)
      NULLIFY (vec_soc_z, pert_cfm, evecs_cfm, work_struct, gsex_struct, col_dist, row_dist)
      NULLIFY (col_blk_size, row_blk_size, row_dist_new, pgrid, dbcsr_sc, dbcsr_sf, dbcsr_work)
      NULLIFY (dbcsr_tmp, dbcsr_ovlp, dbcsr_prod)

      CALL timeset(routineN, handle)

! Initialization
      sc_coeffs => donor_state%sc_coeffs
      sf_coeffs => donor_state%sf_coeffs
      sc_evals => donor_state%sc_evals
      sf_evals => donor_state%sf_evals
      nsc = SIZE(sc_evals)
      nsf = SIZE(sf_evals)
      ntot = 1 + nsc + nsf
      nex = nsc !by contrutciotn nsc == nsf, but keep 2 counts for clarity
      ndo_mo = donor_state%ndo_mo
      ndo_so = 2*ndo_mo
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env, mos=mos, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      orb_soc_x => xas_tdp_env%orb_soc(1)%matrix
      orb_soc_y => xas_tdp_env%orb_soc(2)%matrix
      orb_soc_z => xas_tdp_env%orb_soc(3)%matrix
      do_roks = xas_tdp_control%do_roks
      do_uks = xas_tdp_control%do_uks
      eps_filter = xas_tdp_control%eps_filter

      ! For the GS coeffs, we use the same structure both for ROKS and UKS here => allows us to write
      ! general code later on, and not use IF (do_roks) statements every second line
      IF (do_uks) gs_coeffs => donor_state%gs_coeffs
      IF (do_roks) THEN
         CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=nao, ncol_global=ndo_so)
         CALL cp_fm_create(gs_coeffs, vec_struct)

         ! only alpha donor MOs are stored, need to copy them intoboth the alpha and the beta slot
         CALL cp_fm_to_fm_submat(msource=donor_state%gs_coeffs, mtarget=gs_coeffs, nrow=nao, &
                                 ncol=ndo_mo, s_firstrow=1, s_firstcol=1, t_firstrow=1, &
                                 t_firstcol=1)
         CALL cp_fm_to_fm_submat(msource=donor_state%gs_coeffs, mtarget=gs_coeffs, nrow=nao, &
                                 ncol=ndo_mo, s_firstrow=1, s_firstcol=1, t_firstrow=1, &
                                 t_firstcol=ndo_mo + 1)

         CALL cp_fm_struct_release(vec_struct)
      END IF

! Creating the real and the imaginary part of the SOC perturbation matrix
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      CALL cp_fm_create(real_fm, full_struct)
      CALL cp_fm_create(img_fm, full_struct)

! Put the excitation energies on the diagonal of the real  matrix. Element 1,1 is the ground state
      DO isc = 1, nsc
         CALL cp_fm_set_element(real_fm, 1 + isc, 1 + isc, sc_evals(isc))
      END DO
      DO isf = 1, nsf
         CALL cp_fm_set_element(real_fm, 1 + nsc + isf, 1 + nsc + isf, sf_evals(isf))
      END DO

! Create the bdcsr machinery
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_distribution_get(dbcsr_dist, group=group, row_dist=row_dist, pgrid=pgrid, &
                                  npcols=npcols, nprows=nprows)
      ALLOCATE (col_dist(nex), row_dist_new(nex))
      DO iex = 1, nex
         col_dist(iex) = MODULO(npcols - iex, npcols)
         row_dist_new(iex) = MODULO(nprows - iex, nprows)
      END DO
      ALLOCATE (coeffs_dist, prod_dist)
      CALL dbcsr_distribution_new(coeffs_dist, group=group, pgrid=pgrid, row_dist=row_dist, &
                                  col_dist=col_dist)
      CALL dbcsr_distribution_new(prod_dist, group=group, pgrid=pgrid, row_dist=row_dist_new, &
                                  col_dist=col_dist)

      !Create the matrices
      ALLOCATE (col_blk_size(nex))
      col_blk_size = ndo_so
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size)

      ALLOCATE (dbcsr_sc, dbcsr_sf, dbcsr_work, dbcsr_ovlp, dbcsr_tmp, dbcsr_prod)
      CALL dbcsr_create(matrix=dbcsr_sc, name="SPIN CONS", matrix_type="N", dist=coeffs_dist, &
                        row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_sf, name="SPIN FLIP", matrix_type="N", dist=coeffs_dist, &
                        row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_work, name="WORK", matrix_type="N", dist=coeffs_dist, &
                        row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_prod, name="PROD", matrix_type="N", dist=prod_dist, &
                        row_blk_size=col_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_ovlp, name="OVLP", matrix_type="N", dist=prod_dist, &
                        row_blk_size=col_blk_size, col_blk_size=col_blk_size)

      col_blk_size = 1
      CALL dbcsr_create(matrix=dbcsr_tmp, name="TMP", matrix_type="N", dist=prod_dist, &
                        row_blk_size=col_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_reserve_all_blocks(dbcsr_tmp)

      dbcsr_soc_package%dbcsr_sc => dbcsr_sc
      dbcsr_soc_package%dbcsr_sf => dbcsr_sf
      dbcsr_soc_package%dbcsr_work => dbcsr_work
      dbcsr_soc_package%dbcsr_ovlp => dbcsr_ovlp
      dbcsr_soc_package%dbcsr_prod => dbcsr_prod
      dbcsr_soc_package%dbcsr_tmp => dbcsr_tmp

      !Filling the coeffs matrices by copying from the stored fms
      CALL copy_fm_to_dbcsr(sc_coeffs, dbcsr_sc)
      CALL copy_fm_to_dbcsr(sf_coeffs, dbcsr_sf)

! Precompute what we can before looping over excited states.
      ! Need to compute the scalar: sum_i sum_sigma <phi^0_i,sigma|SOC|phi^0_i,sigma>, where all
      ! occupied MOs are taken into account

      !start with the alpha MOs
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, homo=homo)
      ALLOCATE (diag(homo))
      CALL cp_fm_get_info(mo_coeff, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=homo, ncol_global=homo)
      CALL cp_fm_create(vec_work, vec_struct)
      CALL cp_fm_create(prod_work, prod_struct)

      ! <alpha|SOC_z|alpha> => spin integration yields +1
      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, mo_coeff, vec_work, ncol=homo)
      CALL cp_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, vec_work, 0.0_dp, prod_work)
      CALL cp_fm_get_diag(prod_work, diag)
      gs_sum = SUM(diag)

      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (diag)
      NULLIFY (vec_struct)

      ! Now do the same with the beta gs coeffs
      CALL get_mo_set(mos(2)%mo_set, mo_coeff=mo_coeff, homo=homo)
      ALLOCATE (diag(homo))
      CALL cp_fm_get_info(mo_coeff, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=homo, ncol_global=homo)
      CALL cp_fm_create(vec_work, vec_struct)
      CALL cp_fm_create(prod_work, prod_struct)

      ! <beta|SOC_z|beta> => spin integration yields -1
      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, mo_coeff, vec_work, ncol=homo)
      CALL cp_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, vec_work, 0.0_dp, prod_work)
      CALL cp_fm_get_diag(prod_work, diag)
      gs_sum = gs_sum - SUM(diag) ! -1 because of spin integration

      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (diag)

      ! Need to compute: <phi^0_Isigma|SOC|phi^0_Jtau> for the donor MOs

      CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nao, ncol_global=ndo_so)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_so, ncol_global=ndo_so)
      CALL cp_fm_create(vec_soc_x, vec_struct) ! for SOC_x*gs_coeffs
      CALL cp_fm_create(vec_soc_y, vec_struct) ! for SOC_y*gs_coeffs
      CALL cp_fm_create(vec_soc_z, vec_struct) ! for SOC_z*gs_coeffs
      CALL cp_fm_create(prod_work, prod_struct)
      ALLOCATE (diag(ndo_so))

      ALLOCATE (domo_soc_x(ndo_so, ndo_so), domo_soc_y(ndo_so, ndo_so), domo_soc_z(ndo_so, ndo_so))

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_x, gs_coeffs, vec_soc_x, ncol=ndo_so)
      CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, vec_soc_x, 0.0_dp, prod_work)
      CALL cp_fm_get_submatrix(prod_work, domo_soc_x)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, gs_coeffs, vec_soc_y, ncol=ndo_so)
      CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, vec_soc_y, 0.0_dp, prod_work)
      CALL cp_fm_get_submatrix(prod_work, domo_soc_y)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, gs_coeffs, vec_soc_z, ncol=ndo_so)
      CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, vec_soc_z, 0.0_dp, prod_work)
      CALL cp_fm_get_submatrix(prod_work, domo_soc_z)

      ! some more useful work matrices
      CALL cp_fm_struct_create(work_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nex, ncol_global=nex)
      CALL cp_fm_create(work_fm, work_struct)

!  Looping over the excited states, computing the SOC and filling the perturbation matrix
!  There are 3 loops to do: sc-sc, sc-sf and sf-sf
!  The final perturbation matrix is Hermitian, only the upper diagonal is needed

      !need some work matrices for the GS stuff
      CALL cp_fm_struct_create(gsex_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nex*ndo_so, ncol_global=ndo_so)
      CALL cp_fm_create(gsex_fm, gsex_struct)
      ALLOCATE (gsex_block(ndo_so, ndo_so))

!  Start with ground-state/spin-conserving SOC:
      !  <Psi_0|SOC|Psi_Jsc> = sum_k,sigma <phi^0_k,sigma|SOC|phi^Jsc_k,sigma>

      !compute -sc_coeffs*SOC_Z*gs_coeffs, minus sign because SOC_z antisymmetric
      CALL cp_gemm('T', 'N', nex*ndo_so, ndo_so, nao, -1.0_dp, sc_coeffs, vec_soc_z, 0.0_dp, gsex_fm)

      DO isc = 1, nsc
         CALL cp_fm_get_submatrix(fm=gsex_fm, target_m=gsex_block, start_row=(isc - 1)*ndo_so + 1, &
                                  start_col=1, n_rows=ndo_so, n_cols=ndo_so)
         diag(:) = get_diag(gsex_block)
         soc = SUM(diag(1:ndo_mo)) - SUM(diag(ndo_mo + 1:ndo_so)) !minus sign because of spin integration

         !purely imaginary contribution
         CALL cp_fm_set_element(img_fm, 1, 1 + isc, soc)
      END DO !isc

!  Then ground-state/spin-flip SOC:
      !<Psi_0|SOC|Psi_Jsf> = sum_k,sigma <phi^0_k,sigma|SOC|phi^Jsc_k,tau>   sigma != tau

      !compute  -sc_coeffs*SOC_x*gs_coeffs, imaginary contribution
      CALL cp_gemm('T', 'N', nex*ndo_so, ndo_so, nao, -1.0_dp, sc_coeffs, vec_soc_x, 0.0_dp, gsex_fm)

      DO isf = 1, nsf
         CALL cp_fm_get_submatrix(fm=gsex_fm, target_m=gsex_block, start_row=(isf - 1)*ndo_so + 1, &
                                  start_col=1, n_rows=ndo_so, n_cols=ndo_so)
         diag(:) = get_diag(gsex_block)
         soc = SUM(diag) !alpha and beta parts are simply added due to spin integration
         CALL cp_fm_set_element(img_fm, 1, 1 + nsc + isf, soc)
      END DO !isf

      !compute -sc_coeffs*SOC_y*gs_coeffs, real contribution
      CALL cp_gemm('T', 'N', nex*ndo_so, ndo_so, nao, -1.0_dp, sc_coeffs, vec_soc_y, 0.0_dp, gsex_fm)

      DO isf = 1, nsf
         CALL cp_fm_get_submatrix(fm=gsex_fm, target_m=gsex_block, start_row=(isf - 1)*ndo_so + 1, &
                                  start_col=1, n_rows=ndo_so, n_cols=ndo_so)
         diag(:) = get_diag(gsex_block)
         soc = SUM(diag(1:ndo_mo)) ! alpha-beta
         soc = soc - SUM(diag(ndo_mo + 1:ndo_so)) !beta-alpha
         CALL cp_fm_set_element(real_fm, 1, 1 + nsc + isf, soc)
      END DO !isf

      !ground-state cleanup
      CALL cp_fm_release(gsex_fm)
      CALL cp_fm_release(vec_soc_x)
      CALL cp_fm_release(vec_soc_y)
      CALL cp_fm_release(vec_soc_z)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(gsex_struct)
      CALL cp_fm_struct_release(prod_struct)
      CALL cp_fm_struct_release(vec_struct)
      DEALLOCATE (gsex_block)

!  Then spin-conserving/spin-conserving SOC
!  <Psi_Isc|SOC|Psi_Jsc> =
!  sum_k,sigma [<psi^Isc_k,sigma|SOC|psi^Jsc_k,sigma> + <psi^Isc_k,sigma|psi^Jsc_k,sigma> * gs_sum]
!  - sum_k,l,sigma <psi^0_k,sigma|SOC|psi^0_l,sigma> * <psi^Isc_l,sigma|psi^Jsc_k,sigma>

      !Same spin integration => only SOC_z matters, and the contribution is purely imaginary
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_z, dbcsr_sc, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      !the overlap as well
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_sc, 0.0_dp, dbcsr_work, &
                          filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

      CALL os_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_soc_z, pref_diaga=1.0_dp, &
                                pref_diagb=-1.0_dp, pref_tracea=-1.0_dp, pref_traceb=1.0_dp, &
                                pref_diags=gs_sum, symmetric=.TRUE.)

      CALL copy_dbcsr_to_fm(dbcsr_tmp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm, mtarget=img_fm, nrow=nex, ncol=nex, s_firstrow=1, &
                              s_firstcol=1, t_firstrow=2, t_firstcol=2)

!  Then spin-flip/spin-flip SOC
!  <Psi_Isf|SOC|Psi_Jsf> =
!  sum_k,sigma [<psi^Isf_k,tau|SOC|psi^Jsf_k,tau> + <psi^Isf_k,tau|psi^Jsf_k,tau> * gs_sum]
!  - sum_k,l,sigma <psi^0_k,sigma|SOC|psi^0_l,sigma> * <psi^Isf_l,tau|psi^Jsf_k,tau> , tau != sigma

      !Same spin integration => only SOC_z matters, and the contribution is purely imaginary
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_z, dbcsr_sf, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sf, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      !the overlap as well
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_sf, 0.0_dp, &
                          dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sf, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

      !note: the different prefactors are derived from the fact that because of spin-flip, we have
      !alpha-gs and beta-lr interaction
      CALL os_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_soc_z, pref_diaga=-1.0_dp, &
                                pref_diagb=1.0_dp, pref_tracea=-1.0_dp, pref_traceb=1.0_dp, &
                                pref_diags=gs_sum, symmetric=.TRUE.)

      CALL copy_dbcsr_to_fm(dbcsr_tmp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm, mtarget=img_fm, nrow=nex, ncol=nex, s_firstrow=1, &
                              s_firstcol=1, t_firstrow=1 + nsc + 1, t_firstcol=1 + nsc + 1)

!  Finally the spin-conserving/spin-flip interaction
! <Psi_Isc|SOC|Psi_Jsf> =   sum_k,sigma <psi^Isc_k,sigma|SOC|psi^Isf_k,tau>
!                           - sum_k,l,sigma <psi^0_k,tau|SOC|psi^0_l,sigma

      tas(1) = 1; tbs(1) = ndo_mo + 1
      tas(2) = ndo_mo + 1; tbs(2) = 1

      !the overlap
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_sf, 0.0_dp, &
                          dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

      !start with the imaginary contribution
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_x, dbcsr_sf, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      CALL os_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_soc_x, pref_diaga=1.0_dp, &
                                pref_diagb=1.0_dp, pref_tracea=-1.0_dp, pref_traceb=-1.0_dp, &
                                tracea_start=tas, traceb_start=tbs)

      CALL copy_dbcsr_to_fm(dbcsr_tmp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm, mtarget=img_fm, nrow=nex, ncol=nex, s_firstrow=1, &
                              s_firstcol=1, t_firstrow=2, t_firstcol=1 + nsc + 1)

      !follow up with the real contribution
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_y, dbcsr_sf, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      CALL os_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_soc_y, pref_diaga=1.0_dp, &
                                pref_diagb=-1.0_dp, pref_tracea=1.0_dp, pref_traceb=-1.0_dp, &
                                tracea_start=tas, traceb_start=tbs)

      CALL copy_dbcsr_to_fm(dbcsr_tmp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm, mtarget=real_fm, nrow=nex, ncol=nex, s_firstrow=1, &
                              s_firstcol=1, t_firstrow=2, t_firstcol=1 + nsc + 1)

!  Setting up the complex Hermitian perturbed matrix
      CALL cp_cfm_create(pert_cfm, full_struct)
      CALL cp_fm_to_cfm(real_fm, img_fm, pert_cfm)

      CALL cp_fm_release(real_fm)
      CALL cp_fm_release(img_fm)

!  Diagonalize the perturbed matrix
      ALLOCATE (tmp_evals(ntot))
      CALL cp_cfm_create(evecs_cfm, full_struct)
      CALL cp_cfm_heevd(pert_cfm, evecs_cfm, tmp_evals)

      !shift the energies such that the GS has zero and store all that in soc_evals (\wo the GS)
      ALLOCATE (donor_state%soc_evals(ntot - 1))
      donor_state%soc_evals(:) = tmp_evals(2:ntot) - tmp_evals(1)

!  The SOC dipole oscillator strengths
      CALL compute_soc_dipole_fosc(evecs_cfm, dbcsr_soc_package, donor_state, xas_tdp_env, &
                                   xas_tdp_control, qs_env, gs_coeffs=gs_coeffs)

!  And quadrupole
      IF (xas_tdp_control%do_quad) THEN
         CALL compute_soc_quadrupole_fosc(evecs_cfm, dbcsr_soc_package, donor_state, xas_tdp_env, &
                                          xas_tdp_control, qs_env, gs_coeffs=gs_coeffs)
      END IF

! Clean-up
      CALL cp_cfm_release(pert_cfm)
      CALL cp_cfm_release(evecs_cfm)
      CALL cp_fm_struct_release(full_struct)
      IF (do_roks) CALL cp_fm_release(gs_coeffs)
      CALL dbcsr_distribution_release(coeffs_dist)
      CALL dbcsr_distribution_release(prod_dist)
      CALL dbcsr_release(dbcsr_sc)
      CALL dbcsr_release(dbcsr_sf)
      CALL dbcsr_release(dbcsr_prod)
      CALL dbcsr_release(dbcsr_ovlp)
      CALL dbcsr_release(dbcsr_tmp)
      CALL dbcsr_release(dbcsr_work)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_struct_release(work_struct)
      DEALLOCATE (coeffs_dist, prod_dist, col_dist, col_blk_size, row_dist_new)
      DEALLOCATE (dbcsr_sc, dbcsr_sf, dbcsr_work, dbcsr_prod, dbcsr_ovlp, dbcsr_tmp)

      CALL timestop(handle)

   END SUBROUTINE include_os_soc

! **************************************************************************************************
!> \brief Includes the SOC effects on the precomputed restricted closed-shell singlet and triplet
!>        excitations. This is a perturbative treatmnent
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note Using AMEWs, build an hermitian matrix with all excited states SOC coupling + the
!>       excitation energies on the diagonal. Then diagonalize it to get the new excitation
!>       energies and corresponding linear combinations of lr_coeffs.
!>       The AMEWs are normalized
!>       Only for spin-restricted calculations
!>       The ms=-1,+1 triplets are not explicitely computed in the first place. Assume they have
!>       the same energy as the ms=0 triplets and apply the spin raising and lowering operators
!>       on the latter to get their AMEWs => this is the qusi-degenerate perturbation theory
!>       approach by Neese (QDPT)
! **************************************************************************************************
   SUBROUTINE include_rcs_soc(donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'include_rcs_soc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: group, handle, iex, isg, itp, nao, &
                                                            ndo_mo, nex, npcols, nprows, nsg, &
                                                            ntot, ntp
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_dist, row_blk_size, &
                                                            row_dist, row_dist_new
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      REAL(dp)                                           :: eps_filter, soc_gst, sqrt2
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag, tmp_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: domo_soc_x, domo_soc_y, domo_soc_z, &
                                                            gstp_block
      REAL(dp), DIMENSION(:), POINTER                    :: sg_evals, tp_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type), POINTER                         :: evecs_cfm, hami_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: full_struct, gstp_struct, prod_struct, &
                                                            vec_struct, work_struct
      TYPE(cp_fm_type), POINTER                          :: gs_coeffs, gstp_fm, img_fm, prod_fm, &
                                                            real_fm, sg_coeffs, tmp_fm, tp_coeffs, &
                                                            vec_soc_x, vec_soc_y, vec_soc_z, &
                                                            work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: coeffs_dist, dbcsr_dist, prod_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_soc_package_type)                       :: dbcsr_soc_package
      TYPE(dbcsr_type), POINTER                          :: dbcsr_ovlp, dbcsr_prod, dbcsr_sg, &
                                                            dbcsr_tmp, dbcsr_tp, dbcsr_work, &
                                                            orb_soc_x, orb_soc_y, orb_soc_z

      NULLIFY (sg_coeffs, tp_coeffs, gs_coeffs, sg_evals, tp_evals, real_fm, img_fm, full_struct)
      NULLIFY (para_env, blacs_env, prod_fm, prod_struct, vec_struct, orb_soc_y, orb_soc_z)
      NULLIFY (matrix_s, orb_soc_x, hami_cfm, evecs_cfm, vec_soc_x, vec_soc_y, vec_soc_z)
      NULLIFY (work_struct, work_fm, tmp_fm, dbcsr_dist, coeffs_dist, prod_dist, pgrid)
      NULLIFY (col_dist, row_dist, col_blk_size, row_blk_size, row_dist_new, gstp_struct)
      NULLIFY (dbcsr_tp, dbcsr_sg, dbcsr_prod, dbcsr_work, dbcsr_ovlp, dbcsr_tmp)

      CALL timeset(routineN, handle)

!  Initialization
      CPASSERT(ASSOCIATED(xas_tdp_control))
      gs_coeffs => donor_state%gs_coeffs
      sg_coeffs => donor_state%sg_coeffs
      tp_coeffs => donor_state%tp_coeffs
      sg_evals => donor_state%sg_evals
      tp_evals => donor_state%tp_evals
      nsg = SIZE(sg_evals)
      ntp = SIZE(tp_evals)
      ntot = 1 + nsg + 3*ntp
      ndo_mo = donor_state%ndo_mo
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      orb_soc_x => xas_tdp_env%orb_soc(1)%matrix
      orb_soc_y => xas_tdp_env%orb_soc(2)%matrix
      orb_soc_z => xas_tdp_env%orb_soc(3)%matrix
      !by construction nsg == ntp, keep those separate for more code clarity though
      CPASSERT(nsg == ntp)
      nex = nsg
      eps_filter = xas_tdp_control%eps_filter

!  Creating the real part and imaginary part of the final SOC fm
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      CALL cp_fm_create(real_fm, full_struct)
      CALL cp_fm_create(img_fm, full_struct)

!  Put the excitation energies on the diagonal of the real matrix
      DO isg = 1, nsg
         CALL cp_fm_set_element(real_fm, 1 + isg, 1 + isg, sg_evals(isg))
      END DO
      DO itp = 1, ntp
         ! first T^-1, then T^0, then T^+1
         CALL cp_fm_set_element(real_fm, 1 + itp + nsg, 1 + itp + nsg, tp_evals(itp))
         CALL cp_fm_set_element(real_fm, 1 + itp + ntp + nsg, 1 + itp + ntp + nsg, tp_evals(itp))
         CALL cp_fm_set_element(real_fm, 1 + itp + 2*ntp + nsg, 1 + itp + 2*ntp + nsg, tp_evals(itp))
      END DO

!  Create the dbcsr machinery (for fast MM, the core of this routine)
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_distribution_get(dbcsr_dist, group=group, row_dist=row_dist, pgrid=pgrid, &
                                  npcols=npcols, nprows=nprows)
      ALLOCATE (col_dist(nex), row_dist_new(nex))
      DO iex = 1, nex
         col_dist(iex) = MODULO(npcols - iex, npcols)
         row_dist_new(iex) = MODULO(nprows - iex, nprows)
      END DO
      ALLOCATE (coeffs_dist, prod_dist)
      CALL dbcsr_distribution_new(coeffs_dist, group=group, pgrid=pgrid, row_dist=row_dist, &
                                  col_dist=col_dist)
      CALL dbcsr_distribution_new(prod_dist, group=group, pgrid=pgrid, row_dist=row_dist_new, &
                                  col_dist=col_dist)

      !Create the matrices
      ALLOCATE (col_blk_size(nex))
      col_blk_size = ndo_mo
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size)

      ALLOCATE (dbcsr_sg, dbcsr_tp, dbcsr_work, dbcsr_ovlp, dbcsr_tmp, dbcsr_prod)
      CALL dbcsr_create(matrix=dbcsr_sg, name="SINGLETS", matrix_type="N", dist=coeffs_dist, &
                        row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_tp, name="TRIPLETS", matrix_type="N", dist=coeffs_dist, &
                        row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_work, name="WORK", matrix_type="N", dist=coeffs_dist, &
                        row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_prod, name="PROD", matrix_type="N", dist=prod_dist, &
                        row_blk_size=col_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_create(matrix=dbcsr_ovlp, name="OVLP", matrix_type="N", dist=prod_dist, &
                        row_blk_size=col_blk_size, col_blk_size=col_blk_size)

      col_blk_size = 1
      CALL dbcsr_create(matrix=dbcsr_tmp, name="TMP", matrix_type="N", dist=prod_dist, &
                        row_blk_size=col_blk_size, col_blk_size=col_blk_size)
      CALL dbcsr_reserve_all_blocks(dbcsr_tmp)

      dbcsr_soc_package%dbcsr_sg => dbcsr_sg
      dbcsr_soc_package%dbcsr_tp => dbcsr_tp
      dbcsr_soc_package%dbcsr_work => dbcsr_work
      dbcsr_soc_package%dbcsr_ovlp => dbcsr_ovlp
      dbcsr_soc_package%dbcsr_prod => dbcsr_prod
      dbcsr_soc_package%dbcsr_tmp => dbcsr_tmp

      !Filling the coeffs matrices by copying from the stored fms
      CALL copy_fm_to_dbcsr(sg_coeffs, dbcsr_sg)
      CALL copy_fm_to_dbcsr(tp_coeffs, dbcsr_tp)

!  Create the work and helper fms
      CALL cp_fm_get_info(gs_coeffs, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_mo, ncol_global=ndo_mo)
      CALL cp_fm_create(prod_fm, prod_struct)
      CALL cp_fm_create(vec_soc_x, vec_struct)
      CALL cp_fm_create(vec_soc_y, vec_struct)
      CALL cp_fm_create(vec_soc_z, vec_struct)
      CALL cp_fm_struct_create(work_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nex, ncol_global=nex)
      CALL cp_fm_create(work_fm, work_struct)
      CALL cp_fm_create(tmp_fm, work_struct)
      ALLOCATE (diag(ndo_mo))

!  Precompute everything we can before looping over excited states
      sqrt2 = SQRT(2.0_dp)

      ! The subset of the donor MOs matrix elements: <phi_I^0|Hsoc|phi_J^0> (kept as global array, small)
      ALLOCATE (domo_soc_x(ndo_mo, ndo_mo), domo_soc_y(ndo_mo, ndo_mo), domo_soc_z(ndo_mo, ndo_mo))

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_x, gs_coeffs, vec_soc_x, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_x, 0.0_dp, prod_fm)
      CALL cp_fm_get_submatrix(prod_fm, domo_soc_x)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, gs_coeffs, vec_soc_y, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_y, 0.0_dp, prod_fm)
      CALL cp_fm_get_submatrix(prod_fm, domo_soc_y)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, gs_coeffs, vec_soc_z, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_z, 0.0_dp, prod_fm)
      CALL cp_fm_get_submatrix(prod_fm, domo_soc_z)

!  Only have SOC between singlet-triplet triplet-triplet and ground_state-triplet, the resulting
!  matrix is Hermitian i.e. the real part is symmetric and the imaginary part is anti-symmetric.
!  Can only fill upper half

      !Start with the ground state/triplet SOC, SOC*gs_coeffs already computed above
      !note: we are computing <0|H|T>, but have SOC*gs_coeffs instead of gs_coeffs*SOC in store. Since
      !      the SOC Hamiltonian is anti-symmetric, a - signs pops up in the gemms below

      CALL cp_fm_struct_create(gstp_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntp*ndo_mo, ncol_global=ndo_mo)
      CALL cp_fm_create(gstp_fm, gstp_struct)
      ALLOCATE (gstp_block(ndo_mo, ndo_mo))

      !gs-triplet with Ms=+-1, imaginary part
      CALL cp_gemm('T', 'N', ndo_mo*ntp, ndo_mo, nao, -1.0_dp, tp_coeffs, vec_soc_x, 0.0_dp, gstp_fm)

      DO itp = 1, ntp
         CALL cp_fm_get_submatrix(fm=gstp_fm, target_m=gstp_block, start_row=(itp - 1)*ndo_mo + 1, &
                                  start_col=1, n_rows=ndo_mo, n_cols=ndo_mo)
         diag(:) = get_diag(gstp_block)
         soc_gst = SUM(diag)
         CALL cp_fm_set_element(img_fm, 1, 1 + nsg + itp, -1.0_dp*soc_gst) ! <0|H_x|T^-1>
         CALL cp_fm_set_element(img_fm, 1, 1 + nsg + 2*ntp + itp, soc_gst) ! <0|H_x|T^+1>
      END DO

      !gs-triplet with Ms=+-1, real part
      CALL cp_gemm('T', 'N', ndo_mo*ntp, ndo_mo, nao, -1.0_dp, tp_coeffs, vec_soc_y, 0.0_dp, gstp_fm)

      DO itp = 1, ntp
         CALL cp_fm_get_submatrix(fm=gstp_fm, target_m=gstp_block, start_row=(itp - 1)*ndo_mo + 1, &
                                  start_col=1, n_rows=ndo_mo, n_cols=ndo_mo)
         diag(:) = get_diag(gstp_block)
         soc_gst = SUM(diag)
         CALL cp_fm_set_element(real_fm, 1, 1 + nsg + itp, -1.0_dp*soc_gst) ! <0|H_y|T^-1>
         CALL cp_fm_set_element(real_fm, 1, 1 + nsg + 2*ntp + itp, -1.0_dp*soc_gst) ! <0|H_y|T^+1>
      END DO

      !gs-triplet with Ms=0, purely imaginary
      CALL cp_gemm('T', 'N', ndo_mo*ntp, ndo_mo, nao, -1.0_dp, tp_coeffs, vec_soc_z, 0.0_dp, gstp_fm)

      DO itp = 1, ntp
         CALL cp_fm_get_submatrix(fm=gstp_fm, target_m=gstp_block, start_row=(itp - 1)*ndo_mo + 1, &
                                  start_col=1, n_rows=ndo_mo, n_cols=ndo_mo)
         diag(:) = get_diag(gstp_block)
         soc_gst = sqrt2*SUM(diag)
         CALL cp_fm_set_element(img_fm, 1, 1 + nsg + ntp + itp, soc_gst)
      END DO

      !gs clean-up
      CALL cp_fm_release(prod_fm)
      CALL cp_fm_release(vec_soc_x)
      CALL cp_fm_release(vec_soc_y)
      CALL cp_fm_release(vec_soc_z)
      CALL cp_fm_release(gstp_fm)
      CALL cp_fm_struct_release(gstp_struct)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (gstp_block)

      !Now do the singlet-triplet SOC
      !start by computing the singlet-triplet overlap
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_tp, 0.0_dp, &
                          dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sg, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

      !singlet-triplet with Ms=+-1, imaginary part
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_x, dbcsr_tp, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sg, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_soc_x, &
                                 pref_trace=-1.0_dp, pref_overall=-0.5_dp*sqrt2)

      !<S|H_x|T^-1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=img_fm, nrow=nex, ncol=nex, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=2, &
                              t_firstcol=1 + nsg + 1)

      !<S|H_x|T^+1> takes a minus sign
      CALL cp_fm_scale(-1.0_dp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=img_fm, nrow=nex, ncol=nex, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=2, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      !singlet-triplet with Ms=+-1, real part
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_y, dbcsr_tp, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sg, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_soc_y, &
                                 pref_trace=-1.0_dp, pref_overall=-0.5_dp*sqrt2)

      !<S|H_y|T^-1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=real_fm, nrow=nex, ncol=nex, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=2, &
                              t_firstcol=1 + nsg + 1)

      !<S|H_y|T^+1>
      CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=real_fm, nrow=nex, ncol=nex, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=2, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      !singlet-triplet with Ms=0, purely imaginary
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_z, dbcsr_tp, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sg, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_soc_z, &
                                 pref_trace=-1.0_dp, pref_overall=1.0_dp)

      !<S|H_z|T^0>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=img_fm, nrow=nex, ncol=nex, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=2, &
                              t_firstcol=1 + nsg + ntp + 1)

      !Now the triplet-triplet SOC
      !start by computing the overlap
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_tp, 0.0_dp, &
                          dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_tp, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

      !Ms=0 to Ms=+-1 SOC, imaginary part
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_x, dbcsr_tp, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_tp, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_soc_x, &
                                 pref_trace=1.0_dp, pref_overall=-0.5_dp*sqrt2)

      !<T^0|H_x|T^+1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=img_fm, nrow=nex, ncol=nex, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=1 + nsg + ntp + 1, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      !<T^-1|H_x|T^0>, takes a minus sign and a transpose (because computed <T^0|H_x|T^-1>)
      CALL cp_fm_transpose(tmp_fm, work_fm)
      CALL cp_fm_scale(-1.0_dp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm, mtarget=img_fm, nrow=nex, ncol=nex, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=1 + nsg + 1, &
                              t_firstcol=1 + nsg + ntp + 1)

      !Ms=0 to Ms=+-1 SOC, real part
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_y, dbcsr_tp, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_tp, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_soc_y, &
                                 pref_trace=1.0_dp, pref_overall=0.5_dp*sqrt2)

      !<T^0|H_y|T^+1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=real_fm, nrow=nex, ncol=nex, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=1 + nsg + ntp + 1, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      !<T^-1|H_y|T^0>, takes a minus sign and a transpose
      CALL cp_fm_transpose(tmp_fm, work_fm)
      CALL cp_fm_scale(-1.0_dp, work_fm)
      CALL cp_fm_to_fm_submat(msource=work_fm, mtarget=real_fm, nrow=nex, ncol=nex, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=1 + nsg + 1, &
                              t_firstcol=1 + nsg + ntp + 1)

      !Ms=1 to Ms=1 and Ms=-1 to Ms=-1 SOC, purely imaginary
      CALL dbcsr_multiply('N', 'N', 1.0_dp, orb_soc_z, dbcsr_tp, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_tp, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

      CALL rcs_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_soc_z, &
                                 pref_trace=1.0_dp, pref_overall=1.0_dp)

      !<T^+1|H_z|T^+1>
      CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=img_fm, nrow=nex, ncol=nex, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=1 + nsg + 2*ntp + 1, &
                              t_firstcol=1 + nsg + 2*ntp + 1)

      !<T^-1|H_z|T^-1>, takes a minus sign
      CALL cp_fm_scale(-1.0_dp, tmp_fm)
      CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=img_fm, nrow=nex, ncol=nex, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=1 + nsg + 1, &
                              t_firstcol=1 + nsg + 1)

!  Intermediate clean-up
      CALL cp_fm_struct_release(work_struct)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(tmp_fm)
      DEALLOCATE (diag, domo_soc_x, domo_soc_y, domo_soc_z)

!  Set-up the complex hermitian perturbation matrix
      CALL cp_cfm_create(hami_cfm, full_struct)
      CALL cp_fm_to_cfm(real_fm, img_fm, hami_cfm)

      CALL cp_fm_release(real_fm)
      CALL cp_fm_release(img_fm)

!  Diagonalize the Hamiltonian
      ALLOCATE (tmp_evals(ntot))
      CALL cp_cfm_create(evecs_cfm, full_struct)
      CALL cp_cfm_heevd(hami_cfm, evecs_cfm, tmp_evals)

      !  Adjust the energies so the GS has zero, and store in the donor_state (without the GS)
      ALLOCATE (donor_state%soc_evals(ntot - 1))
      donor_state%soc_evals(:) = tmp_evals(2:ntot) - tmp_evals(1)

!  Compute the dipole oscillator strengths
      CALL compute_soc_dipole_fosc(evecs_cfm, dbcsr_soc_package, donor_state, xas_tdp_env, &
                                   xas_tdp_control, qs_env)

!  And the quadrupole (if needed)
      IF (xas_tdp_control%do_quad) THEN
         CALL compute_soc_quadrupole_fosc(evecs_cfm, dbcsr_soc_package, donor_state, xas_tdp_env, &
                                          xas_tdp_control, qs_env)
      END IF

!  Clean-up
      CALL cp_fm_struct_release(full_struct)
      CALL cp_cfm_release(hami_cfm)
      CALL cp_cfm_release(evecs_cfm)
      CALL dbcsr_distribution_release(coeffs_dist)
      CALL dbcsr_distribution_release(prod_dist)
      CALL dbcsr_release(dbcsr_sg)
      CALL dbcsr_release(dbcsr_tp)
      CALL dbcsr_release(dbcsr_prod)
      CALL dbcsr_release(dbcsr_ovlp)
      CALL dbcsr_release(dbcsr_tmp)
      CALL dbcsr_release(dbcsr_work)
      DEALLOCATE (coeffs_dist, prod_dist, col_dist, col_blk_size, row_dist_new)
      DEALLOCATE (dbcsr_sg, dbcsr_tp, dbcsr_work, dbcsr_prod, dbcsr_ovlp, dbcsr_tmp)

      CALL timestop(handle)

   END SUBROUTINE include_rcs_soc

! **************************************************************************************************
!> \brief Computes the matrix elements of a one-body operator (given wrt AOs) in the basis of the
!>        excited state AMEWs with ground state, for the open-shell case
!> \param amew_op the operator in the basis of the AMEWs (array because could have x,y,z components)
!> \param ao_op the operator in the basis of the atomic orbitals
!> \param gs_coeffs the coefficient of the GS donor MOs. Ecplicitely passed because of special
!>                  format in the ROKS case (see include_os_soc routine)
!> \param dbcsr_soc_package inhertited from the main SOC routine
!> \param donor_state ...
!> \param eps_filter ...
!> \param qs_env ...
!> \note The ordering of the AMEWs is consistent with SOC and is gs, sc, sf
!>       We assume that the operator is spin-independent => only <0|0>, <0|sc>, <sc|sc> and <sf|sf>
!>       yield non-zero matrix elements
!>       Only for open-shell calculations
! **************************************************************************************************
   SUBROUTINE get_os_amew_op(amew_op, ao_op, gs_coeffs, dbcsr_soc_package, donor_state, &
                             eps_filter, qs_env)

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: amew_op
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ao_op
      TYPE(cp_fm_type), POINTER                          :: gs_coeffs
      TYPE(dbcsr_soc_package_type)                       :: dbcsr_soc_package
      TYPE(donor_state_type), POINTER                    :: donor_state
      REAL(dp), INTENT(IN)                               :: eps_filter
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'get_os_amew_op', routineP = moduleN//':'//routineN

      INTEGER                                            :: dim_op, homo, i, isc, nao, ndo_mo, &
                                                            ndo_so, nex, nsc, nsf, ntot
      REAL(dp)                                           :: op
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag, gsgs_op
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: domo_op, gsex_block, tmp
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: full_struct, gsex_struct, prod_struct, &
                                                            tmp_struct, vec_struct
      TYPE(cp_fm_type), POINTER                          :: amew_op_i, gsex_fm, mo_coeff, prod_work, &
                                                            sc_coeffs, sf_coeffs, tmp_fm, &
                                                            vec_work, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: ao_op_i, dbcsr_ovlp, dbcsr_prod, &
                                                            dbcsr_sc, dbcsr_sf, dbcsr_tmp, &
                                                            dbcsr_work
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      NULLIFY (matrix_s, para_env, blacs_env, full_struct, vec_struct, prod_struct, mos, vec_work)
      NULLIFY (prod_work, mo_coeff, ao_op_i, amew_op_i, work_fm, tmp_fm, tmp_struct)
      NULLIFY (dbcsr_sc, dbcsr_sf, dbcsr_ovlp, dbcsr_work, dbcsr_tmp, dbcsr_prod)

!  Iinitialization
      dim_op = SIZE(ao_op)
      sc_coeffs => donor_state%sc_coeffs
      sf_coeffs => donor_state%sf_coeffs
      nsc = SIZE(donor_state%sc_evals)
      nsf = SIZE(donor_state%sf_evals)
      nex = nsc
      ntot = 1 + nsc + nsf
      ndo_mo = donor_state%ndo_mo
      ndo_so = 2*donor_state%ndo_mo !open-shell => nspins = 2
      CALL get_qs_env(qs_env, matrix_s=matrix_s, para_env=para_env, blacs_env=blacs_env, mos=mos)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

      dbcsr_sc => dbcsr_soc_package%dbcsr_sc
      dbcsr_sf => dbcsr_soc_package%dbcsr_sf
      dbcsr_work => dbcsr_soc_package%dbcsr_work
      dbcsr_tmp => dbcsr_soc_package%dbcsr_tmp
      dbcsr_prod => dbcsr_soc_package%dbcsr_prod
      dbcsr_ovlp => dbcsr_soc_package%dbcsr_ovlp

!  Create the amew_op matrix set
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      ALLOCATE (amew_op(dim_op))
      DO i = 1, dim_op
         CALL cp_fm_create(amew_op(i)%matrix, full_struct)
      END DO

!  Before looping, need to evaluate sum_j,sigma <phi^0_j,sgima|op|phi^0_j,sigma>, for each dimension
!  of the operator
      ALLOCATE (gsgs_op(dim_op))

      !start with the alpha MOs
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, homo=homo)
      ALLOCATE (diag(homo))
      CALL cp_fm_get_info(mo_coeff, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=homo, ncol_global=homo)
      CALL cp_fm_create(vec_work, vec_struct)
      CALL cp_fm_create(prod_work, prod_struct)

      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix

         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, mo_coeff, vec_work, ncol=homo)
         CALL cp_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, vec_work, 0.0_dp, prod_work)
         CALL cp_fm_get_diag(prod_work, diag)
         gsgs_op(i) = SUM(diag)

      END DO !i

      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (diag)
      NULLIFY (vec_struct)

      !then beta orbitals
      CALL get_mo_set(mos(2)%mo_set, mo_coeff=mo_coeff, homo=homo)
      ALLOCATE (diag(homo))
      CALL cp_fm_get_info(mo_coeff, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=homo, ncol_global=homo)
      CALL cp_fm_create(vec_work, vec_struct)
      CALL cp_fm_create(prod_work, prod_struct)

      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix

         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, mo_coeff, vec_work, ncol=homo)
         CALL cp_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, vec_work, 0.0_dp, prod_work)
         CALL cp_fm_get_diag(prod_work, diag)
         gsgs_op(i) = gsgs_op(i) + SUM(diag)

      END DO !i

      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (diag)
      NULLIFY (vec_struct)

!  Before looping over excited AMEWs, define some work matrices and structures
      CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nao, ncol_global=ndo_so)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_so, ncol_global=ndo_so)
      CALL cp_fm_struct_create(gsex_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_so*nex, ncol_global=ndo_so)
      CALL cp_fm_struct_create(tmp_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nex, ncol_global=nex)
      CALL cp_fm_create(vec_work, vec_struct) !for op*|phi>
      CALL cp_fm_create(prod_work, prod_struct) !for any <phi|op|phi>
      CALL cp_fm_create(work_fm, full_struct)
      CALL cp_fm_create(gsex_fm, gsex_struct)
      CALL cp_fm_create(tmp_fm, tmp_struct)
      ALLOCATE (diag(ndo_so))
      ALLOCATE (domo_op(ndo_so, ndo_so))
      ALLOCATE (tmp(ndo_so, ndo_so))
      ALLOCATE (gsex_block(ndo_so, ndo_so))

!  Loop over the dimensions of the operator
      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix
         amew_op_i => amew_op(i)%matrix

         !put the gs-gs contribution
         CALL cp_fm_set_element(amew_op_i, 1, 1, gsgs_op(i))

         !  Precompute what we can before looping over excited states
         ! Need the operator in the donor MOs basis <phi^0_I,sigma|op_i|phi^0_J,tau>
         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, gs_coeffs, vec_work, ncol=ndo_so)
         CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, vec_work, 0.0_dp, prod_work)
         CALL cp_fm_get_submatrix(prod_work, domo_op)

         !  Do the ground-state/spin-conserving operator
         CALL cp_gemm('T', 'N', ndo_so*nsc, ndo_so, nao, 1.0_dp, sc_coeffs, vec_work, 0.0_dp, gsex_fm)
         DO isc = 1, nsc
            CALL cp_fm_get_submatrix(fm=gsex_fm, target_m=gsex_block, start_row=(isc - 1)*ndo_so + 1, &
                                     start_col=1, n_rows=ndo_so, n_cols=ndo_so)
            diag(:) = get_diag(gsex_block)
            op = SUM(diag)
            CALL cp_fm_set_element(amew_op_i, 1, 1 + isc, op)
         END DO !isc

         !  The spin-conserving/spin-conserving operator
         !overlap
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_sc, 0.0_dp, &
                             dbcsr_work, filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

         !operator in SC LR-orbital basis
         CALL dbcsr_multiply('N', 'N', 1.0_dp, ao_op_i, dbcsr_sc, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sc, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

         CALL os_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_op, pref_diaga=1.0_dp, &
                                   pref_diagb=1.0_dp, pref_tracea=-1.0_dp, pref_traceb=-1.0_dp, &
                                   pref_diags=gsgs_op(i), symmetric=.TRUE.)

         CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
         CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=amew_op_i, nrow=nex, ncol=nex, &
                                 s_firstrow=1, s_firstcol=1, t_firstrow=2, t_firstcol=2)

         !  The spin-flip/spin-flip operator
         !overlap
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_sf, 0.0_dp, &
                             dbcsr_work, filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sf, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

         !operator in SF LR-orbital basis
         CALL dbcsr_multiply('N', 'N', 1.0_dp, ao_op_i, dbcsr_sf, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sf, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

         !need to reorganize the domo_op array by swapping the alpha-alpha and the beta-beta quarter
         tmp(1:ndo_mo, 1:ndo_mo) = domo_op(ndo_mo + 1:ndo_so, ndo_mo + 1:ndo_so)
         tmp(ndo_mo + 1:ndo_so, ndo_mo + 1:ndo_so) = domo_op(1:ndo_mo, 1:ndo_mo)

         CALL os_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, tmp, pref_diaga=1.0_dp, &
                                   pref_diagb=1.0_dp, pref_tracea=-1.0_dp, pref_traceb=-1.0_dp, &
                                   pref_diags=gsgs_op(i), symmetric=.TRUE.)

         CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
         CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=amew_op_i, nrow=nex, ncol=nex, &
                                 s_firstrow=1, s_firstcol=1, t_firstrow=1 + nsc + 1, t_firstcol=1 + nsc + 1)

         !Symmetry => only upper diag explicitly built
         CALL cp_fm_upper_to_full(amew_op_i, work_fm)

      END DO !i

!  Clean-up
      CALL cp_fm_struct_release(full_struct)
      CALL cp_fm_struct_release(prod_struct)
      CALL cp_fm_struct_release(vec_struct)
      CALL cp_fm_struct_release(tmp_struct)
      CALL cp_fm_struct_release(gsex_struct)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(tmp_fm)
      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_release(gsex_fm)

   END SUBROUTINE get_os_amew_op

! **************************************************************************************************
!> \brief Computes the matrix elements of a one-body operator (given wrt AOs) in the basis of the
!>        excited state AMEWs with ground state, singlet and triplet with Ms = -1,0,+1
!> \param amew_op the operator in the basis of the AMEWs (array because could have x,y,z components)
!> \param ao_op the operator in the basis of the atomic orbitals
!> \param dbcsr_soc_package inherited from the main SOC routine
!> \param donor_state ...
!> \param eps_filter for dbcsr multiplication
!> \param qs_env ...
!> \note The ordering of the AMEWs is consistent with SOC and is gs, sg, tp(-1), tp(0). tp(+1)
!>       We assume that the operator is spin-independent => only <0|0>, <0|S>, <S|S> and <T|T>
!>       yield non-zero matrix elements
!>       Only for spin-restricted calculations
! **************************************************************************************************
   SUBROUTINE get_rcs_amew_op(amew_op, ao_op, dbcsr_soc_package, donor_state, eps_filter, qs_env)

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: amew_op
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ao_op
      TYPE(dbcsr_soc_package_type)                       :: dbcsr_soc_package
      TYPE(donor_state_type), POINTER                    :: donor_state
      REAL(dp), INTENT(IN)                               :: eps_filter
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'get_rcs_amew_op', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: dim_op, homo, i, isg, nao, ndo_mo, nex, &
                                                            nsg, ntot, ntp
      REAL(dp)                                           :: op, sqrt2
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag, gs_diag, gsgs_op
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: domo_op, sggs_block
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: full_struct, gsgs_struct, prod_struct, &
                                                            sggs_struct, std_struct, tmp_struct, &
                                                            vec_struct
      TYPE(cp_fm_type), POINTER                          :: amew_op_i, gs_coeffs, gs_fm, mo_coeff, &
                                                            prod_fm, sg_coeffs, sggs_fm, tmp_fm, &
                                                            vec_op, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: ao_op_i, dbcsr_ovlp, dbcsr_prod, &
                                                            dbcsr_sg, dbcsr_tmp, dbcsr_tp, &
                                                            dbcsr_work
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      NULLIFY (gs_coeffs, sg_coeffs, tmp_fm, matrix_s, full_struct, prod_fm, prod_struct, work_fm)
      NULLIFY (vec_struct, blacs_env, para_env, mo_coeff, mos, gsgs_struct, std_struct)
      NULLIFY (vec_op, gs_fm, tmp_struct, sggs_struct, sggs_fm)
      NULLIFY (ao_op_i, dbcsr_tp, dbcsr_sg, dbcsr_ovlp, dbcsr_work, dbcsr_tmp, dbcsr_prod)

!  Initialization
      gs_coeffs => donor_state%gs_coeffs
      sg_coeffs => donor_state%sg_coeffs
      nsg = SIZE(donor_state%sg_evals)
      ntp = nsg; nex = nsg !all the same by construction, keep them separate for clarity
      ntot = 1 + nsg + 3*ntp
      ndo_mo = donor_state%ndo_mo
      CALL get_qs_env(qs_env, matrix_s=matrix_s, para_env=para_env, blacs_env=blacs_env, mos=mos)
      sqrt2 = SQRT(2.0_dp)
      dim_op = SIZE(ao_op)

      dbcsr_sg => dbcsr_soc_package%dbcsr_sg
      dbcsr_tp => dbcsr_soc_package%dbcsr_tp
      dbcsr_work => dbcsr_soc_package%dbcsr_work
      dbcsr_prod => dbcsr_soc_package%dbcsr_prod
      dbcsr_ovlp => dbcsr_soc_package%dbcsr_ovlp
      dbcsr_tmp => dbcsr_soc_package%dbcsr_tmp

!  Create the amew_op matrix
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      ALLOCATE (amew_op(dim_op))
      DO i = 1, dim_op
         CALL cp_fm_create(amew_op(i)%matrix, full_struct)
      END DO !i

!  Deal with the GS-GS contribution <0|0> = 2*sum_j <phi_j|op|phi_j>
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao, homo=homo)
      CALL cp_fm_struct_create(gsgs_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=homo, ncol_global=homo)
      CALL cp_fm_get_info(mo_coeff, matrix_struct=std_struct)
      CALL cp_fm_create(gs_fm, gsgs_struct)
      CALL cp_fm_create(work_fm, std_struct)
      ALLOCATE (gsgs_op(dim_op))
      ALLOCATE (gs_diag(homo))

      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix

         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, mo_coeff, work_fm, ncol=homo)
         CALL cp_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, work_fm, 0.0_dp, gs_fm)
         CALL cp_fm_get_diag(gs_fm, gs_diag)
         gsgs_op(i) = 2.0_dp*SUM(gs_diag)

      END DO !i

      CALL cp_fm_release(gs_fm)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_struct_release(gsgs_struct)
      DEALLOCATE (gs_diag)

!  Create the work and helper fms
      CALL cp_fm_get_info(gs_coeffs, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_mo, ncol_global=ndo_mo)
      CALL cp_fm_create(prod_fm, prod_struct)
      CALL cp_fm_create(vec_op, vec_struct)
      CALL cp_fm_struct_create(tmp_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nex, ncol_global=nex)
      CALL cp_fm_struct_create(sggs_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_mo*nsg, ncol_global=ndo_mo)
      CALL cp_fm_create(tmp_fm, tmp_struct)
      CALL cp_fm_create(work_fm, full_struct)
      CALL cp_fm_create(sggs_fm, sggs_struct)
      ALLOCATE (diag(ndo_mo))
      ALLOCATE (domo_op(ndo_mo, ndo_mo))
      ALLOCATE (sggs_block(ndo_mo, ndo_mo))

! Iterate over the dimensions of the operator
! Note: operator matrices are asusmed symmetric, can only do upper half
      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix
         amew_op_i => amew_op(i)%matrix

         ! The GS-GS contribution
         CALL cp_fm_set_element(amew_op_i, 1, 1, gsgs_op(i))

         ! Compute the operator in the donor MOs basis
         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, gs_coeffs, vec_op, ncol=ndo_mo)
         CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_op, 0.0_dp, prod_fm)
         CALL cp_fm_get_submatrix(prod_fm, domo_op)

         ! Compute the ground-state/singlet components. ao_op*gs_coeffs already stored in vec_op
         CALL cp_gemm('T', 'N', ndo_mo*nsg, ndo_mo, nao, 1.0_dp, sg_coeffs, vec_op, 0.0_dp, sggs_fm)
         DO isg = 1, nsg
            CALL cp_fm_get_submatrix(fm=sggs_fm, target_m=sggs_block, start_row=(isg - 1)*ndo_mo + 1, &
                                     start_col=1, n_rows=ndo_mo, n_cols=ndo_mo)
            diag(:) = get_diag(sggs_block)
            op = sqrt2*SUM(diag)
            CALL cp_fm_set_element(amew_op_i, 1, 1 + isg, op)
         END DO

         ! do the singlet-singlet components
         !start with the overlap
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_sg, 0.0_dp, &
                             dbcsr_work, filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sg, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

         !then the operator in the LR orbital basis
         CALL dbcsr_multiply('N', 'N', 1.0_dp, ao_op_i, dbcsr_sg, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sg, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

         !use the soc routine, it is compatible
         CALL rcs_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_op, pref_trace=-1.0_dp, &
                                    pref_overall=1.0_dp, pref_diags=gsgs_op(i), symmetric=.TRUE.)

         CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
         CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=amew_op_i, nrow=nex, ncol=nex, &
                                 s_firstrow=1, s_firstcol=1, t_firstrow=2, t_firstcol=2)

         ! compute the triplet-triplet components
         !the overlap
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s(1)%matrix, dbcsr_tp, 0.0_dp, &
                             dbcsr_work, filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_tp, dbcsr_work, 0.0_dp, dbcsr_ovlp, filter_eps=eps_filter)

         !the operator in the LR orbital basis
         CALL dbcsr_multiply('N', 'N', 1.0_dp, ao_op_i, dbcsr_sg, 0.0_dp, dbcsr_work, filter_eps=eps_filter)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, dbcsr_sg, dbcsr_work, 0.0_dp, dbcsr_prod, filter_eps=eps_filter)

         CALL rcs_amew_soc_elements(dbcsr_tmp, dbcsr_prod, dbcsr_ovlp, domo_op, pref_trace=-1.0_dp, &
                                    pref_overall=1.0_dp, pref_diags=gsgs_op(i), symmetric=.TRUE.)

         CALL copy_dbcsr_to_fm(dbcsr_tmp, tmp_fm)
         !<T^-1|op|T^-1>
         CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=amew_op_i, nrow=nex, ncol=nex, &
                                 s_firstrow=1, s_firstcol=1, t_firstrow=1 + nsg + 1, t_firstcol=1 + nsg + 1)
         !<T^0|op|T^0>
         CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=amew_op_i, nrow=nex, ncol=nex, &
                                 s_firstrow=1, s_firstcol=1, t_firstrow=1 + nsg + ntp + 1, &
                                 t_firstcol=1 + nsg + ntp + 1)
         !<T^-1|op|T^-1>
         CALL cp_fm_to_fm_submat(msource=tmp_fm, mtarget=amew_op_i, nrow=nex, ncol=nex, &
                                 s_firstrow=1, s_firstcol=1, t_firstrow=1 + nsg + 2*ntp + 1, &
                                 t_firstcol=1 + nsg + 2*ntp + 1)

         ! Symmetrize the matrix (only upper triangle built)
         CALL cp_fm_upper_to_full(amew_op_i, work_fm)

      END DO !i

!  Clean-up
      CALL cp_fm_release(prod_fm)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(tmp_fm)
      CALL cp_fm_release(vec_op)
      CALL cp_fm_release(sggs_fm)
      CALL cp_fm_struct_release(prod_struct)
      CALL cp_fm_struct_release(full_struct)
      CALL cp_fm_struct_release(tmp_struct)
      CALL cp_fm_struct_release(sggs_struct)

   END SUBROUTINE get_rcs_amew_op

! **************************************************************************************************
!> \brief Computes the os SOC matrix elements between excited states AMEWs based on the LR orbitals
!> \param amew_soc output dbcsr matrix with the SOC in the AMEW basis (needs to be fully resereved)
!> \param lr_soc dbcsr matrix with the SOC wrt the LR orbitals
!> \param lr_overlap dbcsr matrix with the excited states LR orbital overlap
!> \param domo_soc the SOC in the basis of the donor MOs
!> \param pref_diaga ...
!> \param pref_diagb ...
!> \param pref_tracea ...
!> \param pref_traceb ...
!> \param pref_diags see notes
!> \param symmetric if the outcome is known to be symmetric, only elements with iex <= jex are done
!> \param tracea_start the indices where to start in the trace part for alpha
!> \param traceb_start the indices where to start in the trace part for beta
!> \note For an excited states pair i,j, the AMEW SOC matrix element is:
!>       soc_ij =   pref_diaga*SUM(alpha part of diag of lr_soc_ij)
!>                + pref_diagb*SUM(beta part of diag of lr_soc_ij)
!>                + pref_tracea*SUM(alpha part of lr_ovlp_ij*TRANSPOSE(domo_soc))
!>                + pref_traceb*SUM(beta part of lr_ovlp_ij*TRANSPOSE(domo_soc))
!>       optinally, one can add pref_diags*SUM(diag lr_ovlp_ij)
! **************************************************************************************************
   SUBROUTINE os_amew_soc_elements(amew_soc, lr_soc, lr_overlap, domo_soc, pref_diaga, &
                                   pref_diagb, pref_tracea, pref_traceb, pref_diags, &
                                   symmetric, tracea_start, traceb_start)

      TYPE(dbcsr_type)                                   :: amew_soc, lr_soc, lr_overlap
      REAL(dp), DIMENSION(:, :)                          :: domo_soc
      REAL(dp)                                           :: pref_diaga, pref_diagb, pref_tracea, &
                                                            pref_traceb
      REAL(dp), OPTIONAL                                 :: pref_diags
      LOGICAL, OPTIONAL                                  :: symmetric
      INTEGER, DIMENSION(2), OPTIONAL                    :: tracea_start, traceb_start

      CHARACTER(len=*), PARAMETER :: routineN = 'os_amew_soc_elements', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, iex, jex, ndo_mo, ndo_so
      INTEGER, DIMENSION(2)                              :: tas, tbs
      LOGICAL                                            :: do_diags, found, my_symm
      REAL(dp)                                           :: soc_elem
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock
      TYPE(dbcsr_iterator_type)                          :: iter

      ndo_so = SIZE(domo_soc, 1)
      ndo_mo = ndo_so/2
      ALLOCATE (diag(ndo_so))
      my_symm = .FALSE.
      IF (PRESENT(symmetric)) my_symm = symmetric
      do_diags = .FALSE.
      IF (PRESENT(pref_diags)) do_diags = .TRUE.

      !by default, alpha part is (1:ndo_mo,1:ndo_mo) and beta is (ndo_mo+1:ndo_so,ndo_mo+1:ndo_so)
      !note: in some SF cases, that might change, mainly because the spin-flip LR-coeffs have
      !inverse order, that is: the beta-coeffs in the alpha spot and the alpha coeffs in the
      !beta spot
      tas = 1
      tbs = ndo_mo + 1
      IF (PRESENT(tracea_start)) tas = tracea_start
      IF (PRESENT(traceb_start)) tbs = traceb_start

      CALL dbcsr_set(amew_soc, 0.0_dp)
      !loop over the excited states pairs as the block of amew_soc (which are all reserved)
      CALL dbcsr_iterator_start(iter, amew_soc)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iex, column=jex, blk=blk)

         IF (my_symm .AND. iex > jex) CYCLE

         !compute the soc matrix element
         soc_elem = 0.0_dp
         CALL dbcsr_get_block_p(lr_soc, iex, jex, pblock, found)
         IF (found) THEN
            diag(:) = get_diag(pblock)
            soc_elem = soc_elem + pref_diaga*SUM(diag(1:ndo_mo)) + pref_diagb*(SUM(diag(ndo_mo + 1:ndo_so)))
         END IF

         CALL dbcsr_get_block_p(lr_overlap, iex, jex, pblock, found)
         IF (found) THEN
            soc_elem = soc_elem &
                       + pref_tracea*SUM(pblock(tas(1):tas(1) + ndo_mo - 1, tas(2):tas(2) + ndo_mo - 1)* &
                                         TRANSPOSE(domo_soc(tas(1):tas(1) + ndo_mo - 1, tas(2):tas(2) + ndo_mo - 1))) &
                       + pref_traceb*SUM(pblock(tbs(1):tbs(1) + ndo_mo - 1, tbs(2):tbs(2) + ndo_mo - 1)* &
                                         TRANSPOSE(domo_soc(tbs(1):tbs(1) + ndo_mo - 1, tbs(2):tbs(2) + ndo_mo - 1)))

            IF (do_diags) THEN
               diag(:) = get_diag(pblock)
               soc_elem = soc_elem + pref_diags*SUM(diag)
            END IF
         END IF

         CALL dbcsr_get_block_p(amew_soc, iex, jex, pblock, found)
         pblock = soc_elem

      END DO
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE os_amew_soc_elements

! **************************************************************************************************
!> \brief Computes the rcs SOC matrix elements between excited states AMEWs based on the LR orbitals
!> \param amew_soc output dbcsr matrix with the SOC in the AMEW basis (needs to be fully resereved)
!> \param lr_soc dbcsr matrix with the SOC wrt the LR orbitals
!> \param lr_overlap dbcsr matrix with the excited states LR orbital overlap
!> \param domo_soc the SOC in the basis of the donor MOs
!> \param pref_trace see notes
!> \param pref_overall see notes
!> \param pref_diags see notes
!> \param symmetric if the outcome is known to be symmetric, only elements with iex <= jex are done
!> \note For an excited states pair i,j, the AMEW SOC matrix element is:
!>       soc_ij = pref_overall*(SUM(diag(lr_soc_ij)) + pref_trace*SUM(lr_overlap_ij*TRANSPOSE(domo_soc)))
!>       optionally, the value pref_diags*SUM(diag(lr_overlap_ij)) can be added (before pref_overall)
! **************************************************************************************************
   SUBROUTINE rcs_amew_soc_elements(amew_soc, lr_soc, lr_overlap, domo_soc, pref_trace, &
                                    pref_overall, pref_diags, symmetric)

      TYPE(dbcsr_type)                                   :: amew_soc, lr_soc, lr_overlap
      REAL(dp), DIMENSION(:, :)                          :: domo_soc
      REAL(dp)                                           :: pref_trace, pref_overall
      REAL(dp), OPTIONAL                                 :: pref_diags
      LOGICAL, OPTIONAL                                  :: symmetric

      CHARACTER(len=*), PARAMETER :: routineN = 'rcs_amew_soc_elements', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, iex, jex
      LOGICAL                                            :: do_diags, found, my_symm
      REAL(dp)                                           :: soc_elem
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock
      TYPE(dbcsr_iterator_type)                          :: iter

      ALLOCATE (diag(SIZE(domo_soc, 1)))
      my_symm = .FALSE.
      IF (PRESENT(symmetric)) my_symm = symmetric
      do_diags = .FALSE.
      IF (PRESENT(pref_diags)) do_diags = .TRUE.

      CALL dbcsr_set(amew_soc, 0.0_dp)
      !loop over the excited states pairs as the block of amew_soc (which are all reserved)
      CALL dbcsr_iterator_start(iter, amew_soc)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iex, column=jex, blk=blk)

         IF (my_symm .AND. iex > jex) CYCLE

         !compute the soc matrix element
         soc_elem = 0.0_dp
         CALL dbcsr_get_block_p(lr_soc, iex, jex, pblock, found)
         IF (found) THEN
            diag(:) = get_diag(pblock)
            soc_elem = soc_elem + SUM(diag)
         END IF

         CALL dbcsr_get_block_p(lr_overlap, iex, jex, pblock, found)
         IF (found) THEN
            soc_elem = soc_elem + pref_trace*SUM(pblock*TRANSPOSE(domo_soc))

            IF (do_diags) THEN
               diag(:) = get_diag(pblock)
               soc_elem = soc_elem + pref_diags*SUM(diag)
            END IF
         END IF

         CALL dbcsr_get_block_p(amew_soc, iex, jex, pblock, found)
         pblock = pref_overall*soc_elem

      END DO
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE rcs_amew_soc_elements

! **************************************************************************************************
!> \brief Computes the dipole oscillator strengths in the AMEWs basis for SOC
!> \param soc_evecs_cfm the complex AMEWs coefficients
!> \param dbcsr_soc_package ...
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \param gs_coeffs the ground state coefficients, given for open-shell because in ROKS, the gs_coeffs
!>                  are stored slightly differently within SOC for efficiency and code uniquness
! **************************************************************************************************
   SUBROUTINE compute_soc_dipole_fosc(soc_evecs_cfm, dbcsr_soc_package, donor_state, xas_tdp_env, &
                                      xas_tdp_control, qs_env, gs_coeffs)

      TYPE(cp_cfm_type), POINTER                         :: soc_evecs_cfm
      TYPE(dbcsr_soc_package_type)                       :: dbcsr_soc_package
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: gs_coeffs

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_soc_dipole_fosc', &
         routineP = moduleN//':'//routineN

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: transdip
      INTEGER                                            :: handle, i, nosc, ntot
      LOGICAL                                            :: do_os, do_rcs
      REAL(dp), DIMENSION(:), POINTER                    :: osc_str, soc_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type), POINTER                         :: dip_cfm, work1_cfm, work2_cfm
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: amew_dip
      TYPE(cp_fm_struct_type), POINTER                   :: dip_struct, full_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env

      NULLIFY (para_env, blacs_env, dip_struct, full_struct, work1_cfm, work2_cfm, osc_str)
      NULLIFY (soc_evals, amew_dip, dip_cfm)

      CALL timeset(routineN, handle)

      !init
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      do_os = xas_tdp_control%do_spin_cons
      do_rcs = xas_tdp_control%do_singlet
      soc_evals => donor_state%soc_evals
      nosc = SIZE(soc_evals)
      ntot = nosc + 1 !because GS AMEW is in there
      ALLOCATE (donor_state%soc_osc_str(nosc))
      osc_str => donor_state%soc_osc_str
      osc_str(:) = 0.0_dp
      IF (do_os .AND. .NOT. PRESENT(gs_coeffs)) CPABORT("Need to pass gs_coeffs for open-shell")

      !get some work arrays/matrix
      CALL cp_fm_struct_create(dip_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=1)
      CALL cp_cfm_get_info(soc_evecs_cfm, matrix_struct=full_struct)
      CALL cp_cfm_create(dip_cfm, dip_struct)
      CALL cp_cfm_create(work1_cfm, full_struct)
      CALL cp_cfm_create(work2_cfm, full_struct)
      ALLOCATE (transdip(ntot, 1))

      !get the dipole in the AMEW basis
      IF (do_os) THEN
         CALL get_os_amew_op(amew_dip, xas_tdp_env%dipmat, gs_coeffs, dbcsr_soc_package, &
                             donor_state, xas_tdp_control%eps_filter, qs_env)
      ELSE
         CALL get_rcs_amew_op(amew_dip, xas_tdp_env%dipmat, dbcsr_soc_package, donor_state, &
                              xas_tdp_control%eps_filter, qs_env)
      END IF

      DO i = 1, 3 !cartesian coord x, y, z

         !Convert the real dipole into the cfm format for calculations
         CALL cp_fm_to_cfm(msourcer=amew_dip(i)%matrix, mtarget=work1_cfm)

         !compute amew_coeffs^dagger * amew_dip * amew_gs to get the transition moments
         CALL cp_cfm_gemm('C', 'N', ntot, ntot, ntot, (1.0_dp, 0.0_dp), soc_evecs_cfm, work1_cfm, &
                          (0.0_dp, 0.0_dp), work2_cfm)
         CALL cp_cfm_gemm('N', 'N', ntot, 1, ntot, (1.0_dp, 0.0_dp), work2_cfm, soc_evecs_cfm, &
                          (0.0_dp, 0.0_dp), dip_cfm)

         CALL cp_cfm_get_submatrix(dip_cfm, transdip)

         !transition dipoles are real numbers
         osc_str(:) = osc_str(:) + REAL(transdip(2:ntot, 1))**2 + AIMAG(transdip(2:ntot, 1))**2

      END DO !i

      !multiply with appropriate prefac depending in the rep
      IF (xas_tdp_control%dipole_form == xas_dip_len) THEN
         osc_str(:) = 2.0_dp/3.0_dp*soc_evals(:)*osc_str(:)
      ELSE
         osc_str(:) = 2.0_dp/3.0_dp/soc_evals(:)*osc_str(:)
      END IF

      !clean-up
      CALL cp_fm_struct_release(dip_struct)
      CALL cp_cfm_release(work1_cfm)
      CALL cp_cfm_release(work2_cfm)
      CALL cp_cfm_release(dip_cfm)
      DO i = 1, 3
         CALL cp_fm_release(amew_dip(i)%matrix)
      END DO
      DEALLOCATE (amew_dip, transdip)

      CALL timestop(handle)

   END SUBROUTINE compute_soc_dipole_fosc

! **************************************************************************************************
!> \brief Computes the quadrupole oscillator strengths in the AMEWs basis for SOC
!> \param soc_evecs_cfm the complex AMEWs coefficients
!> \param dbcsr_soc_package inherited from the main SOC routine
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \param gs_coeffs the ground state coefficients, given for open-shell because in ROKS, the gs_coeffs
!>                  are stored slightly differently within SOC for efficiency and code uniquness
! **************************************************************************************************
   SUBROUTINE compute_soc_quadrupole_fosc(soc_evecs_cfm, dbcsr_soc_package, donor_state, &
                                          xas_tdp_env, xas_tdp_control, qs_env, gs_coeffs)

      TYPE(cp_cfm_type), POINTER                         :: soc_evecs_cfm
      TYPE(dbcsr_soc_package_type)                       :: dbcsr_soc_package
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: gs_coeffs

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_soc_quadrupole_fosc', &
         routineP = moduleN//':'//routineN

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: trace
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: transquad
      INTEGER                                            :: handle, i, nosc, ntot
      LOGICAL                                            :: do_os, do_rcs
      REAL(dp), DIMENSION(:), POINTER                    :: osc_str, soc_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type), POINTER                         :: quad_cfm, work1_cfm, work2_cfm
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: amew_quad
      TYPE(cp_fm_struct_type), POINTER                   :: full_struct, quad_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env

      NULLIFY (para_env, blacs_env, quad_struct, full_struct, work1_cfm, work2_cfm, osc_str)
      NULLIFY (soc_evals, amew_quad, quad_cfm)

      CALL timeset(routineN, handle)

      !init
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      do_os = xas_tdp_control%do_spin_cons
      do_rcs = xas_tdp_control%do_singlet
      soc_evals => donor_state%soc_evals
      nosc = SIZE(soc_evals)
      ntot = nosc + 1 !because GS AMEW is in there
      ALLOCATE (donor_state%soc_quad_osc_str(nosc))
      osc_str => donor_state%soc_quad_osc_str
      osc_str(:) = 0.0_dp
      IF (do_os .AND. .NOT. PRESENT(gs_coeffs)) CPABORT("Need to pass gs_coeffs for open-shell")

      !get some work arrays/matrix
      CALL cp_fm_struct_create(quad_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=1)
      CALL cp_cfm_get_info(soc_evecs_cfm, matrix_struct=full_struct)
      CALL cp_cfm_create(quad_cfm, quad_struct)
      CALL cp_cfm_create(work1_cfm, full_struct)
      CALL cp_cfm_create(work2_cfm, full_struct)
      ALLOCATE (transquad(ntot, 1))
      ALLOCATE (trace(nosc))
      trace = (0.0_dp, 0.0_dp)

      !get the quadrupole in the AMEWs basis
      IF (do_os) THEN
         CALL get_os_amew_op(amew_quad, xas_tdp_env%quadmat, gs_coeffs, dbcsr_soc_package, &
                             donor_state, xas_tdp_control%eps_filter, qs_env)
      ELSE
         CALL get_rcs_amew_op(amew_quad, xas_tdp_env%quadmat, dbcsr_soc_package, donor_state, &
                              xas_tdp_control%eps_filter, qs_env)
      END IF

      DO i = 1, 6 ! x2, xy, xz, y2, yz, z2

         !Convert the real quadrupole into a cfm for further calculation
         CALL cp_fm_to_cfm(msourcer=amew_quad(i)%matrix, mtarget=work1_cfm)

         !compute amew_coeffs^dagger * amew_quad * amew_gs to get the transition moments
         CALL cp_cfm_gemm('C', 'N', ntot, ntot, ntot, (1.0_dp, 0.0_dp), soc_evecs_cfm, work1_cfm, &
                          (0.0_dp, 0.0_dp), work2_cfm)
         CALL cp_cfm_gemm('N', 'N', ntot, 1, ntot, (1.0_dp, 0.0_dp), work2_cfm, soc_evecs_cfm, &
                          (0.0_dp, 0.0_dp), quad_cfm)

         CALL cp_cfm_get_submatrix(quad_cfm, transquad)

         !if x2, y2 or z2, need to keep track of trace
         IF (i == 1 .OR. i == 4 .OR. i == 6) THEN
            osc_str(:) = osc_str(:) + REAL(transquad(2:ntot, 1))**2 + AIMAG(transquad(2:ntot, 1))**2
            trace(:) = trace(:) + transquad(2:ntot, 1)

            !if xy, xz, or yz, need to count twice (for yx, zx and zy)
         ELSE
            osc_str(:) = osc_str(:) + 2.0_dp*(REAL(transquad(2:ntot, 1))**2 + AIMAG(transquad(2:ntot, 1))**2)
         END IF

      END DO !i

      !remove a third of the trace
      osc_str(:) = osc_str(:) - 1._dp/3._dp*(REAL(trace(:))**2 + AIMAG(trace(:))**2)

      !multiply by the prefactor
      osc_str(:) = osc_str(:)*1._dp/20._dp*a_fine**2*soc_evals(:)**3

      !clean-up
      CALL cp_fm_struct_release(quad_struct)
      CALL cp_cfm_release(work1_cfm)
      CALL cp_cfm_release(work2_cfm)
      CALL cp_cfm_release(quad_cfm)
      DO i = 1, 6
         CALL cp_fm_release(amew_quad(i)%matrix)
      END DO
      DEALLOCATE (amew_quad, transquad, trace)

      CALL timestop(handle)

   END SUBROUTINE compute_soc_quadrupole_fosc

END MODULE xas_tdp_utils

