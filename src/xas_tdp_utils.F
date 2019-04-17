!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for X-ray absorption spectroscopy using TDDFPT
!> \author AB (01.2018)
! **************************************************************************************************

MODULE xas_tdp_utils
   USE ai_contraction_sphi,             ONLY: ab_contract
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_fm_to_cfm, cp_cfm_to_fm, cp_cfm_type, &
                                              cp_cfm_create, cp_cfm_release, cp_cfm_get_submatrix, &
                                              cp_cfm_get_element
   USE cp_eri_mme_interface,            ONLY: cp_eri_mme_param, cp_eri_mme_set_params
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale, cp_fm_trace, cp_fm_upper_to_full
   USE cp_fm_diag,                      ONLY: choose_eigv_solver, cp_fm_geeig, cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create, cp_fm_struct_release, &     
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, &
        cp_fm_release, cp_fm_set_element, cp_fm_to_fm_submat, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_complete_redistribute, dbcsr_copy, dbcsr_create, &
        dbcsr_desymmetrize, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_block_p, &
        dbcsr_get_diag, dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_set, dbcsr_set_diag, dbcsr_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE eri_mme_integrate,               ONLY: eri_mme_2c_integrate
   USE eri_mme_types,                   ONLY: eri_mme_init, eri_mme_release
   USE generic_os_integrals,            ONLY: int_operators_r12_ab_os
   USE input_constants,                 ONLY: tddfpt_singlet, tddfpt_triplet, xas_tdp_uks, &
                                              xas_tdp_roks, tddfpt_spin_cons, tddfpt_spin_flip, &
                                              do_hfx_potential_coulomb, do_hfx_potential_truncated,&
                                              xas_dip_len, xas_dip_vel, do_eri_mme
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat_symm, get_diag
   USE message_passing,                 ONLY: mp_min, mp_max, mp_sync
   USE molecule_types,                  ONLY: molecule_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: qs_environment_type, get_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind, qs_kind_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE qs_o3c_methods,                  ONLY: calculate_o3c_coulomb_integrals
   USE qs_o3c_types,                    ONLY: get_o3c_iterator_info,&
                                              init_o3c_container,&
                                              o3c_container_type,&
                                              o3c_iterate,&
                                              o3c_iterator_create,&
                                              o3c_iterator_release,&
                                              o3c_iterator_type
   USE qs_ot_eigensolver,               ONLY: ot_eigensolver
   USE util,                            ONLY: locate
   USE xas_tdp_kernel,                  ONLY: kernel_coulomb_xc, kernel_exchange
   USE xas_tdp_types,                   ONLY: donor_state_type, xas_tdp_env_type, &
                                              xas_tdp_control_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_utils'

   PUBLIC :: setup_xas_tdp_prob, solve_xas_tdp_prob, build_xas_tdp_3c_nl, &
             compute_ri_coulomb2_int, include_rcs_soc, include_os_soc, &
             compute_ri_o3c_coulomb, compute_ri_o3c_exchange, build_xas_tdp_ovlp_nl

CONTAINS

! **************************************************************************************************
!> \brief Builds the matrix that defines the XAS TDDFPT generalized eigenvalue problem to be solved
!>        for excitation energies omega. The problem has the form omega*G*C = M*C, where C contains
!>        the reponse orbitals coefficients. The matrix M and the metric G are stored in the given
!>        donor_state
!> \param donor_state the donor_state for which the problem is restricted
!> \param qs_env ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \note  the matrix M is symmetric and has the form | M_d   M_o |
!>                                                   | M_o   M_d |,
!>       -In the SPIN-RESTRICTED case:
!>        depending on whther we consider singlet or triplet excitation, the diagonal (M_d) and
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
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: handle
      REAL(dp)                                        :: sx
      LOGICAL                                         :: do_xc, do_hfx, do_sg, do_tp, do_tda, &
                                                         do_os, do_sc, do_sf
      TYPE(dbcsr_type), POINTER                       :: sg_matrix_tdp, tp_matrix_tdp, &
                                                         matrix_d, matrix_e_sc, matrix_c_sf,&
                                                         matrix_c_sg, matrix_c_tp, & 
                                                         sc_matrix_tdp, sf_matrix_tdp, matrix_c_sc
      TYPE(dbcsr_type)                                :: matrix_b, proj_Q, proj_Q_sf, matrix_a, &
                                                         matrix_a_sf, work
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: ex_ker, xc_ker
      TYPE(dbcsr_distribution_type), POINTER          :: submat_dist
      INTEGER, DIMENSION(:), POINTER                  :: submat_blk_size

      NULLIFY(sg_matrix_tdp, tp_matrix_tdp, submat_dist, submat_blk_size, matrix_c_sf)
      NULLIFY(matrix_c_sg, matrix_c_tp, matrix_c_sc, matrix_d, matrix_e_sc)
      NULLIFY(sc_matrix_tdp, sf_matrix_tdp, ex_ker, xc_ker)

      CALL timeset(routineN, handle)

!  Initialization
      do_os = xas_tdp_control%do_uks .OR. xas_tdp_control%do_roks
      do_sc = xas_tdp_control%do_spin_cons
      do_sf = xas_tdp_control%do_spin_flip
      do_sg = xas_tdp_control%do_singlet
      do_tp = xas_tdp_control%do_triplet
      do_xc = xas_tdp_control%do_xc
      do_hfx = xas_tdp_control%do_hfx
      do_tda = xas_tdp_control%tamm_dancoff
      sx = xas_tdp_control%sx
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
      CALL compute_submat_dist_and_blk_size(submat_dist, submat_blk_size, donor_state, do_os, qs_env)

!  Allocate and compute all the matrices A, B, C, etc we will need

      ! The projector(s) on the unoccupied unperturbed ground state 1-SP and associated work matrix
      IF (do_sg .OR. do_tp .OR. do_sc) THEN !spin-conserving
         CALL get_q_projector(proj_Q, submat_dist, submat_blk_size, donor_state, do_os, xas_tdp_env)
      END IF
      IF (do_sf) THEN !spin-flip
         CALL get_q_projector(proj_Q_sf, submat_dist, submat_blk_size, donor_state, do_os, &
                              xas_tdp_env, do_sf=.TRUE.)
      END IF
      CALL dbcsr_create(matrix=work, matrix_type="N", dist=submat_dist, name="WORK", &
                        row_blk_size=submat_blk_size, col_blk_size=submat_blk_size)

      ! The ground state contribution(s)
      IF (do_sg .OR. do_tp .OR. do_sc) THEN !spin-conserving
         CALL build_gs_contribution(matrix_a, submat_dist, submat_blk_size, donor_state, do_os, &
                                    qs_env)
      END IF
      IF (do_sf) THEN !spin-flip
         CALL build_gs_contribution(matrix_a_sf, submat_dist, submat_blk_size, donor_state, do_os, &
                                    qs_env, do_sf=.TRUE.)
      END IF

      ! The Coulomb and XC kernels. Internal analysis to know which matrix to compute
      CALL dbcsr_allocate_matrix_set(xc_ker, 4)
      ALLOCATE (xc_ker(1)%matrix, xc_ker(2)%matrix, xc_ker(3)%matrix, xc_ker(4)%matrix)
      CALL kernel_coulomb_xc(matrix_b, xc_ker, submat_dist, submat_blk_size, donor_state, &
                             xas_tdp_env, xas_tdp_control, qs_env)
      matrix_c_sg => xc_ker(1)%matrix; matrix_c_tp => xc_ker(2)%matrix
      matrix_c_sc => xc_ker(3)%matrix; matrix_c_sf => xc_ker(4)%matrix

      ! The exact exchange. Internal analysis to know which matrices to compute
      CALL dbcsr_allocate_matrix_set(ex_ker, 2)
      ALLOCATE(ex_ker(1)%matrix, ex_ker(2)%matrix)
      CALL kernel_exchange(ex_ker, submat_dist, submat_blk_size, donor_state, xas_tdp_env, &
                           xas_tdp_control, qs_env)
      matrix_d => ex_ker(1)%matrix; matrix_e_sc => ex_ker(2)%matrix

      ! Build the metric G, also need its inverse in case of full-TDDFT
      IF (do_tda) THEN
         ALLOCATE (donor_state%metric(1))
         CALL build_metric(donor_state%metric, submat_dist, submat_blk_size, donor_state, &
                           qs_env, do_os)
      ELSE
         ALLOCATE (donor_state%metric(2))
         CALL build_metric(donor_state%metric, submat_dist, submat_blk_size, donor_state, &
                           qs_env, do_os, do_inv=.TRUE.)
      END IF

!  Build the eigenvalue problem, depending on the case (TDA, singlet, triplet, hfx, etc ...)
      IF (do_tda) THEN

         IF (do_sc) THEN ! open-shell spin-conserving under TDA

            ! The final matrix is M = A + B + C - D
            CALL dbcsr_copy(sc_matrix_tdp, matrix_a, name="OS MATRIX TDP")
            CALL dbcsr_add(sc_matrix_tdp, matrix_b, 1.0_dp, 1.0_dp)

            IF (do_xc) CALL dbcsr_add(sc_matrix_tdp, matrix_c_sc, 1.0_dp, 1.0_dp) !xc kernel
            IF (do_hfx) CALL dbcsr_add(sc_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx) !scaled hfx

            ! The product with the Q projector
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, sc_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, sc_matrix_tdp)

         END IF !do_sc

         IF (do_sf) THEN ! open-shell spin-flip under TDA

            ! The final matrix is M = A + C - D 
            CALL dbcsr_copy(sf_matrix_tdp, matrix_a_sf, name="OS MATRIX TDP")

            IF (do_xc) CALL dbcsr_add(sf_matrix_tdp, matrix_c_sf, 1.0_dp, 1.0_dp) !xc kernel
            IF (do_hfx) CALL dbcsr_add(sf_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx) !scaled hfx

            ! Take the product with the (spin-flip) Q projector
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q_sf, sf_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q_sf, 0.0_dp, sf_matrix_tdp)

         END IF !do_sf

         IF (do_sg) THEN ! singlets under TDA

            ! The final matrix is M = A + 2B + (C_aa + C_ab) - D
            CALL dbcsr_copy(sg_matrix_tdp, matrix_a, name="SINGLET MATRIX TDP")
            CALL dbcsr_add(sg_matrix_tdp, matrix_b, 1.0_dp, 2.0_dp)

            IF (do_xc) CALL dbcsr_add(sg_matrix_tdp, matrix_c_sg, 1.0_dp, 1.0_dp) ! xc kernel
            IF (do_hfx) CALL dbcsr_add(sg_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx) ! scaled hfx

            ! Take the product with the Q projector:
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, sg_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, sg_matrix_tdp)

         END IF !do_sg (TDA)

         IF (do_tp) THEN ! triplets under TDA

            ! The final matrix is M =  A + (C_aa - C_ab) - D
            CALL dbcsr_copy(tp_matrix_tdp, matrix_a, name="TRIPLET MATRIX TDP")

            IF (do_xc) CALL dbcsr_add(tp_matrix_tdp, matrix_c_tp, 1.0_dp, 1.0_dp) ! xc_kernel
            IF (do_hfx) CALL dbcsr_add(tp_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx) ! scaled hfx

            ! Take the product with the Q projector:
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, tp_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, tp_matrix_tdp)

         END IF !do_tp (TDA)

      ELSE ! not TDA

         ! In the case of full-TDDFT, the problem is turned Hermitian with the help of auxiliary
         ! matrices AUX = (A-D+E)^(+-0.5) that are stored in donor_state
         CALL build_aux_matrix(1.0E-8_dp, sx, matrix_a, matrix_d, matrix_e_sc, do_hfx, proj_Q, &
                               work, donor_state, qs_env)

         IF (do_sc) THEN !full-TDDFT open-shell spin-conserving

            ! The final matrix is the sum of the on- and off-diagonal elements as in the description
            ! M = A + 2B + 2C - D - E
            CALL dbcsr_copy(sc_matrix_tdp, matrix_a, name="OS MATRIX TDP")
            CALL dbcsr_add(sc_matrix_tdp, matrix_b, 1.0_dp, 2.0_dp)

            IF (do_hfx) THEN !scaled hfx
               CALL dbcsr_add(sc_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx)
               CALL dbcsr_add(sc_matrix_tdp, matrix_e_sc, 1.0_dp, -1.0_dp*sx)
            END IF
            IF (do_xc) THEN
               CALL dbcsr_add(sc_matrix_tdp, matrix_c_sc, 1.0_dp, 2.0_dp)
            END IF

            ! Take the product with the Q projector
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, sc_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, sc_matrix_tdp)

            ! Take the product with the inverse metric
            ! M <= G^-1 * M * G^-1
            CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, sc_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, work, donor_state%metric(2)%matrix, 0.0_dp, sc_matrix_tdp)

         END IF

         IF (do_sg) THEN ! full-TDDFT singlets

            ! The final matrix is the sum of the on- and off-diagonal elements as in the description
            ! M = A + 4B + 2(C_aa + C_ab) - D - E
            CALL dbcsr_copy(sg_matrix_tdp, matrix_a, name="SINGLET MATRIX TDP")
            CALL dbcsr_add(sg_matrix_tdp, matrix_b, 1.0_dp, 4.0_dp)

            IF (do_hfx) THEN !scaled hfx
               CALL dbcsr_add(sg_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx)
               CALL dbcsr_add(sg_matrix_tdp, matrix_e_sc, 1.0_dp, -1.0_dp*sx)
            END IF
            IF (do_xc) THEN !xc kernel
               CALL dbcsr_add(sg_matrix_tdp, matrix_c_sg, 1.0_dp, 2.0_dp)
            END IF

            ! Take the product with the Q projector
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, sg_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, sg_matrix_tdp)

            ! Take the product with the inverse metric
            ! M <= G^-1 * M * G^-1
            CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, sg_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, work, donor_state%metric(2)%matrix, 0.0_dp, sg_matrix_tdp)

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
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, tp_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, tp_matrix_tdp)

            ! Take the product with the inverse metric
            ! M <= G^-1 * M * G^-1
            CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, tp_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, work, donor_state%metric(2)%matrix, 0.0_dp, tp_matrix_tdp)

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
      CALL dbcsr_distribution_release(submat_dist)
      DEALLOCATE (submat_dist, submat_blk_size)

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

      INTEGER                                            :: handle, i, imo, ispin, nao, ndo_mo, &
                                                            nelectron, nevals, nocc, nrow, nspins
      LOGICAL                                            :: do_os, do_sf, full_diag
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: scaling, tmp_evals
      REAL(dp), DIMENSION(:), POINTER                    :: lr_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: lr_coeffs
      TYPE(cp_fm_struct_type), POINTER                   :: ex_struct, fm_struct, ot_fm_struct
      TYPE(cp_fm_type), POINTER                          :: c_diff, c_sum, lhs_matrix, ot_coeffs, &
                                                            rhs_matrix, work
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type)                                   :: tmp_mat
      TYPE(dbcsr_type), POINTER                          :: matrix_tdp

      CALL timeset(routineN, handle)

      NULLIFY (para_env, blacs_env, fm_struct, rhs_matrix, matrix_tdp, lhs_matrix, work)
      NULLIFY (c_diff, c_sum, ex_struct, lr_evals, lr_coeffs, ot_fm_struct, ot_coeffs)
      CPASSERT(ASSOCIATED(xas_tdp_env))

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
      full_diag = xas_tdp_control%full_diag
      CALL cp_fm_get_info(donor_state%gs_coeffs, nrow_global=nao)
      CALL dbcsr_get_info(matrix_tdp, nfullrows_total=nrow)
      ndo_mo = donor_state%ndo_mo
      nocc = nelectron/2; IF (do_os) nocc = nelectron
      nevals = nspins*nao-nocc
      IF (xas_tdp_control%n_excited > 0 .AND. xas_tdp_control%n_excited < nevals) THEN
         nevals = xas_tdp_control%n_excited
      END IF
      nevals = ndo_mo*nevals !as in input description, multiply by # of donor MOs
      nocc = ndo_mo*nocc

      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nrow, &
                               para_env=para_env, ncol_global=nrow)
      CALL cp_fm_create(c_sum, fm_struct)

!     Allocating space for solutions and oscilaltor strengths
      IF (ex_type == tddfpt_spin_cons) THEN
         ALLOCATE (donor_state%sc_evals(nevals))
         lr_evals => donor_state%sc_evals
         ALLOCATE (donor_state%sc_coeffs(nevals))
         lr_coeffs => donor_state%sc_coeffs
         ALLOCATE (donor_state%sc_osc_str(nevals))
         donor_state%sc_osc_str = 0.0_dp
      ELSE IF (ex_type == tddfpt_spin_flip) THEN
         ALLOCATE (donor_state%sf_evals(nevals))
         lr_evals => donor_state%sf_evals
         ALLOCATE (donor_state%sf_coeffs(nevals))
         lr_coeffs => donor_state%sf_coeffs
         !note: no oscillator strength because spin-forbidden
      ELSE IF (ex_type == tddfpt_singlet) THEN
         ALLOCATE (donor_state%sg_evals(nevals))
         lr_evals => donor_state%sg_evals
         ALLOCATE (donor_state%sg_coeffs(nevals))
         lr_coeffs => donor_state%sg_coeffs
         ALLOCATE (donor_state%sg_osc_str(nevals))
         donor_state%sg_osc_str = 0.0_dp
      ELSE IF (ex_type == tddfpt_triplet) THEN
         ALLOCATE (donor_state%tp_evals(nevals))
         lr_evals => donor_state%tp_evals
         ALLOCATE (donor_state%tp_coeffs(nevals))
         lr_coeffs => donor_state%tp_coeffs
         !note: no oscillator strength for the triplets (by construction)
      END IF

!     Test on the diagonalization method
      IF (full_diag) THEN

         ! create the fm infrastructure
         ALLOCATE (tmp_evals(nrow))
         CALL cp_fm_create(rhs_matrix, fm_struct)
         CALL cp_fm_create(work, fm_struct)

!        Test on TDA
         IF (xas_tdp_control%tamm_dancoff) THEN

!           Get the main matrix_tdp as an fm
            CALL copy_dbcsr_to_fm(matrix_tdp, rhs_matrix)

!           Get the metric as a fm
            CALL cp_fm_create(lhs_matrix, fm_struct)
            CALL copy_dbcsr_to_fm(donor_state%metric(1)%matrix, lhs_matrix)

!           Diagonalisation (Cholesky decomposition). In TDA, c_sum = c^-
            CALL cp_fm_geeig(rhs_matrix, lhs_matrix, c_sum, tmp_evals, work)

!           Take the first nevals eigenvalues that are not projected to zero
            lr_evals(:) = tmp_evals(nocc+1:nocc+nevals)

!           TDA specific clean-up
            CALL cp_fm_release(lhs_matrix)

         ELSE ! not TDA

!           Need to multiply the current matrix_tdp with the auxiliary matrix
!           tmp_mat =  (A-D+E)^0.5 * M * (A-D+E)^0.5
            CALL dbcsr_create(matrix=tmp_mat, template=matrix_tdp, matrix_type="N")
            CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%matrix_aux, matrix_tdp, 0.0_dp, tmp_mat)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, tmp_mat, donor_state%matrix_aux, 0.0_dp, tmp_mat)

!           Get the matrix as a fm
            CALL copy_dbcsr_to_fm(tmp_mat, rhs_matrix)

!           Solve the "turned-Hermitian" eigenvalue problem
            CALL choose_eigv_solver(rhs_matrix, work, tmp_evals)

!           Currently, work = (A-D+E)^0.5 (c^+ - c^-) and tmp_evals = omega^2
!           Put tiny almost zero eigenvalues to zero (corresponding to occupied MOs)
            WHERE (tmp_evals < 1.0E-4_dp) tmp_evals = 0.0_dp

!           Retrive c_diff = (c^+ - c^-) for normalization
!           (c^+ - c^-) = 1/omega^2 * M * (A-D+E)^0.5 * work
            CALL cp_fm_create(c_diff, fm_struct)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_tdp, donor_state%matrix_aux, 0.0_dp, tmp_mat)
            CALL cp_dbcsr_sm_fm_multiply(tmp_mat, work, c_diff, ncol=nrow)

            ALLOCATE (scaling(nrow))
            WHERE (tmp_evals .NE. 0.0_dp) scaling = 1.0_dp/tmp_evals
            CALL cp_fm_column_scale(c_diff, scaling)

!           Normalize with the metric: c_diff * G * c_diff = +- 1
            scaling = 0.0_dp
            CALL get_normal_scaling(scaling, c_diff, donor_state)
            CALL cp_fm_column_scale(c_diff, scaling)

!           Get the actual eigenvalues
            tmp_evals = SQRT(tmp_evals)
            lr_evals(:) = tmp_evals(nocc+1:nocc+nevals)

!           Get c_sum = (c^+ + c^-), which appears in all transition density related expressions
!           c_sum = -1/omega G^-1 * (A-D+E) * (c^+ - c^-)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%matrix_aux, donor_state%matrix_aux, 0.0_dp, tmp_mat)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, tmp_mat, 0.0_dp, tmp_mat)
            CALL cp_dbcsr_sm_fm_multiply(tmp_mat, c_diff, c_sum, ncol=nrow)
            WHERE (tmp_evals .NE. 0) scaling = -1.0_dp/tmp_evals
            CALL cp_fm_column_scale(c_sum, scaling)

!           Full TDDFT specific clean-up
            CALL cp_fm_release(c_diff)
            CALL dbcsr_release(tmp_mat)
            DEALLOCATE (scaling)

         END IF ! TDA

!        Full matrix clean-up
         CALL cp_fm_release(rhs_matrix)
         CALL cp_fm_release(work)

      ELSE ! iterative solver using the OT method

         !only TDA is available for iterative solver
         CPASSERT(xas_tdp_control%tamm_dancoff)

         !Need fm of the correct size
         CALL cp_fm_struct_create(ot_fm_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=nrow, ncol_global=nevals)
         CALL cp_fm_create(ot_coeffs, ot_fm_struct)

         CALL xas_ot_solver(matrix_tdp, donor_state%metric(1)%matrix, ot_coeffs, lr_evals, nevals, &
                            donor_state, xas_tdp_env, xas_tdp_control, qs_env, do_sf)

         !Put the newly obtained coeffs into c_sum
         CALL cp_fm_to_fm_submat(msource=ot_coeffs, mtarget=c_sum, nrow=nrow, ncol=nevals, &
                                 s_firstrow=1, s_firstcol=1, t_firstrow=1, t_firstcol=nocc+1)

         !ot specific clean-up
         CALL cp_fm_release(ot_coeffs)
         CALL cp_fm_struct_release(ot_fm_struct)

      END IF ! diagonalization type

!  Reorganize the eigenvectors in array of cp_fm_p_types so that each element corresponds to an
!  excited state. Makes later calls to those easier and more efficient
!  In case of open-shell, we store the coeffs in the same logic as the matrix => first block where
!  the columns are the c_Ialpha and second block with columns as c_Ibeta
      CALL cp_fm_struct_create(ex_struct, nrow_global=nao, ncol_global=ndo_mo*nspins, &
                               para_env=para_env, context=blacs_env)

      DO i = 1, nevals

         CALL cp_fm_create(lr_coeffs(i)%matrix, ex_struct)

         DO ispin = 1, nspins
            DO imo = 1, ndo_mo

               CALL cp_fm_to_fm_submat(msource=c_sum, mtarget=lr_coeffs(i)%matrix, &
                                       nrow=nao, ncol=1, s_firstrow=((ispin-1)*ndo_mo+imo-1)*nao+1, &
                                       s_firstcol=nocc+i, t_firstrow=1, t_firstcol=(ispin-1)*ndo_mo+imo)
            END DO !imo
         END DO !ispin
      END DO !istate

!  Clean-up
      CALL cp_fm_release(c_sum)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_struct_release(ex_struct)

!  Perform a partial clean-up of the donor_state
      CALL dbcsr_release(matrix_tdp)

      CALL timestop(handle)

   END SUBROUTINE solve_xas_tdp_prob

! **************************************************************************************************
!> \brief An iterative solver based on OT for the TDA generalized eigV problem lambda Sx = Hx
!> \param matrix_tdp the RHS matrix (dbcsr)
!> \param metric the LHS matrix (dbcsr)
!> \param evecs the corresponding eigenvectors (fm)
!> \param evals the corresponding eigenvalues
!> \param neig the number of wanted eigenvalues
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control
!> \param qs_env ...
!> \param do_sf Whtether the problem involves spin-flip => the guesses are done differently
!> \note Requires LUMOs from the ground state calculations for guesses and will converge much faster
!>       if there is a gap in the spectrum right after the last computed excitation energy
!>       All checks on number of LUMOs and appropriate N_EXCITED are done before hand
!>       Note: extensive experimentations have shown that no preconditioner is actually the best
! **************************************************************************************************
   SUBROUTINE xas_ot_solver(matrix_tdp, metric, evecs, evals, neig, donor_state, xas_tdp_env, &
                            xas_tdp_control, qs_env, do_sf)

      TYPE(dbcsr_type), POINTER                       :: matrix_tdp, metric
      TYPE(cp_fm_type), POINTER                       :: evecs
      REAL(dp), DIMENSION(:), INTENT(INOUT)           :: evals
      INTEGER, INTENT(IN)                             :: neig
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                   :: do_sf

      CHARACTER(len=*), PARAMETER :: routineN = "xas_ot_solver", routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, ndo_mo, iex, imo, nao,&
                                                         ido_mo, max_iter, output_unit, &
                                                         nelec_spin(2), minel, offset, ilumo, &
                                                         ispin, nspins, ex_spin, homo(2), nrows, &
                                                         ncols, t_firstcol, t_firstrow
      REAL(dp)                                        :: eps_iter
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      LOGICAL                                         :: do_os, my_dosf
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: lumo_guess
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      TYPE(cp_fm_type), POINTER                       :: mo_coeff, ortho_space
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct

      NULLIFY(para_env, blacs_env, lumo_guess, mos, fm_struct, ortho_space)

      CALL timeset(routineN, handle)

      output_unit = cp_logger_get_default_io_unit()
      IF (output_unit > 0) THEN
         WRITE (output_unit, "(/,T5,A)") &
            "Using OT eigensolver for diagonalization: "
      END IF

!  Initialization
      do_os = xas_tdp_control%do_uks .OR. xas_tdp_control%do_roks
      nspins = 1; IF(do_os) nspins = 2
      ndo_mo = donor_state%ndo_mo
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env, nelectron_spin=nelec_spin)
      CALL cp_fm_get_info(evecs, nrow_global=nao)
      nao = nao/ndo_mo/nspins
      max_iter = xas_tdp_control%max_iter
      eps_iter = xas_tdp_control%eps_iter
      ALLOCATE (lumo_guess(nspins))
      lumo_guess(1)%matrix => xas_tdp_env%lumo_guess(1)%matrix
      IF (do_os) lumo_guess(2)%matrix => xas_tdp_env%lumo_guess(2)%matrix
      my_dosf = .FALSE.
      IF (PRESENT(do_sf)) my_dosf = do_sf !then either roks or uks

!  Find where/if there are partially occupied MOs
      minel = 1; offset = 0
      IF (do_os) THEN
         minel = MINLOC(nelec_spin,1)
         offset = ABS(nelec_spin(2)-nelec_spin(1))
      END IF

!  Build the guesses based on the LUMOs, starting with the singly occupied MOs
!  For each donor MOs, put the LUMO once at different positions in the guess vector
      iex = 0
      DO imo = 1, offset
         DO ido_mo = 1, ndo_mo
            iex = iex+1
            ex_spin = minel; IF (my_dosf) ex_spin = 3-minel
            CALL cp_fm_to_fm_submat(msource=lumo_guess(minel)%matrix, mtarget=evecs, &
                                    nrow=nao, ncol=1, s_firstrow=1, s_firstcol=imo, &
                                    t_firstrow=((ex_spin-1)*ndo_mo+ido_mo-1)*nao+1, t_firstcol=iex)
         END DO !ido_mo
      END DO !imo

      ! Put the unoccupied LUMOs next
      DO imo = 1, neig/ndo_mo
         DO ispin = 1, nspins
            ilumo = imo; IF (ispin == minel) ilumo = imo+offset
            ex_spin = ispin; IF (my_dosf) ex_spin = 3-ispin
            DO ido_mo = 1, ndo_mo
               iex = iex+1
               IF (iex > neig) EXIT
               CALL cp_fm_to_fm_submat(msource=lumo_guess(ispin)%matrix, mtarget=evecs, &
                                       nrow=nao, ncol=1, s_firstrow=1, s_firstcol=ilumo, &
                                       t_firstrow=((ex_spin-1)*ndo_mo+ido_mo-1)*nao+1, t_firstcol=iex)
            END DO !ido_mo
         END DO !ispin
      END DO !imo

!  Creating a matrix containing the orthogonal space to the solutions (occupied MOs, replicated too)
      nrows = nao*ndo_mo*nspins
      homo = 0
      CALL get_qs_env(qs_env, mos=mos)
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin)%mo_set, homo=homo(ispin))
      END DO
      ncols = SUM(homo)*ndo_mo

      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=nrows, ncol_global=ncols)
      CALL cp_fm_create(ortho_space, fm_struct)

      DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
            DO ido_mo = 1, ndo_mo

               t_firstrow = ((ispin-1)*ndo_mo+ido_mo-1)*nao+1
               IF (my_dosf) t_firstrow = ((2-ispin)*ndo_mo+ido_mo-1)*nao+1

               IF (ispin == 1) t_firstcol = (ido_mo-1)*homo(1)+1
               IF (ispin == 2) t_firstcol = ndo_mo*homo(1) + (ido_mo-1)*homo(2)+1

               CALL cp_fm_to_fm_submat(msource=mo_coeff, mtarget=ortho_space, nrow=nao, &
                                       ncol=homo(ispin), s_firstrow=1, s_firstcol=1, &
                                       t_firstrow=t_firstrow, t_firstcol=t_firstcol)
            END DO !ido_mo
      END DO !ispin

!  Actually solving the eigenvalue problem
      CALL ot_eigensolver(matrix_h=matrix_tdp, matrix_s=metric, matrix_c_fm=evecs, &
                          eps_gradient=eps_iter, iter_max=max_iter, silent=.FALSE., &
                          ot_settings=xas_tdp_control%ot_settings, &
                          matrix_orthogonal_space_fm=ortho_space)
      CALL calculate_subspace_eigenvalues(evecs, matrix_tdp, evals_arg=evals)

!  Clean-up
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_release(ortho_space)
      DEALLOCATE(lumo_guess)

      CALL timestop(handle)

   END SUBROUTINE xas_ot_solver

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
      WHERE (diag .NE. 0.0_dp) diag = 1.0_dp/SQRT(ABS(diag))

      nscal = SIZE(scaling)
      scaling(1:nscal) = diag(1:nscal)

!  Clean-up
      CALL cp_fm_release(work)
      CALL cp_fm_release(fm_norm)
      CALL cp_fm_struct_release(norm_struct)

   END SUBROUTINE get_normal_scaling

! **************************************************************************************************
!> \breif This subroutine computes the row/column block structure as well as the dbcsr ditrinution
!>        for the submatrices making up the generalized XAS TDP eigenvalue problem. They all share
!>        the smae properties, which are based on the replication of the KS matrix
!> \param submat_dist the dbcsr distribution of the submatrices
!> \param submat_blk_size the row/col block sizes for the submatrices
!> \param do_os whther this is a open-shell calculation
!> \param donor_state ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_submat_dist_and_blk_size(submat_dist, submat_blk_size, donor_state, do_os, &
                                               qs_env)

      TYPE(dbcsr_distribution_type), POINTER             :: submat_dist
      INTEGER, DIMENSION(:), POINTER                     :: submat_blk_size
      TYPE(donor_state_type), POINTER                    :: donor_state
      LOGICAL, INTENT(IN)                                :: do_os
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_submat_dist_and_blk_size', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: group, i, nao, nblk_row, ndo_mo, nspins, &
                                                            scol_dist, srow_dist
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, col_dist_sub, row_blk_size, &
                                                            row_dist, row_dist_sub
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks

      NULLIFY (matrix_ks, dbcsr_dist, row_blk_size, row_dist, col_dist, pgrid, col_dist_sub)
      NULLIFY (row_dist_sub)

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
         submat_blk_size((i-1)*nblk_row+1:i*nblk_row) = row_blk_size
         row_dist_sub((i-1)*srow_dist+1:i*srow_dist) = row_dist
         col_dist_sub((i-1)*scol_dist+1:i*scol_dist) = col_dist
      END DO

!  Create the submatrix dbcsr distribution
      ALLOCATE (submat_dist)
      CALL dbcsr_distribution_new(submat_dist, group=group, pgrid=pgrid, row_dist=row_dist_sub, &
                                  col_dist=col_dist_sub)

!  Clean-up
      DEALLOCATE (col_dist_sub, row_dist_sub)

   END SUBROUTINE compute_submat_dist_and_blk_size

! **************************************************************************************************
!> \brief Returns the projector on the unperturbed unoccupied ground state Q = 1 - SP on the block
!>        diagonal of a matrix with the stendard size and distribution.
!> \param proj_Q the matrix with the projector
!> \param dist_q the dbcsr distribution
!> \param blk_size_q the block size
!> \param donor_state ...
!> \param do_os whether it is open-shell calculation
!> \param do_sf whether the projector should be prepared for spin-flip excitations
!> \param xas_tdp_env ...
!> \note In the spin-flip case, the alpha spins are sent to beta and vice-versa. The structure of
!>       the projector is changed by swapping the alpha-alpha with the beta-beta block, which
!>       naturally take the spin change into account. Only for open-shell.
! **************************************************************************************************
   SUBROUTINE get_q_projector(proj_Q, dist_q, blk_size_q, donor_state, do_os, xas_tdp_env, do_sf)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: proj_Q
      TYPE(dbcsr_distribution_type), POINTER             :: dist_q
      INTEGER, DIMENSION(:), POINTER                     :: blk_size_q
      TYPE(donor_state_type), POINTER                    :: donor_state
      LOGICAL, INTENT(IN)                                :: do_os
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_sf

      CHARACTER(len=*), PARAMETER :: routineN = 'get_q_projector', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, blk, iblk, iblk, imo, ispin, jblk, &
                                                            jblk, nblk_row, nblk_row, ndo_mo, &
                                                            ndo_mo, nspins
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size, row_blk_size
      LOGICAL                                            :: found_block, found_block, my_dosf
      REAL(dp), DIMENSION(:), POINTER                    :: work_block, work_block
      TYPE(dbcsr_iterator_type)                          :: iter, iter
      TYPE(dbcsr_type), POINTER                          :: one_sp, one_sp

      NULLIFY (work_block, one_sp, row_blk_size)

!  Initialization
      nspins = 1; IF (do_os) nspins = 2
      ndo_mo = donor_state%ndo_mo
      one_sp => xas_tdp_env%q_projector(1)%matrix
      CALL dbcsr_get_info(one_sp, row_blk_size=row_blk_size)
      nblk_row = SIZE(row_blk_size)
      my_dosf = .FALSE.
      IF (PRESENT(do_sf)) my_dosf = do_sf

      ! the projector is not symmetric
      CALL dbcsr_create(matrix=proj_Q, name="PROJ Q", matrix_type="N", dist=dist_q, &
                        row_blk_size=blk_size_q, col_blk_size=blk_size_q)

!  Fill the projector by looping over 1-SP and duplicating blocks. (all on the spin-MO block diagonal)
      DO ispin = 1, nspins
         one_sp => xas_tdp_env%q_projector(ispin)%matrix

         !if spin-flip, swap the alpha-alpha and beta-beta blocks
         IF (my_dosf) one_sp => xas_tdp_env%q_projector(3-ispin)%matrix

         CALL dbcsr_iterator_start(iter, one_sp)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

            ! get the block
            CALL dbcsr_get_block_p(one_sp, iblk, jblk, work_block, found_block)

            IF (found_block) THEN

               DO imo = 1, ndo_mo
                  CALL dbcsr_put_block(proj_Q, ((ispin-1)*ndo_mo+imo-1)*nblk_row+iblk, &
                                       ((ispin-1)*ndo_mo+imo-1)*nblk_row+jblk, work_block)
               END DO

            END IF
            NULLIFY (work_block)

         END DO !iterator
         CALL dbcsr_iterator_stop(iter)
      END DO !ispin

      CALL dbcsr_finalize(proj_Q)

   END SUBROUTINE get_q_projector

! **************************************************************************************************
!> \brief Builds the matrix containing the ground state contribution to the matrix_tdp (aka matrix A)
!>         => A_{pis,qjt} = (F_pq*delta_ij - epsilon_ij*S_pq) delta_st, where:
!>         F is the KS matrix
!>         S is the overlap matrix
!>         epsilon_ij is the donor MO eigenvalue
!>         i,j labels the MOs, p,q the AOs and s,t the spins
!> \param matrix_a  pointer to a DBCSR matrix containing A
!> \param dist_a the dbcsr distribution of matrix A
!> \param blk_size_a the row/col block sizes of matrix A
!> \param donor_state ...
!> \param do_os ...
!> \param qs_env ...
!> \param do_sf whther the ground state contribution should accomodate spin-flip
!> \note Even localized non-canonical MOs are diagonalized in their subsapce => eps_ij = eps_ii*delta_ij
! **************************************************************************************************
   SUBROUTINE build_gs_contribution(matrix_a, dist_a, blk_size_a, donor_state, do_os, qs_env, do_sf)

      TYPE(dbcsr_type), INTENT(INOUT)                 :: matrix_a
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(dbcsr_distribution_type), POINTER          :: dist_a
      INTEGER, DIMENSION(:), POINTER                  :: blk_size_a
      LOGICAL, INTENT(IN)                             :: do_os
      LOGICAL, INTENT(IN), OPTIONAL                   :: do_sf
      
      CHARACTER(len=*), PARAMETER :: routineN = 'build_gs_contribution', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, iblk, imo, ispin, jblk, nblk_row, &
                                                            ndo_mo, nspins
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size
      LOGICAL                                            :: found_block, my_dosf
      REAL(dp), DIMENSION(:), POINTER                    :: energy_evals, work_block
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: m_ks, matrix_ks, matrix_s
      TYPE(dbcsr_type)                                   :: work_matrix

      NULLIFY (matrix_ks, dbcsr_dist, row_blk_size, work_block, energy_evals, matrix_s, m_ks)

      !  Note: matrix A is symmetric and block diagonal.

!  Initialization
      nspins = 1; IF (do_os) nspins = 2
      ndo_mo = donor_state%ndo_mo
      CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size)
      nblk_row = SIZE(row_blk_size)

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
                  CALL dbcsr_put_block(matrix_a, ((ispin-1)*ndo_mo+imo-1)*nblk_row+iblk, &
                                       ((ispin-1)*ndo_mo+imo-1)*nblk_row+jblk, work_block)

               END DO !imo
            END IF !found_block
            NULLIFY (work_block)
         END DO ! iteration on KS blocks
         CALL dbcsr_iterator_stop(iter)

!     Loop over the blocks of S and put them on the block diagonal of work
         energy_evals => donor_state%energy_evals(:,ispin)

         CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

!           Get the block
            CALL dbcsr_get_block_p(matrix_s(1)%matrix, iblk, jblk, work_block, found_block)

            IF (found_block) THEN

!              Add S matrix on block diagonal as epsilon_ii*S_pq
               DO imo = 1, ndo_mo

                  CALL dbcsr_put_block(work_matrix, ((ispin-1)*ndo_mo+imo-1)*nblk_row+iblk, &
                                       ((ispin-1)*ndo_mo+imo-1)*nblk_row+jblk, &
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

   END SUBROUTINE build_gs_contribution

! **************************************************************************************************
!> \brief Creates the metric (aka  matrix G) needed for the generalized eigenvalue problem and inverse
!>         => G_{pis,qjt} = S_pq*delta_ij*delta_st
!> \param matrix_g dbcsr matrix containing G
!> \param dist_g dbcsr distribution of G
!> \param blk_size_g block sizes of G
!> \param donor_state ...
!> \param qs_env ...
!> \param do_inv if the inverse of G should be computed
!> \param do_os if open-shell calculation
! **************************************************************************************************
   SUBROUTINE build_metric(matrix_g, dist_g, blk_size_g, donor_state, qs_env, do_os, do_inv)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_g
      TYPE(dbcsr_distribution_type), POINTER             :: dist_g
      INTEGER, DIMENSION(:), POINTER                     :: blk_size_g
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: do_os
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_inv

      CHARACTER(len=*), PARAMETER :: routineN = 'build_metric', routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, i, iblk, jblk, nao, nblk_row, &
                                                            ndo_mo, nspins
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size
      LOGICAL                                            :: found_block, my_do_inv
      REAL(dp), DIMENSION(:), POINTER                    :: work_block
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: matrix_sinv

      NULLIFY (matrix_s, row_blk_size, work_block, para_env, blacs_env)

!  TODO: this operation of taking a normal size matrix and replicating it on all blocks is quite
!        common, maybe shoudl refactor it ?!?

!  Initilization
      nspins = 1; IF (do_os) nspins = 2
      ndo_mo = donor_state%ndo_mo
      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size, nfullrows_total=nao)
      nblk_row = SIZE(row_blk_size)
      my_do_inv = .FALSE.
      IF (PRESENT(do_inv)) my_do_inv = do_inv

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
               CALL dbcsr_put_block(matrix_g(1)%matrix, (i-1)*nblk_row+iblk, (i-1)*nblk_row+jblk, work_block)
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
                  CALL dbcsr_put_block(matrix_g(2)%matrix, (i-1)*nblk_row+iblk, (i-1)*nblk_row+jblk, work_block)
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

   END SUBROUTINE build_metric

! **************************************************************************************************
! \brief Builds the auxiliary matrix (A-D+E)^+0.5 needed for the transofrmation of the
!>       full-TDDFT problem into an Hermitian one
!> \param threshold a threshold for allowed negative eigenvalues
!> \param sx the amount of exact exchange
!> \param matrix_a the ground state contribution matrix A
!> \param matrix_d the on-diagonal exchange kernel matrix (ab|IJ)
!> \param matrix_e the off-diagonal exchange kernel matrix (aJ|Ib)
!> \param do_hfx if exact exchange is included
!> \param proj_Q ...
!> \param work ...
!> \param donor_state ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE build_aux_matrix(threshold, sx, matrix_a, matrix_d, matrix_e, do_hfx, proj_Q, &
                               work, donor_state, qs_env)

      REAL(dp), INTENT(IN)                               :: threshold, sx
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_a, matrix_d, matrix_e
      LOGICAL, INTENT(IN)                                :: do_hfx
      TYPE(dbcsr_type), INTENT(INOUT)                    :: proj_Q, work
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_aux_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ndep, nrow
      REAL(dp)                                           :: evals(2)
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: struct
      TYPE(cp_fm_type), POINTER                          :: aux_fm, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type)                                   :: tmp_mat

      NULLIFY (aux_fm, work_fm, struct, blacs_env, para_env)

      CALL dbcsr_copy(tmp_mat, matrix_a)
      IF (do_hfx) THEN
         CALL dbcsr_add(tmp_mat, matrix_d, 1.0_dp, -1.0_dp*sx) !scaled hfx
         CALL dbcsr_add(tmp_mat, matrix_e, 1.0_dp, 1.0_dp*sx)
      END IF

      ! Take the product with the Q projector:
      CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, tmp_mat, 0.0_dp, work)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, tmp_mat)

      ! Actually computing and storing the auxiliary matrix
      ALLOCATE (donor_state%matrix_aux)
      CALL dbcsr_create(matrix=donor_state%matrix_aux, template=matrix_a, name="MAT AUX")

      ! Go through fm and cp_fm_power for good quality sqrt
      CALL dbcsr_get_info(matrix_a, nfullrows_total=nrow)
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      CALL cp_fm_struct_create(struct, context=blacs_env, para_env=para_env, &
                               ncol_global=nrow, nrow_global=nrow)
      CALL cp_fm_create(aux_fm, struct)
      CALL cp_fm_create(work_fm, struct)

      CALL copy_dbcsr_to_fm(tmp_mat, aux_fm)
      CALL cp_fm_power(aux_fm, work_fm, 0.5_dp, threshold, ndep, eigvals=evals)

      CALL copy_fm_to_dbcsr(aux_fm, donor_state%matrix_aux)

      ! Warn the user if matrix not positive semi-definite
      IF (evals(1) < 0.0_dp .AND. ABS(evals(1)) > threshold) THEN
         CPWARN("The full TDDFT problem might not have been soundly turned Hermitian. Try TDA.")
      END IF

      ! clean-up
      CALL cp_fm_struct_release(struct)
      CALL cp_fm_release(aux_fm)
      CALL cp_fm_release(work_fm)
      CALL dbcsr_release(tmp_mat)

   END SUBROUTINE build_aux_matrix

! **************************************************************************************************
!> \brief Computes the two-center Coulomb integral needed for the RI in kernel calculation. Stores
!>        the integrals (P|Q)^-1 in the xas_tdp_env as global (small) arrays. Does that for a given
!>        excited kind
!> \param ex_kind ...
!> \param xas_tdp_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_ri_coulomb2_int(ex_kind, xas_tdp_env, xas_tdp_control, qs_env)

      INTEGER, INTENT(IN)                             :: ex_kind
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env 

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_coulomb2_int', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: nsgf
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (ri_basis, qs_kind_set)

!  Initialization
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      IF (ASSOCIATED(xas_tdp_env%ri_inv_coul)) THEN
         DEALLOCATE (xas_tdp_env%ri_inv_coul)
      END IF

!  Get the RI basis of interest and its quantum numbers
      CALL get_qs_kind(qs_kind_set(ex_kind), basis_set=ri_basis, basis_type="RI_XAS")
      CALL get_gto_basis_set(ri_basis, nsgf=nsgf) 
      ALLOCATE(xas_tdp_env%ri_inv_coul(nsgf, nsgf))
      xas_tdp_env%ri_inv_coul = 0.0_dp

      !IF (.NOT. xas_tdp_control%is_periodic) THEN
         CALL int_operators_r12_ab_os(r12_operator=operator_coulomb, vab=xas_tdp_env%ri_inv_coul, &
                                      rab=(/0.0_dp,0.0_dp,0.0_dp/), fba=ri_basis, fbb=ri_basis, &
                                      calculate_forces=.FALSE.)
         CPASSERT(ASSOCIATED(xas_tdp_control))
      !ELSE
         !CALL periodic_ri_coulomb2(xas_tdp_env%ri_inv_coul, ri_basis, qs_env)
      !END IF

!  Inverting 
      CALL invmat_symm(xas_tdp_env%ri_inv_coul)  

   END SUBROUTINE compute_ri_coulomb2_int

! **************************************************************************************************
!> \brief Computes the two-center inverse coulomb integral in the case of PBCs
!> \param ri_coul2 ...
!> \param ri_absis ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE periodic_ri_coulomb2(ri_coul2, ri_basis, qs_env)

      REAL(dp), DIMENSION(:,:), INTENT(INOUT)         :: ri_coul2
      TYPE(gto_basis_set_type), POINTER               :: ri_basis
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "periodic_ri_coulomb2", &
                                     routinep = moduleN//":"//routineN

      INTEGER                                         :: op, oq, ncop, ncoq, nset, maxco, &
                                                         pset, qset, sgfp, sgfq, ppgf, qpgf
      REAL(dp)                                        :: r(3)
      TYPE(cp_eri_mme_param)                          :: mme_param
      TYPE(cell_type), POINTER                        :: cell
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(cp_para_env_type), POINTER                 :: para_env
      INTEGER, DIMENSION(:), POINTER                  :: lmax, lmin, nsgf, npgf
      REAL(dp), DIMENSION(:,:), POINTER               :: zet, sphi
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgf
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: hpq

      NULLIFY(cell, qs_kind_set, lmin, lmax, nsgf, npgf, zet, sphi, first_sgf)

      ! Use eri_mme for this. Don't want to annoy the user with a full input section just for this
      ! tiny bit => initialize our own eri_mme section with the defaults

      CALL get_qs_env(qs_env, cell=cell, qs_kind_set=qs_kind_set, para_env=para_env)

      CALL eri_mme_init(mme_param%par, n_minimax=20, cutoff=300.0_dp, do_calib_cutoff=.TRUE., &
                        cutoff_min=10.0_dp, cutoff_max=10000.0_dp, cutoff_eps=0.01_dp, &
                        cutoff_delta=0.9_dp, sum_precision=1.0E-12_dp, debug=.FALSE., &
                        debug_delta=1.0E-6_dp, debug_nsum=1000000, unit_nr=0, print_calib=.FALSE.)

      CALL cp_eri_mme_set_params(mme_param, cell, qs_kind_set, basis_type_1="RI_XAS", para_env=para_env)

      CALL get_gto_basis_set(ri_basis, lmax=lmax, npgf=npgf, zet=zet, lmin=lmin, nset=nset, &
                             nsgf_set=nsgf, sphi=sphi, first_sgf=first_sgf, maxco=maxco)

      r = 0.0_dp
      ALLOCATE(hpq(nset*maxco, nset*maxco))
      hpq = 0.0_dp

      DO pset = 1, nset
         ncop = npgf(pset)*ncoset(lmax(pset))
         sgfp = first_sgf(1, pset)

         DO qset = 1, nset
            ncoq = npgf(qset)*ncoset(lmax(qset))
            sgfq = first_sgf(1, qset)

            DO ppgf = 1, npgf(pset)
               op = (pset-1)*maxco + (ppgf-1)*ncoset(lmax(pset))
               DO qpgf = 1, npgf(qset)
                  oq = (qset-1)*maxco + (qpgf-1)*ncoset(lmax(qset))

                  CALL eri_mme_2c_integrate(mme_param%par, lmin(pset), lmax(pset), lmin(qset), &
                                            lmax(qset), zet(ppgf,pset), zet(qpgf,qset), r, hpq, &
                                            op, oq)

               END DO !qpgf
            END DO ! ppgf

            !contraction into sgfs
            op = (pset-1)*maxco+1
            oq = (qset-1)*maxco+1

            CALL ab_contract(ri_coul2(sgfp:sgfp+nsgf(pset)-1, sgfq:sgfq+nsgf(qset)-1), &
                             hpq(op:op+ncop-1, oq:oq+ncoq-1), sphi(:,sgfp:), sphi(:,sgfq:), &
                             ncop, ncoq, nsgf(pset), nsgf(qset))

         END DO !qset
      END DO !pset

      !celan-up
      CALL eri_mme_release(mme_param%par)

   END SUBROUTINE periodic_ri_coulomb2

! **************************************************************************************************
!> \brief Computes the RI Coulomb 3-center integrals (ab|c), where c is from the RI_XAS basis and
!>        centered on the excited atoms of xas_tdp_env
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note  The ri_o3c_coul member of xas_tdp_env is defined and allocated here. Only computed once
!>        for the whole system (for optimized load balance). Ok because not too much memory needed
! **************************************************************************************************
   SUBROUTINE compute_ri_o3c_coulomb(xas_tdp_env, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_o3c_coulomb', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ilist, nkind
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_orb, basis_set_ri
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
                                             POINTER  :: ac_list, ab_list
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(cp_para_env_type), POINTER                 :: para_env

      NULLIFY(basis_set_ri, basis_set_orb, ac_list, ab_list, qs_kind_set, para_env)

      CALL timeset(routineN, handle)

!  Take what we need from the qs_env
      CALL get_qs_env(qs_env, nkind=nkind, qs_kind_set=qs_kind_set, para_env=para_env)

!  Build the basis set lists
      ALLOCATE (basis_set_ri(nkind))
      ALLOCATE (basis_set_orb(nkind))
      CALL basis_set_list_setup(basis_set_ri, "RI_XAS", qs_kind_set)
      CALL basis_set_list_setup(basis_set_orb, "ORB", qs_kind_set)

!  Build a neighbor list for the ab centers. Assume (aI|c) = sum_b c_bI (ab|c), with c_bI only
!  non-zero for b centered on the same atom as c => build overlap nl, but only keeping b if centered
!  on an excited atom
      CALL build_xas_tdp_ovlp_nl(ab_list, basis_set_orb, basis_set_orb, &
                                 qs_env, excited_atoms=xas_tdp_env%ex_atom_indices)

!  Build a neighbor list for the ac centers. Since we later contract as (aI|c) and we assume I is 
!  very localized on the same atom as c, we take a,c as neighbors if they overlap
      CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, "OVERLAP", &
                               qs_env, excited_atoms=xas_tdp_env%ex_atom_indices)

!  Allocate, init and compute the integrals
      ALLOCATE (xas_tdp_env%ri_o3c_coul)
      !only one spin needed for copulomb integrals
      CALL init_o3c_container(xas_tdp_env%ri_o3c_coul, 1, basis_set_orb, basis_set_orb, &
                              basis_set_ri, ab_list, ac_list, only_bc_same_center=.TRUE.)
      CALL calculate_o3c_coulomb_integrals(xas_tdp_env%ri_o3c_coul)

! Clean-up
      DO ilist = 1,SIZE(ab_list)
         CALL deallocate_neighbor_list_set(ab_list(ilist)%neighbor_list_set)
      END DO
      DO ilist = 1,SIZE(ac_list)
         CALL deallocate_neighbor_list_set(ac_list(ilist)%neighbor_list_set)
      END DO
      DEALLOCATE(basis_set_ri, basis_set_orb, ab_list, ac_list)

      !not strictly necessary but avoid having any load unbalance here being reported in the 
      !timings for other routines
      CALL mp_sync(para_env%group)

      CALL timestop(handle)

   END SUBROUTINE compute_ri_o3c_coulomb

! **************************************************************************************************
!> \brief Computes the RI exchange 3-center integrals (ab|c), where c is from the RI_XAS basis and
!>        centered on excited atoms and kind
!> \param ex_atom the excited atom on which the third center is
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note  This routine is called once for each excited atom. Because there are many different a,b
!>        pairs involved, load balance is ok. This allows memory saving
! **************************************************************************************************
   SUBROUTINE compute_ri_o3c_exchange(ex_atom, xas_tdp_env, xas_tdp_control, qs_env)

      INTEGER, INTENT(IN)                             :: ex_atom
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_o3c_exchange', &
         routineP = moduleN//':'//routineN

      INTEGER                                         :: handle, nkind, ilist
      INTEGER, DIMENSION(:), ALLOCATABLE              :: exat
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_set_ri, basis_set_orb
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
                                             POINTER  :: ac_list, sab_orb
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(cp_para_env_type), POINTER                 :: para_env

      NULLIFY(basis_set_ri, basis_set_orb, ac_list, sab_orb, qs_kind_set, para_env)

      CALL timeset(routineN, handle)

      !TODO: give the choice to the user to compute all integrals at once ?!?
      !      maybe do that if it turns out that load imbalance gets really bad
      
!  Take what we need from the qs_env
      CALL get_qs_env(qs_env, nkind=nkind, sab_orb=sab_orb, qs_kind_set=qs_kind_set, para_env=para_env)

!  Build the basis set lists
      ALLOCATE (basis_set_ri(nkind))
      ALLOCATE (basis_set_orb(nkind))
      CALL basis_set_list_setup(basis_set_ri, "RI_XAS", qs_kind_set)
      CALL basis_set_list_setup(basis_set_orb, "ORB", qs_kind_set)

!  Build a neighbor list for the ac centers.
      ALLOCATE(exat(1))
      exat(1) = ex_atom
      !test on the exchange operator
      IF (xas_tdp_control%x_op == do_hfx_potential_coulomb) THEN
         CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, "COULOMB", &
                                  qs_env, excited_atoms=exat)
      ELSE IF (xas_tdp_control%x_op == do_hfx_potential_truncated) THEN
         CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, "TRUNCATED COULOMB", &
                                  qs_env, excited_atoms=exat, x_range=xas_tdp_control%x_range)
      END IF

!  Allocate, init and compute the integrals. Use stendard sab_orb because in (ab|IJ), the only
!  restriction is that a and b overlap
      ALLOCATE(xas_tdp_env%ri_o3c_ex)
      CALL init_o3c_container(xas_tdp_env%ri_o3c_ex, 1, basis_set_orb, basis_set_orb, &
                              basis_set_ri, sab_orb, ac_list)
      CALL calculate_o3c_coulomb_integrals(xas_tdp_env%ri_o3c_ex)

! Clean-up
      DO ilist = 1, SIZE(ac_list)
         CALL deallocate_neighbor_list_set(ac_list(ilist)%neighbor_list_set)
      END DO
      DEALLOCATE (basis_set_ri, basis_set_orb, ac_list)

      !not strictly necessary but avoid having any load unbalance here being reported in the 
      !timings for other routines
      CALL mp_sync(para_env%group)

      CALL timestop(handle)

   END SUBROUTINE compute_ri_o3c_exchange

! **************************************************************************************************
!> \brief Builds a neighbor lists set for overlaping 2-center S_ab, where b is restricted on a 
!>        a given list of atoms. Used for Coulomb RI where (aI|P) = sum_b C_bI (ab|P), where 
!>        contraction coeff C_bI is assumed to be non-zero only on excited atoms
!> \param ab_list the neighbor list 
!> \param basis_a basis set list for atom a
!> \param basis_b basis set list for atom b
!> \param excited_atoms the indices of the excited atoms on which b is centered
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE build_xas_tdp_ovlp_nl(ab_list, basis_a, basis_b, qs_env, excited_atoms)

      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
                                             POINTER  :: ab_list
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_a, basis_b
      TYPE(qs_environment_type), POINTER              :: qs_env
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL     :: excited_atoms

      CHARACTER(len=*), PARAMETER :: routineN = "build_xas_tdp_ovlp_nl", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: nkind, ikind
      REAL(dp)                                        :: subcells
      LOGICAL                                         :: my_restrictb
      LOGICAL, DIMENSION(:), ALLOCATABLE              :: a_present, b_present
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: a_radius, b_radius
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
      TYPE(cell_type), POINTER                        :: cell
      TYPE(distribution_1d_type), POINTER             :: distribution_1d
      TYPE(distribution_2d_type), POINTER             :: distribution_2d
      TYPE(local_atoms_type), DIMENSION(:), &
                                          ALLOCATABLE :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER      :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set

      NULLIFY(atomic_kind_set, distribution_1d, distribution_2d, molecule_set, particle_set, cell)

!  Initialization
      CALL get_qs_env(qs_env, nkind=nkind)
      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)

      my_restrictb = .FALSE.
      IF (PRESENT(excited_atoms)) THEN
         my_restrictb = .TRUE.
      END IF

      ALLOCATE(a_present(nkind), b_present(nkind))
      a_present = .FALSE.
      b_present = .FALSE.
      ALLOCATE(a_radius(nkind), b_radius(nkind))
      a_radius = 0.0_dp
      b_radius = 0.0_dp

!  Set up the radii
      DO ikind = 1,nkind
         IF (ASSOCIATED(basis_a(ikind)%gto_basis_set)) THEN
            a_present(ikind) = .TRUE.
            CALL get_gto_basis_set(basis_a(ikind)%gto_basis_set, kind_radius=a_radius(ikind))
         END IF

         IF (ASSOCIATED(basis_b(ikind)%gto_basis_set)) THEN
            b_present(ikind) = .TRUE.
            CALL get_gto_basis_set(basis_b(ikind)%gto_basis_set, kind_radius=b_radius(ikind))
         END IF
      END DO !ikind

      ALLOCATE(pair_radius(nkind, nkind))
      pair_radius = 0.0_dp
      CALL pair_radius_setup(a_present, b_present, a_radius, b_radius, pair_radius)

!  Set up the nl
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                      distribution_2d=distribution_2d, local_particles=distribution_1d, &
                      particle_set=particle_set, molecule_set=molecule_set)

      ALLOCATE(atom2d(nkind))
      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, .FALSE., particle_set)

      IF (my_restrictb) THEN

         CALL build_neighbor_lists(ab_list, particle_set, atom2d, cell, pair_radius, subcells, &
                                   atomb_to_keep=excited_atoms, nlname="XAS_TDP_ovlp_nl")

      ELSE

         CALL build_neighbor_lists(ab_list, particle_set, atom2d, cell, pair_radius, subcells, &
                                   nlname="XAS_TDP_ovlp_nl")

      END IF
!  Clean-up
      CALL atom2d_cleanup(atom2d)

   END SUBROUTINE build_xas_tdp_ovlp_nl
   
! **************************************************************************************************
!> \brief Builds a neighbor lists set taylored for 3-center integral within XAS TDP, such that only
!>        excited atoms are taken into account for the list_c
!> \param ac_list the neighbor list ready for 3-center integrals
!> \param basis_a basis set list for atom a
!> \param basis_c basis set list for atom c
!> \param op_type to indicate whther the list should be built with overlap or Coulomb in mind
!> \param excited_atoms the indices of the excited atoms to consider (if not given, all atoms are taken)
!> \param x_range in case some truncated/screened operator is used, gives its range
!> \param qs_env ...
!> \note Based on setup_neighbor_list with added features
! **************************************************************************************************
   SUBROUTINE build_xas_tdp_3c_nl(ac_list, basis_a, basis_c, op_type, qs_env, excited_atoms, x_range)

      TYPE(neighbor_list_set_p_type), DIMENSION(:), &                                                
                                             POINTER  :: ac_list
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_a, basis_c
      CHARACTER(len=*), INTENT(IN)                    :: op_type
      TYPE(qs_environment_type), POINTER              :: qs_env
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL     :: excited_atoms
      REAL(dp), INTENT(IN), OPTIONAL                  :: x_range

      CHARACTER(len=*), PARAMETER :: routineN =  "build_xas_tdp_3c_nl", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: nkind, ikind
      LOGICAL                                         :: sort_atoms
      REAL(dp)                                        :: subcells
      LOGICAL, DIMENSION(:), ALLOCATABLE              :: a_present, c_present
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: a_radius, c_radius
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
      TYPE(cell_type), POINTER                        :: cell
      TYPE(distribution_1d_type), POINTER             :: distribution_1d
      TYPE(distribution_2d_type), POINTER             :: distribution_2d
      TYPE(local_atoms_type), DIMENSION(:), &
                                          ALLOCATABLE :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER      :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      CHARACTER(2)                                    :: my_op

      NULLIFY(atomic_kind_set, distribution_1d, distribution_2d, molecule_set, particle_set, cell)
      
!  By default, take overlap operator
      my_op = "O "
      IF (op_type == "COULOMB") my_op = "C "
      IF (op_type == "TRUNCATED COULOMB") THEN
         my_op = "TC"
         IF (.NOT. PRESENT(x_range)) CPABORT("Need a range for the truncated Coulomb operator")
      END IF

!  Initialization
      CALL get_qs_env(qs_env, nkind=nkind)
      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)
      sort_atoms = .FALSE.
      IF ((PRESENT(excited_atoms))) sort_atoms = .TRUE.

      ALLOCATE (a_present(nkind), c_present(nkind))
      a_present = .FALSE.
      c_present = .FALSE.
      ALLOCATE (a_radius(nkind), c_radius(nkind))
      a_radius = 0.0_dp
      c_radius = 0.0_dp

!  Set up the radii, depending on the operator type
      IF (my_op == "O ") THEN
         !overlap => use the kind radius for both a and c
         DO ikind = 1,nkind
            !orbital basis set
            IF (ASSOCIATED(basis_a(ikind)%gto_basis_set)) THEN
               a_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_a(ikind)%gto_basis_set, kind_radius=a_radius(ikind))
            END IF
            !RI_XAS basis set
            IF (ASSOCIATED(basis_c(ikind)%gto_basis_set)) THEN
               c_present(ikind) = .TRUE.
               CALL get_gto_basis_set(basis_c(ikind)%gto_basis_set, kind_radius=c_radius(ikind))
            END IF
         END DO !ikind

      ELSE IF (my_op == "C ") THEN
         !Coulomb operator, virtually infinite range => set c_radius to arbitrarily large number
         DO ikind = 1,nkind
            IF (ASSOCIATED(basis_c(ikind)%gto_basis_set)) THEN
               c_present(ikind) = .TRUE.
               c_radius(ikind) = 1000000.0_dp
            END IF
            IF (ASSOCIATED(basis_a(ikind)%gto_basis_set)) a_present(ikind) = .TRUE.
         END DO !ikind

      ELSE
         !Truncated coulomb: set c_radius to x_range and leave a_radius at 0
         DO ikind = 1, nkind
            IF (ASSOCIATED(basis_c(ikind)%gto_basis_set)) THEN
               c_present(ikind) = .TRUE.
               c_radius(ikind) = x_range
            END IF
            IF (ASSOCIATED(basis_a(ikind)%gto_basis_set)) a_present(ikind) = .TRUE.
         END DO !ikind
      END IF

      ALLOCATE (pair_radius(nkind, nkind))
      pair_radius = 0.0_dp
      CALL pair_radius_setup(a_present, c_present, a_radius, c_radius, pair_radius)

!  Actually setup the list
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                      distribution_2d=distribution_2d, local_particles=distribution_1d, &
                      particle_set=particle_set, molecule_set=molecule_set)

      ALLOCATE (atom2d(nkind))
      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, .FALSE., particle_set)

      IF (sort_atoms) THEN
         CALL build_neighbor_lists(ac_list, particle_set, atom2d, cell, pair_radius, subcells, &
                                   operator_type = "ABC", atomb_to_keep=excited_atoms, &
                                   nlname="XAS_TDP_3c_nl")
      ELSE
         CALL build_neighbor_lists(ac_list, particle_set, atom2d, cell, pair_radius, subcells, &
                                   operator_type = "ABC", nlname="XAS_TDP_3c_nl")
      END IF

!  Clean-up
      CALL atom2d_cleanup(atom2d)

   END SUBROUTINE build_xas_tdp_3c_nl

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

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "include_os_soc", routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, nsc, nsf, ntot, isc, isf, ndo_mo, &
                                                         ndo_so, nao, homo, jsc, jsf, i
      REAL(dp)                                        :: gs_sum, soc
      LOGICAL                                         :: do_roks, do_uks
      REAL(dp), DIMENSION(:), POINTER                 :: sc_evals, sf_evals
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: diag, tmp_evals
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: domo_soc_x, domo_soc_y, domo_soc_z, &
                                                         S_ex, SOC_ex
      TYPE(cp_fm_type), POINTER                       :: gs_coeffs, real_fm, img_fm, mo_coeff, &
                                                         vec_work, prod_work, vec_s, vec_soc_x, &
                                                         vec_soc_y, vec_soc_z
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: sc_coeffs, sf_coeffs, amew_dip
      TYPE(cp_fm_struct_type), POINTER                :: full_struct, vec_struct, prod_struct, &
                                                         dip_struct
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(dbcsr_type), POINTER                       :: orb_soc_x, orb_soc_y, orb_soc_z
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      TYPE(cp_cfm_type), POINTER                      :: pert_cfm, evecs_cfm, work_cfm, dip_cfm
      COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE        :: transdip

      NULLIFY(gs_coeffs, sc_coeffs, sf_coeffs, matrix_s, orb_soc_x, orb_soc_y, orb_soc_z, mos)
      NULLIFY(real_fm, img_fm, full_struct, para_env, blacs_env, mo_coeff, vec_s, vec_soc_x)
      NULLIFY(sc_evals, sf_evals, vec_work, prod_work, vec_struct, prod_struct, vec_soc_y)
      NULLIFY(vec_soc_z, pert_cfm, evecs_cfm, work_cfm, dip_cfm, dip_struct, amew_dip)

      CALL timeset(routineN, handle)

! Initialization
      sc_coeffs => donor_state%sc_coeffs
      sf_coeffs => donor_state%sf_coeffs
      sc_evals => donor_state%sc_evals
      sf_evals => donor_state%sf_evals
      nsc = SIZE(sc_evals)
      nsf = SIZE(sf_evals)
      ntot = 1+nsc+nsf
      ndo_mo = donor_state%ndo_mo
      ndo_so = 2*ndo_mo
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env, mos=mos, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      orb_soc_x => xas_tdp_env%orb_soc(1)%matrix
      orb_soc_y => xas_tdp_env%orb_soc(2)%matrix
      orb_soc_z => xas_tdp_env%orb_soc(3)%matrix
      do_roks = xas_tdp_control%do_roks
      do_uks = xas_tdp_control%do_uks

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
                                 t_firstcol=ndo_mo+1)

         CALL cp_fm_struct_release(vec_struct)
      END IF

! Creating the real and the imaginary part of the SOC perturbation matrix
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      CALL cp_fm_create(real_fm, full_struct)
      CALL cp_fm_create(img_fm, full_struct)

! Put the excitation energies on the diagonal of the real  matrix. Element 1,1 is the ground state
      DO isc = 1, nsc
         CALL cp_fm_set_element(real_fm, 1+isc, 1+isc, sc_evals(isc))
      END DO
      DO isf = 1, nsf
         CALL cp_fm_set_element(real_fm, 1+nsc+isf, 1+nsc+isf, sf_evals(isf))
      END DO

! Precompute what we can before looping over excited states.
      ! Need to compute: <phi^0_Isigma|SOC|phi^0_Jtau> for the donor MOs

      CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nao, ncol_global=ndo_so)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_so, ncol_global=ndo_so)
      CALL cp_fm_create(vec_work, vec_struct)
      CALL cp_fm_create(prod_work, prod_struct)

      ALLOCATE (domo_soc_x(ndo_so, ndo_so), domo_soc_y(ndo_so, ndo_so), domo_soc_z(ndo_so, ndo_so))

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_x, gs_coeffs, vec_work, ncol=ndo_so)
      CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, vec_work, 0.0_dp, prod_work)
      CALL cp_fm_get_submatrix(prod_work, domo_soc_x)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, gs_coeffs, vec_work, ncol=ndo_so)
      CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, vec_work, 0.0_dp, prod_work)
      CALL cp_fm_get_submatrix(prod_work, domo_soc_y)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, gs_coeffs, vec_work, ncol=ndo_so)
      CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, vec_work, 0.0_dp, prod_work)
      CALL cp_fm_get_submatrix(prod_work, domo_soc_z)

      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(vec_struct)
      CALL cp_fm_struct_release(prod_struct)

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
      gs_sum = gs_sum-SUM(diag) ! -1 because of spin integration

      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)
      CALL cp_fm_struct_release(prod_struct)
      DEALLOCATE (diag)
      NULLIFY (vec_struct)

      ! Create some work structures and matrices where intermediate results will be stored
      CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nao, ncol_global=ndo_so)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_so, ncol_global=ndo_so)
      CALL cp_fm_create(vec_s, vec_struct) ! for S*coeffs storage
      CALL cp_fm_create(vec_soc_x, vec_struct) ! for SOC_x*coeffs
      CALL cp_fm_create(vec_soc_y, vec_struct) ! for SOC_y*coeffs
      CALL cp_fm_create(vec_soc_z, vec_struct) ! for SOC_z*coeffs
      CALL cp_fm_create(prod_work, prod_struct) ! for anything like coeffs^T * MAT * coeffs
      ALLOCATE (S_ex(ndo_so, ndo_so)) ! Overlap between LR coeffs
      ALLOCATE (SOC_ex(ndo_so, ndo_so)) ! SOC between LR coeffs
      ALLOCATE (diag(ndo_so))

! Looping over the excited states, computing the SOC and filling the perturbation matrix
! There are 3 loops to do: sc-sc, sc-sf and sf-sf

!  Do the spin-conserving spin-conserving coupling loop. Block is on the diagonal and the end
!  perturbation matrix is Hermitian => only loop over the upper diagonal
      DO jsc = 1, nsc
         ! For each SC, need to couple it to the ground state SD: <Psi_0|SOC|Psi_Jsc>
         ! Compute SOC_z|psi^Jc_k,sigma>, which we store as it will be needed afterwards
         CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, sc_coeffs(jsc)%matrix, vec_soc_z, ncol=ndo_so)

         !<Psi_0|SOC|Psi_Jsc> = sum_k,sigma <phi^0_k,sigma|SOC|phi^Jsc_k,sigma>
         CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, vec_soc_z, 0.0_dp, prod_work)
         CALL cp_fm_get_diag(prod_work, diag)
         !same spin integration => <alpha|sigma_z|alpha> = 1, <beta|sigma_z|beta> = -1
         soc = SUM(diag(1:ndo_mo))-SUM(diag(ndo_mo+1:ndo_so))

         !purely imaginary element
         CALL cp_fm_set_element(img_fm, 1, 1+jsc, soc)

         !we also precompute and store S|phi^Jsc_k,sigma>
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, sc_coeffs(jsc)%matrix, vec_s, ncol=ndo_so)

         DO isc = 1, jsc

            ! <Psi_Isc|SOC|Psi_Jsc> =
            ! sum_k,sigma [<psi^Isc_k,sigma|SOC|psi^Jsc_k,sigma> + <psi^Isc_k,sigma|psi^Jsc_k,sigma> * gs_sum]
            ! - sum_k,l,sigma <psi^0_k,sigma|SOC|psi^0_l,sigma> * <psi^Isc_l,sigma|psi^Jsc_k,sigma>

            ! Compute the LR orbitals overlap
            CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, sc_coeffs(isc)%matrix, vec_s, 0.0_dp, prod_work)
            CALL cp_fm_get_submatrix(prod_work, S_ex)

            ! Compute the LR SOC: same-spin integration => only need SOC_z
            CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, sc_coeffs(isc)%matrix, vec_soc_z, 0.0_dp, prod_work)
            CALL cp_fm_get_submatrix(prod_work, SOC_ex)

            ! Actually compute the SOC, with <alpha|sigma_z|alpha> = 1, <beta|sigma_z|beta> = -1
            diag(:) = get_diag(SOC_ex)
            soc = SUM(diag(1:ndo_mo))-SUM(diag(ndo_mo+1:ndo_so))

            diag(:) = get_diag(S_ex)
            soc = soc+SUM(diag)*gs_sum

            soc = soc-SUM(S_ex(1:ndo_mo, 1:ndo_mo)*TRANSPOSE(domo_soc_z(1:ndo_mo, 1:ndo_mo))) !alpha
            soc = soc+SUM(S_ex(ndo_mo+1:ndo_so, ndo_mo+1:ndo_so)* & !beta
                          TRANSPOSE(domo_soc_z(ndo_mo+1:ndo_so, ndo_mo+1:ndo_so)))

            !purely imaginary element
            CALL cp_fm_set_element(img_fm, 1+isc, 1+jsc, soc)

         END DO !isc
      END DO !jsc

!  Do the spin-flip spin-flip loop followed by the spin-conserving spin-flip loop
!  Reminder: for the SF LR orbitals, the alpha-beta transition coeffs sit at the usual alpha spot
!            and the beta-alpha transition coeffs sit the the usual beta spot
      DO jsf = 1, nsf
         ! For each SF, need to couple it to the ground state SD: <Psi_0|SOC|Psi_Jsf>
         ! Need to compute both SOC_x|psi^Jsf_k,sigma> and SOC_y|psi^Jsf_k,sigma>
         CALL cp_dbcsr_sm_fm_multiply(orb_soc_x, sf_coeffs(jsf)%matrix, vec_soc_x, ncol=ndo_so)
         CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, sf_coeffs(jsf)%matrix, vec_soc_y, ncol=ndo_so)

         !<Psi_0|SOC|Psi_Jsf> = sum_k,sigma <phi^0_k,sigma|SOC|phi^Jsc_k,tau>   sigma != tau
         ! the imaginary SOC_x part, both <alpha|sigma_x|beta> and <beta|sigma_x|alpha> = 1
         CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, vec_soc_x, 0.0_dp, prod_work)
         CALL cp_fm_get_diag(prod_work, diag)
         soc = SUM(diag)

         CALL cp_fm_set_element(img_fm, 1, 1+nsc+jsf, soc)

         !the real SOC_y part, <alpha|sigma_y|beta> = -i, <beta|sigma_y|alpha> = i
         CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, vec_soc_y, 0.0_dp, prod_work)
         CALL cp_fm_get_diag(prod_work, diag)
         soc = SUM(diag(1:ndo_mo)) ! alpha-beta
         soc = soc-SUM(diag(ndo_mo+1:ndo_so)) !beta-alpha

         CALL cp_fm_set_element(real_fm, 1, 1+nsc+jsf, soc)

         !Some other stuff we can precompute
         CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, sf_coeffs(jsf)%matrix, vec_soc_z, ncol=ndo_so)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, sf_coeffs(jsf)%matrix, vec_s, ncol=ndo_so)

         ! Do the spin-flip spin-flip loop => block is on diagonal, only cover upper diag
         DO isf = 1, jsf

            !<Psi_Isf|SOC|Psi_Jsf> =
            ! sum_k,sigma [<psi^Isf_k,tau|SOC|psi^Jsf_k,tau> + <psi^Isf_k,tau|psi^Jsf_k,tau> * gs_sum]
            ! - sum_k,l,sigma <psi^0_k,sigma|SOC|psi^0_l,sigma> * <psi^Isf_l,tau|psi^Jsf_k,tau> , tau != sigma

            ! Compute the LR orbitals overlap
            CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, sf_coeffs(isf)%matrix, vec_s, 0.0_dp, prod_work)
            CALL cp_fm_get_submatrix(prod_work, S_ex)

            ! Compute the LR SOC: same-spin integration => only need SOC_z
            CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, sf_coeffs(isf)%matrix, vec_soc_z, 0.0_dp, prod_work)
            CALL cp_fm_get_submatrix(prod_work, SOC_ex)

            ! Actually compute the SOC, with <alpha|sigma_z|alpha> = 1, <beta|sigma_z|beta> = -1
            ! Careful: the alpha-BETA LR are in the alpha spot and the beta-ALPHA in the beta spot
            diag(:) = get_diag(SOC_ex)
            soc = -SUM(diag(1:ndo_mo))+SUM(diag(ndo_mo+1:ndo_so))

            diag(:) = get_diag(S_ex)
            soc = soc+SUM(diag)*gs_sum

            soc = soc-SUM(S_ex(1:ndo_mo, 1:ndo_mo)*TRANSPOSE(domo_soc_z(1:ndo_mo, 1:ndo_mo))) !beta LR, alpha gs
            soc = soc+SUM(S_ex(ndo_mo+1:ndo_so, ndo_mo+1:ndo_so)* & !alpha LR, beta gs
                          TRANSPOSE(domo_soc_z(ndo_mo+1:ndo_so, ndo_mo+1:ndo_so)))

            !purely imaginary element
            CALL cp_fm_set_element(img_fm, 1+nsc+isf, 1+nsc+jsf, soc)

         END DO !isf

         ! Do the spin-conserving spin-flip SOC. Block is off-diagonal => cover all combinations
         DO isc = 1, nsc

            !<Psi_Isc|SOC|Psi_Jsf> =   sum_k,sigma <psi^Isc_k,sigma|SOC|psi^Isf_k,tau>
            !                        - sum_k,l,sigma <psi^0_k,tau|SOC|psi^0_l,sigma

            !Compute the LR orbitals overlap
            CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, sc_coeffs(isc)%matrix, vec_s, 0.0_dp, prod_work)
            CALL cp_fm_get_submatrix(prod_work, S_ex)

            !Start with the fully imaginary SOC_x part => <alpha|sigma_x|beta> = <beta|sigma_x|alpha> = 1
            CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, sc_coeffs(isc)%matrix, vec_soc_x, 0.0_dp, prod_work)
            CALL cp_fm_get_submatrix(prod_work, SOC_ex)

            !TODO: since inly use diag of soc_ex, no need to get the submatrix
            diag(:) = get_diag(SOC_ex)
            soc = SUM(diag)

            soc = soc-SUM(S_ex(1:ndo_mo, ndo_mo+1:ndo_so)*TRANSPOSE(domo_soc_x(1:ndo_mo, ndo_mo+1:ndo_so)))
            soc = soc-SUM(S_ex(ndo_mo+1:ndo_so, 1:ndo_mo)*TRANSPOSE(domo_soc_x(ndo_mo+1:ndo_so, 1:ndo_mo)))

            CALL cp_fm_set_element(img_fm, 1+isc, 1+nsc+jsf, soc)

            !Now do the fully real SOC_y part => <alpha|sigma_y|beta> = -i, <beta|sigma_y|alpha> = i
            CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, sc_coeffs(isc)%matrix, vec_soc_y, 0.0_dp, prod_work)
            CALL cp_fm_get_submatrix(prod_work, SOC_ex)

            diag(:) = get_diag(SOC_ex)
            soc = SUM(diag(1:ndo_mo))-SUM(diag(ndo_mo+1:ndo_so))

            soc = soc+SUM(S_ex(1:ndo_mo, ndo_mo+1:ndo_so)*TRANSPOSE(domo_soc_y(1:ndo_mo, ndo_mo+1:ndo_so)))
            soc = soc-SUM(S_ex(ndo_mo+1:ndo_so, 1:ndo_mo)*TRANSPOSE(domo_soc_y(ndo_mo+1:ndo_so, 1:ndo_mo)))

            CALL cp_fm_set_element(real_fm, 1+isc, 1+nsc+jsf, soc)

         END DO !isc

      END DO !jsf

!  Intermediate clean-up
      CALL cp_fm_release(vec_s)
      CALL cp_fm_release(vec_soc_x)
      CALL cp_fm_release(vec_soc_y)
      CALL cp_fm_release(vec_soc_z)
      CALL cp_fm_release(prod_work)
      IF (do_roks) CALL cp_fm_release(gs_coeffs)
      CALL cp_fm_struct_release(prod_struct)
      CALL cp_fm_struct_release(vec_struct)

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
      ALLOCATE (donor_state%soc_evals(ntot-1))
      donor_state%soc_evals(:) = tmp_evals(2:ntot)-tmp_evals(1)

!  The SOC oscillator strengths
      ALLOCATE(donor_state%soc_osc_str(ntot-1))
      donor_state%soc_osc_str = 0.0_dp

      ! work matrices
      CALL cp_fm_struct_create(dip_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=1)
      CALL cp_cfm_create(dip_cfm, dip_struct)
      CALL cp_cfm_create(work_cfm, full_struct)
      ALLOCATE(transdip(ntot,1))

      ! get the dipole in the AMEW basis
      CALL get_os_amew_op(amew_dip, xas_tdp_env%dipmat, gs_coeffs, donor_state, qs_env)

      DO i = 1, SIZE(amew_dip)

         ! the dipoles are real, but need a cfm for calculations
         CALL cp_fm_to_cfm(msourcer=amew_dip(i)%matrix, mtarget=work_cfm)

         ! store amew_coeffs^c*amew_dip in the pert_cfm
         CALL cp_cfm_gemm('C', 'N', ntot, ntot, ntot, (1.0_dp, 0.0_dp), evecs_cfm, work_cfm, &
                          (0.0_dp, 0.0_dp), pert_cfm)

         ! only interested in transition from the ground-state => take first column of amew_coeffs
         CALL cp_cfm_gemm('N', 'N', ntot, 1, ntot, (1.0_dp, 0.0_dp), pert_cfm, evecs_cfm, &
                          (0.0_dp, 0.0_dp), dip_cfm)

         CALL cp_cfm_get_submatrix(dip_cfm, transdip)

         ! transition dipoles are real numbers
         donor_state%soc_osc_str(:) = donor_state%soc_osc_str(:) + REAL(transdip(2:ntot,1))**2 &
                                      + AIMAG(transdip(2:ntot,1))**2

      END DO !i

      !which representation ?
      IF (xas_tdp_control%dipole_form == xas_dip_len) THEN 
         donor_state%soc_osc_str(:) = 2.0_dp/3.0_dp*donor_state%soc_evals(:)*donor_state%soc_osc_str(:)
      ELSE
         donor_state%soc_osc_str(:) = 2.0_dp/3.0_dp/donor_state%soc_evals(:)*donor_state%soc_osc_str(:)
      END IF

! Clean-up
      DO i = 1,SIZE(amew_dip)
         CALL cp_fm_release(amew_dip(i)%matrix)
      END DO
      DEALLOCATE(amew_dip)
      CALL cp_cfm_release(pert_cfm)
      CALL cp_cfm_release(evecs_cfm)
      CALL cp_cfm_release(work_cfm)
      CALL cp_cfm_release(dip_cfm) 
      CALL cp_fm_struct_release(full_struct)
      CALL cp_fm_struct_release(dip_struct)

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

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "include_rcs_soc", routineP = moduleN//":"//routineN

      INTEGER                                         :: itp, ntp, isg, nsg, ntot, ndo_mo, nao, &
                                                         jtp, i, handle
      REAL(dp), DIMENSION(:), POINTER                 :: sg_evals, tp_evals
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: sg_coeffs, tp_coeffs, amew_dip
      TYPE(cp_fm_type), POINTER                       :: real_fm, img_fm, gs_coeffs, vec_s, &
                                                         prod_fm, vec_soc_x, vec_soc_y, vec_soc_z
      TYPE(cp_fm_struct_type), POINTER                :: full_struct, prod_struct, vec_struct, &
                                                         dip_struct
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(dbcsr_type), POINTER                       :: orb_soc_x, orb_soc_y, orb_soc_z
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: diag, tmp_evals
      REAL(dp)                                        :: soc_st, soc_tt, sqrt2, trace, soc_gst
      TYPE(cp_cfm_type), POINTER                      :: hami_cfm, evecs_cfm, dip_cfm, work_cfm
      COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE        :: transdip
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: domo_soc_x, domo_soc_y, domo_soc_z, s_st, &
                                                         s_tt

      NULLIFY(sg_coeffs, tp_coeffs, gs_coeffs, sg_evals, tp_evals, real_fm, img_fm, full_struct)
      NULLIFY(para_env, blacs_env, vec_s, prod_fm, prod_struct, vec_struct, work_cfm)
      NULLIFY(matrix_s, orb_soc_x, hami_cfm, evecs_cfm, vec_soc_x, vec_soc_y, vec_soc_z)
      NULLIFY(orb_soc_y, orb_soc_z, amew_dip)

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
      ntot = 1+nsg+3*ntp
      ndo_mo = donor_state%ndo_mo
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      orb_soc_x => xas_tdp_env%orb_soc(1)%matrix
      orb_soc_y => xas_tdp_env%orb_soc(2)%matrix
      orb_soc_z => xas_tdp_env%orb_soc(3)%matrix

!  Creating the real part and imaginary part of the final SOC fm
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      CALL cp_fm_create(real_fm, full_struct)
      CALL cp_fm_create(img_fm, full_struct)

!  Put the excitation energies on the diagonal of the real matrix
      DO isg = 1, nsg
         CALL cp_fm_set_element(real_fm, 1+isg, 1+isg, sg_evals(isg))
      END DO
      DO itp = 1, ntp
         ! first T^-1, then T^0, then T^+1
         CALL cp_fm_set_element(real_fm, 1+itp+nsg, 1+itp+nsg, tp_evals(itp))
         CALL cp_fm_set_element(real_fm, 1+itp+ntp+nsg, 1+itp+ntp+nsg, tp_evals(itp))
         CALL cp_fm_set_element(real_fm, 1+itp+2*ntp+nsg, 1+itp+2*ntp+nsg, tp_evals(itp))
      END DO

!  Create the work and helper fms
      CALL cp_fm_get_info(gs_coeffs, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_mo, ncol_global=ndo_mo)
      CALL cp_fm_create(prod_fm, prod_struct)
      CALL cp_fm_create(vec_s, vec_struct)
      CALL cp_fm_create(vec_soc_x, vec_struct)
      CALL cp_fm_create(vec_soc_y, vec_struct)
      CALL cp_fm_create(vec_soc_z, vec_struct)
      ALLOCATE(diag(ndo_mo))
      ALLOCATE(s_st(ndo_mo, ndo_mo), s_tt(ndo_mo, ndo_mo))

!  Precompute everything we can before looping over excited states
      sqrt2 = SQRT(2.0_dp)

   ! The subset of the donor MOs matrix elements: <phi_I^0|Hsoc|phi_J^0>
      ALLOCATE(domo_soc_x(ndo_mo, ndo_mo), domo_soc_y(ndo_mo, ndo_mo), domo_soc_z(ndo_mo, ndo_mo))

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_x, gs_coeffs, vec_soc_x, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_x, 0.0_dp, prod_fm)
      CALL cp_fm_get_submatrix(prod_fm, domo_soc_x)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, gs_coeffs, vec_soc_y, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_y, 0.0_dp, prod_fm)
      CALL cp_fm_get_submatrix(prod_fm, domo_soc_y)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, gs_coeffs, vec_soc_z, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_z, 0.0_dp, prod_fm)
      CALL cp_fm_get_submatrix(prod_fm, domo_soc_z)
      
!  Only have SOC between singlet-triplet triplet-tripelt and ground_state-triplet, the resulting 
!  matrix is Hermitian i.e. the real part is symmetric and the imaginary part is anti-symmetric. 
!  Can only fill upper half

!  Loop over triplet states
      DO jtp = 1, ntp
         ! Compute ground-state-triplet Ms=+-1 SOC
         ! imaginary part
         CALL cp_dbcsr_sm_fm_multiply(orb_soc_x, tp_coeffs(jtp)%matrix, vec_soc_x, ncol=ndo_mo)
         CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_x, 0.0_dp, prod_fm)
         CALL cp_fm_get_diag(prod_fm, diag)
         soc_gst = SUM(diag)
         CALL cp_fm_set_element(img_fm, 1, 1+nsg+jtp, -1.0_dp*soc_gst) ! <0|H_x|T^-1>
         CALL cp_fm_set_element(img_fm, 1, 1+nsg+2*ntp+jtp, soc_gst) ! <0|H_x|T^+1>
         
         ! real part
         CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, tp_coeffs(jtp)%matrix, vec_soc_y, ncol=ndo_mo)
         CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_y, 0.0_dp, prod_fm)
         CALL cp_fm_get_diag(prod_fm, diag)
         soc_gst = SUM(diag)
         CALL cp_fm_set_element(real_fm, 1, 1+nsg+jtp, -1.0_dp*soc_gst) ! <0|H_y|T^-1>
         CALL cp_fm_set_element(real_fm, 1, 1+nsg+2*ntp+jtp, -1.0_dp*soc_gst) ! <0|H_y|T^+1>

         ! Compute the groud_state-triplet SOC: 2* sum_I c^0_I * orb_soc_z * c^T_I
         CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, tp_coeffs(jtp)%matrix, vec_soc_z, ncol=ndo_mo)
         CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_z, 0.0_dp, prod_fm)
         CALL cp_fm_get_diag(prod_fm, diag)
         soc_gst = sqrt2*SUM(diag)
         CALL cp_fm_set_element(img_fm, 1, 1+nsg+ntp+jtp, soc_gst)

         ! Note that the products orb_soc_i*tp_coeffs(jtp) are stored in the vec_soc_i matrices
         ! Also store the matrix_s*tp_coeffs(jtp) in vec_s
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, tp_coeffs(jtp)%matrix, vec_s, ncol=ndo_mo)

         ! Loop over the singlets for singlet-triplet SOC
         DO isg = 1, nsg

            ! compute the singlet-triplet overlap
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, sg_coeffs(isg)%matrix, vec_s, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_submatrix(prod_fm, s_st)

            ! compute the singlet-triplet (Ms=+-1) SOC
            ! imaginary part
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, sg_coeffs(isg)%matrix, vec_soc_x, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)

            soc_st = 0.5_dp*sqrt2*(SUM(diag) - SUM(s_st*TRANSPOSE(domo_soc_x)))

            CALL cp_fm_set_element(img_fm, 1+isg, 1+nsg+jtp, -1.0_dp*soc_st) ! <S|H_x|T^-1>
            CALL cp_fm_set_element(img_fm, 1+isg, 1+nsg+2*ntp+jtp, soc_st) !<S|H_x|T^+1> 

            ! real part
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, sg_coeffs(isg)%matrix, vec_soc_y, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)

            soc_st = 0.5_dp*sqrt2*(SUM(diag) - SUM(s_st*TRANSPOSE(domo_soc_y)))

            CALL cp_fm_set_element(real_fm, 1+isg, 1+nsg+jtp, -1.0_dp*soc_st) ! <S|H_y|T^-1> 
            CALL cp_fm_set_element(real_fm, 1+isg, 1+nsg+2*ntp+jtp, -1.0_dp*soc_st) ! <S|H_y|T^+1> 

            ! compute the singlet-triplet (Ms=0) SOC (purely imaginary)
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, sg_coeffs(isg)%matrix, vec_soc_z, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)

            soc_st = SUM(diag) - SUM(s_st*TRANSPOSE(domo_soc_z))

            CALL cp_fm_set_element(img_fm, 1+isg, 1+nsg+ntp+jtp, soc_st) ! <S|H_x|T^0>

         END DO !isg

         ! Loop over the triplet to do the triplet-triplet SOC
         DO itp = 1, ntp

            ! compute the triplet-triplet overlap
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, tp_coeffs(itp)%matrix, vec_s, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_submatrix(prod_fm, s_tt)

            ! The Ms=0 to Ms=+-1 triplet-triplet SOC
            ! imaginary part
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, tp_coeffs(itp)%matrix, vec_soc_x, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)
            trace = SUM(s_tt*TRANSPOSE(domo_soc_x))

            soc_tt = -0.5_dp*sqrt2*(SUM(diag) + trace)

            CALL cp_fm_set_element(img_fm, 1+nsg+jtp, 1+nsg+ntp+itp, -1.0_dp*soc_tt) ! <T^-1|H_x|T^0>
            CALL cp_fm_set_element(img_fm, 1+nsg+ntp+itp, 1+nsg+2*ntp+jtp, soc_tt) ! <T^0|H_x|T^+1>

            ! real part
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, tp_coeffs(itp)%matrix, vec_soc_y, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)
            trace = SUM(s_tt*TRANSPOSE(domo_soc_y))

            soc_tt = -0.5_dp*sqrt2*(SUM(diag) + trace )

            CALL cp_fm_set_element(real_fm, 1+nsg+jtp, 1+nsg+ntp+itp, soc_tt) ! <T^-1|H_x|T^0>
            CALL cp_fm_set_element(real_fm, 1+nsg+ntp+itp, 1+nsg+2*ntp+jtp, -1.0_dp*soc_tt) ! <T^0|H_x|T^+1>

            ! Do the Ms=+1-Ms=+1 and Ms=-1-Ms=-1 coupling
            IF (itp >= jtp) CYCLE ! SOC purely imaginary, must be zero on the diagonal

            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, tp_coeffs(itp)%matrix, vec_soc_z, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)
            trace = SUM(s_tt*TRANSPOSE(domo_soc_z))

            soc_tt = SUM(diag) + trace 

            CALL cp_fm_set_element(img_fm, 1+nsg+itp, 1+nsg+jtp, -1.0_dp*soc_tt) ! <T^-1|H_x|T^-1>
            CALL cp_fm_set_element(img_fm, 1+nsg+2*ntp+itp, 1+nsg+2*ntp+jtp, soc_tt) !<T^+1|H_x|T^+1>

         END DO !itp
      END DO !jtp

!  Intermediate clean-up
      CALL cp_fm_struct_release(prod_struct)
      CALL cp_fm_release(prod_fm)
      CALL cp_fm_release(vec_s)
      CALL cp_fm_release(vec_soc_x)
      CALL cp_fm_release(vec_soc_y)
      CALL cp_fm_release(vec_soc_z)
      DEALLOCATE (diag)

!  Set-up the complex hermitian matrix
      CALL cp_cfm_create(hami_cfm, full_struct)
      CALL cp_fm_to_cfm(real_fm, img_fm, hami_cfm)

!  Diagonalize the Hamiltonian
      ALLOCATE (tmp_evals(ntot))
      CALL cp_cfm_create(evecs_cfm, full_struct)
      CALL cp_cfm_heevd(hami_cfm, evecs_cfm, tmp_evals)

      !  Adjust the energies so the GS has zero, and store in the donor_state (without the GS)
      ALLOCATE (donor_state%soc_evals(ntot-1))
      donor_state%soc_evals(:) = tmp_evals(2:ntot)-tmp_evals(1)

!  Compute the dipole oscillator strengths
      ALLOCATE (donor_state%soc_osc_str(ntot-1))
      donor_state%soc_osc_str = 0.0_dp

      ! Get some work array/matrices
      CALL cp_fm_struct_create(dip_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=1)
      CALL cp_cfm_create(dip_cfm, dip_struct)
      CALL cp_cfm_create(work_cfm, full_struct)
      ALLOCATE (transdip(ntot, 1))

   ! Get the dipole in the basis of the AMEWs
      CALL get_rcs_amew_op(amew_dip, xas_tdp_env%dipmat, donor_state, qs_env)

      DO i = 1, SIZE(amew_dip)

         ! the dipoles are real, but need a cfm for calculations
         CALL cp_fm_to_cfm(msourcer=amew_dip(i)%matrix, mtarget=work_cfm)

         ! store amew_coeffs^c*amew_dip in the hami_cfm
         CALL cp_cfm_gemm('C', 'N', ntot, ntot, ntot, (1.0_dp, 0.0_dp), evecs_cfm, work_cfm, &
                          (0.0_dp, 0.0_dp), hami_cfm)

         ! only interested in transition from the ground-state => take first column of amew_coeffs
         CALL cp_cfm_gemm('N', 'N', ntot, 1, ntot, (1.0_dp, 0.0_dp), hami_cfm, evecs_cfm, &
                          (0.0_dp, 0.0_dp), dip_cfm)

         CALL cp_cfm_get_submatrix(dip_cfm, transdip)

         ! transition dipoles are real numbers
         donor_state%soc_osc_str(:) = donor_state%soc_osc_str(:)+REAL(transdip(2:ntot, 1))**2 &
                                      +AIMAG(transdip(2:ntot, 1))**2

      END DO

      !which rep ?
      IF (xas_tdp_control%dipole_form == xas_dip_len) THEN
         donor_state%soc_osc_str(:) = 2.0_dp/3.0_dp*donor_state%soc_evals(:)*donor_state%soc_osc_str(:)
      ELSE
         donor_state%soc_osc_str(:) = 2.0_dp/3.0_dp/donor_state%soc_evals(:)*donor_state%soc_osc_str(:)
      END IF

!  Clean-up
      DO i = 1, SIZE(amew_dip)
         CALL cp_fm_release(amew_dip(i)%matrix)
      END DO
      DEALLOCATE (amew_dip)
      CALL cp_fm_struct_release(full_struct)
      CALL cp_fm_struct_release(dip_struct)
      CALL cp_fm_release(real_fm)
      CALL cp_fm_release(img_fm)
      CALL cp_cfm_release(hami_cfm)
      CALL cp_cfm_release(evecs_cfm)
      CALL cp_cfm_release(work_cfm)
      CALL cp_cfm_release(dip_cfm)

      CALL timestop(handle)

   END SUBROUTINE include_rcs_soc

! **************************************************************************************************
!> \brief Computes the matrix elements of a one-body operator (given wrt AOs) in the basis of the
!>        excited state AMEWs with groud state, for the open-shell case
!> \param amew_op the operator in the basis of the AMEWs (array because could have x,y,z components)
!> \param ao_op the operator in the basis of the atomic orbitals
!> \param gs_coeffs the coefficient of the GS donor MOs. Ecplicitely passed because of special 
!>                  format in the ROKS case (see include_os_soc routine)
!> \donor_state ...
!> \qs_env ...
!> \note The ordering of the AMEWs is consistent with SOC and is gs, sc, sf
!>       We assume that the operator is spin-independent => only <0|0>, <0|sc>, <sc|sc> and <sf|sf> 
!>       yield non-zero matrix elements
!>       Only for open-shell calculations
! **************************************************************************************************
   SUBROUTINE get_os_amew_op(amew_op, ao_op, gs_coeffs, donor_state, qs_env)

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: amew_op
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: ao_op
      TYPE(cp_fm_type), POINTER                       :: gs_coeffs
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "get_os_amew_op", routineP = moduleN//":"//routineN

      INTEGER                                         :: i, nsc, nsf, ntot, ndo_so, nao, dim_op, &
                                                         homo, isc, jsc, isf, jsf, ndo_mo
      REAL(dp)                                        :: op
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: sc_coeffs, sf_coeffs
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(dbcsr_type), POINTER                       :: ao_op_i
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                :: full_struct, vec_struct, prod_struct
      TYPE(cp_fm_type), POINTER                       :: vec_work, prod_work, mo_coeff, amew_op_i, &
                                                         vec_s, work
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: domo_op, S_ex
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: gsgs_op, diag

      NULLIFY(matrix_s, para_env, blacs_env, full_struct, vec_struct, prod_struct, mos, vec_work)
      NULLIFY(prod_work, mo_coeff, ao_op_i, amew_op_i, work)

!  Iinitialization
      dim_op = SIZE(ao_op)
      sc_coeffs => donor_state%sc_coeffs
      sf_coeffs => donor_state%sf_coeffs
      nsc = SIZE(sc_coeffs)
      nsf = SIZE(sf_coeffs)
      ntot = 1+nsc+nsf
      ndo_mo = donor_state%ndo_mo
      ndo_so = 2*donor_state%ndo_mo !open-shell => nspins = 2
      CALL get_qs_env(qs_env, matrix_s=matrix_s, para_env=para_env, blacs_env=blacs_env, mos=mos)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

!  Create the amew_op matrix set
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      ALLOCATE(amew_op(dim_op))
      DO i = 1, dim_op
         CALL cp_fm_create(amew_op(i)%matrix, full_struct)
      END DO

!  Before looping, need to evaluate sum_j,sigma <phi^0_j,sgima|op|phi^j,sigma>, for each dimension
!  of the operator
      ALLOCATE(gsgs_op(dim_op))

      !start with the alpha MOs
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, homo=homo)
      ALLOCATE(diag(homo))
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
      DEALLOCATE(diag)
      NULLIFY(vec_struct)

      !then beta orbitals
      CALL get_mo_set(mos(2)%mo_set, mo_coeff=mo_coeff, homo=homo)
      ALLOCATE(diag(homo))
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
      DEALLOCATE(diag)
      NULLIFY(vec_struct)

!  Before looping over excited AMEWs, define some work matrices and structures
      CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nao, ncol_global=ndo_so)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_so, ncol_global=ndo_so)
      CALL cp_fm_create(vec_work, vec_struct) !for op*|phi>
      CALL cp_fm_create(prod_work, prod_struct) !for any <phi|op|phi>
      CALL cp_fm_create(vec_s, vec_struct) !for S*|phi>
      CALL cp_fm_create(work, full_struct)
      ALLOCATE(diag(ndo_so))
      ALLOCATE(domo_op(ndo_so, ndo_so), S_ex(ndo_so, ndo_so))

!  Loop over the dimensions of the operator
      DO  i = 1, dim_op

         ao_op_i => ao_op(i)%matrix
         amew_op_i => amew_op(i)%matrix

         !put the gs-gs contribution
         CALL cp_fm_set_element(amew_op_i, 1, 1, gsgs_op(i))

      !  Precompute what we can before looping over excited states
         ! Need the operator in the donor MOs basis <phi^0_I,sigma|op_i|phi^0_J,tau>
         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, gs_coeffs, vec_work, ncol=ndo_so)
         CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, vec_work, 0.0_dp, prod_work)
         CALL cp_fm_get_submatrix(prod_work, domo_op)


      ! Loop over the spin-conserving excitations
         DO jsc = 1, nsc

            !compute everything ket we can with jsc
            ! ao_op|Psi_Jsc>
            CALL cp_dbcsr_sm_fm_multiply(ao_op_i, sc_coeffs(jsc)%matrix, vec_work, ncol=ndo_so)

            ! S|Psi_Jsc>
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, sc_coeffs(jsc)%matrix, vec_s, ncol=ndo_so)

            !compute the GS-SC contribution: <0|op|Psi_Jsc> = sum_i,sigma <phi^0_i,sigma|op|phi^Jsc_i,sigma>
            CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, gs_coeffs, vec_work, 0.0_dp, prod_work)
            CALL cp_fm_get_diag(prod_work, diag)
            op = SUM(diag)

            CALL cp_fm_set_element(amew_op_i, 1, 1+jsc, op)

            DO isc = 1, jsc

               !<Psi_Isc|op|Psi_Jsc> = 
               ! sum_k,sigma [<psi^Isc_k,sigma|op|psi^Jsc_k,sigma> + <psi^Isc_k,sigma|psi^Jsc_k,sigma> * gsgs_op] 
               ! - sum_k,l,sigma <psi^0_k,sigma|op|psi^0_l,sigma> * <psi^Isc_l,sigma|psi^Jsc_k,sigma>

               ! LR orbital overlap
               CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, sc_coeffs(isc)%matrix, vec_s, 0.0_dp, prod_work)
               CALL cp_fm_get_submatrix(prod_work, S_ex)

               ! LR orbital op
               CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, sc_coeffs(isc)%matrix, vec_work, 0.0_dp, prod_work)
               CALL cp_fm_get_diag(prod_work, diag)

               !the matrix element itself
               op = SUM(diag)

               diag(:) = get_diag(S_ex)
               op = op + SUM(diag)*gsgs_op(i)

               !alpha
               op = op - SUM(S_ex(1:ndo_mo,1:ndo_mo)*TRANSPOSE(domo_op(1:ndo_mo,1:ndo_mo)))
               !beta
               op = op - SUM(S_ex(ndo_mo+1:ndo_so,ndo_mo+1:ndo_so)*TRANSPOSE(domo_op(ndo_mo+1:ndo_so,ndo_mo+1:ndo_so)))

               CALL cp_fm_set_element(amew_op_i, 1+isc, 1+jsc, op)

            END DO !isc
         END DO !jsc

      ! Loop over the spin-flip excitations
         DO jsf = 1, nsf

            !compute everything ket we can with jsf
            ! ao_op|Psi_Jsf>
            CALL cp_dbcsr_sm_fm_multiply(ao_op_i, sf_coeffs(jsf)%matrix, vec_work, ncol=ndo_so)

            ! S|Psi_Jsf>
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, sf_coeffs(jsf)%matrix, vec_s, ncol=ndo_so)

            DO isf = 1, jsf

               !<Psi_Isf|op|Psi_Jsf> = 
               ! sum_k,sigma [<psi^Isf_k,sigma|op|psi^Jsf_k,sigma> + <psi^Isf_k,sigma|psi^Jsf_k,sigma> * gsgs_op] 
               ! - sum_k,l,sigma <psi^0_k,sigma|op|psi^0_l,sigma> * <psi^Isf_l,tau|psi^Jsf_k,tau>, tau != sigma

               ! LR orbital overlap
               CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, sf_coeffs(isf)%matrix, vec_s, 0.0_dp, prod_work)
               CALL cp_fm_get_submatrix(prod_work, S_ex)

               ! LR orbital op
               CALL cp_gemm('T', 'N', ndo_so, ndo_so, nao, 1.0_dp, sf_coeffs(isf)%matrix, vec_work, 0.0_dp, prod_work)
               CALL cp_fm_get_diag(prod_work, diag)

               !the matrix element itself
               op = SUM(diag)

               diag(:) = get_diag(S_ex)
               op = op + SUM(diag)*gsgs_op(i)

               ! alpha-beta
               op = op - SUM(S_ex(1:ndo_mo,1:ndo_mo)*TRANSPOSE(domo_op(ndo_mo+1:ndo_so,ndo_mo+1:ndo_so)))
               ! beta-alpha
               op = op - SUM(S_ex(ndo_mo+1:ndo_so,ndo_mo+1:ndo_so)*TRANSPOSE(domo_op(1:ndo_mo,1:ndo_mo)))

               CALL cp_fm_set_element(amew_op_i, 1+nsc+isf, 1+nsc+jsf, op)

            END DO !isc
         END DO !jsc

         !Symmetry => only upper diag explicitly built
         CALL cp_fm_upper_to_full(amew_op_i, work)

      END DO !i

!  Clean-up
      CALL cp_fm_struct_release(full_struct)
      CALL cp_fm_struct_release(prod_struct)
      CALL cp_fm_struct_release(vec_struct)
      CALL cp_fm_release(work)
      CALL cp_fm_release(vec_s)
      CALL cp_fm_release(vec_work)
      CALL cp_fm_release(prod_work)

   END SUBROUTINE get_os_amew_op

! **************************************************************************************************
!> \brief Computes the matrix elements of a one-body operator (given wrt AOs) in the basis of the
!>        excited state AMEWs with ground state, singlet and triplet with Ms = -1,0,+1
!> \param amew_op the operator in the basis of the AMEWs (array because could have x,y,z components)
!> \param ao_op the operator in the basis of the atomic orbitals
!> \donor_state ...
!> \qs_env ...
!> \note The ordering of the AMEWs is consistent with SOC and is gs, sg, tp(-1), tp(0). tp(+1)
!>       We assume that the operator is spin-independent => only <0|0>, <0|S>, <S|S> and <T|T>
!>       yield non-zero matrix elements
!>       Only for spin-restricted calculations
! **************************************************************************************************
   SUBROUTINE get_rcs_amew_op(amew_op, ao_op, donor_state, qs_env)

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: amew_op
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: ao_op
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN ="get_rcs_amew_op", routineP = moduleN//":"//routineN

      INTEGER                                         :: i, nsg, ntp, ntot, ndo_mo, nao, homo, isg,&
                                                         jsg, itp, jtp, dim_op
      REAL(dp)                                        :: sqrt2, op
      TYPE(cp_fm_type), POINTER                       :: gs_coeffs, prod_fm, amew_op_i, work_fm, &
                                                         vec_op, mo_coeff, vec_s, gs_fm, &
                                                         work
      TYPE(cp_fm_struct_type), POINTER                :: full_struct, prod_struct, vec_struct, &
                                                         std_struct, gsgs_struct
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: sg_coeffs, tp_coeffs
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: diag, gs_diag
      TYPE(dbcsr_type), POINTER                       :: ao_op_i
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: gsgs_op
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: domo_op, S_ex

      NULLIFY(gs_coeffs, sg_coeffs, tp_coeffs, matrix_s, full_struct, prod_fm, prod_struct, work_fm)
      NULLIFY(vec_struct, blacs_env, para_env, mo_coeff, mos, gsgs_struct, std_struct)
      NULLIFY(vec_op, vec_s, gs_fm)

!  Initialization
      gs_coeffs => donor_state%gs_coeffs
      sg_coeffs => donor_state%sg_coeffs
      tp_coeffs => donor_state%tp_coeffs
      nsg = SIZE(sg_coeffs)
      ntp = SIZE(tp_coeffs)
      ntot = 1+nsg+3*ntp
      ndo_mo = donor_state%ndo_mo
      CALL get_qs_env(qs_env, matrix_s=matrix_s, para_env=para_env, blacs_env=blacs_env, mos=mos)
      sqrt2 = SQRT(2.0_dp)
      dim_op = SIZE(ao_op)

!  Create the amew_op matrix
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      ALLOCATE(amew_op(dim_op))
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
      ALLOCATE(gsgs_op(dim_op))
      ALLOCATE(gs_diag(homo))

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
      DEALLOCATE(gs_diag)

!  Create the work and helper fms
      CALL cp_fm_get_info(gs_coeffs, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_mo, ncol_global=ndo_mo)
      CALL cp_fm_create(prod_fm, prod_struct)
      CALL cp_fm_create(vec_op, vec_struct)
      CALL cp_fm_create(vec_s, vec_struct)
      CALL cp_fm_create(work, full_struct)
      ALLOCATE(diag(ndo_mo))
      ALLOCATE(domo_op(ndo_mo, ndo_mo), S_ex(ndo_mo, ndo_mo))

! Iterate over the dimensions of the operator
      DO i = 1, dim_op

         ao_op_i => ao_op(i)%matrix
         amew_op_i => amew_op(i)%matrix

      ! The GS-GS contribution
         CALL cp_fm_set_element(amew_op_i, 1, 1, gsgs_op(i))
         
      ! Compute the operator for the donor MOs 
         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, gs_coeffs, vec_op, ncol=ndo_mo)    
         CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_op, 0.0_dp, prod_fm)
         CALL cp_fm_get_submatrix(prod_fm, domo_op)

         ! Loop over singlet states
         DO jsg = 1, nsg

            ! compute the ground-state/singlet components. ao_op*sg_coeff(jsg) stored in vec_op
            CALL cp_dbcsr_sm_fm_multiply(ao_op_i, sg_coeffs(jsg)%matrix, vec_op, ncol=ndo_mo)
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_op, 0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)
            op = sqrt2*SUM(diag)
            CALL cp_fm_set_element(amew_op_i, 1, 1+jsg, op)

            ! also store matrix_s*sg_coeff(jsg) in vec_s
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, sg_coeffs(jsg)%matrix, vec_s, ncol=ndo_mo)

            ! compute the singlet-singlet components
            DO isg = 1, jsg

               ! isg*ao_op*jsg
               CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, sg_coeffs(isg)%matrix, vec_op, &
                            0.0_dp, prod_fm)
               CALL cp_fm_get_diag(prod_fm, diag)
               op = SUM(diag)

               ! the singlet-singlet overlap
               CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, sg_coeffs(isg)%matrix, vec_s, &
                            0.0_dp, prod_fm)
               CALL cp_fm_get_submatrix(prod_fm, S_ex)

               diag(:) = get_diag(S_ex)
               op = op + SUM(diag)*gsgs_op(i)
               op = op - SUM(S_ex*TRANSPOSE(domo_op))

               CALL cp_fm_set_element(amew_op_i, 1+isg, 1+jsg, op)

            END DO !isg
         END DO !jsg

         ! compute the triplet-triplet components
         DO jtp = 1, ntp

            ! precompute matrix_s*jtp and ao_op*jtp
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, tp_coeffs(jtp)%matrix, vec_s, ncol=ndo_mo)
            CALL cp_dbcsr_sm_fm_multiply(ao_op_i, tp_coeffs(jtp)%matrix, vec_op, ncol=ndo_mo)

            DO itp = 1, jtp

               ! itp*ao_op*jtp
               CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, tp_coeffs(itp)%matrix, vec_op, &
                            0.0_dp, prod_fm)
               CALL cp_fm_get_diag(prod_fm, diag)
               op = SUM(diag)

               ! triplet-triplet overlap
               CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, tp_coeffs(itp)%matrix, vec_s, &
                            0.0_dp, prod_fm)
               CALL cp_fm_get_submatrix(prod_fm, S_ex)

               diag(:) = get_diag(S_ex)
               op = op + SUM(diag)*gsgs_op(i)
               op = op - SUM(S_ex*TRANSPOSE(domo_op))

               ! Same value for the three triplet pairs
               CALL cp_fm_set_element(amew_op_i, 1+nsg+itp, 1+nsg+jtp, op)
               CALL cp_fm_set_element(amew_op_i, 1+nsg+ntp+itp, 1+nsg+ntp+jtp, op)
               CALL cp_fm_set_element(amew_op_i, 1+nsg+2*ntp+itp, 1+nsg+2*ntp+jtp, op)

            END DO !itp
         END DO !jtp

         ! Symmetrize the matrix (only upper triangle built)
         CALL cp_fm_upper_to_full(amew_op_i, work)

      END DO !i

!  Clean-up
      CALL cp_fm_release(prod_fm)
      CALL cp_fm_release(work)
      CALL cp_fm_release(vec_s)
      CALL cp_fm_release(vec_op)
      CALL cp_fm_struct_release(prod_struct)
      CALL cp_fm_struct_release(full_struct)

   END SUBROUTINE get_rcs_amew_op

END MODULE xas_tdp_utils

