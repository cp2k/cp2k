!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for X-ray absorption spectroscopy using TDDFPT
!> \author AB (01.2018)
! **************************************************************************************************

MODULE xas_tdp_utils
   USE arnoldi_api,                     ONLY: arnoldi_data_type,&
                                              arnoldi_ev,&
                                              deallocate_arnoldi_data,&
                                              get_selected_ritz_val,&
                                              setup_arnoldi_data
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose, cp_dbcsr_cholesky_invert
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm, cp_dbcsr_sm_fm_multiply, &
                                              copy_fm_to_dbcsr
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_fm_to_cfm, cp_cfm_to_fm, cp_cfm_type, &
                                              cp_cfm_create, cp_cfm_release, cp_cfm_get_submatrix
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose, cp_fm_cholesky_invert
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale, cp_fm_trace, &
                                              cp_fm_triangular_invert, cp_fm_triangular_multiply, &
                                              cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver, cp_fm_geeig, cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create, cp_fm_struct_release, &     
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, &
        cp_fm_release, cp_fm_set_element, cp_fm_to_fm, cp_fm_to_fm_submat, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_complete_redistribute, dbcsr_copy, &
        dbcsr_copy_into_existing, dbcsr_create, dbcsr_desymmetrize, dbcsr_distribution_get, &
        dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_reserve_block2d, &
        dbcsr_set, dbcsr_transposed, dbcsr_type, dbcsr_type_symmetric
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE generic_os_integrals,            ONLY: int_operators_r12_ab_os
   USE input_constants,                 ONLY: xas_dip_len
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat_symm
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_max,&
                                              mp_min
   USE molecule_types,                  ONLY: molecule_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE qs_o3c_methods,                  ONLY: calculate_o3c_coulomb_integrals,&
                                              contract3_o3c
   USE qs_o3c_types,                    ONLY: &
        get_o3c_container, get_o3c_iterator_info, get_o3c_vec, init_o3c_container, &
        o3c_container_type, o3c_iterate, o3c_iterator_create, o3c_iterator_release, &
        o3c_iterator_type, o3c_vec_create, o3c_vec_release, o3c_vec_type, release_o3c_container
   USE util,                            ONLY: get_limit
   USE xas_tdp_types,                   ONLY: donor_state_type,&
                                              xas_tdp_control_type,&
                                              xas_tdp_env_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_utils'

   PUBLIC :: setup_xas_tdp_prob, solve_xas_tdp_prob, &
             compute_ri_coulomb2_int, include_full_soc, build_xas_tdp_3c_nl, &
             compute_ri_o3c_coulomb, compute_ri_o3c_exchange

CONTAINS

! **************************************************************************************************
!> \brief Builds the matrix that defines the XAS TDDFPT generalized eigenvalue problem to be solved
!>        for excitation energies omega. The problem has the form omega*G*C = M*C, where C contains
!>        the reponse orbitals coefficients. The matrix M and the metric G are stored in the given
!>        donor_state .
!> \param donor_state the donor_state for which the problem is restricted
!> \param qs_env ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \note the matrix M is symmetric and has the form | M_d   M_o |
!>                                                  | M_o   M_d |,
!>       depending on whther we consider singlet or triplet excitation, the diagonal (M_d) and
!>       off-diagonal (M_o) parts of M differ:
!>       - For singlet: M_d = A + 2B + C_aa + C_ab - D
!>                      M_o = 2B + C_aa + C_ab - E
!>       - For triplet: M_d = A + C_aa - C_ab - D
!>                      M_o = C_aa - C_ab - E
!>       where other subroutines computes the matrices A, B, E, D and G, which are:
!>       - A: the ground-state contribution: F_ab*delta_IJ - epsilon_IJ*S_ab
!>       - B: the Coulob kernel ~(aI|Jb)
!>       - C: the xc kernel c_aa (double derivatibe wrt to n_alpha) and C_ab (wrt n_alpha and n_beta)
!>       - D: the on-digonal exact exchange kernel ~(ab|IJ)
!>       - E: the off-diagonal exact exchange kernel ~(aJ|Ib)
!>       - G: the metric  S_ab*delta_IJ
!>       For the xc functionals, C_aa + C_ab or C_aa - C_ab are stored in the same matrix
!>       In the above definitions, I,J label the donnor MOs and a,b the sgfs of the basis
!>       Under TDA, only the on-diagonal elements of M are computed
!>       In the case of non-TDA, one turns the problem Hermitian
!>       The matrix M is always  projected onto the unperturbed unoccupied  ground state:
!>       M <= Q * M * Q^T = (1 - SP) * M * (1 - PS)
!>       All possibilities: TDA or full-TDDFT, singlet or triplet, xc or hybrid, etc are treated
!>       in the same routine to avoid recomputing stuff
! **************************************************************************************************
   SUBROUTINE setup_xas_tdp_prob(donor_state, qs_env, xas_tdp_env, xas_tdp_control)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_xas_tdp_prob', &
                                     routineP = moduleN//":"//routineN

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: submat_blk_size
      LOGICAL                                            :: do_hfx, do_sg, do_tda, do_tp, do_xc
      REAL(dp)                                           :: sx
      TYPE(dbcsr_distribution_type), POINTER             :: submat_dist
      TYPE(dbcsr_type), POINTER                          :: matrix_a, matrix_b, matrix_c_sg, &
                                                            matrix_c_tp, matrix_d, matrix_e, &
                                                            proj_Q, sg_matrix_tdp, tp_matrix_tdp, &
                                                            work

      NULLIFY (sg_matrix_tdp, tp_matrix_tdp, submat_dist, submat_blk_size, matrix_a, matrix_b)
      NULLIFY (matrix_c_sg, matrix_c_tp, matrix_d, matrix_e, proj_Q, work)

      CALL timeset(routineN, handle)

!  Initialization
      do_sg = xas_tdp_control%do_singlet
      do_tp = xas_tdp_control%do_triplet
      do_xc = xas_tdp_control%do_xc
      do_hfx = xas_tdp_control%do_hfx
      do_tda = xas_tdp_control%tamm_dancoff
      sx = xas_tdp_control%sx
      IF (do_sg) THEN
         ALLOCATE (donor_state%sg_matrix_tdp)
         sg_matrix_tdp => donor_state%sg_matrix_tdp
      END IF
      IF (do_tp) THEN
         ALLOCATE (donor_state%tp_matrix_tdp)
         tp_matrix_tdp => donor_state%tp_matrix_tdp
      END IF

!  Get the dist and block size of all matrices A, B, C, etc
      CALL compute_submat_dist_and_blk_size(submat_dist, submat_blk_size, donor_state, qs_env)

!  Allocate and compute all the matrices A, B, C, etc we will need
      ALLOCATE (matrix_a, matrix_b, matrix_c_sg, matrix_c_tp, matrix_d, matrix_e, proj_Q, work)

      ! The projector on the unoccupied unperturbed ground state 1-SP and associated work matrix
      CALL get_q_projector(proj_Q, submat_dist, submat_blk_size, donor_state, xas_tdp_env)
      CALL dbcsr_create(matrix=work, matrix_type="N", dist=submat_dist, name="WORK", &
                        row_blk_size=submat_blk_size, col_blk_size=submat_blk_size)

      ! The ground state contribution
      CALL build_gs_contribution(matrix_a, submat_dist, submat_blk_size, donor_state, qs_env)

      ! The Coulomb and XC kernels (note: in the rare case of CIS triplets, none is needed)
      IF (do_xc .OR. do_sg) THEN
         CALL kernel_coulomb_xc(matrix_b, matrix_c_sg, matrix_c_tp, do_sg, (do_sg .AND. do_xc), &
                                (do_tp .AND. do_xc), submat_dist, submat_blk_size, donor_state, &
                                xas_tdp_env, qs_env)
      END IF

      ! The exact exchange. Off-diagonal kernel only needed in full-TDDFT
      IF (do_hfx) THEN
         CALL kernel_exchange(matrix_d, matrix_e, (.NOT. do_tda), submat_dist, submat_blk_size, &
                              donor_state, xas_tdp_env, qs_env)
      END IF

      ! Build the metric G, also need its inverse in case of full-TDDFT
      IF (do_tda) THEN
         ALLOCATE (donor_state%metric(1))
         CALL build_metric(donor_state%metric, submat_dist, submat_blk_size, donor_state, qs_env)
      ELSE
         ALLOCATE (donor_state%metric(2))
         CALL build_metric(donor_state%metric, submat_dist, submat_blk_size, donor_state, qs_env, &
                           do_inv=.TRUE.)
      END IF

!  Build the eigenvalue problem, denpending on the case (TDA, singlet, triplet, hfx, etc ...)
      IF (do_tda) THEN

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
         CALL build_aux_matrix(1.0E-8_dp, sx, matrix_a, matrix_d, matrix_e, do_hfx, proj_Q, &
                               work, donor_state, qs_env)   

         IF (do_sg) THEN ! full-TDDFT singlets

            ! The final matrix is the sum of the on- and off-diagonal elements as in the description
            ! M = A + 4B + 2(C_aa + C_ab) - D - E
            CALL dbcsr_copy(sg_matrix_tdp, matrix_a, name="SINGLET MATRIX TDP")
            CALL dbcsr_add(sg_matrix_tdp, matrix_b, 1.0_dp, 4.0_dp)

            IF (do_hfx) THEN !scaled hfx
               CALL dbcsr_add(sg_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx)
               CALL dbcsr_add(sg_matrix_tdp, matrix_e, 1.0_dp, -1.0_dp*sx)
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
               CALL dbcsr_add(tp_matrix_tdp, matrix_e, 1.0_dp, -1.0_dp*sx)
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
      CALL dbcsr_release(matrix_b)
      CALL dbcsr_release(matrix_c_sg)
      CALL dbcsr_release(matrix_c_tp)
      CALL dbcsr_release(matrix_d)
      CALL dbcsr_release(matrix_e)
      CALL dbcsr_release(proj_Q)
      CALL dbcsr_release(work)
      CALL dbcsr_distribution_release(submat_dist)
      DEALLOCATE (matrix_a, matrix_b, matrix_c_sg, matrix_c_tp, matrix_d, matrix_e, proj_Q, work)
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
!> \param singlet if the problem has to be solved for singlet (vs triplet)
!> \note The computed eigenvalues and eigenvectors are stored in the donor_state
!>       The eigenvectors are the LR-coefficients. In case of TDA, c^- is stored. In the general
!>       case, the sum c^+ + c^- is stored.
!>       In case both singlets and triplets are considered, this routine must be called twice. This
!>       is the choice that was made because the body of the routine is exactly the same in both cases
!>       Note that for singlet we solve for u = 1/sqrt(2)*(c_alpha + c_beta) = sqrt(2)*c
!>       and that for triplets we solve for v = 1/sqrt(2)*(c_alpha - c_beta) = sqrt(2)*c
!>       Additional speed-up is reached by projecting out the ground state
!>       eigenvectors (with degenerate 0 eigenvalue).
! **************************************************************************************************
   SUBROUTINE solve_xas_tdp_prob(donor_state, xas_tdp_control, xas_tdp_env, qs_env, singlet)  

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(qs_environment_type), POINTER              :: qs_env
      LOGICAL, INTENT(IN)                             :: singlet

      CHARACTER(len=*), PARAMETER :: routineN = "solve_xas_tdp_prob", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: nrow, i, nao, imo, &
                                                         ndo_mo, nelectron, &
                                                         nevals, nocc, handle
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(dbcsr_type), POINTER                       :: matrix_tdp, tmp_mat
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct, gs_struct
      TYPE(cp_fm_type), POINTER                       :: rhs_matrix, lhs_matrix, work, c_diff, c_sum
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: scaling, tmp_evals
      REAL(dp), DIMENSION(:), POINTER                 :: lr_evals
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: lr_coeffs

      CALL timeset(routineN, handle)

      NULLIFY(para_env, blacs_env, fm_struct, rhs_matrix, matrix_tdp, lhs_matrix, work, tmp_mat)
      NULLIFY(c_diff, c_sum, gs_struct, lr_evals, lr_coeffs)
      CPASSERT(ASSOCIATED(xas_tdp_env))

      IF (singlet) THEN
         matrix_tdp => donor_state%sg_matrix_tdp
      ELSE
         matrix_tdp => donor_state%tp_matrix_tdp
      END IF
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env, nelectron_total=nelectron)

!     Initialization
      CALL dbcsr_get_info(matrix_tdp, nfullrows_total=nrow)
      ndo_mo = SIZE(donor_state%mo_indices)
      nocc = ndo_mo*nelectron/2
      nevals = nrow - nocc
      IF (xas_tdp_control%n_excited > 0 .AND. xas_tdp_control%n_excited < nevals) THEN
         nevals = xas_tdp_control%n_excited
      END IF
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nrow, para_env=para_env, &
                               ncol_global=nrow)

!     Need to go to full matrices to diagonalize
      CALL cp_fm_create(rhs_matrix, fm_struct)
         
!     Allocating space for solutions and oscilaltor strengths
      ALLOCATE(tmp_evals(nrow))
      IF (singlet) THEN
         ALLOCATE (donor_state%sg_evals(nevals))
         lr_evals => donor_state%sg_evals
         ALLOCATE(donor_state%sg_coeffs(nevals))
         lr_coeffs => donor_state%sg_coeffs
         ALLOCATE(donor_state%sg_osc_str(nevals))
         donor_state%sg_osc_str = 0.0_dp
      ELSE 
         ALLOCATE(donor_state%tp_evals(nevals))
         lr_evals => donor_state%tp_evals
         ALLOCATE(donor_state%tp_coeffs(nevals))
         lr_coeffs => donor_state%tp_coeffs
         ALLOCATE(donor_state%tp_osc_str(nevals))
         donor_state%tp_osc_str = 0.0_dp
      END IF
      CALL cp_fm_create(c_sum, fm_struct)
      CALL cp_fm_create(work, fm_struct)

!     Test on TDA
      IF (xas_tdp_control%tamm_dancoff) THEN

!        Shift away zero eigenvalues
!   TODO: since we have the LUMOS as guess, we also have guess energies and one could shift
!         the spectrum by that amount (a bit more actually) => more gentle approach
!         Maybe it is also possible that this is not needed at all since the preconditioner
!         uses Q and Q^T, so the occupied MOs are out of the way (but maybe not)
!         Finally, even if it is not strictly needed, might be usefull to have the same treatment
!         of the result, whatever the choice of solver is. Also it might make the matrix dioganlly
!         dominant, which would explain why the preconditioner works well
         !CALL get_gs_filter(filter, matrix_tdp, donor_state, xas_tdp_env)
         !CALL dbcsr_add(matrix_tdp, filter, 1.0_dp, max_ev+10.0_dp)

!        Get the main matrix_tdp as an fm
         CALL copy_dbcsr_to_fm(matrix_tdp, rhs_matrix)

!        Get the metric as a fm
         CALL cp_fm_create(lhs_matrix, fm_struct)
         CALL copy_dbcsr_to_fm(donor_state%metric(1)%matrix, lhs_matrix)

!        Diagonalisation (Cholesky decomposition). In TDA, c_sum = c^-
         CALL cp_fm_geeig(rhs_matrix, lhs_matrix, c_sum, tmp_evals, work)

!        Take the first nevals eigenvalues that are not projected to zero
         lr_evals(:) = tmp_evals(nocc+1:nocc+nevals)

!        TDA specific clean-up
         CALL cp_fm_release(lhs_matrix)

      ELSE ! not TDA

!        Need to multiply the current matrix_tdp with the auxiliary matrix
!        tmp_mat =  (A-D+E)^0.5 * M * (A-D+E)^0.5
         ALLOCATE(tmp_mat)
         CALL dbcsr_create(matrix=tmp_mat, template=matrix_tdp, matrix_type="N")
         CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%matrix_aux, matrix_tdp, 0.0_dp, tmp_mat)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, tmp_mat, donor_state%matrix_aux, 0.0_dp, tmp_mat)

!        Get the matrix as a fm
         CALL copy_dbcsr_to_fm(tmp_mat, rhs_matrix)

!        Solve the "turned-Hermitian" eigenvalue problem
         CALL choose_eigv_solver(rhs_matrix, work, tmp_evals)

!        Currently, work = (A-D+E)^0.5 (c^+ - c^-) and tmp_evals = omega^2
!        Put tiny almost zero eigenvalues to zero (corresponding to occupied MOs)
         WHERE (tmp_evals < 1.0E-4_dp) tmp_evals = 0.0_dp

!        Retrive c_diff = (c^+ - c^-) for normalization
!        (c^+ - c^-) = 1/omega^2 * M * (A-D+E)^0.5 * work
         CALL cp_fm_create(c_diff, fm_struct)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_tdp, donor_state%matrix_aux, 0.0_dp, tmp_mat)
         CALL cp_dbcsr_sm_fm_multiply(tmp_mat, work, c_diff, ncol=nrow)

         ALLOCATE(scaling(nrow))
         WHERE (tmp_evals .NE. 0.0_dp) scaling = 1.0_dp/tmp_evals
         CALL cp_fm_column_scale(c_diff, scaling)

!        Normalize with the metric: c_diff * G * c_diff = +- 1
         scaling = 0.0_dp
         CALL get_normal_scaling(scaling, c_diff, donor_state)
         CALL cp_fm_column_scale(c_diff, scaling)

!        Get the actual eigenvalues
         tmp_evals = SQRT(tmp_evals)
         lr_evals(:) = tmp_evals(nocc+1:nocc+nevals)

!        Get c_sum = (c^+ + c^-), which appears in all transition density related expressions
!        c_sum = -1/omega G^-1 * (A-D+E) * (c^+ - c^-)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%matrix_aux, donor_state%matrix_aux, 0.0_dp, tmp_mat)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, tmp_mat, 0.0_dp, tmp_mat)
         CALL cp_dbcsr_sm_fm_multiply(tmp_mat, c_diff, c_sum, ncol=nrow)
         WHERE (tmp_evals .NE. 0) scaling = -1.0_dp/tmp_evals
         CALL cp_fm_column_scale(c_sum, scaling)

!        Full TDDFT specific clean-up
         CALL cp_fm_release(c_diff)
         CALL dbcsr_release(tmp_mat)
         DEALLOCATE(tmp_mat, scaling)

      END IF ! TDA

!  Reorganize the eigenvectors in array of cp_fm_p_types so that each element corresponds to an
!  excited state. Makes later calls to those easier and more efficient
      CALL cp_fm_get_info(donor_state%gs_coeffs, matrix_struct=gs_struct, nrow_global=nao)      

      DO i = 1, nevals

         CALL cp_fm_create(lr_coeffs(i)%matrix, gs_struct)

         DO imo = 1, ndo_mo

            CALL cp_fm_to_fm_submat(msource=c_sum, mtarget=lr_coeffs(i)%matrix, &
                                    nrow=nao, ncol=1, s_firstrow=(imo-1)*nao+1, s_firstcol=nocc+i, &
                                    t_firstrow=1, t_firstcol=imo)

         END DO !imo
      END DO !istate

!  Clean-up
      CALL cp_fm_release(rhs_matrix)
      CALL cp_fm_release(work)
      CALL cp_fm_release(c_sum)
      CALL cp_fm_struct_release(fm_struct)

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
!>        the properties, which are based on the replication of the KS matrix
!> \param submat_dist the dbcsr distribution of the submatrices
!> \param submat_blk_size the row/col block sizes for the submatrices
!> \param donor_state ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_submat_dist_and_blk_size(submat_dist, submat_blk_size, donor_state, qs_env)

      TYPE(dbcsr_distribution_type), POINTER             :: submat_dist
      INTEGER, DIMENSION(:), POINTER                     :: submat_blk_size
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_submat_dist_and_blk_size', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: group, i, nao, nblk_row, ndo_mo, &
                                                            scol_dist, srow_dist
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, col_dist_sub, row_blk_size, &
                                                            row_dist, row_dist_sub
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks

      NULLIFY (matrix_ks, dbcsr_dist, row_blk_size, row_dist, col_dist, pgrid, col_dist_sub)
      NULLIFY (row_dist_sub)

!  The submatrices are indexed by M_{pi,qj}, where p,q label basis functions and i,j donor MOs
!  We consider them as being made of blocks of size nao x nao and i,j label those blocks
!  Each of these blocks has the distribution and the block sizes of the KS matrix, that we simply
!  replicate ndo_mo*ndo_mo (# donor MOs) times

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices, 1)
      CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_get_info(matrix_ks(1)%matrix, row_blk_size=row_blk_size)
      CALL dbcsr_distribution_get(dbcsr_dist, row_dist=row_dist, col_dist=col_dist, group=group, &
                                  pgrid=pgrid)
      nao = SUM(row_blk_size)
      nblk_row = SIZE(row_blk_size)
      srow_dist = SIZE(row_dist)
      scol_dist = SIZE(col_dist)

!  Creation if submatrix block size and col/row distribution
      ALLOCATE (submat_blk_size(ndo_mo*nblk_row))
      ALLOCATE (row_dist_sub(ndo_mo*srow_dist))
      ALLOCATE (col_dist_sub(ndo_mo*scol_dist))

      DO i = 1, ndo_mo
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
!> \brief Get the projector on the unperturbed occupied ground-state SPS on the MO blok diagonal
!>        of a matrix M, so that when added, degenerate zero eigenvalues are shifted
!> \param filter the filter to be added to M
!> \param matrix_m the matrix M
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \note It is assumed that the matrix M is symmetric
! **************************************************************************************************
   SUBROUTINE get_gs_filter(filter, matrix_m, donor_state, xas_tdp_env)

      TYPE(dbcsr_type), POINTER                          :: filter, matrix_m
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'get_gs_filter', routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, i, iblk, jblk, nblk_row, ndo_mo
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size
      LOGICAL                                            :: found_block
      REAL(dp), DIMENSION(:), POINTER                    :: work_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type), POINTER                          :: sps

      NULLIFY (work_block, sps, row_blk_size)

! TODO: check that it works well when parallel (i.e. that the small and big matrix distribution
!       match). Do the same for build metric, and get_q_projector
!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices)
      sps => xas_tdp_env%p_projector
      CALL dbcsr_get_info(sps, row_blk_size=row_blk_size)
      nblk_row = SIZE(row_blk_size)

      CALL dbcsr_create(matrix=filter, template=matrix_m)

!  Fill the filter by looping over SPS and duplicating blocks. (all on the block diagonal)
      CALL dbcsr_iterator_start(iter, sps)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

         ! get the block
         CALL dbcsr_get_block_p(sps, iblk, jblk, work_block, found_block)

         IF (found_block .AND. jblk >= iblk) THEN

            DO i = 1, ndo_mo
               CALL dbcsr_put_block(filter, (i-1)*nblk_row+iblk, (i-1)*nblk_row+jblk, work_block)
            END DO

         END IF
         NULLIFY (work_block)

      END DO !iterator
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(filter)

   END SUBROUTINE get_gs_filter

! **************************************************************************************************
!> \brief Returns the projector on the unperturbed unoccupied ground state Q = 1 - SP on the block
!>        diagonal of a matrix with the stendard size and distribution.
!> \param proj_Q the matrix with the projector
!> \param dist_q the dbcsr distribution
!> \param blk_size_q the block size
!> \param donor_state ...
!> \param xas_tdp_env ...
! **************************************************************************************************
   SUBROUTINE get_q_projector(proj_Q, dist_q, blk_size_q, donor_state, xas_tdp_env)

      TYPE(dbcsr_type), POINTER                          :: proj_Q
      TYPE(dbcsr_distribution_type), POINTER             :: dist_q
      INTEGER, DIMENSION(:), POINTER                     :: blk_size_q
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'get_q_projector', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, i, iblk, jblk, nblk_row, ndo_mo
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size
      LOGICAL                                            :: found_block
      REAL(dp), DIMENSION(:), POINTER                    :: work_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type), POINTER                          :: one_sp

      NULLIFY (work_block, one_sp, row_blk_size)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices)
      one_sp => xas_tdp_env%q_projector
      CALL dbcsr_get_info(one_sp, row_blk_size=row_blk_size)
      nblk_row = SIZE(row_blk_size)

      ! the projector is not symmetric
      CALL dbcsr_create(matrix=proj_Q, name="PROJ Q", matrix_type="N", dist=dist_q, &
                        row_blk_size=blk_size_q, col_blk_size=blk_size_q)

!  Fill the filter by looping over 1-SP and duplicating blocks. (all on the block diagonal)
      CALL dbcsr_iterator_start(iter, one_sp)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

         ! get the block
         CALL dbcsr_get_block_p(one_sp, iblk, jblk, work_block, found_block)

         IF (found_block) THEN

            DO i = 1, ndo_mo
               CALL dbcsr_put_block(proj_Q, (i-1)*nblk_row+iblk, (i-1)*nblk_row+jblk, work_block)
            END DO

         END IF
         NULLIFY (work_block)

      END DO !iterator
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(proj_Q)

   END SUBROUTINE get_q_projector

! **************************************************************************************************
!> \brief Builds the matrix containing the ground state contribution to the matrix_tdp (aka matrix A)
!>         => A_{pi,qj} = F_pq*delta_ij - epsilon_ij*S_pq, where:
!>         F is the KS matrix
!>         S is the overlap matrix
!> \param matrix_a  pointer to a DBCSR matrix containing A
!> \param dist_a the dbcsr distribution of matrix A
!> \param blk_size_a the row/col block sizes of matrix A
!> \param donor_state ...
!> \param qs_env ...
!> \note Even localized non-canonical MOs are diagonalized in their subsapce => eps_ij = eps_ii*delta_ij
! **************************************************************************************************
   SUBROUTINE build_gs_contribution(matrix_a, dist_a, blk_size_a, donor_state, qs_env)

      TYPE(dbcsr_type), POINTER                          :: matrix_a
      TYPE(dbcsr_distribution_type), POINTER             :: dist_a
      INTEGER, DIMENSION(:), POINTER                     :: blk_size_a
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_gs_contribution', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, i, iblk, jblk, nblk_row, ndo_mo
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size
      LOGICAL                                            :: found_block
      REAL(dp), DIMENSION(:), POINTER                    :: energy_evals, work_block
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dbcsr_type)                                   :: work_matrix

      NULLIFY (matrix_ks, dbcsr_dist, row_blk_size, work_block, energy_evals, matrix_s)

!  Note: matrix A is symmetric. However, build it from non-symmetric basis because dbcsr distributions
!        only match for non-symmetric matrices (that is when the same dist is replicated). So we go
!        over the top triangular part and then redistribute into a symmetric matrix

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices, 1)
      CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_get_info(matrix_ks(1)%matrix, row_blk_size=row_blk_size)
      nblk_row = SIZE(row_blk_size)
      energy_evals => donor_state%energy_evals

!  Creating the symmetric matrix A (and work)
      CALL dbcsr_create(matrix=matrix_a, name="MATRIX A", matrix_type="S", dist=dist_a, &
                        row_blk_size=blk_size_a, col_blk_size=blk_size_a)
      CALL dbcsr_create(matrix=work_matrix, name="WORK MAT", matrix_type="S", dist=dist_a, &
                        row_blk_size=blk_size_a, col_blk_size=blk_size_a)

!  Loop over the blocks of KS and put them on the block diagonal of matrix A
      CALL dbcsr_iterator_start(iter, matrix_ks(1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

!        Get the block
         CALL dbcsr_get_block_p(matrix_ks(1)%matrix, iblk, jblk, work_block, found_block)

         IF (found_block .AND. jblk >= iblk) THEN

!           The KS matrix only appears on diagonal of matrix A => loop over II donor MOs
            DO i = 1, ndo_mo

!              Put the block as it is
               CALL dbcsr_put_block(matrix_a, (i-1)*nblk_row+iblk, (i-1)*nblk_row+jblk, work_block)

            END DO !i
         END IF !found_block
         NULLIFY (work_block)
      END DO ! iteration on KS blocks
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(matrix_a)

!  Loop over the blocks of S and put them on the block diagonal of work
      CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

!        Get the block
         CALL dbcsr_get_block_p(matrix_s(1)%matrix, iblk, jblk, work_block, found_block)

         IF (found_block .AND. jblk >= iblk) THEN

!           Add S matrix on block diagonal as epsilon_ii*S_pq
            DO i = 1, ndo_mo

               CALL dbcsr_put_block(work_matrix, (i-1)*nblk_row+iblk, (i-1)*nblk_row+jblk, &
                                    energy_evals(i)*work_block)
            END DO !i
         END IF !found block
         NULLIFY (work_block)
      END DO ! iteration on S blocks
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(work_matrix)

!  Take matrix_a = matrix_a - work
      CALL dbcsr_add(matrix_a, work_matrix, 1.0_dp, -1.0_dp)
      CALL dbcsr_finalize(matrix_a)

!  Clean-up
      CALL dbcsr_release(work_matrix)

   END SUBROUTINE build_gs_contribution

! **************************************************************************************************
!> \brief Creates the metric (aka  matrix G) needed for the generalized eigenvalue problem and inverse
!>         => G_{pi,qj} = S_pq*delta_ij
!> \param matrix_g dbcsr matrix containing G
!> \param dist_g dbcsr distribution of G
!> \param blk_size_g block sizes of G
!> \param donor_state ...
!> \param qs_env ...
!> \param do_inv if the inverse of G should be computed
! **************************************************************************************************
   SUBROUTINE build_metric(matrix_g, dist_g, blk_size_g, donor_state, qs_env, do_inv)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_g
      TYPE(dbcsr_distribution_type), POINTER          :: dist_g
      INTEGER, DIMENSION(:), POINTER                  :: blk_size_g
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      LOGICAL, OPTIONAL, INTENT(IN)                   :: do_inv

      CHARACTER(len=*), PARAMETER :: routineN = 'build_metric', routineP = moduleN//":"//routineN
      
      INTEGER                                         :: i, iblk, jblk, ndo_mo, nblk_row, blk, nao
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      INTEGER, DIMENSION(:), POINTER                  :: row_blk_size
      REAL(dp), DIMENSION(:), POINTER                 :: work_block
      LOGICAL                                         :: found_block
      TYPE(dbcsr_iterator_type)                       :: iter     
      TYPE(dbcsr_type), POINTER                       :: matrix_sinv
      LOGICAL                                         :: my_do_inv
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env

      NULLIFY(matrix_s, row_blk_size, work_block, para_env, blacs_env)

!  Initilization
      ndo_mo = SIZE(donor_state%mo_indices, 1)
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

         IF (found_block .AND. jblk >= iblk) THEN

!           Go over the diagonal of G => donor MOs ii
            DO i = 1, ndo_mo
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
         ALLOCATE(matrix_sinv)
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

            IF (found_block .AND. jblk >= iblk) THEN

!              Go over the diagonal of G => donor MOs ii
               DO i = 1, ndo_mo
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
         DEALLOCATE (matrix_sinv)
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

      REAL(dp), INTENT(IN)                            :: threshold, sx
      TYPE(dbcsr_type), POINTER                       :: matrix_a, matrix_d, matrix_e, proj_Q, work
      LOGICAL, INTENT(IN)                             :: do_hfx
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "build_aux_matrix", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: nrow, ndep
      REAL(dp)                                        :: evals(2)
      TYPE(dbcsr_type), POINTER                       :: tmp_mat
      TYPE(cp_fm_type), POINTER                       :: aux_fm, work_fm
      TYPE(cp_fm_struct_type), POINTER                :: struct
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(cp_para_env_type), POINTER                 :: para_env

      NULLIFY(tmp_mat, aux_fm, work_fm, struct, blacs_env, para_env)

      ALLOCATE (tmp_mat)
      CALL dbcsr_copy(tmp_mat, matrix_a)
      IF (do_hfx) THEN
         CALL dbcsr_add(tmp_mat, matrix_d, 1.0_dp, -1.0_dp*sx) !scaled hfx
         CALL dbcsr_add(tmp_mat, matrix_e, 1.0_dp, 1.0_dp*sx)
      END IF

      ! Take the product with the Q projector:
      CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, tmp_mat, 0.0_dp, work)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, tmp_mat)

      ! Actually computing and storing the auxiliary matrix
      ALLOCATE(donor_state%matrix_aux)
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
      DEALLOCATE (tmp_mat)

   END SUBROUTINE build_aux_matrix

! **************************************************************************************************
!> \brief Computes the two-center Coulomb integral needed for the RI in kernel calculation. Stores
!>        the integrals (P|Q)^-1 in the xas_tdp_env as global (small) arrays. Does that for a given
!>        excited kind
!> \param ex_kind ...
!> \param xas_tdp_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_ri_coulomb2_int(ex_kind, xas_tdp_env, qs_env)

      INTEGER, INTENT(IN)                                :: ex_kind
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

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
      ALLOCATE (xas_tdp_env%ri_inv_coul(nsgf, nsgf))

      CALL int_operators_r12_ab_os(r12_operator=operator_coulomb, vab=xas_tdp_env%ri_inv_coul, &
                                   rab=(/0.0_dp, 0.0_dp, 0.0_dp/), fba=ri_basis, fbb=ri_basis, &
                                   calculate_forces=.FALSE.)
!  Inverting
      CALL invmat_symm(xas_tdp_env%ri_inv_coul)

   END SUBROUTINE compute_ri_coulomb2_int

! **************************************************************************************************
!> \brief Computes the RI Coulomb 3-center integrals (ab|c), where c is from the RI_XAS basis and
!>        centered on the given excited atom
!> \param ex_atom the excited atom on which the third center is
!> \param ex_kind the kind of the excited atom
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note  The ri_o3c_coul member of xas_tdp_env is defined and allocated here
! **************************************************************************************************
   SUBROUTINE compute_ri_o3c_coulomb(ex_atom, ex_kind, xas_tdp_env, qs_env)

      INTEGER, INTENT(IN)                                :: ex_atom, ex_kind
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_o3c_coulomb', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ilist, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: exat, exki
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_orb, basis_set_ri
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ac_list, sab_orb
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (basis_set_ri, basis_set_orb, ac_list, sab_orb, qs_kind_set, para_env)

      CALL timeset(routineN, handle)

!  Release the previous o3c centered on another atom
      IF (ASSOCIATED(xas_tdp_env%ri_o3c_coul)) THEN
         CALL release_o3c_container(xas_tdp_env%ri_o3c_coul)
         DEALLOCATE (xas_tdp_env%ri_o3c_coul)
      END IF

!  Take what we need from the qs_env
      CALL get_qs_env(qs_env, nkind=nkind, sab_orb=sab_orb, qs_kind_set=qs_kind_set, para_env=para_env)

!  Build the basis set lists
      ALLOCATE (basis_set_ri(nkind))
      ALLOCATE (basis_set_orb(nkind))
      CALL basis_set_list_setup(basis_set_ri, "RI_XAS", qs_kind_set)
      CALL basis_set_list_setup(basis_set_orb, "ORB", qs_kind_set)

!  Build a neighbor list for the ac centers. Since we later contract as (aI|c) and we assume I is
!  very localized on the same atom as c, we take a,c as neighbors if they overlap => the effective
!  range is then second nearest neightbors (relatively weak assumption)
      ALLOCATE (exat(1), exki(1))
      exat(1) = ex_atom; exki(1) = ex_kind; 
      CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, "OVERLAP", &
                               exat, exki, qs_env)

!  Allocate, init and compute the integrals
      ALLOCATE (xas_tdp_env%ri_o3c_coul)
      CALL init_o3c_container(xas_tdp_env%ri_o3c_coul, 1, basis_set_orb, basis_set_orb, &
                              basis_set_ri, sab_orb, ac_list)
      CALL calculate_o3c_coulomb_integrals(xas_tdp_env%ri_o3c_coul)

!  Get th lower and upper bonds of a,b in (ab|P) integrals and sotre in xas_tdp_env
      CALL get_o3c_container(xas_tdp_env%ri_o3c_coul, lbij=xas_tdp_env%lb_o3c_coul, &
                             ubij=xas_tdp_env%ub_o3c_coul)
      CALL mp_min(xas_tdp_env%lb_o3c_coul, para_env%group)
      CALL mp_max(xas_tdp_env%ub_o3c_coul, para_env%group)

! Clean-up
      DO ilist = 1, SIZE(ac_list)
         CALL deallocate_neighbor_list_set(ac_list(ilist)%neighbor_list_set)
      END DO
      DEALLOCATE (basis_set_ri, basis_set_orb, ac_list)

      CALL timestop(handle)

   END SUBROUTINE compute_ri_o3c_coulomb

! **************************************************************************************************
!> \brief Computes the RI exchange 3-center integrals (ab|c), where c is from the RI_XAS basis and
!>        centered on the given excited atom and kind
!> \param ex_atom the excited atom on which the third center is
!> \param ex_kind the kind of the excited atom
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note  The ri_o3c_ex member of xas_tdp_env is defined and allocated here
! **************************************************************************************************
   SUBROUTINE compute_ri_o3c_exchange(ex_atom, ex_kind, xas_tdp_env, qs_env)

      INTEGER, INTENT(IN)                                :: ex_atom, ex_kind
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_o3c_exchange', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ilist, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: exat, exki
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_orb, basis_set_ri
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ac_list, sab_orb
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (basis_set_ri, basis_set_orb, ac_list, sab_orb, qs_kind_set)

      CALL timeset(routineN, handle)

!  Release the previous o3c centered on another atom
      IF (ASSOCIATED(xas_tdp_env%ri_o3c_ex)) THEN
         CALL release_o3c_container(xas_tdp_env%ri_o3c_ex)
         DEALLOCATE (xas_tdp_env%ri_o3c_ex)
      END IF

!  Take what we need from the qs_env
      CALL get_qs_env(qs_env, nkind=nkind, sab_orb=sab_orb, qs_kind_set=qs_kind_set)

!  Build the basis set lists
      ALLOCATE (basis_set_ri(nkind))
      ALLOCATE (basis_set_orb(nkind))
      CALL basis_set_list_setup(basis_set_ri, "RI_XAS", qs_kind_set)
      CALL basis_set_list_setup(basis_set_orb, "ORB", qs_kind_set)

!  Build a neighbor list for the ac centers. Currently only coulomb operator (with "inifinte" range)
      ALLOCATE (exat(1), exki(1))
      exat(1) = ex_atom; exki(1) = ex_kind; 
      CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, "COULOMB", &
                               exat, exki, qs_env)

!  Allocate, init and compute the integrals
      ALLOCATE (xas_tdp_env%ri_o3c_ex)
      CALL init_o3c_container(xas_tdp_env%ri_o3c_ex, 1, basis_set_orb, basis_set_orb, &
                              basis_set_ri, sab_orb, ac_list)
      CALL calculate_o3c_coulomb_integrals(xas_tdp_env%ri_o3c_ex)

! Clean-up
      DO ilist = 1, SIZE(ac_list)
         CALL deallocate_neighbor_list_set(ac_list(ilist)%neighbor_list_set)
      END DO
      DEALLOCATE (basis_set_ri, basis_set_orb, ac_list)

      CALL timestop(handle)

   END SUBROUTINE compute_ri_o3c_exchange

! **************************************************************************************************
!> \brief Contract the ri 3-center integrals contained in o3c with repect to the donor MOs coeffs,
!>        for a given excited atom k => (aI|k) = sum_b c_Ib (ab|k)
!> \param contr_int the contracted integrals as array of dbcsr matrices
!> \param op_type for which operator type we contract (COULOMB or EXCHANGE)
!> \param donor_state ...
!> \param qs_env ...
!> \note  In the output matrices, (aI_b|k) is stored at block a,b where I_b is the partial
!>        contraction that only includes coeffs from atom b. Note that the contracted matrix is
!>        not diagonal, but the blocks (aI_b|P) and (bI_a|P) are strictly identical
! **************************************************************************************************
   SUBROUTINE contract_o3c_int(contr_int, op_type, donor_state, xas_tdp_env, qs_env)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: contr_int
      CHARACTER(len=*), INTENT(IN)                       :: op_type
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'contract_o3c_int', &
         routineP = moduleN//':'//routineN

      CHARACTER                                          :: my_op
      INTEGER                                            :: er, handle, iat, imo, katom, kkind, lb, &
                                                            natom, ndo_mo, nsgfk, sr, ub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nsgf_orb
      INTEGER, DIMENSION(:), POINTER                     :: ri_blk_size, std_blk_size
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: coeffs
      REAL(dp), DIMENSION(:), POINTER                    :: v
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrices, matrix_s
      TYPE(dbcsr_type), POINTER                          :: aI_bP, PI_ab, work1, work2
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(o3c_container_type), POINTER                  :: o3c_coul, o3c_ex
      TYPE(o3c_vec_type), DIMENSION(:), POINTER          :: o3c_vec
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (matrix_s, std_blk_size, ri_blk_size, dbcsr_dist, qs_kind_set, ri_basis, o3c_vec)
      NULLIFY (aI_bP, PI_ab, work1, work2, matrices)

      CALL timeset(routineN, handle)

!  Initialization
      CALL get_qs_env(qs_env, natom=natom, matrix_s=matrix_s, qs_kind_set=qs_kind_set, &
                      dbcsr_dist=dbcsr_dist)
      ndo_mo = SIZE(donor_state%mo_indices)
      kkind = donor_state%kind_index
      katom = donor_state%at_index
      !by default contract for Coulomb
      my_op = "C"
      o3c_coul => xas_tdp_env%ri_o3c_coul
      IF (op_type == "EXCHANGE") THEN
         my_op = "E"
         CPASSERT(ASSOCIATED(xas_tdp_env%ri_o3c_ex))
         o3c_ex => xas_tdp_env%ri_o3c_ex
      END IF

!  contracted integrals block sizes
      CALL dbcsr_get_info(matrix_s(1)%matrix, col_blk_size=std_blk_size)
      ! getting the block dimensions for the RI basis
      CALL get_qs_kind(qs_kind_set(kkind), basis_set=ri_basis, basis_type="RI_XAS")
      CALL get_gto_basis_set(ri_basis, nsgf=nsgfk)
      ALLOCATE (ri_blk_size(natom))
      ri_blk_size = nsgfk

!  Get some info about the o3c_coul. Whether we contract for coulomb or exchange, MOs are assumed
!  to be localized as in the Coulomb integrals
      lb = xas_tdp_env%lb_o3c_coul
      ub = xas_tdp_env%ub_o3c_coul

!  Create  work matrices. They must have the block distribution of a symmetric matric to enter the
!  o3c routines. HACK => create a symmetric matrix with non-symmetric block size (aI_bP, PI_ab)
!  Also create their normal matrix conterparts (work1, work2)
      ALLOCATE (aI_bP, PI_ab, work1, work2, matrices(2))
      CALL dbcsr_create(aI_bP, dist=dbcsr_dist, matrix_type="S", name="(aI_b|P)", &
                        row_blk_size=std_blk_size, col_blk_size=ri_blk_size)
      CALL dbcsr_create(work1, dist=dbcsr_dist, matrix_type="N", name="WORK 1", &
                        row_blk_size=std_blk_size, col_blk_size=ri_blk_size)

      CALL dbcsr_create(PI_ab, dist=dbcsr_dist, matrix_type="S", name="(P|I_ab)", &
                        row_blk_size=ri_blk_size, col_blk_size=std_blk_size)
      CALL dbcsr_create(work2, dist=dbcsr_dist, matrix_type="N", name="WORK 2", &
                        row_blk_size=ri_blk_size, col_blk_size=std_blk_size)

      !reserve the blocks accordingly
      matrices(1)%matrix => aI_bP; matrices(2)%matrix => PI_ab
      CALL reserve_contraction_blocks(matrices, lb, ub, qs_env)

      matrices(1)%matrix => work1; matrices(2)%matrix => work2
      CALL reserve_contraction_blocks(matrices, lb, ub, qs_env, matrix_type="N")
      DEALLOCATE (matrices)

      !  Create the contracted integral matrices
      ALLOCATE (contr_int(ndo_mo))
      DO imo = 1, ndo_mo
         ALLOCATE (contr_int(imo)%matrix)
      END DO

!  Prepare the o3c_vec with the MO coeffs
      ALLOCATE (nsgf_orb(natom), o3c_vec(natom))
      nsgf_orb = 1
      nsgf_orb(lb:ub) = std_blk_size(lb:ub)
      CALL o3c_vec_create(o3c_vec, nsgf_orb)

      ! Only take the coeffs for atoms that are present in o3c_coul
      sr = 1; IF (lb > 1) sr = SUM(std_blk_size(1:lb-1))+1
      er = SUM(std_blk_size(1:ub))

      ALLOCATE (coeffs(er-sr+1, ndo_mo))
      CALL cp_fm_get_submatrix(fm=donor_state%gs_coeffs, target_m=coeffs, start_row=sr, &
                               n_rows=er-sr+1, n_cols=ndo_mo)

!  Loop over the donor MOs, fill the o3c_vec and contract
      DO imo = 1, ndo_mo

         DO iat = lb, ub
            sr = 1; IF (iat > lb) sr = SUM(std_blk_size(lb:iat-1))+1
            er = SUM(std_blk_size(lb:iat))

            CALL get_o3c_vec(o3c_vec, iat, v)
            v(:) = coeffs(sr:er, imo)
         END DO !iat

         ! do the contraction
         CALL dbcsr_set(aI_bP, 0.0_dp); CALL dbcsr_set(PI_ab, 0.0_dp)
         IF (my_op == "C") THEN
            CALL contract_o3c_once(o3c_coul, o3c_vec, aI_bP, PI_ab, katom)
         ELSE
            CALL contract_o3c_once(o3c_ex, o3c_vec, aI_bP, PI_ab, katom)
         END IF

         ! Get the full "normal" (aI_b|P) contracted integrals
         CALL dbcsr_copy_into_existing(work2, PI_ab)
         CALL dbcsr_transposed(contr_int(imo)%matrix, work2)
         CALL dbcsr_copy_into_existing(work1, aI_bP)
         CALL dbcsr_add(contr_int(imo)%matrix, work1, 1.0_dp, 1.0_dp)

      END DO !imo

!  Clean-up
      CALL o3c_vec_release(o3c_vec)
      CALL dbcsr_release(aI_bP)
      CALL dbcsr_release(PI_ab)
      CALL dbcsr_release(work1)
      CALL dbcsr_release(work2)
      DEALLOCATE (aI_bP, PI_ab, o3c_vec, ri_blk_size, work1, work2)

      CALL timestop(handle)

   END SUBROUTINE contract_o3c_int

! **************************************************************************************************
!> \brief Reserves the blocks in of a dbcsr matrix as needed for RI 3-center contraction (aI_b|P)
!> \param matrices the matrices for which blocks are reserved
!> \param lb the lower bond of atoms a,b interaction with basis P
!> \param ub the upper bond of atoms a,b interaction with basis P
!> \param qs_env ...
!> \param matrix_type whether the matrices are symmetric or not
!> \param block_content in the case where all the blocks have the same dimension, can fill them
!>        with the values taken from this array
! **************************************************************************************************
   SUBROUTINE reserve_contraction_blocks(matrices, lb, ub, qs_env, matrix_type, block_content)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrices
      INTEGER, INTENT(IN)                                :: lb, ub
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER, INTENT(IN), OPTIONAL                    :: matrix_type
      REAL(dp), DIMENSION(:, :), INTENT(IN), OPTIONAL    :: block_content

      CHARACTER(len=*), PARAMETER :: routineN = 'reserve_contraction_blocks', &
         routineP = moduleN//':'//routineN

      CHARACTER                                          :: my_type
      INTEGER                                            :: blk, i, iblk, jblk, n, s1, s2
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock_m, pblock_s
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: template

      NULLIFY (matrix_s, pblock_s, pblock_m, template)

!  Initialization
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      IF (PRESENT(block_content)) THEN
         s1 = SIZE(block_content, 1)
         s2 = SIZE(block_content, 2)
      END IF
      n = SIZE(matrices)
      !assume symmetric type
      my_type = dbcsr_type_symmetric
      IF (PRESENT(matrix_type)) my_type = matrix_type

!  If non-symmetric, need to desymmetrize matrix as as a template
      IF (my_type .NE. dbcsr_type_symmetric) THEN
         ALLOCATE (template)
         CALL dbcsr_desymmetrize(matrix_s(1)%matrix, template)
      ELSE
         template => matrix_s(1)%matrix
      END IF

!  Loop over matrix_s as need a,b to overlap
      CALL dbcsr_iterator_start(iter, template)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
         IF (iblk < lb .OR. jblk < lb) CYCLE
         IF (iblk > ub .OR. jblk > ub) CYCLE

         CALL dbcsr_get_block_p(template, iblk, jblk, pblock_s, found)

         IF (found) THEN
            DO i = 1, n
               NULLIFY (pblock_m)
               CALL dbcsr_reserve_block2d(matrices(i)%matrix, iblk, jblk, pblock_m)
               IF (PRESENT(block_content)) THEN
                  CALL dcopy(s1*s2, block_content, 1, pblock_m, 1)
               ELSE
                  pblock_m = 0.0_dp
               END IF
            END DO
         END IF

      END DO !dbcsr iter
      CALL dbcsr_iterator_stop(iter)
      DO i = 1, n
         CALL dbcsr_finalize(matrices(i)%matrix)
      END DO

!  Clean-up
      IF (my_type .NE. dbcsr_type_symmetric) THEN
         CALL dbcsr_release(template)
         DEALLOCATE (template)
      END IF

   END SUBROUTINE reserve_contraction_blocks

! **************************************************************************************************
!> \brief Contraction of the 3-center integrals over index 1 and 2, for a given atom_k. The results
!>        are stored in two matrices, such that (a,b are block indices):
!>        mat_aIb(ab) = mat_aIb(ab) + sum j_b (i_aj_b|k)*v(j_b) and
!>        mat_bIa(ba) = mat_bIa(ba) + sum i_a (i_aj_b|k)*v(i_a)
!>        The block size of the columns of mat_aIb and the rows of mat_bIa are the size of k
!> \param o3c the container for the integrals
!> \param vec the contraction coefficients
!> \param mat_aIb ...
!> \param mat_bIa ...
!> \param atom_k the atom for which we contract
!> \note To get the full (aI|P) matrix, one has to sum mat_aIb + mat_bIa^T (the latter has empty diag)
! **************************************************************************************************
   SUBROUTINE contract_o3c_once(o3c, vec, mat_aIb, mat_bIa, atom_k)
      TYPE(o3c_container_type), POINTER                  :: o3c
      TYPE(o3c_vec_type), DIMENSION(:), POINTER          :: vec
      TYPE(dbcsr_type), POINTER                          :: mat_aIb, mat_bIa
      INTEGER, INTENT(IN)                                :: atom_k

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_o3c_once', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iatom, icol, irow, jatom, &
                                                            katom, mepos, n, nthread, s1, s2
      INTEGER, DIMENSION(:), POINTER                     :: atom_blk_size
      LOGICAL                                            :: found, ijsymmetric, trans
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: work
      REAL(dp), DIMENSION(:), POINTER                    :: v
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock
      REAL(dp), DIMENSION(:, :, :), POINTER              :: iabc
      TYPE(o3c_iterator_type)                            :: o3c_iterator

      CALL timeset(routineN, handle)

      CALL get_o3c_container(o3c, ijsymmetric=ijsymmetric)
      CPASSERT(ijsymmetric)

      CALL dbcsr_get_info(mat_aIb, row_blk_size=atom_blk_size)

      nthread = 1
!$    nthread = omp_get_max_threads()
      CALL o3c_iterator_create(o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,o3c_iterator,vec,mat_aIb,mat_bIa,atom_k,atom_blk_size)&
!$OMP PRIVATE (mepos,iabc,iatom,jatom,katom,irow,icol,trans,pblock,v,found,i,n,s1,s2,work)

      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, iatom=iatom, jatom=jatom, katom=katom, &
                                    integral=iabc)

         IF (katom .NE. atom_k) CYCLE

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
            trans = .FALSE.
         ELSE
            irow = jatom
            icol = iatom
            trans = .TRUE.
         END IF

         ! Deal with mat_aIb
         CALL get_o3c_vec(vec, icol, v)
         n = SIZE(v)
         s1 = atom_blk_size(irow)
         s2 = SIZE(iabc, 3)
         ALLOCATE (work(s1, s2))
         work = 0.0_dp

         IF (trans) THEN
            DO i = 1, n
               CALL daxpy(s1*s2, v(i), iabc(i, :, :), 1, work(:, :), 1)
            END DO
         ELSE
            DO i = 1, n
               CALL daxpy(s1*s2, v(i), iabc(:, i, :), 1, work(:, :), 1)
            END DO
         END IF

!$OMP CRITICAL(mat_aIb)
         CALL dbcsr_get_block_p(matrix=mat_aIb, row=irow, col=icol, BLOCK=pblock, found=found)
         IF (found) THEN
            CALL daxpy(s1*s2, 1.0_dp, work(:, :), 1, pblock(:, :), 1)
         END IF
!$OMP END CRITICAL(mat_aIb)
         DEALLOCATE (work)

         ! Deal with mat_bIa, keep block diagonal empty
         IF (irow == icol) CYCLE

         CALL get_o3c_vec(vec, irow, v)
         n = SIZE(v)
         s1 = SIZE(iabc, 3)
         s2 = atom_blk_size(icol)
         ALLOCATE (work(s1, s2))
         work = 0.0_dp

         IF (trans) THEN
            DO i = 1, n
               CALL daxpy(s1*s2, v(i), TRANSPOSE(iabc(:, i, :)), 1, work(:, :), 1)
            END DO
         ELSE
            DO i = 1, n
               CALL daxpy(s1*s2, v(i), TRANSPOSE(iabc(i, :, :)), 1, work(:, :), 1)
            END DO
         END IF

!$OMP CRITICAL(mat_bIa)
         CALL dbcsr_get_block_p(matrix=mat_bIa, row=irow, col=icol, BLOCK=pblock, found=found)
         IF (found) THEN
            CALL daxpy(s1*s2, 1.0_dp, work(:, :), 1, pblock(:, :), 1)
         END IF
!$OMP END CRITICAL(mat_bIa)
         DEALLOCATE (work)

      END DO !o3c iterator
!$OMP END PARALLEL
      CALL o3c_iterator_release(o3c_iterator)

      CALL timestop(handle)

   END SUBROUTINE contract_o3c_once

! **************************************************************************************************
!> \brief Computes, if asked for it, the Coulomb and XC kernel matrices, in the usuall matrix format
!> \param coulomb_kernel pointer the the Coulomb kernel matrix (can be void pointer)
!> \param xc_kernel_sg pointer the the XC singlet kernel matrix (can be void pointer)
!> \param xc_kernel_tp pointer the the XC triplet kernel matrix (can be void pointer)
!> \param do_coulomb whether the coulomb kernel should actually be computed
!> \param do_xc_sg whether the XC kernel should be computed for singlets
!> \param do_xc_tp whether the XC kernel should be computed for triplets
!> \param dist the dbcsr_dist of the kernel matrices
!> \param blk_size the block size of the kernel matrices
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note Coulomb and xc kernel are put together in the same routine because they use the same RI
!>       Coulomb: (aI|Jb) = (aI|P) (P|Q)^-1 (Q|Jb)
!>       XC : (aI|fxc|Jb) = (aI|P) (P|Q)^-1 (Q|fxc|R) (R|S)^-1 (S|Jb)
! **************************************************************************************************
   SUBROUTINE kernel_coulomb_xc(coulomb_kernel, xc_kernel_sg, xc_kernel_tp, do_coulomb, do_xc_sg, &
                                do_xc_tp, dist, blk_size, donor_state, xas_tdp_env, qs_env)

      TYPE(dbcsr_type), POINTER                          :: coulomb_kernel, xc_kernel_sg, &
                                                            xc_kernel_tp
      LOGICAL, INTENT(IN)                                :: do_coulomb, do_xc_sg, do_xc_tp
      TYPE(dbcsr_distribution_type), POINTER             :: dist
      INTEGER, DIMENSION(:), POINTER                     :: blk_size
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'kernel_coulomb_xc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: bo(2), handle, iex, imo, ip, natom, &
                                                            ndo_mo, nex_atom, nsgfp, ri_atom, &
                                                            source
      INTEGER, DIMENSION(:), POINTER                     :: ri_blk_size
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: fxc_rs
      REAL(dp), DIMENSION(:, :), POINTER                 :: fxc, PQ
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: lhs_int, lhs_tp, matrices, rhs_int
      TYPE(dbcsr_type), POINTER                          :: work_mat

      NULLIFY (lhs_int, rhs_int, PQ, ri_blk_size, dbcsr_dist, lhs_tp, matrices, work_mat, para_env, fxc)

      CALL timeset(routineN, handle)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices)
      ri_atom = donor_state%at_index
      CALL get_qs_env(qs_env, natom=natom, dbcsr_dist=dbcsr_dist)

!  Get the RI inverse coulomb for the excited kind (P|Q)^-1
      PQ => xas_tdp_env%ri_inv_coul
      nsgfp = SIZE(PQ, 1)

!  Get the contracted 3-center integrals in lhs_int and copy them into rhs_int
      CALL contract_o3c_int(lhs_int, "COULOMB", donor_state, xas_tdp_env, qs_env)
      ALLOCATE (rhs_int(ndo_mo))
      DO imo = 1, ndo_mo
         ALLOCATE (rhs_int(imo)%matrix)
         CALL dbcsr_copy(rhs_int(imo)%matrix, lhs_int(imo)%matrix)
      END DO

!  In both Coulomb and XC need to multiply lhs_int by (P|Q)^-1
!  Coming out of contract_o3c_int, the intrgrals are (aI_b|P), so need to sum over b as well
!  => create a symmetric matrix with the same block dist as contr_int and (P|Q)^-1 everywhere
      ALLOCATE (ri_blk_size(natom), work_mat)
      ri_blk_size = nsgfp
      CALL dbcsr_create(work_mat, name="WORK", matrix_type="S", dist=dbcsr_dist, &
                        row_blk_size=ri_blk_size, col_blk_size=ri_blk_size)
      ALLOCATE (matrices(1))
      matrices(1)%matrix => work_mat
      CALL reserve_contraction_blocks(matrices, xas_tdp_env%lb_o3c_coul, xas_tdp_env%ub_o3c_coul, &
                                      qs_env, block_content=PQ)
      DEALLOCATE (matrices)

      ! After this all rows of LHS contain (aI|P) integrals (identitcal blocks)
      DO imo = 1, ndo_mo
         CALL dbcsr_multiply('N', 'N', 1.0_dp, lhs_int(imo)%matrix, work_mat, 0.0_dp, lhs_int(imo)%matrix)
      END DO

      CALL dbcsr_release(work_mat)
      DEALLOCATE (ri_blk_size)

!  Deal with the Coulomb case
      IF (do_coulomb) THEN

         ! Create a normal type work matrix
         CALL dbcsr_create(work_mat, name="WORK", matrix_type="N", dist=dist, &
                           row_blk_size=blk_size, col_blk_size=blk_size)

         !Compute the product (rhs_int is still in (aI_b|P) format, which is what we want)
         CALL ri_int_product(work_mat, lhs_int, rhs_int, qs_env)

         !Create the symmetric kernel matrix and redistribute work_mat into it
         CALL dbcsr_create(coulomb_kernel, name="COULOMB KERNEL", matrix_type="S", dist=dist, &
                           row_blk_size=blk_size, col_blk_size=blk_size)
         CALL dbcsr_complete_redistribute(work_mat, coulomb_kernel)

         CALL dbcsr_release(work_mat)

      END IF !do_coulomb

!  Deal with the XC case
      IF (do_xc_sg .OR. do_xc_tp) THEN

         ! some work array
         ALLOCATE (fxc_rs(nsgfp, nsgfp))

         CALL get_qs_env(qs_env, para_env=para_env)

         !If both singlet AND triplet, need to save LHS before (Q|fxc|R) multiplication
         IF (do_xc_sg .AND. do_xc_tp) THEN
            ALLOCATE (lhs_tp(ndo_mo))
            DO imo = 1, ndo_mo
               ALLOCATE (lhs_tp(imo)%matrix)
               CALL dbcsr_copy(lhs_tp(imo)%matrix, lhs_int(imo)%matrix)
            END DO
         END IF

         ! Find on which processor the integral (Q|fxc|R) is stored
         found = .FALSE.
         nex_atom = SIZE(xas_tdp_env%ex_atom_indices)
         DO ip = 0, para_env%num_pe-1

            bo = get_limit(nex_atom, para_env%num_pe, ip)
            DO iex = bo(1), bo(2)

               IF (xas_tdp_env%ex_atom_indices(iex) == ri_atom) THEN
                  source = ip
                  found = .TRUE.
                  EXIT
               END IF
            END DO !iex
            IF (found) EXIT
         END DO !ip

         ! Case study on excited state multiplicity
         IF (do_xc_sg) THEN

            ! Fetch (Q|fxc|R) and broadcast it to all processors
            IF (.NOT. ASSOCIATED(xas_tdp_env%ri_fxc_sg(ri_atom)%array)) THEN
               ALLOCATE (xas_tdp_env%ri_fxc_sg(ri_atom)%array(nsgfp, nsgfp))
            END IF
            CALL mp_bcast(xas_tdp_env%ri_fxc_sg(ri_atom)%array, source, para_env%group)
            fxc => xas_tdp_env%ri_fxc_sg(ri_atom)%array

            ! Need to further multiply LHS by (Q|fxc|R)*(R|S)^-1
            CALL dgemm('N', 'N', nsgfp, nsgfp, nsgfp, 1.0_dp, fxc, nsgfp, PQ, nsgfp, 0.0_dp, fxc_rs, nsgfp)
            CALL ri_all_blocks_mm(lhs_int, fxc_rs)

            ! Create a normal type work matrix
            CALL dbcsr_create(work_mat, name="WORK", matrix_type="N", dist=dist, &
                              row_blk_size=blk_size, col_blk_size=blk_size)

            ! Compute the product (rhs_int is still in (aI_b|P) format, which is what we want)
            CALL ri_int_product(work_mat, lhs_int, rhs_int, qs_env)

            !Create the symmetric kernel matrix and redistribute work_mat into it
            CALL dbcsr_create(xc_kernel_sg, name="XC SINGLET KERNEL", matrix_type="S", dist=dist, &
                              row_blk_size=blk_size, col_blk_size=blk_size)
            CALL dbcsr_complete_redistribute(work_mat, xc_kernel_sg)

            CALL dbcsr_release(work_mat)

         END IF ! singlet

         IF (do_xc_tp) THEN

            ! Need a fresh set of LHS integrals. If no singlet, simply point to initial LHS
            IF (.NOT. do_xc_sg) THEN
               lhs_tp => lhs_int
            END IF

            ! Fetch (Q|fxc|R) and broadcast it to all processors
            IF (.NOT. ASSOCIATED(xas_tdp_env%ri_fxc_tp(ri_atom)%array)) THEN
               ALLOCATE (xas_tdp_env%ri_fxc_tp(ri_atom)%array(nsgfp, nsgfp))
            END IF
            CALL mp_bcast(xas_tdp_env%ri_fxc_tp(ri_atom)%array, source, para_env%group)
            fxc => xas_tdp_env%ri_fxc_tp(ri_atom)%array

            ! Need to further multiply LHS by (Q|fxc|R)*(R|S)^-1
            CALL dgemm('N', 'N', nsgfp, nsgfp, nsgfp, 1.0_dp, fxc, nsgfp, PQ, nsgfp, 0.0_dp, fxc_rs, nsgfp)
            CALL ri_all_blocks_mm(lhs_tp, fxc)

            ! Create a normal type work matrix
            CALL dbcsr_create(work_mat, name="WORK", matrix_type="N", dist=dist, &
                              row_blk_size=blk_size, col_blk_size=blk_size)

            ! Compute the product (rhs_int is still in (aI_b|P) format, which is what we want)
            CALL ri_int_product(work_mat, lhs_tp, rhs_int, qs_env)

            !Create the symmetric kernel matrix and redistribute work_mat into it
            CALL dbcsr_create(xc_kernel_tp, name="XC TRIPLET KERNEL", matrix_type="S", dist=dist, &
                              row_blk_size=blk_size, col_blk_size=blk_size)
            CALL dbcsr_complete_redistribute(work_mat, xc_kernel_tp)

            CALL dbcsr_release(work_mat)

         END IF ! triplet
      END IF ! do_xc

!  Clean-up
      DO imo = 1, ndo_mo
         CALL dbcsr_release(lhs_int(imo)%matrix)
         CALL dbcsr_release(rhs_int(imo)%matrix)
         DEALLOCATE (lhs_int(imo)%matrix, rhs_int(imo)%matrix)
         IF (do_xc_sg .AND. do_xc_tp) THEN
            CALL dbcsr_release(lhs_tp(imo)%matrix)
            DEALLOCATE (lhs_tp(imo)%matrix)
         END IF
      END DO
      DEALLOCATE (lhs_int, rhs_int, work_mat)
      IF (do_xc_sg .AND. do_xc_tp) DEALLOCATE (lhs_tp)
      ! also clean up the fxc integral, not needed anymore
      IF (do_xc_sg) DEALLOCATE (xas_tdp_env%ri_fxc_sg(ri_atom)%array)
      IF (do_xc_tp) DEALLOCATE (xas_tdp_env%ri_fxc_tp(ri_atom)%array)

      CALL timestop(handle)

   END SUBROUTINE kernel_coulomb_xc

! **************************************************************************************************
!> \brief Computes the exact exchange kernel matrix using RI. Computes the off-diagonal exchange
!>        only is asked for it (for full tddft)
!> \param ex_ker_ondiag the diagonal part of the kernel => (ab|IJ)
!> \param ex_ker_offdiag the off-diagonal part of the kernel => (aJ|Ib)
!> \param do_offdiag whether the offdiag part should be computed
!> \param dist the usuall dbcsr distribution
!> \param blk_size the block size
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note The RI for the diagonal and off-diagonal exchange kernels are respectively
!>       (ab|IJ) = (ab|P) (P|Q)^-1 (Q|IJ)
!>       (aJ|Ib) = (aJ|P) (P|Q)^-1 (Q|Ib)
! **************************************************************************************************
   SUBROUTINE kernel_exchange(ex_ker_ondiag, ex_ker_offdiag, do_offdiag, dist, blk_size, &
                              donor_state, xas_tdp_env, qs_env)

      TYPE(dbcsr_type), POINTER                          :: ex_ker_ondiag, ex_ker_offdiag
      LOGICAL, INTENT(IN)                                :: do_offdiag
      TYPE(dbcsr_distribution_type), POINTER             :: dist
      INTEGER, DIMENSION(:), POINTER                     :: blk_size
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'kernel_exchange', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, er, handle, iblk, imo, jblk, jmo, &
                                                            natom, nblk, ndo_mo, nsgfp, ri_atom, &
                                                            s2, sr
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, ri_blk_size, vec_blk_size
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: coeffs
      REAL(dp), DIMENSION(:), POINTER                    :: v
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock, PQ
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: contr2_int
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: lhs_int, matrices, matrix_s, rhs_int
      TYPE(dbcsr_type)                                   :: lhs_t
      TYPE(dbcsr_type), POINTER                          :: abIJ, abIJ_desymm, work_mat
      TYPE(o3c_vec_type), DIMENSION(:), POINTER          :: o3c_vec

      NULLIFY (PQ, lhs_int, rhs_int, pblock, fm_struct, contr2_int, para_env, blacs_env, col_blk_size)
      NULLIFY (o3c_vec, v, matrix_s, abIJ, abIJ_desymm, dbcsr_dist, work_mat, matrices, vec_blk_size)

      CALL timeset(routineN, handle)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices)
      ri_atom = donor_state%at_index
      CALL get_qs_env(qs_env, natom=natom, dbcsr_dist=dbcsr_dist)

!  Get the RI inverse coulomb for the excited kind (P|Q)^-1
      PQ => xas_tdp_env%ri_inv_coul
      nsgfp = SIZE(PQ, 1)

!  Need the once contracted integrals (aI_b|P) in any case
      CALL contract_o3c_int(lhs_int, "EXCHANGE", donor_state, xas_tdp_env, qs_env)

      ! Also need a copy if off-diag
      IF (do_offdiag) THEN
         ALLOCATE (rhs_int(ndo_mo))
         DO imo = 1, ndo_mo
            ALLOCATE (rhs_int(imo)%matrix)
            CALL dbcsr_copy(rhs_int(imo)%matrix, lhs_int(imo)%matrix)
         END DO
      END IF

!  Multiply lhs by (P|Q)^-1 => Play same trick as in coulomb_xc kernel
      ALLOCATE (ri_blk_size(natom), work_mat)
      ri_blk_size = nsgfp
      CALL dbcsr_create(work_mat, name="WORK", matrix_type="S", dist=dbcsr_dist, &
                        row_blk_size=ri_blk_size, col_blk_size=ri_blk_size)
      ALLOCATE (matrices(1))
      matrices(1)%matrix => work_mat
      CALL reserve_contraction_blocks(matrices, xas_tdp_env%lb_o3c_coul, xas_tdp_env%ub_o3c_coul, &
                                      qs_env, block_content=PQ)
      DEALLOCATE (matrices)

      DO imo = 1, ndo_mo
         CALL dbcsr_multiply('N', 'N', 1.0_dp, lhs_int(imo)%matrix, work_mat, 0.0_dp, lhs_int(imo)%matrix)
      END DO

      CALL dbcsr_release(work_mat)
      DEALLOCATE (ri_blk_size)

!  Deal with the off-diagonal case fisrt
      IF (do_offdiag) THEN

         ! Create a normal type work matrix
         CALL dbcsr_create(work_mat, name="WORK", matrix_type="N", dist=dist, &
                           row_blk_size=blk_size, col_blk_size=blk_size)

         ! Compute the product with transposed MO blocks
         CALL ri_int_product(work_mat, lhs_int, rhs_int, qs_env, mo_transpose=.TRUE.)

         !Create the symmetric kernel matrix and redistribute work_mat into it
         CALL dbcsr_create(ex_ker_offdiag, name="OFFDIAG EX KERNEL", matrix_type="S", dist=dist, &
                           row_blk_size=blk_size, col_blk_size=blk_size)
         CALL dbcsr_complete_redistribute(work_mat, ex_ker_offdiag)

         ! some clean-up
         CALL dbcsr_release(work_mat)
         DO imo = 1, ndo_mo
            CALL dbcsr_release(rhs_int(imo)%matrix)
            DEALLOCATE (rhs_int(imo)%matrix)
         END DO
         DEALLOCATE (rhs_int)

      END IF ! do_offdiag

!  The on-diagonal exchange : (IJ|ab) = (IJ|P) (P|Q)^-1 (Q|ab)
!  Already have (aJ|P) (P|Q)^-1 stored in lhs_int, need to further contract
!  Contract by sm_fm multiplication, need to transpose (aJ|Q) first
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(lhs_int(1)%matrix, nfullcols_total=s2, col_blk_size=col_blk_size)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=s2, ncol_global=ndo_mo)
      CALL cp_fm_create(contr2_int, fm_struct)

      ALLOCATE (coeffs(col_blk_size(ri_atom), ndo_mo))

      nblk = SIZE(col_blk_size)
      ALLOCATE (o3c_vec(nblk), vec_blk_size(nblk))
      vec_blk_size = 1; vec_blk_size(ri_atom) = col_blk_size(ri_atom)
      CALL o3c_vec_create(o3c_vec, vec_blk_size)

      ! a and b need to overlap for non-zero integral => same block structure as S
      ALLOCATE (abIJ, abIJ_desymm)
      CALL dbcsr_copy(abIJ, matrix_s(1)%matrix, name="(ab|IJ)")

      CALL dbcsr_create(work_mat, name="WORK", matrix_type="N", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)

      ! Loop over MOs. End matrix is symmetric => span upper half
      DO imo = 1, ndo_mo

         ! contract the second time
         CALL dbcsr_transposed(lhs_t, lhs_int(imo)%matrix)
         CALL cp_dbcsr_sm_fm_multiply(lhs_t, donor_state%gs_coeffs, contr2_int, ncol=ndo_mo)
         sr = 1
         IF (ri_atom > 1) sr = SUM(col_blk_size(1:ri_atom-1))+1
         er = SUM(col_blk_size(1:ri_atom))

         CALL cp_fm_get_submatrix(contr2_int, coeffs, start_row=sr, start_col=1, n_rows=er-sr+1, &
                                  n_cols=ndo_mo)

         DO jmo = imo, ndo_mo

            ! compute (ab|IJ) = sum (ab|P) * (P|Q)^-1 * (Q|IJ)
            CALL get_o3c_vec(o3c_vec, ri_atom, v)
            v(:) = coeffs(:, jmo)
            CALL dbcsr_set(abIJ, 0.0_dp)

            CALL contract3_o3c(xas_tdp_env%ri_o3c_ex, o3c_vec, abIJ)

            ! (ab|IJ) is symmetric, but need it as normal for distributions to match
            CALL dbcsr_desymmetrize(abIJ, abIJ_desymm)

            CALL dbcsr_iterator_start(iter, abIJ_desymm)
            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
               IF (imo == jmo .AND. jblk < iblk) CYCLE

               CALL dbcsr_get_block_p(abIJ_desymm, iblk, jblk, pblock, found)

               IF (found) THEN
                  CALL dbcsr_put_block(work_mat, (imo-1)*nblk+iblk, (jmo-1)*nblk+jblk, pblock)
               END IF

            END DO ! dbcsr iterator
            CALL dbcsr_iterator_stop(iter)
            CALL dbcsr_release(abIJ_desymm)

         END DO !jmo
         CALL dbcsr_release(lhs_t)
      END DO !imo

      CALL dbcsr_finalize(work_mat)
      CALL dbcsr_create(ex_ker_ondiag, name="ONDIAG EX KERNEL", matrix_type="S", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)
      CALL dbcsr_complete_redistribute(work_mat, ex_ker_ondiag)

!  Clean-up
      CALL cp_fm_release(contr2_int)
      CALL cp_fm_struct_release(fm_struct)
      CALL dbcsr_release(abIJ)
      CALL dbcsr_release(work_mat)
      CALL o3c_vec_release(o3c_vec)
      DO imo = 1, ndo_mo
         CALL dbcsr_release(lhs_int(imo)%matrix)
         DEALLOCATE (lhs_int(imo)%matrix)
      END DO
      DEALLOCATE (lhs_int, o3c_vec, abIJ, abIJ_desymm, work_mat, vec_blk_size)

      CALL timestop(handle)

   END SUBROUTINE kernel_exchange

! **************************************************************************************************
!> \brief Multiply all the blocks of a contracted RI integral (aI|P) by a matrix of type (P|...|Q)
!> \param contr_int the integral array
!> \param PQ the smaller matrix to multiply all blocks
!> \note  It is assumed that all non-zero blocks have the same number of columns
! **************************************************************************************************
   SUBROUTINE ri_all_blocks_mm(contr_int, PQ)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: contr_int
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: PQ

      CHARACTER(len=*), PARAMETER :: routineN = 'ri_all_blocks_mm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, iblk, imo, jblk, ndo_mo, s1, s2
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: work
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock
      TYPE(dbcsr_iterator_type)                          :: iter

      NULLIFY (pblock)

      ndo_mo = SIZE(contr_int)

      DO imo = 1, ndo_mo
         CALL dbcsr_iterator_start(iter, contr_int(imo)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
            CALL dbcsr_get_block_p(contr_int(imo)%matrix, iblk, jblk, pblock, found)

            IF (found) THEN
               s1 = SIZE(pblock, 1)
               s2 = SIZE(pblock, 2)
               ALLOCATE (work(s1, s2))
               CALL dgemm('N', 'N', s1, s2, s2, 1.0_dp, pblock, s1, PQ, s2, 0.0_dp, work, s1)
               CALL dcopy(s1*s2, work, 1, pblock, 1)
               DEALLOCATE (work)
            END IF

         END DO ! dbcsr iterator
         CALL dbcsr_iterator_stop(iter)
      END DO !imo

   END SUBROUTINE ri_all_blocks_mm

! **************************************************************************************************
!> \brief Takes the product of contracted integrals and put them in a kernel matrix
!> \param kernel the matrix where the products are stored
!> \param lhs_int the left-hand side contracted integrals
!> \param rhs_int the right-hand side contracted integrals
!> \param qs_env ...
!> \param mo_transpose whether the MO blocks should be transpose, i.e. (aI|Jb) => (aJ|Ib)
!> \Note It is assumed that the kerenl matrix is not symmetric
! **************************************************************************************************
   SUBROUTINE ri_int_product(kernel, lhs_int, rhs_int, qs_env, mo_transpose)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: kernel
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: lhs_int, rhs_int
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: mo_transpose

      CHARACTER(len=*), PARAMETER :: routineN = 'ri_int_product', routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, i, iblk, imo, j, jblk, jmo, nblk, &
                                                            ndo_mo
      LOGICAL                                            :: found, my_mt
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: prod

      NULLIFY (matrix_s, pblock)

!  Initialization
      ndo_mo = SIZE(lhs_int)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, natom=nblk)
      CALL dbcsr_create(prod, template=matrix_s(1)%matrix, matrix_type="N")
      my_mt = .FALSE.
      IF (PRESENT(mo_transpose)) my_mt = mo_transpose

      ! The kernel matrix is symmetric (even if normal type) => only fill upper half
      DO imo = 1, ndo_mo
         DO jmo = imo, ndo_mo

            i = imo; j = jmo; 
            IF (my_mt) THEN
               i = jmo; j = imo; 
            END IF

            ! Take the product lhs*rhs^T
            CALL dbcsr_multiply('N', 'T', 1.0_dp, lhs_int(i)%matrix, rhs_int(j)%matrix, &
                                0.0_dp, prod)

            ! Loop over blocks of prod and fill kernel matrix => ok cuz same (but replicated) dist
            CALL dbcsr_iterator_start(iter, prod)
            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
               IF (imo == jmo .AND. jblk < iblk) CYCLE

               CALL dbcsr_get_block_p(prod, iblk, jblk, pblock, found)

               IF (found) THEN
                  CALL dbcsr_put_block(kernel, (imo-1)*nblk+iblk, (jmo-1)*nblk+jblk, pblock)
               END IF

            END DO ! dbcsr iterator
            CALL dbcsr_iterator_stop(iter)

         END DO !jmo
      END DO !imo
      CALL dbcsr_finalize(kernel)

!  Clean-up
      CALL dbcsr_release(prod)

   END SUBROUTINE ri_int_product

! **************************************************************************************************
!> \brief Builds a neighbor lists set taylored for 3-center integral within XAS TDP, such that only
!>        excited atoms are taken into account for the list_c
!> \param ac_list the neighbor list ready for 3-center integrals
!> \param basis_a basis set list for atom a
!> \param basis_c basis set list for atom c
!> \param op_type to indicate whther the list should be built with overlap or Coulomb in mind
!> \param excited_atoms the indices of the excited atoms to consider
!> \param excited_kinds the indices of the excited kinds to consider
!> \param qs_env ...
!> \note Based on setup_neighbor_list with added features
! **************************************************************************************************
   SUBROUTINE build_xas_tdp_3c_nl(ac_list, basis_a, basis_c, op_type, excited_atoms, &
                                  excited_kinds, qs_env)

      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ac_list
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_a, basis_c
      CHARACTER(len=*), INTENT(IN)                       :: op_type
      INTEGER, DIMENSION(:), INTENT(IN)                  :: excited_atoms, excited_kinds
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_xas_tdp_3c_nl', &
         routineP = moduleN//':'//routineN

      CHARACTER                                          :: my_op
      INTEGER                                            :: ikind, nkind
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: a_present, c_present
      REAL(dp)                                           :: subcells
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: a_radius, c_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (atomic_kind_set, distribution_1d, distribution_2d, molecule_set, particle_set, cell)

!  By default, take overlap operator
      my_op = "O"
      IF (op_type == "COULOMB ") my_op = "C"

!  Initialization
      CALL get_qs_env(qs_env, nkind=nkind)
      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)

      ALLOCATE (a_present(nkind), c_present(nkind))
      a_present = .FALSE.
      c_present = .FALSE.
      ALLOCATE (a_radius(nkind), c_radius(nkind))
      a_radius = 0.0_dp
      c_radius = 0.0_dp

!  Set up the radii
      DO ikind = 1, nkind
         !orbital basis set
         IF (ASSOCIATED(basis_a(ikind)%gto_basis_set)) THEN
            a_present(ikind) = .TRUE.
            CALL get_gto_basis_set(basis_a(ikind)%gto_basis_set, kind_radius=a_radius(ikind))
         END IF
         !RI_XAS basis set, take only excited kinds
         IF (ASSOCIATED(basis_c(ikind)%gto_basis_set) .AND. ANY(excited_kinds == ikind)) THEN
            c_present(ikind) = .TRUE.
            ! If operator is coulomb, then center c interacts with everything
            IF (my_op == "O") THEN
               CALL get_gto_basis_set(basis_c(ikind)%gto_basis_set, kind_radius=c_radius(ikind))
            ELSE
               c_radius(ikind) = 1000000.0_dp
            END IF
         END IF
      END DO !ikind

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

      CALL build_neighbor_lists(ac_list, particle_set, atom2d, cell, pair_radius, subcells, &
                                operator_type="ABC", atomb_to_keep=excited_atoms, &
                                nlname="XAS_TDP_3c_nl")

!  Clean-up
      CALL atom2d_cleanup(atom2d)

   END SUBROUTINE build_xas_tdp_3c_nl

! **************************************************************************************************
!> \brief Includes the SOC effects on the precomputed singlet and triplet excitations
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note Using AMEWs, build an hermitian matrix with all excited states SOC coupling + the
!>       excitation energies on the diagonal. Then diagonalize it to get the new excitation
!>       energies and corresponding linear combinations of lr_coeffs.
!>       The AMEWs are normalized
! **************************************************************************************************
   SUBROUTINE include_full_soc(donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'include_full_soc', &
         routineP = moduleN//':'//routineN

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: transdip
      INTEGER                                            :: handle, homo, i, isg, itp, jtp, nao, &
                                                            ndo_mo, nsg, ntot, ntp
      REAL(dp)                                           :: gs_sum_x, gs_sum_y, gs_sum_z, &
                                                            s_tt_sumdiag, soc_gst, soc_st, soc_tt, &
                                                            sqrt2, trace
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag, gs_diag, tmp_evals
      REAL(dp), DIMENSION(:), POINTER                    :: sg_evals, tp_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type), POINTER                         :: dip_cfm, evecs_cfm, hami_cfm, work_cfm
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: amew_dip, sg_coeffs, tp_coeffs
      TYPE(cp_fm_struct_type), POINTER                   :: dip_struct, full_struct, gsgs_struct, &
                                                            prod_struct, std_struct, vec_struct
      TYPE(cp_fm_type), POINTER :: domo_soc_x, domo_soc_y, domo_soc_z, gs_coeffs, gs_fm, img_fm, &
         mo_coeff, prod_fm, real_fm, s_st, s_tt, vec_s, vec_soc_x, vec_soc_y, vec_soc_z, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: orb_soc_x, orb_soc_y, orb_soc_z
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      NULLIFY (sg_coeffs, tp_coeffs, gs_coeffs, sg_evals, tp_evals, real_fm, img_fm, full_struct)
      NULLIFY (para_env, blacs_env, vec_s, prod_fm, prod_struct, vec_struct, domo_soc_z, work_cfm)
      NULLIFY (matrix_s, orb_soc_x, hami_cfm, evecs_cfm, s_st, s_tt, domo_soc_x, domo_soc_y)
      NULLIFY (orb_soc_y, orb_soc_z, mo_coeff, mos, std_struct, work_fm, gsgs_struct, amew_dip)
      NULLIFY (vec_soc_x, vec_soc_y, vec_soc_z)

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
      ndo_mo = SIZE(donor_state%mo_indices)
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      orb_soc_x => xas_tdp_env%orb_soc(1)%matrix
      orb_soc_y => xas_tdp_env%orb_soc(2)%matrix
      orb_soc_z => xas_tdp_env%orb_soc(3)%matrix

!  Creating the real part and imaginary part of the final SOC fm
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env, mos=mos)
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      CALL cp_fm_create(real_fm, full_struct)
      CALL cp_fm_create(img_fm, full_struct)
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao, homo=homo)

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
      CALL cp_fm_create(s_tt, prod_struct)
      CALL cp_fm_create(s_st, prod_struct)
      CALL cp_fm_create(vec_s, vec_struct)
      CALL cp_fm_create(vec_soc_x, vec_struct)
      CALL cp_fm_create(vec_soc_y, vec_struct)
      CALL cp_fm_create(vec_soc_z, vec_struct)

!  Some other helper arrays
      ALLOCATE (diag(ndo_mo))

!  Precompute everything we can before looping over excited states
      sqrt2 = SQRT(2.0_dp)

      ! The sum of the ground state KS orbital contributions: sum_i <phi_i^0|Hsoc|phi_i^0>
      ALLOCATE (gs_diag(homo))
      CALL cp_fm_struct_create(gsgs_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=homo, ncol_global=homo)
      CALL cp_fm_get_info(mo_coeff, matrix_struct=std_struct)
      CALL cp_fm_create(gs_fm, gsgs_struct)
      CALL cp_fm_create(work_fm, std_struct)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_x, mo_coeff, work_fm, ncol=homo)
      CALL cp_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, work_fm, 0.0_dp, gs_fm)
      CALL cp_fm_get_diag(gs_fm, gs_diag)
      gs_sum_x = SUM(gs_diag)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, mo_coeff, work_fm, ncol=homo)
      CALL cp_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, work_fm, 0.0_dp, gs_fm)
      CALL cp_fm_get_diag(gs_fm, gs_diag)
      gs_sum_y = SUM(gs_diag)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, mo_coeff, work_fm, ncol=homo)
      CALL cp_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, work_fm, 0.0_dp, gs_fm)
      CALL cp_fm_get_diag(gs_fm, gs_diag)
      gs_sum_z = SUM(gs_diag)

      CALL cp_fm_release(gs_fm)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_struct_release(gsgs_struct)
      DEALLOCATE (gs_diag)

      ! The subset of the donor MOs matrix elements: <phi_I^0|Hsoc|phi_J^0>
      CALL cp_fm_create(domo_soc_x, prod_struct)
      CALL cp_fm_create(domo_soc_y, prod_struct)
      CALL cp_fm_create(domo_soc_z, prod_struct)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_x, gs_coeffs, vec_soc_x, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_x, 0.0_dp, domo_soc_x)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, gs_coeffs, vec_soc_y, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_y, 0.0_dp, domo_soc_y)

      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, gs_coeffs, vec_soc_z, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_z, 0.0_dp, domo_soc_z)

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
         CALL cp_fm_set_element(img_fm, 1, 1+nsg+jtp, soc_gst) ! <0|H_x|T^-1>
         CALL cp_fm_set_element(img_fm, 1, 1+nsg+2*ntp+jtp, -1.0_dp*soc_gst) ! <0|H_x|T^+1>

         ! real part
         CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, tp_coeffs(jtp)%matrix, vec_soc_y, ncol=ndo_mo)
         CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_soc_y, 0.0_dp, prod_fm)
         CALL cp_fm_get_diag(prod_fm, diag)
         soc_gst = SUM(diag)
         CALL cp_fm_set_element(real_fm, 1, 1+nsg+jtp, soc_gst) ! <0|H_y|T^-1>
         CALL cp_fm_set_element(real_fm, 1, 1+nsg+2*ntp+jtp, soc_gst) ! <0|H_y|T^+1>

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
                         0.0_dp, s_st)

            ! compute the singlet-triplet (Ms=+-1) SOC
            ! imaginary part
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, sg_coeffs(isg)%matrix, vec_soc_x, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)
            CALL cp_fm_trace(domo_soc_x, s_st, trace)

            soc_st = 0.5_dp*sqrt2*(SUM(diag)-trace)

            CALL cp_fm_set_element(img_fm, 1+isg, 1+nsg+jtp, soc_st) ! <S|H_x|T^-1>
            CALL cp_fm_set_element(img_fm, 1+isg, 1+nsg+2*ntp+jtp, -1.0_dp*soc_st) !<S|H_x|T^+1>

            ! real part
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, sg_coeffs(isg)%matrix, vec_soc_y, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)
            CALL cp_fm_trace(domo_soc_y, s_st, trace)

            soc_st = 0.5_dp*sqrt2*(SUM(diag)-trace)

            CALL cp_fm_set_element(real_fm, 1+isg, 1+nsg+jtp, soc_st) ! <S|H_x|T^-1>
            CALL cp_fm_set_element(real_fm, 1+isg, 1+nsg+2*ntp+jtp, soc_st) ! <S|H_x|T^+1>

            ! compute the singlet-triplet (Ms=0) SOC (purely imaginary)
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, sg_coeffs(isg)%matrix, vec_soc_z, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)
            CALL cp_fm_trace(domo_soc_z, s_st, trace)

            soc_st = SUM(diag)-trace

            CALL cp_fm_set_element(img_fm, 1+isg, 1+nsg+ntp+jtp, soc_st) ! <S|H_x|T^0>

         END DO !isg

         ! Loop over the triplet to do the triplet-triplet SOC
         DO itp = 1, ntp

            ! compute the triplet-triplet overlap
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, tp_coeffs(itp)%matrix, vec_s, &
                         0.0_dp, s_tt)
            CALL cp_fm_get_diag(s_tt, diag)
            s_tt_sumdiag = SUM(diag)

            ! The Ms=0 to Ms=+-1 triplet-triplet SOC
            ! imaginary part
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, tp_coeffs(itp)%matrix, vec_soc_x, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)
            CALL cp_fm_trace(domo_soc_x, s_tt, trace)

            soc_tt = 0.5_dp*sqrt2*(SUM(diag)-trace+2.0_dp*gs_sum_x*s_tt_sumdiag)

            CALL cp_fm_set_element(img_fm, 1+nsg+jtp, 1+nsg+ntp+itp, -1.0_dp*soc_tt) ! <T^-1|H_x|T^0>
            CALL cp_fm_set_element(img_fm, 1+nsg+ntp+itp, 1+nsg+2*ntp+jtp, soc_tt) ! <T^0|H_x|T^+1>

            ! real part
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, tp_coeffs(itp)%matrix, vec_soc_y, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)
            CALL cp_fm_trace(domo_soc_y, s_tt, trace)

            soc_tt = 0.5_dp*sqrt2*(SUM(diag)-trace+2.0_dp*gs_sum_y*s_tt_sumdiag)

            CALL cp_fm_set_element(real_fm, 1+nsg+jtp, 1+nsg+ntp+itp, soc_tt) ! <T^-1|H_x|T^0>
            CALL cp_fm_set_element(real_fm, 1+nsg+ntp+itp, 1+nsg+2*ntp+jtp, -1.0_dp*soc_tt) ! <T^0|H_x|T^+1>

            ! Do the Ms=+1-Ms=+1 and Ms=-1-Ms=-1 coupling
            IF (itp >= jtp) CYCLE ! SOC purely imaginary, must be zero on the diagonal

            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, tp_coeffs(itp)%matrix, vec_soc_z, &
                         0.0_dp, prod_fm)
            CALL cp_fm_get_diag(prod_fm, diag)
            CALL cp_fm_trace(domo_soc_z, s_tt, trace)

            soc_tt = SUM(diag)-trace+2.0_dp*gs_sum_z*s_tt_sumdiag

            CALL cp_fm_set_element(img_fm, 1+nsg+itp, 1+nsg+jtp, -1.0_dp*soc_tt) ! <T^-1|H_x|T^-1>
            CALL cp_fm_set_element(img_fm, 1+nsg+2*ntp+itp, 1+nsg+2*ntp+jtp, soc_tt) !<T^+1|H_x|T^+1>

         END DO !itp
      END DO !jtp

!  Intermediate clean-up
      CALL cp_fm_struct_release(prod_struct)
      CALL cp_fm_release(prod_fm)
      CALL cp_fm_release(s_tt)
      CALL cp_fm_release(s_st)
      CALL cp_fm_release(domo_soc_x)
      CALL cp_fm_release(domo_soc_y)
      CALL cp_fm_release(domo_soc_z)
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
      CALL get_amew_op(amew_dip, xas_tdp_env%dipmat, donor_state, qs_env)

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

   END SUBROUTINE include_full_soc

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
! **************************************************************************************************
   SUBROUTINE get_amew_op(amew_op, ao_op, donor_state, qs_env)

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: amew_op
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ao_op
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'get_amew_op', routineP = moduleN//':'//routineN

      INTEGER                                            :: homo, i, isg, itp, jsg, jtp, nao, &
                                                            ndo_mo, nsg, ntot, ntp
      REAL(dp)                                           :: gsgs_op, op, sqrt2, trace
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag, gs_diag
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: sg_coeffs, tp_coeffs
      TYPE(cp_fm_struct_type), POINTER                   :: full_struct, gsgs_struct, prod_struct, &
                                                            std_struct, vec_struct
      TYPE(cp_fm_type), POINTER                          :: amew_op_i, domo_op, gs_coeffs, gs_fm, &
                                                            mo_coeff, prod_fm, vec_op, vec_s, &
                                                            work, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: ao_op_i
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      NULLIFY (gs_coeffs, sg_coeffs, tp_coeffs, matrix_s, full_struct, prod_fm, prod_struct, work_fm)
      NULLIFY (vec_struct, blacs_env, para_env, mo_coeff, mos, gsgs_struct, std_struct, domo_op)
      NULLIFY (vec_op, vec_s, gs_fm)

!  Initialization
      gs_coeffs => donor_state%gs_coeffs
      sg_coeffs => donor_state%sg_coeffs
      tp_coeffs => donor_state%tp_coeffs
      nsg = SIZE(sg_coeffs)
      ntp = SIZE(tp_coeffs)
      ntot = 1+nsg+3*ntp
      ndo_mo = SIZE(donor_state%mo_indices)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, para_env=para_env, blacs_env=blacs_env, mos=mos)
      sqrt2 = SQRT(2.0_dp)

!  Create the amew_op matrix
      CALL cp_fm_struct_create(full_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ntot, ncol_global=ntot)
      ALLOCATE (amew_op(SIZE(ao_op)))
      DO i = 1, SIZE(ao_op)
         CALL cp_fm_create(amew_op(i)%matrix, full_struct)
      END DO !i

!  Create the work and helper fms
      CALL cp_fm_get_info(gs_coeffs, matrix_struct=vec_struct)
      CALL cp_fm_struct_create(prod_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_mo, ncol_global=ndo_mo)
      CALL cp_fm_create(prod_fm, prod_struct)
      CALL cp_fm_create(domo_op, prod_struct)
      CALL cp_fm_create(vec_op, vec_struct)
      CALL cp_fm_create(vec_s, vec_struct)
      CALL cp_fm_create(work, full_struct)
      ALLOCATE (diag(ndo_mo))
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao, homo=homo)

! Iterate over the dimensions of the operator
      DO i = 1, SIZE(ao_op)

         ao_op_i => ao_op(i)%matrix
         amew_op_i => amew_op(i)%matrix

         ! compute the ground-state/ground-state matrix element <0|0> = 2*sum_j <phi_j|op|phi_j>
         ALLOCATE (gs_diag(homo))
         CALL cp_fm_struct_create(gsgs_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=homo, ncol_global=homo)
         CALL cp_fm_get_info(mo_coeff, matrix_struct=std_struct)
         CALL cp_fm_create(gs_fm, gsgs_struct)
         CALL cp_fm_create(work_fm, std_struct)

         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, mo_coeff, work_fm, ncol=homo)
         CALL cp_gemm('T', 'N', homo, homo, nao, 1.0_dp, mo_coeff, work_fm, 0.0_dp, gs_fm)
         CALL cp_fm_get_diag(gs_fm, gs_diag)
         gsgs_op = 2.0_dp*SUM(gs_diag)
         CALL cp_fm_set_element(amew_op_i, 1, 1, gsgs_op)

         CALL cp_fm_release(gs_fm)
         CALL cp_fm_release(work_fm)
         CALL cp_fm_struct_release(gsgs_struct)
         DEALLOCATE (gs_diag)

         ! Compute the operator for the donor MOs
         CALL cp_dbcsr_sm_fm_multiply(ao_op_i, gs_coeffs, vec_op, ncol=ndo_mo)
         CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, gs_coeffs, vec_op, 0.0_dp, domo_op)

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
               CALL cp_fm_get_diag(prod_fm, diag)
               op = op+SUM(diag)*gsgs_op

               CALL cp_fm_trace(prod_fm, domo_op, trace)
               op = op-trace

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
               CALL cp_fm_get_diag(prod_fm, diag)
               op = op+SUM(diag)*gsgs_op

               CALL cp_fm_trace(prod_fm, domo_op, trace)
               op = op-trace

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
      CALL cp_fm_release(domo_op)
      CALL cp_fm_struct_release(prod_struct)
      CALL cp_fm_struct_release(full_struct)

   END SUBROUTINE get_amew_op

END MODULE xas_tdp_utils

