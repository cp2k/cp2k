!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for X-ray absorption spectroscopy using TDDFPT
!> \author AB (01.2018)
! **************************************************************************************************

MODULE xas_tdp_utils
   USE ai_contraction_sphi,             ONLY: abc_contract
   USE ai_coulomb,                      ONLY: coulomb3
   USE arnoldi_api,                     ONLY: arnoldi_data_type,&
                                              arnoldi_ev,&
                                              deallocate_arnoldi_data,&
                                              get_selected_ritz_val,&
                                              setup_arnoldi_data
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_trace,&
                                              cp_fm_triangular_invert,&
                                              cp_fm_triangular_multiply,&
                                              cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: cp_fm_syevx
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, &
        cp_fm_release, cp_fm_set_element, cp_fm_to_fm, cp_fm_to_fm_submat, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, &
        dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_type
   USE generic_os_integrals,            ONLY: int_operators_r12_ab_os
   USE input_constants,                 ONLY: xas_dip_len, xas_dip_vel
   USE iterate_matrix,                  ONLY: invert_Hotelling, matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat_symm
   USE orbital_pointers,                ONLY: ncoset
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE util,                            ONLY: locate
   USE xas_tdp_types,                   ONLY: donor_state_type,&
                                              xas_tdp_control_type,&
                                              xas_tdp_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_utils'

   PUBLIC :: setup_xas_tdp_full_prob, solve_xas_tdp_full_prob, compute_dipole_fosc, &
             compute_ri_coulomb3_int, compute_ri_coulomb2_int, include_full_soc

CONTAINS

! **************************************************************************************************
!> \brief Builds the matrix that defines the XAS TDDFPT generalized eigenvalue problem to be solved
!>        for excitation energies omega. The problem has the form omega*C = M*C, where C contains
!>        the reponse orbitals coefficients. The matrix is stored in the given donor_state
!> \param donor_state the donor_state for which the problem is restricted
!> \param qs_env ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \note the matrix M has the form | -G^-1*M_d   -G^-1*M_o |
!>                                 |  G^-1*M_o    G^-1*M_d |,
!>       depending on whther we consider singlet or triplet excitation, the diagonal (M_d) and
!>       off-diagonal (M_o) parts of M differ:
!>       - For singlet: M_d = A + 2B + C_aa + C_ab - D
!>                      M_o = 2B + C_aa + C_ab - E
!>       - For triplet: M_d = A + C_aa - C_ab - D
!>                      M_o = C_aa - C_ab - E
!>       where other subroutines computes the matrices A, B, E, D and G, which are:
!>       - A: the ground-state contribution F_ab*delta_IJ - epsilon_IJ*S_ab
!>       - B: the Coulob kernel ~(aI|Jb)
!>       - C: the xc kernel c_aa (double derivatibe wrt to n_alpha) and C_ab (wrt n_alpha and n_beta)
!>       - D: the on-digonal exact exchange kernel ~(ab|IJ)
!>       - E: the off-diagonal exact exchange kernel ~(aJ|Ib)
!>       - G: the metric  S_ab*delta_IJ
!>       For the xc functionals, C_aa + C_ab or C_aa - C_ab are stored in the same matrix
!>       In the above definitions, I,J label the donnor MOs and a,b the sgfs of the basis
!>       In the case of non-TDA, one turns the problem Hermitian
! **************************************************************************************************
   SUBROUTINE setup_xas_tdp_full_prob(donor_state, qs_env, xas_tdp_env, xas_tdp_control)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_xas_tdp_full_prob', &
                                     routineP = moduleN//":"//routineN
      
      INTEGER                                         :: handle, nblk_sub
      TYPE(dbcsr_type), POINTER                       :: sg_matrix_tdp, tp_matrix_tdp
      TYPE(dbcsr_type), POINTER                       :: matrix_a, matrix_b, matrix_d, matrix_e, &
                                                         tmp_mat, matrix_c, proj_Q, work
      TYPE(dbcsr_distribution_type), POINTER          :: submat_dist
      INTEGER, DIMENSION(:), POINTER                  :: submat_blk_size, col_dist_sub, &
                                                         row_dist_sub
      INTEGER, DIMENSION(:,:), POINTER                :: pgrid
      REAL(dp)                                        :: min_ev, threshold, sx
      TYPE(arnoldi_data_type)                         :: my_arnoldi
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrices
      LOGICAL                                         :: conv, do_xc, do_hfx

      CALL timeset(routineN, handle)

      NULLIFY (matrix_a, submat_dist, submat_blk_size, matrices)
      NULLIFY (sg_matrix_tdp, tp_matrix_tdp, matrix_b, matrix_d, matrix_e, &
               tmp_mat, matrix_c)

      CPASSERT(ASSOCIATED(donor_state))
      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(xas_tdp_env))
      CPASSERT(ASSOCIATED(xas_tdp_control))
      IF (xas_tdp_control%do_singlet) THEN
         ALLOCATE (donor_state%sg_matrix_tdp)
         sg_matrix_tdp => donor_state%sg_matrix_tdp
      END IF
      IF (xas_tdp_control%do_triplet) THEN
         ALLOCATE (donor_state%tp_matrix_tdp)
         tp_matrix_tdp => donor_state%tp_matrix_tdp
      END IF

      do_xc = xas_tdp_control%do_xc
      do_hfx = xas_tdp_control%do_hfx
      sx = xas_tdp_control%sx

!  Get the submatrices distribution and block sizes
      CALL compute_submat_dist_and_blk_size(submat_dist, submat_blk_size, donor_state, qs_env)
      nblk_sub = SIZE(submat_blk_size)

!  Get the projector on the unperturbed unoccupied GS in the correct format
      ALLOCATE (proj_Q, work)
      CALL get_q_projector(proj_Q, submat_dist, submat_blk_size, donor_state, xas_tdp_env)
      CALL dbcsr_create(matrix=work, matrix_type="N", dist=submat_dist, name="WORK", &
                        row_blk_size=submat_blk_size, col_blk_size=submat_blk_size)

!  IF TDA, simplified problem:
      IF (xas_tdp_control%tamm_dancoff) THEN

!     Compute the submatrices, store the metric and ignore the offdiag exchange kernel
!     Only compute needed forms of the kernel according to the functional
         ALLOCATE (matrix_a)
         ALLOCATE (matrix_b)
         ALLOCATE (matrix_c)
         ALLOCATE (matrix_d)
         ALLOCATE (donor_state%metric(1))

         CALL build_gs_contribution(matrix_a, submat_dist, submat_blk_size, donor_state, qs_env)
         IF (do_hfx) THEN
            CALL build_kernel_exchange(ex_ker_diag=matrix_d, dist=submat_dist, blk_size=submat_blk_size, &
                                       donor_state=donor_state, xas_tdp_env=xas_tdp_env, qs_env=qs_env)
         END IF
         CALL build_metric(donor_state%metric, submat_dist, submat_blk_size, donor_state, qs_env)

!     Create the simplfied TDA matrix_tdp. The xc kernel changes with the excited state multiplicity
         IF (xas_tdp_control%do_singlet) THEN
            CALL dbcsr_copy(sg_matrix_tdp, matrix_a, name="SINGLET MATRIX TDP")

            ! Coulomb kernel only needed for singlets
            CALL build_kernel_coulomb(matrix_b, submat_dist, submat_blk_size, donor_state, xas_tdp_env, qs_env)

            CALL dbcsr_add(sg_matrix_tdp, matrix_b, 1.0_dp, 2.0_dp) 
            IF (do_hfx) CALL dbcsr_add(sg_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx) !scaled hfx
            IF (do_xc) THEN
               CALL build_kernel_xc(matrix_c, submat_dist, submat_blk_size, donor_state, &
                                    xas_tdp_env, qs_env, singlet=.TRUE.)
               CALL dbcsr_add(sg_matrix_tdp, matrix_c, 1.0_dp, 1.0_dp) ! xc kernel
            END IF ! do_xc

            ! Take the product with the Q projector:
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, sg_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, sg_matrix_tdp)

         END IF

         IF (xas_tdp_control%do_triplet) THEN

            CALL dbcsr_copy(tp_matrix_tdp, matrix_a, name="TRIPLET MATRIX TDP")
            IF (do_hfx) CALL dbcsr_add(tp_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx) !scaled hfx
            IF (do_xc) THEN
               CALL build_kernel_xc(matrix_c, submat_dist, submat_blk_size, donor_state, &
                                    xas_tdp_env, qs_env, singlet=.FALSE.)
               CALL dbcsr_add(tp_matrix_tdp, matrix_c, 1.0_dp, 1.0_dp) ! xc_kernel
            END IF ! do_xc

            ! Take the product with the Q projector:
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, tp_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, tp_matrix_tdp)

         END IF

!     Some clean-up
         CALL dbcsr_release(matrix_a)
         CALL dbcsr_release(matrix_b)
         CALL dbcsr_release(matrix_c)
         CALL dbcsr_release(matrix_d)
         DEALLOCATE (matrix_a, matrix_b, matrix_c, matrix_d)

      ELSE IF (.NOT. xas_tdp_control%tamm_dancoff) THEN

!     Compute the submatrices A, G and the different kernels (B, C, D and E)
         ALLOCATE (tmp_mat)
         ALLOCATE (matrix_a)
         ALLOCATE (matrix_b)
         ALLOCATE (matrix_c)
         ALLOCATE (matrix_d)
         ALLOCATE (matrix_e)
         ALLOCATE (donor_state%metric(2))
         CALL build_gs_contribution(matrix_a, submat_dist, submat_blk_size, donor_state, qs_env)
         IF (do_hfx) THEN
            CALL build_kernel_exchange(matrix_d, matrix_e, submat_dist, submat_blk_size, &
                                       donor_state, xas_tdp_env, qs_env)
         END IF
         CALL build_metric(donor_state%metric, submat_dist, submat_blk_size, donor_state, qs_env, do_inv=.TRUE.)

!     Build the auxiliary matrix (A-D+E)^0.5
         CALL dbcsr_copy(tmp_mat, matrix_a)
         IF (do_hfx) THEN
            CALL dbcsr_add(tmp_mat, matrix_d, 1.0_dp, -1.0_dp*sx) !scaled hfx
            CALL dbcsr_add(tmp_mat, matrix_e, 1.0_dp, 1.0_dp*sx)
         END IF

!     Take the product with the Q projector
         CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, tmp_mat, 0.0_dp, work)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, tmp_mat)

!     Warn the user if the matrix is not positive semi-definite and the procedure not sound
         threshold = 1.0E-8_dp ! kind of arbitrary TODO: should we ask the user for that ?
         ALLOCATE (matrices(1))
         matrices(1)%matrix => tmp_mat
         CALL setup_arnoldi_data(my_arnoldi, matrices, max_iter=200, threshold=threshold, selection_crit=3, &
                                 nval_request=1, nrestarts=30, generalized_ev=.FALSE., iram=.FALSE.)
         CALL arnoldi_ev(matrices, my_arnoldi)
         min_ev = REAL(get_selected_ritz_val(my_arnoldi, 1), dp)
         IF (min_ev < 0.0_dp .AND. ABS(min_ev) > threshold) THEN
            CPWARN("The full TDDFT problem might not have been soundly turned Hermitian. Try TDA.")
         END IF

         ALLOCATE (donor_state%matrix_aux(2))
         ALLOCATE (donor_state%matrix_aux(1)%matrix, donor_state%matrix_aux(2)%matrix)
         CALL dbcsr_create(matrix=donor_state%matrix_aux(1)%matrix, matrix_type="S", dist=submat_dist, &
                           row_blk_size=submat_blk_size, col_blk_size=submat_blk_size, name="MAT_AUX 1")
         CALL dbcsr_create(matrix=donor_state%matrix_aux(2)%matrix, matrix_type="S", dist=submat_dist, &
                           row_blk_size=submat_blk_size, col_blk_size=submat_blk_size, name="MAT_AUX 2")

         ! at this point, might have tiny negative eigenvalues < threshold. To make sure it works,
         ! add threshold on the diag so that no more negative eV
         CALL dbcsr_add_on_diag(tmp_mat, threshold)

         CALL matrix_sqrt_Newton_Schulz(donor_state%matrix_aux(1)%matrix, donor_state%matrix_aux(2)%matrix, &
                                        tmp_mat, threshold=1.0E-10_dp, order=3, eps_lanczos=1.0E-8_dp, &
                                        max_iter_lanczos=500, converged=conv)
         CPASSERT(conv)
         CALL deallocate_arnoldi_data(my_arnoldi)
         DEALLOCATE (matrices)

!     Create the main matrix_tdp as:
!     M = (A-D+E)^0.5 * G^-1 * (A+4B+2C-E-D) * G^-1 * (A-D+E)^0.5
!     Only pass the middle sum as matrix_tdp, the product is done in the solve routine

         IF (xas_tdp_control%do_singlet) THEN
            CALL dbcsr_create(matrix=sg_matrix_tdp, name="SINGLET XAS TDP MATRIX", matrix_type="S", &
                              dist=submat_dist, row_blk_size=submat_blk_size, &
                              col_blk_size=submat_blk_size)
            CALL build_kernel_coulomb(matrix_b, submat_dist, submat_blk_size, donor_state, &
                                      xas_tdp_env, qs_env)

            CALL dbcsr_copy(sg_matrix_tdp, matrix_a)
            CALL dbcsr_add(sg_matrix_tdp, matrix_b, 1.0_dp, 4.0_dp)
            IF (do_hfx) THEN
               CALL dbcsr_add(sg_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx) !scaled hfx
               CALL dbcsr_add(sg_matrix_tdp, matrix_e, 1.0_dp, -1.0_dp*sx)
            END IF !HFX
            IF (do_xc) THEN
               CALL build_kernel_xc(matrix_c, submat_dist, submat_blk_size, &
                                    donor_state, xas_tdp_env, qs_env, singlet=.TRUE.)
               CALL dbcsr_add(sg_matrix_tdp, matrix_c, 1.0_dp, 2.0_dp)
            END IF

            ! Take the product with the Q projector
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, sg_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, sg_matrix_tdp)

         END IF ! singlet

         IF (xas_tdp_control%do_triplet) THEN
            CALL dbcsr_create(matrix=tp_matrix_tdp, name="TRIPLET XAS TDP MATRIX", matrix_type="S", &
                              dist=submat_dist, row_blk_size=submat_blk_size, &
                              col_blk_size=submat_blk_size)

            CALL dbcsr_copy(tp_matrix_tdp, matrix_a)
            IF (do_hfx) THEN
               CALL dbcsr_add(tp_matrix_tdp, matrix_d, 1.0_dp, -1.0_dp*sx) !scaled hfx
               CALL dbcsr_add(tp_matrix_tdp, matrix_e, 1.0_dp, -1.0_dp*sx)
            END IF !HFX
            IF (do_xc) THEN
               CALL build_kernel_xc(matrix_c, submat_dist, submat_blk_size, &
                                    donor_state, xas_tdp_env, qs_env, singlet=.FALSE.)
               CALL dbcsr_add(tp_matrix_tdp, matrix_c, 1.0_dp, 2.0_dp)
            END IF

            ! Take the product with the Q projector
            CALL dbcsr_multiply('N', 'N', 1.0_dp, proj_Q, tp_matrix_tdp, 0.0_dp, work)
            CALL dbcsr_multiply('N', 'T', 1.0_dp, work, proj_Q, 0.0_dp, tp_matrix_tdp)

         END IF !triplet

!     Some clean-up
         CALL dbcsr_release(matrix_a)
         CALL dbcsr_release(matrix_b)
         CALL dbcsr_release(matrix_c)
         CALL dbcsr_release(matrix_d)
         CALL dbcsr_release(matrix_e)
         CALL dbcsr_release(tmp_mat)
         DEALLOCATE (matrix_a, matrix_b, matrix_c, matrix_d, matrix_e, tmp_mat)

      END IF ! Tamm_Dancoff

!  Finalize
      IF (xas_tdp_control%do_singlet) CALL dbcsr_finalize(sg_matrix_tdp)
      IF (xas_tdp_control%do_triplet) CALL dbcsr_finalize(tp_matrix_tdp)

!  General Clean-up
      CALL dbcsr_release(proj_Q)
      CALL dbcsr_release(work)
      CALL dbcsr_distribution_release(submat_dist)
      DEALLOCATE (submat_dist, submat_blk_size, proj_Q, work)

      CALL timestop(handle)

   END SUBROUTINE setup_xas_tdp_full_prob

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
   SUBROUTINE solve_xas_tdp_full_prob(donor_state, xas_tdp_control, xas_tdp_env, qs_env, singlet)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: singlet

      CHARACTER(len=*), PARAMETER :: routineN = 'solve_xas_tdp_full_prob', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, imo, nao, ndo_mo, nelectron, &
                                                            nevals, nrow
      REAL(dp)                                           :: max_ev
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: scaling, tmp_evals
      REAL(dp), DIMENSION(:), POINTER                    :: lr_evals
      TYPE(arnoldi_data_type)                            :: my_arnoldi
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: lr_coeffs
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, gs_struct
      TYPE(cp_fm_type), POINTER                          :: c_diff, c_sum, lhs_matrix, rhs_matrix, &
                                                            work
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrices
      TYPE(dbcsr_type), POINTER                          :: filter, matrix_tdp, maux, paux, tmp_mat

      CALL timeset(routineN, handle)

      NULLIFY (para_env, blacs_env, fm_struct, rhs_matrix, matrix_tdp, lhs_matrix, work, tmp_mat)
      NULLIFY (c_diff, c_sum, gs_struct, lr_evals, lr_coeffs, filter, paux, maux)

      IF (singlet) THEN
         matrix_tdp => donor_state%sg_matrix_tdp
      ELSE
         matrix_tdp => donor_state%tp_matrix_tdp
      END IF
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env, nelectron_total=nelectron)

!     Initialization
      CALL dbcsr_get_info(matrix_tdp, nfullrows_total=nrow)
      ndo_mo = SIZE(donor_state%mo_indices)
      nevals = nrow-ndo_mo*nelectron/2
      IF (xas_tdp_control%n_excited > 0 .AND. xas_tdp_control%n_excited < nevals) THEN
         nevals = xas_tdp_control%n_excited
      END IF
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nrow, para_env=para_env, &
                               ncol_global=nrow)
      CPASSERT(ASSOCIATED(xas_tdp_env))

!     Need to go to full matrices to diagonalize
      CALL cp_fm_create(rhs_matrix, fm_struct)

!     Allocating space for solutions
      ALLOCATE (tmp_evals(nrow))
      IF (singlet) THEN
         ALLOCATE (donor_state%sg_evals(nevals))
         lr_evals => donor_state%sg_evals
      ELSE
         ALLOCATE (donor_state%tp_evals(nevals))
         lr_evals => donor_state%tp_evals
      END IF
      CALL cp_fm_create(c_sum, fm_struct)
      CALL cp_fm_create(work, fm_struct)

!     Allocate the scaling vector used for normalization
      ALLOCATE (scaling(nevals))

!     The filter to get projected to zero eigenvalues out of the way
      ALLOCATE (filter)

!     Test on TDA
      IF (xas_tdp_control%tamm_dancoff) THEN

!        Get an idea of how much we must shift the gs eigenvalues to have them out of the way
         ALLOCATE (matrices(2))
         matrices(1)%matrix => matrix_tdp
         matrices(2)%matrix => donor_state%metric(1)%matrix
         CALL setup_arnoldi_data(my_arnoldi, matrices, max_iter=50, threshold=1.0E-3_dp, selection_crit=2, &
                                 nval_request=1, nrestarts=30, generalized_ev=.TRUE., iram=.FALSE.)
         CALL arnoldi_ev(matrices, my_arnoldi)
         max_ev = REAL(get_selected_ritz_val(my_arnoldi, 1), dp)
         CALL deallocate_arnoldi_data(my_arnoldi)
         DEALLOCATE (matrices)

!        Shift away zero eigenvalues
         CALL get_gs_filter(filter, matrix_tdp, donor_state, xas_tdp_env)
         CALL dbcsr_add(matrix_tdp, filter, 1.0_dp, max_ev+10.0_dp)

!        Get the main matrix as fm
         CALL copy_dbcsr_to_fm(matrix_tdp, rhs_matrix)

!        Get the metric as a fm
         CALL cp_fm_create(lhs_matrix, fm_struct)
         CALL copy_dbcsr_to_fm(donor_state%metric(1)%matrix, lhs_matrix)

!        Diagonalisation (Cholesky decomposition)
         CALL cp_fm_geeig_syevx(rhs_matrix, lhs_matrix, c_sum, tmp_evals, nevals, work)

!        Take the first nevals eigenvalues and eigenvectors
         lr_evals(:) = tmp_evals(1:nevals)

!        TDA specific clean-up
         CALL cp_fm_release(lhs_matrix)

!        Normalize. c_sum is just c^- for TDA
         scaling = 0.0_dp
         CALL get_normal_scaling(scaling, c_sum, donor_state)
         CALL cp_fm_column_scale(c_sum, scaling)

      ELSE ! not TDA

!        Get the auxiliary matrices (A-D+E)^+-0.5
         paux => donor_state%matrix_aux(1)%matrix
         maux => donor_state%matrix_aux(2)%matrix

!        Make the problem Hermitian by solving M*z = omega^2*z with
!        M = (A-D+E)^0.5 * G^-1 * (A+4B+2C-E-D) * G^-1 * (A-D+E)^0.5
         ALLOCATE (tmp_mat)
         CALL dbcsr_create(matrix=tmp_mat, template=paux, matrix_type="N")

!        First compute G^-1 * (A+4B+2C-E-D) * G^-1 and store it in matrix_tdp
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_tdp, donor_state%metric(2)%matrix, 0.0_dp, tmp_mat)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, tmp_mat, 0.0_dp, matrix_tdp)

!        Compute the full matrix M and store it in matrix_tdp
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_tdp, paux, 0.0_dp, tmp_mat)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, paux, tmp_mat, 0.0_dp, matrix_tdp)

!        Shift the eigenvalues due to the unperturbed orbitals away.
         ALLOCATE (matrices(1))
         matrices(1)%matrix => matrix_tdp
         CALL setup_arnoldi_data(my_arnoldi, matrices, max_iter=50, threshold=1.0E-3_dp, selection_crit=2, &
                                 nval_request=1, nrestarts=30, generalized_ev=.FALSE., iram=.FALSE.)
         CALL arnoldi_ev(matrices, my_arnoldi)
         max_ev = REAL(get_selected_ritz_val(my_arnoldi, 1), dp)
         CALL deallocate_arnoldi_data(my_arnoldi)
         DEALLOCATE (matrices)

!        The filter is: (A-D+E)^-0.5 * SPS * (A-D+E)^-0.5
         CALL get_gs_filter(filter, matrix_tdp, donor_state, xas_tdp_env)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, filter, maux, 0.0_dp, tmp_mat)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, maux, tmp_mat, 0.0_dp, filter)
         CALL dbcsr_add(matrix_tdp, filter, 1.0_dp, max_ev+10.0_dp)

!        Get the matrix M as fm
         CALL copy_dbcsr_to_fm(matrix_tdp, rhs_matrix)

!        The problemed has been turned Hermitian. Call stendard routine
         CALL cp_fm_syevx(rhs_matrix, work, tmp_evals, nevals, work_syevx=1.0_dp)
         lr_evals(:) = tmp_evals(1:nevals)

!        Take the square root since we computed omega^2
         lr_evals = SQRT(lr_evals)

!        Need to normalize the eigenvectors. Need to retrieve the difference of the LR coeffs:
!        (c^+ - c^-) = (A-D+E)^-0.5 * z
         CALL cp_fm_create(c_diff, fm_struct)
         CALL cp_dbcsr_sm_fm_multiply(maux, work, c_diff, ncol=nevals)

!        The normalization is c_diff * G * c_diff = +- 1
         scaling = 0.0_dp
         CALL get_normal_scaling(scaling, c_diff, donor_state)
         CALL cp_fm_column_scale(c_diff, scaling)

!        Get the sum c^+ + c^- = -1/omega G^-1 * (A-D+E)*(c^+ - c^-)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, paux, paux, 0.0_dp, tmp_mat)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, tmp_mat, 0.0_dp, tmp_mat)
         CALL cp_dbcsr_sm_fm_multiply(tmp_mat, c_diff, c_sum, ncol=nevals)
         scaling = 0.0_dp
         WHERE (lr_evals > 0.0_dp) scaling = -1.0_dp/lr_evals
         CALL cp_fm_column_scale(c_sum, scaling)

!        Some clean-up
         CALL dbcsr_release(tmp_mat)
         CALL cp_fm_release(c_diff)
         DEALLOCATE (tmp_mat)

      END IF ! TDA

!  Reorganize the eigenvectors in array of cp_fm_p_types so that each element corresponds to an
!  excited state. Makes later calls to those easier and more efficient
      CALL cp_fm_get_info(donor_state%gs_coeffs, matrix_struct=gs_struct, nrow_global=nao)
      IF (singlet) THEN
         ALLOCATE (donor_state%sg_coeffs(nevals))
         lr_coeffs => donor_state%sg_coeffs
      ELSE
         ALLOCATE (donor_state%tp_coeffs(nevals))
         lr_coeffs => donor_state%tp_coeffs
      END IF

      DO i = 1, nevals

         CALL cp_fm_create(lr_coeffs(i)%matrix, gs_struct)

         DO imo = 1, ndo_mo

            CALL cp_fm_to_fm_submat(msource=c_sum, mtarget=lr_coeffs(i)%matrix, &
                                    nrow=nao, ncol=1, s_firstrow=(imo-1)*nao+1, s_firstcol=i, &
                                    t_firstrow=1, t_firstcol=imo)

         END DO !imo
      END DO !istate

!  Allocate the oscilaltor strength vectors
      IF (singlet) THEN
         ALLOCATE (donor_state%sg_osc_str(nevals))
         donor_state%sg_osc_str = 0.0_dp
      ELSE
         ALLOCATE (donor_state%tp_osc_str(nevals))
         donor_state%tp_osc_str = 0.0_dp
      END IF

!  Clean-up
      CALL cp_fm_release(rhs_matrix)
      CALL cp_fm_release(work)
      CALL cp_fm_release(c_sum)
      CALL cp_fm_struct_release(fm_struct)
      CALL dbcsr_release(filter)
      DEALLOCATE (scaling, filter)

!  Perform a partial clean-up of the donor_state
      CALL dbcsr_release(matrix_tdp)

      CALL timestop(handle)

   END SUBROUTINE solve_xas_tdp_full_prob

! **************************************************************************************************
!> \brief A geeig solver using Cholesky decomposition and syevx so that only a couple of eigenvalues
!>        can be computed.
!> amatrix ...
!> \param bmatrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \param neig the number of eigenvalues to compute
!> \param work ...
!> \note Blatantly copied from cp_fm_geeig
! **************************************************************************************************
   SUBROUTINE cp_fm_geeig_syevx(amatrix, bmatrix, eigenvectors, eigenvalues, neig, work)

      TYPE(cp_fm_type), POINTER                          :: amatrix, bmatrix, eigenvectors
      REAL(dp), DIMENSION(:)                             :: eigenvalues
      INTEGER, INTENT(IN)                                :: neig
      TYPE(cp_fm_type), POINTER                          :: work

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_geeig_syevx', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: nao, nmo

      CALL cp_fm_get_info(amatrix, nrow_global=nao)
      nmo = SIZE(eigenvalues)
      ! Cholesky decompose S=U(T)U
      CALL cp_fm_cholesky_decompose(bmatrix)
      ! Invert to get U^(-1)
      CALL cp_fm_triangular_invert(bmatrix)
      ! Reduce to get U^(-T) * H * U^(-1)
      CALL cp_fm_triangular_multiply(bmatrix, amatrix, side="R")
      CALL cp_fm_triangular_multiply(bmatrix, amatrix, transpose_tr=.TRUE.)
      ! Diagonalize
      CALL cp_fm_syevx(amatrix, work, eigenvalues, neig, work_syevx=1.0_dp)
      ! Restore vectors C = U^(-1) * C*
      CALL cp_fm_triangular_multiply(bmatrix, work)
      CALL cp_fm_to_fm(work, eigenvectors, nmo)

   END SUBROUTINE cp_fm_geeig_syevx

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
!> \brief Computes the oscillator strength based on the dipole moment (velocity or length rep) for
!>        all available excitation energies and store the results in the donor_state. There is no
!>        triplet dipole in the spin-restricted ground state
!> \param donor_state the donor state which is excited
!> \param xas_tdp_control ...
!> \param xas_tdp_env ...
!> \note The oscillator strength is a scalar: osc_str = 2/3*omega*(dipole_r)^2 = 2/(3*omega)*(dipole_v)^2
!>       The formulae for the dipoles come from the trace of the dipole operator with the transition
!>       densities, i.e. what we get from solving the xas_tdp problem. Same procedure with or wo TDA
! **************************************************************************************************
   SUBROUTINE compute_dipole_fosc(donor_state, xas_tdp_control, xas_tdp_env)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_dipole_fosc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iosc, j, nao, ndo_mo, nosc
      REAL(dp)                                           :: dip_mom
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: diag, osc_str
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: col_struct, mat_struct
      TYPE(cp_fm_type), POINTER                          :: col_work, lr_coeffs, mat_work
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipmat

      NULLIFY (dipmat, col_work, mat_work, col_struct, mat_struct, para_env, blacs_env, lr_coeffs)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices)
      nosc = SIZE(donor_state%sg_evals)
      ALLOCATE (osc_str(nosc))
      osc_str = 0.0_dp

      dipmat => xas_tdp_env%dipmat
      !compute the prefactor
      IF (xas_tdp_control%dipole_form == xas_dip_len) THEN
         osc_str(:) = 2.0_dp/3.0_dp*donor_state%sg_evals(:)
      ELSE
         osc_str(:) = 2.0_dp/3.0_dp/donor_state%sg_evals(:)
      END IF

      ! do some work matrix initialization
      CALL cp_fm_get_info(donor_state%gs_coeffs, matrix_struct=col_struct, para_env=para_env, &
                          context=blacs_env, nrow_global=nao)
      CALL cp_fm_struct_create(mat_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=ndo_mo, ncol_global=ndo_mo)
      CALL cp_fm_create(col_work, col_struct)
      CALL cp_fm_create(mat_work, mat_struct)
      ALLOCATE (diag(ndo_mo))

!  Looping over the excited singlets
      DO iosc = 1, nosc

         dip_mom = 0.0_dp
         lr_coeffs => donor_state%sg_coeffs(iosc)%matrix

         !looping over the cartesian coordinates
         DO j = 1, 3

            !computing c^0 * dipmat(j) * c^i
            CALL cp_dbcsr_sm_fm_multiply(dipmat(j)%matrix, lr_coeffs, col_work, ncol=ndo_mo)
            CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, donor_state%gs_coeffs, &
                         col_work, 0.0_dp, mat_work)

            !get the sum of the dipole moments. Factor 2 because of definition of singlet
            CALL cp_fm_get_diag(mat_work, diag)
            dip_mom = dip_mom+2.0_dp*SUM(diag)**2

         END DO !j

         !Taking the excitation energy into account
         osc_str(iosc) = osc_str(iosc)*dip_mom

      END DO !iosc

!  Saving the results in the donor state
      donor_state%sg_osc_str(:) = osc_str

      !clean-up
      CALL cp_fm_release(mat_work)
      CALL cp_fm_release(col_work)
      CALL cp_fm_struct_release(mat_struct)
      DEALLOCATE (diag)

   END SUBROUTINE compute_dipole_fosc

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
         found_block = .FALSE.
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
         found_block = .FALSE.
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
         found_block = .FALSE.
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
         found_block = .FALSE.
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
!> \note Compute the inverse only if the dimension of the incoming array is 1
! **************************************************************************************************
   SUBROUTINE build_metric(matrix_g, dist_g, blk_size_g, donor_state, qs_env, do_inv)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_g
      TYPE(dbcsr_distribution_type), POINTER             :: dist_g
      INTEGER, DIMENSION(:), POINTER                     :: blk_size_g
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_inv

      CHARACTER(len=*), PARAMETER :: routineN = 'build_metric', routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, i, iblk, jblk, nblk_row, ndo_mo
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size
      LOGICAL                                            :: found_block, my_do_inv
      REAL(dp), DIMENSION(:), POINTER                    :: work_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: matrix_sinv

      NULLIFY (matrix_s, row_blk_size, work_block)

!  Initilization
      ndo_mo = SIZE(donor_state%mo_indices, 1)
      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size)
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
         found_block = .FALSE.
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
         ALLOCATE (matrix_sinv)
         CALL dbcsr_create(matrix=matrix_sinv, template=matrix_s(1)%matrix)
         CALL invert_Hotelling(matrix_sinv, matrix_s(1)%matrix, threshold=1.0E-10_dp, use_inv_as_guess=.FALSE., &
                               norm_convergence=1.0E-10_dp, filter_eps=1.0E-12_dp, silent=.TRUE.)
         CALL dbcsr_finalize(matrix_sinv)

!     Fill the matrices G^-1 by looping over the block of S^-1 and putting them on the diagonal
         CALL dbcsr_iterator_start(iter, matrix_sinv)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

!           Get the block
            found_block = .FALSE.
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

! **************************************************************************************************!
!> \brief Computes the two-center Coulomb integral needed for the RI in kernel calculation. Given
!>        the excited kind indices, returns the (for now) cp_fm (P|Q)^-1 for the RI basis P,Q
!> \param ex_kinds the indices of the excited kinds
!> \param nao the number of AOs
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note Stores  (P|Q)^-1 in the xas_tdp_env, use the OS scheme
! **************************************************************************************************
   SUBROUTINE compute_ri_coulomb2_int(xas_tdp_env, qs_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_coulomb2_int', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ikind, nsgf
      INTEGER, DIMENSION(:), POINTER                     :: ex_kinds
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: tmp_int
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (ri_basis, qs_kind_set)
      CPASSERT(ASSOCIATED(xas_tdp_env))

!  Initialization
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      ex_kinds => xas_tdp_env%ex_kind_indices
      ALLOCATE (xas_tdp_env%ri_inv_mats(SIZE(ex_kinds)))

!  Loop over the excited kinds
      DO i = 1, SIZE(ex_kinds)
         ikind = ex_kinds(i)

!        Get the RI basis of interest and its quantum numbers
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_XAS")
         CALL get_gto_basis_set(ri_basis, nsgf=nsgf)

!        The 2D array with the integrals for this kind
         ALLOCATE (tmp_int(nsgf, nsgf))
         tmp_int = 0.0_dp

         CALL int_operators_r12_ab_os(r12_operator=operator_coulomb, vab=tmp_int, &
                                      rab=(/0.0_dp, 0.0_dp, 0.0_dp/), fba=ri_basis, fbb=ri_basis, &
                                      calculate_forces=.FALSE.)

!        Inverting and storing in the xas_tdp_env
         CALL invmat_symm(tmp_int)

         ALLOCATE (xas_tdp_env%ri_inv_mats(i)%array(nsgf, nsgf))
         xas_tdp_env%ri_inv_mats(i)%array(:, :) = tmp_int(:, :)

!        clean-up before next kind
         DEALLOCATE (tmp_int)

      END DO ! loop on kinds

   END SUBROUTINE compute_ri_coulomb2_int

! **************************************************************************************************
!> \brief Computes the three-centers Coulomb integrals needed for RI in kernel calculations and
!>        stores it as 3D arrays. Computes the integrals (ab|P), where a,b are overlaping atomic
!>        basis set and P is the RI basis centered on the given excited atom. The first dimension
!>        of the returned array corresponds to sgf a, the second to sgf b and the third to sgf P.
!>        Keep the result in xas_tdp_env as the integrals for the current excited atom
!> \param exat_index the index of the excited atom, on which the RI basis P is centered
!> \para, exat_kind the kind of the excited atom, to get the RI basis
!> \param xas_tdp_env the xas_tdp_Env to update with the 3-center Coulomb integrals
!> \param qs_env ...
!> \note Currently creates a 3D array, will later switch to DBCSR tensors (TODO). Will be used for
!>       the Coulomb kernel after contraction (aI|P) = sum_b C_b (ab|P) and possibly for the exchange
!>       part too (if non-PBCs and coulomb operator requested). TODO: if exchange operator is not
!>       coulomb, screen the integrals to only compute (ab|P) where C_b is non-zero
! **************************************************************************************************
   SUBROUTINE compute_ri_coulomb3_int(exat_index, exat_kind, xas_tdp_env, qs_env)

      INTEGER, INTENT(IN)                                :: exat_index, exat_kind
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_coulomb3_int', &
         routineP = moduleN//':'//routineN

      INTEGER :: aat, akind, aset, bat, bkind, bset, endp, i, igp, ipgf, nao, natom, ncoa, ncob, &
         ncop, nkind, nseta, nsetb, nsetp, pset, sgfa, sgfb, sgfp, startp
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf_set
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, lp_max, &
                                                            lp_min, npgfa, npgfb, npgfp, nsgfa, &
                                                            nsgfb, nsgfp
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb, first_sgfp
      REAL(dp)                                           :: rab2, rap2, rbp2, rpgfp
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: f, gccc, rpgfa, rpgfb
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: vabc
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: int_abp, int_abp_set, int_array
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :, :)       :: v
      REAL(dp), DIMENSION(3)                             :: ra, rab, rb, rp
      REAL(dp), DIMENSION(:, :), POINTER                 :: sphi_a, sphi_b, sphi_p, zeta, zetb, zetp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: a_basis, b_basis, p_basis
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, p_basis, b_basis, a_basis, particle_set, mos)
      NULLIFY (lp_max, lp_min, npgfp, zetp, sphi_p, first_sgfp)
      NULLIFY (lb_max, lb_min, npgfb, zetb, sphi_b, first_sgfb)
      NULLIFY (la_max, la_min, npgfa, zeta, sphi_a, first_sgfa)
      NULLIFY (sab_orb, nl_iterator, cell, basis_set_list, para_env)

!  Get all we need from the qs_env
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, mos=mos, &
                      sab_orb=sab_orb, cell=cell, natom=natom, para_env=para_env)
      CALL get_mo_set(mos(1)%mo_set, nao=nao)
      ALLOCATE (first_sgf_set(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf_set)

!  Retrieve the basis set P centered on the excited atom
      CALL get_qs_kind(qs_kind_set(exat_kind), basis_set=p_basis, basis_type="RI_XAS")

      CALL get_gto_basis_set(p_basis, lmax=lp_max, npgf=npgfp, zet=zetp, lmin=lp_min, nset=nsetp, &
                             nsgf_set=nsgfp, sphi=sphi_p, first_sgf=first_sgfp)

!  Allocate the 3D array where results will be stored
      ALLOCATE (int_array(nao, nao, SUM(nsgfp)))
      int_array = 0.0_dp

!  Create lists of basis for atoms a and b
      nkind = SIZE(qs_kind_set)
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      !        TODO: do that with arnoldi geeig !!!
!  TODO: be very carfeul when enabeling PBCs
!  Loop over all neighbouring atoms a,b and compute (ab|P)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=akind, jkind=bkind, iatom=aat, jatom=bat, r=rab)

!        Retrieve the basis sets for a and b
         a_basis => basis_set_list(akind)%gto_basis_set
         b_basis => basis_set_list(bkind)%gto_basis_set

!        Retrieve the quantum numbers for both basis sets
         CALL get_gto_basis_set(a_basis, lmax=la_max, npgf=npgfa, zet=zeta, lmin=la_min, &
                                nset=nseta, nsgf_set=nsgfa, sphi=sphi_a, first_sgf=first_sgfa)
         CALL get_gto_basis_set(b_basis, lmax=lb_max, npgf=npgfb, zet=zetb, lmin=lb_min, &
                                nset=nsetb, nsgf_set=nsgfb, sphi=sphi_b, first_sgf=first_sgfb)

!        Get the different positions/distances
         ra = pbc(particle_set(aat)%r, cell)
         rb = pbc(particle_set(bat)%r, cell)
         rp = pbc(particle_set(exat_index)%r, cell)
         rab2 = DOT_PRODUCT(rb-ra, rb-ra)
         rap2 = DOT_PRODUCT(rp-ra, rp-ra)
         rbp2 = DOT_PRODUCT(rp-rb, rp-rb)

!        Loop over sets of the 3 basis
         DO aset = 1, nseta
            DO bset = 1, nsetb
               DO pset = 1, nsetp

!                 Some initialization
                  sgfa = first_sgfa(1, aset)
                  sgfb = first_sgfb(1, bset)
                  sgfp = first_sgfp(1, pset)
                  ncoa = npgfa(aset)*ncoset(la_max(aset))
                  ncob = npgfb(bset)*ncoset(lb_max(bset))
                  ncop = npgfp(pset)*ncoset(lp_max(pset))

!                 Allocate array for the current comnbination of sets
                  ALLOCATE (int_abp_set(nsgfa(aset), nsgfb(bset), nsgfp(pset)))
                  int_abp_set = 0.0_dp

!                 Allocate some work arrays
                  ALLOCATE (f(0:la_max(aset)+lb_max(bset)+lp_max(pset)+2))
                  ALLOCATE (v(ncoa, ncob, ncop, la_max(aset)+lb_max(bset)+lp_max(pset)+1))
                  ALLOCATE (gccc(ncop))
                  ALLOCATE (vabc(ncoa, ncob))
                  ALLOCATE (int_abp(ncoa, ncob, ncop))
                  f = 0.0_dp
                  v = 0.0_dp
                  gccc = 0.0_dp
                  vabc = 0.0_dp
                  int_abp = 0.0_dp

!                 hack: set the rpgf to avoid screening out
                  ALLOCATE (rpgfa(npgfa(aset)))
                  ALLOCATE (rpgfb(npgfb(bset)))
                  rpgfa = 1.0E10_dp
                  rpgfb = 1.0E10_dp
                  rpgfp = 1.0E10_dp

!                 Go one pgf at a time because of how coulomb3 is designed
                  endp = 0
                  DO ipgf = 1, npgfp(pset)
                     startp = endp+ncoset(lp_min(pset)-1)+1
                     endp = endp+ncoset(lp_max(pset))
                     CALL coulomb3(la_max(aset), npgfa(aset), zeta(:, aset), rpgfa(:), la_min(aset), &
                                   lb_max(bset), npgfb(bset), zetb(:, bset), rpgfb(:), lb_min(bset), &
                                   lp_max(pset), zetp(ipgf, pset), rpgfp, lp_min(pset), gccc, &
                                   rb-ra, rab2, rp-ra, rap2, rbp2, vabc, &
                                   int_abp(:, :, startp:endp), v, f)
                  END DO !ipgf

!                 Contract into sgf
                  CALL abc_contract(int_abp_set, int_abp, sphi_a(:, sgfa:), sphi_b(:, sgfb:), &
                                    sphi_p(:, sgfp:), ncoa, ncob, ncop, nsgfa(aset), nsgfb(bset), &
                                    nsgfp(pset))

!                 Put it at the right place in the final array
                  int_array(first_sgf_set(aat)+sgfa-1:first_sgf_set(aat)+sgfa+nsgfa(aset)-2, &
                            first_sgf_set(bat)+sgfb-1:first_sgf_set(bat)+sgfb+nsgfb(bset)-2, &
                            sgfp:sgfp+nsgfp(pset)-1) = int_abp_set

!                 In the neighbor list iterator, only get a,b pair once. But (ab|P) = (ba|P)
                  DO i = 1, nsgfp(pset)
                     igp = sgfp-1+i
                     int_array(first_sgf_set(bat)+sgfb-1:first_sgf_set(bat)+sgfb+nsgfb(bset)-2, &
                               first_sgf_set(aat)+sgfa-1:first_sgf_set(aat)+sgfa+nsgfa(aset)-2, &
                               igp) = TRANSPOSE(int_abp_set(:, :, i))
                  END DO

!                 Deallocate before next loop
                  DEALLOCATE (int_abp_set, int_abp, f, v, gccc, vabc, rpgfa, rpgfb)

               END DO ! pset
            END DO ! bset
         END DO !aset

      END DO !neighbor iteration on a,b
      CALL neighbor_list_iterator_release(nl_iterator)

!  Updating the xas_tdp_env with the new integrals
      IF (ASSOCIATED(xas_tdp_env%ri_coul3)) DEALLOCATE (xas_tdp_env%ri_coul3)
      ALLOCATE (xas_tdp_env%ri_coul3(nao, nao, SUM(nsgfp)))
      xas_tdp_env%ri_coul3 = int_array

!  Clean-up
      DEALLOCATE (basis_set_list)

   END SUBROUTINE compute_ri_coulomb3_int

! **************************************************************************************************
!> \brief Contracts the RI coulomb 3 integrals with the groud state coefficients to get (aI|P) for
!>        coulomb kernel or (IJ|P) for the exchange. For now, returns an array. TODO: switch to
!>        dbcsr tensors
!> \param contr_int the array containing the contracted integrals
!> \param domo_index the indices of the donor MOs for which we contract
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note Currently working with array in a serial manner. Will be improved with the switch to tensors
!>       Will also change: call to the gs_coeffs
! **************************************************************************************************
   SUBROUTINE contract_int_kernel(contr_int, domo_index, xas_tdp_env, qs_env)

      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: contr_int
      INTEGER, DIMENSION(:)                              :: domo_index
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'contract_int_kernel', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, nao, nsgfp
      LOGICAL                                            :: contract_two
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: coeffs, contr_once, contr_twice
      REAL(dp), DIMENSION(:, :, :), POINTER              :: integrals
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      NULLIFY (mos, mo_coeff, integrals)

!  Check whther to contract on 1 or 2 donor MOs
      IF (SIZE(domo_index) == 1) THEN
         contract_two = .FALSE.
      ELSE IF (SIZE(domo_index) == 2) THEN
         contract_two = .TRUE.
      ELSE
         CPABORT("Cannot contract integrals for more than two donor MOs")
      END IF

!  Initialization
      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao)
      integrals => xas_tdp_env%ri_coul3
      nsgfp = SIZE(integrals, 3)
      ALLOCATE (contr_once(nao, nsgfp))
      contr_once = 0.0_dp
      IF (contract_two) ALLOCATE (contr_twice(1, nsgfp)) ! not efficient but not there to stay

!  Get the MO coeffs for the first contraction
      ALLOCATE (coeffs(nao, 1))
      CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=coeffs, start_row=1, start_col=domo_index(1), &
                               n_rows=nao, n_cols=1, transpose=.FALSE.)

      DO i = 1, nao
         CALL daxpy(nao*nsgfp, coeffs(i, 1), integrals(:, i, :), 1, contr_once(:, :), 1)
      END DO

!  Do the second contraction if requested
      IF (contract_two) THEN
         contr_twice = 0.0_dp
         CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=coeffs, start_row=1, &
                                  start_col=domo_index(2), n_rows=nao, n_cols=1, transpose=.FALSE.)

         DO i = 1, nao
            CALL daxpy(nsgfp, coeffs(i, 1), contr_once(i, :), 1, contr_twice(1, :), 1)
         END DO
      END IF

!  Allocate and assign the contr_int array
      IF (.NOT. contract_two) THEN
         ALLOCATE (contr_int(nao, nsgfp))
         contr_int(:, :) = contr_once(:, :)
      ELSE
         ALLOCATE (contr_int(1, nsgfp))
         contr_int(:, :) = contr_twice(:, :)
      END IF

   END SUBROUTINE contract_int_kernel

! **************************************************************************************************
!> \brief Computes the Coulomb part of the kernel using RI.
!> \param coulomb_kernel the dbcsr_matrix where the coulomb kernel is stored
!> \param dist the dbcsr_distribution of the matrix. Same as for matrices A,B or G
!> \param blk_size the block sizes of the matrix
!> \param donor_state the donor_state for which we compute the kernel
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note currently coded with arrays and such low level types. TODO: switch to dbcsr tensors
!>       Everythin is gonna change later on, so don't care about performance
! **************************************************************************************************
   SUBROUTINE build_kernel_coulomb(coulomb_kernel, dist, blk_size, donor_state, xas_tdp_env, qs_env)

      TYPE(dbcsr_type), POINTER                          :: coulomb_kernel
      TYPE(dbcsr_distribution_type), POINTER             :: dist
      INTEGER, DIMENSION(:), POINTER                     :: blk_size
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_kernel_coulomb', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iblk, j, jblk, nao, nblk, &
                                                            ncol_blk, ndo_mo, nrow_blk, nsgfp
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: aIJb_prod, int_aI, int_bJ, work
      REAL(dp), DIMENSION(:, :), POINTER                 :: PQ
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (PQ, qs_kind_set, particle_set)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices, 1)
      CALL get_qs_env(qs_env, natom=nblk, particle_set=particle_set, qs_kind_set=qs_kind_set)
      ALLOCATE (first_sgf(nblk))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf)

!     Get the (P|Q)^-1 from the ri_inv_mats
      i = locate(xas_tdp_env%ex_kind_indices, donor_state%kind_index)
      PQ => xas_tdp_env%ri_inv_mats(i)%array
      nsgfp = SIZE(PQ, 1)
      nao = SIZE(xas_tdp_env%ri_coul3, 1)

!     Build the dbcsr matrix for the kernel
      CALL dbcsr_create(matrix=coulomb_kernel, name="COULOMB KERNEL", matrix_type='S', dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)

!     Allocate the work arrays to store the intermediate matrix products
      ALLOCATE (work(nao, nsgfp))
      ALLOCATE (aIJb_prod(nao, nao))

!  Loop over donor MOs and compute W_ab = (aI|P)* (P|Q)^-1 * (Q|Jb)
!  Matrix is symmetric => only span the upper triangle
      DO i = 1, ndo_mo
         CALL contract_int_kernel(int_aI, (/donor_state%mo_indices(i)/), xas_tdp_env, qs_env)
         DO j = i, ndo_mo
            CALL contract_int_kernel(int_bJ, (/donor_state%mo_indices(j)/), xas_tdp_env, qs_env)

            work = 0.0_dp
            aIJb_prod = 0.0_dp

!           Do the product (aI|P)*(P|Q)^-1
            CALL dgemm('N', 'N', nao, nsgfp, nsgfp, 1.0_dp, int_aI, nao, PQ, nsgfp, 0.0_dp, &
                       work, nao)

!           Mulitply by (Q|Jb)
            CALL dgemm('N', 'T', nao, nao, nsgfp, 1.0_dp, work, nao, int_bJ, nao, 0.0_dp, &
                       aIJb_prod, nao)

!           Put it in the dbcsr matrix, for the right I,J and block by block
            DO iblk = 1, nblk
               DO jblk = 1, nblk

!                 If on the block diagonal, just take jblk >= iblk
                  IF (jblk < iblk .AND. i == j) CYCLE

!                 Size of the current block
                  nrow_blk = blk_size(iblk)
                  ncol_blk = blk_size(jblk)

                  CALL dbcsr_put_block(coulomb_kernel, (i-1)*nblk+iblk, (j-1)*nblk+jblk, &
                                       aIJb_prod(first_sgf(iblk):first_sgf(iblk)+nrow_blk-1, &
                                                 first_sgf(jblk):first_sgf(jblk)+ncol_blk-1))

               END DO ! jblk
            END DO ! iblk

            DEALLOCATE (int_bJ)
         END DO ! donor MO j
         DEALLOCATE (int_aI)
      END DO ! donor MO i

      CALL dbcsr_finalize(coulomb_kernel)

      CALL timestop(handle)

   END SUBROUTINE build_kernel_coulomb

! **************************************************************************************************
!> \brief Computes the exact exchange part of the kernel using RI distinguishing between on- and
!>        off-diagonal blocks of matrix_tdp
!> \param ex_ker_diag the dbcsr_matrix where the  on-diagonal exchange kernel is stored
!> \param ex_ker_offdiag the dbcsr_matrix where the off-diagonal exchange kernel is stored
!> \param dist the dbcsr distribution of the matrix. Same as for matrices A and G
!> \param blk_size the block sizes of the matrix
!> \param donor_state the donor_state for which we compute the kernel
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note currently coded with arrays and such low level types. TODO: switch to dbcsr tensors
!>       Two matrices must be computed for the exact exchange: (ab|IJ) for the digonal of matrix_tdp
!>       and (aJ|Ib) for the off diagonal part. Note that is not the same as the Coulomb (aI|Jb)
! **************************************************************************************************
   SUBROUTINE build_kernel_exchange(ex_ker_diag, ex_ker_offdiag, dist, blk_size, donor_state, &
                                    xas_tdp_env, qs_env)

      TYPE(dbcsr_type), POINTER                          :: ex_ker_diag
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: ex_ker_offdiag
      TYPE(dbcsr_distribution_type), POINTER             :: dist
      INTEGER, DIMENSION(:), POINTER                     :: blk_size
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_kernel_exchange', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iblk, j, jblk, k, nao, nblk, &
                                                            ncol_blk, ndo_mo, nrow_blk, nsgfp
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf
      LOGICAL                                            :: do_offdiag
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: abIJ_prod, aJIb_prod, int_aJ, int_bI, &
                                                            int_IJ, work1, work3
      REAL(dp), DIMENSION(:, :), POINTER                 :: PQ
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (particle_set, qs_kind_set, PQ)

!  Initialization
      do_offdiag = .FALSE.
      IF (PRESENT(ex_ker_offdiag)) do_offdiag = .TRUE.
      ndo_mo = SIZE(donor_state%mo_indices, 1)
      CALL get_qs_env(qs_env, natom=nblk, particle_set=particle_set, qs_kind_set=qs_kind_set)
      ALLOCATE (first_sgf(nblk))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf)

!     Get the (P|Q)^-1 from the xas_tdp_env
      i = locate(xas_tdp_env%ex_kind_indices, donor_state%kind_index)
      PQ => xas_tdp_env%ri_inv_mats(i)%array
      nsgfp = SIZE(PQ, 1)
      nao = SIZE(xas_tdp_env%ri_coul3, 1)

!     Build the dbcsr matrix for the kernels
      CALL dbcsr_create(matrix=ex_ker_diag, name="EX KER DIAG", matrix_type='S', dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)
      IF (do_offdiag) THEN
         CALL dbcsr_create(matrix=ex_ker_offdiag, name="EX KER OFFDIAG", matrix_type='S', &
                           dist=dist, row_blk_size=blk_size, col_blk_size=blk_size)
      END IF

!     Allocate the many work arrays
      IF (do_offdiag) THEN
         ALLOCATE (aJIb_prod(nao, nao))
         ALLOCATE (work1(nao, nsgfp))
      END IF
      ALLOCATE (work3(nsgfp, 1))
      ALLOCATE (abIJ_prod(nao, nao))

!  Loop over the donor MOs and compute both the digonal and off-digonal parts, that is:
!     -Diag_ab = (ab|P)* (P|Q)^-1 *(Q|IJ)
!     -Off_ab = (aJ|P)* (P|Q)^-1 *(Q|Ib)
!  The matrices are symmetric => only loop on upper triangular part
      DO i = 1, ndo_mo
         IF (do_offdiag) CALL contract_int_kernel(int_bI, (/donor_state%mo_indices(i)/), &
                                                  xas_tdp_env, qs_env)
         DO j = i, ndo_mo
            IF (do_offdiag) CALL contract_int_kernel(int_aJ, (/donor_state%mo_indices(j)/), &
                                                     xas_tdp_env, qs_env)
            CALL contract_int_kernel(int_IJ, (/donor_state%mo_indices(i), donor_state%mo_indices(j)/), &
                                     xas_tdp_env, qs_env)

            IF (do_offdiag) THEN
               work1 = 0.0_dp
               aJIb_prod = 0.0_dp
            END IF
            work3 = 0.0_dp
            abIJ_prod = 0.0_dp

!           Start with off-diag product (aJ|P)* (P|Q)^-1 *(Q|Ib)
            IF (do_offdiag) THEN

!              (aJ|P)* (P|Q)^-1
               CALL dgemm('N', 'N', nao, nsgfp, nsgfp, 1.0_dp, int_aJ, nao, PQ, nsgfp, 0.0_dp, &
                          work1, nao)
!              Muliply by (Q|Ib)
               CALL dgemm('N', 'T', nao, nao, nsgfp, 1.0_dp, work1, nao, int_bI, nao, 0.0_dp, &
                          aJIb_prod, nao)
            END IF

!           Then the on-digonal product (ab|P)* (P|Q)^-1 *(Q|IJ)

!           (P|Q)^-1 *(Q|IJ)
            CALL dgemm('N', 'T', nsgfp, 1, nsgfp, 1.0_dp, PQ, nsgfp, int_IJ, 1, 0.0_dp, &
                       work3, nsgfp)

!           Multiply by (ab|P) from the left
            DO k = 1, nsgfp
               CALL daxpy(nao*nao, work3(k, 1), xas_tdp_env%ri_coul3(:, :, k), 1, abIJ_prod(:, :), 1)
            END DO !k

!           Put it all in the dbcsr matrices, at the right I,J and block by block
            DO iblk = 1, nblk
               DO jblk = 1, nblk

!                 If on the block digonal, only take jblk >= iblk
                  IF (jblk < iblk .AND. i == j) CYCLE

!                 Size of the current block
                  nrow_blk = blk_size(iblk)
                  ncol_blk = blk_size(jblk)

!                 Put (ab|IJ) on the digonal exchange kernel block, at position I,J
                  CALL dbcsr_put_block(ex_ker_diag, (i-1)*nblk+iblk, (j-1)*nblk+jblk, &
                                       abIJ_prod(first_sgf(iblk):first_sgf(iblk)+nrow_blk-1, &
                                                 first_sgf(jblk):first_sgf(jblk)+ncol_blk-1))

                  IF (do_offdiag) THEN
!                    Put (aJ|Ib) on the off digonal exchange kernel block, at position I,J
                     CALL dbcsr_put_block(ex_ker_offdiag, (i-1)*nblk+iblk, (j-1)*nblk+jblk, &
                                          aJIb_prod(first_sgf(iblk):first_sgf(iblk)+nrow_blk-1, &
                                                    first_sgf(jblk):first_sgf(jblk)+ncol_blk-1))
                  END IF

               END DO !jblk
            END DO !iblk
            IF (do_offdiag) DEALLOCATE (int_aJ)
            DEALLOCATE (int_IJ)
         END DO !j
         IF (do_offdiag) DEALLOCATE (int_bI)
      END DO !i

      CALL dbcsr_finalize(ex_ker_diag)
      IF (do_offdiag) CALL dbcsr_finalize(ex_ker_offdiag)

      CALL timestop(handle)

   END SUBROUTINE build_kernel_exchange

! **************************************************************************************************
!> \brief Computes the XC part of the kernel base on atomic grid integrals
!> \param xc_kernel the dbcsr_matrix where the xc kernel is stored
!> \param dist the dbcsr_distribution of the matrix. Same as for matrices A,B or G
!> \param blk_size the block sizes of the matrix
!> \param donor_state the donor_state for which we compute the kernel
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \param singlet whether the kernel is built for singlet excitations (vs triplet)
!> \note currently coded with arrays and such low level types. TODO: switch to dbcsr tensors
!>       Everythin is gonna change later on, so don't care about performance
! **************************************************************************************************
   SUBROUTINE build_kernel_xc(xc_kernel, dist, blk_size, donor_state, xas_tdp_env, qs_env, singlet)

      TYPE(dbcsr_type), POINTER                          :: xc_kernel
      TYPE(dbcsr_distribution_type), POINTER             :: dist
      INTEGER, DIMENSION(:), POINTER                     :: blk_size
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: singlet

      CHARACTER(len=*), PARAMETER :: routineN = 'build_kernel_xc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iblk, j, jblk, nao, nblk, &
                                                            ncol_blk, ndo_mo, nrow_blk, nsgfp
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf
      LOGICAL                                            :: found_block
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: aIJb_prod, int_aI, int_bJ, work1, work2
      REAL(dp), DIMENSION(:, :), POINTER                 :: fxc, PQ
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (particle_set, qs_kind_set, PQ)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices, 1)
      CALL get_qs_env(qs_env, natom=nblk, particle_set=particle_set, qs_kind_set=qs_kind_set)
      ALLOCATE (first_sgf(nblk))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf)

!     Get the (P|Q)^-1 as an array
      i = locate(xas_tdp_env%ex_kind_indices, donor_state%kind_index)
      PQ => xas_tdp_env%ri_inv_mats(i)%array
      nsgfp = SIZE(PQ, 1)
      nao = SIZE(xas_tdp_env%ri_coul3, 1)

!     Get the 2-center integral (Q|fxc|R) as an array (a block from a dbcsr matrix)
!     For now assume serial => don't check for which proc the block is on
!     TODO: maybe that storing the atom integrals in a dbcsr matrix is not the best
      found_block = .FALSE.
      IF (singlet) THEN
         CALL dbcsr_get_block_p(xas_tdp_env%ri_fxc_sg, donor_state%at_index, &
                                donor_state%at_index, fxc, found_block)
      ELSE
         CALL dbcsr_get_block_p(xas_tdp_env%ri_fxc_tp, donor_state%at_index, &
                                donor_state%at_index, fxc, found_block)
      END IF
      CPASSERT(found_block)
      CPASSERT(ALL(SHAPE(fxc) == nsgfp))

!     Build the dbcsr matrix for the kernel
      CALL dbcsr_create(matrix=xc_kernel, name="XC KERNEL", matrix_type='S', dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)

!     Allocating work arrays to store the intemediate products
      ALLOCATE (work1(nao, nsgfp))
      ALLOCATE (work2(nao, nsgfp))
      ALLOCATE (aIJb_prod(nao, nao))

!  Loop over the donor MOs and compute the matrix product:
!  W_ab = (aI|P) * (P|Q)^-1 * (Q|fxc|R) * (R|S)^-1 * (S|Jb) ~= (aI|fxc|Jb)
!  Symmetric matrix => only span upper triangle
      DO i = 1, ndo_mo
         CALL contract_int_kernel(int_aI, (/donor_state%mo_indices(i)/), xas_tdp_env, qs_env)
         DO j = i, ndo_mo
            CALL contract_int_kernel(int_bJ, (/donor_state%mo_indices(j)/), xas_tdp_env, qs_env)

            work1 = 0.0_dp
            work2 = 0.0_dp
            aIJb_prod = 0.0_dp

!           Do the product (aI|P)*(P|Q)^-1 and stroe in work1
            CALL dgemm('N', 'N', nao, nsgfp, nsgfp, 1.0_dp, int_aI, nao, PQ, nsgfp, 0.0_dp, &
                       work1, nao)

!           Muliply the above by (Q|fxc|R) and store in work2
            CALL dgemm('N', 'N', nao, nsgfp, nsgfp, 1.0_dp, work1, nao, fxc, nsgfp, 0.0_dp, &
                       work2, nao)
            work1 = 0.0_dp

!           Multiply the above by (R|S)^-1 => (aI|P)*(P|Q)^-1 *(Q|fxc|R)*(R|S)^-1 = work1
            CALL dgemm('N', 'N', nao, nsgfp, nsgfp, 1.0_dp, work2, nao, PQ, nsgfp, 0.0_dp, &
                       work1, nao)

!           Multiply the above by (S|Jb) = (Jb|S)^T and temporarly store in aIJb_prod
            CALL dgemm('N', 'T', nao, nao, nsgfp, 1.0_dp, work1, nao, int_bJ, nao, 0.0_dp, &
                       aIJb_prod, nao)

!           Put the product in the dbcsr_matrix, for the right I,J block
            DO iblk = 1, nblk
               DO jblk = 1, nblk

!                 If on the block diagonal, just take jblk >= iblk
                  IF (jblk < iblk .AND. i == j) CYCLE

!                 Size of the current block
                  nrow_blk = blk_size(iblk)
                  ncol_blk = blk_size(jblk)

                  CALL dbcsr_put_block(xc_kernel, (i-1)*nblk+iblk, (j-1)*nblk+jblk, &
                                       aIJb_prod(first_sgf(iblk):first_sgf(iblk)+nrow_blk-1, &
                                                 first_sgf(jblk):first_sgf(jblk)+ncol_blk-1))

               END DO ! jblk
            END DO ! iblk

            DEALLOCATE (int_bJ)
         END DO ! donor MO j
         DEALLOCATE (int_aI)
      END DO ! donor MO i

      CALL dbcsr_finalize(xc_kernel)

      CALL timestop(handle)

   END SUBROUTINE build_kernel_xc

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
      INTEGER                                            :: homo, i, isg, itp, jtp, nao, ndo_mo, &
                                                            nsg, ntot, ntp
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

! TODO: can probably compute orb_soc_xyz * tp_coeffs(j)%matrix only once for each jtp and store

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

