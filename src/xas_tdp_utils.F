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
   USE arnoldi_api,                     ONLY: arnoldi_extremal
   USE basis_set_types,                 ONLY: gto_basis_set_type, get_gto_basis_set, &
                                              gto_basis_set_p_type
   USE cell_types,                      ONLY: cell_type, pbc
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm, cp_dbcsr_sm_fm_multiply, &
                                              copy_fm_to_dbcsr
   USE cp_cfm_types,                    ONLY: cp_fm_to_cfm, cp_cfm_to_fm, cp_cfm_type, &
                                              cp_cfm_create, cp_cfm_release
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add, cp_fm_gemm, cp_fm_column_scale
   USE cp_fm_diag,                      ONLY: cp_fm_geeig, cp_fm_power, cp_fm_syevd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create, cp_fm_struct_release, &     
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create, cp_fm_release, cp_fm_type, &
                                              cp_fm_to_fm_submat, cp_fm_get_submatrix, &
                                              cp_fm_get_info, cp_fm_set_submatrix, cp_fm_set_all, &
                                              cp_fm_p_type, cp_fm_get_diag
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_type, dbcsr_p_type, dbcsr_distribution_type, &
                                              dbcsr_get_info, dbcsr_distribution_get, &
                                              dbcsr_distribution_new, dbcsr_create, &
                                              dbcsr_get_block_p, dbcsr_print, &
                                              dbcsr_put_block, dbcsr_release, dbcsr_finalize,&
                                              dbcsr_distribution_release, &
                                              dbcsr_copy, dbcsr_set, dbcsr_multiply, &
                                              dbcsr_get_stored_coordinates, dbcsr_desymmetrize, &
                                              dbcsr_complete_redistribute, dbcsr_iterator_type, &
                                              dbcsr_iterator_blocks_left, dbcsr_iterator_start, &
                                              dbcsr_iterator_stop, dbcsr_iterator_next_block, &
                                              dbcsr_add
   USE generic_os_integrals,            ONLY: int_operators_r12_ab_os
   USE input_constants,                 ONLY: xas_dip_len, xas_dip_vel, do_eri_mme, do_eri_os, &
                                              xas_tdp_hfx, xas_tdp_pade
   USE iterate_matrix,                  ONLY: invert_Hotelling
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat_symm
   USE orbital_pointers,                ONLY: ncoset
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: qs_environment_type, get_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind, qs_kind_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_mo_types,                     ONLY: get_mo_set, mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info, neighbor_list_iterate, &
                                              neighbor_list_iterator_create, &
                                              neighbor_list_iterator_p_type, &
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE util,                            ONLY: locate
   USE xas_tdp_types,                   ONLY: donor_state_type, xas_tdp_env_type, &
                                              xas_tdp_control_type

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_utils'

   PUBLIC :: setup_xas_tdp_full_prob, solve_xas_tdp_full_prob, compute_dipole_fosc, &
             compute_ri_coulomb3_int, compute_ri_coulomb2_int

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
!>       For the PADE xc functional, C_aa + C_ab is the same as 2*C (with d2fxc_d2n)
!>                               and C_aa - C_ab simplifies a lot
!>       In the above definitions, I,J label the donnor MOs and a,b the sgfs of the basis
!>       In the case of non-TDA, one turns the problem Hermitian 
! **************************************************************************************************
   SUBROUTINE setup_xas_tdp_full_prob(donor_state, qs_env, xas_tdp_env, xas_tdp_control)

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      
      CHARACTER(len=*), PARAMETER :: routineN = 'setup_xas_tdp_full_prob', &
                                     routineP = moduleN//":"//routineN
      
      INTEGER                                         :: handle , fxc, nrow, n_dependent, nblk_sub
      TYPE(dbcsr_type), POINTER                       :: matrix_tdp
      TYPE(dbcsr_type), POINTER                       :: matrix_a, matrix_b, matrix_d, matrix_e, &
                                                         tmp_mat, matrix_c
      TYPE(dbcsr_distribution_type), POINTER          :: submat_dist
      INTEGER, DIMENSION(:), POINTER                  :: submat_blk_size, col_dist_sub, &
                                                         row_dist_sub
      INTEGER, DIMENSION(:,:), POINTER                :: pgrid
      LOGICAL                                         :: conv
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct
      TYPE(cp_fm_type), POINTER                       :: fm_tmp, work
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      REAL(dp)                                        :: min_ev, max_ev

      CALL timeset(routineN, handle)

      NULLIFY(matrix_a, submat_dist, submat_blk_size, col_dist_sub, row_dist_sub)
      NULLIFY(pgrid, matrix_tdp, matrix_b, matrix_d, matrix_e, tmp_mat, matrix_c) 
      NULLIFY(fm_struct, fm_tmp, para_env, blacs_env, work)

      CPASSERT(ASSOCIATED(donor_state))
      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(xas_tdp_env))
      CPASSERT(ASSOCIATED(xas_tdp_control))
      ALLOCATE(donor_state%matrix_tdp)
      matrix_tdp => donor_state%matrix_tdp
      fxc = xas_tdp_control%fxc_func

!  Get the submatrices distribution and block sizes
      CALL compute_submat_dist_and_blk_size(submat_dist, submat_blk_size, donor_state, qs_env)
      nblk_sub = SIZE(submat_blk_size)

!  IF TDA, simplified problem:
      IF (xas_tdp_control%tamm_dancoff) THEN

!     Compute the submatrices, store the metric and ignore the offdiag exchange kernel
!     Only compute needed forms of the kernel according to the functional
         ALLOCATE(matrix_a)
         ALLOCATE(matrix_b)
         ALLOCATE(matrix_c)
         ALLOCATE(matrix_d)
         ALLOCATE(donor_state%metric(1))
         CALL build_gs_contribution(matrix_a, submat_dist, submat_blk_size, donor_state, xas_tdp_env, qs_env)
         CALL build_kernel_coulomb(matrix_b, submat_dist, submat_blk_size, donor_state, xas_tdp_env, qs_env)
         IF (fxc == xas_tdp_pade) CALL build_kernel_xc(matrix_c, submat_dist, &
                                                       submat_blk_size, donor_state, xas_tdp_env, qs_env)
         IF (fxc == xas_tdp_hfx ) CALL build_kernel_exchange(ex_ker_diag=matrix_d, &
                                                            dist=submat_dist, blk_size=submat_blk_size, &
                                                            donor_state=donor_state, xas_tdp_env=xas_tdp_env, &
                                                            qs_env=qs_env)     
         CALL build_metric(donor_state%metric, submat_dist, submat_blk_size, donor_state, qs_env)

!     Create the simplfied TDA matrix_tdp
         CALL dbcsr_copy(matrix_tdp, matrix_a, name="MATRIX TDP")
         IF (xas_tdp_control%do_singlet) CALL dbcsr_add(matrix_tdp, matrix_b, 1.0_dp, 2.0_dp) 
         IF (fxc == xas_tdp_pade) CALL dbcsr_add(matrix_tdp, matrix_c, 1.0_dp, 1.0_dp) !pure pade
         IF (fxc == xas_tdp_hfx) CALL dbcsr_add(matrix_tdp, matrix_d, 1.0_dp, -1.0_dp) !pure hfx

!     Some clean-up
         CALL dbcsr_release(matrix_a)   
         CALL dbcsr_release(matrix_b) 
         CALL dbcsr_release(matrix_c)
         CALL dbcsr_release(matrix_d)
         DEALLOCATE(matrix_a, matrix_b, matrix_c, matrix_d)

      ELSE IF (.NOT. xas_tdp_control%tamm_dancoff) THEN

!     Compute the submatrices A, G and the different kernels (B, C, D and E)
         ALLOCATE(tmp_mat)
         ALLOCATE(matrix_a)
         ALLOCATE(matrix_b)
         ALLOCATE(matrix_c)
         ALLOCATE(matrix_d)
         ALLOCATE(matrix_e)
         ALLOCATE(donor_state%metric(2))
         CALL build_gs_contribution(matrix_a, submat_dist, submat_blk_size, donor_state, xas_tdp_env, qs_env)
         CALL build_kernel_coulomb(matrix_b, submat_dist, submat_blk_size, donor_state, xas_tdp_env, qs_env)
         IF (fxc == xas_tdp_pade )CALL build_kernel_xc(matrix_c, submat_dist, submat_blk_size, &
                                                       donor_state, xas_tdp_env, qs_env)
         IF (fxc == xas_tdp_hfx) CALL build_kernel_exchange(matrix_d, matrix_e, submat_dist, &
                                                            submat_blk_size, donor_state, xas_tdp_env, qs_env)
         CALL build_metric(donor_state%metric, submat_dist, submat_blk_size, donor_state, qs_env, do_inv=.TRUE.)       

!     Build the auxiliary matrix (A-D+E)^0.5
         CALL dbcsr_copy(tmp_mat, matrix_a)
         IF (fxc == xas_tdp_hfx) THEN
            CALL dbcsr_add(tmp_mat, matrix_d, 1.0_dp, -1.0_dp)
            CALL dbcsr_add(tmp_mat, matrix_e, 1.0_dp, 1.0_dp)
         END IF
         
!     Warn the user if the matrix is not positive semi-definite and the procedure not sound
         CALL arnoldi_extremal(tmp_mat, max_ev, min_ev, threshold=1.0E-4_dp, max_iter=50, converged=conv)         
         IF (min_ev < 0.0_dp .AND. ABS(min_ev) > 1.0E-6_dp) THEN
            CPWARN("The full TDDFT problem might not have been soundly turned Hermitian. Try TDA.")
         END IF

         ALLOCATE(donor_state%matrix_aux)
         CALL dbcsr_create(matrix=donor_state%matrix_aux, matrix_type="S", dist=submat_dist, &
                           row_blk_size=submat_blk_size, col_blk_size=submat_blk_size, name="MAT_AUX")

!        Go through fm to get the square root. Iterative procedures have failed in tests
         CALL get_qs_env(qs_env, blacs_env=blacs_env, para_env=para_env)
         CALL dbcsr_get_info(tmp_mat, nfullrows_total=nrow) 
         CALL cp_fm_struct_create(fm_struct,context=blacs_env, nrow_global=nrow, para_env=para_env,&
                                  ncol_global=nrow)
         CALL cp_fm_create(fm_tmp, fm_struct)
         CALL cp_fm_create(work, fm_struct)

         CALL copy_dbcsr_to_fm(tmp_mat, fm_tmp)
         CALL cp_fm_power(fm_tmp, work, 0.5_dp, 1.0E-8_dp, n_dependent)
         CALL copy_fm_to_dbcsr(fm_tmp, donor_state%matrix_aux)

         CALL dbcsr_finalize(donor_state%matrix_aux)

         CALL cp_fm_release(work)
         CALL cp_fm_release(fm_tmp)
         CALL cp_fm_struct_release(fm_struct)

!     Create the main matrix_tdp as:
!     M = (A-D+E)^0.5 * G^-1 * (A+4B+2C-E-D) * G^-1 * (A-D+E)^0.5
!     Only pass the middle sum as matrix_tdp, the product is done in the solve routine
         CALL dbcsr_create(matrix=matrix_tdp, name="XAS TDP MATRIX",matrix_type="S",dist=submat_dist,&
                           row_blk_size=submat_blk_size, col_blk_size=submat_blk_size) 

         !the central matrix
         CALL dbcsr_copy(matrix_tdp, matrix_a)
         CALL dbcsr_add(matrix_tdp, matrix_b, 1.0_dp, 4.0_dp)
         IF (fxc == xas_tdp_pade) CALL dbcsr_add(matrix_tdp, matrix_c, 1.0_dp, 2.0_dp)
         If (fxc == xas_tdp_hfx) CALL dbcsr_add(matrix_tdp, matrix_d, 1.0_dp, -1.0_dp)
         If (fxc == xas_tdp_hfx) CALL dbcsr_add(matrix_tdp, matrix_e, 1.0_dp, -1.0_dp)


!     Some clean-up
         CALL dbcsr_release(matrix_a)
         CALL dbcsr_release(matrix_b)
         CALL dbcsr_release(matrix_c)
         CALL dbcsr_release(matrix_d)
         CALL dbcsr_release(matrix_e)
         CALL dbcsr_release(tmp_mat)
         DEALLOCATE(matrix_a, matrix_b, matrix_c, matrix_d, matrix_e, tmp_mat)

      END IF ! Tamm_Dancoff

!  Finalize
      CALL dbcsr_finalize(matrix_tdp)

!  General Clean-up
      CALL dbcsr_distribution_release(submat_dist)
      DEALLOCATE(submat_dist, submat_blk_size)

      CALL timestop(handle)

   END SUBROUTINE setup_xas_tdp_full_prob

! **************************************************************************************************
!> \brief Solves the XAS TDP generalized eigenvalue problem omega*C = matrix_tdp*C using standard
!>        full diagonalization methods (for non-hermitian matrices), for a given donor state
!> \param donor_state ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note The computed eigenvalues and eigenvectors are stored in the donor_state
!>       The eigenvectors are the LR-coefficients. In case of TDA, c^- is stored. In the general
!>       case, the sum c^+ + c^- is stored
! **************************************************************************************************
   SUBROUTINE solve_xas_tdp_full_prob(donor_state, xas_tdp_control, qs_env)  

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "solve_xas_tdp_full_prob", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: nrow, i, nao, nstates, imo, &
                                                         to_remove, ndo_mo, nelectron, istate
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(dbcsr_type), POINTER                       :: matrix_tdp, tmp_mat 
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      INTEGER                                         :: handle
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct, gs_struct
      TYPE(cp_fm_type), POINTER                       :: rhs_matrix, lhs_matrix, work, c_diff, &
                                                         c_sum
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: scaling
      REAL(dp), DIMENSION(:), POINTER                 :: tmp_evals

      CALL timeset(routineN, handle)

      NULLIFY(para_env, blacs_env, fm_struct, rhs_matrix, matrix_tdp, lhs_matrix, work, tmp_mat)
      NULLIFY(c_diff, c_sum, gs_struct, tmp_evals)

      matrix_tdp => donor_state%matrix_tdp
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env) 

!     Initialization
      CALL dbcsr_get_info(matrix_tdp, nfullrows_total=nrow)
      CALL cp_fm_struct_create(fm_struct,context=blacs_env, nrow_global=nrow, para_env=para_env,&
                               ncol_global=nrow)

!     Need to go to full matrices to diagonalize
      CALL cp_fm_create(rhs_matrix, fm_struct)
         
!     Allocating space for solutions
      ALLOCATE(donor_state%lr_evals(nrow))
      CALL cp_fm_create(c_sum, fm_struct)
      CALL cp_fm_create(work, fm_struct)

!     Allocate the scaling vector used for normalization
      ALLOCATE(scaling(nrow))

!     Test on TDA
      IF (xas_tdp_control%tamm_dancoff) THEN

!        Get the main matrix as fm
         CALL copy_dbcsr_to_fm(matrix_tdp, rhs_matrix)
         
!        Get the metric as a fm
         CALL cp_fm_create(lhs_matrix, fm_struct)
         CALL copy_dbcsr_to_fm(donor_state%metric(1)%matrix, lhs_matrix)

!        Diagonalisation (Cholesky decomposition)
         CALL cp_fm_geeig(rhs_matrix, lhs_matrix, c_sum, donor_state%lr_evals, work)

!        TDA specific clean-up
         CALL cp_fm_release(lhs_matrix)

!        Normalize. c_sum is just c^- for TDA
         scaling = 0.0_dp
         CALL get_normal_scaling(scaling, c_sum, donor_state)
         CALL cp_fm_column_scale(c_sum, scaling)

      ELSE ! not TDA

!        Make the problem Hermitian by solving M*z = omega^2*z with 
!        M = (A-D+E)^0.5 * G^-1 * (A+4B+2C-E-D) * G^-1 * (A-D+E)^0.5
         ALLOCATE(tmp_mat)
         CALL dbcsr_create(matrix=tmp_mat, template=donor_state%matrix_aux, matrix_type="N")

!        First compute G^-1 * (A+4B+2C-E-D) * G^-1 and store it in matrix_tdp
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_tdp, donor_state%metric(2)%matrix, 0.0_dp, tmp_mat)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, tmp_mat, 0.0_dp, matrix_tdp)

!        Compute the full matrix M and store it in tmp_mat
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_tdp, donor_state%matrix_aux, 0.0_dp, tmp_mat)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%matrix_aux, tmp_mat, 0.0_dp, tmp_mat)

!        Get the matrix M as fm
         CALL copy_dbcsr_to_fm(tmp_mat, rhs_matrix)

!        The problemed has been turned Hermitian. Call stendard routine
         CALL cp_fm_syevd(rhs_matrix, work, donor_state%lr_evals)

!        Filter out residual negative (and tiny) eigenvalues
         WHERE (donor_state%lr_evals < 1.0E-6_dp) donor_state%lr_evals = 0.0_dp

!        Take the square root since we computed omega^2
         donor_state%lr_evals = SQRT(donor_state%lr_evals)

!        Need to normalize the eigenvectors. Need to retrieve the difference of the LR coeffs:
!        (c^+ - c^-) = 1/omega^2 * G^-1 * (A+4B+2C-E-D) * G^-1 * (A-D+E)^0.5 * z
         CALL cp_fm_create(c_diff, fm_struct)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_tdp, donor_state%matrix_aux, 0.0_dp, tmp_mat)
         CALL cp_dbcsr_sm_fm_multiply(tmp_mat, work, c_diff, ncol=nrow)

         scaling = 0.0_dp
         WHERE (donor_state%lr_evals > 0.0_dp) scaling = 1.0_dp/donor_state%lr_evals**2
         CALL cp_fm_column_scale(c_diff, scaling)

!        The normalization is c_diff * G * c_diff = +- 1
         scaling = 0.0_dp
         CALL get_normal_scaling(scaling, c_diff, donor_state)
         CALL cp_fm_column_scale(c_diff, scaling)

!        Get the sum c^+ + c^- = -1/omega G^-1 * (A-D+E)*(c^+ - c^-)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%matrix_aux, donor_state%matrix_aux, &
                             0.0_dp, tmp_mat)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, donor_state%metric(2)%matrix, tmp_mat, 0.0_dp, tmp_mat)
         CALL cp_dbcsr_sm_fm_multiply(tmp_mat, c_diff, c_sum, ncol=nrow)
         scaling = 0.0_dp
         WHERE (donor_state%lr_evals > 0.0_dp) scaling = -1.0_dp/donor_state%lr_evals
         CALL cp_fm_column_scale(c_sum, scaling)

!        Some clean-up
         CALL dbcsr_release(tmp_mat)
         CALL cp_fm_release(c_diff)
         DEALLOCATE(tmp_mat)

      END IF ! TDA

!  Get rid of all the states with zero eigenvalues because of the projection on the unoccupied KS
!  orbitals. Need to remove the number of GS MO time the number of donor states
      CALL cp_fm_get_info(c_sum, ncol_global=nstates)
      ndo_mo = SIZE(donor_state%mo_indices)  
      CALL get_qs_env(qs_env, nelectron_total=nelectron)
      to_remove = ndo_mo*nelectron/2
      ALLOCATE(tmp_evals(nstates-to_remove))

      !evals are ordered
      tmp_evals(:) = donor_state%lr_evals(to_remove+1:nstates)
      DEALLOCATE(donor_state%lr_evals)
      donor_state%lr_evals => tmp_evals

!  Reorganize the eigenvectors in array of cp_fm_p_types so that each element corresponds to an 
!  excited state. Makes later calls to those easier and more efficient
      CALL cp_fm_get_info(donor_state%gs_coeffs, matrix_struct=gs_struct, nrow_global=nao)      
      ALLOCATE(donor_state%lr_coeffs(nstates-to_remove))

      DO i = 1,nstates-to_remove
         istate = to_remove + i

         CALL cp_fm_create(donor_state%lr_coeffs(i)%matrix, gs_struct)

         DO imo = 1, ndo_mo
            
            CALL cp_fm_to_fm_submat(msource=c_sum, mtarget=donor_state%lr_coeffs(i)%matrix, &
                                    nrow=nao, ncol=1, s_firstrow=(imo-1)*nao+1, s_firstcol=istate, &
                                    t_firstrow=1, t_firstcol=imo)

         END DO !imo
      END DO !istate
         
!  Clean-up
      DEALLOCATE(scaling)
      CALL cp_fm_release(rhs_matrix)
      CALL cp_fm_release(work)
      CALL cp_fm_release(c_sum)
      CALL cp_fm_struct_release(fm_struct)

!  Perform a partial clean-up of the donor_state
      CALL dbcsr_release(matrix_tdp)
      CALL dbcsr_release(donor_state%metric(1)%matrix)
      IF (.NOT. xas_tdp_control%tamm_dancoff) THEN
         CALL dbcsr_release(donor_state%matrix_aux)
         CALL dbcsr_release(donor_state%metric(2)%matrix)
      END IF

      CALL timestop(handle)

   END SUBROUTINE solve_xas_tdp_full_prob

! **************************************************************************************************
!> \brief Returns the scaling to apply to normalize the LR eigenvectors. 
!> \param scaling the scaling array to apply 
!> \param lr_coeffs the linear response coefficients as a fm
!> \param donor_state ...
!> \note The LR coeffs are normalized when c^T G c = +- 1, G is the metric, c = c^- for TDA and
!>       c = c^+ - c^- for the full problem
! **************************************************************************************************
   SUBROUTINE get_normal_scaling(scaling, lr_coeffs, donor_state)

      REAL(dp), DIMENSION(:), ALLOCATABLE             :: scaling
      TYPE(cp_fm_type), POINTER                       :: lr_coeffs
      TYPE(donor_state_type), POINTER                 :: donor_state

      CHARACTER(len=*), PARAMETER :: routineN = "get_normal_scaling", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: nvals, nrow
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(cp_fm_type), POINTER                       :: work, fm_norm
      TYPE(cp_fm_struct_type), POINTER                :: norm_struct, work_struct

      NULLIFY(para_env, blacs_env, norm_struct, work, fm_norm, work_struct)

!  Creating the matrix structures and initializing the work matrices
      CALL cp_fm_get_info(lr_coeffs, context=blacs_env, para_env=para_env, &
                          matrix_struct=work_struct, ncol_global=nvals, nrow_global=nrow)
      CALL cp_fm_struct_create(norm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=nvals, ncol_global=nvals)

      CALL cp_fm_create(work, work_struct)
      CALL cp_fm_create(fm_norm, norm_struct)


!  Taking c^T * G * C
      CALL cp_dbcsr_sm_fm_multiply(donor_state%metric(1)%matrix, lr_coeffs, work, ncol=nvals)
      CALL cp_fm_gemm('T', 'N', nvals, nvals, nrow, 1.0_dp, lr_coeffs, work, 0.0_dp, fm_norm)

!  Computing the needed scaling
      CALL cp_fm_get_diag(fm_norm, scaling)
      WHERE(scaling .NE. 0.0_dp) scaling = 1.0_dp/SQRT(ABS(scaling))

!  Clean-up
      CALL cp_fm_release(work)
      CALL cp_fm_release(fm_norm)
      CALL cp_fm_struct_release(norm_struct)

   END SUBROUTINE get_normal_scaling


! **************************************************************************************************
!> \brief Computes the oscillator strength based on the dipole moment (velocity or length rep) for
!>        all available excitation energies and store the results in the donor_state
!> \param donor_state the donor state which is excited
!> \param xas_tdp_control ...
!> \param xas_tdp_env ... 
!> \note The oscillator strength is a scalar: osc_str = 2/3*omega*(dipole_r)^2 = 2/(3*omega)*(dipole_v)^2
!>       The formulae for the dipoles come from the trace of the dipole operator with the transition
!>       densities, i.e. what we get from solving the xas_tdp problem. Same procedure with or wo TDA
! **************************************************************************************************
   SUBROUTINE compute_dipole_fosc(donor_state, xas_tdp_control, xas_tdp_env)

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = "compute_dipole_fosc", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: nosc, iosc, ndo_mo, j, nao
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: dipmat
      TYPE(cp_fm_type), POINTER                       :: col_work, mat_work, lr_coeffs
      TYPE(cp_fm_struct_type), POINTER                :: col_struct, mat_struct
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: diag
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      REAL(dp)                                        :: dip_mom
      REAL(dp), DIMENSION(:), POINTER                 :: osc_str

      NULLIFY(dipmat, col_work, mat_work, col_struct, mat_struct, para_env, blacs_env, lr_coeffs, osc_str)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices)
      nosc = SIZE(donor_state%lr_evals)
      ALLOCATE(donor_state%osc_str(nosc))
      osc_str => donor_state%osc_str
      osc_str = 0.0_dp

!     For spin-unpolarized ground state, only have transition density if singlet
      IF (xas_tdp_control%do_singlet) THEN

         dipmat => xas_tdp_env%dipmat
         !compute the prefactor
         IF (xas_tdp_control%dipole_form == xas_dip_len) THEN
            osc_str = 2.0_dp/3.0_dp*donor_state%lr_evals
         ELSE
            osc_str = 2.0_dp/3.0_dp/donor_state%lr_evals
         END IF

         ! do some work matrix initialization
         CALL cp_fm_get_info(donor_state%gs_coeffs, matrix_struct=col_struct, para_env=para_env, &
                             context=blacs_env, nrow_global=nao)
         CALL cp_fm_struct_create(mat_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=ndo_mo, ncol_global=ndo_mo)
         CALL cp_fm_create(col_work, col_struct)
         CALL cp_fm_create(mat_work, mat_struct)
         ALLOCATE(diag(ndo_mo))

         DO iosc = 1,nosc
            
            dip_mom = 0.0_dp
            !looping over the cartesian coordinates
            DO j = 1,3

               lr_coeffs => donor_state%lr_coeffs(iosc)%matrix
               !computing c^0 * dipmat(j) * c^i
               CALL cp_dbcsr_sm_fm_multiply(dipmat(j)%matrix, lr_coeffs, col_work, ncol=ndo_mo)
               CALL cp_fm_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, donor_state%gs_coeffs, &
                               col_work, 0.0_dp, mat_work)

               !get the sum of the dipole moments
               CALL cp_fm_get_diag(mat_work, diag)
               dip_mom = dip_mom + SUM(diag)**2

            END DO !j
            
            !factor 2 because of singlet definition u = 0.5*(c_alpha + c_beta)
            osc_str(iosc) = 2*osc_str(iosc)*dip_mom

         END DO !iosc

         !clean-up
         CALL cp_fm_release(mat_work)
         CALL cp_fm_release(col_work)
         CALL cp_fm_struct_release(mat_struct)
         DEALLOCATE(diag)

      END IF !singlet

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
!> \brief Builds the matrix containing the ground state contribution to the matrix_tdp (aka matrix A)
!>         => A_{pi,qj} = sum_rs Q_pr* (F_rs*delta_ij - epsilon_ij*S_rs) Q_sq^T, where:
!>         Q = 1 - SP is the projector on the unoccupied unperturbed ground state
!>         F is the KS matrix
!>         S is the overlap matrix
!> \param matrix_a  pointer to a DBCSR matrix containing A
!> \param dist_a the dbcsr distribution of matrix A
!> \param blk_size_a the row/col block sizes of matrix A
!> \param donor_state ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE build_gs_contribution(matrix_a, dist_a, blk_size_a, donor_state, xas_tdp_env, qs_env)

      TYPE(dbcsr_type), POINTER                       :: matrix_a
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(dbcsr_distribution_type), POINTER          :: dist_a
      INTEGER, DIMENSION(:), POINTER                  :: blk_size_a
      
      CHARACTER(len=*), PARAMETER :: routineN = 'build_gs_contribution', &
                                     routineP = moduleN//":"//routineN
      
      INTEGER                                         :: ndo_mo, i, nblk_row, j, iblk, jblk, blk
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks, matrix_s
      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      INTEGER, DIMENSION(:), POINTER                  :: row_blk_size 
      REAL(dp), DIMENSION(:,:), POINTER               :: energy_evals     
      REAL(dp), DIMENSION(:), POINTER                 :: work_block
      LOGICAL                                         :: found_block
      TYPE(dbcsr_type)                                :: desym_matrix_s, desym_matrix_ks, &
                                                         work_matrix, tmp_matrix
      TYPE(dbcsr_iterator_type)                       :: iter

      NULLIFY(matrix_ks, dbcsr_dist, row_blk_size, work_block, energy_evals, matrix_s)

!  Note: matrix A is symmetric. However, build it from non-symmetric basis because dbcsr distributions
!        only match for non-symmetric matrices (that is when the same dist is replicated). So we go
!        over the top triangular part and then redistribute into a symmetric matrix

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices, 1)
      CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_get_info(matrix_ks(1)%matrix, row_blk_size=row_blk_size)
      nblk_row = SIZE(row_blk_size)
      energy_evals => donor_state%energy_evals

!  Creating the work matrix, which is not symmetric (as a type)
      CALL dbcsr_create(matrix=work_matrix, name="WORK MAT", matrix_type="N", dist=dist_a, &
                        row_blk_size=blk_size_a, col_blk_size=blk_size_a)

!  Desymmetrize matrices S and KS
      CALL dbcsr_create(matrix=desym_matrix_s, name="DESYM MATRIX S", matrix_type="N", &
                        dist=dbcsr_dist, row_blk_size=row_blk_size, col_blk_size=row_blk_size)
      CALL dbcsr_desymmetrize(matrix_s(1)%matrix, desym_matrix_s)
   
      CALL dbcsr_create(matrix=desym_matrix_ks, name="DESYM MATRIX KS", matrix_type="N", &
                        dist=dbcsr_dist, row_blk_size=row_blk_size, col_blk_size=row_blk_size)
      CALL dbcsr_desymmetrize(matrix_ks(1)%matrix, desym_matrix_ks)

!  multiply S and KS by proj_Q from the left (keep the names)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, xas_tdp_env%q_projector, desym_matrix_s, 0.0_dp, desym_matrix_s)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, xas_tdp_env%q_projector, desym_matrix_ks, 0.0_dp, desym_matrix_ks)

!  multiply by the transpose of proj_Q from the right
      CALL dbcsr_multiply('N', 'T', 1.0_dp, desym_matrix_s, xas_tdp_env%q_projector, 0.0_dp, desym_matrix_s)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, desym_matrix_ks, xas_tdp_env%q_projector, 0.0_dp, desym_matrix_ks)


!  Loop over the blocks of desymmetrized matrix KS and put them in the work matrix
      CALL dbcsr_iterator_start(iter, desym_matrix_ks)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

!        Get the block
         found_block = .FALSE.
         CALL dbcsr_get_block_p(desym_matrix_ks, iblk, jblk, work_block, found_block)

         IF (found_block) THEN
   
!           The KS matrix only appears on diagonal of matrix A => loop over II donor MOs
            DO i = 1,ndo_mo

!              Put the block as it is
               CALL dbcsr_put_block(work_matrix,(i-1)*nblk_row+iblk, (i-1)*nblk_row+jblk, work_block)

            END DO ! i

         END IF ! found_block
         NULLIFY(work_block)

      END DO ! iteration on KS blocks
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(work_matrix)

!  Put the current work matrix into a temporary storage matrix
      CALL dbcsr_copy(tmp_matrix, work_matrix, name="TMP MATRIX")

!  Reinitialize the work matrix
      CALL dbcsr_release(work_matrix)
      CALL dbcsr_create(matrix=work_matrix, name="WORK MAT", matrix_type="N", dist=dist_a, & 
                        row_blk_size=blk_size_a, col_blk_size=blk_size_a)      

!  Iterate over the desymmetrized matrix S
      CALL dbcsr_iterator_start(iter, desym_matrix_s)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)                                    

!        Get the block
         found_block = .FALSE.
         CALL dbcsr_get_block_p(desym_matrix_s, iblk, jblk, work_block, found_block) 

         IF (found_block) THEN

!           The S matrix appears everywhere, need to loop over all donor MOs
!           Since the end matrix is symmetric, keep to upper triangle
            DO i = 1,ndo_mo
               DO j = i,ndo_mo

!                 Put epsilon_ij*S_pq (if non-zero eval)
                  IF (energy_evals(i,j) .NE. 0.0_dp) THEN
                     CALL dbcsr_put_block(work_matrix,(i-1)*nblk_row+iblk, (j-1)*nblk_row+jblk, & 
                                          energy_evals(i,j)*work_block)   
                  END IF

               END DO ! j
            END DO ! i

         END IF ! found_block
         NULLIFY(work_block)

      END DO ! iteration over S
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(work_matrix)
      CALL dbcsr_release(desym_matrix_s)
      CALL dbcsr_release(desym_matrix_ks) 

!  Substract the current work matrix to the tmp matrix to get the full F_pq*delta_ij - epsilon_ij*S_pq
      CALL dbcsr_add(tmp_matrix, work_matrix, 1.0_dp, -1.0_dp)
      CALL dbcsr_finalize(tmp_matrix)
      CALL dbcsr_release(work_matrix)

!  At this point the content of tmp_matrix is that of matrix_a, but its type is normal. Create a
!  symmetric matrix A and redistribute the tmp matrix in there
      CALL dbcsr_create(matrix=matrix_a, name="MATRIX A", matrix_type="S", dist=dist_a, &         
                        row_blk_size=blk_size_a, col_blk_size=blk_size_a)
      CALL dbcsr_complete_redistribute(tmp_matrix, matrix_a)
      CALL dbcsr_release(tmp_matrix)
      
!  Finalize
      CALL dbcsr_finalize(matrix_a)   
   
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

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_g
      TYPE(dbcsr_distribution_type), POINTER          :: dist_g
      INTEGER, DIMENSION(:), POINTER                  :: blk_size_g
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      LOGICAL, OPTIONAL, INTENT(IN)                   :: do_inv

      CHARACTER(len=*), PARAMETER :: routineN = 'build_metric', routineP = moduleN//":"//routineN
      
      INTEGER                                         :: i, iblk, jblk, ndo_mo, nblk_row, blk
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      INTEGER, DIMENSION(:), POINTER                  :: row_blk_size
      REAL(dp), DIMENSION(:), POINTER                 :: work_block
      LOGICAL                                         :: found_block
      TYPE(dbcsr_iterator_type)                       :: iter     
      TYPE(dbcsr_type), POINTER                       :: matrix_sinv
      LOGICAL                                         :: my_do_inv

      NULLIFY(matrix_s, row_blk_size, work_block)

!  Initilization
      ndo_mo = SIZE(donor_state%mo_indices, 1)
      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size)
      nblk_row = SIZE(row_blk_size)
      my_do_inv = .FALSE.
      IF (PRESENT(do_inv)) my_do_inv = do_inv

!  Creating the symmetric  matrices G and G^-1 with the right size and distribution 
      ALLOCATE(matrix_g(1)%matrix)
      CALL dbcsr_create(matrix=matrix_g(1)%matrix, name="MATRIX G", matrix_type="S", dist=dist_g, &
                        row_blk_size=blk_size_g, col_blk_size=blk_size_g)

!  Fill the matrices G by looping over the block of S and putting them on the diagonal
      CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

!        Get the block
         found_block = .FALSE.
         CALL dbcsr_get_block_p(matrix_s(1)%matrix, iblk, jblk, work_block, found_block)

         IF (found_block .AND. jblk>=iblk) THEN

!           Go over the diagonal of G => donor MOs ii
            DO i = 1,ndo_mo
               CALL dbcsr_put_block(matrix_g(1)%matrix, (i-1)*nblk_row+iblk, (i-1)*nblk_row+jblk, work_block)
            END DO

         END IF
         NULLIFY(work_block)

      END DO ! dbcsr_iterator
      CALL dbcsr_iterator_stop(iter)

!  Finalize
      CALL dbcsr_finalize(matrix_g(1)%matrix)

!  If the inverse of G is required, do the same as above with the inverse
      IF (my_do_inv) THEN
         
         CPASSERT(SIZE(matrix_g) == 2)

         ! Create the matrix
         ALLOCATE(matrix_g(2)%matrix)
         CALL dbcsr_create(matrix=matrix_g(2)%matrix, name="MATRIX GINV", matrix_type="S", &
                           dist=dist_g, row_blk_size=blk_size_g, col_blk_size=blk_size_g)  

         ! Invert the overlap matrix
         ALLOCATE(matrix_sinv)
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

            IF (found_block .AND. jblk>=iblk) THEN

!              Go over the diagonal of G => donor MOs ii
               DO i = 1,ndo_mo
                  CALL dbcsr_put_block(matrix_g(2)%matrix, (i-1)*nblk_row+iblk, (i-1)*nblk_row+jblk, work_block)
               END DO

            END IF
            NULLIFY(work_block)

         END DO ! dbcsr_iterator
         CALL dbcsr_iterator_stop(iter)

         !  Finalize
         CALL dbcsr_finalize(matrix_g(2)%matrix)

         !  Clean-up
         CALL dbcsr_release(matrix_sinv)
         DEALLOCATE(matrix_sinv)
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

      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(qs_environment_type), POINTER              :: qs_env 

      CHARACTER(len=*), PARAMETER :: routineN = "compute_ri_coulomb2_int", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: ikind, i, nsgf
      TYPE(gto_basis_set_type), POINTER               :: ri_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set        
      INTEGER, DIMENSION(:), POINTER                  :: ex_kinds
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: tmp_int

      NULLIFY(ri_basis, qs_kind_set)
      CPASSERT(ASSOCIATED(xas_tdp_env))

!  Initialization
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      ex_kinds => xas_tdp_env%ex_kind_indices
      ALLOCATE(xas_tdp_env%ri_inv_mats(SIZE(ex_kinds)))

!  Loop over the excited kinds
      DO i = 1,SIZE(ex_kinds)
         ikind = ex_kinds(i) 

!        Get the RI basis of interest and its quantum numbers
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_XAS")
         CALL get_gto_basis_set(ri_basis, nsgf=nsgf) 

!        The 2D array with the integrals for this kind
         ALLOCATE(tmp_int(nsgf, nsgf))
         tmp_int = 0.0_dp

         CALL int_operators_r12_ab_os(r12_operator=operator_coulomb, vab=tmp_int, &
                                      rab=(/0.0_dp,0.0_dp,0.0_dp/), fba=ri_basis, fbb=ri_basis, &
                                      calculate_forces=.FALSE.)

!        Inverting and storing in the xas_tdp_env
         CALL invmat_symm(tmp_int)  

         ALLOCATE(xas_tdp_env%ri_inv_mats(i)%array(nsgf,nsgf))
         xas_tdp_env%ri_inv_mats(i)%array(:,:) = tmp_int(:,:)

!        clean-up before next kind
         DEALLOCATE(tmp_int)

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

      INTEGER, INTENT(IN)                             :: exat_index, exat_kind
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_coulomb3_int', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                         :: nsetp, nsetb, natom, &
                                                         nao, nseta, bat, aset, bset, pset, &
                                                         sgfa, sgfp, sgfb, aat, akind, bkind, & 
                                                         nkind, igp, i, ncoa, ncob, ncop, ipgf, &
                                                         startp, endp
      REAL(dp)                                        :: rab2, rap2, rbp2, rpgfp
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER               :: p_basis, b_basis, a_basis
      INTEGER, DIMENSION(:), POINTER                  :: lp_max, lp_min, npgfp, lb_max,lb_min,&
                                                         npgfb, la_max, la_min, npgfa, nsgfa, &
                                                         nsgfb, nsgfp
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgfa, first_sgfb, first_sgfp
      REAL(dp), DIMENSION(:,:), POINTER               :: zetp, zetb, zeta, &
                                                         sphi_p, sphi_a, sphi_b
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos 
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE         :: int_abp_set, int_array, int_abp
      REAL(dp), DIMENSION(3)                          :: ra, rb, rp, rab
      INTEGER, DIMENSION(:), ALLOCATABLE              :: first_sgf_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), & 
                                          POINTER     :: sab_orb
      TYPE(neighbor_list_iterator_p_type), &                                                         
                            DIMENSION(:), POINTER     :: nl_iterator   
      TYPE(cell_type), POINTER                        :: cell
      TYPE(cp_para_env_type), POINTER                 :: para_env
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: gccc, f, rpgfa, rpgfb
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: vabc
      REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE       :: v

      NULLIFY(qs_kind_set, p_basis, b_basis, a_basis, particle_set, mos)
      NULLIFY(lp_max, lp_min, npgfp, zetp, sphi_p, first_sgfp)
      NULLIFY(lb_max, lb_min, npgfb, zetb, sphi_b, first_sgfb)
      NULLIFY(la_max, la_min, npgfa, zeta, sphi_a, first_sgfa)
      NULLIFY(sab_orb, nl_iterator, cell, basis_set_list, para_env)

!  Get all we need from the qs_env
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, mos=mos, &
                      sab_orb=sab_orb, cell=cell, natom=natom, para_env=para_env)
      CALL get_mo_set(mos(1)%mo_set, nao=nao)
      ALLOCATE(first_sgf_set(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf_set)

!  Retrieve the basis set P centered on the excited atom
      CALL get_qs_kind(qs_kind_set(exat_kind), basis_set=p_basis, basis_type="RI_XAS")

      CALL get_gto_basis_set(p_basis, lmax=lp_max, npgf=npgfp, zet=zetp, lmin=lp_min, nset=nsetp, &
                             nsgf_set=nsgfp, sphi=sphi_p, first_sgf=first_sgfp) 

!  Allocate the 3D array where results will be stored
      ALLOCATE(int_array(nao, nao, SUM(nsgfp)))
      int_array = 0.0_dp

!  Create lists of basis for atoms a and b
      nkind = SIZE(qs_kind_set)
      ALLOCATE(basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

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
         rab2 = DOT_PRODUCT(rb-ra,rb-ra)
         rap2 = DOT_PRODUCT(rp-ra,rp-ra)
         rbp2 = DOT_PRODUCT(rp-rb,rp-rb) 

!        Loop over sets of the 3 basis
         DO aset = 1,nseta
            DO bset = 1,nsetb
               DO pset = 1,nsetp

!                 Some initialization
                  sgfa = first_sgfa(1, aset) 
                  sgfb = first_sgfb(1, bset)
                  sgfp = first_sgfp(1, pset)   
                  ncoa = npgfa(aset)*ncoset(la_max(aset))
                  ncob = npgfb(bset)*ncoset(lb_max(bset))
                  ncop = npgfp(pset)*ncoset(lp_max(pset))

!                 Allocate array for the current comnbination of sets
                  ALLOCATE(int_abp_set(nsgfa(aset), nsgfb(bset), nsgfp(pset)))
                  int_abp_set = 0.0_dp

!                 Allocate some work arrays
                  ALLOCATE(f(0:la_max(aset)+lb_max(bset)+lp_max(pset)+2))
                  ALLOCATE(v(ncoa, ncob, ncop, la_max(aset)+lb_max(bset)+lp_max(pset)+1))
                  ALLOCATE(gccc(ncop))
                  ALLOCATE(vabc(ncoa, ncob))
                  ALLOCATE(int_abp(ncoa, ncob, ncop))
                  f = 0.0_dp
                  v = 0.0_dp
                  gccc = 0.0_dp
                  vabc = 0.0_dp
                  int_abp = 0.0_dp

!                 hack: set the rpgf to avoid screening out
                  ALLOCATE(rpgfa(npgfa(aset)))
                  ALLOCATE(rpgfb(npgfb(bset)))
                  rpgfa = 1.0E10_dp
                  rpgfb = 1.0E10_dp
                  rpgfp = 1.0E10_dp

!                 Go one pgf at a time because of how coulomb3 is designed
                  endp = 0
                  DO ipgf = 1, npgfp(pset)
                     startp = endp + ncoset(lp_min(pset)-1) + 1
                     endp = endp + ncoset(lp_max(pset))
                     CALL coulomb3(la_max(aset), npgfa(aset), zeta(:,aset), rpgfa(:), la_min(aset), &
                                   lb_max(bset), npgfb(bset), zetb(:,bset), rpgfb(:), lb_min(bset), &
                                   lp_max(pset), zetp(ipgf,pset), rpgfp, lp_min(pset), gccc, &
                                   rb-ra, rab2, rp-ra, rap2, rbp2, vabc, &
                                   int_abp(:,:,startp:endp), v, f)
                  END DO !ipgf

!                 Contract into sgf
                  CALL abc_contract(int_abp_set, int_abp, sphi_a(:,sgfa:), sphi_b(:,sgfb:), &
                                    sphi_p(:,sgfp:), ncoa, ncob, ncop, nsgfa(aset), nsgfb(bset), &
                                    nsgfp(pset)) 

!                 Put it at the right place in the final array
                  int_array(first_sgf_set(aat)+sgfa-1:first_sgf_set(aat)+sgfa+nsgfa(aset)-2, &
                            first_sgf_set(bat)+sgfb-1:first_sgf_set(bat)+sgfb+nsgfb(bset)-2, &
                            sgfp:sgfp+nsgfp(pset)-1) = int_abp_set

!                 In the neighbor list iterator, only get a,b pair once. But (ab|P) = (ba|P) 
                  DO i = 1,nsgfp(pset)
                     igp = sgfp -1 +i
                     int_array(first_sgf_set(bat)+sgfb-1:first_sgf_set(bat)+sgfb+nsgfb(bset)-2, &
                               first_sgf_set(aat)+sgfa-1:first_sgf_set(aat)+sgfa+nsgfa(aset)-2, &       
                               igp) = TRANSPOSE(int_abp_set(:,:,i))
                  END DO

!                 Deallocate before next loop
                  DEALLOCATE(int_abp_set, int_abp, f, v, gccc, vabc, rpgfa, rpgfb)

               END DO ! pset
            END DO ! bset
         END DO !aset

      END DO !neighbor iteration on a,b
      CALL neighbor_list_iterator_release(nl_iterator)

!  Updating the xas_tdp_env with the new integrals
      IF (ASSOCIATED(xas_tdp_env%ri_coul3)) DEALLOCATE(xas_tdp_env%ri_coul3)
      ALLOCATE(xas_tdp_env%ri_coul3(nao, nao, SUM(nsgfp)))
      xas_tdp_env%ri_coul3 = int_array

!  Clean-up
      DEALLOCATE(basis_set_list)

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
      
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: contr_int
      INTEGER, DIMENSION(:)                           :: domo_index
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'contract_int_kernel', &
                                     routineP = moduleN//":"//routineN

      LOGICAL                                         :: contract_two
      INTEGER                                         :: nao, nsgfp, i
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      TYPE(cp_fm_type), POINTER                       :: mo_coeff 
      REAL(dp), DIMENSION(:,:,:), POINTER             :: integrals
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: contr_once, contr_twice, coeffs

      NULLIFY(mos, mo_coeff, integrals)

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
      nsgfp = SIZE(integrals,3)
      ALLOCATE(contr_once(nao,nsgfp))
      contr_once = 0.0_dp
      IF (contract_two) ALLOCATE(contr_twice(1,nsgfp)) ! not efficient but not there to stay

!  Get the MO coeffs for the first contraction
      ALLOCATE(coeffs(nao,1))
      CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=coeffs, start_row=1, start_col=domo_index(1), &
                               n_rows=nao, n_cols=1, transpose=.FALSE.) 

      DO i = 1,nao
         CALL daxpy(nao*nsgfp, coeffs(i,1), integrals(:,i,:), 1, contr_once(:,:), 1)
      END DO

!  Do the second contraction if requested
      IF (contract_two) THEN
         contr_twice = 0.0_dp
         CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=coeffs, start_row=1, &
                                  start_col=domo_index(2), n_rows=nao, n_cols=1, transpose=.FALSE.)

         DO i = 1,nao
            CALL daxpy(nsgfp, coeffs(i,1), contr_once(i,:), 1, contr_twice(1,:), 1)
         END DO
      END IF

!  Allocate and assign the contr_int array 
      IF (.NOT. contract_two) THEN
         ALLOCATE(contr_int(nao,nsgfp))
         contr_int(:,:) = contr_once(:,:)
      ELSE
         ALLOCATE(contr_int(1,nsgfp))
         contr_int(:,:) = contr_twice(:,:)
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

      TYPE(dbcsr_type), POINTER                       :: coulomb_kernel
      TYPE(dbcsr_distribution_type), POINTER          :: dist 
      INTEGER, DIMENSION(:), POINTER                  :: blk_size  
      TYPE(donor_state_type), POINTER                 :: donor_state 
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(qs_environment_type), POINTER              :: qs_env      

      CHARACTER(len=*), PARAMETER :: routineN = "build_kernel_coulomb", &  
                                     routineP = moduleN//":"//routineN  

      INTEGER                                         :: handle, ndo_mo, i, j, nsgfp, nao, iblk, &
                                                         jblk, nblk,  nrow_blk, ncol_blk
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: int_aI, int_bJ, proj_Q, work1, work2, &
                                                         aIJb_prod
      REAL(dp), DIMENSION(:,:), POINTER               :: PQ
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: context
      TYPE(cp_fm_type), POINTER                       :: work_fm     
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct        
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set      
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      INTEGER, DIMENSION(:), ALLOCATABLE              :: first_sgf

      CALL timeset(routineN, handle)

      NULLIFY(PQ, para_env, context, work_fm, fm_struct, qs_kind_set, particle_set)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices,1)                       
      CALL get_qs_env(qs_env, natom=nblk, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                      blacs_env=context, para_env=para_env)
      ALLOCATE(first_sgf(nblk))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf)

!     Get the (P|Q)^-1 from the ri_inv_mats
      i = locate(xas_tdp_env%ex_kind_indices, donor_state%kind_index)
      PQ => xas_tdp_env%ri_inv_mats(i)%array
      nsgfp = SIZE(PQ,1)
      nao = SIZE(xas_tdp_env%ri_coul3, 1)
      
!     Get the projector Q as an array (dbcsr => fm => array)
      ALLOCATE(proj_Q(nao,nao))
      CALL cp_fm_struct_create(fm_struct, context=context, para_env=para_env, nrow_global=nao, &
                               ncol_global=nao)
      CALL cp_fm_create(work_fm, fm_struct)
      CALL copy_dbcsr_to_fm(xas_tdp_env%q_projector, work_fm)
      CALL cp_fm_get_submatrix(fm=work_fm, target_m=proj_Q)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_release(work_fm)

!     Build the dbcsr matrix for the kernel
      CALL dbcsr_create(matrix=coulomb_kernel, name="COULOMB KERNEL", matrix_type='S', dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)

!     Allocate the work arrays to store the intermediate matrix products
      ALLOCATE(work1(nao, nsgfp))
      ALLOCATE(work2(nao,nao))
      ALLOCATE(aIJb_prod(nao,nao))

!  Loop over donor MOs and compute W_ab = (aI|P)* (P|Q)^-1 * (Q|Jb), then take projQ*W*projQ^T
!  Matrix is symmetric => only span the upper triangle
      DO i = 1,ndo_mo
         CALL contract_int_kernel(int_aI, (/donor_state%mo_indices(i)/), xas_tdp_env, qs_env)
         DO j = i,ndo_mo
            CALL contract_int_kernel(int_bJ, (/donor_state%mo_indices(j)/), xas_tdp_env, qs_env)

            work1 = 0.0_dp
            work2 = 0.0_dp
            aIJb_prod  = 0.0_dp

!           Do the product (aI|P)*(P|Q)^-1
            CALL dgemm('N', 'N', nao, nsgfp, nsgfp, 1.0_dp, int_aI, nao, PQ, nsgfp, 0.0_dp, &
                        work1, nao)

!           Mulitply by (Q|Jb)
            CALL dgemm('N', 'T', nao, nao, nsgfp, 1.0_dp, work1, nao, int_bJ, nao, 0.0_dp, &
                        aIJb_prod, nao)

!           Multiply by the projector
            CALL dgemm('N', 'N', nao, nao, nao,  1.0_dp, proj_Q, nao, aIJb_prod, nao, 0.0_dp, &
                       work2, nao)
            aIJb_prod = 0.0_dp
            CALL dgemm('N', 'T', nao, nao, nao,  1.0_dp, work2, nao, proj_Q, nao, 0.0_dp, & 
                        aIJb_prod, nao)    

!           Put it in the dbcsr matrix, for the right I,J and block by block
            DO iblk = 1,nblk
               DO jblk = 1,nblk

!                 If on the block diagonal, just take jblk >= iblk
                  IF (jblk < iblk .AND. i==j) CYCLE

!                 Size of the current block
                  nrow_blk = blk_size(iblk)
                  ncol_blk = blk_size(jblk)

                  CALL dbcsr_put_block(coulomb_kernel, (i-1)*nblk+iblk, (j-1)*nblk+jblk, &
                                       aIJb_prod(first_sgf(iblk):first_sgf(iblk)+nrow_blk-1, &
                                                  first_sgf(jblk):first_sgf(jblk)+ncol_blk-1))

               END DO ! jblk
            END DO ! iblk
            
            DEALLOCATE(int_bJ)
         END DO ! donor MO j
         DEALLOCATE(int_aI)
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
   SUBROUTINE build_kernel_exchange(ex_ker_diag, ex_ker_offdiag , dist, blk_size, donor_state, &
                                    xas_tdp_env, qs_env)

      TYPE(dbcsr_type), POINTER                       :: ex_ker_diag
      TYPE(dbcsr_type), POINTER, OPTIONAL             :: ex_ker_offdiag  
      TYPE(dbcsr_distribution_type), POINTER          :: dist     
      INTEGER, DIMENSION(:), POINTER                  :: blk_size                
      TYPE(donor_state_type), POINTER                 :: donor_state       
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(qs_environment_type), POINTER              :: qs_env    

      CHARACTER(len=*), PARAMETER :: routineN = "build_kernel_exchange", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, ndo_mo, i, j, nao, nsgfp, nblk, k,&
                                                         iblk, jblk, nrow_blk, ncol_blk
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set     
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      INTEGER, DIMENSION(:), ALLOCATABLE              :: first_sgf          
      TYPE(cp_blacs_env_type), POINTER                :: context
      TYPE(cp_para_env_type), POINTER                 :: para_env
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: proj_Q, int_aJ, int_bI, int_IJ, work1,&
                                                         work2, work3, abIJ_prod, aJIb_prod, work4
      REAL(dp), DIMENSION(:,:), POINTER               :: PQ
      TYPE(cp_fm_type), POINTER                       :: work_fm                      
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct      
      LOGICAL                                         :: do_offdiag

      CALL timeset(routineN, handle)

      NULLIFY(particle_set, qs_kind_set, para_env, context, work_fm, fm_struct, PQ)

!  Initialization
      do_offdiag = .FALSE.
      IF (PRESENT(ex_ker_offdiag)) do_offdiag = .TRUE.
      ndo_mo = SIZE(donor_state%mo_indices,1)
      CALL get_qs_env(qs_env, natom=nblk, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                      blacs_env=context, para_env=para_env)       
      ALLOCATE(first_sgf(nblk))        
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf)

!     Get the (P|Q)^-1 from the xas_tdp_env 
      i = locate(xas_tdp_env%ex_kind_indices, donor_state%kind_index)
      PQ => xas_tdp_env%ri_inv_mats(i)%array
      nsgfp = SIZE(PQ,1)
      nao = SIZE(xas_tdp_env%ri_coul3, 1)

!     Get the projector Q as an array (dbcsr => fm => array)
      ALLOCATE(proj_Q(nao,nao))
      CALL cp_fm_struct_create(fm_struct, context=context, para_env=para_env, nrow_global=nao, &
                               ncol_global=nao)
      CALL cp_fm_create(work_fm, fm_struct) 
      CALL copy_dbcsr_to_fm(xas_tdp_env%q_projector, work_fm)
      CALL cp_fm_get_submatrix(fm=work_fm, target_m=proj_Q) 
      CALL cp_fm_struct_release(fm_struct) 
      CALL cp_fm_release(work_fm)    

!     Build the dbcsr matrix for the kernels
      CALL dbcsr_create(matrix=ex_ker_diag, name="EX KER DIAG", matrix_type='S', dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size) 
      IF (do_offdiag) THEN
         CALL dbcsr_create(matrix=ex_ker_offdiag, name="EX KER OFFDIAG", matrix_type='S',&  
                           dist=dist, row_blk_size=blk_size, col_blk_size=blk_size)
      END IF

!     Allocate the many work arrays
      IF (do_offdiag) THEN
         ALLOCATE(aJIb_prod(nao,nao))
         ALLOCATE(work1(nao, nsgfp))  
         ALLOCATE(work2(nao,nao))
      END IF
      ALLOCATE(work3(nsgfp,1))
      ALLOCATE(work4(nao,nao))
      ALLOCATE(abIJ_prod(nao,nao))

!  Loop over the donor MOs and compute both the digonal and off-digonal parts, that is:
!     -Diag_ab = (ab|P)* (P|Q)^-1 *(Q|IJ), then multiply by proj_Q from the left and proj_Q^T from the right
!     -Off_ab = (aJ|P)* (P|Q)^-1 *(Q|Ib), then multiply by proj_Q from the left and proj_Q^T from the right
!  The matrices are symmetric => only loop on upper triangular part
      DO i = 1,ndo_mo
         IF (do_offdiag) CALL contract_int_kernel(int_bI, (/donor_state%mo_indices(i)/), & 
                                                  xas_tdp_env, qs_env)
         DO j = i,ndo_mo
            IF (do_offdiag) CALL contract_int_kernel(int_aJ, (/donor_state%mo_indices(j)/), &
                                                     xas_tdp_env, qs_env)
            CALL contract_int_kernel(int_IJ, (/donor_state%mo_indices(i), donor_state%mo_indices(j)/), &
                                     xas_tdp_env, qs_env)

            IF (do_offdiag) THEN
               work1 = 0.0_dp
               work2 = 0.0_dp
               aJIb_prod = 0.0_dp
            END IF
            work3 = 0.0_dp
            work4 = 0.0_dp
            abIJ_prod = 0.0_dp

!           Start with off-diag product (aJ|P)* (P|Q)^-1 *(Q|Ib)
            IF (do_offdiag) THEN

!              (aJ|P)* (P|Q)^-1
               CALL dgemm('N', 'N', nao, nsgfp, nsgfp, 1.0_dp, int_aJ, nao, PQ, nsgfp, 0.0_dp, &
                           work1, nao)
!              Muliply by (Q|Ib)  
               CALL dgemm('N', 'T', nao, nao, nsgfp, 1.0_dp, work1, nao, int_bI, nao, 0.0_dp, &
                          aJIb_prod, nao)
!              Multiply that by the Q projector
               CALL dgemm('N', 'N', nao, nao, nao, 1.0_dp, proj_Q, nao, aJIb_prod, nao, 0.0_dp, &
                          work2, nao)
               aJIb_prod = 0.0_dp
               CALL dgemm('N', 'T', nao, nao, nao, 1.0_dp, work2, nao, proj_Q, nao, 0.0_dp, &
                          aJIb_prod, nao)  
            END IF

!           Then the on-digonal product (ab|P)* (P|Q)^-1 *(Q|IJ) 

!           (P|Q)^-1 *(Q|IJ)
            CALL dgemm('N', 'T', nsgfp, 1, nsgfp, 1.0_dp, PQ, nsgfp, int_IJ, 1, 0.0_dp, &
                       work3, nsgfp)

!           Multiply by (ab|P) from the left 
            DO k = 1,nsgfp
               CALL daxpy(nao*nao, work3(k,1), xas_tdp_env%ri_coul3(:,:,k), 1, abIJ_prod(:,:), 1)
            END DO !k

!           Multiply by the Q projector
            CALL dgemm('N', 'N', nao, nao, nao, 1.0_dp, proj_Q, nao, abIJ_prod, nao, 0.0_dp, &
                       work4, nao)
            abIJ_prod = 0.0_dp
            CALL dgemm('N', 'T', nao, nao, nao, 1.0_dp, work4, nao, proj_Q, nao, 0.0_dp, &
                       abIJ_prod, nao)

!           Put it all in the dbcsr matrices, at the right I,J and block by block
            DO iblk = 1,nblk
               DO jblk = 1,nblk

!                 If on the block digonal, only take jblk >= iblk
                  IF (jblk < iblk .AND. i==j) CYCLE

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
            IF (do_offdiag) DEALLOCATE(int_aJ)
            DEALLOCATE(int_IJ) 
         END DO !j
         IF (do_offdiag) DEALLOCATE(int_bI)
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
!> \note currently coded with arrays and such low level types. TODO: switch to dbcsr tensors         
!>       Everythin is gonna change later on, so don't care about performance
! **************************************************************************************************
   SUBROUTINE build_kernel_xc(xc_kernel, dist, blk_size, donor_state, xas_tdp_env, qs_env)

      TYPE(dbcsr_type), POINTER                       :: xc_kernel
      TYPE(dbcsr_distribution_type), POINTER          :: dist
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "build_kernel_xc", &                            
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, ndo_mo, nblk, i, nsgfp, nao, j, &
                                                         iblk, jblk, nrow_blk, ncol_blk
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set 
      INTEGER, DIMENSION(:), ALLOCATABLE              :: first_sgf
      TYPE(cp_blacs_env_type), POINTER                :: context
      TYPE(cp_para_env_type), POINTER                 :: para_env
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: proj_Q, int_aI, int_bJ, work1, &
                                                         work2, aIJb_prod, work3
      REAL(dp), DIMENSION(:,:), POINTER               :: PQ
      REAL(dp), DIMENSION(:,:), POINTER               :: fxc
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct
      TYPE(cp_fm_type), POINTER                       :: work_fm
      LOGICAL                                         :: found_block

      CALL timeset(routineN, handle)

      NULLIFY(particle_set, qs_kind_set, context, para_env, fm_struct, work_fm, PQ)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices,1) 
      CALL get_qs_env(qs_env, natom=nblk, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                      blacs_env=context, para_env=para_env)
      ALLOCATE(first_sgf(nblk))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf)

!     Get the (P|Q)^-1 as an array
      i = locate(xas_tdp_env%ex_kind_indices, donor_state%kind_index)
      PQ => xas_tdp_env%ri_inv_mats(i)%array
      nsgfp = SIZE(PQ, 1)
      nao = SIZE(xas_tdp_env%ri_coul3, 1)

!     Get the projector Q as an array (dbcsr => fm => array)
      ALLOCATE(proj_Q(nao,nao))
      CALL cp_fm_struct_create(fm_struct, context=context, para_env=para_env, nrow_global=nao, &     
                                        ncol_global=nao)
      CALL cp_fm_create(work_fm, fm_struct)
      CALL copy_dbcsr_to_fm(xas_tdp_env%q_projector, work_fm)
      CALL cp_fm_get_submatrix(fm=work_fm, target_m=proj_Q)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_release(work_fm)

!     Get the 2-center integrel (Q|fxc|R) as an array (a block from a dbcsr matrix)
!     For now assume serial => don't check for which proc the block is on
!     TODO: maybe that storing the atom integrals in a dbcsr matrix is not the best
      found_block = .FALSE.
      CALL dbcsr_get_block_p(xas_tdp_env%ri_atom_int%matrix, donor_state%at_index, &
                            donor_state%at_index, fxc, found_block)
      CPASSERT(found_block)
      CPASSERT(ALL(SHAPE(fxc)==nsgfp))

!     Build the dbcsr matrix for the kernel
      CALL dbcsr_create(matrix=xc_kernel, name="XC KERNEL", matrix_type='S', dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)

!     Allocating work arrays to store the intemediate products
      ALLOCATE(work1(nao, nsgfp))
      ALLOCATE(work2(nao, nsgfp))
      ALLOCATE(work3(nao, nao))
      ALLOCATE(aIJb_prod(nao, nao))

!  Loop over the donor MOs and compute the matrix product projQ * W *projQ^T, where W is
!  W_ab = (aI|P) * (P|Q)^-1 * (Q|fxc|R) * (R|S)^-1 * (S|Jb) ~= (aI|fxc|Jb)
!  Symmetric matrix => only span upper triangle
      DO i = 1,ndo_mo
         CALL contract_int_kernel(int_aI, (/donor_state%mo_indices(i)/), xas_tdp_env, qs_env)
         DO j = i,ndo_mo
            CALL contract_int_kernel(int_bJ, (/donor_state%mo_indices(j)/), xas_tdp_env, qs_env)

            work1 = 0.0_dp
            work2 = 0.0_dp
            work3 = 0.0_dp
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

!           Muliply by the projector from both sides => projQ * W *projQ^T
            CALL dgemm('N', 'N', nao, nao, nao,  1.0_dp, proj_Q, nao, aIJb_prod, nao, 0.0_dp, &      
                        work3, nao)
            aIJb_prod = 0.0_dp
            CALL dgemm('N', 'T', nao, nao, nao,  1.0_dp, work3, nao, proj_Q, nao, 0.0_dp, &          
                        aIJb_prod, nao)

!           Put the product in the dbcsr_matrix, for the right I,J block
            DO iblk = 1,nblk
               DO jblk = 1,nblk

!                 If on the block diagonal, just take jblk >= iblk
                  IF (jblk < iblk .AND. i==j) CYCLE

!                 Size of the current block
                  nrow_blk = blk_size(iblk)
                  ncol_blk = blk_size(jblk)

                  CALL dbcsr_put_block(xc_kernel, (i-1)*nblk+iblk, (j-1)*nblk+jblk, & 
                                       aIJb_prod(first_sgf(iblk):first_sgf(iblk)+nrow_blk-1, &
                                                 first_sgf(jblk):first_sgf(jblk)+ncol_blk-1))

               END DO ! jblk
            END DO ! iblk

            DEALLOCATE(int_bJ)
         END DO ! donor MO j
         DEALLOCATE(int_aI)
      END DO ! donor MO i

      CALL dbcsr_finalize(xc_kernel)

      CALL timestop(handle)

   END SUBROUTINE build_kernel_xc

END MODULE xas_tdp_utils

