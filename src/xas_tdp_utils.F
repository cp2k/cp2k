!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for X-ray absorption spectroscopy using TDDFPT
!> \author AB (01.2018)
! **************************************************************************************************

MODULE xas_tdp_utils
   
   USE basis_set_types,                 ONLY: gto_basis_set_type, get_gto_basis_set, &
                                              gto_basis_set_p_type
   USE cell_types,                      ONLY: cell_type, pbc
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm, cp_dbcsr_sm_fm_multiply, &
                                              copy_fm_to_dbcsr
   USE cp_eri_mme_interface,            ONLY: cp_eri_mme_finalize, cp_eri_mme_param, &
                                              cp_eri_mme_set_params, cp_eri_mme_init_read_input, &
                                              cp_eri_mme_update_local_counts
   USE cp_cfm_types,                    ONLY: cp_fm_to_cfm, cp_cfm_to_fm, cp_cfm_type, &
                                              cp_cfm_create, cp_cfm_release
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add, cp_fm_gemm, cp_fm_column_scale, &
                                              cp_fm_invert
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
                                              dbcsr_copy, dbcsr_set, dbcsr_deallocate_matrix_set, &
                                              dbcsr_allocate_matrix_set, dbcsr_multiply, &
                                              dbcsr_get_stored_coordinates, dbcsr_desymmetrize, &
                                              dbcsr_complete_redistribute, dbcsr_iterator_type, &
                                              dbcsr_iterator_blocks_left, dbcsr_iterator_start, &
                                              dbcsr_iterator_stop, dbcsr_iterator_next_block, &
                                              dbcsr_add
   USE eri_mme_types,                   ONLY: eri_mme_init
   USE input_constants,                 ONLY: xas_dip_len, xas_dip_vel, do_eri_mme, do_eri_os
   USE iterate_matrix,                  ONLY: invert_Hotelling
   USE kinds,                           ONLY: dp
   USE mp2_eri,                         ONLY: integrate_set_3c, integrate_set_2c
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
   USE xas_tdp_types,                   ONLY: donor_state_type, xas_tdp_control_type, &
                                              xas_tdp_env_type

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_utils'

   PUBLIC :: setup_xas_tdp_full_prob, solve_xas_tdp_full_prob, compute_oscillator_strength, &
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
!>       - For singlet: M_d = A + 2B - D
!>                      M_o = 2B - E
!>       - For triplet: M_d = A - D
!>                      M_o - E
!>       where other subroutines computes the matrices A, B, E, D and G, which are:
!>       - A: the ground-state contribution F_ab*delta_IJ - epsilon_IJ*S_ab
!>       - B: the Coulob kernel ~(aI|Jb)
!>       - D: the on-digonal exact exchange kernel ~(ab|IJ)
!>       - E: the off-diagonal exact exchange kernel ~(aJ|Ib)
!>       - G: the metric  S_ab*delta_IJ
!>       In the above definitions, I,J label the donnor MOs and a,b the sgfs of the basis
!>       The matrix M (matrix_tdp) has the same (i.e. 2x2 replicated) dbcsr_dist as usual
! **************************************************************************************************
   SUBROUTINE setup_xas_tdp_full_prob(donor_state, qs_env, xas_tdp_env, xas_tdp_control)

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      
      CHARACTER(len=*), PARAMETER :: routineN = 'setup_xas_tdp_full_prob', &
                                     routineP = moduleN//":"//routineN
      
      INTEGER                                         :: handle, group , iblk, jblk, nblk_sub, blk
      
      TYPE(dbcsr_type), POINTER                       :: matrix_tdp
      TYPE(dbcsr_type), POINTER                       :: matrix_a, matrix_g, diag_mat, & 
                                                         matrix_g_inv, offdiag_mat, matrix_b, &
                                                         matrix_d, matrix_e
      TYPE(dbcsr_distribution_type), POINTER          :: submat_dist, supmat_dist
      INTEGER, DIMENSION(:), POINTER                  :: submat_blk_size, supmat_blk_size, &
                                                         col_dist_sub, row_dist_sub, &
                                                         col_dist_sup, row_dist_sup
      INTEGER, DIMENSION(:,:), POINTER                :: pgrid
      LOGICAL                                         :: found_block
      REAL(dp), DIMENSION(:), POINTER                 :: work_block
      TYPE(dbcsr_iterator_type)                       :: iter     

      CALL timeset(routineN, handle)

      NULLIFY(matrix_a, submat_dist, submat_blk_size, supmat_dist, supmat_blk_size, col_dist_sub)
      NULLIFY(row_dist_sub, col_dist_sup, row_dist_sup, work_block, matrix_g, pgrid, diag_mat)
      NULLIFY(matrix_tdp, offdiag_mat, matrix_b, matrix_d, matrix_e) 

      CPASSERT(ASSOCIATED(donor_state))
      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(xas_tdp_env))
      CPASSERT(ASSOCIATED(xas_tdp_control))
      ALLOCATE(donor_state%matrix_tdp)
      matrix_tdp => donor_state%matrix_tdp

!  Get the submatrices distribution and block sizes
      CALL compute_submat_dist_and_blk_size(submat_dist, submat_blk_size, donor_state, qs_env)
      nblk_sub = SIZE(submat_blk_size)

!  Compute the submatrices A, G and the different kernels (B, D and E)
      ALLOCATE(matrix_a)
      ALLOCATE(matrix_b)
      ALLOCATE(matrix_d)
      ALLOCATE(matrix_e)
      ALLOCATE(matrix_g)
      CALL build_gs_contribution(matrix_a, submat_dist, submat_blk_size, donor_state, qs_env)
      CALL build_kernel_coulomb(matrix_b, submat_dist, submat_blk_size, donor_state, xas_tdp_env, &
                                qs_env)
      CALL build_kernel_exchange(matrix_d, matrix_e, submat_dist, submat_blk_size, &
                                 donor_state, xas_tdp_env, qs_env)
      CALL build_metric(matrix_g, submat_dist, submat_blk_size, donor_state, qs_env)       

!  Invert matrix G
      ALLOCATE(matrix_g_inv)
      CALL dbcsr_create(matrix=matrix_g_inv, name="MATRIX G INV", matrix_type="S", &              
                        dist=submat_dist, row_blk_size=submat_blk_size, col_blk_size=submat_blk_size)
!     TODO: impose values of threshold and eps_filter or not?. Current values taken from LRI
      CALL invert_Hotelling(matrix_g_inv, matrix_g, threshold=1.0E-10_dp, use_inv_as_guess=.FALSE., &
                            norm_convergence=1.0E-10_dp, filter_eps=1.0E-12_dp, silent=.TRUE.)
      CALL dbcsr_finalize(matrix_g_inv)

!  Create the matrix that will be on the diagonal of matrix_tdp: G^-1*(A+B-D)
      ALLOCATE(diag_mat)
      CALL dbcsr_copy(diag_mat, matrix_a, name="Diag TDP block")
      IF (xas_tdp_control%do_singlet) CALL dbcsr_add(diag_mat, matrix_b, 1.0_dp, 2.0_dp)
      CALL dbcsr_add(diag_mat, matrix_d, 1.0_dp, -1.0_dp)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_g_inv, diag_mat, 0.0_dp, diag_mat)
      CALL dbcsr_finalize(diag_mat)

!  Create the matrix that will be in the off-diagonal of matrix_tdp: G^-1*(B-E)
      ALLOCATE(offdiag_mat)
      CALL dbcsr_copy(offdiag_mat, matrix_b, name="Offdiag TDP block")
      IF (xas_tdp_control%do_singlet) THEN
         CALL dbcsr_add(offdiag_mat, matrix_e, 2.0_dp, -1.0_dp)
      ELSE
         CALL dbcsr_add(offdiag_mat, matrix_e, 0.0_dp, -1.0_dp)
      END IF
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_g_inv, offdiag_mat, 0.0_dp, offdiag_mat)
      CALL dbcsr_finalize(offdiag_mat)

!  Some clean-up
      CALL dbcsr_release(matrix_g_inv)
      CALL dbcsr_release(matrix_a)
      CALL dbcsr_release(matrix_b)
      CALL dbcsr_release(matrix_d)
      CALL dbcsr_release(matrix_e)
      CALL dbcsr_release(matrix_g)
      DEALLOCATE(matrix_g_inv, matrix_a, matrix_b, matrix_g, matrix_d, matrix_e)

!  Test on Tamm Dancoff approximation 

      IF (xas_tdp_control%tamm_dancoff) THEN

!        If TDA is assumed, then the matrix_tdp is G^-1*(A+2B-D-E)
         CALL dbcsr_copy(matrix_tdp, diag_mat, name="MATRIX TDP")
         CALL dbcsr_add(matrix_tdp, offdiag_mat, 1.0_dp, 1.0_dp)

      ELSE IF (.NOT. xas_tdp_control%tamm_dancoff) THEN

!     The supermatrix matrix_tdp is built with the submatrices as a 2x2 block matrices 
!     Create the dbcsr dist by repeating the submatrices structure
         ALLOCATE (supmat_dist)

         CALL dbcsr_distribution_get(submat_dist, row_dist=row_dist_sub, col_dist=col_dist_sub, &
                                     group=group, pgrid=pgrid)

         ALLOCATE (row_dist_sup(2*SIZE(row_dist_sub)))
         ALLOCATE (col_dist_sup(2*SIZE(col_dist_sub)))

         row_dist_sup(1:SIZE(row_dist_sub)) = row_dist_sub
         col_dist_sup(1:SIZE(col_dist_sub)) = col_dist_sub
         row_dist_sup(1+SIZE(row_dist_sub):2*SIZE(row_dist_sub)) = row_dist_sub
         col_dist_sup(1+SIZE(col_dist_sub):2*SIZE(col_dist_sub)) = col_dist_sub

         CALL dbcsr_distribution_new(supmat_dist, group=group, pgrid=pgrid, &
                                     col_dist=col_dist_sup, row_dist=row_dist_sup)

         DEALLOCATE (row_dist_sup, col_dist_sup)

!        Creating the non-symmetric supermatrix matrix_tdp  with the 2x2 replication of the submatrices
!        block sizes and distribution.
         ALLOCATE(supmat_blk_size(2*nblk_sub))
         supmat_blk_size(1:nblk_sub) = submat_blk_size
         supmat_blk_size(nblk_sub+1:2*nblk_sub) = submat_blk_size 
         CALL dbcsr_create(matrix=matrix_tdp, name="XAS TDP MATRIX",matrix_type="N", dist=supmat_dist,&
                           row_blk_size=supmat_blk_size, col_blk_size=supmat_blk_size)
         DEALLOCATE(supmat_blk_size)

!     Filling the supermatrix matrix_tdp with submatrices +-G^-1*(A+B-D) and +-G^-1*(B-E)

!        For each diagonal block of the super matrix, loop over blocks of G^-1*(A+B-D)
         CALL dbcsr_iterator_start(iter, diag_mat)
         DO WHILE(dbcsr_iterator_blocks_left(iter))

            found_block = .FALSE.
            CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
            CALL dbcsr_get_block_p(diag_mat, iblk, jblk, work_block, found_block)

            IF (found_block) THEN
!              Put -G^-1*(A+B-D) in the upper left block and G^-1*(A+B-D) in the lower right block
               CALL dbcsr_put_block(matrix_tdp, iblk, jblk, -work_block) 
               CALL dbcsr_put_block(matrix_tdp, nblk_sub+iblk, nblk_sub+jblk, work_block)
            END IF
            NULLIFY(work_block)

         END DO ! iterator
         CALL dbcsr_iterator_stop(iter) 

!        For each off-diagonal block of supermatrix, loop over blocks of G^-1*(B-E)
         CALL dbcsr_iterator_start(iter, offdiag_mat)
         DO WHILE(dbcsr_iterator_blocks_left(iter))

            found_block = .FALSE.
            CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
            CALL dbcsr_get_block_p(offdiag_mat, iblk, jblk, work_block, found_block) 

            IF (found_block) THEN 
!              Put -G^-1*B in the upper right and G^-1*B in the lower left
               CALL dbcsr_put_block(matrix_tdp, iblk, nblk_sub+jblk, -work_block)
               CALL dbcsr_put_block(matrix_tdp, nblk_sub+iblk, jblk, work_block)
            END IF
            NULLIFY(work_block)

         END DO ! iterator
         CALL dbcsr_iterator_stop(iter)

!        Some clean-up
         CALL dbcsr_distribution_release(supmat_dist)
         DEALLOCATE(supmat_dist)

      END IF ! Tamm_Dancoff

!  Finalize
      CALL dbcsr_finalize(matrix_tdp)

!  General Clean-up
      CALL dbcsr_distribution_release(submat_dist)
      CALL dbcsr_release(diag_mat)
      CALL dbcsr_release(offdiag_mat)
      DEALLOCATE(submat_dist, submat_blk_size, diag_mat, offdiag_mat)

      CALL timestop(handle)

   END SUBROUTINE setup_xas_tdp_full_prob

! **************************************************************************************************
!> \brief Solves the XAS TDP generalized eigenvalue problem omega*C = matrix_tdp*C using standard
!>        full diagonalization methods (for non-hermitian matrices), for a given donor state
!> \param donor_state ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note The computed eigenvalues and eigenvectors are stored in the donor_state
! **************************************************************************************************
   SUBROUTINE solve_xas_tdp_full_prob(donor_state, xas_tdp_control, qs_env)  

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "solve_xas_tdp_full_prob", &
                                     routineP = moduleN//":"//routineN

      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(dbcsr_type), POINTER                       :: matrix_tdp
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      INTEGER                                         :: handle
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct
      TYPE(cp_fm_type), POINTER                       :: amatrix
      INTEGER                                         :: nrow

      CALL timeset(routineN, handle)

      NULLIFY(para_env, blacs_env, fm_struct, amatrix, matrix_tdp)

      matrix_tdp => donor_state%matrix_tdp
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env) 

!     Initialization
      CALL dbcsr_get_info(matrix_tdp, nfullrows_total=nrow)
      CALL cp_fm_struct_create(fm_struct,context=blacs_env, nrow_global=nrow, para_env=para_env,&
                               ncol_global=nrow)

!     Need to go to full matrices to diagonalize
      CALL cp_fm_create(amatrix, fm_struct)
      CALL copy_dbcsr_to_fm(matrix_tdp, amatrix)
         
!     Allocating space for solutions
      ALLOCATE(donor_state%lr_evals(nrow))
      CALL cp_fm_create(donor_state%lr_coeffs, fm_struct)

!     Diagonalize the matrix
#if defined(__SCALAPACK) 
      CALL pdiag_nonsymm_mat(amatrix, donor_state%lr_evals, donor_state%lr_coeffs)
#else
      CALL diag_nonsymm_mat(amatrix, donor_state%lr_evals, donor_state%lr_coeffs)
#endif
         
!     Clean-up
      CALL cp_fm_release(amatrix)
      CALL cp_fm_struct_release(fm_struct)

!  Making sure that the response orbitals are normalized with respect to matrix_s
!  TODO: RETEST WHEN KERNEL IS AVAILABLE
      CALL normalize_lr_orbitals(donor_state, xas_tdp_control, qs_env)

      CALL timestop(handle)

   END SUBROUTINE solve_xas_tdp_full_prob

! **************************************************************************************************
!> \brief Normalizes the response orbitals so that |C^-T * G * C^- - C^+T * G *C^+| = +-1 
!> \param donor_state ...
!> \param xas_tdp_control ...
!> \param qs_env ... 
! **************************************************************************************************
   SUBROUTINE normalize_lr_orbitals(donor_state, xas_tdp_control, qs_env)

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "normalize_lr_orbitals", &
                                      routineP = moduleN//":"//routineN     

      INTEGER                                         :: handle, nao, nvals, ndo_mo, ido_mo 
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(cp_fm_struct_type), POINTER                :: c_struct, norm_struct
      TYPE(cp_fm_type), POINTER                       :: c_plus, c_minus, work_vect, lr_coeffs, &
                                                         fm_norm
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: scaling, norm
      
      CALL timeset(routineN, handle) 

      NULLIFY(para_env, blacs_env, matrix_s)
      NULLIFY(c_struct, c_plus, c_minus, work_vect, lr_coeffs, fm_norm, norm_struct)
      
      CPASSERT(ASSOCIATED(donor_state%lr_coeffs))

!  Initialization
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      nvals = SIZE(donor_state%lr_evals)

      CALL cp_fm_struct_create(c_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=nao, ncol_global= nvals)
      CALL cp_fm_struct_create(norm_struct, para_env=para_env, context=blacs_env, & 
                               nrow_global=nvals, ncol_global=nvals)   

      CALL cp_fm_create(c_plus, c_struct)
      CALL cp_fm_create(c_minus, c_struct)
      CALL cp_fm_create(work_vect, c_struct)
      CALL cp_fm_create(fm_norm, norm_struct)

      ndo_mo = SIZE(donor_state%mo_indices)
      ALLOCATE(scaling(nvals))
      scaling = 0.0_dp
      ALLOCATE(norm(nvals))
      norm = 0.0_dp

      lr_coeffs => donor_state%lr_coeffs

!  Loop over the different donor MOs, consider all eigenvectors at once
      DO ido_mo = 1,ndo_mo

!        If TDA, only have the C^-
         IF (xas_tdp_control%tamm_dancoff) THEN

!           Retrieve the C^- coeff for current donor_state
            CALL cp_fm_to_fm_submat(msource=lr_coeffs, mtarget=c_minus, nrow=nao, ncol=nvals, &
                                    s_firstrow=(ido_mo-1)*nao+1, s_firstcol=1, t_firstrow=1, &
                                    t_firstcol=1)

         ELSE ! not TDA

!           Retrieve the C^+ coeff for current donor_state
            CALL cp_fm_to_fm_submat(msource=lr_coeffs, mtarget=c_plus, nrow=nao, ncol=nvals, &      
                                    s_firstrow=(ido_mo-1)*nao+1, s_firstcol=1, t_firstrow=1, &
                                    t_firstcol=1)  

!           Retrieve the C^- coeff for current donor_state
            CALL cp_fm_to_fm_submat(msource=lr_coeffs, mtarget=c_minus, nrow=nao, ncol=nvals, &
                                    s_firstrow=(ndo_mo+ido_mo-1)*nao+1, s_firstcol=1, t_firstrow=1, &
                                    t_firstcol=1)

         END IF ! TDA

!        Compute C^-T * S * C^- for current donor state
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, c_minus, work_vect, ncol=nvals)
         CALL cp_fm_gemm('T', 'N', nvals, nvals, nao, 1.0_dp, c_minus, work_vect, 0.0_dp, fm_norm)

!        Add the diagonal (corrsponding to the product for each eigenvector) to the scaling array
         CALL cp_fm_get_diag(fm_norm, norm) 
         scaling(:) = scaling(:) + norm(:)

!        If not TDA, need to substract the C^+T * S * C^+ product
         IF (.NOT. xas_tdp_control%tamm_dancoff) THEN
   
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, c_plus, work_vect, ncol=nvals)
            CALL cp_fm_gemm('T', 'N', nvals, nvals, nao, 1.0_dp, c_plus, work_vect, 0.0_dp, fm_norm)
            CALL cp_fm_get_diag(fm_norm, norm)
            scaling(:) = scaling(:) - norm(:)

         END IF ! not TDA

      END DO ! ido_mo

!  Need to scale each eigenvector 
      scaling = 1.0_dp/SQRT(ABS(scaling))
      CALL cp_fm_column_scale(lr_coeffs, scaling)

!  Clean-up
      CALL cp_fm_struct_release(c_struct)
      CALL cp_fm_struct_release(norm_struct)
      CALL cp_fm_release(c_plus)
      CALL cp_fm_release(c_minus)
      CALL cp_fm_release(work_vect)
      CALL cp_fm_release(fm_norm)
      CALL timestop(handle)
      
   END SUBROUTINE normalize_lr_orbitals

! **************************************************************************************************
!> \brief Computes the oscillator strength based on the dipole moment (velocity or position rep
!>        for all available excitation energies of a donor_state. Stores the result in donor_state
!> \param donor_state the donor state which is excited
!> \param xas_tdp_control ...
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note The oscillator strength is a scalar: osc_str = 2/3*omega*(dipole_r)^2 = 2/(3*omega)*(dipole_v)^2
!>       But because there is a risk of very tiny omega, multiply both cases by omega to avoid division
!>       The formulae for the dipoles come from the trace of the dipole operator with the transition
!>       densities, i.e. what we get from digonalizing the matrix_tdp
! **************************************************************************************************
   SUBROUTINE compute_oscillator_strength(donor_state, xas_tdp_control, xas_tdp_env, qs_env)

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      
      CHARACTER(len=*), PARAMETER :: routineN = "compute_oscillator_strength", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, nosc, nao, i, ndo_mo,j, k, jdo_mo,l
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      TYPE(cp_fm_type), POINTER                       :: mo_coeff
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: dipmat
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      REAL(dp), DIMENSION(3)                          :: xyz_moments
      REAL(dp), DIMENSION(:), POINTER                 :: osc_str
      REAL(dp), DIMENSION(:,:), POINTER               :: gs_coeff, work_array
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(cp_fm_struct_type), POINTER                :: col_vect_struct
      TYPE(cp_fm_type), POINTER                       :: c_vect, work_vect
      REAL(dp)                                        :: contrib

      CALL timeset(routineN, handle)

      NULLIFY(mos, mo_coeff, para_env, blacs_env, dipmat, particle_set, osc_str, matrix_s)
      NULLIFY(col_vect_struct, c_vect, work_vect, gs_coeff, work_array)

!     General initialization
      ndo_mo = SIZE(donor_state%mo_indices)
      nosc = SIZE(donor_state%lr_evals) 
      ALLOCATE(donor_state%osc_str(nosc))
      osc_str => donor_state%osc_str
      osc_str = 0.0_dp
      CALL get_qs_env(qs_env, mos=mos, para_env=para_env, blacs_env=blacs_env, &
                      particle_set=particle_set, matrix_s=matrix_s)
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao)
      dipmat => xas_tdp_env%dipmat

!     TODO: The normalization of the response coeff is kind of abitrary (but consistent) 
!     => factor 2/3 can probably be dropped
      IF (xas_tdp_control%dipole_form == xas_dip_len) THEN
!        Add the dipole form specific prefactor 2/3*omega
         osc_str = 2.0_dp/3.0_dp*donor_state%lr_evals
         
      ELSE IF (xas_tdp_control%dipole_form == xas_dip_vel) THEN
!        Add the dipole form specific prefactor 2/3/omega, no minus sign as p_xyz is real here
         osc_str = 2.0_dp/3.0_dp/donor_state%lr_evals 

      END IF

!     Some more initialization before heavy computing
      CALL cp_fm_struct_create(col_vect_struct, para_env=para_env, context=blacs_env, &
                                nrow_global=nao, ncol_global=1)
      CALL cp_fm_create(c_vect, col_vect_struct)
      CALL cp_fm_create(work_vect, col_vect_struct)
      ALLOCATE(gs_coeff(nao,1))
      ALLOCATE(work_array(nao,1))

!     the dipole moment is <mu> = sum_j C^0_j*M*C^1_j, where j labels donor MOs, M is the dipmat,
!     C^0_i is the ground state coefficient row vector and and C^1 the column vector of lr coeffs
      DO i = 1,nosc ! loop on excitation energies
          xyz_moments = 0.0_dp
          DO j = 1,ndo_mo ! loop on donor MOs

            jdo_mo = donor_state%mo_indices(j)

!           Retrieve column vector C^1. Test for Tamm Dancoff
            IF (xas_tdp_control%tamm_dancoff) THEN

!              For TDA, C^+ = 0, hence C^1 = C^- and lr_coeffs = C^-
               CALL cp_fm_to_fm_submat(msource=donor_state%lr_coeffs, mtarget=c_vect, nrow=nao, & 
                                       ncol=1,s_firstrow=(j-1)*nao+1, s_firstcol=i, t_firstrow=1, &
                                       t_firstcol=1)

            ELSE IF (.NOT. xas_tdp_control%tamm_dancoff) THEN

!              In the general case, C^1 = C^+ + C^-, which are both stored in lr_coeffs
               CALL cp_fm_to_fm_submat(msource=donor_state%lr_coeffs, mtarget=c_vect, nrow=nao, &
                                       ncol=1, s_firstrow=(j-1)*nao+1, s_firstcol=i, t_firstrow=1, &
                                       t_firstcol=1) ! C^+

               CALL cp_fm_to_fm_submat(msource=donor_state%lr_coeffs, mtarget=work_vect, nrow=nao, &
                                       ncol=1, s_firstrow=(ndo_mo+j-1)*nao+1, s_firstcol=i, & 
                                       t_firstrow=1, t_firstcol=1) ! C^-

!              store the sum in c_vect
               CALL cp_fm_scale_and_add(1.0_dp, c_vect, 1.0_dp, work_vect)

            END IF

            gs_coeff = 0.0_dp
!           Retrieve array of ground state MO coeff
            CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=gs_coeff, start_row=1, & 
                                     start_col=jdo_mo, n_rows=nao, n_cols=1, transpose=.FALSE.)

            DO k =1,3 ! loop on cartesian coordinates
!              Computing M*C^-_j
               CALL cp_dbcsr_sm_fm_multiply(dipmat(k)%matrix, c_vect, work_vect, ncol=1)
               work_array = 0.0_dp
               CALL cp_fm_get_submatrix(work_vect, work_array)

!              Computing C^0_j*M*C^-_j
               DO l = 1,nao 
                  xyz_moments(k) = xyz_moments(k) + gs_coeff(l,1)*work_array(l,1) 
               END DO ! l
            END DO ! k
         END DO ! j

!        Compute the oscilaltor strength for this specific energy
         contrib = 0.0_dp
         DO k = 1,3
            contrib = contrib +  osc_str(i)*xyz_moments(k)**2
         END DO ! k
         osc_str(i) = contrib
      END DO ! i
         
!     Clean-up
      DEALLOCATE(gs_coeff, work_array)
      CALL cp_fm_release(c_vect)
      CALL cp_fm_release(work_vect)
      CALL cp_fm_struct_release(col_vect_struct)

      CALL timestop(handle)

   END SUBROUTINE compute_oscillator_strength

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
!>         => A_{pi,qj} = F_pq*delta_ij - epsilon_ij*S_pq
!> \param matrix_a  pointer to a DBCSR matrix containing A
!> \param dist_a the dbcsr distribution of matrix A
!> \param blk_size_a the row/col block sizes of matrix A
!> \param donor_state ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE build_gs_contribution(matrix_a, dist_a, blk_size_a, donor_state, qs_env)

      TYPE(dbcsr_type), POINTER                       :: matrix_a
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(dbcsr_distribution_type), POINTER          :: dist_a
      INTEGER, DIMENSION(:), POINTER                  :: blk_size_a
      
      CHARACTER(len=*), PARAMETER :: routineN = 'build_gs_contribution', &
                                     routineP = moduleN//":"//routineN
      
      INTEGER                                         :: handle, ndo_mo, i, nblk_row,&
                                                         j, iblk, jblk, blk
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks, matrix_s
      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      INTEGER, DIMENSION(:), POINTER                  :: row_blk_size 
      REAL(dp), DIMENSION(:,:), POINTER               :: energy_evals     
      REAL(dp), DIMENSION(:), POINTER                 :: work_block
      LOGICAL                                         :: found_block
      TYPE(dbcsr_type)                                :: desym_matrix_s, desym_matrix_ks, work_matrix
      TYPE(dbcsr_iterator_type)                       :: iter

      CALL timeset(routineN, handle)

      NULLIFY(matrix_ks, dbcsr_dist, row_blk_size, work_block, energy_evals, matrix_s)

!  Note: build matrix A as a non-symmetric matrix. Ok since need to desymmetrixe matrix S and KS and
!        will be later added to kernel elements, which are non-symmetric themselves

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

!  Put the current work matrix into the matrix A
      CALL dbcsr_copy(matrix_a, work_matrix, name="MATRIX A")

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
            DO i = 1,ndo_mo
               DO j = 1,ndo_mo

!                 Put epsilon_ij*S_pq
                  CALL dbcsr_put_block(work_matrix,(i-1)*nblk_row+iblk, (j-1)*nblk_row+jblk, & 
                                       energy_evals(i,j)*work_block)   

               END DO ! j
            END DO ! i

         END IF ! found_block
         NULLIFY(work_block)

      END DO ! iteration over S
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(work_matrix)

!  Substract the current work matrix to matrix A to get the full F_pq*delta_ij - epsilon_ij*S_pq
      CALL dbcsr_add(matrix_a, work_matrix, 1.0_dp, -1.0_dp)
      
!  Finalize
      CALL dbcsr_finalize(matrix_a)   

!  Cleanu-up
      CALL dbcsr_release(desym_matrix_s)
      CALL dbcsr_release(desym_matrix_ks)
      CALL dbcsr_release(work_matrix)
   
      CALL timestop(handle)

   END SUBROUTINE build_gs_contribution

! **************************************************************************************************
!> \brief Creates the metric (aka  matrix G) needed for the generalized eigenvalue problem
!>         => G_{pi,qj} = S_pq*delta_ij
!> \param matrix_g dbcsr matrix containing G
!> \param dist_g dbcsr distribution of G
!> \param blk_size_g block sizes of G
!> \param donor_state ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE build_metric(matrix_g, dist_g, blk_size_g, donor_state, qs_env)

      TYPE(dbcsr_type), POINTER                          :: matrix_g
      TYPE(dbcsr_distribution_type), POINTER             :: dist_g
      INTEGER, DIMENSION(:), POINTER                     :: blk_size_g
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_metric', routineP = moduleN//":"//routineN
      
      INTEGER                                         :: handle, i, iblk, jblk, ndo_mo, nblk_row, &
                                                         blk
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      INTEGER, DIMENSION(:), POINTER                  :: row_blk_size
      REAL(dp), DIMENSION(:), POINTER                 :: work_block
      LOGICAL                                         :: found_block
      TYPE(dbcsr_iterator_type)                       :: iter     

      CALL timeset(routineN, handle)

      NULLIFY (matrix_s, row_blk_size, work_block)

!  Initilization
      ndo_mo = SIZE(donor_state%mo_indices, 1)
      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size)
      nblk_row = SIZE(row_blk_size)

!  Creating the symmetric  matrix G with the right size and distribution
      CALL dbcsr_create(matrix=matrix_g, name="MATRIX G", matrix_type="S", dist=dist_g, &
                        row_blk_size=blk_size_g, col_blk_size=blk_size_g)

!  Fill the matrix G by looping over the block of S and putting them on the diagonal
      CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

!        Get the block
         found_block = .FALSE.
         CALL dbcsr_get_block_p(matrix_s(1)%matrix, iblk, jblk, work_block, found_block)

         IF (found_block) THEN

!           Go over the diagonal of G => donor MOs ii
            DO i = 1,ndo_mo
               CALL dbcsr_put_block(matrix_g, (i-1)*nblk_row+iblk, (i-1)*nblk_row+jblk, work_block)
            END DO

         END IF
         NULLIFY(work_block)

      END DO ! dbcsr_iterator
      CALL dbcsr_iterator_stop(iter)

!  Finalize
      CALL dbcsr_finalize(matrix_g)

      CALL timestop(handle)

   END SUBROUTINE build_metric

! **************************************************************************************************
!> \brief Diagonalises a non-symmetric real square matrix. Returns the the assumed real eigenvectors
!>        and eigenvalues. This is base on the LAPACK library and is not paralellized
!> \param matrix the full matrix to diagonalize
!> \param evals the array containing the computed eigenvalues
!> \param evecs the square matrix where the eigenvectors are stored
!> \note the input matrix is changed during the subroutine
! **************************************************************************************************
   SUBROUTINE diag_nonsymm_mat(matrix, evals, evecs)

      TYPE(cp_fm_type), POINTER                       :: matrix
      REAL(dp), DIMENSION(:), POINTER                 :: evals
      TYPE(cp_fm_type), POINTER                       :: evecs

      CHARACTER(len=*), PARAMETER :: routineN = "diag_nonsymm_mat", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: order, lda, ilo, ihi, info, lwork, m, &
                                                         handle
      REAL(dp), DIMENSION(:,:), POINTER               :: qmat, vr, wmat
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: scal, work, tau, wi
      LOGICAL, DIMENSION(:), ALLOCATABLE              :: sel

      NULLIFY(qmat, vr, wmat)

      CALL timeset(routineN, handle)

!  Initialization
      CPASSERT(ASSOCIATED(matrix))
      CPASSERT(ASSOCIATED(evals))
      CPASSERT(ASSOCIATED(evecs))

!  Saqure matrix => order = leading dimension (lda)
      CALL cp_fm_get_info(matrix, nrow_global=order)
      lda = order
      wmat => matrix%local_data
!     Not implemented for parallel run yet !
      CPASSERT(SIZE(wmat,1) == order)
      CPASSERT(SIZE(wmat,2) == order)

!  Balancing the matrix, to make it block upper diagonal: gebal subroutine
      ALLOCATE(scal(lda))
      CALL dgebal("B", order, wmat, lda, ilo, ihi, scal, info)
      IF (info .NE. 0 ) CPABORT("gebal failed in non-symmetric matrix diagonalization")

!  Reduction to upper-Hessenberg form: gehrd
      lwork = -1
      ALLOCATE(work(1))
      ALLOCATE(tau(order-1))

!     First call to determine optimal size of work array
      CALL dgehrd(order, ilo, ihi, wmat, lda, tau, work, lwork, info)
      IF (info .NE. 0 ) CPABORT("gehrd failed in non-symmetric matrix diagonalization")

      lwork = INT(work(1))
      DEALLOCATE(work)
      ALLOCATE(work(lwork))

!     Actual matrix reduction
      CALL dgehrd(order, ilo, ihi, wmat, lda, tau, work, lwork, info) 
      IF (info .NE. 0 ) CPABORT("gehrd failed in non-symmetric matrix diagonalization") 

!  Compute the matrix Q that reduces the matrix to upper-Hessenberg: orghr
      ALLOCATE(qmat(order,order)) 
      qmat = wmat
      CALL dorghr(order, ilo, ihi, qmat, lda, tau, work, lwork, info)
      IF (info .NE. 0 ) CPABORT("orghr failed in non-symmetric matrix diagonalization")

!  Compute eigenvalues and Schur factorization: hseqr
      CPASSERT(SIZE(evals) == order)
      ALLOCATE(wi(order))
      CALL dhseqr('S', 'V', order, ilo, ihi, wmat, lda, evals, wi, qmat, lda, work, lwork, info)
      IF (info .NE. 0 ) CPABORT("hseqr failed in non-symmetric matrix diagonalization")

!  Compute eigenvectors: trevc
      ALLOCATE(vr(order,order))
      ALLOCATE(sel(1))
      CALL dtrevc('R', 'A', sel, order, wmat, lda, vr, lda,  vr, lda, order, m, work, info)
      IF (info .NE. 0 ) CPABORT("trevc failed in non-symmetric matrix diagonalization")

!  Transform eigenvectors to those of the original matrix: gebak
      vr = MATMUL(qmat,vr)
      CALL dgebak("B", "R", order, ilo, ihi, scal, order, vr, lda, info)
      IF (info .NE. 0 ) CPABORT("gebak failed in non-symmetric matrix diagonalization") 

!  Return the eigenvectors
      CPASSERT(SIZE(evecs%local_data,1) == order)
      CPASSERT(SIZE(evecs%local_data,2) == order)
      evecs%local_data = vr

! Clean-up
      DEALLOCATE(qmat, vr)

      CALL timestop(handle)

   END SUBROUTINE diag_nonsymm_mat


! **************************************************************************************************
!> \brief Solves the non-symmetric eigenvalue problem in parallel using ScaLapack. It is assumed
!>        that the eigenvectors and eigenvalues are real
!> \param matrix the full matrix to be diagonalised
!> \param evals the array the will contain the real eigenvalues as output
!> \param evecs the full matrix containing the real eigenvectors
!> \note  The input matrix is modified during execution
! **************************************************************************************************
   SUBROUTINE pdiag_nonsymm_mat(matrix, evals, evecs)

      TYPE(cp_fm_type), POINTER                       :: matrix 
      REAL(dp), DIMENSION(:), POINTER                 :: evals    
      TYPE(cp_fm_type), POINTER                       :: evecs       

      CHARACTER(len=*), PARAMETER :: routineN = "pdiag_nonsymm_mat", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, order, ilo, ihi, ia, ja, info, &
                                                         lwork, m, liwork
      REAL(dp), DIMENSION(:,:), POINTER               :: amat, wmat
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: work, tau, wi
      COMPLEX(dp), DIMENSION(:), ALLOCATABLE          :: cwork
      COMPLEX(dp), DIMENSION(:,:), POINTER            :: cmat, cevecs
      LOGICAL, DIMENSION(:), ALLOCATABLE              :: selec
      INTEGER, DIMENSION(:), ALLOCATABLE              :: iwork
      INTEGER, DIMENSION(9)                           :: desca, descw, descvr
      TYPE(cp_fm_type), POINTER                       :: work_mat, work_evecs
      TYPE(cp_cfm_type), POINTER                      :: comp_mat, comp_evecs

      NULLIFY(amat, work_mat, wmat, comp_mat, comp_evecs, cmat, cevecs, work_evecs)

!  We assume SCALAPACK is available
      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(matrix))
      CPASSERT(ASSOCIATED(evals))
      CPASSERT(ASSOCIATED(evecs))

!  Initialization
      order = matrix%matrix_struct%nrow_global 
      amat => matrix%local_data
      desca = matrix%matrix_struct%descriptor
!     intereseted in the whole matrix
      ilo = 1
      ihi = order
      ia = 1
      ja = 1

!  Workspace querry
      info = -1
      lwork = -1
      ALLOCATE(work(1))
      ALLOCATE(tau(ja+order-2))

#if defined(__SCALAPACK)
      CALL pdgehrd(order, ilo, ihi, amat(1,1), ia, ja, desca, tau, work(1), lwork, info)
#endif

      IF (info .NE. 0) CPABORT("Routine pdgehrd of ScaLAPACK failed for non-symmetric matrix diag")
      lwork = INT(work(1))
      DEALLOCATE(work)
      ALLOCATE(work(lwork))

!  Actual call to pdgehrd => Reduction to upper Hessenberg form
#if defined(__SCALAPACK)                                                                             
      CALL pdgehrd(order, ilo, ihi, amat(1,1), ia, ja, desca, tau(1), work(1), lwork, info)                     
#endif

      IF (info .NE. 0) CPABORT("Routine pdgehrd of ScaLAPACK failed for non-symmetric matrix diag")

!  Need to get the orthogonal transormation matrix Q, do that by multiplying the ouput of pdgehrd by
!  identity (hack since output of pdgehrd is in a very weird format)
      CALL cp_fm_create(work_mat, matrix%matrix_struct)
      CALL cp_fm_set_all(work_mat, alpha=0.0_dp, beta=1.0_dp)
      wmat => work_mat%local_data
      descw = work_mat%matrix_struct%descriptor

!  First work space querry
      DEALLOCATE(work)
      ALLOCATE(work(1))
      lwork = -1

#if defined(__SCALAPACK) 
      CALL pdormhr('R', 'N', order, order, ilo, ihi, amat(1,1), ia, ja, desca, tau(1), wmat(1,1), ia, ja, descw,&
                   work(1), lwork, info)
#endif

      IF (info .NE. 0 ) CPABORT("Routine pdormhr of ScaLAPACK failed for non-symmetric matrix diag")

      lwork = INT(work(1))
      DEALLOCATE(work)
      ALLOCATE(work(lwork))

!  Actual call of pdormhr

#if defined(__SCALAPACK)                                                                             
      CALL pdormhr('R', 'N', order, order, ilo, ihi, amat(1,1), ia, ja, desca, tau(1), wmat(1,1), ia, ja, descw,& 
                   work(1), lwork, info)                                                                
#endif 

      IF (info .NE. 0 ) CPABORT("Routine pdormhr of ScaLAPACK failed for non-symmetric matrix diag")
      DEALLOCATE(work)
      DEALLOCATE(tau)

!  Compute eigenvalues and Schur factorization using pdlaqr1
      CPASSERT(SIZE(evals) == order)

!  workspace querry
      lwork = -1
      liwork = -1
      ALLOCATE(work(1))
      ALLOCATE(iwork(1)) 
      ALLOCATE(wi(order))

#if defined(__SCALAPACK)  
      CALL pdlaqr1(.TRUE., .TRUE., order, ilo, ihi, amat, desca, evals, wi, ilo, ihi, wmat, &
                   descw, work, lwork, iwork, liwork, info)
#endif
      
      IF (info .NE. 0 ) CPABORT("Routine pdlaqr1 of ScaLAPACK failed for non-symmetric matrix diag")

      lwork = INT(work(1))
      liwork = INT(iwork(1)) 
      DEALLOCATE(work, iwork)         
      ALLOCATE(work(lwork))
      ALLOCATE(iwork(liwork))

!  Actual call to pdlaqr1

#if defined(__SCALAPACK)                                                                             
      CALL pdlaqr1(.TRUE., .TRUE., order, ilo, ihi, amat, desca, evals, wi, ilo, ihi, wmat, &
                   descw, work, lwork, iwork, liwork, info)                                    
#endif
      IF (info .NE. 0 ) CPABORT("Routine pdlaqr1 of ScaLAPACK failed for non-symmetric matrix diag")
            
      DEALLOCATE(iwork)
      DEALLOCATE(wi)

!  Now get the eigenvectors. Note: we have a real matrix and expect real eigenvectors. However,
!  ScalLAPACK only has such a routine for complex cases => need to convert
      
!     real to complex conversion, keep same matrix structure
      CALL cp_cfm_create(comp_mat, matrix%matrix_struct)
      CALL cp_cfm_create(comp_evecs, evecs%matrix_struct)
      CALL cp_fm_to_cfm(msourcer=matrix, mtarget=comp_mat)
      descvr = comp_evecs%matrix_struct%descriptor
      cmat => comp_mat%local_data
      cevecs => comp_evecs%local_data

!     other initializations
      ALLOCATE(selec(1))
      selec(1) = .TRUE.
      DEALLOCATE(work)
      ALLOCATE(work(desca(9))) !local leading dimension
      ALLOCATE(cwork(2*desca(9)))
      m = 0

!  todo: need to properly check the performances of this guy.. Apparently slow
!        also, can only compute a subset of the eigenvectors. Maybe do that to

#if defined(__SCALAPACK) 
      CALL pztrevc('R', 'A', selec(1), order, cmat(1,1), desca, cevecs(1,1), descvr, cevecs(1,1), descvr, order, m, &
                   cwork(1), work(1), info)
#endif

      IF (info .NE. 0 ) CPABORT("Routine pztrevc of ScaLAPACK failed for non-symmetric matrix diag")

!  Go back to real matrices (assume no imaginary part)
      CALL cp_fm_create(work_evecs, evecs%matrix_struct)
      CALL cp_cfm_to_fm(msource=comp_mat, mtargetr=matrix)
      CALL cp_cfm_to_fm(msource=comp_evecs, mtargetr=work_evecs)

!  To obtain the proper eigenvectors, need to compute evecs = Q*evecs, Q stored in work_mat
      CALL cp_fm_gemm('N', 'N', order, order, order, 1.0_dp, work_mat, work_evecs, 0.0_dp, evecs)

!  Clean-up
      CALL cp_fm_release(work_mat)
      CALL cp_fm_release(work_evecs)
      CALL cp_cfm_release(comp_mat)
      CALL cp_cfm_release(comp_evecs)

      CALL timestop(handle)

   END SUBROUTINE pdiag_nonsymm_mat

! **************************************************************************************************!
!> \brief Computes the two-center Coulomb integral needed for the RI in kernel calculation. Given
!>        the excited kind indices, returns the (for now) cp_fm (P|Q)^-1 for the RI basis P,Q
!> \param ex_kinds the indices of the excited kinds
!> \param nao the number of AOs
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \param xas_tdp_control ...
!> \note Stores  (P|Q)^-1 in the xas_tdp_env, for now, uses eri_mme and only works for PBCs
! **************************************************************************************************
   SUBROUTINE compute_ri_coulomb2_int(xas_tdp_env, qs_env, xas_tdp_control)

      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(qs_environment_type), POINTER              :: qs_env 
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control  

      CHARACTER(len=*), PARAMETER :: routineN = "compute_ri_coulomb2_int", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, nset, pset, qset, G_count, & 
                                                         R_count, ikind, i, offset_p_start, &
                                                         offset_p_end, offset_q_start, &
                                                         offset_q_end, sgfp, sgfq, eri_method
      TYPE(gto_basis_set_type), POINTER               :: ri_basis
      REAL(dp), DIMENSION(:,:), POINTER               :: zet, sphi
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set        
      INTEGER, DIMENSION(:), POINTER                  :: lmax, lmin, ex_kinds, nsgf, npgf
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgf
      TYPE(cp_eri_mme_param), POINTER                 :: eri_param
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cell_type), POINTER                        :: cell 
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: tmp_int
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(cp_fm_type), POINTER                       :: work_fm  
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct 
      REAL(dp), DIMENSION(3)                          :: r

      CALL timeset(routineN, handle)
      
      NULLIFY(ri_basis, zet, sphi, qs_kind_set, lmax, lmin, first_sgf, eri_param, para_env)
      NULLIFY(work_fm, fm_struct, blacs_env)
      CPASSERT(ASSOCIATED(xas_tdp_env))

!  Initialization
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, cell=cell, para_env=para_env, blacs_env=blacs_env)
      ex_kinds => xas_tdp_env%ex_atom_indices
      ALLOCATE(xas_tdp_env%ri_inv_mats(SIZE(ex_kinds)))

!  Check on the periodicity to choose the eri method
      IF (ALL(cell%perd == 1)) eri_method = do_eri_mme
      IF (ALL(cell%perd == 0)) eri_method = do_eri_os
      

!  Reading the ERI_MME params and calibrate them (not centralized to avoid circular dependency)
      ALLOCATE(eri_param)
      CALL cp_eri_mme_init_read_input(xas_tdp_control%eri_mme_subsection, eri_param)
      CALL cp_eri_mme_set_params(eri_param, cell, qs_kind_set, basis_type_1="RI_AUX", &
                                 basis_type_2="RI_AUX", para_env=para_env)      
      G_count = 0
      R_count = 0

!  Loop over the excited kinds
      DO i = 1,SIZE(ex_kinds)
         ikind = ex_kinds(i) 

!        Get the RI basis of interest and its quantum numbers
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_AUX")
         CALL get_gto_basis_set(ri_basis, lmax=lmax, npgf=npgf, zet=zet, lmin=lmin, nset=nset, &
                                nsgf_set=nsgf, sphi=sphi, first_sgf=first_sgf) 

!        The 2D array with the integrals for this kind
         ALLOCATE(tmp_int(SUM(nsgf,1), SUM(nsgf,1)))
         tmp_int = 0.0_dp

!        Loop over the sets of the basis
!        TODO: can probably avoid to go over all sets by using some symmetry

         offset_p_end = 0

         DO pset = 1,nset
            offset_p_start = offset_p_end
            offset_p_end = offset_p_end+nsgf(pset)
            offset_q_end = 0
            sgfp = first_sgf(1, pset)
            
            DO qset = 1,nset
               offset_q_start = offset_q_end 
               offset_q_end = offset_q_end+nsgf(qset) 
               sgfq = first_sgf(1, qset)

               r = 0.0_dp
               CALL integrate_set_2c(eri_param%par, lmin(pset), lmax(pset), lmin(qset), lmax(qset),&
                                     npgf(pset), npgf(qset), zet(:,pset), zet(:,qset), r, r, tmp_int, &
                                     nsgf(pset), nsgf(qset), offset_p_start, offset_q_start, 0, 0,&
                                     sphi, sphi, sgfp, sgfq, nsgf(pset), nsgf(qset), eri_method, &
                                     G_count=G_count, R_count=R_count)

            END DO ! qset
         END DO !pset

!        Storing it into the xas_tdp_env
         CALL cp_fm_struct_create(fm_struct, context=blacs_env, para_env=para_env, & 
                                  nrow_global=SUM(nsgf,1), ncol_global=SUM(nsgf,1))

         CALL cp_fm_create(work_fm, fm_struct)
         CALL cp_fm_set_submatrix(work_fm, tmp_int)

         CALL cp_fm_create(xas_tdp_env%ri_inv_mats(i)%matrix, fm_struct)
         CALL cp_fm_invert(work_fm, xas_tdp_env%ri_inv_mats(i)%matrix) 

!        clean-up before next kind
         CALL cp_fm_struct_release(fm_struct)   
         CALL cp_fm_release(work_fm)       
         DEALLOCATE(tmp_int)

      END DO ! loop on kinds

!  Clean-up
      CALL cp_eri_mme_update_local_counts(eri_param, para_env, G_count_2c=G_count, &
                                          R_count_2c=R_count)
      CALL cp_eri_mme_finalize(eri_param)
      DEALLOCATE(eri_param)

      CALL timestop(handle)

   END SUBROUTINE compute_ri_coulomb2_int

! **************************************************************************************************
!> \brief Computes the three-centers Coulomb integrals needed for RI in kernel calculations and 
!>        stores it as 3D arrays. Computes the integrals (ab|P), where a,b are overlaping atomic
!>        basis set and P is the RI basis centered on the given excited atom. The first dimension
!>        of the returned array corresponds to sgf a, the second to sgf b and the third to sgf P.
!>        Keep the result in xas_tdp_env as the integrals for the current excited atom
!> \param exat_index the index of the excited atom, on which the RI basis P is centered
!> \param xas_tdp_env the xas_tdp_Env to update with the 3-center Coulomb integrals
!> \param qs_env ...
!> \param xas_tdp_control ...
!> \note Currently creates a 3D array, will later switch to DBCSR tensors (TODO). Will be used for 
!>       the Coulomb kernel after contraction (aI|P) = sum_b C_b (ab|P) and possibly for the exchange
!>       part too (if non-PBCs and coulomb operator requested). TODO: if exchange operator is not
!>       coulomb, screen the integrals to only compute (ab|P) where C_b is non-zero
!>       Based on the 3-center ERI_MME integrals
! **************************************************************************************************
   SUBROUTINE compute_ri_coulomb3_int(exat_index, xas_tdp_env, qs_env, xas_tdp_control)

      INTEGER, INTENT(IN)                             :: exat_index
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_coulomb3_int', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                         :: handle, nsetp, nsetb, natom, &
                                                         nao, nseta, bat, aset, bset, pset, &
                                                         sgfa, sgfp, sgfb, aat, akind, bkind, & 
                                                         nkind, exat_kind, idx, igp, i, GG_count, &
                                                         GR_count, RR_count, inode, last_bat, &
                                                         eri_method
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER               :: p_basis, b_basis, a_basis
      INTEGER, DIMENSION(:), POINTER                  :: lp_max, lp_min, npgfp, lb_max,lb_min,&
                                                         npgfb, la_max, la_min, npgfa, nsgfa, &
                                                         nsgfb, nsgfp
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgfa, first_sgfb, first_sgfp
      REAL(dp), DIMENSION(:,:), POINTER               :: zetp, rpgfp, zetb, rpgfb, zeta, rpgfa,&
                                                         sphi_p, sphi_a, sphi_b
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos 
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE         :: int_abp_set, int_array
      REAL(dp), DIMENSION(3)                          :: ra, rb, rp, rab
      INTEGER, DIMENSION(:), ALLOCATABLE              :: first_sgf_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), & 
                                          POINTER     :: sab_orb
      TYPE(neighbor_list_iterator_p_type), &                                                         
                            DIMENSION(:), POINTER     :: nl_iterator   
      TYPE(cell_type), POINTER                        :: cell
      TYPE(cp_eri_mme_param), POINTER                 :: eri_param   
      TYPE(cp_para_env_type), POINTER                 :: para_env
      REAL(dp)                                        :: dap, dbp

      CALL timeset(routineN, handle)

      NULLIFY(qs_kind_set, p_basis, b_basis, a_basis, particle_set, mos)
      NULLIFY(lp_max, lp_min, npgfp, zetp, rpgfp, sphi_p, first_sgfp)
      NULLIFY(lb_max, lb_min, npgfb, zetb, rpgfb, sphi_b, first_sgfb)
      NULLIFY(la_max, la_min, npgfa, zeta, rpgfa, sphi_a, first_sgfa)
      NULLIFY(sab_orb, nl_iterator, cell, basis_set_list, eri_param, para_env)

!  Get all we need from the qs_env
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, mos=mos, &
                      sab_orb=sab_orb, cell=cell, natom=natom, para_env=para_env)
      CALL get_mo_set(mos(1)%mo_set, nao=nao)
      ALLOCATE(first_sgf_set(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf_set)

!  Retrieve the basis set P centered on the excited atom
      idx = MINLOC(ABS(xas_tdp_env%ex_atom_indices-exat_index), 1) 
      exat_kind = xas_tdp_env%ex_kind_indices(idx)
      CALL get_qs_kind(qs_kind_set(exat_kind), basis_set=p_basis, basis_type="RI_AUX")

      CALL get_gto_basis_set(p_basis, lmax=lp_max, npgf=npgfp, zet=zetp, lmin=lp_min, nset=nsetp, &
                             nsgf_set=nsgfp, sphi=sphi_p, first_sgf=first_sgfp, pgf_radius=rpgfp) 

!  Check on periodicity to impose ERI method
      IF (ALL(cell%perd == 1)) eri_method = do_eri_mme                                               
      IF (ALL(cell%perd == 0)) eri_method = do_eri_os

!  Allocate the 3D array where results will be stored
      ALLOCATE(int_array(nao, nao, SUM(nsgfp)))
      int_array = 0.0_dp

!  Create lists of basis for atoms a and b
      nkind = SIZE(qs_kind_set)
      ALLOCATE(basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

!  Get the mme params form the ERI_MME subsection. Do that here to avoid circular dependency
!  TODO: is that expensive, is there a way around ?
      ALLOCATE(eri_param)
      CALL cp_eri_mme_init_read_input(xas_tdp_control%eri_mme_subsection, eri_param)
      CALL cp_eri_mme_set_params(eri_param, cell, qs_kind_set, basis_type_1="ORB", & 
                                              basis_type_2="RI_AUX", para_env=para_env)
      GG_count = 0
      GR_count = 0
      RR_count = 0

!  Loop over all neighbouring atoms a,b and compute (ab|P)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=akind, jkind=bkind, iatom=aat, jatom=bat, r=rab,&
                                inode=inode)

!        Discard periodic images because ERI_MME intrinsically periodic
         IF (inode == 1) last_bat = 0

         IF (bat /= last_bat) THEN
            last_bat = bat
         ELSE
            CYCLE
         END IF

!        Retrieve the basis sets for a and b
         a_basis => basis_set_list(akind)%gto_basis_set
         b_basis => basis_set_list(bkind)%gto_basis_set
         
!        Retrieve the quantum numbers for both basis sets
         CALL get_gto_basis_set(a_basis, lmax=la_max, npgf=npgfa, zet=zeta, lmin=la_min, &
                                nset=nseta, nsgf_set=nsgfa, sphi=sphi_a, first_sgf=first_sgfa, &
                                pgf_radius=rpgfa)   
         CALL get_gto_basis_set(b_basis, lmax=lb_max, npgf=npgfb, zet=zetb, lmin=lb_min, &
                                nset=nsetb, nsgf_set=nsgfb, sphi=sphi_b, first_sgf=first_sgfb, &
                                pgf_radius=rpgfb)         

!        Get the different positions/distances
         ra = pbc(particle_set(aat)%r, cell) 
         rb = pbc(particle_set(bat)%r, cell)         
         rp = pbc(particle_set(exat_index)%r, cell)

!        Loop over sets of the 3 basis
         DO aset = 1,nseta
            DO bset = 1,nsetb
               DO pset = 1,nsetp

!                 Abort if coulomb3 subroutine screens stuff out, abort if more than one pgf per 
!                 set in RI basis
                  IF (eri_method == do_eri_os) THEN
                     IF (SIZE(rpgfp,1) .GT. 1) THEN
                        CPABORT("For OS scheme in NPBCs, only one pgf per set is allowed in the RI basis")
                     END IF

                     dap = SQRT((rp(1)-ra(1))**2+(rp(2)-ra(2))**2+(rp(3)-ra(3))**2)
                     dbp = SQRT((rp(1)-rb(1))**2+(rp(2)-rb(2))**2+(rp(3)-rb(3))**2) 

                     IF (dap .GT. MINVAL(rpgfa(:,aset))+rpgfp(1,pset) .OR. &
                         dbp .GT. MINVAL(rpgfb(:,bset))+rpgfp(1,pset)) THEN
                        CPABORT("Coulomb3 subroutine screens some if the ERI integrals...")
                     END IF

                  END IF

!                 Some initialization
                  sgfa = first_sgfa(1, aset) 
                  sgfb = first_sgfb(1, bset)
                  sgfp = first_sgfp(1, pset)   

!                 Allocate array for the current comnbination of sets
                  ALLOCATE(int_abp_set(nsgfa(aset), nsgfb(bset), nsgfp(pset)))
                  int_abp_set = 0.0_dp

!                 Took inspiration form mp2. Integrals are already contracted as sgf
!                 the do_symmetric argument only matters for forces
                  CALL integrate_set_3c(eri_param%par, la_min(aset), la_max(aset), lb_min(bset), &
                                        lb_max(bset), lp_min(pset), lp_max(pset), npgfa(aset), &
                                        npgfb(bset), npgfp(pset), zeta(:, aset), zetb(:, bset), &
                                        zetp(:, pset), ra, rb, rp, int_abp_set, nsgfa(aset), &
                                        nsgfb(bset), nsgfp(pset), 0, 0, 0, 0, 0, 0, sphi_a, &
                                        sphi_b, sphi_p, sgfa, sgfb, sgfp, nsgfa(aset), &
                                        nsgfb(bset), nsgfp(pset), eri_method, do_symmetric=.FALSE.,&
                                        GG_count=GG_count, GR_count=GR_count, RR_count=RR_count)

!                 Put it at the right place in the final array
!                 TODO: in case of PBCs, might need to do += since pair a,b might involve images
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
                  DEALLOCATE(int_abp_set)

               END DO ! pset
            END DO ! bset
         END DO !aset
      
      END DO !neighbor iteration on a,b
      CALL neighbor_list_iterator_release(nl_iterator)

!  Clean the eri params 
      CALL cp_eri_mme_update_local_counts(eri_param, para_env, GG_count_3c=GG_count, &
                                          GR_count_3c=GR_count, RR_count_3c=RR_count)
      CALL cp_eri_mme_finalize(eri_param)
      DEALLOCATE(eri_param)

!  Updating the xas_tdp_env with the new integrals
      IF (ASSOCIATED(xas_tdp_env%ri_coul3)) DEALLOCATE(xas_tdp_env%ri_coul3)
      ALLOCATE(xas_tdp_env%ri_coul3(nao, nao, SUM(nsgfp)))
      xas_tdp_env%ri_coul3 = int_array

!  Clean-up
      DEALLOCATE(basis_set_list)

      CALL timestop(handle)

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
         contr_once(:,:) = contr_once(:,:) + coeffs(i,1)*integrals(:,i,:) 
      END DO

!  Do the second contraction if requested
      IF (contract_two) THEN
         contr_twice = 0.0_dp
         CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=coeffs, start_row=1, &
                                  start_col=domo_index(2), n_rows=nao, n_cols=1, transpose=.FALSE.)

         DO i = 1,nao
            contr_twice(1,:) = contr_twice(1,:) + coeffs(i,1)*contr_once(i,:)
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
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: int_aI, int_bJ, PQ, proj_Q, work1, work2, &
                                                         aIJb_prod
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: context
      TYPE(cp_fm_type), POINTER                       :: work_fm     
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct        
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set      
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      INTEGER, DIMENSION(:), ALLOCATABLE              :: first_sgf

      CALL timeset(routineN, handle)

      NULLIFY(para_env, context, work_fm, fm_struct, qs_kind_set, particle_set)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices,1)                       
      CALL get_qs_env(qs_env, natom=nblk, particle_set=particle_set, qs_kind_set=qs_kind_set)
      ALLOCATE(first_sgf(nblk))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf)

!     Get the (P|Q)^-1 as an array
      i = MINLOC(ABS(xas_tdp_env%ex_kind_indices-donor_state%kind_index), 1)
      CALL cp_fm_get_info(xas_tdp_env%ri_inv_mats(i)%matrix, ncol_global=nsgfp, &
                          context=context, para_env=para_env) 
      nao = SIZE(xas_tdp_env%ri_coul3, 1)
      ALLOCATE(PQ(nsgfp, nsgfp))
      CALL cp_fm_get_submatrix(fm=xas_tdp_env%ri_inv_mats(i)%matrix, target_m=PQ)

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
      CALL dbcsr_create(matrix=coulomb_kernel, name="COULOMB KERNEL", matrix_type='N', dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)

!     Allocate the work arrays to store the intermediate matrix products
      ALLOCATE(work1(nao, nsgfp))
      ALLOCATE(work2(nao,nao))
      ALLOCATE(aIJb_prod(nao,nao))

!  Loop over donor MOs and compute W_ab = (aI|P)* (P|Q)^-1 * (Q|Jb), then take projQ*W
!  In principle IJ and JI are the same. Here we span everything for safety though
      DO i = 1,ndo_mo
         CALL contract_int_kernel(int_aI, (/donor_state%mo_indices(i)/), xas_tdp_env, qs_env)
         DO j = 1,ndo_mo
            CALL contract_int_kernel(int_bJ, (/donor_state%mo_indices(j)/), xas_tdp_env, qs_env)

            work1 = 0.0_dp
            work2 = 0.0_dp
            aIJb_prod  = 0.0_dp

!           Do the product (aI|P)*(P|Q)^-1
            CALL dgemm('N', 'N', nao, nsgfp, nsgfp, 1.0_dp, int_aI, nao, PQ, nsgfp, 0.0_dp, &
                        work1, nao)

!           Mulitply by (Q|Jb)
            CALL dgemm('N', 'T', nao, nao, nsgfp, 1.0_dp, work1, nao, int_bJ, nao, 0.0_dp, &
                        work2, nao)

!           Multiply by the projector
            CALL dgemm('N', 'N', nao, nao, nao,  1.0_dp, proj_Q, nao, work2, nao, 0.0_dp, &
                       aIJb_prod, nao)

!           Put it in the dbcsr matrix, for the right I,J and block by block
            DO iblk = 1,nblk
               DO jblk = 1,nblk

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
!> \param ex_ker_offdiag the dbcsr_matrix where the  off-diagonal exchange kernel is stored
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

      TYPE(dbcsr_type), POINTER                       :: ex_ker_diag, ex_ker_offdiag  
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
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: PQ, proj_Q, int_aJ, int_bI, int_IJ, work1,&
                                                         work2, work3, abIJ_prod, aJIb_prod, work4
      TYPE(cp_fm_type), POINTER                       :: work_fm                      
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct      

      CALL timeset(routineN, handle)

      NULLIFY(particle_set, qs_kind_set, para_env, context, work_fm, fm_struct)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices,1)
      CALL get_qs_env(qs_env, natom=nblk, particle_set=particle_set, qs_kind_set=qs_kind_set)       
      ALLOCATE(first_sgf(nblk))        
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf)

!     Get the (P|Q)^-1 as an array      
      i = MINLOC(ABS(xas_tdp_env%ex_kind_indices-donor_state%kind_index), 1)          
      CALL cp_fm_get_info(xas_tdp_env%ri_inv_mats(i)%matrix, ncol_global=nsgfp, &
                          context=context, para_env=para_env)                                 
      nao = SIZE(xas_tdp_env%ri_coul3, 1)    
      ALLOCATE(PQ(nsgfp, nsgfp))
      CALL cp_fm_get_submatrix(fm=xas_tdp_env%ri_inv_mats(i)%matrix, target_m=PQ)               

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
      CALL dbcsr_create(matrix=ex_ker_diag, name="EX KER DIAG", matrix_type='N', dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size) 
      CALL dbcsr_create(matrix=ex_ker_offdiag, name="EX KER OFFDIAG", matrix_type='N', dist=dist, &  
                        row_blk_size=blk_size, col_blk_size=blk_size)

!     Allocate the many work arrays
      ALLOCATE(work1(nao, nsgfp)) 
      ALLOCATE(work2(nao,nao))
      ALLOCATE(work3(nsgfp,1))
      ALLOCATE(work4(nao,nao))
      ALLOCATE(abIJ_prod(nao,nao))
      ALLOCATE(aJIb_prod(nao,nao))

!  Loop over the donor MOs and compute both the digonal and off-digonal parts, that is:
!     -Diag_ab = (ab|P)* (P|Q)^-1 *(Q|IJ), then multiply by proj_Q from the left
!     -Off_ab = (aJ|P)* (P|Q)^-1 *(Q|Ib), then multiply by proj_Q from the left
      DO i = 1,ndo_mo
         CALL contract_int_kernel(int_bI, (/donor_state%mo_indices(i)/), xas_tdp_env, qs_env)
         DO j = 1,ndo_mo
            CALL contract_int_kernel(int_aJ, (/donor_state%mo_indices(j)/), xas_tdp_env, qs_env)
            CALL contract_int_kernel(int_IJ, (/donor_state%mo_indices(i), donor_state%mo_indices(j)/), &
                                     xas_tdp_env, qs_env)

            work1 = 0.0_dp
            work2 = 0.0_dp
            work3 = 0.0_dp
            work4 = 0.0_dp
            abIJ_prod = 0.0_dp
            aJIb_prod = 0.0_dp

!           Start with off-diag product (aJ|P)* (P|Q)^-1 *(Q|Ib)

!           (aJ|P)* (P|Q)^-1
            CALL dgemm('N', 'N', nao, nsgfp, nsgfp, 1.0_dp, int_aJ, nao, PQ, nsgfp, 0.0_dp, &
                        work1, nao)
!           Muliply by (Q|Ib)  
            CALL dgemm('N', 'T', nao, nao, nsgfp, 1.0_dp, work1, nao, int_bI, nao, 0.0_dp, &
                       work2, nao)
!           Multiply that by the Q projector
            CALL dgemm('N', 'N', nao, nao, nao, 1.0_dp, proj_Q, nao, work2, nao, 0.0_dp, &
                       aJIb_prod, nao)

!           Then the on-digonal product (ab|P)* (P|Q)^-1 *(Q|IJ) 

!           (P|Q)^-1 *(Q|IJ)
            CALL dgemm('N', 'T', nsgfp, 1, nsgfp, 1.0_dp, PQ, nsgfp, int_IJ, 1, 0.0_dp, &
                       work3, nsgfp)

!           Multiply by (ab|P) from the left 
            DO k = 1,nsgfp
               work4(:,:) = work4(:,:) + xas_tdp_env%ri_coul3(:,:,k)*work3(k,1)
            END DO !k

!           Multiply by the Q projector from the left
            CALL dgemm('N', 'N', nao, nao, nao, 1.0_dp, proj_Q, nao, work4, nao, 0.0_dp, &
                       abIJ_prod, nao)

!           Put it all in the dbcsr matrices, at the right I,J and block by block
            DO iblk = 1,nblk
               DO jblk = 1,nblk

!                 Size of the current block
                  nrow_blk = blk_size(iblk)
                  ncol_blk = blk_size(jblk)

!                 Put (ab|IJ) on the digonal exchange kernel block, at position I,J
                  CALL dbcsr_put_block(ex_ker_diag, (i-1)*nblk+iblk, (j-1)*nblk+jblk, & 
                                       abIJ_prod(first_sgf(iblk):first_sgf(iblk)+nrow_blk-1, &
                                       first_sgf(jblk):first_sgf(jblk)+ncol_blk-1))

!                 Put (aJ|Ib) on the off digonal exchange kernel block, at position I,J
                  CALL dbcsr_put_block(ex_ker_offdiag, (i-1)*nblk+iblk, (j-1)*nblk+jblk, &
                                       aJIb_prod(first_sgf(iblk):first_sgf(iblk)+nrow_blk-1, &
                                       first_sgf(jblk):first_sgf(jblk)+ncol_blk-1))     
                  
               END DO !jblk
            END DO !iblk
            
         END DO !j
      END DO !i

      CALL dbcsr_finalize(ex_ker_diag)
      CALL dbcsr_finalize(ex_ker_offdiag)

      CALL timestop(handle)

   END SUBROUTINE build_kernel_exchange

END MODULE xas_tdp_utils

