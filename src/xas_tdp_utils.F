!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for X-ray absorption spectroscopy using TDDFPT
!> \author AB (01.2018)
! **************************************************************************************************

MODULE xas_tdp_utils
   
   USE ai_coulomb,                      ONLY: coulomb3
   USE ai_contraction_sphi,             ONLY: abc_contract
   USE basis_set_types,                 ONLY: gto_basis_set_type, get_gto_basis_set, &
                                              gto_basis_set_p_type
   USE cell_types,                      ONLY: cell_type, pbc
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose, &
                                              cp_dbcsr_cholesky_invert, &
                                              cp_dbcsr_cholesky_restore
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm, cp_dbcsr_sm_fm_multiply, &
                                              copy_fm_to_dbcsr
   USE cp_cfm_types,                    ONLY: cp_fm_to_cfm, cp_cfm_to_fm, cp_cfm_type, &
                                              cp_cfm_create, cp_cfm_release
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add, cp_fm_gemm, cp_fm_column_scale
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
   USE input_constants,                 ONLY: xas_dip_len, xas_dip_vel
   USE iterate_matrix,                  ONLY: invert_Hotelling
   USE kinds,                           ONLY: dp
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
   USE qs_operators_ao,                 ONLY: p_xyz_ao, rRc_xyz_ao
   USE xas_tdp_types,                   ONLY: donor_state_type, xas_tdp_control_type, &
                                              xas_tdp_env_type

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_utils'

   PUBLIC :: setup_xas_tdp_full_prob, solve_xas_tdp_full_prob, compute_oscillator_strength

CONTAINS

! **************************************************************************************************
!> \brief Builds the matrix that define the XAS TDDFPT generalized eigenvalue problem to be solved
!>        for excitation energies omega. The problem has the form omega*C = M*C, where C contains
!>        the reponse orbitals coefficients. The matrix is stored in the given donor_state 
!> \param donor_state the donor_state for which the problem is restricted
!> \param qs_env ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \note the matrix M has the form |-G^-1*(A+B)  -G^-1*B  | where other subroutines builds A, B and G.
!>                                 |   G^-1*B   G^-1*(A+B)|                                       
!>       The matrix M (matrix_tdp) has the same (i.e. 2x2 replicated) dbcsr_dist as usual
!>       If the Tamm-Dancoff approximation is assumed, then C^+ = 0, M = G^-1*(A+B) 
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
                                                         matrix_b, matrix_g_inv, offdiag_mat
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
      NULLIFY(matrix_b, matrix_tdp, offdiag_mat)

      CPASSERT(ASSOCIATED(donor_state))
      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(xas_tdp_env))
      CPASSERT(ASSOCIATED(xas_tdp_control))
      ALLOCATE(donor_state%matrix_tdp)
      matrix_tdp => donor_state%matrix_tdp

!  Get the submatrices distribution and block sizes
      CALL compute_submat_dist_and_blk_size(submat_dist, submat_blk_size, donor_state, qs_env)
      nblk_sub = SIZE(submat_blk_size)

!  Compute the submatrices A, B  and G
      ALLOCATE(matrix_a)
      ALLOCATE(matrix_g)
      ALLOCATE(matrix_b)
      CALL build_matrix_a(matrix_a, submat_dist, submat_blk_size, donor_state, qs_env)
      CALL build_matrix_b(matrix_b, submat_dist, submat_blk_size, donor_state, xas_tdp_env, qs_env)
      CALL build_matrix_g(matrix_g, submat_dist, submat_blk_size, donor_state, qs_env)       

!  Invert matrix G
      ALLOCATE(matrix_g_inv)
      CALL dbcsr_create(matrix=matrix_g_inv, name="MATRIX G INV", matrix_type="S", &              
                        dist=submat_dist, row_blk_size=submat_blk_size, col_blk_size=submat_blk_size)
!     TODO: impose values of threshold and eps_filter or not?. Current values taken from LRI
      CALL invert_Hotelling(matrix_g_inv, matrix_g, threshold=1.0E-10_dp, use_inv_as_guess=.FALSE., &
                            norm_convergence=1.0E-10_dp, filter_eps=1.0E-12_dp, silent=.TRUE.)
      CALL dbcsr_finalize(matrix_g_inv)

!  Create the matrix that will be on the diagonal of matrix_tdp: G^-1*(A+B)
      ALLOCATE(diag_mat)
      CALL dbcsr_copy(diag_mat, matrix_a, name="Diag TDP block")
      CALL dbcsr_add(diag_mat, matrix_b, 1.0_dp, 1.0_dp)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_g_inv, diag_mat, 0.0_dp, diag_mat)
      CALL dbcsr_finalize(diag_mat)

!  Create the matrix that will be in the off-diagonal of matrix_tdp: G^-1*B
      ALLOCATE(offdiag_mat)
      CALL dbcsr_copy(offdiag_mat, matrix_b, name="Offdiag TDP block")
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_g_inv, offdiag_mat, 0.0_dp, offdiag_mat)
      CALL dbcsr_finalize(offdiag_mat)

!  Some clean-up
      CALL dbcsr_release(matrix_g_inv)
      CALL dbcsr_release(matrix_a)
      CALL dbcsr_release(matrix_b)
      CALL dbcsr_release(matrix_g)
      DEALLOCATE(matrix_g_inv, matrix_a, matrix_b, matrix_g)

!  Test on Tamm Dancoff approximation 

      IF (xas_tdp_control%tamm_dancoff) THEN

!        If TDA is assumed, then the matrix_tdp is G^-1*(A+2B)
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

!     Filling the supermatrix matrix_tdp with submatrices +-G^-1*(A+B) and +-G^-1*B

!        For each diagonal block of the super matrix, loop over blocks of G^-1*(A+B)
         CALL dbcsr_iterator_start(iter, diag_mat)
         DO WHILE(dbcsr_iterator_blocks_left(iter))

            found_block = .FALSE.
            CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
            CALL dbcsr_get_block_p(diag_mat, iblk, jblk, work_block, found_block)

            IF (found_block) THEN
!              Put -G^-1*(A+B) in the upper left block and G^-1*(A+B) in the lower right block
               CALL dbcsr_put_block(matrix_tdp, iblk, jblk, -work_block) 
               CALL dbcsr_put_block(matrix_tdp, nblk_sub+iblk, nblk_sub+jblk, work_block)
            END IF
            NULLIFY(work_block)

         END DO ! iterator
         CALL dbcsr_iterator_stop(iter) 

!        For each off-diagonal block of supermatrix, loop over blocks of G^-1*B
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
!> \brief Normalizes the response orbitals so that |C^-T * G * C^- - C^+T * G *C^+| = 1 
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
!> \brief Builds the so called matrix A, which is needed for the full xas_tdp matrix
!>        A_{pi,qj} = F_pq*delta_ij - epsilon_ij*S_pq
!> \param matrix_a  pointer to a DBCSR matrix containing A
!> \param dist_a the dbcsr distribution of matrix A
!> \param blk_size_a the row/col block sizes of matrix A
!> \param donor_state ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE build_matrix_a(matrix_a, dist_a, blk_size_a, donor_state, qs_env)

      TYPE(dbcsr_type), POINTER                       :: matrix_a
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(dbcsr_distribution_type), POINTER          :: dist_a
      INTEGER, DIMENSION(:), POINTER                  :: blk_size_a
      
      CHARACTER(len=*), PARAMETER :: routineN = 'build_matrix_a', routineP = moduleN//":"//routineN
      
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
!        will be later added to the matrix B (kernel), which is itslef non-symmetric

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

   END SUBROUTINE build_matrix_a

! **************************************************************************************************
!> \brief Creates the so-called matrix G needed for the LHS of the generalized eigenvalue problem
!>        omega*|-G  0| |C+ | = M |C+|, G_{pi,qj} = S_pq*delta_ij
!>              | 0  G| |C- |     |C-|
!> \param matrix_g dbcsr matrix containing G
!> \param dist_g dbcsr distribution of G
!> \param blk_size_g block sizes of G
!> \param donor_state ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE build_matrix_g(matrix_g, dist_g, blk_size_g, donor_state, qs_env)

      TYPE(dbcsr_type), POINTER                          :: matrix_g
      TYPE(dbcsr_distribution_type), POINTER             :: dist_g
      INTEGER, DIMENSION(:), POINTER                     :: blk_size_g
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_matrix_g', routineP = moduleN//":"//routineN
      
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

   END SUBROUTINE build_matrix_g

! **************************************************************************************************
!> \brief Builds the so-called matrix B, which contains the kernel parts (Coulomb, xc, exact 
!>        exchange) with appropriate weights. Same block structure and distribution as A and G
!> \param matrix_b pointer to a DBCSR matrix containing B
!> \param dist_b the dbcsr distribution of b
!> \param blk_size_b the row/col block sizes of matrix B
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE build_matrix_b(matrix_b, dist_b, blk_size_b, donor_state, xas_tdp_env, qs_env)

      TYPE(dbcsr_type), POINTER                       :: matrix_b
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(dbcsr_distribution_type), POINTER          :: dist_b
      INTEGER, DIMENSION(:), POINTER                  :: blk_size_b
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_matrix_b', routineP = moduleN//":"//routineN

      INTEGER                                         :: handle
      TYPE(dbcsr_type), POINTER                       :: coulomb_kernel

      CALL timeset(routineN, handle)

!  Creating the kernel matrix B
      CALL dbcsr_create(matrix=matrix_b, name="MATRIX B", matrix_type='N', dist=dist_b, &
                        row_blk_size=blk_size_b, col_blk_size=blk_size_b)

!  Compute the Coulomb part of the kernel
      ALLOCATE(coulomb_kernel) 
      CALL build_kernel_coulomb(coulomb_kernel, dist_b, blk_size_b, donor_state, xas_tdp_env, qs_env)
      CALL dbcsr_add(matrix_b, coulomb_kernel, 0.0_dp, 1.0_dp)
      CALL dbcsr_release(coulomb_kernel)
      DEALLOCATE(coulomb_kernel)

!  Finalize 
      CALL dbcsr_finalize(matrix_b)

      CALL timestop(handle)

   END SUBROUTINE build_matrix_b

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

! **************************************************************************************************
!> \brief Computes the three-centers Coulomb integrals needed for RI in kernel calculations and 
!>        returns it in appropriate array of full matrices. Each entry of the output array corresponds
!>        to the sgfs centered on atom a of: sum_i c_i*(ai|p), {p} is the RI basis centered on the 
!>        excited atom.
!> \param int_array array containing the integrals in full matrix form
!)> \param nonzero_int array containing the indices of the atoms, for which there is a nonzero integral 
!> \param donor_state ...
!> \param domo_index the index of donor_state's MO to consider (1 for 1s or 2s, 1,2,3 for 2p)
!> \param qs_env ...
!> \note  Computes integrals of the form: sum_i (ai|p)*c_i, where c_i are the ground state coefficient
!>        for the donor MO. TODO: there is probably a high potential for parallelization since we have
!>        many loop that do their own stuff...
! **************************************************************************************************
SUBROUTINE compute_ri_coulomb3_int(int_array, nonzero_int, donor_state, domo_index, qs_env)

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: int_array
      TYPE(donor_state_type), POINTER                 :: donor_state
      INTEGER, DIMENSION(:), ALLOCATABLE              :: nonzero_int
      INTEGER, INTENT(IN)                             :: domo_index
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_coulomb3_int', &
                                     routinP=moduleN//':'//routineN

      INTEGER                                         :: handle, nsetp, nseti, natom, &
                                                         nao, nseta, iat, aset, iset, pset, &
                                                         ncoa, ncoi, ncop, igp, sgfa, sgfp, &
                                                         sgfi, aat, iao, isgf, akind, ikind, n_int,&
                                                         asgf, aao, i, nsgf, nkind
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_set_list_a, basis_set_list_i
      TYPE(gto_basis_set_type), POINTER               :: p_basis, i_basis, a_basis
      INTEGER, DIMENSION(:), POINTER                  :: lp_max, lp_min, npgfp, li_max,li_min,&
                                                         npgfi, la_max, la_min, npgfa, nsgfa, &
                                                         nsgfi, nsgfp
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgfa, first_sgfi, first_sgfp
      REAL(dp), DIMENSION(:,:), POINTER               :: zetp, rpgfp, zeti, rpgfi, zeta, rpgfa,&
                                                         sphi_p, sphi_a, sphi_i, coeffi
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos 
      TYPE(cp_fm_type), POINTER                       :: mo_coeff
      REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE       :: v
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE         :: int_aip, int_aip_set, int_aip_contr, &
                                                         int_aip_atom
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: vaip, int_coeff_contr
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: f, gcc
      REAL(dp), DIMENSION(3)                          :: ra, ri, rp, rai, rap, rip
      REAL(dp)                                        :: rai2, rap2, rip2, dap, dip, rada, radi, &
                                                         radp
      INTEGER, DIMENSION(:), ALLOCATABLE              :: first_sgf_set, last_sgf_set, wnonzero_int
      TYPE(neighbor_list_set_p_type), DIMENSION(:), & 
                                          POINTER     :: sab_orb
      TYPE(neighbor_list_iterator_p_type), &                                                         
                            DIMENSION(:), POINTER     :: nl_iterator   
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(cell_type), POINTER                        :: cell

      CALL timeset(routineN, handle)

      NULLIFY(qs_kind_set, p_basis, i_basis, a_basis, particle_set, mos, mo_coeff)
      NULLIFY(lp_max, lp_min, npgfp, zetp, rpgfp, sphi_p, first_sgfp)
      NULLIFY(li_max, li_min, npgfi, zeti, rpgfi, sphi_i, first_sgfi)
      NULLIFY(la_max, la_min, npgfa, zeta, rpgfa, sphi_a, first_sgfa)
      NULLIFY(coeffi, int_array, sab_orb, nl_iterator, fm_struct, para_env, blacs_env, cell)
      NULLIFY(basis_set_list_a, basis_set_list_i)

!  Get all we need from the qs_env
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, mos=mos, & 
                      natom=natom, sab_orb=sab_orb, para_env=para_env, blacs_env=blacs_env, &
                      cell=cell)
      ALLOCATE(first_sgf_set(natom))
      ALLOCATE(last_sgf_set(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf_set, last_sgf=last_sgf_set)
      
!  Get the RI basis set for the excited atom
      CALL get_qs_kind(qs_kind_set(donor_state%kind_index), basis_set=p_basis, basis_type="RI_AUX")

      CALL get_gto_basis_set(p_basis, lmax=lp_max, npgf=npgfp, zet=zetp, pgf_radius=rpgfp, &
                             lmin=lp_min, nset=nsetp, nsgf_set=nsgfp, sphi=sphi_p, &
                             first_sgf=first_sgfp, kind_radius=radp)

!  Create lists of basis for atoms a and i
      nkind = SIZE(qs_kind_set)
      ALLOCATE(basis_set_list_a(nkind))
      ALLOCATE(basis_Set_list_i(nkind))
      CALL basis_set_list_setup(basis_set_list_a, "ORB", qs_kind_set)
      CALL basis_set_list_setup(basis_set_list_i, "ORB", qs_kind_set)
      

!  We also need the coefficients of the donor_state MO
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao)
      ALLOCATE(coeffi(nao,1))
      CALL cp_fm_get_submatrix(fm=mo_coeff,target_m=coeffi,start_row=1, &
         start_col=donor_state%mo_indices(domo_index), n_rows=nao, n_cols=1 )

!  Get the atoms aat and iat that will give a non-zero integral (ai|p). They all need to be neighbor
!  or they are screened away in the coulomb3 subroutine. 
      n_int = 0 
      ALLOCATE(wnonzero_int(natom))
      wnonzero_int = 0
      rp = pbc(particle_set(donor_state%at_index)%r, cell) 
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0) 
         CALL get_iterator_info(nl_iterator, iatom=aat, jatom=iat, ikind=akind, jkind=ikind, r=rai)

!        Check if atoms a, i and p are all neighbors by checking sum of radii versus distance
         ra = pbc(particle_set(aat)%r, cell)
         ri = rai + ra
         rap = pbc(ra, rp, cell)
         rip = rap - rai
         dap = SQRT(rap(1)**2+rap(2)**2+rap(3)**2)
         dip = SQRT(rip(1)**2+rip(2)**2+rai(3)**2)
         a_basis => basis_set_list_a(akind)%gto_basis_set
         CALL get_gto_basis_set(a_basis, kind_radius=rada)
         i_basis => basis_set_list_i(ikind)%gto_basis_set
         CALL get_gto_basis_set(i_basis, kind_radius=radi)

!        Add both aat and iat as nonzero integrals if not there already
         IF (rada+radp >= dap .AND. radi+radi >= dip) THEN
            IF (.NOT. MINVAL(ABS(wnonzero_int - aat))==0 ) THEN
               n_int = n_int+1
               wnonzero_int(n_int) = aat
            END IF
            IF (.NOT. MINVAL(ABS(wnonzero_int - iat))==0) THEN
               n_int = n_int+1
               wnonzero_int(n_int) = iat
            END IF
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ALLOCATE(nonzero_int(n_int)) 
      nonzero_int(:) = wnonzero_int(1:n_int)
      DEALLOCATE(wnonzero_int)
      ALLOCATE(int_array(n_int))

!  Create fm as elements of the int_array
!  TODO: many very small non-square fm. Should we just store arrays and then parallelized with open MP ?
!  TODO: causes same issue in parallel as with multiple donor states in a row: it looks like the 
!        cp_fm_struct get mixed up because processors do not go at same speed
      DO i = 1,n_int
         iat = nonzero_int(i)
         nsgf = last_sgf_set(iat) - first_sgf_set(iat) + 1
         
         CALL cp_fm_struct_create(fm_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=nsgf, ncol_global=SUM(nsgfp)) 

         CALL cp_fm_create(int_array(i)%matrix, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END DO

!  Loop over all neighbouring atoms, compute (ai|p) only if all three atoms are neighbors.
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=akind, jkind=ikind, iatom=aat, jatom=iat, &
                                r=rai)

!        the iterator does not do double counting => compute (ai|p) and (ia|p) in the same loop
!        that is: compute one integral and contract with MO coeff once on i and once on a

!        Retrieve the basis sets for both atoms
         a_basis => basis_set_list_a(akind)%gto_basis_set
         i_basis => basis_set_list_i(ikind)%gto_basis_Set


!        Retrieve the quantum numbers for both basis sets
         CALL get_gto_basis_set(a_basis, lmax=la_max, npgf=npgfa, zet=zeta, pgf_radius=rpgfa, &
                                lmin=la_min, nset=nseta, nsgf_set=nsgfa, sphi=sphi_a, & 
                                first_sgf=first_sgfa, kind_radius=rada) 
         CALL get_gto_basis_set(i_basis, lmax=li_max, npgf=npgfi, zet=zeti, pgf_radius=rpgfi, &
                                lmin=li_min, nset=nseti, nsgf_set=nsgfi, sphi=sphi_i, &
                                first_sgf=first_sgfi, kind_radius=radi)
            
!        Get the different postition/distances
         ra = pbc(particle_set(aat)%r, cell)                                                                 
         ri = rai + ra
         rp = pbc(particle_set(donor_state%at_index)%r, cell) 
         rap = pbc(ra, rp, cell)
         rip = rap - rai
         rai2 = rai(1)**2+rai(2)**2+rai(3)**2                                
         rap2 = rap(1)**2+rap(2)**2+rap(3)**2                                
         rip2 = rip(1)**2+rip(2)**2+rai(3)**2
         dap = SQRT(rap2)
         dip = SQRT(rip2)

!        Check for nonzero integral 
         IF (rada+radp < dap .OR. radi+radp < dip) CYCLE

!        Allocate an array where the integrals for the current atoms is stored (in contracted sgf form)
         ALLOCATE(int_aip_atom(SUM(nsgfa), SUM(nsgfi), SUM(nsgfp)))
         int_aip_atom = 0.0_dp 
         
!        Loop over the sets of the three basis
         DO aset = 1,nseta
            DO iset = 1,nseti
               DO pset = 1,nsetp

!                 Some initialization
                  ncoa = npgfa(aset)*ncoset(la_max(aset))
                  ncoi = npgfi(iset)*ncoset(li_max(iset))
                  ncop = npgfp(pset)*ncoset(lp_max(pset))
                  sgfa = first_sgfa(1, aset)
                  sgfi = first_sgfi(1, iset)
                  sgfp = first_sgfp(1, pset) 

!                 Subroutine coulomb3 only takes one Gaussian at a time for the c of (ab|c)       
!                 Need to loop on the primitive gaussians of the RI basis                         
!                 Still, save the integral in part of a big matrix that spans all GTOs of the RI set
                  ALLOCATE(int_aip_set(ncoa, ncoi, ncop))
                  int_aip_set = 0.0_dp

                  DO igp = 1,npgfp(pset)    

!                    Allocate work array memory
                     ALLOCATE(f(0:la_max(aset)+li_max(iset)+lp_max(pset)+2))    
                     ALLOCATE(v(ncoa, ncoi, ncoset(lp_max(pset)), la_max(aset)+li_max(iset)+lp_max(pset)+1))
                     ALLOCATE(vaip(ncoa, ncoi)) !this one does not play a role here
                     ALLOCATE(gcc(ncoset(lp_max(pset)))) !this one does not play a role here
                     ALLOCATE(int_aip(ncoa, ncoi, ncoset(lp_max(pset))))

                     f = 0.0_dp
                     v = 0.0_dp
                     vaip = 0.0_dp
                     gcc = 0.0_dp
                     int_aip = 0.0_dp

!                    Compute the integral 
                     CALL coulomb3(la_max(aset), npgfa(aset), zeta(:,aset), rpgfa(:,aset), la_min(aset), &
                                   li_max(iset), npgfi(iset), zeti(:,iset), rpgfi(:,iset), li_min(iset), &
                                   lp_max(pset), zetp(igp,pset), rpgfp(igp,pset), lp_min(pset), &
                                   gcc, rai, rai2, rap, rap2, rip2, vaip, int_aip, v, f)

!                    Add this as a block in an array that includes all contributions of the GTOs
!                    of the current set of the RI basis
                     int_aip_set(:,:,(igp-1)*ncoset(lp_max(pset))+1:igp*ncoset(lp_max(pset))) = int_aip

!                    Deallocate arrays, might have different size in next loop
                     DEALLOCATE(f, v, vaip, gcc, int_aip)

                  END DO ! igp

!                 Perform contraction to get spherical GTOs
                  ALLOCATE(int_aip_contr(nsgfa(aset), nsgfi(iset), nsgfp(pset)))
                  CALL abc_contract(int_aip_contr, int_aip_set, sphi_a(:,sgfa:), sphi_i(:,sgfi:),&
                                    sphi_p(:,sgfp:), ncoa, ncoi, ncop, nsgfa(aset), &
                                    nsgfi(iset), nsgfp(pset))

!                 Store the set specific combination of integrals
                  int_aip_atom(sgfa:sgfa+nsgfa(aset)-1,sgfi:sgfi+nsgfi(iset)-1, &
                               sgfp:sgfp+nsgfp(pset)-1) = int_aip_contr

!                 Deallocate before next loop
                  DEALLOCATE(int_aip_set, int_aip_contr)

               END DO ! pset
            END DO !iset
         END DO ! aset

!        Need to contract with the donor MO coefficient. 3 possible cases: aat == iat and do
!        contraction once on either a or i. aat /= iat and need to do contraction on each
!        TODO: can we make those contraction product more efficient using LAPACK gemm instead ?
         IF (aat == iat) THEN

            ALLOCATE(int_coeff_contr(SUM(nsgfa),SUM(nsgfp)))
            int_coeff_contr = 0.0_dp
            DO isgf = 1,SUM(nsgfi)
               iao = first_sgf_set(iat) + isgf -1
               int_coeff_contr(:,:) = int_coeff_contr(:,:) + int_aip_atom(:,isgf,:)*coeffi(iao,1)
            END DO
!           Add to existing full matrix of int_array, at the right spot
            i = MINLOC(ABS(nonzero_int-aat), 1)
            CALL cp_fm_set_submatrix(int_array(i)%matrix, int_coeff_contr, alpha=1.0_dp, beta=1.0_dp)
            DEALLOCATE(int_coeff_contr)

         ELSE ! iat /= aat

!           First contract with coeff c_i as entry aat in int_array
            ALLOCATE(int_coeff_contr(SUM(nsgfa),SUM(nsgfp)))
            int_coeff_contr = 0.0_dp
            DO isgf = 1,SUM(nsgfi)
               iao = first_sgf_set(iat) + isgf -1
               int_coeff_contr(:,:) = int_coeff_contr(:,:) + int_aip_atom(:,isgf,:)*coeffi(iao,1)
            END DO
!           Add to existing full matrix in int_array at the right spot
            i = MINLOC(ABS(nonzero_int-aat), 1)
            CALL cp_fm_set_submatrix(int_array(i)%matrix, int_coeff_contr, alpha=1.0_dp, beta=1.0_dp)
            DEALLOCATE(int_coeff_contr)

!           Then contract with c_a as entry iat in int_array
            ALLOCATE(int_coeff_contr(SUM(nsgfi),SUM(nsgfp)))
            int_coeff_contr = 0.0_dp
            DO asgf = 1,SUM(nsgfa)
               aao = first_sgf_set(aat) + asgf -1
               int_coeff_contr(:,:) = int_coeff_contr(:,:) + int_aip_atom(asgf,:,:)*coeffi(aao,1)
            END DO
!           Add to existing fm in int_array at the right spot
            i = MINLOC(ABS(nonzero_int-iat), 1)
            CALL cp_fm_set_submatrix(int_array(i)%matrix, int_coeff_contr, alpha=1.0_dp, beta=1.0_dp)
            DEALLOCATE(int_coeff_contr)
         END IF ! iat==aat
            
         DEALLOCATE(int_aip_atom)

      END DO ! nl_iterator
      CALL neighbor_list_iterator_release(nl_iterator)

!  Clean-up
      DEALLOCATE(coeffi, basis_set_list_a, basis_set_list_i)

      CALL timestop(handle)

   END SUBROUTINE compute_ri_coulomb3_int

! **************************************************************************************************
!> \brief Creates a matrix Wab = sum_k Q_ak*(kI|P)*(P|Q)^-1(Q|Jb) based on the RI 3 centers integrals
!>        from compute_ri_coulomb3 and the projector Q on the unperturbed, unoccupied ground-state. 
!>        Returns a dbcsr matrix the size of the KS matrix. 
!> \param ker_coulomb_block the dbcsr block containing the RI Coulomb integrals for I, J
!> \param int_arri the ri integral array for donor MO I
!> \param int_arrj the ri integral array for donor MO J
!> \param nonzero_inti the indices of non-zero integrals for donor MO I
!> \param nonzero_intj the indices of non-zero integrals for donor MO J
!> \param donor_state the core donor state from which we excite an electron 
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \notes for excitation from a s state, I=J=1, for a p state I,J = 1,2,3
! **************************************************************************************************
   SUBROUTINE build_ri_coulomb_block(ker_coulomb_blk, int_arri, int_arrj, nonzero_inti, &
                                     nonzero_intj, donor_state, xas_tdp_env,qs_env)

      TYPE(dbcsr_type), POINTER                       :: ker_coulomb_blk
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: int_arri, int_arrj
      INTEGER, DIMENSION(:), ALLOCATABLE              :: nonzero_inti, nonzero_intj
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_ri_coulomb_block', &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, iat, i, jat, j, ikind, mynode, &
                                                         nrowi, ncoli, ncolj, node_holds_blk
      TYPE(cp_fm_type), POINTER                       :: pq_inv, worki, workj
      TYPE(cp_fm_struct_type), POINTER                :: fm_structi, fm_structj
      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(cp_para_env_type), POINTER                 :: para_env 
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: block_data
      TYPE(dbcsr_type), POINTER                       :: q_projector, work_coulomb

      CALL timeset(routineN, handle)
      
      NULLIFY(pq_inv, dbcsr_dist, matrix_s, worki, workj, fm_structi, fm_structj, para_env, &
              blacs_env, q_projector, blk_size)

!  Get what we need from the qs_env
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist, matrix_s=matrix_s, blacs_env=blacs_env, &
                      para_env=para_env)

!  Create the dbcsr_matrix that will hold everything
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=blk_size)
      CALL dbcsr_create(matrix=ker_coulomb_blk, name="RI COULOMB BLK", dist=dbcsr_dist, matrix_type="N", &
                        row_blk_size=blk_size, col_blk_size=blk_size)
      ALLOCATE(work_coulomb)
      CALL dbcsr_create(matrix=work_coulomb, name="WORK COULOMB", dist=dbcsr_dist, matrix_type="N", &
                        row_blk_size=blk_size, col_blk_size=blk_size)

!  Retrieve the (P|Q)^-1 matrix from xas_tdp_env, and the Q projector 
      ikind = MINLOC(ABS(xas_tdp_env%ex_kind_indices-donor_state%kind_index),1)
      pq_inv => xas_tdp_env%ri_inv_mats(ikind)%matrix
      q_projector => xas_tdp_env%q_projector
      
!  Loop over the nonzero coulomb integrals and fill work_coulomb
      CALL dbcsr_distribution_get(dbcsr_dist, mynode=mynode)
      DO i = 1,SIZE(nonzero_inti)
         DO j = 1,SIZE(nonzero_intj)

            iat = nonzero_inti(i)
            jat = nonzero_intj(j)

!           Create work structures
            CALL cp_fm_get_info(int_arri(i)%matrix, nrow_global=nrowi, ncol_global=ncoli)
            CALL cp_fm_get_info(int_arrj(j)%matrix, nrow_global=ncolj)
            CALL cp_fm_struct_create(fm_structi, context=blacs_env, para_env=para_env, &
                                     nrow_global=nrowi, ncol_global=ncoli)
            CALL cp_fm_struct_create(fm_structj, context=blacs_env, para_env=para_env, &
                                     nrow_global=ncoli, ncol_global=ncolj)
            CALL cp_fm_create(worki, fm_structi)
            CALL cp_fm_create(workj, fm_structj)

!           Compute the matrix product (aI|P)*(P|Q)^-1
            CALL cp_fm_gemm('N', 'N', nrowi, ncoli, ncoli, 1.0_dp, int_arri(i)%matrix, pq_inv, &
                             0.0_dp, worki)

!           Compute the product [(aI|P)*(P|Q)^-1]*(Q|Jb)
            CALL cp_fm_gemm('N', 'T', nrowi, ncolj, ncoli, 1.0_dp, worki, int_arrj(j)%matrix, &
                            0.0_dp, workj)

!           Get the whole data block
            ALLOCATE(block_data(nrowi, ncolj))
            CALL cp_fm_get_submatrix(workj, block_data)

!           Put it in the dbcsr_matrix
            CALL dbcsr_get_stored_coordinates(work_coulomb, iat, jat, node_holds_blk)
            IF (node_holds_blk == mynode) CALL dbcsr_put_block(work_coulomb, iat, jat, block_data) 

!           Clean-up
            DEALLOCATE(block_data)
            CALL cp_fm_struct_release(fm_structi)
            CALL cp_fm_struct_release(fm_structj)
            CALL cp_fm_release(worki)
            CALL cp_fm_release(workj)

         END DO ! j
      END DO !i
      CALL dbcsr_finalize(work_coulomb)

!  Multiply by the projector Q
      CALL dbcsr_multiply('N','N',1.0_dp, q_projector, work_coulomb, 0.0_dp, ker_coulomb_blk)
      CALL dbcsr_finalize(ker_coulomb_blk)
       
!  Clean-up
      CALL dbcsr_release(work_coulomb)
      DEALLOCATE(work_coulomb)

      CALL timestop(handle)

   END SUBROUTINE build_ri_coulomb_block

! **************************************************************************************************
!> \brief Computes the Coulomb part of the kernel matrix using RI
!> \param coulomb_kernel the dbcsr matrix containing the Coulomb part of the kernel, based on RI
!> \param dist the dbcsr ditribution of the coulomb_kernel matrix
!> \param blk_size the block sizes of the coulomb_kernel matrix
!> \param donor_state the donor_state from which we excite a core electron
!> \param xas_tdp_env ...
!> \param qs_env ...
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

      INTEGER                                         :: handle, i, j, ndo_mo, k, iblk, jblk, blk, &
                                                         nblk_row
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER       :: int_arri, int_arrj
      INTEGER, DIMENSION(:), ALLOCATABLE              :: nonzero_inti, nonzero_intj
      TYPE(dbcsr_type), POINTER                       :: coulomb_blk
      TYPE(dbcsr_iterator_type)                       :: iter
      REAL(dp), DIMENSION(:), POINTER                 :: work_block
      LOGICAL                                         :: found_block

      CALL timeset(routineN, handle)

      NULLIFY(int_arri, int_arrj, coulomb_blk, work_block)

!  Initilaization
      ndo_mo = SIZE(donor_state%mo_indices,1)
      nblk_row = SIZE(blk_size)/ndo_mo

!  Create the dbcsr matrix for the coulomb kernel
      CPASSERT(ASSOCIATED(coulomb_kernel))
      CALL dbcsr_create(matrix=coulomb_kernel, name="COULOMB KERNEL", matrix_type="N", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)

!  Loop over donor MOs  TODO: compute int_array only once for each I and then loop over that) 
      DO i = 1,ndo_mo
!        Compute the RI int_array for donor MO I
         CALL compute_ri_coulomb3_int(int_arri, nonzero_inti, donor_state, i, qs_env)

         DO j = 1,ndo_mo
!           Compute the RI int_array for donor MO J
            CALL compute_ri_coulomb3_int(int_arrj, nonzero_intj, donor_state, j, qs_env)

!           Build the IJ block
            ALLOCATE(coulomb_blk)
            CALL build_ri_coulomb_block(coulomb_blk, int_arri, int_arrj, nonzero_inti, nonzero_intj, &
                                        donor_state, xas_tdp_env, qs_env)

!           And put it in the coulomb kernel at the right spot
            CALL dbcsr_iterator_start(iter, coulomb_blk)
            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

               found_block = .FALSE.
               CALL dbcsr_get_block_p(coulomb_blk, iblk, jblk, work_block, found_block)

               IF (found_block) THEN

                  CALL dbcsr_put_block(coulomb_kernel, (i-1)*nblk_row+iblk, (j-1)*nblk_row+jblk, &
                                       work_block)
                  
               END IF
               NULLIFY(work_block)

            END DO ! dbcsr_iterator
            CALL dbcsr_iterator_stop(iter)

!           Clean-up for J
            CALL dbcsr_release(coulomb_blk)
            DO k = 1,SIZE(nonzero_intj)
               CALL cp_fm_release(int_arrj(k)%matrix)
            END DO ! k
            DEALLOCATE(int_arrj, nonzero_intj, coulomb_blk)
         END DO !j

!        Clean-up for I
         DO k = 1,SIZE(nonzero_inti)
            CALL cp_fm_release(int_arri(k)%matrix)
         END DO ! k
         DEALLOCATE(int_arri, nonzero_inti)
      END DO !i

!  Finalize
      CALL dbcsr_finalize(coulomb_kernel)

      CALL timestop(handle)

   END SUBROUTINE build_kernel_coulomb

END MODULE xas_tdp_utils

