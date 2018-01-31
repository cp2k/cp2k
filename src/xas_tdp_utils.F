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
   USE basis_set_types,                 ONLY: gto_basis_set_type, get_gto_basis_set
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose, &
                                              cp_dbcsr_cholesky_invert, &
                                              cp_dbcsr_cholesky_restore
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm, cp_dbcsr_sm_fm_multiply, &
                                              copy_fm_to_dbcsr
   USE cp_cfm_types,                    ONLY: cp_fm_to_cfm, cp_cfm_to_fm, cp_cfm_type, &
                                              cp_cfm_create, cp_cfm_release
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert, cp_fm_scale_and_add, cp_fm_gemm
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create, cp_fm_release, cp_fm_type, &
                                              cp_fm_to_fm_submat, cp_fm_get_submatrix, &
                                              cp_fm_get_info, cp_fm_set_submatrix, cp_fm_set_all
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
                                              dbcsr_complete_redistribute
   USE input_constants,                 ONLY: xas_dip_len, xas_dip_vel
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: ncoset
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: qs_environment_type, get_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind, qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set, mo_set_p_type
   USE qs_operators_ao,                 ONLY: p_xyz_ao, rRc_xyz_ao
   USE xas_tdp_types,                   ONLY: donor_state_type, xas_tdp_control_type

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_utils'

   PUBLIC :: setup_xas_tdp_full_prob, solve_xas_tdp_full_prob, compute_oscillator_strength

CONTAINS

! **************************************************************************************************
!> \brief Builds the matrices that define the XAS TDDFPT generalized eigenvalue problem to be solved
!>        for excitation energies omega. The problem has the form omega*S*C = M*C, where C contains
!>        the reponse orbitals coefficients. This is a full build, so that conventional diag
!>        methods can be applied.
!> \param matrix_tdp pointer to the dbcsr matrix M (right hand side)
!> \parma matrix_lhs pointer to the dbcsr matrix S (left hand side)
!> \param donor_state donor_state type for which the problem is restricted
!> \param qs_env ...
!> \param xas_tdp_control ...
!> \note At this point in time (01.2018), the Kernel is not taken into account, i.e. matrix M  has
!>       the form | A  0 | and another subroutine builds A. Matrix S has the form | -G  0 |, where
!>                | 0  A |                                                        |  0  G |
!>       G is also built elsewhere. Both tdp and lhs matrices have the same block sizes and dbcsr
!>       distribution, which are refered to as supmat_dist and supmat_blk_size (for supermatrix)
!>       If the Tamm-Dancoff approximation is assumed, then C^+ = 0, M = A and S = G
!>       If TDA not required then need to solve the non-hermitian problem omega*C = S^-1*M*C, where
!>       S^-1*M*C  = | -G^-1*A   0    |
!>                   |   0    G^-1*A  |  
! **************************************************************************************************
   SUBROUTINE setup_xas_tdp_full_prob(matrix_tdp, matrix_lhs, donor_state, qs_env, xas_tdp_control)

      TYPE(dbcsr_type), POINTER                       :: matrix_tdp
      TYPE(dbcsr_type), POINTER                       :: matrix_lhs
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      
      CHARACTER(len=*), PARAMETER :: routineN = 'setup_xas_tdp_full_prob', &
                                     routineP = moduleN//":"//routineN
      
      INTEGER                                         :: handle, group , iblk, jblk, nblk_sub, &
                                                         mynode, node_holds_blk
      TYPE(dbcsr_type), POINTER                       :: matrix_a, matrix_g, work_matrix
      TYPE(dbcsr_distribution_type), POINTER          :: submat_dist, supmat_dist
      INTEGER, DIMENSION(:), POINTER                  :: submat_blk_size, supmat_blk_size, &
                                                         col_dist_sub, row_dist_sub, &
                                                         col_dist_sup, row_dist_sup
      INTEGER, DIMENSION(:,:), POINTER                :: pgrid
      LOGICAL                                         :: found_block
      REAL(dp), DIMENSION(:), POINTER                 :: work_block
      TYPE(cp_fm_type), POINTER                       :: fm_matrix_g, fm_matrix_inv_g
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env

      CALL timeset(routineN, handle)

      NULLIFY(matrix_a, submat_dist, submat_blk_size, supmat_dist, supmat_blk_size, col_dist_sub)
      NULLIFY(row_dist_sub, col_dist_sup, row_dist_sup, work_block, matrix_g, pgrid, work_matrix)
      NULLIFY(fm_matrix_g, fm_matrix_inv_g, fm_struct, para_env, blacs_env)

      CPASSERT(ASSOCIATED(donor_state))
      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(xas_tdp_control))

!  Get the submatrices distribution and block sizes
      CALL compute_submat_dist_and_blk_size(submat_dist, submat_blk_size, donor_state, qs_env)
      nblk_sub = SIZE(submat_blk_size)

!  Compute the submatrices A and G
      ALLOCATE (matrix_a)
      ALLOCATE (matrix_g)
      CALL build_matrix_a(matrix_a, submat_dist, submat_blk_size, donor_state, qs_env)
      CALL build_matrix_g(matrix_g, submat_dist, submat_blk_size, donor_state, qs_env)

!  Test on Tamm Dancoff approximation

      IF (xas_tdp_control%tamm_dancoff) THEN
!        If TDA is assumed, then the matrix_tdp and the matrix_lhs are basically A and G
         CALL dbcsr_copy(matrix_tdp, matrix_a, name="MATRIX TDP")
         CALL dbcsr_copy(matrix_lhs, matrix_g, name="MATRIX LHS")

      ELSE IF (.NOT. xas_tdp_control%tamm_dancoff) THEN

!     The supermatrices matrix_tdp and matrix_lhs are built with the submatrices as a 2x2 block matrices
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
!        block sizes and distribution. In this case, no need of matrix_lhs
         ALLOCATE(supmat_blk_size(2*nblk_sub))
         supmat_blk_size(1:nblk_sub) = submat_blk_size
         supmat_blk_size(nblk_sub+1:2*nblk_sub) = submat_blk_size 
         CALL dbcsr_create(matrix=matrix_tdp, name="XAS TDP MATRIX",matrix_type="N", dist=supmat_dist,&
                           row_blk_size=supmat_blk_size, col_blk_size=supmat_blk_size)
         DEALLOCATE(supmat_blk_size)

!     Filling the supermatrix matrix_tdp with submatrix with +- G^-1*A

!        Creating a work matrix where the product G^-1*A will be stored. Same dist and blk_size as
!        A and G but non-symmetric
         ALLOCATE(work_matrix) 
         CALL dbcsr_create(matrix=work_matrix, name="WORK MATRIX", matrix_type="N",dist=submat_dist,&
                           row_blk_size=submat_blk_size, col_blk_size=submat_blk_size)

!        Inversion of matrix G, need to go through full matrix
         CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=SUM(submat_blk_size),ncol_global=SUM(submat_blk_size))
         CALL cp_fm_create(fm_matrix_g, fm_struct)
         CALL cp_fm_create(fm_matrix_inv_g, fm_struct)

         CALL copy_dbcsr_to_fm(matrix_g, fm_matrix_g)
         CALL cp_fm_invert(fm_matrix_g, fm_matrix_inv_g)

!        Store the inverse of G in G itself (inverse is symmetric too)
         CALL copy_fm_to_dbcsr(fm_matrix_inv_g, matrix_g)

         CALL cp_fm_release(fm_matrix_g)
         CALL cp_fm_release(fm_matrix_inv_g)
         CALL cp_fm_struct_release(fm_struct)

!        Compute the product G^-1*A
         CALL dbcsr_multiply('N','N',1.0_dp, matrix_g, matrix_a, 0.0_dp, work_matrix)
         CALL dbcsr_finalize(work_matrix)

!        For each diagonal block of the super matrix, loop over blocks of G^-1*A
!        TODO: replace for loops by dbcsr_iterator
         CALL dbcsr_distribution_get(submat_dist, mynode=mynode)
         DO iblk = 1,nblk_sub
            DO jblk = 1,nblk_sub
           
!              Checking that block is stored on current node
               CALL dbcsr_get_stored_coordinates(work_matrix, iblk, jblk, node_holds_blk)

               IF(node_holds_blk == mynode) THEN
                  CALL dbcsr_get_block_p(work_matrix, iblk, jblk, work_block, found_block)

                  IF (found_block) THEN
                     CALL dbcsr_put_block(matrix_tdp, iblk, jblk, -work_block)
                     CALL dbcsr_put_block(matrix_tdp, nblk_sub+iblk, nblk_sub+jblk, work_block)
                  END IF
            
                  NULLIFY(work_block)
               END IF

            END DO
         END DO

!        Some clean-up
         CALL dbcsr_distribution_release(supmat_dist)
         CALL dbcsr_release(work_matrix)
         DEALLOCATE(supmat_dist, work_matrix)

      END IF ! Tamm_Dancoff

!  Finalize
      CALL dbcsr_finalize(matrix_tdp)
      IF(xas_tdp_control%tamm_dancoff) CALL dbcsr_finalize(matrix_lhs)

!  General Clean-up
      CALL dbcsr_distribution_release(submat_dist)
      CALL dbcsr_release(matrix_a)
      CALL dbcsr_release(matrix_g)
      DEALLOCATE (submat_dist, submat_blk_size, matrix_a, matrix_g)

      CALL timestop(handle)

   END SUBROUTINE setup_xas_tdp_full_prob

! **************************************************************************************************
!> \brief Solves the XAS TDP generalized eigenvalue problem omega*matrix_lhs*C = matrix_tdp*C
!>        using standard full diagonalization methods
!> \param matrix_tdp ...
!> \param matrix_lhs ...
!> \param donor_state ...
!> \param qs_env ...
!> \param xas_tdp_control ...
!> \note The computed eigenvalues and eigenvectors are stored in the donor_state
! **************************************************************************************************
   SUBROUTINE solve_xas_tdp_full_prob(matrix_tdp, matrix_lhs, donor_state, qs_env, xas_tdp_control)  

      TYPE(dbcsr_type), POINTER                          :: matrix_tdp, matrix_lhs
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = "solve_xas_tdp_full_prob", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                            :: handle, nrow
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: amatrix, bmatrix, work
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (para_env, blacs_env)
      NULLIFY (fm_struct, amatrix, bmatrix, work)

      CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)

      IF (xas_tdp_control%tamm_dancoff) THEN
!        In this case matrix_tdp is symmetric and so is matrix_lhs. The latter is also positive
!        definite, which allows us to use Cholesky decomposition to solve the geeig
!        We convert to full matrices to make use of the existing machinery

!        Initialization
         CALL dbcsr_get_info(matrix_tdp, nfullrows_total=nrow)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nrow, para_env=para_env, &
                                  ncol_global=nrow)

         CALL cp_fm_create(amatrix, fm_struct)
         CALL cp_fm_create(bmatrix, fm_struct)
         CALL cp_fm_create(work, fm_struct)

!        Copy geeig in full matrix format
         CALL copy_dbcsr_to_fm(matrix_tdp, amatrix)
         CALL copy_dbcsr_to_fm(matrix_lhs, bmatrix)

!        Allcoating space for solutions
         ALLOCATE (donor_state%lr_evals(nrow))
         CALL cp_fm_create(donor_state%lr_coeffs, fm_struct)

!        Solve the geeig
         CALL cp_fm_geeig(amatrix, bmatrix, donor_state%lr_coeffs, donor_state%lr_evals, work)

!        Clean-up
         CALL cp_fm_release(amatrix)
         CALL cp_fm_release(bmatrix)
         CALL cp_fm_release(work)
         CALL cp_fm_struct_release(fm_struct)

      ELSE IF (.NOT. xas_tdp_control%tamm_dancoff) THEN
!        In this case matrix_lhs is irrelevant and matrix_tdp is non-symmetric. We do a direct
!        diagonalization of the latter using inefficient, non-parallelized methods

!        Initialization
         CALL dbcsr_get_info(matrix_tdp, nfullrows_total=nrow)
         CALL cp_fm_struct_create(fm_struct,context=blacs_env, nrow_global=nrow, para_env=para_env,&
                                  ncol_global=nrow)

!        Need to go to full matrices to diagonalize
         CALL cp_fm_create(amatrix, fm_struct)
         CALL copy_dbcsr_to_fm(matrix_tdp, amatrix)
         
!        Allocating space for solutions
         ALLOCATE(donor_state%lr_evals(nrow))
         CALL cp_fm_create(donor_state%lr_coeffs, fm_struct)

!        Diagonalize the matrix
#if defined(__SCALAPACK) 
         CALL pdiag_nonsymm_mat(amatrix, donor_state%lr_evals, donor_state%lr_coeffs)
#else
         CALL diag_nonsymm_mat(amatrix, donor_state%lr_evals, donor_state%lr_coeffs)
#endif
         
!        Clean-up
         CALL cp_fm_release(amatrix)
         CALL cp_fm_struct_release(fm_struct)

      END IF

!  Making sure that the response orbitals are normalized with respect to matrix_s
!  TODO: RETEST WHEN KERNEL IS AVAILABLE
      CALL normalize_lr_orbitals(donor_state, qs_env, xas_tdp_control)

      CALL timestop(handle)

   END SUBROUTINE solve_xas_tdp_full_prob

! **************************************************************************************************
!> \brief Normalizes the response orbitals so that <phi^+|phi^+>=1 and <phi^-|phi^->=1
!> \param donor_state ...
!> \param qs_env ... 
!> \param xas_tdp_control ...
! **************************************************************************************************
   SUBROUTINE normalize_lr_orbitals(donor_state, qs_env, xas_tdp_control)

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control               

      CHARACTER(len=*), PARAMETER :: routineN = "normalize_lr_orbitals", &
                                      routineP = moduleN//":"//routineN     

      INTEGER                                         :: handle, nao, iex, j, ipm, imo, npm, &
                                                         s_firstrow, ndo_mo
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos 
      TYPE(cp_fm_type), POINTER                       :: mo_coeff
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      REAL(dp), DIMENSION(:,:), POINTER               :: coeff_array, work_array, full_col
      TYPE(cp_fm_struct_type), POINTER                :: col_vect_struct
      TYPE(cp_fm_type), POINTER                       :: c_vect, work_vect, lr_coeffs
      REAL(dp)                                        :: norm
      LOGICAL                                         :: found_nonzero
      
      CALL timeset(routineN, handle) 

      NULLIFY(mos, mo_coeff, para_env, blacs_env, matrix_s, coeff_array, work_array)
      NULLIFY(col_vect_struct, c_vect, work_vect, lr_coeffs, full_col)

!  Notes: we could simply take one coefficient vector C^+_j from lr_coeffs, compute the norm
!  C^+_j*S*C^+_j and divide the whole column of lr_coeffs by SQRT(norm). The problem is that
!  there is a chance that a linear response orbital has zero coefficient, which could be a
!  catastrophy when dividing. So instead we span all vectors C^+-_j until we find a non-zero one

!  Initialization
      CALL get_qs_env(qs_env, mos=mos, para_env=para_env, blacs_env=blacs_env, matrix_s=matrix_s)
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao)

      CALL cp_fm_struct_create(col_vect_struct, para_env=para_env, context=blacs_env, &              
                               nrow_global=nao, ncol_global=1)                                      
      CALL cp_fm_create(c_vect, col_vect_struct)                                                     
      CALL cp_fm_create(work_vect, col_vect_struct)                                                  
      ALLOCATE(coeff_array(nao,1))                                                                      
      ALLOCATE(work_array(nao,1))  
      lr_coeffs => donor_state%lr_coeffs

      npm = 2
      IF (xas_tdp_control%tamm_dancoff) npm=1
      ndo_mo = SIZE(donor_state%mo_indices)
      ALLOCATE(full_col(npm*ndo_mo*nao,1))

!  Loop over the response orbital for different excitation energies
      DO iex = 1,SIZE(donor_state%lr_evals)
         norm = 0.0_dp
         found_nonzero = .FALSE.
!        Loop over the coeffs for different donor MOs
         DO imo = 1,ndo_mo
            IF (.NOT. found_nonzero) THEN
!              Loop over C^+/C^-
               DO ipm = 1,npm
                  
                  IF (.NOT. found_nonzero) THEN
         
                     s_firstrow = (ipm-1)*ndo_mo*nao + (imo-1)*nao + 1
!                    Retrieving the C^+ (C^- if TDA) coefficients corresponding to the current donor MO
                     CALL cp_fm_to_fm_submat(msource=lr_coeffs, mtarget=c_vect, nrow=nao, ncol=1, &
                                             s_firstrow=s_firstrow, s_firstcol=iex, &
                                             t_firstrow=1, t_firstcol=1)
                     CALL cp_fm_get_submatrix(c_vect, coeff_array)

!                    Computing matrix_s*c_vect
                     CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, c_vect, work_vect, ncol=1)
                     CALL cp_fm_get_submatrix(work_vect, work_array)

                     DO j = 1,nao
                        norm = norm + coeff_array(j,1)*work_array(j,1)
                     END DO
                     If (norm .NE. 0.0_dp) found_nonzero = .TRUE.

                  END IF ! found_nonzero
               END DO ! ipm
            END IF ! found_nonzero
         END DO ! imo

!        Scaling the lr_coeff for this energy
         full_col = 0.0_dp
         CALL cp_fm_get_submatrix(fm=lr_coeffs, target_m=full_col, start_row=1, start_col = iex, &
                                  n_rows=nao*ndo_mo*npm, n_cols=1, transpose=.FALSE.)
         full_col  = 1.0_dp/SQRT(norm) * full_col
         CALL cp_fm_set_submatrix(fm=lr_coeffs, new_values=full_col, start_row=1, start_col = iex,&
                                  n_rows=nao*ndo_mo*npm, n_cols=1, alpha=1.0_dp, beta=0.0_dp, &
                                  transpose=.FALSE.)

      END DO ! iex

!  Clean-up
      CALL cp_fm_struct_release(col_vect_struct)
      CALL cp_fm_release(c_vect)
      CALL cp_fm_release(work_vect)
      DEALLOCATE(coeff_array, work_array, full_col)

      CALL timestop(handle)
      
   END SUBROUTINE normalize_lr_orbitals

! **************************************************************************************************
!> \brief Computes the oscillator strength based on the dipole moment (velocity or position rep
!>        for all available excitation energies of a donor_state. Stores the result in donor_state
!> \param donor_state the donor state which is excited
!> \param qs_env ...
!> \param xas_tdp_control ...
!> \note The oscillator strength is a scalar: osc_str = 2/3*omega*(dipole_r)^2 = -2/(3*omega)*(dipole_v)^2
!>       But because there is a risk of very tiny omega, multiply both cases by omega to avoid division
! **************************************************************************************************
   SUBROUTINE compute_oscillator_strength(donor_state, qs_env, xas_tdp_control)

      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      
      CHARACTER(len=*), PARAMETER :: routineN = "compute_oscillator_strength", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, nosc, nao, i, ndo_mo,j, k, jdo_mo,l
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos
      TYPE(cp_fm_type), POINTER                       :: mo_coeff
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: dipmat
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      REAL(dp), DIMENSION(3)                          :: rc, xyz_moments
      REAL(dp), DIMENSION(:), POINTER                 :: osc_str
      REAL(dp), DIMENSION(:,:), POINTER               :: gs_coeff, work_array
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(cp_fm_struct_type), POINTER                :: col_vect_struct
      TYPE(cp_fm_type), POINTER                       :: c_vect, work_vect
      REAL(dp)                                        :: contrib
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: int_array

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

      CALL dbcsr_allocate_matrix_set(dipmat,3)
      DO i = 1,3
         ALLOCATE(dipmat(i)%matrix)
         CALL dbcsr_copy(dipmat(i)%matrix, matrix_s(1)%matrix, name="XAS TDP dipole matrix")
         CALL dbcsr_set(dipmat(i)%matrix, 0.0_dp)
      END DO

      IF (xas_tdp_control%dipole_form == xas_dip_len) THEN
!        Take the origin of r on the excited atom
         rc = particle_set(donor_state%at_index)%r
         CALL rRc_xyz_ao(dipmat, qs_env, rc, order=1)
!        Add the dipole form specific prefactor 2/3*omega**2
         osc_str = 2.0_dp/3.0_dp*donor_state%lr_evals**2

      ELSE IF (xas_tdp_control%dipole_form == xas_dip_vel) THEN
!        Simply compute the moment in the basis set
         CALL p_xyz_ao(dipmat, qs_env)
!        Add the dipole form specific prefactor 2/3, no minus sign as p_xyz is real here
         osc_str = 2.0_dp/3.0_dp
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
      CALL dbcsr_deallocate_matrix_set(dipmat)

      CALL timestop(handle)

!  Testing three-center integrals
      NULLIFY(int_array)
!      CALL compute_ri_coulomb3_int(int_array, donor_state, qs_env) 
!      CALL dbcsr_deallocate_matrix_set(int_array)

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
                                                         j, iblk, jblk, mynode, node_holds_blk
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks, matrix_s
      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      INTEGER, DIMENSION(:), POINTER                  :: row_blk_size 
      REAL(dp), DIMENSION(:,:), POINTER               :: energy_evals     
      REAL(dp), DIMENSION(:), POINTER                 :: work_block, work_block_ks, work_block_s 
      LOGICAL                                         :: found_block, found_block_ks 
      TYPE(dbcsr_type)                                :: desym_matrix_s, desym_matrix_ks, work_matrix

      CALL timeset(routineN, handle)

      NULLIFY (matrix_ks, dbcsr_dist, row_blk_size, work_block, energy_evals, matrix_s)
      NULLIFY (work_block_ks, work_block_s)

!  Note: All matrices are based on the basic CP2K dbcsr_distribution. Larger matrices simply replicate
!        the distribution. For non-symmetric matrices, the outcome is straight forward, not so much
!        otherwise (pretty much impossible to predict which processors holds which block)
!  Strategy: Use a non-symmetric work matrix and desymmetrize matrices S and KS, hence the blocks
!            are kept on the same processors for all matrices. After filling the upper block 
!            triangle part of the work matrix, redistribute it into a symmetric matrix_a.

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
   
!  Loop over the blocks of desymmetrized matrices S and KS and put them in the work matrix
      CALL dbcsr_distribution_get(dist_a, mynode=mynode)
      DO iblk = 1,nblk_row
         DO jblk = 1,nblk_row

!           Check if block is on current processor
            CALL dbcsr_get_stored_coordinates(desym_matrix_s, iblk, jblk, node_holds_blk)
            IF (node_holds_blk == mynode) THEN

               found_block = .FALSE.                                                           
               found_block_ks = .FALSE.                                                        
               ALLOCATE(work_block(row_blk_size(iblk)*row_blk_size(jblk)))

!              Retrieve S_pq
               CALL dbcsr_get_block_p(desym_matrix_s, iblk, jblk, work_block_s, found_block)

!              Retrieve the upper triangular part of F_pq (all we need since only on diag of mat A)
               IF (jblk .GE. iblk) CALL dbcsr_get_block_p(desym_matrix_ks, iblk, jblk, &
                                                          work_block_ks, found_block_ks)

!              Now loop over the upper triangular part of the work matrix.
               DO i = 1,ndo_mo
                  DO j = i,ndo_mo

!                    Take -epsilon_ij*Spq
                     IF (found_block .AND. (energy_evals(i,j) .NE. 0.0_dp)) THEN
                        work_block = -(energy_evals(i,j))*work_block_s

                     ELSE IF (energy_evals(i,j) ==  0.0_dp) THEN
!                       If multiplying block by zero, then equivalent to no block at all
                        found_block = .FALSE.
                     END IF

!                    If not on the block digonal of work mat, can simply put the block now
                     IF ((i .NE. j) .AND. found_block) THEN
                        CALL dbcsr_put_block(work_matrix,(i-1)*nblk_row+iblk, (j-1)*nblk_row+jblk, & 
                                             work_block)
                     END IF ! i != j

!                    If on the block diagonal, take F_pq - epsilon_ii*S_pq and only consider upper part
                     IF ((i==j) .AND. (jblk .GE. iblk)) THEN

                        IF (found_block .AND. found_block_ks) THEN
                           work_block = work_block + work_block_ks
                           CALL dbcsr_put_block(work_matrix,(i-1)*nblk_row+iblk, &
                                                (j-1)*nblk_row+jblk, work_block)

                        ELSE IF ((.NOT. found_block) .AND. found_block_ks) THEN
                           work_block = work_block_ks
                           CALL dbcsr_put_block(work_matrix,(i-1)*nblk_row+iblk, &
                                                (j-1)*nblk_row+jblk, work_block)
                        END IF
                        
                     END IF ! i==j

                  END DO ! j
               END DO ! i

               DEALLOCATE(work_block)
               NULLIFY(work_block_s, work_block_ks)

            END IF ! node_holds_blk==mynode
         END DO ! jblk
      END DO ! iblk

!  Create symmetric matrix A
      CALL dbcsr_create(matrix=matrix_a, name="MATRIX A", matrix_type="S", dist=dist_a, &
                        row_blk_size=blk_size_a, col_blk_size=blk_size_a)    

!  Redistribute the work matrix into A so that it is stored as symmetric
      CALL dbcsr_finalize(work_matrix)
      CALL dbcsr_complete_redistribute(work_matrix, matrix_a)

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
                                                         mynode, node_holds_blk
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      INTEGER, DIMENSION(:), POINTER                  :: row_blk_size
      REAL(dp), DIMENSION(:), POINTER                 :: work_block
      LOGICAL                                         :: found_block

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

!  Fill the matrix G. First loop over the superblocks of size nao x nao (only on the diagonal), then
!  their respective blocks which have the standard cp2k structure. They also have the standard dbcsr
!  distribution. Since G is block diagonal, no trap like for build_matrix_a
      CALL dbcsr_distribution_get(dist_g, mynode=mynode)
      DO i = 1,ndo_mo
!        loop over blocks of matrix_s (upper triangular part, symmetric anyway)
!        TO DO: replace for loops by dbcsr_iterator
         DO iblk = 1,nblk_row
            DO jblk = iblk,nblk_row

!              Check that the block is stored on current processor.
               CALL dbcsr_get_stored_coordinates(matrix_s(1)%matrix, iblk, jblk, node_holds_blk)

               IF (node_holds_blk == mynode) THEN
                  CALL dbcsr_get_block_p(matrix_s(1)%matrix, iblk, jblk, work_block, found_block)

!                 if non zero-block, copy its content in G
                  IF (found_block) THEN
                     CALL dbcsr_put_block(matrix_g, (i-1)*nblk_row+iblk, (i-1)*nblk_row+jblk, &
                                          work_block)
                  END IF

                  NULLIFY(work_block)

               END IF 

            END DO
         END DO
      END DO

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
!   SUBROUTINE build_matrix_b(matrix_b, dist_b, blk_size_b, donor_state, xas_tdp_env, qs_env)
!
!      TYPE(dbcsr_type), POINTER                       :: matrix_b
!      TYPE(donor_state_type), POINTER                 :: donor_state
!      TYPE(qs_environment_type), POINTER              :: qs_env
!      TYPE(dbcsr_distribution_type), POINTER          :: dist_b
!      INTEGER, DIMENSION(:), POINTER                  :: blk_size_b
!      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env

!      CHARACTER(len=*), PARAMETER :: routineN = 'build_matrix_b', routineP = moduleN//":"//routineN

!      INTEGER                                         :: handle

!      CALL timeset(routineN, handle)

!      CALL timestop(handle)

!   END SUBROUTINE build_matrix_b

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
!>        returns it in appropriate array of full matrices
!> \param int_array array containing the integrals in full matrix form
!> \param donor_state ...
!> \param qs_env ...
!> \note  Computes integrals of the form (ai|p)*c_i, where c_i is the ground state coefficient for 
!>        the donor MO
! **************************************************************************************************
   SUBROUTINE compute_ri_coulomb3_int(int_array, donor_state, qs_env)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: int_array
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_ri_coulomb3_int', &
                                     routinP=moduleN//':'//routineN

      INTEGER                                         :: handle, nsetp, nseti, natom, &
                                                         nao, nseta, iat, aset, iset, pset, &
                                                         ncoa, ncoi, ncop, igp, sgfa, sgfp, &
                                                         sgfi, group, aat, iao, isgf
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(gto_basis_set_type), POINTER               :: p_basis, i_basis, a_basis
      INTEGER, DIMENSION(:), POINTER                  :: lp_max, lp_min, npgfp, li_max,li_min,&
                                                         npgfi, la_max, la_min, npgfa, nsgfa, &
                                                         nsgfi, nsgfp, row_dist, col_dist, &
                                                         row_blk_size, col_blk_size
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgfa, first_sgfi, first_sgfp, pgrid
      REAL(dp), DIMENSION(:,:), POINTER               :: zetp, rpgfp, zeti, rpgfi, zeta, rpgfa,&
                                                         sphi_p, sphi_a, sphi_i, coeffi
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mos 
      TYPE(cp_fm_type), POINTER                       :: mo_coeff
      REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE       :: v
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE         :: int_aip, int_aip_set, int_aip_contr, &
                                                         int_aip_iat
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: vaip
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: f, gcc
      REAL(dp), DIMENSION(3)                          :: ra, rb, rc
      REAL(dp)                                        :: rab2, rac2, rbc2
      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      TYPE(dbcsr_distribution_type)                   :: int_aip_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      INTEGER, DIMENSION(:), ALLOCATABLE              :: first_sgf_set

      CALL timeset(routineN, handle)

!  We need to screen the which of the atom pairs iat/aat actually have overlaping basis, so that we
!  don't compute integrals that we now are zero. This way, can also go to full matrices and only
!  store a short number, along with a system to which matrix is what

!      IF (SIZE(donor_state%energy_evals,1) .NE. 1) &
!         CPABORT("RI Coulomb integrals not implemented yet for cases with multiple donor MOs")

      NULLIFY(qs_kind_set, p_basis, i_basis, a_basis, particle_set, mos, mo_coeff)
      NULLIFY(lp_max, lp_min, npgfp, zetp, rpgfp, sphi_p, first_sgfp)
      NULLIFY(li_max, li_min, npgfi, zeti, rpgfi, sphi_i, first_sgfi)
      NULLIFY(la_max, la_min, npgfa, zeta, rpgfa, sphi_a, first_sgfa)
      NULLIFY(dbcsr_dist, pgrid, row_blk_size, col_blk_size, coeffi)

!  Get all we need from the qs_env
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, mos=mos, & 
                      natom=natom, dbcsr_dist=dbcsr_dist, matrix_s=matrix_s)
      ALLOCATE(first_sgf_set(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf_set)
      
!  Get the RI basis set for the excited atom
      CALL get_qs_kind(qs_kind_set(donor_state%kind_index), basis_set=p_basis, basis_type="RI_AUX")

      CALL get_gto_basis_set(p_basis, lmax=lp_max, npgf=npgfp, zet=zetp, pgf_radius=rpgfp, &
                             lmin=lp_min, nset=nsetp, nsgf_set=nsgfp, sphi=sphi_p, &
                             first_sgf=first_sgfp)

!  We also need the coefficients of the donor_state MOs
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao)
      ALLOCATE(coeffi(nao,1))
      CALL cp_fm_get_submatrix(fm=mo_coeff,target_m=coeffi,start_row=1, &
                               start_col=donor_state%mo_indices(1), n_rows=nao, n_cols=1 )


!  Build the array of integral matrices. It has an entry for each atomic orbital, which contains
!  the matrix (ai|p), where a is the atomic orbital (each entry of array has a different a), i spans
!  all the atomic orbitals (row index of matrix) and p spans the elements of the RI basis 
!  (column index). 
      
      CALL dbcsr_allocate_matrix_set(int_array, nao)

!  Each of the element of the array will be a super-sparse matrix. For now, don't care about the 
!  performance, we distribute all blocks on the same processor... We put one block for each atom
      
      CALL dbcsr_distribution_get(dbcsr_dist, group=group, pgrid=pgrid)      
      ALLOCATE(row_dist(natom))
      ALLOCATE(col_dist(1))
!     TODO: do not use reuse array, it detroys pgrid on release
      CALL dbcsr_distribution_new(int_aip_dist, group=group, pgrid=pgrid, row_dist=row_dist, &
                                  col_dist=col_dist, reuse_arrays=.TRUE.)

      ALLOCATE(col_blk_size(1))
      col_blk_size(1) = SUM(nsgfp)
      CALL dbcsr_get_info(matrix_s(1)%matrix, row_blk_size=row_blk_size)

      DO iao = 1,nao
         ALLOCATE(int_array(iao)%matrix)
         CALL dbcsr_create(matrix=int_array(iao)%matrix, name="int_array(i)%matrix", &
                           matrix_type="N", dist=int_aip_dist, col_blk_size=col_blk_size, &
                           row_blk_size=row_blk_size)
      END DO

      DEALLOCATE(col_blk_size)
      CALL dbcsr_distribution_release(int_aip_dist)      

!  Loop over all atoms and compute the three-center integrals
      DO aat = 1,natom
!        Retrieve the orbital basis of the current atom (for the a of (ai|p))
         CALL get_qs_kind(qs_kind_set(particle_set(aat)%atomic_kind%kind_number), &
                          basis_set=a_basis, basis_type="ORB")
         CALL get_gto_basis_set(a_basis, lmax=la_max, npgf=npgfa, zet=zeta, pgf_radius=rpgfa, &
                               lmin=la_min, nset=nseta, nsgf_set=nsgfa, sphi=sphi_a, &
                               first_sgf=first_sgfa)

         DO iat = 1,natom

!           Retrieve the basis of the second atom (for the i of (ai|p))
            CALL get_qs_kind(qs_kind_set(particle_set(iat)%atomic_kind%kind_number), &
                                         basis_set=i_basis, basis_type="ORB")

            CALL get_gto_basis_set(i_basis, lmax=li_max, npgf=npgfi, zet=zeti, pgf_radius=rpgfi, &
                                   lmin=li_min, nset=nseti, nsgf_set=nsgfi, sphi=sphi_i, &
                                   first_sgf=first_sgfi)

!           Create an aray where all integrals for the atom iat are strored (contracted sgf form)
            ALLOCATE(int_aip_iat(SUM(nsgfa), SUM(nsgfi), SUM(nsgfp)))
            int_aip_iat = 0.0_dp


!           Get set the different position/distances 
            ra = particle_set(aat)%r
            rb = particle_set(iat)%r
            rc = particle_set(donor_state%at_index)%r 
            rab2 = (rb(1)-ra(1))**2+(rb(2)-ra(2))**2+(rb(3)-ra(3))**2                                   
            rac2 = (rc(1)-ra(1))**2+(rc(2)-ra(2))**2+(rc(3)-ra(3))**2                                   
            rbc2 = (rc(1)-rb(1))**2+(rc(2)-rb(2))**2+(rc(3)-rb(3))**2

!           Loop over the sets of the three basis
            DO aset = 1,nseta
               DO iset = 1,nseti
                  DO pset = 1,nsetp 

!                    Subroutine coulomb3 only takes one Gaussian at a time for the c of (ab|c)
!                    Need to loop on the primitive gaussians of the RI basis
!                    Still, save the integral in part of a big matrix that spans all GTOs of the RI set
                     ncoa = npgfa(aset)*ncoset(la_max(aset))
                     ncoi = npgfi(iset)*ncoset(li_max(iset))
                     ncop = npgfp(pset)*ncoset(lp_max(pset))
                     ALLOCATE(int_aip_set(ncoa, ncoi, ncop))
                     int_aip_set = 0.0_dp

                     sgfa = first_sgfa(1, aset)
                     sgfi = first_sgfi(1, iset)
                     sgfp = first_sgfp(1, pset)

                     DO igp = 1,npgfp(pset)

!                       Allocate memory for the different work arrays
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

!                       Compute the integral 
                        CALL coulomb3(la_max(aset), npgfa(aset), zeta(:,aset), rpgfa(:,aset), la_min(aset), &
                                      li_max(iset), npgfi(iset), zeti(:,iset), rpgfi(:,iset), li_min(iset), &
                                      lp_max(pset), zetp(igp,pset), rpgfp(igp,pset), lp_min(pset), &
                                      gcc, rb-ra, rab2, rc-ra, rac2, rbc2, vaip, int_aip, v, f)

!                       Add this as a block in an array that includes all contributions of the GTOs
!                       of the current set of the RI basis
                        int_aip_set(:,:,(igp-1)*ncoset(lp_max(pset))+1:igp*ncoset(lp_max(pset)))&
                                   = int_aip

!                       Deallocate arrays, might have different size in next loop
                        DEALLOCATE(f, v, vaip, gcc, int_aip)
                     
                     END DO ! igp

!                    Perform contraction to get spherical GTOs
                     ALLOCATE(int_aip_contr(nsgfa(aset), nsgfi(iset), nsgfp(pset)))
                     CALL abc_contract(int_aip_contr, int_aip_set, sphi_a(:,sgfa:), sphi_i(:,sgfi:),&
                                       sphi_p(:,sgfp:), ncoa, ncoi, ncop, nsgfa(aset), &
                                       nsgfi(iset), nsgfp(pset))

!                    Store the specific set combination of integrals
                     int_aip_iat(sgfa:sgfa+nsgfa(aset)-1,sgfi:sgfi+nsgfi(iset)-1, &
                                 sgfp:sgfp+nsgfp(pset)-1) = int_aip_contr

!                    Deallocate before next loop
                     DEALLOCATE(int_aip_set, int_aip_contr)

                  END DO ! pset
               END DO ! iset
            END DO !aset

!           We compute the integral for the donor MO (i in (ai|p)), so  we need to consider the
!           MO coefficients too
            
            DO isgf = 1,SUM(nsgfi)
               iao = first_sgf_set(iat)+isgf-1
               int_aip_iat(:,isgf,:) = int_aip_iat(:,isgf,:)*coeffi(iao,1)
            END DO

!           TODO: make sure we compute the right thing... double and even triple check
!           If so then we have our blocks

!           Deallocate before next loop
            DEALLOCATE(int_aip_iat)

         END DO !iat
      END DO ! aat

!  Clean-up
      DEALLOCATE(coeffi)

      CALL timestop(handle)

   END SUBROUTINE compute_ri_coulomb3_int

END MODULE xas_tdp_utils

