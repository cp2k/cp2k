!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for X-ray absorption spectroscopy using TDDFPT
!> \author AB (01.2018)
! **************************************************************************************************

MODULE xas_tdp_utils

   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_allocate_matrix_set, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix_set, &
        dbcsr_distribution_get, dbcsr_distribution_new, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, dbcsr_p_type, &
        dbcsr_put_block, dbcsr_release, dbcsr_set, dbcsr_type
   USE input_constants,                 ONLY: xas_dip_len,&
                                              xas_dip_vel
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_operators_ao,                 ONLY: p_xyz_ao,&
                                              rRc_xyz_ao
   USE xas_tdp_types,                   ONLY: donor_state_type,&
                                              xas_tdp_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_utils'

   PUBLIC :: setup_xas_tdp_geeig_prob, solve_xas_tdp_geeig_prob, compute_oscillator_strength

CONTAINS

! **************************************************************************************************
!> \brief Builds the matrices that define the XAS TDDFPT generalized eigenvalue problem to be solved
!>        for excitation energies omega. The problem has the form omega*S*C = M*C, where C contains
!>        the reponse orbitals coefficients.
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
! **************************************************************************************************
   SUBROUTINE setup_xas_tdp_geeig_prob(matrix_tdp, matrix_lhs, donor_state, qs_env, xas_tdp_control)

      TYPE(dbcsr_type), POINTER                          :: matrix_tdp, matrix_lhs
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_xas_tdp_geeig_prob', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: group, handle, iblk, jblk, nblk_sub
      INTEGER, DIMENSION(:), POINTER                     :: col_dist_sub, col_dist_sup, &
                                                            row_dist_sub, row_dist_sup, &
                                                            submat_blk_size, supmat_blk_size
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      LOGICAL                                            :: found_block
      REAL(dp), DIMENSION(:), POINTER                    :: work_block
      TYPE(dbcsr_distribution_type), POINTER             :: submat_dist, supmat_dist
      TYPE(dbcsr_type), POINTER                          :: matrix_a, matrix_g

      CALL timeset(routineN, handle)

      NULLIFY (matrix_a, submat_dist, submat_blk_size, supmat_dist, supmat_blk_size, col_dist_sub)
      NULLIFY (row_dist_sub, col_dist_sup, row_dist_sup, work_block, matrix_g, pgrid)

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

!        Creating the symmetric supermatrices matrix_tdp and matrix_lhs with the 2x2 replication of the
!        submatrices block sizes and distribution
         ALLOCATE (supmat_blk_size(2*nblk_sub))
         supmat_blk_size(1:nblk_sub) = submat_blk_size
         supmat_blk_size(nblk_sub+1:2*nblk_sub) = submat_blk_size
         CALL dbcsr_create(matrix=matrix_tdp, name="XAS TDP MATRIX", matrix_type="S", dist=supmat_dist, &
                           row_blk_size=supmat_blk_size, col_blk_size=supmat_blk_size)
         CALL dbcsr_create(matrix=matrix_lhs, name="LHS MATRIX", matrix_type="S", dist=supmat_dist, &
                           row_blk_size=supmat_blk_size, col_blk_size=supmat_blk_size)
         DEALLOCATE (supmat_blk_size)

!     Filling the supermatrix matrix_tdp with submatrix A

!        For each diagonal block of the super matrix, loop over super-diagonal blocks of submat A
         DO iblk = 1, nblk_sub
            DO jblk = iblk, nblk_sub

               CALL dbcsr_get_block_p(matrix_a, iblk, jblk, work_block, found_block)

               IF (found_block) THEN
                  CALL dbcsr_put_block(matrix_tdp, iblk, jblk, work_block)
                  CALL dbcsr_put_block(matrix_tdp, nblk_sub+iblk, nblk_sub+jblk, work_block)
               END IF

               NULLIFY (work_block)

            END DO
         END DO

!     Filling the supermatrix matrix_lhs with the submatrix G

!        Loop over super-diagonal blocks of submatrix G
         DO iblk = 1, nblk_sub
            DO jblk = iblk, nblk_sub

               CALL dbcsr_get_block_p(matrix_g, iblk, jblk, work_block, found_block)

               IF (found_block) THEN
                  CALL dbcsr_put_block(matrix_lhs, iblk, jblk, -work_block)
                  CALL dbcsr_put_block(matrix_lhs, nblk_sub+iblk, nblk_sub+jblk, work_block)
               END IF

               NULLIFY (work_block)

            END DO
         END DO

!        Some clean-up
         CALL dbcsr_distribution_release(supmat_dist)
         DEALLOCATE (supmat_dist)

      END IF ! Tamm_Dancoff

!  Finalize
      CALL dbcsr_finalize(matrix_tdp)
      CALL dbcsr_finalize(matrix_lhs)

!  General Clean-up
      CALL dbcsr_distribution_release(submat_dist)
      CALL dbcsr_release(matrix_a)
      CALL dbcsr_release(matrix_g)
      DEALLOCATE (submat_dist, submat_blk_size, matrix_a, matrix_g)

      CALL timestop(handle)

   END SUBROUTINE setup_xas_tdp_geeig_prob

! **************************************************************************************************
!> \brief Solves the XAS TDP generalized eigenvalue problem omega*matrix_lhs*C = matrix_tdp*C
!> \param matrix_tdp ...
!> \param matrix_lhs ...
!> \param donor_state ...
!> \param qs_env ...
!> \param xas_tdp_control ...
!> \note The computed eigenvalues and eigenvectors are stored in the donor_state
! **************************************************************************************************
   SUBROUTINE solve_xas_tdp_geeig_prob(matrix_tdp, matrix_lhs, donor_state, qs_env, xas_tdp_control)

      TYPE(dbcsr_type), POINTER                          :: matrix_tdp, matrix_lhs
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = 'solve_xas_tdp_geeig_prob', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nrow
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: amatrix, bmatrix, work
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (para_env, blacs_env)
      NULLIFY (fm_struct, amatrix, bmatrix, work)

      CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)

!  Full TDA diagonalization

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

      END IF

      CALL timestop(handle)

   END SUBROUTINE solve_xas_tdp_geeig_prob

! **************************************************************************************************
!> \brief Computes the oscillator strength based on the dipole moment (velocity or position rep
!>        for all available excitation energies of a donor_state. Stores the result in donor_state
!> \param donor_state the donor state which is excited
!> \param qs_env ...
!> \param xas_tdp_control ...
!> \note The oscillator strength is a scalar: osc_str = 2/3*omega*(dipole_r)^2
!>                                                    = -2/(3*omega)*(dipole_v)^2
! **************************************************************************************************
   SUBROUTINE compute_oscillator_strength(donor_state, qs_env, xas_tdp_control)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_oscillator_strength', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, jdo_mo, k, l, nao, ndo_mo, &
                                                            nosc
      REAL(dp)                                           :: contrib
      REAL(dp), DIMENSION(3)                             :: rc, xyz_moments
      REAL(dp), DIMENSION(:), POINTER                    :: osc_str
      REAL(dp), DIMENSION(:, :), POINTER                 :: gs_coeff, work_array
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: col_vect_struct
      TYPE(cp_fm_type), POINTER                          :: c_minus, mo_coeff, work_vect
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipmat, matrix_s
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      NULLIFY (mos, mo_coeff, para_env, blacs_env, dipmat, particle_set, osc_str, matrix_s)
      NULLIFY (col_vect_struct, c_minus, work_vect, gs_coeff, work_array)

!  Test for TDA because of different data structure (should be unified in the future)
      IF (xas_tdp_control%tamm_dancoff) THEN
!        General initialization
         ndo_mo = SIZE(donor_state%mo_indices)
         nosc = SIZE(donor_state%lr_evals)
         ALLOCATE (donor_state%osc_str(nosc))
         osc_str => donor_state%osc_str
         CALL get_qs_env(qs_env, mos=mos, para_env=para_env, blacs_env=blacs_env, &
                         particle_set=particle_set, matrix_s=matrix_s)
         CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao)

         CALL dbcsr_allocate_matrix_set(dipmat, 3)
         DO i = 1, 3
            ALLOCATE (dipmat(i)%matrix)
            CALL dbcsr_copy(dipmat(i)%matrix, matrix_s(1)%matrix, name="XAS TDP dipole matrix")
            CALL dbcsr_set(dipmat(i)%matrix, 0.0_dp)
         END DO

         IF (xas_tdp_control%dipole_form == xas_dip_len) THEN
!           Take the origin of r on the excited atom
            rc = particle_set(donor_state%at_index)%r
            CALL rRc_xyz_ao(dipmat, qs_env, rc, order=1)
!           Add the dipole form specific prefactor 2/3*omega
            osc_str = 2.0_dp/3.0_dp*donor_state%lr_evals

         ELSE IF (xas_tdp_control%dipole_form == xas_dip_vel) THEN
!           Simply compute the moment in the basis set
            CALL p_xyz_ao(dipmat, qs_env)
!           Add the dipole form specific prefactor 2/(3*omega), no minus sign as p_xyz is real here
            osc_str = 2.0_dp/3.0_dp*1.0_dp/donor_state%lr_evals
         END IF

!        Some more initialization before heavy computing
         CALL cp_fm_struct_create(col_vect_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=nao, ncol_global=1)
         CALL cp_fm_create(c_minus, col_vect_struct)
         CALL cp_fm_create(work_vect, col_vect_struct)
         ALLOCATE (gs_coeff(nao, 1))
         ALLOCATE (work_array(nao, 1))

!        the dipole moment is <mu> = sum_j C^0_j*M*C^-_j, where j labels donor MOs, M is the dipmat,
!        C^0_i is the ground state coefficient row vector and and C^- the column vector of lr coeffs
         DO i = 1, nosc ! loop on excitation energies
            xyz_moments = 0.0_dp
            DO j = 1, ndo_mo ! loop on donor MOs

               jdo_mo = donor_state%mo_indices(j)

!              Retrieve column vector C^-
               CALL cp_fm_to_fm_submat(msource=donor_state%lr_coeffs, mtarget=c_minus, nrow=nao, &
                                       ncol=1, s_firstrow=(j-1)*nao+1, s_firstcol=i, t_firstrow=1, &
                                       t_firstcol=1)

               gs_coeff = 0.0_dp
!              Retrieve array of ground state MO coeff
               CALL cp_fm_get_submatrix(fm=mo_coeff, target_m=gs_coeff, start_row=1, &
                                        start_col=jdo_mo, n_rows=nao, n_cols=1, transpose=.FALSE.)

               DO k = 1, 3 ! loop on cartesian coordinates
!                 Computing M*C^-_j
                  CALL cp_dbcsr_sm_fm_multiply(dipmat(k)%matrix, c_minus, work_vect, ncol=1)
                  work_array = 0.0_dp
                  CALL cp_fm_get_submatrix(work_vect, work_array)

!                 Computing C^0_j*M*C^-_j
                  DO l = 1, nao
                     xyz_moments(k) = xyz_moments(k)+gs_coeff(l, 1)*work_array(l, 1)
                  END DO ! l
               END DO ! k
            END DO ! j

!           Compute the oscilaltor strength for this specific energy
            contrib = 0.0_dp
            DO k = 1, 3
               contrib = contrib+osc_str(i)*xyz_moments(k)**2
            END DO ! k
            osc_str(i) = contrib
         END DO ! i

!        Clean-up
         DEALLOCATE (gs_coeff, work_array)
         CALL cp_fm_release(c_minus)
         CALL cp_fm_release(work_vect)
         CALL cp_fm_struct_release(col_vect_struct)
         CALL dbcsr_deallocate_matrix_set(dipmat)

      ELSE IF (.NOT. xas_tdp_control%tamm_dancoff) THEN

      END IF

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

      TYPE(dbcsr_type), POINTER                          :: matrix_a
      TYPE(dbcsr_distribution_type), POINTER             :: dist_a
      INTEGER, DIMENSION(:), POINTER                     :: blk_size_a
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_matrix_a', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iblk, j, jblk, nblk_row, &
                                                            ndo_mo
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size
      LOGICAL                                            :: found_block, found_block_ks
      REAL(dp), DIMENSION(:), POINTER                    :: work_block, work_block_ks, work_block_s
      REAL(dp), DIMENSION(:, :), POINTER                 :: energy_evals
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s

      CALL timeset(routineN, handle)

      NULLIFY (matrix_ks, dbcsr_dist, row_blk_size, work_block, energy_evals, matrix_s)
      NULLIFY (work_block_ks, work_block_s)

!  Initialization
      ndo_mo = SIZE(donor_state%mo_indices, 1)
      CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks, matrix_s=matrix_s, dbcsr_dist=dbcsr_dist)
      CALL dbcsr_get_info(matrix_ks(1)%matrix, row_blk_size=row_blk_size)
      nblk_row = SIZE(row_blk_size)
      energy_evals => donor_state%energy_evals

!  Creating the matrix A with the right size and distribution (it is symmetric)
      CALL dbcsr_create(matrix=matrix_a, name="MATRIX A", matrix_type="S", dist=dist_a, &
                        row_blk_size=blk_size_a, col_blk_size=blk_size_a)

!  Fill the matrix A, looping on blocks the size of matrix_ks, only iterate on upper triangle since
!  the matrix is symmetric
      DO i = 1, ndo_mo
         DO j = i, ndo_mo

!           Loop over the blocks of the symmetric submatrices KS and S  and copy one by one into A
            DO iblk = 1, nblk_row
               DO jblk = iblk, nblk_row

                  found_block = .FALSE.
                  found_block_ks = .FALSE.
                  ALLOCATE (work_block(row_blk_size(iblk)*row_blk_size(jblk)))

!                 Retrieve S_pq
                  CALL dbcsr_get_block_p(matrix_s(1)%matrix, iblk, jblk, work_block_s, found_block)

!                 Take -epsilon_ij*S_pq
                  IF (found_block) THEN
                     work_block = -(energy_evals(i, j))*work_block_s

                     IF (energy_evals(i, j) == 0.0_dp) THEN
!                       If the whole block is zero, then equivalent to no block at all
                        found_block = .FALSE.
                     END IF
                  END IF

!                 If on the block diagonal, take F_pq - epsilon_ii*S_pq
                  IF (i == j) THEN
                     CALL dbcsr_get_block_p(matrix_ks(1)%matrix, iblk, jblk, work_block_ks, &
                                            found_block_ks)

                     IF (found_block .AND. found_block_ks) THEN
                        work_block = work_block+work_block_ks
                     ELSE IF ((.NOT. found_block) .AND. found_block_ks) THEN
                        work_block = work_block_ks
                     END IF
                  END IF

!                 Put the block in A if either found in KS,S or both
                  IF (found_block .OR. found_block_ks) THEN

                     CALL dbcsr_put_block(matrix_a, (i-1)*nblk_row+iblk, (j-1)*nblk_row+jblk, &
                                          work_block)
!                    Need to fill the whole upper part of A, by making sure that non-diagonal
!                    sub-blocks are coiped twice in the upper diagonal part of A
                     IF ((iblk .NE. jblk) .AND. (i .NE. j)) THEN
                        CALL dbcsr_put_block(matrix_a, (i-1)*nblk_row+jblk, (j-1)*nblk_row+iblk, &
                                             work_block)
                     END IF
                  END IF

                  DEALLOCATE (work_block)
                  NULLIFY (work_block_ks, work_block_s)
               END DO ! jblk
            END DO ! iblk

         END DO ! j
      END DO ! i

!  Finalize
      CALL dbcsr_finalize(matrix_a)

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

      CHARACTER(len=*), PARAMETER :: routineN = 'build_matrix_g', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iblk, jblk, nblk_row, ndo_mo
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size
      LOGICAL                                            :: found_block
      REAL(dp), DIMENSION(:), POINTER                    :: work_block
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s

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
!  their repective blocks which have the standard cp2k structure.
      DO i = 1, ndo_mo
!        loop over blocks of matrix_s
         DO iblk = 1, nblk_row
            DO jblk = 1, nblk_row

               CALL dbcsr_get_block_p(matrix_s(1)%matrix, iblk, jblk, work_block, found_block)

!              if non zero-block, copy its content in G
               IF (found_block) THEN
                  CALL dbcsr_put_block(matrix_g, (i-1)*nblk_row+iblk, (i-1)*nblk_row+jblk, &
                                       work_block)
               END IF

               NULLIFY (work_block)

            END DO
         END DO
      END DO

!  Finalize
      CALL dbcsr_finalize(matrix_g)

      CALL timestop(handle)

   END SUBROUTINE build_matrix_g

END MODULE xas_tdp_utils

