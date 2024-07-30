!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the full diagonalization of GW + Bethe-Salpeter for computing
!> electronic excitations
!> \par History
!>      10.2023 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_full_diag

   USE bse_util,                        ONLY: comp_eigvec_coeff_BSE,&
                                              compute_absorption_spectrum,&
                                              filter_eigvec_contrib,&
                                              fm_general_add_bse,&
                                              get_oscillator_strengths
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                              cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE input_constants,                 ONLY: bse_singlet,&
                                              bse_triplet
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_full_diag'

   PUBLIC :: create_A, diagonalize_A, create_B, create_hermitian_form_of_ABBA, &
             diagonalize_C

CONTAINS

! **************************************************************************************************
!> \brief Matrix A constructed from GW energies and 3c-B-matrices (cf. subroutine mult_B_with_W)
!>   A_ia,jb = (ε_a-ε_i) δ_ij δ_ab + α * v_ia,jb - W_ij,ab
!>   ε_a, ε_i are GW singleparticle energies from Eigenval_reduced
!>   α is a spin-dependent factor
!>   v_ia,jb = \sum_P B^P_ia B^P_jb (unscreened Coulomb interaction)
!>   W_ij,ab = \sum_P \bar{B}^P_ij B^P_ab (screened Coulomb interaction)
!> \param fm_mat_S_ia_bse ...
!> \param fm_mat_S_bar_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_A ...
!> \param Eigenval ...
!> \param unit_nr ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param mp2_env ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE create_A(fm_mat_S_ia_bse, fm_mat_S_bar_ij_bse, fm_mat_S_ab_bse, &
                       fm_A, Eigenval, unit_nr, &
                       homo, virtual, dimen_RI, mp2_env, &
                       para_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S_ia_bse, fm_mat_S_bar_ij_bse, &
                                                            fm_mat_S_ab_bse
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_A
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER, INTENT(IN)                                :: unit_nr, homo, virtual, dimen_RI
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      TYPE(mp_para_env_type), INTENT(INOUT)              :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_A'

      INTEGER                                            :: a_virt_row, handle, i_occ_row, &
                                                            i_row_global, ii, j_col_global, jj, &
                                                            ncol_local_A, nrow_local_A
      INTEGER, DIMENSION(4)                              :: reordering
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_A, row_indices_A
      REAL(KIND=dp)                                      :: alpha, eigen_diff
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_A, fm_struct_W
      TYPE(cp_fm_type)                                   :: fm_W

      CALL timeset(routineN, handle)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A10)') 'BSE|DEBUG|', 'Creating A'
      END IF

      !Determines factor of exchange term, depending on requested spin configuration (cf. input_constants.F)
      SELECT CASE (mp2_env%ri_g0w0%bse_spin_config)
      CASE (bse_singlet)
         alpha = 2.0_dp
      CASE (bse_triplet)
         alpha = 0.0_dp
      END SELECT

      ! create the blacs env for ij matrices (NOT fm_mat_S_ia_bse%matrix_struct related parallel_gemms!)
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

      ! We have to use the same blacs_env for A as for the matrices fm_mat_S_ia_bse from RPA
      ! Logic: A_ia,jb = (ε_a-ε_i) δ_ij δ_ab + α * v_ia,jb - W_ij,ab
      ! We create v_ia,jb and W_ij,ab, then we communicate entries from local W_ij,ab
      ! to the full matrix v_ia,jb. By adding these and the energy diffenences: v_ia,jb -> A_ia,jb
      ! We use the A matrix already from the start instead of v
      CALL cp_fm_struct_create(fm_struct_A, context=fm_mat_S_ia_bse%matrix_struct%context, nrow_global=homo*virtual, &
                               ncol_global=homo*virtual, para_env=fm_mat_S_ia_bse%matrix_struct%para_env)
      CALL cp_fm_create(fm_A, fm_struct_A, name="fm_A_iajb")
      CALL cp_fm_set_all(fm_A, 0.0_dp)

      CALL cp_fm_struct_create(fm_struct_W, context=fm_mat_S_ab_bse%matrix_struct%context, nrow_global=homo**2, &
                               ncol_global=virtual**2, para_env=fm_mat_S_ab_bse%matrix_struct%para_env)
      CALL cp_fm_create(fm_W, fm_struct_W, name="fm_W_ijab")
      CALL cp_fm_set_all(fm_W, 0.0_dp)

      ! Create A matrix from GW Energies, v_ia,jb and W_ij,ab (different blacs_env!)
      ! v_ia,jb, which is directly initialized in A (with a factor of alpha)
      ! v_ia,jb = \sum_P B^P_ia B^P_jb
      CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=alpha, &
                         matrix_a=fm_mat_S_ia_bse, matrix_b=fm_mat_S_ia_bse, beta=0.0_dp, &
                         matrix_c=fm_A)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A16)') 'BSE|DEBUG|', 'Allocated A_iajb'
      END IF

      !W_ij,ab = \sum_P \bar{B}^P_ij B^P_ab
      CALL parallel_gemm(transa="T", transb="N", m=homo**2, n=virtual**2, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_S_bar_ij_bse, matrix_b=fm_mat_S_ab_bse, beta=0.0_dp, &
                         matrix_c=fm_W)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A16)') 'BSE|DEBUG|', 'Allocated W_ijab'
      END IF

      ! We start by moving data from local parts of W_ij,ab to the full matrix A_ia,jb using buffers
      CALL cp_fm_get_info(matrix=fm_A, &
                          nrow_local=nrow_local_A, &
                          ncol_local=ncol_local_A, &
                          row_indices=row_indices_A, &
                          col_indices=col_indices_A)
      ! Writing -1.0_dp * W_ij,ab to A_ia,jb, i.e. beta = -1.0_dp,
      ! W_ij,ab: nrow_secidx_in  = homo,    ncol_secidx_in  = virtual
      ! A_ia,jb: nrow_secidx_out = virtual, ncol_secidx_out = virtual
      reordering = (/1, 3, 2, 4/)
      CALL fm_general_add_bse(fm_A, fm_W, -1.0_dp, homo, virtual, &
                              virtual, virtual, unit_nr, reordering, mp2_env)
      !full matrix W is not needed anymore, release it to save memory
      CALL cp_fm_struct_release(fm_struct_W)
      CALL cp_fm_release(fm_W)
      CALL cp_fm_struct_release(fm_struct_A)

      !Now add the energy differences (ε_a-ε_i) on the diagonal (i.e. δ_ij δ_ab) of A_ia,jb
      DO ii = 1, nrow_local_A

         i_row_global = row_indices_A(ii)

         DO jj = 1, ncol_local_A

            j_col_global = col_indices_A(jj)

            IF (i_row_global == j_col_global) THEN
               i_occ_row = (i_row_global - 1)/virtual + 1
               a_virt_row = MOD(i_row_global - 1, virtual) + 1
               eigen_diff = Eigenval(a_virt_row + homo) - Eigenval(i_occ_row)
               fm_A%local_data(ii, jj) = fm_A%local_data(ii, jj) + eigen_diff

            END IF
         END DO
      END DO

      CALL cp_fm_struct_release(fm_struct_A)
      CALL cp_fm_struct_release(fm_struct_W)

      CALL cp_blacs_env_release(blacs_env)

      CALL timestop(handle)

   END SUBROUTINE create_A

! **************************************************************************************************
!> \brief Matrix B constructed from 3c-B-matrices (cf. subroutine mult_B_with_W)
!>   B_ia,jb = α * v_ia,jb - W_ib,aj
!>   α is a spin-dependent factor
!>   v_ia,jb = \sum_P B^P_ia B^P_jb (unscreened Coulomb interaction)
!>   W_ib,aj = \sum_P \bar{B}^P_ib B^P_aj (screened Coulomb interaction)
!> \param fm_mat_S_ia_bse ...
!> \param fm_mat_S_bar_ia_bse ...
!> \param fm_B ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE create_B(fm_mat_S_ia_bse, fm_mat_S_bar_ia_bse, fm_B, &
                       homo, virtual, dimen_RI, unit_nr, mp2_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S_ia_bse, fm_mat_S_bar_ia_bse
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_B
      INTEGER, INTENT(IN)                                :: homo, virtual, dimen_RI, unit_nr
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_B'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(4)                              :: reordering
      REAL(KIND=dp)                                      :: alpha
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_v
      TYPE(cp_fm_type)                                   :: fm_W

      CALL timeset(routineN, handle)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A10)') 'BSE|DEBUG|', 'Creating B'
      END IF

      ! Determines factor of exchange term, depending on requested spin configuration (cf. input_constants.F)
      SELECT CASE (mp2_env%ri_g0w0%bse_spin_config)
      CASE (bse_singlet)
         alpha = 2.0_dp
      CASE (bse_triplet)
         alpha = 0.0_dp
      END SELECT

      CALL cp_fm_struct_create(fm_struct_v, context=fm_mat_S_ia_bse%matrix_struct%context, nrow_global=homo*virtual, &
                               ncol_global=homo*virtual, para_env=fm_mat_S_ia_bse%matrix_struct%para_env)
      CALL cp_fm_create(fm_B, fm_struct_v, name="fm_B_iajb")
      CALL cp_fm_set_all(fm_B, 0.0_dp)

      CALL cp_fm_create(fm_W, fm_struct_v, name="fm_W_ibaj")
      CALL cp_fm_set_all(fm_W, 0.0_dp)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A16)') 'BSE|DEBUG|', 'Allocated B_iajb'
      END IF
      ! v_ia,jb = \sum_P B^P_ia B^P_jb
      CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=alpha, &
                         matrix_a=fm_mat_S_ia_bse, matrix_b=fm_mat_S_ia_bse, beta=0.0_dp, &
                         matrix_c=fm_B)

      ! W_ib,aj = \sum_P \bar{B}^P_ib B^P_aj
      CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_S_bar_ia_bse, matrix_b=fm_mat_S_ia_bse, beta=0.0_dp, &
                         matrix_c=fm_W)
      ! from W_ib,ja to A_ia,jb (formally: W_ib,aj, but our internal indexorder is different)
      ! Writing -1.0_dp * W_ib,ja to A_ia,jb, i.e. beta = -1.0_dp,
      ! W_ib,ja: nrow_secidx_in  = virtual,    ncol_secidx_in  = virtual
      ! A_ia,jb: nrow_secidx_out = virtual, ncol_secidx_out = virtual
      reordering = (/1, 4, 3, 2/)
      CALL fm_general_add_bse(fm_B, fm_W, -1.0_dp, virtual, virtual, &
                              virtual, virtual, unit_nr, reordering, mp2_env)

      CALL cp_fm_release(fm_W)
      CALL cp_fm_struct_release(fm_struct_v)
      CALL timestop(handle)

   END SUBROUTINE create_B

   ! **************************************************************************************************
!> \brief Construct Matrix C=(A-B)^0.5 (A+B) (A-B)^0.5 to solve full BSE matrix as a hermitian problem
!>   (cf. Eq. (A7) in F. Furche J. Chem. Phys., Vol. 114, No. 14, (2001)).
!>   We keep fm_sqrt_A_minus_B and fm_inv_sqrt_A_minus_B for print of singleparticle transitions
!>   of ABBA as described in Eq. (A10) in F. Furche J. Chem. Phys., Vol. 114, No. 14, (2001).
!> \param fm_A ...
!> \param fm_B ...
!> \param fm_C ...
!> \param fm_sqrt_A_minus_B ...
!> \param fm_inv_sqrt_A_minus_B ...
!> \param homo ...
!> \param virtual ...
!> \param unit_nr ...
!> \param mp2_env ...
!> \param diag_est ...
! **************************************************************************************************
   SUBROUTINE create_hermitian_form_of_ABBA(fm_A, fm_B, fm_C, &
                                            fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                                            homo, virtual, unit_nr, mp2_env, diag_est)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_A, fm_B
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_C, fm_sqrt_A_minus_B, &
                                                            fm_inv_sqrt_A_minus_B
      INTEGER, INTENT(IN)                                :: homo, virtual, unit_nr
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      REAL(KIND=dp), INTENT(IN)                          :: diag_est

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_hermitian_form_of_ABBA'

      INTEGER                                            :: dim_mat, handle, n_dependent
      REAL(KIND=dp), DIMENSION(2)                        :: eigvals_AB_diff
      TYPE(cp_fm_type)                                   :: fm_A_minus_B, fm_A_plus_B, fm_dummy, &
                                                            fm_work_product

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A25,A39,ES6.0,A3)') 'BSE|', 'Diagonalizing aux. matrix', &
            ' with size of A. This will take around ', diag_est, " s."
      END IF

      ! Create work matrices, which will hold A+B and A-B and their powers
      ! C is created afterwards to save memory
      ! Final result: C = (A-B)^0.5             (A+B)              (A-B)^0.5              EQ.I
      !                   \_______/             \___/              \______/
      !               fm_sqrt_A_minus_B      fm_A_plus_B     fm_sqrt_A_minus_B
      !                    (EQ.Ia)             (EQ.Ib)              (EQ.Ia)
      ! Intermediate work matrices:
      ! fm_inv_sqrt_A_minus_B: (A-B)^-0.5                                                 EQ.II
      ! fm_A_minus_B: (A-B)                                                               EQ.III
      ! fm_work_product: (A-B)^0.5 (A+B) from (EQ.Ia) and (EQ.Ib)                         EQ.IV
      CALL cp_fm_create(fm_A_plus_B, fm_A%matrix_struct)
      CALL cp_fm_to_fm(fm_A, fm_A_plus_B)
      CALL cp_fm_create(fm_A_minus_B, fm_A%matrix_struct)
      CALL cp_fm_to_fm(fm_A, fm_A_minus_B)
      CALL cp_fm_create(fm_sqrt_A_minus_B, fm_A%matrix_struct)
      CALL cp_fm_set_all(fm_sqrt_A_minus_B, 0.0_dp)
      CALL cp_fm_create(fm_inv_sqrt_A_minus_B, fm_A%matrix_struct)
      CALL cp_fm_set_all(fm_inv_sqrt_A_minus_B, 0.0_dp)

      CALL cp_fm_create(fm_work_product, fm_A%matrix_struct)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A19)') 'BSE|DEBUG|', 'Created work arrays'
      END IF

      ! Add/Substract B (cf. EQs. Ib and III)
      CALL cp_fm_scale_and_add(1.0_dp, fm_A_plus_B, 1.0_dp, fm_B)
      CALL cp_fm_scale_and_add(1.0_dp, fm_A_minus_B, -1.0_dp, fm_B)

      ! cp_fm_power will overwrite matrix, therefore we create copies
      CALL cp_fm_to_fm(fm_A_minus_B, fm_inv_sqrt_A_minus_B)

      ! In order to avoid a second diagonalization (cp_fm_power), we create (A-B)^0.5 (EQ.Ia)
      ! from (A-B)^-0.5 (EQ.II) by multiplication with (A-B) (EQ.III) afterwards.

      ! Raise A-B to -0.5_dp, no quenching of eigenvectors, hence threshold=0.0_dp
      CALL cp_fm_create(fm_dummy, fm_A%matrix_struct)
      ! Create (A-B)^-0.5 (cf. EQ.II)
      CALL cp_fm_power(fm_inv_sqrt_A_minus_B, fm_dummy, -0.5_dp, 0.0_dp, n_dependent, eigvals=eigvals_AB_diff)
      CALL cp_fm_release(fm_dummy)
      ! Raise an error in case the the matrix A-B is not positive definite (i.e. negative eigenvalues)
      ! In this case, the procedure for hermitian form of ABBA is not applicable
      IF (eigvals_AB_diff(1) < 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "Matrix (A-B) is not positive definite. "// &
                       "Hermitian diagonalization of full ABBA matrix is ill-defined.")
      END IF

      ! We keep fm_inv_sqrt_A_minus_B for print of singleparticle transitions of ABBA
      ! We further create (A-B)^0.5 for the singleparticle transitions of ABBA
      ! Create (A-B)^0.5= (A-B)^-0.5 * (A-B) (EQ.Ia)
      dim_mat = homo*virtual
      CALL parallel_gemm("N", "N", dim_mat, dim_mat, dim_mat, 1.0_dp, fm_inv_sqrt_A_minus_B, fm_A_minus_B, 0.0_dp, &
                         fm_sqrt_A_minus_B)

      ! Compute and store LHS of C, i.e. (A-B)^0.5 (A+B) (EQ.IV)
      CALL parallel_gemm("N", "N", dim_mat, dim_mat, dim_mat, 1.0_dp, fm_sqrt_A_minus_B, fm_A_plus_B, 0.0_dp, &
                         fm_work_product)

      ! Release to save memory
      CALL cp_fm_release(fm_A_plus_B)
      CALL cp_fm_release(fm_A_minus_B)

      ! Now create full
      CALL cp_fm_create(fm_C, fm_A%matrix_struct)
      CALL cp_fm_set_all(fm_C, 0.0_dp)
      ! Compute C=(A-B)^0.5 (A+B) (A-B)^0.5 (EQ.I)
      CALL parallel_gemm("N", "N", dim_mat, dim_mat, dim_mat, 1.0_dp, fm_work_product, fm_sqrt_A_minus_B, 0.0_dp, &
                         fm_C)
      CALL cp_fm_release(fm_work_product)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A36)') 'BSE|DEBUG|', 'Filled C=(A-B)^0.5 (A+B) (A-B)^0.5'
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Solving eigenvalue equation C Z^n = (Ω^n)^2 Z^n .
!>   Here, the eigenvectors Z^n relate to X^n via
!>   Eq. (A10) in F. Furche J. Chem. Phys., Vol. 114, No. 14, (2001).
!> \param fm_C ...
!> \param homo ...
!> \param virtual ...
!> \param homo_irred ...
!> \param fm_sqrt_A_minus_B ...
!> \param fm_inv_sqrt_A_minus_B ...
!> \param unit_nr ...
!> \param diag_est ...
!> \param mp2_env ...
!> \param qs_env ...
!> \param mo_coeff ...
! **************************************************************************************************
   SUBROUTINE diagonalize_C(fm_C, homo, virtual, homo_irred, &
                            fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B, &
                            unit_nr, diag_est, mp2_env, qs_env, mo_coeff)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_C
      INTEGER, INTENT(IN)                                :: homo, virtual, homo_irred
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_sqrt_A_minus_B, fm_inv_sqrt_A_minus_B
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), INTENT(IN)                          :: diag_est
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'diagonalize_C'

      INTEGER                                            :: diag_info, handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens
      TYPE(cp_fm_type)                                   :: fm_eigvec_X, fm_eigvec_Y, fm_eigvec_Z, &
                                                            fm_mat_eigvec_transform_diff, &
                                                            fm_mat_eigvec_transform_sum

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A17,A22,ES6.0,A3)') 'BSE|', 'Diagonalizing C. ', &
            'This will take around ', diag_est, ' s.'
      END IF

      !We have now the full matrix C=(A-B)^0.5 (A+B) (A-B)^0.5
      !Now: Diagonalize it
      CALL cp_fm_create(fm_eigvec_Z, fm_C%matrix_struct)

      ALLOCATE (Exc_ens(homo*virtual))

      CALL choose_eigv_solver(fm_C, fm_eigvec_Z, Exc_ens, diag_info)
      CPASSERT(diag_info == 0)
      ! C could have negative eigenvalues, since we do not explicitly check A+B
      ! for positive definiteness (would make another O(N^6) Diagon. necessary)
      ! Instead, we include a check here
      IF (Exc_ens(1) < 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "Matrix C=(A-B)^0.5 (A+B) (A-B)^0.5 has negative eigenvalues, i.e. "// &
                       "(A+B) is not positive definite.")
      END IF
      Exc_ens = SQRT(Exc_ens)

      ! Prepare eigenvector for interpretation of singleparticle transitions
      ! Compare: F. Furche J. Chem. Phys., Vol. 114, No. 14, (2001)
      ! We aim for the upper part of the vector (X,Y) for a direct comparison with the TDA result

      ! Following Furche, we basically use Eqs. (A10): First, we multiply
      ! the (A-B)^+-0.5 with eigenvectors and then the eigenvalues
      ! One has to be careful about the index structure, since the eigenvector matrix is not symmetric anymore!

      ! First, Eq. I from (A10) from Furche: (X+Y)_n = (Ω_n)^-0.5 (A-B)^0.5 T_n
      CALL cp_fm_create(fm_mat_eigvec_transform_sum, fm_C%matrix_struct)
      CALL cp_fm_set_all(fm_mat_eigvec_transform_sum, 0.0_dp)
      CALL parallel_gemm(transa="N", transb="N", m=homo*virtual, n=homo*virtual, k=homo*virtual, alpha=1.0_dp, &
                         matrix_a=fm_sqrt_A_minus_B, matrix_b=fm_eigvec_Z, beta=0.0_dp, &
                         matrix_c=fm_mat_eigvec_transform_sum)
      CALL cp_fm_release(fm_sqrt_A_minus_B)
      CALL comp_eigvec_coeff_BSE(fm_mat_eigvec_transform_sum, Exc_ens, -0.5_dp, gamma=2.0_dp, do_transpose=.TRUE.)

      ! Second, Eq. II from (A10) from Furche: (X-Y)_n = (Ω_n)^0.5 (A-B)^-0.5 T_n
      CALL cp_fm_create(fm_mat_eigvec_transform_diff, fm_C%matrix_struct)
      CALL cp_fm_set_all(fm_mat_eigvec_transform_diff, 0.0_dp)
      CALL parallel_gemm(transa="N", transb="N", m=homo*virtual, n=homo*virtual, k=homo*virtual, alpha=1.0_dp, &
                         matrix_a=fm_inv_sqrt_A_minus_B, matrix_b=fm_eigvec_Z, beta=0.0_dp, &
                         matrix_c=fm_mat_eigvec_transform_diff)
      CALL cp_fm_release(fm_inv_sqrt_A_minus_B)
      CALL cp_fm_release(fm_eigvec_Z)
      CALL comp_eigvec_coeff_BSE(fm_mat_eigvec_transform_diff, Exc_ens, 0.5_dp, gamma=2.0_dp, do_transpose=.TRUE.)

      ! Now, we add the two equations to obtain X_n
      ! Add overwrites the first argument, therefore we copy it beforehand
      CALL cp_fm_create(fm_eigvec_X, fm_C%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_eigvec_transform_sum, fm_eigvec_X)
      CALL cp_fm_scale_and_add(1.0_dp, fm_eigvec_X, 1.0_dp, fm_mat_eigvec_transform_diff)

      ! Now, we subtract the two equations to obtain Y_n
      ! Add overwrites the first argument, therefore we copy it beforehand
      CALL cp_fm_create(fm_eigvec_Y, fm_C%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_eigvec_transform_sum, fm_eigvec_Y)
      CALL cp_fm_scale_and_add(1.0_dp, fm_eigvec_Y, -1.0_dp, fm_mat_eigvec_transform_diff)

      !Cleanup
      CALL cp_fm_release(fm_mat_eigvec_transform_diff)
      CALL cp_fm_release(fm_mat_eigvec_transform_sum)

      CALL success_message(Exc_ens, fm_eigvec_X, mp2_env, qs_env, mo_coeff, &
                           homo, virtual, homo_irred, unit_nr, &
                           .FALSE., fm_eigvec_Y)

      DEALLOCATE (Exc_ens)
      CALL cp_fm_release(fm_eigvec_X)
      CALL cp_fm_release(fm_eigvec_Y)

      CALL timestop(handle)

   END SUBROUTINE diagonalize_C

! **************************************************************************************************
!> \brief Solving hermitian eigenvalue equation A X^n = Ω^n X^n
!> \param fm_A ...
!> \param homo ...
!> \param virtual ...
!> \param homo_irred ...
!> \param unit_nr ...
!> \param diag_est ...
!> \param mp2_env ...
!> \param qs_env ...
!> \param mo_coeff ...
! **************************************************************************************************
   SUBROUTINE diagonalize_A(fm_A, homo, virtual, homo_irred, &
                            unit_nr, diag_est, mp2_env, qs_env, mo_coeff)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_A
      INTEGER, INTENT(IN)                                :: homo, virtual, homo_irred, unit_nr
      REAL(KIND=dp), INTENT(IN)                          :: diag_est
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'diagonalize_A'

      INTEGER                                            :: diag_info, handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens
      TYPE(cp_fm_type)                                   :: fm_eigvec

      CALL timeset(routineN, handle)

      !Continue with formatting of subroutine create_A
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A17,A22,ES6.0,A3)') 'BSE|', 'Diagonalizing A. ', &
            'This will take around ', diag_est, ' s.'
      END IF

      !We have now the full matrix A_iajb, distributed over all ranks
      !Now: Diagonalize it
      CALL cp_fm_create(fm_eigvec, fm_A%matrix_struct)

      ALLOCATE (Exc_ens(homo*virtual))

      CALL choose_eigv_solver(fm_A, fm_eigvec, Exc_ens, diag_info)
      CPASSERT(diag_info == 0)
      CALL success_message(Exc_ens, fm_eigvec, mp2_env, qs_env, mo_coeff, &
                           homo, virtual, homo_irred, unit_nr, .TRUE.)

      CALL cp_fm_release(fm_eigvec)
      DEALLOCATE (Exc_ens)

      CALL timestop(handle)

   END SUBROUTINE diagonalize_A

! **************************************************************************************************
!> \brief Prints the success message (incl. energies) for full diag of BSE (TDA/full ABBA via flag)
!> \param Exc_ens ...
!> \param fm_eigvec_X ...
!> \param mp2_env ...
!> \param qs_env ...
!> \param mo_coeff ...
!> \param homo ...
!> \param virtual ...
!> \param homo_irred ...
!> \param unit_nr ...
!> \param flag_TDA ...
!> \param fm_eigvec_Y ...
! **************************************************************************************************
   SUBROUTINE success_message(Exc_ens, fm_eigvec_X, mp2_env, qs_env, mo_coeff, &
                              homo, virtual, homo_irred, unit_nr, &
                              flag_TDA, fm_eigvec_Y)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_eigvec_X
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      INTEGER                                            :: homo, virtual, homo_irred, unit_nr
      LOGICAL, OPTIONAL                                  :: flag_TDA
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: fm_eigvec_Y

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'success_message'

      CHARACTER(LEN=10)                                  :: info_approximation, multiplet
      INTEGER                                            :: handle, i_exc, k, num_entries
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: idx_homo, idx_virt
      REAL(KIND=dp)                                      :: alpha
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigvec_entries, oscill_str
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: abs_spectrum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dipole_exc

      CALL timeset(routineN, handle)

      !Prepare variables for printing
      IF (mp2_env%ri_g0w0%bse_spin_config == 0) THEN
         multiplet = "Singlet"
         alpha = 2.0_dp
      ELSE
         multiplet = "Triplet"
         alpha = 0.0_dp
      END IF
      IF (.NOT. PRESENT(flag_TDA)) THEN
         flag_TDA = .FALSE.
      END IF
      IF (flag_TDA) THEN
         info_approximation = " -TDA- "
      ELSE
         info_approximation = "-ABBA-"
      END IF

      !Get information about eigenvector

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         IF (flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '**************************************************************************'
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '*   Bethe Salpeter equation (BSE) with Tamm Dancoff approximation (TDA)  *'
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '**************************************************************************'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A48,A23)') 'BSE|', 'The excitations are calculated by diagonalizing ', &
               'the BSE within the TDA:'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T29,A16)') 'BSE|', 'A X^n = Ω^n X^n'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A23)') 'BSE|', 'i.e. in index notation:'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A41)') 'BSE|', 'sum_jb ( A_ia,jb   X_jb^n ) = Ω^n X_ia^n'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A30)') 'BSE|', 'prelim Ref.: Eq. (36) with B=0'
            WRITE (unit_nr, '(T2,A4,T7,A71)') 'BSE|', 'in PRB 92,045209 (2015); http://dx.doi.org/10.1103/PhysRevB.92.045209 .'
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '**************************************************************************'
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '*      Full ("ABBA") Bethe Salpeter equation (BSE) (i.e. without TDA)    *'
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '**************************************************************************'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A48,A24)') 'BSE|', 'The excitations are calculated by diagonalizing ', &
               'the BSE without the TDA:'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T22,A30)') 'BSE|', '|A B| |X^n|       |1  0| |X^n|'
            WRITE (unit_nr, '(T2,A4,T22,A31)') 'BSE|', '|B A| |Y^n| = Ω^n |0 -1| |Y^n|'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A23)') 'BSE|', 'i.e. in index notation:'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A62)') 'BSE|', '  sum_jb ( A_ia,jb   X_jb^n + B_ia,jb   Y_jb^n ) = Ω^n X_ia^n'
            WRITE (unit_nr, '(T2,A4,T7,A62)') 'BSE|', '- sum_jb ( B_ia,jb   X_jb^n + A_ia,jb   Y_jb^n ) = Ω^n Y_ia^n'
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A4,T18,A42,T70,A1,I4,A1,I4,A1)') 'BSE|', 'i,j:', &
            'occupied molecular orbitals, i.e. state in', '[', homo_irred - homo + 1, ',', homo_irred, ']'
         WRITE (unit_nr, '(T2,A4,T7,A4,T18,A44,T70,A1,I4,A1,I4,A1)') 'BSE|', 'a,b:', &
            'unoccupied molecular orbitals, i.e. state in', '[', homo_irred + 1, ',', homo_irred + virtual, ']'
         WRITE (unit_nr, '(T2,A4,T7,A2,T18,A16)') 'BSE|', 'n:', 'Excitation index'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A58)') 'BSE|', 'A_ia,jb = (ε_a-ε_i) δ_ij δ_ab + α * v_ia,jb - W_ij,ab'
         IF (.NOT. flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4,T7,A32)') 'BSE|', 'B_ia,jb = α * v_ia,jb - W_ib,aj'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A35)') 'BSE|', 'prelim Ref.: Eqs. (24-27),(30),(35)'
            WRITE (unit_nr, '(T2,A4,T7,A71)') 'BSE|', 'in PRB 92,045209 (2015); http://dx.doi.org/10.1103/PhysRevB.92.045209 .'
         END IF
         IF (.NOT. flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', 'The BSE is solved for Ω^n and X_ia^n as a hermitian problem, e.g. Eq.(42)'
            WRITE (unit_nr, '(T2,A4,T7,A71)') 'BSE|', 'in PRB 92,045209 (2015); http://dx.doi.org/10.1103/PhysRevB.92.045209 .'
            ! WRITE (unit_nr, '(T2,A4,T7,A69)') 'BSE|', 'C_ia,jb = sum_kc,ld ((A-B)^0.5)_ia,kc (A+B)_kc,ld ((A-B)^0.5)_ld,jb'
            ! WRITE (unit_nr, '(T2,A4)') 'BSE|'
            ! WRITE (unit_nr, '(T2,A4,T7,A58)') 'BSE|', '(X+Y)_ia,n = sum_jb,m  ((A-B)^0.5)ia,jb Z_jb,m (Ω_m)^-0.5'
            ! WRITE (unit_nr, '(T2,A4,T7,A58)') 'BSE|', '(X-Y)_ia,n = sum_jb,m  ((A-B)^-0.5)ia,jb Z_jb,m (Ω_m)^0.5'
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A7,T19,A23)') 'BSE|', 'ε_...:', 'GW quasiparticle energy'
         WRITE (unit_nr, '(T2,A4,T7,A7,T19,A15)') 'BSE|', 'δ_...:', 'Kronecker delta'
         WRITE (unit_nr, '(T2,A4,T7,A3,T19,A21)') 'BSE|', 'α:', 'spin-dependent factor (Singlet/Triplet)'
         WRITE (unit_nr, '(T2,A4,T7,A6,T18,A34)') 'BSE|', 'v_...:', 'Electron-hole exchange interaction'
         WRITE (unit_nr, '(T2,A4,T7,A6,T18,A27)') 'BSE|', 'W_...:', 'Screened direct interaction'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A47,A7,A9,F3.1)') 'BSE|', &
            'The spin-dependent factor is for the requested ', multiplet, " is α = ", alpha
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         IF (flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4,T7,A56)') 'BSE|', 'Excitation energies from solving the BSE within the TDA:'
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A57)') 'BSE|', 'Excitation energies from solving the BSE without the TDA:'
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T11,A12,T26,A11,T44,A8,T55,A27)') 'BSE|', &
            'Excitation n', "Spin Config", 'TDA/ABBA', 'Excitation energy Ω^n (eV)'
      END IF
      !prints actual energies values
      IF (unit_nr > 0) THEN
         DO i_exc = 1, MIN(homo*virtual, mp2_env%ri_g0w0%num_print_exc)
            WRITE (unit_nr, '(T2,A4,T7,I16,T30,A7,T46,A6,T59,F22.4)') &
               'BSE|', i_exc, multiplet, info_approximation, Exc_ens(i_exc)*evolt
         END DO
      END IF
      !prints single particle transitions
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'

         WRITE (unit_nr, '(T2,A4,T7,A70)') &
            'BSE|', "Excitations are built up by the following single-particle transitions,"
         WRITE (unit_nr, '(T2,A4,T7,A42,F5.2,A2)') &
            'BSE|', "neglecting contributions where |X_ia^n| < ", mp2_env%ri_g0w0%eps_x, " :"

         WRITE (unit_nr, '(T2,A4,T15,A27,I5,A13,I5,A3)') 'BSE|', '-- Quick reminder: HOMO i =', &
            homo_irred, ' and LUMO a =', homo_irred + 1, " --"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A12,T25,A1,T27,A2,T34,A1,T54,A8,T71,A10)') &
            "BSE|", "Excitation n", "i", "=>", "a", 'TDA/ABBA', "|X_ia^n|"
      END IF
      DO i_exc = 1, MIN(homo*virtual, mp2_env%ri_g0w0%num_print_exc)
         !Iterate through eigenvector and print values above threshold
         CALL filter_eigvec_contrib(fm_eigvec_X, idx_homo, idx_virt, eigvec_entries, &
                                    i_exc, virtual, num_entries, mp2_env)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            DO k = 1, num_entries
               WRITE (unit_nr, '(T2,A4,T14,I5,T21,I5,T27,A2,T30,I5,T56,A6,T65,F16.4)') &
                  "BSE|", i_exc, homo_irred - homo + idx_homo(k), "=>", &
                  homo_irred + idx_virt(k), info_approximation, ABS(eigvec_entries(k))
            END DO
         END IF

         DEALLOCATE (idx_homo)
         DEALLOCATE (idx_virt)
         DEALLOCATE (eigvec_entries)
      END DO
      ! Compute BSE dipoles and oscillator strengths
      CALL get_oscillator_strengths(fm_eigvec_X, Exc_ens, &
                                    dipole_exc, oscill_str, &
                                    mp2_env, qs_env, mo_coeff, &
                                    homo, virtual, unit_nr, &
                                    fm_eigvec_Y)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A63)') &
            'BSE|', "Dipoles d_r^n and oscillator strength f^n of excitation level n"
         WRITE (unit_nr, '(T2,A4,T7,A17)') &
            'BSE|', "are obtained from"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         IF (flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4,T12,A51)') &
               'BSE|', "d_r^n = sqrt(2) sum_ia < ψ_i | µ | ψ_a >  X_ia^n"
         ELSE
            WRITE (unit_nr, '(T2,A4,T12,A63)') &
               'BSE|', "d_r^n = sqrt(2) sum_ia < ψ_i | µ | ψ_a > ( X_ia^n + Y_ia^n )"
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T20,A44)') &
            'BSE|', "f^n = 2/3 * Ω^n sum_r∈(x,y,z) ( d_r^n )^2"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A19)') &
            'BSE|', "where we introduced"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A5,T15,A28)') &
            'BSE|', "ψ_i:", "occupied molecular orbitals,"
         WRITE (unit_nr, '(T2,A4,T7,A5,T15,A28)') &
            'BSE|', "ψ_a:", "empty molecular orbitals and"
         WRITE (unit_nr, '(T2,A4,T7,A3,T15,A16)') &
            'BSE|', "µ:", "dipole operator."
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A28)') &
            'BSE|', "prelim Ref.: Eqs. (23), (24)"
         WRITE (unit_nr, '(T2,A4,T7,A71)') &
            'BSE|', "in J. Chem. Phys. 152, 044105 (2020); https://doi.org/10.1063/1.5123290"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T8,A12,T22,A8,T35,A8,T45,A8,T55,A8,T65,A16)') 'BSE|', &
            'Excitation n', "TDA/ABBA", "Dipole X", "Dipole Y", "Dipole Z", 'Osc. strength'
         DO i_exc = 1, MIN(homo*virtual, mp2_env%ri_g0w0%num_print_exc)
            WRITE (unit_nr, '(T2,A4,T8,I12,T24,A6,T35,F8.3,T45,F8.3,T55,F8.3,T65,F16.3)') &
               'BSE|', i_exc, info_approximation, dipole_exc(1, 1, i_exc), dipole_exc(2, 1, i_exc), &
               dipole_exc(3, 1, i_exc), oscill_str(i_exc)
         END DO
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
      END IF

      ! Compute and print the absorption spectrum to external file
      IF (mp2_env%ri_g0w0%bse_print_spectrum) THEN
         CALL compute_absorption_spectrum(abs_spectrum, oscill_str, Exc_ens, &
                                          info_approximation, mp2_env)
         DEALLOCATE (abs_spectrum)
      END IF
      DEALLOCATE (oscill_str, dipole_exc)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
      END IF

      CALL timestop(handle)

   END SUBROUTINE success_message

END MODULE bse_full_diag
