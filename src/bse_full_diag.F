!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
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

   USE bse_print,                       ONLY: print_excitation_energies,&
                                              print_exciton_descriptors,&
                                              print_optical_properties,&
                                              print_output_header,&
                                              print_transition_amplitudes
   USE bse_properties,                  ONLY: calculate_NTOs,&
                                              exciton_descr_type,&
                                              get_exciton_descriptors,&
                                              get_oscillator_strengths
   USE bse_util,                        ONLY: comp_eigvec_coeff_BSE,&
                                              fm_general_add_bse,&
                                              get_multipoles_mo,&
                                              reshuffle_eigvec
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
   USE input_constants,                 ONLY: bse_screening_alpha,&
                                              bse_screening_rpa,&
                                              bse_singlet,&
                                              bse_triplet
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
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
      REAL(KIND=dp)                                      :: alpha, alpha_screening, eigen_diff
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_A, fm_struct_W
      TYPE(cp_fm_type)                                   :: fm_W

      CALL timeset(routineN, handle)

      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A10)') 'BSE|DEBUG|', 'Creating A'
      END IF

      !Determines factor of exchange term, depending on requested spin configuration (cf. input_constants.F)
      SELECT CASE (mp2_env%bse%bse_spin_config)
      CASE (bse_singlet)
         alpha = 2.0_dp
      CASE (bse_triplet)
         alpha = 0.0_dp
      END SELECT

      IF (mp2_env%bse%screening_method == bse_screening_alpha) THEN
         alpha_screening = mp2_env%bse%screening_factor
      ELSE
         alpha_screening = 1.0_dp
      END IF

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

      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A16)') 'BSE|DEBUG|', 'Allocated A_iajb'
      END IF

      ! If infinite screening is applied, fm_W is simply 0 - Otherwise it needs to be computed from 3c integrals
      IF (mp2_env%bse%screening_method /= bse_screening_rpa) THEN
         !W_ij,ab = \sum_P \bar{B}^P_ij B^P_ab
         CALL parallel_gemm(transa="T", transb="N", m=homo**2, n=virtual**2, k=dimen_RI, alpha=alpha_screening, &
                            matrix_a=fm_mat_S_bar_ij_bse, matrix_b=fm_mat_S_ab_bse, beta=0.0_dp, &
                            matrix_c=fm_W)
      END IF

      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
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

      ! If infinite screening is applied, fm_W is simply 0 - Otherwise it needs to be computed from 3c integrals
      IF (mp2_env%bse%screening_method /= bse_screening_rpa) THEN
         reordering = (/1, 3, 2, 4/)
         CALL fm_general_add_bse(fm_A, fm_W, -1.0_dp, homo, virtual, &
                                 virtual, virtual, unit_nr, reordering, mp2_env)
      END IF
      !full matrix W is not needed anymore, release it to save memory
      CALL cp_fm_release(fm_W)

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
      REAL(KIND=dp)                                      :: alpha, alpha_screening
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_v
      TYPE(cp_fm_type)                                   :: fm_W

      CALL timeset(routineN, handle)

      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A10)') 'BSE|DEBUG|', 'Creating B'
      END IF

      ! Determines factor of exchange term, depending on requested spin configuration (cf. input_constants.F)
      SELECT CASE (mp2_env%bse%bse_spin_config)
      CASE (bse_singlet)
         alpha = 2.0_dp
      CASE (bse_triplet)
         alpha = 0.0_dp
      END SELECT

      IF (mp2_env%bse%screening_method == bse_screening_alpha) THEN
         alpha_screening = mp2_env%bse%screening_factor
      ELSE
         alpha_screening = 1.0_dp
      END IF

      CALL cp_fm_struct_create(fm_struct_v, context=fm_mat_S_ia_bse%matrix_struct%context, nrow_global=homo*virtual, &
                               ncol_global=homo*virtual, para_env=fm_mat_S_ia_bse%matrix_struct%para_env)
      CALL cp_fm_create(fm_B, fm_struct_v, name="fm_B_iajb")
      CALL cp_fm_set_all(fm_B, 0.0_dp)

      CALL cp_fm_create(fm_W, fm_struct_v, name="fm_W_ibaj")
      CALL cp_fm_set_all(fm_W, 0.0_dp)

      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A16)') 'BSE|DEBUG|', 'Allocated B_iajb'
      END IF
      ! v_ia,jb = \sum_P B^P_ia B^P_jb
      CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=alpha, &
                         matrix_a=fm_mat_S_ia_bse, matrix_b=fm_mat_S_ia_bse, beta=0.0_dp, &
                         matrix_c=fm_B)

      ! If infinite screening is applied, fm_W is simply 0 - Otherwise it needs to be computed from 3c integrals
      IF (mp2_env%bse%screening_method /= bse_screening_rpa) THEN
         ! W_ib,aj = \sum_P \bar{B}^P_ib B^P_aj
         CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=alpha_screening, &
                            matrix_a=fm_mat_S_bar_ia_bse, matrix_b=fm_mat_S_ia_bse, beta=0.0_dp, &
                            matrix_c=fm_W)

         ! from W_ib,ja to A_ia,jb (formally: W_ib,aj, but our internal indexorder is different)
         ! Writing -1.0_dp * W_ib,ja to A_ia,jb, i.e. beta = -1.0_dp,
         ! W_ib,ja: nrow_secidx_in  = virtual,    ncol_secidx_in  = virtual
         ! A_ia,jb: nrow_secidx_out = virtual, ncol_secidx_out = virtual
         reordering = (/1, 4, 3, 2/)
         CALL fm_general_add_bse(fm_B, fm_W, -1.0_dp, virtual, virtual, &
                                 virtual, virtual, unit_nr, reordering, mp2_env)
      END IF

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

      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
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

      IF (unit_nr > 0 .AND. mp2_env%bse%bse_debug_print) THEN
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

      IF (diag_info /= 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "Diagonalization of C=(A-B)^0.5 (A+B) (A-B)^0.5 failed in BSE")
      END IF

      ! C could have negative eigenvalues, since we do not explicitly check A+B
      ! for positive definiteness (would make another O(N^6) Diagon. necessary)
      ! Instead, we include a check here
      IF (Exc_ens(1) < 0) THEN
         IF (unit_nr > 0) THEN
            CALL cp_abort(__LOCATION__, &
                          "Matrix C=(A-B)^0.5 (A+B) (A-B)^0.5 has negative eigenvalues, i.e. "// &
                          "(A+B) is not positive definite.")
         END IF
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
      ! This normalizes the eigenvectors
      CALL comp_eigvec_coeff_BSE(fm_mat_eigvec_transform_sum, Exc_ens, -0.5_dp, gamma=2.0_dp, do_transpose=.TRUE.)

      ! Second, Eq. II from (A10) from Furche: (X-Y)_n = (Ω_n)^0.5 (A-B)^-0.5 T_n
      CALL cp_fm_create(fm_mat_eigvec_transform_diff, fm_C%matrix_struct)
      CALL cp_fm_set_all(fm_mat_eigvec_transform_diff, 0.0_dp)
      CALL parallel_gemm(transa="N", transb="N", m=homo*virtual, n=homo*virtual, k=homo*virtual, alpha=1.0_dp, &
                         matrix_a=fm_inv_sqrt_A_minus_B, matrix_b=fm_eigvec_Z, beta=0.0_dp, &
                         matrix_c=fm_mat_eigvec_transform_diff)
      CALL cp_fm_release(fm_inv_sqrt_A_minus_B)
      CALL cp_fm_release(fm_eigvec_Z)

      ! This normalizes the eigenvectors
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

      CALL postprocess_bse(Exc_ens, fm_eigvec_X, mp2_env, qs_env, mo_coeff, &
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

      IF (diag_info /= 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "Diagonalization of A failed in TDA-BSE")
      END IF

      CALL postprocess_bse(Exc_ens, fm_eigvec, mp2_env, qs_env, mo_coeff, &
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
   SUBROUTINE postprocess_bse(Exc_ens, fm_eigvec_X, mp2_env, qs_env, mo_coeff, &
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

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'postprocess_bse'

      CHARACTER(LEN=10)                                  :: info_approximation, multiplet
      INTEGER                                            :: handle, i_exc, idir, n_moments_di, &
                                                            n_moments_quad
      REAL(KIND=dp)                                      :: alpha
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: oscill_str, ref_point_multipole
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: polarizability_residues, trans_mom_bse
      TYPE(cp_fm_type)                                   :: fm_X_ia, fm_Y_ia
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: fm_dipole_ab_trunc, fm_dipole_ai_trunc, &
         fm_dipole_ij_trunc, fm_quadpole_ab_trunc, fm_quadpole_ai_trunc, fm_quadpole_ij_trunc
      TYPE(exciton_descr_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: exc_descr

      CALL timeset(routineN, handle)

      !Prepare variables for printing
      IF (mp2_env%bse%bse_spin_config == 0) THEN
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

      n_moments_di = 3
      n_moments_quad = 9
      ! Compute BSE dipoles and oscillator strengths - Keep in memory for later usage
      ! Need dipoles also for spatial expectation values, which are well-defined also for triplets
      ALLOCATE (fm_dipole_ij_trunc(n_moments_di))
      ALLOCATE (fm_dipole_ab_trunc(n_moments_di))
      ALLOCATE (fm_dipole_ai_trunc(n_moments_di))
      ALLOCATE (ref_point_multipole(3))
      ! Obtain dipoles in MO basis
      CALL get_multipoles_mo(fm_dipole_ai_trunc, fm_dipole_ij_trunc, fm_dipole_ab_trunc, &
                             qs_env, mo_coeff, ref_point_multipole, 1, &
                             homo, virtual, fm_eigvec_X%matrix_struct%context)
      ! Compute exciton descriptors from these multipoles
      IF (mp2_env%bse%num_print_exc_descr > 0) THEN
         ! Obtain quadrupoles in MO basis
         ALLOCATE (fm_quadpole_ij_trunc(n_moments_quad))
         ALLOCATE (fm_quadpole_ab_trunc(n_moments_quad))
         ALLOCATE (fm_quadpole_ai_trunc(n_moments_quad))
         CALL get_multipoles_mo(fm_quadpole_ai_trunc, fm_quadpole_ij_trunc, fm_quadpole_ab_trunc, &
                                qs_env, mo_coeff, ref_point_multipole, 2, &
                                homo, virtual, fm_eigvec_X%matrix_struct%context)
         ! Iterate over excitation index outside of routine to make it compatible with tddft module
         ALLOCATE (exc_descr(mp2_env%bse%num_print_exc_descr))
         DO i_exc = 1, mp2_env%bse%num_print_exc_descr
            CALL reshuffle_eigvec(fm_eigvec_X, fm_X_ia, homo, virtual, i_exc, &
                                  .FALSE., unit_nr, mp2_env)
            IF (.NOT. flag_TDA) THEN
               CALL reshuffle_eigvec(fm_eigvec_Y, fm_Y_ia, homo, virtual, i_exc, &
                                     .FALSE., unit_nr, mp2_env)

               CALL get_exciton_descriptors(exc_descr, fm_X_ia, &
                                            fm_quadpole_ij_trunc, fm_quadpole_ab_trunc, &
                                            fm_quadpole_ai_trunc, &
                                            i_exc, homo, virtual, &
                                            fm_Y_ia)
            ELSE
               CALL get_exciton_descriptors(exc_descr, fm_X_ia, &
                                            fm_quadpole_ij_trunc, fm_quadpole_ab_trunc, &
                                            fm_quadpole_ai_trunc, &
                                            i_exc, homo, virtual)
            END IF
            CALL cp_fm_release(fm_X_ia)
            IF (.NOT. flag_TDA) THEN
               CALL cp_fm_release(fm_Y_ia)
            END IF
         END DO
      END IF

      IF (mp2_env%bse%bse_spin_config == 0) THEN
         CALL get_oscillator_strengths(fm_eigvec_X, Exc_ens, fm_dipole_ai_trunc, &
                                       trans_mom_bse, oscill_str, polarizability_residues, &
                                       mp2_env, homo, virtual, unit_nr, &
                                       fm_eigvec_Y)
      END IF

      ! Prints basic definitions used in BSE calculation
      CALL print_output_header(homo, virtual, homo_irred, flag_TDA, &
                               multiplet, alpha, mp2_env, unit_nr)

      ! Prints excitation energies up to user-specified number
      CALL print_excitation_energies(Exc_ens, homo, virtual, flag_TDA, multiplet, &
                                     info_approximation, mp2_env, unit_nr)

      ! Print single particle transition amplitudes, i.e. components of eigenvectors X and Y
      CALL print_transition_amplitudes(fm_eigvec_X, homo, virtual, homo_irred, &
                                       info_approximation, mp2_env, unit_nr, fm_eigvec_Y)

      ! Prints optical properties, if state is a singlet
      CALL print_optical_properties(Exc_ens, oscill_str, trans_mom_bse, polarizability_residues, &
                                    homo, virtual, homo_irred, flag_TDA, &
                                    info_approximation, mp2_env, unit_nr)
      ! Print exciton descriptors if keyword is invoked
      IF (mp2_env%bse%num_print_exc_descr > 0) THEN
         CALL print_exciton_descriptors(exc_descr, ref_point_multipole, unit_nr, &
                                        mp2_env%bse%num_print_exc_descr, mp2_env%bse%bse_debug_print, &
                                        mp2_env%bse%print_directional_exc_descr, &
                                        'BSE|', qs_env)
      END IF

      ! Compute and print excitation wavefunctions
      IF (mp2_env%bse%do_nto_analysis) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A47)') &
               'BSE|', "Calculating Natural Transition Orbitals (NTOs)."
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
         END IF
         CALL calculate_NTOs(fm_eigvec_X, fm_eigvec_Y, &
                             mo_coeff, homo, virtual, &
                             info_approximation, &
                             oscill_str, &
                             qs_env, unit_nr, mp2_env)
      END IF

      DO idir = 1, n_moments_di
         CALL cp_fm_release(fm_dipole_ai_trunc(idir))
         CALL cp_fm_release(fm_dipole_ij_trunc(idir))
         CALL cp_fm_release(fm_dipole_ab_trunc(idir))
      END DO
      IF (mp2_env%bse%num_print_exc_descr > 0) THEN
         DO idir = 1, n_moments_quad
            CALL cp_fm_release(fm_quadpole_ai_trunc(idir))
            CALL cp_fm_release(fm_quadpole_ij_trunc(idir))
            CALL cp_fm_release(fm_quadpole_ab_trunc(idir))
         END DO
         DEALLOCATE (fm_quadpole_ai_trunc, fm_quadpole_ij_trunc, fm_quadpole_ab_trunc)
         DEALLOCATE (exc_descr)
      END IF
      DEALLOCATE (fm_dipole_ai_trunc, fm_dipole_ij_trunc, fm_dipole_ab_trunc)
      DEALLOCATE (ref_point_multipole)
      IF (mp2_env%bse%bse_spin_config == 0) THEN
         DEALLOCATE (oscill_str, trans_mom_bse, polarizability_residues)
      END IF
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
      END IF

      CALL timestop(handle)

   END SUBROUTINE postprocess_bse

END MODULE bse_full_diag
