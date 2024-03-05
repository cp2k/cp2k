!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the full diagonalization of GW + Bethe-Salpeter for computing
!> \brief electronic excitations
!> \par History
!>      10.2023 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_full_diag

   USE bse_util,                        ONLY: comp_eigvec_coeff_BSE,&
                                              redist_fm,&
                                              sort_excitations,&
                                              truncate_fm
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
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
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE input_constants,                 ONLY: bse_singlet,&
                                              bse_triplet
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_full_diag'

   PUBLIC :: create_A, diagonalize_A, create_B, create_hermitian_form_of_ABBA, &
             diagonalize_C, truncate_BSE_matrices

CONTAINS

! **************************************************************************************************
!> \brief Create A Matrix: A_iajb = (E_a-E_i) delta_ij delta_ab + alpha * v_iajb - W_ijab
!> \param fm_mat_S ...
!> \param fm_mat_S_bar_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_mat_A ...
!> \param Eigenval ...
!> \param unit_nr ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param mp2_env ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE create_A(fm_mat_S, fm_mat_S_bar_ij_bse, fm_mat_S_ab_bse, &
                       fm_mat_A, Eigenval, unit_nr, &
                       homo, virtual, dimen_RI, mp2_env, &
                       para_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S, fm_mat_S_bar_ij_bse, &
                                                            fm_mat_S_ab_bse
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_A
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
      TYPE(cp_fm_type)                                   :: fm_mat_W

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

      ! create the blacs env for ij matrices (NOT fm_mat_S%matrix_struct related parallel_gemms!)
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

      !We have to use the same blacs_env for A as for the matrices fm_mat_S from RPA
      !Logic: A_iajb = (eps_a-eps_i) delta_ij delta _ab + alpha * v_iajb - W_ijab
      ! We create v_iajb and W_ijab, then we communicate entries from local W_ijab
      ! to the full matrix v_iajb. By adding these and the energy diffenences: v_iajb -> A_iajb
      ! We use the A matrix already from the start instead of v
      CALL cp_fm_struct_create(fm_struct_a, context=fm_mat_S%matrix_struct%context, nrow_global=homo*virtual, &
                               ncol_global=homo*virtual, para_env=fm_mat_S%matrix_struct%para_env)
      CALL cp_fm_create(fm_mat_A, fm_struct_a, name="fm_A_iajb")
      CALL cp_fm_set_all(fm_mat_A, 0.0_dp)

      CALL cp_fm_struct_create(fm_struct_W, context=fm_mat_S_ab_bse%matrix_struct%context, nrow_global=homo**2, &
                               ncol_global=virtual**2, para_env=fm_mat_S_ab_bse%matrix_struct%para_env)
      CALL cp_fm_create(fm_mat_W, fm_struct_W, name="fm_W_ijab")
      CALL cp_fm_set_all(fm_mat_W, 0.0_dp)

      !Create A matrix from GW Energies, v_iajb and W_ijab (different blacs_env!)
      !v_iajb, which is directly initialized in A (with a factor of alpha)
      CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=alpha, &
                         matrix_a=fm_mat_S, matrix_b=fm_mat_S, beta=0.0_dp, &
                         matrix_c=fm_mat_A)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A16)') 'BSE|DEBUG|', 'Allocated A_iajb'
      END IF

      !W_ijab
      CALL parallel_gemm(transa="T", transb="N", m=homo**2, n=virtual**2, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_S_bar_ij_bse, matrix_b=fm_mat_S_ab_bse, beta=0.0_dp, &
                         matrix_c=fm_mat_W)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A16)') 'BSE|DEBUG|', 'Allocated W_ijab'
      END IF

      ! We start by moving data from local parts of W_ijab to the full matrix A_iajb using buffers
      CALL cp_fm_get_info(matrix=fm_mat_A, &
                          nrow_local=nrow_local_A, &
                          ncol_local=ncol_local_A, &
                          row_indices=row_indices_A, &
                          col_indices=col_indices_A)
      ! Writing -1.0_dp * W_ijab to A_iajb, i.e. beta = -1.0_dp,
      ! W_ijab: nrow_secidx_in  = homo,    ncol_secidx_in  = virtual
      ! A_iajb: nrow_secidx_out = virtual, ncol_secidx_out = virtual
      reordering = (/1, 3, 2, 4/)
      CALL redist_fm(fm_mat_A, fm_mat_W, -1.0_dp, homo, virtual, &
                     virtual, virtual, unit_nr, reordering, mp2_env)
      !full matrix W is not needed anymore, release it to save memory
      CALL cp_fm_struct_release(fm_struct_W)
      CALL cp_fm_release(fm_mat_W)
      CALL cp_fm_struct_release(fm_struct_A)

      !Now add the energy differences to A
      DO ii = 1, nrow_local_A

         i_row_global = row_indices_A(ii)

         DO jj = 1, ncol_local_A

            j_col_global = col_indices_A(jj)

            IF (i_row_global == j_col_global) THEN
               i_occ_row = (i_row_global - 1)/virtual + 1
               a_virt_row = MOD(i_row_global - 1, virtual) + 1
               eigen_diff = Eigenval(a_virt_row + homo) - Eigenval(i_occ_row)
               fm_mat_A%local_data(ii, jj) = fm_mat_A%local_data(ii, jj) + eigen_diff

            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE create_A

! **************************************************************************************************
!> \brief Create B matrix for full BSE: B_iajb   =  alpha * v_iajb   -  W_ibaj
!> \param fm_mat_S ...
!> \param fm_mat_S_bar_ia_bse ...
!> \param fm_mat_B ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE create_B(fm_mat_S, fm_mat_S_bar_ia_bse, fm_mat_B, &
                       homo, virtual, dimen_RI, unit_nr, mp2_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S, fm_mat_S_bar_ia_bse
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_B
      INTEGER, INTENT(IN)                                :: homo, virtual, dimen_RI, unit_nr
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_B'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(4)                              :: reordering
      REAL(KIND=dp)                                      :: alpha
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_v
      TYPE(cp_fm_type)                                   :: fm_mat_W

      CALL timeset(routineN, handle)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A10)') 'BSE|DEBUG|', 'Creating B'
      END IF

      !Determines factor of exchange term, depending on requested spin configuration (cf. input_constants.F)
      SELECT CASE (mp2_env%ri_g0w0%bse_spin_config)
      CASE (bse_singlet)
         alpha = 2.0_dp
      CASE (bse_triplet)
         alpha = 0.0_dp
      END SELECT

      CALL cp_fm_struct_create(fm_struct_v, context=fm_mat_S%matrix_struct%context, nrow_global=homo*virtual, &
                               ncol_global=homo*virtual, para_env=fm_mat_S%matrix_struct%para_env)
      CALL cp_fm_create(fm_mat_B, fm_struct_v, name="fm_B_iajb")
      CALL cp_fm_set_all(fm_mat_B, 0.0_dp)

      CALL cp_fm_create(fm_mat_W, fm_struct_v, name="fm_W_ibaj")
      CALL cp_fm_set_all(fm_mat_W, 0.0_dp)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A16)') 'BSE|DEBUG|', 'Allocated B_iajb'
      END IF

      CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=alpha, &
                         matrix_a=fm_mat_S, matrix_b=fm_mat_S, beta=0.0_dp, &
                         matrix_c=fm_mat_B)

      !W_ibaj
      CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_S_bar_ia_bse, matrix_b=fm_mat_S, beta=0.0_dp, &
                         matrix_c=fm_mat_W)
      !from W_iajb to W_ibja(formally: W_ibaj, but our internal indexorder is different)
      reordering = (/1, 4, 3, 2/)
      CALL redist_fm(fm_mat_B, fm_mat_W, -1.0_dp, virtual, virtual, &
                     virtual, virtual, unit_nr, reordering, mp2_env)

      CALL cp_fm_release(fm_mat_W)
      CALL timestop(handle)

   END SUBROUTINE create_B

   ! **************************************************************************************************
!> \brief Create the hermitian form of the ABBA matrix, i.e. C = (A-B)^-0.5 (A+B) (A-B)^-0.5
!> \param fm_mat_A ...
!> \param fm_mat_B ...
!> \param fm_mat_C ...
!> \param fm_mat_work_diff_raised_pos ...
!> \param fm_mat_work_diff_raised_neg ...
!> \param homo ...
!> \param virtual ...
!> \param unit_nr ...
!> \param mp2_env ...
!> \param diag_est ...
! **************************************************************************************************
   SUBROUTINE create_hermitian_form_of_ABBA(fm_mat_A, fm_mat_B, fm_mat_C, &
                                            fm_mat_work_diff_raised_pos, fm_mat_work_diff_raised_neg, &
                                            homo, virtual, unit_nr, mp2_env, diag_est)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_A, fm_mat_B
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_C, fm_mat_work_diff_raised_pos, &
                                                            fm_mat_work_diff_raised_neg
      INTEGER, INTENT(IN)                                :: homo, virtual, unit_nr
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      REAL(KIND=dp), INTENT(IN)                          :: diag_est

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_hermitian_form_of_ABBA'

      INTEGER                                            :: dim_mat, handle, n_dependent
      REAL(KIND=dp), DIMENSION(2)                        :: eigvals_AB_diff
      TYPE(cp_fm_type)                                   :: fm_mat_work_diff, &
                                                            fm_mat_work_diff_eigvec, &
                                                            fm_mat_work_product, fm_mat_work_sum

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A25,A39,ES6.0,A3)') 'BSE|', 'Diagonalizing aux. matrix', &
            ' with size of A. This will take around ', diag_est, " s."
      END IF

      ! Create work matrices, which will hold A+B and A-B and their powers
      ! C is created afterwards to save memory
      CALL cp_fm_create(fm_mat_work_sum, fm_mat_A%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_A, fm_mat_work_sum)
      CALL cp_fm_create(fm_mat_work_diff, fm_mat_A%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_A, fm_mat_work_diff)
      CALL cp_fm_create(fm_mat_work_diff_raised_pos, fm_mat_A%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work_diff_raised_pos, 0.0_dp)
      CALL cp_fm_create(fm_mat_work_diff_raised_neg, fm_mat_A%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work_diff_raised_neg, 0.0_dp)

      CALL cp_fm_create(fm_mat_work_diff_eigvec, fm_mat_A%matrix_struct)
      CALL cp_fm_create(fm_mat_work_product, fm_mat_A%matrix_struct)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A19)') 'BSE|DEBUG|', 'Created work arrays'
      END IF

      !Add/Substract B
      CALL cp_fm_scale_and_add(1.0_dp, fm_mat_work_sum, 1.0_dp, fm_mat_B)
      CALL cp_fm_scale_and_add(1.0_dp, fm_mat_work_diff, -1.0_dp, fm_mat_B)

      !cp_fm_power will overwrite matrix, therefore we copy it
      CALL cp_fm_to_fm(fm_mat_work_diff, fm_mat_work_diff_raised_neg)
      CALL cp_fm_to_fm(fm_mat_work_diff, fm_mat_work_diff_raised_pos)

      ! In the output, we want to print singleparticle transitions of ABBA, where we need
      ! (A-B)**-0.5 as well as (A-B)**0.5. In order to avoid a second diagonalization (cp_fm_power),
      ! we create (A-B)**0.5 from (A-B)**-0.5 by multiplication with (A-B).

      !Raise A-B to 1/2, no quenching of eigenvectors, hence threshold=0.0_dp
      CALL cp_fm_power(fm_mat_work_diff_raised_neg, fm_mat_work_diff_eigvec, -0.5_dp, 0.0_dp, n_dependent, eigvals=eigvals_AB_diff)
      !Raise an error in case the the matrix A-B is not positive definite (i.e. negative eigenvalues)
      !In this case, the procedure for hermitian form of ABBA is not applicable
      IF (eigvals_AB_diff(1) < 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "Matrix (A-B) is not positive definite. "// &
                       "Hermitian diagonalization of full BSE matrix is ill-defined.")
      END IF
      !(A-B)**-0.5 is written to the eigvec - work matrix
      CALL cp_fm_release(fm_mat_work_diff_eigvec)

      ! We keep diff_eigvec for print of singleparticle transitions of ABBA
      ! We further create (A-B)**0.5 for the singleparticle transitions of ABBA
      dim_mat = homo*virtual
      CALL parallel_gemm("N", "N", dim_mat, dim_mat, dim_mat, 1.0_dp, fm_mat_work_diff_raised_neg, fm_mat_work_diff, 0.0_dp, &
                         fm_mat_work_diff_raised_pos)

      CALL cp_fm_release(fm_mat_work_diff)

      !Create C = (A-B)^0.5 (A+B) (A-B)^0.5
      CALL parallel_gemm("N", "N", dim_mat, dim_mat, dim_mat, 1.0_dp, fm_mat_work_diff_raised_pos, fm_mat_work_sum, 0.0_dp, &
                         fm_mat_work_product)
      CALL cp_fm_release(fm_mat_work_sum)

      !Now create C itself
      CALL cp_fm_create(fm_mat_C, fm_mat_A%matrix_struct)
      CALL cp_fm_set_all(fm_mat_C, 0.0_dp)

      CALL parallel_gemm("N", "N", dim_mat, dim_mat, dim_mat, 1.0_dp, fm_mat_work_product, fm_mat_work_diff_raised_pos, 0.0_dp, &
                         fm_mat_C)
      CALL cp_fm_release(fm_mat_work_product)

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A10,T13,A36)') 'BSE|DEBUG|', 'Filled C=(A-B)^-0.5 (A+B) (A-B)^-0.5'
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Diagonalize the hermitian form of the ABBA matrix, i.e. C = (A-B)^-0.5 (A+B) (A-B)^-0.5
!> \param fm_mat_C ...
!> \param homo ...
!> \param virtual ...
!> \param homo_irred ...
!> \param virtual_irred ...
!> \param fm_mat_work_diff_raised_pos ...
!> \param fm_mat_work_diff_raised_neg ...
!> \param unit_nr ...
!> \param diag_est ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE diagonalize_C(fm_mat_C, homo, virtual, homo_irred, virtual_irred, &
                            fm_mat_work_diff_raised_pos, fm_mat_work_diff_raised_neg, &
                            unit_nr, diag_est, mp2_env)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_C
      INTEGER, INTENT(IN)                                :: homo, virtual, homo_irred, virtual_irred
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_work_diff_raised_pos, &
                                                            fm_mat_work_diff_raised_neg
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), INTENT(IN)                          :: diag_est
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'diagonalize_C'

      INTEGER                                            :: diag_info, handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens
      TYPE(cp_fm_type)                                   :: fm_mat_eigvec, fm_mat_eigvec_transform, &
                                                            fm_mat_eigvec_transform_neg

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A17,A22,ES6.0,A3)') 'BSE|', 'Diagonalizing C. ', &
            'This will take around ', diag_est, ' s.'
      END IF

      !We have now the full matrix C=(A-B)^-0.5 (A+B) (A-B)^-0.5
      !Now: Diagonalize it
      CALL cp_fm_create(fm_mat_eigvec, fm_mat_C%matrix_struct)

      ALLOCATE (Exc_ens(homo*virtual))

      CALL choose_eigv_solver(fm_mat_C, fm_mat_eigvec, Exc_ens, diag_info)
      CPASSERT(diag_info == 0)
      Exc_ens = SQRT(Exc_ens)

      ! Prepare eigenvector for interpretation of singleparticle transitions
      ! Compare: F. Furche J. Chem. Phys., Vol. 114, No. 14, (2001)
      ! We aim for the upper part of the vector (X,Y) for a direct comparison with the TDA result

      ! Following Furche, we basically use Eqs. (A10): First, we multiply
      ! the (A-B)^+-0.5 with eigenvectors and then the eigenvalues
      ! One has to be careful about the index structure, since the eigenvector matrix is not symmetric anymore!

      ! First, Eq. I from (A10) from Furche: (X+Y)_n = (Ω_n)^0.5 (A-B)^-0.5 T_n
      CALL cp_fm_create(fm_mat_eigvec_transform, fm_mat_C%matrix_struct)
      CALL cp_fm_set_all(fm_mat_eigvec_transform, 0.0_dp)
      CALL parallel_gemm(transa="N", transb="N", m=homo*virtual, n=homo*virtual, k=homo*virtual, alpha=1.0_dp, &
                         matrix_a=fm_mat_work_diff_raised_pos, matrix_b=fm_mat_eigvec, beta=0.0_dp, &
                         matrix_c=fm_mat_eigvec_transform)
      CALL cp_fm_release(fm_mat_work_diff_raised_pos)
      CALL comp_eigvec_coeff_BSE(fm_mat_eigvec_transform, Exc_ens, -0.5_dp, gamma=2.0_dp, do_transpose=.TRUE.)

      ! Second, Eq. II from (A10) from Furche: (X-Y)_n = (Ω_n)^0.5 (A-B)^-0.5 T_n
      CALL cp_fm_create(fm_mat_eigvec_transform_neg, fm_mat_C%matrix_struct)
      CALL cp_fm_set_all(fm_mat_eigvec_transform_neg, 0.0_dp)
      CALL parallel_gemm(transa="N", transb="N", m=homo*virtual, n=homo*virtual, k=homo*virtual, alpha=1.0_dp, &
                         matrix_a=fm_mat_work_diff_raised_neg, matrix_b=fm_mat_eigvec, beta=0.0_dp, &
                         matrix_c=fm_mat_eigvec_transform_neg)
      CALL cp_fm_release(fm_mat_work_diff_raised_neg)
      CALL cp_fm_release(fm_mat_eigvec)
      CALL comp_eigvec_coeff_BSE(fm_mat_eigvec_transform_neg, Exc_ens, 0.5_dp, gamma=2.0_dp, do_transpose=.TRUE.)

      ! Now, we add the two equations to obtain X_n
      CALL cp_fm_scale_and_add(1.0_dp, fm_mat_eigvec_transform, 1.0_dp, fm_mat_eigvec_transform_neg)

      !Cleanup
      CALL cp_fm_release(fm_mat_eigvec_transform_neg)

      CALL success_message(Exc_ens, fm_mat_eigvec_transform, mp2_env, &
                           homo, virtual, homo_irred, virtual_irred, unit_nr, .FALSE.)

      DEALLOCATE (Exc_ens)
      CALL cp_fm_release(fm_mat_eigvec_transform)

      CALL timestop(handle)

   END SUBROUTINE diagonalize_C

! **************************************************************************************************
!> \brief Diagonalizes the A matrix, i.e. TDA of BSE
!> \param fm_mat_A ...
!> \param homo ...
!> \param virtual ...
!> \param homo_irred ...
!> \param virtual_irred ...
!> \param unit_nr ...
!> \param diag_est ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE diagonalize_A(fm_mat_A, homo, virtual, homo_irred, virtual_irred, &
                            unit_nr, diag_est, mp2_env)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_A
      INTEGER, INTENT(IN)                                :: homo, virtual, homo_irred, &
                                                            virtual_irred, unit_nr
      REAL(KIND=dp), INTENT(IN)                          :: diag_est
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'diagonalize_A'

      INTEGER                                            :: diag_info, handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens
      TYPE(cp_fm_type)                                   :: fm_mat_eigvec

      CALL timeset(routineN, handle)

      !Continue with formatting of subroutine create_A
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A17,A22,ES6.0,A3)') 'BSE|', 'Diagonalizing A. ', &
            'This will take around ', diag_est, ' s.'
      END IF

      !We have now the full matrix A_iajb, distributed over all ranks
      !Now: Diagonalize it
      CALL cp_fm_create(fm_mat_eigvec, fm_mat_A%matrix_struct)

      ALLOCATE (Exc_ens(homo*virtual))

      CALL choose_eigv_solver(fm_mat_A, fm_mat_eigvec, Exc_ens, diag_info)
      CPASSERT(diag_info == 0)
      CALL success_message(Exc_ens, fm_mat_eigvec, mp2_env, &
                           homo, virtual, homo_irred, virtual_irred, unit_nr, .TRUE.)

      CALL cp_fm_release(fm_mat_eigvec)
      DEALLOCATE (Exc_ens)

      CALL timestop(handle)

   END SUBROUTINE diagonalize_A

! **************************************************************************************************
!> \brief Prints the success message (incl. energies) for full diag of BSE (TDA/full ABBA via flag)
!> \param Exc_ens ...
!> \param fm_mat_eigvec ...
!> \param mp2_env ...
!> \param homo ...
!> \param virtual ...
!> \param homo_irred ...
!> \param virtual_irred ...
!> \param unit_nr ...
!> \param flag_TDA ...
! **************************************************************************************************
   SUBROUTINE success_message(Exc_ens, fm_mat_eigvec, mp2_env, &
                              homo, virtual, homo_irred, virtual_irred, unit_nr, flag_TDA)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_eigvec
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      INTEGER                                            :: homo, virtual, homo_irred, &
                                                            virtual_irred, unit_nr
      LOGICAL, OPTIONAL                                  :: flag_TDA

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'success_message'

      CHARACTER(LEN=10)                                  :: info_approximation, multiplet
      INTEGER                                            :: handle, i_exc, k, num_entries
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: idx_homo, idx_virt
      REAL(KIND=dp)                                      :: alpha
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigvec_entries

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
         info_approximation = "-full-"
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
            WRITE (unit_nr, '(T2,A4,T7,A65)') 'BSE|', 'prelim.: cf. Eq. (XX) in  Ref. XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '**************************************************************************'
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '*          Full Bethe Salpeter equation (BSE) (i.e. without TDA)         *'
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
            WRITE (unit_nr, '(T2,A4,T7,A65)') 'BSE|', 'prelim.: cf. Eq. (XX) in  Ref. XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
         END IF
         IF (.NOT. flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A60)') 'BSE|', 'The BSE is solved for Ω^n and X_ia^n as a hermitian problem,'
            WRITE (unit_nr, '(T2,A4,T7,A70)') 'BSE|', 'e.g., Eq. (A8) in F. Furche, J. Chem. Phys., Vol. 114, No. 14, (2001).'
            ! WRITE (unit_nr, '(T2,A4,T7,A69)') 'BSE|', 'C_ia,jb = sum_kc,ld ((A-B)^-0.5)_ia,kc (A+B)_kc,ld ((A-B)^-0.5)_ld,jb'
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
         WRITE (unit_nr, '(T2,A4,T7,A47,A7,A9,F3.1))') 'BSE|', &
            'The spin-dependent factor is for the requested ', multiplet, " is α = ", alpha
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         IF (flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4,T7,A56)') 'BSE|', 'Excitation energies from solving the BSE within the TDA:'
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A57)') 'BSE|', 'Excitation energies from solving the BSE without the TDA:'
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T13,A10,T27,A13,T42,A12,T59,A22)') 'BSE|', &
            'Excitation', "Multiplet", 'TDA/full BSE', 'Excitation energy (eV)'
      END IF
      !prints actual energies values
      IF (unit_nr > 0) THEN
         DO i_exc = 1, MIN(homo*virtual, mp2_env%ri_g0w0%num_print_exc)
            WRITE (unit_nr, '(T2,A4,T7,I16,T27,A7,A6,T48,A6,T59,F22.4)') &
               'BSE|', i_exc, multiplet, " State", info_approximation, Exc_ens(i_exc)*evolt
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
         WRITE (unit_nr, '(T2,A4,T7,A9,T20,A1,T22,A2,T29,A1,T42,A12,T71,A10)') &
            "BSE|", "n-th exc.", "i", "=>", "a", 'TDA/full BSE', "|X_ia^n|"
      END IF
      DO i_exc = 1, MIN(homo*virtual, mp2_env%ri_g0w0%num_print_exc)
         !Iterate through eigenvector and print values above threshold
         CALL filter_eigvec_contrib(fm_mat_eigvec, idx_homo, idx_virt, eigvec_entries, &
                                    i_exc, homo, virtual, num_entries, mp2_env)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            DO k = 1, num_entries
               WRITE (unit_nr, '(T2,A4,T11,I5,T16,I5,T22,A2,T25,I5,T48,A6,T59,F22.4)') &
                  "BSE|", i_exc, homo_irred - homo + idx_homo(k), "=>", &
                  homo_irred + idx_virt(k), info_approximation, ABS(eigvec_entries(k))
            END DO
         END IF

         DEALLOCATE (idx_homo)
         DEALLOCATE (idx_virt)
         DEALLOCATE (eigvec_entries)
      END DO
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A53)') 'BSE|', 'The BSE was successfully calculated. Have a nice day!'
      END IF

      CALL timestop(handle)
   END SUBROUTINE success_message

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_eigvec ...
!> \param idx_homo ...
!> \param idx_virt ...
!> \param eigvec_entries ...
!> \param i_exc ...
!> \param homo ...
!> \param virtual ...
!> \param num_entries ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE filter_eigvec_contrib(fm_mat_eigvec, idx_homo, idx_virt, eigvec_entries, &
                                    i_exc, homo, virtual, num_entries, mp2_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_eigvec
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: idx_homo, idx_virt
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: eigvec_entries
      INTEGER                                            :: i_exc, homo, virtual, num_entries
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'filter_eigvec_contrib'

      INTEGER                                            :: eigvec_idx, handle, ii, iproc, jj, kk, &
                                                            ncol_local, nrow_local, &
                                                            num_entries_local
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_entries_to_comm
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: eigvec_entry
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_entries
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      para_env => fm_mat_eigvec%matrix_struct%para_env

      CALL cp_fm_get_info(matrix=fm_mat_eigvec, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ALLOCATE (num_entries_to_comm(0:para_env%num_pe - 1))
      num_entries_to_comm(:) = 0

      DO jj = 1, ncol_local
         !First check if i is localized on this proc
         IF (col_indices(jj) /= i_exc) THEN
            CYCLE
         END IF
         DO ii = 1, nrow_local
            eigvec_idx = row_indices(ii)
            eigvec_entry = fm_mat_eigvec%local_data(ii, jj)/SQRT(2.0_dp)
            IF (ABS(eigvec_entry) > mp2_env%ri_g0w0%eps_x) THEN
               num_entries_to_comm(para_env%mepos) = num_entries_to_comm(para_env%mepos) + 1
            END IF
         END DO
      END DO

      !Gather number of entries of other processes
      CALL para_env%sum(num_entries_to_comm)

      num_entries_local = num_entries_to_comm(para_env%mepos)

      ALLOCATE (buffer_entries(0:para_env%num_pe - 1))

      DO iproc = 0, para_env%num_pe - 1
         ALLOCATE (buffer_entries(iproc)%msg(num_entries_to_comm(iproc)))
         ALLOCATE (buffer_entries(iproc)%indx(num_entries_to_comm(iproc), 2))
         buffer_entries(iproc)%msg = 0.0_dp
         buffer_entries(iproc)%indx = 0
      END DO

      kk = 1
      DO jj = 1, ncol_local
         !First check if i is localized on this proc
         IF (col_indices(jj) /= i_exc) THEN
            CYCLE
         END IF
         DO ii = 1, nrow_local
            eigvec_idx = row_indices(ii)
            eigvec_entry = fm_mat_eigvec%local_data(ii, jj)/SQRT(2.0_dp)
            IF (ABS(eigvec_entry) > mp2_env%ri_g0w0%eps_x) THEN
               buffer_entries(para_env%mepos)%indx(kk, 1) = (eigvec_idx - 1)/virtual + 1
               buffer_entries(para_env%mepos)%indx(kk, 2) = MOD(eigvec_idx - 1, virtual) + 1
               buffer_entries(para_env%mepos)%msg(kk) = eigvec_entry
               kk = kk + 1
            END IF
         END DO
      END DO

      DO iproc = 0, para_env%num_pe - 1
         CALL para_env%sum(buffer_entries(iproc)%msg)
         CALL para_env%sum(buffer_entries(iproc)%indx)
      END DO

      !Now sum up gathered information
      num_entries = SUM(num_entries_to_comm)
      ALLOCATE (idx_homo(num_entries))
      ALLOCATE (idx_virt(num_entries))
      ALLOCATE (eigvec_entries(num_entries))

      kk = 1
      DO iproc = 0, para_env%num_pe - 1
         IF (num_entries_to_comm(iproc) /= 0) THEN
            DO ii = 1, num_entries_to_comm(iproc)
               idx_homo(kk) = buffer_entries(iproc)%indx(ii, 1)
               idx_virt(kk) = buffer_entries(iproc)%indx(ii, 2)
               eigvec_entries(kk) = buffer_entries(iproc)%msg(ii)
               kk = kk + 1
            END DO
         END IF
      END DO

      !Deallocate all the used arrays
      DO iproc = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_entries(iproc)%msg)
         DEALLOCATE (buffer_entries(iproc)%indx)
      END DO
      DEALLOCATE (buffer_entries)
      DEALLOCATE (num_entries_to_comm)
      NULLIFY (row_indices)
      NULLIFY (col_indices)

      !Now sort the results according to the involved singleparticle orbitals
      ! (homo first, then virtual)
      CALL sort_excitations(idx_homo, idx_virt, eigvec_entries)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S ...
!> \param fm_mat_S_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_mat_S_trunc ...
!> \param fm_mat_S_ij_trunc ...
!> \param fm_mat_S_ab_trunc ...
!> \param Eigenval ...
!> \param Eigenval_reduced ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param homo_red ...
!> \param virt_red ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE truncate_BSE_matrices(fm_mat_S, fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                    fm_mat_S_trunc, fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, &
                                    Eigenval, Eigenval_reduced, &
                                    homo, virtual, dimen_RI, unit_nr, &
                                    homo_red, virt_red, &
                                    mp2_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S, fm_mat_S_ij_bse, &
                                                            fm_mat_S_ab_bse
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_S_trunc, fm_mat_S_ij_trunc, &
                                                            fm_mat_S_ab_trunc
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Eigenval_reduced
      INTEGER, INTENT(IN)                                :: homo, virtual, dimen_RI, unit_nr
      INTEGER, INTENT(OUT)                               :: homo_red, virt_red
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'truncate_BSE_matrices'

      INTEGER                                            :: handle, homo_incl, i_homo, j_virt, &
                                                            virt_incl
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_ab, fm_struct_ia, fm_struct_ij
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      ! Determine index in homo and virtual for truncation
      ! Uses indices of outermost orbitals within energy range (-mp2_env%ri_g0w0%bse_cutoff_occ,mp2_env%ri_g0w0%bse_cutoff_virt)
      IF (mp2_env%ri_g0w0%bse_cutoff_occ > 0 .OR. mp2_env%ri_g0w0%bse_cutoff_virt > 0) THEN
         IF (-mp2_env%ri_g0w0%bse_cutoff_occ .LT. Eigenval(1) .OR. mp2_env%ri_g0w0%bse_cutoff_occ < 0) THEN
            homo_red = homo
            homo_incl = 1
         ELSE
            homo_incl = 1
            DO i_homo = 1, homo
               IF (Eigenval(i_homo) .GT. -mp2_env%ri_g0w0%bse_cutoff_occ) THEN
                  homo_incl = i_homo
                  EXIT
               END IF
            END DO
            homo_red = homo - homo_incl + 1
         END IF

         IF (mp2_env%ri_g0w0%bse_cutoff_virt .GT. Eigenval(homo + virtual) .OR. mp2_env%ri_g0w0%bse_cutoff_virt < 0) THEN
            virt_red = virtual
            virt_incl = virtual
         ELSE
            virt_incl = homo + 1
            DO j_virt = 1, virtual
               IF (Eigenval(homo + j_virt) .GT. mp2_env%ri_g0w0%bse_cutoff_virt) THEN
                  virt_incl = j_virt - 1
                  EXIT
               END IF
            END DO
            virt_red = virt_incl
         END IF
      ELSE
         homo_red = homo
         virt_red = virtual
         homo_incl = 1
         virt_incl = virtual
      END IF
      IF (unit_nr > 0) THEN
         IF (mp2_env%ri_g0w0%bse_cutoff_occ > 0) THEN
            WRITE (unit_nr, '(T2,A4,T7,A29,T71,F10.3)') 'BSE|', 'Cutoff occupied orbitals [eV]', &
               mp2_env%ri_g0w0%bse_cutoff_occ*evolt
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A37)') 'BSE|', 'No cutoff given for occupied orbitals'
         END IF
         IF (mp2_env%ri_g0w0%bse_cutoff_virt > 0) THEN
            WRITE (unit_nr, '(T2,A4,T7,A28,T71,F10.3)') 'BSE|', 'Cutoff virtual orbitals [eV]', &
               mp2_env%ri_g0w0%bse_cutoff_virt*evolt
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A36)') 'BSE|', 'No cutoff given for virtual orbitals'
         END IF
         WRITE (unit_nr, '(T2,A4,T7,A20,T71,I10)') 'BSE|', 'First occupied index', homo_incl
         WRITE (unit_nr, '(T2,A4,T7,A18,T71,I10)') 'BSE|', 'Last virtual index', virt_incl
         WRITE (unit_nr, '(T2,A4,T7,A26,T71,F10.3)') 'BSE|', 'First occupied energy [eV]', Eigenval(homo_incl)*evolt
         WRITE (unit_nr, '(T2,A4,T7,A24,T71,F10.3)') 'BSE|', 'Last virtual energy [eV]', Eigenval(homo + virt_incl)*evolt
         WRITE (unit_nr, '(T2,A4,T7,A35,T71,I10)') 'BSE|', 'Number of GW-corrected occupied MOs', mp2_env%ri_g0w0%corr_mos_occ
         WRITE (unit_nr, '(T2,A4,T7,A34,T71,I10)') 'BSE|', 'Number of GW-corrected virtual MOs', mp2_env%ri_g0w0%corr_mos_virt
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
      END IF
      IF (unit_nr > 0) THEN
         IF (homo - homo_incl + 1 > mp2_env%ri_g0w0%corr_mos_occ) THEN
            CPABORT("Number of GW-corrected occupied MOs too small for chosen BSE cutoff")
         END IF
         IF (virt_incl > mp2_env%ri_g0w0%corr_mos_virt) THEN
            CPABORT("Number of GW-corrected virtual MOs too small for chosen BSE cutoff")
         END IF
      END IF
      !Truncate full fm_S matrices
      !Allocate new truncated matrices of proper size
      para_env => fm_mat_S%matrix_struct%para_env
      context => fm_mat_S%matrix_struct%context

      CALL cp_fm_struct_create(fm_struct_ia, para_env, context, dimen_RI, homo_red*virt_red)
      CALL cp_fm_struct_create(fm_struct_ij, para_env, context, dimen_RI, homo_red*homo_red)
      CALL cp_fm_struct_create(fm_struct_ab, para_env, context, dimen_RI, virt_red*virt_red)

      CALL cp_fm_create(fm_mat_S_trunc, fm_struct_ia, "fm_S_trunc")
      CALL cp_fm_create(fm_mat_S_ij_trunc, fm_struct_ij, "fm_S_ij_trunc")
      CALL cp_fm_create(fm_mat_S_ab_trunc, fm_struct_ab, "fm_S_ab_trunc")

      !Copy parts of original matrices to truncated ones
      IF (mp2_env%ri_g0w0%bse_cutoff_occ > 0 .OR. mp2_env%ri_g0w0%bse_cutoff_virt > 0) THEN
         !Truncate eigenvals
         ALLOCATE (Eigenval_reduced(homo_red + virt_red))
         Eigenval_reduced(:) = Eigenval(homo_incl:homo + virt_incl)

         CALL truncate_fm(fm_mat_S_trunc, fm_mat_S, homo, virtual, &
                          homo_red, virt_red, unit_nr, mp2_env, &
                          nrow_offset=homo_incl)
         CALL truncate_fm(fm_mat_S_ij_trunc, fm_mat_S_ij_bse, homo, homo, &
                          homo_red, homo_red, unit_nr, mp2_env, &
                          homo_incl, homo_incl)
         CALL truncate_fm(fm_mat_S_ab_trunc, fm_mat_S_ab_bse, virtual, virtual, &
                          virt_red, virt_red, unit_nr, mp2_env)

      ELSE
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A4,T7,A37)') 'BSE|', 'No truncation of BSE matrices applied'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
         END IF
         Eigenval_reduced = Eigenval
         CALL cp_fm_to_fm_submat_general(fm_mat_S, fm_mat_S_trunc, dimen_RI, homo_red*virt_red, &
                                         1, 1, 1, 1, context)
         CALL cp_fm_to_fm_submat_general(fm_mat_S_ij_bse, fm_mat_S_ij_trunc, dimen_RI, homo_red*homo_red, &
                                         1, 1, 1, 1, context)
         CALL cp_fm_to_fm_submat_general(fm_mat_S_ab_bse, fm_mat_S_ab_trunc, dimen_RI, virt_red*virt_red, &
                                         1, 1, 1, 1, context)
      END IF

      CALL timestop(handle)

   END SUBROUTINE truncate_BSE_matrices

END MODULE bse_full_diag
