!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief evaluations of colvar for electronic properties
!> \par History
!>      09.2024 created [Rangsiman Ketkaew]
!> \author Rangsiman Ketkaew - University of Zurich (2024)
! **************************************************************************************************
MODULE colvar_elec_utils
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_linres_types,                 ONLY: dcdr_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: mult_fm_sm_fm, &
             mult_fm_fm_fm, &
             charge_cv_linres_init, &
             charge_cv_cleanup

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colvar_elec_utils'

CONTAINS

! **************************************************************************************************
!> \brief Matrix multiplication of fm, sm, and fm for the first and second term in dqdr contribution
!> \param mat_A Full matrix
!> \param mat_B DBCSR sparse matrix
!> \param mat_C Full matrix
!> \param mat_tmp Temporary full matrix
!> \param product Product matrix
!> \param nao Number of AOs
!> \param nmo Number of MOs
!> \date  09.2024
!> \author Rangsiman Ketkaew
! **************************************************************************************************

   SUBROUTINE mult_fm_sm_fm(mat_A, mat_B, mat_C, mat_tmp, product, nao, nmo)

      TYPE(cp_fm_type)                                   :: mat_A
      TYPE(dbcsr_type)                                   :: mat_B
      TYPE(cp_fm_type)                                   :: mat_C, mat_tmp, product
      INTEGER                                            :: nao, nmo

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mult_fm_sm_fm'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! A*B*C = A*(B*C)
      !       = A*(tmp)

      CALL cp_dbcsr_sm_fm_multiply(mat_B, mat_C, &
                                   mat_tmp, nmo)
      CALL parallel_gemm("T", "N", nmo, nmo, nao, &
                         1.0_dp, mat_A, mat_tmp, &
                         0.0_dp, product)

      CALL timestop(handle)

   END SUBROUTINE mult_fm_sm_fm

! **************************************************************************************************
!> \brief Matrix multiplication of fm, fm, and fm for the third term in dqdr contribution
!> \param mat_A Full matrix
!> \param mat_B Full matrix
!> \param mat_C Full matrix
!> \param mat_tmp Temporary full matrix
!> \param product Product matrix
!> \param nao Number of AOs
!> \param nmo Number of MOs
!> \date  09.2024
!> \author Rangsiman Ketkaew
! **************************************************************************************************

   SUBROUTINE mult_fm_fm_fm(mat_A, mat_B, mat_C, mat_tmp, product, nao, nmo)

      TYPE(cp_fm_type)                                   :: mat_A, mat_B, mat_C, mat_tmp, product
      INTEGER                                            :: nao, nmo

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mult_fm_fm_fm'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! A*B*C = A*(B*C)
      !       = A*(tmp)

      CALL parallel_gemm("N", "N", nao, nmo, nao, &
                         1.0_dp, mat_B, mat_C, &
                         0.0_dp, mat_tmp)
      CALL parallel_gemm("T", "N", nmo, nmo, nao, &
                         1.0_dp, mat_A, mat_tmp, &
                         0.0_dp, product)

      ! We can (probably) also convert fm to sm and use mult_fm_sm_fm:
      !  CALL dbcsr_set(dsdr_dbcsr, 0.0_dp)
      !  CALL copy_fm_to_dbcsr(dsdr(i,ii)%matrix, dsdr_dbcsr(i,ii)%matrix)
      !  CALL cp_dbcsr_sm_fm_multiply(dsdr_dbcsr(i,ii)%matrix, mo_coeff, tmp, nmo)

      CALL timestop(handle)

   END SUBROUTINE mult_fm_fm_fm

! **************************************************************************************************
!> \brief Prepare matrices for storing in linres
!> \param qs_env Quickstep environment
!> \param dcdr_env ...
!> \date  09.2024
!> \author Rangsiman Ketkaew
! **************************************************************************************************

   SUBROUTINE charge_cv_linres_init(qs_env, dcdr_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dcdr_env_type)                                :: dcdr_env

      INTEGER                                            :: ispin, m, n, nspins
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control)

      nspins = dft_control%nspins

      ! Store the perturbed wavefunction (dCdR) (no. of atoms x 3 directions x no. of spins)
      ALLOCATE (qs_env%linres_control%dcdr(SIZE(dcdr_env%list_of_atoms), 3, nspins))
      CPASSERT(ASSOCIATED(dcdr_env%likemos_fm_struct))

      ! Create full matrices for derivatives for all atoms
      DO m = 1, SIZE(dcdr_env%list_of_atoms)
         DO n = 1, 3
            ! dCdR
            DO ispin = 1, nspins
               CALL cp_fm_create(matrix=qs_env%linres_control%dcdr(m, n, ispin), &
                                 matrix_struct=dcdr_env%likemos_fm_struct(ispin)%struct, &
                                 name="dCdR MATRIX")
               CALL cp_fm_set_all(qs_env%linres_control%dcdr(m, n, ispin), 0.0_dp)
            END DO
            ! dSdR
            ! CALL cp_fm_create(matrix=qs_env%linres_control%matrix_s1(m, n), &
            !                   matrix_struct=dcdr_env%nsgf_fm_struct, &
            !                   name="dSdR MATRIX")
            ! CALL cp_fm_set_all(qs_env%linres_control%matrix_s1(m, n), 0.0_dp)
         END DO
      END DO

      ! Store overlap matrix derivatives (dSdR) (nsgf x nsgf)
      ! where nsgf = contracted spherical Gaussian basis functions
      NULLIFY (qs_env%linres_control%matrix_s1)
      CALL dbcsr_allocate_matrix_set(qs_env%linres_control%matrix_s1, 3)

      DO n = 1, 3
         CALL dbcsr_init_p(qs_env%linres_control%matrix_s1(n)%matrix)
         CALL dbcsr_create(qs_env%linres_control%matrix_s1(n)%matrix, &
                           template=dcdr_env%matrix_s1(1 + n)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
      END DO

   END SUBROUTINE charge_cv_linres_init

! **************************************************************************************************
!> \brief Clean the charge colvar from the Quickstep environment
!> \param qs_env Quickstep environment
!> \date  09.2024
!> \author Rangsiman Ketkaew
! **************************************************************************************************
   SUBROUTINE charge_cv_cleanup(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: i, ii, iii

      IF (ASSOCIATED(qs_env%linres_control%gapt_charge)) THEN
         DEALLOCATE (qs_env%linres_control%gapt_charge)
      END IF
      IF (ASSOCIATED(qs_env%linres_control%dcdr)) THEN
         DO iii = 1, SIZE(qs_env%linres_control%dcdr, 3)
            DO ii = 1, SIZE(qs_env%linres_control%dcdr, 2)
               DO i = 1, SIZE(qs_env%linres_control%dcdr, 1)
                  CALL cp_fm_release(qs_env%linres_control%dcdr(i, ii, iii))
               END DO
            END DO
         END DO
         DEALLOCATE (qs_env%linres_control%dcdr)
      END IF
      IF (ASSOCIATED(qs_env%linres_control%matrix_s1)) THEN
         ! DO ii = 1, SIZE(qs_env%linres_control%matrix_s1, 2)
         !    DO i = 1, SIZE(qs_env%linres_control%matrix_s1, 1)
         !       CALL cp_fm_release(qs_env%linres_control%matrix_s1(i, ii))
         !    END DO
         ! END DO
         ! DEALLOCATE (qs_env%linres_control%matrix_s1)
         CALL dbcsr_deallocate_matrix_set(qs_env%linres_control%matrix_s1)
      END IF

   END SUBROUTINE charge_cv_cleanup

END MODULE colvar_elec_utils
