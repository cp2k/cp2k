!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate MP2 energy with laplace approach
!> \par History
!>      11.2012 created [Mauro Del Ben]
! **************************************************************************************************
MODULE mp2_laplace
!
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_laplace'

   PUBLIC :: calc_fm_mat_S_laplace, SOS_MP2_postprocessing

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S ...
!> \param first_cycle ...
!> \param homo ...
!> \param virtual ...
!> \param Eigenval ...
!> \param ajquad ...
!> \param alpha_old ...
! **************************************************************************************************
   SUBROUTINE calc_fm_mat_S_laplace(fm_mat_S, first_cycle, homo, virtual, Eigenval, ajquad, alpha_old)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      LOGICAL, INTENT(IN)                                :: first_cycle
      INTEGER, INTENT(IN)                                :: homo, virtual
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: ajquad, alpha_old

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_fm_mat_S_laplace'

      INTEGER                                            :: avirt, handle, i_global, iiB, iocc, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      REAL(KIND=dp)                                      :: laplace_transf

      CALL timeset(routineN, handle)

      ! get info of fm_mat_S
      CALL cp_fm_get_info(matrix=fm_mat_S, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices)

      ! update G matrix with the new value of w and a
      IF (first_cycle) THEN
         DO jjB = 1, ncol_local
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)

               iocc = MAX(1, i_global - 1)/virtual + 1
               avirt = i_global - (iocc - 1)*virtual

               laplace_transf = EXP(Eigenval(iocc)*ajquad)*EXP(-Eigenval(avirt + homo)*ajquad)
               laplace_transf = SQRT(laplace_transf)

               fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)*laplace_transf

            END DO
         END DO
      ELSE
         DO jjB = 1, ncol_local
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)

               iocc = MAX(1, i_global - 1)/virtual + 1
               avirt = i_global - (iocc - 1)*virtual

               laplace_transf = EXP(Eigenval(iocc)*(ajquad - alpha_old))*EXP(-Eigenval(avirt + homo)*(ajquad - alpha_old))
               laplace_transf = SQRT(laplace_transf)

               fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)*laplace_transf

            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE calc_fm_mat_S_laplace

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_Q ...
!> \param Erpa ...
!> \param tau_wjquad ...
! **************************************************************************************************
   SUBROUTINE SOS_MP2_postprocessing(fm_mat_Q, Erpa, tau_wjquad)
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_Q
      REAL(KIND=dp), INTENT(INOUT)                       :: Erpa
      REAL(KIND=dp), INTENT(IN)                          :: tau_wjquad

      CHARACTER(LEN=*), PARAMETER :: routineN = 'SOS_MP2_postprocessing'

      INTEGER                                            :: handle, jjB, ncol_local
      REAL(KIND=dp)                                      :: trace_XX

      CALL timeset(routineN, handle)

      ! get info of fm_mat_Q
      CALL cp_fm_get_info(matrix=fm_mat_Q(1)%matrix, &
                          ncol_local=ncol_local)

      ! calculate the trace of the product Q*Q
      trace_XX = 0.0_dp
      DO jjB = 1, ncol_local
         trace_XX = trace_XX + DOT_PRODUCT(fm_mat_Q(1)%matrix%local_data(:, jjB), &
                                           fm_mat_Q(SIZE(fm_mat_Q))%matrix%local_data(:, jjB))
      END DO

      Erpa = Erpa - trace_XX*tau_wjquad

      CALL timestop(handle)

   END SUBROUTINE SOS_MP2_postprocessing

END MODULE mp2_laplace
