!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of commutator [H,r] matrices
!> \par History
!>      JGH: [7.2016]
!> \author Juerg Hutter
! **************************************************************************************************
MODULE qs_commutators
   USE commutator_rkinetic,             ONLY: build_com_tr_matrix
   USE commutator_rpnl,                 ONLY: build_com_rpnl
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE dbcsr_api,                       ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_commutators'

! *** Public subroutines ***

   PUBLIC :: build_com_hr_matrix

CONTAINS

! **************************************************************************************************
!> \brief   Calculation of the [H,r] commutators matrices over Cartesian Gaussian functions.
!> \param qs_env ...
!> \param matrix_hr ...
!> \date    26.07.2016
!> \par     History
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE build_com_hr_matrix(qs_env, matrix_hr)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_hr

      CHARACTER(len=*), PARAMETER :: routineN = 'build_com_hr_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ir
      REAL(KIND=dp)                                      :: eps_ppnl
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (sab_orb, sap_ppnl)
      CALL get_qs_env(qs_env=qs_env, sab_orb=sab_orb, sap_ppnl=sap_ppnl)
      !
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, dft_control=dft_control)
      eps_ppnl = dft_control%qs_control%eps_ppnl
      !
      CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrix_s)
      CPASSERT(.NOT. ASSOCIATED(matrix_hr))
      CALL dbcsr_allocate_matrix_set(matrix_hr, 3)
      DO ir = 1, 3
         ALLOCATE (matrix_hr(ir)%matrix)
         CALL dbcsr_create(matrix_hr(ir)%matrix, template=matrix_s(1, 1)%matrix, &
                           name="COMMUTATOR")
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_hr(ir)%matrix, sab_orb)
         CALL dbcsr_set(matrix_hr(ir)%matrix, 0.0_dp)
      END DO

      CALL build_com_tr_matrix(matrix_hr, qs_kind_set, "ORB", sab_orb)
      CALL build_com_rpnl(matrix_hr, qs_kind_set, sab_orb, sap_ppnl, eps_ppnl)

      CALL timestop(handle)

   END SUBROUTINE build_com_hr_matrix

END MODULE qs_commutators

