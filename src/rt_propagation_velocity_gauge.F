!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to perform the RTP in the velocity gauge
! **************************************************************************************************

MODULE rt_propagation_velocity_gauge

   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_print
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_velocity_gauge'

   PUBLIC :: velocity_gauge_ks_matrix

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE velocity_gauge_ks_matrix(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'velocity_gauge_ks_matrix'

      INTEGER                                            :: handle
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_im

      CALL timeset(routineN, handle)

      NULLIFY (matrix_ks, matrix_ks_im)

      CALL get_qs_env(qs_env, &
                      matrix_ks=matrix_ks, &
                      matrix_ks_im=matrix_ks_im)

      CALL dbcsr_print(matrix_ks(1)%matrix)
      CALL dbcsr_print(matrix_ks_im(1)%matrix)

      CALL timestop(handle)

   END SUBROUTINE velocity_gauge_ks_matrix

END MODULE rt_propagation_velocity_gauge
