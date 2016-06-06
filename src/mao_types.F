!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The types needed for the calculation of modified atomic orbitals (MAO)
!> \par History
!>      03.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE mao_types

#include "./base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mao_types'

! **************************************************************************************************
!> \brief Quantities needed for MAO's
!> \author JGH
! **************************************************************************************************
   TYPE mao_type
      INTEGER                       :: ref_basis
   END TYPE mao_type

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mao_env ...
! **************************************************************************************************
   SUBROUTINE create_mao_type(mao_env)
      TYPE(mao_type), POINTER                            :: mao_env

      CHARACTER(len=*), PARAMETER :: routineN = 'create_mao_type', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(mao_env)) THEN
         CALL release_mao_type(mao_env)
      END IF

      ALLOCATE (mao_env)

   END SUBROUTINE create_mao_type

! **************************************************************************************************
!> \brief ...
!> \param mao_env ...
! **************************************************************************************************
   SUBROUTINE release_mao_type(mao_env)
      TYPE(mao_type), POINTER                            :: mao_env

      CHARACTER(len=*), PARAMETER :: routineN = 'release_mao_type', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(mao_env)) THEN
         DEALLOCATE (mao_env)
      END IF

   END SUBROUTINE release_mao_type

! **************************************************************************************************

END MODULE mao_types
