!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief dielectric constant data type
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
MODULE dielectric_types

   USE kinds,                           ONLY: dp
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dielectric_types'

   INTEGER, PARAMETER, PUBLIC           :: derivative_cd3 = 1, &
                                           derivative_cd5 = 2, &
                                           derivative_cd7 = 3, &
                                           derivative_fft = 100, &
                                           derivative_fft_use_deps = 101, &
                                           derivative_fft_use_drho = 102

   INTEGER, PARAMETER, PUBLIC           :: rho_dependent = 10, &
                                           spatially_dependent = 20, &
                                           spatially_rho_dependent = 30

   TYPE dielectric_parameters
      LOGICAL                                :: dielec_core_correction = .FALSE.
      INTEGER                                :: dielec_functiontype = rho_dependent
      REAL(KIND=dp)                          :: eps0 = 0.0_dp
      REAL(KIND=dp)                          :: rho_min = 0.0_dp
      REAL(KIND=dp)                          :: rho_max = 0.0_dp
      INTEGER                                :: derivative_method = derivative_cd3
      INTEGER                                :: times_called = 0

      INTEGER                                :: n_aa_cuboidal = 0
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: aa_cuboidal_xxtnt => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: aa_cuboidal_yxtnt => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: aa_cuboidal_zxtnt => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER   :: aa_cuboidal_eps => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER   :: aa_cuboidal_zeta => NULL()

      INTEGER                                :: n_xaa_annular = 0
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: xaa_annular_xxtnt => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: xaa_annular_bctr => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: xaa_annular_brad => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER   :: xaa_annular_eps => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER   :: xaa_annular_zeta => NULL()
   END TYPE dielectric_parameters

   TYPE dielectric_type
      TYPE(dielectric_parameters)       :: params = dielectric_parameters()
      TYPE(pw_type), POINTER            :: eps => NULL()
      TYPE(pw_type), POINTER            :: deps_drho => NULL()
      TYPE(pw_type), DIMENSION(3)     :: dln_eps = pw_type()
   END TYPE dielectric_type

   PUBLIC :: dielectric_type, dielectric_parameters
   PUBLIC :: dielectric_release

CONTAINS

! **************************************************************************************************
!> \brief   deallocates dielectric data type
!> \param dielectric the dielectric data type to be released
!> \param pw_pool pool of the plane wave grid
!> \par History
!>       06.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE dielectric_release(dielectric, pw_pool)
      TYPE(dielectric_type), POINTER                     :: dielectric
      TYPE(pw_pool_type), INTENT(IN), OPTIONAL, POINTER  :: pw_pool

      CHARACTER(len=*), PARAMETER :: routineN = 'dielectric_release'

      INTEGER                                            :: handle, i
      LOGICAL                                            :: can_give_back

      CALL timeset(routineN, handle)
      IF (ASSOCIATED(dielectric)) THEN
         can_give_back = PRESENT(pw_pool)
         IF (can_give_back) can_give_back = ASSOCIATED(pw_pool)
         IF (can_give_back) THEN
            CALL pw_pool%give_back_pw(dielectric%eps)
            CALL pw_pool%give_back_pw(dielectric%deps_drho)
            DEALLOCATE (dielectric%eps, dielectric%deps_drho)
            DO i = 1, 3
               CALL pw_pool%give_back_pw(dielectric%dln_eps(i))
            END DO
         ELSE
            CALL dielectric%eps%release()
            CALL dielectric%deps_drho%release()
            DEALLOCATE (dielectric%eps, dielectric%deps_drho)
            DO i = 1, 3
               CALL dielectric%dln_eps(i)%release()
            END DO
         END IF
         CALL dielectric_parameters_dealloc(dielectric%params)
         DEALLOCATE (dielectric)
      END IF

      CALL timestop(handle)

   END SUBROUTINE dielectric_release

! **************************************************************************************************
!> \brief   deallocates dielectric_parameters type
!> \param dielec_params dielectric parameters
!> \par History
!>       07.2015 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE dielectric_parameters_dealloc(dielec_params)

      TYPE(dielectric_parameters), INTENT(INOUT)         :: dielec_params

      IF (ASSOCIATED(dielec_params%aa_cuboidal_xxtnt)) DEALLOCATE (dielec_params%aa_cuboidal_xxtnt)
      IF (ASSOCIATED(dielec_params%aa_cuboidal_yxtnt)) DEALLOCATE (dielec_params%aa_cuboidal_yxtnt)
      IF (ASSOCIATED(dielec_params%aa_cuboidal_zxtnt)) DEALLOCATE (dielec_params%aa_cuboidal_zxtnt)
      IF (ASSOCIATED(dielec_params%aa_cuboidal_eps)) DEALLOCATE (dielec_params%aa_cuboidal_eps)
      IF (ASSOCIATED(dielec_params%aa_cuboidal_zeta)) DEALLOCATE (dielec_params%aa_cuboidal_zeta)

      IF (ASSOCIATED(dielec_params%xaa_annular_xxtnt)) DEALLOCATE (dielec_params%xaa_annular_xxtnt)
      IF (ASSOCIATED(dielec_params%xaa_annular_bctr)) DEALLOCATE (dielec_params%xaa_annular_bctr)
      IF (ASSOCIATED(dielec_params%xaa_annular_brad)) DEALLOCATE (dielec_params%xaa_annular_brad)
      IF (ASSOCIATED(dielec_params%xaa_annular_eps)) DEALLOCATE (dielec_params%xaa_annular_eps)
      IF (ASSOCIATED(dielec_params%xaa_annular_zeta)) DEALLOCATE (dielec_params%xaa_annular_zeta)

   END SUBROUTINE dielectric_parameters_dealloc

END MODULE dielectric_types
