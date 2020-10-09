!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!>
!> \par History
!>     init 17.03.2020
!> \author JGH
! **************************************************************************************************
MODULE qs_fxc

   USE cp_control_types,                ONLY: dft_control_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_p_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_methods,                  ONLY: qs_rho_copy,&
                                              qs_rho_scale_and_add
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_release,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: qs_fxc_create, qs_fxc_release

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fxc'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param ks_env ...
!> \param rho0_struct ...
!> \param rho1_struct ...
!> \param xc_section ...
!> \param do_analytic ...
!> \param do_numeric ...
!> \param is_triplet ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param gxc_rho ...
!> \param gxc_tau ...
! **************************************************************************************************
   SUBROUTINE qs_fxc_create(ks_env, rho0_struct, rho1_struct, xc_section, do_analytic, do_numeric, &
                            is_triplet, fxc_rho, fxc_tau, gxc_rho, gxc_tau)

      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho0_struct, rho1_struct
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: do_analytic, do_numeric, is_triplet
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: fxc_rho, fxc_tau, gxc_rho, gxc_tau

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_fxc_create', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: epsrho = 5.e-4_dp

      INTEGER                                            :: handle, ispin, istep, nspins
      REAL(KIND=dp)                                      :: alpha, beta, exc, oeps1, oeps2
      REAL(KIND=dp), DIMENSION(-2:2)                     :: ak, bl
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_tau_rspace, vxc00
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rhoin

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(fxc_rho))
      CPASSERT(.NOT. ASSOCIATED(fxc_tau))
      CPASSERT(.NOT. ASSOCIATED(gxc_rho))
      CPASSERT(.NOT. ASSOCIATED(gxc_tau))

      IF (do_analytic .AND. .NOT. do_numeric) THEN
         CPABORT("TDA Kernel derivatives could not be calculated: Analytic option not implemented")
         IF (is_triplet) THEN
            CPABORT("TDA Kernel derivatives could not be calculated: Triplet not implemented")
         END IF
      END IF

      IF (do_numeric) THEN

         CPASSERT(ASSOCIATED(rho0_struct))
         CPASSERT(ASSOCIATED(rho1_struct))

         ak(-2:2) = (/1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp/)/12.0_dp
         bl(-2:2) = (/-1.0_dp, 16.0_dp, -30.0_dp, 16.0_dp, -1.0_dp/)/12.0_dp

         CALL get_ks_env(ks_env, dft_control=dft_control, pw_env=pw_env)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

         nspins = dft_control%nspins
         exc = 0.0_dp

         DO istep = -2, 2

            alpha = 1.0_dp
            beta = REAL(istep, KIND=dp)*epsrho
            NULLIFY (rhoin)
            CALL qs_rho_create(rhoin)
            NULLIFY (vxc00, v_tau_rspace)
            IF (is_triplet) THEN
               CPASSERT(nspins == 1)
               CALL qs_rho_copy(rho0_struct, rhoin, auxbas_pw_pool, 2)
               CALL qs_rho_scale_and_add(rhoin, rho1_struct, alpha, 0.5_dp*beta)
               CALL qs_vxc_create(ks_env=ks_env, rho_struct=rhoin, xc_section=xc_section, &
                                  vxc_rho=vxc00, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE.)
               CALL pw_axpy(vxc00(2)%pw, vxc00(1)%pw, -1.0_dp)
            ELSE
               CALL qs_rho_copy(rho0_struct, rhoin, auxbas_pw_pool, nspins)
               CALL qs_rho_scale_and_add(rhoin, rho1_struct, alpha, beta)
               CALL qs_vxc_create(ks_env=ks_env, rho_struct=rhoin, xc_section=xc_section, &
                                  vxc_rho=vxc00, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE.)
            END IF
            CALL qs_rho_release(rhoin)
            IF (.NOT. ASSOCIATED(fxc_rho)) THEN
               ALLOCATE (fxc_rho(nspins))
               DO ispin = 1, nspins
                  CALL pw_pool_create_pw(auxbas_pw_pool, fxc_rho(ispin)%pw, &
                                         in_space=REALSPACE, use_data=REALDATA3D)
                  CALL pw_zero(fxc_rho(ispin)%pw)
               END DO
            END IF
            IF (.NOT. ASSOCIATED(gxc_rho)) THEN
               ALLOCATE (gxc_rho(nspins))
               DO ispin = 1, nspins
                  CALL pw_pool_create_pw(auxbas_pw_pool, gxc_rho(ispin)%pw, &
                                         in_space=REALSPACE, use_data=REALDATA3D)
                  CALL pw_zero(gxc_rho(ispin)%pw)
               END DO
            END IF
            CPASSERT(.NOT. ASSOCIATED(v_tau_rspace))
            DO ispin = 1, nspins
               IF (ak(istep) /= 0.0_dp) THEN
                  CALL pw_axpy(vxc00(ispin)%pw, fxc_rho(ispin)%pw, ak(istep))
               END IF
               IF (bl(istep) /= 0.0_dp) THEN
                  CALL pw_axpy(vxc00(ispin)%pw, gxc_rho(ispin)%pw, bl(istep))
               END IF
            END DO
            DO ispin = 1, SIZE(vxc00)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc00(ispin)%pw)
            END DO
            DEALLOCATE (vxc00)

         END DO

         oeps1 = 1.0_dp/epsrho
         oeps2 = 1.0_dp/(epsrho**2)
         DO ispin = 1, nspins
            CALL pw_scale(fxc_rho(ispin)%pw, oeps1)
            CALL pw_scale(gxc_rho(ispin)%pw, oeps2)
         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_fxc_create

! **************************************************************************************************
!> \brief ...
!> \param ks_env ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param gxc_rho ...
!> \param gxc_tau ...
! **************************************************************************************************
   SUBROUTINE qs_fxc_release(ks_env, fxc_rho, fxc_tau, gxc_rho, gxc_tau)

      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: fxc_rho, fxc_tau, gxc_rho, gxc_tau

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_fxc_release', routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL get_ks_env(ks_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      IF (ASSOCIATED(fxc_rho)) THEN
         DO ispin = 1, SIZE(fxc_rho)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, fxc_rho(ispin)%pw)
         END DO
         DEALLOCATE (fxc_rho)
      END IF
      IF (ASSOCIATED(fxc_tau)) THEN
         DO ispin = 1, SIZE(fxc_tau)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, fxc_tau(ispin)%pw)
         END DO
         DEALLOCATE (fxc_tau)
      END IF
      IF (ASSOCIATED(gxc_rho)) THEN
         DO ispin = 1, SIZE(gxc_rho)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, gxc_rho(ispin)%pw)
         END DO
         DEALLOCATE (gxc_rho)
      END IF
      IF (ASSOCIATED(gxc_tau)) THEN
         DO ispin = 1, SIZE(gxc_tau)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, gxc_tau(ispin)%pw)
         END DO
         DEALLOCATE (gxc_tau)
      END IF

   END SUBROUTINE qs_fxc_release

END MODULE qs_fxc
