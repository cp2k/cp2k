!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Does all kind of post scf calculations for GPW/GAPW
!> \par History
!>      Taken out from qs_scf_post_gpw
!> \author JGH
! **************************************************************************************************
MODULE qs_elf_methods
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_derive,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Global parameters
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_elf_methods'

   PUBLIC :: qs_elf_calc

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param elf_r ...
!> \param rho_cutoff ...
! **************************************************************************************************
   SUBROUTINE qs_elf_calc(qs_env, elf_r, rho_cutoff)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_type), DIMENSION(:), INTENT(IN)            :: elf_r
      REAL(kind=dp), INTENT(IN)                          :: rho_cutoff

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_elf_calc'
      INTEGER, DIMENSION(3, 3), PARAMETER :: nd = RESHAPE((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3, 3/))
      REAL(KIND=dp), PARAMETER                           :: ELFCUT = 0.0001_dp, &
                                                            f18 = (1.0_dp/8.0_dp), &
                                                            f23 = (2.0_dp/3.0_dp), &
                                                            f53 = (5.0_dp/3.0_dp)

      INTEGER                                            :: handle, i, idir, ispin, j, k, nspin
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: deriv_pw, drho_r_valid, tau_r_valid
      REAL(kind=dp)                                      :: cfermi, elf_kernel, norm_drho, rho_53, &
                                                            udvol
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_struct_ao
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: tmp_g
      TYPE(pw_type), DIMENSION(3)                        :: drho_r
      TYPE(pw_type), DIMENSION(:), POINTER               :: rho_struct_r, tau_struct_r
      TYPE(pw_type), DIMENSION(:, :), POINTER            :: drho_struct_r
      TYPE(pw_type), POINTER                             :: rho_r, tau_r
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_struct

      CALL timeset(routineN, handle)

      NULLIFY (rho_struct, rho_r, tau_r, pw_env, auxbas_pw_pool, pw_pools, ks_env)
      NULLIFY (rho_struct_ao, rho_struct_r, tau_struct_r, drho_struct_r)

      CALL get_qs_env(qs_env, ks_env=ks_env, pw_env=pw_env, rho=rho_struct)

      CALL qs_rho_get(rho_struct, &
                      rho_ao_kp=rho_struct_ao, &
                      rho_r=rho_struct_r, &
                      tau_r=tau_struct_r, &
                      drho_r=drho_struct_r, &
                      tau_r_valid=tau_r_valid, &
                      drho_r_valid=drho_r_valid)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      pw_pools=pw_pools)
      nspin = SIZE(rho_struct_r)
      bo = rho_struct_r(1)%pw_grid%bounds_local
      cfermi = (3.0_dp/10.0_dp)*(pi*pi*3.0_dp)**f23

      ! In this case, we need a work matrix containing tau in g space
      ! We will not have further use for it, so we will need only one
      IF (.NOT. tau_r_valid) THEN
         ALLOCATE (tau_r)
         CALL auxbas_pw_pool%create_pw(tau_r, &
                                       use_data=REALDATA3D, in_space=REALSPACE)
      END IF
      IF (.NOT. tau_r_valid .OR. .NOT. drho_r_valid) THEN
         CALL auxbas_pw_pool%create_pw(tmp_g, &
                                       use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      END IF
      IF (.NOT. drho_r_valid) THEN
         DO idir = 1, 3
            CALL auxbas_pw_pool%create_pw(drho_r(idir), &
                                          use_data=REALDATA3D, in_space=REALSPACE)
         END DO
      END IF

      DO ispin = 1, nspin
         rho_r => rho_struct_r(ispin)
         IF (tau_r_valid) THEN
            tau_r => tau_struct_r(ispin)
         ELSE
            rho_ao => rho_struct_ao(ispin, :)
            CALL pw_zero(tau_r)
            CALL calculate_rho_elec(matrix_p_kp=rho_ao, &
                                    rho=tau_r, &
                                    rho_gspace=tmp_g, &
                                    ks_env=ks_env, soft_valid=.FALSE., &
                                    compute_tau=.TRUE.)
         END IF

         IF (drho_r_valid) THEN
            drho_r(:) = drho_struct_r(:, ispin)
         ELSE
            deriv_pw = .FALSE.
            IF (deriv_pw) THEN
               udvol = 1.0_dp/rho_r%pw_grid%dvol
               DO idir = 1, 3
                  CALL pw_transfer(rho_r, tmp_g)
                  CALL pw_derive(tmp_g, nd(:, idir))
                  CALL pw_transfer(tmp_g, drho_r(idir))
               END DO

            ELSE
               DO idir = 1, 3
                  rho_ao => rho_struct_ao(ispin, :)
                  CALL calculate_rho_elec(matrix_p_kp=rho_ao, &
                                          rho=drho_r(idir), &
                                          rho_gspace=tmp_g, &
                                          ks_env=ks_env, soft_valid=.FALSE., &
                                          compute_tau=.FALSE., compute_grad=.TRUE., idir=idir)

               END DO
            END IF
         END IF

         ! Calculate elf_r
!$OMP        PARALLEL DO DEFAULT(NONE) SHARED(bo,elf_r, ispin, drho_r,rho_r, tau_r, cfermi, rho_cutoff)&
!$OMP                    PRIVATE(k,j,i, norm_drho, rho_53, elf_kernel)
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  norm_drho = drho_r(1)%cr3d(i, j, k)**2 + &
                              drho_r(2)%cr3d(i, j, k)**2 + &
                              drho_r(3)%cr3d(i, j, k)**2
                  norm_drho = norm_drho/MAX(rho_r%cr3d(i, j, k), rho_cutoff)
                  rho_53 = cfermi*MAX(rho_r%cr3d(i, j, k), rho_cutoff)**f53
                  elf_kernel = (tau_r%cr3d(i, j, k) - f18*norm_drho) + 2.87E-5_dp
                  elf_kernel = (elf_kernel/rho_53)**2
                  elf_r(ispin)%cr3d(i, j, k) = 1.0_dp/(1.0_dp + elf_kernel)
                  IF (elf_r(ispin)%cr3d(i, j, k) < ELFCUT) elf_r(ispin)%cr3d(i, j, k) = 0.0_dp
               END DO
            END DO
         END DO
      END DO

      IF (.NOT. drho_r_valid) THEN
         DO idir = 1, 3
            CALL auxbas_pw_pool%give_back_pw(drho_r(idir))
         END DO
      END IF
      IF (.NOT. tau_r_valid) THEN
         CALL auxbas_pw_pool%give_back_pw(tau_r)
         DEALLOCATE (tau_r)
      END IF
      IF (.NOT. tau_r_valid .OR. .NOT. drho_r_valid) THEN
         CALL auxbas_pw_pool%give_back_pw(tmp_g)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_elf_calc

END MODULE qs_elf_methods
