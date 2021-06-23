!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  https://en.wikipedia.org/wiki/Finite_difference_coefficient
!---------------------------------------------------------------------------------------------------
!Derivative    Accuracy            4       3       2       1       0       1       2      3       4
!---------------------------------------------------------------------------------------------------
!    1             2                                    -1/2       0     1/2
!                  4                            1/12    -2/3       0     2/3   -1/12
!                  6                   -1/60    3/20    -3/4       0     3/4   -3/20   1/60
!                  8           1/280  -4/105     1/5    -4/5       0     4/5    -1/5  4/105  -1/280
!---------------------------------------------------------------------------------------------------
!    2             2                                       1      -2       1
!                  4                           -1/12     4/3    -5/2     4/3   -1/12
!                  6                    1/90   -3/20     3/2  -49/18     3/2   -3/20   1/90
!                  8          -1/560   8/315    -1/5     8/5 -205/72     8/5    -1/5  8/315  -1/560
!---------------------------------------------------------------------------------------------------
!> \par History
!>     init 17.03.2020
!> \author JGH
! **************************************************************************************************
MODULE qs_fxc

   USE cp_control_types,                ONLY: dft_control_type
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
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
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: qs_fxc_fdiff, qs_fxc_analytic, qs_fgxc_create, qs_fgxc_release

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fxc'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param rho0 ...
!> \param rho1_r ...
!> \param xc_section ...
!> \param auxbas_pw_pool ...
!> \param is_triplet ...
!> \param v_xc ...
! **************************************************************************************************
   SUBROUTINE qs_fxc_analytic(rho0, rho1_r, xc_section, auxbas_pw_pool, is_triplet, v_xc)

      TYPE(qs_rho_type), POINTER                         :: rho0
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_r
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      LOGICAL, INTENT(IN)                                :: is_triplet
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_xc

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_fxc_analytic', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: lsd
      REAL(KIND=dp)                                      :: fac
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho0_g, rho0_r, rho1_g, tau0_r, tau1_r
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho0_set, rho1_set

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(v_xc))

      CALL qs_rho_get(rho0, rho_r=rho0_r, rho_g=rho0_g, tau_r=tau0_r)
      nspins = SIZE(rho0_r)

      ALLOCATE (v_xc(nspins))
      DO ispin = 1, nspins
         NULLIFY (v_xc(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(v_xc(ispin)%pw)
      END DO

      lsd = (nspins == 2)
      fac = 0._dp
      IF (is_triplet .AND. nspins == 1) fac = -1.0_dp

      NULLIFY (deriv_set, rho0_set, rho1_set)
      NULLIFY (rho1_g, tau1_r)
      bo = rho1_r(1)%pw%pw_grid%bounds_local
      CALL xc_rho_set_create(rho1_set, bo, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))
      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)
      ! calculate the arguments needed by the functionals
      CALL xc_rho_set_update(rho1_set, rho1_r, rho1_g, tau1_r, needs, &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                             auxbas_pw_pool)
      CALL xc_prep_2nd_deriv(deriv_set, rho0_set, rho0_r, auxbas_pw_pool, xc_section=xc_section)
      CALL xc_calc_2nd_deriv(v_xc, deriv_set, rho0_set, &
                             rho1_set, auxbas_pw_pool, xc_section=xc_section, tddfpt_fac=fac)
      CALL xc_dset_release(deriv_set)
      CALL xc_rho_set_release(rho0_set)
      CALL xc_rho_set_release(rho1_set)

      CALL timestop(handle)

   END SUBROUTINE qs_fxc_analytic

! **************************************************************************************************
!> \brief ...
!> \param ks_env ...
!> \param rho0_struct ...
!> \param rho1_struct ...
!> \param xc_section ...
!> \param accuracy ...
!> \param is_triplet ...
!> \param fxc_rho ...
!> \param fxc_tau ...
! **************************************************************************************************
   SUBROUTINE qs_fxc_fdiff(ks_env, rho0_struct, rho1_struct, xc_section, accuracy, is_triplet, &
                           fxc_rho, fxc_tau)

      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho0_struct, rho1_struct
      TYPE(section_vals_type), POINTER                   :: xc_section
      INTEGER, INTENT(IN)                                :: accuracy
      LOGICAL, INTENT(IN)                                :: is_triplet
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: fxc_rho, fxc_tau

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_fxc_fdiff', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: epsrho = 5.e-4_dp

      INTEGER                                            :: handle, ispin, istep, nspins, nstep
      REAL(KIND=dp)                                      :: alpha, beta, exc, oeps1
      REAL(KIND=dp), DIMENSION(-4:4)                     :: ak
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_tau_rspace, vxc00
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rhoin

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(fxc_rho))
      CPASSERT(.NOT. ASSOCIATED(fxc_tau))
      CPASSERT(ASSOCIATED(rho0_struct))
      CPASSERT(ASSOCIATED(rho1_struct))

      ak = 0.0_dp
      SELECT CASE (accuracy)
      CASE (:4)
         nstep = 2
         ak(-2:2) = (/1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp/)/12.0_dp
      CASE (5:7)
         nstep = 3
         ak(-3:3) = (/-1.0_dp, 9.0_dp, -45.0_dp, 0.0_dp, 45.0_dp, -9.0_dp, 1.0_dp/)/60.0_dp
      CASE (8:)
         nstep = 4
         ak(-4:4) = (/1.0_dp, -32.0_dp/3.0_dp, 56.0_dp, -224.0_dp, 0.0_dp, &
                      224.0_dp, -56.0_dp, 32.0_dp/3.0_dp, -1.0_dp/)/280.0_dp
      END SELECT

      CALL get_ks_env(ks_env, dft_control=dft_control, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      nspins = dft_control%nspins
      exc = 0.0_dp

      DO istep = -nstep, nstep

         IF (ak(istep) /= 0.0_dp) THEN
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
!deb        CPASSERT(.NOT. ASSOCIATED(v_tau_rspace))
            DO ispin = 1, nspins
               CALL pw_axpy(vxc00(ispin)%pw, fxc_rho(ispin)%pw, ak(istep))
            END DO
            DO ispin = 1, SIZE(vxc00)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc00(ispin)%pw)
            END DO
            DEALLOCATE (vxc00)
         END IF

      END DO

      oeps1 = 1.0_dp/epsrho
      DO ispin = 1, nspins
         CALL pw_scale(fxc_rho(ispin)%pw, oeps1)
      END DO

      CALL timestop(handle)

   END SUBROUTINE qs_fxc_fdiff

! **************************************************************************************************
!> \brief ...
!> \param ks_env ...
!> \param rho0_struct ...
!> \param rho1_struct ...
!> \param xc_section ...
!> \param accuracy ...
!> \param is_triplet ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param gxc_rho ...
!> \param gxc_tau ...
! **************************************************************************************************
   SUBROUTINE qs_fgxc_create(ks_env, rho0_struct, rho1_struct, xc_section, accuracy, is_triplet, &
                             fxc_rho, fxc_tau, gxc_rho, gxc_tau)

      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho0_struct, rho1_struct
      TYPE(section_vals_type), POINTER                   :: xc_section
      INTEGER, INTENT(IN)                                :: accuracy
      LOGICAL, INTENT(IN)                                :: is_triplet
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: fxc_rho, fxc_tau, gxc_rho, gxc_tau

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_fgxc_create', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: epsrho = 5.e-4_dp

      INTEGER                                            :: handle, ispin, istep, nspins, nstep
      REAL(KIND=dp)                                      :: alpha, beta, exc, oeps1, oeps2
      REAL(KIND=dp), DIMENSION(-4:4)                     :: ak, bl
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
      CPASSERT(ASSOCIATED(rho0_struct))
      CPASSERT(ASSOCIATED(rho1_struct))

      ak = 0.0_dp
      bl = 0.0_dp
      SELECT CASE (accuracy)
      CASE (:4)
         nstep = 2
         ak(-2:2) = (/1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp/)/12.0_dp
         bl(-2:2) = (/-1.0_dp, 16.0_dp, -30.0_dp, 16.0_dp, -1.0_dp/)/12.0_dp
      CASE (5:7)
         nstep = 3
         ak(-3:3) = (/-1.0_dp, 9.0_dp, -45.0_dp, 0.0_dp, 45.0_dp, -9.0_dp, 1.0_dp/)/60.0_dp
         bl(-3:3) = (/2.0_dp, -27.0_dp, 270.0_dp, -490.0_dp, 270.0_dp, -27.0_dp, 2.0_dp/)/180.0_dp
      CASE (8:)
         nstep = 4
         ak(-4:4) = (/1.0_dp, -32.0_dp/3.0_dp, 56.0_dp, -224.0_dp, 0.0_dp, &
                      224.0_dp, -56.0_dp, 32.0_dp/3.0_dp, -1.0_dp/)/280.0_dp
         bl(-4:4) = (/-1.0_dp, 128.0_dp/9.0_dp, -112.0_dp, 896.0_dp, -14350.0_dp/9.0_dp, &
                      896.0_dp, -112.0_dp, 128.0_dp/9.0_dp, -1.0_dp/)/560.0_dp
      END SELECT

      CALL get_ks_env(ks_env, dft_control=dft_control, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      nspins = dft_control%nspins
      exc = 0.0_dp

      DO istep = -nstep, nstep

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

      CALL timestop(handle)

   END SUBROUTINE qs_fgxc_create

! **************************************************************************************************
!> \brief ...
!> \param ks_env ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param gxc_rho ...
!> \param gxc_tau ...
! **************************************************************************************************
   SUBROUTINE qs_fgxc_release(ks_env, fxc_rho, fxc_tau, gxc_rho, gxc_tau)

      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: fxc_rho, fxc_tau, gxc_rho, gxc_tau

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_fgxc_release', &
         routineP = moduleN//':'//routineN

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

   END SUBROUTINE qs_fgxc_release

END MODULE qs_fxc
