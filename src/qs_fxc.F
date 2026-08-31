!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Setup Routine for Fxc Potentials
!         https://en.wikipedia.org/wiki/Finite_difference_coefficient
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
!>     complete refactoring 08.2026
!> \author JGH
! **************************************************************************************************
MODULE qs_fxc

   USE cp_control_types,                ONLY: dft_control_type
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_lval,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grids,                        ONLY: pw_grid_compare
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_fxc_atom,                     ONLY: fxc_atom_calc
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_rho_methods,                  ONLY: qs_rho_copy,&
                                              qs_rho_scale_and_add,&
                                              qs_rho_transfer
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE xc,                              ONLY: xc_calc_2nd_deriv_analytical,&
                                              xc_calc_2nd_deriv_numerical,&
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
   PUBLIC :: qs_fxc_create, qs_fxc_prep, qs_fxc_apply
   PUBLIC :: qs_fxc_fdiff

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fxc'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho0_struct ...
!> \param rho1_struct ...
!> \param rho0_atom_set ...
!> \param xc_section ...
!> \param do_onecenter ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param rho1_atom_set ...
!> \param do_scale ...
!> \param is_triplet ...
!> \param spinflip ...
!> \param no_weights ...
!> \param uf_grid_results ...
!> \param pw_env_ext ...
!> \param kind_set_external ...
!> \param para_env_external ...
!> \param compute_virial ...
!> \param virial_xc ...
! **************************************************************************************************
   SUBROUTINE qs_fxc_create(qs_env, rho0_struct, rho1_struct, rho0_atom_set, &
                            xc_section, do_onecenter, &
                            fxc_rho, fxc_tau, rho1_atom_set, &
                            do_scale, is_triplet, spinflip, no_weights, uf_grid_results, &
                            pw_env_ext, kind_set_external, para_env_external, &
                            compute_virial, virial_xc)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: rho0_struct, rho1_struct
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho0_atom_set
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: do_onecenter
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho, fxc_tau
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho1_atom_set
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_scale, is_triplet, spinflip, &
                                                            no_weights, uf_grid_results
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_ext
      TYPE(qs_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: kind_set_external
      TYPE(mp_para_env_type), INTENT(IN), OPTIONAL       :: para_env_external
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial_xc

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_fxc_create'

      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: do_virial, do_w, ret_uf, uf_grid
      REAL(KIND=dp)                                      :: factor
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_nlcc_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool, xc_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho_lo, fxc_rho_uf, fxc_tau_lo, &
                                                            fxc_tau_uf, rho_r
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_nlcc, weights, weights_uf
      TYPE(qs_rho_type), POINTER                         :: rho0_uf, rho1_uf

      CALL timeset(routineN, handle)

      do_virial = .FALSE.
      IF (PRESENT(compute_virial)) do_virial = compute_virial

      CALL get_qs_env(qs_env, dft_control=dft_control)

      IF (PRESENT(pw_env_ext)) THEN
         pw_env => pw_env_ext
      ELSE
         CALL get_qs_env(qs_env, pw_env=pw_env)
      END IF
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, xc_pw_pool=xc_pw_pool)
      uf_grid = .NOT. pw_grid_compare(auxbas_pw_pool%pw_grid, xc_pw_pool%pw_grid)

      nspins = dft_control%nspins
      IF (ASSOCIATED(fxc_rho)) THEN
         CPASSERT(nspins == SIZE(fxc_rho))
      END IF
      IF (ASSOCIATED(fxc_tau)) THEN
         CPASSERT(nspins == SIZE(fxc_tau))
      END IF

      NULLIFY (rho_nlcc, rho_nlcc_g)
      CALL get_qs_env(qs_env, rho_nlcc=rho_nlcc, rho_nlcc_g=rho_nlcc_g)
      IF (ASSOCIATED(rho_nlcc)) THEN
         NULLIFY (rho_r, rho_g)
         CALL qs_rho_get(rho0_struct, rho_r=rho_r, rho_g=rho_g)
         factor = 1.0_dp
         DO ispin = 1, nspins
            CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
            CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
         END DO
      END IF

      do_w = .TRUE.
      IF (PRESENT(no_weights)) do_w = .NOT. no_weights
      IF (do_w) THEN
         CALL get_qs_env(qs_env, xcint_weights=weights)
      ELSE
         NULLIFY (weights)
      END IF

      NULLIFY (fxc_rho_lo, fxc_tau_lo)
      IF (uf_grid) THEN
         IF (PRESENT(uf_grid_results)) THEN
            ret_uf = uf_grid_results
         ELSE
            ret_uf = .FALSE.
         END IF
         NULLIFY (weights_uf)
         IF (ASSOCIATED(weights)) THEN
            ALLOCATE (weights_uf)
            CALL xc_pw_pool%create_pw(weights_uf)
            BLOCK
               TYPE(pw_c1d_gs_type) :: weights_g, weights_g_uf
               CALL auxbas_pw_pool%create_pw(weights_g)
               CALL xc_pw_pool%create_pw(weights_g_uf)
               CALL pw_transfer(weights, weights_g)
               CALL pw_transfer(weights_g, weights_g_uf)
               CALL pw_transfer(weights_g_uf, weights_uf)
               CALL xc_pw_pool%give_back_pw(weights_g_uf)
               CALL auxbas_pw_pool%give_back_pw(weights_g)
            END BLOCK
         END IF
         !
         ALLOCATE (rho0_uf, rho1_uf)
         CALL qs_rho_create(rho0_uf)
         CALL qs_rho_create(rho1_uf)
         CALL qs_rho_transfer(rho0_struct, rho0_uf, auxbas_pw_pool, xc_pw_pool)
         CALL qs_rho_transfer(rho1_struct, rho1_uf, auxbas_pw_pool, xc_pw_pool)
         !
         NULLIFY (fxc_rho_uf, fxc_tau_uf)
         CALL qs_fxc_calculate(rho0_uf, rho1_uf, xc_section, weights_uf, xc_pw_pool, &
                               fxc_rho_uf, fxc_tau_uf, &
                               is_triplet=is_triplet, spinflip=spinflip, &
                               compute_virial=do_virial, virial_xc=virial_xc)
         !
         CALL qs_rho_release(rho0_uf)
         CALL qs_rho_release(rho1_uf)
         DEALLOCATE (rho0_uf, rho1_uf)
         IF (ASSOCIATED(weights_uf)) THEN
            CALL xc_pw_pool%give_back_pw(weights_uf)
            DEALLOCATE (weights_uf)
         END IF
         IF (ret_uf) THEN
            fxc_rho_lo => fxc_rho_uf
            fxc_tau_lo => fxc_tau_uf
         ELSE
            IF (ASSOCIATED(fxc_rho_uf)) THEN
               ALLOCATE (fxc_rho_lo(nspins))
               DO ispin = 1, nspins
                  CALL auxbas_pw_pool%create_pw(fxc_rho_lo(ispin))
                  BLOCK
                     TYPE(pw_c1d_gs_type) :: fxc_g, fxc_g_uf
                     CALL auxbas_pw_pool%create_pw(fxc_g)
                     CALL xc_pw_pool%create_pw(fxc_g_uf)
                     CALL pw_transfer(fxc_rho_uf(ispin), fxc_g_uf)
                     CALL pw_transfer(fxc_g_uf, fxc_g)
                     CALL pw_transfer(fxc_g, fxc_rho_lo(ispin))
                     CALL xc_pw_pool%give_back_pw(fxc_g_uf)
                     CALL auxbas_pw_pool%give_back_pw(fxc_g)
                  END BLOCK
                  CALL xc_pw_pool%give_back_pw(fxc_rho_uf(ispin))
               END DO
               DEALLOCATE (fxc_rho_uf)
            END IF
            IF (ASSOCIATED(fxc_tau_uf)) THEN
               ALLOCATE (fxc_tau_lo(nspins))
               DO ispin = 1, nspins
                  CALL auxbas_pw_pool%create_pw(fxc_tau_lo(ispin))
                  BLOCK
                     TYPE(pw_c1d_gs_type) :: fxc_g, fxc_g_uf
                     CALL auxbas_pw_pool%create_pw(fxc_g)
                     CALL xc_pw_pool%create_pw(fxc_g_uf)
                     CALL pw_transfer(fxc_tau_uf(ispin), fxc_g_uf)
                     CALL pw_transfer(fxc_g_uf, fxc_g)
                     CALL pw_transfer(fxc_g, fxc_tau_lo(ispin))
                     CALL xc_pw_pool%give_back_pw(fxc_g_uf)
                     CALL auxbas_pw_pool%give_back_pw(fxc_g)
                  END BLOCK
                  CALL xc_pw_pool%give_back_pw(fxc_tau_uf(ispin))
               END DO
            END IF
         END IF

      ELSE
         CALL qs_fxc_calculate(rho0_struct, rho1_struct, xc_section, weights, auxbas_pw_pool, &
                               fxc_rho_lo, fxc_tau_lo, &
                               is_triplet=is_triplet, spinflip=spinflip, &
                               compute_virial=do_virial, virial_xc=virial_xc)
      END IF

      ! de-apply NLCC density
      IF (ASSOCIATED(rho_nlcc)) THEN
         factor = -1.0_dp
         DO ispin = 1, nspins
            CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
            CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
         END DO
      END IF

      ! return potentials
      IF (ASSOCIATED(fxc_rho)) THEN
         DO ispin = 1, MIN(SIZE(fxc_rho_lo), SIZE(fxc_rho))
            CALL pw_transfer(fxc_rho_lo(ispin), fxc_rho(ispin))
         END DO
         DO ispin = 1, SIZE(fxc_rho_lo)
            CALL auxbas_pw_pool%give_back_pw(fxc_rho_lo(ispin))
         END DO
         DEALLOCATE (fxc_rho_lo)
      ELSE
         fxc_rho => fxc_rho_lo
      END IF
      IF (ASSOCIATED(fxc_tau)) THEN
         IF (ASSOCIATED(fxc_tau_lo)) THEN
            DO ispin = 1, MIN(SIZE(fxc_tau_lo), SIZE(fxc_tau))
               CALL pw_transfer(fxc_tau_lo(ispin), fxc_tau(ispin))
            END DO
            DO ispin = 1, SIZE(fxc_tau_lo)
               CALL auxbas_pw_pool%give_back_pw(fxc_tau_lo(ispin))
            END DO
            DEALLOCATE (fxc_tau_lo)
         ELSE
            DO ispin = 1, nspins
               CALL pw_zero(fxc_tau(ispin))
            END DO
         END IF
      ELSE
         fxc_tau => fxc_tau_lo
      END IF

      IF (do_onecenter) THEN
         CALL fxc_atom_calc(qs_env, rho0_atom_set, rho1_atom_set, xc_section, &
                            do_scale=do_scale, do_triplet=is_triplet, do_sf=spinflip, &
                            para_env_ext=para_env_external, &
                            kind_set_external=kind_set_external)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_fxc_create

! **************************************************************************************************
!> \brief ...
!> \param rho0 ...
!> \param rho1 ...
!> \param xc_section ...
!> \param weights ...
!> \param auxbas_pw_pool ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param is_triplet ...
!> \param spinflip ...
!> \param compute_virial ...
!> \param virial_xc ...
! **************************************************************************************************
   SUBROUTINE qs_fxc_calculate(rho0, rho1, xc_section, weights, auxbas_pw_pool, &
                               fxc_rho, fxc_tau, is_triplet, spinflip, &
                               compute_virial, virial_xc)

      TYPE(qs_rho_type), POINTER                         :: rho0, rho1
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho, fxc_tau
      LOGICAL, INTENT(IN), OPTIONAL                      :: is_triplet, spinflip, compute_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial_xc

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_fxc_calculate'

      INTEGER                                            :: handle, ispin, mspins, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: do_analytic, do_sf, do_triplet, &
                                                            do_virial, lsd
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: vxg
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho0_g, rho1_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho0_r, rho1_r, tau0_r, tau1_r
      TYPE(qs_rho_type), POINTER                         :: rhot0, rhot1
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_derivative_set_type)                       :: xc_deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho0_set, rho1_set

      CALL timeset(routineN, handle)

      do_triplet = .FALSE.
      IF (PRESENT(is_triplet)) do_triplet = is_triplet

      do_sf = .FALSE.
      IF (PRESENT(spinflip)) do_sf = spinflip

      do_virial = .FALSE.
      IF (PRESENT(compute_virial)) do_virial = compute_virial

      do_analytic = section_get_lval(xc_section, "2ND_DERIV_ANALYTICAL")

      CALL qs_rho_get(rho0, rho_r=rho0_r, rho_g=rho0_g, tau_r=tau0_r)
      CALL qs_rho_get(rho1, rho_r=rho1_r, tau_r=tau1_r)
      NULLIFY (rho1_g)

      mspins = SIZE(rho0_r)
      nspins = SIZE(rho0_r)
      lsd = (nspins == 2)
      IF (nspins == 1 .AND. do_triplet) THEN
         nspins = 2
         lsd = .TRUE.
      ELSE IF (do_sf) THEN
         nspins = 1
         mspins = 1
         lsd = .TRUE.
      END IF

      CPASSERT(.NOT. ASSOCIATED(fxc_rho))
      CPASSERT(.NOT. ASSOCIATED(fxc_tau))
      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)
      ALLOCATE (fxc_rho(mspins))
      DO ispin = 1, mspins
         CALL auxbas_pw_pool%create_pw(fxc_rho(ispin))
         CALL pw_zero(fxc_rho(ispin))
      END DO
      IF (needs%tau .OR. needs%tau_spin) THEN
         IF (.NOT. ASSOCIATED(tau1_r)) THEN
            CPABORT("Tau-dependent functionals requires allocated kinetic energy density grid")
         END IF
         ALLOCATE (fxc_tau(mspins))
         DO ispin = 1, mspins
            CALL auxbas_pw_pool%create_pw(fxc_tau(ispin))
            CALL pw_zero(fxc_tau(ispin))
         END DO
      END IF

      IF (mspins == 1 .AND. do_triplet) THEN
         ! split the density and response density arrays for triplet calculation
         ALLOCATE (rhot0)
         CALL qs_rho_create(rhot0)
         CALL qs_rho_copy(rho0, rhot0, auxbas_pw_pool, 2, factor=2.0_dp)
         !
         ALLOCATE (rhot1)
         CALL qs_rho_create(rhot1)
         CALL qs_rho_copy(rho1, rhot1, auxbas_pw_pool, 2, factor=2.0_dp)
         !
         CALL qs_rho_get(rhot0, rho_r=rho0_r, rho_g=rho0_g, tau_r=tau0_r)
         CALL qs_rho_get(rhot1, rho_r=rho1_r, tau_r=tau1_r)

         CALL xc_prep_2nd_deriv(xc_deriv_set, rho0_set, rho0_r, auxbas_pw_pool, weights, &
                                xc_section=xc_section, tau_r=tau0_r)
         bo = rho1_r(1)%pw_grid%bounds_local
         ! create the place where to store the argument for the functionals
         CALL xc_rho_set_create(rho1_set, bo, &
                                rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                                drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                                tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

         ! calculate the arguments needed by the functionals
         CALL xc_rho_set_update(rho1_set, rho1_r, rho1_g, tau1_r, needs, &
                                section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                                section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                                auxbas_pw_pool, spinflip=do_sf)
      ELSE
         CALL xc_prep_2nd_deriv(xc_deriv_set, rho0_set, rho0_r, auxbas_pw_pool, weights, &
                                xc_section=xc_section, tau_r=tau0_r)
         bo = rho1_r(1)%pw_grid%bounds_local
         ! create the place where to store the argument for the functionals
         CALL xc_rho_set_create(rho1_set, bo, &
                                rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                                drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                                tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

         ! calculate the arguments needed by the functionals
         CALL xc_rho_set_update(rho1_set, rho1_r, rho1_g, tau1_r, needs, &
                                section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                                section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                                auxbas_pw_pool, spinflip=do_sf)
      END IF

      IF (mspins == 1 .AND. do_triplet .AND. do_analytic) THEN

         CALL xc_calc_2nd_deriv_analytical(fxc_rho, fxc_tau, xc_deriv_set, rho0_set, &
                                           rho1_set, auxbas_pw_pool, xc_section, &
                                           gapw=.FALSE., vxg=vxg, tddfpt_fac=-1.0_dp, spinflip=do_sf, &
                                           compute_virial=compute_virial, virial_xc=virial_xc)

      ELSE IF (do_analytic) THEN

         CALL xc_calc_2nd_deriv_analytical(fxc_rho, fxc_tau, xc_deriv_set, rho0_set, &
                                           rho1_set, auxbas_pw_pool, xc_section, &
                                           gapw=.FALSE., vxg=vxg, spinflip=do_sf, &
                                           compute_virial=compute_virial, virial_xc=virial_xc)

      ELSE

         CALL xc_calc_2nd_deriv_numerical(fxc_rho, fxc_tau, rho0_set, rho1_r, rho1_g, tau1_r, &
                                          auxbas_pw_pool, weights, xc_section, &
                                          do_triplet, compute_virial, virial_xc, xc_deriv_set)

      END IF

      IF (mspins == 1 .AND. do_triplet) THEN
         CALL qs_rho_release(rhot0)
         DEALLOCATE (rhot0)
         CALL qs_rho_release(rhot1)
         DEALLOCATE (rhot1)
      END IF

      CALL xc_dset_release(xc_deriv_set)
      CALL xc_rho_set_release(rho0_set)
      CALL xc_rho_set_release(rho1_set)

      CALL timestop(handle)

   END SUBROUTINE qs_fxc_calculate

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho0_struct ...
!> \param xc_rho_set ...
!> \param xc_deriv_set ...
!> \param xc_section ...
!> \param pw_env_ext ...
!> \param is_triplet ...
! **************************************************************************************************
   SUBROUTINE qs_fxc_prep(qs_env, rho0_struct, xc_rho_set, xc_deriv_set, &
                          xc_section, pw_env_ext, is_triplet)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: rho0_struct
      TYPE(xc_rho_set_type)                              :: xc_rho_set
      TYPE(xc_derivative_set_type)                       :: xc_deriv_set
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_ext
      LOGICAL, INTENT(IN), OPTIONAL                      :: is_triplet

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_fxc_prep'

      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: uf_grid
      REAL(KIND=dp)                                      :: factor
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_nlcc_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool, xc_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_nlcc, weights, weights_uf
      TYPE(qs_rho_type), POINTER                         :: rho0_uf

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)

      IF (PRESENT(pw_env_ext)) THEN
         pw_env => pw_env_ext
      ELSE
         CALL get_qs_env(qs_env, pw_env=pw_env)
      END IF
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, xc_pw_pool=xc_pw_pool)
      uf_grid = .NOT. pw_grid_compare(auxbas_pw_pool%pw_grid, xc_pw_pool%pw_grid)

      nspins = dft_control%nspins

      NULLIFY (rho_nlcc, rho_nlcc_g)
      CALL get_qs_env(qs_env, rho_nlcc=rho_nlcc, rho_nlcc_g=rho_nlcc_g)
      IF (ASSOCIATED(rho_nlcc)) THEN
         NULLIFY (rho_r, rho_g)
         CALL qs_rho_get(rho0_struct, rho_r=rho_r, rho_g=rho_g)
         factor = 1.0_dp
         DO ispin = 1, nspins
            CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
            CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
         END DO
      END IF

      NULLIFY (weights)
      CALL get_qs_env(qs_env, xcint_weights=weights)

      IF (uf_grid) THEN
         NULLIFY (weights_uf)
         IF (ASSOCIATED(weights)) THEN
            ALLOCATE (weights_uf)
            CALL xc_pw_pool%create_pw(weights_uf)
            BLOCK
               TYPE(pw_c1d_gs_type) :: weights_g, weights_g_uf
               CALL auxbas_pw_pool%create_pw(weights_g)
               CALL xc_pw_pool%create_pw(weights_g_uf)
               CALL pw_transfer(weights, weights_g)
               CALL pw_transfer(weights_g, weights_g_uf)
               CALL pw_transfer(weights_g_uf, weights_uf)
               CALL xc_pw_pool%give_back_pw(weights_g_uf)
               CALL auxbas_pw_pool%give_back_pw(weights_g)
            END BLOCK
         END IF
         !
         ALLOCATE (rho0_uf)
         CALL qs_rho_create(rho0_uf)
         CALL qs_rho_transfer(rho0_struct, rho0_uf, auxbas_pw_pool, xc_pw_pool)
         !
         CALL qs_fxc_deriv(rho0_uf, xc_rho_set, xc_deriv_set, &
                           xc_section, weights_uf, xc_pw_pool, is_triplet)
         !
         CALL qs_rho_release(rho0_uf)
         DEALLOCATE (rho0_uf)
         IF (ASSOCIATED(weights_uf)) THEN
            CALL xc_pw_pool%give_back_pw(weights_uf)
            DEALLOCATE (weights_uf)
         END IF
      ELSE
         CALL qs_fxc_deriv(rho0_struct, xc_rho_set, xc_deriv_set, &
                           xc_section, weights, auxbas_pw_pool, is_triplet)
      END IF

      ! de-apply NLCC density
      IF (ASSOCIATED(rho_nlcc)) THEN
         factor = -1.0_dp
         DO ispin = 1, nspins
            CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
            CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_fxc_prep

! **************************************************************************************************
!> \brief ...
!> \param rho0 ...
!> \param xc_rho_set ...
!> \param xc_deriv_set ...
!> \param xc_section ...
!> \param weights ...
!> \param auxbas_pw_pool ...
!> \param is_triplet ...
! **************************************************************************************************
   SUBROUTINE qs_fxc_deriv(rho0, xc_rho_set, xc_deriv_set, xc_section, weights, auxbas_pw_pool, &
                           is_triplet)

      TYPE(qs_rho_type), POINTER                         :: rho0
      TYPE(xc_rho_set_type)                              :: xc_rho_set
      TYPE(xc_derivative_set_type)                       :: xc_deriv_set
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      LOGICAL, INTENT(IN)                                :: is_triplet

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_fxc_deriv'

      INTEGER                                            :: handle
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho0_r, tau0_r
      TYPE(qs_rho_type), POINTER                         :: rhot0

      CALL timeset(routineN, handle)

      NULLIFY (rho0_r, tau0_r)
      IF (is_triplet) THEN
         ! split the density and response density arrays for triplet calculation
         ALLOCATE (rhot0)
         CALL qs_rho_create(rhot0)
         CALL qs_rho_copy(rho0, rhot0, auxbas_pw_pool, 2, factor=2.0_dp)
         CALL qs_rho_get(rhot0, rho_r=rho0_r, tau_r=tau0_r)
         CALL xc_prep_2nd_deriv(xc_deriv_set, xc_rho_set, rho0_r, auxbas_pw_pool, weights, &
                                xc_section=xc_section, tau_r=tau0_r)
         CALL qs_rho_release(rhot0)
         DEALLOCATE (rhot0)
      ELSE
         CALL qs_rho_get(rho0, rho_r=rho0_r, tau_r=tau0_r)
         CALL xc_prep_2nd_deriv(xc_deriv_set, xc_rho_set, rho0_r, auxbas_pw_pool, weights, &
                                xc_section=xc_section, tau_r=tau0_r)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_fxc_deriv

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param xc_deriv_set ...
!> \param xc_rho_set ...
!> \param rho1_struct ...
!> \param rho0_atom_set ...
!> \param xc_section ...
!> \param do_onecenter ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param rho1_atom_set ...
!> \param do_scale ...
!> \param is_triplet ...
!> \param spinflip ...
!> \param pw_env_ext ...
!> \param kind_set_external ...
!> \param para_env_external ...
!> \param compute_virial ...
!> \param virial_xc ...
! **************************************************************************************************
   SUBROUTINE qs_fxc_apply(qs_env, xc_deriv_set, xc_rho_set, rho1_struct, rho0_atom_set, &
                           xc_section, do_onecenter, fxc_rho, fxc_tau, rho1_atom_set, &
                           do_scale, is_triplet, spinflip, pw_env_ext, &
                           kind_set_external, para_env_external, compute_virial, virial_xc)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(xc_derivative_set_type)                       :: xc_deriv_set
      TYPE(xc_rho_set_type)                              :: xc_rho_set
      TYPE(qs_rho_type), POINTER                         :: rho1_struct
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho0_atom_set
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: do_onecenter
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho, fxc_tau
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho1_atom_set
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_scale, is_triplet, spinflip
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_ext
      TYPE(qs_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: kind_set_external
      TYPE(mp_para_env_type), INTENT(IN), OPTIONAL       :: para_env_external
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial_xc

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_fxc_apply'

      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: do_virial, uf_grid
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool, xc_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho_lo, fxc_rho_uf, fxc_tau_lo, &
                                                            fxc_tau_uf
      TYPE(pw_r3d_rs_type), POINTER                      :: weights, weights_uf
      TYPE(qs_rho_type), POINTER                         :: rho1_uf

      CALL timeset(routineN, handle)

      do_virial = .FALSE.
      IF (PRESENT(compute_virial)) do_virial = compute_virial

      CALL get_qs_env(qs_env, dft_control=dft_control)

      IF (PRESENT(pw_env_ext)) THEN
         pw_env => pw_env_ext
      ELSE
         CALL get_qs_env(qs_env, pw_env=pw_env)
      END IF
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, xc_pw_pool=xc_pw_pool)
      uf_grid = .NOT. pw_grid_compare(auxbas_pw_pool%pw_grid, xc_pw_pool%pw_grid)

      nspins = dft_control%nspins
      IF (ASSOCIATED(fxc_rho)) THEN
         CPASSERT(nspins == SIZE(fxc_rho))
      END IF
      IF (ASSOCIATED(fxc_tau)) THEN
         CPASSERT(nspins == SIZE(fxc_tau))
      END IF

      CALL get_qs_env(qs_env, xcint_weights=weights)

      NULLIFY (fxc_rho_lo, fxc_tau_lo)
      IF (uf_grid) THEN
         NULLIFY (weights_uf)
         IF (ASSOCIATED(weights)) THEN
            ALLOCATE (weights_uf)
            CALL xc_pw_pool%create_pw(weights_uf)
            BLOCK
               TYPE(pw_c1d_gs_type) :: weights_g, weights_g_uf
               CALL auxbas_pw_pool%create_pw(weights_g)
               CALL xc_pw_pool%create_pw(weights_g_uf)
               CALL pw_transfer(weights, weights_g)
               CALL pw_transfer(weights_g, weights_g_uf)
               CALL pw_transfer(weights_g_uf, weights_uf)
               CALL xc_pw_pool%give_back_pw(weights_g_uf)
               CALL auxbas_pw_pool%give_back_pw(weights_g)
            END BLOCK
         END IF
         !
         ALLOCATE (rho1_uf)
         CALL qs_rho_create(rho1_uf)
         CALL qs_rho_transfer(rho1_struct, rho1_uf, auxbas_pw_pool, xc_pw_pool)
         !
         NULLIFY (fxc_rho_uf, fxc_tau_uf)
         CALL qs_fxc_eval(xc_deriv_set, xc_rho_set, rho1_uf, xc_section, &
                          weights_uf, xc_pw_pool, fxc_rho_uf, fxc_tau_uf, &
                          is_triplet=is_triplet, spinflip=spinflip, &
                          compute_virial=do_virial, virial_xc=virial_xc)
         !
         CALL qs_rho_release(rho1_uf)
         DEALLOCATE (rho1_uf)
         IF (ASSOCIATED(weights_uf)) THEN
            CALL xc_pw_pool%give_back_pw(weights_uf)
            DEALLOCATE (weights_uf)
         END IF
         IF (ASSOCIATED(fxc_rho_uf)) THEN
            ALLOCATE (fxc_rho_lo(nspins))
            DO ispin = 1, nspins
               CALL auxbas_pw_pool%create_pw(fxc_rho_lo(ispin))
               BLOCK
                  TYPE(pw_c1d_gs_type) :: fxc_g, fxc_g_uf
                  CALL auxbas_pw_pool%create_pw(fxc_g)
                  CALL xc_pw_pool%create_pw(fxc_g_uf)
                  CALL pw_transfer(fxc_rho_uf(ispin), fxc_g_uf)
                  CALL pw_transfer(fxc_g_uf, fxc_g)
                  CALL pw_transfer(fxc_g, fxc_rho_lo(ispin))
                  CALL xc_pw_pool%give_back_pw(fxc_g_uf)
                  CALL auxbas_pw_pool%give_back_pw(fxc_g)
               END BLOCK
               CALL xc_pw_pool%give_back_pw(fxc_rho_uf(ispin))
            END DO
            DEALLOCATE (fxc_rho_uf)
         END IF
         IF (ASSOCIATED(fxc_tau_uf)) THEN
            ALLOCATE (fxc_tau_lo(nspins))
            DO ispin = 1, nspins
               CALL auxbas_pw_pool%create_pw(fxc_tau_lo(ispin))
               BLOCK
                  TYPE(pw_c1d_gs_type) :: fxc_g, fxc_g_uf
                  CALL auxbas_pw_pool%create_pw(fxc_g)
                  CALL xc_pw_pool%create_pw(fxc_g_uf)
                  CALL pw_transfer(fxc_tau_uf(ispin), fxc_g_uf)
                  CALL pw_transfer(fxc_g_uf, fxc_g)
                  CALL pw_transfer(fxc_g, fxc_tau_lo(ispin))
                  CALL xc_pw_pool%give_back_pw(fxc_g_uf)
                  CALL auxbas_pw_pool%give_back_pw(fxc_g)
               END BLOCK
               CALL xc_pw_pool%give_back_pw(fxc_tau_uf(ispin))
            END DO
         END IF

      ELSE
         CALL qs_fxc_eval(xc_deriv_set, xc_rho_set, rho1_struct, xc_section, &
                          weights, auxbas_pw_pool, fxc_rho_lo, fxc_tau_lo, &
                          is_triplet=is_triplet, spinflip=spinflip, &
                          compute_virial=do_virial, virial_xc=virial_xc)
      END IF

      ! return potentials
      IF (ASSOCIATED(fxc_rho)) THEN
         DO ispin = 1, MIN(SIZE(fxc_rho_lo), SIZE(fxc_rho))
            CALL pw_transfer(fxc_rho_lo(ispin), fxc_rho(ispin))
         END DO
         DO ispin = 1, SIZE(fxc_rho_lo)
            CALL auxbas_pw_pool%give_back_pw(fxc_rho_lo(ispin))
         END DO
         DEALLOCATE (fxc_rho_lo)
      ELSE
         fxc_rho => fxc_rho_lo
      END IF
      IF (ASSOCIATED(fxc_tau)) THEN
         IF (ASSOCIATED(fxc_tau_lo)) THEN
            DO ispin = 1, MIN(SIZE(fxc_tau_lo), SIZE(fxc_tau))
               CALL pw_transfer(fxc_tau_lo(ispin), fxc_tau(ispin))
            END DO
            DO ispin = 1, SIZE(fxc_tau_lo)
               CALL auxbas_pw_pool%give_back_pw(fxc_tau_lo(ispin))
            END DO
            DEALLOCATE (fxc_tau_lo)
         ELSE
            DO ispin = 1, nspins
               CALL pw_zero(fxc_tau(ispin))
            END DO
         END IF
      ELSE
         fxc_tau => fxc_tau_lo
      END IF

      IF (do_onecenter) THEN
         CALL fxc_atom_calc(qs_env, rho0_atom_set, rho1_atom_set, xc_section, &
                            do_scale=do_scale, do_triplet=is_triplet, do_sf=spinflip, &
                            para_env_ext=para_env_external, &
                            kind_set_external=kind_set_external)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_fxc_apply

! **************************************************************************************************
!> \brief ...
!> \param xc_deriv_set ...
!> \param xc_rho_set ...
!> \param rho1 ...
!> \param xc_section ...
!> \param weights ...
!> \param auxbas_pw_pool ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param is_triplet ...
!> \param spinflip ...
!> \param compute_virial ...
!> \param virial_xc ...
! **************************************************************************************************
   SUBROUTINE qs_fxc_eval(xc_deriv_set, xc_rho_set, rho1, xc_section, weights, auxbas_pw_pool, &
                          fxc_rho, fxc_tau, is_triplet, spinflip, &
                          compute_virial, virial_xc)

      TYPE(xc_derivative_set_type)                       :: xc_deriv_set
      TYPE(xc_rho_set_type)                              :: xc_rho_set
      TYPE(qs_rho_type), POINTER                         :: rho1
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho, fxc_tau
      LOGICAL, INTENT(IN), OPTIONAL                      :: is_triplet, spinflip, compute_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial_xc

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_fxc_eval'

      INTEGER                                            :: handle, ispin, mspins, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: do_analytic, do_sf, do_triplet, &
                                                            do_virial, lsd
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: vxg
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho1_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho1_r, tau1_r
      TYPE(qs_rho_type), POINTER                         :: rhot1
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho1_set

      CALL timeset(routineN, handle)

      do_triplet = .FALSE.
      IF (PRESENT(is_triplet)) do_triplet = is_triplet

      do_sf = .FALSE.
      IF (PRESENT(spinflip)) do_sf = spinflip

      do_virial = .FALSE.
      IF (PRESENT(compute_virial)) do_virial = compute_virial

      do_analytic = section_get_lval(xc_section, "2ND_DERIV_ANALYTICAL")

      CALL qs_rho_get(rho1, rho_r=rho1_r, tau_r=tau1_r)
      NULLIFY (rho1_g)

      mspins = SIZE(rho1_r)
      nspins = SIZE(rho1_r)
      lsd = (nspins == 2)
      IF (nspins == 1 .AND. do_triplet) THEN
         nspins = 2
         lsd = .TRUE.
      ELSE IF (do_sf) THEN
         nspins = 1
         mspins = 1
         lsd = .TRUE.
      END IF

      CPASSERT(.NOT. ASSOCIATED(fxc_rho))
      CPASSERT(.NOT. ASSOCIATED(fxc_tau))
      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)
      ALLOCATE (fxc_rho(mspins))
      DO ispin = 1, mspins
         CALL auxbas_pw_pool%create_pw(fxc_rho(ispin))
         CALL pw_zero(fxc_rho(ispin))
      END DO
      IF (needs%tau .OR. needs%tau_spin) THEN
         IF (.NOT. ASSOCIATED(tau1_r)) THEN
            CPABORT("Tau-dependent functionals requires allocated kinetic energy density grid")
         END IF
         ALLOCATE (fxc_tau(mspins))
         DO ispin = 1, mspins
            CALL auxbas_pw_pool%create_pw(fxc_tau(ispin))
            CALL pw_zero(fxc_tau(ispin))
         END DO
      END IF

      IF (mspins == 1 .AND. do_triplet) THEN
         ! split the density and response density arrays for triplet calculation
         ALLOCATE (rhot1)
         CALL qs_rho_create(rhot1)
         CALL qs_rho_copy(rho1, rhot1, auxbas_pw_pool, 2, factor=2.0_dp)
         !
         CALL qs_rho_get(rhot1, rho_r=rho1_r, tau_r=tau1_r)
      END IF

      bo = rho1_r(1)%pw_grid%bounds_local
      ! create the place where to store the argument for the functionals
      CALL xc_rho_set_create(rho1_set, bo, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      ! calculate the arguments needed by the functionals
      CALL xc_rho_set_update(rho1_set, rho1_r, rho1_g, tau1_r, needs, &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                             auxbas_pw_pool, spinflip=do_sf)

      IF (mspins == 1 .AND. do_triplet .AND. do_analytic) THEN

         CALL xc_calc_2nd_deriv_analytical(fxc_rho, fxc_tau, xc_deriv_set, xc_rho_set, &
                                           rho1_set, auxbas_pw_pool, xc_section, &
                                           gapw=.FALSE., vxg=vxg, tddfpt_fac=-1.0_dp, spinflip=do_sf, &
                                           compute_virial=compute_virial, virial_xc=virial_xc)

      ELSE IF (do_analytic) THEN

         CALL xc_calc_2nd_deriv_analytical(fxc_rho, fxc_tau, xc_deriv_set, xc_rho_set, &
                                           rho1_set, auxbas_pw_pool, xc_section, &
                                           gapw=.FALSE., vxg=vxg, spinflip=do_sf, &
                                           compute_virial=compute_virial, virial_xc=virial_xc)

      ELSE

         CALL xc_calc_2nd_deriv_numerical(fxc_rho, fxc_tau, xc_rho_set, rho1_r, rho1_g, tau1_r, &
                                          auxbas_pw_pool, weights, xc_section, &
                                          do_triplet, compute_virial, virial_xc, xc_deriv_set)

      END IF

      IF (mspins == 1 .AND. do_triplet) THEN
         CALL qs_rho_release(rhot1)
         DEALLOCATE (rhot1)
      END IF
      CALL xc_rho_set_release(rho1_set)

      CALL timestop(handle)

   END SUBROUTINE qs_fxc_eval

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho0_struct ...
!> \param rho1_struct ...
!> \param xc_section ...
!> \param accuracy ...
!> \param fxc_rho ...
!> \param fxc_tau ...
!> \param is_triplet ...
!> \param spinflip ...
! **************************************************************************************************
   SUBROUTINE qs_fxc_fdiff(qs_env, rho0_struct, rho1_struct, xc_section, accuracy, &
                           fxc_rho, fxc_tau, is_triplet, spinflip)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: rho0_struct, rho1_struct
      TYPE(section_vals_type), POINTER                   :: xc_section
      INTEGER, INTENT(IN)                                :: accuracy
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho, fxc_tau
      LOGICAL, INTENT(IN), OPTIONAL                      :: is_triplet, spinflip

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_fxc_fdiff'
      REAL(KIND=dp), PARAMETER                           :: epsrho = 5.e-4_dp

      INTEGER                                            :: handle, ispin, istep, nspins, nstep
      LOGICAL                                            :: do_sf, do_triplet
      REAL(KIND=dp)                                      :: alpha, beta, exc, oeps1
      REAL(KIND=dp), DIMENSION(-4:4)                     :: ak
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: v_tau_rspace, vxc00
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rhoin

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(fxc_rho))
      CPASSERT(.NOT. ASSOCIATED(fxc_tau))
      CPASSERT(ASSOCIATED(rho0_struct))
      CPASSERT(ASSOCIATED(rho1_struct))

      do_triplet = .FALSE.
      IF (PRESENT(is_triplet)) do_triplet = is_triplet

      do_sf = .FALSE.
      IF (PRESENT(spinflip)) do_sf = spinflip
      IF (do_sf) THEN
         CPABORT("Spin Flip TDDFT only available with analytic 2nd xc derivatives")
      END IF

      ak = 0.0_dp
      SELECT CASE (accuracy)
      CASE (:4)
         nstep = 2
         ak(-2:2) = [1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp]/12.0_dp
      CASE (5:7)
         nstep = 3
         ak(-3:3) = [-1.0_dp, 9.0_dp, -45.0_dp, 0.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]/60.0_dp
      CASE (8:)
         nstep = 4
         ak(-4:4) = [1.0_dp, -32.0_dp/3.0_dp, 56.0_dp, -224.0_dp, 0.0_dp, &
                     224.0_dp, -56.0_dp, 32.0_dp/3.0_dp, -1.0_dp]/280.0_dp
      END SELECT

      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      nspins = dft_control%nspins
      exc = 0.0_dp

      DO istep = -nstep, nstep

         IF (ak(istep) /= 0.0_dp) THEN
            alpha = 1.0_dp
            beta = REAL(istep, KIND=dp)*epsrho
            NULLIFY (rhoin)
            ALLOCATE (rhoin)
            CALL qs_rho_create(rhoin)
            NULLIFY (vxc00, v_tau_rspace)
            IF (do_triplet) THEN
               CPASSERT(nspins == 1)
               ! rhoin = (0.5 rho0, 0.5 rho0)
               CALL qs_rho_copy(rho0_struct, rhoin, auxbas_pw_pool, 2)
               ! rhoin = (0.5 rho0 + 0.5 rho1, 0.5 rho0)
               CALL qs_rho_scale_and_add(rhoin, rho1_struct, alpha, 0.5_dp*beta)
               CALL qs_vxc_create(ks_env, rhoin, xc_section, vxc00, v_tau_rspace, exc)
               CALL pw_axpy(vxc00(2), vxc00(1), -1.0_dp)
               IF (ASSOCIATED(v_tau_rspace)) CALL pw_axpy(v_tau_rspace(2), v_tau_rspace(1), -1.0_dp)
            ELSE
               CALL qs_rho_copy(rho0_struct, rhoin, auxbas_pw_pool, nspins)
               CALL qs_rho_scale_and_add(rhoin, rho1_struct, alpha, beta)
               CALL qs_vxc_create(ks_env, rhoin, xc_section, vxc00, v_tau_rspace, exc)
            END IF
            CALL qs_rho_release(rhoin)
            DEALLOCATE (rhoin)
            IF (.NOT. ASSOCIATED(fxc_rho)) THEN
               ALLOCATE (fxc_rho(nspins))
               DO ispin = 1, nspins
                  CALL auxbas_pw_pool%create_pw(fxc_rho(ispin))
                  CALL pw_zero(fxc_rho(ispin))
               END DO
            END IF
            DO ispin = 1, nspins
               CALL pw_axpy(vxc00(ispin), fxc_rho(ispin), ak(istep))
            END DO
            DO ispin = 1, SIZE(vxc00)
               CALL auxbas_pw_pool%give_back_pw(vxc00(ispin))
            END DO
            DEALLOCATE (vxc00)
            IF (ASSOCIATED(v_tau_rspace)) THEN
               IF (.NOT. ASSOCIATED(fxc_tau)) THEN
                  ALLOCATE (fxc_tau(nspins))
                  DO ispin = 1, nspins
                     CALL auxbas_pw_pool%create_pw(fxc_tau(ispin))
                     CALL pw_zero(fxc_tau(ispin))
                  END DO
               END IF
               DO ispin = 1, nspins
                  CALL pw_axpy(v_tau_rspace(ispin), fxc_tau(ispin), ak(istep))
               END DO
               DO ispin = 1, SIZE(v_tau_rspace)
                  CALL auxbas_pw_pool%give_back_pw(v_tau_rspace(ispin))
               END DO
               DEALLOCATE (v_tau_rspace)
            END IF
         END IF

      END DO

      oeps1 = 1.0_dp/epsrho
      DO ispin = 1, nspins
         CALL pw_scale(fxc_rho(ispin), oeps1)
      END DO
      IF (ASSOCIATED(fxc_tau)) THEN
         DO ispin = 1, nspins
            CALL pw_scale(fxc_tau(ispin), oeps1)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_fxc_fdiff

END MODULE qs_fxc
