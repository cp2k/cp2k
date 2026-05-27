!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Experimental CP2K-native GPW real-space-grid path for SKALA TorchScript models.
! **************************************************************************************************
MODULE skala_gpw_functional
   USE cell_types,                      ONLY: cell_type
   USE input_section_types,             ONLY: section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE particle_types,                  ONLY: particle_type
   USE pw_methods,                      ONLY: pw_scale,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE skala_gpw_features,              ONLY: skala_gpw_feature_build,&
                                              skala_gpw_feature_release,&
                                              skala_gpw_feature_type
   USE skala_torch_api,                 ONLY: skala_torch_model_get_exc,&
                                              skala_torch_model_load,&
                                              skala_torch_model_release,&
                                              skala_torch_model_type
   USE string_utilities,                ONLY: uppercase
   USE torch_api,                       ONLY: torch_tensor_backward_scalar,&
                                              torch_tensor_data_ptr,&
                                              torch_tensor_grad,&
                                              torch_tensor_release,&
                                              torch_tensor_type,&
                                              torch_use_cuda
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
   USE xc_util,                         ONLY: xc_pw_divergence,&
                                              xc_requires_tmp_g
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_gpw_functional'

   PUBLIC :: skala_gpw_eval, xc_section_uses_native_skala_grid

   TYPE(skala_torch_model_type), SAVE                  :: cached_model
   CHARACTER(len=default_path_length), SAVE            :: cached_model_path = ""
   LOGICAL, SAVE                                       :: cached_model_loaded = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief Return true if the GAUXC subsection requests the CP2K-native GPW grid path.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_section_uses_native_skala_grid(xc_section) RESULT(uses_native_grid)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: uses_native_grid

      TYPE(section_vals_type), POINTER                   :: gauxc_section

      uses_native_grid = .FALSE.
      gauxc_section => get_gauxc_section(xc_section)
      IF (ASSOCIATED(gauxc_section)) THEN
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID", l_val=uses_native_grid)
      END IF

   END FUNCTION xc_section_uses_native_skala_grid

! **************************************************************************************************
!> \brief Evaluate SKALA energy and first derivatives on a CP2K GPW grid.
!> \param vxc_rho ...
!> \param vxc_tau ...
!> \param exc ...
!> \param rho_r ...
!> \param rho_g ...
!> \param tau ...
!> \param xc_section ...
!> \param weights ...
!> \param pw_pool ...
!> \param particle_set ...
!> \param cell ...
!> \param compute_virial ...
!> \param virial_xc ...
!> \param just_energy ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_eval(vxc_rho, vxc_tau, exc, rho_r, rho_g, tau, xc_section, &
                             weights, pw_pool, particle_set, cell, compute_virial, virial_xc, just_energy)
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: vxc_rho, vxc_tau
      REAL(KIND=dp), INTENT(OUT)                         :: exc
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: tau
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN)                                :: compute_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: virial_xc
      LOGICAL, INTENT(IN), OPTIONAL                      :: just_energy

      CHARACTER(len=default_path_length)                 :: model_path
      INTEGER                                            :: nspins, xc_deriv_method_id, &
                                                            xc_rho_smooth_id
      LOGICAL                                            :: lsd, my_just_energy
      TYPE(skala_gpw_feature_type)                       :: features
      TYPE(torch_tensor_type)                            :: density_grad_t, exc_tensor, grad_grad_t, &
                                                            kin_grad_t
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      virial_xc = 0.0_dp
      exc = 0.0_dp
      my_just_energy = .FALSE.
      IF (PRESENT(just_energy)) my_just_energy = just_energy

      IF (compute_virial) THEN
         CALL cp_abort(__LOCATION__, &
                       "Native SKALA GPW stress/virial is not implemented yet.")
      END IF
      IF (.NOT. ASSOCIATED(rho_g)) THEN
         CALL cp_abort(__LOCATION__, &
                       "Native SKALA GPW requires the reciprocal-space density to form density gradients.")
      END IF
      IF (.NOT. ASSOCIATED(tau)) THEN
         CALL cp_abort(__LOCATION__, &
                       "Native SKALA GPW requires the kinetic-energy density.")
      END IF

      nspins = SIZE(rho_r)
      lsd = (nspins /= 1)
      CALL get_skala_model_path(xc_section, model_path)
      ! The current SKALA 1.1 TorchScript model creates internal ragged-index tensors on CPU.
      CALL torch_use_cuda(.FALSE.)
      CALL ensure_model_loaded(model_path)

      IF (lsd) THEN
         needs%rho_spin = .TRUE.
         needs%drho_spin = .TRUE.
         needs%tau_spin = .TRUE.
      ELSE
         needs%rho = .TRUE.
         needs%drho = .TRUE.
         needs%tau = .TRUE.
      END IF

      CALL section_vals_val_get(xc_section, "XC_GRID%XC_DERIV", i_val=xc_deriv_method_id)
      CALL section_vals_val_get(xc_section, "XC_GRID%XC_SMOOTH_RHO", i_val=xc_rho_smooth_id)

      CALL xc_rho_set_create(rho_set, &
                             rho_r(1)%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "density_cutoff"), &
                             drho_cutoff=section_get_rval(xc_section, "gradient_cutoff"), &
                             tau_cutoff=section_get_rval(xc_section, "tau_cutoff"))
      CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs, &
                             xc_deriv_method_id, xc_rho_smooth_id, pw_pool)

      CALL skala_gpw_feature_build(features, rho_set, rho_r, particle_set, cell, &
                                   requires_grad=(.NOT. my_just_energy), weights=weights)
      CALL skala_torch_model_get_exc(cached_model, features%inputs, &
                                     features%grid_weights_t, exc_tensor, exc)

      IF (.NOT. my_just_energy) THEN
         CALL torch_tensor_backward_scalar(exc_tensor)
         CALL torch_tensor_grad(features%density_t, density_grad_t)
         CALL torch_tensor_grad(features%grad_t, grad_grad_t)
         CALL torch_tensor_grad(features%kin_t, kin_grad_t)
         CALL build_vxc_from_feature_grads(vxc_rho, vxc_tau, rho_r, pw_pool, features, &
                                           density_grad_t, grad_grad_t, kin_grad_t, &
                                           xc_deriv_method_id)
         CALL torch_tensor_release(density_grad_t)
         CALL torch_tensor_release(grad_grad_t)
         CALL torch_tensor_release(kin_grad_t)
      END IF

      CALL torch_tensor_release(exc_tensor)
      CALL skala_gpw_feature_release(features)
      CALL xc_rho_set_release(rho_set, pw_pool=pw_pool)
      CALL torch_use_cuda(.TRUE.)

   END SUBROUTINE skala_gpw_eval

! **************************************************************************************************
!> \brief Fill CP2K VXC real-space arrays from Torch feature gradients.
!> \param vxc_rho ...
!> \param vxc_tau ...
!> \param rho_r ...
!> \param pw_pool ...
!> \param features ...
!> \param density_grad_t ...
!> \param grad_grad_t ...
!> \param kin_grad_t ...
!> \param xc_deriv_method_id ...
! **************************************************************************************************
   SUBROUTINE build_vxc_from_feature_grads(vxc_rho, vxc_tau, rho_r, pw_pool, features, &
                                           density_grad_t, grad_grad_t, kin_grad_t, &
                                           xc_deriv_method_id)
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: vxc_rho, vxc_tau, rho_r
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(torch_tensor_type), INTENT(IN)                :: density_grad_t, grad_grad_t, kin_grad_t
      INTEGER, INTENT(IN)                                :: xc_deriv_method_id

      INTEGER                                            :: i, ipt, ispin, j, k, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp)                                      :: dvol_inv
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: density_grad, kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grad_grad
      TYPE(pw_c1d_gs_type)                               :: tmp_g, vxc_g
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: grad_pw

      NULLIFY (density_grad, grad_grad, kin_grad)
      CALL torch_tensor_data_ptr(density_grad_t, density_grad)
      CALL torch_tensor_data_ptr(grad_grad_t, grad_grad)
      CALL torch_tensor_data_ptr(kin_grad_t, kin_grad)

      nspins = SIZE(rho_r)
      bo = rho_r(1)%pw_grid%bounds_local
      dvol_inv = 1.0_dp/rho_r(1)%pw_grid%dvol

      ALLOCATE (vxc_rho(nspins), vxc_tau(nspins))
      DO ispin = 1, nspins
         CALL pw_pool%create_pw(vxc_rho(ispin))
         CALL pw_pool%create_pw(vxc_tau(ispin))
         CALL pw_zero(vxc_rho(ispin))
         CALL pw_zero(vxc_tau(ispin))
      END DO

      IF (xc_requires_tmp_g(xc_deriv_method_id) .OR. rho_r(1)%pw_grid%spherical) THEN
         CALL pw_pool%create_pw(vxc_g)
         IF (.NOT. rho_r(1)%pw_grid%spherical) CALL pw_pool%create_pw(tmp_g)
      END IF

      DO ispin = 1, nspins
         DO i = 1, 3
            CALL pw_pool%create_pw(grad_pw(i))
            CALL pw_zero(grad_pw(i))
         END DO

         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  ipt = features%feature_index(i, j, k)
                  IF (nspins == 1) THEN
                     vxc_rho(1)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                 (density_grad(ipt, 1) + density_grad(ipt, 2))
                     vxc_tau(1)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                 (kin_grad(ipt, 1) + kin_grad(ipt, 2))
                     grad_pw(1)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                 (grad_grad(ipt, 1, 1) + grad_grad(ipt, 1, 2))
                     grad_pw(2)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                 (grad_grad(ipt, 2, 1) + grad_grad(ipt, 2, 2))
                     grad_pw(3)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                 (grad_grad(ipt, 3, 1) + grad_grad(ipt, 3, 2))
                  ELSE
                     vxc_rho(ispin)%array(i, j, k) = dvol_inv*density_grad(ipt, ispin)
                     vxc_tau(ispin)%array(i, j, k) = dvol_inv*kin_grad(ipt, ispin)
                     grad_pw(1)%array(i, j, k) = dvol_inv*grad_grad(ipt, 1, ispin)
                     grad_pw(2)%array(i, j, k) = dvol_inv*grad_grad(ipt, 2, ispin)
                     grad_pw(3)%array(i, j, k) = dvol_inv*grad_grad(ipt, 3, ispin)
                  END IF
               END DO
            END DO
         END DO

         DO i = 1, 3
            CALL pw_scale(grad_pw(i), -1.0_dp)
         END DO
         CALL xc_pw_divergence(xc_deriv_method_id, grad_pw, tmp_g, vxc_g, vxc_rho(ispin))

         DO i = 1, 3
            CALL pw_pool%give_back_pw(grad_pw(i))
         END DO
      END DO

      IF (ASSOCIATED(vxc_g%pw_grid)) CALL pw_pool%give_back_pw(vxc_g)
      IF (ASSOCIATED(tmp_g%pw_grid)) CALL pw_pool%give_back_pw(tmp_g)

   END SUBROUTINE build_vxc_from_feature_grads

! **************************************************************************************************
!> \brief Load and cache the TorchScript SKALA model.
!> \param model_path ...
! **************************************************************************************************
   SUBROUTINE ensure_model_loaded(model_path)
      CHARACTER(len=*), INTENT(IN)                       :: model_path

      IF (cached_model_loaded) THEN
         IF (TRIM(cached_model_path) == TRIM(model_path)) RETURN
         CALL skala_torch_model_release(cached_model)
         cached_model_loaded = .FALSE.
      END IF

      CALL skala_torch_model_load(cached_model, TRIM(model_path))
      cached_model_path = model_path
      cached_model_loaded = .TRUE.

   END SUBROUTINE ensure_model_loaded

! **************************************************************************************************
!> \brief Resolve the SKALA TorchScript model path from the GAUXC subsection.
!> \param xc_section ...
!> \param model_path ...
! **************************************************************************************************
   SUBROUTINE get_skala_model_path(xc_section, model_path)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      CHARACTER(len=default_path_length), INTENT(OUT)    :: model_path

      CHARACTER(len=default_path_length)                 :: model_key
      INTEGER                                            :: env_status
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      gauxc_section => get_gauxc_section(xc_section)
      IF (.NOT. ASSOCIATED(gauxc_section)) THEN
         CPABORT("Native SKALA GPW requires an XC_FUNCTIONAL%GAUXC section")
      END IF

      CALL section_vals_val_get(gauxc_section, "MODEL", c_val=model_path)
      model_key = ADJUSTL(model_path)
      CALL uppercase(model_key)
      IF (TRIM(model_key) == "NONE" .OR. TRIM(model_key) == "") THEN
         CPABORT("Native SKALA GPW requires GAUXC%MODEL SKALA or a TorchScript model path")
      ELSE IF (TRIM(model_key) == "SKALA") THEN
         CALL GET_ENVIRONMENT_VARIABLE("GAUXC_SKALA_MODEL", model_path, STATUS=env_status)
         IF (env_status /= 0 .OR. LEN_TRIM(model_path) == 0) THEN
            CPABORT("MODEL SKALA requires the GAUXC_SKALA_MODEL environment variable")
         END IF
      END IF

   END SUBROUTINE get_skala_model_path

! **************************************************************************************************
!> \brief Return the first GAUXC functional subsection, if present.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_gauxc_section(xc_section) RESULT(gauxc_section)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      INTEGER                                            :: ifun
      TYPE(section_vals_type), POINTER                   :: functionals, xc_fun

      NULLIFY (gauxc_section)
      IF (.NOT. ASSOCIATED(xc_section)) RETURN

      functionals => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      IF (.NOT. ASSOCIATED(functionals)) RETURN

      ifun = 0
      DO
         ifun = ifun + 1
         xc_fun => section_vals_get_subs_vals2(functionals, i_section=ifun)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         IF (xc_fun%section%name == "GAUXC") THEN
            gauxc_section => xc_fun
            EXIT
         END IF
      END DO

   END FUNCTION get_gauxc_section

END MODULE skala_gpw_functional
