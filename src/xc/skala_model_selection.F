!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared SKALA/GauXC section and model selection helpers
! **************************************************************************************************
MODULE skala_model_selection

   USE input_constants,                 ONLY: do_gpu_auto,&
                                              do_gpu_device,&
                                              do_gpu_host
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length
   USE offload_api,                     ONLY: offload_get_chosen_device,&
                                              offload_get_device_count
   USE skala_gpw_features,              ONLY: skala_execution_target_type,&
                                              skala_state
   USE string_utilities,                ONLY: uppercase
   USE torch_api,                       ONLY: torch_cuda_device_count,&
                                              torch_cuda_is_available,&
                                              torch_use_cuda
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: find_single_functional_subsection, &
             gauxc_model_none_selected, get_gauxc_section, get_skala_native_section, &
             model_key_is_none, &
             resolve_skala_execution_target, resolve_skala_model_path, &
             skala_execution_target_type, &
             xc_section_uses_gauxc, xc_section_uses_gauxc_model, &
             xc_section_uses_native_skala_grid, xc_section_uses_skala_model

CONTAINS

! **************************************************************************************************
!> \brief Return the XC_FUNCTIONAL subsection with the given name, if present.
!> \param xc_section ...
!> \param section_name ...
!> \return ...
! **************************************************************************************************
   FUNCTION find_functional_subsection(xc_section, section_name) RESULT(functional_section)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      CHARACTER(len=*), INTENT(IN)                       :: section_name
      TYPE(section_vals_type), POINTER                   :: functional_section

      INTEGER                                            :: ifun
      TYPE(section_vals_type), POINTER                   :: functionals, xc_fun

      NULLIFY (functional_section)
      IF (.NOT. ASSOCIATED(xc_section)) RETURN

      functionals => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      IF (.NOT. ASSOCIATED(functionals)) RETURN

      ifun = 0
      DO
         ifun = ifun + 1
         xc_fun => section_vals_get_subs_vals2(functionals, i_section=ifun)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         IF (TRIM(xc_fun%section%name) == TRIM(section_name)) THEN
            functional_section => xc_fun
            EXIT
         END IF
      END DO

   END FUNCTION find_functional_subsection

! **************************************************************************************************
!> \brief Return the named functional subsection and the number of explicit XC_FUNCTIONAL subsections.
!> \param xc_section ...
!> \param section_name ...
!> \param functional_section ...
!> \param n_subsections ...
! **************************************************************************************************
   SUBROUTINE find_single_functional_subsection(xc_section, section_name, functional_section, &
                                                n_subsections)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      CHARACTER(len=*), INTENT(IN)                       :: section_name
      TYPE(section_vals_type), POINTER                   :: functional_section
      INTEGER, INTENT(OUT)                               :: n_subsections

      INTEGER                                            :: ifun
      TYPE(section_vals_type), POINTER                   :: functionals, xc_fun

      NULLIFY (functional_section)
      n_subsections = 0
      IF (.NOT. ASSOCIATED(xc_section)) RETURN

      functionals => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      IF (.NOT. ASSOCIATED(functionals)) RETURN

      ifun = 0
      DO
         ifun = ifun + 1
         xc_fun => section_vals_get_subs_vals2(functionals, i_section=ifun)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         n_subsections = n_subsections + 1
         IF (TRIM(xc_fun%section%name) == TRIM(section_name)) functional_section => xc_fun
      END DO

   END SUBROUTINE find_single_functional_subsection

! **************************************************************************************************
!> \brief Whether a functional subsection is present and not switched off with
!>        its _SECTION_PARAMETERS_ keyword (e.g. "&GAUXC OFF").
!> \param section ...
!> \return ...
! **************************************************************************************************
   FUNCTION section_is_active(section) RESULT(active)
      TYPE(section_vals_type), POINTER                   :: section
      LOGICAL                                            :: active

      active = .FALSE.
      IF (.NOT. ASSOCIATED(section)) RETURN
      CALL section_vals_val_get(section, "_SECTION_PARAMETERS_", l_val=active)

   END FUNCTION section_is_active

! **************************************************************************************************
!> \brief Return the active CP2K-native SKALA functional subsection, if any.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_skala_native_section(xc_section) RESULT(native_section)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      TYPE(section_vals_type), POINTER                   :: native_section

      native_section => find_functional_subsection(xc_section, "SKALA_NATIVE")
      IF (.NOT. section_is_active(native_section)) NULLIFY (native_section)

   END FUNCTION get_skala_native_section

! **************************************************************************************************
!> \brief Return the first active GAUXC functional subsection, if any.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_gauxc_section(xc_section) RESULT(gauxc_section)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      gauxc_section => find_functional_subsection(xc_section, "GAUXC")
      IF (.NOT. section_is_active(gauxc_section)) NULLIFY (gauxc_section)

   END FUNCTION get_gauxc_section

! **************************************************************************************************
!> \brief Return true if an active GAUXC or SKALA_NATIVE subsection is present.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_section_uses_gauxc(xc_section) RESULT(uses_gauxc)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: uses_gauxc

      uses_gauxc = ASSOCIATED(get_gauxc_section(xc_section)) .OR. &
                   ASSOCIATED(get_skala_native_section(xc_section))

   END FUNCTION xc_section_uses_gauxc

! **************************************************************************************************
!> \brief Return true if an active SKALA_NATIVE subsection is present.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_section_uses_native_skala_grid(xc_section) RESULT(uses_native_grid)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: uses_native_grid

      uses_native_grid = ASSOCIATED(get_skala_native_section(xc_section))

   END FUNCTION xc_section_uses_native_skala_grid

! **************************************************************************************************
!> \brief Return true if an active GAUXC or SKALA_NATIVE subsection requests a model evaluation.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_section_uses_gauxc_model(xc_section) RESULT(uses_gauxc_model)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: uses_gauxc_model

      CHARACTER(len=default_path_length)                 :: model_name, xc_name
      TYPE(section_vals_type), POINTER                   :: gauxc_section, native_section

      uses_gauxc_model = .FALSE.
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         uses_gauxc_model = .TRUE.
         RETURN
      END IF
      gauxc_section => get_gauxc_section(xc_section)
      IF (.NOT. ASSOCIATED(gauxc_section)) RETURN

      CALL section_vals_val_get(gauxc_section, "MODEL", c_val=model_name)
      CALL section_vals_val_get(gauxc_section, "FUNCTIONAL", c_val=xc_name)
      uses_gauxc_model = .NOT. model_key_is_none(model_name, xc_name)

   END FUNCTION xc_section_uses_gauxc_model

! **************************************************************************************************
!> \brief Return true if the XC section requests a SKALA model (GauXC or native/FTorch).
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_section_uses_skala_model(xc_section) RESULT(uses_skala_model)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: uses_skala_model

      CHARACTER(len=default_path_length)                 :: model_key, model_name
      TYPE(section_vals_type), POINTER                   :: gauxc_section, native_section

      uses_skala_model = .FALSE.
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         uses_skala_model = .TRUE.
         RETURN
      END IF

      gauxc_section => get_gauxc_section(xc_section)
      IF (.NOT. ASSOCIATED(gauxc_section)) RETURN

      CALL section_vals_val_get(gauxc_section, "MODEL", c_val=model_name)
      model_key = ADJUSTL(model_name)
      CALL uppercase(model_key)
      uses_skala_model = (TRIM(model_key) == "SKALA")

   END FUNCTION xc_section_uses_skala_model

! **************************************************************************************************
!> \brief Return true when a model key selects NONE or the functional itself.
!> \param model_name ...
!> \param xc_name ...
!> \return ...
! **************************************************************************************************
   FUNCTION model_key_is_none(model_name, xc_name) RESULT(none_selected)
      CHARACTER(len=*), INTENT(IN)                       :: model_name, xc_name
      LOGICAL                                            :: none_selected

      CHARACTER(len=default_path_length)                 :: model_key, xc_key

      model_key = ADJUSTL(model_name)
      xc_key = ADJUSTL(xc_name)
      CALL uppercase(model_key)
      CALL uppercase(xc_key)
      none_selected = (TRIM(model_key) == "" .OR. TRIM(model_key) == "NONE" .OR. &
                       TRIM(model_key) == TRIM(xc_key))

   END FUNCTION model_key_is_none

! **************************************************************************************************
!> \brief Return true when MODEL NONE (or no model distinct from FUNCTIONAL) is selected.
!> \param functional the GAUXC section
!> \return true when no distinct model is selected
! **************************************************************************************************
   FUNCTION gauxc_model_none_selected(functional)
      TYPE(section_vals_type), POINTER                   :: functional
      LOGICAL                                            :: gauxc_model_none_selected

      CHARACTER(len=default_path_length)                 :: model_name, xc_fun_name

      CALL section_vals_val_get(functional, "MODEL", c_val=model_name)
      CALL section_vals_val_get(functional, "FUNCTIONAL", c_val=xc_fun_name)
      gauxc_model_none_selected = model_key_is_none(model_name, xc_fun_name)

   END FUNCTION gauxc_model_none_selected

! **************************************************************************************************
!> \brief Resolve the CPU/GPU execution target once for the whole calculation.
!>        CUDA is opt-out: on a CUDA-capable backend AUTO/DEVICE select the GPU.
!>        HOST always selects the CPU. A CUDA-capable backend that cannot reach a
!>        usable GPU aborts with a hint to set USE_GPU HOST. The decision is
!>        memoized so later calls return the same target, and the Torch wrapper is
!>        configured once.
!> \param gpu_mode the USE_GPU selection
!> \param backend_has_gpu whether this backend was built with CUDA support
!> \param backend_name backend name for diagnostics
!> \param exec_target the resolved execution target
! **************************************************************************************************
   SUBROUTINE resolve_skala_execution_target(gpu_mode, backend_has_gpu, backend_name, exec_target)
      INTEGER, INTENT(IN)                                :: gpu_mode
      LOGICAL, INTENT(IN)                                :: backend_has_gpu
      CHARACTER(len=*), INTENT(IN)                       :: backend_name
      TYPE(skala_execution_target_type), INTENT(OUT)     :: exec_target

      CHARACTER(len=256)                                 :: msg

      exec_target%gpu_mode = gpu_mode
      exec_target%use_cuda = .FALSE.
      exec_target%device = -1

      SELECT CASE (gpu_mode)
      CASE (do_gpu_host)
         exec_target%use_cuda = .FALSE.
      CASE (do_gpu_device)
         IF (.NOT. backend_has_gpu) THEN
            msg = TRIM(backend_name)//": USE_GPU GPU requested, but this backend "// &
                  "was not built with CUDA support. Set USE_GPU HOST to run on the CPU."
            CPABORT(TRIM(msg))
         END IF
         CALL select_cuda_target(backend_name, exec_target)
      CASE (do_gpu_auto)
         IF (backend_has_gpu) THEN
            CALL select_cuda_target(backend_name, exec_target)
         END IF
      CASE DEFAULT
         CPABORT("Unknown USE_GPU execution-space value.")
      END SELECT

      IF (skala_state%execution_resolved) THEN
         IF (skala_state%execution_target%use_cuda .NEQV. exec_target%use_cuda) THEN
            msg = "Conflicting CPU/GPU execution targets across backends. "// &
                  "All backends must agree; set USE_GPU HOST consistently."
            CPABORT(TRIM(msg))
         END IF
         exec_target = skala_state%execution_target
      ELSE
         skala_state%execution_target = exec_target
         skala_state%execution_resolved = .TRUE.
         CALL torch_use_cuda(exec_target%use_cuda)
      END IF
   END SUBROUTINE resolve_skala_execution_target

! **************************************************************************************************
!> \brief Validate that a usable CUDA device exists and record it in the target.
!> \param backend_name backend name for diagnostics
!> \param exec_target execution target to update
! **************************************************************************************************
   SUBROUTINE select_cuda_target(backend_name, exec_target)
      CHARACTER(len=*), INTENT(IN)                       :: backend_name
      TYPE(skala_execution_target_type), INTENT(INOUT)   :: exec_target

      CHARACTER(len=256)                                 :: msg
      INTEGER                                            :: device_count

      IF (.NOT. torch_cuda_is_available()) THEN
         msg = TRIM(backend_name)//": GPU execution is the default for this CUDA "// &
               "build, but no usable CUDA device was found. Set USE_GPU HOST to run on the CPU."
         CPABORT(TRIM(msg))
      END IF
      device_count = torch_cuda_device_count()
      IF (device_count <= 0) THEN
         msg = TRIM(backend_name)//": GPU execution is the default for this CUDA "// &
               "build, but no usable CUDA device was found. Set USE_GPU HOST to run on the CPU."
         CPABORT(TRIM(msg))
      END IF
      IF (offload_get_device_count() <= 0) THEN
         msg = TRIM(backend_name)//": GPU execution was requested, but the offload "// &
               "library reports no usable device. Check CUDA_VISIBLE_DEVICES and the "// &
               "MPI rank-to-device mapping, or set USE_GPU HOST."
         CPABORT(TRIM(msg))
      END IF
      exec_target%device = offload_get_chosen_device()
      IF (exec_target%device < 0 .OR. exec_target%device >= device_count) THEN
         msg = TRIM(backend_name)//": GPU execution was requested, but the offload "// &
               "library did not select a device visible to Torch. Check "// &
               "CUDA_VISIBLE_DEVICES and the MPI rank-to-device mapping, or set USE_GPU HOST."
         CPABORT(TRIM(msg))
      END IF
      exec_target%use_cuda = .TRUE.
   END SUBROUTINE select_cuda_target

! **************************************************************************************************
!> \brief Resolve MODEL SKALA to the host or CUDA model path, or return a custom path unchanged.
!> \param model_name ...
!> \param use_cuda ...
!> \param model_path ...
! **************************************************************************************************
   SUBROUTINE resolve_skala_model_path(model_name, use_cuda, model_path)
      CHARACTER(len=*), INTENT(IN)                       :: model_name
      LOGICAL, INTENT(IN)                                :: use_cuda
      CHARACTER(len=*), INTENT(OUT)                      :: model_path

      CHARACTER(len=default_path_length)                 :: model_key
      INTEGER                                            :: env_status

      model_key = ADJUSTL(model_name)
      CALL uppercase(model_key)
      IF (TRIM(model_key) == "NONE" .OR. TRIM(model_key) == "") THEN
         CPABORT("SKALA requires MODEL SKALA or a TorchScript model path")
      ELSE IF (TRIM(model_key) == "SKALA") THEN
         IF (use_cuda) THEN
            CALL GET_ENVIRONMENT_VARIABLE("SKALA_CUDA_MODEL", model_path, STATUS=env_status)
            IF (env_status /= 0 .OR. LEN_TRIM(model_path) == 0) THEN
               CPABORT("GPU execution of SKALA requires SKALA_CUDA_MODEL")
            END IF
         ELSE
            CALL GET_ENVIRONMENT_VARIABLE("SKALA_MODEL", model_path, STATUS=env_status)
            IF (env_status /= 0 .OR. LEN_TRIM(model_path) == 0) THEN
               CPABORT("SKALA requires the SKALA_MODEL environment variable")
            END IF
         END IF
      ELSE
         model_path = model_name
      END IF

   END SUBROUTINE resolve_skala_model_path

END MODULE skala_model_selection
