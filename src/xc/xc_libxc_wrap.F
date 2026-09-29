!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Includes all necessary routines, functions and parameters from
!>        libxc. Provides CP2K routines/functions where the LibXC calling list
!>        is version dependent (>=4.0.3). The naming convention for such
!>        routines/functions is xc_f03_XXX --> 'xc_libxc_wrap_XXX'. All version
!>        independent routines/functions are just bypassed to higher level
!>        module file 'xc_libxc'.
!>
!> \par History
!>      08.2015 created [A. Gloess (agloess)]
!>      01.2018 refactoring [A. Gloess (agloess)]
!>      10.2018/04.2019 added hyb_mgga [S. Simko, included by F. Stein]
!> \author A. Gloess (agloess)
! **************************************************************************************************
#:include 'xc_gamma.fypp'

#! Every staging buffer of a worker, so that the component count and the
#! release cannot fall out of step with the allocation.
#:set worker_1d = ["exc", "nd", "nda", "ndb"]
#:set worker_2d = ["rho", "sigma", "lapl", "tau", "vrho", "vsigma", "vlapl", "vtau", &
                 & "v2rho2", "v2rhosigma", "v2sigma2", "v2rholapl", "v2rhotau", &
                 & "v2sigmalapl", "v2sigmatau", "v2lapl2", "v2lapltau", "v2tau2", &
                 & "v3rho3", "v3rho2sigma", "v3rhosigma2", "v3sigma3"] + v3_new

MODULE xc_libxc_wrap
#if defined (__LIBXC)
#include <xc_version.h>
! check for LibXC version - updated for LibXC 7.x
#if (XC_MAJOR_VERSION < 7)
#error "This version of CP2K requires libxc version 7.x or above for full functionality. Some features may not be available with older versions."
#else
   USE kinds, ONLY: dp
   USE input_constants, ONLY: do_libxc_cpu, do_libxc_gpu
   USE xc_f03_lib_m, ONLY: xc_f03_func_end, &
                           xc_f03_func_init, &
                           xc_f03_functional_get_name, &
                           xc_f03_func_set_ext_params, &
                           xc_f03_func_set_dens_threshold, &
                           xc_f03_func_set_tau_threshold, &
                           xc_f03_version_string, &
                           xc_f03_reference, &
                           xc_f03_reference_doi, &
                           xc_f03_functional_get_number, &
                           xc_f03_available_functional_numbers, &
                           xc_f03_available_functional_names, &
                           xc_f03_maximum_name_length, &
                           xc_f03_number_of_functionals, &
                           !
                           xc_f03_gga_exc, &
                           xc_f03_gga_exc_vxc, &
                           xc_f03_gga_exc_vxc_fxc, &
                           xc_f03_gga_exc_vxc_fxc_kxc, &
                           xc_f03_mgga_exc_vxc_fxc_kxc, &
                           xc_f03_mgga_vxc_fxc_kxc, &
                           xc_f03_gga_fxc, &
                           xc_f03_gga_vxc, &
                           xc_f03_gga_vxc_fxc, &
                           xc_f03_gga_vxc_fxc_kxc, &
                           !
                           xc_f03_func_get_info, &
                           xc_f03_func_info_get_family, &
                           xc_f03_func_info_get_kind, &
                           xc_f03_func_info_get_name, &
                           xc_f03_func_info_get_references, &
                           xc_f03_func_info_get_flags, &
                           xc_f03_func_info_get_n_ext_params, &
                           xc_f03_func_info_get_ext_params_name, &
                           xc_f03_func_info_get_ext_params_default_value, &
                           xc_f03_func_info_get_ext_params_description, &
                           !
                           xc_f03_func_reference_get_ref, &
                           xc_f03_func_reference_get_doi, &
                           !
                           xc_f03_lda_exc_vxc_fxc_kxc, &
                           xc_f03_lda_exc, &
                           xc_f03_lda_exc_vxc, &
                           xc_f03_lda_exc_vxc_fxc, &
                           xc_f03_lda_fxc, &
                           xc_f03_lda_kxc, &
                           xc_f03_lda_vxc, &
                           xc_f03_lda_vxc_fxc, &
                           xc_f03_lda_vxc_fxc_kxc, &
                           !
                           xc_f03_mgga_exc_vxc_fxc, &
                           xc_f03_mgga_exc, &
                           xc_f03_mgga_exc_vxc, &
                           xc_f03_mgga_fxc, &
                           xc_f03_mgga_vxc, &
                           xc_f03_mgga_vxc_fxc, &
                           !
                           xc_f03_func_t, &
                           xc_f03_func_info_t, &
                           xc_f03_func_reference_t, &
                           !
                           XC_FAMILY_LDA, &
                           XC_FAMILY_GGA, &
                           XC_FAMILY_MGGA, &
                           XC_FAMILY_HYB_LDA, &
                           XC_FAMILY_HYB_GGA, &
                           XC_FAMILY_HYB_MGGA, &
                           !
                           XC_UNPOLARIZED, &
                           XC_POLARIZED, &
                           !
                           XC_EXCHANGE, &
                           XC_CORRELATION, &
                           XC_EXCHANGE_CORRELATION, &
                           XC_KINETIC, &
                           !
                           XC_FLAGS_NEEDS_LAPLACIAN, &
                           XC_FLAGS_HAVE_EXC, &
                           XC_FLAGS_DEVELOPMENT, &
                           XC_FLAGS_ON_HOST, &
                           XC_FLAGS_ON_DEVICE, &
                           xc_f03_func_info_set_default_flags

#include "../base/base_uses.f90"
   USE ISO_C_BINDING, ONLY: C_PTR, C_NULL_PTR, C_SIZE_T, C_ASSOCIATED, C_LOC
#endif
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_LIBXC)
   USE offload_api, ONLY: offload_buffer_type
#endif
   IMPLICIT NONE
   PRIVATE
#if defined (__LIBXC)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_libxc_wrap'

   CHARACTER(LEN=*), PARAMETER, PUBLIC :: libxc_version = XC_VERSION

   PUBLIC :: xc_f03_func_t, xc_f03_func_info_t
   PUBLIC :: xc_f03_func_init, xc_f03_func_end, xc_f03_func_set_ext_params
   PUBLIC :: xc_f03_functional_get_name, xc_f03_available_functional_numbers, xc_f03_maximum_name_length, &
             xc_f03_number_of_functionals, xc_f03_available_functional_names
   PUBLIC :: xc_f03_func_get_info, xc_f03_func_info_get_family, xc_f03_func_info_get_kind, &
             xc_f03_func_info_get_name, xc_f03_func_info_get_ext_params_name, &
             xc_f03_func_info_get_ext_params_description, xc_f03_func_info_get_ext_params_default_value, &
             xc_f03_func_info_get_n_ext_params
   PUBLIC :: xc_f03_gga_exc, xc_f03_gga_exc_vxc, xc_f03_gga_exc_vxc_fxc, xc_f03_gga_fxc, &
             xc_f03_gga_vxc, xc_f03_gga_vxc_fxc, &
             xc_f03_gga_exc_vxc_fxc_kxc, xc_f03_gga_vxc_fxc_kxc
   PUBLIC :: xc_f03_mgga_exc_vxc_fxc_kxc, xc_f03_mgga_vxc_fxc_kxc, &
             xc_f03_mgga_exc, xc_f03_mgga_exc_vxc, xc_f03_mgga_exc_vxc_fxc, &
             xc_f03_mgga_fxc, xc_f03_mgga_vxc, xc_f03_mgga_vxc_fxc

   PUBLIC :: XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, &
             XC_FAMILY_HYB_LDA, XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA

   PUBLIC :: XC_UNPOLARIZED, XC_POLARIZED

   PUBLIC :: XC_EXCHANGE, XC_CORRELATION, XC_EXCHANGE_CORRELATION, XC_KINETIC
   PUBLIC :: XC_FLAGS_ON_HOST, XC_FLAGS_ON_DEVICE

   ! Create aliases for backward compatibility with the LibXC 6.x interface,
   ! which referenced xc_f03_lda / xc_f03_mgga as all-in-one evaluation entries.
   PUBLIC :: xc_f03_lda, xc_f03_mgga
   PUBLIC :: xc_f03_lda_exc, xc_f03_lda_exc_vxc, xc_f03_lda_exc_vxc_fxc, &
             xc_f03_lda_fxc, xc_f03_lda_kxc, xc_f03_lda_vxc

   ! Runtime backend selection. Whether a device backend can ever be selected is
   ! a compile-time property (__OFFLOAD and !__NO_OFFLOAD_LIBXC), but the flag
   ! itself is always compiled so that callers keep a single runtime branch; on
   ! CPU-only builds it can never become the GPU.
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_LIBXC)
   LOGICAL, PARAMETER, PRIVATE :: libxc_have_device = .TRUE.
#else
   LOGICAL, PARAMETER, PRIVATE :: libxc_have_device = .FALSE.
#endif
   INTEGER, SAVE, PRIVATE :: libxc_backend = MERGE(do_libxc_gpu, do_libxc_cpu, libxc_have_device)

   ! Debug switch: test the freshly evaluated block for NaNs. Kept here because
   ! the per-block evaluation now runs inside this module.
   LOGICAL, PARAMETER, PRIVATE :: libxc_check_nan = .FALSE.

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_LIBXC)
   ! Device buffer type for GPU memory management. It wraps an offload buffer so
   ! that the host-side pinned staging memory and the device allocation share a
   ! lifetime and are released exactly once, together.
   TYPE, PRIVATE :: libxc_device_buffer_type
      TYPE(offload_buffer_type) :: buffer
      INTEGER(KIND=C_SIZE_T) :: size = 0
   END TYPE libxc_device_buffer_type
#endif

! wrappers for routines
   PUBLIC :: xc_libxc_wrap_info_refs, &
             xc_libxc_wrap_version, &
             xc_libxc_wrap_library_reference, &
             xc_libxc_wrap_info_needs_laplace, &
             xc_libxc_wrap_info_no_exc, &
             xc_libxc_wrap_set_thresholds, &
             xc_libxc_wrap_functional_get_number, &
             xc_libxc_wrap_is_under_development, &
             xc_libxc_get_reference_length, &
             xc_libxc_check_functional, &
             xc_libxc_wrap_init_safe, &
             libxc_set_backend, &
             libxc_func_init_wrapper, &
             libxc_gga_exc_vxc_wrapper, libxc_gga_exc_wrapper, &
             libxc_lda_exc_vxc_wrapper, libxc_lda_exc_wrapper, &
             libxc_mgga_exc_vxc_wrapper, libxc_mgga_exc_wrapper, &
             libxc_mgga_vxc_wrapper
   PUBLIC :: libxc_is_gpu

   ! Worker type: one LibXC functional object plus the host staging buffers used
   ! to hand it a block of grid points, together with the device buffers mirroring
   ! them on the GPU backend. It lives in the wrapper because it owns the device
   ! buffers; CPU-only builds compile an empty device section but the same type.
   ! One worker is kept per OpenMP thread, so neither the functional object nor
   ! the buffers are rebuilt between batches.
   INTEGER, PARAMETER, PRIVATE :: libxc_max_dev_buffers = 64

   TYPE, PUBLIC :: libxc_worker_type
      TYPE(xc_f03_func_t)                                :: func = xc_f03_func_t()
      TYPE(xc_f03_func_info_t)                           :: info = xc_f03_func_info_t()
      LOGICAL                                            :: is_init = .FALSE.
      ! number of points the staging buffers below were allocated for
      INTEGER                                            :: nb = 0
      ! staging buffers, dimensioned (component, point)
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: rho, sigma, lapl, tau
      ! the floored gradient norms the spin-polarized chain rule is written in
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: nd, nda, ndb
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: exc
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: vrho, vsigma, vlapl, vtau
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: v2rho2, v2rhosigma, v2sigma2, &
                                                            v2rholapl, v2rhotau, v2sigmalapl, &
                                                            v2sigmatau, v2lapl2, v2lapltau, v2tau2
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: v3rho3, v3rho2sigma, &
                                                            v3rhosigma2, v3sigma3
      #:for arr in v3_new
         REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: ${arr}$
      #:endfor
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_LIBXC)
      ! GPU device buffers for the active backend, indexed by the caller with the
      ! libxc_dev_* slot constants defined in this module.
      TYPE(libxc_device_buffer_type), DIMENSION(libxc_max_dev_buffers) :: dev
#endif
   END TYPE libxc_worker_type

   ! Slot integers naming the device buffers, used only by the host/device copy
   ! routines above; the generated third-derivative arrays continue after the
   ! hand-written ones.
   INTEGER, PARAMETER, PRIVATE :: libxc_dev_rho = 1, libxc_dev_sigma = 2, &
                                  libxc_dev_lapl = 3, libxc_dev_tau = 4, &
                                  libxc_dev_exc = 5, libxc_dev_vrho = 6, &
                                  libxc_dev_vsigma = 7, libxc_dev_vlapl = 8, &
                                  libxc_dev_vtau = 9, libxc_dev_v2rho2 = 10, &
                                  libxc_dev_v2rhosigma = 11, libxc_dev_v2sigma2 = 12, &
                                  libxc_dev_v2rholapl = 13, libxc_dev_v2rhotau = 14, &
                                  libxc_dev_v2sigmalapl = 15, libxc_dev_v2sigmatau = 16, &
                                  libxc_dev_v2lapl2 = 17, libxc_dev_v2lapltau = 18, &
                                  libxc_dev_v2tau2 = 19, libxc_dev_v3rho3 = 20, &
                                  libxc_dev_v3rho2sigma = 21, libxc_dev_v3rhosigma2 = 22, &
                                  libxc_dev_v3sigma3 = 23
   INTEGER, PARAMETER, PRIVATE :: libxc_dev_first_new3 = 24
   #:for i, arr in enumerate(v3_new)
      INTEGER, PARAMETER, PRIVATE :: libxc_dev_${arr}$ = libxc_dev_first_new3 + ${i}$
   #:endfor

   PUBLIC :: libxc_worker_alloc, libxc_worker_free_buffers, &
             libxc_worker_resize_cuda, libxc_worker_eval, libxc_family_flags

   ! Alias definitions for backward compatibility
   PROCEDURE(xc_f03_lda_exc_vxc_fxc_kxc), POINTER :: xc_f03_lda => xc_f03_lda_exc_vxc_fxc_kxc
   PROCEDURE(xc_f03_mgga_exc_vxc_fxc), POINTER :: xc_f03_mgga => xc_f03_mgga_exc_vxc_fxc

CONTAINS

   ! Backend selection: the runtime libxc_backend switch, fed by
   ! &XC/LIBXC_GPU_BACKEND (see libxc_set_backend). Whether a device backend can
   ! ever be selected is a compile-time property (__OFFLOAD and
   ! !defined(__NO_OFFLOAD_LIBXC)).
! **************************************************************************************************
!> \brief Runtime switch for the LibXC backend, fed by &XC/LIBXC_GPU_BACKEND.
!> \param backend one of do_libxc_cpu / do_libxc_gpu
!> \note Whether a device backend can ever be selected is a compile-time property
!>       (__OFFLOAD and !defined(__NO_OFFLOAD_LIBXC)). A DEVICE request on a build
!>       without a device path is an error.
! **************************************************************************************************
   SUBROUTINE libxc_set_backend(backend)
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_LIBXC)
      USE offload_api, ONLY: offload_get_device_count
#endif
      INTEGER, INTENT(IN)                                :: backend

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_LIBXC)
      ! Honor the request only when an accelerator is usable at runtime;
      ! a DEVICE request without one is an error.
      SELECT CASE (backend)
      CASE (do_libxc_gpu)
         IF (offload_get_device_count() > 0) THEN
            libxc_backend = do_libxc_gpu
         ELSE
            CALL cp__b(__SHORT_FILE__, __LINE__, &
                       "&XC/LIBXC_GPU_BACKEND requests DEVICE, "// &
                       "but no usable device/GPU was found at runtime.")
         END IF
      CASE default
         libxc_backend = do_libxc_cpu
      END SELECT
#else
      ! Without a device path a DEVICE request is an error, HOST runs.
      IF (backend == do_libxc_gpu) THEN
         CALL cp__b(__SHORT_FILE__, __LINE__, &
                    "LibXC was compiled without a device/GPU backend, "// &
                    "but &XC/LIBXC_GPU_BACKEND requests DEVICE.")
      END IF
      libxc_backend = do_libxc_cpu
#endif
   END SUBROUTINE libxc_set_backend

! **************************************************************************************************
!> \brief Set the default LibXC functional flags to host.
! **************************************************************************************************
   SUBROUTINE xc_libxc_wrap_init_safe()
      CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_wrap_init_safe'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Set default flags
      CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_HOST)

      CALL timestop(handle)
   END SUBROUTINE xc_libxc_wrap_init_safe

! **************************************************************************************************
!> \brief Provides the reference(s) for this functional.
!> \param xc_info func_info object of the functional
!> \return upper bound for the length of the reference string
!> \author F. Stein
! **************************************************************************************************
   FUNCTION xc_libxc_get_reference_length(xc_info) RESULT(length)

      TYPE(xc_f03_func_info_t), INTENT(IN)               :: xc_info
      INTEGER                                            :: length

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_get_reference_length'
      INTEGER, PARAMETER                                 :: maxlen = 67

      CHARACTER(LEN=1024)                                :: doi_string, ref_string
      CHARACTER(LEN=128)                                 :: descr_string
      INTEGER                                            :: handle, i, i_ref, i_ref_old, n_params
      TYPE(xc_f03_func_reference_t)                      :: xc_ref

      CALL timeset(routineN, handle)

      ! We are counting the number of necessary lines by carrying out a dry run of xc_libxc_wrap_info_refs
      i_ref = 0
      i_ref_old = -1
      length = 0
      DO WHILE (i_ref >= 0 .AND. i_ref /= i_ref_old)
         ! information about functional references
         xc_ref = xc_f03_func_info_get_references(xc_info, i_ref)
         ref_string = xc_f03_func_reference_get_ref(xc_ref)
         doi_string = xc_f03_func_reference_get_doi(xc_ref)
         length = length + LEN_TRIM(ref_string) + LEN_TRIM(doi_string) + 11
         IF (MOD(length, maxlen) /= 0) length = length + maxlen - MOD(length, maxlen)
         ! information about (optional) external parameters
         n_params = xc_f03_func_info_get_n_ext_params(xc_info)
         IF (n_params > 0) THEN
            length = length + maxlen
         END IF
         DO i = 1, n_params
            descr_string = xc_f03_func_info_get_ext_params_description(xc_info, i - 1)
            length = length + LEN_TRIM(descr_string) + 3
            IF (MOD(length, maxlen) /= 0) length = length + maxlen - MOD(length, maxlen)
         END DO
         i_ref_old = i_ref
      END DO
      ! two additional lines for spin polarization, scaling factor and buffer
      length = length + 2*maxlen

      CALL timestop(handle)

   END FUNCTION xc_libxc_get_reference_length

! **************************************************************************************************
!> \brief Provides the reference(s) for this functional.
!> \param xc_info ...
!> \param polarized ...
!> \param sc ...
!> \param reference ...
!>
!> \author A. Gloess (agloess)
! **************************************************************************************************
   SUBROUTINE xc_libxc_wrap_info_refs(xc_info, polarized, sc, reference)
      TYPE(xc_f03_func_info_t), INTENT(IN)               :: xc_info
      INTEGER, INTENT(IN)                                :: polarized
      REAL(KIND=dp), INTENT(IN)                          :: sc
      CHARACTER(LEN=*), INTENT(OUT)                      :: reference

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_wrap_info_refs'
      INTEGER, PARAMETER                                 :: maxlen = 67

      CHARACTER(LEN=1028)                                :: doi_string, ref_string
      CHARACTER(LEN=128)                                 :: descr_string
      CHARACTER(LEN=2067)                                :: tmp_string
      INTEGER                                            :: empty, first, handle, i, i_ref, &
                                                            i_ref_old, idx, last, n_params
      TYPE(xc_f03_func_reference_t)                      :: xc_ref

! conservative estimate of the necessary length: 2*1028+11=2067

      CALL timeset(routineN, handle)

      i_ref = 0
      i_ref_old = -1
      idx = 1
      first = 1
      DO WHILE (i_ref >= 0 .AND. i_ref /= i_ref_old)
         ! information about functional references
         xc_ref = xc_f03_func_info_get_references(xc_info, i_ref)
         ref_string = xc_f03_func_reference_get_ref(xc_ref)
         doi_string = xc_f03_func_reference_get_doi(xc_ref)
         WRITE (tmp_string, '(a1,i1,a2,a,a7,a)') '[', idx, '] ', &
            TRIM(ref_string), ', doi: ', TRIM(doi_string)
         last = first + LEN_TRIM(tmp_string) - 1
         reference(first:last) = TRIM(tmp_string)
         first = last + 1
         empty = last + (maxlen - 1) - MOD(last - 1, maxlen)
         ! fill up line with 'spaces'
         IF (empty /= last) THEN
            reference(first:empty) = ' '
            first = empty + 1
         END IF
         ! information about (optional) external parameters
         n_params = xc_f03_func_info_get_n_ext_params(xc_info)
         IF (n_params > 0) THEN
            reference(first:first + maxlen - 1) = 'Optional external parameters:'//REPEAT(' ', maxlen - 28)
            first = first + maxlen
         END IF
         DO i = 1, n_params
            descr_string = xc_f03_func_info_get_ext_params_description(xc_info, i - 1)
            last = first + LEN_TRIM(descr_string) - 1 + 3
            reference(first:last) = ' * '//TRIM(descr_string)
            first = last + 1
            empty = last + (maxlen - 1) - MOD(last - 1, maxlen)
            ! fill up line with 'spaces'

            IF (empty /= last) THEN
               reference(first:empty) = ' '
               first = empty + 1
            END IF
         END DO
         idx = idx + 1
         i_ref_old = i_ref
      END DO
      SELECT CASE (polarized)
      CASE (XC_UNPOLARIZED)
         WRITE (tmp_string, "('{scale=',f5.3,', spin-unpolarized}')") sc
      CASE (XC_POLARIZED)
         WRITE (tmp_string, "('{scale=',f5.3,', spin-polarized}')") sc
      CASE default
         CPABORT("Unsupported value for variable 'polarized'.")
      END SELECT
      last = first + LEN_TRIM(tmp_string) - 1
      reference(first:last) = TRIM(tmp_string)
      first = last + 1
      ! fill with 'spaces'
      reference(first:LEN(reference)) = ' '

      IF (last > LEN(reference)) &
         CPABORT("Faulty reference length.")

      CALL timestop(handle)

   END SUBROUTINE xc_libxc_wrap_info_refs

! **************************************************************************************************
!> \brief Provides the version of the LibXC library CP2K is running against.
!> \param version version string reported by the library at run time
!> \param compiled_version version of the headers CP2K was compiled against
!> \author A. Gloess (agloess)
!> \note The run-time and compile-time versions differ whenever CP2K is executed
!>       against a shared LibXC other than the one it was built with. Reporting
!>       the run-time version is what actually describes the numbers produced.
!>
! **************************************************************************************************
   SUBROUTINE xc_libxc_wrap_version(version, compiled_version)
      CHARACTER(LEN=*), INTENT(OUT)                      :: version
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: compiled_version

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_wrap_version'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL xc_f03_version_string(version)
      IF (PRESENT(compiled_version)) compiled_version = TRIM(libxc_version)

      CALL timestop(handle)

   END SUBROUTINE xc_libxc_wrap_version

! **************************************************************************************************
!> \brief Provides the citation LibXC asks to be given for the library itself.
!> \param reference bibliographic reference
!> \param doi digital object identifier of that reference
!> \author S. Lehtola
!> \note Taken from the library at run time rather than hard-coded, so that the
!>       citation stays correct when LibXC updates its own self-citation.
!>
! **************************************************************************************************
   SUBROUTINE xc_libxc_wrap_library_reference(reference, doi)
      CHARACTER(LEN=*), INTENT(OUT)                      :: reference, doi

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_wrap_library_reference'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL xc_f03_reference(reference)
      CALL xc_f03_reference_doi(doi)

      CALL timestop(handle)

   END SUBROUTINE xc_libxc_wrap_library_reference

! **************************************************************************************************
!> \brief Applies the CP2K screening cutoffs to a LibXC functional object.
!> \param xc_func LibXC functional object
!> \param xc_info info object belonging to xc_func
!> \param epsilon_rho density cutoff (XC%DENSITY_CUTOFF)
!> \param epsilon_tau kinetic energy density cutoff (XC%TAU_CUTOFF)
!> \author S. Lehtola
!> \note LibXC screens every point against the density threshold internally, which
!>       is what allows CP2K to hand it whole blocks of points instead of filtering
!>       point by point. The sigma threshold is deliberately left at LibXC's own
!>       default: CP2K has no cutoff of its own on the gradient here, and raising
!>       it in step with the density cutoff would floor sigma for points with a
!>       small gradient that used to be evaluated as given.
!>
! **************************************************************************************************
   SUBROUTINE xc_libxc_wrap_set_thresholds(xc_func, xc_info, epsilon_rho, epsilon_tau)
      TYPE(xc_f03_func_t), INTENT(INOUT)                 :: xc_func
      TYPE(xc_f03_func_info_t), INTENT(IN)               :: xc_info
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_rho, epsilon_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_wrap_set_thresholds'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL xc_f03_func_set_dens_threshold(xc_func, epsilon_rho)

      SELECT CASE (xc_f03_func_info_get_family(xc_info))
      CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
         CALL xc_f03_func_set_tau_threshold(xc_func, epsilon_tau)
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE xc_libxc_wrap_set_thresholds

! **************************************************************************************************
!> \brief Checks existence of functional in LibXC
!> \param func_string ...
!> \return ...
!> \author F. Stein
!> \note Remove prefix to keep compatibility, functionals can be specified (in
!>       LIBXC section) as:
!>       GGA_X_...  or  XC_GGA_X_...
!>       Starting from version 2.2.0 both name conventions are allowed, before
!>       the 'XC_' prefix was necessary.
!>
! **************************************************************************************************
   LOGICAL FUNCTION xc_libxc_check_functional(func_string) RESULT(exists)
      CHARACTER(LEN=*), INTENT(IN)                       :: func_string

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_check_functional'

      INTEGER                                            :: func_id, handle

      CALL timeset(routineN, handle)

      IF (func_string(1:3) == "XC_") THEN
         func_id = xc_f03_functional_get_number(func_string(4:LEN_TRIM(func_string)))
      ELSE
         func_id = xc_f03_functional_get_number(func_string(1:LEN_TRIM(func_string)))
      END IF

      exists = .TRUE.
      IF (func_id == -1) exists = .FALSE.

      CALL timestop(handle)

   END FUNCTION xc_libxc_check_functional

! **************************************************************************************************
!> \brief Provides the functional ID.
!> \param func_string ...
!> \return ...
!> \author A. Gloess (agloess)
!> \note Remove prefix to keep compatibility, functionals can be specified (in
!>       LIBXC section) as:
!>       GGA_X_...  or  XC_GGA_X_...
!>       Starting from version 2.2.0 both name conventions are allowed, before
!>       the 'XC_' prefix was necessary.
!>
! **************************************************************************************************
   INTEGER FUNCTION xc_libxc_wrap_functional_get_number(func_string) RESULT(func_id)
      CHARACTER(LEN=*), INTENT(IN)                       :: func_string

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_wrap_functional_get_number'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (func_string(1:3) == "XC_") THEN
         func_id = xc_f03_functional_get_number(func_string(4:LEN_TRIM(func_string)))
      ELSE
         func_id = xc_f03_functional_get_number(func_string(1:LEN_TRIM(func_string)))
      END IF
      IF (func_id == -1) THEN
         CPABORT(TRIM(func_string)//": wrong functional name")
      END IF

      CALL timestop(handle)

   END FUNCTION xc_libxc_wrap_functional_get_number

! **************************************************************************************************
!> \brief Wrapper to test wether functional is considered under development in Libxc
!> \param xc_info ...
!>
!> \return ...
!> \author F. Stein (fstein93)
! **************************************************************************************************
   LOGICAL FUNCTION xc_libxc_wrap_is_under_development(xc_info)
      TYPE(xc_f03_func_info_t)                           :: xc_info

      IF (IAND(xc_f03_func_info_get_flags(xc_info), XC_FLAGS_DEVELOPMENT) == XC_FLAGS_DEVELOPMENT) THEN
         xc_libxc_wrap_is_under_development = .TRUE.
      ELSE
         xc_libxc_wrap_is_under_development = .FALSE.
      END IF

   END FUNCTION xc_libxc_wrap_is_under_development

! **************************************************************************************************
!> \brief Wrapper for functionals that need the Laplacian, all others can use
!>        a dummy array.
!> \param xc_info info object of an already initialized functional
!>
!> \return ...
!> \author A. Gloess (agloess)
!> \note Reads the flag off an info object the caller already holds. Creating a
!>       throw-away functional object just to look at one flag used to cost an
!>       extra xc_f03_func_init/xc_f03_func_end pair per evaluated batch.
! **************************************************************************************************
   LOGICAL FUNCTION xc_libxc_wrap_info_needs_laplace(xc_info)
      ! Only some MGGA functionals needs the laplacian
      TYPE(xc_f03_func_info_t), INTENT(IN)               :: xc_info

      IF (IAND(xc_f03_func_info_get_flags(xc_info), XC_FLAGS_NEEDS_LAPLACIAN) == XC_FLAGS_NEEDS_LAPLACIAN) THEN
         xc_libxc_wrap_info_needs_laplace = .TRUE.
      ELSE
         xc_libxc_wrap_info_needs_laplace = .FALSE.
      END IF

   END FUNCTION xc_libxc_wrap_info_needs_laplace

! **************************************************************************************************
!> \brief Whether the functional does not provide the energy density.
!> \param xc_info info object of an already initialized functional
!>
!> \return .TRUE. when LibXC cannot supply Exc for this functional
!> \author S. Lehtola
! **************************************************************************************************
   LOGICAL FUNCTION xc_libxc_wrap_info_no_exc(xc_info)
      TYPE(xc_f03_func_info_t), INTENT(IN)               :: xc_info

      IF (IAND(xc_f03_func_info_get_flags(xc_info), XC_FLAGS_HAVE_EXC) == XC_FLAGS_HAVE_EXC) THEN
         xc_libxc_wrap_info_no_exc = .FALSE.
      ELSE
         xc_libxc_wrap_info_no_exc = .TRUE.
      END IF

   END FUNCTION xc_libxc_wrap_info_no_exc

! **************************************************************************************************
!> \brief Initialize a LibXC functional honoring the runtime backend.
!> \param func the functional object to initialize
!> \param functional LibXC functional number
!> \param nspin XC_UNPOLARIZED or XC_POLARIZED (defaults to XC_UNPOLARIZED)
!> \param force_host initialize on the host even when the runtime backend is the
!>        device (used for query/metadata functional objects)
!> \note The functional's default flags are set to XC_FLAGS_ON_DEVICE when the
!>       device backend is active (and the initialization is not forced to the
!>       host), so that later evaluation calls run on the device.
! **************************************************************************************************
   SUBROUTINE libxc_func_init_wrapper(func, functional, nspin, force_host)
      TYPE(xc_f03_func_t), INTENT(OUT)                   :: func
      INTEGER, INTENT(IN)                                :: functional
      INTEGER, INTENT(IN), OPTIONAL                      :: nspin
      LOGICAL, INTENT(IN), OPTIONAL                      :: force_host

      LOGICAL                                            :: my_force_host

      my_force_host = .FALSE.
      IF (PRESENT(force_host)) my_force_host = force_host

      ! Query/metadata initialization always runs on the host (force_host).
      IF (libxc_backend == do_libxc_gpu .AND. .NOT. my_force_host) THEN
         CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_DEVICE)
         IF (PRESENT(nspin)) THEN
            CALL xc_f03_func_init(func, functional, nspin)
         ELSE
            CALL xc_f03_func_init(func, functional, XC_UNPOLARIZED)
         END IF
      ELSE
         CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_HOST)
         IF (PRESENT(nspin)) THEN
            CALL xc_f03_func_init(func, functional, nspin)
         ELSE
            CALL xc_f03_func_init(func, functional, XC_UNPOLARIZED)
         END IF
      END IF
   END SUBROUTINE libxc_func_init_wrapper

#endif

! **************************************************************************************************
!> \brief Whether the runtime LibXC backend is the GPU.
!> \return ...
!> \note Always available under __LIBXC: CPU-only/__NO_OFFLOAD_LIBXC builds
!>       compile the flag but can never select the device backend at runtime.
! **************************************************************************************************
   FUNCTION libxc_is_gpu() RESULT(is_gpu)
      LOGICAL                                            :: is_gpu

      is_gpu = (libxc_backend == do_libxc_gpu)
   END FUNCTION libxc_is_gpu

! **************************************************************************************************
!> \brief Allocates the staging buffers of a single worker, sizing the block so
!>        that they stay within the cache budget.
!> \param worker the worker
!> \param family LibXC functional family
!> \param nspin XC_UNPOLARIZED or XC_POLARIZED
!> \param grad_deriv highest derivative order the buffers must hold
!> \author S. Lehtola
! **************************************************************************************************
   SUBROUTINE libxc_worker_alloc(worker, family, nspin, grad_deriv)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker
      INTEGER, INTENT(IN)                                :: family, nspin, grad_deriv

      INTEGER, PARAMETER                                 :: libxc_block_bytes = 256*1024
      INTEGER, PARAMETER                                 :: libxc_block_min = 64
      INTEGER, PARAMETER                                 :: libxc_block_max = 512

      INTEGER                                            :: nb, ncomp

      ! How many doubles a point occupies spans a factor of twenty over the
      ! families, spin cases and derivative orders, so a block length fixed in
      ! points would size the buffers anywhere between 49 KiB and 950 KiB. Derive
      ! it from the footprint instead: build a one-point worker, ask it how wide
      ! a point is, and cut the real block to the byte budget. Taking the count
      ! from the buffers themselves is what stops it drifting away from the
      ! allocation it is meant to describe.
      BLOCK
         TYPE(libxc_worker_type)                         :: probe
         ! The probe only counts host components, so skip the device buffers:
         ! they would otherwise be allocated and leaked for a worker that is
         ! discarded immediately.
         CALL libxc_worker_alloc_buffers(probe, family, nspin, grad_deriv, 1, &
                                         alloc_device=.FALSE.)
         ncomp = libxc_worker_ncomp(probe)
         CALL libxc_worker_free_buffers(probe)
      END BLOCK
      nb = MIN(libxc_block_max, &
               MAX(libxc_block_min, libxc_block_bytes/(ncomp*(STORAGE_SIZE(0.0_dp)/8))))

      CALL libxc_worker_alloc_buffers(worker, family, nspin, grad_deriv, nb)

   END SUBROUTINE libxc_worker_alloc

! **************************************************************************************************
!> \brief Number of buffer components a single grid point occupies in a worker.
!> \param worker the worker
!> \return components per point, summed over the buffers it actually carries
!> \author S. Lehtola
! **************************************************************************************************
   FUNCTION libxc_worker_ncomp(worker) RESULT(ncomp)
      TYPE(libxc_worker_type), INTENT(IN)                :: worker
      INTEGER                                            :: ncomp

      ncomp = 0
      #:for arr in worker_1d
         IF (ALLOCATED(worker%${arr}$)) ncomp = ncomp + 1
      #:endfor
      #:for arr in worker_2d
         IF (ALLOCATED(worker%${arr}$)) ncomp = ncomp + SIZE(worker%${arr}$, 1)
      #:endfor

   END FUNCTION libxc_worker_ncomp

! **************************************************************************************************
!> \brief Allocates the staging buffers of a single worker for a given block length.
!> \param worker the worker
!> \param family LibXC functional family
!> \param nspin XC_UNPOLARIZED or XC_POLARIZED
!> \param grad_deriv highest derivative order the buffers must hold
!> \param nb number of grid points per block
!> \param alloc_device whether to allocate the matching GPU device buffers
!> \author S. Lehtola
!> \note The buffers follow LibXC's own layout, i.e. the leading dimension is the
!>       number of components of the quantity and the points run along the second
!>       dimension. They are allocated once per worker and reused for every batch.
!>       Only the orders this evaluation asks for are allocated: the third
!>       derivatives alone are 130 components per point for a spin-polarized
!>       meta-GGA, which a ground-state calculation would never fill.
! **************************************************************************************************
   SUBROUTINE libxc_worker_alloc_buffers(worker, family, nspin, grad_deriv, nb, alloc_device)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker
      INTEGER, INTENT(IN)                                :: family, nspin, grad_deriv, nb
      LOGICAL, INTENT(IN), OPTIONAL                      :: alloc_device

      INTEGER                                            :: nl, np, ns
      LOGICAL                                            :: pol, want_device

      pol = (nspin == XC_POLARIZED)
      want_device = .TRUE.
      IF (PRESENT(alloc_device)) want_device = alloc_device
      CALL libxc_worker_dims(nspin, np, ns, nl)

      worker%nb = nb

      ! The energy density and the first derivatives are always staged: LibXC
      ! writes them whatever higher order was asked for.
      ALLOCATE (worker%rho(np, nb))
      ALLOCATE (worker%exc(nb))
      ALLOCATE (worker%vrho(np, nb))
      IF (grad_deriv >= 2) ALLOCATE (worker%v2rho2(ns, nb))
      IF (grad_deriv >= 3) ALLOCATE (worker%v3rho3(nl, nb))

      SELECT CASE (family)
      CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
         ALLOCATE (worker%sigma(ns, nb))
         ALLOCATE (worker%vsigma(ns, nb))
         IF (pol) THEN
            ALLOCATE (worker%nd(nb), worker%nda(nb), worker%ndb(nb))
         END IF
         IF (grad_deriv >= 2) THEN
            ALLOCATE (worker%v2rhosigma(MERGE(6, 1, pol), nb))
            ALLOCATE (worker%v2sigma2(MERGE(6, 1, pol), nb))
         END IF
         IF (grad_deriv >= 3) THEN
            ALLOCATE (worker%v3rho2sigma(MERGE(9, 1, pol), nb))
            ALLOCATE (worker%v3rhosigma2(MERGE(12, 1, pol), nb))
            ALLOCATE (worker%v3sigma3(MERGE(10, 1, pol), nb))
         END IF
      END SELECT

      SELECT CASE (family)
      CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
         ALLOCATE (worker%lapl(np, nb))
         ALLOCATE (worker%tau(np, nb))
         ALLOCATE (worker%vlapl(np, nb))
         ALLOCATE (worker%vtau(np, nb))
         IF (grad_deriv >= 2) THEN
            ALLOCATE (worker%v2rholapl(nl, nb))
            ALLOCATE (worker%v2rhotau(nl, nb))
            ALLOCATE (worker%v2sigmalapl(MERGE(6, 1, pol), nb))
            ALLOCATE (worker%v2sigmatau(MERGE(6, 1, pol), nb))
            ALLOCATE (worker%v2lapl2(ns, nb))
            ALLOCATE (worker%v2lapltau(nl, nb))
            ALLOCATE (worker%v2tau2(ns, nb))
         END IF
         IF (grad_deriv >= 3) THEN
            #:for arr in v3_new
               ALLOCATE (worker%${arr}$ (MERGE(${mgga_ncomp[arr]}$, 1, pol), nb))
            #:endfor
         END IF
         ! LibXC leaves vlapl untouched for functionals that do not use the
         ! Laplacian, and CP2K then ignores it; zero it once so it never holds
         ! uninitialized memory.
         worker%vlapl = 0.0_dp
      END SELECT

      ! Allocate GPU device buffers when the GPU backend is active. They follow
      ! the host staging buffers, so the same derivative-order guards apply.
      IF (want_device .AND. libxc_is_gpu()) THEN
         CALL libxc_worker_device_alloc(worker, family, pol, nl, np, ns, nb)
      END IF

   END SUBROUTINE libxc_worker_alloc_buffers

! **************************************************************************************************
!> \brief Releases the staging buffers of a worker, keeping its functional object.
!> \param worker the worker
!> \author S. Lehtola
!> \note Frees the host buffers and the matching device buffers, but leaves the
!>       LibXC functional object untouched so that it can be reused.
! **************************************************************************************************
   SUBROUTINE libxc_worker_free_buffers(worker)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker

      #:for arr in worker_1d + worker_2d
         IF (ALLOCATED(worker%${arr}$)) DEALLOCATE (worker%${arr}$)
      #:endfor
      worker%nb = 0

      ! Freeing a buffer that was never allocated is a no-op, so this is done
      ! unconditionally rather than gating on the currently active backend.
      CALL libxc_worker_device_free(worker)

   END SUBROUTINE libxc_worker_free_buffers

! **************************************************************************************************
!> \brief Activates the rank's chosen device on the calling thread.
!> \note Needed because CUDA's current device is thread-local; a no-op without
!>       a device, so xc_libxc can call it unconditionally.
! **************************************************************************************************
   SUBROUTINE libxc_activate_device()
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_LIBXC)
      USE offload_api, ONLY: offload_activate_chosen_device

      CALL offload_activate_chosen_device()
#endif
   END SUBROUTINE libxc_activate_device

! **************************************************************************************************
!> \brief Waits for all outstanding device-to-host copies to complete.
! **************************************************************************************************
   SUBROUTINE libxc_sync_device()
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_LIBXC)
      USE offload_api, ONLY: offload_device_synchronize

      CALL offload_device_synchronize()
#endif
   END SUBROUTINE libxc_sync_device

! **************************************************************************************************
!> \brief Copies host staging data into a worker's device buffer.
!> \param worker the worker
!> \param slot device buffer to copy into (one of libxc_dev_*)
!> \param host host array to copy from
!> \param n number of elements
! **************************************************************************************************
   SUBROUTINE libxc_worker_h2d(worker, slot, host, n)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker
      INTEGER, INTENT(IN)                                :: slot
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: host
      INTEGER, INTENT(IN)                                :: n
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_LIBXC)
      CALL libxc_host_to_device(worker%dev(slot), host, n)
#else
      MARK_USED(worker)
      MARK_USED(slot)
      MARK_USED(n)
      CPABORT("libxc_worker_h2d: device backend not compiled in.")
#endif
   END SUBROUTINE libxc_worker_h2d

! **************************************************************************************************
!> \brief Issues an async device-to-host copy of one of a worker's buffers.
! **************************************************************************************************
   SUBROUTINE libxc_worker_d2h_issue(worker, slot, n)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker
      INTEGER, INTENT(IN)                                :: slot
      INTEGER, INTENT(IN)                                :: n
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_LIBXC)
      CALL libxc_d2h_issue(worker%dev(slot), n)
#else
      MARK_USED(worker)
      MARK_USED(slot)
      MARK_USED(n)
      CPABORT("libxc_worker_d2h_issue: device backend not compiled in.")
#endif
   END SUBROUTINE libxc_worker_d2h_issue

! **************************************************************************************************
!> \brief Reads one of a worker's device buffers back into host memory.
! **************************************************************************************************
   SUBROUTINE libxc_worker_d2h_read(worker, slot, host, n)
      TYPE(libxc_worker_type), INTENT(IN)                :: worker
      INTEGER, INTENT(IN)                                :: slot
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT)           :: host
      INTEGER, INTENT(IN)                                :: n
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_LIBXC)
      CALL libxc_d2h_read(host, worker%dev(slot), n)
#else
      MARK_USED(worker)
      MARK_USED(slot)
      MARK_USED(n)
      CPABORT("libxc_worker_d2h_read: device backend not compiled in.")
#endif
   END SUBROUTINE libxc_worker_d2h_read

#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_LIBXC)
! **************************************************************************************************
!> \brief Allocates (or, with grow_only, grows) one device buffer.
!> \param buffer the device buffer
!> \param n required element count
!> \param grow_only when true, leave the buffer alone if it is already large enough
! **************************************************************************************************
   SUBROUTINE libxc_dev_alloc_maybe(buffer, n, grow_only)
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER, INTENT(IN)                                :: n
      LOGICAL, INTENT(IN)                                :: grow_only

      IF (grow_only .AND. INT(buffer%size) >= n) RETURN
      CALL libxc_device_buffer_alloc(buffer, n)
   END SUBROUTINE libxc_dev_alloc_maybe

! **************************************************************************************************
!> \brief Allocates the device buffers of a worker, sizing them for n points.
!> \param worker the worker
!> \param family LibXC functional family
!> \param pol whether the functional is spin-polarized
!> \param nl, np, ns leading dimensions from libxc_worker_dims
!> \param n number of grid points the buffers must hold
!> \param grow_only when true, only grow buffers that are too small (used on the
!>        GPU path, where the buffers are regrown to the whole grid)
!> \note Only the buffers whose host counterpart is allocated are touched, so the
!>       derivative-order and family guards of libxc_worker_alloc_buffers hold.
! **************************************************************************************************
   SUBROUTINE libxc_worker_device_alloc(worker, family, pol, nl, np, ns, n, grow_only)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker
      INTEGER, INTENT(IN)                                :: family, nl, np, ns, n
      LOGICAL, INTENT(IN)                                :: pol
      LOGICAL, INTENT(IN), OPTIONAL                      :: grow_only

      INTEGER                                            :: base, i
      LOGICAL                                            :: g

      g = .FALSE.
      IF (PRESENT(grow_only)) g = grow_only

      CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_rho), np*n, g)
      CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_exc), n, g)
      CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_vrho), np*n, g)
      IF (ALLOCATED(worker%v2rho2)) CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v2rho2), ns*n, g)
      IF (ALLOCATED(worker%v3rho3)) CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v3rho3), nl*n, g)

      SELECT CASE (family)
      CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
         CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_sigma), ns*n, g)
         CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_vsigma), ns*n, g)
         IF (ALLOCATED(worker%v2rhosigma)) &
            CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v2rhosigma), MERGE(6, 1, pol)*n, g)
         IF (ALLOCATED(worker%v2sigma2)) &
            CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v2sigma2), MERGE(6, 1, pol)*n, g)
         IF (ALLOCATED(worker%v3rho2sigma)) &
            CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v3rho2sigma), MERGE(9, 1, pol)*n, g)
         IF (ALLOCATED(worker%v3rhosigma2)) &
            CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v3rhosigma2), MERGE(12, 1, pol)*n, g)
         IF (ALLOCATED(worker%v3sigma3)) &
            CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v3sigma3), MERGE(10, 1, pol)*n, g)
      END SELECT

      SELECT CASE (family)
      CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
         CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_lapl), np*n, g)
         CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_tau), np*n, g)
         CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_vlapl), np*n, g)
         CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_vtau), np*n, g)
         IF (ALLOCATED(worker%v2rholapl)) &
            CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v2rholapl), nl*n, g)
         IF (ALLOCATED(worker%v2rhotau)) &
            CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v2rhotau), nl*n, g)
         IF (ALLOCATED(worker%v2sigmalapl)) &
            CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v2sigmalapl), MERGE(6, 1, pol)*n, g)
         IF (ALLOCATED(worker%v2sigmatau)) &
            CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v2sigmatau), MERGE(6, 1, pol)*n, g)
         IF (ALLOCATED(worker%v2lapl2)) &
            CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v2lapl2), ns*n, g)
         IF (ALLOCATED(worker%v2lapltau)) &
            CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v2lapltau), nl*n, g)
         IF (ALLOCATED(worker%v2tau2)) &
            CALL libxc_dev_alloc_maybe(worker%dev(libxc_dev_v2tau2), ns*n, g)
         base = libxc_dev_first_new3
         #:for idx, arr in enumerate(v3_new)
            IF (ALLOCATED(worker%${arr}$)) THEN
               i = base + ${idx}$
               CALL libxc_dev_alloc_maybe(worker%dev(i), MERGE(${mgga_ncomp[arr]}$, 1, pol)*n, g)
            END IF
         #:endfor
      END SELECT

   END SUBROUTINE libxc_worker_device_alloc

! **************************************************************************************************
!> \brief Frees every device buffer of a worker.
! **************************************************************************************************
   SUBROUTINE libxc_worker_device_free(worker)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker

      INTEGER                                            :: i

      DO i = 1, libxc_max_dev_buffers
         CALL libxc_device_buffer_free(worker%dev(i))
      END DO

   END SUBROUTINE libxc_worker_device_free

! **************************************************************************************************
!> \brief Runs the device evaluation matching the functional family and the
!>        requested derivative order on a worker's device buffers.
!> \param worker the worker owning the device buffers
!> \param family LibXC functional family
!> \param grad_deriv derivative order requested (0..3)
!> \param no_exc functional does not provide the energy density
!> \param np number of grid points in this block
!> \note This is the only copy of the device dispatch; xc_libxc only has to
!>       stage its host buffers and then read the results back.
! **************************************************************************************************
   SUBROUTINE libxc_worker_eval_device(worker, family, grad_deriv, no_exc, np)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker
      INTEGER, INTENT(IN)                                :: family, grad_deriv
      LOGICAL, INTENT(IN)                                :: no_exc
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      LOGICAL                                            :: is_gga, is_mgga

      CALL libxc_family_flags(family, is_gga, is_mgga, "libxc_worker_eval_device")

      IF (is_mgga) THEN
         SELECT CASE (grad_deriv)
         CASE (0)
            CALL libxc_mgga_exc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_sigma), &
                                       worker%dev(libxc_dev_lapl), worker%dev(libxc_dev_tau), &
                                       worker%dev(libxc_dev_exc), np)
         CASE (1)
            IF (no_exc) THEN
               CALL libxc_mgga_vxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_sigma), &
                                          worker%dev(libxc_dev_lapl), worker%dev(libxc_dev_tau), &
                                          worker%dev(libxc_dev_vrho), worker%dev(libxc_dev_vsigma), &
                                          worker%dev(libxc_dev_vlapl), worker%dev(libxc_dev_vtau), np)
               worker%exc = 0.0_dp
            ELSE
               CALL libxc_mgga_exc_vxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_sigma), &
                                              worker%dev(libxc_dev_lapl), worker%dev(libxc_dev_tau), &
                                              worker%dev(libxc_dev_exc), worker%dev(libxc_dev_vrho), &
                                              worker%dev(libxc_dev_vsigma), worker%dev(libxc_dev_vlapl), &
                                              worker%dev(libxc_dev_vtau), np)
            END IF
         CASE (2)
            IF (no_exc) THEN
               CALL libxc_mgga_vxc_fxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_sigma), &
                                              worker%dev(libxc_dev_lapl), worker%dev(libxc_dev_tau), &
                                              worker%dev(libxc_dev_vrho), worker%dev(libxc_dev_vsigma), &
                                              worker%dev(libxc_dev_vlapl), worker%dev(libxc_dev_vtau), &
                                              worker%dev(libxc_dev_v2rho2), worker%dev(libxc_dev_v2rhosigma), &
                                              worker%dev(libxc_dev_v2rholapl), worker%dev(libxc_dev_v2rhotau), &
                                              worker%dev(libxc_dev_v2sigma2), worker%dev(libxc_dev_v2sigmalapl), &
                                              worker%dev(libxc_dev_v2sigmatau), worker%dev(libxc_dev_v2lapl2), &
                                              worker%dev(libxc_dev_v2lapltau), worker%dev(libxc_dev_v2tau2), np)
               worker%exc = 0.0_dp
            ELSE
               CALL libxc_mgga_exc_vxc_fxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_sigma), &
                                                  worker%dev(libxc_dev_lapl), worker%dev(libxc_dev_tau), &
                                                  worker%dev(libxc_dev_exc), worker%dev(libxc_dev_vrho), &
                                                  worker%dev(libxc_dev_vsigma), worker%dev(libxc_dev_vlapl), &
                                                  worker%dev(libxc_dev_vtau), worker%dev(libxc_dev_v2rho2), &
                                                  worker%dev(libxc_dev_v2rhosigma), worker%dev(libxc_dev_v2rholapl), &
                                                  worker%dev(libxc_dev_v2rhotau), worker%dev(libxc_dev_v2sigma2), &
                                                  worker%dev(libxc_dev_v2sigmalapl), worker%dev(libxc_dev_v2sigmatau), &
                                                  worker%dev(libxc_dev_v2lapl2), worker%dev(libxc_dev_v2lapltau), &
                                                  worker%dev(libxc_dev_v2tau2), np)
            END IF
         CASE (3)
            IF (no_exc) THEN
               CALL libxc_mgga_vxc_fxc_kxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_sigma), &
                                                  worker%dev(libxc_dev_lapl), worker%dev(libxc_dev_tau), &
                                                  worker%dev(libxc_dev_vrho), worker%dev(libxc_dev_vsigma), &
                                                  worker%dev(libxc_dev_vlapl), worker%dev(libxc_dev_vtau), &
                                                  worker%dev(libxc_dev_v2rho2), worker%dev(libxc_dev_v2rhosigma), &
                                                  worker%dev(libxc_dev_v2rholapl), worker%dev(libxc_dev_v2rhotau), &
                                                  worker%dev(libxc_dev_v2sigma2), worker%dev(libxc_dev_v2sigmalapl), &
                                                  worker%dev(libxc_dev_v2sigmatau), worker%dev(libxc_dev_v2lapl2), &
                                                  worker%dev(libxc_dev_v2lapltau), worker%dev(libxc_dev_v2tau2), &
                                                  worker%dev(libxc_dev_v3rho3), worker%dev(libxc_dev_v3rho2sigma), &
                                                  worker%dev(libxc_dev_v3rho2lapl), worker%dev(libxc_dev_v3rho2tau), &
                                                  worker%dev(libxc_dev_v3rhosigma2), worker%dev(libxc_dev_v3rhosigmalapl), &
                                                  worker%dev(libxc_dev_v3rhosigmatau), worker%dev(libxc_dev_v3rholapl2), &
                                                  worker%dev(libxc_dev_v3rholapltau), worker%dev(libxc_dev_v3rhotau2), &
                                                  worker%dev(libxc_dev_v3sigma3), worker%dev(libxc_dev_v3sigma2lapl), &
                                                  worker%dev(libxc_dev_v3sigma2tau), worker%dev(libxc_dev_v3sigmalapl2), &
                                                  worker%dev(libxc_dev_v3sigmalapltau), worker%dev(libxc_dev_v3sigmatau2), &
                                                  worker%dev(libxc_dev_v3lapl3), worker%dev(libxc_dev_v3lapl2tau), &
                                                  worker%dev(libxc_dev_v3lapltau2), worker%dev(libxc_dev_v3tau3), np)
               worker%exc = 0.0_dp
            ELSE
               CALL libxc_mgga_exc_vxc_fxc_kxc_device(worker%func, worker%dev(libxc_dev_rho), &
                                                      worker%dev(libxc_dev_sigma), worker%dev(libxc_dev_lapl), &
                                                      worker%dev(libxc_dev_tau), worker%dev(libxc_dev_exc), &
                                                      worker%dev(libxc_dev_vrho), worker%dev(libxc_dev_vsigma), &
                                                      worker%dev(libxc_dev_vlapl), worker%dev(libxc_dev_vtau), &
                                                      worker%dev(libxc_dev_v2rho2), worker%dev(libxc_dev_v2rhosigma), &
                                                      worker%dev(libxc_dev_v2rholapl), worker%dev(libxc_dev_v2rhotau), &
                                                      worker%dev(libxc_dev_v2sigma2), worker%dev(libxc_dev_v2sigmalapl), &
                                                      worker%dev(libxc_dev_v2sigmatau), worker%dev(libxc_dev_v2lapl2), &
                                                      worker%dev(libxc_dev_v2lapltau), worker%dev(libxc_dev_v2tau2), &
                                                      worker%dev(libxc_dev_v3rho3), worker%dev(libxc_dev_v3rho2sigma), &
                                                      worker%dev(libxc_dev_v3rho2lapl), worker%dev(libxc_dev_v3rho2tau), &
                                                      worker%dev(libxc_dev_v3rhosigma2), worker%dev(libxc_dev_v3rhosigmalapl), &
                                                      worker%dev(libxc_dev_v3rhosigmatau), worker%dev(libxc_dev_v3rholapl2), &
                                                      worker%dev(libxc_dev_v3rholapltau), worker%dev(libxc_dev_v3rhotau2), &
                                                      worker%dev(libxc_dev_v3sigma3), worker%dev(libxc_dev_v3sigma2lapl), &
                                                      worker%dev(libxc_dev_v3sigma2tau), worker%dev(libxc_dev_v3sigmalapl2), &
                                                      worker%dev(libxc_dev_v3sigmalapltau), worker%dev(libxc_dev_v3sigmatau2), &
                                                      worker%dev(libxc_dev_v3lapl3), worker%dev(libxc_dev_v3lapl2tau), &
                                                      worker%dev(libxc_dev_v3lapltau2), worker%dev(libxc_dev_v3tau3), np)
            END IF
         END SELECT
      ELSE IF (is_gga) THEN
         SELECT CASE (grad_deriv)
         CASE (0)
            CALL libxc_gga_exc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_sigma), &
                                      worker%dev(libxc_dev_exc), np)
         CASE (1)
            IF (no_exc) THEN
               CALL libxc_gga_vxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_sigma), &
                                         worker%dev(libxc_dev_vrho), worker%dev(libxc_dev_vsigma), np)
               worker%exc = 0.0_dp
            ELSE
               CALL libxc_gga_exc_vxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_sigma), &
                                             worker%dev(libxc_dev_exc), worker%dev(libxc_dev_vrho), &
                                             worker%dev(libxc_dev_vsigma), np)
            END IF
         CASE (2)
            IF (no_exc) THEN
               CALL libxc_gga_vxc_fxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_sigma), &
                                             worker%dev(libxc_dev_vrho), worker%dev(libxc_dev_vsigma), &
                                             worker%dev(libxc_dev_v2rho2), worker%dev(libxc_dev_v2rhosigma), &
                                             worker%dev(libxc_dev_v2sigma2), np)
               worker%exc = 0.0_dp
            ELSE
               CALL libxc_gga_exc_vxc_fxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_sigma), &
                                                 worker%dev(libxc_dev_exc), worker%dev(libxc_dev_vrho), &
                                                 worker%dev(libxc_dev_vsigma), worker%dev(libxc_dev_v2rho2), &
                                                 worker%dev(libxc_dev_v2rhosigma), worker%dev(libxc_dev_v2sigma2), np)
            END IF
         CASE (3)
            IF (no_exc) THEN
               CALL libxc_gga_vxc_fxc_kxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_sigma), &
                                                 worker%dev(libxc_dev_vrho), worker%dev(libxc_dev_vsigma), &
                                                 worker%dev(libxc_dev_v2rho2), worker%dev(libxc_dev_v2rhosigma), &
                                                 worker%dev(libxc_dev_v2sigma2), worker%dev(libxc_dev_v3rho3), &
                                                 worker%dev(libxc_dev_v3rho2sigma), worker%dev(libxc_dev_v3rhosigma2), &
                                                 worker%dev(libxc_dev_v3sigma3), np)
               worker%exc = 0.0_dp
            ELSE
               CALL libxc_gga_exc_vxc_fxc_kxc_device(worker%func, worker%dev(libxc_dev_rho), &
                                                     worker%dev(libxc_dev_sigma), worker%dev(libxc_dev_exc), &
                                                     worker%dev(libxc_dev_vrho), worker%dev(libxc_dev_vsigma), &
                                                     worker%dev(libxc_dev_v2rho2), worker%dev(libxc_dev_v2rhosigma), &
                                                     worker%dev(libxc_dev_v2sigma2), worker%dev(libxc_dev_v3rho3), &
                                                     worker%dev(libxc_dev_v3rho2sigma), worker%dev(libxc_dev_v3rhosigma2), &
                                                     worker%dev(libxc_dev_v3sigma3), np)
            END IF
         END SELECT
      ELSE
         SELECT CASE (grad_deriv)
         CASE (0)
            CALL libxc_lda_exc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_exc), np)
         CASE (1)
            CALL libxc_lda_exc_vxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_exc), &
                                          worker%dev(libxc_dev_vrho), np)
         CASE (2)
            CALL libxc_lda_exc_vxc_fxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_exc), &
                                              worker%dev(libxc_dev_vrho), worker%dev(libxc_dev_v2rho2), np)
         CASE (3)
            CALL libxc_lda_exc_vxc_fxc_kxc_device(worker%func, worker%dev(libxc_dev_rho), worker%dev(libxc_dev_exc), &
                                                  worker%dev(libxc_dev_vrho), worker%dev(libxc_dev_v2rho2), &
                                                  worker%dev(libxc_dev_v3rho3), np)
         END SELECT
      END IF

   END SUBROUTINE libxc_worker_eval_device

! **************************************************************************************************
!> \brief Grows a worker's staging and device buffers to hold at least npoints
!>        grid points.
!> \param worker the worker
!> \param family LibXC functional family
!> \param nspin XC_UNPOLARIZED or XC_POLARIZED
!> \param npoints number of grid points the buffers must hold
!> \author S. Lehtola
!> \note Only used on the GPU path, where the whole grid is handed to LibXC in
!>       a single call (no block cap), so the staging and device buffers have
!>       to hold the entire grid. They grow to the largest point count seen and
!>       are reused (never shrunk) for the run. Only the arrays the worker
!>       actually carries are touched, mirroring libxc_worker_alloc_buffers.
! **************************************************************************************************
   SUBROUTINE libxc_worker_resize_cuda(worker, family, nspin, npoints)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker
      INTEGER, INTENT(IN)                                :: family, nspin, npoints

      INTEGER                                            :: nl, np, ns
      LOGICAL                                            :: pol

      pol = (nspin == XC_POLARIZED)
      CALL libxc_worker_dims(nspin, np, ns, nl)

      ! grow the host staging buffers that are actually allocated
      IF (ALLOCATED(worker%rho)) THEN
         IF (SIZE(worker%rho, 2) < npoints) THEN
            DEALLOCATE (worker%rho)
            ALLOCATE (worker%rho(np, npoints))
         END IF
      END IF
      IF (ALLOCATED(worker%exc)) THEN
         IF (SIZE(worker%exc) < npoints) THEN
            DEALLOCATE (worker%exc)
            ALLOCATE (worker%exc(npoints))
         END IF
      END IF
      IF (ALLOCATED(worker%vrho)) THEN
         IF (SIZE(worker%vrho, 2) < npoints) THEN
            DEALLOCATE (worker%vrho)
            ALLOCATE (worker%vrho(np, npoints))
         END IF
      END IF
      IF (ALLOCATED(worker%v2rho2)) THEN
         IF (SIZE(worker%v2rho2, 2) < npoints) THEN
            DEALLOCATE (worker%v2rho2)
            ALLOCATE (worker%v2rho2(ns, npoints))
         END IF
      END IF
      IF (ALLOCATED(worker%v3rho3)) THEN
         IF (SIZE(worker%v3rho3, 2) < npoints) THEN
            DEALLOCATE (worker%v3rho3)
            ALLOCATE (worker%v3rho3(nl, npoints))
         END IF
      END IF

      SELECT CASE (family)
      CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
         IF (ALLOCATED(worker%sigma)) THEN
            IF (SIZE(worker%sigma, 2) < npoints) THEN
               DEALLOCATE (worker%sigma)
               ALLOCATE (worker%sigma(ns, npoints))
            END IF
         END IF
         IF (ALLOCATED(worker%vsigma)) THEN
            IF (SIZE(worker%vsigma, 2) < npoints) THEN
               DEALLOCATE (worker%vsigma)
               ALLOCATE (worker%vsigma(ns, npoints))
            END IF
         END IF
         IF (pol) THEN
            IF (ALLOCATED(worker%nd) .AND. SIZE(worker%nd) < npoints) THEN
               DEALLOCATE (worker%nd, worker%nda, worker%ndb)
               ALLOCATE (worker%nd(npoints), worker%nda(npoints), worker%ndb(npoints))
            END IF
         END IF
         IF (ALLOCATED(worker%v2rhosigma)) THEN
            IF (SIZE(worker%v2rhosigma, 2) < npoints) THEN
               DEALLOCATE (worker%v2rhosigma)
               ALLOCATE (worker%v2rhosigma(MERGE(6, 1, pol), npoints))
            END IF
         END IF
         IF (ALLOCATED(worker%v2sigma2)) THEN
            IF (SIZE(worker%v2sigma2, 2) < npoints) THEN
               DEALLOCATE (worker%v2sigma2)
               ALLOCATE (worker%v2sigma2(MERGE(6, 1, pol), npoints))
            END IF
         END IF
         #:for arr, dim in [("v3rho2sigma", "MERGE(9, 1, pol)"), &
            ("v3rhosigma2", "MERGE(12, 1, pol)"), &
            ("v3sigma3", "MERGE(10, 1, pol)")]
            IF (ALLOCATED(worker%${arr}$)) THEN
               IF (SIZE(worker%${arr}$, 2) < npoints) THEN
                  DEALLOCATE (worker%${arr}$)
                  ALLOCATE (worker%${arr}$ (${dim}$, npoints))
               END IF
            END IF
         #:endfor
      END SELECT

      SELECT CASE (family)
      CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
         #:for arr, dim in [("lapl", "np"), ("tau", "np"), ("vlapl", "np"), ("vtau", "np")]
            IF (ALLOCATED(worker%${arr}$)) THEN
               IF (SIZE(worker%${arr}$, 2) < npoints) THEN
                  DEALLOCATE (worker%${arr}$)
                  ALLOCATE (worker%${arr}$ (${dim}$, npoints))
               END IF
            END IF
         #:endfor
         ! LibXC leaves vlapl untouched for functionals that do not use the
         ! Laplacian; re-zero a freshly grown buffer as in the initial alloc.
         IF (ALLOCATED(worker%vlapl)) worker%vlapl = 0.0_dp
         #:for arr, dim in [("v2rholapl", "nl"), ("v2rhotau", "nl"), &
            ("v2sigmalapl", "MERGE(6, 1, pol)"), &
            ("v2sigmatau", "MERGE(6, 1, pol)"), &
            ("v2lapl2", "ns"), ("v2lapltau", "nl"), ("v2tau2", "ns")]
            IF (ALLOCATED(worker%${arr}$)) THEN
               IF (SIZE(worker%${arr}$, 2) < npoints) THEN
                  DEALLOCATE (worker%${arr}$)
                  ALLOCATE (worker%${arr}$ (${dim}$, npoints))
               END IF
            END IF
         #:endfor
         #:for arr in v3_new
            IF (ALLOCATED(worker%${arr}$)) THEN
               IF (SIZE(worker%${arr}$, 2) < npoints) THEN
                  DEALLOCATE (worker%${arr}$)
                  ALLOCATE (worker%${arr}$ (MERGE(${mgga_ncomp[arr]}$, 1, pol), npoints))
               END IF
            END IF
         #:endfor
      END SELECT

      ! grow only the CUDA device buffers that are already large enough to be
      ! reused; the grow-only guards keep the derivative-order and family
      ! logic of libxc_worker_alloc_buffers intact.
      CALL libxc_worker_device_alloc(worker, family, pol, nl, np, ns, npoints, grow_only=.TRUE.)

   END SUBROUTINE libxc_worker_resize_cuda

! **************************************************************************************************
!> \brief Allocate a device buffer for LibXC.
!> \param buffer ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE libxc_device_buffer_alloc(buffer, n)
      USE ISO_C_BINDING, ONLY: C_INT64_T
      USE offload_api, ONLY: offload_create_buffer
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER, INTENT(IN)                                :: n

      ! Release any previous allocation first: the GPU path regrows the device
      ! buffers to the whole grid (libxc_worker_resize_cuda), and
      ! offload_create_buffer does not free the old device pointer on its own.
      ! Freeing a buffer that was never allocated is a no-op.
      CALL libxc_device_buffer_free(buffer)

      buffer%size = INT(n, KIND=C_INT64_T)
      CALL offload_create_buffer(n, buffer%buffer)
   END SUBROUTINE libxc_device_buffer_alloc

! **************************************************************************************************
!> \brief Free a LibXC device buffer.
!> \param buffer ...
! **************************************************************************************************
   SUBROUTINE libxc_device_buffer_free(buffer)
      USE offload_api, ONLY: offload_free_buffer
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: buffer

      CALL offload_free_buffer(buffer%buffer)
      buffer%size = 0
   END SUBROUTINE libxc_device_buffer_free

! **************************************************************************************************
!> \brief Copy host data into a LibXC device buffer.
!> \param buffer ...
!> \param host ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE libxc_host_to_device(buffer, host, n)
      USE ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER
      USE offload_api, ONLY: offload_buffer_h2d, offload_get_buffer_host_pointer
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: buffer
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: host
      INTEGER, INTENT(IN)                                :: n

      REAL(KIND=dp), DIMENSION(:), POINTER               :: host_p
      TYPE(C_PTR)                                        :: hptr

      hptr = offload_get_buffer_host_pointer(buffer%buffer)
      CALL C_F_POINTER(hptr, host_p, [n])
      host_p(1:n) = host(1:n)
      CALL offload_buffer_h2d(buffer%buffer, n)
   END SUBROUTINE libxc_host_to_device
! **************************************************************************************************
!> \brief Issue a device-to-host copy of a LibXC buffer.
!> \param buffer ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE libxc_d2h_issue(buffer, n)
      USE offload_api, ONLY: offload_buffer_d2h
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER, INTENT(IN)                                :: n

      CALL offload_buffer_d2h(buffer%buffer, n)
   END SUBROUTINE libxc_d2h_issue
! **************************************************************************************************
!> \brief Read the host side of a LibXC device buffer after a copy.
!> \param host ...
!> \param buffer ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE libxc_d2h_read(host, buffer, n)
      USE ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER
      USE offload_api, ONLY: offload_get_buffer_host_pointer
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT)           :: host
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: buffer
      INTEGER, INTENT(IN)                                :: n

      REAL(KIND=dp), DIMENSION(:), POINTER               :: host_p
      TYPE(C_PTR)                                        :: hptr

      hptr = offload_get_buffer_host_pointer(buffer%buffer)
      CALL C_F_POINTER(hptr, host_p, [n])
      host(1:n) = host_p(1:n)
   END SUBROUTINE libxc_d2h_read

! **************************************************************************************************
!> \brief LDA exc and vxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param exc_d ...
!> \param vrho_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_lda_exc_vxc_device(func, rho_d, exc_d, vrho_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: exc_d, vrho_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER               :: exc_p, rho_p, vrho_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(exc_d, INT(exc_d%size), exc_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)

      CALL xc_f03_lda_exc_vxc(func, np, rho_p, exc_p, vrho_p)
   END SUBROUTINE libxc_lda_exc_vxc_device

! **************************************************************************************************
!> \brief LDA exc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param exc_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_lda_exc_device(func, rho_d, exc_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: exc_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER               :: exc_p, rho_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(exc_d, INT(exc_d%size), exc_p)

      CALL xc_f03_lda_exc(func, np, rho_p, exc_p)
   END SUBROUTINE libxc_lda_exc_device

! **************************************************************************************************
!> \brief LDA exc, vxc and fxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param exc_d ...
!> \param vrho_d ...
!> \param v2rho2_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_lda_exc_vxc_fxc_device(func, rho_d, exc_d, vrho_d, v2rho2_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: exc_d, vrho_d, v2rho2_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER               :: exc_p, rho_p, v2rho2_p, vrho_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(exc_d, INT(exc_d%size), exc_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(v2rho2_d, INT(v2rho2_d%size), v2rho2_p)

      CALL xc_f03_lda_exc_vxc_fxc(func, np, rho_p, exc_p, vrho_p, v2rho2_p)
   END SUBROUTINE libxc_lda_exc_vxc_fxc_device

! **************************************************************************************************
!> \brief LDA exc, vxc, fxc and kxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param exc_d ...
!> \param vrho_d ...
!> \param v2rho2_d ...
!> \param v3rho3_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_lda_exc_vxc_fxc_kxc_device(func, rho_d, exc_d, vrho_d, v2rho2_d, v3rho3_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: exc_d, vrho_d, v2rho2_d, v3rho3_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER               :: exc_p, rho_p, v2rho2_p, v3rho3_p, vrho_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(exc_d, INT(exc_d%size), exc_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(v2rho2_d, INT(v2rho2_d%size), v2rho2_p)
      CALL libxc_device_buffer_get_ptr(v3rho3_d, INT(v3rho3_d%size), v3rho3_p)

      CALL xc_f03_lda_exc_vxc_fxc_kxc(func, np, rho_p, exc_p, vrho_p, v2rho2_p, v3rho3_p)
   END SUBROUTINE libxc_lda_exc_vxc_fxc_kxc_device

! **************************************************************************************************
!> \brief GGA exc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param exc_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_gga_exc_device(func, rho_d, sigma_d, exc_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: exc_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER               :: exc_p, rho_p, sigma_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(exc_d, INT(exc_d%size), exc_p)

      CALL xc_f03_gga_exc(func, np, rho_p, sigma_p, exc_p)
   END SUBROUTINE libxc_gga_exc_device

! **************************************************************************************************
!> \brief GGA exc and vxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param exc_d ...
!> \param vrho_d ...
!> \param vsigma_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_gga_exc_vxc_device(func, rho_d, sigma_d, exc_d, vrho_d, vsigma_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: exc_d, vrho_d, vsigma_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER               :: exc_p, rho_p, sigma_p, vrho_p, vsigma_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(exc_d, INT(exc_d%size), exc_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(vsigma_d, INT(vsigma_d%size), vsigma_p)

      CALL xc_f03_gga_exc_vxc(func, np, rho_p, sigma_p, exc_p, vrho_p, vsigma_p)
   END SUBROUTINE libxc_gga_exc_vxc_device

! **************************************************************************************************
!> \brief GGA exc, vxc and fxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param exc_d ...
!> \param vrho_d ...
!> \param vsigma_d ...
!> \param v2rho2_d ...
!> \param v2rhosigma_d ...
!> \param v2sigma2_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_gga_exc_vxc_fxc_device(func, rho_d, sigma_d, exc_d, vrho_d, vsigma_d, &
                                           v2rho2_d, v2rhosigma_d, v2sigma2_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: exc_d, vrho_d, vsigma_d, v2rho2_d, &
                                                            v2rhosigma_d, v2sigma2_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER               :: exc_p, rho_p, sigma_p, v2rho2_p, &
                                                            v2rhosigma_p, v2sigma2_p, vrho_p, &
                                                            vsigma_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(exc_d, INT(exc_d%size), exc_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(vsigma_d, INT(vsigma_d%size), vsigma_p)
      CALL libxc_device_buffer_get_ptr(v2rho2_d, INT(v2rho2_d%size), v2rho2_p)
      CALL libxc_device_buffer_get_ptr(v2rhosigma_d, INT(v2rhosigma_d%size), v2rhosigma_p)
      CALL libxc_device_buffer_get_ptr(v2sigma2_d, INT(v2sigma2_d%size), v2sigma2_p)

      CALL xc_f03_gga_exc_vxc_fxc(func, np, rho_p, sigma_p, exc_p, vrho_p, vsigma_p, &
                                  v2rho2_p, v2rhosigma_p, v2sigma2_p)
   END SUBROUTINE libxc_gga_exc_vxc_fxc_device

! **************************************************************************************************
!> \brief GGA vxc and fxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param vrho_d ...
!> \param vsigma_d ...
!> \param v2rho2_d ...
!> \param v2rhosigma_d ...
!> \param v2sigma2_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_gga_vxc_fxc_device(func, rho_d, sigma_d, vrho_d, vsigma_d, &
                                       v2rho2_d, v2rhosigma_d, v2sigma2_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: vrho_d, vsigma_d, v2rho2_d, &
                                                            v2rhosigma_d, v2sigma2_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER               :: rho_p, sigma_p, v2rho2_p, v2rhosigma_p, &
                                                            v2sigma2_p, vrho_p, vsigma_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(vsigma_d, INT(vsigma_d%size), vsigma_p)
      CALL libxc_device_buffer_get_ptr(v2rho2_d, INT(v2rho2_d%size), v2rho2_p)
      CALL libxc_device_buffer_get_ptr(v2rhosigma_d, INT(v2rhosigma_d%size), v2rhosigma_p)
      CALL libxc_device_buffer_get_ptr(v2sigma2_d, INT(v2sigma2_d%size), v2sigma2_p)

      CALL xc_f03_gga_vxc_fxc(func, np, rho_p, sigma_p, vrho_p, vsigma_p, &
                              v2rho2_p, v2rhosigma_p, v2sigma2_p)
   END SUBROUTINE libxc_gga_vxc_fxc_device

! **************************************************************************************************
!> \brief GGA exc, vxc, fxc and kxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param exc_d ...
!> \param vrho_d ...
!> \param vsigma_d ...
!> \param v2rho2_d ...
!> \param v2rhosigma_d ...
!> \param v2sigma2_d ...
!> \param v3rho3_d ...
!> \param v3rho2sigma_d ...
!> \param v3rhosigma2_d ...
!> \param v3sigma3_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_gga_exc_vxc_fxc_kxc_device(func, rho_d, sigma_d, exc_d, vrho_d, vsigma_d, &
                                               v2rho2_d, v2rhosigma_d, v2sigma2_d, &
                                               v3rho3_d, v3rho2sigma_d, v3rhosigma2_d, v3sigma3_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: exc_d, vrho_d, vsigma_d, v2rho2_d, &
                                                            v2rhosigma_d, v2sigma2_d, v3rho3_d, &
                                                            v3rho2sigma_d, v3rhosigma2_d, &
                                                            v3sigma3_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER :: exc_p, rho_p, sigma_p, v2rho2_p, v2rhosigma_p, &
                                              v2sigma2_p, v3rho2sigma_p, v3rho3_p, v3rhosigma2_p, v3sigma3_p, vrho_p, vsigma_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(exc_d, INT(exc_d%size), exc_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(vsigma_d, INT(vsigma_d%size), vsigma_p)
      CALL libxc_device_buffer_get_ptr(v2rho2_d, INT(v2rho2_d%size), v2rho2_p)
      CALL libxc_device_buffer_get_ptr(v2rhosigma_d, INT(v2rhosigma_d%size), v2rhosigma_p)
      CALL libxc_device_buffer_get_ptr(v2sigma2_d, INT(v2sigma2_d%size), v2sigma2_p)
      CALL libxc_device_buffer_get_ptr(v3rho3_d, INT(v3rho3_d%size), v3rho3_p)
      CALL libxc_device_buffer_get_ptr(v3rho2sigma_d, INT(v3rho2sigma_d%size), v3rho2sigma_p)
      CALL libxc_device_buffer_get_ptr(v3rhosigma2_d, INT(v3rhosigma2_d%size), v3rhosigma2_p)
      CALL libxc_device_buffer_get_ptr(v3sigma3_d, INT(v3sigma3_d%size), v3sigma3_p)

      CALL xc_f03_gga_exc_vxc_fxc_kxc(func, np, rho_p, sigma_p, exc_p, vrho_p, vsigma_p, &
                                      v2rho2_p, v2rhosigma_p, v2sigma2_p, v3rho3_p, &
                                      v3rho2sigma_p, v3rhosigma2_p, v3sigma3_p)
   END SUBROUTINE libxc_gga_exc_vxc_fxc_kxc_device

! **************************************************************************************************
!> \brief GGA vxc, fxc and kxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param vrho_d ...
!> \param vsigma_d ...
!> \param v2rho2_d ...
!> \param v2rhosigma_d ...
!> \param v2sigma2_d ...
!> \param v3rho3_d ...
!> \param v3rho2sigma_d ...
!> \param v3rhosigma2_d ...
!> \param v3sigma3_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_gga_vxc_fxc_kxc_device(func, rho_d, sigma_d, vrho_d, vsigma_d, &
                                           v2rho2_d, v2rhosigma_d, v2sigma2_d, &
                                           v3rho3_d, v3rho2sigma_d, v3rhosigma2_d, v3sigma3_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: vrho_d, vsigma_d, v2rho2_d, &
                                                            v2rhosigma_d, v2sigma2_d, v3rho3_d, &
                                                            v3rho2sigma_d, v3rhosigma2_d, &
                                                            v3sigma3_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER :: rho_p, sigma_p, v2rho2_p, v2rhosigma_p, v2sigma2_p, &
                                              v3rho2sigma_p, v3rho3_p, v3rhosigma2_p, v3sigma3_p, vrho_p, vsigma_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(vsigma_d, INT(vsigma_d%size), vsigma_p)
      CALL libxc_device_buffer_get_ptr(v2rho2_d, INT(v2rho2_d%size), v2rho2_p)
      CALL libxc_device_buffer_get_ptr(v2rhosigma_d, INT(v2rhosigma_d%size), v2rhosigma_p)
      CALL libxc_device_buffer_get_ptr(v2sigma2_d, INT(v2sigma2_d%size), v2sigma2_p)
      CALL libxc_device_buffer_get_ptr(v3rho3_d, INT(v3rho3_d%size), v3rho3_p)
      CALL libxc_device_buffer_get_ptr(v3rho2sigma_d, INT(v3rho2sigma_d%size), v3rho2sigma_p)
      CALL libxc_device_buffer_get_ptr(v3rhosigma2_d, INT(v3rhosigma2_d%size), v3rhosigma2_p)
      CALL libxc_device_buffer_get_ptr(v3sigma3_d, INT(v3sigma3_d%size), v3sigma3_p)

      CALL xc_f03_gga_vxc_fxc_kxc(func, np, rho_p, sigma_p, vrho_p, vsigma_p, &
                                  v2rho2_p, v2rhosigma_p, v2sigma2_p, v3rho3_p, &
                                  v3rho2sigma_p, v3rhosigma2_p, v3sigma3_p)
   END SUBROUTINE libxc_gga_vxc_fxc_kxc_device

! **************************************************************************************************
!> \brief GGA vxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param vrho_d ...
!> \param vsigma_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_gga_vxc_device(func, rho_d, sigma_d, vrho_d, vsigma_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: vrho_d, vsigma_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER               :: rho_p, sigma_p, vrho_p, vsigma_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(vsigma_d, INT(vsigma_d%size), vsigma_p)

      CALL xc_f03_gga_vxc(func, np, rho_p, sigma_p, vrho_p, vsigma_p)
   END SUBROUTINE libxc_gga_vxc_device

! **************************************************************************************************
!> \brief MGGA exc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param lapl_d ...
!> \param tau_d ...
!> \param exc_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_mgga_exc_device(func, rho_d, sigma_d, lapl_d, tau_d, exc_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d, lapl_d, tau_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: exc_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER               :: exc_p, lapl_p, rho_p, sigma_p, tau_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(lapl_d, INT(lapl_d%size), lapl_p)
      CALL libxc_device_buffer_get_ptr(tau_d, INT(tau_d%size), tau_p)
      CALL libxc_device_buffer_get_ptr(exc_d, INT(exc_d%size), exc_p)

      CALL xc_f03_mgga_exc(func, np, rho_p, sigma_p, lapl_p, tau_p, exc_p)
   END SUBROUTINE libxc_mgga_exc_device

! **************************************************************************************************
!> \brief MGGA exc and vxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param lapl_d ...
!> \param tau_d ...
!> \param exc_d ...
!> \param vrho_d ...
!> \param vsigma_d ...
!> \param vlapl_d ...
!> \param vtau_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_mgga_exc_vxc_device(func, rho_d, sigma_d, lapl_d, tau_d, exc_d, vrho_d, vsigma_d, vlapl_d, vtau_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d, lapl_d, tau_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: exc_d, vrho_d, vsigma_d, vlapl_d, vtau_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER               :: exc_p, lapl_p, rho_p, sigma_p, tau_p, &
                                                            vlapl_p, vrho_p, vsigma_p, vtau_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(lapl_d, INT(lapl_d%size), lapl_p)
      CALL libxc_device_buffer_get_ptr(tau_d, INT(tau_d%size), tau_p)
      CALL libxc_device_buffer_get_ptr(exc_d, INT(exc_d%size), exc_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(vsigma_d, INT(vsigma_d%size), vsigma_p)
      CALL libxc_device_buffer_get_ptr(vlapl_d, INT(vlapl_d%size), vlapl_p)
      CALL libxc_device_buffer_get_ptr(vtau_d, INT(vtau_d%size), vtau_p)

      CALL xc_f03_mgga_exc_vxc(func, np, rho_p, sigma_p, lapl_p, tau_p, &
                               exc_p, vrho_p, vsigma_p, vlapl_p, vtau_p)
   END SUBROUTINE libxc_mgga_exc_vxc_device

! **************************************************************************************************
!> \brief MGGA exc, vxc and fxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param lapl_d ...
!> \param tau_d ...
!> \param exc_d ...
!> \param vrho_d ...
!> \param vsigma_d ...
!> \param vlapl_d ...
!> \param vtau_d ...
!> \param v2rho2_d ...
!> \param v2rhosigma_d ...
!> \param v2rholapl_d ...
!> \param v2rhotau_d ...
!> \param v2sigma2_d ...
!> \param v2sigmalapl_d ...
!> \param v2sigmatau_d ...
!> \param v2lapl2_d ...
!> \param v2lapltau_d ...
!> \param v2tau2_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_mgga_exc_vxc_fxc_device(func, rho_d, sigma_d, lapl_d, tau_d, exc_d, &
                                            vrho_d, vsigma_d, vlapl_d, vtau_d, v2rho2_d, &
                                            v2rhosigma_d, v2rholapl_d, v2rhotau_d, v2sigma2_d, &
                                            v2sigmalapl_d, v2sigmatau_d, v2lapl2_d, v2lapltau_d, &
                                            v2tau2_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d, lapl_d, tau_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT) :: exc_d, vrho_d, vsigma_d, vlapl_d, vtau_d, &
                                         v2rho2_d, v2rhosigma_d, v2rholapl_d, v2rhotau_d, v2sigma2_d, v2sigmalapl_d, v2sigmatau_d, &
                                                       v2lapl2_d, v2lapltau_d, v2tau2_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER :: exc_p, lapl_p, rho_p, sigma_p, tau_p, v2lapl2_p, &
                                          v2lapltau_p, v2rho2_p, v2rholapl_p, v2rhosigma_p, v2rhotau_p, v2sigma2_p, v2sigmalapl_p, &
                                              v2sigmatau_p, v2tau2_p, vlapl_p, vrho_p, vsigma_p, vtau_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(lapl_d, INT(lapl_d%size), lapl_p)
      CALL libxc_device_buffer_get_ptr(tau_d, INT(tau_d%size), tau_p)
      CALL libxc_device_buffer_get_ptr(exc_d, INT(exc_d%size), exc_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(vsigma_d, INT(vsigma_d%size), vsigma_p)
      CALL libxc_device_buffer_get_ptr(vlapl_d, INT(vlapl_d%size), vlapl_p)
      CALL libxc_device_buffer_get_ptr(vtau_d, INT(vtau_d%size), vtau_p)
      CALL libxc_device_buffer_get_ptr(v2rho2_d, INT(v2rho2_d%size), v2rho2_p)
      CALL libxc_device_buffer_get_ptr(v2rhosigma_d, INT(v2rhosigma_d%size), v2rhosigma_p)
      CALL libxc_device_buffer_get_ptr(v2rholapl_d, INT(v2rholapl_d%size), v2rholapl_p)
      CALL libxc_device_buffer_get_ptr(v2rhotau_d, INT(v2rhotau_d%size), v2rhotau_p)
      CALL libxc_device_buffer_get_ptr(v2sigma2_d, INT(v2sigma2_d%size), v2sigma2_p)
      CALL libxc_device_buffer_get_ptr(v2sigmalapl_d, INT(v2sigmalapl_d%size), v2sigmalapl_p)
      CALL libxc_device_buffer_get_ptr(v2sigmatau_d, INT(v2sigmatau_d%size), v2sigmatau_p)
      CALL libxc_device_buffer_get_ptr(v2lapl2_d, INT(v2lapl2_d%size), v2lapl2_p)
      CALL libxc_device_buffer_get_ptr(v2lapltau_d, INT(v2lapltau_d%size), v2lapltau_p)
      CALL libxc_device_buffer_get_ptr(v2tau2_d, INT(v2tau2_d%size), v2tau2_p)

      CALL xc_f03_mgga_exc_vxc_fxc(func, np, rho_p, sigma_p, lapl_p, tau_p, &
                                   exc_p, vrho_p, vsigma_p, vlapl_p, vtau_p, &
                                   v2rho2_p, v2rhosigma_p, v2rholapl_p, &
                                   v2rhotau_p, v2sigma2_p, v2sigmalapl_p, &
                                   v2sigmatau_p, v2lapl2_p, v2lapltau_p, v2tau2_p)
   END SUBROUTINE libxc_mgga_exc_vxc_fxc_device

! **************************************************************************************************
!> \brief MGGA vxc and fxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param lapl_d ...
!> \param tau_d ...
!> \param vrho_d ...
!> \param vsigma_d ...
!> \param vlapl_d ...
!> \param vtau_d ...
!> \param v2rho2_d ...
!> \param v2rhosigma_d ...
!> \param v2rholapl_d ...
!> \param v2rhotau_d ...
!> \param v2sigma2_d ...
!> \param v2sigmalapl_d ...
!> \param v2sigmatau_d ...
!> \param v2lapl2_d ...
!> \param v2lapltau_d ...
!> \param v2tau2_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_mgga_vxc_fxc_device(func, rho_d, sigma_d, lapl_d, tau_d, vrho_d, &
                                        vsigma_d, vlapl_d, vtau_d, v2rho2_d, &
                                        v2rhosigma_d, v2rholapl_d, v2rhotau_d, v2sigma2_d, &
                                        v2sigmalapl_d, v2sigmatau_d, v2lapl2_d, v2lapltau_d, &
                                        v2tau2_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d, lapl_d, tau_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT) :: vrho_d, vsigma_d, vlapl_d, vtau_d, &
                                         v2rho2_d, v2rhosigma_d, v2rholapl_d, v2rhotau_d, v2sigma2_d, v2sigmalapl_d, v2sigmatau_d, &
                                                       v2lapl2_d, v2lapltau_d, v2tau2_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER :: lapl_p, rho_p, sigma_p, tau_p, v2lapl2_p, &
                                          v2lapltau_p, v2rho2_p, v2rholapl_p, v2rhosigma_p, v2rhotau_p, v2sigma2_p, v2sigmalapl_p, &
                                              v2sigmatau_p, v2tau2_p, vlapl_p, vrho_p, vsigma_p, vtau_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(lapl_d, INT(lapl_d%size), lapl_p)
      CALL libxc_device_buffer_get_ptr(tau_d, INT(tau_d%size), tau_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(vsigma_d, INT(vsigma_d%size), vsigma_p)
      CALL libxc_device_buffer_get_ptr(vlapl_d, INT(vlapl_d%size), vlapl_p)
      CALL libxc_device_buffer_get_ptr(vtau_d, INT(vtau_d%size), vtau_p)
      CALL libxc_device_buffer_get_ptr(v2rho2_d, INT(v2rho2_d%size), v2rho2_p)
      CALL libxc_device_buffer_get_ptr(v2rhosigma_d, INT(v2rhosigma_d%size), v2rhosigma_p)
      CALL libxc_device_buffer_get_ptr(v2rholapl_d, INT(v2rholapl_d%size), v2rholapl_p)
      CALL libxc_device_buffer_get_ptr(v2rhotau_d, INT(v2rhotau_d%size), v2rhotau_p)
      CALL libxc_device_buffer_get_ptr(v2sigma2_d, INT(v2sigma2_d%size), v2sigma2_p)
      CALL libxc_device_buffer_get_ptr(v2sigmalapl_d, INT(v2sigmalapl_d%size), v2sigmalapl_p)
      CALL libxc_device_buffer_get_ptr(v2sigmatau_d, INT(v2sigmatau_d%size), v2sigmatau_p)
      CALL libxc_device_buffer_get_ptr(v2lapl2_d, INT(v2lapl2_d%size), v2lapl2_p)
      CALL libxc_device_buffer_get_ptr(v2lapltau_d, INT(v2lapltau_d%size), v2lapltau_p)
      CALL libxc_device_buffer_get_ptr(v2tau2_d, INT(v2tau2_d%size), v2tau2_p)

      CALL xc_f03_mgga_vxc_fxc(func, np, rho_p, sigma_p, lapl_p, tau_p, &
                               vrho_p, vsigma_p, vlapl_p, vtau_p, &
                               v2rho2_p, v2rhosigma_p, v2rholapl_p, &
                               v2rhotau_p, v2sigma2_p, v2sigmalapl_p, &
                               v2sigmatau_p, v2lapl2_p, v2lapltau_p, v2tau2_p)
   END SUBROUTINE libxc_mgga_vxc_fxc_device

! **************************************************************************************************
!> \brief MGGA exc, vxc, fxc and kxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param lapl_d ...
!> \param tau_d ...
!> \param exc_d ...
!> \param vrho_d ...
!> \param vsigma_d ...
!> \param vlapl_d ...
!> \param vtau_d ...
!> \param v2rho2_d ...
!> \param v2rhosigma_d ...
!> \param v2rholapl_d ...
!> \param v2rhotau_d ...
!> \param v2sigma2_d ...
!> \param v2sigmalapl_d ...
!> \param v2sigmatau_d ...
!> \param v2lapl2_d ...
!> \param v2lapltau_d ...
!> \param v2tau2_d ...
!> \param v3rho3_d ...
!> \param v3rho2sigma_d ...
!> \param v3rho2lapl_d ...
!> \param v3rho2tau_d ...
!> \param v3rhosigma2_d ...
!> \param v3rhosigmalapl_d ...
!> \param v3rhosigmatau_d ...
!> \param v3rholapl2_d ...
!> \param v3rholapltau_d ...
!> \param v3rhotau2_d ...
!> \param v3sigma3_d ...
!> \param v3sigma2lapl_d ...
!> \param v3sigma2tau_d ...
!> \param v3sigmalapl2_d ...
!> \param v3sigmalapltau_d ...
!> \param v3sigmatau2_d ...
!> \param v3lapl3_d ...
!> \param v3lapl2tau_d ...
!> \param v3lapltau2_d ...
!> \param v3tau3_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_mgga_exc_vxc_fxc_kxc_device(func, rho_d, sigma_d, lapl_d, tau_d, exc_d, &
                                                vrho_d, vsigma_d, vlapl_d, vtau_d, v2rho2_d, &
                                                v2rhosigma_d, v2rholapl_d, v2rhotau_d, v2sigma2_d, &
                                                v2sigmalapl_d, v2sigmatau_d, v2lapl2_d, v2lapltau_d, &
                                                v2tau2_d, v3rho3_d, v3rho2sigma_d, v3rho2lapl_d, &
                                                v3rho2tau_d, v3rhosigma2_d, v3rhosigmalapl_d, &
                                                v3rhosigmatau_d, v3rholapl2_d, v3rholapltau_d, &
                                                v3rhotau2_d, v3sigma3_d, v3sigma2lapl_d, v3sigma2tau_d, &
                                                v3sigmalapl2_d, v3sigmalapltau_d, v3sigmatau2_d, &
                                                v3lapl3_d, v3lapl2tau_d, v3lapltau2_d, v3tau3_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d, lapl_d, tau_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT) :: exc_d, vrho_d, vsigma_d, vlapl_d, vtau_d, &
                                         v2rho2_d, v2rhosigma_d, v2rholapl_d, v2rhotau_d, v2sigma2_d, v2sigmalapl_d, v2sigmatau_d, &
                                             v2lapl2_d, v2lapltau_d, v2tau2_d, v3rho3_d, v3rho2sigma_d, v3rho2lapl_d, v3rho2tau_d, &
                                                   v3rhosigma2_d, v3rhosigmalapl_d, v3rhosigmatau_d, v3rholapl2_d, v3rholapltau_d, &
                                         v3rhotau2_d, v3sigma3_d, v3sigma2lapl_d, v3sigma2tau_d, v3sigmalapl2_d, v3sigmalapltau_d, &
                                                       v3sigmatau2_d, v3lapl3_d, v3lapl2tau_d, v3lapltau2_d, v3tau3_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER :: exc_p, lapl_p, rho_p, sigma_p, tau_p, v2lapl2_p, &
                                          v2lapltau_p, v2rho2_p, v2rholapl_p, v2rhosigma_p, v2rhotau_p, v2sigma2_p, v2sigmalapl_p, &
                                              v2sigmatau_p, v2tau2_p, v3lapl2tau_p, v3lapl3_p, v3lapltau2_p, v3rho2lapl_p, &
                                              v3rho2sigma_p, v3rho2tau_p, v3rho3_p, v3rholapl2_p, v3rholapltau_p, v3rhosigma2_p, &
                                              v3rhosigmalapl_p, v3rhosigmatau_p, v3rhotau2_p, v3sigma2lapl_p, v3sigma2tau_p, &
                                           v3sigma3_p, v3sigmalapl2_p, v3sigmalapltau_p, v3sigmatau2_p, v3tau3_p, vlapl_p, vrho_p, &
                                              vsigma_p, vtau_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(lapl_d, INT(lapl_d%size), lapl_p)
      CALL libxc_device_buffer_get_ptr(tau_d, INT(tau_d%size), tau_p)
      CALL libxc_device_buffer_get_ptr(exc_d, INT(exc_d%size), exc_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(vsigma_d, INT(vsigma_d%size), vsigma_p)
      CALL libxc_device_buffer_get_ptr(vlapl_d, INT(vlapl_d%size), vlapl_p)
      CALL libxc_device_buffer_get_ptr(vtau_d, INT(vtau_d%size), vtau_p)
      CALL libxc_device_buffer_get_ptr(v2rho2_d, INT(v2rho2_d%size), v2rho2_p)
      CALL libxc_device_buffer_get_ptr(v2rhosigma_d, INT(v2rhosigma_d%size), v2rhosigma_p)
      CALL libxc_device_buffer_get_ptr(v2rholapl_d, INT(v2rholapl_d%size), v2rholapl_p)
      CALL libxc_device_buffer_get_ptr(v2rhotau_d, INT(v2rhotau_d%size), v2rhotau_p)
      CALL libxc_device_buffer_get_ptr(v2sigma2_d, INT(v2sigma2_d%size), v2sigma2_p)
      CALL libxc_device_buffer_get_ptr(v2sigmalapl_d, INT(v2sigmalapl_d%size), v2sigmalapl_p)
      CALL libxc_device_buffer_get_ptr(v2sigmatau_d, INT(v2sigmatau_d%size), v2sigmatau_p)
      CALL libxc_device_buffer_get_ptr(v2lapl2_d, INT(v2lapl2_d%size), v2lapl2_p)
      CALL libxc_device_buffer_get_ptr(v2lapltau_d, INT(v2lapltau_d%size), v2lapltau_p)
      CALL libxc_device_buffer_get_ptr(v2tau2_d, INT(v2tau2_d%size), v2tau2_p)
      CALL libxc_device_buffer_get_ptr(v3rho3_d, INT(v3rho3_d%size), v3rho3_p)
      CALL libxc_device_buffer_get_ptr(v3rho2sigma_d, INT(v3rho2sigma_d%size), v3rho2sigma_p)
      CALL libxc_device_buffer_get_ptr(v3rho2lapl_d, INT(v3rho2lapl_d%size), v3rho2lapl_p)
      CALL libxc_device_buffer_get_ptr(v3rho2tau_d, INT(v3rho2tau_d%size), v3rho2tau_p)
      CALL libxc_device_buffer_get_ptr(v3rhosigma2_d, INT(v3rhosigma2_d%size), v3rhosigma2_p)
      CALL libxc_device_buffer_get_ptr(v3rhosigmalapl_d, INT(v3rhosigmalapl_d%size), v3rhosigmalapl_p)
      CALL libxc_device_buffer_get_ptr(v3rhosigmatau_d, INT(v3rhosigmatau_d%size), v3rhosigmatau_p)
      CALL libxc_device_buffer_get_ptr(v3rholapl2_d, INT(v3rholapl2_d%size), v3rholapl2_p)
      CALL libxc_device_buffer_get_ptr(v3rholapltau_d, INT(v3rholapltau_d%size), v3rholapltau_p)
      CALL libxc_device_buffer_get_ptr(v3rhotau2_d, INT(v3rhotau2_d%size), v3rhotau2_p)
      CALL libxc_device_buffer_get_ptr(v3sigma3_d, INT(v3sigma3_d%size), v3sigma3_p)
      CALL libxc_device_buffer_get_ptr(v3sigma2lapl_d, INT(v3sigma2lapl_d%size), v3sigma2lapl_p)
      CALL libxc_device_buffer_get_ptr(v3sigma2tau_d, INT(v3sigma2tau_d%size), v3sigma2tau_p)
      CALL libxc_device_buffer_get_ptr(v3sigmalapl2_d, INT(v3sigmalapl2_d%size), v3sigmalapl2_p)
      CALL libxc_device_buffer_get_ptr(v3sigmalapltau_d, INT(v3sigmalapltau_d%size), v3sigmalapltau_p)
      CALL libxc_device_buffer_get_ptr(v3sigmatau2_d, INT(v3sigmatau2_d%size), v3sigmatau2_p)
      CALL libxc_device_buffer_get_ptr(v3lapl3_d, INT(v3lapl3_d%size), v3lapl3_p)
      CALL libxc_device_buffer_get_ptr(v3lapl2tau_d, INT(v3lapl2tau_d%size), v3lapl2tau_p)
      CALL libxc_device_buffer_get_ptr(v3lapltau2_d, INT(v3lapltau2_d%size), v3lapltau2_p)
      CALL libxc_device_buffer_get_ptr(v3tau3_d, INT(v3tau3_d%size), v3tau3_p)

      CALL xc_f03_mgga_exc_vxc_fxc_kxc(func, np, rho_p, sigma_p, lapl_p, tau_p, &
                                       exc_p, vrho_p, vsigma_p, vlapl_p, vtau_p, &
                                       v2rho2_p, v2rhosigma_p, v2rholapl_p, &
                                       v2rhotau_p, v2sigma2_p, v2sigmalapl_p, &
                                       v2sigmatau_p, v2lapl2_p, v2lapltau_p, v2tau2_p, &
                                       v3rho3_p, v3rho2sigma_p, v3rho2lapl_p, v3rho2tau_p, &
                                       v3rhosigma2_p, v3rhosigmalapl_p, v3rhosigmatau_p, &
                                       v3rholapl2_p, v3rholapltau_p, v3rhotau2_p, &
                                       v3sigma3_p, v3sigma2lapl_p, v3sigma2tau_p, &
                                       v3sigmalapl2_p, v3sigmalapltau_p, v3sigmatau2_p, &
                                       v3lapl3_p, v3lapl2tau_p, v3lapltau2_p, v3tau3_p)
   END SUBROUTINE libxc_mgga_exc_vxc_fxc_kxc_device

! **************************************************************************************************
!> \brief MGGA vxc, fxc and kxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param lapl_d ...
!> \param tau_d ...
!> \param vrho_d ...
!> \param vsigma_d ...
!> \param vlapl_d ...
!> \param vtau_d ...
!> \param v2rho2_d ...
!> \param v2rhosigma_d ...
!> \param v2rholapl_d ...
!> \param v2rhotau_d ...
!> \param v2sigma2_d ...
!> \param v2sigmalapl_d ...
!> \param v2sigmatau_d ...
!> \param v2lapl2_d ...
!> \param v2lapltau_d ...
!> \param v2tau2_d ...
!> \param v3rho3_d ...
!> \param v3rho2sigma_d ...
!> \param v3rho2lapl_d ...
!> \param v3rho2tau_d ...
!> \param v3rhosigma2_d ...
!> \param v3rhosigmalapl_d ...
!> \param v3rhosigmatau_d ...
!> \param v3rholapl2_d ...
!> \param v3rholapltau_d ...
!> \param v3rhotau2_d ...
!> \param v3sigma3_d ...
!> \param v3sigma2lapl_d ...
!> \param v3sigma2tau_d ...
!> \param v3sigmalapl2_d ...
!> \param v3sigmalapltau_d ...
!> \param v3sigmatau2_d ...
!> \param v3lapl3_d ...
!> \param v3lapl2tau_d ...
!> \param v3lapltau2_d ...
!> \param v3tau3_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_mgga_vxc_fxc_kxc_device(func, rho_d, sigma_d, lapl_d, tau_d, vrho_d, &
                                            vsigma_d, vlapl_d, vtau_d, v2rho2_d, &
                                            v2rhosigma_d, v2rholapl_d, v2rhotau_d, v2sigma2_d, &
                                            v2sigmalapl_d, v2sigmatau_d, v2lapl2_d, v2lapltau_d, &
                                            v2tau2_d, v3rho3_d, v3rho2sigma_d, v3rho2lapl_d, &
                                            v3rho2tau_d, v3rhosigma2_d, v3rhosigmalapl_d, &
                                            v3rhosigmatau_d, v3rholapl2_d, v3rholapltau_d, &
                                            v3rhotau2_d, v3sigma3_d, v3sigma2lapl_d, v3sigma2tau_d, &
                                            v3sigmalapl2_d, v3sigmalapltau_d, v3sigmatau2_d, &
                                            v3lapl3_d, v3lapl2tau_d, v3lapltau2_d, v3tau3_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d, lapl_d, tau_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT) :: vrho_d, vsigma_d, vlapl_d, vtau_d, &
                                         v2rho2_d, v2rhosigma_d, v2rholapl_d, v2rhotau_d, v2sigma2_d, v2sigmalapl_d, v2sigmatau_d, &
                                             v2lapl2_d, v2lapltau_d, v2tau2_d, v3rho3_d, v3rho2sigma_d, v3rho2lapl_d, v3rho2tau_d, &
                                                   v3rhosigma2_d, v3rhosigmalapl_d, v3rhosigmatau_d, v3rholapl2_d, v3rholapltau_d, &
                                         v3rhotau2_d, v3sigma3_d, v3sigma2lapl_d, v3sigma2tau_d, v3sigmalapl2_d, v3sigmalapltau_d, &
                                                       v3sigmatau2_d, v3lapl3_d, v3lapl2tau_d, v3lapltau2_d, v3tau3_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER :: lapl_p, rho_p, sigma_p, tau_p, v2lapl2_p, &
                                          v2lapltau_p, v2rho2_p, v2rholapl_p, v2rhosigma_p, v2rhotau_p, v2sigma2_p, v2sigmalapl_p, &
                                              v2sigmatau_p, v2tau2_p, v3lapl2tau_p, v3lapl3_p, v3lapltau2_p, v3rho2lapl_p, &
                                              v3rho2sigma_p, v3rho2tau_p, v3rho3_p, v3rholapl2_p, v3rholapltau_p, v3rhosigma2_p, &
                                              v3rhosigmalapl_p, v3rhosigmatau_p, v3rhotau2_p, v3sigma2lapl_p, v3sigma2tau_p, &
                                           v3sigma3_p, v3sigmalapl2_p, v3sigmalapltau_p, v3sigmatau2_p, v3tau3_p, vlapl_p, vrho_p, &
                                              vsigma_p, vtau_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(lapl_d, INT(lapl_d%size), lapl_p)
      CALL libxc_device_buffer_get_ptr(tau_d, INT(tau_d%size), tau_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(vsigma_d, INT(vsigma_d%size), vsigma_p)
      CALL libxc_device_buffer_get_ptr(vlapl_d, INT(vlapl_d%size), vlapl_p)
      CALL libxc_device_buffer_get_ptr(vtau_d, INT(vtau_d%size), vtau_p)
      CALL libxc_device_buffer_get_ptr(v2rho2_d, INT(v2rho2_d%size), v2rho2_p)
      CALL libxc_device_buffer_get_ptr(v2rhosigma_d, INT(v2rhosigma_d%size), v2rhosigma_p)
      CALL libxc_device_buffer_get_ptr(v2rholapl_d, INT(v2rholapl_d%size), v2rholapl_p)
      CALL libxc_device_buffer_get_ptr(v2rhotau_d, INT(v2rhotau_d%size), v2rhotau_p)
      CALL libxc_device_buffer_get_ptr(v2sigma2_d, INT(v2sigma2_d%size), v2sigma2_p)
      CALL libxc_device_buffer_get_ptr(v2sigmalapl_d, INT(v2sigmalapl_d%size), v2sigmalapl_p)
      CALL libxc_device_buffer_get_ptr(v2sigmatau_d, INT(v2sigmatau_d%size), v2sigmatau_p)
      CALL libxc_device_buffer_get_ptr(v2lapl2_d, INT(v2lapl2_d%size), v2lapl2_p)
      CALL libxc_device_buffer_get_ptr(v2lapltau_d, INT(v2lapltau_d%size), v2lapltau_p)
      CALL libxc_device_buffer_get_ptr(v2tau2_d, INT(v2tau2_d%size), v2tau2_p)
      CALL libxc_device_buffer_get_ptr(v3rho3_d, INT(v3rho3_d%size), v3rho3_p)
      CALL libxc_device_buffer_get_ptr(v3rho2sigma_d, INT(v3rho2sigma_d%size), v3rho2sigma_p)
      CALL libxc_device_buffer_get_ptr(v3rho2lapl_d, INT(v3rho2lapl_d%size), v3rho2lapl_p)
      CALL libxc_device_buffer_get_ptr(v3rho2tau_d, INT(v3rho2tau_d%size), v3rho2tau_p)
      CALL libxc_device_buffer_get_ptr(v3rhosigma2_d, INT(v3rhosigma2_d%size), v3rhosigma2_p)
      CALL libxc_device_buffer_get_ptr(v3rhosigmalapl_d, INT(v3rhosigmalapl_d%size), v3rhosigmalapl_p)
      CALL libxc_device_buffer_get_ptr(v3rhosigmatau_d, INT(v3rhosigmatau_d%size), v3rhosigmatau_p)
      CALL libxc_device_buffer_get_ptr(v3rholapl2_d, INT(v3rholapl2_d%size), v3rholapl2_p)
      CALL libxc_device_buffer_get_ptr(v3rholapltau_d, INT(v3rholapltau_d%size), v3rholapltau_p)
      CALL libxc_device_buffer_get_ptr(v3rhotau2_d, INT(v3rhotau2_d%size), v3rhotau2_p)
      CALL libxc_device_buffer_get_ptr(v3sigma3_d, INT(v3sigma3_d%size), v3sigma3_p)
      CALL libxc_device_buffer_get_ptr(v3sigma2lapl_d, INT(v3sigma2lapl_d%size), v3sigma2lapl_p)
      CALL libxc_device_buffer_get_ptr(v3sigma2tau_d, INT(v3sigma2tau_d%size), v3sigma2tau_p)
      CALL libxc_device_buffer_get_ptr(v3sigmalapl2_d, INT(v3sigmalapl2_d%size), v3sigmalapl2_p)
      CALL libxc_device_buffer_get_ptr(v3sigmalapltau_d, INT(v3sigmalapltau_d%size), v3sigmalapltau_p)
      CALL libxc_device_buffer_get_ptr(v3sigmatau2_d, INT(v3sigmatau2_d%size), v3sigmatau2_p)
      CALL libxc_device_buffer_get_ptr(v3lapl3_d, INT(v3lapl3_d%size), v3lapl3_p)
      CALL libxc_device_buffer_get_ptr(v3lapl2tau_d, INT(v3lapl2tau_d%size), v3lapl2tau_p)
      CALL libxc_device_buffer_get_ptr(v3lapltau2_d, INT(v3lapltau2_d%size), v3lapltau2_p)
      CALL libxc_device_buffer_get_ptr(v3tau3_d, INT(v3tau3_d%size), v3tau3_p)

      CALL xc_f03_mgga_vxc_fxc_kxc(func, np, rho_p, sigma_p, lapl_p, tau_p, &
                                   vrho_p, vsigma_p, vlapl_p, vtau_p, &
                                   v2rho2_p, v2rhosigma_p, v2rholapl_p, &
                                   v2rhotau_p, v2sigma2_p, v2sigmalapl_p, &
                                   v2sigmatau_p, v2lapl2_p, v2lapltau_p, v2tau2_p, &
                                   v3rho3_p, v3rho2sigma_p, v3rho2lapl_p, v3rho2tau_p, &
                                   v3rhosigma2_p, v3rhosigmalapl_p, v3rhosigmatau_p, &
                                   v3rholapl2_p, v3rholapltau_p, v3rhotau2_p, &
                                   v3sigma3_p, v3sigma2lapl_p, v3sigma2tau_p, &
                                   v3sigmalapl2_p, v3sigmalapltau_p, v3sigmatau2_p, &
                                   v3lapl3_p, v3lapl2tau_p, v3lapltau2_p, v3tau3_p)
   END SUBROUTINE libxc_mgga_vxc_fxc_kxc_device

! **************************************************************************************************
!> \brief MGGA vxc on device buffers.
!> \param func ...
!> \param rho_d ...
!> \param sigma_d ...
!> \param lapl_d ...
!> \param tau_d ...
!> \param vrho_d ...
!> \param vsigma_d ...
!> \param vlapl_d ...
!> \param vtau_d ...
!> \param np ...
! **************************************************************************************************
   SUBROUTINE libxc_mgga_vxc_device(func, rho_d, sigma_d, lapl_d, tau_d, vrho_d, vsigma_d, vlapl_d, vtau_d, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      TYPE(libxc_device_buffer_type), INTENT(IN)         :: rho_d, sigma_d, lapl_d, tau_d
      TYPE(libxc_device_buffer_type), INTENT(INOUT)      :: vrho_d, vsigma_d, vlapl_d, vtau_d
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      REAL(KIND=dp), DIMENSION(:), POINTER               :: lapl_p, rho_p, sigma_p, tau_p, vlapl_p, &
                                                            vrho_p, vsigma_p, vtau_p

      CALL libxc_device_buffer_get_ptr(rho_d, INT(rho_d%size), rho_p)
      CALL libxc_device_buffer_get_ptr(sigma_d, INT(sigma_d%size), sigma_p)
      CALL libxc_device_buffer_get_ptr(lapl_d, INT(lapl_d%size), lapl_p)
      CALL libxc_device_buffer_get_ptr(tau_d, INT(tau_d%size), tau_p)
      CALL libxc_device_buffer_get_ptr(vrho_d, INT(vrho_d%size), vrho_p)
      CALL libxc_device_buffer_get_ptr(vsigma_d, INT(vsigma_d%size), vsigma_p)
      CALL libxc_device_buffer_get_ptr(vlapl_d, INT(vlapl_d%size), vlapl_p)
      CALL libxc_device_buffer_get_ptr(vtau_d, INT(vtau_d%size), vtau_p)

      CALL xc_f03_mgga_vxc(func, np, rho_p, sigma_p, lapl_p, tau_p, vrho_p, vsigma_p, vlapl_p, vtau_p)
   END SUBROUTINE libxc_mgga_vxc_device

#else
! **************************************************************************************************
!> \brief CPU-only build: a worker never carries device buffers. The device
!>        entry points libxc_worker_eval calls are compiled as no-ops/aborts so
!>        that no __OFFLOAD guards are needed above; the runtime backend
!>        selection prevents them from ever being reached.
! **************************************************************************************************
   SUBROUTINE libxc_worker_device_alloc(worker, family, pol, nl, np, ns, n, grow_only)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker
      INTEGER, INTENT(IN)                                :: family, nl, np, ns, n
      LOGICAL, INTENT(IN)                                :: pol
      LOGICAL, INTENT(IN), OPTIONAL                      :: grow_only

      MARK_USED(worker)
      MARK_USED(family)
      MARK_USED(pol)
      MARK_USED(nl)
      MARK_USED(np)
      MARK_USED(ns)
      MARK_USED(n)
      MARK_USED(grow_only)
   END SUBROUTINE libxc_worker_device_alloc

   SUBROUTINE libxc_worker_device_free(worker)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker

      MARK_USED(worker)
   END SUBROUTINE libxc_worker_device_free

! **************************************************************************************************
!> \brief CPU-only build: the whole dispatch is a single abort. The runtime
!>        backend selection can never route here.
! **************************************************************************************************
   SUBROUTINE libxc_worker_eval_device(worker, family, grad_deriv, no_exc, np)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker
      INTEGER, INTENT(IN)                                :: family, grad_deriv
      LOGICAL, INTENT(IN)                                :: no_exc
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np

      MARK_USED(worker)
      MARK_USED(family)
      MARK_USED(grad_deriv)
      MARK_USED(no_exc)
      MARK_USED(np)
      CPABORT("libxc_worker_eval_device: device backend not compiled in.")
   END SUBROUTINE libxc_worker_eval_device

! **************************************************************************************************
!> \brief CPU-only build: growing the device buffers is a no-op.
! **************************************************************************************************
   SUBROUTINE libxc_worker_resize_cuda(worker, family, nspin, npoints)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker
      INTEGER, INTENT(IN)                                :: family, nspin, npoints

      MARK_USED(worker)
      MARK_USED(family)
      MARK_USED(nspin)
      MARK_USED(npoints)
   END SUBROUTINE libxc_worker_resize_cuda
#endif
! **************************************************************************************************
!> \brief Natural family classification: whether a family carries the gradient
!>        (GGA) and the kinetic-energy density/Laplacian (MGGA) variables.
!> \param family LibXC functional family
!> \param is_gga the family has sigma/gradient components
!> \param is_mgga the family has tau/laplacian components
!> \param where optional routine name used in the abort message
! **************************************************************************************************
   SUBROUTINE libxc_family_flags(family, is_gga, is_mgga, where)
      INTEGER, INTENT(IN)                                :: family
      LOGICAL, INTENT(OUT)                               :: is_gga, is_mgga
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: where

      SELECT CASE (family)
      CASE (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
         is_gga = .FALSE.
         is_mgga = .FALSE.
      CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
         is_gga = .TRUE.
         is_mgga = .FALSE.
      CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
         is_gga = .TRUE.
         is_mgga = .TRUE.
      CASE default
         is_gga = .FALSE.
         is_mgga = .FALSE.
         IF (PRESENT(where)) THEN
            CPABORT(TRIM(where)//": this XC_FAMILY is currently not supported.")
         ELSE
            CPABORT("libxc_family_flags: unsupported XC_FAMILY.")
         END IF
      END SELECT
   END SUBROUTINE libxc_family_flags

! **************************************************************************************************
!> \brief Leading dimension (number of components per point) of the worker
!>        buffers for a spin polarization, cf. the LibXC manual.
!> \param nspin XC_UNPOLARIZED or XC_POLARIZED
!> \param np number of rho-like components (rho, vrho, lapl, tau, vlapl, vtau)
!> \param ns number of sigma-like components (sigma, vsigma, v2rho2, ...)
!> \param nl number of rho-lapl/tau-like components (v2rholapl, v3rho3, ...)
! **************************************************************************************************
   SUBROUTINE libxc_worker_dims(nspin, np, ns, nl)
      INTEGER, INTENT(IN)                                :: nspin
      INTEGER, INTENT(OUT)                               :: np, ns, nl

      IF (nspin == XC_POLARIZED) THEN
         np = 2
         ns = 3
         nl = 4
      ELSE
         np = 1
         ns = 1
         nl = 1
      END IF
   END SUBROUTINE libxc_worker_dims

! **************************************************************************************************
!> \brief Aborts when a freshly evaluated block contains a NaN, if checking is on.
!> \param worker the worker
!> \param npoints number of points in this block
!> \param grad_deriv derivative order requested
!> \param no_exc functional does not provide the energy density
!> \param family LibXC functional family
! **************************************************************************************************
   SUBROUTINE libxc_worker_check_nan(worker, npoints, grad_deriv, no_exc, family)
      TYPE(libxc_worker_type), INTENT(IN)                :: worker
      INTEGER, INTENT(IN)                                :: npoints, grad_deriv
      LOGICAL, INTENT(IN)                                :: no_exc
      INTEGER, INTENT(IN)                                :: family

      IF (.NOT. libxc_check_nan) RETURN
      IF (ANY(worker%exc(1:npoints) /= worker%exc(1:npoints)) .OR. &
          (grad_deriv >= 1 .AND. ANY(worker%vrho(:, 1:npoints) /= worker%vrho(:, 1:npoints)))) THEN
         WRITE (*, *) "NaN DETECTED: family=", family, " grad_deriv=", grad_deriv, &
            " no_exc=", no_exc, " nb=", npoints
         WRITE (*, *) "exc has NaN:", ANY(worker%exc(1:npoints) /= worker%exc(1:npoints))
         IF (grad_deriv >= 1) WRITE (*, *) "vrho has NaN:", ANY(worker%vrho(:, 1:npoints) /= worker%vrho(:, 1:npoints))
         CPABORT("NaN found right after device_to_host")
      END IF
   END SUBROUTINE libxc_worker_check_nan

! **************************************************************************************************
!> \brief Evaluates one block of points for a worker, on the device when the GPU
!>        backend is active and on the host otherwise.
!> \param worker the worker holding the functional and the staging buffers
!> \param family LibXC functional family
!> \param nspin XC_UNPOLARIZED or XC_POLARIZED
!> \param grad_deriv derivative order requested
!> \param no_exc functional does not provide the energy density
!> \param npoints number of points in this block; on the GPU path the whole grid
!>        is one block, so this is also the full grid size
!> \note The staged inputs (worker%rho/sigma/lapl/tau) must already be filled by
!>       the caller; the outputs (worker%exc/vrho/...) are filled on return.
! **************************************************************************************************
   SUBROUTINE libxc_worker_eval(worker, family, nspin, grad_deriv, no_exc, npoints)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker
      INTEGER, INTENT(IN)                                :: family, nspin, grad_deriv, npoints
      LOGICAL, INTENT(IN)                                :: no_exc

      INTEGER                                            :: nl, ns, np
      INTEGER(KIND=C_SIZE_T)                             :: npts
      LOGICAL                                            :: is_gga, is_mgga, pol

      CALL libxc_family_flags(family, is_gga, is_mgga)
      CALL libxc_worker_dims(nspin, np, ns, nl)
      pol = (nspin == XC_POLARIZED)
      npts = INT(npoints, KIND=C_SIZE_T)

      ! CUDA's current device is thread-local; make sure every thread that runs
      ! this block uses the rank's selected device. No-op without a device.
      CALL libxc_activate_device()

      IF (libxc_is_gpu()) THEN
         ! Copy the staged block to the worker's device buffers.
         CALL libxc_worker_h2d(worker, libxc_dev_rho, worker%rho(1, 1), INT(np*npts))
         IF (is_gga) THEN
            CALL libxc_worker_h2d(worker, libxc_dev_sigma, worker%sigma(1, 1), INT(ns*npts))
         END IF
         IF (is_mgga) THEN
            CALL libxc_worker_h2d(worker, libxc_dev_lapl, worker%lapl(1, 1), INT(np*npts))
            CALL libxc_worker_h2d(worker, libxc_dev_tau, worker%tau(1, 1), INT(np*npts))
         END IF

         ! One device call for the whole block; worker%exc is zeroed there if the
         ! functional does not provide it.
         CALL libxc_worker_eval_device(worker, family, grad_deriv, no_exc, npts)

         ! Issue all D2H copies first (async), then a single synchronize() covers
         ! all of them, then read every result.
         IF (.NOT. no_exc) THEN
            CALL libxc_worker_d2h_issue(worker, libxc_dev_exc, npoints)
         END IF
         IF (grad_deriv >= 1) THEN
            CALL libxc_worker_d2h_issue(worker, libxc_dev_vrho, np*npoints)
            IF (is_gga) THEN
               CALL libxc_worker_d2h_issue(worker, libxc_dev_vsigma, ns*npoints)
            END IF
            IF (is_mgga) THEN
               CALL libxc_worker_d2h_issue(worker, libxc_dev_vlapl, np*npoints)
               CALL libxc_worker_d2h_issue(worker, libxc_dev_vtau, np*npoints)
            END IF
         END IF
         IF (grad_deriv >= 2) THEN
            CALL libxc_worker_d2h_issue(worker, libxc_dev_v2rho2, ns*npoints)
            IF (is_gga) THEN
               CALL libxc_worker_d2h_issue(worker, libxc_dev_v2rhosigma, MERGE(6, 1, pol)*npoints)
               CALL libxc_worker_d2h_issue(worker, libxc_dev_v2sigma2, MERGE(6, 1, pol)*npoints)
            END IF
            IF (is_mgga) THEN
               CALL libxc_worker_d2h_issue(worker, libxc_dev_v2rholapl, nl*npoints)
               CALL libxc_worker_d2h_issue(worker, libxc_dev_v2rhotau, nl*npoints)
               CALL libxc_worker_d2h_issue(worker, libxc_dev_v2sigmalapl, MERGE(6, 1, pol)*npoints)
               CALL libxc_worker_d2h_issue(worker, libxc_dev_v2sigmatau, MERGE(6, 1, pol)*npoints)
               CALL libxc_worker_d2h_issue(worker, libxc_dev_v2lapl2, ns*npoints)
               CALL libxc_worker_d2h_issue(worker, libxc_dev_v2lapltau, nl*npoints)
               CALL libxc_worker_d2h_issue(worker, libxc_dev_v2tau2, ns*npoints)
            END IF
         END IF
         IF (grad_deriv >= 3) THEN
            CALL libxc_worker_d2h_issue(worker, libxc_dev_v3rho3, nl*npoints)
            IF (is_gga) THEN
               CALL libxc_worker_d2h_issue(worker, libxc_dev_v3rho2sigma, MERGE(9, 1, pol)*npoints)
               CALL libxc_worker_d2h_issue(worker, libxc_dev_v3rhosigma2, MERGE(12, 1, pol)*npoints)
               CALL libxc_worker_d2h_issue(worker, libxc_dev_v3sigma3, MERGE(10, 1, pol)*npoints)
            END IF
            #:for arr in v3_new
               IF (is_mgga) THEN
                  CALL libxc_worker_d2h_issue(worker, libxc_dev_${arr}$, MERGE(${mgga_ncomp[arr]}$, 1, pol)*npoints)
               END IF
            #:endfor
         END IF

         CALL libxc_sync_device()

         IF (.NOT. no_exc) THEN
            CALL libxc_worker_d2h_read(worker, libxc_dev_exc, worker%exc(1:npoints), npoints)
         END IF
         IF (grad_deriv >= 1) THEN
            CALL libxc_worker_d2h_read(worker, libxc_dev_vrho, worker%vrho(1, 1), np*npoints)
            IF (is_gga) THEN
               CALL libxc_worker_d2h_read(worker, libxc_dev_vsigma, worker%vsigma(1, 1), ns*npoints)
            END IF
            IF (is_mgga) THEN
               CALL libxc_worker_d2h_read(worker, libxc_dev_vlapl, worker%vlapl(1, 1), np*npoints)
               CALL libxc_worker_d2h_read(worker, libxc_dev_vtau, worker%vtau(1, 1), np*npoints)
            END IF
         END IF
         IF (grad_deriv >= 2) THEN
            CALL libxc_worker_d2h_read(worker, libxc_dev_v2rho2, worker%v2rho2(1, 1), ns*npoints)
            IF (is_gga) THEN
               CALL libxc_worker_d2h_read(worker, libxc_dev_v2rhosigma, worker%v2rhosigma(1, 1), MERGE(6, 1, pol)*npoints)
               CALL libxc_worker_d2h_read(worker, libxc_dev_v2sigma2, worker%v2sigma2(1, 1), MERGE(6, 1, pol)*npoints)
            END IF
            IF (is_mgga) THEN
               CALL libxc_worker_d2h_read(worker, libxc_dev_v2rholapl, worker%v2rholapl(1, 1), nl*npoints)
               CALL libxc_worker_d2h_read(worker, libxc_dev_v2rhotau, worker%v2rhotau(1, 1), nl*npoints)
               CALL libxc_worker_d2h_read(worker, libxc_dev_v2sigmalapl, worker%v2sigmalapl(1, 1), &
                                          MERGE(6, 1, pol)*npoints)
               CALL libxc_worker_d2h_read(worker, libxc_dev_v2sigmatau, worker%v2sigmatau(1, 1), &
                                          MERGE(6, 1, pol)*npoints)
               CALL libxc_worker_d2h_read(worker, libxc_dev_v2lapl2, worker%v2lapl2(1, 1), ns*npoints)
               CALL libxc_worker_d2h_read(worker, libxc_dev_v2lapltau, worker%v2lapltau(1, 1), nl*npoints)
               CALL libxc_worker_d2h_read(worker, libxc_dev_v2tau2, worker%v2tau2(1, 1), ns*npoints)
            END IF
         END IF
         IF (grad_deriv >= 3) THEN
            CALL libxc_worker_d2h_read(worker, libxc_dev_v3rho3, worker%v3rho3(1, 1), nl*npoints)
            IF (is_gga) THEN
               CALL libxc_worker_d2h_read(worker, libxc_dev_v3rho2sigma, worker%v3rho2sigma(1, 1), &
                                          MERGE(9, 1, pol)*npoints)
               CALL libxc_worker_d2h_read(worker, libxc_dev_v3rhosigma2, worker%v3rhosigma2(1, 1), &
                                          MERGE(12, 1, pol)*npoints)
               CALL libxc_worker_d2h_read(worker, libxc_dev_v3sigma3, worker%v3sigma3(1, 1), &
                                          MERGE(10, 1, pol)*npoints)
            END IF
            #:for arr in v3_new
               IF (is_mgga) THEN
                  CALL libxc_worker_d2h_read(worker, libxc_dev_${arr}$, worker%${arr}$ (1, 1), &
                                             MERGE(${mgga_ncomp[arr]}$, 1, pol)*npoints)
               END IF
            #:endfor
         END IF
      ELSE
         ! Host evaluation. The staging buffers already hold the block; these
         ! wrappers write the requested derivatives in place.
         IF (is_mgga) THEN
            SELECT CASE (grad_deriv)
            CASE (0)
               CALL libxc_mgga_exc_wrapper(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                           worker%lapl(1, 1), worker%tau(1, 1), worker%exc(1))
            CASE (1)
               IF (no_exc) THEN
                  CALL libxc_mgga_vxc_wrapper(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                              worker%lapl(1, 1), worker%tau(1, 1), worker%vrho(1, 1), &
                                              worker%vsigma(1, 1), worker%vlapl(1, 1), worker%vtau(1, 1))
                  worker%exc(1:npoints) = 0.0_dp
               ELSE
                  CALL libxc_mgga_exc_vxc_wrapper(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                                  worker%lapl(1, 1), worker%tau(1, 1), worker%exc(1), &
                                                  worker%vrho(1, 1), worker%vsigma(1, 1), &
                                                  worker%vlapl(1, 1), worker%vtau(1, 1))
               END IF
            CASE (2)
               IF (no_exc) THEN
                  CALL xc_f03_mgga_vxc_fxc(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                           worker%lapl(1, 1), worker%tau(1, 1), worker%vrho(1, 1), &
                                           worker%vsigma(1, 1), worker%vlapl(1, 1), worker%vtau(1, 1), &
                                           worker%v2rho2(1, 1), worker%v2rhosigma(1, 1), worker%v2rholapl(1, 1), &
                                           worker%v2rhotau(1, 1), worker%v2sigma2(1, 1), worker%v2sigmalapl(1, 1), &
                                           worker%v2sigmatau(1, 1), worker%v2lapl2(1, 1), worker%v2lapltau(1, 1), &
                                           worker%v2tau2(1, 1))
                  worker%exc(1:npoints) = 0.0_dp
               ELSE
                  CALL xc_f03_mgga(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                   worker%lapl(1, 1), worker%tau(1, 1), worker%exc(1), &
                                   worker%vrho(1, 1), worker%vsigma(1, 1), worker%vlapl(1, 1), worker%vtau(1, 1), &
                                   worker%v2rho2(1, 1), worker%v2rhosigma(1, 1), worker%v2rholapl(1, 1), &
                                   worker%v2rhotau(1, 1), worker%v2sigma2(1, 1), worker%v2sigmalapl(1, 1), &
                                   worker%v2sigmatau(1, 1), worker%v2lapl2(1, 1), worker%v2lapltau(1, 1), &
                                   worker%v2tau2(1, 1))
               END IF
            CASE (3)
               IF (no_exc) THEN
                  CALL xc_f03_mgga_vxc_fxc_kxc(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                               worker%lapl(1, 1), worker%tau(1, 1), worker%vrho(1, 1), &
                                               worker%vsigma(1, 1), worker%vlapl(1, 1), worker%vtau(1, 1), &
                                               worker%v2rho2(1, 1), worker%v2rhosigma(1, 1), worker%v2rholapl(1, 1), &
                                               worker%v2rhotau(1, 1), worker%v2sigma2(1, 1), worker%v2sigmalapl(1, 1), &
                                               worker%v2sigmatau(1, 1), worker%v2lapl2(1, 1), worker%v2lapltau(1, 1), &
                                               worker%v2tau2(1, 1), worker%v3rho3(1, 1), worker%v3rho2sigma(1, 1), &
                                               worker%v3rho2lapl(1, 1), worker%v3rho2tau(1, 1), &
                                               worker%v3rhosigma2(1, 1), worker%v3rhosigmalapl(1, 1), &
                                               worker%v3rhosigmatau(1, 1), worker%v3rholapl2(1, 1), &
                                               worker%v3rholapltau(1, 1), worker%v3rhotau2(1, 1), &
                                               worker%v3sigma3(1, 1), worker%v3sigma2lapl(1, 1), &
                                               worker%v3sigma2tau(1, 1), worker%v3sigmalapl2(1, 1), &
                                               worker%v3sigmalapltau(1, 1), worker%v3sigmatau2(1, 1), &
                                               worker%v3lapl3(1, 1), worker%v3lapl2tau(1, 1), &
                                               worker%v3lapltau2(1, 1), worker%v3tau3(1, 1))
                  worker%exc(1:npoints) = 0.0_dp
               ELSE
                  CALL xc_f03_mgga_exc_vxc_fxc_kxc(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                                   worker%lapl(1, 1), worker%tau(1, 1), worker%exc(1), &
                                                   worker%vrho(1, 1), worker%vsigma(1, 1), worker%vlapl(1, 1), &
                                                   worker%vtau(1, 1), worker%v2rho2(1, 1), worker%v2rhosigma(1, 1), &
                                                   worker%v2rholapl(1, 1), worker%v2rhotau(1, 1), &
                                                   worker%v2sigma2(1, 1), worker%v2sigmalapl(1, 1), &
                                                   worker%v2sigmatau(1, 1), worker%v2lapl2(1, 1), &
                                                   worker%v2lapltau(1, 1), worker%v2tau2(1, 1), &
                                                   worker%v3rho3(1, 1), worker%v3rho2sigma(1, 1), &
                                                   worker%v3rho2lapl(1, 1), worker%v3rho2tau(1, 1), &
                                                   worker%v3rhosigma2(1, 1), worker%v3rhosigmalapl(1, 1), &
                                                   worker%v3rhosigmatau(1, 1), worker%v3rholapl2(1, 1), &
                                                   worker%v3rholapltau(1, 1), worker%v3rhotau2(1, 1), &
                                                   worker%v3sigma3(1, 1), worker%v3sigma2lapl(1, 1), &
                                                   worker%v3sigma2tau(1, 1), worker%v3sigmalapl2(1, 1), &
                                                   worker%v3sigmalapltau(1, 1), worker%v3sigmatau2(1, 1), &
                                                   worker%v3lapl3(1, 1), worker%v3lapl2tau(1, 1), &
                                                   worker%v3lapltau2(1, 1), worker%v3tau3(1, 1))
               END IF
            END SELECT
         ELSE IF (is_gga) THEN
            SELECT CASE (grad_deriv)
            CASE (0)
               CALL libxc_gga_exc_wrapper(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), worker%exc(1))
            CASE (1)
               IF (no_exc) THEN
                  CALL xc_f03_gga_vxc(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                      worker%vrho(1, 1), worker%vsigma(1, 1))
                  worker%exc(1:npoints) = 0.0_dp
               ELSE
                  CALL libxc_gga_exc_vxc_wrapper(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                                 worker%exc(1), worker%vrho(1, 1), worker%vsigma(1, 1))
               END IF
            CASE (2)
               IF (no_exc) THEN
                  CALL xc_f03_gga_vxc_fxc(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                          worker%vrho(1, 1), worker%vsigma(1, 1), &
                                          worker%v2rho2(1, 1), worker%v2rhosigma(1, 1), worker%v2sigma2(1, 1))
                  worker%exc(1:npoints) = 0.0_dp
               ELSE
                  CALL xc_f03_gga_exc_vxc_fxc(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                              worker%exc(1), worker%vrho(1, 1), worker%vsigma(1, 1), &
                                              worker%v2rho2(1, 1), worker%v2rhosigma(1, 1), worker%v2sigma2(1, 1))
               END IF
            CASE (3)
               IF (no_exc) THEN
                  CALL xc_f03_gga_vxc_fxc_kxc(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                              worker%vrho(1, 1), worker%vsigma(1, 1), &
                                              worker%v2rho2(1, 1), worker%v2rhosigma(1, 1), worker%v2sigma2(1, 1), &
                                              worker%v3rho3(1, 1), worker%v3rho2sigma(1, 1), &
                                              worker%v3rhosigma2(1, 1), worker%v3sigma3(1, 1))
                  worker%exc(1:npoints) = 0.0_dp
               ELSE
                  CALL xc_f03_gga_exc_vxc_fxc_kxc(worker%func, npts, worker%rho(1, 1), worker%sigma(1, 1), &
                                                  worker%exc(1), worker%vrho(1, 1), worker%vsigma(1, 1), &
                                                  worker%v2rho2(1, 1), worker%v2rhosigma(1, 1), worker%v2sigma2(1, 1), &
                                                  worker%v3rho3(1, 1), worker%v3rho2sigma(1, 1), &
                                                  worker%v3rhosigma2(1, 1), worker%v3sigma3(1, 1))
               END IF
            END SELECT
         ELSE
            SELECT CASE (grad_deriv)
            CASE (0)
               CALL libxc_lda_exc_wrapper(worker%func, npts, worker%rho(1, 1), worker%exc(1))
            CASE (1)
               CALL libxc_lda_exc_vxc_wrapper(worker%func, npts, worker%rho(1, 1), worker%exc(1), worker%vrho(1, 1))
            CASE (2)
               CALL xc_f03_lda_exc_vxc_fxc(worker%func, npts, worker%rho(1, 1), &
                                           worker%exc(1), worker%vrho(1, 1), worker%v2rho2(1, 1))
            CASE (3)
               CALL xc_f03_lda(worker%func, npts, worker%rho(1, 1), &
                               worker%exc(1), worker%vrho(1, 1), worker%v2rho2(1, 1), worker%v3rho3(1, 1))
            END SELECT
         END IF
      END IF

      CALL libxc_worker_check_nan(worker, npoints, grad_deriv, no_exc, family)

   END SUBROUTINE libxc_worker_eval

! **************************************************************************************************
!> \brief GGA exc and vxc on host arrays.
!> \param func ...
!> \param np ...
!> \param rho ...
!> \param sigma ...
!> \param exc ...
!> \param vrho ...
!> \param vsigma ...
! **************************************************************************************************
   SUBROUTINE libxc_gga_exc_vxc_wrapper(func, np, rho, sigma, exc, vrho, vsigma)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, sigma
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT)           :: exc, vrho, vsigma

      CALL xc_f03_gga_exc_vxc(func, np, rho, sigma, exc, vrho, vsigma)
   END SUBROUTINE libxc_gga_exc_vxc_wrapper

! **************************************************************************************************
!> \brief GGA exc on host arrays.
!> \param func ...
!> \param np ...
!> \param rho ...
!> \param sigma ...
!> \param exc ...
! **************************************************************************************************
   SUBROUTINE libxc_gga_exc_wrapper(func, np, rho, sigma, exc)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, sigma
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT)           :: exc

      CALL xc_f03_gga_exc(func, np, rho, sigma, exc)
   END SUBROUTINE libxc_gga_exc_wrapper

! **************************************************************************************************
!> \brief LDA exc and vxc on host arrays.
!> \param func ...
!> \param np ...
!> \param rho ...
!> \param exc ...
!> \param vrho ...
! **************************************************************************************************
   SUBROUTINE libxc_lda_exc_vxc_wrapper(func, np, rho, exc, vrho)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT)           :: exc, vrho

      CALL xc_f03_lda_exc_vxc(func, np, rho, exc, vrho)
   END SUBROUTINE libxc_lda_exc_vxc_wrapper

! **************************************************************************************************
!> \brief LDA exc on host arrays.
!> \param func ...
!> \param np ...
!> \param rho ...
!> \param exc ...
! **************************************************************************************************
   SUBROUTINE libxc_lda_exc_wrapper(func, np, rho, exc)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT)           :: exc

      CALL xc_f03_lda_exc(func, np, rho, exc)
   END SUBROUTINE libxc_lda_exc_wrapper

! **************************************************************************************************
!> \brief MGGA exc and vxc on host arrays.
!> \param func ...
!> \param np ...
!> \param rho ...
!> \param sigma ...
!> \param lapl ...
!> \param tau ...
!> \param exc ...
!> \param vrho ...
!> \param vsigma ...
!> \param vlapl ...
!> \param vtau ...
! **************************************************************************************************
   SUBROUTINE libxc_mgga_exc_vxc_wrapper(func, np, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, sigma, lapl, tau
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT)           :: exc, vrho, vsigma, vlapl, vtau

      CALL xc_f03_mgga_exc_vxc(func, np, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau)
   END SUBROUTINE libxc_mgga_exc_vxc_wrapper

! **************************************************************************************************
!> \brief MGGA exc on host arrays.
!> \param func ...
!> \param np ...
!> \param rho ...
!> \param sigma ...
!> \param lapl ...
!> \param tau ...
!> \param exc ...
! **************************************************************************************************
   SUBROUTINE libxc_mgga_exc_wrapper(func, np, rho, sigma, lapl, tau, exc)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, sigma, lapl, tau
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT)           :: exc

      CALL xc_f03_mgga_exc(func, np, rho, sigma, lapl, tau, exc)
   END SUBROUTINE libxc_mgga_exc_wrapper

! **************************************************************************************************
!> \brief MGGA vxc on host arrays.
!> \param func ...
!> \param np ...
!> \param rho ...
!> \param sigma ...
!> \param lapl ...
!> \param tau ...
!> \param vrho ...
!> \param vsigma ...
!> \param vlapl ...
!> \param vtau ...
! **************************************************************************************************
   SUBROUTINE libxc_mgga_vxc_wrapper(func, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN)                    :: func
      INTEGER(KIND=C_SIZE_T), INTENT(IN)                 :: np
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, sigma, lapl, tau
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT)           :: vrho, vsigma, vlapl, vtau

      CALL xc_f03_mgga_vxc(func, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
   END SUBROUTINE libxc_mgga_vxc_wrapper

#endif

END MODULE xc_libxc_wrap
