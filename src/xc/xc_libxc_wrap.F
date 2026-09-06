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
MODULE xc_libxc_wrap
#if defined (__LIBXC)
#include <xc_version.h>
! check for LibXC version - updated for LibXC 7.x
#if (XC_MAJOR_VERSION < 7)
   This version of CP2K requires libxc version 7.x or above for full functionality.
   Some features may not be available with older versions.
#else
   USE kinds, ONLY: dp
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
   PUBLIC :: xc_f03_mgga_exc_vxc_fxc_kxc, xc_f03_mgga_vxc_fxc_kxc
   PUBLIC :: xc_f03_lda, &

   PUBLIC :: XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, &
             XC_FAMILY_HYB_LDA, XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA

   PUBLIC :: XC_UNPOLARIZED, XC_POLARIZED

    PUBLIC :: XC_EXCHANGE, XC_CORRELATION, XC_EXCHANGE_CORRELATION, XC_KINETIC
    PUBLIC :: XC_FLAGS_ON_HOST, XC_FLAGS_ON_DEVICE

    ! Create aliases for backward compatibility with LibXC 5.x interface
    PUBLIC :: xc_f03_lda, xc_f03_mgga

    ! Backend selection constants (similar to SPLA)
    INTEGER, PARAMETER, PUBLIC :: do_libxc_cpu = 1, &
                               do_libxc_gpu = 2

     ! Device buffer type for CUDA memory management
     TYPE, PUBLIC :: libxc_device_buffer_type
        TYPE(C_PTR) :: ptr = C_NULL_PTR
        INTEGER(KIND=C_SIZE_T) :: size = 0
     END TYPE libxc_device_buffer_type

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
             xc_libxc_wrap_init_cuda_safe, &
             libxc_set_backend, libxc_get_backend, &
             libxc_func_init_wrapper, libxc_func_end_wrapper, &
             libxc_set_library, &
             libxc_is_cuda, &
             libxc_device_buffer_alloc, libxc_device_buffer_free, &
             libxc_host_to_device, libxc_device_to_host, &
             libxc_assert_device_buffer, &
             libxc_buffer_create, libxc_buffer_free, &
             libxc_buffer_h2d, libxc_buffer_d2h, &
             libxc_gga_exc_vxc_wrapper, libxc_gga_exc_wrapper, &
             libxc_lda_exc_vxc_wrapper, libxc_lda_exc_wrapper, &
             libxc_mgga_exc_vxc_wrapper, libxc_mgga_exc_wrapper, &
             libxc_mgga_vxc_wrapper, &
              libxc_lda_exc_vxc_device, libxc_lda_exc_device, &
              libxc_gga_exc_device, libxc_gga_exc_vxc_device, &
              libxc_gga_vxc_device, &
              libxc_mgga_exc_device, libxc_mgga_exc_vxc_device, &
              libxc_mgga_vxc_device

   ! Backend selection for CPU vs GPU Libxc
   INTEGER, PARAMETER :: LIBXC_BACKEND_CPU = 1
   INTEGER, PARAMETER :: LIBXC_BACKEND_GPU = 2
    INTEGER :: libxc_current_backend = LIBXC_BACKEND_CPU
    INTEGER :: libxc_selected_library = do_libxc_cpu

   ! Alias definitions for backward compatibility
   PROCEDURE(xc_f03_lda_exc_vxc_fxc_kxc), POINTER :: xc_f03_lda => xc_f03_lda_exc_vxc_fxc_kxc
   PROCEDURE(xc_f03_mgga_exc_vxc_fxc), POINTER :: xc_f03_mgga => xc_f03_mgga_exc_vxc_fxc

CONTAINS

   ! Backend selection functions
! **************************************************************************************************
!> \brief ...
!> \param backend ...
! **************************************************************************************************
    SUBROUTINE libxc_set_backend(backend)
       INTEGER, INTENT(IN)                                :: backend

       SELECT CASE (backend)
       CASE (LIBXC_BACKEND_CPU)
          libxc_current_backend = LIBXC_BACKEND_CPU
       CASE (LIBXC_BACKEND_GPU)
#if defined (__LIBXC_CUDA)
          libxc_current_backend = LIBXC_BACKEND_GPU
#else
          CPABORT("GPU backend requested but CUDA Libxc not available")
#endif
       CASE DEFAULT
          CPABORT("Invalid Libxc backend specified")
       END SELECT
    END SUBROUTINE libxc_set_backend

! **************************************************************************************************
!> \brief Set the libxc library to use (CPU or GPU), similar to SPLA backend selection
!> \param library The library to use (do_libxc_cpu or do_libxc_gpu)
! **************************************************************************************************
    SUBROUTINE libxc_set_library(library)
       INTEGER, INTENT(IN)                                :: library

       SELECT CASE (library)
       CASE (do_libxc_cpu)
          libxc_selected_library = do_libxc_cpu
#if defined (__LIBXC_CUDA)
       CASE (do_libxc_gpu)
          libxc_selected_library = do_libxc_gpu
#else
       CASE (do_libxc_gpu)
          CPABORT("GPU library requested but CUDA Libxc not available")
#endif
       CASE DEFAULT
          CPABORT("Invalid Libxc library specified")
       END SELECT

       ! Automatically set the backend based on library selection
       IF (libxc_selected_library == do_libxc_gpu) THEN
          libxc_current_backend = LIBXC_BACKEND_GPU
       ELSE
          libxc_current_backend = LIBXC_BACKEND_CPU
       END IF
    END SUBROUTINE libxc_set_library

! **************************************************************************************************
!> \brief ...
!> \return ...
! **************************************************************************************************
   FUNCTION libxc_get_backend() RESULT(backend)
      INTEGER                                            :: backend

      backend = libxc_current_backend
   END FUNCTION libxc_get_backend

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
    SUBROUTINE xc_libxc_wrap_init_cuda_safe()
       CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_wrap_init_cuda_safe', &
           routineP = moduleN//':'//routineN

       INTEGER                                            :: handle

       CALL timeset(routineN, handle)

       ! Set default flags based on CUDA availability
#if defined (__LIBXC_CUDA)
       CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_DEVICE)
       libxc_selected_library = do_libxc_gpu
       libxc_current_backend = LIBXC_BACKEND_GPU
#else
       CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_HOST)
       libxc_selected_library = do_libxc_cpu
       libxc_current_backend = LIBXC_BACKEND_CPU
#endif

       CALL timestop(handle)
    END SUBROUTINE xc_libxc_wrap_init_cuda_safe

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

   ! Backend-aware wrapper functions for key Libxc operations
   ! These functions will dispatch to either CPU or CUDA implementations
   ! based on the current backend setting

! **************************************************************************************************
!> \brief ...
!> \param func ...
!> \param functional ...
!> \param nspin ...
! **************************************************************************************************
    SUBROUTINE libxc_func_init_wrapper(func, functional, nspin)
       TYPE(xc_f03_func_t), INTENT(OUT)                   :: func
       INTEGER, INTENT(IN)                                :: functional
       INTEGER, INTENT(IN), OPTIONAL                      :: nspin

       SELECT CASE (libxc_current_backend)
       CASE (LIBXC_BACKEND_CPU)
          CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_HOST)
          IF (PRESENT(nspin)) THEN
             CALL xc_f03_func_init(func, functional, nspin)
          ELSE
             CALL xc_f03_func_init(func, functional, XC_UNPOLARIZED)
          END IF
       CASE (LIBXC_BACKEND_GPU)
#if defined (__LIBXC_CUDA)
          CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_DEVICE)
          IF (PRESENT(nspin)) THEN
             CALL xc_f03_func_init(func, functional, nspin)
          ELSE
             CALL xc_f03_func_init(func, functional, XC_UNPOLARIZED)
          END IF
#else
          CPABORT("GPU backend requested but CUDA Libxc not available")
#endif
       END SELECT
    END SUBROUTINE libxc_func_init_wrapper

! **************************************************************************************************
!> \brief ...
!> \param func ...
! **************************************************************************************************
   SUBROUTINE libxc_func_end_wrapper(func)
      TYPE(xc_f03_func_t), INTENT(INOUT)                 :: func

      SELECT CASE (libxc_current_backend)
      CASE (LIBXC_BACKEND_CPU)
         CALL xc_f03_func_end(func)
      CASE (LIBXC_BACKEND_GPU)
#if defined (__LIBXC_CUDA)
         ! For now, use CPU version until proper CUDA symbols are available
         CALL xc_f03_func_end(func)
         ! Directly call CUDA-enabled functions (same symbols, GPU-accelerated internally)
#else
         CPABORT("GPU backend requested but CUDA Libxc not available")
#endif
      END SELECT
   END SUBROUTINE libxc_func_end_wrapper

#endif
#endif

    ! Device buffer management functions
    FUNCTION libxc_is_cuda() RESULT(is_cuda)
       LOGICAL :: is_cuda
       is_cuda = (libxc_current_backend == LIBXC_BACKEND_GPU)
    END FUNCTION libxc_is_cuda

     SUBROUTINE libxc_device_buffer_alloc(buffer, n)
        USE ISO_C_BINDING, ONLY: C_PTR, C_NULL_PTR, C_INT64_T
        USE offload_api, ONLY: offload_create_buffer, offload_buffer_type
        TYPE(libxc_device_buffer_type), INTENT(OUT) :: buffer
        INTEGER, INTENT(IN) :: n
        TYPE(offload_buffer_type) :: offload_buffer

        buffer%size = INT(n, KIND=C_INT64_T)
        CALL offload_create_buffer(n, offload_buffer)
        ! Store the device pointer from the offload buffer
        buffer%ptr = offload_buffer%c_ptr
     END SUBROUTINE libxc_device_buffer_alloc

     SUBROUTINE libxc_device_buffer_free(buffer)
        USE ISO_C_BINDING, ONLY: C_PTR, C_NULL_PTR
        USE offload_api, ONLY: offload_free_buffer, offload_buffer_type
        TYPE(libxc_device_buffer_type), INTENT(INOUT) :: buffer
        TYPE(offload_buffer_type) :: offload_buffer

        ! Create a temporary offload buffer with the device pointer
        offload_buffer%c_ptr = buffer%ptr
        offload_buffer%host_buffer => NULL()
        
        CALL offload_free_buffer(offload_buffer)
        buffer%ptr = C_NULL_PTR
        buffer%size = 0
     END SUBROUTINE libxc_device_buffer_free

      SUBROUTINE libxc_host_to_device(buffer, host, n)
        USE ISO_C_BINDING, ONLY: C_PTR, C_LOC
        USE offload_api, ONLY: offload_buffer_type, offload_buffer_h2d
        TYPE(libxc_device_buffer_type), INTENT(IN) :: buffer
        REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: host
        INTEGER, INTENT(IN) :: n
        TYPE(offload_buffer_type) :: offload_buffer
        REAL(KIND=dp), DIMENSION(n), TARGET :: temp_host

        ! Copy host data to a temporary array with proper bounds
        temp_host(:) = host(1:n)
        
        ! Create a temporary offload buffer with the device pointer
        offload_buffer%c_ptr = buffer%ptr
        offload_buffer%host_buffer => temp_host
        
        CALL offload_buffer_h2d(offload_buffer, n)
     END SUBROUTINE libxc_host_to_device

      SUBROUTINE libxc_device_to_host(host, buffer, n)
        USE ISO_C_BINDING, ONLY: C_PTR, C_LOC
        USE offload_api, ONLY: offload_buffer_type, offload_buffer_d2h
        REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: host
        TYPE(libxc_device_buffer_type), INTENT(IN) :: buffer
        INTEGER, INTENT(IN) :: n
        TYPE(offload_buffer_type) :: offload_buffer
        REAL(KIND=dp), DIMENSION(n), TARGET :: temp_host

        ! Create a temporary offload buffer with the device pointer
        offload_buffer%c_ptr = buffer%ptr
        offload_buffer%host_buffer => temp_host
        
        CALL offload_buffer_d2h(offload_buffer, n)
        
        ! Copy data from temporary array to output
        host(1:n) = temp_host(:)
     END SUBROUTINE libxc_device_to_host

     LOGICAL FUNCTION libxc_assert_device_buffer(buffer)
        USE ISO_C_BINDING, ONLY: C_PTR, C_NULL_PTR, C_ASSOCIATED
        TYPE(libxc_device_buffer_type), INTENT(IN) :: buffer

        libxc_assert_device_buffer = .TRUE.
        IF (.NOT. C_ASSOCIATED(buffer%ptr)) THEN
           CPABORT("Device buffer pointer is not associated")
           libxc_assert_device_buffer = .FALSE.
        END IF
        IF (buffer%size <= 0) THEN
           CPABORT("Device buffer size is invalid")
           libxc_assert_device_buffer = .FALSE.
        END IF
     END FUNCTION libxc_assert_device_buffer

    SUBROUTINE libxc_buffer_create(buffer, length)
       USE offload_api, ONLY: offload_buffer_type, offload_create_buffer
       TYPE(offload_buffer_type), INTENT(OUT) :: buffer
       INTEGER, INTENT(IN) :: length
       CALL offload_create_buffer(length, buffer)
    END SUBROUTINE libxc_buffer_create

   SUBROUTINE libxc_buffer_free(buffer)
      USE offload_api, ONLY: offload_buffer_type, offload_free_buffer
      TYPE(offload_buffer_type), INTENT(INOUT) :: buffer
      CALL offload_free_buffer(buffer)
   END SUBROUTINE libxc_buffer_free

   SUBROUTINE libxc_buffer_h2d(buffer, length)
      USE offload_api, ONLY: offload_buffer_type, offload_buffer_h2d
      TYPE(offload_buffer_type), INTENT(INOUT) :: buffer
      INTEGER, INTENT(IN) :: length
      CALL offload_buffer_h2d(buffer, length)
   END SUBROUTINE libxc_buffer_h2d

   SUBROUTINE libxc_buffer_d2h(buffer, length)
      USE offload_api, ONLY: offload_buffer_type, offload_buffer_d2h
      TYPE(offload_buffer_type), INTENT(INOUT) :: buffer
      INTEGER, INTENT(IN) :: length
      CALL offload_buffer_d2h(buffer, length)
   END SUBROUTINE libxc_buffer_d2h

   SUBROUTINE libxc_gga_exc_vxc_wrapper(func, np, rho, sigma, exc, vrho, vsigma)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN) :: func
      INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np
      REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: rho, sigma
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: exc, vrho, vsigma

      IF (libxc_current_backend == LIBXC_BACKEND_GPU) THEN
         CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_DEVICE)
         CALL xc_f03_gga_exc_vxc(func, np, rho, sigma, exc, vrho, vsigma)
         CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_HOST)
      ELSE
         CALL xc_f03_gga_exc_vxc(func, np, rho, sigma, exc, vrho, vsigma)
      END IF
   END SUBROUTINE libxc_gga_exc_vxc_wrapper

   SUBROUTINE libxc_gga_exc_wrapper(func, np, rho, sigma, exc)
      USE ISO_C_BINDING, ONLY: C_SIZE_T
      TYPE(xc_f03_func_t), INTENT(IN) :: func
      INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np
      REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: rho, sigma
      REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: exc

      IF (libxc_current_backend == LIBXC_BACKEND_GPU) THEN
         CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_DEVICE)
         CALL xc_f03_gga_exc(func, np, rho, sigma, exc)
         CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_HOST)
      ELSE
         CALL xc_f03_gga_exc(func, np, rho, sigma, exc)
      END IF
   END SUBROUTINE libxc_gga_exc_wrapper

    SUBROUTINE libxc_lda_exc_vxc_wrapper(func, np, rho, exc, vrho)
       USE ISO_C_BINDING, ONLY: C_SIZE_T
       TYPE(xc_f03_func_t), INTENT(IN) :: func
       INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np
       REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: rho
       REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: exc, vrho

       IF (libxc_current_backend == LIBXC_BACKEND_GPU) THEN
          CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_DEVICE)
          CALL xc_f03_lda_exc_vxc(func, np, rho, exc, vrho)
          CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_HOST)
       ELSE
          CALL xc_f03_lda_exc_vxc(func, np, rho, exc, vrho)
       END IF
    END SUBROUTINE libxc_lda_exc_vxc_wrapper

     SUBROUTINE libxc_lda_exc_wrapper(func, np, rho, exc)
        USE ISO_C_BINDING, ONLY: C_SIZE_T
        TYPE(xc_f03_func_t), INTENT(IN) :: func
        INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np
        REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: rho
        REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: exc

        IF (libxc_current_backend == LIBXC_BACKEND_GPU) THEN
           CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_DEVICE)
           CALL xc_f03_lda_exc(func, np, rho, exc)
           CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_HOST)
        ELSE
           CALL xc_f03_lda_exc(func, np, rho, exc)
        END IF
     END SUBROUTINE libxc_lda_exc_wrapper

     ! Device-aware LDA wrapper functions
     ! These functions accept libxc_device_buffer_type and handle the H2D/D2H transfers
     SUBROUTINE libxc_lda_exc_vxc_device(func, rho_d, exc_d, vrho_d, np)
        USE ISO_C_BINDING, ONLY: C_SIZE_T
        TYPE(xc_f03_func_t), INTENT(IN) :: func
        TYPE(libxc_device_buffer_type), INTENT(IN) :: rho_d
        TYPE(libxc_device_buffer_type), INTENT(OUT) :: exc_d, vrho_d
        INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np

        REAL(KIND=dp), DIMENSION(INT(np)) :: rho_host, exc_host, vrho_host

        ! Copy input from device to host
        CALL libxc_device_to_host(rho_host, rho_d, INT(np))
        
        ! Call regular Libxc function on host
        CALL xc_f03_lda_exc_vxc(func, np, rho_host, exc_host, vrho_host)
        
        ! Copy results back to device
        CALL libxc_host_to_device(exc_d, exc_host, INT(np))
        CALL libxc_host_to_device(vrho_d, vrho_host, INT(np))
     END SUBROUTINE libxc_lda_exc_vxc_device

     SUBROUTINE libxc_lda_exc_device(func, rho_d, exc_d, np)
        USE ISO_C_BINDING, ONLY: C_SIZE_T
        TYPE(xc_f03_func_t), INTENT(IN) :: func
        TYPE(libxc_device_buffer_type), INTENT(IN) :: rho_d
        TYPE(libxc_device_buffer_type), INTENT(OUT) :: exc_d
        INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np

        REAL(KIND=dp), DIMENSION(INT(np)) :: rho_host, exc_host

        ! Copy input from device to host
        CALL libxc_device_to_host(rho_host, rho_d, INT(np))
        
        ! Call regular Libxc function on host
        CALL xc_f03_lda_exc(func, np, rho_host, exc_host)
        
         ! Copy results back to device
         CALL libxc_host_to_device(exc_d, exc_host, INT(np))
      END SUBROUTINE libxc_lda_exc_device

     ! Device-aware GGA wrapper functions
     SUBROUTINE libxc_gga_exc_device(func, rho_d, sigma_d, exc_d, np)
        USE ISO_C_BINDING, ONLY: C_SIZE_T
        TYPE(xc_f03_func_t), INTENT(IN) :: func
        TYPE(libxc_device_buffer_type), INTENT(IN) :: rho_d, sigma_d
        TYPE(libxc_device_buffer_type), INTENT(OUT) :: exc_d
        INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np

        REAL(KIND=dp), DIMENSION(INT(np)) :: rho_host, sigma_host, exc_host

        ! Copy input from device to host
        CALL libxc_device_to_host(rho_host, rho_d, INT(np))
        CALL libxc_device_to_host(sigma_host, sigma_d, INT(np))
        
        ! Call regular Libxc function on host
        CALL xc_f03_gga_exc(func, np, rho_host, sigma_host, exc_host)
        
        ! Copy results back to device
        CALL libxc_host_to_device(exc_d, exc_host, INT(np))
     END SUBROUTINE libxc_gga_exc_device

      SUBROUTINE libxc_gga_exc_vxc_device(func, rho_d, sigma_d, exc_d, vrho_d, vsigma_d, np)
         USE ISO_C_BINDING, ONLY: C_SIZE_T
         TYPE(xc_f03_func_t), INTENT(IN) :: func
         TYPE(libxc_device_buffer_type), INTENT(IN) :: rho_d, sigma_d
         TYPE(libxc_device_buffer_type), INTENT(OUT) :: exc_d, vrho_d, vsigma_d
         INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np

         REAL(KIND=dp), DIMENSION(INT(np)) :: rho_host, sigma_host, exc_host, vrho_host, vsigma_host

         ! Copy input from device to host
         CALL libxc_device_to_host(rho_host, rho_d, INT(np))
         CALL libxc_device_to_host(sigma_host, sigma_d, INT(np))
         
         ! Call regular Libxc function on host
         CALL xc_f03_gga_exc_vxc(func, np, rho_host, sigma_host, exc_host, vrho_host, vsigma_host)
         
         ! Copy results back to device
         CALL libxc_host_to_device(exc_d, exc_host, INT(np))
         CALL libxc_host_to_device(vrho_d, vrho_host, INT(np))
         CALL libxc_host_to_device(vsigma_d, vsigma_host, INT(np))
      END SUBROUTINE libxc_gga_exc_vxc_device

      ! Device-aware GGA vxc wrapper function (without exc)
      SUBROUTINE libxc_gga_vxc_device(func, rho_d, sigma_d, vrho_d, vsigma_d, np)
         USE ISO_C_BINDING, ONLY: C_SIZE_T
         TYPE(xc_f03_func_t), INTENT(IN) :: func
         TYPE(libxc_device_buffer_type), INTENT(IN) :: rho_d, sigma_d
         TYPE(libxc_device_buffer_type), INTENT(OUT) :: vrho_d, vsigma_d
         INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np

         REAL(KIND=dp), DIMENSION(INT(np)) :: rho_host, sigma_host, vrho_host, vsigma_host

         ! Copy input from device to host
         CALL libxc_device_to_host(rho_host, rho_d, INT(np))
         CALL libxc_device_to_host(sigma_host, sigma_d, INT(np))
         
         ! Call regular Libxc function on host
         CALL xc_f03_gga_vxc(func, np, rho_host, sigma_host, vrho_host, vsigma_host)
         
         ! Copy results back to device
         CALL libxc_host_to_device(vrho_d, vrho_host, INT(np))
         CALL libxc_host_to_device(vsigma_d, vsigma_host, INT(np))
      END SUBROUTINE libxc_gga_vxc_device

     ! Device-aware MGGA wrapper functions
     SUBROUTINE libxc_mgga_exc_device(func, rho_d, sigma_d, lapl_d, tau_d, exc_d, np)
        USE ISO_C_BINDING, ONLY: C_SIZE_T
        TYPE(xc_f03_func_t), INTENT(IN) :: func
        TYPE(libxc_device_buffer_type), INTENT(IN) :: rho_d, sigma_d, lapl_d, tau_d
        TYPE(libxc_device_buffer_type), INTENT(OUT) :: exc_d
        INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np

        REAL(KIND=dp), DIMENSION(INT(np)) :: rho_host, sigma_host, lapl_host, tau_host, exc_host

        ! Copy input from device to host
        CALL libxc_device_to_host(rho_host, rho_d, INT(np))
        CALL libxc_device_to_host(sigma_host, sigma_d, INT(np))
        CALL libxc_device_to_host(lapl_host, lapl_d, INT(np))
        CALL libxc_device_to_host(tau_host, tau_d, INT(np))
        
        ! Call regular Libxc function on host
        CALL xc_f03_mgga_exc(func, np, rho_host, sigma_host, lapl_host, tau_host, exc_host)
        
        ! Copy results back to device
        CALL libxc_host_to_device(exc_d, exc_host, INT(np))
     END SUBROUTINE libxc_mgga_exc_device

      SUBROUTINE libxc_mgga_exc_vxc_device(func, rho_d, sigma_d, lapl_d, tau_d, exc_d, vrho_d, vsigma_d, vlapl_d, vtau_d, np)
         USE ISO_C_BINDING, ONLY: C_SIZE_T
         TYPE(xc_f03_func_t), INTENT(IN) :: func
         TYPE(libxc_device_buffer_type), INTENT(IN) :: rho_d, sigma_d, lapl_d, tau_d
         TYPE(libxc_device_buffer_type), INTENT(OUT) :: exc_d, vrho_d, vsigma_d, vlapl_d, vtau_d
         INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np

         REAL(KIND=dp), DIMENSION(INT(np)) :: rho_host, sigma_host, lapl_host, tau_host, exc_host
         REAL(KIND=dp), DIMENSION(INT(np)) :: vrho_host, vsigma_host, vlapl_host, vtau_host

         ! Copy input from device to host
         CALL libxc_device_to_host(rho_host, rho_d, INT(np))
         CALL libxc_device_to_host(sigma_host, sigma_d, INT(np))
         CALL libxc_device_to_host(lapl_host, lapl_d, INT(np))
         CALL libxc_device_to_host(tau_host, tau_d, INT(np))
         
         ! Call regular Libxc function on host
         CALL xc_f03_mgga_exc_vxc(func, np, rho_host, sigma_host, lapl_host, tau_host, &
                                    exc_host, vrho_host, vsigma_host, vlapl_host, vtau_host)
         
         ! Copy results back to device
         CALL libxc_host_to_device(exc_d, exc_host, INT(np))
         CALL libxc_host_to_device(vrho_d, vrho_host, INT(np))
         CALL libxc_host_to_device(vsigma_d, vsigma_host, INT(np))
         CALL libxc_host_to_device(vlapl_d, vlapl_host, INT(np))
         CALL libxc_host_to_device(vtau_d, vtau_host, INT(np))
      END SUBROUTINE libxc_mgga_exc_vxc_device

      ! Device-aware MGGA vxc wrapper function (without exc)
      SUBROUTINE libxc_mgga_vxc_device(func, rho_d, sigma_d, lapl_d, tau_d, vrho_d, vsigma_d, vlapl_d, vtau_d, np)
         USE ISO_C_BINDING, ONLY: C_SIZE_T
         TYPE(xc_f03_func_t), INTENT(IN) :: func
         TYPE(libxc_device_buffer_type), INTENT(IN) :: rho_d, sigma_d, lapl_d, tau_d
         TYPE(libxc_device_buffer_type), INTENT(OUT) :: vrho_d, vsigma_d, vlapl_d, vtau_d
         INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np

         REAL(KIND=dp), DIMENSION(INT(np)) :: rho_host, sigma_host, lapl_host, tau_host
         REAL(KIND=dp), DIMENSION(INT(np)) :: vrho_host, vsigma_host, vlapl_host, vtau_host

         ! Copy input from device to host
         CALL libxc_device_to_host(rho_host, rho_d, INT(np))
         CALL libxc_device_to_host(sigma_host, sigma_d, INT(np))
         CALL libxc_device_to_host(lapl_host, lapl_d, INT(np))
         CALL libxc_device_to_host(tau_host, tau_d, INT(np))
         
         ! Call regular Libxc function on host
         CALL xc_f03_mgga_vxc(func, np, rho_host, sigma_host, lapl_host, tau_host, &
                               vrho_host, vsigma_host, vlapl_host, vtau_host)
         
         ! Copy results back to device
         CALL libxc_host_to_device(vrho_d, vrho_host, INT(np))
         CALL libxc_host_to_device(vsigma_d, vsigma_host, INT(np))
         CALL libxc_host_to_device(vlapl_d, vlapl_host, INT(np))
         CALL libxc_host_to_device(vtau_d, vtau_host, INT(np))
      END SUBROUTINE libxc_mgga_vxc_device

    SUBROUTINE libxc_mgga_exc_vxc_wrapper(func, np, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau)
       USE ISO_C_BINDING, ONLY: C_SIZE_T
       TYPE(xc_f03_func_t), INTENT(IN) :: func
       INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np
       REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: rho, sigma, lapl, tau
       REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: exc, vrho, vsigma, vlapl, vtau

       IF (libxc_current_backend == LIBXC_BACKEND_GPU) THEN
          WRITE(*,*) "MGGA exc_vxc: Using GPU backend, np =", np
          CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_DEVICE)
          CALL xc_f03_mgga_exc_vxc(func, np, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau)
          CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_HOST)
       ELSE
          WRITE(*,*) "MGGA exc_vxc: Using CPU backend, np =", np
          CALL xc_f03_mgga_exc_vxc(func, np, rho, sigma, lapl, tau, exc, vrho, vsigma, vlapl, vtau)
       END IF
    END SUBROUTINE libxc_mgga_exc_vxc_wrapper

    SUBROUTINE libxc_mgga_exc_wrapper(func, np, rho, sigma, lapl, tau, exc)
       USE ISO_C_BINDING, ONLY: C_SIZE_T
       TYPE(xc_f03_func_t), INTENT(IN) :: func
       INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np
       REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: rho, sigma, lapl, tau
       REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: exc

       IF (libxc_current_backend == LIBXC_BACKEND_GPU) THEN
          WRITE(*,*) "MGGA exc: Using GPU backend, np =", np
          CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_DEVICE)
          CALL xc_f03_mgga_exc(func, np, rho, sigma, lapl, tau, exc)
          CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_HOST)
       ELSE
          WRITE(*,*) "MGGA exc: Using CPU backend, np =", np
          CALL xc_f03_mgga_exc(func, np, rho, sigma, lapl, tau, exc)
       END IF
    END SUBROUTINE libxc_mgga_exc_wrapper

    SUBROUTINE libxc_mgga_vxc_wrapper(func, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
       USE ISO_C_BINDING, ONLY: C_SIZE_T
       TYPE(xc_f03_func_t), INTENT(IN) :: func
       INTEGER(KIND=C_SIZE_T), INTENT(IN) :: np
       REAL(KIND=dp), DIMENSION(*), INTENT(IN) :: rho, sigma, lapl, tau
       REAL(KIND=dp), DIMENSION(*), INTENT(OUT) :: vrho, vsigma, vlapl, vtau

       IF (libxc_current_backend == LIBXC_BACKEND_GPU) THEN
          WRITE(*,*) "MGGA vxc: Using GPU backend, np =", np
          CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_DEVICE)
          CALL xc_f03_mgga_vxc(func, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
          CALL xc_f03_func_info_set_default_flags(XC_FLAGS_ON_HOST)
       ELSE
          WRITE(*,*) "MGGA vxc: Using CPU backend, np =", np
          CALL xc_f03_mgga_vxc(func, np, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau)
       END IF
    END SUBROUTINE libxc_mgga_vxc_wrapper

END MODULE xc_libxc_wrap
