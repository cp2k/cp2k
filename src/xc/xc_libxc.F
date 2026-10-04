!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief calculates a functional from libxc and its derivatives
!> \note
!>      LibXC:
!>      (Marques, Oliveira, Burnus, CPC 183, 2272 (2012)).
!>
!>      WARNING: In the subroutine libxc_spin_polarized_calc, it could be that the
!>      ordering for the 1st index of v2lapltau, v2rholapl, v2rhotau,
!>      v2sigmalapl and v2sigmatau is not correct. For the moment it does not
!>      matter since the calculation of the 2nd derivatives for meta-GGA
!>      functionals is not implemented in CP2K.
!>
!> \par History
!>      01.2013 created [F. Tran]
!>      07.2014 updates to versions 2.1 [JGH]
!>      08.2015 refactoring [A. Gloess (agloess)]
!>      01.2018 refactoring [A. Gloess (agloess)]
!>      10.2018/04.2019 added hyb_mgga [S. Simko, included by F. Stein]
!> \author F. Tran
! **************************************************************************************************
#:include 'xc_gamma.fypp'

MODULE xc_libxc
   USE bibliography, ONLY: Lehtola2018, &
                           Marques2012, &
                           cite_reference
   USE input_section_types, ONLY: section_add_keyword, &
                                  section_add_subsection, &
                                  section_create, &
                                  section_release, &
                                  section_type, &
                                  section_vals_type, &
                                  section_vals_val_get
   USE kinds, ONLY: default_string_length, &
                    dp
   USE xc_derivative_set_types, ONLY: xc_derivative_set_type, &
                                      xc_dset_get_derivative
   USE xc_derivative_types, ONLY: xc_derivative_get, &
                                  xc_derivative_type
   USE xc_rho_cflags_types, ONLY: xc_rho_cflags_type
   USE xc_rho_set_types, ONLY: xc_rho_set_get, &
                               xc_rho_set_type
#include "../base/base_uses.f90"
#if defined (__LIBXC)
   USE input_keyword_types, ONLY: keyword_create, &
                                  keyword_release, &
                                  keyword_type
   USE iso_c_binding, ONLY: C_INT, C_DOUBLE
   USE xc_derivative_desc, ONLY: &
      deriv_rho, deriv_rhoa, deriv_rhob, &
      deriv_norm_drhoa, deriv_norm_drhob, deriv_norm_drho, deriv_tau_a, deriv_tau_b, deriv_tau, &
      deriv_laplace_rho, deriv_laplace_rhoa, deriv_laplace_rhob, &
      deriv_gamma, deriv_gamma_aa, deriv_gamma_ab, deriv_gamma_bb
   USE xc_libxc_wrap, ONLY: xc_f03_func_t, &
                            xc_f03_func_init, &
                            xc_f03_func_end, &
                            xc_f03_func_info_t, &
                            xc_f03_functional_get_name, &
                            xc_f03_func_get_info, &
                            xc_f03_func_info_get_family, &
                            xc_f03_func_info_get_kind, &
                            xc_f03_func_info_get_n_ext_params, &
                            xc_f03_func_info_get_name, &
                            xc_f03_available_functional_numbers, &
                            xc_f03_available_functional_names, &
                            xc_f03_maximum_name_length, &
                            xc_f03_number_of_functionals, &
                            xc_f03_func_info_get_ext_params_name, &
                            xc_f03_func_info_get_ext_params_description, &
                            xc_f03_func_info_get_ext_params_default_value, &
                            xc_f03_gga_exc, &
                            xc_f03_gga_exc_vxc, &
                            xc_f03_gga_exc_vxc_fxc, &
                            xc_f03_gga_exc_vxc_fxc_kxc, &
                            xc_f03_gga_fxc, &
                            xc_f03_gga_vxc, &
                            xc_f03_gga_vxc_fxc, &
                            xc_f03_gga_vxc_fxc_kxc, &
                            xc_f03_lda, &
                            xc_f03_lda_exc, &
                            xc_f03_lda_exc_vxc, &
                            xc_f03_lda_exc_vxc_fxc, &
                            xc_f03_lda_fxc, &
                            xc_f03_lda_kxc, &
                            xc_f03_lda_vxc, &
                            xc_f03_mgga, &
                            xc_f03_mgga_exc, &
                            xc_f03_mgga_exc_vxc, &
                            xc_f03_mgga_exc_vxc_fxc, &
                            xc_f03_mgga_exc_vxc_fxc_kxc, &
                            xc_f03_mgga_fxc, &
                            xc_f03_mgga_vxc, &
                            xc_f03_mgga_vxc_fxc, &
                            xc_f03_mgga_vxc_fxc_kxc, &
                            XC_POLARIZED, &
                            XC_UNPOLARIZED, &
                            XC_FAMILY_LDA, &
                            XC_FAMILY_GGA, &
                            XC_FAMILY_MGGA, &
                            XC_FAMILY_HYB_LDA, &
                            XC_FAMILY_HYB_GGA, &
                            XC_FAMILY_HYB_MGGA, &
                            libxc_set_backend, &
                            libxc_func_init_wrapper, &
                            libxc_gga_exc_vxc_wrapper, libxc_gga_exc_wrapper, &
                            libxc_lda_exc_vxc_wrapper, libxc_lda_exc_wrapper, &
                            libxc_mgga_exc_vxc_wrapper, libxc_mgga_exc_wrapper, &
                            libxc_mgga_vxc_wrapper, &
                            xc_libxc_wrap_functional_get_number, &
                            xc_libxc_wrap_info_needs_laplace, &
                            xc_libxc_wrap_info_no_exc, &
                            xc_libxc_wrap_is_under_development, &
                            xc_libxc_check_functional, &
                            xc_libxc_get_reference_length, &
                            xc_libxc_wrap_init_safe, &
                            xc_libxc_wrap_library_reference, &
                            xc_libxc_wrap_version, &
                            xc_libxc_wrap_set_thresholds, &
                            xc_libxc_wrap_info_refs, &
                            xc_f03_func_set_ext_params, &
                            XC_EXCHANGE, XC_CORRELATION, XC_EXCHANGE_CORRELATION, XC_KINETIC
   USE xc_libxc_wrap, ONLY: libxc_is_gpu
   USE xc_libxc_wrap, ONLY: libxc_worker_type, &
                            libxc_worker_alloc, libxc_worker_free_buffers, &
                            libxc_worker_resize_cuda, libxc_worker_eval, libxc_family_flags
#endif
!$ USE OMP_LIB, ONLY: omp_get_max_threads, &
!$                    omp_get_num_threads, &
!$                    omp_get_thread_num

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_libxc'

   PUBLIC :: libxc_spin_unpolarized_info, libxc_spin_unpolarized_eval, &
             libxc_spin_polarized_info, libxc_spin_polarized_eval, &
             libxc_version_info, libxc_library_reference, &
             libxc_get_reference_length, libxc_add_sections, &
             libxc_check_existence_in_libxc, libxc_release_workers

#if defined (__LIBXC)
   PUBLIC :: libxc_set_backend

! **************************************************************************************************
!> \brief All threads' workers for one distinct functional setup, together with
!>        the key that identifies that setup. Each worker (functional object plus
!>        host staging buffers and, on GPU builds, the mirroring device buffers)
!>        is implemented in xc_libxc_wrap; here they are only aggregated.
! **************************************************************************************************
   TYPE libxc_worker_set_type
      INTEGER                                            :: func_id = -1
      INTEGER                                            :: nspin = -1
      INTEGER                                            :: family = -1
      ! highest derivative order the workers' buffers were allocated for
      INTEGER                                            :: grad_deriv = -1
      LOGICAL                                            :: has_laplace = .FALSE.
      LOGICAL                                            :: no_exc = .FALSE.
      LOGICAL                                            :: has_params = .FALSE.
      REAL(KIND=dp)                                      :: epsilon_rho = -1.0_dp
      REAL(KIND=dp)                                      :: epsilon_tau = -1.0_dp
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: params
      CHARACTER(LEN=128), DIMENSION(:), ALLOCATABLE      :: param_names
      TYPE(libxc_worker_type), DIMENSION(:), ALLOCATABLE :: worker
   END TYPE libxc_worker_set_type

   ! Workers are acquired from a serial region only (see libxc_get_workers), so
   ! the cache itself needs no locking.
   TYPE(libxc_worker_set_type), DIMENSION(:), ALLOCATABLE, TARGET, PRIVATE, SAVE :: libxc_cache

   ! Stands in for the derivative arrays a given evaluation does not produce.
   ! Those arguments still have to be associated with something, but they are
   ! never read or written; pointing them at the density instead would make an
   ! output argument alias a read-only input. Grown on demand and kept for the
   ! run rather than allocated per evaluation, and like the cache above it is
   ! only ever touched from a serial region.
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, TARGET, PRIVATE, SAVE :: libxc_unused
#endif

   ! Backend initialization bookkeeping
   LOGICAL, PRIVATE, SAVE :: libxc_initialized = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief This function checks whether a functional name belongs to LibXC
!> \param libxc_params (possible) LibXC input section
!> \return exists whether the functional exists in LibXC
! **************************************************************************************************
   FUNCTION libxc_check_existence_in_libxc(libxc_params) RESULT(exists)

      TYPE(section_vals_type), POINTER, INTENT(IN)         :: libxc_params
      LOGICAL                                  :: exists

#if defined (__LIBXC)
      ! Initialize once on first call
      IF (.NOT. libxc_initialized) THEN
         CALL xc_libxc_wrap_init_safe()
         libxc_initialized = .TRUE.
      END IF

      exists = xc_libxc_check_functional(libxc_params%section%name)
#else
      MARK_USED(libxc_params)
      exists = .FALSE.
#endif

   END FUNCTION libxc_check_existence_in_libxc

! **************************************************************************************************
!> \brief This function returns the maximum length of the reference string for a given LibXC functional
!> \param libxc_params LibXC input section
!> \param lsd spin polarized calculation
!> \return maximum length of the string
! **************************************************************************************************
   FUNCTION libxc_get_reference_length(libxc_params, lsd) RESULT(length)

      TYPE(section_vals_type), POINTER, INTENT(IN)         :: libxc_params
      LOGICAL, INTENT(IN)                      :: lsd
      INTEGER                                  :: length

#if defined (__LIBXC)
      CHARACTER(len=*), PARAMETER :: routineN = 'libxc_get_reference_length'

      CHARACTER(LEN=default_string_length)     :: func_name
      INTEGER                                  :: func_id, handle
      TYPE(xc_f03_func_t)                      :: xc_func
      TYPE(xc_f03_func_info_t)                 :: xc_info

      CALL timeset(routineN, handle)

      func_name = libxc_params%section%name

      func_id = xc_libxc_wrap_functional_get_number(func_name)
      IF (lsd) THEN
         CALL libxc_func_init_wrapper(xc_func, func_id, XC_POLARIZED, force_host=.TRUE.)
      ELSE
         CALL libxc_func_init_wrapper(xc_func, func_id, XC_UNPOLARIZED, force_host=.TRUE.)
      END IF
      xc_info = xc_f03_func_get_info(xc_func)

      length = xc_libxc_get_reference_length(xc_info)

      CALL xc_f03_func_end(xc_func)

      CALL timestop(handle)
#else
      MARK_USED(libxc_params)
      MARK_USED(lsd)
      length = 0
      CPABORT("In order to use LibXC you have to download and install it!")
#endif

   END FUNCTION libxc_get_reference_length

! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE libxc_add_sections(section)

      TYPE(section_type), POINTER, INTENT(IN) :: section

#if defined (__LIBXC)
      CHARACTER(len=*), PARAMETER :: routineN = 'libxc_add_sections'

      TYPE(section_type), POINTER :: subsection
      TYPE(keyword_type), POINTER :: keyword
      INTEGER :: handle, no_func, len_name, ii, func_id, n_param, iparam
      REAL(KIND=C_DOUBLE) :: default_val
      CHARACTER(LEN=128) :: func_name, param_name, param_descr, description
      CHARACTER(LEN=2*default_string_length) :: warning
      INTEGER(KIND=C_INT), DIMENSION(:), ALLOCATABLE :: func_ids
      TYPE(xc_f03_func_t)                      :: xc_func
      TYPE(xc_f03_func_info_t)                 :: xc_info

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(section))
      NULLIFY (subsection, keyword)

      no_func = xc_f03_number_of_functionals()
      len_name = xc_f03_maximum_name_length()

      ALLOCATE (func_ids(no_func))

      CALL xc_f03_available_functional_numbers(func_ids)

      DO ii = 1, no_func

         func_id = func_ids(ii)
         IF (ii > 1) THEN
            IF (func_id == func_ids(ii - 1)) CYCLE
         END IF
         CALL libxc_func_init_wrapper(xc_func, func_id, XC_UNPOLARIZED, force_host=.TRUE.)
         xc_info = xc_f03_func_get_info(xc_func)

         func_name = xc_f03_functional_get_name(func_id)
         description = xc_f03_func_info_get_name(xc_info)
         n_param = xc_f03_func_info_get_n_ext_params(xc_info)

         NULLIFY (subsection)
         CALL section_create(subsection, __LOCATION__, name=TRIM(func_name), description=TRIM(description), &
                             n_keywords=2 + n_param, n_subsections=0, repeats=.FALSE.)

         IF (description(1:1) == "_") THEN
            warning = " This parameter is an internal parameter of the functional. Changing this "// &
                      "parameter effectively changes the functional."
         ELSE
            warning = " "
         END IF

         NULLIFY (keyword)
         CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                             description="Activates the functional."//TRIM(warning), &
                             lone_keyword_l_val=.TRUE., default_l_val=.FALSE.)
         CALL section_add_keyword(subsection, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="SCALE", description="Scales this functional", &
                             default_r_val=1.0_dp)
         CALL section_add_keyword(subsection, keyword)
         CALL keyword_release(keyword)

         DO iparam = 1, n_param
            param_name = xc_f03_func_info_get_ext_params_name(xc_info, iparam - 1)
            param_descr = xc_f03_func_info_get_ext_params_description(xc_info, iparam - 1)
            default_val = xc_f03_func_info_get_ext_params_default_value(xc_info, iparam - 1)
            NULLIFY (keyword)
            CALL keyword_create(keyword, __LOCATION__, name=TRIM(param_name), &
                                description=TRIM(param_descr), default_r_val=default_val)
            CALL section_add_keyword(subsection, keyword)
            CALL keyword_release(keyword)
         END DO

         CALL section_add_subsection(section, subsection)
         CALL section_release(subsection)

         CALL xc_f03_func_end(xc_func)

      END DO

      DEALLOCATE (func_ids)

      CALL timestop(handle)
#else
      MARK_USED(section)

#endif

   END SUBROUTINE libxc_add_sections

! **************************************************************************************************
!> \brief info about the functional from libxc
!> \param libxc_params input parameter (functional name, scaling and parameters)
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv maximum implemented derivative of the xc functional
!> \param print_warn whether to print warning about development status of a functional
!> \param func_name_override optional LibXC functional name overriding the section name
!> \author F. Tran
! **************************************************************************************************
   SUBROUTINE libxc_spin_unpolarized_info(libxc_params, reference, shortform, needs, max_deriv, print_warn, &
                                          func_name_override)

      TYPE(section_vals_type), POINTER         :: libxc_params
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
      TYPE(xc_rho_cflags_type), &
         INTENT(inout), OPTIONAL               :: needs
      INTEGER, INTENT(out), OPTIONAL           :: max_deriv
      LOGICAL, INTENT(IN), OPTIONAL            :: print_warn
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: func_name_override

#if defined (__LIBXC)
      CHARACTER(LEN=128)                       :: s1, s2
      CHARACTER(LEN=default_string_length)     :: func_name
      INTEGER                                  :: func_id
      REAL(KIND=dp)                            :: func_scale
      TYPE(xc_f03_func_t)                      :: xc_func
      TYPE(xc_f03_func_info_t)                 :: xc_info

      IF (PRESENT(func_name_override)) THEN
         func_name = func_name_override
         func_scale = 1.0_dp
      ELSE
         func_name = libxc_params%section%name
         CALL section_vals_val_get(libxc_params, "scale", r_val=func_scale)
      END IF

      CALL cite_reference(Marques2012)
      CALL cite_reference(Lehtola2018)

      IF (ABS(func_scale - 1.0_dp) < 1.0e-10_dp) func_scale = 1.0_dp

      func_id = xc_libxc_wrap_functional_get_number(func_name)
      CALL libxc_func_init_wrapper(xc_func, func_id, XC_UNPOLARIZED, force_host=.TRUE.)
      xc_info = xc_f03_func_get_info(xc_func)

      s1 = xc_f03_func_info_get_name(xc_info)
      SELECT CASE (xc_f03_func_info_get_kind(xc_info))
      CASE (XC_EXCHANGE); WRITE (s2, '(a)') "exchange"
      CASE (XC_CORRELATION); WRITE (s2, '(a)') "correlation"
      CASE (XC_EXCHANGE_CORRELATION); WRITE (s2, '(a)') "exchange-correlation"
      CASE (XC_KINETIC); WRITE (s2, '(a)') "kinetic"
      CASE default
         CPABORT(TRIM(func_name)//": this XC_KIND is currently not supported.")
      END SELECT
      IF (PRESENT(shortform)) THEN
         shortform = TRIM(s1)//' ('//TRIM(s2)//')'
      END IF
      IF (PRESENT(reference)) THEN
         CALL xc_libxc_wrap_info_refs(xc_info, XC_UNPOLARIZED, func_scale, reference)
      END IF
      IF (PRESENT(needs)) THEN
         SELECT CASE (xc_f03_func_info_get_family(xc_info))
         CASE (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
            needs%rho = .TRUE.
         CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            needs%rho = .TRUE.
            needs%norm_drho = .TRUE.
         CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            needs%rho = .TRUE.
            needs%norm_drho = .TRUE.
            needs%tau = .TRUE.
            needs%laplace_rho = xc_libxc_wrap_info_needs_laplace(xc_info)
         CASE default
            CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
         END SELECT
      END IF
      IF (PRESENT(max_deriv)) THEN
         SELECT CASE (xc_f03_func_info_get_family(xc_info))
         CASE (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
            max_deriv = 3
         CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            max_deriv = 2
         CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            max_deriv = 2
         CASE default
            CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
         END SELECT
      END IF
      IF (PRESENT(print_warn)) THEN
         IF (print_warn .AND. xc_libxc_wrap_is_under_development(xc_info)) THEN
            CPWARN(TRIM(func_name)//" is under development. Use with caution.")
         END IF
      END IF

      CALL xc_f03_func_end(xc_func)
#else
      MARK_USED(libxc_params)
      MARK_USED(reference)
      MARK_USED(shortform)
      MARK_USED(needs)
      MARK_USED(max_deriv)
      MARK_USED(print_warn)
      MARK_USED(func_name_override)

      CALL cp_abort(__LOCATION__, "Unknown functional! If you are asking "// &
                    "for a functional of the LibXC library, "// &
                    "you have to download and install the library!")
#endif

   END SUBROUTINE libxc_spin_unpolarized_info

! **************************************************************************************************
!> \brief info about the functional from libxc
!> \param libxc_params input parameter (functional name, scaling and parameters)
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv maximum implemented derivative of the xc functional
!> \param print_warn whether to print warning about development status of a functional
!> \param func_name_override optional LibXC functional name overriding the section name
!> \author F. Tran
! **************************************************************************************************
   SUBROUTINE libxc_spin_polarized_info(libxc_params, reference, shortform, needs, max_deriv, print_warn, &
                                        func_name_override)

      TYPE(section_vals_type), POINTER         :: libxc_params
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
      TYPE(xc_rho_cflags_type), &
         INTENT(inout), OPTIONAL               :: needs
      INTEGER, INTENT(out), OPTIONAL           :: max_deriv
      LOGICAL, INTENT(IN), OPTIONAL            :: print_warn
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: func_name_override

#if defined (__LIBXC)
      CHARACTER(LEN=128)                       :: s1, s2
      CHARACTER(LEN=default_string_length)     :: func_name
      INTEGER                                  :: func_id
      REAL(KIND=dp)                            :: func_scale
      TYPE(xc_f03_func_t)                      :: xc_func
      TYPE(xc_f03_func_info_t)                 :: xc_info

      IF (PRESENT(func_name_override)) THEN
         func_name = func_name_override
         func_scale = 1.0_dp
      ELSE
         func_name = libxc_params%section%name
         CALL section_vals_val_get(libxc_params, "scale", r_val=func_scale)
      END IF

      CALL cite_reference(Marques2012)
      CALL cite_reference(Lehtola2018)

      IF (ABS(func_scale - 1.0_dp) < 1.0e-10_dp) func_scale = 1.0_dp

      func_id = xc_libxc_wrap_functional_get_number(func_name)
      CALL libxc_func_init_wrapper(xc_func, func_id, XC_POLARIZED, force_host=.TRUE.)
      xc_info = xc_f03_func_get_info(xc_func)

      s1 = xc_f03_func_info_get_name(xc_info)
      SELECT CASE (xc_f03_func_info_get_kind(xc_info))
      CASE (XC_EXCHANGE); WRITE (s2, '(a)') "exchange"
      CASE (XC_CORRELATION); WRITE (s2, '(a)') "correlation"
      CASE (XC_EXCHANGE_CORRELATION); WRITE (s2, '(a)') "exchange-correlation"
      CASE (XC_KINETIC); WRITE (s2, '(a)') "kinetic"
      CASE default
         CPABORT(TRIM(func_name)//": this XC_KIND is currently not supported.")
      END SELECT
      IF (PRESENT(shortform)) THEN
         shortform = TRIM(s1)//' ('//TRIM(s2)//')'
      END IF
      IF (PRESENT(reference)) THEN
         CALL xc_libxc_wrap_info_refs(xc_info, XC_POLARIZED, func_scale, reference)
      END IF
      IF (PRESENT(needs)) THEN
         SELECT CASE (xc_f03_func_info_get_family(xc_info))
         CASE (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
            needs%rho_spin = .TRUE.
         CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            needs%rho_spin = .TRUE.
            needs%norm_drho = .TRUE.
            needs%norm_drho_spin = .TRUE.
         CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            needs%rho_spin = .TRUE.
            needs%norm_drho = .TRUE.
            needs%norm_drho_spin = .TRUE.
            needs%tau_spin = .TRUE.
            needs%laplace_rho_spin = xc_libxc_wrap_info_needs_laplace(xc_info)
         CASE default
            CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
         END SELECT
      END IF
      IF (PRESENT(max_deriv)) THEN
         SELECT CASE (xc_f03_func_info_get_family(xc_info))
         CASE (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
            max_deriv = 3
         CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            max_deriv = 2
         CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            max_deriv = 2
         CASE default
            CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
         END SELECT
      END IF
      IF (PRESENT(print_warn)) THEN
         IF (print_warn .AND. xc_libxc_wrap_is_under_development(xc_info)) THEN
            CPWARN(TRIM(func_name)//" is under development. Use with caution.")
         END IF
      END IF

      CALL xc_f03_func_end(xc_func)
#else
      MARK_USED(libxc_params)
      MARK_USED(reference)
      MARK_USED(shortform)
      MARK_USED(needs)
      MARK_USED(max_deriv)
      MARK_USED(print_warn)
      MARK_USED(func_name_override)

      CALL cp_abort(__LOCATION__, "Unknown functional! If you are "// &
                    "asking for a functional of the LibXC library, "// &
                    "you have to download and install the library!")
#endif

   END SUBROUTINE libxc_spin_polarized_info

! **************************************************************************************************
!> \brief info about the LibXC version
!> \param version version string reported by the library at run time
!> \param compiled_version version of the LibXC headers CP2K was compiled against
!> \author A. Gloess (agloess)
! **************************************************************************************************
   SUBROUTINE libxc_version_info(version, compiled_version)
      CHARACTER(LEN=*), INTENT(OUT)      :: version ! the string that is output
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: compiled_version

#if defined (__LIBXC)
      CALL xc_libxc_wrap_version(version, compiled_version)
#else
      version = "none"
      IF (PRESENT(compiled_version)) compiled_version = "none"
      CPABORT("In order to use libxc you need to download and install it")
#endif

   END SUBROUTINE libxc_version_info

! **************************************************************************************************
!> \brief Returns the citation LibXC asks for the library itself.
!> \param reference bibliographic reference
!> \param doi digital object identifier of that reference
!> \author S. Lehtola
! **************************************************************************************************
   SUBROUTINE libxc_library_reference(reference, doi)
      CHARACTER(LEN=*), INTENT(OUT)      :: reference, doi

#if defined (__LIBXC)
      CALL xc_libxc_wrap_library_reference(reference, doi)
#else
      reference = "none"
      doi = "none"
      CPABORT("In order to use libxc you need to download and install it")
#endif

   END SUBROUTINE libxc_library_reference

! **************************************************************************************************
!> \brief Reads the external parameters of a LibXC functional from its input section.
!> \param libxc_params LibXC input section
!> \param param_names names of the external parameters, as reported by LibXC
!> \param params values read from the input
!> \author S. Lehtola
! **************************************************************************************************
#if defined (__LIBXC)
   SUBROUTINE libxc_read_params(libxc_params, param_names, params)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: libxc_params
      CHARACTER(LEN=128), DIMENSION(:), INTENT(IN)       :: param_names
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: params

      INTEGER                                            :: i

      DO i = 1, SIZE(param_names)
         CALL section_vals_val_get(libxc_params, TRIM(param_names(i)), r_val=params(i))
      END DO

   END SUBROUTINE libxc_read_params

! **************************************************************************************************
!> \brief Makes sure the placeholder for unproduced derivatives holds npoints values.
!> \param npoints number of grid points in this evaluation
!> \author S. Lehtola
! **************************************************************************************************
   SUBROUTINE libxc_reserve_unused(npoints)
      INTEGER, INTENT(IN)                                :: npoints

      IF (ALLOCATED(libxc_unused)) THEN
         IF (SIZE(libxc_unused) >= npoints) RETURN
         DEALLOCATE (libxc_unused)
      END IF
      ALLOCATE (libxc_unused(npoints))

   END SUBROUTINE libxc_reserve_unused

! **************************************************************************************************
!> \brief Returns the cache slot holding the per-thread workers for a functional.
!> \param func_name LibXC functional name
!> \param nspin XC_UNPOLARIZED or XC_POLARIZED
!> \param libxc_params LibXC input section, or null when no parameters are to be set
!> \param epsilon_rho density cutoff
!> \param epsilon_tau kinetic energy density cutoff
!> \param grad_deriv highest derivative order this evaluation needs
!> \return index into libxc_cache
!> \author S. Lehtola
!> \note Must be called from a serial region: it is the only writer of the cache,
!>       and keeping it serial is what makes the cache lock-free. The workers it
!>       returns are then used concurrently, one per thread.
! **************************************************************************************************
   FUNCTION libxc_get_workers(func_name, nspin, libxc_params, epsilon_rho, epsilon_tau, &
                              grad_deriv) RESULT(idx)
      CHARACTER(LEN=*), INTENT(IN)                       :: func_name
      INTEGER, INTENT(IN)                                :: nspin
      TYPE(section_vals_type), INTENT(IN), POINTER       :: libxc_params
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_rho, epsilon_tau
      INTEGER, INTENT(IN)                                :: grad_deriv
      INTEGER                                            :: idx

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'libxc_get_workers'

      INTEGER                                            :: func_id, handle, i, ithread, n_params, &
                                                            nthreads
      LOGICAL                                            :: has_params
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: params
      TYPE(libxc_worker_set_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: tmp_cache

      CALL timeset(routineN, handle)

      func_id = xc_libxc_wrap_functional_get_number(func_name)
      has_params = ASSOCIATED(libxc_params)

      nthreads = 1
!$    nthreads = omp_get_max_threads()

      IF (.NOT. ALLOCATED(libxc_cache)) ALLOCATE (libxc_cache(0))

      idx = 0
      DO i = 1, SIZE(libxc_cache)
         IF (libxc_cache(i)%func_id == func_id .AND. &
             libxc_cache(i)%nspin == nspin .AND. &
             (libxc_cache(i)%has_params .EQV. has_params)) THEN
            IF (libxc_cache(i)%epsilon_rho == epsilon_rho .AND. &
                libxc_cache(i)%epsilon_tau == epsilon_tau .AND. &
                SIZE(libxc_cache(i)%worker) >= nthreads) THEN
               idx = i
               EXIT
            END IF
         END IF
      END DO

      ! A set built for a lower derivative order carries buffers that are too
      ! narrow for this call. Widening them in place keeps one set per functional
      ! setup: the SCF fills the cache at first order and a later response
      ! calculation grows that same slot to third order.
      IF (idx /= 0) THEN
         IF (libxc_cache(idx)%grad_deriv < grad_deriv) THEN
            DO ithread = 1, SIZE(libxc_cache(idx)%worker)
               CALL libxc_worker_free_buffers(libxc_cache(idx)%worker(ithread))
               CALL libxc_worker_alloc(libxc_cache(idx)%worker(ithread), &
                                       libxc_cache(idx)%family, nspin, grad_deriv)
            END DO
            libxc_cache(idx)%grad_deriv = grad_deriv
         END IF
      END IF

      IF (idx == 0) THEN
         ! Grow the cache by one slot and build the workers for it.
         CALL MOVE_ALLOC(libxc_cache, tmp_cache)
         ALLOCATE (libxc_cache(SIZE(tmp_cache) + 1))
         DO i = 1, SIZE(tmp_cache)
            CALL libxc_move_worker_set(tmp_cache(i), libxc_cache(i))
         END DO
         DEALLOCATE (tmp_cache)
         idx = SIZE(libxc_cache)

         ASSOCIATE (ws => libxc_cache(idx))
            ws%func_id = func_id
            ws%nspin = nspin
            ws%has_params = has_params
            ws%epsilon_rho = epsilon_rho
            ws%epsilon_tau = epsilon_tau
            ws%grad_deriv = grad_deriv
            ALLOCATE (ws%worker(nthreads))

            DO ithread = 1, nthreads
               CALL libxc_func_init_wrapper(ws%worker(ithread)%func, func_id, nspin)
               ws%worker(ithread)%info = xc_f03_func_get_info(ws%worker(ithread)%func)
               ws%worker(ithread)%is_init = .TRUE.
            END DO

            ws%family = xc_f03_func_info_get_family(ws%worker(1)%info)
            ws%has_laplace = xc_libxc_wrap_info_needs_laplace(ws%worker(1)%info)
            ws%no_exc = xc_libxc_wrap_info_no_exc(ws%worker(1)%info)

            ! Remember the parameter names so that later evaluations can re-read
            ! the input without another functional object to ask for them.
            n_params = 0
            IF (has_params) n_params = xc_f03_func_info_get_n_ext_params(ws%worker(1)%info)
            ALLOCATE (ws%params(n_params))
            ALLOCATE (ws%param_names(n_params))
            DO i = 1, n_params
               ws%param_names(i) = xc_f03_func_info_get_ext_params_name(ws%worker(1)%info, i - 1)
            END DO
            IF (n_params > 0) THEN
               CALL libxc_read_params(libxc_params, ws%param_names, ws%params)
               DO ithread = 1, nthreads
                  CALL xc_f03_func_set_ext_params(ws%worker(ithread)%func, ws%params)
               END DO
            END IF

            DO ithread = 1, nthreads
               CALL xc_libxc_wrap_set_thresholds(ws%worker(ithread)%func, ws%worker(ithread)%info, &
                                                 epsilon_rho, epsilon_tau)
               CALL libxc_worker_alloc(ws%worker(ithread), ws%family, nspin, grad_deriv)
            END DO
         END ASSOCIATE
      ELSE
         ! Cached slot: the input section may in principle carry different
         ! parameter values than the ones the workers were last set up with, so
         ! re-read them and push them through only when they actually changed.
         ASSOCIATE (ws => libxc_cache(idx))
            n_params = SIZE(ws%params)
            IF (n_params > 0) THEN
               ALLOCATE (params(n_params))
               CALL libxc_read_params(libxc_params, ws%param_names, params)
               IF (ANY(params /= ws%params)) THEN
                  ws%params(:) = params(:)
                  DO ithread = 1, SIZE(ws%worker)
                     CALL xc_f03_func_set_ext_params(ws%worker(ithread)%func, ws%params)
                     ! re-assert the cutoffs: setting parameters runs the
                     ! functional's own callback, which is free to rebuild the
                     ! component functionals a mixed functional is made of
                     CALL xc_libxc_wrap_set_thresholds(ws%worker(ithread)%func, &
                                                       ws%worker(ithread)%info, &
                                                       ws%epsilon_rho, ws%epsilon_tau)
                  END DO
               END IF
               DEALLOCATE (params)
            END IF
         END ASSOCIATE
      END IF

      CALL timestop(handle)

   END FUNCTION libxc_get_workers

! **************************************************************************************************
!> \brief Moves a worker set to a new cache slot without touching LibXC.
!> \param from source slot, left empty
!> \param to destination slot
!> \author S. Lehtola
!> \note Used when the cache array grows. The LibXC objects are opaque handles, so
!>       they survive being moved; only the Fortran allocatables need transferring.
! **************************************************************************************************
   SUBROUTINE libxc_move_worker_set(from, to)
      TYPE(libxc_worker_set_type), INTENT(INOUT)         :: from, to

      to%func_id = from%func_id
      to%nspin = from%nspin
      to%family = from%family
      to%grad_deriv = from%grad_deriv
      to%has_laplace = from%has_laplace
      to%no_exc = from%no_exc
      to%has_params = from%has_params
      to%epsilon_rho = from%epsilon_rho
      to%epsilon_tau = from%epsilon_tau
      IF (ALLOCATED(from%params)) CALL MOVE_ALLOC(from%params, to%params)
      IF (ALLOCATED(from%param_names)) CALL MOVE_ALLOC(from%param_names, to%param_names)
      IF (ALLOCATED(from%worker)) CALL MOVE_ALLOC(from%worker, to%worker)

   END SUBROUTINE libxc_move_worker_set
#endif

! **************************************************************************************************
!> \brief Destroys all cached LibXC functional objects.
!> \author S. Lehtola
!> \note Call once at the end of a run. The workers are deliberately kept alive
!>       across SCF steps and across calculations, so nothing else releases them.
! **************************************************************************************************
   SUBROUTINE libxc_release_workers()

#if defined (__LIBXC)
      INTEGER                                            :: i, ithread

      IF (.NOT. ALLOCATED(libxc_cache)) RETURN

      DO i = 1, SIZE(libxc_cache)
         IF (.NOT. ALLOCATED(libxc_cache(i)%worker)) CYCLE
         DO ithread = 1, SIZE(libxc_cache(i)%worker)
            IF (libxc_cache(i)%worker(ithread)%is_init) THEN
               CALL xc_f03_func_end(libxc_cache(i)%worker(ithread)%func)
               libxc_cache(i)%worker(ithread)%is_init = .FALSE.
            END IF
            ! Free the staging buffers and any device buffers they mirror.
            CALL libxc_worker_free_buffers(libxc_cache(i)%worker(ithread))
         END DO
      END DO
      DEALLOCATE (libxc_cache)
      IF (ALLOCATED(libxc_unused)) DEALLOCATE (libxc_unused)
#endif

   END SUBROUTINE libxc_release_workers

! **************************************************************************************************
!> \brief evaluates the functional from libxc
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv maximum degree of the derivative that should be evaluated;
!>        all derivatives up to the given degree are evaluated, in a single
!>        LibXC call per block of grid points
!> \param libxc_params input parameter (functional name, scaling and parameters)
!> \param func_name_override optional LibXC functional name overriding the section name
!> \author F. Tran
! **************************************************************************************************
   SUBROUTINE libxc_spin_unpolarized_eval(rho_set, deriv_set, grad_deriv, libxc_params, func_name_override)

      TYPE(xc_rho_set_type), INTENT(IN)        :: rho_set
      TYPE(xc_derivative_set_type), INTENT(IN) :: deriv_set
      INTEGER, INTENT(in)                      :: grad_deriv
      TYPE(section_vals_type), POINTER         :: libxc_params
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: func_name_override

#if defined (__LIBXC)
      CHARACTER(len=*), PARAMETER :: routineN = 'libxc_spin_unpolarized_eval'

      CHARACTER(LEN=default_string_length)     :: func_name
      INTEGER                                  :: handle, iw, npoints
      INTEGER, DIMENSION(2, 3)                 :: bo
      LOGICAL                                  :: has_laplace
      REAL(KIND=dp)                            :: epsilon_rho, epsilon_tau, func_scale
      TYPE(libxc_worker_set_type), POINTER     :: workers
      TYPE(section_vals_type), POINTER         :: no_params
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), POINTER :: dummy, e_0, e_laplace_rho, &
                                                                e_laplace_rho_laplace_rho, e_laplace_rho_tau, e_ndrho, &
                                                              e_ndrho_laplace_rho, e_ndrho_ndrho, e_ndrho_rho, e_ndrho_tau, e_rho, &
                                                                e_rho_laplace_rho, e_rho_rho, e_rho_rho_rho, e_rho_tau, e_tau, &
                                                                e_tau_tau, laplace_rho, norm_drho, rho, tau
      ! derivatives w.r.t. the reduced gradient gamma = |grad rho|^2
      #:for descs, arr, idx in umgga_all
         REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), POINTER :: e_${'_'.join(descs)}$
      #:endfor
      TYPE(xc_derivative_type), POINTER        :: deriv
      TYPE(xc_f03_func_info_t)                 :: xc_info

      CALL timeset(routineN, handle)

      ! Only "everything up to grad_deriv" is supported. Asking for a single
      ! derivative order in isolation is rejected here rather than further down,
      ! where it would quietly evaluate nothing at all.
      IF (grad_deriv > 3) THEN
         CPABORT("derivatives larger than 3 not implemented")
      END IF

      has_laplace = .FALSE.
      NULLIFY (dummy)
      NULLIFY (workers)
      NULLIFY (rho, norm_drho, laplace_rho, tau)

      IF (PRESENT(func_name_override)) THEN
         func_name = func_name_override
         func_scale = 1.0_dp
      ELSE
         func_name = libxc_params%section%name
         CALL section_vals_val_get(libxc_params, "scale", r_val=func_scale)
      END IF

      IF (ABS(func_scale - 1.0_dp) < 1.0e-10_dp) func_scale = 1.0_dp

      CALL xc_rho_set_get(rho_set, can_return_null=.TRUE., &
                          rho=rho, norm_drho=norm_drho, laplace_rho=laplace_rho, &
                          rho_cutoff=epsilon_rho, tau_cutoff=epsilon_tau, &
                          tau=tau, local_bounds=bo)

      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)

      ! One functional object per thread, kept across calls. The cutoffs are part
      ! of the worker setup because LibXC does the screening itself.
      IF (PRESENT(func_name_override)) THEN
         ! an overriding name carries no input section, so no parameters are set
         NULLIFY (no_params)
         iw = libxc_get_workers(func_name, XC_UNPOLARIZED, no_params, epsilon_rho, epsilon_tau, &
                                grad_deriv)
      ELSE
         iw = libxc_get_workers(func_name, XC_UNPOLARIZED, libxc_params, epsilon_rho, epsilon_tau, &
                                grad_deriv)
      END IF
      workers => libxc_cache(iw)
      xc_info = workers%worker(1)%info
      has_laplace = workers%has_laplace

      ! see libxc_unused: the arguments the requested order does not produce are
      ! never touched, but they still have to point somewhere
      CALL libxc_reserve_unused(npoints)
      dummy(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)) => libxc_unused(1:npoints)

      ! due to assumed shape array usage in next routine
      IF (.NOT. ASSOCIATED(norm_drho)) norm_drho => rho
      IF (.NOT. ASSOCIATED(tau)) tau => rho

      ! only some MGGA functionals really need the Laplacian,
      ! all others can work with rho (read-only) as dummy
      IF (.NOT. has_laplace) laplace_rho => rho

      e_0 => dummy
      e_rho => dummy
      e_ndrho => dummy
      e_laplace_rho => dummy
      e_tau => dummy
      e_rho_rho => dummy
      e_ndrho_rho => dummy
      e_ndrho_ndrho => dummy
      e_rho_laplace_rho => dummy
      e_rho_tau => dummy
      e_ndrho_laplace_rho => dummy
      e_ndrho_tau => dummy
      e_laplace_rho_laplace_rho => dummy
      e_laplace_rho_tau => dummy
      e_tau_tau => dummy
      e_rho_rho_rho => dummy
      #:for descs, arr, idx in umgga_all
         e_${'_'.join(descs)}$ => dummy
      #:endfor

      deriv => xc_dset_get_derivative(deriv_set, [INTEGER::], &
                                      allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv, deriv_data=e_0)
      IF (grad_deriv >= 1) THEN
         SELECT CASE (xc_f03_func_info_get_family(xc_info))
         CASE (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho)
         CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho)
            #:for descs, arr, idx in ugga_fill_1
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
         CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho)
            #:for descs, arr, idx in ugga_fill_1
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
            deriv => xc_dset_get_derivative(deriv_set, [deriv_tau], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_tau)
            IF (has_laplace) THEN
               deriv => xc_dset_get_derivative(deriv_set, [deriv_laplace_rho], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_laplace_rho)
            END IF
         CASE default
            CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
         END SELECT
      END IF
      IF (grad_deriv >= 2) THEN
         SELECT CASE (xc_f03_func_info_get_family(xc_info))
         CASE (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho)
         CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_rho)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_norm_drho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_ndrho)
            #:for descs, arr, idx in ugga_fill_2
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
         CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            ! not implemented ...

            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_rho)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_norm_drho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_ndrho)
            #:for descs, arr, idx in ugga_fill_2
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_tau], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_tau)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_tau], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_tau)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_tau, deriv_tau], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_tau_tau)
            #:for descs, arr, idx in umgga_extra_2_nl
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
            IF (has_laplace) THEN
               deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_laplace_rho], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_rho_laplace_rho)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_laplace_rho], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_ndrho_laplace_rho)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_laplace_rho, deriv_laplace_rho], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_laplace_rho_laplace_rho)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_laplace_rho, deriv_tau], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_laplace_rho_tau)
            END IF
         CASE default
            CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
         END SELECT
      END IF
      IF (grad_deriv >= 3) THEN
         SELECT CASE (xc_f03_func_info_get_family(xc_info))
         CASE (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
            IF (has_laplace) THEN
               #:for descs, arr, idx in umgga_fill_3_l
                  deriv => xc_dset_get_derivative(deriv_set, &
                                                  [${', '.join('deriv_' + d for d in descs)}$], &
                                                  allocate_deriv=.TRUE.)
                  CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
               #:endfor
            END IF
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_rho, deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho_rho)
         CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            IF (has_laplace) THEN
               #:for descs, arr, idx in umgga_fill_3_l
                  deriv => xc_dset_get_derivative(deriv_set, &
                                                  [${', '.join('deriv_' + d for d in descs)}$], &
                                                  allocate_deriv=.TRUE.)
                  CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
               #:endfor
            END IF
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_rho, deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho_rho)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_rho, deriv_gamma], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho_gamma)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_gamma, deriv_gamma], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_gamma_gamma)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_gamma, deriv_gamma, deriv_gamma], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_gamma_gamma_gamma)
         CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            #:for descs, arr, idx in umgga_fill_3_nl
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
            IF (has_laplace) THEN
               #:for descs, arr, idx in umgga_fill_3_l
                  deriv => xc_dset_get_derivative(deriv_set, &
                                                  [${', '.join('deriv_' + d for d in descs)}$], &
                                                  allocate_deriv=.TRUE.)
                  CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
               #:endfor
            END IF
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_rho, deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho_rho)
         CASE default
            CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
         END SELECT
      END IF

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP SHARED(rho,norm_drho,laplace_rho,tau,e_0,e_rho,e_ndrho,e_laplace_rho),&
!$OMP SHARED(e_tau,e_rho_rho,e_ndrho_rho,e_ndrho_ndrho,e_rho_laplace_rho),&
!$OMP SHARED(e_rho_tau,e_ndrho_laplace_rho,e_ndrho_tau,e_laplace_rho_laplace_rho),&
!$OMP SHARED(e_laplace_rho_tau,e_tau_tau,e_rho_rho_rho),&
      #:for descs, arr, idx in umgga_all
!$OMP SHARED(e_${'_'.join(descs)}$),&
      #:endfor
!$OMP SHARED(grad_deriv,npoints),&
!$OMP SHARED(func_name,func_scale,workers)

      CALL libxc_spin_unpolarized_calc(rho=rho, norm_drho=norm_drho, &
                                       laplace_rho=laplace_rho, tau=tau, &
                                       e_0=e_0, e_rho=e_rho, e_ndrho=e_ndrho, e_laplace_rho=e_laplace_rho, &
                                       e_tau=e_tau, e_rho_rho=e_rho_rho, e_ndrho_rho=e_ndrho_rho, &
                                       e_ndrho_ndrho=e_ndrho_ndrho, e_rho_laplace_rho=e_rho_laplace_rho, &
                                       e_rho_tau=e_rho_tau, e_ndrho_laplace_rho=e_ndrho_laplace_rho, &
                                       e_ndrho_tau=e_ndrho_tau, e_laplace_rho_laplace_rho=e_laplace_rho_laplace_rho, &
                                       e_laplace_rho_tau=e_laplace_rho_tau, e_tau_tau=e_tau_tau, &
                                       e_rho_rho_rho=e_rho_rho_rho, &
                                       #:for descs, arr, idx in umgga_all
                                       e_${'_'.join(descs)}$=e_${'_'.join(descs)}$, &
                                       #:endfor
                                       grad_deriv=grad_deriv, npoints=npoints, &
                                       func_name=func_name, sc=func_scale, workers=workers)

!$OMP END PARALLEL

      NULLIFY (dummy)
      NULLIFY (workers)

      CALL timestop(handle)
#else
      MARK_USED(rho_set)
      MARK_USED(deriv_set)
      MARK_USED(grad_deriv)
      MARK_USED(libxc_params)
      MARK_USED(func_name_override)
      CALL cp_abort(__LOCATION__, "Unknown functional! If you are asking "// &
                    "for a functional of the LibXC library, "// &
                    "you have to download and install the library!")
#endif
   END SUBROUTINE libxc_spin_unpolarized_eval

! **************************************************************************************************
!> \brief evaluates the functional from libxc
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param grad_deriv degree of the derivative that should be evaluated;
!>        all derivatives up to the given degree are evaluated, in a single
!>        LibXC call per block of grid points
!> \param libxc_params input parameter (functional name, scaling and parameters)
!> \param func_name_override optional LibXC functional name overriding the section name
!> \author F. Tran
! **************************************************************************************************
   SUBROUTINE libxc_spin_polarized_eval(rho_set, deriv_set, grad_deriv, libxc_params, func_name_override)

      TYPE(xc_rho_set_type), INTENT(IN)        :: rho_set
      TYPE(xc_derivative_set_type), INTENT(IN) :: deriv_set
      INTEGER, INTENT(in)                      :: grad_deriv
      TYPE(section_vals_type), POINTER         :: libxc_params
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: func_name_override

#if defined (__LIBXC)
      CHARACTER(len=*), PARAMETER :: routineN = 'libxc_spin_polarized_eval'

      CHARACTER(LEN=default_string_length)     :: func_name
      INTEGER                                  :: handle, iw, npoints
      INTEGER, DIMENSION(2, 3)                 :: bo
      LOGICAL                                  :: has_laplace
      REAL(KIND=dp)                            :: epsilon_rho, epsilon_tau, func_scale
      TYPE(libxc_worker_set_type), POINTER     :: workers
      TYPE(section_vals_type), POINTER         :: no_params
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), POINTER :: dummy, e_0, e_laplace_rhoa, &
                                                                e_laplace_rhoa_laplace_rhoa, e_laplace_rhoa_laplace_rhob, &
                                                                e_laplace_rhoa_tau_a, e_laplace_rhoa_tau_b, e_laplace_rhob, &
                                                                e_laplace_rhob_laplace_rhob, e_laplace_rhob_tau_a, &
                                                                e_laplace_rhob_tau_b, e_ndrho, e_ndrho_laplace_rhoa, &
                                                              e_ndrho_laplace_rhob, e_ndrho_ndrho, e_ndrho_ndrhoa, e_ndrho_ndrhob, &
                                                               e_ndrho_rhoa, e_ndrho_rhob, e_ndrho_tau_a, e_ndrho_tau_b, e_ndrhoa, &
                                                                e_ndrhoa_laplace_rhoa, e_ndrhoa_laplace_rhob, e_ndrhoa_ndrhoa, &
                                                                e_ndrhoa_ndrhob, e_ndrhoa_rhoa, e_ndrhoa_rhob, e_ndrhoa_tau_a, &
                                                                e_ndrhoa_tau_b, e_ndrhob
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), POINTER :: e_ndrhob_laplace_rhoa, &
                                                             e_ndrhob_laplace_rhob, e_ndrhob_ndrhob, e_ndrhob_rhoa, e_ndrhob_rhob, &
                                                                e_ndrhob_tau_a, e_ndrhob_tau_b, e_rhoa, e_rhoa_laplace_rhoa, &
                                                             e_rhoa_laplace_rhob, e_rhoa_rhoa, e_rhoa_rhoa_rhoa, e_rhoa_rhoa_rhob, &
                                                                e_rhoa_rhob, e_rhoa_rhob_rhob, e_rhoa_tau_a, e_rhoa_tau_b, e_rhob, &
                                                                e_rhob_laplace_rhoa, e_rhob_laplace_rhob, e_rhob_rhob, &
                                                             e_rhob_rhob_rhob, e_rhob_tau_a, e_rhob_tau_b, e_tau_a, e_tau_a_tau_a, &
                                                                e_tau_a_tau_b, e_tau_b, e_tau_b_tau_b, laplace_rhoa, laplace_rhob, &
                                                                norm_drho, norm_drhoa, norm_drhob, rhoa, rhob, tau_a, tau_b
      ! derivatives w.r.t. the reduced gradients gamma_ij = grad rho_i . grad rho_j
      #:for descs, arr, idx in mgga_all
         REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), POINTER :: e_${'_'.join(descs)}$
      #:endfor
      TYPE(xc_derivative_type), POINTER        :: deriv
      TYPE(xc_f03_func_info_t)                 :: xc_info

      CALL timeset(routineN, handle)

      ! Only "everything up to grad_deriv" is supported. Asking for a single
      ! derivative order in isolation is rejected here rather than further down,
      ! where it would quietly evaluate nothing at all.
      IF (grad_deriv > 3) THEN
         CPABORT("derivatives larger than 3 not implemented")
      END IF

      NULLIFY (dummy)
      NULLIFY (workers)
      NULLIFY (rhoa, rhob, norm_drho, norm_drhoa, norm_drhob, laplace_rhoa, &
               laplace_rhob, tau_a, tau_b)

      IF (PRESENT(func_name_override)) THEN
         func_name = func_name_override
         func_scale = 1.0_dp
      ELSE
         func_name = libxc_params%section%name
         CALL section_vals_val_get(libxc_params, "scale", r_val=func_scale)
      END IF

      IF (ABS(func_scale - 1.0_dp) < 1.0e-10_dp) func_scale = 1.0_dp

      CALL xc_rho_set_get(rho_set, can_return_null=.TRUE., &
                          rhoa=rhoa, rhob=rhob, norm_drho=norm_drho, &
                          norm_drhoa=norm_drhoa, norm_drhob=norm_drhob, &
                          laplace_rhoa=laplace_rhoa, laplace_rhob=laplace_rhob, &
                          rho_cutoff=epsilon_rho, tau_cutoff=epsilon_tau, &
                          tau_a=tau_a, tau_b=tau_b, local_bounds=bo)

      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)

      ! One functional object per thread, kept across calls. The cutoffs are part
      ! of the worker setup because LibXC does the screening itself.
      IF (PRESENT(func_name_override)) THEN
         ! an overriding name carries no input section, so no parameters are set
         NULLIFY (no_params)
         iw = libxc_get_workers(func_name, XC_POLARIZED, no_params, epsilon_rho, epsilon_tau, &
                                grad_deriv)
      ELSE
         iw = libxc_get_workers(func_name, XC_POLARIZED, libxc_params, epsilon_rho, epsilon_tau, &
                                grad_deriv)
      END IF
      workers => libxc_cache(iw)
      xc_info = workers%worker(1)%info
      has_laplace = workers%has_laplace

      ! see libxc_unused: the arguments the requested order does not produce are
      ! never touched, but they still have to point somewhere
      CALL libxc_reserve_unused(npoints)
      dummy(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)) => libxc_unused(1:npoints)

      ! due to assumed shape array usage in next routine
      IF (.NOT. ASSOCIATED(norm_drho)) norm_drho => rhoa
      IF (.NOT. ASSOCIATED(norm_drhoa)) norm_drhoa => rhoa
      IF (.NOT. ASSOCIATED(norm_drhob)) norm_drhob => rhoa
      IF (.NOT. ASSOCIATED(tau_a)) tau_a => rhoa
      IF (.NOT. ASSOCIATED(tau_b)) tau_b => rhoa

      ! only some MGGA functionals really need the Laplacian,
      ! all others can work with rhoa (read-only) as dummy
      IF (.NOT. has_laplace) laplace_rhoa => rhoa
      IF (.NOT. has_laplace) laplace_rhob => rhoa

      e_0 => dummy
      e_rhoa => dummy
      e_rhob => dummy
      e_ndrho => dummy
      e_ndrhoa => dummy
      e_ndrhob => dummy
      e_laplace_rhoa => dummy
      e_laplace_rhob => dummy
      e_tau_a => dummy
      e_tau_b => dummy
      e_rhoa_rhoa => dummy
      e_rhoa_rhob => dummy
      e_rhob_rhob => dummy
      e_ndrho_rhoa => dummy
      e_ndrho_rhob => dummy
      e_ndrhoa_rhoa => dummy
      e_ndrhoa_rhob => dummy
      e_ndrhob_rhoa => dummy
      e_ndrhob_rhob => dummy
      e_ndrho_ndrho => dummy
      e_ndrho_ndrhoa => dummy
      e_ndrho_ndrhob => dummy
      e_ndrhoa_ndrhoa => dummy
      e_ndrhoa_ndrhob => dummy
      e_ndrhob_ndrhob => dummy
      e_rhoa_laplace_rhoa => dummy
      e_rhoa_laplace_rhob => dummy
      e_rhob_laplace_rhoa => dummy
      e_rhob_laplace_rhob => dummy
      e_rhoa_tau_a => dummy
      e_rhoa_tau_b => dummy
      e_rhob_tau_a => dummy
      e_rhob_tau_b => dummy
      e_ndrho_laplace_rhoa => dummy
      e_ndrho_laplace_rhob => dummy
      e_ndrhoa_laplace_rhoa => dummy
      e_ndrhoa_laplace_rhob => dummy
      e_ndrhob_laplace_rhoa => dummy
      e_ndrhob_laplace_rhob => dummy
      e_ndrho_tau_a => dummy
      e_ndrho_tau_b => dummy
      e_ndrhoa_tau_a => dummy
      e_ndrhoa_tau_b => dummy
      e_ndrhob_tau_a => dummy
      e_ndrhob_tau_b => dummy
      e_laplace_rhoa_laplace_rhoa => dummy
      e_laplace_rhoa_laplace_rhob => dummy
      e_laplace_rhob_laplace_rhob => dummy
      e_laplace_rhoa_tau_a => dummy
      e_laplace_rhoa_tau_b => dummy
      e_laplace_rhob_tau_a => dummy
      e_laplace_rhob_tau_b => dummy
      e_tau_a_tau_a => dummy
      e_tau_a_tau_b => dummy
      e_tau_b_tau_b => dummy
      e_rhoa_rhoa_rhoa => dummy
      e_rhoa_rhoa_rhob => dummy
      e_rhoa_rhob_rhob => dummy
      e_rhob_rhob_rhob => dummy
      #:for descs, arr, idx in mgga_all
         e_${'_'.join(descs)}$ => dummy
      #:endfor

      deriv => xc_dset_get_derivative(deriv_set, [INTEGER::], &
                                      allocate_deriv=.TRUE.)
      CALL xc_derivative_get(deriv, deriv_data=e_0)
      IF (grad_deriv >= 1) THEN
         SELECT CASE (xc_f03_func_info_get_family(xc_info))
         CASE (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhob)
         CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhob)
            #:for descs, arr, idx in gamma_only_1
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
         CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhob)
            #:for descs, arr, idx in gamma_only_1
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
            deriv => xc_dset_get_derivative(deriv_set, [deriv_tau_a], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_tau_a)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_tau_b], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_tau_b)
            IF (has_laplace) THEN
               deriv => xc_dset_get_derivative(deriv_set, [deriv_laplace_rhoa], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_laplace_rhoa)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_laplace_rhob], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_laplace_rhob)
            END IF
         CASE default
            CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
         END SELECT
      END IF
      IF (grad_deriv >= 2) THEN
         SELECT CASE (xc_f03_func_info_get_family(xc_info))
         CASE (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhob_rhob)
         CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhob_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa, deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob, deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhob_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhob_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_norm_drho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_ndrho)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_norm_drhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_ndrhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_norm_drhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_ndrhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa, deriv_norm_drhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa_ndrhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa, deriv_norm_drhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa_ndrhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob, deriv_norm_drhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhob_ndrhob)
            #:for descs, arr, idx in gamma_only_2
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
         CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)

            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhob_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa, deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob, deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhob_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhob_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_norm_drho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_ndrho)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_norm_drhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_ndrhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_norm_drhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_ndrhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa, deriv_norm_drhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa_ndrhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa, deriv_norm_drhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa_ndrhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob, deriv_norm_drhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhob_ndrhob)
            #:for descs, arr, idx in mgga_fill_2_nl
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
            IF (has_laplace) THEN
               #:for descs, arr, idx in mgga_fill_2_l
                  deriv => xc_dset_get_derivative(deriv_set, &
                                                  [${', '.join('deriv_' + d for d in descs)}$], &
                                                  allocate_deriv=.TRUE.)
                  CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
               #:endfor
            END IF
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_tau_a], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa_tau_a)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_tau_b], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa_tau_b)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob, deriv_tau_a], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhob_tau_a)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob, deriv_tau_b], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhob_tau_b)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_tau_a], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_tau_a)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_tau_b], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_tau_b)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa, deriv_tau_a], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa_tau_a)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa, deriv_tau_b], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa_tau_b)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob, deriv_tau_a], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhob_tau_a)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob, deriv_tau_b], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrhob_tau_b)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_tau_a, deriv_tau_a], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_tau_a_tau_a)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_tau_a, deriv_tau_b], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_tau_a_tau_b)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_tau_b, deriv_tau_b], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_tau_b_tau_b)
            IF (has_laplace) THEN
               deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_laplace_rhoa], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_rhoa_laplace_rhoa)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_laplace_rhob], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_rhoa_laplace_rhob)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob, deriv_laplace_rhoa], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_rhob_laplace_rhoa)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob, deriv_laplace_rhob], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_rhob_laplace_rhob)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_laplace_rhoa], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_ndrho_laplace_rhoa)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_laplace_rhob], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_ndrho_laplace_rhob)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa, deriv_laplace_rhoa], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa_laplace_rhoa)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhoa, deriv_laplace_rhob], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_ndrhoa_laplace_rhob)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob, deriv_laplace_rhoa], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_ndrhob_laplace_rhoa)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drhob, deriv_laplace_rhob], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_ndrhob_laplace_rhob)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_laplace_rhoa, deriv_laplace_rhoa], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_laplace_rhoa_laplace_rhoa)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_laplace_rhoa, deriv_laplace_rhob], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_laplace_rhoa_laplace_rhob)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_laplace_rhob, deriv_laplace_rhob], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_laplace_rhob_laplace_rhob)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_laplace_rhoa, deriv_tau_a], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_laplace_rhoa_tau_a)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_laplace_rhoa, deriv_tau_b], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_laplace_rhoa_tau_b)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_laplace_rhob, deriv_tau_a], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_laplace_rhob_tau_a)
               deriv => xc_dset_get_derivative(deriv_set, [deriv_laplace_rhob, deriv_tau_b], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_laplace_rhob_tau_b)
            END IF
         CASE default
            CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
         END SELECT
      END IF
      IF (grad_deriv >= 3) THEN
         SELECT CASE (xc_f03_func_info_get_family(xc_info))
         CASE (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_rhoa, deriv_rhoa], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa_rhoa_rhoa)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_rhoa, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa_rhoa_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhoa, deriv_rhob, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhoa_rhob_rhob)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rhob, deriv_rhob, deriv_rhob], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rhob_rhob_rhob)
         CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            #:for descs, arr, idx in gamma_only_3 + rho3_entries
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
         CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            #:for descs, arr, idx in mgga_fill_3_nl
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
            IF (has_laplace) THEN
               #:for descs, arr, idx in mgga_fill_3_l
                  deriv => xc_dset_get_derivative(deriv_set, &
                                                  [${', '.join('deriv_' + d for d in descs)}$], &
                                                  allocate_deriv=.TRUE.)
                  CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
               #:endfor
            END IF
            #:for descs, arr, idx in rho3_entries
               deriv => xc_dset_get_derivative(deriv_set, &
                                               [${', '.join('deriv_' + d for d in descs)}$], &
                                               allocate_deriv=.TRUE.)
               CALL xc_derivative_get(deriv, deriv_data=e_${'_'.join(descs)}$)
            #:endfor
         CASE default
            CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
         END SELECT
      END IF

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP SHARED(rhoa,rhob,norm_drho,norm_drhoa,norm_drhob),&
!$OMP SHARED(laplace_rhoa,laplace_rhob,tau_a,tau_b),&
!$OMP SHARED(e_0,e_rhoa,e_rhob,e_ndrho,e_ndrhoa,e_ndrhob),&
!$OMP SHARED(e_laplace_rhoa,e_laplace_rhob,e_tau_a,e_tau_b),&
!$OMP SHARED(e_rhoa_rhoa,e_rhoa_rhob,e_rhob_rhob),&
!$OMP SHARED(e_ndrho_rhoa,e_ndrho_rhob),&
!$OMP SHARED(e_ndrhoa_rhoa,e_ndrhoa_rhob,e_ndrhob_rhoa,e_ndrhob_rhob),&
!$OMP SHARED(e_ndrho_ndrho,e_ndrho_ndrhoa,e_ndrho_ndrhob),&
!$OMP SHARED(e_ndrhoa_ndrhoa,e_ndrhoa_ndrhob,e_ndrhob_ndrhob),&
!$OMP SHARED(e_rhoa_laplace_rhoa,e_rhoa_laplace_rhob,e_rhob_laplace_rhoa,e_rhob_laplace_rhob),&
!$OMP SHARED(e_rhoa_tau_a,e_rhoa_tau_b,e_rhob_tau_a,e_rhob_tau_b),&
!$OMP SHARED(e_ndrho_laplace_rhoa,e_ndrho_laplace_rhob),&
!$OMP SHARED(e_ndrhoa_laplace_rhoa,e_ndrhoa_laplace_rhob,e_ndrhob_laplace_rhoa,e_ndrhob_laplace_rhob),&
!$OMP SHARED(e_ndrho_tau_a,e_ndrho_tau_b),&
!$OMP SHARED(e_ndrhoa_tau_a,e_ndrhoa_tau_b,e_ndrhob_tau_a,e_ndrhob_tau_b),&
!$OMP SHARED(e_laplace_rhoa_laplace_rhoa,e_laplace_rhoa_laplace_rhob,e_laplace_rhob_laplace_rhob),&
!$OMP SHARED(e_laplace_rhoa_tau_a,e_laplace_rhoa_tau_b,e_laplace_rhob_tau_a,e_laplace_rhob_tau_b),&
!$OMP SHARED(e_tau_a_tau_a,e_tau_a_tau_b,e_tau_b_tau_b),&
!$OMP SHARED(e_rhoa_rhoa_rhoa,e_rhoa_rhoa_rhob,e_rhoa_rhob_rhob,e_rhob_rhob_rhob),&
      #:for descs, arr, idx in mgga_all
!$OMP SHARED(e_${'_'.join(descs)}$),&
      #:endfor
!$OMP SHARED(grad_deriv,npoints),&
!$OMP SHARED(func_name,func_scale,workers)

      CALL libxc_spin_polarized_calc(rhoa=rhoa, rhob=rhob, norm_drho=norm_drho, &
                                     norm_drhoa=norm_drhoa, norm_drhob=norm_drhob, laplace_rhoa=laplace_rhoa, &
                                     laplace_rhob=laplace_rhob, tau_a=tau_a, tau_b=tau_b, &
                                     e_0=e_0, e_rhoa=e_rhoa, e_rhob=e_rhob, e_ndrho=e_ndrho, &
                                     e_ndrhoa=e_ndrhoa, e_ndrhob=e_ndrhob, e_laplace_rhoa=e_laplace_rhoa, &
                                     e_laplace_rhob=e_laplace_rhob, e_tau_a=e_tau_a, e_tau_b=e_tau_b, &
                                     e_rhoa_rhoa=e_rhoa_rhoa, e_rhoa_rhob=e_rhoa_rhob, e_rhob_rhob=e_rhob_rhob, &
                                     e_ndrho_rhoa=e_ndrho_rhoa, e_ndrho_rhob=e_ndrho_rhob, &
                                     e_ndrhoa_rhoa=e_ndrhoa_rhoa, e_ndrhoa_rhob=e_ndrhoa_rhob, &
                                     e_ndrhob_rhoa=e_ndrhob_rhoa, e_ndrhob_rhob=e_ndrhob_rhob, &
                                     e_ndrho_ndrho=e_ndrho_ndrho, e_ndrho_ndrhoa=e_ndrho_ndrhoa, &
                                     e_ndrho_ndrhob=e_ndrho_ndrhob, e_ndrhoa_ndrhoa=e_ndrhoa_ndrhoa, &
                                     e_ndrhoa_ndrhob=e_ndrhoa_ndrhob, e_ndrhob_ndrhob=e_ndrhob_ndrhob, &
                                     e_rhoa_laplace_rhoa=e_rhoa_laplace_rhoa, &
                                     e_rhoa_laplace_rhob=e_rhoa_laplace_rhob, &
                                     e_rhob_laplace_rhoa=e_rhob_laplace_rhoa, &
                                     e_rhob_laplace_rhob=e_rhob_laplace_rhob, &
                                     e_rhoa_tau_a=e_rhoa_tau_a, e_rhoa_tau_b=e_rhoa_tau_b, &
                                     e_rhob_tau_a=e_rhob_tau_a, e_rhob_tau_b=e_rhob_tau_b, &
                                     e_ndrho_laplace_rhoa=e_ndrho_laplace_rhoa, &
                                     e_ndrho_laplace_rhob=e_ndrho_laplace_rhob, &
                                     e_ndrhoa_laplace_rhoa=e_ndrhoa_laplace_rhoa, &
                                     e_ndrhoa_laplace_rhob=e_ndrhoa_laplace_rhob, &
                                     e_ndrhob_laplace_rhoa=e_ndrhob_laplace_rhoa, &
                                     e_ndrhob_laplace_rhob=e_ndrhob_laplace_rhob, &
                                     e_ndrho_tau_a=e_ndrho_tau_a, e_ndrho_tau_b=e_ndrho_tau_b, &
                                     e_ndrhoa_tau_a=e_ndrhoa_tau_a, e_ndrhoa_tau_b=e_ndrhoa_tau_b, &
                                     e_ndrhob_tau_a=e_ndrhob_tau_a, e_ndrhob_tau_b=e_ndrhob_tau_b, &
                                     e_laplace_rhoa_laplace_rhoa=e_laplace_rhoa_laplace_rhoa, &
                                     e_laplace_rhoa_laplace_rhob=e_laplace_rhoa_laplace_rhob, &
                                     e_laplace_rhob_laplace_rhob=e_laplace_rhob_laplace_rhob, &
                                     e_laplace_rhoa_tau_a=e_laplace_rhoa_tau_a, &
                                     e_laplace_rhoa_tau_b=e_laplace_rhoa_tau_b, &
                                     e_laplace_rhob_tau_a=e_laplace_rhob_tau_a, &
                                     e_laplace_rhob_tau_b=e_laplace_rhob_tau_b, &
                                     e_tau_a_tau_a=e_tau_a_tau_a, &
                                     e_tau_a_tau_b=e_tau_a_tau_b, &
                                     e_tau_b_tau_b=e_tau_b_tau_b, &
                                     e_rhoa_rhoa_rhoa=e_rhoa_rhoa_rhoa, &
                                     e_rhoa_rhoa_rhob=e_rhoa_rhoa_rhob, &
                                     e_rhoa_rhob_rhob=e_rhoa_rhob_rhob, &
                                     e_rhob_rhob_rhob=e_rhob_rhob_rhob, &
                                      #:for descs, arr, idx in mgga_all
                                     e_${'_'.join(descs)}$=e_${'_'.join(descs)}$, &
                                      #:endfor
                                     grad_deriv=grad_deriv, npoints=npoints, &
                                     func_name=func_name, sc=func_scale, workers=workers)

!$OMP END PARALLEL

      NULLIFY (dummy)
      NULLIFY (workers)

      CALL timestop(handle)
#else
      MARK_USED(rho_set)
      MARK_USED(deriv_set)
      MARK_USED(grad_deriv)
      MARK_USED(libxc_params)
      MARK_USED(func_name_override)

      CALL cp_abort(__LOCATION__, "Unknown functional! If you are asking "// &
                    "for a functional of the LibXC library, "// &
                    "you have to download and install the library!")
#endif
   END SUBROUTINE libxc_spin_polarized_eval

#if defined (__LIBXC)

! **************************************************************************************************
!> \brief Compact C2 maximum, with a transition width relative to its positive lower argument.
!> \param a First argument.
!> \param b Positive second argument.
!> \param value Projected maximum.
!> \param first First derivatives in a and b.
!> \param second Second derivatives in a and b.
!> \param active Optional indication that the first argument is not on the identity branch.
! **************************************************************************************************
   PURE SUBROUTINE libxc_c2_maximum(a, b, value, first, second, active)
      REAL(dp), INTENT(IN) :: a, b
      REAL(dp), INTENT(OUT) :: value, first(2), second(2, 2)
      LOGICAL, INTENT(OUT), OPTIONAL :: active
      REAL(dp), PARAMETER :: relative_width = 64.0_dp*SQRT(EPSILON(0.0_dp))
      REAL(dp) :: width, u, polynomial, slope, curvature, ratio

      width = relative_width*b
      IF (PRESENT(active)) active = a - b < width
      second = 0.0_dp
      IF (a - b >= width) THEN
         value = a
         first = [1.0_dp, 0.0_dp]
      ELSE IF (a - b <= -width) THEN
         value = b
         first = [0.0_dp, 1.0_dp]
      ELSE
         u = (a - b)/width
         polynomial = (1.0_dp + u)**3*(3.0_dp - u)/16.0_dp
         slope = (1.0_dp + u)**2*(2.0_dp - u)/4.0_dp
         curvature = 0.75_dp*(1.0_dp - u*u)/width
         ratio = a/b
         value = b + width*polynomial
         first = [slope, 1.0_dp + relative_width*polynomial - ratio*slope]
         second(1, 1) = curvature
         second(1, 2) = -ratio*curvature
         second(2, 1) = second(1, 2)
         second(2, 2) = ratio*ratio*curvature
      END IF
   END SUBROUTINE libxc_c2_maximum

! **************************************************************************************************
!> \brief C2 kinetic projection and its full derivatives in rho, sigma and raw tau.
!> \param rho Density used in the bounded denominator.
!> \param sigma Same-spin gradient contraction.
!> \param raw_tau Unprojected kinetic density.
!> \param epsilon_rho Density cutoff.
!> \param epsilon_tau Positive kinetic floor.
!> \param first Optional first derivatives.
!> \param second Optional second derivatives.
!> \param active Optional projection activity, independent of rounding of the value.
!> \return Projected kinetic density.
! **************************************************************************************************
   FUNCTION libxc_c2_project_tau(rho, sigma, raw_tau, epsilon_rho, epsilon_tau, first, second, active) RESULT(value)
      REAL(dp), INTENT(IN) :: rho, sigma, raw_tau, epsilon_rho, epsilon_tau
      REAL(dp), INTENT(OUT), OPTIONAL :: first(3), second(3, 3)
      LOGICAL, INTENT(OUT), OPTIONAL :: active
      LOGICAL :: projection_active
      REAL(dp) :: value, denominator, tau_w, lower, lower_first(2), lower_second(2, 2), &
                  final_first(2), final_second(2, 2), w_first(2), w_second(2, 2), outer(2, 2)
      INTEGER :: i, j

      denominator = MAX(rho, epsilon_rho, TINY(0.0_dp))
      tau_w = sigma/(8.0_dp*denominator)
      CALL libxc_c2_maximum(tau_w, MAX(epsilon_tau, EPSILON(0.0_dp)), lower, lower_first, lower_second)
      CALL libxc_c2_maximum(raw_tau, lower, value, final_first, final_second, projection_active)
      IF (PRESENT(active)) active = projection_active
      IF (.NOT. projection_active) THEN
         IF (PRESENT(first)) first = [0.0_dp, 0.0_dp, 1.0_dp]
         IF (PRESENT(second)) second = 0.0_dp
         RETURN
      END IF
      IF (.NOT. PRESENT(first) .AND. .NOT. PRESENT(second)) RETURN
      w_first = [0.0_dp, 1.0_dp/(8.0_dp*denominator)]
      w_second = 0.0_dp
      IF (rho > epsilon_rho) THEN
         w_first(1) = -tau_w/denominator
         w_second(1, 1) = 2.0_dp*tau_w/denominator**2
         w_second(1, 2) = -w_first(2)/denominator
         w_second(2, 1) = w_second(1, 2)
      END IF
      IF (PRESENT(first)) THEN
         first(1:2) = final_first(2)*lower_first(1)*w_first
         first(3) = final_first(1)
      END IF
      IF (PRESENT(second)) THEN
         DO j = 1, 2
            DO i = 1, 2
               outer(i, j) = w_first(i)*w_first(j)
            END DO
         END DO
         second = 0.0_dp
         second(1:2, 1:2) = final_first(2)*(lower_first(1)*w_second + lower_second(1, 1)*outer) + &
                            final_second(2, 2)*lower_first(1)**2*outer
         second(1:2, 3) = final_second(1, 2)*lower_first(1)*w_first
         second(3, 1:2) = second(1:2, 3)
         second(3, 3) = final_second(1, 1)
      END IF
   END FUNCTION libxc_c2_project_tau

! **************************************************************************************************
!> \brief Pull back Libxc meta-GGA derivatives through the kinetic-density projection.
!> \param w Functional worker with derivatives in Libxc variables.
!> \param point Point within the worker's current block.
!> \param raw_tau Kinetic densities before projection.
!> \param nspin Number of spin channels.
!> \param grad_deriv Requested derivative order.
!> \param has_laplace Whether the functional uses the density Laplacian.
!> \param epsilon_rho Density cutoff used in the projection denominator.
!> \param epsilon_tau Kinetic-density cutoff used as a positive floor.
! **************************************************************************************************
   SUBROUTINE libxc_mgga_tau_pullback(w, point, raw_tau, nspin, grad_deriv, has_laplace, &
                                      epsilon_rho, epsilon_tau)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: w
      INTEGER, INTENT(IN)                                :: point, nspin, grad_deriv
      REAL(dp), INTENT(IN)                               :: raw_tau(nspin), epsilon_rho, epsilon_tau
      LOGICAL, INTENT(IN)                                :: has_laplace

      INTEGER                                            :: channel, i, irho, isigma, itau, nvar
      REAL(dp)                                           :: projected, tau_first(3), tau_second(3, 3)
      INTEGER                                            :: indices(3), j
      LOGICAL                                            :: active(2)
      REAL(dp)                                           :: gradient(9), hessian(9, 9), jacobian(9, 9), &
                                                            curvature(9, 9, 2), pulled_gradient(9), &
                                                            pulled_hessian(9, 9)

      IF (grad_deriv < 1) RETURN
      CPASSERT(nspin == 1 .OR. nspin == 2)
      nvar = 4
      IF (nspin == 2) nvar = 9
      gradient = 0.0_dp
      hessian = 0.0_dp
      jacobian = 0.0_dp
      curvature = 0.0_dp
      DO i = 1, nvar
         jacobian(i, i) = 1.0_dp
      END DO

      ! The transition can change curvature even when its value rounds to raw tau.
      active = .FALSE.
      DO channel = 1, nspin
         irho = channel
         isigma = 2
         itau = 4
         IF (nspin == 2) THEN
            isigma = 2*channel + 1
            itau = 7 + channel
         END IF
         projected = libxc_c2_project_tau(w%rho(channel, point), w%sigma(isigma - nspin, point), &
                                          raw_tau(channel), epsilon_rho, epsilon_tau, tau_first, tau_second, active(channel))
         IF (.NOT. active(channel)) CYCLE
         jacobian(itau, :) = 0.0_dp
         indices = [irho, isigma, itau]
         DO i = 1, 3
            jacobian(itau, indices(i)) = tau_first(i)
            DO j = 1, 3
               curvature(indices(i), indices(j), channel) = tau_second(i, j)
            END DO
         END DO
      END DO
      IF (.NOT. ANY(active(1:nspin))) RETURN
      IF (grad_deriv > 2) CPABORT("Libxc meta-GGA tau projection supports derivatives through second order")

      IF (nspin == 1) THEN
         #:for descs, arr, idx in umgga_deriv_1
            #:if has_lap((descs, arr, idx))
            IF (has_laplace) THEN
               gradient(${umgga_vars.index(descs[0]) + 1}$) = w%${arr}$ (${idx}$, point)
            END IF
            #:else
            gradient(${umgga_vars.index(descs[0]) + 1}$) = w%${arr}$ (${idx}$, point)
            #:endif
         #:endfor
         IF (grad_deriv >= 2) THEN
            #:for descs, arr, idx in umgga_deriv_2
               #:set a, b = [umgga_vars.index(d) + 1 for d in descs]
               #:if has_lap((descs, arr, idx))
               IF (has_laplace) THEN
                  hessian(${a}$, ${b}$) = w%${arr}$ (${idx}$, point)
               END IF
               #:else
               hessian(${a}$, ${b}$) = w%${arr}$ (${idx}$, point)
               #:endif
               hessian(${b}$, ${a}$) = hessian(${a}$, ${b}$)
            #:endfor
         END IF
      ELSE
         #:for descs, arr, idx in mgga_deriv_1
            #:if has_lap((descs, arr, idx))
            IF (has_laplace) THEN
               gradient(${mgga_vars.index(descs[0]) + 1}$) = w%${arr}$ (${idx}$, point)
            END IF
            #:else
            gradient(${mgga_vars.index(descs[0]) + 1}$) = w%${arr}$ (${idx}$, point)
            #:endif
         #:endfor
         IF (grad_deriv >= 2) THEN
            #:for descs, arr, idx in mgga_deriv_2
               #:set a, b = [mgga_vars.index(d) + 1 for d in descs]
               #:if has_lap((descs, arr, idx))
               IF (has_laplace) THEN
                  hessian(${a}$, ${b}$) = w%${arr}$ (${idx}$, point)
               END IF
               #:else
               hessian(${a}$, ${b}$) = w%${arr}$ (${idx}$, point)
               #:endif
               hessian(${b}$, ${a}$) = hessian(${a}$, ${b}$)
            #:endfor
         END IF
      END IF

      pulled_gradient(1:nvar) = MATMUL(TRANSPOSE(jacobian(1:nvar, 1:nvar)), gradient(1:nvar))
      IF (grad_deriv >= 2) THEN
         pulled_hessian(1:nvar, 1:nvar) = &
            MATMUL(TRANSPOSE(jacobian(1:nvar, 1:nvar)), &
                   MATMUL(hessian(1:nvar, 1:nvar), jacobian(1:nvar, 1:nvar)))
         DO channel = 1, nspin
            itau = 4
            IF (nspin == 2) itau = 7 + channel
            pulled_hessian(1:nvar, 1:nvar) = pulled_hessian(1:nvar, 1:nvar) + &
                                             gradient(itau)*curvature(1:nvar, 1:nvar, channel)
         END DO
      END IF
      IF (nspin == 1) THEN
         #:for descs, arr, idx in umgga_deriv_1
            #:if has_lap((descs, arr, idx))
            IF (has_laplace) THEN
               w%${arr}$ (${idx}$, point) = pulled_gradient(${umgga_vars.index(descs[0]) + 1}$)
            END IF
            #:else
            w%${arr}$ (${idx}$, point) = pulled_gradient(${umgga_vars.index(descs[0]) + 1}$)
            #:endif
         #:endfor
         IF (grad_deriv >= 2) THEN
            #:for descs, arr, idx in umgga_deriv_2
               #:set a, b = [umgga_vars.index(d) + 1 for d in descs]
               #:if has_lap((descs, arr, idx))
               IF (has_laplace) THEN
                  w%${arr}$ (${idx}$, point) = pulled_hessian(${a}$, ${b}$)
               END IF
               #:else
               w%${arr}$ (${idx}$, point) = pulled_hessian(${a}$, ${b}$)
               #:endif
            #:endfor
         END IF
      ELSE
         #:for descs, arr, idx in mgga_deriv_1
            #:if has_lap((descs, arr, idx))
            IF (has_laplace) THEN
               w%${arr}$ (${idx}$, point) = pulled_gradient(${mgga_vars.index(descs[0]) + 1}$)
            END IF
            #:else
            w%${arr}$ (${idx}$, point) = pulled_gradient(${mgga_vars.index(descs[0]) + 1}$)
            #:endif
         #:endfor
         IF (grad_deriv >= 2) THEN
            #:for descs, arr, idx in mgga_deriv_2
               #:set a, b = [mgga_vars.index(d) + 1 for d in descs]
               #:if has_lap((descs, arr, idx))
               IF (has_laplace) THEN
                  w%${arr}$ (${idx}$, point) = pulled_hessian(${a}$, ${b}$)
               END IF
               #:else
               w%${arr}$ (${idx}$, point) = pulled_hessian(${a}$, ${b}$)
               #:endif
            #:endfor
         END IF
      END IF
   END SUBROUTINE libxc_mgga_tau_pullback

! **************************************************************************************************
!> \brief Differentiate exc times the original density through Libxc's individual spin floors.
!> \param w Evaluated worker buffers before the kinetic-density and negative-density pullbacks.
!> \param point Grid point in the worker block.
!> \param grad_deriv Highest evaluated derivative order.
!> \param has_laplace Whether the functional uses the density Laplacian.
!> \param epsilon_rho Libxc density cutoff.
!> \note Libxc returns derivatives of F at m_s=MAX(rho_s,epsilon_rho), whereas CP2K
!>       accumulates E=n*exc=(n/m)*F with n=SUM(rho_s) and m=SUM(m_s). This pullback
!>       leaves that energy unchanged, including its original-density normalization.
! **************************************************************************************************
   SUBROUTINE libxc_spin_floor_pullback(w, point, grad_deriv, has_laplace, epsilon_rho)
      TYPE(libxc_worker_type), INTENT(INOUT) :: w
      INTEGER, INTENT(IN) :: point, grad_deriv
      LOGICAL, INTENT(IN) :: has_laplace
      REAL(dp), INTENT(IN) :: epsilon_rho

      INTEGER :: i, j
      REAL(dp) :: n, m, ratio, functional, third_ratio
      REAL(dp) :: diagonal(9), dn(9), dm(9), first_ratio(9), second_ratio(9, 9), &
                  gradient(9), hessian(9, 9)

      IF (grad_deriv < 1 .OR. ALL(w%rho(:, point) >= epsilon_rho)) RETURN
      n = SUM(w%rho(:, point))
      IF (n < epsilon_rho) RETURN ! Libxc screened this point to zero.
      m = SUM(MAX(w%rho(:, point), epsilon_rho))
      ratio = n/m
      functional = m*w%exc(point)
      diagonal = 1.0_dp
      diagonal(1:2) = MERGE(1.0_dp, 0.0_dp, w%rho(:, point) >= epsilon_rho)
      dn = 0.0_dp
      dn(1:2) = 1.0_dp
      dm = 0.0_dp
      dm(1:2) = diagonal(1:2)
      first_ratio = (dn - ratio*dm)/m
      second_ratio = 0.0_dp
      DO i = 1, 2
         DO j = 1, 2
            second_ratio(i, j) = (-dn(i)*dm(j) - dn(j)*dm(i) + 2.0_dp*ratio*dm(i)*dm(j))/m**2
         END DO
      END DO
      gradient = 0.0_dp
      hessian = 0.0_dp
      #:for descs, arr, idx in mgga_deriv_1
         #:set a = mgga_vars.index(descs[0]) + 1
         IF (ALLOCATED(w%${arr}$) ${' .AND. has_laplace' if has_lap((descs,arr,idx)) else ''}$) THEN
            gradient(${a}$) = w%${arr}$ (${idx}$, point)
         END IF
      #:endfor
      IF (grad_deriv >= 2) THEN
         #:for descs, arr, idx in mgga_deriv_2
            #:set a,b = [mgga_vars.index(d)+1 for d in descs]
            IF (ALLOCATED(w%${arr}$) ${' .AND. has_laplace' if has_lap((descs,arr,idx)) else ''}$) THEN
               hessian(${a}$, ${b}$) = w%${arr}$ (${idx}$, point)
            END IF
            hessian(${b}$, ${a}$) = hessian(${a}$, ${b}$)
         #:endfor
      END IF
      ! Use unmodified first and second derivatives in every product-rule term.
      IF (grad_deriv >= 3) THEN
         #:for descs, arr, idx in mgga_deriv_3
            #:set a,b,c = [mgga_vars.index(d)+1 for d in descs]
            IF (ALLOCATED(w%${arr}$) ${' .AND. has_laplace' if has_lap((descs,arr,idx)) else ''}$) THEN
               third_ratio = (2.0_dp*(dn(${a}$)*dm(${b}$)*dm(${c}$) + &
                                      dn(${b}$)*dm(${a}$)*dm(${c}$) + dn(${c}$)*dm(${a}$)*dm(${b}$)) - &
                              6.0_dp*ratio*dm(${a}$)*dm(${b}$)*dm(${c}$))/m**3
               w%${arr}$ (${idx}$, point) = ratio*diagonal(${a}$)*diagonal(${b}$)*diagonal(${c}$)* &
                           w%${arr}$ (${idx}$, point) + first_ratio(${a}$)*diagonal(${b}$)*diagonal(${c}$)*hessian(${b}$, ${c}$) + &
                                            first_ratio(${b}$)*diagonal(${a}$)*diagonal(${c}$)*hessian(${a}$, ${c}$) + &
                                            first_ratio(${c}$)*diagonal(${a}$)*diagonal(${b}$)*hessian(${a}$, ${b}$) + &
                                            second_ratio(${a}$, ${b}$)*diagonal(${c}$)*gradient(${c}$) + &
                                            second_ratio(${a}$, ${c}$)*diagonal(${b}$)*gradient(${b}$) + &
                                            second_ratio(${b}$, ${c}$)*diagonal(${a}$)*gradient(${a}$) + third_ratio*functional
            END IF
         #:endfor
      END IF
      IF (grad_deriv >= 2) THEN
         #:for descs, arr, idx in mgga_deriv_2
            #:set a,b = [mgga_vars.index(d)+1 for d in descs]
            IF (ALLOCATED(w%${arr}$) ${' .AND. has_laplace' if has_lap((descs,arr,idx)) else ''}$) w%${arr}$ (${idx}$, point) = &
               ratio*diagonal(${a}$)*diagonal(${b}$)*hessian(${a}$, ${b}$) + &
               first_ratio(${a}$)*diagonal(${b}$)*gradient(${b}$) + &
               first_ratio(${b}$)*diagonal(${a}$)*gradient(${a}$) + second_ratio(${a}$, ${b}$)*functional

         #:endfor
      END IF
      #:for descs, arr, idx in mgga_deriv_1
         #:set a = mgga_vars.index(descs[0])+1
         IF (ALLOCATED(w%${arr}$) ${' .AND. has_laplace' if has_lap((descs,arr,idx)) else ''}$) w%${arr}$ (${idx}$, point) = &
            ratio*diagonal(${a}$)*gradient(${a}$) + first_ratio(${a}$)*functional

      #:endfor
   END SUBROUTINE libxc_spin_floor_pullback

! **************************************************************************************************
!> \brief Pull Libxc derivatives back through the negative spin-density clamp.
!> \param w Evaluated worker buffers, updated before functional contributions are accumulated.
!> \param point Grid point in the worker block.
!> \param raw_rho Original alpha and beta densities before MAX(rho,0).
!> \param grad_deriv Highest evaluated derivative order.
! **************************************************************************************************
   SUBROUTINE libxc_spin_density_pullback(w, point, raw_rho, grad_deriv)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: w
      INTEGER, INTENT(IN)                                :: point, grad_deriv
      REAL(dp), INTENT(IN)                               :: raw_rho(2)

      INTEGER                                            :: channel

      IF (grad_deriv < 1 .OR. ALL(raw_rho >= 0.0_dp)) RETURN
      DO channel = 1, 2
         IF (raw_rho(channel) >= 0.0_dp) CYCLE
         #:for order, entries in [(1, mgga_deriv_1), (2, mgga_deriv_2), (3, mgga_deriv_3)]
         IF (grad_deriv >= ${order}$) THEN
            #:for descs, arr, idx in entries
               #:if 'rhoa' in descs or 'rhob' in descs
                  #:set channels = [i + 1 for i, name in enumerate(['rhoa', 'rhob']) if name in descs]
                  IF (${ ' .OR. '.join('channel == ' + str(i) for i in channels) }$) THEN
                     IF (ALLOCATED(w%${arr}$)) w%${arr}$ (${idx}$, point) = 0.0_dp
                  END IF
               #:endif
            #:endfor
         END IF
         #:endfor
      END DO
   END SUBROUTINE libxc_spin_density_pullback
#endif

! **************************************************************************************************
!> \brief libxc exchange-correlation functionals
!> \param rho density
!> \param norm_drho norm of the gradient of the density
!> \param laplace_rho laplacian of the density
!> \param tau kinetic-energy density
!> \param e_0 energy density
!> \param e_rho derivative of the energy density with respect to rho
!> \param e_ndrho derivative of the energy density with respect to ndrho
!> \param e_laplace_rho derivative of the energy density with respect to laplace_rho
!> \param e_tau derivative of the energy density with respect to tau
!> \param e_rho_rho derivative of the energy density with respect to rho_rho
!> \param e_ndrho_rho derivative of the energy density with respect to ndrho_rho
!> \param e_ndrho_ndrho derivative of the energy density with respect to ndrho_ndrho
!> \param e_rho_laplace_rho derivative of the energy density with respect to rho_laplace_rho
!> \param e_rho_tau derivative of the energy density with respect to rho_tau
!> \param e_ndrho_laplace_rho derivative of the energy density with respect to ndrho_laplace_rho
!> \param e_ndrho_tau derivative of the energy density with respect to ndrho_tau
!> \param e_laplace_rho_laplace_rho derivative of the energy density with respect to laplace_rho_laplace_rho
!> \param e_laplace_rho_tau derivative of the energy density with respect to laplace_rho_tau
!> \param e_tau_tau derivative of the energy density with respect to tau_tau
!> \param e_rho_rho_rho derivative of the energy density with respect to rho_rho_rho
!> \param grad_deriv maximum degree of the derivative that should be evaluated;
!>        all derivatives up to the given degree are evaluated, in a single
!>        LibXC call per block of grid points
!> \param npoints number of points on the grid
!> \param func_name name of the functional
!> \param sc scaling factor of the functional
!> \param workers cached LibXC functional objects and staging buffers, one per thread
!> \author F. Tran
! **************************************************************************************************
#if defined (__LIBXC)
   SUBROUTINE libxc_spin_unpolarized_calc(rho, norm_drho, laplace_rho, tau, &
                                          e_0, e_rho, e_ndrho, e_laplace_rho, e_tau, e_rho_rho, e_ndrho_rho, &
                                          e_ndrho_ndrho, e_rho_laplace_rho, e_rho_tau, e_ndrho_laplace_rho, &
                                          e_ndrho_tau, e_laplace_rho_laplace_rho, e_laplace_rho_tau, &
                                          e_tau_tau, e_rho_rho_rho, &
#:for descs, arr, idx in umgga_all
                                          e_${'_'.join(descs)}$, &
#:endfor
                                          grad_deriv, npoints, func_name, sc, workers)

   REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, norm_drho, laplace_rho, tau
   REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_0, e_rho, e_ndrho, e_laplace_rho, e_tau, &
                                         e_rho_rho, e_ndrho_rho, e_ndrho_ndrho, e_rho_laplace_rho, e_rho_tau, e_ndrho_laplace_rho, &
                                                 e_ndrho_tau, e_laplace_rho_laplace_rho, e_laplace_rho_tau, e_tau_tau, e_rho_rho_rho
   ! Derivatives with respect to the reduced gradient gamma = |grad rho|^2.
   ! These are LibXC's own sigma derivatives, copied through unchanged: unlike
   ! the norm_drho ones they need no chain rule and carry no 1/|grad rho|.
   #:for descs, arr, idx in umgga_all
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_${'_'.join(descs)}$
   #:endfor
   INTEGER, INTENT(in)                                :: grad_deriv, npoints
   CHARACTER(LEN=default_string_length), INTENT(IN)   :: func_name
   REAL(KIND=dp), INTENT(in)                          :: sc
   TYPE(libxc_worker_set_type), INTENT(INOUT)         :: workers

   INTEGER                                            :: bsize, family, i, i0, ib, ii, ithread, &
                                                         nb, nblocks, nthreads
   LOGICAL                                            :: is_gga, is_mgga

   ithread = 0
   nthreads = 1
!$ ithread = omp_get_thread_num()
!$ nthreads = omp_get_num_threads()
   CPASSERT(ithread < SIZE(workers%worker))

   ! On the GPU path the whole grid is handed to LibXC in a single call (no
   ! block cap), so the per-thread staging and device buffers are grown to
   ! npoints on demand; one large kernel and one synchronize beat hundreds of
   ! cache-sized launches. The CPU/OpenMP blocking is left untouched: cap the
   ! block at what the staging buffers were allocated for, but shrink it when
   ! there are too few points to give every thread a block of its own, so small
   ! atomic grids are not all evaluated by thread 0.
   IF (libxc_is_gpu()) THEN
      bsize = npoints
   ELSE
      bsize = MAX(1, MIN(workers%worker(1)%nb, (npoints + nthreads - 1)/nthreads))
   END IF
   nblocks = (npoints + bsize - 1)/bsize

   family = workers%family
   CALL libxc_family_flags(family, is_gga, is_mgga, TRIM(func_name))

   ! Below the cutoff LibXC leaves outputs at zero, so no test is needed here.
   ASSOCIATE (w => workers%worker(ithread + 1), has_laplace => workers%has_laplace, &
              no_exc => workers%no_exc, eps_rho => workers%epsilon_rho, &
              eps_tau => workers%epsilon_tau)
      !$OMP    DO
      DO ib = 1, nblocks
         IF (libxc_is_gpu()) CALL libxc_worker_resize_cuda(w, workers%family, workers%nspin, npoints)
         i0 = (ib - 1)*bsize
         nb = MIN(bsize, npoints - i0)

         ! stage this block's inputs in LibXC's layout
         DO i = 1, nb
            w%rho(1, i) = rho(i0 + i)
         END DO
         IF (is_gga) THEN
            DO i = 1, nb
               w%sigma(1, i) = norm_drho(i0 + i)**2
            END DO
         END IF
         IF (is_mgga) THEN
            DO i = 1, nb
               ii = i0 + i
               w%lapl(1, i) = laplace_rho(ii)
               w%tau(1, i) = tau(ii)
            END DO
            ! Continue at nonpositive tau without discarding finite density.
            ! The projection's Jacobian and curvature are applied after evaluation.
            DO i = 1, nb
               w%tau(1, i) = libxc_c2_project_tau(w%rho(1, i), w%sigma(1, i), tau(i0 + i), eps_rho, eps_tau)
            END DO
         END IF

         ! one LibXC call for the whole block, on the host or the device
         CALL libxc_worker_eval(w, workers%family, workers%nspin, grad_deriv, no_exc, nb)
         IF (is_mgga .AND. grad_deriv >= 1) THEN
            DO i = 1, nb
               CALL libxc_mgga_tau_pullback(w, i, [tau(i0 + i)], 1, grad_deriv, has_laplace, eps_rho, eps_tau)
            END DO
         END IF

         ! accumulate; each derivative is added by exactly one loop
         DO i = 1, nb
            ii = i0 + i
            e_0(ii) = e_0(ii) + sc*w%exc(i)*rho(ii)
         END DO
         IF (grad_deriv >= 1) THEN
            DO i = 1, nb
               ii = i0 + i
               e_rho(ii) = e_rho(ii) + sc*w%vrho(1, i)
            END DO
            IF (is_gga) THEN
               DO i = 1, nb
                  ii = i0 + i
                  e_ndrho(ii) = e_ndrho(ii) + sc*2.0_dp*w%vsigma(1, i)*norm_drho(ii)
               END DO
               DO i = 1, nb
                  ii = i0 + i
                  e_gamma(ii) = e_gamma(ii) + sc*w%vsigma(1, i)
               END DO
            END IF
            IF (is_mgga) THEN
               DO i = 1, nb
                  ii = i0 + i
                  e_tau(ii) = e_tau(ii) + sc*w%vtau(1, i)
               END DO
            END IF
            IF (is_mgga .AND. has_laplace) THEN
               DO i = 1, nb
                  ii = i0 + i
                  e_laplace_rho(ii) = e_laplace_rho(ii) + sc*w%vlapl(1, i)
               END DO
            END IF
         END IF
         IF (grad_deriv >= 2) THEN
            DO i = 1, nb
               ii = i0 + i
               e_rho_rho(ii) = e_rho_rho(ii) + sc*w%v2rho2(1, i)
            END DO
            IF (is_gga) THEN
               DO i = 1, nb
                  ii = i0 + i
                  e_ndrho_rho(ii) = e_ndrho_rho(ii) + sc*2.0_dp*w%v2rhosigma(1, i)*norm_drho(ii)
                  e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii) + &
                                      sc*2.0_dp*(2.0_dp*w%sigma(1, i)*w%v2sigma2(1, i) + w%vsigma(1, i))
               END DO
               DO i = 1, nb
                  ii = i0 + i
                  #:for descs, arr, idx in ugga_fill_2
                     e_${'_'.join(descs)}$ (ii) = e_${'_'.join(descs)}$ (ii) + sc*w%${arr}$ (${idx}$, i)
                  #:endfor
               END DO
            END IF
            IF (is_mgga) THEN
               DO i = 1, nb
                  ii = i0 + i
                  #:for descs, arr, idx in umgga_extra_2
                     e_${'_'.join(descs)}$ (ii) = e_${'_'.join(descs)}$ (ii) + sc*w%${arr}$ (${idx}$, i)
                  #:endfor
                  e_rho_tau(ii) = e_rho_tau(ii) + sc*w%v2rhotau(1, i)
                  e_ndrho_tau(ii) = e_ndrho_tau(ii) + sc*2.0_dp*w%v2sigmatau(1, i)*norm_drho(ii)
                  e_tau_tau(ii) = e_tau_tau(ii) + sc*w%v2tau2(1, i)
               END DO
            END IF
            IF (is_mgga .AND. has_laplace) THEN
               DO i = 1, nb
                  ii = i0 + i
                  e_rho_laplace_rho(ii) = e_rho_laplace_rho(ii) + sc*w%v2rholapl(1, i)
                  e_ndrho_laplace_rho(ii) = e_ndrho_laplace_rho(ii) + &
                                            sc*2.0_dp*w%v2sigmalapl(1, i)*norm_drho(ii)
                  e_laplace_rho_laplace_rho(ii) = e_laplace_rho_laplace_rho(ii) + sc*w%v2lapl2(1, i)
                  e_laplace_rho_tau(ii) = e_laplace_rho_tau(ii) + sc*w%v2lapltau(1, i)
               END DO
            END IF
         END IF
         IF (grad_deriv >= 3) THEN
            DO i = 1, nb
               ii = i0 + i
               e_rho_rho_rho(ii) = e_rho_rho_rho(ii) + sc*w%v3rho3(1, i)
            END DO
            IF (is_mgga) THEN
               DO i = 1, nb
                  ii = i0 + i
                  #:for descs, arr, idx in umgga_fill_3
                     e_${'_'.join(descs)}$ (ii) = e_${'_'.join(descs)}$ (ii) + sc*w%${arr}$ (${idx}$, i)
                  #:endfor
               END DO
            ELSE IF (is_gga) THEN
               DO i = 1, nb
                  ii = i0 + i
                  #:for descs, arr, idx in ugga_fill_3
                     e_${'_'.join(descs)}$ (ii) = e_${'_'.join(descs)}$ (ii) + sc*w%${arr}$ (${idx}$, i)
                  #:endfor
               END DO
            END IF
         END IF
      END DO
!$OMP    END DO
   END ASSOCIATE

END SUBROUTINE libxc_spin_unpolarized_calc
#endif

! **************************************************************************************************
!> \brief libxc exchange-correlation functionals
!> \param rhoa alpha density
!> \param rhob beta density
!> \param norm_drho ...
!> \param norm_drhoa norm of the gradient of the alpha density
!> \param norm_drhob norm of the gradient of the beta density
!> \param laplace_rhoa laplacian of the alpha density
!> \param laplace_rhob laplacian of the beta density
!> \param tau_a alpha kinetic-energy density
!> \param tau_b beta kinetic-energy density
!> \param e_0 energy density
!> \param e_rhoa derivative of the energy density with respect to rhoa
!> \param e_rhob derivative of the energy density with respect to rhob
!> \param e_ndrho derivative of the energy density with respect to ndrho
!> \param e_ndrhoa derivative of the energy density with respect to ndrhoa
!> \param e_ndrhob derivative of the energy density with respect to ndrhob
!> \param e_laplace_rhoa derivative of the energy density with respect to laplace_rhoa
!> \param e_laplace_rhob derivative of the energy density with respect to laplace_rhob
!> \param e_tau_a derivative of the energy density with respect to tau_a
!> \param e_tau_b derivative of the energy density with respect to tau_b
!> \param e_rhoa_rhoa derivative of the energy density with respect to rhoa_rhoa
!> \param e_rhoa_rhob derivative of the energy density with respect to rhoa_rhob
!> \param e_rhob_rhob derivative of the energy density with respect to rhob_rhob
!> \param e_ndrho_rhoa derivative of the energy density with respect to ndrho_rhoa
!> \param e_ndrho_rhob derivative of the energy density with respect to ndrho_rhob
!> \param e_ndrhoa_rhoa derivative of the energy density with respect to ndrhoa_rhoa
!> \param e_ndrhoa_rhob derivative of the energy density with respect to ndrhoa_rhob
!> \param e_ndrhob_rhoa derivative of the energy density with respect to ndrhob_rhoa
!> \param e_ndrhob_rhob derivative of the energy density with respect to ndrhob_rhob
!> \param e_ndrho_ndrho derivative of the energy density with respect to ndrho_ndrho
!> \param e_ndrho_ndrhoa derivative of the energy density with respect to ndrho_ndrhoa
!> \param e_ndrho_ndrhob derivative of the energy density with respect to ndrho_ndrhob
!> \param e_ndrhoa_ndrhoa derivative of the energy density with respect to ndrhoa_ndrhoa
!> \param e_ndrhoa_ndrhob derivative of the energy density with respect to ndrhoa_ndrhob
!> \param e_ndrhob_ndrhob derivative of the energy density with respect to ndrhob_ndrhob
!> \param e_rhoa_laplace_rhoa derivative of the energy density with respect to rhoa_laplace_rhoa
!> \param e_rhoa_laplace_rhob derivative of the energy density with respect to rhoa_laplace_rhob
!> \param e_rhob_laplace_rhoa derivative of the energy density with respect to rhob_laplace_rhoa
!> \param e_rhob_laplace_rhob derivative of the energy density with respect to rhob_laplace_rhob
!> \param e_rhoa_tau_a derivative of the energy density with respect to rhoa_tau_a
!> \param e_rhoa_tau_b derivative of the energy density with respect to rhoa_tau_b
!> \param e_rhob_tau_a derivative of the energy density with respect to rhob_tau_a
!> \param e_rhob_tau_b derivative of the energy density with respect to rhob_tau_b
!> \param e_ndrho_laplace_rhoa derivative of the energy density with respect to ndrho_laplace_rhoa
!> \param e_ndrho_laplace_rhob derivative of the energy density with respect to ndrho_laplace_rhob
!> \param e_ndrhoa_laplace_rhoa derivative of the energy density with respect to ndrhoa_laplace_rhoa
!> \param e_ndrhoa_laplace_rhob derivative of the energy density with respect to ndrhoa_laplace_rhob
!> \param e_ndrhob_laplace_rhoa derivative of the energy density with respect to ndrhob_laplace_rhoa
!> \param e_ndrhob_laplace_rhob derivative of the energy density with respect to ndrhob_laplace_rhob
!> \param e_ndrho_tau_a derivative of the energy density with respect to ndrho_tau_a
!> \param e_ndrho_tau_b derivative of the energy density with respect to ndrho_tau_b
!> \param e_ndrhoa_tau_a derivative of the energy density with respect to ndrhoa_tau_a
!> \param e_ndrhoa_tau_b derivative of the energy density with respect to ndrhoa_tau_b
!> \param e_ndrhob_tau_a derivative of the energy density with respect to ndrhob_tau_a
!> \param e_ndrhob_tau_b derivative of the energy density with respect to ndrhob_tau_b
!> \param e_laplace_rhoa_laplace_rhoa derivative of the energy density with respect to laplace_rhoa_laplace_rhoa
!> \param e_laplace_rhoa_laplace_rhob derivative of the energy density with respect to laplace_rhoa_laplace_rhob
!> \param e_laplace_rhob_laplace_rhob derivative of the energy density with respect to laplace_rhob_laplace_rhob
!> \param e_laplace_rhoa_tau_a derivative of the energy density with respect to laplace_rhoa_tau_a
!> \param e_laplace_rhoa_tau_b derivative of the energy density with respect to laplace_rhoa_tau_b
!> \param e_laplace_rhob_tau_a derivative of the energy density with respect to laplace_rhob_tau_a
!> \param e_laplace_rhob_tau_b derivative of the energy density with respect to laplace_rhob_tau_b
!> \param e_tau_a_tau_a derivative of the energy density with respect to tau_a_tau_a
!> \param e_tau_a_tau_b derivative of the energy density with respect to tau_a_tau_b
!> \param e_tau_b_tau_b derivative of the energy density with respect to tau_b_tau_b
!> \param e_rhoa_rhoa_rhoa derivative of the energy density with respect to rhoa_rhoa_rhoa
!> \param e_rhoa_rhoa_rhob derivative of the energy density with respect to rhoa_rhoa_rhob
!> \param e_rhoa_rhob_rhob derivative of the energy density with respect to rhoa_rhob_rhob
!> \param e_rhob_rhob_rhob derivative of the energy density with respect to rhob_rhob_rhob
!> \param grad_deriv maximum degree of the derivative that should be evaluated;
!>        all derivatives up to the given degree are evaluated, in a single
!>        LibXC call per block of grid points
!> \param npoints number of points on the grid
!> \param func_name name of the functional
!> \param sc scaling factor of the functional
!> \param workers cached LibXC functional objects and staging buffers, one per thread
!> \author F. Tran
! **************************************************************************************************
#if defined (__LIBXC)
SUBROUTINE libxc_spin_polarized_calc(rhoa, rhob, norm_drho, norm_drhoa, &
                                     norm_drhob, laplace_rhoa, laplace_rhob, tau_a, tau_b, &
                                     e_0, e_rhoa, e_rhob, e_ndrho, e_ndrhoa, e_ndrhob, &
                                     e_laplace_rhoa, e_laplace_rhob, e_tau_a, e_tau_b, &
                                     e_rhoa_rhoa, e_rhoa_rhob, e_rhob_rhob, &
                                     e_ndrho_rhoa, e_ndrho_rhob, e_ndrhoa_rhoa, &
                                     e_ndrhoa_rhob, e_ndrhob_rhoa, e_ndrhob_rhob, &
                                     e_ndrho_ndrho, e_ndrho_ndrhoa, e_ndrho_ndrhob, &
                                     e_ndrhoa_ndrhoa, e_ndrhoa_ndrhob, e_ndrhob_ndrhob, &
                                     e_rhoa_laplace_rhoa, e_rhoa_laplace_rhob, &
                                     e_rhob_laplace_rhoa, e_rhob_laplace_rhob, &
                                     e_rhoa_tau_a, e_rhoa_tau_b, e_rhob_tau_a, e_rhob_tau_b, &
                                     e_ndrho_laplace_rhoa, e_ndrho_laplace_rhob, &
                                     e_ndrhoa_laplace_rhoa, e_ndrhoa_laplace_rhob, &
                                     e_ndrhob_laplace_rhoa, e_ndrhob_laplace_rhob, &
                                     e_ndrho_tau_a, e_ndrho_tau_b, &
                                     e_ndrhoa_tau_a, e_ndrhoa_tau_b, &
                                     e_ndrhob_tau_a, e_ndrhob_tau_b, &
                                     e_laplace_rhoa_laplace_rhoa, &
                                     e_laplace_rhoa_laplace_rhob, &
                                     e_laplace_rhob_laplace_rhob, &
                                     e_laplace_rhoa_tau_a, e_laplace_rhoa_tau_b, &
                                     e_laplace_rhob_tau_a, e_laplace_rhob_tau_b, &
                                     e_tau_a_tau_a, e_tau_a_tau_b, e_tau_b_tau_b, &
                                     e_rhoa_rhoa_rhoa, e_rhoa_rhoa_rhob, &
                                     e_rhoa_rhob_rhob, e_rhob_rhob_rhob, &
                                        #:for descs, arr, idx in mgga_all
                                     e_${'_'.join(descs)}$, &
                                        #:endfor
                                     grad_deriv, npoints, func_name, sc, workers)

REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rhoa, rhob, norm_drho, norm_drhoa, &
                                                      norm_drhob, laplace_rhoa, &
                                                      laplace_rhob, tau_a, tau_b
REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_0, e_rhoa, e_rhob, e_ndrho, e_ndrhoa, &
                                             e_ndrhob, e_laplace_rhoa, e_laplace_rhob, e_tau_a, e_tau_b, e_rhoa_rhoa, e_rhoa_rhob, &
                                             e_rhob_rhob, e_ndrho_rhoa, e_ndrho_rhob, e_ndrhoa_rhoa, e_ndrhoa_rhob, e_ndrhob_rhoa, &
                                              e_ndrhob_rhob, e_ndrho_ndrho, e_ndrho_ndrhoa, e_ndrho_ndrhob, e_ndrhoa_ndrhoa, &
                                              e_ndrhoa_ndrhob, e_ndrhob_ndrhob, e_rhoa_laplace_rhoa, e_rhoa_laplace_rhob, &
                                              e_rhob_laplace_rhoa, e_rhob_laplace_rhob, e_rhoa_tau_a, e_rhoa_tau_b, e_rhob_tau_a, &
                                              e_rhob_tau_b, e_ndrho_laplace_rhoa, e_ndrho_laplace_rhob, e_ndrhoa_laplace_rhoa
REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_ndrhoa_laplace_rhob, e_ndrhob_laplace_rhoa, &
                                              e_ndrhob_laplace_rhob, e_ndrho_tau_a, e_ndrho_tau_b, e_ndrhoa_tau_a, e_ndrhoa_tau_b, &
                                         e_ndrhob_tau_a, e_ndrhob_tau_b, e_laplace_rhoa_laplace_rhoa, e_laplace_rhoa_laplace_rhob, &
                                              e_laplace_rhob_laplace_rhob, e_laplace_rhoa_tau_a, e_laplace_rhoa_tau_b, &
                                          e_laplace_rhob_tau_a, e_laplace_rhob_tau_b, e_tau_a_tau_a, e_tau_a_tau_b, e_tau_b_tau_b, &
                                              e_rhoa_rhoa_rhoa, e_rhoa_rhoa_rhob, e_rhoa_rhob_rhob, e_rhob_rhob_rhob
#:for descs, arr, idx in mgga_all
   REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_${'_'.join(descs)}$
#:endfor
INTEGER, INTENT(in)                                :: grad_deriv, npoints
CHARACTER(LEN=default_string_length), INTENT(IN)   :: func_name
REAL(KIND=dp), INTENT(in)                          :: sc
TYPE(libxc_worker_set_type), INTENT(INOUT)         :: workers

INTEGER                                            :: bsize, family, i, i0, ib, ii, ithread, &
                                                      nb, nblocks, nthreads
LOGICAL                                            :: is_gga, is_mgga

ithread = 0
nthreads = 1
!$ ithread = omp_get_thread_num()
!$ nthreads = omp_get_num_threads()
CPASSERT(ithread < SIZE(workers%worker))

! See the spin-unpolarized calc: the device is activated per evaluation.

! See the spin-unpolarized calc: on the GPU path the whole grid is evaluated
! in a single LibXC call, while the CPU/OpenMP path keeps its per-footprint
! block cap (shrunk so every thread gets a block on small grids).
IF (libxc_is_gpu()) THEN
   bsize = npoints
ELSE
   bsize = MAX(1, MIN(workers%worker(1)%nb, (npoints + nthreads - 1)/nthreads))
END IF
nblocks = (npoints + bsize - 1)/bsize

family = workers%family
CALL libxc_family_flags(family, is_gga, is_mgga, TRIM(func_name))

! Density and tau cutoffs are applied by LibXC itself, as in the
! spin-unpolarized case.
ASSOCIATE (w => workers%worker(ithread + 1), has_laplace => workers%has_laplace, &
           no_exc => workers%no_exc, eps_rho => workers%epsilon_rho, &
           eps_tau => workers%epsilon_tau)
!$OMP    DO
   DO ib = 1, nblocks
      IF (libxc_is_gpu()) CALL libxc_worker_resize_cuda(w, workers%family, workers%nspin, npoints)
      i0 = (ib - 1)*bsize
      nb = MIN(bsize, npoints - i0)

      ! stage this block's inputs in LibXC's layout
      DO i = 1, nb
         ii = i0 + i
         w%rho(1, i) = MAX(rhoa(ii), 0.0_dp)
         w%rho(2, i) = MAX(rhob(ii), 0.0_dp)
      END DO
      IF (is_gga) THEN
         ! CP2K works with gradient norms, LibXC with their contractions.
         DO i = 1, nb
            ii = i0 + i
            w%nda(i) = MAX(norm_drhoa(ii), EPSILON(0.0_dp)*1.e4_dp)
            w%ndb(i) = MAX(norm_drhob(ii), EPSILON(0.0_dp)*1.e4_dp)
            w%nd(i) = MAX(norm_drho(ii), EPSILON(0.0_dp)*1.e4_dp)
            w%sigma(1, i) = w%nda(i)**2
            w%sigma(3, i) = w%ndb(i)**2
            w%sigma(2, i) = 0.5_dp*(w%nd(i)**2 - w%sigma(1, i) - w%sigma(3, i))
         END DO
      END IF
      IF (is_mgga) THEN
         DO i = 1, nb
            ii = i0 + i
            w%lapl(1, i) = laplace_rhoa(ii)
            w%lapl(2, i) = laplace_rhob(ii)
            w%tau(1, i) = MAX(tau_a(ii), 0.0_dp)
            w%tau(2, i) = MAX(tau_b(ii), 0.0_dp)
         END DO
         ! Apply the same continuous kinetic-density projection per spin channel.
         DO i = 1, nb
            w%tau(1, i) = libxc_c2_project_tau(w%rho(1, i), w%sigma(1, i), tau_a(i0 + i), eps_rho, eps_tau)
            w%tau(2, i) = libxc_c2_project_tau(w%rho(2, i), w%sigma(3, i), tau_b(i0 + i), eps_rho, eps_tau)
         END DO
      END IF

      ! one LibXC call for the whole block, on the host or the device
      CALL libxc_worker_eval(w, workers%family, workers%nspin, grad_deriv, no_exc, nb)
      IF (grad_deriv >= 1 .AND. .NOT. no_exc) THEN
         DO i = 1, nb
            CALL libxc_spin_floor_pullback(w, i, grad_deriv, has_laplace, eps_rho)
         END DO
      END IF
      IF (is_mgga .AND. grad_deriv >= 1) THEN
         DO i = 1, nb
            CALL libxc_mgga_tau_pullback(w, i, [tau_a(i0 + i), tau_b(i0 + i)], 2, &
                                         grad_deriv, has_laplace, eps_rho, eps_tau)
         END DO
      END IF
      IF (grad_deriv >= 1) THEN
         DO i = 1, nb
            CALL libxc_spin_density_pullback(w, i, [rhoa(i0 + i), rhob(i0 + i)], grad_deriv)
         END DO
      END IF
      ! accumulate; one loop per group of derivatives, selected by the
      ! order requested and the family
      DO i = 1, nb
         ii = i0 + i
         e_0(ii) = e_0(ii) + sc*w%exc(i)*(w%rho(1, i) + w%rho(2, i))
      END DO
      IF (grad_deriv >= 1) THEN
         DO i = 1, nb
            ii = i0 + i
            e_rhoa(ii) = e_rhoa(ii) + sc*w%vrho(1, i)
            e_rhob(ii) = e_rhob(ii) + sc*w%vrho(2, i)
         END DO
         IF (is_gga) THEN
            DO i = 1, nb
               ii = i0 + i
               e_ndrho(ii) = e_ndrho(ii) + sc*w%vsigma(2, i)*w%nd(i)
               e_ndrhoa(ii) = e_ndrhoa(ii) + &
                              sc*(2.0_dp*w%vsigma(1, i) - w%vsigma(2, i))*w%nda(i)
               e_ndrhob(ii) = e_ndrhob(ii) + &
                              sc*(2.0_dp*w%vsigma(3, i) - w%vsigma(2, i))*w%ndb(i)
            END DO
         END IF
         IF (is_mgga) THEN
            DO i = 1, nb
               ii = i0 + i
               e_tau_a(ii) = e_tau_a(ii) + sc*w%vtau(1, i)
               e_tau_b(ii) = e_tau_b(ii) + sc*w%vtau(2, i)
            END DO
         END IF
         IF (is_mgga .AND. has_laplace) THEN
            DO i = 1, nb
               ii = i0 + i
               e_laplace_rhoa(ii) = e_laplace_rhoa(ii) + sc*w%vlapl(1, i)
               e_laplace_rhob(ii) = e_laplace_rhob(ii) + sc*w%vlapl(2, i)
            END DO
         END IF
      END IF
      IF (grad_deriv >= 2) THEN
         DO i = 1, nb
            ii = i0 + i
            e_rhoa_rhoa(ii) = e_rhoa_rhoa(ii) + sc*w%v2rho2(1, i)
            e_rhoa_rhob(ii) = e_rhoa_rhob(ii) + sc*w%v2rho2(2, i)
            e_rhob_rhob(ii) = e_rhob_rhob(ii) + sc*w%v2rho2(3, i)
         END DO
         IF (is_gga) THEN
            DO i = 1, nb
               ii = i0 + i
               e_ndrho_rhoa(ii) = e_ndrho_rhoa(ii) + sc*w%v2rhosigma(2, i)*w%nd(i)
               e_ndrho_rhob(ii) = e_ndrho_rhob(ii) + sc*w%v2rhosigma(5, i)*w%nd(i)
               e_ndrhoa_rhoa(ii) = e_ndrhoa_rhoa(ii) + &
                                   sc*(2.0_dp*w%v2rhosigma(1, i) - w%v2rhosigma(2, i))*w%nda(i)
               e_ndrhoa_rhob(ii) = e_ndrhoa_rhob(ii) + &
                                   sc*(2.0_dp*w%v2rhosigma(4, i) - w%v2rhosigma(5, i))*w%nda(i)
               e_ndrhob_rhoa(ii) = e_ndrhob_rhoa(ii) + &
                                   sc*(2.0_dp*w%v2rhosigma(3, i) - w%v2rhosigma(2, i))*w%ndb(i)
               e_ndrhob_rhob(ii) = e_ndrhob_rhob(ii) + &
                                   sc*(2.0_dp*w%v2rhosigma(6, i) - w%v2rhosigma(5, i))*w%ndb(i)
               e_ndrho_ndrho(ii) = e_ndrho_ndrho(ii) + &
                                   sc*(w%vsigma(2, i) + w%nd(i)**2*w%v2sigma2(4, i))
               e_ndrho_ndrhoa(ii) = e_ndrho_ndrhoa(ii) + &
                                    sc*(2.0_dp*w%v2sigma2(2, i) - w%v2sigma2(4, i))*w%nd(i)*w%nda(i)
               e_ndrho_ndrhob(ii) = e_ndrho_ndrhob(ii) + &
                                    sc*(2.0_dp*w%v2sigma2(5, i) - w%v2sigma2(4, i))*w%nd(i)*w%ndb(i)
               e_ndrhoa_ndrhoa(ii) = e_ndrhoa_ndrhoa(ii) + &
                                     sc*(2.0_dp*w%vsigma(1, i) - w%vsigma(2, i) + w%nda(i)**2*( &
                                         4.0_dp*w%v2sigma2(1, i) - 4.0_dp*w%v2sigma2(2, i) + w%v2sigma2(4, i)))
               e_ndrhoa_ndrhob(ii) = e_ndrhoa_ndrhob(ii) + &
                                     sc*(4.0_dp*w%v2sigma2(3, i) - 2.0_dp*w%v2sigma2(2, i) - &
                                         2.0_dp*w%v2sigma2(5, i) + w%v2sigma2(4, i))*w%nda(i)*w%ndb(i)
               e_ndrhob_ndrhob(ii) = e_ndrhob_ndrhob(ii) + &
                                     sc*(2.0_dp*w%vsigma(3, i) - w%vsigma(2, i) + w%ndb(i)**2*( &
                                         4.0_dp*w%v2sigma2(6, i) - 4.0_dp*w%v2sigma2(5, i) + w%v2sigma2(4, i)))
            END DO
         END IF
         IF (is_mgga) THEN
            DO i = 1, nb
               ii = i0 + i
               e_rhoa_tau_a(ii) = e_rhoa_tau_a(ii) + sc*w%v2rhotau(1, i)
               e_rhoa_tau_b(ii) = e_rhoa_tau_b(ii) + sc*w%v2rhotau(2, i)
               e_rhob_tau_a(ii) = e_rhob_tau_a(ii) + sc*w%v2rhotau(3, i)
               e_rhob_tau_b(ii) = e_rhob_tau_b(ii) + sc*w%v2rhotau(4, i)
               e_ndrho_tau_a(ii) = e_ndrho_tau_a(ii) + sc*w%v2sigmatau(3, i)*w%nd(i)
               e_ndrho_tau_b(ii) = e_ndrho_tau_b(ii) + sc*w%v2sigmatau(4, i)*w%nd(i)
               e_ndrhoa_tau_a(ii) = e_ndrhoa_tau_a(ii) + &
                                    sc*(2.0_dp*w%v2sigmatau(1, i) - w%v2sigmatau(3, i))*w%nda(i)
               e_ndrhoa_tau_b(ii) = e_ndrhoa_tau_b(ii) + &
                                    sc*(2.0_dp*w%v2sigmatau(2, i) - w%v2sigmatau(4, i))*w%nda(i)
               e_ndrhob_tau_a(ii) = e_ndrhob_tau_a(ii) + &
                                    sc*(2.0_dp*w%v2sigmatau(5, i) - w%v2sigmatau(3, i))*w%ndb(i)
               e_ndrhob_tau_b(ii) = e_ndrhob_tau_b(ii) + &
                                    sc*(2.0_dp*w%v2sigmatau(6, i) - w%v2sigmatau(4, i))*w%ndb(i)
               e_tau_a_tau_a(ii) = e_tau_a_tau_a(ii) + sc*w%v2tau2(1, i)
               e_tau_a_tau_b(ii) = e_tau_a_tau_b(ii) + sc*w%v2tau2(2, i)
               e_tau_b_tau_b(ii) = e_tau_b_tau_b(ii) + sc*w%v2tau2(3, i)
            END DO
         END IF
         IF (is_mgga .AND. has_laplace) THEN
            DO i = 1, nb
               ii = i0 + i
               e_rhoa_laplace_rhoa(ii) = e_rhoa_laplace_rhoa(ii) + sc*w%v2rholapl(1, i)
               e_rhoa_laplace_rhob(ii) = e_rhoa_laplace_rhob(ii) + sc*w%v2rholapl(2, i)
               e_rhob_laplace_rhoa(ii) = e_rhob_laplace_rhoa(ii) + sc*w%v2rholapl(3, i)
               e_rhob_laplace_rhob(ii) = e_rhob_laplace_rhob(ii) + sc*w%v2rholapl(4, i)
               e_ndrho_laplace_rhoa(ii) = e_ndrho_laplace_rhoa(ii) + sc*w%v2sigmalapl(3, i)*w%nd(i)
               e_ndrho_laplace_rhob(ii) = e_ndrho_laplace_rhob(ii) + sc*w%v2sigmalapl(4, i)*w%nd(i)
               e_ndrhoa_laplace_rhoa(ii) = e_ndrhoa_laplace_rhoa(ii) + &
                                           sc*(2.0_dp*w%v2sigmalapl(1, i) - w%v2sigmalapl(3, i))*w%nda(i)
               e_ndrhoa_laplace_rhob(ii) = e_ndrhoa_laplace_rhob(ii) + &
                                           sc*(2.0_dp*w%v2sigmalapl(2, i) - w%v2sigmalapl(4, i))*w%nda(i)
               e_ndrhob_laplace_rhoa(ii) = e_ndrhob_laplace_rhoa(ii) + &
                                           sc*(2.0_dp*w%v2sigmalapl(5, i) - w%v2sigmalapl(3, i))*w%ndb(i)
               e_ndrhob_laplace_rhob(ii) = e_ndrhob_laplace_rhob(ii) + &
                                           sc*(2.0_dp*w%v2sigmalapl(6, i) - w%v2sigmalapl(4, i))*w%ndb(i)
               e_laplace_rhoa_laplace_rhoa(ii) = e_laplace_rhoa_laplace_rhoa(ii) + sc*w%v2lapl2(1, i)
               e_laplace_rhoa_laplace_rhob(ii) = e_laplace_rhoa_laplace_rhob(ii) + sc*w%v2lapl2(2, i)
               e_laplace_rhob_laplace_rhob(ii) = e_laplace_rhob_laplace_rhob(ii) + sc*w%v2lapl2(3, i)
               e_laplace_rhoa_tau_a(ii) = e_laplace_rhoa_tau_a(ii) + sc*w%v2lapltau(1, i)
               e_laplace_rhoa_tau_b(ii) = e_laplace_rhoa_tau_b(ii) + sc*w%v2lapltau(2, i)
               e_laplace_rhob_tau_a(ii) = e_laplace_rhob_tau_a(ii) + sc*w%v2lapltau(3, i)
               e_laplace_rhob_tau_b(ii) = e_laplace_rhob_tau_b(ii) + sc*w%v2lapltau(4, i)
            END DO
         END IF
      END IF
      IF (grad_deriv >= 3) THEN
         DO i = 1, nb
            ii = i0 + i
            e_rhoa_rhoa_rhoa(ii) = e_rhoa_rhoa_rhoa(ii) + sc*w%v3rho3(1, i)
            e_rhoa_rhoa_rhob(ii) = e_rhoa_rhoa_rhob(ii) + sc*w%v3rho3(2, i)
            e_rhoa_rhob_rhob(ii) = e_rhoa_rhob_rhob(ii) + sc*w%v3rho3(3, i)
            e_rhob_rhob_rhob(ii) = e_rhob_rhob_rhob(ii) + sc*w%v3rho3(4, i)
         END DO
      END IF
      IF (grad_deriv >= 1 .AND. is_gga) THEN
         DO i = 1, nb
            ii = i0 + i
            #:for descs, arr, idx in gamma_only_1
               e_${'_'.join(descs)}$ (ii) = e_${'_'.join(descs)}$ (ii) + sc*w%${arr}$ (${idx}$, i)
            #:endfor
         END DO
      END IF
      IF (grad_deriv >= 2 .AND. is_gga) THEN
         DO i = 1, nb
            ii = i0 + i
            #:for descs, arr, idx in gamma_only_2
               e_${'_'.join(descs)}$ (ii) = e_${'_'.join(descs)}$ (ii) + sc*w%${arr}$ (${idx}$, i)
            #:endfor
         END DO
      END IF
      IF (grad_deriv >= 2 .AND. is_mgga) THEN
         DO i = 1, nb
            ii = i0 + i
            #:for descs, arr, idx in [e for e in mgga_fill_2 if e not in gga_fill_2]
               e_${'_'.join(descs)}$ (ii) = e_${'_'.join(descs)}$ (ii) + sc*w%${arr}$ (${idx}$, i)
            #:endfor
         END DO
      END IF
      IF (grad_deriv >= 3 .AND. is_mgga) THEN
         DO i = 1, nb
            ii = i0 + i
            #:for descs, arr, idx in mgga_fill_3
               e_${'_'.join(descs)}$ (ii) = e_${'_'.join(descs)}$ (ii) + sc*w%${arr}$ (${idx}$, i)
            #:endfor
         END DO
      END IF
      IF (grad_deriv >= 3 .AND. is_gga .AND. .NOT. is_mgga) THEN
         DO i = 1, nb
            ii = i0 + i
            #:for descs, arr, idx in gamma_only_3
               e_${'_'.join(descs)}$ (ii) = e_${'_'.join(descs)}$ (ii) + sc*w%${arr}$ (${idx}$, i)
            #:endfor
         END DO
      END IF
   END DO
!$OMP    END DO
END ASSOCIATE

END SUBROUTINE libxc_spin_polarized_calc
#endif

END MODULE xc_libxc
