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
#if defined (__LIBXC)
   USE input_keyword_types, ONLY: keyword_create, &
                                  keyword_release, &
                                  keyword_type
   USE iso_c_binding, ONLY: C_SIZE_T, C_INT, C_DOUBLE
   USE xc_derivative_desc, ONLY: &
      deriv_rho, deriv_rhoa, deriv_rhob, &
      deriv_norm_drhoa, deriv_norm_drhob, deriv_norm_drho, deriv_tau_a, deriv_tau_b, deriv_tau, &
      deriv_laplace_rho, deriv_laplace_rhoa, deriv_laplace_rhob
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
                            xc_f03_gga_fxc, &
                            xc_f03_gga_vxc, &
                            xc_f03_gga_vxc_fxc, &
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
                            xc_f03_mgga_fxc, &
                            xc_f03_mgga_vxc, &
                            xc_f03_mgga_vxc_fxc, &
                            XC_POLARIZED, &
                            XC_UNPOLARIZED, &
                            XC_FAMILY_LDA, &
                            XC_FAMILY_GGA, &
                            XC_FAMILY_MGGA, &
                            XC_FAMILY_HYB_LDA, &
                            XC_FAMILY_HYB_GGA, &
                            XC_FAMILY_HYB_MGGA, &
                            XC_CORRELATION, &
                            XC_EXCHANGE, &
                            XC_EXCHANGE_CORRELATION, &
                            XC_KINETIC, &
                            xc_libxc_wrap_info_refs, &
                            xc_libxc_wrap_version, &
                            xc_libxc_wrap_library_reference, &
                            xc_libxc_wrap_functional_get_number, &
                            xc_libxc_wrap_info_needs_laplace, &
                            xc_libxc_wrap_info_no_exc, &
                            xc_libxc_wrap_set_thresholds, &
                            xc_f03_func_set_ext_params, &
                            xc_libxc_wrap_is_under_development, &
                            xc_libxc_get_reference_length, &
                            xc_libxc_check_functional
#endif
!$ USE OMP_LIB, ONLY: omp_get_max_threads, &
!$                    omp_get_num_threads, &
!$                    omp_get_thread_num

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_libxc'

   PUBLIC :: libxc_spin_unpolarized_info, libxc_spin_unpolarized_eval, &
             libxc_spin_polarized_info, libxc_spin_polarized_eval, &
             libxc_version_info, libxc_library_reference, &
             libxc_get_reference_length, libxc_add_sections, &
             libxc_check_existence_in_libxc, libxc_release_workers

#if defined (__LIBXC)
   ! Number of grid points handed to LibXC in a single call. LibXC amortizes its
   ! per-call work (argument checks, output initialization and, for functionals
   ! built by mixing components, a full allocate/free cycle of the component
   ! buffers) over the points in the call, so evaluating point by point paid that
   ! cost once per point. Blocking rather than passing the whole grid keeps the
   ! staging buffers inside the cache; the timing is flat between roughly 128 and
   ! 2048 points, so the exact value is not critical.
   INTEGER, PARAMETER, PRIVATE :: libxc_block_size = 512

! **************************************************************************************************
!> \brief One LibXC functional object plus the buffers used to hand it a block of
!>        grid points. One of these is kept per OpenMP thread so that neither the
!>        functional object nor the buffers are rebuilt between batches.
! **************************************************************************************************
   TYPE libxc_worker_type
      TYPE(xc_f03_func_t)                                :: func = xc_f03_func_t()
      TYPE(xc_f03_func_info_t)                           :: info = xc_f03_func_info_t()
      LOGICAL                                            :: is_init = .FALSE.
      ! staging buffers, dimensioned (component, point)
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: rho, sigma, lapl, tau
      ! the floored gradient norms the spin-polarized chain rule is written in
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: nd, nda, ndb
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: exc
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: vrho, vsigma, vlapl, vtau
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: v2rho2, v2rhosigma, v2sigma2, &
                                                            v2rholapl, v2rhotau, v2sigmalapl, &
                                                            v2sigmatau, v2lapl2, v2lapltau, v2tau2
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: v3rho3
   END TYPE libxc_worker_type

! **************************************************************************************************
!> \brief All threads' workers for one distinct functional setup, together with
!>        the key that identifies that setup.
! **************************************************************************************************
   TYPE libxc_worker_set_type
      INTEGER                                            :: func_id = -1
      INTEGER                                            :: nspin = -1
      INTEGER                                            :: family = -1
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
         CALL xc_f03_func_init(xc_func, func_id, XC_POLARIZED)
      ELSE
         CALL xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
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
         CALL xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
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
      CALL xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
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
      CALL xc_f03_func_init(xc_func, func_id, XC_POLARIZED)
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
!> \brief Allocates the staging buffers of a single worker.
!> \param worker the worker
!> \param family LibXC functional family
!> \param nspin XC_UNPOLARIZED or XC_POLARIZED
!> \author S. Lehtola
!> \note The buffers follow LibXC's own layout, i.e. the leading dimension is the
!>       number of components of the quantity and the points run along the second
!>       dimension. They are allocated once per worker and reused for every batch.
! **************************************************************************************************
   SUBROUTINE libxc_worker_alloc(worker, family, nspin)
      TYPE(libxc_worker_type), INTENT(INOUT)             :: worker
      INTEGER, INTENT(IN)                                :: family, nspin

      INTEGER                                            :: nb, nl, np, ns

      nb = libxc_block_size
      ! number of components per point, cf. the LibXC manual
      IF (nspin == XC_POLARIZED) THEN
         np = 2   ! rho, vrho, lapl, tau, vlapl, vtau
         ns = 3   ! sigma, vsigma, v2rho2, v2lapl2, v2tau2
         nl = 4   ! v2rholapl, v2rhotau, v2lapltau, v3rho3
      ELSE
         np = 1
         ns = 1
         nl = 1
      END IF

      ALLOCATE (worker%rho(np, nb))
      ALLOCATE (worker%exc(nb))
      ALLOCATE (worker%vrho(np, nb))
      ALLOCATE (worker%v2rho2(ns, nb))
      ALLOCATE (worker%v3rho3(nl, nb))

      SELECT CASE (family)
      CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
         ALLOCATE (worker%sigma(ns, nb))
         ALLOCATE (worker%vsigma(ns, nb))
         IF (nspin == XC_POLARIZED) THEN
            ALLOCATE (worker%nd(nb), worker%nda(nb), worker%ndb(nb))
         END IF
         ALLOCATE (worker%v2rhosigma(MERGE(6, 1, nspin == XC_POLARIZED), nb))
         ALLOCATE (worker%v2sigma2(MERGE(6, 1, nspin == XC_POLARIZED), nb))
      END SELECT

      SELECT CASE (family)
      CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
         ALLOCATE (worker%lapl(np, nb))
         ALLOCATE (worker%tau(np, nb))
         ALLOCATE (worker%vlapl(np, nb))
         ALLOCATE (worker%vtau(np, nb))
         ALLOCATE (worker%v2rholapl(nl, nb))
         ALLOCATE (worker%v2rhotau(nl, nb))
         ALLOCATE (worker%v2sigmalapl(MERGE(6, 1, nspin == XC_POLARIZED), nb))
         ALLOCATE (worker%v2sigmatau(MERGE(6, 1, nspin == XC_POLARIZED), nb))
         ALLOCATE (worker%v2lapl2(ns, nb))
         ALLOCATE (worker%v2lapltau(nl, nb))
         ALLOCATE (worker%v2tau2(ns, nb))
         ! LibXC leaves vlapl untouched for functionals that do not use the
         ! Laplacian, and CP2K then ignores it; zero it once so it never holds
         ! uninitialized memory.
         worker%vlapl = 0.0_dp
      END SELECT

   END SUBROUTINE libxc_worker_alloc

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
!> \return index into libxc_cache
!> \author S. Lehtola
!> \note Must be called from a serial region: it is the only writer of the cache,
!>       and keeping it serial is what makes the cache lock-free. The workers it
!>       returns are then used concurrently, one per thread.
! **************************************************************************************************
   FUNCTION libxc_get_workers(func_name, nspin, libxc_params, epsilon_rho, epsilon_tau) RESULT(idx)
      CHARACTER(LEN=*), INTENT(IN)                       :: func_name
      INTEGER, INTENT(IN)                                :: nspin
      TYPE(section_vals_type), INTENT(IN), POINTER       :: libxc_params
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_rho, epsilon_tau
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
            ALLOCATE (ws%worker(nthreads))

            DO ithread = 1, nthreads
               CALL xc_f03_func_init(ws%worker(ithread)%func, func_id, nspin)
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
               CALL libxc_worker_alloc(ws%worker(ithread), ws%family, nspin)
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
!> \param grad_deriv degree of the derivative that should be evaluated;
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
      TYPE(xc_derivative_type), POINTER        :: deriv
      TYPE(xc_f03_func_info_t)                 :: xc_info

      CALL timeset(routineN, handle)

      ! Only "everything up to grad_deriv" is supported. Asking for a single
      ! derivative order in isolation is rejected here rather than further down,
      ! where it would quietly evaluate nothing at all.
      IF (grad_deriv < 0) THEN
         CPABORT("Evaluating a single derivative order is not supported.")
      END IF
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
         iw = libxc_get_workers(func_name, XC_UNPOLARIZED, no_params, epsilon_rho, epsilon_tau)
      ELSE
         iw = libxc_get_workers(func_name, XC_UNPOLARIZED, libxc_params, epsilon_rho, epsilon_tau)
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

      IF (grad_deriv >= 0) THEN
         deriv => xc_dset_get_derivative(deriv_set, [INTEGER::], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_0)
      END IF
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
         CASE (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho)
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
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_tau], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_tau)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_norm_drho, deriv_tau], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_ndrho_tau)
            deriv => xc_dset_get_derivative(deriv_set, [deriv_tau, deriv_tau], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_tau_tau)
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
            deriv => xc_dset_get_derivative(deriv_set, [deriv_rho, deriv_rho, deriv_rho], &
                                            allocate_deriv=.TRUE.)
            CALL xc_derivative_get(deriv, deriv_data=e_rho_rho_rho)
         CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            CPABORT("derivatives larger than 2 not implemented")
         CASE default
            CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
         END SELECT
      END IF

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP SHARED(rho,norm_drho,laplace_rho,tau,e_0,e_rho,e_ndrho,e_laplace_rho),&
!$OMP SHARED(e_tau,e_rho_rho,e_ndrho_rho,e_ndrho_ndrho,e_rho_laplace_rho),&
!$OMP SHARED(e_rho_tau,e_ndrho_laplace_rho,e_ndrho_tau,e_laplace_rho_laplace_rho),&
!$OMP SHARED(e_laplace_rho_tau,e_tau_tau,e_rho_rho_rho),&
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
      TYPE(xc_derivative_type), POINTER        :: deriv
      TYPE(xc_f03_func_info_t)                 :: xc_info

      CALL timeset(routineN, handle)

      ! Only "everything up to grad_deriv" is supported. Asking for a single
      ! derivative order in isolation is rejected here rather than further down,
      ! where it would quietly evaluate nothing at all.
      IF (grad_deriv < 0) THEN
         CPABORT("Evaluating a single derivative order is not supported.")
      END IF
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
         iw = libxc_get_workers(func_name, XC_POLARIZED, no_params, epsilon_rho, epsilon_tau)
      ELSE
         iw = libxc_get_workers(func_name, XC_POLARIZED, libxc_params, epsilon_rho, epsilon_tau)
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

      IF (grad_deriv >= 0) THEN
         deriv => xc_dset_get_derivative(deriv_set, [INTEGER::], &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_0)
      END IF
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
         CASE (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            CPABORT("derivatives larger than 2 not implemented")
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
!> \param grad_deriv degree of the derivative that should be evaluated;
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
                                          grad_deriv, npoints, func_name, sc, workers)

      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho, norm_drho, laplace_rho, tau
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT) :: e_0, e_rho, e_ndrho, e_laplace_rho, e_tau, &
         e_rho_rho, e_ndrho_rho, e_ndrho_ndrho, e_rho_laplace_rho, e_rho_tau, e_ndrho_laplace_rho, &
         e_ndrho_tau, e_laplace_rho_laplace_rho, e_laplace_rho_tau, e_tau_tau, e_rho_rho_rho
      INTEGER, INTENT(in)                                :: grad_deriv, npoints
      CHARACTER(LEN=default_string_length), INTENT(IN)   :: func_name
      REAL(KIND=dp), INTENT(in)                          :: sc
      TYPE(libxc_worker_set_type), INTENT(INOUT)         :: workers

      INTEGER                                            :: bsize, family, i, i0, ib, ii, ithread, &
                                                            nb, nblocks, nthreads
      INTEGER(C_SIZE_T)                                  :: np
      LOGICAL                                            :: is_gga, is_mgga

      ithread = 0
      nthreads = 1
!$    ithread = omp_get_thread_num()
!$    nthreads = omp_get_num_threads()
      CPASSERT(ithread < SIZE(workers%worker))

      ! Blocks are capped by libxc_block_size, which is what the staging buffers
      ! were allocated for, but shrunk when there are too few points to give every
      ! thread a block of its own: small atomic grids would otherwise all be
      ! evaluated by thread 0.
      bsize = MAX(1, MIN(libxc_block_size, (npoints + nthreads - 1)/nthreads))
      nblocks = (npoints + bsize - 1)/bsize

      family = workers%family
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
         CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
      END SELECT

      ! Points below the density cutoff need no test here. The cutoffs were handed
      ! to the functional object in xc_libxc_wrap_set_thresholds, and LibXC leaves
      ! the outputs of the points it screens at zero, so accumulating them is a
      ! no-op.
      ASSOCIATE (w => workers%worker(ithread + 1), has_laplace => workers%has_laplace, &
                 no_exc => workers%no_exc, eps_rho => workers%epsilon_rho, &
                 eps_tau => workers%epsilon_tau)
!$OMP    DO
         DO ib = 1, nblocks
            i0 = (ib - 1)*bsize
            nb = MIN(bsize, npoints - i0)
            np = INT(nb, KIND=C_SIZE_T)

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
               ! Meta-GGAs are screened on the kinetic energy density as well as on the
               ! density, and LibXC only screens on the latter. A point failing the tau
               ! test is handed over with zero density, which makes LibXC screen it out
               ! and leave its outputs at zero: what skipping it used to achieve. The
               ! Fermi hole curvature bound tau >= sigma/(8*rho) follows in the same pass,
               ! with the denominator floored at the cutoff so the division stays defined
               ! for the points that are about to be screened.
               DO i = 1, nb
                  IF (w%tau(1, i) <= eps_tau) w%rho(1, i) = 0.0_dp
                  w%tau(1, i) = MAX(w%tau(1, i), &
                                    w%sigma(1, i)/(8.0_dp*MAX(w%rho(1, i), eps_rho)))
               END DO
            END IF

            ! one LibXC call for the whole block
            IF (is_mgga) THEN
               SELECT CASE (grad_deriv)
               CASE (0)
                  CALL xc_f03_mgga_exc(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%lapl(1, 1), w%tau(1, 1), w%exc(1))
               CASE (1)
                  IF (no_exc) THEN
                     CALL xc_f03_mgga_vxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%lapl(1, 1), w%tau(1, 1), &
                                          w%vrho(1, 1), w%vsigma(1, 1), w%vlapl(1, 1), w%vtau(1, 1))
                     w%exc(1:nb) = 0.0_dp
                  ELSE
                     CALL xc_f03_mgga_exc_vxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%lapl(1, 1), w%tau(1, 1), &
                                              w%exc(1), w%vrho(1, 1), w%vsigma(1, 1), w%vlapl(1, 1), w%vtau(1, 1))
                  END IF
               CASE (2)
                  IF (no_exc) THEN
                     CALL xc_f03_mgga_vxc_fxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%lapl(1, 1), w%tau(1, 1), &
                                              w%vrho(1, 1), w%vsigma(1, 1), w%vlapl(1, 1), w%vtau(1, 1), &
                                              w%v2rho2(1, 1), w%v2rhosigma(1, 1), w%v2rholapl(1, 1), &
                                              w%v2rhotau(1, 1), w%v2sigma2(1, 1), w%v2sigmalapl(1, 1), &
                                              w%v2sigmatau(1, 1), w%v2lapl2(1, 1), w%v2lapltau(1, 1), w%v2tau2(1, 1))
                     w%exc(1:nb) = 0.0_dp
                  ELSE
                     CALL xc_f03_mgga(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%lapl(1, 1), w%tau(1, 1), &
                                      w%exc(1), w%vrho(1, 1), w%vsigma(1, 1), w%vlapl(1, 1), w%vtau(1, 1), &
                                      w%v2rho2(1, 1), w%v2rhosigma(1, 1), w%v2rholapl(1, 1), &
                                      w%v2rhotau(1, 1), w%v2sigma2(1, 1), w%v2sigmalapl(1, 1), &
                                      w%v2sigmatau(1, 1), w%v2lapl2(1, 1), w%v2lapltau(1, 1), w%v2tau2(1, 1))
                  END IF
               END SELECT
            ELSE IF (is_gga) THEN
               SELECT CASE (grad_deriv)
               CASE (0)
                  CALL xc_f03_gga_exc(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%exc(1))
               CASE (1)
                  IF (no_exc) THEN
                     CALL xc_f03_gga_vxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%vrho(1, 1), w%vsigma(1, 1))
                     w%exc(1:nb) = 0.0_dp
                  ELSE
                     CALL xc_f03_gga_exc_vxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), &
                                             w%exc(1), w%vrho(1, 1), w%vsigma(1, 1))
                  END IF
               CASE (2)
                  IF (no_exc) THEN
                     CALL xc_f03_gga_vxc_fxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), &
                                             w%vrho(1, 1), w%vsigma(1, 1), &
                                             w%v2rho2(1, 1), w%v2rhosigma(1, 1), w%v2sigma2(1, 1))
                     w%exc(1:nb) = 0.0_dp
                  ELSE
                     CALL xc_f03_gga_exc_vxc_fxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), &
                                                 w%exc(1), w%vrho(1, 1), w%vsigma(1, 1), &
                                                 w%v2rho2(1, 1), w%v2rhosigma(1, 1), w%v2sigma2(1, 1))
                  END IF
               END SELECT
            ELSE
               SELECT CASE (grad_deriv)
               CASE (0)
                  CALL xc_f03_lda_exc(w%func, np, w%rho(1, 1), w%exc(1))
               CASE (1)
                  CALL xc_f03_lda_exc_vxc(w%func, np, w%rho(1, 1), w%exc(1), w%vrho(1, 1))
               CASE (2)
                  CALL xc_f03_lda_exc_vxc_fxc(w%func, np, w%rho(1, 1), &
                                              w%exc(1), w%vrho(1, 1), w%v2rho2(1, 1))
               CASE (3)
                  CALL xc_f03_lda(w%func, np, w%rho(1, 1), &
                                  w%exc(1), w%vrho(1, 1), w%v2rho2(1, 1), w%v3rho3(1, 1))
               END SELECT
            END IF

            ! accumulate; each derivative is added by exactly one loop, selected by
            ! the order requested and the family, not by a per-branch copy of the whole
            ! block scaffolding
            IF (grad_deriv >= 0) THEN
               DO i = 1, nb
                  ii = i0 + i
                  e_0(ii) = e_0(ii) + sc*w%exc(i)*rho(ii)
               END DO
            END IF
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
               END IF
               IF (is_mgga) THEN
                  DO i = 1, nb
                     ii = i0 + i
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
!> \param grad_deriv degree of the derivative that should be evaluated;
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
      INTEGER, INTENT(in)                                :: grad_deriv, npoints
      CHARACTER(LEN=default_string_length), INTENT(IN)   :: func_name
      REAL(KIND=dp), INTENT(in)                          :: sc
      TYPE(libxc_worker_set_type), INTENT(INOUT)         :: workers

      INTEGER                                            :: bsize, family, i, i0, ib, ii, ithread, &
                                                            nb, nblocks, nthreads
      INTEGER(C_SIZE_T)                                  :: np
      LOGICAL                                            :: is_gga, is_mgga

      ithread = 0
      nthreads = 1
!$    ithread = omp_get_thread_num()
!$    nthreads = omp_get_num_threads()
      CPASSERT(ithread < SIZE(workers%worker))

      bsize = MAX(1, MIN(libxc_block_size, (npoints + nthreads - 1)/nthreads))
      nblocks = (npoints + bsize - 1)/bsize

      family = workers%family
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
         CPABORT(TRIM(func_name)//": this XC_FAMILY is currently not supported.")
      END SELECT

      ! As in the spin-unpolarized case the density and kinetic energy density
      ! cutoffs are applied by LibXC itself.
      ASSOCIATE (w => workers%worker(ithread + 1), has_laplace => workers%has_laplace, &
                 no_exc => workers%no_exc, eps_rho => workers%epsilon_rho, &
                 eps_tau => workers%epsilon_tau)
!$OMP    DO
         DO ib = 1, nblocks
            i0 = (ib - 1)*bsize
            nb = MIN(bsize, npoints - i0)
            np = INT(nb, KIND=C_SIZE_T)

            ! stage this block's inputs in LibXC's layout
            DO i = 1, nb
               ii = i0 + i
               w%rho(1, i) = MAX(rhoa(ii), 0.0_dp)
               w%rho(2, i) = MAX(rhob(ii), 0.0_dp)
            END DO
            IF (is_gga) THEN
               ! CP2K works with the norms of the gradients, LibXC with their contractions
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
               ! screening on tau and the Fermi hole curvature bound, per spin channel;
               ! see the spin-unpolarized routine. The tau test is on the sum of the two
               ! channels, as the density test is.
               DO i = 1, nb
                  IF (w%tau(1, i) + w%tau(2, i) <= eps_tau) THEN
                     w%rho(1, i) = 0.0_dp
                     w%rho(2, i) = 0.0_dp
                  END IF
                  w%tau(1, i) = MAX(w%tau(1, i), &
                                    w%sigma(1, i)/(8.0_dp*MAX(w%rho(1, i), eps_rho)))
                  w%tau(2, i) = MAX(w%tau(2, i), &
                                    w%sigma(3, i)/(8.0_dp*MAX(w%rho(2, i), eps_rho)))
               END DO
            END IF

            ! one LibXC call for the whole block
            IF (is_mgga) THEN
               SELECT CASE (grad_deriv)
               CASE (0)
                  CALL xc_f03_mgga_exc(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%lapl(1, 1), w%tau(1, 1), w%exc(1))
               CASE (1)
                  IF (no_exc) THEN
                     CALL xc_f03_mgga_vxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%lapl(1, 1), w%tau(1, 1), &
                                          w%vrho(1, 1), w%vsigma(1, 1), w%vlapl(1, 1), w%vtau(1, 1))
                     w%exc(1:nb) = 0.0_dp
                  ELSE
                     CALL xc_f03_mgga_exc_vxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%lapl(1, 1), w%tau(1, 1), &
                                              w%exc(1), w%vrho(1, 1), w%vsigma(1, 1), w%vlapl(1, 1), w%vtau(1, 1))
                  END IF
               CASE (2)
                  IF (no_exc) THEN
                     CALL xc_f03_mgga_vxc_fxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%lapl(1, 1), w%tau(1, 1), &
                                              w%vrho(1, 1), w%vsigma(1, 1), w%vlapl(1, 1), w%vtau(1, 1), &
                                              w%v2rho2(1, 1), w%v2rhosigma(1, 1), w%v2rholapl(1, 1), &
                                              w%v2rhotau(1, 1), w%v2sigma2(1, 1), w%v2sigmalapl(1, 1), &
                                              w%v2sigmatau(1, 1), w%v2lapl2(1, 1), w%v2lapltau(1, 1), w%v2tau2(1, 1))
                     w%exc(1:nb) = 0.0_dp
                  ELSE
                     CALL xc_f03_mgga(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%lapl(1, 1), w%tau(1, 1), &
                                      w%exc(1), w%vrho(1, 1), w%vsigma(1, 1), w%vlapl(1, 1), w%vtau(1, 1), &
                                      w%v2rho2(1, 1), w%v2rhosigma(1, 1), w%v2rholapl(1, 1), &
                                      w%v2rhotau(1, 1), w%v2sigma2(1, 1), w%v2sigmalapl(1, 1), &
                                      w%v2sigmatau(1, 1), w%v2lapl2(1, 1), w%v2lapltau(1, 1), w%v2tau2(1, 1))
                  END IF
               END SELECT
            ELSE IF (is_gga) THEN
               SELECT CASE (grad_deriv)
               CASE (0)
                  CALL xc_f03_gga_exc(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%exc(1))
               CASE (1)
                  IF (no_exc) THEN
                     CALL xc_f03_gga_vxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), w%vrho(1, 1), w%vsigma(1, 1))
                     w%exc(1:nb) = 0.0_dp
                  ELSE
                     CALL xc_f03_gga_exc_vxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), &
                                             w%exc(1), w%vrho(1, 1), w%vsigma(1, 1))
                  END IF
               CASE (2)
                  IF (no_exc) THEN
                     CALL xc_f03_gga_vxc_fxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), &
                                             w%vrho(1, 1), w%vsigma(1, 1), &
                                             w%v2rho2(1, 1), w%v2rhosigma(1, 1), w%v2sigma2(1, 1))
                     w%exc(1:nb) = 0.0_dp
                  ELSE
                     CALL xc_f03_gga_exc_vxc_fxc(w%func, np, w%rho(1, 1), w%sigma(1, 1), &
                                                 w%exc(1), w%vrho(1, 1), w%vsigma(1, 1), &
                                                 w%v2rho2(1, 1), w%v2rhosigma(1, 1), w%v2sigma2(1, 1))
                  END IF
               END SELECT
            ELSE
               SELECT CASE (grad_deriv)
               CASE (0)
                  CALL xc_f03_lda_exc(w%func, np, w%rho(1, 1), w%exc(1))
               CASE (1)
                  CALL xc_f03_lda_exc_vxc(w%func, np, w%rho(1, 1), w%exc(1), w%vrho(1, 1))
               CASE (2)
                  CALL xc_f03_lda_exc_vxc_fxc(w%func, np, w%rho(1, 1), &
                                              w%exc(1), w%vrho(1, 1), w%v2rho2(1, 1))
               CASE (3)
                  CALL xc_f03_lda(w%func, np, w%rho(1, 1), &
                                  w%exc(1), w%vrho(1, 1), w%v2rho2(1, 1), w%v3rho3(1, 1))
               END SELECT
            END IF

            ! accumulate; one loop per group of derivatives, selected by the order
            ! requested and the family
            IF (grad_deriv >= 0) THEN
               DO i = 1, nb
                  ii = i0 + i
                  e_0(ii) = e_0(ii) + sc*w%exc(i)*(w%rho(1, i) + w%rho(2, i))
               END DO
            END IF
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
         END DO
!$OMP    END DO
      END ASSOCIATE

   END SUBROUTINE libxc_spin_polarized_calc
#endif

END MODULE xc_libxc
