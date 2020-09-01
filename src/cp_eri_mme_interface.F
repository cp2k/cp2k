!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to Minimax-Ewald method for periodic ERI's to be used in CP2K.
! **************************************************************************************************

MODULE cp_eri_mme_interface
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_create,&
                                              cell_release,&
                                              cell_type,&
                                              init_cell,&
                                              pbc
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_section_create,&
                                              cp_print_key_unit_nr,&
                                              medium_print_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE eri_mme_test,                    ONLY: eri_mme_2c_perf_acc_test,&
                                              eri_mme_3c_perf_acc_test
   USE eri_mme_types,                   ONLY: eri_mme_coulomb,&
                                              eri_mme_init,&
                                              eri_mme_longrange,&
                                              eri_mme_param,&
                                              eri_mme_print_grid_info,&
                                              eri_mme_release,&
                                              eri_mme_set_params,&
                                              eri_mme_yukawa
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE input_val_types,                 ONLY: real_t
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_sum
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_eri_mme_interface'

   PUBLIC :: &
      cp_eri_mme_finalize, &
      cp_eri_mme_init_read_input, &
      cp_eri_mme_param, &
      cp_eri_mme_perf_acc_test, &
      cp_eri_mme_set_params, &
      cp_eri_mme_update_local_counts, &
      create_eri_mme_section, &
      create_eri_mme_test_section

   INTERFACE cp_eri_mme_set_params
      MODULE PROCEDURE eri_mme_set_params_from_basis
      MODULE PROCEDURE eri_mme_set_params_custom
   END INTERFACE

   TYPE cp_eri_mme_param
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(eri_mme_param)                                :: par
      TYPE(section_vals_type), &
         POINTER                                         :: mme_section => NULL()
      INTEGER                                            :: G_count_2c, R_count_2c
      INTEGER                                            :: GG_count_3c, GR_count_3c, RR_count_3c
      LOGICAL                                            :: do_calib
   END TYPE cp_eri_mme_param

CONTAINS
! **************************************************************************************************
!> \brief Create main input section
!> \param section ...
!> \param default_n_minimax ...
! **************************************************************************************************
   SUBROUTINE create_eri_mme_section(section, default_n_minimax)
      TYPE(section_type), POINTER                        :: section
      INTEGER, INTENT(IN), OPTIONAL                      :: default_n_minimax

      INTEGER                                            :: my_default_n_minimax
      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      NULLIFY (keyword, print_key, subsection)
      CPASSERT(.NOT. ASSOCIATED(section))

      IF (PRESENT(default_n_minimax)) THEN
         my_default_n_minimax = default_n_minimax
      ELSE
         my_default_n_minimax = 20
      ENDIF

      CALL section_create(section, __LOCATION__, name="ERI_MME", &
                          description="Parameters for the calculation of periodic electron repulsion "// &
                          "integrals (ERI) using the Minimax-Ewald (MME) method. "// &
                          "Note: N_MINIMAX is the only parameter to be tuned for accuracy, "// &
                          "all other parameters can be left to default. MME method is faster "// &
                          "than numerical GPW.", &
                          n_keywords=5, n_subsections=1)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="N_MINIMAX", &
                          description="Number of terms in minimax approximation of "// &
                          "reciprocal space potential. ", &
                          default_i_val=my_default_n_minimax)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="CUTOFF", &
                          description="User-defined energy cutoff to be used only if "// &
                          "DO_CALIBRATE_CUTOFF is set to .FALSE. ", &
                          default_r_val=300.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="SUM_PRECISION", &
                          description="Terms in lattice sums are ignored if absolute value smaller than this value.", &
                          default_r_val=1.0E-16_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DO_CALIBRATE_CUTOFF", &
                          description="Whether the energy cutoff shall be calibrated to "// &
                          "minimize upper bound error estimate. ", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DO_ERROR_ESTIMATE", &
                          description="Whether the error due to minimax approx. and cutoff shall be estimated", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "ERI_MME_INFO", &
                                       description="Controls the printing info.", &
                                       print_level=medium_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="PRINT_CALIB", &
                          description="Print detailed info on calibration. ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DEBUG", &
                          description="debug mode (consistency of summation methods is checked).", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DEBUG_TOLERANCE", &
                          description="tolerance for rel. numerical error in debug mode.", &
                          default_r_val=1.0E-06_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DEBUG_NSUM_MAX", &
                          description="restrict debug mode for non-ortho cells to this number of summands. "// &
                          "Sums with more terms are not checked.", &
                          default_i_val=1000000)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, __LOCATION__, name="CUTOFF_CALIB", &
                          description="Parameters for the calibration of the energy cutoff by "// &
                          "minimizing the errors due to finite cutoff and minimax approximation. "// &
                          "Implemented as bisection of error(minimax) - error(cutoff). Not "// &
                          "implemented for non-orthorhombic cells. ", &
                          n_keywords=5, n_subsections=0)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="MIN", &
                          description="Initial guess of lower bound for cutoff. ", &
                          default_r_val=10.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="MAX", &
                          description="Initial guess of upper bound for cutoff. ", &
                          default_r_val=10000.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="DELTA", &
                          description="Relative widening of cutoff interval in case starting "// &
                          "values are not valid. ", &
                          default_r_val=0.9_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="EPS", &
                          description="Relative cutoff precision required to stop calibration. ", &
                          default_r_val=0.01_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_eri_mme_section

! **************************************************************************************************
!> \brief Read input and initialize parameter type
!> \param mme_section ...
!> \param param ...
! **************************************************************************************************
   SUBROUTINE cp_eri_mme_init_read_input(mme_section, param)
      TYPE(section_vals_type), POINTER                   :: mme_section
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param

      INTEGER                                            :: debug_nsum, n_minimax, unit_nr
      LOGICAL                                            :: debug, do_calib_cutoff, do_error_est, &
                                                            print_calib
      REAL(KIND=dp)                                      :: cutoff, cutoff_delta, cutoff_eps, &
                                                            cutoff_max, cutoff_min, debug_delta, &
                                                            sum_precision
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: subsection

      logger => cp_get_default_logger()
      unit_nr = cp_print_key_unit_nr(logger, mme_section, "ERI_MME_INFO", &
                                     extension=".eri_mme")

      NULLIFY (subsection)
      CALL section_vals_val_get(mme_section, "N_MINIMAX", i_val=n_minimax)

      CALL section_vals_val_get(mme_section, "CUTOFF", r_val=cutoff)
      CALL section_vals_val_get(mme_section, "SUM_PRECISION", r_val=sum_precision)
      CALL section_vals_val_get(mme_section, "DO_CALIBRATE_CUTOFF", l_val=do_calib_cutoff)
      CALL section_vals_val_get(mme_section, "DO_ERROR_ESTIMATE", l_val=do_error_est)
      CALL section_vals_val_get(mme_section, "PRINT_CALIB", l_val=print_calib)
      subsection => section_vals_get_subs_vals(mme_section, "CUTOFF_CALIB")
      CALL section_vals_val_get(subsection, "MIN", r_val=cutoff_min)
      CALL section_vals_val_get(subsection, "MAX", r_val=cutoff_max)
      CALL section_vals_val_get(subsection, "EPS", r_val=cutoff_eps)
      CALL section_vals_val_get(subsection, "DELTA", r_val=cutoff_delta)
      CALL section_vals_val_get(mme_section, "DEBUG", l_val=debug)
      CALL section_vals_val_get(mme_section, "DEBUG_TOLERANCE", r_val=debug_delta)
      CALL section_vals_val_get(mme_section, "DEBUG_NSUM_MAX", i_val=debug_nsum)

      param%mme_section => mme_section

      CALL eri_mme_init(param%par, n_minimax, &
                        cutoff, do_calib_cutoff, do_error_est, cutoff_min, cutoff_max, cutoff_eps, cutoff_delta, &
                        sum_precision, debug, debug_delta, debug_nsum, unit_nr, print_calib)

      param%do_calib = do_calib_cutoff

      param%G_count_2c = 0
      param%R_count_2c = 0
      param%GG_count_3c = 0
      param%GR_count_3c = 0
      param%RR_count_3c = 0

      param%logger => logger
   END SUBROUTINE cp_eri_mme_init_read_input

! **************************************************************************************************
!> \brief Release eri mme data. Prints some statistics on summation methods chosen.
!> \param param ...
! **************************************************************************************************
   SUBROUTINE cp_eri_mme_finalize(param)
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param

      INTEGER                                            :: count_2c, count_3c, unit_nr

      count_2c = param%G_count_2c + param%R_count_2c
      count_3c = param%GG_count_3c + param%GR_count_3c + param%RR_count_3c

      unit_nr = param%par%unit_nr

      IF (unit_nr > 0) THEN
         IF (count_2c .GT. 0) THEN
            WRITE (unit_nr, '(/T2, A)') "ERI_MME| Percentage of 2-center integrals evaluated in"
            WRITE (unit_nr, '(T2, A, T76, F5.1)') "ERI_MME|   G space:", &
               100.0_dp*param%G_count_2c/count_2c
            WRITE (unit_nr, '(T2, A, T76, F5.1/)') "ERI_MME|   R space:", &
               100.0_dp*param%R_count_2c/count_2c
         ENDIF
         IF (count_3c .GT. 0) THEN
            WRITE (unit_nr, '(/T2, A)') "ERI_MME| Percentage of 3-center integrals evaluated in"
            WRITE (unit_nr, '(T2, A, T76, F5.1)') "ERI_MME|   G/G space:", &
               100.0_dp*param%GG_count_3c/count_3c
            WRITE (unit_nr, '(T2, A, T76, F5.1)') "ERI_MME|   G/R space:", &
               100.0_dp*param%GR_count_3c/count_3c
            WRITE (unit_nr, '(T2, A, T76, F5.1/)') "ERI_MME|   R/R space:", &
               100.0_dp*param%RR_count_3c/count_3c
         ENDIF
      ENDIF

      CALL eri_mme_release(param%par)
      CALL cp_print_key_finished_output(unit_nr, param%logger, param%mme_section, "ERI_MME_INFO")
   END SUBROUTINE cp_eri_mme_finalize

! **************************************************************************************************
!> \brief Set parameters for MME method by deriving basis info from basis set.
!>        Cutoff can be auto-calibrated to minimize total error.
!> \param param ...
!> \param cell ...
!> \param qs_kind_set ...
!> \param basis_type_1 ...
!> \param basis_type_2 ...
!> \param para_env ...
!> \param potential ...
!> \param pot_par ...
! **************************************************************************************************
   SUBROUTINE eri_mme_set_params_from_basis(param, cell, qs_kind_set, basis_type_1, basis_type_2, para_env, &
                                            potential, pot_par)
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param
      TYPE(cell_type), INTENT(IN)                        :: cell
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      CHARACTER(len=*), INTENT(IN)                       :: basis_type_1
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type_2
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(IN), OPTIONAL                      :: potential
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: pot_par

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_set_params_from_basis'

      INTEGER                                            :: handle, l_max, l_max_zet
      REAL(KIND=dp)                                      :: zet_max, zet_min

      CALL timeset(routineN, handle)

      CALL error_est_pgf_params_from_basis(qs_kind_set, basis_type_1, basis_type_2, &
                                           zet_min, zet_max, l_max_zet, l_max)

      CALL eri_mme_set_params_custom(param, cell%hmat, cell%orthorhombic, &
                                     zet_min, zet_max, l_max_zet, &
                                     l_max, para_env, &
                                     potential, pot_par)

      CALL timestop(handle)
   END SUBROUTINE eri_mme_set_params_from_basis

! **************************************************************************************************
!> \brief Wrapper for eri_mme_set_params
!> \param param ...
!> \param hmat ...
!> \param is_ortho ...
!> \param zet_min ...
!> \param zet_max ...
!> \param l_max_zet ...
!> \param l_max ...
!> \param para_env ...
!> \param potential ...
!> \param pot_par ...
! **************************************************************************************************
   SUBROUTINE eri_mme_set_params_custom(param, hmat, is_ortho, zet_min, zet_max, l_max_zet, l_max, para_env, &
                                        potential, pot_par)
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat
      LOGICAL, INTENT(IN)                                :: is_ortho
      REAL(KIND=dp), INTENT(IN)                          :: zet_min, zet_max
      INTEGER, INTENT(IN)                                :: l_max_zet, l_max
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(IN), OPTIONAL                      :: potential
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: pot_par

      REAL(KIND=dp), PARAMETER                           :: eps_changed = 1.0E-14_dp

      IF (param%do_calib) THEN
         IF (.NOT. param%par%is_valid) THEN
            param%par%do_calib_cutoff = .TRUE.
         ELSE
            ! only calibrate cutoff if parameters (cell, basis coefficients) have changed
            IF (ALL(ABS(param%par%hmat - hmat) < eps_changed) .AND. &
                ABS(param%par%zet_min - zet_min) < eps_changed .AND. &
                ABS(param%par%zet_max - zet_max) < eps_changed .AND. &
                param%par%l_max_zet == l_max_zet) THEN
               param%par%do_calib_cutoff = .FALSE.
            ELSE
               param%par%do_calib_cutoff = .TRUE.
            ENDIF
         ENDIF
      ELSE
         param%par%do_calib_cutoff = .FALSE.
      ENDIF

      CALL eri_mme_set_params(param%par, hmat, is_ortho, zet_min, zet_max, l_max_zet, l_max, para_env, &
                              potential, pot_par)

      CALL eri_mme_print_info(param)
   END SUBROUTINE eri_mme_set_params_custom

! **************************************************************************************************
!> \brief Get basis parameters for estimating cutoff and minimax error from cp2k basis
!> \param qs_kind_set ...
!> \param basis_type_1 ...
!> \param basis_type_2 ...
!> \param zet_min    Smallest exponent, used to estimate error due to minimax approx.
!> \param zet_max    contains max. exponent,
!>                   used to estimate cutoff error
!> \param l_max_zet  contains the largest l for max. exponent,
!>                   used to estimate cutoff error
!> \param l_max ...
! **************************************************************************************************
   SUBROUTINE error_est_pgf_params_from_basis(qs_kind_set, basis_type_1, basis_type_2, zet_min, zet_max, l_max_zet, l_max)
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      CHARACTER(len=*), INTENT(IN)                       :: basis_type_1
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type_2
      REAL(KIND=dp), INTENT(OUT)                         :: zet_min, zet_max
      INTEGER, INTENT(OUT)                               :: l_max_zet, l_max

      CHARACTER(LEN=*), PARAMETER :: routineN = 'error_est_pgf_params_from_basis'

      CHARACTER(len=default_string_length)               :: basis_type
      INTEGER                                            :: handle, ibasis, ikind, ipgf, iset, l_m, &
                                                            l_zet, nbasis, nkind, nset
      INTEGER, DIMENSION(:), POINTER                     :: npgf
      REAL(KIND=dp)                                      :: zet_m
      TYPE(gto_basis_set_type), POINTER                  :: basis_set

      CALL timeset(routineN, handle)

      l_m = 0
      zet_m = 0.0_dp
      l_zet = -1
      zet_min = -1.0_dp

      nkind = SIZE(qs_kind_set)
      nbasis = MERGE(2, 1, PRESENT(basis_type_2))

      ! 1) get global max l and max zet
      ! (and min zet for minimax error)
      DO ikind = 1, nkind
         DO ibasis = 1, nbasis
            IF (ibasis .EQ. 1) THEN
               basis_type = basis_type_1
            ELSE
               basis_type = basis_type_2
            ENDIF
            CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set, &
                             basis_type=basis_type)
            CPASSERT(ASSOCIATED(basis_set))
            npgf => basis_set%npgf
            nset = basis_set%nset
            l_m = MAX(l_m, MAXVAL(basis_set%lmax(:)))
            DO iset = 1, nset
               zet_m = MAX(zet_m, MAXVAL(basis_set%zet(1:npgf(iset), iset)))
               IF (zet_min .LT. 0.0_dp) THEN
                  zet_min = MINVAL(basis_set%zet(1:npgf(iset), iset))
               ELSE
                  zet_min = MIN(zet_min, MINVAL(basis_set%zet(1:npgf(iset), iset)))
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ! 2) get largest zet for max l and largest l for max zet
      DO ikind = 1, nkind
         DO ibasis = 1, nbasis
            IF (ibasis .EQ. 1) THEN
               basis_type = basis_type_1
            ELSE
               basis_type = basis_type_2
            ENDIF
            CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set, &
                             basis_type=basis_type)
            DO iset = 1, basis_set%nset
               DO ipgf = 1, basis_set%npgf(iset)
                  IF (ABS(zet_m - basis_set%zet(ipgf, iset)) .LE. (zet_m*1.0E-12_dp) &
                      .AND. (basis_set%lmax(iset) .GT. l_zet)) THEN
                     l_zet = basis_set%lmax(iset)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CPASSERT(l_zet .GE. 0)

      zet_max = zet_m

      ! l + 1 because we may calculate forces
      ! this is probably a safe choice also for the case that forces are not needed
      l_max_zet = l_zet + 1
      l_max = l_m + 1

      CALL timestop(handle)
   END SUBROUTINE error_est_pgf_params_from_basis

! **************************************************************************************************
!> \brief ...
!> \param param ...
! **************************************************************************************************
   SUBROUTINE eri_mme_print_info(param)
      TYPE(cp_eri_mme_param)                             :: param

      INTEGER                                            :: igrid, unit_nr
      LOGICAL                                            :: print_multigrids
      TYPE(cp_logger_type), POINTER                      :: logger

      print_multigrids = .FALSE.

      logger => param%logger

      unit_nr = param%par%unit_nr

      IF (unit_nr > 0) THEN
         SELECT CASE (param%par%potential)
         CASE (eri_mme_coulomb)
            WRITE (unit_nr, '(/T2, A)') "ERI_MME| Potential: Coulomb"
         CASE (eri_mme_yukawa)
            WRITE (unit_nr, '(/T2, A, ES9.2)') "ERI_MME| Potential: Yukawa with a=", param%par%pot_par
         CASE (eri_mme_longrange)
            WRITE (unit_nr, '(/T2, A, ES9.2)') "ERI_MME| Potential: long-range Coulomb with a=", param%par%pot_par
         END SELECT
      ENDIF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/T2, A, T71, F10.1)') "ERI_MME| Cutoff for ERIs [a.u.]:", param%par%cutoff
         WRITE (unit_nr, '(/T2, A, T78, I3/)') "ERI_MME| Number of terms in minimax approximation:", param%par%n_minimax
      ENDIF
      IF (param%par%is_ortho) THEN
         IF (unit_nr > 0) THEN
            IF (param%par%do_error_est) THEN
               WRITE (unit_nr, '(T2, A)') "ERI_MME| Estimated absolute error for normalized Hermite-Gaussian basis"
               WRITE (unit_nr, '(T2, A, T72, ES9.2)') "ERI_MME|   Minimax error:", param%par%err_mm
               WRITE (unit_nr, '(T2, A, T72, ES9.2)') "ERI_MME|   Cutoff error:", param%par%err_c
             WRITE (unit_nr, '(T2, A, T72, ES9.2)') "ERI_MME|   Total error (minimax + cutoff):", param%par%err_mm + param%par%err_c
            ENDIF
            IF (param%par%print_calib) &
               WRITE (unit_nr, '(T2, A, T68, F13.10)') "ERI_MME|   Minimax scaling constant in AM-GM estimate:", param%par%C_mm
         ENDIF
      ENDIF

      IF (print_multigrids) THEN
         DO igrid = 1, param%par%n_grids
            CALL eri_mme_print_grid_info(param%par%minimax_grid(igrid), igrid, unit_nr)
         ENDDO
      ENDIF

      IF (unit_nr > 0) WRITE (unit_nr, *)

   END SUBROUTINE eri_mme_print_info

! **************************************************************************************************
!> \brief Create input section for unit testing
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_eri_mme_test_section(section)
      TYPE(section_type), INTENT(INOUT), POINTER         :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      NULLIFY (keyword, subsection)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="ERI_MME_TEST", &
                          description="Parameters for testing the ERI_MME method for electron repulsion integrals. "// &
                          "Testing w.r.t. performance and accuracy. ", &
                          n_keywords=5, n_subsections=1)

      CALL create_eri_mme_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, __LOCATION__, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation the ERI_MME test. ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TEST_3C", &
                          description="Whether to test 3-center integrals.", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TEST_2C", &
                          description="Whether to test 2-center integrals.", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ABC", &
                          description="Specify the lengths of the cell vectors A, B, and C. ", &
                          usage="ABC 10.000 10.000 10.000", unit_str="angstrom", &
                          n_var=3, type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MIN_NPOS", &
                          description="Minimum number of atomic distances to consider. ", &
                          default_i_val=8)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NREP", &
                          description="Number of repeated calculation of each integral. ", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="CHECK_2C_ACCURACY", &
                          description="Whether integrals should be compared to reference values, "// &
                          "created on the fly by exact method (G-space sum on grid without "// &
                          "minimax approximation). Note: only feasible for not too many "// &
                          "integrals and maximum exponent around 10.0. ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)

      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LMAX", &
                          description="Maximum total angular momentum quantum number. ", &
                          default_i_val=6)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ZET_MIN", &
                          description="Minimum exponent. ", &
                          default_r_val=0.001_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ZET_MAX", &
                          description="Maximum exponent. ", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NZET", &
                          description="Number of exponents (logarithmic partition between ZET_MIN and ZET_MAX). ", &
                          default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NSAMPLE_3C", &
                          description="If NSAMPLE_3C = N, only calculate every Nth 3-center integral.", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POTENTIAL", &
                          description="Operator to test", &
                          default_i_val=eri_mme_coulomb, &
                          enum_i_vals=(/eri_mme_coulomb, eri_mme_yukawa, eri_mme_longrange/), &
                          enum_c_vals=s2a("COULOMB", "YUKAWA", "LONGRANGE"), &
                          enum_desc=s2a("1/r", "exp(-a*r)/r", "erf(a*r)/r"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POTENTIAL_PARAM", &
                          description="Parameter 'a' for chosen potential, ignored for Coulomb", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_eri_mme_test_section

! **************************************************************************************************
!> \brief Update local counters to gather statistics on different paths taken in MME algorithm
!> (each Ewald sum can be performed over direct or reciprocal lattice vectors)
!> \param param ...
!> \param para_env ...
!> \param G_count_2c ...
!> \param R_count_2c ...
!> \param GG_count_3c ...
!> \param GR_count_3c ...
!> \param RR_count_3c ...
! **************************************************************************************************
   SUBROUTINE cp_eri_mme_update_local_counts(param, para_env, G_count_2c, R_count_2c, GG_count_3c, GR_count_3c, RR_count_3c)
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(INOUT), OPTIONAL                   :: G_count_2c, R_count_2c, GG_count_3c, &
                                                            GR_count_3c, RR_count_3c

      IF (PRESENT(G_count_2c)) THEN
         CALL mp_sum(G_count_2c, para_env%group)
         param%G_count_2c = param%G_count_2c + G_count_2c
      ENDIF

      IF (PRESENT(R_count_2c)) THEN
         CALL mp_sum(R_count_2c, para_env%group)
         param%R_count_2c = param%R_count_2c + R_count_2c
      ENDIF

      IF (PRESENT(GG_count_3c)) THEN
         CALL mp_sum(GG_count_3c, para_env%group)
         param%GG_count_3c = param%GG_count_3c + GG_count_3c
      ENDIF

      IF (PRESENT(GR_count_3c)) THEN
         CALL mp_sum(GR_count_3c, para_env%group)
         param%GR_count_3c = param%GR_count_3c + GR_count_3c
      ENDIF

      IF (PRESENT(RR_count_3c)) THEN
         CALL mp_sum(RR_count_3c, para_env%group)
         param%RR_count_3c = param%RR_count_3c + RR_count_3c
      ENDIF

   END SUBROUTINE cp_eri_mme_update_local_counts

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param iw ...
!> \param eri_mme_test_section ...
! **************************************************************************************************
   SUBROUTINE cp_eri_mme_perf_acc_test(para_env, iw, eri_mme_test_section)
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(IN)                                :: iw
      TYPE(section_vals_type), POINTER                   :: eri_mme_test_section

      INTEGER                                            :: count_r, G_count, GG_count, GR_count, i, &
                                                            ix, iy, iz, l_max, min_nR, nR, nR_xyz, &
                                                            nrep, nsample, nzet, potential, &
                                                            R_count, RR_count
      LOGICAL                                            :: test_2c, test_3c, test_accuracy
      REAL(KIND=dp)                                      :: pot_par, zet_fac, zetmax, zetmin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rabc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cell_par
      TYPE(cell_type), POINTER                           :: box
      TYPE(cp_eri_mme_param)                             :: param
      TYPE(section_vals_type), POINTER                   :: eri_mme_section

      NULLIFY (box, eri_mme_section, cell_par)

      eri_mme_section => section_vals_get_subs_vals(eri_mme_test_section, "ERI_MME")
      CALL cp_eri_mme_init_read_input(eri_mme_section, param)
      CALL section_vals_val_get(eri_mme_test_section, "TEST_3C", l_val=test_3c)
      CALL section_vals_val_get(eri_mme_test_section, "TEST_2C", l_val=test_2c)

      CALL section_vals_val_get(eri_mme_test_section, "ABC", r_vals=cell_par)
      CALL section_vals_val_get(eri_mme_test_section, "MIN_NPOS", i_val=min_nR)
      CALL section_vals_val_get(eri_mme_test_section, "NREP", i_val=nrep)
      CALL section_vals_val_get(eri_mme_test_section, "CHECK_2C_ACCURACY", l_val=test_accuracy)
      CALL section_vals_val_get(eri_mme_test_section, "LMAX", i_val=l_max)
      CALL section_vals_val_get(eri_mme_test_section, "ZET_MIN", r_val=zetmin)
      CALL section_vals_val_get(eri_mme_test_section, "ZET_MAX", r_val=zetmax)
      CALL section_vals_val_get(eri_mme_test_section, "NZET", i_val=nzet)
      CALL section_vals_val_get(eri_mme_test_section, "NSAMPLE_3C", i_val=nsample)
      CALL section_vals_val_get(eri_mme_test_section, "POTENTIAL", i_val=potential)
      CALL section_vals_val_get(eri_mme_test_section, "POTENTIAL_PARAM", r_val=pot_par)

      IF (nzet .LE. 0) &
         CPABORT("Number of exponents NZET must be greater than 0.")

      CALL init_orbital_pointers(l_max)

      ! Create ranges of zet to be tested
      ALLOCATE (zet(nzet))

      zet(1) = zetmin
      IF (nzet .GT. 1) THEN
         zet_fac = (zetmax/zetmin)**(1.0_dp/(nzet - 1))
         DO i = 1, nzet - 1
            zet(i + 1) = zet(i)*zet_fac
         ENDDO
      ENDIF

      ! initialize cell
      CALL cell_create(box)
      box%hmat = 0.0_dp
      box%hmat(1, 1) = cell_par(1)
      box%hmat(2, 2) = cell_par(2)
      box%hmat(3, 3) = cell_par(3)
      CALL init_cell(box)

      ! Create range of rab (atomic distances) to be tested
      nR_xyz = CEILING(REAL(min_nR, KIND=dp)**(1.0_dp/3.0_dp) - 1.0E-06)
      nR = nR_xyz**3

      ALLOCATE (rabc(3, nR))
      count_r = 0
      DO ix = 1, nR_xyz
      DO iy = 1, nR_xyz
      DO iz = 1, nR_xyz
         count_r = count_r + 1
         ! adding 10% of cell size to positions to avoid atoms exactly at boundary or center of a cell
         rabc(:, count_r) = pbc([ix*ABS(cell_par(1)), &
                                 iy*ABS(cell_par(2)), &
                                 iz*ABS(cell_par(3))]/nR_xyz + &
                                0.1_dp*ABS(cell_par(:)), box)
      ENDDO
      ENDDO
      ENDDO

      ! initialize MME method
      CALL cp_eri_mme_set_params(param, box%hmat, box%orthorhombic, MINVAL(zet), MAXVAL(zet), l_max, l_max, para_env, &
                                 potential, pot_par)

      IF (iw > 0) WRITE (iw, '(T2, A, T61, I20)') "ERI_MME| Number of atomic distances:", nR

      G_count = 0; R_count = 0
      GG_count = 0; GR_count = 0; RR_count = 0

      IF (test_2c) CALL eri_mme_2c_perf_acc_test(param%par, l_max, zet, rabc, nrep, test_accuracy, para_env, iw, &
                                                 potential=potential, pot_par=pot_par, G_count=G_count, R_count=R_count)
      IF (test_3c) CALL eri_mme_3c_perf_acc_test(param%par, l_max, zet, rabc, nrep, nsample, &
                                                 para_env, iw, potential=potential, pot_par=pot_par, &
                                                 GG_count=GG_count, GR_count=GR_count, RR_count=RR_count)
      CALL cp_eri_mme_update_local_counts(param, para_env, G_count, R_count, GG_count, GR_count, RR_count)
      CALL cp_eri_mme_finalize(param)
      CALL cell_release(box)
   END SUBROUTINE cp_eri_mme_perf_acc_test

END MODULE cp_eri_mme_interface
