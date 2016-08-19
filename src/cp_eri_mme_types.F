!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to Minimax-Ewald method for periodic ERI's to be used in CP2K.
!>        Types and initialization / finalization methods.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE cp_eri_mme_types
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_section_create,&
                                              cp_print_key_unit_nr,&
                                              medium_print_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE eri_mme_types,                   ONLY: eri_mme_init,&
                                              eri_mme_param,&
                                              eri_mme_release,&
                                              eri_mme_set_params
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
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_eri_mme_types'

   PUBLIC :: cp_eri_mme_param, &
             cp_eri_mme_set_params, &
             cp_eri_mme_init_read_input, &
             cp_eri_mme_finalize, &
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
   END TYPE cp_eri_mme_param

CONTAINS

! **************************************************************************************************
!> \brief Create main input section
!> \param section ...
!> \param default_n_minimax ...
! **************************************************************************************************
   SUBROUTINE create_eri_mme_section(section, default_n_minimax)
      TYPE(section_type), INTENT(INOUT), POINTER         :: section
      INTEGER, INTENT(IN), OPTIONAL                      :: default_n_minimax

      CHARACTER(len=*), PARAMETER :: routineN = 'create_eri_mme_section', &
         routineP = moduleN//':'//routineN

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

      CALL section_create(section, "ERI_MME", &
                          description="Parameters for the calculation of periodic electron repulsion "// &
                          "integrals (ERI) using the Minimax-Ewald (MME) method. "// &
                          "Note: N_MINIMAX is the only parameter to be tuned for accuracy, "// &
                          "all other parameters can be left to default. MME method is faster "// &
                          "than numerical GPW.", &
                          n_keywords=5, n_subsections=1)

      CALL keyword_create(keyword, &
                          name="N_MINIMAX", &
                          description="Number of terms in minimax approximation of "// &
                          "reciprocal space potential. ", &
                          default_i_val=my_default_n_minimax)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, &
                          name="CUTOFF", &
                          description="User-defined energy cutoff to be used only if "// &
                          "DO_CALIBRATE_CUTOFF is set to .FALSE. ", &
                          default_r_val=300.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, &
                          name="DO_CALIBRATE_CUTOFF", &
                          description="Whether the energy cutoff shall be calibrated to "// &
                          "minimize upper bound error estimate. ", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL cp_print_key_section_create(print_key, "ERI_MME_INFO", &
                                       description="Controls the printing info.", &
                                       print_level=medium_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL keyword_create(keyword, &
                          name="PRINT_CALIB", &
                          description="Print detailed info on calibration. ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, &
                          name="DEBUG", &
                          description="debug mode (consistency of summation methods is checked).", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, &
                          name="DEBUG_TOLERANCE", &
                          description="tolerance for rel. numerical error in debug mode.", &
                          default_r_val=1.0E-06_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL section_create(subsection, "CUTOFF_CALIB", &
                          description="Parameters for the calibration of the energy cutoff by "// &
                          "minimizing the errors due to finite cutoff and minimax approximation. "// &
                          "Implemented as bisection of error(minimax) - error(cutoff). Not "// &
                          "implemented for non-orthorhombic cells. ", &
                          n_keywords=5, n_subsections=0)

      CALL keyword_create(keyword, &
                          name="MIN", &
                          description="Initial guess of lower bound for cutoff. ", &
                          default_r_val=10.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, &
                          name="MAX", &
                          description="Initial guess of upper bound for cutoff. ", &
                          default_r_val=10000.0_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, &
                          name="DELTA", &
                          description="Relative widening of cutoff interval in case starting "// &
                          "values are not valid. ", &
                          default_r_val=0.9_dp)
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, &
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
      TYPE(section_vals_type), INTENT(INOUT), POINTER    :: mme_section
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_eri_mme_init_read_input', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: n_minimax, unit_nr
      LOGICAL                                            :: debug, do_calib_cutoff, print_calib
      REAL(KIND=dp)                                      :: cutoff, cutoff_delta, cutoff_eps, &
                                                            cutoff_max, cutoff_min, debug_delta
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: subsection

      logger => cp_get_default_logger()
      unit_nr = cp_print_key_unit_nr(logger, mme_section, "ERI_MME_INFO", &
                                     extension=".eri_mme")

      NULLIFY (subsection)
      CALL section_vals_val_get(mme_section, "N_MINIMAX", i_val=n_minimax)

      CALL section_vals_val_get(mme_section, "CUTOFF", r_val=cutoff)
      CALL section_vals_val_get(mme_section, "DO_CALIBRATE_CUTOFF", l_val=do_calib_cutoff)
      CALL section_vals_val_get(mme_section, "PRINT_CALIB", l_val=print_calib)
      subsection => section_vals_get_subs_vals(mme_section, "CUTOFF_CALIB")
      CALL section_vals_val_get(subsection, "MIN", r_val=cutoff_min)
      CALL section_vals_val_get(subsection, "MAX", r_val=cutoff_max)
      CALL section_vals_val_get(subsection, "EPS", r_val=cutoff_eps)
      CALL section_vals_val_get(subsection, "DELTA", r_val=cutoff_delta)
      CALL section_vals_val_get(mme_section, "DEBUG", l_val=debug)
      CALL section_vals_val_get(mme_section, "DEBUG_TOLERANCE", r_val=debug_delta)

      param%mme_section => mme_section

      CALL eri_mme_init(param%par, n_minimax, &
                        cutoff, do_calib_cutoff, cutoff_min, cutoff_max, cutoff_eps, cutoff_delta, &
                        debug, debug_delta, unit_nr, print_calib)

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_eri_mme_finalize', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: count_2c, count_3c, unit_nr

      count_2c = param%G_count_2c+param%R_count_2c
      count_3c = param%GG_count_3c+param%GR_count_3c+param%RR_count_3c

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
!>        Switch between G / R space representation can be auto-calibrated to
!>        minimize run time.
!> \param param ...
!> \param cell ...
!> \param qs_kind_set ...
!> \param basis_type_1 ...
!> \param basis_type_2 ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE eri_mme_set_params_from_basis(param, cell, qs_kind_set, basis_type_1, basis_type_2, para_env)
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param
      TYPE(cell_type), INTENT(IN)                        :: cell
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      CHARACTER(len=*), INTENT(IN)                       :: basis_type_1
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type_2
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_set_params_from_basis', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, l_max
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: l_err_cutoff
      REAL(KIND=dp)                                      :: zet_err_minimax
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: zet_err_cutoff

      CALL timeset(routineN, handle)

      CALL error_est_pgf_params_from_basis(qs_kind_set, basis_type_1, basis_type_2, &
                                           zet_err_minimax, zet_err_cutoff, l_err_cutoff)

      l_max = MAXVAL(l_err_cutoff)
      CALL eri_mme_set_params_custom(param, cell%hmat, cell%orthorhombic, &
                                     zet_err_minimax, zet_err_cutoff, l_err_cutoff, &
                                     l_max, para_env)

      CALL timestop(handle)
   END SUBROUTINE eri_mme_set_params_from_basis

! **************************************************************************************************
!> \brief Wrapper for eri_mme_set_params
!> \param param ...
!> \param hmat ...
!> \param is_ortho ...
!> \param zet_err_minimax ...
!> \param zet_err_cutoff ...
!> \param l_err_cutoff ...
!> \param l_max ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE eri_mme_set_params_custom(param, hmat, is_ortho, zet_err_minimax, zet_err_cutoff, l_err_cutoff, l_max, para_env)
      TYPE(cp_eri_mme_param), INTENT(INOUT)              :: param
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat
      LOGICAL, INTENT(IN)                                :: is_ortho
      REAL(KIND=dp), INTENT(IN)                          :: zet_err_minimax
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zet_err_cutoff
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: l_err_cutoff
      INTEGER, INTENT(IN)                                :: l_max
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env

      CALL eri_mme_set_params(param%par, hmat, is_ortho, zet_err_minimax, zet_err_cutoff, l_err_cutoff, l_max, para_env)

      CALL eri_mme_print_info(param)
   END SUBROUTINE eri_mme_set_params_custom

! **************************************************************************************************
!> \brief Get basis parameters for estimating cutoff and minimax error from cp2k basis
!> \param qs_kind_set ...
!> \param basis_type_1 ...
!> \param basis_type_2 ...
!> \param zet_mm Smallest exponent, used to estimate error due to minimax approx.
!> \param zet_c  contains max. exponent and the largest exponent for max. l,
!>               used to estimate cutoff error
!> \param l_c    contains the largest l for max. exponent and max. l,
!>               used to estimate cutoff error
! **************************************************************************************************
   SUBROUTINE error_est_pgf_params_from_basis(qs_kind_set, basis_type_1, basis_type_2, zet_mm, zet_c, l_c)
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      CHARACTER(len=*), INTENT(IN)                       :: basis_type_1
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type_2
      REAL(KIND=dp), INTENT(OUT)                         :: zet_mm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: zet_c
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: l_c

      CHARACTER(LEN=*), PARAMETER :: routineN = 'error_est_pgf_params_from_basis', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: basis_type
      INTEGER                                            :: handle, i, ibasis, ikind, iset, j, l_m, &
                                                            l_zet, nbasis, nkind, nset
      INTEGER, DIMENSION(:), POINTER                     :: npgf
      REAL(KIND=dp)                                      :: zet_l, zet_m
      TYPE(gto_basis_set_type), POINTER                  :: basis_set

      CALL timeset(routineN, handle)

      l_m = 0
      zet_m = 0.0_dp
      l_zet = -1
      zet_l = -1.0_dp
      zet_mm = -1.0_dp

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
            l_m = MAX(l_m, MAXVAL(basis_set%l(:, :)))
            DO iset = 1, nset
               zet_m = MAX(zet_m, MAXVAL(basis_set%zet(1:npgf(iset), iset)))
               IF (zet_mm .LT. 0.0_dp) THEN
                  zet_mm = MINVAL(basis_set%zet(1:npgf(iset), iset))
               ELSE
                  zet_mm = MIN(zet_mm, MINVAL(basis_set%zet(1:npgf(iset), iset)))
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
            DO i = LBOUND(basis_set%l, 1), UBOUND(basis_set%l, 1)
               DO j = LBOUND(basis_set%l, 2), UBOUND(basis_set%l, 2)
                  IF (ABS(zet_m-basis_set%zet(i, j)) .LE. (zet_m*1.0E-12_dp) .AND. (basis_set%l(i, j) .GT. l_zet)) THEN
                     l_zet = basis_set%l(i, j)
                  ENDIF
                  IF ((basis_set%l(i, j) .EQ. l_m) .AND. (basis_set%zet(i, j) .GT. zet_l)) THEN
                     zet_l = basis_set%zet(i, j)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      ALLOCATE (l_c(2))
      ALLOCATE (zet_c(2))

      CPASSERT(zet_l .GT. 0.0_dp .AND. l_zet .GE. 0)

      zet_c(:) = [zet_m, zet_l]
      l_c(:) = [l_zet, l_m]

      CALL timestop(handle)
   END SUBROUTINE error_est_pgf_params_from_basis

! **************************************************************************************************
!> \brief ...
!> \param param ...
! **************************************************************************************************
   SUBROUTINE eri_mme_print_info(param)
      TYPE(cp_eri_mme_param)                             :: param

      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => param%logger

      unit_nr = param%par%unit_nr
      IF (unit_nr > 0) WRITE (unit_nr, '(/T2, A, T71, F10.1)') "ERI_MME| Cutoff for ERIs [a.u.]:", param%par%cutoff
      IF (param%par%is_ortho) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(/T2, A, T78, I3/)') "ERI_MME| Number of terms in minimax approximation:", param%par%n_minimax
            WRITE (unit_nr, '(T2, A)') "ERI_MME| Estimated absolute error for normalized Hermite-Gaussian basis"
            WRITE (unit_nr, '(T2, A, T72, ES9.2)') "ERI_MME|   Minimax error:", param%par%err_mm
            WRITE (unit_nr, '(T2, A, T72, ES9.2)') "ERI_MME|   Cutoff error:", param%par%err_c
            WRITE (unit_nr, '(T2, A, T72, ES9.2)') "ERI_MME|   Total error (minimax + cutoff):", param%par%err_mm+param%par%err_c
            IF (param%par%print_calib) &
               WRITE (unit_nr, '(T2, A, T68, F13.10)') "ERI_MME|   Minimax scaling constant in AM-GM estimate:", param%par%C_mm
         ENDIF
      ENDIF
      IF (unit_nr > 0) WRITE (unit_nr, *)

   END SUBROUTINE eri_mme_print_info

! **************************************************************************************************
!> \brief Create input section for unit testing
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_eri_mme_test_section(section)
      TYPE(section_type), INTENT(INOUT), POINTER         :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_eri_mme_test_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      NULLIFY (keyword, subsection)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, "ERI_MME_TEST", &
                          description="Parameters for testing the ERI_MME method for electron repulsion integrals. "// &
                          "Testing w.r.t. performance and accuracy. ", &
                          n_keywords=5, n_subsections=1)

      CALL create_eri_mme_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation the ERI_MME test. ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TEST_3C", &
                          description="Whether to test 3-center integrals.", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TEST_2C", &
                          description="Whether to test 2-center integrals.", &
                          default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="ABC", &
                          description="Specify the lengths of the cell vectors A, B, and C. ", &
                          usage="ABC 10.000 10.000 10.000", unit_str="angstrom", &
                          n_var=3, type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="MIN_NPOS", &
                          description="Minimum number of atomic distances to consider. ", &
                          default_i_val=8)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="NREP", &
                          description="Number of repeated calculation of each integral. ", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="CHECK_2C_ACCURACY", &
                          description="Whether integrals should be compared to reference values, "// &
                          "created on the fly by exact method (G-space sum on grid without "// &
                          "minimax approximation). Note: only feasible for not too many "// &
                          "integrals and maximum exponent around 10.0. ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)

      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="LMAX", &
                          description="Maximum total angular momentum quantum number. ", &
                          default_i_val=6)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="ZET_MIN", &
                          description="Minimum exponent. ", &
                          default_r_val=0.001_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="ZET_MAX", &
                          description="Maximum exponent. ", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="NZET", &
                          description="Number of exponents (logarithmic partition between ZET_MIN and ZET_MAX). ", &
                          default_i_val=4)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="NSAMPLE_3C", &
                          description="If NSAMPLE_3C = N, only calculate every Nth 3-center integral.", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
   END SUBROUTINE create_eri_mme_test_section

END MODULE cp_eri_mme_types
