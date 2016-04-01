!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types and initialization routines for Minimax-Ewald method for electron
!>        repulsion integrals.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE eri_mme_types

   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_section_create,&
                                              cp_print_key_unit_nr,&
                                              medium_print_level
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
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_types'

   INTEGER, PARAMETER, PUBLIC :: lmax = 6
   REAL(KIND=dp), PARAMETER, PUBLIC :: lattice_sum_prec = 1.0E-20_dp ! truncation threshold for summing PGFs
   REAL(KIND=dp), PARAMETER, PUBLIC :: min_Rc = 2.0E+01_dp ! minimum allowed value for minimax upper bound Rc
   INTEGER, PARAMETER, PUBLIC :: dynamic_lattice_space = 0, rec_lattice_space = 1, dir_lattice_space = 2
   INTEGER, PARAMETER, PUBLIC :: n_minimax_max = 15

   PUBLIC :: eri_mme_param, &
             eri_mme_init_read_input, &
             eri_mme_release, &
             create_eri_mme_section

   TYPE eri_mme_param
      TYPE(section_vals_type), &
         POINTER                    :: mme_section => NULL()
      INTEGER                          :: n_minimax
      REAL(KIND=dp), POINTER, &
         DIMENSION(:)               :: minimax_aw => NULL()
      REAL(KIND=dp), DIMENSION(3, 3) :: hmat, h_inv
      REAL(KIND=dp)                  :: vol
      LOGICAL                          :: is_ortho
      REAL(KIND=dp)                  :: cutoff
      LOGICAL                          :: do_calib_cutoff, do_calib_alpha
      LOGICAL                          :: print_calib
      REAL(KIND=dp)                  :: cutoff_min, cutoff_max, cutoff_delta, &
                                        cutoff_eps, error_tot
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :)            :: alpha_c => NULL()

      INTEGER                          :: G_count, R_count
      REAL(KIND=dp)                  :: G_min, R_min
      LOGICAL                          :: is_valid
      INTEGER                          :: lmax
   END TYPE eri_mme_param

CONTAINS

! **************************************************************************************************
!> \brief Create main input section
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_eri_mme_section(section)
      TYPE(section_type), INTENT(INOUT), POINTER         :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_eri_mme_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key, subsection

      NULLIFY (keyword, print_key, subsection)
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, "ERI_MME", &
                          description="Parameters for the calculation of electron repulsion "// &
                          "integrals (ERI) using the Minimax-Ewald (MME) method. "// &
                          "Note: N_MINIMAX is the only parameter to be tuned for accuracy, "// &
                          "all other parameters can be left to default. "// &
                          "Method may not yet be sufficiently accurate since the number "// &
                          "of terms in minimax approximation is limited to 15 in the current "// &
                          "implementation.", &
                          n_keywords=5, n_subsections=1)

      CALL keyword_create(keyword, &
                          name="N_MINIMAX", &
                          description="Number of terms in minimax approximation of "// &
                          "reciprocal space potential. ", &
                          default_i_val=15)
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

      CALL keyword_create(keyword, &
                          name="DO_CALIBRATE_GR_SWITCH", &
                          description="Whether the optimal compound exponent for switching between G and R space "// &
                          "summation shall be calibrated w.r.t. CPU time. Not implemented "// &
                          "for non-orthorhombic cells. Note: this feature is experimental. ", &
                          default_l_val=.FALSE., &
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
   SUBROUTINE eri_mme_init_read_input(mme_section, param)
      TYPE(section_vals_type), INTENT(INOUT), POINTER    :: mme_section
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param

      CHARACTER(len=*), PARAMETER :: routineN = 'eri_mme_init_read_input', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: string
      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: subsection

      logger => cp_get_default_logger()
      unit_nr = cp_print_key_unit_nr(logger, mme_section, "ERI_MME_INFO", &
                                     extension=".eri_mme")

      NULLIFY (subsection)
      CALL section_vals_val_get(mme_section, "N_MINIMAX", i_val=param%n_minimax)
      WRITE (string, '(I2)') n_minimax_max
      IF (param%n_minimax .GT. n_minimax_max) &
         CPABORT("The maximum allowed number of minimax points N_MINIMAX is "//TRIM(string))
      CALL section_vals_val_get(mme_section, "CUTOFF", r_val=param%cutoff)
      CALL section_vals_val_get(mme_section, "DO_CALIBRATE_GR_SWITCH", l_val=param%do_calib_alpha)
      CALL section_vals_val_get(mme_section, "DO_CALIBRATE_CUTOFF", l_val=param%do_calib_cutoff)
      CALL section_vals_val_get(mme_section, "PRINT_CALIB", l_val=param%print_calib)
      subsection => section_vals_get_subs_vals(mme_section, "CUTOFF_CALIB")
      CALL section_vals_val_get(subsection, "MIN", r_val=param%cutoff_min)
      CALL section_vals_val_get(subsection, "MAX", r_val=param%cutoff_max)
      CALL section_vals_val_get(subsection, "EPS", r_val=param%cutoff_eps)
      CALL section_vals_val_get(subsection, "DELTA", r_val=param%cutoff_delta)
      param%mme_section => mme_section

      ALLOCATE (param%minimax_aw(2*param%n_minimax))
      IF (param%do_calib_alpha) ALLOCATE (param%alpha_c(0:2*lmax, 3))

      param%is_valid = .FALSE.

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/T2, A, T78, I3/)') "ERI_MME| Number of terms in minimax approximation:", param%n_minimax
      ENDIF

      CALL cp_print_key_finished_output(unit_nr, logger, param%mme_section, "ERI_MME_INFO")

   END SUBROUTINE eri_mme_init_read_input

! **************************************************************************************************
!> \brief ...
!> \param param ...
! **************************************************************************************************
   SUBROUTINE eri_mme_release(param)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(param%minimax_aw)) DEALLOCATE (param%minimax_aw)
      IF (ASSOCIATED(param%alpha_c)) DEALLOCATE (param%alpha_c)
   END SUBROUTINE eri_mme_release

END MODULE eri_mme_types
