!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types and initialization / release routines for Minimax-Ewald method for electron
!>        repulsion integrals.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE eri_mme_types

   USE cp_para_types,                   ONLY: cp_para_env_type
   USE eri_mme_error_control,           ONLY: calibrate_cutoff,&
                                              cutoff_minimax_error
   USE eri_mme_gaussian,                ONLY: get_minimax_coeff_v_gspace
   USE eri_mme_util,                    ONLY: G_abs_min,&
                                              R_abs_min
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3
   USE orbital_pointers,                ONLY: init_orbital_pointers
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_types'

   INTEGER, PARAMETER, PUBLIC :: n_minimax_max = 53

   PUBLIC :: eri_mme_param, &
             eri_mme_init, &
             eri_mme_release, &
             eri_mme_set_params

   TYPE eri_mme_param
      INTEGER                          :: n_minimax
      REAL(KIND=dp), POINTER, &
         DIMENSION(:)                  :: minimax_aw => NULL()
      REAL(KIND=dp), DIMENSION(3, 3)   :: hmat, h_inv
      REAL(KIND=dp)                    :: vol
      LOGICAL                          :: is_ortho
      REAL(KIND=dp)                    :: cutoff
      LOGICAL                          :: do_calib_cutoff
      LOGICAL                          :: print_calib
      REAL(KIND=dp)                    :: cutoff_min, cutoff_max, cutoff_delta, &
                                          cutoff_eps, err_mm, err_c
      REAL(KIND=dp)                    :: G_min, R_min
      LOGICAL                          :: is_valid
      LOGICAL                          :: debug
      REAL(KIND=dp)                    :: debug_delta
      REAL(KIND=dp)                    :: C_mm
      INTEGER                          :: unit_nr
   END TYPE eri_mme_param

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param param ...
!> \param n_minimax ...
!> \param cutoff ...
!> \param do_calib_cutoff ...
!> \param cutoff_min ...
!> \param cutoff_max ...
!> \param cutoff_eps ...
!> \param cutoff_delta ...
!> \param debug ...
!> \param debug_delta ...
!> \param unit_nr ...
!> \param print_calib ...
! **************************************************************************************************
   SUBROUTINE eri_mme_init(param, n_minimax, cutoff, do_calib_cutoff, &
                           cutoff_min, cutoff_max, cutoff_eps, cutoff_delta, &
                           debug, debug_delta, unit_nr, print_calib)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      INTEGER, INTENT(IN)                                :: n_minimax
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      LOGICAL, INTENT(IN)                                :: do_calib_cutoff
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_min, cutoff_max, cutoff_eps, &
                                                            cutoff_delta
      LOGICAL, INTENT(IN)                                :: debug
      REAL(KIND=dp), INTENT(IN)                          :: debug_delta
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: print_calib

      CHARACTER(len=2)                                   :: string

      WRITE (string, '(I2)') n_minimax_max
      IF (n_minimax .GT. n_minimax_max) &
         CPABORT("The maximum allowed number of minimax points N_MINIMAX is "//TRIM(string))

      param%n_minimax = n_minimax
      param%cutoff = cutoff
      param%do_calib_cutoff = do_calib_cutoff
      param%cutoff_min = cutoff_min
      param%cutoff_max = cutoff_max
      param%cutoff_eps = cutoff_eps
      param%cutoff_delta = cutoff_delta
      param%debug = debug
      param%debug_delta = debug_delta
      param%print_calib = print_calib
      param%unit_nr = unit_nr

      param%is_valid = .FALSE.
      ALLOCATE (param%minimax_aw(2*param%n_minimax))
   END SUBROUTINE eri_mme_init

! **************************************************************************************************
!> \brief Set parameters for MME method with manual specification of basis parameters.
!>        Takes care of cutoff calibration if requested.
!> \param param ...
!> \param hmat ...
!> \param is_ortho ...
!> \param zet_err_minimax Exponent used to estimate error of minimax approximation.
!> \param zet_err_cutoff  Exponents used to estimate error of finite cutoff.
!> \param l_err_cutoff    Total ang. mom. quantum numbers l to be combined with exponents in
!>                        zet_err_cutoff.
!> \param l_max           Maximum total angular momentum quantum number
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE eri_mme_set_params(param, hmat, is_ortho, zet_err_minimax, zet_err_cutoff, l_err_cutoff, l_max, para_env)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat
      LOGICAL, INTENT(IN)                                :: is_ortho
      REAL(KIND=dp), INTENT(IN)                          :: zet_err_minimax
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zet_err_cutoff
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: l_err_cutoff
      INTEGER, INTENT(IN)                                :: l_max
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_set_params', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, l_mm
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: l_c
      LOGICAL                                            :: s_only
      REAL(KIND=dp)                                      :: cutoff, E_mm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: zet_c

      CALL timeset(routineN, handle)

      ! Note: in MP2 default logger hacked and does not use global default print level
      s_only = l_max .EQ. 0

      CPASSERT(PRESENT(l_err_cutoff) .OR. s_only)

      IF (PRESENT(l_err_cutoff)) THEN
         CPASSERT(SIZE(zet_err_cutoff) .EQ. SIZE(l_err_cutoff))
      ENDIF

      CALL init_orbital_pointers(l_max)

      ! l values for minimax error estimate (l_mm) and for cutoff error estimate (l_c)
      l_mm = MERGE(0, 1, s_only)
      ALLOCATE (l_c(SIZE(zet_err_cutoff)), zet_c(SIZE(zet_err_cutoff)))
      IF (s_only) THEN
         l_c(:) = 0
      ELSE
         l_c(:) = l_err_cutoff(:)
      ENDIF

      ! cell info
      ! Note: we recompute basic quantities from hmat to avoid dependency on cp2k cell type
      param%hmat = hmat
      param%h_inv = inv_3x3(hmat)
      param%vol = ABS(det_3x3(hmat))
      param%is_ortho = is_ortho

      ! Minimum lattice vectors
      param%G_min = G_abs_min(param%h_inv)
      param%R_min = R_abs_min(param%hmat)

      ! cutoff calibration not yet implemented for general cell
      IF (.NOT. param%is_ortho) THEN
         param%do_calib_cutoff = .FALSE.
      ENDIF

      ! Cutoff calibration and error estimate for orthorhombic cell
      IF (param%do_calib_cutoff) THEN
         CALL calibrate_cutoff(param%hmat, param%h_inv, param%G_min, param%vol, &
                               zet_err_minimax, l_mm, zet_err_cutoff, l_c, param%n_minimax, &
                               param%cutoff_min, param%cutoff_max, param%cutoff_eps, &
                               param%cutoff_delta, cutoff, param%err_mm, param%err_c, &
                               param%C_mm, para_env, param%print_calib, param%unit_nr)

         param%cutoff = cutoff
      ELSE IF (param%is_ortho) THEN
         CALL cutoff_minimax_error(param%cutoff, param%hmat, param%h_inv, param%vol, param%G_min, &
                                   zet_err_minimax, l_mm, zet_err_cutoff, l_c, param%n_minimax, &
                                   param%minimax_aw, param%err_mm, param%err_c, param%C_mm, para_env)
      ENDIF

      CALL get_minimax_coeff_v_gspace(param%n_minimax, param%cutoff, param%G_min, param%minimax_aw, &
                                      E_mm)

      param%is_valid = .TRUE.

      CALL timestop(handle)
   END SUBROUTINE eri_mme_set_params

! **************************************************************************************************
!> \brief ...
!> \param param ...
! **************************************************************************************************
   SUBROUTINE eri_mme_release(param)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param

      IF (ASSOCIATED(param%minimax_aw)) DEALLOCATE (param%minimax_aw)
   END SUBROUTINE eri_mme_release

END MODULE eri_mme_types
