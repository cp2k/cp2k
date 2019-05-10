!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
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
             eri_mme_set_params, &
             eri_mme_print_grid_info, &
             get_minimax_from_cutoff

   TYPE minimax_grid
      REAL(KIND=dp)                    :: cutoff
      INTEGER                          :: n_minimax
      REAL(KIND=dp), POINTER, &
         DIMENSION(:)                  :: minimax_aw => NULL()
      REAL(KIND=dp)                    :: error
   END TYPE

   TYPE eri_mme_param
      INTEGER                          :: n_minimax
      REAL(KIND=dp), DIMENSION(3, 3)   :: hmat, h_inv
      REAL(KIND=dp)                    :: vol
      LOGICAL                          :: is_ortho
      REAL(KIND=dp)                    :: cutoff
      LOGICAL                          :: do_calib_cutoff
      LOGICAL                          :: do_error_est
      LOGICAL                          :: print_calib
      REAL(KIND=dp)                    :: cutoff_min, cutoff_max, cutoff_delta, &
                                          cutoff_eps
      REAL(KIND=dp)                    :: err_mm, err_c
      REAL(KIND=dp)                    :: G_min, R_min
      LOGICAL                          :: is_valid
      LOGICAL                          :: debug
      REAL(KIND=dp)                    :: debug_delta
      INTEGER                          :: debug_nsum
      REAL(KIND=dp)                    :: C_mm
      INTEGER                          :: unit_nr
      REAL(KIND=dp)                    :: sum_precision
      INTEGER                          :: n_grids
      TYPE(minimax_grid), DIMENSION(:), &
         ALLOCATABLE                   :: minimax_grid
   END TYPE eri_mme_param

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param param ...
!> \param n_minimax ...
!> \param cutoff ...
!> \param do_calib_cutoff ...
!> \param do_error_est ...
!> \param cutoff_min ...
!> \param cutoff_max ...
!> \param cutoff_eps ...
!> \param cutoff_delta ...
!> \param sum_precision ...
!> \param debug ...
!> \param debug_delta ...
!> \param debug_nsum ...
!> \param unit_nr ...
!> \param print_calib ...
! **************************************************************************************************
   SUBROUTINE eri_mme_init(param, n_minimax, cutoff, do_calib_cutoff, do_error_est, &
                           cutoff_min, cutoff_max, cutoff_eps, cutoff_delta, sum_precision, &
                           debug, debug_delta, debug_nsum, unit_nr, print_calib)
      TYPE(eri_mme_param), INTENT(OUT)                   :: param
      INTEGER, INTENT(IN)                                :: n_minimax
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      LOGICAL, INTENT(IN)                                :: do_calib_cutoff, do_error_est
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_min, cutoff_max, cutoff_eps, &
                                                            cutoff_delta, sum_precision
      LOGICAL, INTENT(IN)                                :: debug
      REAL(KIND=dp), INTENT(IN)                          :: debug_delta
      INTEGER, INTENT(IN)                                :: debug_nsum, unit_nr
      LOGICAL, INTENT(IN)                                :: print_calib

      CHARACTER(len=2)                                   :: string

      WRITE (string, '(I2)') n_minimax_max
      IF (n_minimax .GT. n_minimax_max) &
         CPABORT("The maximum allowed number of minimax points N_MINIMAX is "//TRIM(string))

      param%n_minimax = n_minimax
      param%n_grids = 1
      param%cutoff = cutoff
      param%do_calib_cutoff = do_calib_cutoff
      param%do_error_est = do_error_est
      param%cutoff_min = cutoff_min
      param%cutoff_max = cutoff_max
      param%cutoff_eps = cutoff_eps
      param%cutoff_delta = cutoff_delta
      param%sum_precision = sum_precision
      param%debug = debug
      param%debug_delta = debug_delta
      param%debug_nsum = debug_nsum
      param%print_calib = print_calib
      param%unit_nr = unit_nr

      param%is_valid = .FALSE.
      ALLOCATE (param%minimax_grid(param%n_grids))
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
      REAL(KIND=dp), INTENT(IN)                          :: zet_err_minimax, zet_err_cutoff
      INTEGER, INTENT(IN)                                :: l_err_cutoff, l_max
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_set_params', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_grid, l_c, l_mm, n_grids, &
                                                            n_minimax, n_mm
      LOGICAL                                            :: s_only
      REAL(KIND=dp) :: cutoff, cutoff_delta, cutoff_max, cutoff_min, cutoff_rel, E_mm, err_c, &
         err_c_prev, err_mm, err_mm_prev, err_target, exp_max, exp_min
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: minimax_aw, minimax_aw_prev

      !REAL(KIND=dp)                                     :: zet_c

      CALL timeset(routineN, handle)

      ! Note: in MP2 default logger hacked and does not use global default print level
      s_only = l_max .EQ. 0

      CALL init_orbital_pointers(3*l_max) ! allow for orbital pointers of combined index

      ! l values for minimax error estimate (l_mm) and for cutoff error estimate (l_c)
      l_mm = MERGE(0, 1, s_only)
      l_c = l_err_cutoff

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
         param%do_error_est = .FALSE.
      ENDIF

      n_grids = param%n_grids

      ! Cutoff calibration and error estimate for orthorhombic cell
      IF (param%do_calib_cutoff) THEN
         CALL calibrate_cutoff(param%hmat, param%h_inv, param%G_min, param%vol, &
                               zet_err_minimax, l_mm, zet_err_cutoff, l_c, param%n_minimax, &
                               param%cutoff_min, param%cutoff_max, param%cutoff_eps, &
                               param%cutoff_delta, cutoff, param%err_mm, param%err_c, &
                               param%C_mm, para_env, param%print_calib, param%unit_nr)

         param%cutoff = cutoff
      ELSE IF (param%do_error_est) THEN
         ALLOCATE (minimax_aw(2*param%n_minimax))
         CALL cutoff_minimax_error(param%cutoff, param%hmat, param%h_inv, param%vol, param%G_min, &
                                   zet_err_minimax, l_mm, zet_err_cutoff, l_c, param%n_minimax, &
                                   minimax_aw, param%err_mm, param%err_c, param%C_mm, para_env)
         DEALLOCATE (minimax_aw)
      ELSE
         param%err_mm = -1.0_dp
         param%err_c = -1.0_dp
      ENDIF

      ALLOCATE (minimax_aw(2*param%n_minimax))
      CALL get_minimax_coeff_v_gspace(param%n_minimax, param%cutoff, param%G_min, minimax_aw, E_mm)

      DEALLOCATE (minimax_aw)

      err_target = E_mm

      exp_max = zet_err_cutoff
      exp_min = zet_err_minimax

      CPASSERT(exp_max+1.0E-12 .GT. exp_min)
      CPASSERT(n_grids .GE. 1)

      cutoff_max = param%cutoff
      cutoff_rel = param%cutoff/exp_max
      cutoff_min = exp_min*cutoff_rel
      cutoff_delta = (cutoff_max/cutoff_min)**(1.0_dp/(n_grids))

      n_minimax = param%n_minimax
      cutoff = cutoff_max
      DO i_grid = n_grids, 1, -1
         DO n_mm = n_minimax, 1, -1
            ALLOCATE (minimax_aw(2*n_mm))
            CALL get_minimax_coeff_v_gspace(n_mm, cutoff, param%G_min, minimax_aw, err_mm)
            IF (err_mm .GT. 1.1_dp*err_target) THEN
               CPASSERT(n_mm .NE. n_minimax)
               CALL create_minimax_grid(param%minimax_grid(i_grid), cutoff, n_mm+1, minimax_aw_prev, err_mm_prev)
               DEALLOCATE (minimax_aw)
               EXIT
            ENDIF
            IF (ALLOCATED(minimax_aw_prev)) DEALLOCATE (minimax_aw_prev)
            ALLOCATE (minimax_aw_prev(2*n_mm))
            minimax_aw_prev(:) = minimax_aw(:)
            DEALLOCATE (minimax_aw)
            err_mm_prev = err_mm
            err_c_prev = err_c
         ENDDO
         cutoff = cutoff/cutoff_delta
      ENDDO

      param%is_valid = .TRUE.

      CALL timestop(handle)
   END SUBROUTINE eri_mme_set_params

! **************************************************************************************************
!> \brief ...
!> \param grid ...
!> \param cutoff ...
!> \param n_minimax ...
!> \param minimax_aw ...
!> \param error ...
! **************************************************************************************************
   SUBROUTINE create_minimax_grid(grid, cutoff, n_minimax, minimax_aw, error)
      TYPE(minimax_grid), INTENT(OUT)                    :: grid
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      INTEGER, INTENT(IN)                                :: n_minimax
      REAL(KIND=dp), DIMENSION(2*n_minimax), INTENT(IN)  :: minimax_aw
      REAL(KIND=dp), INTENT(IN)                          :: error

      grid%cutoff = cutoff
      grid%n_minimax = n_minimax
      ALLOCATE (grid%minimax_aw(2*n_minimax))
      grid%minimax_aw(:) = minimax_aw(:)
      grid%error = error

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param grids ...
!> \param cutoff ...
!> \param n_minimax ...
!> \param minimax_aw ...
!> \param grid_no ...
! **************************************************************************************************
   SUBROUTINE get_minimax_from_cutoff(grids, cutoff, n_minimax, minimax_aw, grid_no)
      TYPE(minimax_grid), DIMENSION(:), INTENT(IN)       :: grids
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      INTEGER, INTENT(OUT)                               :: n_minimax
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT), POINTER  :: minimax_aw
      INTEGER, INTENT(OUT)                               :: grid_no

      INTEGER                                            :: igrid

      grid_no = 0
      DO igrid = 1, SIZE(grids)
         IF (grids(igrid)%cutoff .GE. cutoff/2) THEN
            n_minimax = grids(igrid)%n_minimax
            minimax_aw => grids(igrid)%minimax_aw
            grid_no = igrid
            EXIT
         ENDIF
      ENDDO
      IF (grid_no == 0) THEN
         igrid = SIZE(grids)
         n_minimax = grids(igrid)%n_minimax
         minimax_aw => grids(igrid)%minimax_aw
         grid_no = igrid
      ENDIF
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param grid ...
! **************************************************************************************************
   SUBROUTINE destroy_minimax_grid(grid)
      TYPE(minimax_grid), INTENT(INOUT)                  :: grid

      DEALLOCATE (grid%minimax_aw)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param grid ...
!> \param grid_no ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE eri_mme_print_grid_info(grid, grid_no, unit_nr)
      TYPE(minimax_grid), INTENT(IN)                     :: grid
      INTEGER, INTENT(IN)                                :: grid_no, unit_nr

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2, A, 1X, I2)') "ERI_MME | Info for grid no.", grid_no
         WRITE (unit_nr, '(T2, A, 1X, ES9.2)') "ERI_MME | Cutoff", grid%cutoff
         WRITE (unit_nr, '(T2, A, 1X, I2)') "ERI_MME | Number of minimax points", grid%n_minimax
         WRITE (unit_nr, '(T2, A, 1X, 2ES9.2)') "ERI_MME | minimax error", grid%error
         WRITE (unit_nr, *)
      ENDIF

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param param ...
! **************************************************************************************************
   SUBROUTINE eri_mme_release(param)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param

      INTEGER                                            :: igrid

      IF (ALLOCATED(param%minimax_grid)) THEN
         DO igrid = 1, param%n_grids
            CALL destroy_minimax_grid(param%minimax_grid(igrid))
         ENDDO
         DEALLOCATE (param%minimax_grid)
      ENDIF
   END SUBROUTINE eri_mme_release

END MODULE eri_mme_types
