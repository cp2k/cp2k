!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types and basic routines needed for a kpoint calculation
!> \par History
!>       2014.07 created [JGH]
!> \author JGH
! *****************************************************************************
MODULE kpoint_types
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_env,                     ONLY: cp_cart_release,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_cart_type,&
                                             cp_para_env_type
  USE f77_blas
  USE input_cp2k_kpoints,              ONLY: use_complex_wfn,&
                                             use_real_wfn
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE qs_matrix_pools,                 ONLY: mpools_release,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE string_utilities,                ONLY: uppercase
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_types'

  PUBLIC :: kpoint_type
  PUBLIC :: kpoint_create, kpoint_release, get_kpoint_info, set_kpoint_info
  PUBLIC :: read_kpoint_section, write_kpoint_info
  PUBLIC :: kpoint_env_type, kpoint_env_p_type
  PUBLIC :: kpoint_env_create, kpoint_env_release, get_kpoint_env
  PUBLIC :: opmat_type, kp_opmat_create, kp_opmat_release

! *****************************************************************************  
!> \brief Keeps operator matrices in real space form
!> \param index_to_cell   given a real space index, returns a cell triple
!> \param cell_to_index   given a cell triple, returns the real space index
!> \param sab_nl          neighbourlist that defines real space matrices
!> \param matrix_ks       Kohn-Sham matrices (spin,rs index)
!> \param matrix_s        Overlap matrices (spin,rs index)
!> \param matrix_t        Kinetic energy matrices (spin,rs index)
!> \param matrix_vhxc     Potential matrices (spin,rs index)
!> \param matrix_h        Core matrices (spin,rs index)
!> \author JGH   
! *****************************************************************************
  TYPE opmat_type
     INTEGER, DIMENSION(:,:), POINTER                  :: index_to_cell
     INTEGER, DIMENSION(:,:,:), POINTER                :: cell_to_index
     TYPE(neighbor_list_set_p_type), DIMENSION(:),&
        POINTER                                        :: sab_nl
     TYPE(cp_dbcsr_p_type), DIMENSION(:,:), POINTER    :: matrix_ks
     TYPE(cp_dbcsr_p_type), DIMENSION(:,:), POINTER    :: matrix_s
     TYPE(cp_dbcsr_p_type), DIMENSION(:,:), POINTER    :: matrix_t
     TYPE(cp_dbcsr_p_type), DIMENSION(:,:), POINTER    :: matrix_vhxc
     TYPE(cp_dbcsr_p_type), DIMENSION(:,:), POINTER    :: matrix_h
  END TYPE opmat_type

! *****************************************************************************  
!> \brief Keeps information about a specific k-point
!> \param nkpoint   the kpoint index
!> \param wkp       weight of this kpoint
!> \param xkp       kpoint coordinates
!> \param is_local  if this kpoint is calculated on a single thread
!> \param mos       associated MOs (r/i,spin)
!> \param pmat      associated density matrix (r/i,spin)
!> \author JGH   
! *****************************************************************************
  TYPE kpoint_env_type
     INTEGER                                           :: nkpoint
     REAL(KIND=dp)                                     :: wkp
     REAL(KIND=dp), DIMENSION(3)                       :: xkp
     LOGICAL                                           :: is_local
     TYPE(mo_set_p_type), DIMENSION(:,:), POINTER      :: mos
     TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER       :: pmat
  END TYPE kpoint_env_type
 
  TYPE kpoint_env_p_type
     TYPE(kpoint_env_type), POINTER                    :: kpoint_env
  END TYPE kpoint_env_p_type

! *****************************************************************************
!> \brief Contains information about kpoints
!> \par History
!>       2014.07 created [JGH]
!> \param kp_scheme           [input] Type of kpoint grid
!> \param nkp_grid            [input] Grid points
!> \param kp_shift            [input] Shift of the grid
!> \param use_real_wfn        [input] real/complex wfn
!> \param symmetry            [input] use symmetry (atoms) to reduce kpoints
!> \param full_grid           [input] don't reduce kpoints at all
!> \param verbose             [input] more output information
!> \param eps_geo             [input] accuracy of atom symmetry detection
!> \param parallel_group_size [input] kpoint group size
!> \param nkp     number of kpoints
!> \param xkp     kpoint coordinates
!> \param wkp     kpoint weights
!> \param para_env 'global' parallel environment
!> \param cart 2d distribution of the processors for the kpoints,
!>        a column work together on the same kpoint
!> \param para_env_full the global para env that contains all the kpoints,
!>        this is just the cart as para_env
!> \param para_env_kp parallel environment of the kpoint calculation
!> \param para_env_inter_kp parallel environment between kpoints
!> \param iogrp  this kpoint group has the IO processor
!> \param nkp_groups   number of kpoint groups
!> \param kp_dist      kpoints distribution on groups
!> \param kp_range     kpoints distribution for local processor
!> \param blacs_env    BLACS env for the kpoint group
!> \param opmats       Operator matrices
!> \param kp_env       Information for each kpoint
!> \param mpools       FM matrix pools for kpoint groups
!> \author JGH
! *****************************************************************************  
  TYPE kpoint_type
     CHARACTER (LEN=default_string_length)   :: kp_scheme
     INTEGER, DIMENSION(3)                   :: nkp_grid
     REAL(KIND=dp), DIMENSION(3)             :: kp_shift
     LOGICAL                                 :: use_real_wfn
     LOGICAL                                 :: symmetry
     LOGICAL                                 :: full_grid
     LOGICAL                                 :: verbose
     REAL(KIND=dp)                           :: eps_geo
     INTEGER                                 :: parallel_group_size
     INTEGER                                 :: nkp
     REAL(KIND=dp), DIMENSION(:,:), POINTER  :: xkp
     REAL(KIND=dp), DIMENSION(:), POINTER    :: wkp
     ! parallel environment
     TYPE(cp_para_env_type), POINTER         :: para_env
     TYPE(cp_para_cart_type), POINTER        :: cart
     TYPE(cp_para_env_type), POINTER         :: para_env_full,&
                                                para_env_kp,&
                                                para_env_inter_kp
     LOGICAL                                 :: iogrp
     INTEGER                                 :: nkp_groups
     INTEGER, DIMENSION(:,:), POINTER        :: kp_dist
     INTEGER, DIMENSION(2)                   :: kp_range
     TYPE(cp_blacs_env_type), POINTER        :: blacs_env
     ! operators
     TYPE(opmat_type), POINTER               :: opmats
     ! environment
     TYPE(kpoint_env_p_type), DIMENSION(:),&
       POINTER                               :: kp_env
     TYPE(qs_matrix_pools_type), POINTER     :: mpools
  END TYPE kpoint_type

! *****************************************************************************  

CONTAINS

! *****************************************************************************
!> \brief Create a kpoint environment
!> \param kpoint  All the kpoint information
!> \param error   CP2K error handling
!> \author JGH
! *****************************************************************************
  SUBROUTINE kpoint_create(kpoint, error)
    TYPE(kpoint_type), POINTER               :: kpoint
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPAssert(.NOT.ASSOCIATED(kpoint),cp_failure_level,routineP,error,failure)

    ALLOCATE(kpoint,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    kpoint%kp_scheme = ""
    kpoint%nkp_grid = 0
    kpoint%kp_shift = 0.0_dp
    kpoint%symmetry = .FALSE.
    kpoint%verbose  = .FALSE.
    kpoint%full_grid= .FALSE.
    kpoint%use_real_wfn = .FALSE.
    kpoint%eps_geo = 1.0e-6_dp
    kpoint%parallel_group_size = -1

    kpoint%nkp = 0

    NULLIFY(kpoint%xkp,kpoint%wkp)

    NULLIFY(kpoint%para_env)
    NULLIFY(kpoint%cart)
    NULLIFY(kpoint%para_env_full,kpoint%para_env_kp,kpoint%para_env_inter_kp)
    NULLIFY(kpoint%blacs_env)
    kpoint%nkp_groups = 0
    kpoint%iogrp = .FALSE.
    kpoint%kp_range = 0

    NULLIFY(kpoint%opmats)

    NULLIFY(kpoint%kp_env)
    NULLIFY(kpoint%mpools)

  END SUBROUTINE kpoint_create

! *****************************************************************************
!> \brief  Release a kpoint environment, deallocate all data
!> \param kpoint  The kpoint environment
!> \param error   CP2K error handling
!> \author JGH
! *****************************************************************************
  SUBROUTINE kpoint_release(kpoint, error)
    TYPE(kpoint_type), POINTER               :: kpoint
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, ik
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF(ASSOCIATED(kpoint)) THEN

       IF(ASSOCIATED(kpoint%xkp)) THEN
         DEALLOCATE(kpoint%xkp,STAT=ierr)
         CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
         NULLIFY(kpoint%xkp)
       END IF
       IF(ASSOCIATED(kpoint%wkp)) THEN
         DEALLOCATE(kpoint%wkp,STAT=ierr)
         CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
         NULLIFY(kpoint%wkp)
       END IF

       CALL mpools_release(kpoint%mpools, error=error)
       CALL cp_blacs_env_release(kpoint%blacs_env,error=error)

       CALL cp_cart_release(kpoint%cart,error=error)
       CALL cp_para_env_release(kpoint%para_env,error=error)
       CALL cp_para_env_release(kpoint%para_env_full,error=error)
       CALL cp_para_env_release(kpoint%para_env_kp,error=error)
       CALL cp_para_env_release(kpoint%para_env_inter_kp,error=error)

       CALL kp_opmat_release(kpoint%opmats,error=error)

       IF(ASSOCIATED(kpoint%kp_env)) THEN
          DO ik=1,SIZE(kpoint%kp_env)
             CALL kpoint_env_release(kpoint%kp_env(ik)%kpoint_env,error)
          END DO
          DEALLOCATE(kpoint%kp_env,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
       END IF

       DEALLOCATE(kpoint,STAT=ierr)
       CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
       NULLIFY(kpoint)

    END IF

  END SUBROUTINE kpoint_release

! *****************************************************************************
!> \brief Retrieve information from a kpoint environment
!> \param kpoint        The kpoint environment
!> \param kp_scheme     Type of kpoint grid
!> \param nkp_grid      Grid points
!> \param kp_shift      Shift of the grid
!> \param symmetry      use symmetry (atoms) to reduce kpoints
!> \param verbose       more output information
!> \param full_grid     don't reduce kpoints at all
!> \param use_real_wfn  real/complex wfn
!> \param eps_geo       accuracy of atom symmetry detection
!> \param parallel_group_size kpoint group size
!> \param kp_range      kpoints distribution for local processor
!> \param nkp           number of kpoints
!> \param xkp           kpoint coordinates
!> \param wkp           kpoint weights
!> \param para_env      'global' parallel environment
!> \param cart          2d distribution of the processors for the kpoints,
!>        a column work together on the same kpoint
!> \param para_env_full the global para env that contains all the kpoints,
!>        this is just the cart as para_env
!> \param para_env_kp   parallel environment of the kpoint calculation
!> \param para_env_inter_kp   parallel environment between kpoints
!> \param blacs_env     BLACS env for the kpoint group
!> \param kp_env        Information for each kpoint
!> \param mpools        FM matrix pools for kpoint groups
!> \param iogrp         this kpoint group has the IO processor
!> \param nkp_groups    number of kpoint groups
!> \param kp_dist       kpoints distribution on groups
!> \param index_to_cell given a real space index, returns a cell triple
!> \param cell_to_index given a cell triple, returns the real space index
!> \param ncount        number of real space matrices
!> \param sab_nl        neighbourlist that defines real space matrices
!> \param matrix_ks     Kohn-Sham matrices (spin,rs index)
!> \param matrix_s      Overlap matrices (spin,rs index)
!> \param matrix_t      Kinetic energy matrices (spin,rs index)
!> \param matrix_vhxc   Potential matrices (spin,rs index)
!> \param matrix_h      Core matrices (spin,rs index)
!> \param error         CP2K error handling
!> \author JGH
! *****************************************************************************
  SUBROUTINE get_kpoint_info(kpoint, kp_scheme, nkp_grid, kp_shift, symmetry, verbose,&
      full_grid, use_real_wfn, eps_geo, parallel_group_size, kp_range, nkp, xkp, wkp,&
      para_env, cart, para_env_full, para_env_kp, para_env_inter_kp, blacs_env,&
      kp_env, mpools, iogrp, nkp_groups, kp_dist, index_to_cell, cell_to_index, ncount, sab_nl,&
      matrix_ks, matrix_s, matrix_t, matrix_vhxc, matrix_h,&
      error)
    TYPE(kpoint_type), POINTER               :: kpoint
    CHARACTER(LEN=*), OPTIONAL               :: kp_scheme
    INTEGER, DIMENSION(3), OPTIONAL          :: nkp_grid
    REAL(KIND=dp), DIMENSION(3), OPTIONAL    :: kp_shift
    LOGICAL, OPTIONAL                        :: symmetry, verbose, full_grid, &
                                                use_real_wfn
    REAL(KIND=dp), OPTIONAL                  :: eps_geo
    INTEGER, OPTIONAL                        :: parallel_group_size
    INTEGER, DIMENSION(2), OPTIONAL          :: kp_range
    INTEGER, OPTIONAL                        :: nkp
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: xkp
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: wkp
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_para_cart_type), OPTIONAL, &
      POINTER                                :: cart
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env_full, para_env_kp, &
                                                para_env_inter_kp
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env
    TYPE(kpoint_env_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: kp_env
    TYPE(qs_matrix_pools_type), OPTIONAL, &
      POINTER                                :: mpools
    LOGICAL, OPTIONAL                        :: iogrp
    INTEGER, OPTIONAL                        :: nkp_groups
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: kp_dist, index_to_cell
    INTEGER, DIMENSION(:, :, :), OPTIONAL, &
      POINTER                                :: cell_to_index
    INTEGER, OPTIONAL                        :: ncount
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sab_nl
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: matrix_ks, matrix_s, &
                                                matrix_t, matrix_vhxc, &
                                                matrix_h
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kpoint_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPAssert(ASSOCIATED(kpoint),cp_failure_level,routineP,error,failure)

    IF(PRESENT(kp_scheme)) kp_scheme = kpoint%kp_scheme
    IF(PRESENT(nkp_grid)) nkp_grid = kpoint%nkp_grid
    IF(PRESENT(kp_shift)) kp_shift = kpoint%kp_shift
    IF(PRESENT(symmetry)) symmetry = kpoint%symmetry
    IF(PRESENT(verbose)) verbose = kpoint%verbose
    IF(PRESENT(full_grid)) full_grid = kpoint%full_grid
    IF(PRESENT(use_real_wfn)) use_real_wfn = kpoint%use_real_wfn
    IF(PRESENT(eps_geo)) eps_geo = kpoint%eps_geo
    IF(PRESENT(parallel_group_size)) parallel_group_size = kpoint%parallel_group_size

    IF(PRESENT(nkp)) nkp = kpoint%nkp
    IF(PRESENT(wkp)) wkp => kpoint%wkp
    IF(PRESENT(xkp)) xkp => kpoint%xkp

    IF(PRESENT(para_env)) para_env => kpoint%para_env
    IF(PRESENT(cart)) cart => kpoint%cart
    IF(PRESENT(para_env_full)) para_env_full => kpoint%para_env_full
    IF(PRESENT(para_env_kp)) para_env_kp => kpoint%para_env_kp
    IF(PRESENT(para_env_inter_kp)) para_env_inter_kp => kpoint%para_env_inter_kp
    IF(PRESENT(blacs_env)) blacs_env => kpoint%blacs_env

    IF(PRESENT(iogrp)) iogrp = kpoint%iogrp
    IF(PRESENT(kp_range)) kp_range = kpoint%kp_range
    IF(PRESENT(nkp_groups)) nkp_groups = kpoint%nkp_groups
    IF(PRESENT(kp_dist)) kp_dist => kpoint%kp_dist

    IF(PRESENT(kp_env)) kp_env => kpoint%kp_env
    IF(PRESENT(mpools)) mpools => kpoint%mpools

    IF(PRESENT(ncount)) ncount = 0
    IF(ASSOCIATED(kpoint%opmats)) THEN
       IF(PRESENT(index_to_cell)) index_to_cell => kpoint%opmats%index_to_cell
       IF(PRESENT(cell_to_index)) cell_to_index => kpoint%opmats%cell_to_index
       IF(PRESENT(ncount)) ncount = SIZE(kpoint%opmats%index_to_cell,2)
       IF(PRESENT(sab_nl)) sab_nl => kpoint%opmats%sab_nl
       IF(PRESENT(matrix_ks)) matrix_ks => kpoint%opmats%matrix_ks
       IF(PRESENT(matrix_s)) matrix_s => kpoint%opmats%matrix_s
       IF(PRESENT(matrix_t)) matrix_t => kpoint%opmats%matrix_t
       IF(PRESENT(matrix_vhxc)) matrix_vhxc => kpoint%opmats%matrix_vhxc
       IF(PRESENT(matrix_h)) matrix_h => kpoint%opmats%matrix_h
    ELSE
       IF(PRESENT(index_to_cell) .OR.  PRESENT(cell_to_index) .OR. PRESENT(ncount) &
          .OR. PRESENT(sab_nl) .OR. PRESENT(matrix_ks) .OR. PRESENT(matrix_s) &
          .OR. PRESENT(matrix_t) .OR. PRESENT(matrix_vhxc) .OR. PRESENT(matrix_h)) THEN
          CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
       END IF
    END IF

  END SUBROUTINE get_kpoint_info

! *****************************************************************************
!> \brief Set information in a kpoint environment
!> \param kpoint        The kpoint environment
!> \param kp_scheme     Type of kpoint grid
!> \param nkp_grid      Grid points
!> \param kp_shift      Shift of the grid
!> \param symmetry      use symmetry (atoms) to reduce kpoints
!> \param verbose       more output information
!> \param full_grid     don't reduce kpoints at all
!> \param use_real_wfn  real/complex wfn
!> \param eps_geo       accuracy of atom symmetry detection
!> \param parallel_group_size kpoint group size
!> \param kp_range      kpoints distribution for local processor
!> \param nkp           number of kpoints
!> \param xkp           kpoint coordinates
!> \param wkp           kpoint weights
!> \param para_env      'global' parallel environment
!> \param cart          2d distribution of the processors for the kpoints,
!>        a column work together on the same kpoint
!> \param para_env_full the global para env that contains all the kpoints,
!>        this is just the cart as para_env
!> \param para_env_kp   parallel environment of the kpoint calculation
!> \param para_env_inter_kp   parallel environment between kpoints
!> \param blacs_env     BLACS env for the kpoint group
!> \param kp_env        Information for each kpoint
!> \param mpools        FM matrix pools for kpoint groups
!> \param iogrp         this kpoint group has the IO processor
!> \param nkp_groups    number of kpoint groups
!> \param kp_dist       kpoints distribution on groups
!> \param index_to_cell given a real space index, returns a cell triple
!> \param cell_to_index given a cell triple, returns the real space index
!> \param sab_nl        neighbourlist that defines real space matrices
!> \param matrix_ks     Kohn-Sham matrices (spin,rs index)
!> \param matrix_s      Overlap matrices (spin,rs index)
!> \param matrix_t      Kinetic energy matrices (spin,rs index)
!> \param matrix_vhxc   Potential matrices (spin,rs index)
!> \param matrix_h      Core matrices (spin,rs index)
!> \param error         CP2K error handling
!> \param ncount        number of real space matrices UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author JGH
! *****************************************************************************
  SUBROUTINE set_kpoint_info(kpoint, kp_scheme, nkp_grid, kp_shift, symmetry, verbose,&
      full_grid, use_real_wfn, eps_geo, parallel_group_size, kp_range, nkp, xkp, wkp,&
      para_env, cart, para_env_full, para_env_kp, para_env_inter_kp, blacs_env,&
      kp_env, mpools, iogrp, nkp_groups, kp_dist, index_to_cell, cell_to_index, sab_nl,&
      matrix_ks, matrix_s, matrix_t, matrix_vhxc, matrix_h,&
      error)
    TYPE(kpoint_type), POINTER               :: kpoint
    CHARACTER(LEN=*), OPTIONAL               :: kp_scheme
    INTEGER, DIMENSION(3), OPTIONAL          :: nkp_grid
    REAL(KIND=dp), DIMENSION(3), OPTIONAL    :: kp_shift
    LOGICAL, OPTIONAL                        :: symmetry, verbose, full_grid, &
                                                use_real_wfn
    REAL(KIND=dp), OPTIONAL                  :: eps_geo
    INTEGER, OPTIONAL                        :: parallel_group_size
    INTEGER, DIMENSION(2), OPTIONAL          :: kp_range
    INTEGER, OPTIONAL                        :: nkp
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: xkp
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: wkp
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_para_cart_type), OPTIONAL, &
      POINTER                                :: cart
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env_full, para_env_kp, &
                                                para_env_inter_kp
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env
    TYPE(kpoint_env_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: kp_env
    TYPE(qs_matrix_pools_type), OPTIONAL, &
      POINTER                                :: mpools
    LOGICAL, OPTIONAL                        :: iogrp
    INTEGER, OPTIONAL                        :: nkp_groups
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: kp_dist, index_to_cell
    INTEGER, DIMENSION(:, :, :), OPTIONAL, &
      POINTER                                :: cell_to_index
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sab_nl
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: matrix_ks, matrix_s, &
                                                matrix_t, matrix_vhxc, &
                                                matrix_h
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_kpoint_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPAssert(ASSOCIATED(kpoint),cp_failure_level,routineP,error,failure)
 
    IF(PRESENT(kp_scheme)) kpoint%kp_scheme = kp_scheme
    IF(PRESENT(nkp_grid)) kpoint%nkp_grid = nkp_grid
    IF(PRESENT(kp_shift)) kpoint%kp_shift = kp_shift
    IF(PRESENT(symmetry)) kpoint%symmetry = symmetry
    IF(PRESENT(verbose)) kpoint%verbose = verbose
    IF(PRESENT(full_grid)) kpoint%full_grid = full_grid
    IF(PRESENT(use_real_wfn)) kpoint%use_real_wfn = use_real_wfn
    IF(PRESENT(eps_geo)) kpoint%eps_geo = eps_geo
    IF(PRESENT(parallel_group_size)) kpoint%parallel_group_size = parallel_group_size

    IF(PRESENT(nkp)) kpoint%nkp = nkp
    IF(PRESENT(wkp)) kpoint%wkp => wkp
    IF(PRESENT(xkp)) kpoint%xkp => xkp

    IF(PRESENT(para_env)) kpoint%para_env => para_env
    IF(PRESENT(cart)) kpoint%cart => cart
    IF(PRESENT(para_env_full)) kpoint%para_env_full => para_env_full
    IF(PRESENT(para_env_kp)) kpoint%para_env_kp => para_env_kp
    IF(PRESENT(para_env_inter_kp)) kpoint%para_env_inter_kp => para_env_inter_kp
    IF(PRESENT(blacs_env)) kpoint%blacs_env => blacs_env

    IF(PRESENT(iogrp)) kpoint%iogrp = iogrp
    IF(PRESENT(kp_range)) kpoint%kp_range = kp_range
    IF(PRESENT(nkp_groups)) kpoint%nkp_groups = nkp_groups
    IF(PRESENT(kp_dist)) kpoint%kp_dist => kp_dist

    IF(PRESENT(kp_env)) kpoint%kp_env => kp_env
    IF(PRESENT(mpools)) kpoint%mpools => mpools

    IF(ASSOCIATED(kpoint%opmats)) THEN
       IF(PRESENT(index_to_cell)) kpoint%opmats%index_to_cell => index_to_cell
       IF(PRESENT(cell_to_index)) kpoint%opmats%cell_to_index => cell_to_index
       IF(PRESENT(sab_nl)) kpoint%opmats%sab_nl => sab_nl
       IF(PRESENT(matrix_ks)) kpoint%opmats%matrix_ks => matrix_ks
       IF(PRESENT(matrix_s)) kpoint%opmats%matrix_s => matrix_s
       IF(PRESENT(matrix_t)) kpoint%opmats%matrix_t => matrix_t
       IF(PRESENT(matrix_vhxc)) kpoint%opmats%matrix_vhxc => matrix_vhxc
       IF(PRESENT(matrix_h)) kpoint%opmats%matrix_h => matrix_h
    ELSE
       CPAssert(.NOT.PRESENT(index_to_cell),cp_failure_level,routineP,error,failure)
       CPAssert(.NOT.PRESENT(cell_to_index),cp_failure_level,routineP,error,failure)
       CPAssert(.NOT.PRESENT(sab_nl),cp_failure_level,routineP,error,failure)
       CPAssert(.NOT.PRESENT(matrix_ks),cp_failure_level,routineP,error,failure)
       CPAssert(.NOT.PRESENT(matrix_s),cp_failure_level,routineP,error,failure)
       CPAssert(.NOT.PRESENT(matrix_t),cp_failure_level,routineP,error,failure)
       CPAssert(.NOT.PRESENT(matrix_h),cp_failure_level,routineP,error,failure)
       CPAssert(.NOT.PRESENT(matrix_vhxc),cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE set_kpoint_info

! *****************************************************************************
!> \brief  Creates an operator matrix type
!> \param opmats  The operator matrix type
!> \param error   CP2K error handling
!> \author JGH
! *****************************************************************************
  SUBROUTINE kp_opmat_create(opmats, error)
    TYPE(opmat_type), POINTER                :: opmats
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kp_opmat_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPAssert(.NOT.ASSOCIATED(opmats),cp_failure_level,routineP,error,failure)

    ALLOCATE(opmats,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    NULLIFY(opmats%index_to_cell,opmats%cell_to_index)
    NULLIFY(opmats%sab_nl)
    NULLIFY(opmats%matrix_ks,opmats%matrix_s,opmats%matrix_t)
    NULLIFY(opmats%matrix_vhxc,opmats%matrix_h)

  END SUBROUTINE kp_opmat_create

! *****************************************************************************  
!> \brief Releases an operator matrix type
!> \param opmats  The operator matrix type
!> \param error   CP2K error handling
!> \author JGH
! *****************************************************************************
  SUBROUTINE kp_opmat_release(opmats, error)
    TYPE(opmat_type), POINTER                :: opmats
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kp_opmat_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF(ASSOCIATED(opmats)) THEN
       IF(ASSOCIATED(opmats%index_to_cell)) THEN
          DEALLOCATE(opmats%index_to_cell,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          NULLIFY(opmats%index_to_cell)
       END IF
       IF(ASSOCIATED(opmats%cell_to_index)) THEN
          DEALLOCATE(opmats%cell_to_index,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          NULLIFY(opmats%cell_to_index)
       END IF

       NULLIFY(opmats%sab_nl)

       CALL cp_dbcsr_deallocate_matrix_set(opmats%matrix_ks,error)
       CALL cp_dbcsr_deallocate_matrix_set(opmats%matrix_s,error)
       CALL cp_dbcsr_deallocate_matrix_set(opmats%matrix_t,error)
       CALL cp_dbcsr_deallocate_matrix_set(opmats%matrix_vhxc,error)
       CALL cp_dbcsr_deallocate_matrix_set(opmats%matrix_h,error)

       DEALLOCATE(opmats,STAT=ierr)
       CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
       NULLIFY(opmats)
    END IF

  END SUBROUTINE kp_opmat_release

! *****************************************************************************
!> \brief Read the kpoint input section
!> \param kpoint  The kpoint environment
!> \param kpoint_section The input section
!> \param error   CP2K error handling
!> \author JGH
! *****************************************************************************
  SUBROUTINE read_kpoint_section(kpoint,kpoint_section,error)
    TYPE(kpoint_type), POINTER               :: kpoint
    TYPE(section_vals_type), POINTER         :: kpoint_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_kpoint_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER                                  :: i, ierr, n_rep, nval, wfntype
    LOGICAL                                  :: available, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: reallist

    failure = .FALSE.
    CPAssert(ASSOCIATED(kpoint),cp_failure_level,routineP,error,failure)

    CALL section_vals_get(kpoint_section,explicit=available, error=error)
    
    IF (available) THEN
       CALL section_vals_val_get(kpoint_section,"SCHEME",c_vals=tmpstringlist,error=error)
       nval = SIZE(tmpstringlist)
       CPPrecondition(nval >= 1, cp_failure_level, routineP, error, failure)
       kpoint%kp_scheme = tmpstringlist(1)
       CALL uppercase(kpoint%kp_scheme)

       ! SCHEME [None, Gamma, Monkhorst-Pack, MacDonald, General]
       SELECT CASE (kpoint%kp_scheme)
          CASE ("NONE")
             ! do nothing
          CASE ("GAMMA")
             ! do nothing
          CASE ("MONKHORST-PACK")
             CPPrecondition(nval >= 4, cp_failure_level, routineP, error, failure)
             DO i=2,4
                READ(tmpstringlist(i),*) kpoint%nkp_grid(i-1)
             END DO
          CASE ("MACDONALD")
             CPPrecondition(nval >= 7, cp_failure_level, routineP, error, failure)
             DO i=2,4
                READ(tmpstringlist(i),*) kpoint%nkp_grid(i-1)
             END DO
             DO i=5,7
                READ(tmpstringlist(i),*) kpoint%kp_shift(i-4)
             END DO
          CASE ("GENERAL")
             CALL section_vals_val_get(kpoint_section,"KPOINT", n_rep_val=n_rep, error=error)
             kpoint%nkp = n_rep
             ALLOCATE(kpoint%xkp(3,n_rep),kpoint%wkp(n_rep),STAT=ierr)
             CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
             DO i=1,n_rep
                CALL section_vals_val_get(kpoint_section,"KPOINT", i_rep_val=i,&
                                  r_vals=reallist, error=error)
                nval = SIZE(reallist)
                CPPostcondition(nval>=4, cp_failure_level, routineP, error, failure)
                kpoint%xkp(1:3,i) = reallist(1:3)
                kpoint%wkp(i) = reallist(4)
             END DO
          CASE DEFAULT
            CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
       END SELECT

       CALL section_vals_val_get(kpoint_section,"SYMMETRY",l_val=kpoint%symmetry,error=error)
       CALL section_vals_val_get(kpoint_section,"WAVEFUNCTIONS",i_val=wfntype,error=error)
       CALL section_vals_val_get(kpoint_section,"VERBOSE",l_val=kpoint%verbose,error=error)
       CALL section_vals_val_get(kpoint_section,"FULL_GRID",l_val=kpoint%full_grid,error=error)
       CALL section_vals_val_get(kpoint_section,"EPS_GEO",r_val=kpoint%eps_geo,error=error)
       CALL section_vals_val_get(kpoint_section,"PARALLEL_GROUP_SIZE",&
                                 i_val=kpoint%parallel_group_size,error=error)
       SELECT CASE (wfntype)
          CASE (use_real_wfn)
             kpoint%use_real_wfn = .TRUE.
          CASE (use_complex_wfn)
             kpoint%use_real_wfn = .FALSE.
          CASE DEFAULT
            CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
       END SELECT

    ELSE
       kpoint%kp_scheme = "NONE"
    END IF

  END SUBROUTINE read_kpoint_section

! *****************************************************************************
!> \brief Write information on the kpoints to output
!> \param kpoint  The kpoint environment
!> \param dft_section  DFT section information
!> \param error   CP2K error handling
!> \author JGH
! *****************************************************************************
  SUBROUTINE write_kpoint_info(kpoint,dft_section,error)
    TYPE(kpoint_type), POINTER               :: kpoint
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_kpoint_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, punit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    CPAssert(ASSOCIATED(kpoint),cp_failure_level,routineP,error,failure)

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    punit = cp_print_key_unit_nr(logger,dft_section,"PRINT%KPOINTS",extension=".Log",error=error)
    IF (punit>0) THEN

       IF(kpoint%kp_scheme /= "NONE") THEN
          WRITE(punit,'(/," ",79("*"),/,T37,A,/," ",79("*"))') "Kpoints"
       END IF
       SELECT CASE (kpoint%kp_scheme)
          CASE ("NONE")
             ! be silent
          CASE ("GAMMA")
             WRITE (punit,'(A,T57,A)') ' BRILLOUIN|', ' Gamma-point calculation'
          CASE ("MONKHORST-PACK")
             WRITE (punit,'(A,T61,A20)') ' BRILLOUIN| K-point scheme ','      Monkhorst-Pack'
             WRITE (punit,'(A,T66,3I5)') ' BRILLOUIN| K-Point grid', kpoint%nkp_grid
             WRITE (punit,'(A,T66,G15.6)') ' BRILLOUIN| Accuracy in Symmetry determination', kpoint%eps_geo
          CASE ("MACDONALD")
             WRITE (punit,'(A,T71,A10)') ' BRILLOUIN| K-point scheme ',' MacDonald'
             WRITE (punit,'(A,T66,3I5)') ' BRILLOUIN| K-Point grid', kpoint%nkp_grid
             WRITE (punit,'(A,T51,3F10.4)') ' BRILLOUIN| K-Point shift', kpoint%kp_shift
             WRITE (punit,'(A,T66,G15.6)') ' BRILLOUIN| Accuracy in Symmetry determination', kpoint%eps_geo
          CASE ("GENERAL")
             WRITE (punit,'(A,T71,A10)') ' BRILLOUIN| K-point scheme ','   General'
          CASE DEFAULT
            CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
       END SELECT
       IF(kpoint%kp_scheme /= "NONE") THEN
          IF (kpoint%symmetry) THEN
             WRITE (punit,'(A,T76,A)') ' BRILLOUIN| K-Point point group symmetrization', '   ON'
          ELSE
             WRITE (punit,'(A,T76,A)') ' BRILLOUIN| K-Point point group symmetrization', '  OFF'
          END IF
          IF(kpoint%use_real_wfn) THEN
             WRITE (punit,'(A,T76,A)') ' BRILLOUIN| Wavefunction type', ' REAL'
          ELSE
             WRITE (punit,'(A,T73,A)') ' BRILLOUIN| Wavefunction type', ' COMPLEX'
          END IF
          IF(kpoint%full_grid) THEN
             WRITE (punit,'(A,T76,A)') ' BRILLOUIN| Use full k-point grid     '
          END IF
          IF(kpoint%kp_scheme /= "GAMMA") THEN
             WRITE (punit,'(A,T71,I10)') ' BRILLOUIN| List of Kpoints',kpoint%nkp
             WRITE (punit,'(A,T30,A,T48,A,T63,A,T78,A)') &
                   ' BRILLOUIN| Number ', 'Weight', 'X', 'Y', 'Z'
             DO i = 1, kpoint%nkp
                WRITE (punit,'(A,I5,3X,4F15.5)') ' BRILLOUIN| ', i, kpoint%wkp(i), &
                      kpoint%xkp(1,i), kpoint%xkp(2,i), kpoint%xkp(3,i)
             END DO
          END IF
          WRITE(punit,'(" ",79("*"))')
       END IF

    END IF
    CALL cp_print_key_finished_output(punit,logger,dft_section,"PRINT%KPOINTS",error=error)

  END SUBROUTINE write_kpoint_info

! *****************************************************************************
!> \brief Create a single kpoint environment
!> \param kp_env  Single kpoint environment
!> \param error   CP2K error handling
!> \author JGH
! *****************************************************************************
  SUBROUTINE kpoint_env_create(kp_env, error)
    TYPE(kpoint_env_type), POINTER           :: kp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPAssert(.NOT.ASSOCIATED(kp_env),cp_failure_level,routineP,error,failure)

    ALLOCATE(kp_env,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    kp_env%nkpoint = 0
    kp_env%wkp     = 0.0_dp
    kp_env%xkp     = 0.0_dp
    kp_env%is_local= .FALSE.
  
    NULLIFY(kp_env%mos)
    NULLIFY(kp_env%pmat)

  END SUBROUTINE kpoint_env_create

! *****************************************************************************
!> \brief Release a single kpoint environment
!> \param kp_env  Single kpoint environment
!> \param error   CP2K error handling
!> \author JGH
! *****************************************************************************
  SUBROUTINE kpoint_env_release(kp_env, error)
    TYPE(kpoint_env_type), POINTER           :: kp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic, ierr, is
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF(ASSOCIATED(kp_env)) THEN

       IF(ASSOCIATED(kp_env%mos)) THEN
          DO is=1,SIZE(kp_env%mos,2)
             DO ic=1,SIZE(kp_env%mos,1)
                CALL deallocate_mo_set(kp_env%mos(ic,is)%mo_set,error)
             END DO
          END DO
          DEALLOCATE(kp_env%mos,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
       END IF
      
       IF(ASSOCIATED(kp_env%pmat)) THEN
          DO is=1,SIZE(kp_env%pmat,2)
             DO ic=1,SIZE(kp_env%pmat,1)
                CALL cp_fm_release(kp_env%pmat(ic,is)%matrix,error)
             END DO
          END DO
          DEALLOCATE(kp_env%pmat,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
       END IF

       DEALLOCATE(kp_env,STAT=ierr)
       CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
       NULLIFY(kp_env)

    END IF

  END SUBROUTINE kpoint_env_release

! *****************************************************************************
!> \brief Get information from a single kpoint environment
!> \param kpoint_env Single kpoint environment
!> \param nkpoint    Index of kpoint
!> \param wkp        Weight of kpoint
!> \param xkp        Coordinates of kpoint
!> \param is_local   Is this kpoint local (single cpu group)
!> \param mos        MOs of this kpoint
!> \param error      CP2K error handling
!> \author JGH
! *****************************************************************************
  SUBROUTINE get_kpoint_env(kpoint_env, nkpoint, wkp, xkp, is_local, mos, error)
    TYPE(kpoint_env_type), POINTER           :: kpoint_env
    INTEGER, OPTIONAL                        :: nkpoint
    REAL(KIND=dp), OPTIONAL                  :: wkp
    REAL(KIND=dp), DIMENSION(3), OPTIONAL    :: xkp
    LOGICAL, OPTIONAL                        :: is_local
    TYPE(mo_set_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: mos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kpoint_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPAssert(ASSOCIATED(kpoint_env),cp_failure_level,routineP,error,failure)

    IF(PRESENT(nkpoint)) nkpoint = kpoint_env%nkpoint
    IF(PRESENT(wkp)) wkp = kpoint_env%wkp
    IF(PRESENT(xkp)) xkp = kpoint_env%xkp
    IF(PRESENT(is_local)) is_local = kpoint_env%is_local
    IF(PRESENT(mos)) mos => kpoint_env%mos

  END SUBROUTINE get_kpoint_env

! *****************************************************************************

END MODULE kpoint_types
