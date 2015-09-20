!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types and basic routines needed for a kpoint calculation
!> \par History
!>       2014.07 created [JGH]
!>       2014.11 unified k-point and gamma-point code [Ole Schuett]
!> \author JGH
! *****************************************************************************
MODULE kpoint_types
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_env,                     ONLY: cp_cart_release,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_cart_type,&
                                             cp_para_env_type
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
  PUBLIC :: kpoint_create, kpoint_release, get_kpoint_info
  PUBLIC :: read_kpoint_section, write_kpoint_info
  PUBLIC :: kpoint_env_type, kpoint_env_p_type
  PUBLIC :: kpoint_env_create


! *****************************************************************************  
!> \brief Keeps information about a specific k-point
!> \param nkpoint   the kpoint index
!> \param wkp       weight of this kpoint
!> \param xkp       kpoint coordinates
!> \param is_local  if this kpoint is calculated on a single thread
!> \param mos       associated MOs (r/i,spin)
!> \param pmat      associated density matrix (r/i,spin)
!> \param wmat      associated energy weighted density matrix (r/i,spin)
!> \author JGH
! *****************************************************************************
  TYPE kpoint_env_type
     INTEGER                                           :: nkpoint
     REAL(KIND=dp)                                     :: wkp
     REAL(KIND=dp), DIMENSION(3)                       :: xkp
     LOGICAL                                           :: is_local
     TYPE(mo_set_p_type), DIMENSION(:,:), POINTER      :: mos
     TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER       :: pmat
     TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER       :: wmat
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
     REAL(KIND=dp), DIMENSION(:,:), POINTER  :: xkp => Null()
     REAL(KIND=dp), DIMENSION(:), POINTER    :: wkp => Null()
     ! parallel environment
     TYPE(cp_para_env_type), POINTER         :: para_env => Null()
     TYPE(cp_blacs_env_type), POINTER        :: blacs_env_all => Null()
     TYPE(cp_para_cart_type), POINTER        :: cart => Null()
     TYPE(cp_para_env_type), POINTER         :: para_env_full => Null(),&
                                                para_env_kp => Null(),&
                                                para_env_inter_kp => Null()
     LOGICAL                                 :: iogrp
     INTEGER                                 :: nkp_groups
     INTEGER, DIMENSION(:,:), POINTER        :: kp_dist => Null()
     INTEGER, DIMENSION(2)                   :: kp_range
     TYPE(cp_blacs_env_type), POINTER        :: blacs_env => Null()
     INTEGER, DIMENSION(:,:,:), POINTER      :: cell_to_index => Null()
     TYPE(neighbor_list_set_p_type), &
       DIMENSION(:), POINTER                 :: sab_nl => Null()
     ! environment
     TYPE(kpoint_env_p_type), DIMENSION(:),&
       POINTER                               :: kp_env => Null()
     TYPE(qs_matrix_pools_type), POINTER     :: mpools => Null()
  END TYPE kpoint_type

! *****************************************************************************  

CONTAINS

! *****************************************************************************
!> \brief Create a kpoint environment
!> \param kpoint  All the kpoint information
!> \author JGH
! *****************************************************************************
  SUBROUTINE kpoint_create(kpoint)
    TYPE(kpoint_type), POINTER               :: kpoint

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(.NOT.ASSOCIATED(kpoint))

    ALLOCATE(kpoint,STAT=ierr)
    CPASSERT(ierr==0)

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
    NULLIFY(kpoint%kp_dist)

    NULLIFY(kpoint%para_env)
    NULLIFY(kpoint%blacs_env_all)
    NULLIFY(kpoint%cart)
    NULLIFY(kpoint%para_env_full,kpoint%para_env_kp,kpoint%para_env_inter_kp)
    NULLIFY(kpoint%blacs_env)
    kpoint%nkp_groups = 0
    kpoint%iogrp = .FALSE.
    kpoint%kp_range = 0


    NULLIFY(kpoint%kp_env)
    NULLIFY(kpoint%mpools)

    ALLOCATE(kpoint%cell_to_index(0:0,0:0,0:0))
    kpoint%cell_to_index(:,:,:) = 1
  END SUBROUTINE kpoint_create

! *****************************************************************************
!> \brief  Release a kpoint environment, deallocate all data
!> \param kpoint  The kpoint environment
!> \author JGH
! *****************************************************************************
  SUBROUTINE kpoint_release(kpoint)
    TYPE(kpoint_type), POINTER               :: kpoint

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, ik
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF(ASSOCIATED(kpoint)) THEN

       IF(ASSOCIATED(kpoint%xkp)) THEN
         DEALLOCATE(kpoint%xkp,STAT=ierr)
         CPASSERT(ierr==0)
         NULLIFY(kpoint%xkp)
       END IF
       IF(ASSOCIATED(kpoint%wkp)) THEN
         DEALLOCATE(kpoint%wkp,STAT=ierr)
         CPASSERT(ierr==0)
         NULLIFY(kpoint%wkp)
       END IF
       IF(ASSOCIATED(kpoint%kp_dist)) THEN
         DEALLOCATE(kpoint%kp_dist,STAT=ierr)
         CPASSERT(ierr==0)
       ENDIF

       CALL mpools_release(kpoint%mpools)
       CALL cp_blacs_env_release(kpoint%blacs_env)
       CALL cp_blacs_env_release(kpoint%blacs_env_all)

       CALL cp_cart_release(kpoint%cart)
       CALL cp_para_env_release(kpoint%para_env)
       CALL cp_para_env_release(kpoint%para_env_full)
       CALL cp_para_env_release(kpoint%para_env_kp)
       CALL cp_para_env_release(kpoint%para_env_inter_kp)

       IF(ASSOCIATED(kpoint%cell_to_index)) DEALLOCATE(kpoint%cell_to_index)

       IF(ASSOCIATED(kpoint%kp_env)) THEN
          DO ik=1,SIZE(kpoint%kp_env)
             CALL kpoint_env_release(kpoint%kp_env(ik)%kpoint_env)
          END DO
          DEALLOCATE(kpoint%kp_env,STAT=ierr)
          CPASSERT(ierr==0)
       END IF

       DEALLOCATE(kpoint,STAT=ierr)
       CPASSERT(ierr==0)
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
!> \param blacs_env_all BLACS env for the total environment
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
!> \param cell_to_index given a cell triple, returns the real space index
!> \param sab_nl        neighbourlist that defines real space matrices
!> \author JGH
! *****************************************************************************
  SUBROUTINE get_kpoint_info(kpoint, kp_scheme, nkp_grid, kp_shift, symmetry, verbose,&
      full_grid, use_real_wfn, eps_geo, parallel_group_size, kp_range, nkp, xkp, wkp,&
      para_env, blacs_env_all, cart, para_env_full, para_env_kp, para_env_inter_kp, blacs_env,&
      kp_env, mpools, iogrp, nkp_groups, kp_dist, cell_to_index, sab_nl)
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
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env_all
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
      POINTER                                :: kp_dist
    INTEGER, DIMENSION(:, :, :), OPTIONAL, &
      POINTER                                :: cell_to_index
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sab_nl

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kpoint_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(ASSOCIATED(kpoint))

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
    IF(PRESENT(blacs_env_all)) blacs_env_all => kpoint%blacs_env_all
    IF(PRESENT(blacs_env)) blacs_env => kpoint%blacs_env

    IF(PRESENT(iogrp)) iogrp = kpoint%iogrp
    IF(PRESENT(kp_range)) kp_range = kpoint%kp_range
    IF(PRESENT(nkp_groups)) nkp_groups = kpoint%nkp_groups
    IF(PRESENT(kp_dist)) kp_dist => kpoint%kp_dist

    IF(PRESENT(kp_env)) kp_env => kpoint%kp_env
    IF(PRESENT(mpools)) mpools => kpoint%mpools


    IF(PRESENT(cell_to_index)) cell_to_index => kpoint%cell_to_index
    IF(PRESENT(sab_nl)) sab_nl => kpoint%sab_nl

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
!> \param blacs_env_all BLACS env for the total environment
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
!> \param cell_to_index given a cell triple, returns the real space index
!> \param sab_nl        neighbourlist that defines real space matrices
!> \author JGH
! *****************************************************************************
  SUBROUTINE set_kpoint_info(kpoint, kp_scheme, nkp_grid, kp_shift, symmetry, verbose,&
      full_grid, use_real_wfn, eps_geo, parallel_group_size, kp_range, nkp, xkp, wkp,&
      para_env, blacs_env_all, cart, para_env_full, para_env_kp, para_env_inter_kp, blacs_env,&
      kp_env, mpools, iogrp, nkp_groups, kp_dist, cell_to_index, sab_nl)
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
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env_all
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
      POINTER                                :: kp_dist
    INTEGER, DIMENSION(:, :, :), OPTIONAL, &
      POINTER                                :: cell_to_index
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sab_nl

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_kpoint_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(ASSOCIATED(kpoint))
 
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
    IF(PRESENT(blacs_env_all)) kpoint%blacs_env_all => blacs_env_all
    IF(PRESENT(blacs_env)) kpoint%blacs_env => blacs_env

    IF(PRESENT(iogrp)) kpoint%iogrp = iogrp
    IF(PRESENT(kp_range)) kpoint%kp_range = kp_range
    IF(PRESENT(nkp_groups)) kpoint%nkp_groups = nkp_groups
    IF(PRESENT(kp_dist)) kpoint%kp_dist => kp_dist

    IF(PRESENT(kp_env)) kpoint%kp_env => kp_env
    IF(PRESENT(mpools)) kpoint%mpools => mpools
    IF(PRESENT(sab_nl)) kpoint%sab_nl => sab_nl

    IF(PRESENT(cell_to_index)) THEN
       IF(ASSOCIATED(kpoint%cell_to_index)) DEALLOCATE(kpoint%cell_to_index)
       kpoint%cell_to_index => cell_to_index
    ENDIF

  END SUBROUTINE set_kpoint_info

! *****************************************************************************
!> \brief Read the kpoint input section
!> \param kpoint  The kpoint environment
!> \param kpoint_section The input section
!> \author JGH
! *****************************************************************************
  SUBROUTINE read_kpoint_section(kpoint,kpoint_section)
    TYPE(kpoint_type), POINTER               :: kpoint
    TYPE(section_vals_type), POINTER         :: kpoint_section

    CHARACTER(len=*), PARAMETER :: routineN = 'read_kpoint_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER                                  :: i, ierr, n_rep, nval, wfntype
    LOGICAL                                  :: available, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: reallist

    failure = .FALSE.
    CPASSERT(ASSOCIATED(kpoint))

    CALL section_vals_get(kpoint_section,explicit=available)

    IF (available) THEN
       CALL section_vals_val_get(kpoint_section,"SCHEME",c_vals=tmpstringlist)
       nval = SIZE(tmpstringlist)
       CPASSERT(nval >= 1)
       kpoint%kp_scheme = tmpstringlist(1)
       CALL uppercase(kpoint%kp_scheme)

       ! SCHEME [None, Gamma, Monkhorst-Pack, MacDonald, General]
       SELECT CASE (kpoint%kp_scheme)
          CASE ("NONE")
             ! do nothing
          CASE ("GAMMA")
             ! do nothing
          CASE ("MONKHORST-PACK")
             CPASSERT(nval >= 4)
             DO i=2,4
                READ(tmpstringlist(i),*) kpoint%nkp_grid(i-1)
             END DO
          CASE ("MACDONALD")
             CPASSERT(nval >= 7)
             DO i=2,4
                READ(tmpstringlist(i),*) kpoint%nkp_grid(i-1)
             END DO
             DO i=5,7
                READ(tmpstringlist(i),*) kpoint%kp_shift(i-4)
             END DO
          CASE ("GENERAL")
             CALL section_vals_val_get(kpoint_section,"KPOINT", n_rep_val=n_rep)
             kpoint%nkp = n_rep
             ALLOCATE(kpoint%xkp(3,n_rep),kpoint%wkp(n_rep),STAT=ierr)
             CPASSERT(ierr==0)
             DO i=1,n_rep
                CALL section_vals_val_get(kpoint_section,"KPOINT", i_rep_val=i,&
                                  r_vals=reallist)
                nval = SIZE(reallist)
                CPASSERT(nval>=4)
                kpoint%xkp(1:3,i) = reallist(1:3)
                kpoint%wkp(i) = reallist(4)
             END DO
          CASE DEFAULT
            CPASSERT(.FALSE.)
       END SELECT

       CALL section_vals_val_get(kpoint_section,"SYMMETRY",l_val=kpoint%symmetry)
       CALL section_vals_val_get(kpoint_section,"WAVEFUNCTIONS",i_val=wfntype)
       CALL section_vals_val_get(kpoint_section,"VERBOSE",l_val=kpoint%verbose)
       CALL section_vals_val_get(kpoint_section,"FULL_GRID",l_val=kpoint%full_grid)
       CALL section_vals_val_get(kpoint_section,"EPS_GEO",r_val=kpoint%eps_geo)
       CALL section_vals_val_get(kpoint_section,"PARALLEL_GROUP_SIZE",&
                                 i_val=kpoint%parallel_group_size)
       SELECT CASE (wfntype)
          CASE (use_real_wfn)
             kpoint%use_real_wfn = .TRUE.
          CASE (use_complex_wfn)
             kpoint%use_real_wfn = .FALSE.
          CASE DEFAULT
            CPASSERT(.FALSE.)
       END SELECT

    ELSE
       kpoint%kp_scheme = "NONE"
    END IF

  END SUBROUTINE read_kpoint_section

! *****************************************************************************
!> \brief Write information on the kpoints to output
!> \param kpoint  The kpoint environment
!> \param dft_section  DFT section information
!> \author JGH
! *****************************************************************************
  SUBROUTINE write_kpoint_info(kpoint,dft_section)
    TYPE(kpoint_type), POINTER               :: kpoint
    TYPE(section_vals_type), POINTER         :: dft_section

    CHARACTER(len=*), PARAMETER :: routineN = 'write_kpoint_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, punit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    CPASSERT(ASSOCIATED(kpoint))

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_get_default_logger()

    punit = cp_print_key_unit_nr(logger,dft_section,"PRINT%KPOINTS",extension=".Log")
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
            CPASSERT(.FALSE.)
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
    CALL cp_print_key_finished_output(punit,logger,dft_section,"PRINT%KPOINTS")

  END SUBROUTINE write_kpoint_info

! *****************************************************************************
!> \brief Create a single kpoint environment
!> \param kp_env  Single kpoint environment
!> \author JGH
! *****************************************************************************
  SUBROUTINE kpoint_env_create(kp_env)
    TYPE(kpoint_env_type), POINTER           :: kp_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(.NOT.ASSOCIATED(kp_env))

    ALLOCATE(kp_env,STAT=ierr)
    CPASSERT(ierr==0)

    kp_env%nkpoint = 0
    kp_env%wkp     = 0.0_dp
    kp_env%xkp     = 0.0_dp
    kp_env%is_local= .FALSE.
  
    NULLIFY(kp_env%mos)
    NULLIFY(kp_env%pmat)
    NULLIFY(kp_env%wmat)

  END SUBROUTINE kpoint_env_create

! *****************************************************************************
!> \brief Release a single kpoint environment
!> \param kp_env  Single kpoint environment
!> \author JGH
! *****************************************************************************
  SUBROUTINE kpoint_env_release(kp_env)
    TYPE(kpoint_env_type), POINTER           :: kp_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic, ierr, is
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF(ASSOCIATED(kp_env)) THEN

       IF(ASSOCIATED(kp_env%mos)) THEN
          DO is=1,SIZE(kp_env%mos,2)
             DO ic=1,SIZE(kp_env%mos,1)
                CALL deallocate_mo_set(kp_env%mos(ic,is)%mo_set)
             END DO
          END DO
          DEALLOCATE(kp_env%mos,STAT=ierr)
          CPASSERT(ierr==0)
       END IF

       IF(ASSOCIATED(kp_env%pmat)) THEN
          DO is=1,SIZE(kp_env%pmat,2)
             DO ic=1,SIZE(kp_env%pmat,1)
                CALL cp_fm_release(kp_env%pmat(ic,is)%matrix)
             END DO
          END DO
          DEALLOCATE(kp_env%pmat,STAT=ierr)
          CPASSERT(ierr==0)
       END IF

       IF(ASSOCIATED(kp_env%wmat)) THEN
          DO is=1,SIZE(kp_env%wmat,2)
             DO ic=1,SIZE(kp_env%wmat,1)
                CALL cp_fm_release(kp_env%wmat(ic,is)%matrix)
             END DO
          END DO
          DEALLOCATE(kp_env%wmat,STAT=ierr)
          CPASSERT(ierr==0)
       END IF

       DEALLOCATE(kp_env,STAT=ierr)
       CPASSERT(ierr==0)
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
!> \author JGH
! *****************************************************************************
  SUBROUTINE get_kpoint_env(kpoint_env, nkpoint, wkp, xkp, is_local, mos)
    TYPE(kpoint_env_type), POINTER           :: kpoint_env
    INTEGER, OPTIONAL                        :: nkpoint
    REAL(KIND=dp), OPTIONAL                  :: wkp
    REAL(KIND=dp), DIMENSION(3), OPTIONAL    :: xkp
    LOGICAL, OPTIONAL                        :: is_local
    TYPE(mo_set_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: mos

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kpoint_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(ASSOCIATED(kpoint_env))

    IF(PRESENT(nkpoint)) nkpoint = kpoint_env%nkpoint
    IF(PRESENT(wkp)) wkp = kpoint_env%wkp
    IF(PRESENT(xkp)) xkp = kpoint_env%xkp
    IF(PRESENT(is_local)) is_local = kpoint_env%is_local
    IF(PRESENT(mos)) mos => kpoint_env%mos

  END SUBROUTINE get_kpoint_env

! *****************************************************************************

END MODULE kpoint_types
