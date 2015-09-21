!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief container for various plainwaves related things
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE pw_env_types
  USE cube_utils,                      ONLY: cube_info_type,&
                                             destroy_cube_info
  USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                             gridlevel_info_type
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE lgrid_types,                     ONLY: lgrid_release,&
                                             lgrid_type
  USE pw_poisson_types,                ONLY: pw_poisson_release,&
                                             pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                             pw_pool_release,&
                                             pw_pool_type,&
                                             pw_pools_dealloc,&
                                             pw_pools_flush_cache
  USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                             realspace_grid_desc_type,&
                                             realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs_grid_release,&
                                             rs_grid_release_descriptor
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_env_types'

  PUBLIC :: pw_env_type
  PUBLIC :: pw_env_retain, pw_env_release, pw_env_get
!***

! *****************************************************************************
!> \brief contained for different pw related things
!> \param pw_pools pools for multigrids
!> \param rs_descs descriptors for the realspace grids
!> \param rs_grids the realspace grids
!> \param lgrid scratch space for OpenMP threads rs_grid operations
!> \param xc_pw_pool the pool used for xc calculation on grid (might be the
!>        same as the finest, or finer still).
!> \param vdw_pw_pool the pool used for vdw-nl calculation on grid (might be the
!>        same or a grid with a new cutoff
!> \param gridlevel_info gaussian gridlevel info
!> \param cube_info info on cube and max sphere size for multigrids
!> \param aux_bas which multigrid is the auxiliar basis
!>
!>      readonly attributes
!> \param auxbas_pw_pool a pool that allocates grids in the auxiliary basis
!> \param auxbas_rs_desc real space grid in the auxiliary basis, be careful
!>        in parallel nsmax is choosen with multigrid in mind!
!> \note
!>      be careful in parallel nsmax for the rs_grids is choosen with
!>      multigrid in mind! (well not as of this writing but...)
!> \par History
!>      10.2002 created [fawzi]
!>      04.2003 added rs grid pools [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE pw_env_type
     INTEGER                                     :: ref_count
     INTEGER                                     :: auxbas_grid
     TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pw_pools
     TYPE(pw_pool_type), POINTER                 :: xc_pw_pool
     TYPE(pw_pool_type), POINTER                 :: vdw_pw_pool
     TYPE(realspace_grid_desc_p_type), &
        DIMENSION(:), POINTER                    :: rs_descs
     TYPE(realspace_grid_p_type), &
        DIMENSION(:), POINTER                    :: rs_grids
     TYPE(lgrid_type), POINTER                   :: lgrid
     TYPE(gridlevel_info_type), POINTER          :: gridlevel_info
     TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
     TYPE(pw_poisson_type), POINTER              :: poisson_env
     TYPE(section_vals_type), POINTER            :: interp_section
     ! store current cell, used to rebuild lazily.
     REAL(KIND = dp), DIMENSION(3,3)          :: cell_hmat = -1.0_dp
  END TYPE pw_env_type

CONTAINS

! *****************************************************************************
!> \brief returns the various attributes of the pw env
!> \param pw_env the einvironment you want info about
!> \param pw_pools ...
!> \param cube_info ...
!> \param gridlevel_info ...
!> \param auxbas_pw_pool ...
!> \param auxbas_grid ...
!> \param auxbas_rs_desc ...
!> \param auxbas_rs_grid ...
!> \param rs_descs ...
!> \param rs_grids ...
!> \param lgrid ...
!> \param xc_pw_pool ...
!> \param vdw_pw_pool ...
!> \param poisson_env ...
!> \param interp_section ...
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE pw_env_get(pw_env, pw_pools, cube_info, gridlevel_info,&
     auxbas_pw_pool,auxbas_grid,auxbas_rs_desc,auxbas_rs_grid,rs_descs,rs_grids,lgrid,&
     xc_pw_pool, vdw_pw_pool, poisson_env, interp_section)
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: pw_pools
    TYPE(cube_info_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: cube_info
    TYPE(gridlevel_info_type), OPTIONAL, &
      POINTER                                :: gridlevel_info
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: auxbas_pw_pool
    INTEGER, INTENT(out), OPTIONAL           :: auxbas_grid
    TYPE(realspace_grid_desc_type), &
      OPTIONAL, POINTER                      :: auxbas_rs_desc
    TYPE(realspace_grid_type), OPTIONAL, &
      POINTER                                :: auxbas_rs_grid
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: rs_grids
    TYPE(lgrid_type), OPTIONAL, POINTER      :: lgrid
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: xc_pw_pool, vdw_pw_pool
    TYPE(pw_poisson_type), OPTIONAL, POINTER :: poisson_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: interp_section

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(pw_env))
  CPASSERT(pw_env%ref_count>0)
  IF (PRESENT(pw_pools)) pw_pools => pw_env%pw_pools
  IF (PRESENT(rs_descs)) rs_descs => pw_env%rs_descs
  IF (PRESENT(rs_grids)) rs_grids => pw_env%rs_grids
  IF (PRESENT(lgrid)) lgrid => pw_env%lgrid
  IF (PRESENT(cube_info)) cube_info => pw_env%cube_info
  IF (PRESENT(gridlevel_info)) gridlevel_info => pw_env%gridlevel_info
  IF (PRESENT(auxbas_pw_pool)) THEN
     auxbas_pw_pool => pw_env%pw_pools(pw_env%auxbas_grid)%pool
  END IF
  IF (PRESENT(auxbas_rs_desc)) THEN
     auxbas_rs_desc => pw_env%rs_descs(pw_env%auxbas_grid)%rs_desc
  END IF
  IF (PRESENT(auxbas_rs_grid)) THEN
     auxbas_rs_grid => pw_env%rs_grids(pw_env%auxbas_grid)%rs_grid
  END IF
  IF (PRESENT(auxbas_grid)) auxbas_grid = pw_env%auxbas_grid
  IF (PRESENT(xc_pw_pool)) xc_pw_pool => pw_env%xc_pw_pool
  IF (PRESENT(vdw_pw_pool)) vdw_pw_pool => pw_env%vdw_pw_pool
  IF (PRESENT(poisson_env)) poisson_env => pw_env%poisson_env
  IF (PRESENT(interp_section)) interp_section => pw_env%interp_section
END SUBROUTINE pw_env_get

! *****************************************************************************
!> \brief retains the pw_env (see doc/ReferenceCounting.html)
!> \param pw_env the pw_env to retain
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE pw_env_retain(pw_env)
    TYPE(pw_env_type), POINTER               :: pw_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(pw_env))
  CPASSERT(pw_env%ref_count>0)
  pw_env%ref_count=pw_env%ref_count+1
END SUBROUTINE pw_env_retain

! *****************************************************************************
!> \brief releases the given pw_env (see doc/ReferenceCounting.html)
!> \param pw_env the pw_env to release
!> \param kg ...
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE pw_env_release(pw_env, kg)
    TYPE(pw_env_type), POINTER               :: pw_env
    LOGICAL, OPTIONAL                        :: kg

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, igrid_level
    LOGICAL                                  :: failure, my_kg

  failure=.FALSE.
  my_kg = .FALSE.
  IF(PRESENT(kg)) my_kg = kg
  IF (ASSOCIATED(pw_env)) THEN
     CPASSERT(pw_env%ref_count>0)
     pw_env%ref_count=pw_env%ref_count-1
     IF (pw_env%ref_count<1) THEN
        CALL pw_poisson_release(pw_env%poisson_env)
        CALL pw_pools_dealloc(pw_env%pw_pools)
        IF (ASSOCIATED(pw_env%gridlevel_info)) THEN
           IF(my_kg) THEN
             DEALLOCATE (pw_env%gridlevel_info%cutoff)
             CALL section_vals_release(pw_env%gridlevel_info%print_section)
             DEALLOCATE (pw_env%gridlevel_info%count)
           ELSE
             CALL destroy_gaussian_gridlevel(pw_env%gridlevel_info)
           END IF
           DEALLOCATE(pw_env%gridlevel_info)
        END IF
        IF (ASSOCIATED(pw_env%cube_info)) THEN
           DO igrid_level=1,SIZE(pw_env%cube_info)
              CALL destroy_cube_info(pw_env%cube_info(igrid_level))
           END DO
           DEALLOCATE(pw_env%cube_info)
        END IF
        CALL pw_pool_release(pw_env%xc_pw_pool)
        CALL pw_pool_release(pw_env%vdw_pw_pool)
        IF (ASSOCIATED(pw_env%rs_descs)) THEN
          DO i=1, SIZE(pw_env%rs_descs)
            CALL rs_grid_release_descriptor(pw_env%rs_descs(i)%rs_desc)
          END DO
          DEALLOCATE(pw_env%rs_descs)
        END IF
        IF (ASSOCIATED(pw_env%rs_grids)) THEN
          DO i=1, SIZE(pw_env%rs_grids)
            CALL rs_grid_release(pw_env%rs_grids(i)%rs_grid)
          END DO
          DEALLOCATE(pw_env%rs_grids)
        END IF
        IF (ASSOCIATED(pw_env%lgrid)) THEN
          CALL lgrid_release(pw_env%lgrid)
        END IF
        DEALLOCATE(pw_env)
    END IF
  END IF
  NULLIFY(pw_env)
END SUBROUTINE pw_env_release

! *****************************************************************************
!> \brief flushes the cached pws
!> \param pw_env the pw_env to be flushed
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE pw_env_flush_cache(pw_env)
    TYPE(pw_env_type), POINTER               :: pw_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_flush_cache', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(pw_env))
  IF (ASSOCIATED(pw_env%pw_pools)) THEN
     CALL pw_pools_flush_cache(pw_env%pw_pools)
  END IF
END SUBROUTINE pw_env_flush_cache

END MODULE pw_env_types
