!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief container for various plainwaves related things
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE pw_env_types
   USE cube_utils,                      ONLY: cube_info_type,&
                                              destroy_cube_info
   USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                              gridlevel_info_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pool_release,&
                                              pw_pool_type,&
                                              pw_pools_dealloc
   USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                              realspace_grid_desc_type,&
                                              realspace_grid_type,&
                                              rs_grid_release,&
                                              rs_grid_release_descriptor
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_env_types'

   PUBLIC :: pw_env_type
   PUBLIC :: pw_env_retain, pw_env_release, pw_env_get

!***

! **************************************************************************************************
!> \brief contained for different pw related things
!> \param pw_pools pools for multigrids
!> \param rs_descs descriptors for the realspace grids
!> \param rs_grids the realspace grids
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
!>        in parallel nsmax is chosen with multigrid in mind!
!> \note
!>      be careful in parallel nsmax for the rs_grids is chosen with
!>      multigrid in mind! (well not as of this writing but...)
!> \par History
!>      10.2002 created [fawzi]
!>      04.2003 added rs grid pools [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE pw_env_type
      INTEGER                                     :: ref_count = -1
      INTEGER                                     :: auxbas_grid = 0
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pw_pools => NULL()
      TYPE(pw_pool_type), POINTER                 :: xc_pw_pool => NULL()
      TYPE(pw_pool_type), POINTER                 :: vdw_pw_pool => NULL()
      TYPE(realspace_grid_desc_p_type), &
         DIMENSION(:), POINTER                    :: rs_descs => NULL()
      TYPE(realspace_grid_type), &
         DIMENSION(:), POINTER                    :: rs_grids => NULL()
      TYPE(gridlevel_info_type), POINTER          :: gridlevel_info => NULL()
      TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info => NULL()
      TYPE(pw_poisson_type), POINTER              :: poisson_env => NULL()
      TYPE(section_vals_type), POINTER            :: interp_section => NULL()
      ! store current cell, used to rebuild lazily.
      REAL(KIND=dp), DIMENSION(3, 3)          :: cell_hmat = -1.0_dp
   END TYPE pw_env_type

CONTAINS

! **************************************************************************************************
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
!> \param xc_pw_pool ...
!> \param vdw_pw_pool ...
!> \param poisson_env ...
!> \param interp_section ...
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_env_get(pw_env, pw_pools, cube_info, gridlevel_info, &
                         auxbas_pw_pool, auxbas_grid, auxbas_rs_desc, auxbas_rs_grid, rs_descs, rs_grids, &
                         xc_pw_pool, vdw_pw_pool, poisson_env, interp_section)
      TYPE(pw_env_type), INTENT(IN)                      :: pw_env
      TYPE(pw_pool_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: pw_pools
      TYPE(cube_info_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: cube_info
      TYPE(gridlevel_info_type), OPTIONAL, POINTER       :: gridlevel_info
      TYPE(pw_pool_type), OPTIONAL, POINTER              :: auxbas_pw_pool
      INTEGER, INTENT(out), OPTIONAL                     :: auxbas_grid
      TYPE(realspace_grid_desc_type), OPTIONAL, POINTER  :: auxbas_rs_desc
      TYPE(realspace_grid_type), OPTIONAL, POINTER       :: auxbas_rs_grid
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: rs_descs
      TYPE(realspace_grid_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: rs_grids
      TYPE(pw_pool_type), OPTIONAL, POINTER              :: xc_pw_pool, vdw_pw_pool
      TYPE(pw_poisson_type), OPTIONAL, POINTER           :: poisson_env
      TYPE(section_vals_type), OPTIONAL, POINTER         :: interp_section

      CPASSERT(pw_env%ref_count > 0)
      IF (PRESENT(pw_pools)) pw_pools => pw_env%pw_pools
      IF (PRESENT(rs_descs)) rs_descs => pw_env%rs_descs
      IF (PRESENT(rs_grids)) rs_grids => pw_env%rs_grids
      IF (PRESENT(cube_info)) cube_info => pw_env%cube_info
      IF (PRESENT(gridlevel_info)) gridlevel_info => pw_env%gridlevel_info
      IF (PRESENT(auxbas_pw_pool)) THEN
         auxbas_pw_pool => pw_env%pw_pools(pw_env%auxbas_grid)%pool
      END IF
      IF (PRESENT(auxbas_rs_desc)) THEN
         auxbas_rs_desc => pw_env%rs_descs(pw_env%auxbas_grid)%rs_desc
      END IF
      IF (PRESENT(auxbas_rs_grid)) THEN
         auxbas_rs_grid => pw_env%rs_grids(pw_env%auxbas_grid)
      END IF
      IF (PRESENT(auxbas_grid)) auxbas_grid = pw_env%auxbas_grid
      IF (PRESENT(xc_pw_pool)) xc_pw_pool => pw_env%xc_pw_pool
      IF (PRESENT(vdw_pw_pool)) vdw_pw_pool => pw_env%vdw_pw_pool
      IF (PRESENT(poisson_env)) poisson_env => pw_env%poisson_env
      IF (PRESENT(interp_section)) interp_section => pw_env%interp_section
   END SUBROUTINE pw_env_get

! **************************************************************************************************
!> \brief retains the pw_env (see doc/ReferenceCounting.html)
!> \param pw_env the pw_env to retain
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_env_retain(pw_env)
      TYPE(pw_env_type), INTENT(INOUT)                   :: pw_env

      CPASSERT(pw_env%ref_count > 0)
      pw_env%ref_count = pw_env%ref_count + 1
   END SUBROUTINE pw_env_retain

! **************************************************************************************************
!> \brief releases the given pw_env (see doc/ReferenceCounting.html)
!> \param pw_env the pw_env to release
!> \param para_env ...
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_env_release(pw_env, para_env)
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(mp_para_env_type), INTENT(IN), OPTIONAL       :: para_env

      INTEGER                                            :: i, igrid_level

      IF (ASSOCIATED(pw_env)) THEN
         CPASSERT(pw_env%ref_count > 0)
         pw_env%ref_count = pw_env%ref_count - 1
         IF (pw_env%ref_count < 1) THEN
            IF (ASSOCIATED(pw_env%poisson_env)) THEN
               CALL pw_env%poisson_env%release()
               DEALLOCATE (pw_env%poisson_env)
            END IF
            CALL pw_pools_dealloc(pw_env%pw_pools)
            IF (ASSOCIATED(pw_env%gridlevel_info)) THEN
               CALL destroy_gaussian_gridlevel(pw_env%gridlevel_info, para_env)
               DEALLOCATE (pw_env%gridlevel_info)
            END IF
            IF (ASSOCIATED(pw_env%cube_info)) THEN
               DO igrid_level = 1, SIZE(pw_env%cube_info)
                  CALL destroy_cube_info(pw_env%cube_info(igrid_level))
               END DO
               DEALLOCATE (pw_env%cube_info)
            END IF
            CALL pw_pool_release(pw_env%xc_pw_pool)
            CALL pw_pool_release(pw_env%vdw_pw_pool)
            IF (ASSOCIATED(pw_env%rs_descs)) THEN
               DO i = 1, SIZE(pw_env%rs_descs)
                  CALL rs_grid_release_descriptor(pw_env%rs_descs(i)%rs_desc)
               END DO
               DEALLOCATE (pw_env%rs_descs)
            END IF
            IF (ASSOCIATED(pw_env%rs_grids)) THEN
               DO i = 1, SIZE(pw_env%rs_grids)
                  CALL rs_grid_release(pw_env%rs_grids(i))
               END DO
               DEALLOCATE (pw_env%rs_grids)
            END IF
            DEALLOCATE (pw_env)
         END IF
      END IF
      NULLIFY (pw_env)
   END SUBROUTINE pw_env_release

END MODULE pw_env_types
