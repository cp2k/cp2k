!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Transfers densities from PW to RS grids and potentials from PW to RS  
!> \par History
!>      - Copied from qs_coolocate_Density and qs_integrate_potenntial
!> \author JGH (04.2014)
! *****************************************************************************
MODULE rs_pw_interface
  USE cp_spline_utils,                 ONLY: pw_prolongate_s3,&
                                             pw_restrict_s3
  USE gaussian_gridlevels,             ONLY: gridlevel_info_type
  USE input_constants,                 ONLY: pw_interp,&
                                             spline3_pbc_interp
  USE input_section_types,             ONLY: section_vals_val_get
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                             pw_pools_create_pws,&
                                             pw_pools_give_back_pws
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE realspace_grid_types,            ONLY: pw2rs,&
                                             realspace_grid_desc_p_type,&
                                             realspace_grid_p_type,&
                                             rs2pw,&
                                             rs_grid_release,&
                                             rs_pw_transfer
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rs_pw_interface'
! *** Public subroutines ***

  PUBLIC :: density_rs2pw,&
            density_rs2pw_basic,&
            potential_pw2rs

CONTAINS

! *****************************************************************************
!> \brief given partial densities on the realspace multigrids,
!>      computes the full density on the plane wave grids, both in real and
!>      gspace
!> \param pw_env ...
!> \param rs_rho ...
!> \param rho ...
!> \param rho_gspace ...
!> \param error ...
!> \note
!>      should contain all communication in the collocation of the density
!>      in the case of replicated grids
! *****************************************************************************
  SUBROUTINE density_rs2pw(pw_env,rs_rho,rho,rho_gspace,error)

    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    TYPE(pw_p_type), INTENT(INOUT)           :: rho, rho_gspace
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'density_rs2pw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, igrid_level, &
                                                interp_kind
    LOGICAL                                  :: failure
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: mgrid_gspace, mgrid_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(gridlevel_info, mgrid_gspace, mgrid_rspace, rs_descs, pw_pools)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_descs=rs_descs, pw_pools=pw_pools, error=error)

    gridlevel_info=>pw_env%gridlevel_info

    CALL section_vals_val_get(pw_env%interp_section,"KIND",i_val=interp_kind,error=error)

    CALL pw_pools_create_pws(pw_pools,mgrid_rspace,&
                              use_data = REALDATA3D,&
                              in_space = REALSPACE, error=error)

    CALL pw_pools_create_pws(pw_pools,mgrid_gspace,&
                              use_data = COMPLEXDATA1D,&
                              in_space = RECIPROCALSPACE, error=error)

   IF (gridlevel_info%ngrid_levels==1) THEN
       CALL rs_pw_transfer(rs_rho(1)%rs_grid,rho%pw,rs2pw,error=error)
       CALL rs_grid_release(rs_rho(1)%rs_grid, error=error)
       CALL pw_transfer(rho%pw,rho_gspace%pw,error=error)
       IF (rho%pw%pw_grid%spherical) THEN ! rho_gspace = rho
          CALL pw_transfer(rho_gspace%pw,rho%pw,error=error)
       ENDIF
    ELSE
       DO igrid_level=1,gridlevel_info%ngrid_levels
          CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid,&
               mgrid_rspace(igrid_level)%pw,rs2pw,error=error)
          CALL rs_grid_release(rs_rho(igrid_level)%rs_grid, error=error)
       ENDDO

       ! we want both rho and rho_gspace, the latter for Hartree and co-workers.
       SELECT CASE(interp_kind)
       CASE(pw_interp)
          CALL pw_zero(rho_gspace%pw,error=error)
          DO igrid_level=1,gridlevel_info%ngrid_levels
             CALL pw_transfer(mgrid_rspace(igrid_level)%pw,&
                  mgrid_gspace(igrid_level)%pw,error=error)
             CALL pw_axpy(mgrid_gspace(igrid_level)%pw,rho_gspace%pw,&
                  error=error)
          END DO
          CALL pw_transfer(rho_gspace%pw,rho%pw,error=error)
       CASE(spline3_pbc_interp)
          DO igrid_level=gridlevel_info%ngrid_levels,2,-1
             CALL pw_prolongate_s3(mgrid_rspace(igrid_level)%pw,&
                  mgrid_rspace(igrid_level-1)%pw,pw_pools(igrid_level)%pool,&
                  pw_env%interp_section,error=error)
          END DO
          CALL pw_copy(mgrid_rspace(1)%pw,rho%pw,error=error)
          CALL pw_transfer(rho%pw,rho_gspace%pw,error=error)
       CASE default
          CALL cp_unimplemented_error(routineN,"interpolator "//&
               cp_to_string(interp_kind),error=error)
       END SELECT
    END IF

    ! *** give back the pw multi-grids
    CALL pw_pools_give_back_pws(pw_pools,mgrid_gspace,error=error)
    CALL pw_pools_give_back_pws(pw_pools,mgrid_rspace,error=error)
    CALL timestop(handle)

  END SUBROUTINE density_rs2pw

! *****************************************************************************
!> \brief given partial densities on the realspace multigrids,
!>      computes the full density on the plane wave grids
!> \param pw_env ...
!> \param rs_rho ...
!> \param rho ...
!> \param rho_gspace ...
!> \param error ...
!> \note
!>      should contain the all communication in the collocation of the density
!>      in the case of replicated grids
! *****************************************************************************
  SUBROUTINE density_rs2pw_basic(pw_env,rs_rho,rho,rho_gspace,error)

    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    TYPE(pw_p_type), INTENT(INOUT)           :: rho, rho_gspace
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'density_rs2pw_basic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, igrid_level, &
                                                interp_kind
    LOGICAL                                  :: failure
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: mgrid_gspace, mgrid_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(gridlevel_info, mgrid_gspace, mgrid_rspace, rs_descs, pw_pools)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_descs=rs_descs, pw_pools=pw_pools, error=error)


    gridlevel_info=>pw_env%gridlevel_info

    CALL section_vals_val_get(pw_env%interp_section,"KIND",i_val=interp_kind,error=error)

    CALL pw_pools_create_pws(pw_pools,mgrid_rspace,&
                              use_data = REALDATA3D,&
                              in_space = REALSPACE, error=error)

    CALL pw_pools_create_pws(pw_pools,mgrid_gspace,&
                              use_data = COMPLEXDATA1D,&
                              in_space = RECIPROCALSPACE, error=error)

    IF (gridlevel_info%ngrid_levels==1) THEN
       CALL rs_pw_transfer(rs_rho(1)%rs_grid,rho%pw,rs2pw,error=error)
       CALL pw_transfer(rho%pw,rho_gspace%pw,error=error)
    ELSE
       DO igrid_level=1,gridlevel_info%ngrid_levels
         CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid,&
               mgrid_rspace(igrid_level)%pw,rs2pw,error=error)
       ENDDO

       ! we want both rho and rho_gspace, the latter for Hartree and co-workers.
       SELECT CASE(interp_kind)
       CASE(pw_interp)
          DO igrid_level=1,gridlevel_info%ngrid_levels
             CALL pw_transfer(mgrid_rspace(igrid_level)%pw,&
                  mgrid_gspace(igrid_level)%pw,error=error)
             IF (igrid_level/=1) THEN
               CALL pw_axpy(mgrid_gspace(igrid_level)%pw,mgrid_gspace(1)%pw,&
                    error=error)
             END IF
          END DO
          CALL pw_transfer(mgrid_gspace(1)%pw,rho%pw,error=error)
          CALL pw_transfer(mgrid_rspace(1)%pw,rho_gspace%pw,error=error)
       CASE(spline3_pbc_interp)
          DO igrid_level=gridlevel_info%ngrid_levels,2,-1
             CALL pw_prolongate_s3(mgrid_rspace(igrid_level)%pw,&
                  mgrid_rspace(igrid_level-1)%pw,pw_pools(igrid_level)%pool,&
                  pw_env%interp_section,error=error)
          END DO
          CALL pw_copy(mgrid_rspace(1)%pw,rho%pw,error=error)
          CALL pw_transfer(rho%pw,rho_gspace%pw,error=error)
       CASE default
          CALL cp_unimplemented_error(routineN,"interpolator "//&
               cp_to_string(interp_kind),error=error)
       END SELECT
    END IF

    ! *** give back the pw multi-grids
    CALL pw_pools_give_back_pws(pw_pools,mgrid_gspace,error=error)
    CALL pw_pools_give_back_pws(pw_pools,mgrid_rspace,error=error)
    CALL timestop(handle)

  END SUBROUTINE density_rs2pw_basic


! *****************************************************************************
!> \brief transfers a potential from a pw_grid to a vector of
!>      realspace multigrids
!> \param rs_v OUTPUT: the potential on the realspace multigrids
!> \param v_rspace INPUT : the potential on a planewave grid in Rspace
!> \param pw_env ...
!> \param error ...
!> \par History
!>      09.2006 created [Joost VandeVondele]
!> \note
!>      extracted from integrate_v_rspace
!>      should contain all parallel communication of integrate_v_rspace in the
!>      case of replicated grids.
! *****************************************************************************
  SUBROUTINE potential_pw2rs(rs_v,v_rspace,pw_env,error)

    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_v
    TYPE(pw_p_type), INTENT(IN)              :: v_rspace
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'potential_pw2rs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: auxbas_grid, handle, &
                                                igrid_level, interp_kind
    REAL(KIND=dp)                            :: scale
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: mgrid_gspace, mgrid_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools

    CALL timeset(routineN,handle)

    ! *** set up of the potential on the multigrids
    CALL pw_env_get(pw_env, pw_pools=pw_pools, gridlevel_info=gridlevel_info, &
           auxbas_grid = auxbas_grid, error=error)

    CALL pw_pools_create_pws(pw_pools,mgrid_rspace,&
                use_data = REALDATA3D,&
                in_space = REALSPACE, error=error)

    ! use either realspace or fft techniques to get the potential on the rs multigrids
    CALL section_vals_val_get(pw_env%interp_section,"KIND",i_val=interp_kind,error=error)
    SELECT CASE(interp_kind)
    CASE (pw_interp)
       CALL pw_pools_create_pws(pw_pools,mgrid_gspace,&
                                 use_data = COMPLEXDATA1D,&
                                 in_space = RECIPROCALSPACE, error=error)
       CALL pw_transfer(v_rspace%pw,mgrid_gspace(auxbas_grid)%pw,error=error)
       DO igrid_level=1,gridlevel_info%ngrid_levels
         IF ( igrid_level /= auxbas_grid ) THEN
              CALL pw_copy(mgrid_gspace(auxbas_grid)%pw,mgrid_gspace(igrid_level)%pw,&
                   error=error)
              CALL pw_transfer(mgrid_gspace(igrid_level)%pw,mgrid_rspace(igrid_level)%pw,&
                   error=error)
         ELSE
              IF (mgrid_gspace(auxbas_grid)%pw%pw_grid%spherical) THEN
                  CALL pw_transfer(mgrid_gspace(auxbas_grid)%pw,mgrid_rspace(auxbas_grid)%pw,&
                       error=error)
              ELSE ! fft forward + backward should be identical
                  CALL pw_copy(v_rspace%pw,mgrid_rspace(auxbas_grid)%pw,error=error)
              ENDIF
         ENDIF
         ! *** Multiply by the grid volume element ratio ***
         IF ( igrid_level /= auxbas_grid ) THEN
            scale = mgrid_rspace(igrid_level)%pw%pw_grid%dvol/&
                    mgrid_rspace(auxbas_grid)%pw%pw_grid%dvol
            mgrid_rspace(igrid_level)%pw%cr3d = &
                                      scale*mgrid_rspace(igrid_level)%pw%cr3d
         END IF
       END DO
       CALL pw_pools_give_back_pws(pw_pools,mgrid_gspace,error=error)
    CASE(spline3_pbc_interp)
       CALL pw_copy(v_rspace%pw,mgrid_rspace(1)%pw,error=error)
       DO igrid_level=1,gridlevel_info%ngrid_levels-1
          CALL pw_zero(mgrid_rspace(igrid_level+1)%pw,error=error)
          CALL pw_restrict_s3(mgrid_rspace(igrid_level)%pw,&
               mgrid_rspace(igrid_level+1)%pw,pw_pools(igrid_level+1)%pool,&
               pw_env%interp_section,error=error)
          ! *** Multiply by the grid volume element ratio
          mgrid_rspace(igrid_level+1) % pw % cr3d = &
                 mgrid_rspace(igrid_level+1) % pw % cr3d * 8._dp
       END DO
    CASE default
       CALL cp_unimplemented_error(routineN,"interpolation not supported "//&
            cp_to_string(interp_kind),error=error)
    END SELECT

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_pw_transfer(rs_v(igrid_level)%rs_grid,&
                           mgrid_rspace(igrid_level)%pw,pw2rs,error=error)
    ENDDO
    ! *** give back the pw multi-grids
    CALL pw_pools_give_back_pws(pw_pools,mgrid_rspace,error=error)

    CALL timestop(handle)

  END SUBROUTINE potential_pw2rs

END MODULE rs_pw_interface
