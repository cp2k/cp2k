!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/kg_rho_methods [1.0] *
!!
!!   NAME
!!     kg_rspw_methods
!!
!!   FUNCTION
!!     rho_methods    
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     JGH (22-Feb-03) PW grid options added
!!     gt 16-nov-03 moved initialization in this new module
!!
!!   SOURCE
!******************************************************************************

MODULE kg_rspw_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&   
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             rs_pool_create,&
                                             rs_pools_dealloc
  USE cube_utils,                      ONLY: destroy_cube_info,&
                                             init_cube_info
  USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                             gaussian_gridlevel,&
                                             init_gaussian_gridlevel
  USE green_methods,                   ONLY: pw_green_fn_rebuild
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_rspw_types,                   ONLY: kg_rspw_type,&
                                             kg_rspw_get
  USE kinds,                           ONLY: dp
  USE l_utils,                         ONLY: init_l_info,&
                                             l_info_retain,&
                                             l_info_release
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_construct,&
                                             pw_grid_release,&
                                             pw_grid_setup
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_p_type,&
                                             pw_pools_dealloc
  USE qs_util,                         ONLY: exp_radius
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE green_types,                     ONLY: greens_function_type
  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='kg_rspw_methods'
  PUBLIC :: kg_rspw_create, kg_rspw_rebuild

!-----------------------------------------------------------------------------!

CONTAINS

!!*** **********************************************************************
SUBROUTINE kg_rspw_create(kg_rspw,kg_env,error)
    TYPE(kg_rspw_type), POINTER               :: kg_rspw
    TYPE(kg_environment_type), &
      OPTIONAL, POINTER                      :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rspw_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
 
  ALLOCATE(kg_rspw, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     NULLIFY( kg_rspw%gridlevel_info, kg_rspw%pw_pools,kg_rspw%rs_pools,&
          kg_rspw%l_info, kg_rspw%cube_info)
     kg_rspw%auxbas_grid=-1
     kg_rspw%ref_count=1
     IF (PRESENT(kg_env)) CALL kg_rspw_rebuild(kg_rspw,kg_env=kg_env,error=error)
  END IF
END SUBROUTINE kg_rspw_create
!***************************************************************************
!!**** kg_rspw_methods/kg_rspw_rebuild        [1.0] *
!!
!!   NAME
!!     kg_rspw_rebuild               
!!
!!   FUNCTION
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kg_rspw: the density/potential env to be initialized
!!     - kg_env: the kg environment
!!     - error: error
!!
!!   AUTHOR
!!    gloria 
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE kg_rspw_rebuild(kg_rspw, kg_env, error )

    TYPE(kg_rspw_type), POINTER              :: kg_rspw
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rspw_rebuild ', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: safety_factor = 1.2

!   INTEGER, PARAMETER                       :: maxgridpoints = 100

    INTEGER                                  :: cmax, cmaxl, handle, &
                                                i, ikind, ipgf, iset, ishell, &
                                                igrid_level, la, lgrid_level, maxlgto, &
                                                ncommensurate, ngrid_level,&
                                                nkind, nseta, nsmax, stat
    INTEGER, DIMENSION(:), POINTER           :: npgfa, nshella
    INTEGER, DIMENSION(:, :), POINTER        :: lshella 
    LOGICAL                                  :: failure,pw_grid_information
    REAL(dp)                                 :: alpha, cutilev, drmin, maxradius,&
                                                rel_cutoff, zetp
    REAL(dp), DIMENSION(:), POINTER          :: cutoff
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(greens_function_type), POINTER      :: green
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(pw_grid_type), POINTER              :: old_pw_grid, pw_grid
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.

    NULLIFY(cutoff,cell,pw_grid,dft_control,para_env,rs_pools,old_pw_grid,&
            atomic_kind_set,green)
    logger => cp_error_get_logger(error)
    pw_grid_information=logger%print_keys%pw_grid_information

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    green=green,&
                    dft_control=dft_control,&
                    para_env=para_env)

    rel_cutoff = dft_control%qs_control%relative_cutoff
    cutoff => dft_control%qs_control%e_cutoff

    IF ( dft_control%qs_control%method == "KG_GPW" ) THEN
      ngrid_level = 1
      CPPrecondition(MAXVAL(cutoff)==cutoff(1),cp_failure_level,routineP,error,failure)
    ELSE
      ngrid_level = SIZE(cutoff)
    END IF

    CPPrecondition(ASSOCIATED(kg_rspw),cp_failure_level,routineP,error,failure)
    CPPrecondition(kg_rspw%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL pw_pools_dealloc(kg_rspw%pw_pools,error=error)
       CALL rs_pools_dealloc(kg_rspw%rs_pools,error=error)
       IF (ASSOCIATED(kg_rspw%gridlevel_info)) THEN
          CALL destroy_gaussian_gridlevel(kg_rspw%gridlevel_info,error=error)
       ELSE
          ALLOCATE(kg_rspw%gridlevel_info,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(kg_rspw%l_info)) THEN
          CALL l_info_release(kg_rspw%l_info)
       ELSE
          ALLOCATE(kg_rspw%l_info,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(kg_rspw%cube_info)) THEN
          DO igrid_level=1,SIZE(kg_rspw%cube_info)
             CALL destroy_cube_info(kg_rspw%cube_info(igrid_level))
          END DO
          DEALLOCATE(kg_rspw%cube_info,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       NULLIFY(kg_rspw%pw_pools, kg_rspw%cube_info,kg_rspw%rs_pools)
    END IF

    IF (.not.failure) THEN
      CALL init_gaussian_gridlevel(kg_rspw%gridlevel_info,&
           ngrid_level,cutoff,rel_cutoff)

!!      ALLOCATE(radius(ngrid_level),stat=stat)
!!      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!!      CALL get_radius ( kg_env, cutoff, radius )

     ! init pw_grids and pools
      ALLOCATE(pw_pools(ngrid_level),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(rs_pools(ngrid_level),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      IF (dft_control % qs_control % commensurate_mgrids) THEN
        ncommensurate=ngrid_level
      ELSE
        ncommensurate=0
      ENDIF

      DO igrid_level=1,ngrid_level 
        ALLOCATE(pw_grid,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        CALL pw_grid_construct(pw_grid)
        pw_grid%para%rs_dims ( 1 ) = para_env%num_pe
        pw_grid%para%rs_dims ( 2 ) = 1

        IF ( dft_control % qs_control % pw_grid_opt % spherical ) THEN
           cutilev = cutoff(igrid_level)
           pw_grid%grid_span = HALFSPACE
        ELSE IF ( dft_control % qs_control % pw_grid_opt % fullspace ) THEN
           cutilev = -cutoff(igrid_level)
           pw_grid%grid_span = FULLSPACE
        ELSE
           cutilev = -cutoff(igrid_level)
           pw_grid%grid_span = HALFSPACE
        END IF

        IF (igrid_level == 1) THEN
           IF (pw_grid_information) THEN
              CALL pw_grid_setup(cell,pw_grid,cutilev,&
                   pe_group=para_env%group,&
                   info= cp_logger_get_default_unit_nr(logger,&
                   local=.FALSE.), blocked=.FALSE.,ncommensurate=ncommensurate,icommensurate=igrid_level)
           ELSE
              CALL pw_grid_setup(cell,pw_grid,cutilev,&
                   pe_group=para_env%group, blocked=.FALSE.,ncommensurate=ncommensurate,icommensurate=igrid_level)
           END IF
           old_pw_grid => pw_grid
        ELSE
           IF (pw_grid_information) THEN
              CALL pw_grid_setup(cell,pw_grid,cutilev,&
                   pe_group=para_env%group,&
                   info=cp_logger_get_default_unit_nr(logger,&
                   local=.FALSE.), blocked=.FALSE., &
                   ref_grid=old_pw_grid,ncommensurate=ncommensurate,icommensurate=igrid_level)
           ELSE
              CALL pw_grid_setup(cell,pw_grid,cutilev,&
                   pe_group=para_env%group, blocked=.FALSE., &
                   ref_grid=old_pw_grid,ncommensurate=ncommensurate,icommensurate=igrid_level)
           END IF
        END IF

      ! init pw_pools
        NULLIFY(pw_pools(igrid_level)%pool)
        CALL pw_pool_create(pw_pools(igrid_level)%pool,&
             pw_grid=pw_grid,error=error)

      ! init rs_pools
        drmin = MINVAL ( pw_grid%dr )

!!        IF ( dft_control % qs_control % rs_distributed ) THEN
!!         nsmax = 2*MAX(1,CEILING((radius(igrid_level))/drmin))+1
!!      ELSE
           nsmax = -1
!!      ENDIF

        NULLIFY(rs_pools(igrid_level)%pool)
        CALL rs_pool_create(rs_pools(igrid_level)%pool,&
             el_struct=pw_grid,nsmax=nsmax,error=error)
        CALL pw_grid_release(pw_grid,error=error)
      END DO

!!     DEALLOCATE(radius,stat=stat)
!!     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      kg_rspw%pw_pools => pw_pools
      kg_rspw%rs_pools => rs_pools

      ALLOCATE (kg_rspw%cube_info(ngrid_level),STAT=stat)
      IF (stat /= 0) CALL stop_memory(routineP,"cube_info",0)

      CALL get_kg_env(kg_env=kg_env, atomic_kind_set=atomic_kind_set)
      nkind=SIZE(atomic_kind_set)

      cmax = 0
      maxradius=0.0_dp
      DO igrid_level=1,ngrid_level
        IF (.NOT. dft_control%qs_control%map_paa) maxradius=0.0_dp

        DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)

          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               orb_basis_set=orb_basis_set,&
                               alpha_core_charge=alpha)
           ! this is to be sure that the core charge is mapped ok
           ! right now, the core is mapped on the auxiliary basis,
           ! this should, at a give point be changed
           ! so that also for the core a multigrid is used
           IF (alpha > 0.0_dp)&
             maxradius=MAX(maxradius,exp_radius( 0, alpha, &
                           dft_control%qs_control%eps_rho_rspace, 10.0_dp))
           CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                npgf=npgfa, nset=nseta, zet=zeta,l=lshella,nshell=nshella)
           DO iset=1,nseta
              DO ipgf=1,npgfa(iset)
                 DO ishell=1,nshella(iset)
                   zetp = zeta(ipgf,iset)
                   la = lshella(ishell,iset)
                   IF (dft_control%qs_control%map_paa) THEN
                      lgrid_level = 1
                   ELSE
                      lgrid_level = gaussian_gridlevel(kg_rspw%gridlevel_info,zetp)
                   ENDIF
                   IF (lgrid_level .EQ. igrid_level) THEN
                     !density
                     maxradius=MAX(maxradius,exp_radius( la, zetp, &
                                  dft_control%qs_control%eps_rho_rspace, 1.0_dp))
                     !potential
                     maxradius=MAX(maxradius,exp_radius( la, zetp, &
                                  dft_control%qs_control%eps_gvg_rspace, 1.0_dp))
                   ENDIF
                 END DO
              END DO
           END DO
         END DO
!!oldcode         maxradius = maxgridpoints*MAXVAL(kg_rspw%pw_grid(igrid_level)%dr)
         maxradius = maxradius * safety_factor
      !  write(6,*) igrid_level,maxradius, CEILING(maxradius/MINVAL(pw_pools(igrid_level)%pool%pw_grid%dr(:)))
         CALL init_cube_info(kg_rspw%cube_info(igrid_level),&
               pw_pools(igrid_level)%pool%pw_grid%dr(:),&
               maxradius,cmaxl)
         cmax = MAX(cmax,cmaxl)
       END DO

       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxlgto=maxlgto)

      ! maxl+1 to allow for the forces
       CALL init_l_info(kg_rspw%l_info,maxlgto+1,cmax)

     ! init auxbas_grid
       DO i=1,ngrid_level
          IF (cutoff(i) == dft_control%qs_control%cutoff) kg_rspw%auxbas_grid=i
       END DO

!      poisson solver initialized only for the reference cutoff grid
!      Total density and potentials are allocated on the reference grid

       CALL pw_green_fn_rebuild(green,&
          kg_rspw%pw_pools(kg_rspw%auxbas_grid)%pool)

    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE kg_rspw_rebuild              

! **************************************************************************

END MODULE kg_rspw_methods

!******************************************************************************
