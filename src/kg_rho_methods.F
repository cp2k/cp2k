!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/kg_rho_methods [1.0] *
!!
!!   NAME
!!     kg_rho_methods
!!
!!   FUNCTION
!!     rho_methods    
!!
!!   AUTHOR
!!     gloria,30.09.2002, previously in kg_force 
!!
!!   MODIFICATION HISTORY
!!     JGH (22-Feb-03) PW grid options added
!!     gt 16-nov-03 moved initialization in this new module
!!
!!   SOURCE
!******************************************************************************

MODULE kg_rho_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
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
  USE cube_utils,                      ONLY: destroy_cube_info,&
                                             init_cube_info
  USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                             init_gaussian_gridlevel
  USE greens_fns,                      ONLY: pw_green_fn_rebuild
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_rho_types,                    ONLY: allocate_kg_densities,&
                                             kg_rho_type
  USE kinds,                           ONLY: dp
  USE l_utils,                         ONLY: destroy_l_info,&
                                             init_l_info
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_construct,&
                                             pw_grid_setup
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_type
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='kg_rho_methods'
  PUBLIC :: kg_rho_create, kg_rho_setup

!-----------------------------------------------------------------------------!

CONTAINS

!!*** **********************************************************************
SUBROUTINE kg_rho_create(kg_rho,kg_env,error)
    TYPE(kg_rho_type), POINTER               :: kg_rho
    TYPE(kg_environment_type), INTENT(in), &
      OPTIONAL, TARGET                       :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rho_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
 
  ALLOCATE(kg_rho, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     NULLIFY( kg_rho%gridlevel_info, kg_rho%pw_pool,&
          kg_rho%l_info, kg_rho%cube_info)
     kg_rho%auxbas_grid=-1
     kg_rho%ref_count=1
     IF (PRESENT(kg_env)) CALL kg_rho_setup(kg_rho,kg_env=kg_env,error=error)
  END IF
END SUBROUTINE kg_rho_create
!***************************************************************************
!!**** kg_rho_methods/kg_rho_setup        [1.0] *
!!
!!   NAME
!!     kg_rho_setup               
!!
!!   FUNCTION
!!     initializes the environment needed for calculation of kg forces
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kg_rho: the density/potential env to be initialized
!!     - kg_env: the kg environment
!!     - error: error
!!
!!   AUTHOR
!!    gloria 
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE kg_rho_setup(kg_rho, kg_env, error )



    TYPE(kg_rho_type), POINTER               :: kg_rho
    TYPE(kg_environment_type), INTENT(IN), &
      TARGET                                 :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rho_setup ', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: maxgridpoints = 100

    INTEGER                                  :: cmax, cmaxl, handle, &
                                                igrid_level, istat, maxlgto, &
                                                ngrid_level, stat
    LOGICAL                                  :: failure
    REAL(dp)                                 :: maxradius, rel_cutoff
    REAL(dp), DIMENSION(:), POINTER          :: cutoff
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dft_control_type), POINTER          :: dft_control

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.

    NULLIFY(cutoff,dft_control,atomic_kind_set)
    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control)

    kg_rho%kgpol = dft_control % qs_control % polarization
    rel_cutoff = dft_control%qs_control%relative_cutoff
    cutoff => dft_control%qs_control%e_cutoff
    IF ( dft_control%qs_control%method == "KG_GPW" ) THEN
      ngrid_level = 1
      CPPrecondition(MAXVAL(cutoff)==cutoff(1),cp_failure_level,routineP,error,failure)
    ELSE
      ngrid_level = SIZE(cutoff)
    END IF

    CPPrecondition(ASSOCIATED(kg_rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(kg_rho%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (ASSOCIATED(kg_rho%gridlevel_info)) THEN
          CALL destroy_gaussian_gridlevel(kg_rho%gridlevel_info,error=error)
       ELSE
          ALLOCATE(kg_rho%gridlevel_info,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(kg_rho%l_info)) THEN
          CALL destroy_l_info(kg_rho%l_info)
       ELSE
          ALLOCATE(kg_rho%l_info,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(kg_rho%cube_info)) THEN
          DO igrid_level=1,SIZE(kg_rho%cube_info)
             CALL destroy_cube_info(kg_rho%cube_info(igrid_level))
          END DO
          DEALLOCATE(kg_rho%cube_info,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       NULLIFY(kg_rho%pw_pool, kg_rho%cube_info)
    END IF

    IF (.not.failure) THEN
      CALL init_gaussian_gridlevel(kg_rho%gridlevel_info,&
           ngrid_level,cutoff,rel_cutoff)

      CALL init_grid (kg_rho,kg_env,&
           ngrid_level=ngrid_level,error=error)

      ALLOCATE (kg_rho%cube_info(ngrid_level),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"cube_info",0)

      cmax = 0
      DO igrid_level=1,ngrid_level
         maxradius = maxgridpoints*MAXVAL(kg_rho%pw_grid(igrid_level)%dr)
         CALL init_cube_info(kg_rho%cube_info(igrid_level),&
              kg_rho%pw_grid(igrid_level)%dr(:),&
              maxradius,cmaxl)
          cmax = MAX(cmax,cmaxl)
      END DO

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxlgto=maxlgto)

      ! maxl+1 to allow for the forces
      CALL init_l_info(kg_rho%l_info,maxlgto+1,cmax)

      kg_rho%total_rho_core_rspace=0.0_dp
      kg_rho%total_rho_rspace=0.0_dp
      kg_rho%total_rho_gspace=0.0_dp
    
      CALL allocate_kg_densities(kg_rho)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE kg_rho_setup              

! **************************************************************************

  SUBROUTINE init_grid ( kg_rho, kg_env, ngrid_level,error )

    TYPE(kg_rho_type), INTENT(INOUT)         :: kg_rho
    TYPE(kg_environment_type), INTENT(IN)    :: kg_env
    INTEGER, INTENT(IN)                      :: ngrid_level
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_grid' , &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: auxbas, handle, i, &
                                                igrid_level, istat
    LOGICAL                                  :: failure, pw_grid_information
    REAL(dp)                                 :: cutilev
    REAL(dp), DIMENSION(:), POINTER          :: cutoff
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_grid_type), DIMENSION(:), &
      POINTER                                :: pw_grid
    TYPE(pw_grid_type), POINTER              :: pw_grid_loc
    TYPE(pw_pool_type), POINTER              :: pw_pool

!   ---------------------------------------------------------------------------

    CALL timeset("init_grid","I","",handle)
    failure=.FALSE.
    NULLIFY(cutoff,cell,pw_grid,dft_control,para_env)              

    logger => cp_error_get_logger(error)
    pw_grid_information=logger%print_keys%pw_grid_information

    IF (.NOT.failure) THEN
      CALL get_kg_env(kg_env=kg_env,&
          dft_control=dft_control, cell=cell, para_env=para_env)    

      cutoff => dft_control % qs_control % e_cutoff

    ! the density grid coincides with the multi-grid having the highest cutoff
    ! auxbas_grid is the index of the auxiliary (REFERENCE) grid

      DO i=1,ngrid_level
        IF ( cutoff(i) == dft_control%qs_control%cutoff ) kg_rho % auxbas_grid = i
      END DO

!   *** Allocate the multi-grid pointer pw_grid  ***

      ALLOCATE (kg_rho%pw_grid(ngrid_level),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"pw_grid",0)

!   *** Construct the various grids ***

      pw_grid => kg_rho%pw_grid
      DO igrid_level=1,ngrid_level
         CALL pw_grid_construct(pw_grid(igrid_level))
         pw_grid(igrid_level)%grid_span = HALFSPACE
         pw_grid(igrid_level)%para%rs_dims ( 1 ) = para_env % num_pe
         pw_grid(igrid_level)%para%rs_dims ( 2 ) = 1

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
             CALL pw_grid_setup(cell,pw_grid(igrid_level),cutilev,&
                                pe_group=para_env%group,&
                                info=cp_logger_get_default_unit_nr(logger,&
                                local=.FALSE.))
           ELSE
             CALL pw_grid_setup(cell,pw_grid(igrid_level),cutilev,&
                                pe_group=para_env%group)
           END IF
         ELSE
           IF (pw_grid_information) THEN
             CALL pw_grid_setup(cell,pw_grid(igrid_level),cutilev,&
                                pe_group=para_env%group,&
                                info=cp_logger_get_default_unit_nr(logger,&
                                 local=.FALSE.),&
                                ref_grid=pw_grid(igrid_level-1))
           ELSE
             CALL pw_grid_setup(cell,pw_grid(igrid_level),cutilev,&
                                pe_group=para_env%group,&
                                ref_grid=pw_grid(igrid_level-1))
           END IF
         END IF

      END DO

      auxbas = kg_rho % auxbas_grid

!   poisson solver initialized only for the reference cutoff grid 
!   Total density and potentials are allocated on the reference grid

      pw_grid_loc=>pw_grid(auxbas)
      !   pw pools initialized
      NULLIFY (pw_pool)
      CALL pw_pool_create(pw_pool, pw_grid=pw_grid_loc)     
      CALL pw_green_fn_rebuild(cell%green,pw_pool)

      kg_rho%pw_pool => pw_pool
    END IF 
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE init_grid

! *****************************************************************************

! SUBROUTINE get_grids ( pw_grids, pw_pool, gridlevel_info, l_info, &
!                         cube_info, ref_grid )

!  IMPLICIT NONE

! TYPE ( pw_grid_type ), DIMENSION (:), POINTER, OPTIONAL :: pw_grids
! TYPE ( pw_pool_type ), POINTER, OPTIONAL :: pw_pool
! TYPE ( gridlevel_info_type ), POINTER, OPTIONAL :: gridlevel_info
! TYPE ( l_info_type ), POINTER, OPTIONAL         :: l_info
! TYPE ( cube_info_type ), DIMENSION(:), POINTER, OPTIONAL :: cube_info
! INTEGER, OPTIONAL, INTENT(OUT) :: ref_grid


! IF ( PRESENT ( pw_grids ) ) pw_grids => kg_rho % pw_grid
! IF ( PRESENT ( pw_pool ) ) pw_pool => kg_rho % pw_pool
! IF ( PRESENT ( gridlevel_info ) ) gridlevel_info => kg_rho % gridlevel_info
! IF ( PRESENT ( l_info ) ) l_info => kg_rho % l_info
! IF ( PRESENT ( cube_info ) ) cube_info => kg_rho % cube_info
! IF ( PRESENT ( ref_grid ) ) ref_grid =  kg_rho % auxbas_grid   

!END SUBROUTINE get_grids
!******************************************************************************

END MODULE kg_rho_methods

!******************************************************************************
