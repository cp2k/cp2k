!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
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
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             rs_pool_create,&
                                             rs_pools_dealloc
  USE cube_utils,                      ONLY: destroy_cube_info,&
                                             init_cube_info
  USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                             gaussian_gridlevel,&
                                             init_gaussian_gridlevel
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_rspw_types,                   ONLY: kg_rspw_type
  USE kinds,                           ONLY: dp
  USE l_utils,                         ONLY: init_l_info,&
                                             l_info_release,&
                                             l_info_retain
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_find_cutoff,&
                                             pw_grid_change,&
                                             pw_grid_create,&
                                             pw_grid_release,&
                                             pw_grid_setup
  USE pw_poisson_methods,              ONLY: pw_poisson_set
  USE pw_poisson_types,                ONLY: pw_poisson_create
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_p_type,&
                                             pw_pools_dealloc
  USE qs_util,                         ONLY: exp_radius
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_rspw_methods'
  PUBLIC :: kg_rspw_create, kg_rspw_rebuild

!-----------------------------------------------------------------------------!

CONTAINS

!!*** **********************************************************************
SUBROUTINE kg_rspw_create(kg_rspw,kg_env,error)
    TYPE(kg_rspw_type), POINTER              :: kg_rspw
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
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
          kg_rspw%l_info, kg_rspw%cube_info, kg_rspw%poisson_env)
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

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rspw_rebuild', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: safety_factor = 1.2

    INTEGER :: cmax, cmaxl, handle, i, igrid_level, ikind, ipgf, iset, &
      ishell, la, lgrid_level, maxlgto, ncommensurate, ngrid_level, nkind, &
      nseta, nsmax, stat
    INTEGER, DIMENSION(:), POINTER           :: npgfa, nshella
    INTEGER, DIMENSION(:, :), POINTER        :: lshella
    LOGICAL                                  :: failure, use_ref_cell
    REAL(dp)                                 :: alpha, cutilev, drmin, &
                                                maxradius, rel_cutoff, zetp, &
                                                my_cut
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cutoff, my_cutoff
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell, cell_ref, my_cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(pw_grid_type), POINTER              :: old_pw_grid, pw_grid
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(section_vals_type), POINTER         :: poisson_section, print_section

!   INTEGER, PARAMETER                       :: maxgridpoints = 100
!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)
!
!
! Part one, deallocate old data if needed
!
!

    failure=.FALSE.
    NULLIFY(cutoff,cell,pw_grid,dft_control,para_env,rs_pools,old_pw_grid,&
            atomic_kind_set)

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    cell_ref=cell_ref,&
                    dft_control=dft_control,&
                    input=poisson_section,&
                    para_env=para_env)

    CPPrecondition(ASSOCIATED(kg_rspw),cp_failure_level,routineP,error,failure)
    CPPrecondition(kg_rspw%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL pw_pools_dealloc(kg_rspw%pw_pools,error=error)
       CALL rs_pools_dealloc(kg_rspw%rs_pools,error=error)
       IF (ASSOCIATED(kg_rspw%gridlevel_info)) THEN
          CALL destroy_gaussian_gridlevel(kg_rspw%gridlevel_info, error=error)
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


!
!
! Part two, setup the pw_grids
!
!

  IF (.NOT.failure) THEN
     CALL get_kg_env ( kg_env, use_ref_cell = use_ref_cell )
     IF (use_ref_cell) THEN
       my_cell => cell_ref
     ELSE
       my_cell => cell
     END IF

     rel_cutoff = dft_control%qs_control%relative_cutoff
     cutoff => dft_control%qs_control%e_cutoff
     IF ( dft_control%qs_control%method == "KG_GPW" ) THEN
       ngrid_level = 1
       CPPrecondition(MAXVAL(cutoff)==cutoff(1),cp_failure_level,routineP,error,failure)
     ELSE
       ngrid_level = SIZE(cutoff)
     END IF
     print_section=>section_vals_get_subs_vals(kg_env%input, &
                       "PRINT%GRID_INFORMATION",error=error)
     CALL init_gaussian_gridlevel(kg_rspw%gridlevel_info,&
           ngrid_level,cutoff,rel_cutoff,print_section=print_section)
     ! init pw_grids and pools
     ALLOCATE(pw_pools(ngrid_level),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     IF (dft_control % qs_control % commensurate_mgrids) THEN
       ncommensurate=ngrid_level
     ELSE
       ncommensurate=0
     ENDIF

     DO igrid_level=1,ngrid_level 
        CALL pw_grid_create(pw_grid,para_env%group)

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
           IF (ASSOCIATED(old_pw_grid)) THEN
             CALL pw_grid_setup(my_cell,pw_grid, print_section, cutilev,&
                  ref_grid=old_pw_grid, &
                  blocked=.FALSE., ncommensurate=ncommensurate, &
                  icommensurate=igrid_level,rs_dims=(/para_env%num_pe,1/))
              old_pw_grid => pw_grid
           ELSE
              CALL pw_grid_setup(my_cell,pw_grid,print_section, cutilev,&
                   blocked=.FALSE.,ncommensurate=ncommensurate,&
                   icommensurate=igrid_level, rs_dims=(/para_env%num_pe,1/))
              old_pw_grid => pw_grid
           END IF
        ELSE
           CALL pw_grid_setup(my_cell,pw_grid, print_section, cutilev,&
                blocked=.FALSE.,ref_grid=old_pw_grid,ncommensurate=ncommensurate,&
                icommensurate=igrid_level, rs_dims=(/para_env%num_pe,1/))
        ENDIF

      ! init pw_pools
         NULLIFY(pw_pools(igrid_level)%pool)
         CALL pw_pool_create(pw_pools(igrid_level)%pool,&
              pw_grid=pw_grid,error=error)
         CALL pw_grid_release(pw_grid,error=error)

       END DO

       kg_rspw%pw_pools => pw_pools

     ! init auxbas_grid
       DO i=1,ngrid_level
         IF (cutoff(i) == dft_control%qs_control%cutoff) kg_rspw%auxbas_grid=i
       END DO
!      poisson solver initialized only for the reference cutoff grid
!      Total density and potentials are allocated on the reference grid

        IF (.NOT.ASSOCIATED(kg_rspw%poisson_env)) THEN
          CALL pw_poisson_create(kg_rspw%poisson_env,error=error)
        END IF
        poisson_section => section_vals_get_subs_vals(poisson_section,"DFT%POISSON",&
                          error=error)
        CALL pw_poisson_set(kg_rspw%poisson_env, cell=my_cell, &
            parameters=poisson_section, pw_pools=kg_rspw%pw_pools,&
            use_level=kg_rspw%auxbas_grid, error=error )
!
! If reference cell is present, then use pw_grid_change to keep bounds constant...
!
!
        IF ( use_ref_cell ) THEN
           ALLOCATE( my_cutoff ( SIZE ( cutoff ) ),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO igrid_level = 1, SIZE ( pw_pools )
             CALL pw_grid_change ( cell, pw_pools ( igrid_level ) % pool % pw_grid )
             CALL pw_find_cutoff ( pw_pools ( igrid_level ) % pool % pw_grid % bounds,  &
                                   cell%h_inv, my_cut )
             my_cutoff (igrid_level) = my_cut
           ENDDO
           CALL destroy_gaussian_gridlevel (kg_rspw%gridlevel_info, error=error )
           print_section=>section_vals_get_subs_vals(kg_env%input, &
                               "PRINT%GRID_INFORMATION",error=error)
           CALL init_gaussian_gridlevel(kg_rspw%gridlevel_info,&
                ngrid_levels=ngrid_level,cutoff=my_cutoff,rel_cutoff=rel_cutoff, &
                print_section=print_section)
           CALL pw_poisson_set(kg_rspw%poisson_env,cell=cell,pw_pools=kg_rspw%pw_pools,&
           parameters=poisson_section, use_level=kg_rspw%auxbas_grid,error=error)
           DEALLOCATE( my_cutoff ,stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ENDIF
!
!
!
!    determine the maximum radii for mapped gaussians, needed to
!    set up distributed rs grids
!
!


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
        ! safety first, not very efficient, nor general
        ! one could possibly decide that this *is* the maximum allowed radius
         maxradius = maxradius * safety_factor
      END DO
!
!
!    set up the rs_grids, 
!
!
      ALLOCATE(rs_pools(ngrid_level),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO igrid_level=1,ngrid_level 
        pw_grid => pw_pools(igrid_level)%pool%pw_grid
! init rs_pools
        drmin = MINVAL ( pw_grid%dr )
!!        IF ( dft_control % qs_control % rs_distributed ) THEN
!!           nsmax = 2*MAX(1,CEILING((radius(igrid_level))/drmin))+1
!!      ELSE
            nsmax = -1
!!        ENDIF
        NULLIFY(rs_pools(igrid_level)%pool)
        CALL rs_pool_create(rs_pools(igrid_level)%pool,&
            el_struct=pw_grid,nsmax=nsmax,error=error)
      ENDDO
      kg_rspw%rs_pools => rs_pools

!
!
! set up the cube info and the cmaxl array
!
!
      DO igrid_level=1,ngrid_level 
         CALL init_cube_info(kg_rspw%cube_info(igrid_level),&
               pw_pools(igrid_level)%pool%pw_grid%dr(:),&
               pw_pools(igrid_level)%pool%pw_grid%dh(:,:),&
               pw_pools(igrid_level)%pool%pw_grid%orthorhombic,&
               maxradius,cmaxl)
         cmax = MAX(cmax,cmaxl)
      END DO
!
!
!  l_info set up here, requires cmaxl
!
!

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxlgto=maxlgto)

! maxl+1 to allow for the forces
      CALL init_l_info(kg_rspw%l_info,maxlgto+1,cmax)



   END IF

   CALL timestop(0.0_dp,handle)

  END SUBROUTINE kg_rspw_rebuild              

! **************************************************************************

END MODULE kg_rspw_methods

!******************************************************************************
