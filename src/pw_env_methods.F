!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002 - 2004 CP2K developers group                           !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pw_env_methods [1.0] *
!!
!!   NAME
!!     pw_env_methods
!!
!!   FUNCTION
!!     methods of pw_env that have dependence on qs_env
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!     JGH (22-Feb-03) PW grid options added
!!     04.2003 added rs grid pools [fawzi]
!!     02.2004 added commensurate grids
!!
!!   SOURCE
!****************************************************************************
MODULE pw_env_methods
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
  USE kinds,                           ONLY: dp
  USE l_utils,                         ONLY: init_l_info,&
                                             l_info_retain,&
                                             l_info_release
  USE orbital_pointers,                ONLY: init_orbital_pointers
  USE pw_env_types,                    ONLY: pw_env_type
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_change,&
                                             pw_grid_construct,&
                                             pw_grid_release,&
                                             pw_grid_setup
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_p_type,&
                                             pw_pools_dealloc
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_interactions,                 ONLY: exp_radius_very_extended
  USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                             rho0_mpole_type
  USE qs_util,                         ONLY: exp_radius
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE green_types,                     ONLY: greens_function_type

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pw_env_methods'

  PUBLIC :: pw_env_create, pw_env_rebuild
!***
!****************************************************************************
CONTAINS

!!****f* pw_env_methods/pw_env_create [1.0] *
!!
!!   NAME
!!     pw_env_create
!!
!!   FUNCTION
!!     creates a pw_env, if qs_env is given calls pw_env_rebuild
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pw_env: the pw_env that gets created
!!     - qs_env: if given it is used to initialize the pw_env
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_env_create(pw_env,qs_env,error)
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_environment_type), OPTIONAL, POINTER :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  ALLOCATE(pw_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     NULLIFY(pw_env%pw_pools, pw_env%gridlevel_info,&
          pw_env%l_info, pw_env%l_info_rho0, pw_env%cube_info, pw_env%rs_pools)
     pw_env%auxbas_grid=-1
     pw_env%ref_count=1
     IF (PRESENT(qs_env)) CALL pw_env_rebuild(pw_env,qs_env=qs_env,error=error)
  END IF
END SUBROUTINE pw_env_create

!***************************************************************************

!!****f* pw_env_methods/pw_env_rebuild [1.0] *
!!
!!   NAME
!!     pw_env_rebuild
!!
!!   FUNCTION
!!     rebuilds the pw_env data (necessary if cell or cutoffs change)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pw_env: the environment to rebuild
!!     - qs_env: the qs_env where to get the cell, cutoffs,...
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_env_rebuild(pw_env, qs_env, error)
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_rebuild', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: safety_factor = 1.2 

    INTEGER :: cmax, cmaxl, handle, i, igrid_level, ikind, ipgf, iset, &
      ishell, istat, jkind, jpgf, jset, jshell, la, lb, lgrid_level, igrid_zet0_s,&
      lmax_rho0, maxlgto, ncommensurate, ngrid_level, nkind, nseta, nsetb, nsmax, stat
    INTEGER, DIMENSION(:), POINTER           :: npgfa, npgfb, nshella, nshellb
    INTEGER, DIMENSION(:, :), POINTER        :: lshella, lshellb
    LOGICAL                                  :: failure, pw_grid_information, &
                                                use_ref_cell
    REAL(KIND=dp)                            :: alpha, cutilev, drmin, &
                                                maxradius, rel_cutoff, &
                                                max_rpgf0_s, zet0_s, zetp
    TYPE(greens_function_type), POINTER      :: green
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: radius
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cutoff
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(pw_grid_type), POINTER              :: old_pw_grid, pw_grid
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole, nddo_mpole

! should be be avoided (see later)

  CALL timeset(routineN,"I","",handle)
  failure=.FALSE.
  NULLIFY(cutoff,cell,pw_grid,old_pw_grid,dft_control,atomic_kind_set,&
          pw_pools,rho0_mpole,rs_pools, para_env,cell_ref, green)
  logger => cp_error_get_logger(error)
  pw_grid_information=logger%print_keys%pw_grid_information
 
  CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL pw_pools_dealloc(pw_env%pw_pools,error=error)
     CALL rs_pools_dealloc(pw_env%rs_pools,error=error)
     IF (ASSOCIATED(pw_env%gridlevel_info)) THEN
        CALL destroy_gaussian_gridlevel(pw_env%gridlevel_info,error=error)
     ELSE
        ALLOCATE(pw_env%gridlevel_info,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF (ASSOCIATED(pw_env%l_info)) THEN
        CALL l_info_release(pw_env%l_info,error)
     ELSE
        ALLOCATE(pw_env%l_info,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF

     IF(qs_env%dft_control%qs_control%gapw) THEN
        IF (ASSOCIATED(pw_env%l_info_rho0)) THEN
           CALL l_info_release(pw_env%l_info_rho0,error)
        ELSE
           ALLOCATE(pw_env%l_info_rho0,stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
        CALL get_qs_env(qs_env=qs_env,rho0_mpole=rho0_mpole)
        CPPostcondition(ASSOCIATED(rho0_mpole),cp_failure_level,routineP,error,failure)
        CALL get_rho0_mpole(rho0_mpole=rho0_mpole,&
                         zet0_s=zet0_s,max_rpgf0_s=max_rpgf0_s)
     ELSEIF(qs_env%dft_control%qs_control%se_control%ewald) THEN
        IF (ASSOCIATED(pw_env%l_info_rho0)) THEN
           CALL l_info_release(pw_env%l_info_rho0,error)
        ELSE
           ALLOCATE(pw_env%l_info_rho0,stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
        CALL get_qs_env(qs_env=qs_env,nddo_mpole=nddo_mpole)
        CPPostcondition(ASSOCIATED(nddo_mpole),cp_failure_level,routineP,error,failure)
        CALL get_rho0_mpole(rho0_mpole=nddo_mpole,&
                         zet0_s=zet0_s,max_rpgf0_s=max_rpgf0_s)
     ELSE
        NULLIFY(pw_env%l_info_rho0)
     END IF

     IF (ASSOCIATED(pw_env%cube_info)) THEN
        DO igrid_level=1,SIZE(pw_env%cube_info)
           CALL destroy_cube_info(pw_env%cube_info(igrid_level))
        END DO
        DEALLOCATE(pw_env%cube_info,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
     NULLIFY(pw_env%pw_pools, pw_env%cube_info)
  END IF
  IF (.NOT.failure) THEN
     CALL get_qs_env(qs_env=qs_env,&
          dft_control=dft_control,atomic_kind_set=atomic_kind_set,&
          cell_ref=cell_ref, green=green, cell=cell,para_env=para_env)
     rel_cutoff = dft_control%qs_control%relative_cutoff
     cutoff => dft_control%qs_control%e_cutoff
     ngrid_level = SIZE(cutoff)
     
     ! init gridlevel_info XXXXXXXXX setup mapping to the effective cutoff ?
     !                     XXXXXXXXX the cutoff array here is more a 'wish-list'
     !                     XXXXXXXXX same holds for radius
     CALL init_gaussian_gridlevel(pw_env%gridlevel_info,&
          ngrid_levels=ngrid_level,cutoff=cutoff,rel_cutoff=rel_cutoff)
     ALLOCATE(radius(ngrid_level),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     CALL get_subgrid_size ( qs_env, cutoff, radius ) 

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
     DO igrid_level=1,ngrid_level ! keep the pw_grids in a separate array?
! Start from scratch
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
              CALL pw_grid_setup(cell_ref,pw_grid,cutilev,&
                   pe_group=para_env%group,&
                   info= cp_logger_get_default_unit_nr(logger,&
                   local=.FALSE.), blocked=.FALSE.,ncommensurate=ncommensurate,icommensurate=igrid_level)
           ELSE
              CALL pw_grid_setup(cell_ref,pw_grid,cutilev,&
                   pe_group=para_env%group, blocked=.FALSE.,ncommensurate=ncommensurate,icommensurate=igrid_level)
           END IF
           old_pw_grid => pw_grid
        ELSE
           IF (pw_grid_information) THEN
              CALL pw_grid_setup(cell_ref,pw_grid,cutilev,&
                   pe_group=para_env%group,&
                   info=cp_logger_get_default_unit_nr(logger,&
                   local=.FALSE.), blocked=.FALSE., &
                   ref_grid=old_pw_grid,ncommensurate=ncommensurate,icommensurate=igrid_level)
           ELSE
              CALL pw_grid_setup(cell_ref,pw_grid,cutilev,&
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

        IF ( dft_control % qs_control % rs_distributed ) THEN
           nsmax = 2*MAX(1,CEILING((radius(igrid_level))/drmin))+1
        ELSE
           nsmax = -1
        ENDIF

        NULLIFY(rs_pools(igrid_level)%pool)
        CALL rs_pool_create(rs_pools(igrid_level)%pool,&
             el_struct=pw_grid,nsmax=nsmax,error=error)
        CALL pw_grid_release(pw_grid,error=error)
     END DO

     ! Find the grid level suitable for rho0_soft
     IF(qs_env%dft_control%qs_control%gapw) THEN
        igrid_zet0_s = gaussian_gridlevel(pw_env%gridlevel_info,2.0_dp*zet0_s)
        rho0_mpole%igrid_zet0_s = igrid_zet0_s
     ELSEIF(qs_env%dft_control%qs_control%se_control%ewald) THEN
        igrid_zet0_s = gaussian_gridlevel(pw_env%gridlevel_info,2.0_dp*zet0_s)
        nddo_mpole%igrid_zet0_s = igrid_zet0_s
     END IF

     DEALLOCATE(radius,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     pw_env%pw_pools => pw_pools
     pw_env%rs_pools => rs_pools

     ! init cube info
     ALLOCATE (pw_env%cube_info(ngrid_level),STAT=istat)
     IF (istat /= 0) CALL stop_memory(routineP,"cube_info",0)

     CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
     nkind=SIZE(atomic_kind_set)

     cmax = 0
     ! XXXXXXXXXXXXXX this duplicates the functionality of get_subgrid_size  XXXXXXXXXXXXXXXXXX
     ! try to predict the maximum radius of the gaussians to be mapped on the grid
     ! up to now, it is not yet very good 
     maxradius=0.0_dp
     DO igrid_level=1,ngrid_level
        IF (.NOT. dft_control%qs_control%map_paa) maxradius=0.0_dp

        ! Take into account the radius of the soft compensation charge rho0_soft1
        IF(qs_env%dft_control%qs_control%gapw .OR. &
           qs_env%dft_control%qs_control%se_control%ewald) THEN
           IF(igrid_zet0_s == igrid_level ) maxradius=MAX(maxradius,max_rpgf0_s)
        END IF

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

           DO jkind=1,nkind
              atomic_kind => atomic_kind_set(jkind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)
              CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                     npgf=npgfb, nset=nsetb, zet=zetb, l=lshellb, nshell=nshellb)
              DO iset=1,nseta
              DO ipgf=1,npgfa(iset)
              DO ishell=1,nshella(iset)
                 ! for calculate wavefunction
                 zetp = zeta(ipgf,iset)
                 la = lshella(ishell,iset)
                 IF (dft_control%qs_control%map_paa) THEN
                    lgrid_level = 1
                 ELSE
                    lgrid_level = gaussian_gridlevel(pw_env%gridlevel_info,zetp)
                 ENDIF
                 IF (lgrid_level .EQ. igrid_level) THEN
                    maxradius=MAX(maxradius,exp_radius( la, zetp, &
                                  dft_control%qs_control%eps_rho_rspace, 1.0_dp))
                 ENDIF
                 DO jset=1,nsetb
                 DO jpgf=1,npgfb(jset)
                 DO jshell=1,nshellb(jset)
                    zetp = zeta(ipgf,iset) + zetb(jpgf,jset)
                    lb = lshellb(jshell,jset)+la
                    lgrid_level = gaussian_gridlevel(pw_env%gridlevel_info,zetp)
                    IF (lgrid_level .EQ. igrid_level) THEN
                       ! density
                       maxradius=MAX(maxradius,exp_radius( lb, zetp, &
                                     dft_control%qs_control%eps_rho_rspace, 1.0_dp))
                       ! potential
                       maxradius=MAX(maxradius,exp_radius( lb, zetp, &
                                     dft_control%qs_control%eps_gvg_rspace, 1.0_dp))
                    ENDIF
                 ENDDO
                 ENDDO
                 ENDDO
              ENDDO
              ENDDO
              ENDDO
         ENDDO
        ENDDO
        ! safety first, not very efficient, nor general
        ! one could possibly decide that this *is* the maximum allowed radius
        maxradius = maxradius * safety_factor 
        pw_env%cube_info%max_rad_ga = maxradius
        ! write(6,*) igrid_level,maxradius, CEILING(maxradius/MINVAL(pw_pools(igrid_level)%pool%pw_grid%dr(:)))
        CALL init_cube_info(pw_env%cube_info(igrid_level),&
             pw_pools(igrid_level)%pool%pw_grid%dr(:),&
             maxradius,cmaxl)
        cmax = MAX(cmax,cmaxl)
     END DO

     ! init l_info
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
          maxlgto=maxlgto)
     CALL init_l_info(pw_env%l_info,maxlgto+2,cmax)

     IF( dft_control%qs_control%gapw .OR. &
         dft_control%qs_control%se_control%ewald ) THEN
        CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
             lmax_rho0 = lmax_rho0)
        IF (lmax_rho0+1 <= maxlgto+2) THEN
           DEALLOCATE(pw_env%l_info_rho0,stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           pw_env%l_info_rho0 => pw_env%l_info
           CALL l_info_retain(pw_env%l_info,error)
        ELSE
           CALL init_l_info(pw_env%l_info_rho0,lmax_rho0+1,cmax,lrho0=.TRUE.)
        END IF
     END IF
 
     ! init auxbas_grid
     DO i=1,ngrid_level
        IF (cutoff(i) == dft_control%qs_control%cutoff) pw_env%auxbas_grid=i
     END DO

     ! complete init of green's function
     CALL pw_green_fn_rebuild(green,&
          pw_env%pw_pools(pw_env%auxbas_grid)%pool, cell=cell)
  END IF

! If reference cell is present, then use pw_grid_change to keep bounds constant...
! MJM is looking at a way to phase this out, because it doesn't seem like everything gets initialized properly
  CALL get_qs_env ( qs_env, use_ref_cell = use_ref_cell )
  IF ( use_ref_cell ) THEN
      DO igrid_level = 1, SIZE ( pw_pools )
        CALL pw_grid_change ( cell, pw_pools ( igrid_level ) % pool % pw_grid )
      ENDDO
      CALL pw_green_fn_rebuild(green,&
          pw_env%pw_pools(pw_env%auxbas_grid)%pool, cell=cell)
  ENDIF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE pw_env_rebuild

!***************************************************************************

SUBROUTINE get_subgrid_size ( qs_env, cutoff, radius )
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: cutoff
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: radius

    INTEGER                                  :: i, igrid_level, ikind, j, &
                                                maxlgto, numgrids, stat
    REAL(KIND=dp)                            :: alpha, charge, &
                                                eps_rho_rspace, lcut, &
                                                maxcharge, mincoreval, &
                                                minexpval, rel_cutoff
    REAL(KIND=dp), DIMENSION(3)              :: ra, rb, rp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p, zet
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

  radius = 1.E10_dp

  numgrids = SIZE ( cutoff )
  ra = 0._dp
  rb = 0._dp
  rp = 0._dp

  NULLIFY(atomic_kind_set,atomic_kind,dft_control,orb_basis_set,p,zet)
  ALLOCATE (p(1,1),STAT=stat)
  IF ( stat /= 0 ) CALL stop_memory ("get_subgrid_size","p",1)
  p = 1._dp

  CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,atomic_kind_set=atomic_kind_set)

  rel_cutoff = dft_control%qs_control%relative_cutoff
  eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

  CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxlgto=maxlgto)
  ! derivatives need l+1
  maxlgto = maxlgto + 1
  ! be sure coset is set as it should for exp_radius_very_extended
  CALL init_orbital_pointers(maxlgto)

  minexpval = 1.E10_dp
  mincoreval = 1.E10_dp
  maxcharge = 0._dp
  DO ikind=1,SIZE(atomic_kind_set)
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set,&
                           alpha_core_charge=alpha,&
                           ccore_charge=charge)
      mincoreval = MIN ( alpha, mincoreval )
      maxcharge = MAX ( charge, maxcharge )

      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             zet=zet)
      DO i=1,SIZE(zet,2)
        DO j=1,SIZE(zet,1)
          IF ( zet(j,i) > 0._dp) THEN
            minexpval = MIN ( zet(j,i), minexpval )
          END IF
        END DO
      END DO
  END DO
  minexpval = 2._dp * minexpval

  DO igrid_level = 1, numgrids

    IF ( igrid_level == 1 ) THEN
      ! this should be the grid level with the highest cutoff
      IF ( MAXVAL ( ABS ( cutoff ) ) /= ABS ( cutoff (1) ) ) &
          CALL stop_program("get_subgrid_size","cutoff not correctly ordered")
      IF ( numgrids == 1 ) THEN
        lcut = minexpval
      ELSE
        lcut = ABS ( cutoff ( igrid_level+1 ) ) / rel_cutoff
        lcut = MAX ( lcut, minexpval )
      END IF

    ELSE IF ( igrid_level == numgrids ) THEN
      ! check that grid-1 has a higher cutoff
      IF ( ABS(cutoff(igrid_level)) > ABS(cutoff(igrid_level-1)) ) &
          CALL stop_program("get_subgrid_size","cutoff not correctly ordered")
      lcut = minexpval

    ELSE
      ! check that grid-1 has a higher cutoff
      IF ( ABS(cutoff(igrid_level)) > ABS(cutoff(igrid_level-1)) ) &
          CALL stop_program("get_subgrid_size","cutoff not correctly ordered")
      ! check that grid+1 has a lower cutoff
      IF ( ABS(cutoff(igrid_level)) < ABS(cutoff(igrid_level+1)) ) &
          CALL stop_program("get_subgrid_size","cutoff not correctly ordered")
      lcut = ABS ( cutoff ( igrid_level+1 ) ) / rel_cutoff
      lcut = MAX ( lcut, minexpval )

    END IF
    radius(igrid_level)=exp_radius_very_extended(0,maxlgto,0,maxlgto,ra=ra,rb=rb,rp=rp,&
                             zetp=lcut,eps=eps_rho_rspace,prefactor=1._dp,cutoff=0._dp)
    IF ( igrid_level == 1 ) THEN
      p=maxcharge
      IF (mincoreval > 0.0_dp)&
        radius(1)=MAX(radius(1), exp_radius_very_extended(0,0,0,0,ra=ra,rb=rb,rp=rp,&
                      zetp=mincoreval,eps=eps_rho_rspace,prefactor=1._dp,cutoff=0._dp) )
      p = 1._dp
    END IF

  END DO

  DEALLOCATE (p,STAT=stat)
  IF ( stat /= 0 ) CALL stop_memory ("get_subgrid_size","p")

END SUBROUTINE get_subgrid_size 

!***************************************************************************

END MODULE pw_env_methods
