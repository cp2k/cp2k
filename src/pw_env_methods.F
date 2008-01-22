!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!


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
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             rs_pool_create,&
                                             rs_pools_create_rs_vect,&
                                             rs_pools_dealloc,&
                                             rs_pools_give_back_rs_vect
  USE cube_utils,                      ONLY: destroy_cube_info,&
                                             init_cube_info,&
                                             return_cube_max_iradius
  USE f77_blas
  USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                             gaussian_gridlevel,&
                                             init_gaussian_gridlevel
  USE input_constants,                 ONLY: use_mt
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_type
  USE pw_grid_info,                    ONLY: pw_find_cutoff,&
                                             pw_grid_init_setup
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_change,&
                                             pw_grid_create,&
                                             pw_grid_release,&
                                             pw_grid_setup
  USE pw_poisson_methods,              ONLY: pw_poisson_set
  USE pw_poisson_types,                ONLY: pw_poisson_create
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_p_type,&
                                             pw_pool_release,&
                                             pw_pool_retain,&
                                             pw_pools_dealloc
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                             rho0_mpole_type
  USE qs_util,                         ONLY: exp_radius
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             rs_grid_print,&
                                             realspace_grid_input_type,&
                                             init_input_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_env_methods'

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
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(.NOT.ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(pw_env, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     NULLIFY(pw_env%pw_pools, pw_env%gridlevel_info,pw_env%poisson_env, &
          pw_env%cube_info, pw_env%rs_pools,&
          pw_env%xc_pw_pool)
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
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_env_rebuild', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: safety_factor = 1.01

    INTEGER :: grid_span, handle, i, igrid_level, igrid_zet0_s, ikind, &
      iounit, ipgf, iset, ishell, jkind, jpgf, jset, jshell, la, lb, &
      lgrid_level, ncommensurate, ngrid_level, nkind, nseta, nsetb, nsmax, &
      stat
    INTEGER, DIMENSION(:), POINTER           :: npgfa, npgfb, nshella, nshellb
    INTEGER, DIMENSION(:, :), POINTER        :: lshella, lshellb
    LOGICAL                                  :: failure, odd, should_output, &
                                                spherical, uf_grid, &
                                                use_ref_cell
    REAL(KIND=dp)                            :: alpha, core_charge, cutilev, &
                                                max_rpgf0_s, maxradius, &
                                                my_cut, rel_cutoff, zet0_s, &
                                                zetp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: radius
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cutoff, my_cutoff
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell, cell_ref, my_cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(pw_grid_type), POINTER              :: mt_super_ref_grid, &
                                                old_pw_grid, pw_grid, &
                                                super_ref_grid, &
                                                xc_super_ref_grid
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    TYPE(section_vals_type), POINTER         :: input, poisson_section, &
                                                print_section, rs_grid_section
    TYPE(realspace_grid_input_type)          :: input_settings

! a very small safety factor might be needed for roundoff issues
! e.g. radius being computed here as 12.998 (13) and 13.002 (14) during the collocation
! te latter can happen due to the lower precision in the computation of the radius in collocate
! parallel cost of rs_pw_transfer goes as safety_factor**3 so it is worthwhile keeping it tight

  CALL timeset(routineN,"I","",handle)

!
!
! Part one, deallocate old data if needed
!
!
  failure=.FALSE.
  NULLIFY(my_cutoff,cutoff,cell,pw_grid,old_pw_grid,dft_control,atomic_kind_set,&
          pw_pools,rho0_mpole,rs_pools, para_env,cell_ref, mt_super_ref_grid,&
          input,poisson_section,xc_super_ref_grid,super_ref_grid,my_cell)

  CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL pw_pool_release(pw_env%xc_pw_pool,error=error)
     CALL pw_pools_dealloc(pw_env%pw_pools,error=error)
     CALL rs_pools_dealloc(pw_env%rs_pools,error=error)
     IF (ASSOCIATED(pw_env%gridlevel_info)) THEN
        CALL destroy_gaussian_gridlevel(pw_env%gridlevel_info,error=error)
     ELSE
        ALLOCATE(pw_env%gridlevel_info,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF(qs_env%dft_control%qs_control%gapw) THEN
        CALL get_qs_env(qs_env=qs_env,rho0_mpole=rho0_mpole,error=error)
        CPPostcondition(ASSOCIATED(rho0_mpole),cp_failure_level,routineP,error,failure)
        CALL get_rho0_mpole(rho0_mpole=rho0_mpole,&
                         zet0_s=zet0_s,max_rpgf0_s=max_rpgf0_s)
     ELSE
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


!
!
! Part two, setup the pw_grids
!
!

  IF (.NOT.failure) THEN
     CALL get_qs_env(qs_env=qs_env,&
          dft_control=dft_control,atomic_kind_set=atomic_kind_set,&
          cell_ref=cell_ref,cell=cell,para_env=para_env, input=input,&
          error=error)
     CALL get_qs_env ( qs_env, use_ref_cell = use_ref_cell ,error=error)
     IF (use_ref_cell) THEN
       my_cell => cell_ref
     ELSE
       my_cell => cell
     END IF
     rel_cutoff = dft_control%qs_control%relative_cutoff
     cutoff => dft_control%qs_control%e_cutoff
     CALL section_vals_val_get(input,"DFT%XC%XC_GRID%USE_FINER_GRID",&
          l_val=uf_grid,error=error)
     ngrid_level = SIZE(cutoff)


     ! init gridlevel_info XXXXXXXXX setup mapping to the effective cutoff ?
     !                     XXXXXXXXX the cutoff array here is more a 'wish-list'
     !                     XXXXXXXXX same holds for radius
     print_section=>section_vals_get_subs_vals(input, &
                        "PRINT%GRID_INFORMATION",error=error)
     CALL init_gaussian_gridlevel(pw_env%gridlevel_info,&
          ngrid_levels=ngrid_level,cutoff=cutoff,rel_cutoff=rel_cutoff, &
          print_section=print_section,error=error)
     ! init pw_grids and pools
     ALLOCATE(pw_pools(ngrid_level),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     IF (dft_control % qs_control % commensurate_mgrids) THEN
        ncommensurate=ngrid_level
     ELSE
        ncommensurate=0
     ENDIF
     !
     ! If Tuckerman is present let's perform the set-up of the super-reference-grid
     !
     cutilev = cutoff(1)
     IF ( dft_control % qs_control % pw_grid_opt % spherical ) THEN
        grid_span = HALFSPACE
        spherical = .TRUE.
     ELSE IF ( dft_control % qs_control % pw_grid_opt % fullspace ) THEN
        grid_span = FULLSPACE
        spherical = .FALSE.
     ELSE
        grid_span = HALFSPACE
        spherical = .FALSE.
     END IF

     CALL setup_super_ref_grid(super_ref_grid,mt_super_ref_grid,&
          xc_super_ref_grid, cutilev, grid_span, spherical, my_cell, para_env, &
          qs_env%input, ncommensurate, uf_grid=uf_grid,&
          print_section=print_section,error=error)
     old_pw_grid => super_ref_grid
     !
     ! Setup of the multi-grid pw_grid and pw_pools
     !
     logger => cp_error_get_logger(error)
     iounit = cp_print_key_unit_nr(logger,print_section,"",&
              extension=".Log",error=error)

     IF ( dft_control % qs_control % pw_grid_opt % spherical ) THEN
        grid_span = HALFSPACE
        spherical = .TRUE.
        odd = .TRUE.
     ELSE IF ( dft_control % qs_control % pw_grid_opt % fullspace ) THEN
        grid_span = FULLSPACE
        spherical = .FALSE.
        odd = .FALSE.
     ELSE
        grid_span = HALFSPACE
        spherical = .FALSE.
        odd = .TRUE.
     END IF

     DO igrid_level=1,ngrid_level
        CALL pw_grid_create(pw_grid,para_env%group,error=error)

        cutilev = cutoff(igrid_level)

        IF (igrid_level == 1) THEN
           IF (ASSOCIATED(old_pw_grid)) THEN
              CALL pw_grid_setup(my_cell,pw_grid,grid_span=grid_span,&
                   cutoff=cutilev,&
                   spherical=spherical,odd=odd,fft_usage=.TRUE.,&
                   ncommensurate=ncommensurate,icommensurate=igrid_level,&
                   blocked=.FALSE.,&
                   ref_grid=old_pw_grid,&
                   rs_dims=(/para_env%num_pe,1/),&
                   iounit=iounit,error=error)
              old_pw_grid => pw_grid
           ELSE
              CALL pw_grid_setup(my_cell,pw_grid,grid_span=grid_span,&
                   cutoff=cutilev,&
                   spherical=spherical,odd=odd,fft_usage=.TRUE.,&
                   ncommensurate=ncommensurate,icommensurate=igrid_level,&
                   blocked=.FALSE.,&
                   rs_dims=(/para_env%num_pe,1/),&
                   iounit=iounit,error=error)
              old_pw_grid => pw_grid
           END IF
        ELSE
           CALL pw_grid_setup(my_cell,pw_grid,grid_span=grid_span,&
                cutoff=cutilev,&
                spherical=spherical,odd=odd,fft_usage=.TRUE.,&
                ncommensurate=ncommensurate,icommensurate=igrid_level,&
                blocked=.FALSE.,&
                ref_grid=old_pw_grid,&
                rs_dims=(/para_env%num_pe,1/),&
                iounit=iounit,error=error)
        END IF

      ! init pw_pools
        NULLIFY(pw_pools(igrid_level)%pool)
        CALL pw_pool_create(pw_pools(igrid_level)%pool,pw_grid=pw_grid,error=error)

        CALL pw_grid_release(pw_grid,error=error)

     END DO

     pw_env%pw_pools => pw_pools

     CALL cp_print_key_finished_output(iounit,logger,print_section,&
            "",error=error)

     ! init auxbas_grid
     DO i=1,ngrid_level
        IF (cutoff(i) == dft_control%qs_control%cutoff) pw_env%auxbas_grid=i
     END DO

     ! init xc_pool
     IF (ASSOCIATED(xc_super_ref_grid)) THEN
        CALL pw_pool_create(pw_env%xc_pw_pool,&
             pw_grid=xc_super_ref_grid,error=error)
        CALL pw_grid_release(xc_super_ref_grid,error=error)
     ELSE
        pw_env%xc_pw_pool => pw_pools(pw_env%auxbas_grid)%pool
        CALL pw_pool_retain(pw_env%xc_pw_pool,error=error)
     END IF

     ! complete init of the poisson_env
     IF (.NOT.ASSOCIATED(pw_env%poisson_env)) THEN
        CALL pw_poisson_create(pw_env%poisson_env,error=error)
     END IF
     poisson_section => section_vals_get_subs_vals(input,"DFT%POISSON",&
          error=error)
     CALL pw_poisson_set(pw_env%poisson_env,cell=my_cell,pw_pools=pw_env%pw_pools,&
          parameters=poisson_section,mt_super_ref_pw_grid=mt_super_ref_grid,&
          use_level=pw_env%auxbas_grid,error=error)
     CALL pw_grid_release(mt_super_ref_grid,error=error)
!
! If reference cell is present, then use pw_grid_change to keep bounds constant...
!
!
    IF ( use_ref_cell ) THEN
         ALLOCATE( my_cutoff ( SIZE ( cutoff ) ),stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         DO igrid_level = 1, SIZE ( pw_pools )
           CALL pw_grid_change ( cell, pw_pools ( igrid_level ) % pool % pw_grid )
           my_cut = pw_find_cutoff ( pw_pools ( igrid_level ) % pool % pw_grid % npts,  &
                                 cell%h_inv, error=error)
           my_cutoff (igrid_level) = 0.5_dp * my_cut*my_cut
         ENDDO
         CALL destroy_gaussian_gridlevel (pw_env%gridlevel_info, error=error )
         print_section=>section_vals_get_subs_vals(qs_env%input, &
                             "PRINT%GRID_INFORMATION",error=error)
         CALL init_gaussian_gridlevel(pw_env%gridlevel_info,&
              ngrid_levels=ngrid_level,cutoff=my_cutoff,rel_cutoff=rel_cutoff, &
              print_section=print_section,error=error)
         CALL pw_poisson_set(pw_env%poisson_env,cell=cell,pw_pools=pw_env%pw_pools,&
          parameters=poisson_section,mt_super_ref_pw_grid=mt_super_ref_grid,&
          use_level=pw_env%auxbas_grid,error=error)
         DEALLOCATE( my_cutoff ,stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
!
!
!    determine the maximum radii for mapped gaussians, needed to
!    set up distributed rs grids
!
!
     ALLOCATE(radius(ngrid_level),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     ! Find the grid level suitable for rho0_soft
     IF(qs_env%dft_control%qs_control%gapw) THEN
        igrid_zet0_s = gaussian_gridlevel(pw_env%gridlevel_info,2.0_dp*zet0_s)
        rho0_mpole%igrid_zet0_s = igrid_zet0_s
     END IF

     CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,error=error)
     nkind=SIZE(atomic_kind_set)

     ! XXXXXXXXXXXXXX this duplicates the functionality of get_subgrid_size  XXXXXXXXXXXXXXXXXX
     ! try to predict the maximum radius of the gaussians to be mapped on the grid
     ! up to now, it is not yet very good
     maxradius=0.0_dp
     DO igrid_level=1,ngrid_level
        IF (.NOT. dft_control%qs_control%map_paa) maxradius=0.0_dp

        ! Take into account the radius of the soft compensation charge rho0_soft1
        IF(qs_env%dft_control%qs_control%gapw) THEN
           IF(igrid_zet0_s == igrid_level ) maxradius=MAX(maxradius,max_rpgf0_s)
        END IF

        DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)

          IF (qs_env%dft_control%qs_control%semi_empirical .OR. &
              qs_env%dft_control%qs_control%dftb) CYCLE

          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               orb_basis_set=orb_basis_set,&
                               alpha_core_charge=alpha,ccore_charge=core_charge)
           ! this is to be sure that the core charge is mapped ok
           ! right now, the core is mapped on the auxiliary basis,
           ! this should, at a give point be changed
           ! so that also for the core a multigrid is used
           IF (alpha > 0.0_dp .AND. core_charge.NE.0.0_dp) THEN
             maxradius=MAX(maxradius,exp_radius( 0, alpha, &
                           dft_control%qs_control%eps_rho_rspace, core_charge))
             ! forces
             maxradius=MAX(maxradius,exp_radius( 1, alpha, &
                           dft_control%qs_control%eps_rho_rspace, core_charge))
           ENDIF
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
                       ! density (scale is at most 2)
                       maxradius=MAX(maxradius,exp_radius( lb, zetp, &
                                     dft_control%qs_control%eps_rho_rspace, 2.0_dp))

                       ! should do tau (if tau will be used)

                       ! potential
                       maxradius=MAX(maxradius,exp_radius( lb, zetp, &
                                     dft_control%qs_control%eps_gvg_rspace, 2.0_dp))
                       ! forces
                       maxradius=MAX(maxradius,exp_radius( lb+1, zetp, &
                                     dft_control%qs_control%eps_gvg_rspace, 2.0_dp))

                       ! should do stresses (if stresses will be used)
                    ENDIF
                 ENDDO
                 ENDDO
                 ENDDO
              ENDDO
              ENDDO
              ENDDO
         ENDDO
        ENDDO
        ! this is a bit of hack, but takes into account numerics and rounding
        maxradius = maxradius * safety_factor
        radius(igrid_level)=maxradius
     END DO

!
!
!    set up the rs_grids and the cubes, requires 'radius' to be set up correctly
!
!
     ALLOCATE(rs_pools(ngrid_level),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     ALLOCATE (pw_env%cube_info(ngrid_level),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     DO igrid_level=1,ngrid_level
        pw_grid => pw_pools(igrid_level)%pool%pw_grid
         
        CALL init_cube_info(pw_env%cube_info(igrid_level),&
             pw_grid%dr(:), pw_grid%dh(:,:), pw_grid%dh_inv(:,:), pw_grid%orthorhombic,&
             radius(igrid_level))

        rs_grid_section=>section_vals_get_subs_vals(input,"DFT%MGRID%RS_GRID",error=error)

        CALL init_input_type(input_settings,nsmax=2*MAX(1,return_cube_max_iradius(pw_env%cube_info(igrid_level)))+1,&
                                 rs_grid_section=rs_grid_section,ilevel=igrid_level,error=error)

        NULLIFY(rs_pools(igrid_level)%pool)
        CALL rs_pool_create(rs_pools(igrid_level)%pool,&
             el_struct=pw_grid,input_settings=input_settings,error=error)
     END DO
     pw_env%rs_pools => rs_pools

     DEALLOCATE(radius,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

!
!    Print grid information
!
     logger => cp_error_get_logger(error)
     iounit = cp_print_key_unit_nr(logger,print_section,"",&
          extension=".Log",error=error)
     should_output=BTEST(cp_print_key_should_output(logger%iter_info,&
          print_section,"",error=error),cp_p_file)

     IF ( should_output ) THEN
       CALL rs_pools_create_rs_vect(rs_pools,rs,error=error)
       DO igrid_level=1,ngrid_level
          CALL rs_grid_print(rs(igrid_level)%rs_grid,iounit,error=error)
       END DO
       CALL rs_pools_give_back_rs_vect(rs_pools,rs,error=error)
     END IF
     CALL cp_print_key_finished_output(iounit,logger,print_section,&
         "",error=error)

  END IF

  CALL timestop(0.0_dp,handle)

END SUBROUTINE pw_env_rebuild

!***************************************************************************
  !!****** pw_poisson_types/setup_super_ref_grid [1.0] *
  !!
  !!   NAME
  !!     setup_super_ref_grid
  !!
  !!   FUNCTION
  !!     Initialize the super-reference grid for Tuckerman or xc
  !!
  !!   AUTHOR
  !!     03-2005 Teodoro Laino [teo]
  !!
  !!   NOTE
  !!     move somewere else?
  !!
  !!*** ***********************************************************************
  SUBROUTINE setup_super_ref_grid(super_ref_pw_grid,mt_super_ref_pw_grid,&
       xc_super_ref_pw_grid, cutilev, grid_span, spherical,&
       cell_ref, para_env, input, my_ncommensurate, uf_grid,print_section, &
       error)
    TYPE(pw_grid_type), POINTER              :: super_ref_pw_grid, &
                                                mt_super_ref_pw_grid, &
                                                xc_super_ref_pw_grid
    REAL(KIND=dp), INTENT(IN)                :: cutilev
    INTEGER, INTENT(IN)                      :: grid_span
    LOGICAL, INTENT(in)                      :: spherical
    TYPE(cell_type), POINTER                 :: cell_ref
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    INTEGER, INTENT(IN)                      :: my_ncommensurate
    LOGICAL, INTENT(in)                      :: uf_grid
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_super_ref_grid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iounit, my_val, nn(3), no(3)
    LOGICAL                                  :: failure, mt_s_grid
    REAL(KIND=dp)                            :: mt_rel_cutoff, my_cutilev
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: poisson_section

    failure = .FALSE.
    NULLIFY(poisson_section)
    CPPrecondition(.NOT.ASSOCIATED(mt_super_ref_pw_grid),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(xc_super_ref_pw_grid),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(super_ref_pw_grid),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       poisson_section => section_vals_get_subs_vals(input,"DFT%POISSON",error=error)
       CALL section_vals_val_get(poisson_section,"POISSON_SOLVER",i_val=my_val,error=error)
       !
       ! Check if grids will be the same... In this case we don't use a super-reference grid
       !
       mt_s_grid=.FALSE.
       IF (my_val==use_mt) THEN
          CALL section_vals_val_get(poisson_section,"MT%REL_CUTOFF",&
               r_val=mt_rel_cutoff,error=error)
          IF (mt_rel_cutoff>1._dp) mt_s_grid=.TRUE.
       END IF

       logger => cp_error_get_logger(error)
       iounit = cp_print_key_unit_nr(logger,print_section,"",&
            extension=".Log",error=error)

       IF (uf_grid) THEN
          CALL pw_grid_create(xc_super_ref_pw_grid,para_env%group,error=error)
          CALL pw_grid_setup(cell_ref, xc_super_ref_pw_grid, grid_span=grid_span,&
               cutoff=4._dp*cutilev,spherical=spherical,odd=.FALSE., fft_usage=.TRUE., &
               ncommensurate=my_ncommensurate,icommensurate=1,&
               blocked=.FALSE., rs_dims=(/para_env%num_pe,1/),&
               iounit=iounit,error=error)
          super_ref_pw_grid => xc_super_ref_pw_grid
       END IF
       IF (mt_s_grid) THEN
          CALL pw_grid_create(mt_super_ref_pw_grid,para_env%group,error=error)

          IF (ASSOCIATED(xc_super_ref_pw_grid)) THEN
             CALL cp_unimplemented_error(routineP,&
                  "special grid for mt and fine xc grid not compatible",&
                  error=error)
          ELSE
             my_cutilev=cutilev*mt_rel_cutoff

             no = pw_grid_init_setup(cell_ref%hmat,cutoff=cutilev,spherical=spherical,&
                  odd=.FALSE.,fft_usage=.TRUE.,ncommensurate=0,icommensurate=1,error=error)
             nn = pw_grid_init_setup(cell_ref%hmat,cutoff=my_cutilev,spherical=spherical,&
                  odd=.FALSE.,fft_usage=.TRUE.,ncommensurate=0,icommensurate=1,error=error)

             ! bug appears for nn==no, also in old versions
             CPPrecondition(ALL(nn>no),cp_failure_level,routineP,error,failure)
             CALL pw_grid_setup(cell_ref, mt_super_ref_pw_grid, &
                  cutoff=my_cutilev,spherical=spherical,fft_usage=.TRUE., &
                  blocked=.FALSE., rs_dims=(/para_env%num_pe,1/),&
                  iounit=iounit,error=error)
             super_ref_pw_grid => mt_super_ref_pw_grid
          END IF
       END IF
       CALL cp_print_key_finished_output(iounit,logger,print_section,&
         "",error=error)
    END IF
  END SUBROUTINE setup_super_ref_grid
!***************************************************************************

END MODULE pw_env_methods

!***************************************************************************
