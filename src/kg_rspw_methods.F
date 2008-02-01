!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief rho_methods
!> \par History
!>      JGH (22-Feb-03) PW grid options added
!>      gt 16-nov-03 moved initialization in this new module
! *****************************************************************************
MODULE kg_rspw_methods
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
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_rspw_types,                   ONLY: kg_rspw_type
  USE kinds,                           ONLY: dp
  USE pw_grid_info,                    ONLY: pw_find_cutoff
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
                                             pw_pools_dealloc
  USE qs_util,                         ONLY: exp_radius
  USE realspace_grid_types,            ONLY: init_input_type,&
                                             realspace_grid_input_type,&
                                             realspace_grid_p_type,&
                                             rs_grid_print
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

! *****************************************************************************
SUBROUTINE kg_rspw_create(kg_rspw,kg_env,error)
    TYPE(kg_rspw_type), POINTER              :: kg_rspw
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rspw_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  ALLOCATE(kg_rspw, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     NULLIFY( kg_rspw%gridlevel_info, kg_rspw%pw_pools,kg_rspw%rs_pools,&
          kg_rspw%cube_info, kg_rspw%poisson_env)
     kg_rspw%auxbas_grid=-1
     kg_rspw%ref_count=1
     IF (PRESENT(kg_env)) CALL kg_rspw_rebuild(kg_rspw,kg_env=kg_env,error=error)
  END IF
END SUBROUTINE kg_rspw_create

! *****************************************************************************
!> \param kg_rspw the density/potential env to be initialized
!> \param kg_env the kg environment
!> \param error error
!> \author gloria
! *****************************************************************************
  SUBROUTINE kg_rspw_rebuild(kg_rspw, kg_env, error )

    TYPE(kg_rspw_type), POINTER              :: kg_rspw
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_rspw_rebuild', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: safety_factor = 1.2

    INTEGER :: grid_span, handle, i, igrid_level, ikind, iounit, ipgf, iset, &
      ishell, la, lgrid_level, ncommensurate, ngrid_level, nkind, nseta, &
      nsmax, stat
    INTEGER, DIMENSION(:), POINTER           :: npgfa, nshella
    INTEGER, DIMENSION(:, :), POINTER        :: lshella
    LOGICAL                                  :: failure, odd, should_output, &
                                                spherical, use_ref_cell
    REAL(dp)                                 :: alpha, core_charge, cutilev, &
                                                maxradius, my_cut, &
                                                rel_cutoff, zetp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: radius
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cutoff, my_cutoff
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zeta
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
    TYPE(pw_grid_type), POINTER              :: old_pw_grid, pw_grid
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(realspace_grid_input_type)          :: input_settings
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs
    TYPE(section_vals_type), POINTER         :: poisson_section, &
                                                print_section, rs_grid_section

! KG doesn't setup the radius correctly. It appears that not all primitives that are actually mapped
! are being tested in the radius calculation. This safety_factor will lead to poorer parallel performance
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
                    para_env=para_env,error=error)

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
     CALL get_kg_env ( kg_env, use_ref_cell = use_ref_cell ,error=error)
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
           ngrid_level,cutoff,rel_cutoff,print_section=print_section,error=error)
     ! init pw_grids and pools
     ALLOCATE(pw_pools(ngrid_level),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

     IF (dft_control % qs_control % commensurate_mgrids) THEN
       ncommensurate=ngrid_level
     ELSE
       ncommensurate=0
     ENDIF

     logger => cp_error_get_logger(error)
     iounit = cp_print_key_unit_nr(logger,print_section,"",&
          extension=".Log",error=error)
     DO igrid_level=1,ngrid_level
        CALL pw_grid_create(pw_grid,para_env%group,error=error)

        cutilev = cutoff(igrid_level)
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
        ENDIF

      ! init pw_pools
         NULLIFY(pw_pools(igrid_level)%pool)
         CALL pw_pool_create(pw_pools(igrid_level)%pool,&
              pw_grid=pw_grid,error=error)
         CALL pw_grid_release(pw_grid,error=error)

     END DO
     CALL cp_print_key_finished_output(iounit,logger,print_section,&
         "",error=error)

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
             my_cut = pw_find_cutoff ( pw_pools ( igrid_level ) % pool % pw_grid % npts,  &
                                   cell%h_inv, error=error)
             my_cutoff (igrid_level) = 0.5_dp * my_cut * my_cut
           ENDDO
           CALL destroy_gaussian_gridlevel (kg_rspw%gridlevel_info, error=error )
           print_section=>section_vals_get_subs_vals(kg_env%input, &
                               "PRINT%GRID_INFORMATION",error=error)
           CALL init_gaussian_gridlevel(kg_rspw%gridlevel_info,&
                ngrid_levels=ngrid_level,cutoff=my_cutoff,rel_cutoff=rel_cutoff, &
                print_section=print_section,error=error)
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
      ALLOCATE(radius(ngrid_level),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      ALLOCATE (kg_rspw%cube_info(ngrid_level),STAT=stat)
      IF (stat /= 0) CALL stop_memory(routineP,"cube_info",0)

      CALL get_kg_env(kg_env=kg_env, atomic_kind_set=atomic_kind_set,error=error)
      nkind=SIZE(atomic_kind_set)

      maxradius=0.0_dp
      DO igrid_level=1,ngrid_level
        IF (.NOT. dft_control%qs_control%map_paa) maxradius=0.0_dp

        DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)

          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               orb_basis_set=orb_basis_set,&
                               alpha_core_charge=alpha, core_charge=core_charge)
           ! this is to be sure that the core charge is mapped ok
           ! right now, the core is mapped on the auxiliary basis,
           ! this should, at a give point be changed
           ! so that also for the core a multigrid is used
           IF (alpha > 0.0_dp  .AND. core_charge.NE.0.0_dp) THEN
             maxradius=MAX(maxradius,exp_radius( 0, alpha, &
                           dft_control%qs_control%eps_rho_rspace, 10.0_dp))
             ! forces
             maxradius=MAX(maxradius,exp_radius( 1, alpha, &
                           dft_control%qs_control%eps_rho_rspace, 10.0_dp))
           ENDIF
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

                     !forces
                     maxradius=MAX(maxradius,exp_radius( la+1, zetp, &
                                  dft_control%qs_control%eps_gvg_rspace, 1.0_dp))
                   ENDIF
                 END DO
              END DO
           END DO
         END DO
        ! safety first, not very efficient, nor general
        ! one could possibly decide that this *is* the maximum allowed radius
         maxradius = maxradius * safety_factor
         radius(igrid_level)=maxradius
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
        NULLIFY(rs_pools(igrid_level)%pool)

        CALL init_cube_info(kg_rspw%cube_info(igrid_level),&
               pw_grid%dr(:), pw_grid%dh(:,:), pw_grid%dh_inv(:,:), pw_grid%orthorhombic,&
               radius(igrid_level))

        rs_grid_section=>section_vals_get_subs_vals(kg_env%input,"DFT%MGRID%RS_GRID",error=error)

        CALL init_input_type(input_settings,nsmax=2*MAX(1,return_cube_max_iradius(kg_rspw%cube_info(igrid_level)))+1,&
                                 rs_grid_section=rs_grid_section,ilevel=igrid_level,error=error)

        CALL rs_pool_create(rs_pools(igrid_level)%pool,&
            el_struct=pw_grid,input_settings=input_settings,error=error)
      ENDDO
      kg_rspw%rs_pools => rs_pools

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

  END SUBROUTINE kg_rspw_rebuild

END MODULE kg_rspw_methods

