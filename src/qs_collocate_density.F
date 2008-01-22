 !-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_collocate_density [1.0] *
!!
!!   NAME
!!     qs_collocate_density
!!
!!   FUNCTION
!!     Calculate the plane wave density by collocating the primitive Gaussian
!!     functions (pgf).
!!
!!   AUTHOR
!!     Matthias Krack (03.04.2001)
!!     1) Joost VandeVondele (01.2002)
!!     Thomas D. Kuehne (04.08.2005)
!!
!!   MODIFICATION HISTORY
!!     - rewrote collocate for increased accuracy and speed
!!     - introduced the PGI hack for increased speed with that compiler
!!       (22.02.02)
!!     - Added Multiple Grid feature
!!     - new way to get over the grid (01.03.02)
!!     - removed timing calls since they were getting expensive
!!     - Updated with the new QS data structures (09.04.02,MK)
!!     - introduction of the real space grid type ( prelim. version JVdV 05.02)
!!     - parallel FFT (JGH 22.05.02)
!!     - multigrid arrays independent from density (JGH 30.08.02)
!!     - old density stored in g space (JGH 30.08.02)
!!     - distributed real space code (JGH 17.07.03)
!!     - refactoring and new loop ordering (JGH 23.11.03)
!!     - OpenMP parallelization (JGH 03.12.03)
!!     - Modified to compute tau (Joost 12.03)
!!     - removed the incremental density rebuild (Joost 01.04)
!!     - introduced realspace multigridding (Joost 02.04)
!!     - introduced map_consistent (Joost 02.04)
!!     - Addition of the subroutine calculate_atomic_charge_density (TdK, 08.05)
!!     - rewrite of the collocate/integrate kernels (Joost VandeVondele, 03.07)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_collocate_density
! *****************************************************************************
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs,&
                                             rs_pools_create_rs_vect,&
                                             rs_pools_give_back_rs_vect
  USE cube_utils,                      ONLY: cube_info_type,&
                                             return_cube,&
                                             return_cube_nonortho
  USE distribution_2d_types,           ONLY: distribution_2d_release
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE input_constants,                 ONLY: pw_interp,&
                                             spline3_pbc_interp
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp, int_8,&
                                             dp_size
  USE mathconstants,                   ONLY: pi,&
                                             twopi, &
                                             fac
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum, &
                                             mp_sync
  USE orbital_pointers,                ONLY: coset,&
                                             indco,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_create_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_give_back_pws,&
                                             pw_pools_create_pws
  USE pw_spline_utils,                 ONLY: pw_prolongate_s3
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_copy,&
                                             pw_integrate_function,&
                                             pw_axpy,&
                                             pw_transfer,&
                                             pw_p_type,&
                                             pw_type,&
                                             pw_zero
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_interactions,                 ONLY: exp_radius_very_extended
  USE qs_modify_pab_block,             ONLY: FUNC_AB,&
                                             FUNC_ADBmDAB,&
                                             FUNC_ARDBmDARB,&
                                             FUNC_DADB,&
                                             prepare_adb_m_dab,&
                                             prepare_ardb_m_darb,&
                                             prepare_dadb
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs2pw,&
                                             rs_grid_zero,&
                                             rs_pw_transfer,&
                                             rs_find_node
  USE realspace_task_selection,        ONLY: rs_get_my_tasks,&
                                             distribute_matrix,&
                                             int2pair
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             deallocate_matrix,&
                                             get_block_node,&
                                             replicate_matrix,&
                                             real_matrix_type
  USE task_list_types,                 ONLY: task_list_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_collocate_density'
! *** Public subroutines ***

  PUBLIC :: calculate_rho_core,&
            calculate_rho_elec,&
            calculate_wavefunction,&
            collocate_pgf_product_gspace,&
            collocate_pgf_product_rspace,&
            collocate_atomic_charge_density,&
            density_rs2pw

! *** Public type ***

  TYPE lgrid_type
     INTEGER :: ldim
     REAL(dp), DIMENSION(:), POINTER :: r
  END TYPE lgrid_type

  PUBLIC :: lgrid_type

  INTEGER :: debug_count=0

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****f* qs_collocate_density/collocate_atomic_charge_density [1.0] *
!!
!!   NAME
!!     collocate_atomic_charge_density
!!
!!   FUNCTION
!!     Collocates an arbitrary density from the aux_basis_set onto a grid.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rho: The PW-Grid onto which the density is collocated
!!     - rho_g: The realspace-grid onto which the density is collocated
!!     - rho_r: The G-space-grid onto which the density is collocated
!!     - total_rho: Gives back the integral of the collocated density
!!     - qs_env: The QS environment of matter
!!     - error: variable to control error logging, stopping,...
!!              see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     08.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE collocate_atomic_charge_density(total_rho, qs_env, error)

    REAL(KIND=dp), INTENT(OUT)               :: total_rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'collocate_atomic_charge_density', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, ierr, igrid_level, ikind, ipgf, iset, &
      ithread, maxco, na1, natom, ncoa, nkind, nseta, sgfa, unit_nr
    INTEGER, DIMENSION(:), POINTER           :: atom_list, la_max, la_min, &
                                                npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: failure, map_consistent
    REAL(KIND=dp)                            :: eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab, sphi_a, work, zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: mgrid_gspace, mgrid_rspace
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho

!TYPE(realspace_grid_type), POINTER       :: rs_rho
!TYPE(cp_rs_pool_type), POINTER           :: auxbas_rs_pool
!TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
!TYPE(pw_p_type)                          :: rhoc_r
!TYPE(cube_info_type), POINTER            :: cube_info
!   ---------------------------------------------------------------------------

    NULLIFY(aux_basis_set, atomic_kind_set, atomic_kind, npgfa, cell, particle_set, &
            sphi_a, rs_rho, pw_env, cube_info, dft_control, rho_struct, atom_list, &
            first_sgfa, pab, work)
    NULLIFY(gridlevel_info, rs_pools, pw_pools, mgrid_gspace, mgrid_rspace)

    CALL timeset("collocate_atomic_charge_density","I"," ",handle)

    !ALLOCATE(pab(1,1),STAT=ierr)
    !IF ( ierr /= 0 ) CALL stop_memory ( routineN, "pab", 1 )

    logger => cp_error_get_logger(error)

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                    particle_set=particle_set, pw_env=pw_env, rho=rho_struct, &
                    dft_control=dft_control, error=error)

    !CALL pw_env_get(pw_env=pw_env, auxbas_rs_pool=auxbas_rs_pool, &
    !                auxbas_pw_pool=auxbas_pw_pool, error=error)

    !CALL rs_pool_create_rs(auxbas_rs_pool, rs_rho)
    !CALL rs_grid_zero(rs_rho)

    !cube_info => pw_env%cube_info(1)
    cube_info => pw_env%cube_info
    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    ithread = 0
    gridlevel_info=>pw_env%gridlevel_info

    ! *** set up the pw multi-grids *** !
    CPPrecondition(ASSOCIATED(pw_env), cp_failure_level, routineN, error, failure)
    CALL pw_env_get(pw_env=pw_env, rs_pools=rs_pools, pw_pools=pw_pools, error=error)

    CALL pw_pools_create_pws(pw_pools, mgrid_rspace, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)

    CALL pw_pools_create_pws(pw_pools, mgrid_gspace, &
                              use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                              error=error)

    ! *** set up the rs multi-grids *** !
    CALL rs_pools_create_rs_vect(rs_pools,rs_rho, error=error)
    DO igrid_level = 1,gridlevel_info%ngrid_levels
      CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
    END DO

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxco=maxco)
    ALLOCATE(pab(maxco,1), STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN, "pab", maxco*1*dp_size)
    ALLOCATE(work(maxco,1), STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN, "work", maxco*1*dp_size)

    nkind = SIZE(atomic_kind_set)

    DO ikind = 1,nkind
      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind, aux_basis_set=aux_basis_set, natom=natom, &
                           atom_list=atom_list)

      CALL get_gto_basis_set(gto_basis_set=aux_basis_set, lmax=la_max, lmin=la_min, zet=zeta, &
                             nset=nseta, npgf=npgfa, sphi=sphi_a, first_sgf=first_sgfa, nsgf_set=nsgfa)

      DO iatom = 1,natom
        ! ra(:) = pbc(particle_set(iatom)%r, cell)
        ra(:) = pbc(particle_set(atom_list(iatom))%r, cell)

        DO iset = 1,nseta

          sgfa = first_sgfa(1,iset)
          ncoa = npgfa(iset)*ncoset(la_max(iset))

          DO i = 1,nsgfa(iset)
            work(i,1) = 1.0_dp
          END DO
          CALL dgemm("N","N",ncoa,1,nsgfa(iset),1.0_dp, sphi_a(1,sgfa),SIZE(sphi_a,1), &
                     work(1,1),SIZE(work,1),0.0_dp,pab(1,1),SIZE(pab,1))

          DO ipgf = 1,npgfa(iset)

            na1 = (ipgf-1)*ncoset(la_max(iset))
            igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf,iset))

            CALL collocate_pgf_product_rspace(la_max=la_max(iset), zeta=zeta(ipgf,iset), &
                                              la_min=la_min(iset),&
                                              lb_max=0, zetb=0.0_dp, lb_min=0,&
                                              ra=ra,rab=(/0.0_dp,0.0_dp,0.0_dp/),rab2=0.0_dp,&
                                              scale=1.0_dp, pab=pab, o1=na1, o2=0,&
                                              rsgrid=rs_rho(igrid_level)%rs_grid,cell=cell,&
                                              cube_info=cube_info(igrid_level),&
                                              eps_rho_rspace=eps_rho_rspace,&
                                              ga_gb_function=FUNC_AB, ithread=ithread, &
                                              map_consistent=map_consistent,error=error)
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE(pab, stat=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN, "pab")

    DEALLOCATE(work, stat=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN, "work")

    IF (gridlevel_info%ngrid_levels==1) THEN
      CALL rs_pw_transfer(rs=rs_rho(1)%rs_grid, pw=qs_env%rho%rho_r(1)%pw, dir=rs2pw)
      CALL rs_pools_give_back_rs_vect(pools=rs_pools, elements=rs_rho, error=error)
      CALL pw_transfer(qs_env%rho%rho_r(1)%pw, qs_env%rho%rho_g(1)%pw)
      IF (qs_env%rho%rho_r(1)%pw%pw_grid%spherical) THEN
        CALL pw_transfer(qs_env%rho%rho_g(1)%pw, qs_env%rho%rho_r(1)%pw)
      END IF
    ELSE
      DO igrid_level = 1,gridlevel_info%ngrid_levels
         CALL rs_pw_transfer(rs=rs_rho(igrid_level)%rs_grid, &
                            pw=mgrid_rspace(igrid_level)%pw, dir=rs2pw)
      END DO
      CALL rs_pools_give_back_rs_vect(pools=rs_pools, elements=rs_rho, error=error)

      CALL pw_zero(qs_env%rho%rho_g(1)%pw)
      DO igrid_level=1, gridlevel_info%ngrid_levels
        CALL pw_transfer(mgrid_rspace(igrid_level)%pw, &
                                   mgrid_gspace(igrid_level)%pw)
        CALL pw_axpy(mgrid_gspace(igrid_level)%pw, qs_env%rho%rho_g(1)%pw)
      END DO
      CALL pw_transfer(qs_env%rho%rho_g(1)%pw, qs_env%rho%rho_r(1)%pw)
    END IF

    total_rho = pw_integrate_function(qs_env%rho%rho_r(1)%pw,isign=-1)
    qs_env%rho%tot_rho_r = total_rho
    IF (logger%para_env%source==logger%para_env%mepos) THEN
      unit_nr=cp_logger_get_default_unit_nr(logger,local=.FALSE.)
      WRITE (unit_nr,*) "Total Rho =", total_rho
    END IF

    ! *** give back the multi-grids *** !
    CALL pw_pools_give_back_pws(pw_pools, mgrid_gspace, error=error)
    CALL pw_pools_give_back_pws(pw_pools, mgrid_rspace, error=error)

    qs_env%rho%rho_r_valid = .TRUE.
    qs_env%rho%rho_g_valid = .TRUE.

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE collocate_atomic_charge_density

!!****f* qs_collocate_density/calculate_rho_core *
!!
!!   NAME
!!     calculate_rho_core
!!
!!   FUNCTION
!!     computes the density of the core charges on the grid
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE calculate_rho_core(rho_core,total_rho,qs_env,error)

    TYPE(pw_p_type), INTENT(INOUT)           :: rho_core
    REAL(KIND=dp), INTENT(OUT)               :: total_rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_core', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, dest, dir, handle, &
                                                iatom, ierr, ikind, ithread, &
                                                j, natom, npme, &
                                                nthread,i,k,n, &
                                                tp(3)
    INTEGER, DIMENSION(:), POINTER           :: atom_list,ub,lb
    INTEGER, DIMENSION(:), POINTER           :: cores
    REAL(KIND=dp)                            :: alpha, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_p_type)                          :: rhoc_r
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_type), POINTER           :: auxbas_rs_pool
    TYPE(cube_info_type)                     :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_type), POINTER       :: rs_rho
    TYPE(lgrid_type)                         :: lgrid

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_rho_core","I"," ",handle)
    NULLIFY(atomic_kind,cell,dft_control,pab,atomic_kind_set,particle_set,&
         atom_list,pw_env,rs_rho,auxbas_rs_pool,auxbas_pw_pool,cores)
    NULLIFY(lgrid%r)
    ALLOCATE(pab(1,1),STAT=ierr)
    IF ( ierr /= 0 ) CALL stop_memory ( "calculate_rho_core", "pab", 1 )

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    para_env=para_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env,auxbas_rs_pool=auxbas_rs_pool,&
         auxbas_pw_pool=auxbas_pw_pool,error=error)
    cube_info=pw_env%cube_info(1)
    ! be careful in parallel nsmax is choosen with multigrid in mind!
    CALL rs_pool_create_rs(auxbas_rs_pool,rs_rho, error=error)
    CALL rs_grid_zero(rs_rho)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           atom_list=atom_list,&
                           alpha_core_charge=alpha,&
                           ccore_charge=pab(1,1))

      IF (alpha == 0.0_dp) CYCLE

      nthread = 1
      ithread=0

      CALL reallocate ( cores, 1, natom )
      npme = 0
      cores = 0

      DO iatom = 1, natom
         atom_a = atom_list(iatom)
         ra(:) = pbc(particle_set(atom_a)%r,cell)
         DO dir = 1,3 
            tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),ra)*rs_rho%npts(dir))
            tp(dir) = MODULO ( tp(dir), rs_rho%npts(dir) )
            tp(dir) = tp(dir) + rs_rho%lb(dir)
         END DO
         IF(rs_rho%distributed)THEN
            IF( ALL (tp.GE.rs_rho%lb_real) .AND. ALL (tp.LE.rs_rho%ub_real) ) THEN
               ! local atom
               CALL rs_find_node(rs_rho,tp,dest)
               npme = npme + 1
               cores (npme) = iatom
            ENDIF
         ELSEIF(rs_rho%parallel) THEN
            ! replicated realspace grid, split the atoms up between procs
            IF (MODULO(natom,rs_rho%group_size) == rs_rho % my_pos ) THEN
               npme = npme + 1
               cores (npme) = iatom
            ENDIF
         ELSE
            ! serial job
            npme = npme + 1
            cores (npme) = iatom            
         ENDIF
      END DO

      IF(npme.GT.0)THEN
         DO j=1,npme

            iatom = cores(j)
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r,cell)

            CALL collocate_pgf_product_rspace(0,alpha,0,0,0.0_dp,0,ra,&
                  (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,-1.0_dp,pab,0,0,rs_rho,&
                  cell,cube_info,eps_rho_rspace,ga_gb_function=FUNC_AB,&
                  ithread=ithread,error=error)

         END DO
      ENDIF

    END DO

    IF (ASSOCIATED(cores)) THEN
      DEALLOCATE (cores,STAT=ierr)
      IF ( ierr /= 0 ) CALL stop_memory ( "calculate_rho_core", "cores" )
    END IF
    DEALLOCATE ( pab, STAT=ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "calculate_rho_core", "pab" )

    CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
         use_data=REALDATA3D,in_space=REALSPACE, error=error)

    CALL rs_pw_transfer(rs_rho,rhoc_r%pw,rs2pw)
    CALL rs_pool_give_back_rs(auxbas_rs_pool,rs_rho,error=error)

    total_rho = pw_integrate_function(rhoc_r%pw,isign=-1)

    CALL pw_transfer(rhoc_r%pw,rho_core%pw)

    CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw, error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_rho_core
!!****f* qs_collocate_density/calculate_rho_elec *
!!
!!   NAME
!!     calculate_rho_elec
!!
!!   FUNCTION
!!     computes the density corresponding to a given density matrix on the grid
!!
!!   NOTES
!!     both rho and rho_gspace contain the new rho
!!     (in real and g-space respectively)
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE calculate_rho_elec(matrix_p, rho,rho_gspace, total_rho,&
       qs_env, soft_valid, compute_tau, error)

    TYPE(real_matrix_type), POINTER          :: matrix_p
    TYPE(pw_p_type), INTENT(INOUT)           :: rho, rho_gspace
    REAL(KIND=dp), INTENT(OUT)               :: total_rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: soft_valid, compute_tau
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
         routine = "SUBROUTINE calculate_rho_elec (MODULE qs_collocate_density)"
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
         max_tasks = 2000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: ab, bcol, brow, curr_tasks, dest, dir, ilevel, &
         ga_gb_function, handle, i, iatom, igrid_level, ijatoms, ijsets, ikind, ikind_old, &
         ilist, inode, ipgf, iset, istat, itask, ithread, j, jatom, jkind, jkind_old, jpgf, &
         jset, k, maxco, maxsgf, maxsgf_set, n, na1, na2, natoms, natom_pairs, nb1, nb2, &
         maxset, maxpgf, ncoa, ncob, nkind, nlist, nnode, nseta, nsetb, nthread, &
         sgfa, sgfb, stat, tp(3), ntasks, ntasks_rep, ntasks_dist,&
         iatom_old, jatom_old, iset_old, jset_old
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb, ub
    INTEGER, DIMENSION(:), POINTER           :: atom_pair_send, atom_pair_recv

    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, &
                                                first_sgfb
    INTEGER(kind=int_8), DIMENSION(:, :),&
         POINTER                             :: tasks
    LOGICAL                                  :: distributed_rs_grids, &
         failure, map_consistent, &
         my_compute_tau, my_soft, &
         atom_pair_changed
    REAL(KIND=dp)                            :: dab, eps_rho_rspace, &
         kind_radius_b, rab2, scale, &
         zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb, rp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block, pab, &
         pblock, rpgfa, rpgfb, sphi_a, &
         sphi_b, work, zeta, zetb, dist_ab
    REAL(KIND=dp), DIMENSION(:, :, :), &
         POINTER                             :: pabt, workt
    TYPE(atomic_kind_type), DIMENSION(:), &
         POINTER                             :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
         POINTER                             :: rs_pools
    TYPE(cube_info_type), DIMENSION(:), &
         POINTER                             :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(lgrid_type)                         :: lgrid
    TYPE(neighbor_list_set_p_type), &
         DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
         POINTER                             :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(real_matrix_type), POINTER          :: deltap
    TYPE(realspace_grid_p_type), &
         DIMENSION(:), POINTER               :: rs_rho
    TYPE(section_vals_type), POINTER         :: input, interp_section
    TYPE(task_list_type), POINTER            :: task_list
    TYPE(task_list_type), POINTER            :: task_list_soft

    !   ---------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(atomic_kind,cell,dft_control,orb_basis_set,sab_orb_neighbor_list,&
         sab_orb_neighbor_node,deltap,atomic_kind_set,sab_orb,particle_set,&
         rs_rho,pw_env,rs_pools,para_env,pblock,dist_ab,&
         set_radius_a,set_radius_b,la_max,la_min,&
         lb_max,lb_min,npgfa,npgfb,nsgfa,nsgfb,p_block,&
         rpgfa,rpgfb,sphi_a,sphi_b,zeta,zetb,first_sgfa,first_sgfb,&
         tasks,pabt,workt)
    NULLIFY(lgrid%r)

    debug_count=debug_count+1

    ! by default, the full density is calculated
    my_soft=.FALSE.
    IF (PRESENT(soft_valid)) my_soft = soft_valid

    ! by default, do not compute the kinetic energy density (tau)
    ! if compute_tau, all grids referening to rho are actually tau
    IF (PRESENT(compute_tau)) THEN
       my_compute_tau = compute_tau
    ELSE
       my_compute_tau = .FALSE.
    ENDIF

    IF (my_compute_tau) THEN
       CALL timeset("calculate_rho_tau","I"," ",handle)
       ga_gb_function = FUNC_DADB
    ELSE
       CALL timeset("calculate_rho_elec","I"," ",handle)
       ga_gb_function = FUNC_AB
    ENDIF

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         particle_set=particle_set,&
         sab_orb=sab_orb,&
         para_env=para_env,&
         task_list=task_list,&
         task_list_soft=task_list_soft,&
         input=input,&
         pw_env=pw_env,error=error)

    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info

    !   *** Allocate work storage ***
    nthread = 1
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxco=maxco,&
         maxsgf=maxsgf,&
         maxsgf_set=maxsgf_set)
    nkind = SIZE(atomic_kind_set)
    CALL reallocate(pabt,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread-1)

    ! find maximum numbers
    nkind = SIZE(atomic_kind_set)
    natoms = SIZE( particle_set )
    maxset=0
    maxpgf=0
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            softb = my_soft, &
            orb_basis_set=orb_basis_set)

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            npgf=npgfa, nset=nseta )

       maxset=MAX(nseta,maxset)
       maxpgf=MAX(MAXVAL(npgfa),maxpgf)
    END DO

    ! get the task lists
    IF (my_soft) task_list=>task_list_soft
    CPPrecondition(ASSOCIATED(task_list),cp_failure_level,routine,error,failure)
    tasks  =>task_list%tasks
    dist_ab=>task_list%dist_ab
    atom_pair_send=>task_list%atom_pair_send
    atom_pair_recv=>task_list%atom_pair_recv
    ntasks=task_list%ntasks
    
    ! *** set up the pw multi-grids
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routine,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, error=error)
    ! *** set up the rs multi-grids
    distributed_rs_grids=.FALSE.
    CALL rs_pools_create_rs_vect(rs_pools, rs_rho, error=error)
    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
       IF ( rs_rho(igrid_level)%rs_grid%distributed ) THEN
          distributed_rs_grids=.TRUE.
       ENDIF
    END DO

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    !   *** Initialize working density matrix ***
    ! distributed rs grids require a matrix that will be changed 
    ! whereas this is not the case for replicated grids
    IF (distributed_rs_grids) THEN
       NULLIFY(deltap)
       CALL replicate_matrix(matrix_p,deltap,target_name="DeltaP",error=error)
       ! this matrix has no strict sparsity pattern in parallel
       deltap%sparsity_id=-1
       CALL distribution_2d_release(deltap%distribution_2d,error=error)
    ELSE
       deltap=>matrix_p
    ENDIF

    ! distribute the matrix
    IF (distributed_rs_grids) THEN
       CALL distribute_matrix (rs_rho(1)%rs_grid, deltap, atom_pair_send, atom_pair_recv, natoms, scatter=.TRUE., error=error)
    ENDIF

    ! map all tasks on the grids

    ithread = 0
    pab => pabt(:,:,ithread)
    work => workt(:,:,ithread)

    iatom_old = -1 ; jatom_old = -1 ; iset_old = -1 ; jset_old = -1 
    ikind_old = -1 ; jkind_old = -1

    loop_tasks: DO itask=1,ntasks

       !decode the atom pair and basis info
       CALL int2pair(tasks(3,itask),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)

       ikind = particle_set(iatom)%atomic_kind%kind_number
       jkind = particle_set(jatom)%atomic_kind%kind_number

       IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old) THEN 

          IF (iatom .NE. iatom_old) ra(:) = pbc(particle_set(iatom)%r,cell)

          IF (iatom <= jatom) THEN
             brow = iatom
             bcol = jatom
          ELSE
             brow = jatom
             bcol = iatom
          END IF

          IF(ikind .NE. ikind_old ) THEN
             CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                  softb = my_soft, &
                  orb_basis_set=orb_basis_set)

             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  first_sgf=first_sgfa,&
                  lmax=la_max,&
                  lmin=la_min,&
                  npgf=npgfa,&
                  nset=nseta,&
                  nsgf_set=nsgfa,&
                  pgf_radius=rpgfa,&
                  set_radius=set_radius_a,&
                  sphi=sphi_a,&
                  zet=zeta)
          ENDIF

          IF (jkind .NE. jkind_old ) THEN
             CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind,&
                  softb = my_soft, &
                  orb_basis_set=orb_basis_set)

             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  first_sgf=first_sgfb,&
                  kind_radius=kind_radius_b,&
                  lmax=lb_max,&
                  lmin=lb_min,&
                  npgf=npgfb,&
                  nset=nsetb,&
                  nsgf_set=nsgfb,&
                  pgf_radius=rpgfb,&
                  set_radius=set_radius_b,&
                  sphi=sphi_b,&
                  zet=zetb)          
          ENDIF

          CALL get_block_node(matrix=deltap,&
               block_row=brow,&
               block_col=bcol,&
               BLOCK=p_block)
          IF (.NOT.ASSOCIATED(p_block)) &
               CALL stop_program(routine,"p_block not associated in deltap")

          iatom_old = iatom
          jatom_old = jatom
          ikind_old = ikind
          jkind_old = jkind
          atom_pair_changed = .TRUE.

       ELSE

          atom_pair_changed = .FALSE.

       ENDIF


       IF (atom_pair_changed .OR. iset_old .NE. iset .OR. jset_old .NE. jset) THEN

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)
          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)

          IF (iatom <= jatom) THEN
             CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                  1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                  p_block(sgfa,sgfb),SIZE(p_block,1),&
                  0.0_dp,work(1,1),maxco)
             CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                  1.0_dp,work(1,1),maxco,&
                  sphi_b(1,sgfb),SIZE(sphi_b,1),&
                  0.0_dp,pab(1,1),maxco)
          ELSE
             CALL dgemm("N","N",ncob,nsgfa(iset),nsgfb(jset),&
                  1.0_dp,sphi_b(1,sgfb),SIZE(sphi_b,1),&
                  p_block(sgfb,sgfa),SIZE(p_block,1),&
                  0.0_dp,work(1,1),maxco)
             CALL dgemm("N","T",ncob,ncoa,nsgfa(iset),&
                  1.0_dp,work(1,1),maxco,&
                  sphi_a(1,sgfa),SIZE(sphi_a,1),&
                  0.0_dp,pab(1,1),maxco)
          END IF

          iset_old = iset
          jset_old = jset

       ENDIF

       rab(:) = dist_ab (:,itask)
       rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       rb(:) = ra(:) + rab(:)
       zetp = zeta(ipgf,iset) + zetb(jpgf,jset)
       dab=sqrt(rab2)
       IF (dab.lt.0.1E0_dp .AND. dft_control%qs_control%map_paa) THEN
          igrid_level = 1
       ELSE
          igrid_level = gaussian_gridlevel(gridlevel_info,zetp)
       ENDIF

       na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
       na2 = ipgf*ncoset(la_max(iset))
       nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
       nb2 = jpgf*ncoset(lb_max(jset))

       ! takes the density matrix symmetry in account, i.e. off-diagonal blocks need to be mapped 'twice'
       IF (iatom == jatom) THEN
          scale = 1.0_dp
       ELSE
          scale = 2.0_dp
       END IF

       IF (iatom <= jatom) THEN
          CALL collocate_pgf_product_rspace(&
               la_max(iset),zeta(ipgf,iset),la_min(iset),&
               lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
               ra,rab,rab2,scale,pab,na1-1,nb1-1,&
               rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
               eps_rho_rspace,&
               ga_gb_function=ga_gb_function, &
               map_consistent=map_consistent,error=error)
       ELSE
          CALL collocate_pgf_product_rspace(&
               lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
               la_max(iset),zeta(ipgf,iset),la_min(iset),&
               rb,-rab,rab2,scale,pab,nb1-1,na1-1,&
               rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
               eps_rho_rspace,&
               ga_gb_function=ga_gb_function, &
               map_consistent=map_consistent,error=error)
       END IF

    END DO loop_tasks


    !   *** Release work storage ***

    IF (distributed_rs_grids) THEN
      CALL deallocate_matrix ( deltap ,error=error)
    ENDIF

    DEALLOCATE (pabt,workt,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pabt,workt")

    CALL mp_sync(para_env%group)

    interp_section => section_vals_get_subs_vals(input,"DFT%MGRID%INTERPOLATOR",error=error)
    CALL density_rs2pw(pw_env,rs_rho,rho,rho_gspace,interp_section=interp_section,error=error)

    total_rho = pw_integrate_function(rho%pw,isign=-1)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_rho_elec

!!*************************************************************************


!!****f* qs_collocate_density/calculate_wavefunction *
!!
!!   NAME
!!     calculate_wavefunction
!!
!!   FUNCTION
!!     maps a given wavefunction on the grid
!!
!!   NOTES
!!     modified calculate_rho_elec, should write the wavefunction represented by
!!     it's presumably dominated by the FFT and the rs->pw and back routines
!!
!!   BUGS ???
!!     does it take correctly the periodic images of the basis functions into account
!!     i.e. is only correct if the basis functions are localised enough to be just in 1 cell ?
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [Joost VandeVondele]
!!     03.2006 made independent of qs_env [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE calculate_wavefunction(mo_vectors,ivector,rho,rho_gspace, &
                   atomic_kind_set,cell,dft_control,particle_set, &
                   pw_env,error)

    TYPE(cp_fm_type), POINTER                :: mo_vectors
    INTEGER                                  :: ivector
    TYPE(pw_p_type), INTENT(INOUT)           :: rho, rho_gspace
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_wavefunction', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, igrid_level, ipgf, iset, istat, maxco, &
      maxsgf_set, na1, na2, nao, natom, ncoa, ncol_global, nseta, offset, &
      sgfa, stat, group, my_pos, group_size, dir
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: failure, local
    REAL(KIND=dp)                            :: dab, eps_rho_rspace, rab2, &
                                                scale, zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvector
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab, sphi_a, work, zeta
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: mgrid_gspace, mgrid_rspace
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(REALSPACE_GRID_P_TYPE), &
      DIMENSION(:), POINTER                  :: rs_rho
    LOGICAL                                  :: map_it_here
    INTEGER, DIMENSION(:), ALLOCATABLE       :: where_is_the_point
    INTEGER, DIMENSION(3)                    :: location, lb, ub, tp

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CALL cp_fm_get_info(matrix=mo_vectors,nrow_global=nao,ncol_global=ncol_global,error=error)

    NULLIFY(eigenvector,  orb_basis_set,&
         pab,work,la_max, la_min,&
         npgfa, nsgfa, &
         sphi_a, zeta, first_sgfa,&
         rs_rho,rs_pools,pw_pools,mgrid_rspace,mgrid_gspace)

    ALLOCATE(eigenvector(nao), stat=istat)
    IF (istat.NE.0) CALL stop_program("calculate_wavefunction","eigenvector")

    ! *** set up the pw multi-grids
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineN,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, pw_pools=pw_pools, &
                            cube_info=cube_info, gridlevel_info=gridlevel_info, error=error)
    CALL pw_pools_create_pws(pw_pools,mgrid_gspace,&
                use_data = COMPLEXDATA1D,&
                in_space = RECIPROCALSPACE, error=error)
    CALL pw_pools_create_pws(pw_pools,mgrid_rspace,&
                use_data = REALDATA3D,&
                in_space = REALSPACE, error=error)

    ! *** set up rs multi-grids
    CALL rs_pools_create_rs_vect(rs_pools, rs_rho, error=error)

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
    END DO

    DO i=1,nao
       CALL cp_fm_get_element(mo_vectors,i,ivector,eigenvector(i),local)
    ENDDO

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

!   *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             natom=natom,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,"pab",maxco*1*dp_size)
    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,"work",maxco*1*dp_size)

    offset=0
    group=mgrid_rspace(1)%pw%pw_grid%para%group
    my_pos=mgrid_rspace(1)%pw%pw_grid%para%my_pos
    group_size=mgrid_rspace(1)%pw%pw_grid%para%group_size
    ALLOCATE(where_is_the_point(0:group_size-1))

    DO iatom=1,natom

      CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               orb_basis_set=orb_basis_set)

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             sphi=sphi_a,&
                             zet=zeta)
      ra(:) = pbc(particle_set(iatom)%r,cell)
      rab(:) = 0.0_dp
      rab2  = 0.0_dp
      dab   = 0.0_dp

      DO iset=1,nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

         DO i=1,nsgfa(iset)
            work(i,1)=eigenvector(offset+i)
         ENDDO

         CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                    1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                    work(1,1),SIZE(work,1),&
                    0.0_dp,pab(1,1),SIZE(pab,1))

         DO ipgf=1,npgfa(iset)

            na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
            na2 = ipgf*ncoset(la_max(iset))

            scale = 1.0_dp
            zetp = zeta(ipgf,iset)
            igrid_level = gaussian_gridlevel(gridlevel_info,zetp)

            map_it_here=.FALSE.

            IF (.NOT. ALL (rs_rho(igrid_level)%rs_grid%perd == 1)) THEN
               DO dir = 1,3
                     ! bounds of local grid (i.e. removing the 'wings')
                     tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),ra)*rs_rho(igrid_level)%rs_grid%npts(dir))
                     tp(dir) = MODULO ( tp(dir), rs_rho(igrid_level)%rs_grid%npts(dir) )
                     lb(dir) = rs_rho(igrid_level)%rs_grid%lb_local ( dir ) + rs_rho(igrid_level)%rs_grid%border
                     ub(dir) = rs_rho(igrid_level)%rs_grid%ub_local ( dir ) - rs_rho(igrid_level)%rs_grid%border
                     ! distributed grid, only map if it is local to the grid
                     location(dir)=tp(dir)+rs_rho(igrid_level)%rs_grid%lb(dir)
               ENDDO
               IF  (lb(1)<=location(1) .AND. location(1)<=ub(1) .AND. &
                    lb(2)<=location(2) .AND. location(2)<=ub(2) .AND. &
                    lb(3)<=location(3) .AND. location(3)<=ub(3)) THEN
                  map_it_here=.TRUE.
                  ! this slice can still be on many CPUs. only one should do the point.
                  ! this should be changed as group_dim is going to!
                  !IF (MODULO(offset,rs_rho(igrid_level)%rs_grid%group_dim ( 2 )) == &
                  !     rs_rho(igrid_level)%rs_grid%group_coor (2 )) map_it_here=.TRUE.
               ENDIF
            ELSE
               ! not distributed, just a round-robin distribution over the full set of CPUs
               IF (MODULO(offset,group_size)==my_pos) map_it_here=.TRUE.
            ENDIF
            
            
            IF (map_it_here) CALL collocate_pgf_product_rspace(&
                 la_max(iset),zeta(ipgf,iset),la_min(iset),&
                 0,0.0_dp,0,&
                 ra,rab,rab2,scale,pab,na1-1,0,&
                 rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
                 eps_rho_rspace,ga_gb_function=FUNC_AB,error=error)

         END DO

         offset=offset+nsgfa(iset)

      END DO

    END DO

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid,&
            mgrid_rspace(igrid_level)%pw,rs2pw)
    ENDDO

    CALL rs_pools_give_back_rs_vect(rs_pools,rs_rho,error=error)

    CALL pw_zero(rho_gspace%pw)
    DO igrid_level=1,gridlevel_info%ngrid_levels
      CALL pw_transfer(mgrid_rspace(igrid_level)%pw,&
           mgrid_gspace(igrid_level)%pw)
      CALL pw_axpy(mgrid_gspace(igrid_level)%pw,rho_gspace%pw)
    END DO

    CALL pw_transfer(rho_gspace%pw,rho%pw)

!   *** Release work storage ***
    DEALLOCATE(eigenvector)

    DEALLOCATE (pab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,"pab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,"work")

    ! *** give back the pw multi-grids
    CALL pw_pools_give_back_pws(pw_pools,mgrid_gspace,error=error)
    CALL pw_pools_give_back_pws(pw_pools,mgrid_rspace,error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_wavefunction

!!****f* qs_collocate_density/collocate_pgf_product_rspace *
!!
!!   NAME
!!     collocate_pgf_product_rspace
!!
!!   FUNCTION
!!     low level collocation of primitive gaussian functions
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE collocate_pgf_product_rspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2,scale,pab,o1,o2,&
                                          rsgrid,cell,cube_info,&
                                          eps_rho_rspace,ga_gb_function,&
                                          lgrid,ithread,&
                                          map_consistent,&
                                          collocate_rho0,&
                                          rpgf0_s,idir,ir,error)

    INTEGER, INTENT(IN)                      :: la_max
    REAL(KIND=dp), INTENT(IN)                :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(KIND=dp), INTENT(IN)                :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rab
    REAL(KIND=dp), INTENT(IN)                :: rab2, scale
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    INTEGER, INTENT(IN)                      :: o1, o2
    TYPE(realspace_grid_type), POINTER       :: rsgrid
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    REAL(KIND=dp), INTENT(IN)                :: eps_rho_rspace
    INTEGER, INTENT(IN)                      :: ga_gb_function
    TYPE(lgrid_type), OPTIONAL               :: lgrid
    INTEGER, INTENT(IN), OPTIONAL            :: ithread
    LOGICAL, INTENT(IN), OPTIONAL            :: map_consistent, collocate_rho0
    REAL(dp), INTENT(IN), OPTIONAL           :: rpgf0_s
    INTEGER, INTENT(IN), OPTIONAL            :: idir, ir
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'collocate_pgf_product_rspace', &
      routineP = moduleN//':'//routineN

    INTEGER :: cmax, gridbounds(2,3), i, ico, icoef, ig, ithread_l, &
      jco, k, l, la_max_local, la_min_local, lb_max_local, lb_min_local, &
      length, lx, lx_max, lxa, lxb, lxy, lxy_max, lxyz, lxyz_max, lya, lyb, &
      lza, lzb, o1_local, o2_local, offset, start
    INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, ng, &
                                                ub_cube
    INTEGER, DIMENSION(:), POINTER           :: ly_max, lz_max, sphere_bounds
    LOGICAL                                  :: failure, my_collocate_rho0, &
                                                my_map_consistent
    REAL(KIND=dp) :: a, b, binomial_k_lxa, binomial_l_lxb, cutoff, f, pg, &
      prefactor, radius, rpg, ya, yap, yb, ybp, za, zap, zb, zbp, zetp
    REAL(KIND=dp), DIMENSION(3)              :: dr, rap, rb, rbp, roffset, rp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab_local
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid

    INTEGER :: lxp,lyp,lzp,lp,lxpm,lypm,iaxis
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: map
    REAL(kind=dp) :: p_ele,ax,ay,az
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: alpha
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coef_xyz
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coef_xyt
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coef_xtt

    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: pol_z
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: pol_y
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: pol_x
    REAL(KIND=dp) :: t_exp_1,t_exp_2,t_exp_min_1,t_exp_min_2,t_exp_plus_1,t_exp_plus_2

    failure = .FALSE.

    IF (PRESENT(ithread)) THEN
       ithread_l=ithread
    ELSE
       ithread_l=0
    ENDIF

    ! use identical radii for integrate and collocate ?
    IF (PRESENT(map_consistent)) THEN
       my_map_consistent=map_consistent
    ELSE
       my_map_consistent=.FALSE.
    ENDIF


    IF (PRESENT(collocate_rho0).AND.PRESENT(rpgf0_s)) THEN
       my_collocate_rho0=collocate_rho0
    ELSE
       my_collocate_rho0=.FALSE.
    END IF

    zetp      = zeta + zetb
    f         = zetb/zetp
    rap(:)    = f*rab(:)
    rbp(:)    = rap(:) - rab(:)
    rp(:)     = ra(:) + rap(:)
    rb(:)     = ra(:)+rab(:)

    IF (my_map_consistent) THEN
       cutoff    = 1.0_dp
       prefactor = EXP(-zeta*f*rab2)
       radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,ra=ra,rb=rb,rp=rp,&
                                       zetp=zetp,eps=eps_rho_rspace,&
                                       prefactor=prefactor,cutoff=cutoff)
       prefactor = scale*EXP(-zeta*f*rab2)
    ELSE IF (my_collocate_rho0) THEN
       cutoff    = 0.0_dp
       prefactor = 1.0_dp
       radius = rpgf0_s
    ELSE
       cutoff    = 0.0_dp
       prefactor = scale*EXP(-zeta*f*rab2)
       radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,o1,o2,ra,rb,rp,&
                                       zetp,eps_rho_rspace,prefactor,cutoff)
    ENDIF

    IF (radius .EQ. 0.0_dp ) THEN
      RETURN
    END IF

    ! it's a choice to compute lX_min/max, pab here,
    ! this way we get the same radius as we use for the corresponding density
    SELECT CASE (ga_gb_function)
    CASE(FUNC_DADB)
        la_max_local=la_max+1
        la_min_local=MAX(la_min-1,0)
        lb_max_local=lb_max+1
        lb_min_local=MAX(lb_min-1,0)
        ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
        ! is equivalent to mapping pab with 0.5 * (nabla pgf_a) . (nabla pgf_b)
        ! (ddx pgf_a ) (ddx pgf_b) = (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x})*(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x})
        ! cleaner would possibly be to touch pzyx directly (avoiding the following allocate)
        ALLOCATE(pab_local(ncoset(la_max_local),ncoset(lb_max_local)))
        pab_local = 0.0_dp
        DO lxa=0,la_max
        DO lxb=0,lb_max
           DO lya=0,la_max-lxa
           DO lyb=0,lb_max-lxb
              DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
              DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb

                 ! this element of pab results in 12 elements of pab_local
                 CALL prepare_dadb(pab_local,pab,lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)

              ENDDO
              ENDDO
           ENDDO
           ENDDO
        ENDDO
        ENDDO
        o1_local=0
        o2_local=0
        pab_local=pab_local * 0.5_dp
    CASE(FUNC_ADBmDAB)
        CPPrecondition(PRESENT(idir),cp_failure_level,routineP,error,failure)
        la_max_local=la_max+1
        la_min_local=MAX(la_min-1,0)
        lb_max_local=lb_max+1
        lb_min_local=MAX(lb_min-1,0)
        ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
        ! is equivalent to mapping pab with
        !    pgf_a (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) pgf_b
        ! ( pgf_a ) (ddx pgf_b) - (ddx pgf_a)( pgf_b ) =
        !          pgf_a *(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x}) -
        !                   (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_b

        ALLOCATE(pab_local(ncoset(la_max_local),ncoset(lb_max_local)))
        pab_local = 0.0_dp
        DO lxa=0,la_max
        DO lxb=0,lb_max
           DO lya=0,la_max-lxa
           DO lyb=0,lb_max-lxb
              DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
              DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
                 ! this element of pab results in 4 elements of pab_local
                 CALL prepare_adb_m_dab(pab_local,pab,idir,&
                      lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)
              END DO
              END DO
           END DO
           END DO
        END DO
        END DO
        o1_local=0
        o2_local=0
    CASE(FUNC_ARDBmDARB)
        CPPrecondition(PRESENT(idir),cp_failure_level,routineP,error,failure)
        CPPrecondition(PRESENT(ir),cp_failure_level,routineP,error,failure)
        la_max_local=la_max+1
        la_min_local=MAX(la_min-1,0)
        lb_max_local=lb_max+2
        lb_min_local=MAX(lb_min-1,0)
        ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
        ! is equivalent to mapping pab with
        ! pgf_a (r-Rb)_{ir} (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) (r-Rb)_{ir}  pgf_b
        ! ( pgf_a )(r-Rb)_{ir} (ddx pgf_b) - (ddx pgf_a) (r-Rb)_{ir} ( pgf_b ) =
        !                        pgf_a *(lbx pgf_{b-1x+1ir} - 2*zetb*pgf_{b+1x+1ir}) -
        !                       (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_{b+1ir}

        ALLOCATE(pab_local(ncoset(la_max_local),ncoset(lb_max_local)))
        pab_local = 0.0_dp
        DO lxa=0,la_max
        DO lxb=0,lb_max
           DO lya=0,la_max-lxa
           DO lyb=0,lb_max-lxb
              DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
              DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb

                 ! this element of pab results in 4 elements of pab_local
                 CALL prepare_ardb_m_darb(pab_local,pab,idir,ir,&
                      lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)
              END DO
              END DO
           END DO
           END DO
        END DO
        END DO
        o1_local=0
        o2_local=0
    CASE(FUNC_AB)
        la_max_local=la_max
        la_min_local=la_min
        lb_max_local=lb_max
        lb_min_local=lb_min
        pab_local => pab
        o1_local=o1
        o2_local=o2
    CASE DEFAULT
        CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    ng(:) = rsgrid%npts(:)
    grid => rsgrid%r(:,:,:)
    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
    lp=la_max_local+lb_max_local
    ALLOCATE(coef_xyz(((lp+1)*(lp+2)*(lp+3))/6))
    ALLOCATE(coef_xyt(((lp+1)*(lp+2))/2))
    ALLOCATE(coef_xtt(0:lp))
    ALLOCATE(alpha(0:lp,0:la_max_local,0:lb_max_local,3))

!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
    alpha(:,:,:,:)=0.0_dp
    DO iaxis=1,3
    DO lxa=0,la_max_local
    DO lxb=0,lb_max_local
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k,lxa,lxb,iaxis)=alpha(lxa-l+lxb-k,lxa,lxb,iaxis)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(iaxis)-(ra(iaxis)+rab(iaxis)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(iaxis)+rp(iaxis))
       ENDDO
    ENDDO
    ENDDO
    ENDDO

!
!   compute P_{lxp,lyp,lzp} given P_{lxa,lya,lza,lxb,lyb,lzb} and alpha(ls,lxa,lxb,1)
!   use a three step procedure
!   we don't store zeros, so counting is done using lxyz,lxy in order to have contiguous memory access in collocate_fast.F
!
    lxyz=0
    DO lzp=0,lp
    DO lyp=0,lp-lzp
    DO lxp=0,lp-lzp-lyp
        lxyz=lxyz+1
        coef_xyz(lxyz)=0.0_dp
    ENDDO
    ENDDO
    ENDDO
    DO lzb=0,lb_max_local
    DO lza=0,la_max_local
       lxy=0
       DO lyp=0,lp-lza-lzb
          DO lxp=0,lp-lza-lzb-lyp
             lxy=lxy+1
             coef_xyt(lxy)=0.0_dp
          ENDDO
          lxy=lxy+lza+lzb
       ENDDO
       DO lyb=0,lb_max_local-lzb
       DO lya=0,la_max_local-lza
          lxpm=(lb_max_local-lzb-lyb)+(la_max_local-lza-lya) 
          coef_xtt(0:lxpm)=0.0_dp
          DO lxb=MAX(lb_min_local-lzb-lyb,0),lb_max_local-lzb-lyb
          DO lxa=MAX(la_min_local-lza-lya,0),la_max_local-lza-lya
             ico=coset(lxa,lya,lza)
             jco=coset(lxb,lyb,lzb)
             p_ele=prefactor*pab_local(o1_local+ico,o2_local+jco)
             DO lxp=0,lxa+lxb
                coef_xtt(lxp)=coef_xtt(lxp)+p_ele*alpha(lxp,lxa,lxb,1)
             ENDDO
          ENDDO
          ENDDO
          lxy=0
          DO lyp=0,lya+lyb
             DO lxp=0,lp-lza-lzb-lya-lyb
               lxy=lxy+1
               coef_xyt(lxy)=coef_xyt(lxy)+alpha(lyp,lya,lyb,2)*coef_xtt(lxp)
             ENDDO
             lxy=lxy+lza+lzb+lya+lyb-lyp
          ENDDO
       ENDDO
       ENDDO
       lxyz=0
       DO lzp=0,lza+lzb
          lxy=0
          DO lyp=0,lp-lza-lzb
             DO lxp=0,lp-lza-lzb-lyp
                    lxy=lxy+1 ; lxyz=lxyz+1
                    coef_xyz(lxyz)=coef_xyz(lxyz)+alpha(lzp,lza,lzb,3)*coef_xyt(lxy)
             ENDDO
             lxy=lxy+lza+lzb ; lxyz=lxyz+lza+lzb-lzp
          ENDDO
          DO lyp=lp-lza-lzb+1,lp-lzp
             DO lxp=0,lp-lyp-lzp
                lxyz=lxyz+1
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDDO

    IF (rsgrid%orthorhombic ) THEN
       CALL collocate_ortho()
    ELSE
       CALL collocate_general_opt()
    END IF

    IF (ga_gb_function /= FUNC_AB) THEN
       DEALLOCATE(pab_local)
    ENDIF
    ! deallocation needed to pass around a pgi bug.. 
    DEALLOCATE(alpha)
    DEALLOCATE(coef_xtt)
    DEALLOCATE(coef_xyt)
    DEALLOCATE(coef_xyz)

  CONTAINS

    !
    ! this treats efficiently the orthogonal case
    !
    SUBROUTINE collocate_ortho()

!   *** properties of the grid ***

    ! notice we're in the ortho case
    dr(1) = rsgrid%dh(1,1)
    dr(2) = rsgrid%dh(2,2)
    dr(3) = rsgrid%dh(3,3)

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)
    cmax=MAXVAL(ub_cube)

!   *** position of the gaussian product
!
!   this is the actual definition of the position on the grid
!   i.e. a point rp(:) gets here grid coordinates
!   MODULO(rp(:)/dr(:),ng(:))+1
!   hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid ((1,1,1) on grid)
!

    ALLOCATE(map(-cmax:cmax,3))
    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:)    = rp(:) - REAL(cubecenter(:),dp)*dr(:)
!   *** a mapping so that the ig corresponds to the right grid point 
    DO i=1,3
      IF ( rsgrid % perd ( i ) == 1 ) THEN
        start=lb_cube(i)
        DO
         offset=MODULO(cubecenter(i)+start,ng(i))+1-start
         length=MIN(ub_cube(i),ng(i)-offset)-start
         DO ig=start,start+length
            map(ig,i) = ig+offset
         END DO
         IF (start+length.GE.ub_cube(i)) EXIT
         start=start+length+1
        END DO
      ELSE
        ! this takes partial grid + border regions into account
        offset=MODULO(cubecenter(i),ng(i))+1+rsgrid%lb(i)-rsgrid%lb_local(i)
        DO ig=lb_cube(i),ub_cube(i)
           map(ig,i) = ig+offset
        END DO
      END IF
    ENDDO
    ALLOCATE(pol_z(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_y(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_x(0:lp,-cmax:cmax))

#include "prep.f90"

    IF ( PRESENT ( lgrid ) ) THEN
       ig = lgrid%ldim * ithread_l + 1
#include "call_collocate_omp.f90"
    ELSE

#include "call_collocate.f90"

    END IF
    ! deallocation needed to pass around a pgi bug..
    DEALLOCATE(pol_z)
    DEALLOCATE(pol_y)
    DEALLOCATE(pol_x)
    DEALLOCATE(map)

    END SUBROUTINE collocate_ortho

!
!   this is a general 'optimized' routine to do the collocation
!
    SUBROUTINE collocate_general_opt()

      INTEGER :: index_min(3),index_max(3),i,j,k,i_index,j_index,k_index
      REAL(KIND=dp) :: point(3),radius2,fval,hmatgrid(3,3)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: coef_ijk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: hmatgridp
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: coef_map
      INTEGER, ALLOCATABLE, DIMENSION(:) :: grid_map
      INTEGER :: ilx,jlx,klx,ily,jly,kly,ilz,jlz,klz,lpx,lpy,lpz,lx,ly,lz,il,jl,kl,ismin,ismax,dir
      INTEGER :: offset(3)
      REAL(KIND=dp) :: res,pointk(3),pointj(3),pointi(3), &
                       rsq,a,b,c,v(3),d,exp2i,exp1i,exp0i,&
                       dk,dj,di,dkp,djp,dip,gp(3)
      ! 
      ! transform P_{lxp,lyp,lzp} into a P_{lip,ljp,lkp} such that
      ! sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-x_p)**lxp (y-y_p)**lyp (z-z_p)**lzp =
      ! sum_{lip,ljp,lkp} P_{lip,ljp,lkp} (i-i_p)**lip (j-j_p)**ljp (k-k_p)**lkp
      ! 
      ALLOCATE(coef_ijk(((lp+1)*(lp+2)*(lp+3))/6))

      ! aux mapping array to simplify life
      ALLOCATE(coef_map(0:lp,0:lp,0:lp))
      coef_map=HUGE(coef_map) 
      lxyz=0
      DO lzp=0,lp
      DO lyp=0,lp-lzp
      DO lxp=0,lp-lzp-lyp
          lxyz=lxyz+1
          coef_ijk(lxyz)=0.0_dp
          coef_map(lxp,lyp,lzp)=lxyz
      ENDDO
      ENDDO
      ENDDO

      ! cell hmat in grid points
      hmatgrid=rsgrid%dh

      ! center in grid coords
      gp=MATMUL(rsgrid%dh_inv,rp)
      cubecenter(:) = FLOOR(gp)

      ! transform using multinomials
      ALLOCATE(hmatgridp(3,3,0:lp))
      hmatgridp(:,:,0)=1.0_dp
      DO k=1,lp
         hmatgridp(:,:,k)=hmatgridp(:,:,k-1)*hmatgrid(:,:)
      ENDDO

      lpx=lp
      DO klx=0,lpx
      DO jlx=0,lpx-klx
      DO ilx=0,lpx-klx-jlx
         lx=ilx+jlx+klx
         lpy=lp-lx
         DO kly=0,lpy
         DO jly=0,lpy-kly
         DO ily=0,lpy-kly-jly
            ly=ily+jly+kly
            lpz=lp-lx-ly
            DO klz=0,lpz
            DO jlz=0,lpz-klz
            DO ilz=0,lpz-klz-jlz
               lz=ilz+jlz+klz

               il=ilx+ily+ilz
               jl=jlx+jly+jlz
               kl=klx+kly+klz
               coef_ijk(coef_map(il,jl,kl))=coef_ijk(coef_map(il,jl,kl))+ coef_xyz(coef_map(lx,ly,lz))* &
                                            hmatgridp(1,1,ilx) * hmatgridp(1,2,jlx) * hmatgridp(1,3,klx) * &
                                            hmatgridp(2,1,ily) * hmatgridp(2,2,jly) * hmatgridp(2,3,kly) * &
                                            hmatgridp(3,1,ilz) * hmatgridp(3,2,jlz) * hmatgridp(3,3,klz) * &
                                            fac(lx)*fac(ly)*fac(lz)/ &
                        (fac(ilx)*fac(ily)*fac(ilz)*fac(jlx)*fac(jly)*fac(jlz)*fac(klx)*fac(kly)*fac(klz))
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL return_cube_nonortho(cube_info,radius,index_min,index_max,rp)

      DO dir=1,3
         offset(dir)=MODULO(cubecenter(dir),ng(dir))+rsgrid%lb(dir)
         offset(dir)=offset(dir)-rsgrid%lb_local(dir)+1
      ENDDO

      ALLOCATE(grid_map(index_min(1):index_max(1)))
      DO i=index_min(1),index_max(1)
        grid_map(i)=MODULO(i,ng(1))+1
        IF (rsgrid % perd ( 1 )==1) THEN
           grid_map(i)=MODULO(i,ng(1))+1
        ELSE
           grid_map(i)=i-cubecenter(1)+offset(1)
        ENDIF
      ENDDO


      ! go over the grid, but cycle if the point is not within the radius
      DO k=index_min(3),index_max(3)
        dk=k-gp(3)
        pointk=hmatgrid(:,3)*dk

        IF (rsgrid % perd ( 3 )==1) THEN
           k_index=MODULO(k,ng(3))+1
        ELSE
           k_index=k-cubecenter(3)+offset(3)
        ENDIF

        coef_xyt=0.0_dp
        lxyz = 0
        dkp=1.0_dp
        DO kl=0,lp
           lxy=0
           DO jl=0,lp-kl
              DO il=0,lp-kl-jl
                 lxyz=lxyz+1 ; lxy=lxy+1
                 coef_xyt(lxy)=coef_xyt(lxy)+coef_ijk(lxyz)*dkp
              ENDDO
              lxy=lxy+kl
           ENDDO
           dkp=dkp*dk
        ENDDO

        DO j=index_min(2),index_max(2)
          dj=j-gp(2)
          pointj=pointk+hmatgrid(:,2)*dj
          IF (rsgrid % perd ( 2 )==1) THEN
             j_index=MODULO(j,ng(2))+1
          ELSE
             j_index=j-cubecenter(2)+offset(2)
          ENDIF

          coef_xtt=0.0_dp
          lxy=0
          djp=1.0_dp
          DO jl=0,lp
            DO il=0,lp-jl
               lxy=lxy+1
               coef_xtt(il)=coef_xtt(il)+coef_xyt(lxy)*djp
            ENDDO
            djp=djp*dj
          ENDDO

          ! find bounds for the inner loop
          ! based on a quadratic equation in i
          ! a*i**2+b*i+c=radius**2
          v=pointj-gp(1)*hmatgrid(:,1)
          a=DOT_PRODUCT(hmatgrid(:,1),hmatgrid(:,1))
          b=2*DOT_PRODUCT(v,hmatgrid(:,1))
          c=DOT_PRODUCT(v,v)
          d=b*b-4*a*(c-radius**2)

          IF (d<0) THEN
              CYCLE
          ELSE
              d=sqrt(d)
              ismin=CEILING((-b-d)/(2*a))
              ismax=FLOOR((-b+d)/(2*a))
          ENDIF
          ! prepare for computing -zetp*rsq
          a=-zetp*a
          b=-zetp*b
          c=-zetp*c
          i=ismin-1

          ! the recursion relation might have to be done
          ! from the center of the gaussian (in both directions)
          ! instead as the current implementation from an edge
          exp2i=exp((a*i+b)*i+c)
          exp1i=exp(2*a*i+a+b)
          exp0i=exp(2*a)

          DO i=ismin,ismax
             di=i-gp(1)

             ! polynomial terms
             res=0.0_dp
             dip=1.0_dp
             DO il=0,lp
                res=res+coef_xtt(il)*dip
                dip=dip*di
             ENDDO

             ! the exponential recursion
             exp2i=exp2i*exp1i
             exp1i=exp1i*exp0i
             res=res*exp2i

             i_index=grid_map(i)
             IF ( PRESENT ( lgrid ) ) THEN
                ig = lgrid%ldim * ithread_l + (k_index-1) * ng(2) * ng(1) + (j_index-1) * ng(1) + (i_index-1) + 1
                lgrid%r(ig)=lgrid%r(ig)+res
             ELSE
                grid(i_index,j_index,k_index)=grid(i_index,j_index,k_index)+res
             ENDIF
          ENDDO
        ENDDO
      ENDDO
      !t2=nanotime_ia32()
      !write(6,*) t2-t1
      ! deallocation needed to pass around a pgi bug..
      DEALLOCATE(coef_ijk)
      DEALLOCATE(coef_map)
      DEALLOCATE(hmatgridp)
      DEALLOCATE(grid_map)

    END SUBROUTINE collocate_general_opt

!
!   this is a general 'reference' routine to do the collocation
!
    SUBROUTINE collocate_general()
 
      INTEGER :: index_min(3),index_max(3),i,j,k,ipoint(3)
      REAL(KIND=dp) :: point(3),radius2,fval

      ! still hard-wired (see MODULO)
      CPPostcondition(ALL(rsgrid % perd==1),cp_failure_level,routineP,error,failure)

      CALL return_cube_nonortho(cube_info,radius,index_min,index_max,rp)

      ! go over the grid, but cycle if the point is not within the radius
      DO k=index_min(3),index_max(3)
      DO j=index_min(2),index_max(2)
      DO i=index_min(1),index_max(1)
         ! point in real space
         point=MATMUL(cell%hmat,REAL((/i,j,k/),KIND=dp)/ng)
         ! skip if outside of the sphere
         if (SUM((point-rp)**2)>radius**2) CYCLE
         ! point on the grid (including pbc)
         ipoint=MODULO((/i,j,k/),ng)+1
         ! add to grid
         IF ( PRESENT ( lgrid ) ) THEN
            ig = lgrid%ldim * ithread_l + ipoint(3) * ng(2) * ng(1) + ipoint(2) * ng(1) + ipoint(1) + 1
            lgrid%r(ig)=lgrid%r(ig)+primitive_value(point)
         ELSE
            grid(ipoint(1),ipoint(2),ipoint(3))=grid(ipoint(1),ipoint(2),ipoint(3))+primitive_value(point)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
       
    END SUBROUTINE collocate_general

    FUNCTION primitive_value(point) RESULT(res)
       REAL(KIND=dp) :: point(3),dra(3),drb(3), res,myexp,drap(3),drbp(3)
       res=0.0_dp

        myexp=exp(-zetp*SUM((point-rp)**2))*prefactor
        dra=point-ra
        drb=point-rb
        drap(1)=1.0_dp
        DO lxa=0,la_max_local
        drbp(1)=1.0_dp
        DO lxb=0,lb_max_local
           drap(2)=1.0_dp
           DO lya=0,la_max_local-lxa
           drbp(2)=1.0_dp
           DO lyb=0,lb_max_local-lxb
              drap(3)=1.0_dp
              DO lza=1,MAX(la_min_local-lxa-lya,0) 
                 drap(3)=drap(3)*dra(3) 
              ENDDO
              DO lza=MAX(la_min_local-lxa-lya,0),la_max_local-lxa-lya
              drbp(3)=1.0_dp
              DO lzb=1,MAX(lb_min_local-lxb-lyb,0) 
                 drbp(3)=drbp(3)*drb(3) 
              ENDDO
              DO lzb=MAX(lb_min_local-lxb-lyb,0),lb_max_local-lxb-lyb
                ico=coset(lxa,lya,lza)
                jco=coset(lxb,lyb,lzb)
                res=res+pab_local(ico+o1_local,jco+o2_local)*myexp*PRODUCT(drap)*PRODUCT(drbp)
                drbp(3)=drbp(3)*drb(3)
              ENDDO
              drap(3)=drap(3)*dra(3)
              ENDDO
           drbp(2)=drbp(2)*drb(2)
           ENDDO
           drap(2)=drap(2)*dra(2)
           ENDDO
        drbp(1)=drbp(1)*drb(1)
        ENDDO
        drap(1)=drap(1)*dra(1)
        ENDDO

    END FUNCTION primitive_value

  END SUBROUTINE collocate_pgf_product_rspace

!!****f* qs_collocate_density/collocate_pgf_product_gspace *
!!
!!   NAME
!!     collocate_pgf_product_gspace
!!
!!   FUNCTION
!!     low level collcation of primitive gaussian functions in g-space
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE collocate_pgf_product_gspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2,scale,pab,na,nb,&
                                          eps_rho_gspace,gsq_max,pw)

    ! NOTE: this routine is much slower than collocate_pgf_product_rspace

    INTEGER, INTENT(IN)                      :: la_max
    REAL(dp), INTENT(IN)                     :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(dp), INTENT(IN)                     :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(dp), DIMENSION(3), INTENT(IN)       :: ra, rab
    REAL(dp), INTENT(IN)                     :: rab2, scale
    REAL(dp), DIMENSION(:, :), POINTER       :: pab
    INTEGER, INTENT(IN)                      :: na, nb
    REAL(dp), INTENT(IN)                     :: eps_rho_gspace, gsq_max
    TYPE(pw_type), POINTER                   :: pw

    CHARACTER(LEN=*), PARAMETER :: routineN = 'collocate_pgf_product_gspace', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp), DIMENSION(3)                :: phasefactor
    COMPLEX(dp), DIMENSION(:), POINTER       :: rag, rbg
    COMPLEX(dp), DIMENSION(:, :, :, :), &
      POINTER                                :: cubeaxis
    INTEGER :: ax, ay, az, bx, by, bz, handle, i, ico, ig, ig2, istat, jco, &
      jg, kg, la, lb, lb_cube_min, lb_grid, ub_cube_max, ub_grid
    INTEGER, DIMENSION(3)                    :: lb_cube, ub_cube
    REAL(dp)                                 :: f, fa, fb, pij, prefactor, &
                                                rzetp, twozetp, zetp
    REAL(dp), DIMENSION(3)                   :: dg, expfactor, fap, fbp, rap, &
                                                rbp, rp
    REAL(dp), DIMENSION(:), POINTER          :: g

    CALL timeset(routineN,"I","",handle)

    dg(:) = twopi/(pw%pw_grid%npts(:)*pw%pw_grid%dr(:))

    zetp = zeta + zetb
    rzetp = 1.0_dp/zetp
    f = zetb*rzetp
    rap(:) = f*rab(:)
    rbp(:) = rap(:) - rab(:)
    rp(:) = ra(:) + rap(:)
    twozetp = 2.0_dp*zetp
    fap(:) = twozetp*rap(:)
    fbp(:) = twozetp*rbp(:)

    prefactor = scale*SQRT((pi*rzetp)**3)*EXP(-zeta*f*rab2)
    phasefactor(:) = EXP(CMPLX(0.0_dp,-rp(:)*dg(:),KIND=dp))
    expfactor(:) = EXP(-0.25*rzetp*dg(:)*dg(:))

    lb_cube(:) = pw%pw_grid%bounds(1,:)
    ub_cube(:) = pw%pw_grid%bounds(2,:)

    lb_cube_min = MINVAL(lb_cube(:))
    ub_cube_max = MAXVAL(ub_cube(:))

    NULLIFY (cubeaxis,g,rag,rbg)

    CALL reallocate(cubeaxis,lb_cube_min,ub_cube_max,1,3,0,la_max,0,lb_max)
    CALL reallocate(g,lb_cube_min,ub_cube_max)
    CALL reallocate(rag,lb_cube_min,ub_cube_max)
    CALL reallocate(rbg,lb_cube_min,ub_cube_max)

    lb_grid = LBOUND(pw%cc,1)
    ub_grid = UBOUND(pw%cc,1)

    DO i=1,3

      DO ig=lb_cube(i),ub_cube(i)
        ig2 = ig*ig
        cubeaxis(ig,i,0,0) = expfactor(i)**ig2*phasefactor(i)**ig
      END DO

      IF (la_max > 0) THEN
        DO ig=lb_cube(i),ub_cube(i)
          g(ig) = REAL(ig,dp)*dg(i)
          rag(ig) = CMPLX(fap(i),-g(ig),KIND=dp)
          cubeaxis(ig,i,1,0) = rag(ig)*cubeaxis(ig,i,0,0)
        END DO
        DO la=2,la_max
          fa = REAL(la-1,dp)*twozetp
          DO ig=lb_cube(i),ub_cube(i)
            cubeaxis(ig,i,la,0) = rag(ig)*cubeaxis(ig,i,la-1,0) +&
                                  fa*cubeaxis(ig,i,la-2,0)
          END DO
        END DO
        IF (lb_max > 0) THEN
          fa = twozetp
          DO ig=lb_cube(i),ub_cube(i)
            rbg(ig) = CMPLX(fbp(i),-g(ig),KIND=dp)
            cubeaxis(ig,i,0,1) = rbg(ig)*cubeaxis(ig,i,0,0)
            cubeaxis(ig,i,1,1) = rbg(ig)*cubeaxis(ig,i,1,0) +&
                                 fa*cubeaxis(ig,i,0,0)
          END DO
          DO lb=2,lb_max
            fb = REAL(lb-1,dp)*twozetp
            DO ig=lb_cube(i),ub_cube(i)
              cubeaxis(ig,i,0,lb) = rbg(ig)*cubeaxis(ig,i,0,lb-1) +&
                                    fb*cubeaxis(ig,i,0,lb-2)
              cubeaxis(ig,i,1,lb) = rbg(ig)*cubeaxis(ig,i,1,lb-1) +&
                                    fb*cubeaxis(ig,i,1,lb-2) +&
                                    fa*cubeaxis(ig,i,0,lb-1)
            END DO
          END DO
          DO la=2,la_max
            fa = REAL(la,dp)*twozetp
            DO ig=lb_cube(i),ub_cube(i)
              cubeaxis(ig,i,la,1) = rbg(ig)*cubeaxis(ig,i,la,0) +&
                                    fa*cubeaxis(ig,i,la-1,0)
            END DO
            DO lb=2,lb_max
              fb = REAL(lb-1,dp)*twozetp
              DO ig=lb_cube(i),ub_cube(i)
                cubeaxis(ig,i,la,lb) = rbg(ig)*cubeaxis(ig,i,la,lb-1) +&
                                       fb*cubeaxis(ig,i,la,lb-2) +&
                                       fa*cubeaxis(ig,i,la-1,lb-1)
              END DO
            END DO
          END DO
        END IF
      ELSE
        IF (lb_max > 0) THEN
          DO ig=lb_cube(i),ub_cube(i)
            g(ig) = REAL(ig,dp)*dg(i)
            rbg(ig) = CMPLX(fbp(i),-g(ig),KIND=dp)
            cubeaxis(ig,i,0,1) = rbg(ig)*cubeaxis(ig,i,0,0)
          END DO
          DO lb=2,lb_max
            fb = REAL(lb-1,dp)*twozetp
            DO ig=lb_cube(i),ub_cube(i)
              cubeaxis(ig,i,0,lb) = rbg(ig)*cubeaxis(ig,i,0,lb-1) +&
                                    fb*cubeaxis(ig,i,0,lb-2)
            END DO
          END DO
        END IF
      END IF

    END DO

    DO la=0,la_max
      DO lb=0,lb_max
        IF (la + lb == 0) CYCLE
        fa = (1.0_dp/twozetp)**(la + lb)
        DO i=1,3
          DO ig=lb_cube(i),ub_cube(i)
            cubeaxis(ig,i,la,lb) = fa*cubeaxis(ig,i,la,lb)
          END DO
        END DO
      END DO
    END DO

    ! Add the current primitive Gaussian function product to grid

    DO ico=ncoset(la_min-1)+1,ncoset(la_max)

      ax = indco(1,ico)
      ay = indco(2,ico)
      az = indco(3,ico)

      DO jco=ncoset(lb_min-1)+1,ncoset(lb_max)

        pij = prefactor*pab(na+ico,nb+jco)

        IF (ABS(pij) < eps_rho_gspace) CYCLE

        bx = indco(1,jco)
        by = indco(2,jco)
        bz = indco(3,jco)

        DO i=lb_grid,ub_grid
          IF (pw%pw_grid%gsq(i) > gsq_max) CYCLE
          ig = pw%pw_grid%g_hat(1,i)
          jg = pw%pw_grid%g_hat(2,i)
          kg = pw%pw_grid%g_hat(3,i)
          pw%cc(i) = pw%cc(i) + pij*cubeaxis(ig,1,ax,bx)*&
                                    cubeaxis(jg,2,ay,by)*&
                                    cubeaxis(kg,3,az,bz)
        END DO

      END DO

    END DO

    DEALLOCATE (cubeaxis,g,rag,rbg,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,"cubeaxis,g,rag,rbg")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE collocate_pgf_product_gspace

!!****f* qs_collocate_density/density_rs2pw *
!!
!!   NAME
!!     density_rs2pw
!!
!!   FUNCTION
!!     given partial densities on the realspace multigrids, 
!!     computes the full density on the plane wave grids, both in real and 
!!     gspace
!!
!!   NOTES
!!     should contain the all communication in the collocation of the density
!!     in the case of replicated grids
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE density_rs2pw(pw_env,rs_rho,rho,rho_gspace,interp_section,error)

    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    TYPE(pw_p_type), INTENT(INOUT)           :: rho, rho_gspace
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: interp_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'density_rs2pw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, igrid_level, &
                                                interp_kind, stat
    LOGICAL                                  :: failure
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: mgrid_gspace, mgrid_rspace, &
                                                mgrid_temp_rspace
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(gridlevel_info,mgrid_gspace, mgrid_rspace,mgrid_temp_rspace,rs_pools,pw_pools)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineN,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, pw_pools=pw_pools, error=error)

    gridlevel_info=>pw_env%gridlevel_info
    IF(PRESENT(interp_section)) THEN
      CALL section_vals_val_get(interp_section,"KIND",i_val=interp_kind,error=error)
    ELSE
      interp_kind = pw_interp
    END IF

    CALL pw_pools_create_pws(pw_pools,mgrid_rspace,&
                              use_data = REALDATA3D,&
                              in_space = REALSPACE, error=error)

    CALL pw_pools_create_pws(pw_pools,mgrid_gspace,&
                              use_data = COMPLEXDATA1D,&
                              in_space = RECIPROCALSPACE, error=error)

   IF (gridlevel_info%ngrid_levels==1) THEN
       CALL rs_pw_transfer(rs_rho(1)%rs_grid,rho%pw,rs2pw)
       CALL rs_pools_give_back_rs_vect(rs_pools, rs_rho, error=error)
       CALL pw_transfer(rho%pw,rho_gspace%pw)
       IF (rho%pw%pw_grid%spherical) THEN ! rho_gspace = rho
          CALL pw_transfer(rho_gspace%pw,rho%pw)
       ENDIF
    ELSE
       DO igrid_level=1,gridlevel_info%ngrid_levels
          CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid,&
               mgrid_rspace(igrid_level)%pw,rs2pw)
       ENDDO
       CALL rs_pools_give_back_rs_vect(rs_pools, rs_rho, error=error)

       ! we want both rho and rho_gspace, the latter for Hartree and co-workers.
       SELECT CASE(interp_kind)
       CASE(pw_interp)
          CALL pw_zero(rho_gspace%pw)
          DO igrid_level=1,gridlevel_info%ngrid_levels
             CALL pw_transfer(mgrid_rspace(igrid_level)%pw,&
                  mgrid_gspace(igrid_level)%pw)
             CALL pw_axpy(mgrid_gspace(igrid_level)%pw,rho_gspace%pw)
          END DO
          CALL pw_transfer(rho_gspace%pw,rho%pw)
       CASE(spline3_pbc_interp)
          DO igrid_level=gridlevel_info%ngrid_levels,2,-1
             CALL pw_prolongate_s3(mgrid_rspace(igrid_level)%pw,&
                  mgrid_rspace(igrid_level-1)%pw,pw_pools(igrid_level)%pool,&
                  interp_section,error=error)
          END DO
          CALL pw_copy(mgrid_rspace(1)%pw,rho%pw)
          CALL pw_transfer(rho%pw,rho_gspace%pw)
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

END MODULE qs_collocate_density
