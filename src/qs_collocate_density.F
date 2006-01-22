!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003, 2004 CP2K developers group                           !
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
  USE coefficient_types,               ONLY: coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
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
                                             return_cube
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE kinds,                           ONLY: dp, mp,&
                                             dp_size
  USE l_utils,                         ONLY: l_info_type,&
                                             return_l_info
  USE mathconstants,                   ONLY: pi,&
                                             twopi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: coset,&
                                             indco,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_give_back_coeffs,&
                                             pw_pools_init_coeffs
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_copy,&
                                             pw_prolongate_l,&
                                             pw_sumup,&
                                             pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_interactions,                 ONLY: exp_radius_very_extended
  USE qs_neighbor_list_types,          ONLY: first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs_get_loop_vars,&
                                             rs_get_my_tasks,&
                                             rs_grid_zero,&
                                             rs_pw_transfer,&
                                             rs2pw
  USE cell_types,                 ONLY: cell_type,&
                                             pbc
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             deallocate_matrix,&
                                             get_block_node,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
  USE pw_spline_utils, only:pw_prolongate_s3
  USE input_constants, ONLY: pw_interp, linear_interp, spline3_pbc_interp
  USE input_section_types, ONLY: section_vals_type, section_vals_get_subs_vals,&
       section_vals_val_get
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN='qs_collocate_density'
! *** Public subroutines ***

  PUBLIC :: calculate_rho_core,&
            calculate_rho_elec,&
            calculate_wavefunction,&
            collocate_pgf_product_gspace,&
            collocate_pgf_product_rspace,&
            collocate_atomic_charge_density

! *** Public functions ***

  PUBLIC :: calculate_total_rho,calculate_total_abs_rho

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
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    INTEGER                                  :: handle, &
                                                ikind, &
                                                nkind, &
                                                nseta, &
                                                iset, &
                                                ipgf, &
                                                iatom, &
                                                natom, &
                                                ierr, &
                                                ithread, &
                                                sgfa, &
                                                igrid_level
    INTEGER, DIMENSION(:), POINTER           :: la_max, &
                                                la_min, &
                                                atom_list, &
                                                npgfa, &
                                                nsgfa
    INTEGER, DIMENSION(:,:), POINTER         :: first_sgfa
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: zeta, &
                                                sphi_a, &
                                                pab, &
                                                work
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    !TYPE(realspace_grid_type), POINTER       :: rs_rho
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    !TYPE(cp_rs_pool_type), POINTER           :: auxbas_rs_pool
    TYPE(pw_env_type), POINTER               :: pw_env
    !TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    !TYPE(coeff_type)                         :: rhoc_r
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    !TYPE(cube_info_type), POINTER            :: cube_info
    TYPE(l_info_type)                        :: l_info
    REAL(KIND=dp)                            :: eps_rho_rspace
    TYPE(dft_control_type), POINTER          :: dft_control
    LOGICAL                                  :: map_consistent, &
                                                failure
    TYPE(qs_rho_type), POINTER               :: rho_struct
    INTEGER                                  :: i, &
                                                maxco, &
                                                ncoa, &
                                                na1, &
                                                unit_nr
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_gspace, &
                                                mgrid_rspace
    TYPE(cp_logger_type), POINTER            :: logger

    CHARACTER(len=*), PARAMETER :: routineN = 'collocate_atomic_charge_density', &
      routineP = moduleN//':'//routineN
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
    l_info = pw_env%l_info
    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    ithread = 0
    gridlevel_info=>pw_env%gridlevel_info

    ! *** set up the pw multi-grids *** !
    CPPrecondition(ASSOCIATED(pw_env), cp_failure_level, routineN, error, failure)
    CALL pw_env_get(pw_env=pw_env, rs_pools=rs_pools, pw_pools=pw_pools, error=error)

    ALLOCATE(mgrid_rspace(SIZE(pw_pools)), stat=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineN, error, failure)
    CALL pw_pools_init_coeffs(pools=pw_pools, coeffs=mgrid_rspace, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)

    ALLOCATE(mgrid_gspace(SIZE(pw_pools)), stat=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineN, error, failure)
    CALL pw_pools_init_coeffs(pools=pw_pools, coeffs=mgrid_gspace, &
                              use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                              error=error)

    ! *** set up the rs multi-grids *** !
    CALL rs_pools_create_rs_vect(rs_pools,rs_rho, force_env_section=qs_env%input,&
         error=error)
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
                                              la_min=la_min(iset), lb_max=0, zetb=0.0_dp, lb_min=0, &
                                              ra=ra, rab=(/0.0_dp,0.0_dp,0.0_dp/), rab2=0.0_dp, &
                                              scale=1.0_dp, pab=pab, o1=na1, o2=0, &
                                              rsgrid=rs_rho(igrid_level)%rs_grid, &
                                              cube_info=cube_info(igrid_level), l_info=l_info, &
                                              eps_rho_rspace=eps_rho_rspace, ithread=ithread, &
                                              compute_tau=.FALSE., map_consistent=map_consistent)
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
      CALL coeff_transform_space(qs_env%rho%rho_r(1), qs_env%rho%rho_g(1))
      IF (qs_env%rho%rho_r(1)%pw%pw_grid%spherical) THEN
        CALL coeff_transform_space(qs_env%rho%rho_g(1), qs_env%rho%rho_r(1))
      END IF
    ELSE
      DO igrid_level = 1,gridlevel_info%ngrid_levels
        CALL rs_pw_transfer(rs=rs_rho(igrid_level)%rs_grid, &
                            pw=mgrid_rspace(igrid_level)%pw, dir=rs2pw)
      END DO
      CALL rs_pools_give_back_rs_vect(pools=rs_pools, elements=rs_rho, error=error)
      
      !CALL coeff_zero(rho_gspace)
      CALL coeff_zero(qs_env%rho%rho_g(1))
      DO igrid_level=1, gridlevel_info%ngrid_levels
        CALL coeff_transform_space(mgrid_rspace(igrid_level), &
                                   mgrid_gspace(igrid_level))
        CALL coeff_sumup(mgrid_gspace(igrid_level), qs_env%rho%rho_g(1))
      END DO
      CALL coeff_transform_space(qs_env%rho%rho_g(1), qs_env%rho%rho_r(1))
    END IF

  !  CALL rs_pw_transfer(rs=rs_rho, pw=qs_env%rho%rho_r(1)%pw, dir=rs2pw)

  !  CALL rs_pool_give_back_rs(pool=auxbas_rs_pool, element=rs_rho, error=error)

    !CALL coeff_transform_space(rhoc_r, qs_env%rho%rho_g(1))
  !  CALL coeff_transform_space(qs_env%rho%rho_r(1), qs_env%rho%rho_g(1))

  !  IF (qs_env%rho%rho_r(1)%pw%pw_grid%spherical) THEN
  !    CALL coeff_transform_space(qs_env%rho%rho_g(1), qs_env%rho%rho_r(1))
  !  END IF

    total_rho = calculate_total_rho(qs_env%rho%rho_r(1))
    qs_env%rho%tot_rho_r = total_rho
    IF (logger%para_env%source==logger%para_env%mepos) THEN
      unit_nr=cp_logger_get_default_unit_nr(logger,local=.FALSE.)
      WRITE (unit_nr,*) "Total Rho =", total_rho
    END IF

    ! *** give back the multi-grids *** !
    CALL pw_pools_give_back_coeffs(pw_pools, mgrid_gspace, error=error)
    DEALLOCATE(mgrid_gspace, stat=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineN, error, failure)
    CALL pw_pools_give_back_coeffs(pw_pools, mgrid_rspace, error=error)
    DEALLOCATE(mgrid_rspace, stat=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineN, error, failure)

    qs_env%rho%rho_r_valid = .TRUE.
    qs_env%rho%rho_g_valid = .TRUE.

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE collocate_atomic_charge_density

! *****************************************************************************

  SUBROUTINE calculate_rho_core(rho_core,total_rho,qs_env,error)

    TYPE(coeff_type), INTENT(INOUT)          :: rho_core
    REAL(KIND=dp), INTENT(OUT)               :: total_rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_core', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, bo(2), dir, handle, &
                                                iatom, ierr, ikind, ithread, &
                                                j, natom, ncurr, npme, &
                                                omp_get_thread_num
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    INTEGER, DIMENSION(:, :), POINTER        :: tasks
    REAL(KIND=dp)                            :: alpha, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: rhoc_r
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_type), POINTER           :: auxbas_rs_pool
    TYPE(cube_info_type)                     :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(l_info_type)                        :: l_info
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_type), POINTER       :: rs_rho

    SAVE pab                ! Temporary addition for Intel 8.0 OpenMP to work
!   ---------------------------------------------------------------------------

    CALL timeset("calculate_rho_core","I"," ",handle)
    NULLIFY(atomic_kind,cell,dft_control,pab,atomic_kind_set,particle_set,&
         atom_list,pw_env,rs_rho,auxbas_rs_pool,auxbas_pw_pool,tasks)
    ALLOCATE(pab(1,1),STAT=ierr)
    IF ( ierr /= 0 ) CALL stop_memory ( "calculate_rho_core", "pab", 1 )

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    para_env=para_env,pw_env=pw_env)
    CALL pw_env_get(pw_env,auxbas_rs_pool=auxbas_rs_pool,&
         auxbas_pw_pool=auxbas_pw_pool,error=error)
    cube_info=pw_env%cube_info(1)
    l_info=pw_env%l_info
    ! be careful in parallel nsmax is choosen with multigrid in mind!
    CALL rs_pool_create_rs(auxbas_rs_pool,rs_rho,force_env_section=qs_env%input,&
         error=error)
    CALL rs_grid_zero(rs_rho)

    ncurr = -1

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           atom_list=atom_list,&
                           alpha_core_charge=alpha,&
                           ccore_charge=pab(1,1))

      IF (alpha == 0.0_dp) CYCLE

      bo = get_limit ( natom, para_env%num_pe, para_env%mepos )

      npme = bo(2) - bo(1) + 1
      IF ( ncurr < npme ) THEN
        CALL reallocate ( tasks, 1, 2, 1, npme )
        ncurr = npme
      ELSE
        tasks = 0
      END IF

      dir = rs_rho%direction

!$OMP parallel do &
!$OMP default(none) &
!$OMP private(iatom,j,atom_a,ra) &
!$OMP shared(bo,tasks,dir,rs_rho,atom_list,particle_set,cell)
      DO iatom=bo(1),bo(2)

        j = iatom - bo(1) + 1
        tasks ( 1, j ) = iatom
        IF ( dir /= 0) THEN
          atom_a = atom_list(iatom)
          ra(:) = pbc(particle_set(atom_a)%r,cell)
          tasks ( 2, j ) = FLOOR(ra(dir)/rs_rho%dr(dir))
          tasks ( 2, j ) = MODULO ( tasks ( 2, j ), rs_rho%npts(dir) )
          tasks ( 2, j ) = tasks ( 2, j ) + rs_rho%lb(dir)
        END IF

      END DO

      CALL rs_get_my_tasks ( rs_rho, tasks, npme )

      ithread=0
!$OMP parallel do &
!$OMP default(none) &
!$OMP private(j,iatom,atom_a,ra,ithread) &
!$OMP shared(npme,tasks,atom_list,particle_set,cell,pab,rs_rho) &
!$OMP shared(cube_info,l_info,eps_rho_rspace,alpha)
      DO j=1,npme

        iatom = tasks(1,j)
        atom_a = atom_list(iatom)
        ra(:) = pbc(particle_set(atom_a)%r,cell)

!$      ithread=omp_get_thread_num()
        CALL collocate_pgf_product_rspace(0,alpha,0,0,0.0_dp,0,ra,&
             (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,-1.0_dp,pab,0,0,rs_rho,&
             cube_info,l_info,eps_rho_rspace,ithread=ithread)

      END DO
    
    END DO

    IF (ASSOCIATED(tasks)) THEN
      DEALLOCATE (tasks,STAT=ierr)
      IF ( ierr /= 0 ) CALL stop_memory ( "calculate_rho_core", "tasks" )
    END IF
    DEALLOCATE ( pab, STAT=ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "calculate_rho_core", "pab" )

    CALL pw_pool_init_coeff(pool=auxbas_pw_pool, coeff=rhoc_r, &
         use_data=REALDATA3D,in_space=REALSPACE, error=error)

    CALL rs_pw_transfer(rs_rho,rhoc_r%pw,rs2pw)
    CALL rs_pool_give_back_rs(auxbas_rs_pool,rs_rho,error=error)

    total_rho = calculate_total_rho(rhoc_r)

    CALL coeff_transform_space(rhoc_r,rho_core)

    CALL pw_pool_give_back_coeff(auxbas_pw_pool, rhoc_r, error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_rho_core

! *****************************************************************************
! both rho,rho_gspace are the new rho,
! *****************************************************************************

  SUBROUTINE calculate_rho_elec(matrix_p,rho,rho_gspace, total_rho,&
                                qs_env, soft_valid,  compute_tau, error)

    TYPE(real_matrix_type), POINTER          :: matrix_p
    TYPE(coeff_type), INTENT(INOUT)          :: rho, rho_gspace
    REAL(KIND=dp), INTENT(OUT)               :: total_rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: soft_valid, compute_tau
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE calculate_rho_elec (MODULE qs_collocate_density)"
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 2000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: ab, bcol, brow, curr_tasks, dir, first_pgfb, first_setb, &
      handle, i, iatom, igrid_level, igridlevel, ijatoms, ijsets, ikind, &
      ilist, inode, ipgf, iset, istat, itask, ithread, j, jatom, jkind, jpgf, &
      jset, k, maxco, maxsgf, maxsgf_set, n, na1, na2,  &
      natom_pairs, nb1, nb2, ncoa, ncob, nkind, nlist, nnode, npme, nseta, &
      nsetb, nthread, omp_get_max_threads, omp_get_thread_num, &
      sgfa, sgfb, stat, tp
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb, ntasks, ub
    INTEGER, DIMENSION(:, :), POINTER        :: asets, atasks, first_sgfa, &
                                                first_sgfb, ival, latom, &
                                                tasks_local
    INTEGER, DIMENSION(:, :, :), POINTER     :: tasks
    LOGICAL                                  :: failure, map_consistent, &
                                                my_compute_tau, my_soft
    REAL(KIND=dp)                            :: dab, eps_rho_rspace, &
                                                kind_radius_b, rab2, scale, &
                                                zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb, rp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dab_local, p_block, pab, &
                                                pblock, rpgfa, rpgfb, sphi_a, &
                                                sphi_b, work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dist_ab, pabt, workt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_gspace, mgrid_rspace, &
                                                mgrid_temp_rspace
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(l_info_type)                        :: l_info
    TYPE(lgrid_type)                         :: lgrid
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(real_matrix_type), POINTER          :: deltap
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    LOGICAL                                  :: distributed_rs_grids
    TYPE(section_vals_type), POINTER :: input, interp_section
    integer :: interp_kind
!   ---------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(atomic_kind,cell,dft_control,orb_basis_set,sab_orb_neighbor_list,&
         sab_orb_neighbor_node,deltap,atomic_kind_set,sab_orb,particle_set,&
         rs_rho,pw_env,rs_pools,para_env,pblock,dist_ab,dab_local,&
         set_radius_a,set_radius_b,la_max,la_min,&
         lb_max,lb_min,npgfa,npgfb,nsgfa,nsgfb,p_block,&
         rpgfa,rpgfb,sphi_a,sphi_b,zeta,zetb,first_sgfa,first_sgfb,&
         tasks,tasks_local,ival,latom,ntasks,asets,atasks,pabt,workt,mgrid_rspace,&
         mgrid_gspace,pw_pools)
    NULLIFY(mgrid_temp_rspace)
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
    ELSE
       CALL timeset("calculate_rho_elec","I"," ",handle)
    ENDIF

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb,&
                    para_env=para_env,&
                    input=input,&
                    pw_env=pw_env)

    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info
    l_info=pw_env%l_info

    interp_section => section_vals_get_subs_vals(input,"DFT%MGRID%INTERPOLATOR",&
         error=error)
    CALL section_vals_val_get(interp_section,"KIND",i_val=interp_kind,error=error)

    ! *** set up the pw multi-grids 
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routine,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, pw_pools=pw_pools, error=error)

    ALLOCATE(mgrid_rspace(SIZE(pw_pools)) ,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
    CALL pw_pools_init_coeffs(pw_pools,mgrid_rspace,&
                              use_data = REALDATA3D,&
                              in_space = REALSPACE, error=error)
 
    IF (interp_kind==linear_interp) THEN
        ALLOCATE(mgrid_temp_rspace(SIZE(pw_pools)) ,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
        CALL pw_pools_init_coeffs(pw_pools,mgrid_temp_rspace,&
                                  use_data = REALDATA3D,&
                                  in_space = REALSPACE, error=error)
    ELSE
        ALLOCATE(mgrid_gspace(SIZE(pw_pools)) ,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
        CALL pw_pools_init_coeffs(pw_pools,mgrid_gspace,&
                                  use_data = COMPLEXDATA1D,&
                                  in_space = RECIPROCALSPACE, error=error)
    ENDIF

    ! *** set up the rs multi-grids
    distributed_rs_grids=.FALSE.
    CALL rs_pools_create_rs_vect(rs_pools, rs_rho, force_env_section=input,&
         error=error)
    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
       IF ( rs_rho(igrid_level)%rs_grid%direction /= 0 ) THEN
          distributed_rs_grids=.TRUE.
       ENDIF
    END DO

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    nthread = 1
!$  nthread = omp_get_max_threads()

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf=maxsgf,&
                             maxsgf_set=maxsgf_set)

    IF ( nthread > 1 ) THEN
      n=0
      DO igrid_level = 1,gridlevel_info%ngrid_levels
        n = MAX(n,rs_rho(igrid_level)%rs_grid%ngpts_local)
      END DO
      n = n*nthread
      CALL reallocate(lgrid%r,1,n)
    END IF

    nkind = SIZE(atomic_kind_set)

    CALL reallocate(pabt,1,maxco,1,maxco,0,nthread-1)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread-1)
    CALL reallocate(ntasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(tasks,1,8,1,max_tasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(dist_ab,1,3,1,max_tasks,1,gridlevel_info%ngrid_levels)
    CALL reallocate(tasks_local,1,2,1,max_tasks)
    CALL reallocate(ival,1,6,1,max_tasks)
    CALL reallocate(latom,1,2,1,max_tasks)
    CALL reallocate(dab_local,1,3,1,max_tasks)
    CALL reallocate(atasks,1,2,1,max_tasks)
    CALL reallocate(asets,1,2,1,max_tasks)
    curr_tasks = max_tasks

!   *** Initialize working density matrix ***


    ! distributed rs grids require a matrix that will be changed (rs_get_my_tasks)
    ! whereas this is not the case for replicated grids
    IF (distributed_rs_grids) THEN
        CALL allocate_matrix(matrix=deltap,&
                         nrow=matrix_p%nrow,&
                         ncol=matrix_p%ncol,&
                         nblock_row=matrix_p%nblock_row,&
                         nblock_col=matrix_p%nblock_col,&
                         first_row=matrix_p%first_row(:),&
                         last_row=matrix_p%last_row(:),&
                         first_col=matrix_p%first_col(:),&
                         last_col=matrix_p%last_col(:),&
                         matrix_name="DeltaP",&
                         sparsity_id=-1, &        ! basically unknown sparsity in parallel
                         matrix_symmetry=matrix_p%symmetry)
    ELSE
        deltap=>matrix_p
    ENDIF

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           softb = my_soft, &
                           orb_basis_set=orb_basis_set)

      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

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

      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           softb = my_soft, &
                           orb_basis_set=orb_basis_set)

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

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

        ab = ikind + nkind*(jkind - 1)

        IF (ASSOCIATED(sab_orb(ab)%neighbor_list_set)) THEN

           CALL get_neighbor_list_set(neighbor_list_set=&
                                      sab_orb(ab)%neighbor_list_set,&
                                      nlist=nlist)
           sab_orb_neighbor_list => first_list(sab_orb(ab)%neighbor_list_set)
        ELSE
           nlist=0
        END IF

        ntasks = 0
        tasks = 0

        DO ilist = 1, nlist

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,&
                                 atom=iatom,nnode=nnode)

          ra(:) = pbc(particle_set(iatom)%r,cell)

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

          DO inode = 1, nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab(:))

             IF (iatom <= jatom) THEN
                brow = iatom
                bcol = jatom
             ELSE
                brow = jatom
                bcol = iatom
             END IF

             ! bad, should do better loop ordering XXXXXXXXXX
             CALL get_block_node(matrix=matrix_p,&
                                 block_row=brow,&
                                 block_col=bcol,&
                                 BLOCK=p_block)

!           *** Check, if the atomic block has to be ***
!           *** calculated by the current processor  ***

             IF (.NOT.ASSOCIATED(p_block)) THEN
               sab_orb_neighbor_node => next(sab_orb_neighbor_node)
               CYCLE
             END IF

             IF (distributed_rs_grids) THEN
                 NULLIFY ( pblock )
                 CALL add_block_node ( deltap, brow, bcol, pblock )
                 pblock = p_block
             ELSE
                 pblock => p_block
             ENDIF

  
             IF (.NOT. map_consistent) THEN
                IF ( ALL ( 100.0_dp*ABS(pblock) < eps_rho_rspace) ) THEN
                  sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                  CYCLE
                END IF
             END IF

             rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
             dab = SQRT(rab2)

             DO iset=1,nseta

               IF (set_radius_a(iset) + kind_radius_b < dab) CYCLE
   
               IF (iatom == jatom) THEN
                 first_setb = iset
               ELSE
                 first_setb = 1
               END IF
 
               DO jset=first_setb,nsetb
 
                 IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
 
                 DO ipgf=1,npgfa(iset)

                   IF (rpgfa(ipgf,iset) + set_radius_b(jset) < dab) CYCLE

                   IF ((iatom == jatom).AND.(iset == jset)) THEN
                     first_pgfb = ipgf
                   ELSE
                     first_pgfb = 1
                   END IF

                   DO jpgf=first_pgfb,npgfb(jset)

                     IF (rpgfa(ipgf,iset) + rpgfb(jpgf,jset) < dab) CYCLE

                     zetp = zeta(ipgf,iset) + zetb(jpgf,jset)

                     IF (dab.lt.0.1E0_dp .AND. dft_control%qs_control%map_paa) THEN
                         igrid_level = 1
                     ELSE
                         igrid_level = gaussian_gridlevel(gridlevel_info,zetp)
                     ENDIF

                     ntasks(igrid_level) = ntasks(igrid_level) + 1
                     n = ntasks(igrid_level)
                     IF ( n > curr_tasks ) THEN
!                       curr_tasks = curr_tasks + add_tasks
                       curr_tasks = curr_tasks*mult_tasks
                       CALL reallocate(tasks,1,8,1,curr_tasks,&
                                       1,gridlevel_info%ngrid_levels)
                       CALL reallocate(dist_ab,1,3,1,curr_tasks,&
                                       1,gridlevel_info%ngrid_levels)
                     END IF
                     dir = rs_rho(igrid_level)%rs_grid%direction
                     tasks (1,n,igrid_level) = n
                     IF ( dir /= 0) THEN
                       rp(:) = ra(:) + zetb(jpgf,jset)/zetp*rab(:)
                       rp(:) = pbc(rp,cell)
                       tp = FLOOR(rp(dir)/rs_rho(igrid_level)%rs_grid%dr(dir))
                       tp = MODULO ( tp, rs_rho(igrid_level)%rs_grid%npts(dir) )
                       tasks (2,n,igrid_level) = tp + rs_rho(igrid_level)%rs_grid%lb(dir)
                     END IF
                     tasks (3,n,igrid_level) = iatom
                     tasks (4,n,igrid_level) = jatom
                     tasks (5,n,igrid_level) = iset
                     tasks (6,n,igrid_level) = jset
                     tasks (7,n,igrid_level) = ipgf
                     tasks (8,n,igrid_level) = jpgf
                     dist_ab (1,n,igrid_level) = rab(1)
                     dist_ab (2,n,igrid_level) = rab(2)
                     dist_ab (3,n,igrid_level) = rab(3)

                   END DO

                 END DO

               END DO

             END DO

             sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO

          sab_orb_neighbor_list => next(sab_orb_neighbor_list)

        END DO

        DO igrid_level = 1, gridlevel_info%ngrid_levels
          n = ntasks ( igrid_level )
          IF ( n > SIZE ( tasks_local, 2 ) ) &
            CALL reallocate(tasks_local,1,2,1,n)
          IF ( n > SIZE ( ival, 2 ) ) &
            CALL reallocate(ival,1,6,1,n)
          IF ( n > SIZE ( latom, 2 ) ) &
            CALL reallocate(latom,1,2,1,n)
          IF ( n > SIZE ( dab_local, 2 ) ) &
            CALL reallocate(dab_local,1,3,1,n)

!$OMP parallel do private(i)
          DO i=1,n
            tasks_local(1,i) = tasks(1,i,igrid_level)
            tasks_local(2,i) = tasks(2,i,igrid_level)
            latom(1,i) = tasks(3,i,igrid_level)
            latom(2,i) = tasks(4,i,igrid_level)
            ival(1,i) = tasks(3,i,igrid_level)
            ival(2,i) = tasks(4,i,igrid_level)
            ival(3,i) = tasks(5,i,igrid_level)
            ival(4,i) = tasks(6,i,igrid_level)
            ival(5,i) = tasks(7,i,igrid_level)
            ival(6,i) = tasks(8,i,igrid_level)
            dab_local(1,i) = dist_ab(1,i,igrid_level)
            dab_local(2,i) = dist_ab(2,i,igrid_level)
            dab_local(3,i) = dist_ab(3,i,igrid_level)
          END DO
!$OMP parallel do private(i)
          DO i=n+1,SIZE(tasks_local,2)
            tasks_local(1,i) = 0
            tasks_local(2,i) = 0
          END DO
          npme = 0
          ! modifies deltap if distributed_rs_grids
          CALL rs_get_my_tasks ( rs_rho(igrid_level)%rs_grid, tasks_local, &
               npme, ival=ival, rval=dab_local, pmat=deltap, pcor=latom )
          CALL rs_get_loop_vars ( npme, ival, natom_pairs, asets, atasks )
 
          IF ( nthread > 1 .AND. natom_pairs > 0 ) THEN
            lb => rs_rho(igrid_level)%rs_grid%lb_local
            ub => rs_rho(igrid_level)%rs_grid%ub_local
            lgrid%ldim = rs_rho(igrid_level)%rs_grid%ngpts_local
!$OMP parallel private(ithread,n)
!$          ithread = omp_get_thread_num()
            n = ithread*lgrid%ldim + 1
            CALL dcopy(lgrid%ldim,0._dp,0,lgrid%r(n),1)
!$OMP end parallel
          END IF
!$OMP parallel &
!$OMP default(none) &
!$OMP private(ijatoms,ithread,itask,iatom,jatom,ra,brow,bcol,p_block) &
!$OMP private(ijsets,iset,jset,ncoa,ncob,sgfa,sgfb) &
!$OMP private(work,pab,istat) &
!$OMP private(rb,rab,rab2,ipgf,jpgf,na1,na2,nb1,nb2,scale) &
!$OMP shared(maxco,maxsgf_set,natom_pairs,atasks,asets,ival,particle_set,cell) &
!$OMP shared(deltap,npgfa,npgfb,ncoset,la_max,lb_max,first_sgfa,first_sgfb) &
!$OMP shared(nsgfa,nsgfb,sphi_a,sphi_b,dab_local,la_min,lb_min,zeta,zetb) &
!$OMP shared(rs_rho,igrid_level,cube_info,l_info,eps_rho_rspace,lgrid,nthread) &
!$OMP shared(workt,pabt,my_compute_tau,map_consistent)
          ithread = 0
!$        ithread = omp_get_thread_num()
          pab => pabt(:,:,ithread)
          work => workt(:,:,ithread)
!$OMP do
          DO ijatoms = 1,natom_pairs
            itask = atasks(1,asets(1,ijatoms))
            iatom  = ival (1,itask)
            jatom  = ival (2,itask)
            ra(:) = pbc(particle_set(iatom)%r,cell)
            IF (iatom <= jatom) THEN
              brow = iatom
              bcol = jatom
            ELSE
              brow = jatom
              bcol = iatom
            END IF
            CALL get_block_node(matrix=deltap,&
                                block_row=brow,&
                                block_col=bcol,&
                                BLOCK=p_block)
            IF (.NOT.ASSOCIATED(p_block)) &
               CALL stop_program(routine,"p_block not associated in deltap")
            DO ijsets = asets(1,ijatoms), asets(2,ijatoms)
              itask = atasks(1,ijsets)
              iset   = ival (3,itask)
              jset   = ival (4,itask)
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
              DO itask = atasks(1,ijsets),atasks(2,ijsets)
                rab(:) = dab_local (:,itask)
                rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                rb(:) = ra(:) + rab(:)
                ipgf   = ival (5,itask)
                jpgf   = ival (6,itask)
                na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
                na2 = ipgf*ncoset(la_max(iset))
                nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
                nb2 = jpgf*ncoset(lb_max(jset))
 
                IF ((iatom == jatom).AND.&
                    (iset == jset).AND.&
                    (ipgf == jpgf)) THEN
                  scale = 1.0_dp
                ELSE
                  scale = 2.0_dp
                END IF

                IF ( nthread > 1 ) THEN
                  IF (iatom <= jatom) THEN
                    CALL collocate_pgf_product_rspace(&
                      la_max(iset),zeta(ipgf,iset),la_min(iset),&
                      lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                      ra,rab,rab2,scale,pab,na1-1,nb1-1,&
                      rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                      l_info,eps_rho_rspace,lgrid=lgrid,ithread=ithread, &
                      compute_tau=my_compute_tau,map_consistent=map_consistent)
                  ELSE
                    CALL collocate_pgf_product_rspace(&
                      lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                      la_max(iset),zeta(ipgf,iset),la_min(iset),&
                      rb,-rab,rab2,scale,pab,nb1-1,na1-1,&
                      rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                      l_info,eps_rho_rspace,lgrid=lgrid,ithread=ithread, &
                      compute_tau=my_compute_tau,map_consistent=map_consistent)
                  END IF
                ELSE
                  IF (iatom <= jatom) THEN
                    CALL collocate_pgf_product_rspace(&
                      la_max(iset),zeta(ipgf,iset),la_min(iset),&
                      lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                      ra,rab,rab2,scale,pab,na1-1,nb1-1,&
                      rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                      l_info,eps_rho_rspace,compute_tau=my_compute_tau, &
                      map_consistent=map_consistent)
                  ELSE
                    CALL collocate_pgf_product_rspace(&
                      lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                      la_max(iset),zeta(ipgf,iset),la_min(iset),&
                      rb,-rab,rab2,scale,pab,nb1-1,na1-1,&
                      rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                      l_info,eps_rho_rspace,compute_tau=my_compute_tau, &
                      map_consistent=map_consistent)
                  END IF
                END IF
 
              END DO
 
            END DO
 
          END DO
!$OMP end parallel
          IF ( nthread > 1 .AND. natom_pairs > 0 ) THEN
            n = (ub(1)-lb(1)+1)*(ub(2)-lb(2)+1)
            DO i=1,nthread
!$OMP parallel do &
!$OMP default(none) &
!$OMP private(j,k) &
!$OMP shared(i,lb,ub,lgrid,rs_rho,n,igrid_level)
              DO j=lb(3),ub(3)
                k = lgrid%ldim*(i-1) + n*(j-lb(3)) + 1
                CALL daxpy (n,1._dp,lgrid%r(k),1,&
                  rs_rho(igrid_level)%rs_grid%r(lb(1),lb(2),j),1)
              END DO
            END DO
          END IF

        END DO

      END DO

    END DO

!   *** Release work storage ***

    IF (distributed_rs_grids) CALL deallocate_matrix ( deltap )

    IF ( nthread > 1 ) THEN
      DEALLOCATE (lgrid%r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"lgrid%r")
    END IF

    DEALLOCATE (pabt,workt,ntasks,tasks,tasks_local,ival,latom,&
        dist_ab,dab_local,asets,atasks,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pabt,workt,ntasks,"//&
        "tasks,tasks_local,ival,latom,dist_ab,dab_local,asets,atasks")

    IF (gridlevel_info%ngrid_levels==1) THEN
       CALL rs_pw_transfer(rs_rho(1)%rs_grid,rho%pw,rs2pw)
       CALL rs_pools_give_back_rs_vect(rs_pools, rs_rho, error=error)
       CALL coeff_transform_space(rho,rho_gspace)
       IF (rho%pw%pw_grid%spherical) THEN ! rho_gspace = rho
          CALL coeff_transform_space(rho_gspace,rho)
       ENDIF
    ELSE
       DO igrid_level=1,gridlevel_info%ngrid_levels
          CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid,&
               mgrid_rspace(igrid_level)%pw,rs2pw)
       ENDDO
       CALL rs_pools_give_back_rs_vect(rs_pools, rs_rho, error=error)

       ! we want both rho and rho_gspace, the latter for Hartree and co-workers.
       SELECT CASE(interp_kind)
       CASE(linear_interp)
          CALL pw_copy(mgrid_rspace(gridlevel_info%ngrid_levels)%pw,mgrid_temp_rspace(gridlevel_info%ngrid_levels)%pw)
          DO igridlevel=gridlevel_info%ngrid_levels,2,-1
             ! prologate to the next grid and addup
             CALL pw_prolongate_l(mgrid_temp_rspace(igridlevel)%pw,mgrid_temp_rspace(igridlevel-1)%pw)
             ! add the next grid to the prolongated grid
             CALL pw_sumup(mgrid_rspace(igridlevel-1)%pw,mgrid_temp_rspace(igridlevel-1)%pw)
          ENDDO
          CALL pw_copy(mgrid_temp_rspace(1)%pw,rho%pw)
          CALL coeff_transform_space(rho,rho_gspace)
       CASE(pw_interp)
          CALL coeff_zero(rho_gspace)
          DO igrid_level=1,gridlevel_info%ngrid_levels
             CALL coeff_transform_space(mgrid_rspace(igrid_level),&
                  mgrid_gspace(igrid_level))
             CALL coeff_sumup(mgrid_gspace(igrid_level),rho_gspace)
          END DO
          CALL coeff_transform_space(rho_gspace,rho)
       CASE(spline3_pbc_interp)
          DO igrid_level=gridlevel_info%ngrid_levels,2,-1
             CALL pw_prolongate_s3(mgrid_rspace(igrid_level)%pw,&
                  mgrid_rspace(igrid_level-1)%pw,pw_pools(igrid_level)%pool,&
                  interp_section,error=error)
          END DO
          CALL pw_copy(mgrid_rspace(1)%pw,rho%pw)
          CALL coeff_transform_space(rho,rho_gspace)
       CASE default
          CALL cp_unimplemented_error(routine,"interpolator "//&
               cp_to_string(interp_kind),error=error)
       END SELECT
    END IF
    
    total_rho = calculate_total_rho(rho)

    ! *** give back the pw multi-grids
    IF (dft_control % qs_control % realspace_mgrids) THEN
       CALL pw_pools_give_back_coeffs(pw_pools,mgrid_temp_rspace,&
                                      error=error)
       DEALLOCATE(mgrid_temp_rspace,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
    ELSE
       CALL pw_pools_give_back_coeffs(pw_pools,mgrid_gspace,&
                                      error=error)
       DEALLOCATE(mgrid_gspace,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
    ENDIF
    CALL pw_pools_give_back_coeffs(pw_pools,mgrid_rspace,&
         error=error)
    DEALLOCATE(mgrid_rspace,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routine,error,failure)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_rho_elec

! *****************************************************************************
! modified calculate_rho_elec, should write the wavefunction represented by 
! the vector eigenvector on the grid. Is not efficient, in the sence that 
! it's presumably dominated by the FFT and the rs->pw and back routines
! especially in parrallel things are not efficient. 
! Currently it's only meant to provide a way to plot an MO
! *****************************************************************************
  SUBROUTINE calculate_wavefunction(mo_vectors,ivector,rho,rho_gspace,&
                                    qs_env,error)

    TYPE(cp_fm_type), POINTER                :: mo_vectors
    INTEGER                                  :: ivector
    TYPE(coeff_type), INTENT(INOUT)          :: rho, rho_gspace
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE calculate_wavefunction (MODULE qs_collocate_density)"

    INTEGER :: handle, i, iatom, igrid_level, ipgf, iset, istat, maxco, &
      maxsgf_set, na1, na2, nao, natom, ncoa, ncol_global, nseta, offset, &
      sgfa, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: failure, local
    REAL(KIND=dp)                            :: dab, eps_rho_rspace, rab2, &
                                                scale, zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvector
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab, sphi_a, work, zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_gspace, mgrid_rspace
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(l_info_type)                        :: l_info
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(REALSPACE_GRID_P_TYPE), &
      DIMENSION(:), POINTER                  :: rs_rho

    failure=.FALSE.

    CALL timeset("calculate_wavefunction","I"," ",handle)

    CALL cp_fm_get_info(matrix=mo_vectors,nrow_global=nao,ncol_global=ncol_global) 

    NULLIFY(eigenvector, cell,dft_control, orb_basis_set,&
         atomic_kind_set, sab_orb, particle_set,&
         pab,work,la_max, la_min,&
         npgfa, nsgfa, &
         sphi_a, zeta, first_sgfa,&
         rs_rho,pw_env,rs_pools,pw_pools,mgrid_rspace,mgrid_gspace)

    ALLOCATE(eigenvector(nao), stat=istat)
    IF (istat.NE.0) THEN
       CALL stop_program("calculate_wavefunction","eigenvector")
    ENDIF
    ! only some processors have non-zero elements. The sum of eigenvectors on 
    ! all cpus is the real eigenvector. combined with the rs->pw this will give the
    ! real density
    DO i=1,nao 
       CALL cp_fm_get_element(mo_vectors,i,ivector,eigenvector(i),local)
       IF (.NOT. local) eigenvector(i)=0.0_dp ! some more intelligent scheme is needed
    ENDDO

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb,&
                    pw_env=pw_env)

    ! *** set up the pw multi-grids
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routine,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, pw_pools=pw_pools, error=error)
    ALLOCATE(mgrid_rspace(SIZE(pw_pools)), mgrid_gspace(SIZE(pw_pools)),&
             stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routine,error,failure)
    CALL pw_pools_init_coeffs(pw_pools,mgrid_gspace,&
                use_data = COMPLEXDATA1D,&
                in_space = RECIPROCALSPACE, error=error)
    CALL pw_pools_init_coeffs(pw_pools,mgrid_rspace,&
                use_data = REALDATA3D,&
                in_space = REALSPACE, error=error)

    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info
    l_info=pw_env%l_info
    ! *** set up rs multi-grids
    CALL rs_pools_create_rs_vect(rs_pools, rs_rho, force_env_section=qs_env%input,&
         error=error)

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
    END DO

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

!   *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             natom=natom,&
                             maxsgf_set=maxsgf_set)

    ALLOCATE (pab(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab",maxco*1*dp_size)
    ALLOCATE (work(maxco,1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*1*dp_size)

    offset=0

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

            CALL collocate_pgf_product_rspace(&
                        la_max(iset),zeta(ipgf,iset),la_min(iset),&
                        0,0.0_dp,0,&
                        ra,rab,rab2,scale,pab,na1-1,0,&
                        rs_rho(igrid_level)%rs_grid,cube_info(igrid_level),&
                        l_info,eps_rho_rspace)
         END DO

         offset=offset+nsgfa(iset)

      END DO

    END DO

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid,&
            mgrid_rspace(igrid_level)%pw,rs2pw)
    ENDDO
    
    CALL rs_pools_give_back_rs_vect(rs_pools,rs_rho,error=error)

    CALL coeff_zero(rho_gspace)
    DO igrid_level=1,gridlevel_info%ngrid_levels
      CALL coeff_transform_space(mgrid_rspace(igrid_level),&
           mgrid_gspace(igrid_level))
      CALL coeff_sumup(mgrid_gspace(igrid_level),rho_gspace)
    END DO

    CALL coeff_transform_space(rho_gspace,rho)

!   *** Release work storage ***
    DEALLOCATE(eigenvector)

    DEALLOCATE (pab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    ! *** give back the pw multi-grids
    CALL pw_pools_give_back_coeffs(pw_pools,mgrid_gspace,&
         error=error)
    CALL pw_pools_give_back_coeffs(pw_pools,mgrid_rspace,&
         error=error)

    DEALLOCATE(mgrid_gspace, mgrid_rspace, stat=stat)
    IF (istat /= 0) CALL stop_memory(routine,"mgrid_(r/g)space")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_wavefunction

! *****************************************************************************

  FUNCTION calculate_total_rho(rho) RESULT(total_rho)

    TYPE(coeff_type), INTENT(IN)             :: rho
    REAL(KIND=dp)                            :: total_rho

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "FUNCTION calculate_total_rho (MODULE qs_collocate_density)"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(rho%pw%cc3d)) THEN
      total_rho = -rho%pw%pw_grid%vol*rho%pw%cc3d(      &
                            rho%pw%pw_grid%bounds(1,1), &
                            rho%pw%pw_grid%bounds(1,2), &
                            rho%pw%pw_grid%bounds(1,3))
    ELSE IF (ASSOCIATED(rho%pw%cr3d)) THEN
      ! do reduction using maximum accuracy
      total_rho = -rho%pw%pw_grid%dvol*SUM(REAL(rho%pw%cr3d,KIND=mp))
    ELSE IF (ASSOCIATED(rho%pw%cc)) THEN
      IF ( rho%pw%pw_grid%have_g0 ) THEN
        total_rho = -rho%pw%pw_grid%vol*rho%pw%cc(1)
      ELSE
        total_rho = 0._dp
      END IF
    ELSE
      CALL stop_program(routine,"No density coefficients available")
    END IF
    IF (rho%pw%pw_grid%para%mode.NE.PW_MODE_LOCAL) THEN
       CALL mp_sum(total_rho,rho%pw%pw_grid%para%group)
    END IF

  END FUNCTION calculate_total_rho

! *****************************************************************************

  FUNCTION calculate_total_abs_rho(rho) RESULT(total_abs_rho)

    TYPE(coeff_type), INTENT(IN)             :: rho
    REAL(KIND=dp)                            :: total_abs_rho

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "FUNCTION calculate_total_rho (MODULE qs_collocate_density)"

!   ---------------------------------------------------------------------------

    total_abs_rho = 0.0_dp
    IF (ASSOCIATED(rho%pw%cr3d)) THEN
      total_abs_rho = rho%pw%pw_grid%dvol*SUM(ABS(rho%pw%cr3d))
    ELSE
      CALL stop_program(routine,"Need density coefficients in real space !")
    END IF
    IF (rho%pw%pw_grid%para%mode.NE.PW_MODE_LOCAL) THEN
       CALL mp_sum(total_abs_rho,rho%pw%pw_grid%para%group)
    END IF

  END FUNCTION calculate_total_abs_rho

! *****************************************************************************

  SUBROUTINE collocate_pgf_product_rspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2,scale,pab,o1,o2,&
                                          rsgrid,cube_info,l_info,&
                                          eps_rho_rspace,lgrid,ithread,&
                                          compute_tau,map_consistent,&
                                          collocate_rho0,rpgf0_s)

    INTEGER, INTENT(IN)                      :: la_max
    REAL(KIND=dp), INTENT(IN)                :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(KIND=dp), INTENT(IN)                :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rab
    REAL(KIND=dp), INTENT(IN)                :: rab2, scale
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    INTEGER                                  :: o1, o2
    TYPE(realspace_grid_type),POINTER        :: rsgrid
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    TYPE(l_info_type), INTENT(IN)            :: l_info
    REAL(KIND=dp), INTENT(IN)                :: eps_rho_rspace
    TYPE(lgrid_type), OPTIONAL               :: lgrid
    INTEGER, INTENT(IN), OPTIONAL            :: ithread
    LOGICAL, INTENT(IN), OPTIONAL            :: compute_tau, map_consistent
    LOGICAL, INTENT(IN), OPTIONAL            :: collocate_rho0
    REAL(dp), INTENT(IN), OPTIONAL           :: rpgf0_s

    INTEGER :: cmax, coef_max, gridbounds(2,3), i, ico, ico_l, icoef, ig, &
      ithread_l, jco, jco_l, k, l, la_max_local, la_min_local, lb_max_local, &
      lb_min_local, length, lx, lx_max, lxa, lxb, lxy, lxy_max, lxyz, &
      lxyz_max, lya, lyb, lza, lzb, o1_local, o2_local, offset, start
    INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, ng, &
                                                ub_cube
    INTEGER, DIMENSION(:), POINTER           :: ly_max, lz_max, sphere_bounds
    INTEGER, DIMENSION(:, :), POINTER        :: map
    INTEGER, POINTER                         :: ipzyx(:,:,:,:,:,:)
    LOGICAL                                  :: my_collocate_rho0, &
                                                my_compute_tau, &
                                                my_map_consistent
    REAL(KIND=dp) :: a, b, binomial_k_lxa, binomial_l_lxb, cutoff, f, pg, &
      prefactor, radius, rpg, ya, yap, yb, ybp, za, zap, zb, zbp, zetp
    REAL(KIND=dp), DIMENSION(3)              :: dr, rap, rb, rbp, roffset, rp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab_local
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid
    REAL(KIND=dp), POINTER                   :: alpha(:,:), dpy(:,:), &
                                                dpz(:,:), polx(:,:), &
                                                poly(:,:), polz(:,:), pzyx(:)

!   ---------------------------------------------------------------------------

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


    IF (PRESENT(compute_tau)) THEN
       my_compute_tau=compute_tau
    ELSE
       my_compute_tau=.FALSE.
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
                                       zetp=zetp,eps=eps_rho_rspace,prefactor=prefactor,cutoff=cutoff)
       prefactor = scale*EXP(-zeta*f*rab2)
    ELSE IF (my_collocate_rho0) THEN
       cutoff    = 0.0_dp
       prefactor = 1.0_dp
       radius = rpgf0_s
!       radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,ra=ra,rb=rb,rp=rp,&
!                     zetp=zetp,eps=eps_rho_rspace,prefactor=prefactor,cutoff=cutoff)

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
    IF (my_compute_tau) THEN
        la_max_local=la_max+1
        la_min_local=MAX(la_min-1,0)
        lb_max_local=lb_max+1
        lb_min_local=MAX(lb_min-1,0)
        ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
        ! is equivalent to mapping pab with 0.5 * (nabla pgf_a) . (nabla pgf_b)
        ! (ddx pgf_a ) (ddx pgf_b) = (a pgf_{a-1} - 2*zeta*pgf_{a+1}*(b pgf_{x-1} - 2 * zetb * pgf_{b+1})
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
                     ico=coset(lxa,lya,lza)
                     jco=coset(lxb,lyb,lzb)
                     ! x  (all safe if lxa = 0, as the spurious added terms have zero prefactor)

                     ico_l=coset(MAX(lxa-1,0),lya,lza)
                     jco_l=coset(MAX(lxb-1,0),lyb,lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lxa* lxb*pab(o1+ico,o2+jco)
                     ico_l=coset(MAX(lxa-1,0),lya,lza)
                     jco_l=coset(   (lxb+1  ),lyb,lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp* lxa*zetb*pab(o1+ico,o2+jco)
                     ico_l=coset(   (lxa+1  ),lya,lza)
                     jco_l=coset(MAX(lxb-1,0),lyb,lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta* lxb*pab(o1+ico,o2+jco)
                     ico_l=coset(   (lxa+1  ),lya,lza)
                     jco_l=coset(   (lxb+1  ),lyb,lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zeta*zetb*pab(o1+ico,o2+jco)

                     ! y 

                     ico_l=coset(lxa,MAX(lya-1,0),lza)
                     jco_l=coset(lxb,MAX(lyb-1,0),lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lya* lyb*pab(o1+ico,o2+jco)
                     ico_l=coset(lxa,MAX(lya-1,0),lza)
                     jco_l=coset(lxb,   (lyb+1  ),lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp* lya*zetb*pab(o1+ico,o2+jco)
                     ico_l=coset(lxa,   (lya+1  ),lza)
                     jco_l=coset(lxb,MAX(lyb-1,0),lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta* lyb*pab(o1+ico,o2+jco)
                     ico_l=coset(lxa,   (lya+1  ),lza)
                     jco_l=coset(lxb,   (lyb+1  ),lzb)
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zeta*zetb*pab(o1+ico,o2+jco)

                     ! z 

                     ico_l=coset(lxa,lya,MAX(lza-1,0))
                     jco_l=coset(lxb,lyb,MAX(lzb-1,0))
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lza* lzb*pab(o1+ico,o2+jco)
                     ico_l=coset(lxa,lya,MAX(lza-1,0))
                     jco_l=coset(lxb,lyb,   (lzb+1  ))
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp* lza*zetb*pab(o1+ico,o2+jco)
                     ico_l=coset(lxa,lya,   (lza+1  ))
                     jco_l=coset(lxb,lyb,MAX(lzb-1,0))
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta* lzb*pab(o1+ico,o2+jco)
                     ico_l=coset(lxa,lya,   (lza+1  ))
                     jco_l=coset(lxb,lyb,   (lzb+1  ))
                     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zeta*zetb*pab(o1+ico,o2+jco)

              ENDDO
              ENDDO
           ENDDO
           ENDDO
        ENDDO
        ENDDO
        o1_local=0
        o2_local=0
        pab_local=pab_local * 0.5_dp
    ELSE
        la_max_local=la_max
        la_min_local=la_min
        lb_max_local=lb_max
        lb_min_local=lb_min
        pab_local => pab
        o1_local=o1
        o2_local=o2
    ENDIF

    coef_max=la_max_local+lb_max_local+1
!   *** properties of the grid ***
    dr(:) = rsgrid%dr(:)
    ng(:) = rsgrid%npts(:)
    
! WARNING: this resets the lower bounds of grid to 1 (unlike grid => rsgrid%r)
    grid => rsgrid%r(:,:,:)
    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)

!   *** get the l_info logic and arrays ***
    CALL return_l_info(l_info,la_min_local,la_max_local,lb_min_local,lb_max_local,ithread_l,lx_max, &
                       lxy_max,lxyz_max,ly_max,lz_max, &
                       map,polx,poly,polz,dpy,dpz,alpha,pzyx,ipzyx,cmax)

!   *** position of the gaussian product
!
!   this is the actual definition of the position on the grid
!   i.e. a point rp(:) gets here grid coordinates
!   MODULO(rp(:)/dr(:),ng(:))+1
!   hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid ((1,1,1) on grid)
!

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
        offset=MODULO(cubecenter(i),ng(i))+rsgrid%lb(i)
        offset=offset-rsgrid%lb_local(i)+1
        DO ig=lb_cube(i),ub_cube(i)
           map(ig,i) = ig+offset
        END DO
      END IF
    ENDDO

!   *** initialise the p terms and loop logic
    lxyz=0
    DO lxa=0,la_max_local
    DO lxb=0,lb_max_local
       DO lya=0,la_max_local-lxa
       DO lyb=0,lb_max_local-lxb
          DO lza=MAX(la_min_local-lxa-lya,0),la_max_local-lxa-lya
          DO lzb=MAX(lb_min_local-lxb-lyb,0),lb_max_local-lxb-lyb
             lxyz=lxyz+1
             ico=coset(lxa,lya,lza)
             jco=coset(lxb,lyb,lzb)
             pzyx(lxyz)=prefactor*pab_local(o1_local+ico,o2_local+jco)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
 
!   *** initialise the pol x,y,z terms
    DO ig=lb_cube(3),ub_cube(3)
      lxyz=0
      rpg = REAL(ig,dp)*dr(3) - roffset(3)
      zap = EXP(-zetp*rpg**2)
      za  = rpg + rap(3)
      zb  = za  - rab(3)
      DO lza=0,la_max_local
       zbp=1.0_dp
       DO lzb=0,lb_max_local
          dpz(lzb,lza)=zap*zbp
          zbp=zbp*zb
       ENDDO
       zap=zap*za
      ENDDO
      DO lxa=0,la_max_local
      DO lxb=0,lb_max_local
       DO lya=0,la_max_local-lxa
       DO lyb=0,lb_max_local-lxb
          DO lza=MAX(la_min_local-lxa-lya,0),la_max_local-lxa-lya
          DO lzb=MAX(lb_min_local-lxb-lyb,0),lb_max_local-lxb-lyb
             lxyz=lxyz+1
             polz(lxyz,ig)=dpz(lzb,lza)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
      ENDDO
      ENDDO 
    ENDDO

    DO ig=lb_cube(2),ub_cube(2)
      rpg = REAL(ig,dp)*dr(2) - roffset(2)
      yap = EXP(-zetp*rpg**2)
      ya  = rpg + rap(2)
      yb  = ya  - rab(2)
      DO lya=0,la_max_local
       ybp=1.0_dp
       DO lyb=0,lb_max_local
          dpy(lyb,lya)=yap*ybp
          ybp=ybp*yb
       ENDDO
       yap=yap*ya
      ENDDO

      lxy=0
      DO lxa=0,la_max_local
      DO lxb=0,lb_max_local
       DO lya=0,la_max_local-lxa
       DO lyb=0,lb_max_local-lxb
             lxy=lxy+1
             poly(lxy,ig)=dpy(lyb,lya)
       ENDDO
       ENDDO
      ENDDO
      ENDDO
    ENDDO
 
!   *** make the alpha matrix ***
    alpha(:,:)=0.0_dp
    lx=0
    DO lxa=0,la_max_local
    DO lxb=0,lb_max_local
       lx=lx+1
       binomial_k_lxa=1.0_dp
       a=1.0_dp
       DO k=0,lxa
        binomial_l_lxb=1.0_dp
        b=1.0_dp
        DO l=0,lxb
           alpha(lxa-l+lxb-k+1,lx)=alpha(lxa-l+lxb-k+1,lx)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*REAL(lxb-l,dp)/REAL(l+1,dp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*REAL(lxa-k,dp)/REAL(k+1,dp)
        a=a*(-ra(1)+rp(1))
       ENDDO
    ENDDO
    ENDDO

    DO ig=lb_cube(1),ub_cube(1)
      rpg = REAL(ig,dp)*dr(1) - roffset(1)
      pg  = EXP(-zetp*rpg**2)
      DO icoef=1,coef_max
         polx(icoef,ig)=pg
         pg=pg*(rpg)
      ENDDO
    ENDDO

    IF ( PRESENT ( lgrid ) ) THEN
      ig = lgrid%ldim * ithread_l + 1
      CALL collocate_core(pzyx(1),polx(1,-cmax),poly(1,-cmax),&
              polz(1,-cmax),lgrid%r(ig),alpha(1,1),lx_max,lxy_max,&
              lxyz_max,coef_max,cmax,ly_max(1),lz_max(1),&
              gridbounds(1,1),map(-cmax,1),sphere_bounds(1))
    ELSE
!   *** do the loop over the grid
!   notice this is not the same as critical or so, since we may have 
!   several different rsgrids that all use the same function
!   I guess we need a flush (of the grid) to guarantee that we are working 
!   on an uptodate copy of the grid
      CALL collocate_core(pzyx(1),polx(1,-cmax),poly(1,-cmax),&
              polz(1,-cmax),grid(1,1,1),alpha(1,1),lx_max,lxy_max,&
              lxyz_max,coef_max,cmax,ly_max(1),lz_max(1),&
              gridbounds(1,1),map(-cmax,1),sphere_bounds(1))
    END IF

    IF (my_compute_tau) THEN
       DEALLOCATE(pab_local)
    ENDIF

  END SUBROUTINE collocate_pgf_product_rspace

! *****************************************************************************

  SUBROUTINE collocate_pgf_product_gspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2,scale,pab,na,nb,&
                                          eps_rho_gspace,gsq_max,pw)

    ! NOTE: this routine is much slower than collocate_pgf_product_rspace

    INTEGER, INTENT(IN)                   :: la_max,la_min,lb_max,lb_min,na,nb
    REAL(dp), INTENT(IN)                  :: eps_rho_gspace,gsq_max,rab2,scale,&
                                             zeta,zetb
    REAL(dp), DIMENSION(3), INTENT(IN)    :: ra,rab
    TYPE(pw_type), POINTER                :: pw
    REAL(dp), DIMENSION(:,:), POINTER     :: pab

    ! Local parameter

    CHARACTER(LEN=*), PARAMETER :: routineN = "collocate_pgf_product_gspace",&
                                   routineP = moduleN//":"//routineN

    ! Local variables

    REAL(dp)    :: f,fa,fb,pij,prefactor,rzetp,twozetp,zetp
    INTEGER     :: ax,ay,az,bx,by,bz,handle,ico,i,ig,ig2,istat,jco,jg,kg,la,&
                   lb,lb_cube_min,lb_grid,ub_cube_max,ub_grid

    COMPLEX(dp), DIMENSION(3) :: phasefactor
    REAL(dp), DIMENSION(3)    :: dg,expfactor,fap,fbp,rap,rbp,rp
    INTEGER, DIMENSION(3)     :: lb_cube,ub_cube

    COMPLEX(dp), DIMENSION(:), POINTER :: rag,rbg
    REAL(dp), DIMENSION(:), POINTER    :: g

    COMPLEX(dp), DIMENSION(:,:,:,:), POINTER :: cubeaxis

    ! -------------------------------------------------------------------------

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

! *****************************************************************************

END MODULE qs_collocate_density
