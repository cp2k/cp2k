!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_rho0_ggrid
!!
!!   NAME
!!     qs_rho0_ggrid
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   SOURCE
!******************************************************************************

#include "cp_prep_globals.h"

MODULE qs_rho0_ggrid

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE coefficient_types,               ONLY: coeff_copy,&
                                             coeff_deallocate,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_sumup,coeff_zero
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_type,&
                                             cp_rs_pool_p_type,&
                                             rs_pool_give_back_rs,&
                                             rs_pool_create_rs
  USE cube_utils,                      ONLY: cube_info_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE hartree_local_types,             ONLY: ecoul_1center_type,&
                                             get_ecoul_1c
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE l_utils,                         ONLY: l_info_type
  USE mathconstants,                   ONLY: dfac,fourpi,pi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indco, &
                                             nco,ncoset,&
                                             nso,nsoset
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type, &
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE
  USE qs_collocate_density,            ONLY: calculate_total_rho,&
                                             collocate_pgf_product_rspace
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_integrate_potential,          ONLY: integrate_pgf_product_rspace
  USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                             rho0_mpole_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                             rs_get_my_tasks,&
                                             rs_grid_zero,&
                                             rs_pw_transfer
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
 
  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_rho0_ggrid"

! *** Public subroutines ***

  PUBLIC :: put_rho0_s_on_grid,rho0_s_grid_create,integrate_vhg0_rspace

!******************************************************************************

 CONTAINS

! *****************************************************************************

  SUBROUTINE put_rho0_s_on_grid(qs_env, atom_set, rho0, tot_rs_int, error)

    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(qs_environment_type), INTENT(INOUT) :: qs_env
    TYPE(rho0_mpole_type),      POINTER      :: rho0
    TYPE(atomic_kind_type), DIMENSION(:), &
                                 POINTER     :: atom_set
    REAL(dp), INTENT(OUT)                    :: tot_rs_int

    CHARACTER(LEN=*), PARAMETER         :: routine_name = "put_rho0_s_on_grid"

    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), POINTER                :: rho0_s_rs,rho0_s_gs
    TYPE(coeff_type), POINTER                :: coeff_rspace,& 
                                                coeff_gspace
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_rs_pool_type), POINTER           :: pool
    TYPE(cube_info_type)                     :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(l_info_type)                        :: l_info
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    REAL(KIND=dp)                            :: eps_rho_rspace,zet0_s
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Qlm_c
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab

    INTEGER                                  :: auxbas_grid,bo(2), dir,&
                                                handle, igrid, iat,&
                                                iatom, ikind, ithread, istat, j,&
                                                l0_ikind, lmax0,  nat, nch_ik,&
                                                nch_max, ncurr, nkind, npme,&
                                                omp_get_thread_num
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    INTEGER, DIMENSION(:, :), POINTER        :: tasks

!   ---------------------------------------------------------------------------

    CALL timeset("put_rho0_s_on_grid","I","",handle)

    NULLIFY(atom_kind,atom_list,tasks,pab,Qlm_c)

    NULLIFY(dft_control,pw_env,particle_set,para_env,cell,rho0_s_gs,rho0_s_rs)
    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    particle_set=particle_set,&
                    para_env=para_env,&
                    pw_env=pw_env,cell=cell,& 
                    rho0_s_gs=rho0_s_gs, &
                    rho0_s_rs=rho0_s_rs)
    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    NULLIFY(pools,pw_pools)
    CALL pw_env_get(pw_env=pw_env,rs_pools=pools,pw_pools=pw_pools,error=error)
    cube_info=pw_env%cube_info(1)
    l_info=pw_env%l_info
    auxbas_grid=pw_env%auxbas_grid

    CALL get_rho0_mpole(rho0_mpole=rho0,lmax_0=lmax0,&
                        zet0_s=zet0_s,igrid_zet0_s=igrid)

    ! *** set up the rs grid at level igrid
    NULLIFY(rs_grid,pool,pw_pool,coeff_rspace,coeff_gspace)
    pool => pools(igrid)%pool
    pw_pool => pw_pools(igrid)%pool

    IF(igrid /= auxbas_grid)  THEN
      CALL pw_pool_init_coeff(pw_pool,coeff_rspace,use_data=REALDATA3D,&
                              in_space=REALSPACE,error=error)
      CALL pw_pool_init_coeff(pw_pool,coeff_gspace,&
                              use_data = COMPLEXDATA1D,&
                              in_space = RECIPROCALSPACE, error=error)
    END IF
    CALL rs_pool_create_rs(pool,rs_grid,error=error)
    CALL rs_grid_zero(rs_grid)

    nch_max  = ncoset(lmax0)

    ALLOCATE(pab(nch_max,1), STAT=istat)
    IF(istat /= 0) CALL stop_memory("put_rho0_s_on_grid","pab",nch_max)

    nkind=SIZE(atom_set,1)

    ncurr = -1

    DO ikind = 1,nkind
      atom_kind => atom_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,natom=nat)

      CALL get_rho0_mpole(rho0_mpole=rho0, ikind=ikind, l0_ikind=l0_ikind)

      nch_ik = ncoset(l0_ikind)
      pab = 0.0_dp

      bo = get_limit ( nat, para_env%num_pe, para_env%mepos )

      npme = bo(2) - bo(1) + 1
      IF ( ncurr < npme ) THEN
        CALL reallocate ( tasks, 1, 2, 1, npme )
        ncurr = npme
      ELSE
        tasks = 0
      END IF

      dir = rs_grid%direction

!$OMP parallel do &
!$OMP default(none) &
!$OMP private(iat,j,iatom,ra) &
!$OMP shared(bo,tasks,dir,rs_grid,atom_list,particle_set,cell)
      DO iat=bo(1),bo(2)

        j = iat - bo(1) + 1
        tasks ( 1, j ) = iat
        IF ( dir /= 0) THEN
          iatom = atom_list(iat)
          ra(:) = pbc(particle_set(iatom)%r,cell)
          tasks ( 2, j ) = FLOOR(ra(dir)/rs_grid%dr(dir))
          tasks ( 2, j ) = MODULO ( tasks ( 2, j ), rs_grid%npts(dir) )
          tasks ( 2, j ) = tasks ( 2, j ) + rs_grid%lb(dir)
        END IF

      END DO

      CALL rs_get_my_tasks ( rs_grid, tasks, npme )

      ithread=0
!$OMP parallel do &
!$OMP default(none) &
!$OMP private(j,iat,iatom,ra,ithread) &
!$OMP shared(npme,tasks,atom_list,particle_set,cell,pab,rs_grid) &
!$OMP shared(cube_info,l_info,eps_rho_rspace,zet0_s)
      DO j=1,npme

        iat =  tasks(1,j)
        iatom = atom_list(iat)

        CALL  get_rho0_mpole(rho0_mpole=rho0, iat=iatom, Qlm_car_s=Qlm_c)

        pab(1:nch_ik,1)=Qlm_c(1:nch_ik)

        ra(:) = pbc(particle_set(iatom)%r,cell)

!$      ithread=omp_get_thread_num()
        CALL collocate_pgf_product_rspace(&
           l0_ikind,zet0_s,0,0,0.0_dp,0,&
           ra, (/0.0_dp,0.0_dp,0.0_dp/), 0.0_dp, 1.0_dp, pab,0,0,&
           rs_grid, cube_info,l_info,eps_rho_rspace,ithread=ithread)

      END DO  !  j

    END DO  ! ikind

    IF (ASSOCIATED(tasks)) THEN
      DEALLOCATE (tasks,STAT=istat)
      IF ( istat /= 0 ) CALL stop_memory ( "put_rho0_s_on_grid", "tasks" )
    END IF

    DEALLOCATE(pab,STAT=istat)
    IF( istat/= 0 ) CALL stop_memory ("put_rho0_s_on_grid","pab")


    IF(igrid /= auxbas_grid) THEN
      CALL rs_pw_transfer(rs_grid,coeff_rspace%pw,"FORWARD")
      CALL rs_pool_give_back_rs(pool,rs_grid,error=error)

      CALL coeff_zero(rho0_s_gs)
      CALL coeff_transform_space(coeff_rspace,coeff_gspace)
      CALL coeff_sumup(coeff_gspace,rho0_s_gs)

      tot_rs_int = calculate_total_rho(coeff_rspace)

      CALL pw_pool_give_back_coeff(pw_pool,coeff_rspace,error=error)
      CALL pw_pool_give_back_coeff(pw_pool,coeff_gspace,error=error)

    ELSE
      CALL rs_pw_transfer(rs_grid,rho0_s_rs%pw,"FORWARD")
      CALL rs_pool_give_back_rs(pool,rs_grid,error=error)

      tot_rs_int = calculate_total_rho(rho0_s_rs)
    
      CALL coeff_transform_space(rho0_s_rs,rho0_s_gs)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE  put_rho0_s_on_grid

! *****************************************************************************

  SUBROUTINE rho0_s_grid_create(qs_env,error)

    TYPE(qs_environment_type), &
      INTENT(inout), TARGET                  :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rho0_s_grid_create', &
      routineP = module_name//':'//routineN

    INTEGER                                  :: stat, igrid
    LOGICAL                                  :: failure
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    TYPE(coeff_type), POINTER                :: rho_rs,rho_gs
    TYPE(pw_env_type), POINTER               :: new_pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    REAL(dp)                                 :: zet

!   ---------------------------------------------------------------------------

  failure=.FALSE.
  NULLIFY(new_pw_env,auxbas_pw_pool,rho0_mpole,rho_rs,rho_gs)

! reallocate rho0 on the global grid in real and reciprocal space
  CALL get_qs_env(qs_env, pw_env=new_pw_env, rho0_mpole=rho0_mpole,&
                  rho0_s_rs=rho_rs, rho0_s_gs=rho_gs)
  CPPrecondition(ASSOCIATED(new_pw_env),cp_failure_level,routineP,error,failure)

  CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)

! rho0 density in real space
  IF (ASSOCIATED(rho_rs)) THEN
     CALL coeff_deallocate(rho_rs)
     DEALLOCATE(rho_rs,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  ELSE
     ALLOCATE(rho_rs,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     CALL pw_pool_init_coeff(pool=auxbas_pw_pool, coeff=rho_rs,&
              use_data=REALDATA3D,in_space=REALSPACE,error=error)
     rho0_mpole%rho0_s_rs => rho_rs
  END IF 

! rho0 density in reciprocal space
  IF (ASSOCIATED(rho_gs)) THEN
     CALL coeff_deallocate(rho_gs)
     DEALLOCATE(rho_gs,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  ELSE
     ALLOCATE(rho_gs,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     CALL pw_pool_init_coeff(pool=auxbas_pw_pool, coeff=rho_gs, &
          use_data=COMPLEXDATA1D, error=error)
     rho_gs%pw%in_space=RECIPROCALSPACE
     rho0_mpole%rho0_s_gs => rho_gs
  END IF

  zet = rho0_mpole%zet0_s*2.0_dp
  igrid =  gaussian_gridlevel(new_pw_env%gridlevel_info,zet)
  rho0_mpole%igrid_zet0_s = igrid

  END SUBROUTINE rho0_s_grid_create

!******************************************************************************

  SUBROUTINE integrate_vhg0_rspace(qs_env,v_rspace, calculate_forces, error)

    TYPE(qs_environment_type), INTENT(IN)    :: qs_env
    TYPE(coeff_type)                         :: v_rspace
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "integrate_vhg0_rspace",&
                     routineP = module_name//':'//routineN

    TYPE(atomic_kind_type), DIMENSION(:),&
                            POINTER          :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), POINTER                :: coeff_rspace,&
                                                coeff_gaux,&
                                                coeff_gspace,&
                                                coeff_raux
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cp_rs_pool_type), POINTER           :: rs_pool
    TYPE(cube_info_type),  &
                         POINTER             :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ecoul_1center_type), DIMENSION(:), &
                            POINTER          :: ecoul_1c
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(l_info_type)                        :: l_info
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(realspace_grid_type), POINTER       :: rs_v
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: pw_pool, pw_aux
    TYPE(rho0_mpole_type),      POINTER      :: rho0_mpole
    TYPE(rho_atom_coeff), DIMENSION(:),&
                                    POINTER  :: cpc_h, cpc_s
    TYPE(rho_atom_type), DIMENSION(:), &
                          POINTER            :: rho_atom_set
    TYPE(rho_atom_type), POINTER             :: rho_atom

    INTEGER :: auxbas_grid, bo(2), dir, handle, iat, iatom, ic, icg, ico, ig1, ig2,&
               igrid, ikind, ipgf1, ipgf2, is, iset1, iset2, iso, iso1, iso2, ispin,&
               istat, ithread, j, lmax0, l0_ikind, lshell, lx, ly, lz, m1, m2, &
               maxso, n1, n2, nat, natom, nch_ik, nch_max, ncurr, npme, nset, nspins, &
               nthread, omp_get_thread_num
    INTEGER, DIMENSION(:), POINTER           :: atom_list,&
                                                lmax,npgf
    INTEGER, DIMENSION(:,:), POINTER         :: tasks
    REAL(dp)  :: eps_rho_rspace, force_tmp(3), ra(3), scale, zet0_s
    REAL(dp), DIMENSION(:), POINTER          :: norm_l
    REAL(dp), DIMENSION(:), POINTER          :: hab_sph
    REAL(dp), DIMENSION(:,:), POINTER        :: hab_spht
    REAL(dp), DIMENSION(:,:),    POINTER     :: aVhsb_0s
    REAL(dp), DIMENSION(:,:), POINTER        :: hab, hdab_sph, pab
    REAL(dp), DIMENSION(:,:,:), POINTER      :: habt, hdab_spht
    REAL(dp), DIMENSION(:,:,:),  POINTER     :: Qlm_gg
    REAL(dp), DIMENSION(:,:,:), POINTER      :: hdab, hadb
    REAL(dp), DIMENSION(:,:,:,:), POINTER    :: hdabt, hadbt
    LOGICAL                                  :: failure, paw_atom
!   ---------------------------------------------------------------------------
    failure=.FALSE.

    CALL timeset("integrate_vhg0_rspace","I","",handle)

    NULLIFY(atomic_kind_set,dft_control,para_env,particle_set)
    NULLIFY(cell,force,pw_env,rho0_mpole,rho_atom_set,ecoul_1c,tasks)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    force=force,pw_env=pw_env,&
                    rho0_mpole=rho0_mpole,&
                    rho_atom_set=rho_atom_set,&
                    ecoul_1c=ecoul_1c,&
                    particle_set=particle_set)

    CALL get_rho0_mpole(rho0_mpole=rho0_mpole,lmax_0=lmax0,&
                        zet0_s=zet0_s,igrid_zet0_s=igrid,&
                        norm_g0l_s=norm_l)

    ! *** set up of the potential on the multigrids
    NULLIFY(rs_pools, pw_pools)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineN,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, pw_pools=pw_pools, error=error)

    ! *** assign from pw_env
    auxbas_grid=pw_env%auxbas_grid
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info(igrid)
    l_info=pw_env%l_info

    ! *** Get the potential on the right grid
    NULLIFY(rs_v,rs_pool,pw_pool,pw_aux,coeff_rspace,coeff_gspace)
    rs_pool => rs_pools(igrid)%pool
    pw_pool => pw_pools(igrid)%pool

    ALLOCATE(coeff_gspace, STAT=istat)
    ALLOCATE(coeff_rspace, STAT=istat)
    IF ( istat /= 0 ) CALL stop_memory ( "integrate_vhg0_rspace", &
                                         "coeff_gspace,coeff_rspace") 
    CALL pw_pool_init_coeff(pw_pool,coeff_gspace,&
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)

    CALL pw_pool_init_coeff(pw_pool,coeff_rspace,use_data=REALDATA3D,&
                            in_space=REALSPACE,error=error)

    IF (igrid /= auxbas_grid) THEN
       pw_aux => pw_pools(auxbas_grid)%pool
       ALLOCATE(coeff_gaux,STAT=istat)
       IF ( istat /= 0 ) CALL stop_memory ("integrate_vhg0_rspace",& 
                                           "coeff_gaux")
       CALL pw_pool_init_coeff(pw_aux, coeff_gaux,&
                               use_data = COMPLEXDATA1D,&
                               in_space = RECIPROCALSPACE, error=error)
       CALL coeff_transform_space(v_rspace,coeff_gaux)
       CALL coeff_copy(coeff_gaux,coeff_gspace)
       CALL coeff_transform_space(coeff_gspace,coeff_rspace)
       CALL pw_pool_give_back_coeff(pw_aux,coeff_gaux,error=error)
       DEALLOCATE(coeff_gaux, STAT=istat)
       IF ( istat /= 0 ) CALL stop_memory ("integrate_vhg0_rspace",&
                                           "coeff_gaux")
       ALLOCATE(coeff_raux,STAT=istat)
       IF ( istat /= 0 ) CALL stop_memory ("integrate_vhg0_rspace",&
                                           "coeff_raux")
       CALL pw_pool_init_coeff(pw_aux, coeff_raux,use_data=REALDATA3D,&
                               in_space=REALSPACE,error=error)
       scale = coeff_rspace%pw%pw_grid%dvol/coeff_raux%pw%pw_grid%dvol
       coeff_rspace%pw%cr3d = scale*coeff_rspace%pw%cr3d
       CALL pw_pool_give_back_coeff(pw_aux,coeff_raux,error=error)
       DEALLOCATE(coeff_raux, STAT=istat)
       IF ( istat /= 0 ) CALL stop_memory ("integrate_vhg0_rspace",&
                                           "coeff_raux")
    ELSE

       IF(coeff_gspace%pw%pw_grid%spherical) THEN
         CALL coeff_transform_space(v_rspace,coeff_gspace)
         CALL coeff_transform_space(coeff_gspace,coeff_rspace)
       ELSE
         CALL coeff_copy(v_rspace,coeff_rspace)
       END IF
    END IF
    CALL pw_pool_give_back_coeff(pw_pool,coeff_gspace,error=error)
    DEALLOCATE(coeff_gspace, STAT=istat)
    IF ( istat /= 0 ) CALL stop_memory ("integrate_vhg0_rspace",&
                                        "coeff_gspace")

    ! *** set up the rs grid at level igrid
    CALL rs_pool_create_rs(rs_pool,rs_v,error=error)
    CALL rs_grid_zero(rs_v)
    CALL rs_pw_transfer(rs_v,coeff_rspace%pw,"BACKWARD")

    CALL pw_pool_give_back_coeff(pw_pool,coeff_rspace,error=error)
    DEALLOCATE(coeff_rspace, STAT=istat)
    IF ( istat /= 0 ) CALL stop_memory ("integrate_vhg0_rspace",&
                                        "coeff_rspace")

    ! *** Now the potential is on the right grid => integration

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    nthread = 1
!$  nthread = omp_get_max_threads()

!   *** Allocate work storage ***

    NULLIFY(habt,hdabt,hdab_spht,hadbt,pab)
    nch_max = ncoset(lmax0)
    CALL reallocate(habt,1,nch_max,1,1,0,nthread)
    CALL reallocate(hab_spht,1,nch_max,0,nthread)
    CALL reallocate(hdabt,1,3,1,nch_max,1,1,0,nthread)
    CALL reallocate(hadbt,1,3,1,nch_max,1,1,0,nthread)
    CALL reallocate(hdab_spht,1,3,1,nch_max,0,nthread)
    CALL reallocate(pab,1,nch_max,1,1)

    ncurr = -1

    DO ikind = 1,SIZE(atomic_kind_set,1)
      atom_kind => atomic_kind_set(ikind)
      NULLIFY(orb_basis_set,atom_list,harmonics) 
      CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
                           natom=nat,paw_atom=paw_atom,&
                           orb_basis_set=orb_basis_set,&
                           harmonics=harmonics)

      IF(.NOT. paw_atom) CYCLE

      NULLIFY(Qlm_gg,lmax,npgf)
      CALL get_rho0_mpole(rho0_mpole=rho0_mpole,ikind=ikind,&
                          l0_ikind=l0_ikind, Qlm_gg=Qlm_gg)

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             lmax=lmax, maxso=maxso,&
                             nset=nset, npgf=npgf)
                             
      nch_ik = ncoset(l0_ikind)
      pab = 1.0_dp

      bo = get_limit ( nat, para_env%num_pe, para_env%mepos )

      npme = bo(2) - bo(1) + 1
      IF ( ncurr < npme ) THEN
        CALL reallocate ( tasks, 1, 2, 1, npme )
        ncurr = npme
      ELSE
        tasks = 0
      END IF

      dir = rs_v%direction

!$OMP parallel do &
!$OMP default(none) &
!$OMP private(iat,j,iatom,ra) &
!$OMP shared(bo,tasks,dir,rs_v,atom_list,particle_set,cell)
      DO iat=bo(1),bo(2)

        j = iat - bo(1) + 1
        tasks ( 1, j ) = iat
        IF ( dir /= 0) THEN
          iatom = atom_list(iat)
          ra(:) = pbc(particle_set(iatom)%r,cell)
          tasks ( 2, j ) = FLOOR(ra(dir)/rs_v%dr(dir))
          tasks ( 2, j ) = MODULO ( tasks ( 2, j ), rs_v%npts(dir) )
          tasks ( 2, j ) = tasks ( 2, j ) + rs_v%lb(dir)
        END IF

      END DO

      CALL rs_get_my_tasks ( rs_v, tasks, npme )

      ithread=0
!$OMP parallel do &
!$OMP default(none) &
!$OMP private(j,iat,iatom,ra,ithread) &
!$OMP private(lshell,is,iso,hab,hdab,hadb,hab_sph,hdab_sph)&
!$OMP private(lx,ly,lz,ic,ico,rho_atom,cpc_h,cpc_s)&
!$OMP private(iset1,iset2,ipgf1,ipgf2,n1,n2,m1,m2,ig1,ig2,iso1,iso2)&
!$OMP private(force_tmp,aVhsb_0s)&
!$OMP shared(npme,tasks,atom_list,particle_set,cell,pab,rs_v) &
!$OMP shared(habt,hdabt,hadbt,hab_spht,hdab_spht)&
!$OMP shared(cube_info,l_info,eps_rho_rspace,zet0_s,l0_ikind)
!$OMP shared(norm_l,calculate_forces,fourpi,rho_atom_set,ecoul_1c)&
!$OMP shared(nset,npgf,lmax,maxso,harmonics,Olm_gg,force)
      DO j=1,npme

!$      ithread=omp_get_thread_num()
        iat =  tasks(1,j)
        iatom = atom_list(iat)
        ra(:) = pbc(particle_set(iatom)%r,cell) 

        hab => habt(:,:,ithread)
        hdab => hdabt(:,:,:,ithread)
        hadb => hadbt(:,:,:,ithread)
        hab_sph => hab_spht(:,ithread)
        hdab_sph => hdab_spht(:,:,ithread)
        hab = 0.0_dp
        hdab = 0.0_dp

        CALL integrate_pgf_product_rspace(&
             l0_ikind,zet0_s,0,0,0.0_dp,0,&
             ra, (/0.0_dp,0.0_dp,0.0_dp/), 0.0_dp,rs_v,&
             cube_info,l_info,hab,pab,o1=0,o2=0,&
             eps_gvg_rspace=eps_rho_rspace,&
             calculate_forces=calculate_forces,&
             hdab=hdab,hadb=hadb,ithread=ithread)

! Convert from cartesian to spherical
        IF (calculate_forces) THEN
          DO lshell = 0,l0_ikind
            DO is = 1,nso(lshell)
              iso = is + nsoset(lshell-1) 
              hab_sph(iso) = 0.0_dp
              hdab_sph(1:3,iso) = 0.0_dp
              DO ic = 1,nco(lshell)
                ico = ic + ncoset(lshell-1)
                lx = indco(1,ico)
                ly = indco(2,ico)
                lz = indco(3,ico)
  
                hab_sph(iso) = hab_sph(iso) + &
                    orbtramat(lshell)%c2s(is,ic)*hab(ico,1)*&
                    norm_l(lshell)/&
                    SQRT(fourpi*dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)/&
                    dfac(2*lshell+1))

                hdab_sph(1:3,iso) = hdab_sph(1:3,iso) +&
                    orbtramat(lshell)%c2s(is,ic)*hdab(1:3,ico,1)*&
                    norm_l(lshell)/&
                    SQRT(fourpi*dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)/&
                    dfac(2*lshell+1))

              END DO  ! ic
            END DO  ! is
          END DO  ! lshell
        ELSE
          DO lshell = 0,l0_ikind
            DO is = 1,nso(lshell)
              iso = is + nsoset(lshell-1) 
              hab_sph(iso) = 0.0_dp
              DO ic = 1,nco(lshell)
                ico = ic + ncoset(lshell-1)
                lx = indco(1,ico)
                ly = indco(2,ico)
                lz = indco(3,ico)
  
                hab_sph(iso) = hab_sph(iso) + &
                    orbtramat(lshell)%c2s(is,ic)*hab(ico,1)*&
                    norm_l(lshell)/&
                    SQRT(fourpi*dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)/&
                    dfac(2*lshell+1))

              END DO  ! ic
            END DO  ! is
          END DO  ! lshell
        END IF

! Calculate Integrals <aVb> and forces on iatom from the gradient og g0_s 
        IF (calculate_forces) THEN
          NULLIFY(rho_atom,cpc_h,cpc_s,aVhsb_0s)
          rho_atom => rho_atom_set(iatom)
          CALL get_rho_atom(rho_atom=rho_atom,cpc_h=cpc_h,cpc_s=cpc_s)
          nspins = SIZE(cpc_h,1)

          CALL  get_ecoul_1c(ecoul_1c,iatom=iatom,aVhsb_0s=aVhsb_0s)

          aVhsb_0s = 0.0_dp
          force_tmp(1:3) = 0.0_dp
          m1 = 0
          DO iset1 = 1,nset
            n1 = nsoset(lmax(iset1))
            DO ipgf1 = 1,npgf(iset1)
              m2 = 0
              DO iset2 = 1,nset
                n2 = nsoset(lmax(iset2))
                DO ipgf2 = 1,npgf(iset2)

                  DO iso = 1,nsoset(l0_ikind)
                    DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                      iso1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                      iso2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)
                      ig1 = iso1 + n1*(ipgf1-1) + m1
                      ig2 = iso2 + n2*(ipgf2-1) + m2

                      ! sum up the contribution
                      aVhsb_0s(ig1,ig2) = aVhsb_0s(ig1,ig2) +&
                                     Qlm_gg(ig1,ig2,iso)*hab_sph(iso)

                      force_tmp(1:3) = force_tmp(1:3) +&
                             (cpc_h(1)%r_coef(ig1,ig2)-cpc_s(1)%r_coef(ig1,ig2))*&
                              Qlm_gg(ig1,ig2,iso)*hdab_sph(1:3,iso)
                      DO ispin = 2,nspins
                        force_tmp(1:3) = force_tmp(1:3) +&
                             (cpc_h(ispin)%r_coef(ig1,ig2)-cpc_s(ispin)%r_coef(ig1,ig2))*&
                              Qlm_gg(ig1,ig2,iso)*hdab_sph(1:3,iso)
                      END DO
                    END DO  ! icg
                  END DO  ! iso

                END DO  ! ipgf2 
                m2 = m2 + maxso
              END DO  ! iset2
            END DO ! ipgf1
            m1 = m1 + maxso
          END DO  ! iset1

!$OMP critical (qs_core_force_g0s_Vh_elec)
          force(ikind)%g0s_Vh_elec(1:3,iat) =  force(ikind)%g0s_Vh_elec(1:3,iat) + &
                                            force_tmp(1:3)
!$OMP end critical (qs_core_force_g0s_Vh_elec)

        ELSE

          CALL  get_ecoul_1c(ecoul_1c,iatom=iatom,aVhsb_0s=aVhsb_0s)

          aVhsb_0s = 0.0_dp
          m1 = 0
          DO iset1 = 1,nset
            n1 = nsoset(lmax(iset1))
            DO ipgf1 = 1,npgf(iset1)
              m2 = 0
              DO iset2 = 1,nset
                n2 = nsoset(lmax(iset2))
                DO ipgf2 = 1,npgf(iset2)

                  DO iso = 1,nsoset(l0_ikind)
                    DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                      iso1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                      iso2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)
                      ig1 = iso1 + n1*(ipgf1-1) + m1
                      ig2 = iso2 + n2*(ipgf2-1) + m2
                      ! sum up the contribution
                      aVhsb_0s(ig1,ig2) = aVhsb_0s(ig1,ig2) +&
                                     Qlm_gg(ig1,ig2,iso)*hab_sph(iso)
                    END DO  ! icg
                  END DO  ! iso

                END DO  ! ipgf2 
                m2 = m2 + maxso
              END DO  ! iset2
            END DO ! ipgf1
            m1 = m1 + maxso
          END DO  ! iset1

        END IF

      END DO  ! j
    END DO  ! ikind  

    CALL rs_pool_give_back_rs(rs_pool,rs_v,error=error)

    IF (ASSOCIATED(tasks)) THEN
      DEALLOCATE (tasks,STAT=istat)
      IF ( istat /= 0 ) CALL stop_memory ( "integrate_vhg0_rspace", "tasks" )
    END IF
    DEALLOCATE (habt,hdabt,hadbt,hab_spht,hdab_spht, STAT= istat)
     IF ( istat /= 0 ) CALL stop_memory ( "integrate_vhg0_rspace", "hab")
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE integrate_vhg0_rspace
 
!******************************************************************************
END MODULE qs_rho0_ggrid
!******************************************************************************
