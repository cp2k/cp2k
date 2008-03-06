!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE qs_rho0_ggrid

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_p_type,&
                                             cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs
  USE cube_utils,                      ONLY: cube_info_type
  USE f77_blas
  USE input_constants,                 ONLY: tddfpt_singlet
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: dfac,&
                                             fourpi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indco,&
                                             nco,&
                                             ncoset,&
                                             nso,&
                                             nsoset
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_integrate_function,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release
  USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_integrate_potential,          ONLY: integrate_pgf_product_rspace
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                             rho0_mpole_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE realspace_grid_types,            ONLY: pw2rs,&
                                             realspace_grid_type,&
                                             rs2pw,&
                                             rs_grid_zero,&
                                             rs_pw_transfer
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho0_ggrid'

! *** Public subroutines ***

  PUBLIC :: put_rho0_s_on_grid,rho0_s_grid_create,integrate_vhg0_rspace

 CONTAINS

! *****************************************************************************
  SUBROUTINE put_rho0_s_on_grid(qs_env, atom_set, rho0, tot_rs_int, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atom_set
    TYPE(rho0_mpole_type), POINTER           :: rho0
    REAL(dp), INTENT(OUT)                    :: tot_rs_int
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'put_rho0_s_on_grid', &
      routineP = moduleN//':'//routineN

    INTEGER :: auxbas_grid, dir, handle, iat, iatom, igrid, ikind, istat, &
      ithread, j, l0_ikind, lmax0, nat, nch_ik, nch_max, nkind, npme, tp(3)
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: eps_rho_rspace, r_V0, &
                                                rpgf0_s, zet0_s
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Qlm_c
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_rs_pool_type), POINTER           :: pool
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), POINTER                 :: coeff_gspace, coeff_rspace, &
                                                rho0_s_gs, rho0_s_rs
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(realspace_grid_type), POINTER       :: rs_grid

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure = .FALSE.
    NULLIFY(atom_kind,atom_list,cores,pab,Qlm_c)

    NULLIFY(dft_control,pw_env,particle_set,para_env,cell,rho0_s_gs,rho0_s_rs)
    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    particle_set=particle_set,&
                    para_env=para_env,&
                    pw_env=pw_env,cell=cell,error=error)
    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    NULLIFY(pools,pw_pools)
    CALL pw_env_get(pw_env=pw_env,rs_pools=pools,pw_pools=pw_pools,error=error)
    cube_info => pw_env%cube_info
    auxbas_grid=pw_env%auxbas_grid

    NULLIFY(rho0_s_gs,rho0_s_rs)
    CALL get_rho0_mpole(rho0_mpole=rho0,lmax_0=lmax0,&
                        zet0_s=zet0_s,igrid_zet0_s=igrid,&
                        rho0_s_gs=rho0_s_gs, &
                        rho0_s_rs=rho0_s_rs)

    ! *** set up the rs grid at level igrid
    NULLIFY(rs_grid,pool,pw_pool,coeff_rspace,coeff_gspace)
    pool => pools(igrid)%pool
    pw_pool => pw_pools(igrid)%pool

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)

    IF(igrid /= auxbas_grid)  THEN
      ALLOCATE(coeff_rspace,coeff_gspace ,stat=istat)
      IF(istat /= 0) CALL stop_memory("put_rho0_s_on_grid","coeff_rspace,coeff_gspace")

      CALL pw_pool_create_pw(pw_pool,coeff_rspace%pw,use_data=REALDATA3D,&
                              in_space=REALSPACE,error=error)
      CALL pw_pool_create_pw(pw_pool,coeff_gspace%pw,&
                              use_data = COMPLEXDATA1D,&
                              in_space = RECIPROCALSPACE, error=error)
    END IF
    CALL rs_pool_create_rs(pool,rs_grid,error=error)
    CALL rs_grid_zero(rs_grid)

    nch_max  = ncoset(lmax0)

    ALLOCATE(pab(nch_max,1), STAT=istat)
    IF(istat /= 0) CALL stop_memory("put_rho0_s_on_grid","pab",nch_max)

    nkind=SIZE(atom_set,1)

    DO ikind = 1,nkind
      atom_kind => atom_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,natom=nat)

      CALL get_rho0_mpole(rho0_mpole=rho0, ikind=ikind, l0_ikind=l0_ikind,&
                           rpgf0_s=rpgf0_s,r_V0hmV0s=r_V0)

      nch_ik = ncoset(l0_ikind)
      pab = 0.0_dp

      CALL reallocate ( cores, 1, nat )
      npme = 0
      cores = 0
                  
      DO iat=1,nat
         iatom = atom_list(iat)
         ra(:) = pbc(particle_set(iatom)%r,cell)
         DO dir = 1,3
            tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),ra)*rs_grid%npts(dir))
            tp(dir) = MODULO ( tp(dir), rs_grid%npts(dir) )
            tp(dir) = tp(dir) + rs_grid%lb(dir)
         END DO
         IF(rs_grid%distributed)THEN
            IF( ALL (tp.GE.rs_grid%lb_real) .AND. ALL (tp.LE.rs_grid%ub_real) ) THEN
               ! local atom
               npme = npme + 1
               cores (npme) = iat
            ENDIF
         ELSEIF(rs_grid%parallel) THEN
            ! replicated realspace grid, split the atoms up between procs
            IF (MODULO(nat,rs_grid%group_size) == rs_grid % my_pos ) THEN
               npme = npme + 1
               cores (npme) = iat
            ENDIF
         ELSE
            ! serial job
            npme = npme + 1
            cores (npme) = iat           
         ENDIF
      END DO

      ithread=0
      DO j=1,npme

        iat =  cores(j)
        iatom = atom_list(iat)

        CALL  get_rho0_mpole(rho0_mpole=rho0, iat=iatom, Qlm_car_s=Qlm_c)

        pab(1:nch_ik,1)=Qlm_c(1:nch_ik)

        ra(:) = pbc(particle_set(iatom)%r,cell)

        CALL collocate_pgf_product_rspace(&
           l0_ikind,zet0_s,0,0,0.0_dp,0,&
           ra, (/0.0_dp,0.0_dp,0.0_dp/), 0.0_dp, 1.0_dp, pab,0,0,&
           rs_grid,cell, cube_info(igrid),eps_rho_rspace,ga_gb_function=401,&
           ithread=ithread,collocate_rho0=.TRUE.,rpgf0_s=rpgf0_s,error=error)

      END DO  !  j

    END DO  ! ikind

    IF (ASSOCIATED(cores)) THEN
      DEALLOCATE (cores,STAT=istat)
      IF ( istat /= 0 ) CALL stop_memory ( "put_rho0_s_on_grid", "cores" )
    END IF

    DEALLOCATE(pab,STAT=istat)
    IF( istat/= 0 ) CALL stop_memory ("put_rho0_s_on_grid","pab")

    IF(igrid /= auxbas_grid) THEN
      CALL rs_pw_transfer(rs_grid,coeff_rspace%pw,rs2pw)
      CALL rs_pool_give_back_rs(pool,rs_grid,error=error)

      CALL pw_zero(rho0_s_gs%pw)
      CALL pw_transfer(coeff_rspace%pw,coeff_gspace%pw)
      CALL pw_axpy(coeff_gspace%pw,rho0_s_gs%pw)

      tot_rs_int = pw_integrate_function(coeff_rspace%pw,isign=-1)

      CALL pw_pool_give_back_pw(pw_pool,coeff_rspace%pw,error=error)
      CALL pw_pool_give_back_pw(pw_pool,coeff_gspace%pw,error=error)

      DEALLOCATE(coeff_rspace,coeff_gspace,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ELSE
      CALL rs_pw_transfer(rs_grid,rho0_s_rs%pw,rs2pw)
      CALL rs_pool_give_back_rs(pool,rs_grid,error=error)

      tot_rs_int = pw_integrate_function(rho0_s_rs%pw,isign=-1)

      CALL pw_transfer(rho0_s_rs%pw,rho0_s_gs%pw)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE  put_rho0_s_on_grid

! *****************************************************************************
  SUBROUTINE rho0_s_grid_create(qs_env, rho0_mpole, tddft, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    LOGICAL, INTENT(IN), OPTIONAL            :: tddft
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rho0_s_grid_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure, my_tddft
    TYPE(pw_env_type), POINTER               :: new_pw_env
    TYPE(pw_p_type), POINTER                 :: rho_gs, rho_rs
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

!   ---------------------------------------------------------------------------

  failure=.FALSE.
  my_tddft = .FALSE.
  IF(PRESENT(tddft)) my_tddft = tddft

  NULLIFY(new_pw_env)
  CALL get_qs_env(qs_env, pw_env=new_pw_env,error=error)
  CPPrecondition(ASSOCIATED(new_pw_env),cp_failure_level,routineP,error,failure)

  NULLIFY(auxbas_pw_pool)
  CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
  CPPrecondition(ASSOCIATED(auxbas_pw_pool),cp_failure_level,routineP,error,failure)

  ! reallocate rho0 on the global grid in real and reciprocal space
  NULLIFY(rho_rs,rho_gs)
  CPPrecondition(ASSOCIATED(rho0_mpole),cp_failure_level,routineP,error,failure)
  rho_rs => rho0_mpole%rho0_s_rs
  rho_gs => rho0_mpole%rho0_s_gs

  ! rho0 density in real space
  IF (ASSOCIATED(rho_rs)) THEN
     CALL pw_release(rho_rs%pw,error=error)
     DEALLOCATE(rho_rs,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  ELSE
     ALLOCATE(rho_rs,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     CALL pw_pool_create_pw(auxbas_pw_pool, rho_rs%pw,&
              use_data=REALDATA3D,in_space=REALSPACE,error=error)
     rho0_mpole%rho0_s_rs => rho_rs
  END IF

  ! rho0 density in reciprocal space
  IF (ASSOCIATED(rho_gs)) THEN
     CALL pw_release(rho_gs%pw,error=error)
     DEALLOCATE(rho_gs,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  ELSE
     ALLOCATE(rho_gs,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     CALL pw_pool_create_pw(auxbas_pw_pool, rho_gs%pw, &
          use_data=COMPLEXDATA1D, error=error)
     rho_gs%pw%in_space=RECIPROCALSPACE
     rho0_mpole%rho0_s_gs => rho_gs
  END IF

  IF(my_tddft) THEN
    rho0_mpole%igrid_zet0_s = qs_env%local_rho_set%rho0_mpole%igrid_zet0_s
  END IF

  END SUBROUTINE rho0_s_grid_create

! *****************************************************************************
  SUBROUTINE integrate_vhg0_rspace(qs_env, v_rspace, calculate_forces, &
                                   tddft, p_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type)                          :: v_rspace
    LOGICAL, INTENT(IN)                      :: calculate_forces
    LOGICAL, INTENT(IN), OPTIONAL            :: tddft
    TYPE(qs_p_env_type), OPTIONAL, POINTER   :: p_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_vhg0_rspace', &
      routineP = moduleN//':'//routineN

    INTEGER :: auxbas_grid, dir, handle, iat, iatom, ic, icg, ico, ig1, ig2, &
      igrid, ikind, ipgf1, ipgf2, is, iset1, iset2, iso, iso1, iso2, ispin, &
      istat, ithread, j, l0_ikind, lmax0, lshell, lx, ly, lz, m1, m2, maxso, &
      n1, n2, nat, nch_ik, nch_max, ncurr, npme, nset, nspins, nthread, tp(3)
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores, lmax, npgf
    LOGICAL                                  :: failure, my_tddft, paw_atom
    REAL(dp)                                 :: c4pi, eps_rho_rspace, &
                                                force_tmp(3), r_V0, ra(3), &
                                                rpgf0_s, scale, zet0_s
    REAL(dp), DIMENSION(:), POINTER          :: hab_sph, norm_l, Qlm
    REAL(dp), DIMENSION(:, :), POINTER       :: hab, hab_spht, hdab_sph, pab
    REAL(dp), DIMENSION(:, :, :), POINTER    :: habt, hadb, hdab, hdab_spht, &
                                                Qlm_gg
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: hadbt, hdabt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: rs_pools
    TYPE(cp_rs_pool_type), POINTER           :: rs_pool
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), POINTER                 :: coeff_gaux, coeff_gspace, &
                                                coeff_raux, coeff_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: pw_aux, pw_pool
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(realspace_grid_type), POINTER       :: rs_v
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_local_h, int_local_s
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(rho_atom_type), POINTER             :: rho_atom

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    c4pi = fourpi

    CALL timeset(routineN,"I","",handle)

    NULLIFY(atomic_kind_set,dft_control,para_env,particle_set)
    NULLIFY(cell,force,pw_env,rho0_mpole,rho_atom_set,cores)

    my_tddft = .FALSE.
    IF (PRESENT(tddft)) my_tddft = tddft

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    force=force,pw_env=pw_env,&
                    rho0_mpole=rho0_mpole,&
                    rho_atom_set=rho_atom_set,&
                    particle_set=particle_set,error=error)

    nspins = dft_control%nspins

    IF (my_tddft) THEN
       IF (nspins==1.AND.dft_control%tddfpt_control%res_etype/=tddfpt_singlet) RETURN
    END IF

    IF (my_tddft) THEN
       rho0_mpole   => p_env%local_rho_set%rho0_mpole
       rho_atom_set   => p_env%local_rho_set%rho_atom_set
    END IF

    CALL get_rho0_mpole(rho0_mpole=rho0_mpole,lmax_0=lmax0,&
                        zet0_s=zet0_s,igrid_zet0_s=igrid,&
                        norm_g0l_s=norm_l)

    ! *** set up of the potential on the multigrids
    NULLIFY(rs_pools, pw_pools)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_pools=rs_pools, pw_pools=pw_pools, error=error)

    ! *** assign from pw_env
    auxbas_grid=pw_env%auxbas_grid
    cube_info => pw_env%cube_info

    ! *** Get the potential on the right grid
    NULLIFY(rs_v,rs_pool,pw_pool,pw_aux,coeff_rspace,coeff_gspace)
    rs_pool => rs_pools(igrid)%pool
    pw_pool => pw_pools(igrid)%pool

    ALLOCATE(coeff_gspace, STAT=istat)
    ALLOCATE(coeff_rspace, STAT=istat)
    IF ( istat /= 0 ) CALL stop_memory ( "integrate_vhg0_rspace", &
                                         "coeff_gspace,coeff_rspace")
    CALL pw_pool_create_pw(pw_pool,coeff_gspace%pw,&
                            use_data = COMPLEXDATA1D,&
                            in_space = RECIPROCALSPACE, error=error)

    CALL pw_pool_create_pw(pw_pool,coeff_rspace%pw,use_data=REALDATA3D,&
                            in_space=REALSPACE,error=error)

    IF (igrid /= auxbas_grid) THEN
       pw_aux => pw_pools(auxbas_grid)%pool
       ALLOCATE(coeff_gaux,STAT=istat)
       IF ( istat /= 0 ) CALL stop_memory ("integrate_vhg0_rspace",&
                                           "coeff_gaux")
       CALL pw_pool_create_pw(pw_aux, coeff_gaux%pw,&
                               use_data = COMPLEXDATA1D,&
                               in_space = RECIPROCALSPACE, error=error)
       CALL pw_transfer(v_rspace%pw,coeff_gaux%pw)
       CALL pw_copy(coeff_gaux%pw,coeff_gspace%pw)
       CALL pw_transfer(coeff_gspace%pw,coeff_rspace%pw)
       CALL pw_pool_give_back_pw(pw_aux,coeff_gaux%pw,error=error)
       DEALLOCATE(coeff_gaux, STAT=istat)
       IF ( istat /= 0 ) CALL stop_memory ("integrate_vhg0_rspace",&
                                           "coeff_gaux")
       ALLOCATE(coeff_raux,STAT=istat)
       IF ( istat /= 0 ) CALL stop_memory ("integrate_vhg0_rspace",&
                                           "coeff_raux")
       CALL pw_pool_create_pw(pw_aux, coeff_raux%pw,use_data=REALDATA3D,&
                               in_space=REALSPACE,error=error)
       scale = coeff_rspace%pw%pw_grid%dvol/coeff_raux%pw%pw_grid%dvol
       coeff_rspace%pw%cr3d = scale*coeff_rspace%pw%cr3d
       CALL pw_pool_give_back_pw(pw_aux,coeff_raux%pw,error=error)
       DEALLOCATE(coeff_raux, STAT=istat)
       IF ( istat /= 0 ) CALL stop_memory ("integrate_vhg0_rspace",&
                                           "coeff_raux")
    ELSE

       IF(coeff_gspace%pw%pw_grid%spherical) THEN
         CALL pw_transfer(v_rspace%pw,coeff_gspace%pw)
         CALL pw_transfer(coeff_gspace%pw,coeff_rspace%pw)
       ELSE
         CALL pw_copy(v_rspace%pw,coeff_rspace%pw)
       END IF
    END IF
    CALL pw_pool_give_back_pw(pw_pool,coeff_gspace%pw,error=error)
    DEALLOCATE(coeff_gspace, STAT=istat)
    IF ( istat /= 0 ) CALL stop_memory ("integrate_vhg0_rspace",&
                                        "coeff_gspace")

    ! *** set up the rs grid at level igrid
    CALL rs_pool_create_rs(rs_pool,rs_v,error=error)
    CALL rs_grid_zero(rs_v)
    CALL rs_pw_transfer(rs_v,coeff_rspace%pw,pw2rs)

    CALL pw_pool_give_back_pw(pw_pool,coeff_rspace%pw,error=error)
    DEALLOCATE(coeff_rspace, STAT=istat)
    IF ( istat /= 0 ) CALL stop_memory ("integrate_vhg0_rspace",&
                                        "coeff_rspace")

    ! *** Now the potential is on the right grid => integration

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    nthread = 1
!   *** Allocate work storage ***

    NULLIFY(habt,hab_spht,hdabt,hdab_spht,hadbt,pab)
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
                          l0_ikind=l0_ikind, Qlm_gg=Qlm_gg,&
                          rpgf0_s=rpgf0_s,r_V0hmV0s=r_V0)

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             lmax=lmax, maxso=maxso,&
                             nset=nset, npgf=npgf)

!     Initialize the local KS integrals

      nch_ik = ncoset(l0_ikind)
      pab = 1.0_dp

      CALL reallocate ( cores, 1, nat )
      npme = 0
      cores = 0

      DO iat=1,nat
         iatom = atom_list(iat)
         ra(:) = pbc(particle_set(iatom)%r,cell)
         DO dir = 1,3
            tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),ra)*rs_v%npts(dir))
            tp(dir) = MODULO ( tp(dir), rs_v%npts(dir) )
            tp(dir) = tp(dir) + rs_v%lb(dir)
         END DO
         IF(rs_v%distributed)THEN
            IF( ALL (tp.GE.rs_v%lb_real) .AND. ALL (tp.LE.rs_v%ub_real) ) THEN
               ! local atom
               npme = npme + 1
               cores (npme) = iat
            ENDIF
         ELSEIF(rs_v%parallel) THEN
            ! replicated realspace grid, split the atoms up between procs
            IF (MODULO(nat,rs_v%group_size) == rs_v % my_pos ) THEN
               npme = npme + 1
               cores (npme) = iat
            ENDIF
         ELSE
            ! serial job
            npme = npme + 1
            cores (npme) = iat           
         ENDIF
      END DO

      ithread=0
      DO j=1,npme

        iat =  cores(j)
        iatom = atom_list(iat)
        ra(:) = pbc(particle_set(iatom)%r,cell)

        NULLIFY(Qlm)
        CALL  get_rho0_mpole(rho0_mpole=rho0_mpole, iat=iatom,&
                             Qlm_tot=Qlm)

        hab => habt(:,:,ithread)
        hdab => hdabt(:,:,:,ithread)
        hadb => hadbt(:,:,:,ithread)
        hab_sph => hab_spht(:,ithread)
        hdab_sph => hdab_spht(:,:,ithread)
        hab = 0.0_dp
        hdab = 0.0_dp

        CALL integrate_pgf_product_rspace(&
             l0_ikind,zet0_s,0,0,0.0_dp,0,&
             ra, (/0.0_dp,0.0_dp,0.0_dp/), 0.0_dp,rs_v,cell,&
             cube_info(igrid),hab,pab,o1=0,o2=0,&
             eps_gvg_rspace=eps_rho_rspace,&
             calculate_forces=calculate_forces,&
             hdab=hdab,hadb=hadb,ithread=ithread,&
             collocate_rho0=.TRUE.,rpgf0_s=rpgf0_s,error=error)

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
                    SQRT(c4pi*dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)/&
                    dfac(2*lshell+1))

                hdab_sph(1:3,iso) = hdab_sph(1:3,iso) +&
                    orbtramat(lshell)%c2s(is,ic)*hdab(1:3,ico,1)*&
                    norm_l(lshell)/&
                    SQRT(c4pi*dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)/&
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
                    SQRT(c4pi*dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)/&
                    dfac(2*lshell+1))

              END DO  ! ic
            END DO  ! is
          END DO  ! lshell
        END IF

! Calculate Integrals <aVb> and forces on iatom from the gradient og g0_s
        rho_atom => rho_atom_set(iatom)
        CALL get_rho_atom(rho_atom=rho_atom,&
              ga_Vlocal_gb_h=int_local_h,ga_Vlocal_gb_s=int_local_s)

        IF (calculate_forces) THEN

          force_tmp(1:3) = 0.0_dp

          DO iso = 1,nsoset(l0_ikind)
            force_tmp(1) = force_tmp(1) + Qlm(iso)*hdab_sph(1,iso)
            force_tmp(2) = force_tmp(2) + Qlm(iso)*hdab_sph(2,iso)
            force_tmp(3) = force_tmp(3) + Qlm(iso)*hdab_sph(3,iso)
          END DO

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

                      DO ispin = 1,nspins

                        int_local_h(ispin)%r_coef(ig1,ig2) = &
                                  int_local_h(ispin)%r_coef(ig1,ig2) + &
                                  Qlm_gg(ig1,ig2,iso)*hab_sph(iso)

                        int_local_s(ispin)%r_coef(ig1,ig2) = &
                                  int_local_s(ispin)%r_coef(ig1,ig2) + &
                                  Qlm_gg(ig1,ig2,iso)*hab_sph(iso)
                      END DO

                    END DO  ! icg
                  END DO  ! iso

                END DO  ! ipgf2
                m2 = m2 + maxso
              END DO  ! iset2
            END DO ! ipgf1
            m1 = m1 + maxso
          END DO  ! iset1

          force(ikind)%g0s_Vh_elec(1:3,iat) =  force(ikind)%g0s_Vh_elec(1:3,iat) + &
                                            force_tmp(1:3)

        ELSE

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

                      DO ispin = 1,nspins

                        int_local_h(ispin)%r_coef(ig1,ig2) = &
                                  int_local_h(ispin)%r_coef(ig1,ig2) + &
                                  Qlm_gg(ig1,ig2,iso)*hab_sph(iso)

                        int_local_s(ispin)%r_coef(ig1,ig2) = &
                                  int_local_s(ispin)%r_coef(ig1,ig2) + &
                                  Qlm_gg(ig1,ig2,iso)*hab_sph(iso)
                      END DO

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

    IF (ASSOCIATED(cores)) THEN
      DEALLOCATE (cores,STAT=istat)
      IF ( istat /= 0 ) CALL stop_memory ( "integrate_vhg0_rspace", "cores" )
    END IF
    DEALLOCATE (habt,hdabt,hadbt,hab_spht,hdab_spht,pab, STAT= istat)
    IF ( istat /= 0 ) CALL stop_memory ( "integrate_vhg0_rspace", &
       "habt,hdabt,hadbt,hab_spht,hdab_spht,pab")
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE integrate_vhg0_rspace

END MODULE qs_rho0_ggrid
