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
!DBG
  USE ai_overlap_new,                  ONLY: overlap
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_transform_space,&
                                             coeff_type
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
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel
  USE kinds,                           ONLY: dp
  USE l_utils,                         ONLY: l_info_type
  USE mathconstants,                   ONLY: pi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE
  USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                             rs_get_my_tasks,&
                                             rs_grid_zero,&
                                             rs_pw_transfer
  USE qs_collocate_density,            ONLY: calculate_total_rho,&
                                             collocate_pgf_product_rspace
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                             rho0_mpole_type
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

  PUBLIC :: put_rho0_s_on_grid,rho0_s_grid_create

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

    TYPE(cube_info_type)                     :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(l_info_type)                        :: l_info
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(coeff_type), POINTER                :: rho0_s_rs,rho0_s_gs
    TYPE(cp_rs_pool_type), POINTER           :: pool
    REAL(KIND=dp)                            :: eps_rho_rspace,zet0_s
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Qlm_c
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab

    INTEGER                                  :: bo(2), dir, handle, igrid, iat,&
                                                iatom, ikind, ithread, istat, j,&
                                                l0_ikind, lmax0,  nat, nch_ik,&
                                                nch_max, ncurr, nkind, npme,&
                                                omp_get_thread_num, ico
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

    NULLIFY(pools)
    CALL pw_env_get(pw_env=pw_env,rs_pools=pools,error=error)
    cube_info=pw_env%cube_info(1)
    l_info=pw_env%l_info

    CALL get_rho0_mpole(rho0_mpole=rho0,lmax_0=lmax0,&
                        zet0_s=zet0_s,igrid_zet0_s=igrid)

    ! *** set up the rs grid at level igrid
    NULLIFY(rs_grid,pool)
    pool => pools(igrid)%pool

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
!$OMP shared(cube_info,l_info,eps_rho_rspace,zet0_s,rho0,Qlm_c,nch_ik,l0_ikind)
      DO j=1,npme

        iat =  tasks(1,j)
        iatom = atom_list(iat)

        CALL  get_rho0_mpole(rho0_mpole=rho0, iat=iatom, Qlm_tot_c=Qlm_c)

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

    CALL rs_pw_transfer(rs_grid,rho0_s_rs%pw,"FORWARD")
    CALL rs_pool_give_back_rs(pool,rs_grid,error=error)

    tot_rs_int = calculate_total_rho(rho0_s_rs)
    WRITE(*,*) rho0%total_rho0_s
     stop 'back'      

    CALL coeff_transform_space(rho0_s_rs,rho0_s_gs)

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


 
!******************************************************************************
END MODULE qs_rho0_ggrid
!******************************************************************************
