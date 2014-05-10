!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the plane wave density by collocating the primitive Gaussian
!>      functions (pgf).
!> \par History
!>      - rewrote collocate for increased accuracy and speed
!>      - introduced the PGI hack for increased speed with that compiler
!>        (22.02.02)
!>      - Added Multiple Grid feature
!>      - new way to get over the grid (01.03.02)
!>      - removed timing calls since they were getting expensive
!>      - Updated with the new QS data structures (09.04.02,MK)
!>      - introduction of the real space grid type ( prelim. version JVdV 05.02)
!>      - parallel FFT (JGH 22.05.02)
!>      - multigrid arrays independent from density (JGH 30.08.02)
!>      - old density stored in g space (JGH 30.08.02)
!>      - distributed real space code (JGH 17.07.03)
!>      - refactoring and new loop ordering (JGH 23.11.03)
!>      - OpenMP parallelization (JGH 03.12.03)
!>      - Modified to compute tau (Joost 12.03)
!>      - removed the incremental density rebuild (Joost 01.04)
!>      - introduced realspace multigridding (Joost 02.04)
!>      - introduced map_consistent (Joost 02.04)
!>      - Addition of the subroutine calculate_atomic_charge_density (TdK, 08.05)
!>      - rewrite of the collocate/integrate kernels (Joost VandeVondele, 03.07)
!> \author Matthias Krack (03.04.2001)
!>      1) Joost VandeVondele (01.2002)
!>      Thomas D. Kuehne (04.08.2005)
! *****************************************************************************
MODULE qs_collocate_density
  USE ao_util,                         ONLY: exp_radius_very_extended
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_deallocate_matrix
  USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cube_utils,                      ONLY: compute_cube_center,&
                                             cube_info_type,&
                                             return_cube,&
                                             return_cube_nonortho
  USE d3_poly,                         ONLY: poly_cp2k2d3
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE gauss_colloc,                    ONLY: collocGauss
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE input_constants,                 ONLY: use_aux_basis_set,&
                                             use_aux_fit_basis_set,&
                                             use_lri_basis_set,&
                                             use_orb_basis_set,&
                                             use_ri_aux_basis_set
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_8
  USE lgrid_types,                     ONLY: lgrid_allocate_grid,&
                                             lgrid_type
  USE lri_environment_types,           ONLY: lri_environment_type,&
                                             lri_kind_type
  USE mathconstants,                   ONLY: dfac,&
                                             fac,&
                                             pi,&
                                             twopi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: coset,&
                                             indco,&
                                             nco,&
                                             ncoset,&
                                             nso,&
                                             nsoset
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_integrate_function,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_create_pws,&
                                             pw_pools_give_back_pws
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_modify_pab_block,             ONLY: &
       FUNC_AB, FUNC_ADBmDAB, FUNC_ARB, FUNC_ARDBmDARB, FUNC_DABpADB, &
       FUNC_DADB, FUNC_DER, FUNC_DX, FUNC_DXDX, FUNC_DXDY, FUNC_DY, &
       FUNC_DYDY, FUNC_DYDZ, FUNC_DZ, FUNC_DZDX, FUNC_DZDZ, &
       prepare_adb_m_dab, prepare_arb, prepare_ardb_m_darb, &
       prepare_dab_p_adb, prepare_dadb, prepare_diadib, prepare_diiadiib, &
       prepare_dijadijb
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE realspace_grid_types,            ONLY: &
       realspace_grid_desc_p_type, realspace_grid_desc_type, &
       realspace_grid_p_type, realspace_grid_type, rs2pw, rs_grid_release, &
       rs_grid_retain, rs_grid_zero, rs_pw_transfer
  USE rs_pw_interface,                 ONLY: density_rs2pw,&
                                             density_rs2pw_basic
  USE scptb_types,                     ONLY: get_scptb_parameter,&
                                             scp_vector_type,&
                                             scptb_parameter_type
  USE task_list_methods,               ONLY: int2pair,&
                                             rs_distribute_matrix
  USE task_list_types,                 ONLY: task_list_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_collocate_density'
! *** Public subroutines ***

  PUBLIC :: calculate_ppl_grid,&
            calculate_rho_core,&
            calculate_lri_rho_elec,&
            calculate_rho_single_gaussian,&
            calculate_rho_metal,&
            calculate_rho_resp_single,&
            calculate_rho_resp_all,&
            calculate_rho_elec,&
            calculate_drho_elec,&
            calculate_scp_charge,&
            calculate_wavefunction,&
            collocate_pgf_product_gspace,&
            collocate_pgf_product_rspace,&
            collocate_atomic_charge_density,&
            calculate_rho_nlcc

  INTEGER :: debug_count=0

CONTAINS

! *****************************************************************************
!> \brief computes the SCP charge distribution on the grid
! *****************************************************************************
  SUBROUTINE calculate_scp_charge(scp,qs_env,scpv,error)

    TYPE(pw_p_type), INTENT(INOUT)           :: scp
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(scp_vector_type), POINTER           :: scpv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_scp_charge', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, i, iatom, &
                                                ierr, ii, ikind, ithread, j, &
                                                jj, l, lmaxscp, natom, ni, &
                                                nj, npme, nthread
    INTEGER(KIND=int_8)                      :: subpatch_pattern
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores
    LOGICAL                                  :: defined, failure
    REAL(KIND=dp)                            :: alpha, eps_rho_rspace, norm, &
                                                pp
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type)                     :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_rho
    TYPE(scptb_parameter_type), POINTER      :: scptb_kind

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(atomic_kind,cell,dft_control,pab,atomic_kind_set,particle_set,&
         atom_list,pw_env,rs_rho,auxbas_rs_desc,auxbas_pw_pool,cores)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    cell=cell,&
                    particle_set=particle_set,&
                    para_env=para_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env,auxbas_rs_desc=auxbas_rs_desc,auxbas_rs_grid=rs_rho,&
         auxbas_pw_pool=auxbas_pw_pool,error=error)
    cube_info=pw_env%cube_info(1)
    CALL rs_grid_retain(rs_rho, error=error)
    CALL rs_grid_zero(rs_rho)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,atom_list=atom_list,&
                           scptb_parameter=scptb_kind)
      CALL get_scptb_parameter(scptb_kind,defined=defined,lmaxscp=lmaxscp,ag=alpha)

      IF (.NOT.defined) CYCLE

      ni = ncoset(lmaxscp)
      ALLOCATE(pab(ni,1),STAT=ierr)
      CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)

      nthread = 1
      ithread=0

      ALLOCATE(cores(natom),STAT=ierr)
      CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
      cores = 0
      npme = 0

      DO iatom = 1, natom
         atom_a = atom_list(iatom)
         ra(:) = pbc(particle_set(atom_a)%r,cell)
         IF(rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
            ! replicated realspace grid, split the atoms up between procs
            IF (MODULO(iatom,rs_rho%desc%group_size) == rs_rho % desc % my_pos ) THEN
               npme = npme + 1
               cores (npme) = iatom
            ENDIF
         ELSE
            npme = npme + 1
            cores (npme) = iatom
         ENDIF
      END DO

      DO nj=1,npme

         iatom = cores(nj)
         atom_a = atom_list(iatom)
         ra(:) = pbc(particle_set(atom_a)%r,cell)
         subpatch_pattern=0
         pab = 0._dp
         DO l=0,lmaxscp
            pp = (2._dp*l+3._dp)/2._dp
            norm = 2._dp**(l+2)/SQRT(pi)/dfac(2*l+1)
            norm = SQRT(0.25_dp*dfac(2*l+1)/pi) * norm*alpha**pp
            DO jj=1,nco(l)
               j = ncoset(l-1) + jj
               DO ii=1,nso(l)
                  i = nsoset(l-1) + ii
                  pab(j,1) = pab(j,1) + scpv%vector(ikind)%vmat(i,iatom)*&
                                        orbtramat(l)%s2c(ii,jj)*norm
               END DO
            END DO
         END DO
         CALL collocate_pgf_product_rspace(lmaxscp,alpha,0,0,0.0_dp,0,ra,&
            (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,1.0_dp,pab,0,0,rs_rho,&
            cell,cube_info,eps_rho_rspace,ga_gb_function=FUNC_AB,&
            ithread=ithread,use_subpatch=.TRUE.,subpatch_pattern=subpatch_pattern,error=error)

      END DO

      DEALLOCATE ( pab, cores, STAT=ierr )
      CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)

    END DO

    CALL rs_pw_transfer(rs_rho,scp%pw,rs2pw,error=error)
    CALL rs_grid_release(rs_rho, error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_scp_charge
! *****************************************************************************
!> \brief computes the density of the non-linear core correction on the grid
! *****************************************************************************
  SUBROUTINE calculate_rho_nlcc(rho_nlcc,qs_env,error)

    TYPE(pw_p_type), INTENT(INOUT)           :: rho_nlcc
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_nlcc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, ierr, &
                                                iexp_nlcc, ikind, ithread, j, &
                                                n, natom, nc, nexp_nlcc, ni, &
                                                npme, nthread
    INTEGER(KIND=int_8)                      :: subpatch_pattern
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores, nct_nlcc
    LOGICAL                                  :: failure, nlcc
    REAL(KIND=dp)                            :: alpha, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha_nlcc
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cval_nlcc, pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type)                     :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_rho

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(atomic_kind,cell,dft_control,pab,atomic_kind_set,particle_set,&
         atom_list,pw_env,rs_rho,auxbas_rs_desc,auxbas_pw_pool,cores)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    para_env=para_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env,auxbas_rs_desc=auxbas_rs_desc,auxbas_rs_grid=rs_rho,&
         auxbas_pw_pool=auxbas_pw_pool,error=error)
    cube_info=pw_env%cube_info(1)
    ! be careful in parallel nsmax is choosen with multigrid in mind!
    CALL rs_grid_retain(rs_rho, error=error)
    CALL rs_grid_zero(rs_rho)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           atom_list=atom_list,&
                           gth_potential=gth_potential)

      IF (.NOT.ASSOCIATED(gth_potential)) CYCLE
      CALL get_potential(potential=gth_potential,nlcc_present=nlcc,nexp_nlcc=nexp_nlcc,&
                         alpha_nlcc=alpha_nlcc,nct_nlcc=nct_nlcc,cval_nlcc=cval_nlcc)

      IF ( .NOT. nlcc ) CYCLE

      DO iexp_nlcc=1,nexp_nlcc

         alpha=alpha_nlcc(iexp_nlcc)
         nc=nct_nlcc(iexp_nlcc)

         ni = ncoset(2*nc-2)
         ALLOCATE(pab(ni,1),STAT=ierr)
         CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
         pab = 0._dp

         nthread = 1
         ithread=0

         CALL reallocate ( cores, 1, natom )
         npme = 0
         cores = 0

         ! prepare core function
         DO j=1,nc
           SELECT CASE (j)
             CASE (1)
               pab(1,1) = cval_nlcc(1,iexp_nlcc)
             CASE (2)
               n = coset(2,0,0)
               pab(n,1) = cval_nlcc(2,iexp_nlcc)/alpha**2
               n = coset(0,2,0)
               pab(n,1) = cval_nlcc(2,iexp_nlcc)/alpha**2
               n = coset(0,0,2)
               pab(n,1) = cval_nlcc(2,iexp_nlcc)/alpha**2
             CASE (3)
               n = coset(4,0,0)
               pab(n,1) = cval_nlcc(3,iexp_nlcc)/alpha**4
               n = coset(0,4,0)
               pab(n,1) = cval_nlcc(3,iexp_nlcc)/alpha**4
               n = coset(0,0,4)
               pab(n,1) = cval_nlcc(3,iexp_nlcc)/alpha**4
               n = coset(2,2,0)
               pab(n,1) = 2._dp*cval_nlcc(3,iexp_nlcc)/alpha**4
               n = coset(2,0,2)
               pab(n,1) = 2._dp*cval_nlcc(3,iexp_nlcc)/alpha**4
               n = coset(0,2,2)
               pab(n,1) = 2._dp*cval_nlcc(3,iexp_nlcc)/alpha**4
             CASE (4)
               n = coset(6,0,0)
               pab(n,1) = cval_nlcc(4,iexp_nlcc)/alpha**6
               n = coset(0,6,0)
               pab(n,1) = cval_nlcc(4,iexp_nlcc)/alpha**6
               n = coset(0,0,6)
               pab(n,1) = cval_nlcc(4,iexp_nlcc)/alpha**6
               n = coset(4,2,0)
               pab(n,1) = 3._dp*cval_nlcc(4,iexp_nlcc)/alpha**6
               n = coset(4,0,2)
               pab(n,1) = 3._dp*cval_nlcc(4,iexp_nlcc)/alpha**6
               n = coset(2,4,0)
               pab(n,1) = 3._dp*cval_nlcc(4,iexp_nlcc)/alpha**6
               n = coset(2,0,4)
               pab(n,1) = 3._dp*cval_nlcc(4,iexp_nlcc)/alpha**6
               n = coset(0,4,2)
               pab(n,1) = 3._dp*cval_nlcc(4,iexp_nlcc)/alpha**6
               n = coset(0,2,4)
               pab(n,1) = 3._dp*cval_nlcc(4,iexp_nlcc)/alpha**6
               n = coset(2,2,2)
               pab(n,1) = 6._dp*cval_nlcc(4,iexp_nlcc)/alpha**6
             CASE DEFAULT
               CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
           END SELECT
         END DO
         IF (dft_control%nspins==2) pab=pab*0.5_dp

         DO iatom = 1, natom
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r,cell)
            IF(rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
               ! replicated realspace grid, split the atoms up between procs
               IF (MODULO(iatom,rs_rho%desc%group_size) == rs_rho % desc % my_pos ) THEN
                  npme = npme + 1
                  cores (npme) = iatom
               ENDIF
            ELSE
               npme = npme + 1
               cores (npme) = iatom
            ENDIF
         END DO

         DO j=1,npme

            iatom = cores(j)
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r,cell)
            subpatch_pattern=0
            ni = 2*nc-2
            CALL collocate_pgf_product_rspace(ni,1/(2*alpha**2),0,0,0.0_dp,0,ra,&
               (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,1.0_dp,pab,0,0,rs_rho,&
               cell,cube_info,eps_rho_rspace,ga_gb_function=FUNC_AB,&
               ithread=ithread,use_subpatch=.TRUE.,subpatch_pattern=subpatch_pattern,error=error)

         END DO

         DEALLOCATE ( pab, STAT=ierr )
         CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)

      END DO

    END DO

    IF (ASSOCIATED(cores)) THEN
      DEALLOCATE (cores,STAT=ierr)
      CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL rs_pw_transfer(rs_rho,rho_nlcc%pw,rs2pw,error=error)
    CALL rs_grid_release(rs_rho, error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_rho_nlcc

! *****************************************************************************
!> \brief computes the local pseudopotential (without erf term) on the grid
! *****************************************************************************
  SUBROUTINE calculate_ppl_grid(vppl,qs_env,error)

    TYPE(pw_p_type), INTENT(INOUT)           :: vppl
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_ppl_grid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, ierr, &
                                                ikind, ithread, j, lppl, n, &
                                                natom, ni, npme, nthread
    INTEGER(KIND=int_8)                      :: subpatch_pattern
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cexp_ppl
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type)                     :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_rho

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(atomic_kind,cell,dft_control,pab,atomic_kind_set,particle_set,&
         atom_list,pw_env,rs_rho,auxbas_rs_desc,auxbas_pw_pool,cores)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    para_env=para_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env,auxbas_rs_desc=auxbas_rs_desc,auxbas_rs_grid=rs_rho,&
         auxbas_pw_pool=auxbas_pw_pool,error=error)
    cube_info=pw_env%cube_info(1)
    ! be careful in parallel nsmax is choosen with multigrid in mind!
    CALL rs_grid_retain(rs_rho, error=error)
    CALL rs_grid_zero(rs_rho)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           atom_list=atom_list,&
                           gth_potential=gth_potential)

      IF (.NOT.ASSOCIATED(gth_potential)) CYCLE
      CALL get_potential(potential=gth_potential,alpha_ppl=alpha,nexp_ppl=lppl,cexp_ppl=cexp_ppl)

      IF ( lppl <= 0 ) CYCLE

      ni = ncoset(2*lppl-2)
      ALLOCATE(pab(ni,1),STAT=ierr)
      CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
      pab = 0._dp

      nthread = 1
      ithread=0

      CALL reallocate ( cores, 1, natom )
      npme = 0
      cores = 0

      ! prepare core function
      DO j=1,lppl
        SELECT CASE (j)
          CASE (1)
            pab(1,1) = cexp_ppl(1)
          CASE (2)
            n = coset(2,0,0)
            pab(n,1) = cexp_ppl(2)
            n = coset(0,2,0)
            pab(n,1) = cexp_ppl(2)
            n = coset(0,0,2)
            pab(n,1) = cexp_ppl(2)
          CASE (3)
            n = coset(4,0,0)
            pab(n,1) = cexp_ppl(3)
            n = coset(0,4,0)
            pab(n,1) = cexp_ppl(3)
            n = coset(0,0,4)
            pab(n,1) = cexp_ppl(3)
            n = coset(2,2,0)
            pab(n,1) = 2._dp*cexp_ppl(3)
            n = coset(2,0,2)
            pab(n,1) = 2._dp*cexp_ppl(3)
            n = coset(0,2,2)
            pab(n,1) = 2._dp*cexp_ppl(3)
          CASE (4)
            n = coset(6,0,0)
            pab(n,1) = cexp_ppl(4)
            n = coset(0,6,0)
            pab(n,1) = cexp_ppl(4)
            n = coset(0,0,6)
            pab(n,1) = cexp_ppl(4)
            n = coset(4,2,0)
            pab(n,1) = 3._dp*cexp_ppl(4)
            n = coset(4,0,2)
            pab(n,1) = 3._dp*cexp_ppl(4)
            n = coset(2,4,0)
            pab(n,1) = 3._dp*cexp_ppl(4)
            n = coset(2,0,4)
            pab(n,1) = 3._dp*cexp_ppl(4)
            n = coset(0,4,2)
            pab(n,1) = 3._dp*cexp_ppl(4)
            n = coset(0,2,4)
            pab(n,1) = 3._dp*cexp_ppl(4)
            n = coset(2,2,2)
            pab(n,1) = 6._dp*cexp_ppl(4)
          CASE DEFAULT
            CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
        END SELECT
      END DO

      DO iatom = 1, natom
         atom_a = atom_list(iatom)
         ra(:) = pbc(particle_set(atom_a)%r,cell)
         IF(rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
            ! replicated realspace grid, split the atoms up between procs
            IF (MODULO(iatom,rs_rho%desc%group_size) == rs_rho % desc % my_pos ) THEN
               npme = npme + 1
               cores (npme) = iatom
            ENDIF
         ELSE
            npme = npme + 1
            cores (npme) = iatom
         ENDIF
      END DO

      IF(npme.GT.0)THEN
         DO j=1,npme

            iatom = cores(j)
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r,cell)
            subpatch_pattern=0
            ni = 2*lppl-2
            CALL collocate_pgf_product_rspace(ni,alpha,0,0,0.0_dp,0,ra,&
               (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,1.0_dp,pab,0,0,rs_rho,&
               cell,cube_info,eps_rho_rspace,ga_gb_function=FUNC_AB,&
               ithread=ithread,use_subpatch=.TRUE.,subpatch_pattern=subpatch_pattern,error=error)

         END DO
      ENDIF

      DEALLOCATE ( pab, STAT=ierr )
      CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)

    END DO

    IF (ASSOCIATED(cores)) THEN
      DEALLOCATE (cores,STAT=ierr)
      CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL rs_pw_transfer(rs_rho,vppl%pw,rs2pw,error=error)
    CALL rs_grid_release(rs_rho, error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_ppl_grid


! *****************************************************************************
!> \brief Collocates an arbitrary density from the aux_basis_set onto a grid.
!> \param rho The PW-Grid onto which the density is collocated
!> \param rho_g The realspace-grid onto which the density is collocated
!> \param rho_r The G-space-grid onto which the density is collocated
!> \param total_rho Gives back the integral of the collocated density
!> \param qs_env The QS environment of matter
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \par History
!>      08.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
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
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: mgrid_gspace, mgrid_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho

    NULLIFY(aux_basis_set, atomic_kind_set, atomic_kind, npgfa, cell, particle_set, &
            sphi_a, rs_rho, pw_env, cube_info, dft_control, rho_struct, atom_list, &
            first_sgfa, pab, work)
    NULLIFY(gridlevel_info, rs_descs, pw_pools, mgrid_gspace, mgrid_rspace)

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                    particle_set=particle_set, pw_env=pw_env, rho=rho_struct, &
                    dft_control=dft_control, error=error)

    cube_info => pw_env%cube_info
    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    ithread = 0
    gridlevel_info=>pw_env%gridlevel_info

    ! *** set up the pw multi-grids *** !
    CPPrecondition(ASSOCIATED(pw_env), cp_failure_level, routineP, error, failure)
    CALL pw_env_get(pw_env=pw_env, rs_descs=rs_descs, rs_grids=rs_rho, &
                    pw_pools=pw_pools, error=error)

    CALL pw_pools_create_pws(pw_pools, mgrid_rspace, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              error=error)

    CALL pw_pools_create_pws(pw_pools, mgrid_gspace, &
                              use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                              error=error)

    ! *** set up the rs multi-grids *** !
    DO igrid_level = 1,gridlevel_info%ngrid_levels
      CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid, error=error)
      CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
    END DO

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxco=maxco)
    ALLOCATE (pab(maxco,1),STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "pab",maxco*dp_size)
    ALLOCATE (work(maxco,1),STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "work",maxco*dp_size)

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
                                              map_consistent=map_consistent,use_subpatch=.FALSE.,error=error)
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE (pab,STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"pab")

    DEALLOCATE (work,STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"work")

    IF (gridlevel_info%ngrid_levels==1) THEN
      CALL rs_pw_transfer(rs=rs_rho(1)%rs_grid, pw=qs_env%rho%rho_r(1)%pw, dir=rs2pw,error=error)
      CALL rs_grid_release(rs_rho(1)%rs_grid, error=error)
      CALL pw_transfer(qs_env%rho%rho_r(1)%pw, qs_env%rho%rho_g(1)%pw,error=error)
      IF (qs_env%rho%rho_r(1)%pw%pw_grid%spherical) THEN
        CALL pw_transfer(qs_env%rho%rho_g(1)%pw, qs_env%rho%rho_r(1)%pw,error=error)
      END IF
    ELSE
      DO igrid_level = 1,gridlevel_info%ngrid_levels
         CALL rs_pw_transfer(rs=rs_rho(igrid_level)%rs_grid, &
                            pw=mgrid_rspace(igrid_level)%pw, dir=rs2pw, error=error)
         CALL rs_grid_release(rs_rho(igrid_level)%rs_grid, error=error)
      END DO

      CALL pw_zero(qs_env%rho%rho_g(1)%pw,error=error)
      DO igrid_level=1, gridlevel_info%ngrid_levels
        CALL pw_transfer(mgrid_rspace(igrid_level)%pw, &
                                   mgrid_gspace(igrid_level)%pw,error=error)
        CALL pw_axpy(mgrid_gspace(igrid_level)%pw, qs_env%rho%rho_g(1)%pw,&
             error=error)
      END DO
      CALL pw_transfer(qs_env%rho%rho_g(1)%pw, qs_env%rho%rho_r(1)%pw,error=error)
    END IF

    total_rho = pw_integrate_function(qs_env%rho%rho_r(1)%pw,isign=-1,&
                                      error=error)
    qs_env%rho%tot_rho_r = total_rho
    unit_nr=cp_logger_get_default_io_unit(logger)
    IF (unit_nr>0) THEN
      WRITE (unit_nr,*) "Total Rho =", total_rho
    END IF

    ! *** give back the multi-grids *** !
    CALL pw_pools_give_back_pws(pw_pools, mgrid_gspace, error=error)
    CALL pw_pools_give_back_pws(pw_pools, mgrid_rspace, error=error)

    qs_env%rho%rho_r_valid = .TRUE.
    qs_env%rho%rho_g_valid = .TRUE.

    CALL timestop(handle)

  END SUBROUTINE collocate_atomic_charge_density

! *****************************************************************************
!> \brief Collocates the fitted lri density on a grid.
!> \par History
!>      04.2013
!> \author Dorothea Golze
! *****************************************************************************
  SUBROUTINE calculate_lri_rho_elec(lri_rho_g, lri_rho_r, qs_env, lri_env,&
                                    lri_coef, total_rho, error)

    TYPE(pw_p_type), INTENT(INOUT)           :: lri_rho_g, lri_rho_r
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(lri_kind_type), DIMENSION(:), &
      POINTER                                :: lri_coef
    REAL(KIND=dp), INTENT(OUT)               :: total_rho
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_lri_rho_elec', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, first_sgf, handle, iatom, ierr, igrid_level, ikind, &
      ipgf, iset, j, m1, maxco, n, na1, natom, ncoa, nkind, npme, nseta, sgfa
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores, la_max, &
                                                la_min, npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: distributed_grids, failure, &
                                                map_consistent, parallel
    REAL(KIND=dp)                            :: eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:), POINTER     :: aci
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab, sphi_a, work, zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: lri_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: mgrid_gspace, mgrid_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho

    NULLIFY(aci, atomic_kind, atomic_kind_set, atom_list, cell, cores, cube_info,&
            dft_control, first_sgfa, gridlevel_info, la_max, la_min, lri_basis_set,&
            logger, mgrid_gspace, mgrid_rspace, npgfa, nsgfa, pab, para_env,&
            particle_set, pw_env, pw_pools, rs_descs, rs_rho, sphi_a, work, zeta)
            
    CALL timeset(routineN,handle)
    
    logger => cp_error_get_logger(error)

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,&
                    cell=cell, particle_set=particle_set, &
                    para_env=para_env, pw_env=pw_env,&
                    dft_control=dft_control, error=error)

    cube_info => pw_env%cube_info
    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    gridlevel_info=>pw_env%gridlevel_info

    ! *** set up the pw multi-grids *** !
    CPPrecondition(ASSOCIATED(pw_env), cp_failure_level, routineP, error, failure)
    CALL pw_env_get(pw_env=pw_env, rs_descs=rs_descs, rs_grids=rs_rho, &
                    pw_pools=pw_pools, error=error)

    CALL pw_pools_create_pws(pw_pools, mgrid_rspace, &
                             use_data=REALDATA3D,&
                             in_space=REALSPACE, &
                             error=error)

    CALL pw_pools_create_pws(pw_pools, mgrid_gspace, &
                             use_data=COMPLEXDATA1D,&
                             in_space=RECIPROCALSPACE,&
                             error=error)

    ! *** set up the rs multi-grids *** !
    DO igrid_level = 1,gridlevel_info%ngrid_levels
      CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid, error=error)
      CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
    END DO

    !take maxco from the LRI basis set!
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,basis_set_id=use_lri_basis_set)

    ALLOCATE (pab(maxco,1),STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "pab",maxco*dp_size)
    nkind = SIZE(atomic_kind_set)
    DO ikind = 1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind, lri_basis_set=lri_basis_set, natom=natom, &
                           atom_list=atom_list)
      
      !Take the lri basis set here!
      CALL get_gto_basis_set(gto_basis_set=lri_basis_set, lmax=la_max,&
                             lmin=la_min, zet=zeta, nset=nseta, npgf=npgfa,&
                             sphi=sphi_a, first_sgf=first_sgfa, nsgf_set=nsgfa)

      CALL reallocate (cores, 1, natom)
      npme = 0
      cores = 0
      
      distributed_grids = .FALSE.
      parallel = .FALSE.
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         IF (rs_descs(igrid_level)%rs_desc%distributed) distributed_grids = .TRUE.
         IF (rs_descs(igrid_level)%rs_desc%parallel) parallel = .TRUE.
      ENDDO

      DO iatom = 1, natom
         IF(parallel .AND. .NOT. distributed_grids) THEN
           ! replicated realspace grid, split the atoms up between procs
            IF (MODULO(iatom,para_env%num_pe) == para_env%mepos ) THEN
               npme = npme + 1
               cores (npme) = iatom
            ENDIF
         ELSE
              npme = npme + 1
              cores (npme) = iatom
         ENDIF
      END DO

      DO j = 1,npme
         iatom = cores(j)
         atom_a=atom_list(iatom)
         ra(:) = pbc(particle_set(atom_a)%r, cell)
         aci => lri_coef(ikind)%acoef(iatom,:) 

         DO iset = 1,nseta
           sgfa = first_sgfa(1,iset)
           ncoa = npgfa(iset)*ncoset(la_max(iset))
           m1=sgfa+nsgfa(iset)-1
           CALL reallocate (work,1,nsgfa(iset),1,1)
           work(1:nsgfa(iset),1)=aci(sgfa:m1)
           pab = 0._dp
 
           CALL dgemm("N","N",ncoa,1,nsgfa(iset),1.0_dp,lri_basis_set%sphi(1,sgfa),&
                SIZE(lri_basis_set%sphi,1),work(1,1),SIZE(work,1),0.0_dp,pab(1,1),SIZE(pab,1))
          
           DO ipgf = 1,npgfa(iset)
             na1 = (ipgf-1)*ncoset(la_max(iset))
             igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf,iset))

             CALL collocate_pgf_product_rspace(la_max=la_max(iset), zeta=zeta(ipgf,iset), &
                                               la_min=la_min(iset),&
                                               lb_max=0, zetb=0.0_dp, lb_min=0,&
                                               ra=ra,rab=(/0.0_dp,0.0_dp,0.0_dp/),rab2=0.0_dp,&
                                               scale=1._dp, pab=pab, o1=na1, o2=0,&     
                                               rsgrid=rs_rho(igrid_level)%rs_grid,cell=cell,&
                                               cube_info=cube_info(igrid_level),&
                                               eps_rho_rspace=eps_rho_rspace,&
                                               ga_gb_function=FUNC_AB, &
                                               map_consistent=map_consistent,&
                                               error=error)
           END DO
         END DO

         DEALLOCATE (work,STAT=ierr)
         IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"work")

      END DO
    END DO

    DEALLOCATE (pab,STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"pab")

    DEALLOCATE (cores,STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cores")
    

    CALL pw_zero(lri_rho_g%pw,error=error)
    CALL pw_zero(lri_rho_r%pw,error=error)

    DO igrid_level = 1,gridlevel_info%ngrid_levels
       CALL pw_zero(mgrid_rspace(igrid_level)%pw,error=error)
       CALL rs_pw_transfer(rs=rs_rho(igrid_level)%rs_grid, &
                          pw=mgrid_rspace(igrid_level)%pw, dir=rs2pw, error=error)
       CALL rs_grid_release(rs_rho(igrid_level)%rs_grid, error=error)
    END DO

    DO igrid_level=1, gridlevel_info%ngrid_levels
      CALL pw_zero(mgrid_gspace(igrid_level)%pw,error=error)
      CALL pw_transfer(mgrid_rspace(igrid_level)%pw, &
                                 mgrid_gspace(igrid_level)%pw,error=error)
      CALL pw_axpy(mgrid_gspace(igrid_level)%pw, lri_rho_g%pw,&
           error=error)
    END DO
    CALL pw_transfer(lri_rho_g%pw,lri_rho_r%pw,error=error)
    total_rho = pw_integrate_function(lri_rho_r%pw,isign=-1,error=error)

    ! *** give back the multi-grids *** !
    CALL pw_pools_give_back_pws(pw_pools, mgrid_gspace, error=error)
    CALL pw_pools_give_back_pws(pw_pools, mgrid_rspace, error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_lri_rho_elec  

! *****************************************************************************
!> \brief computes the density of the core charges on the grid
! *****************************************************************************
  SUBROUTINE calculate_rho_core(rho_core,total_rho,qs_env,only_nopaw,error)

    TYPE(pw_p_type), INTENT(INOUT)           :: rho_core
    REAL(KIND=dp), INTENT(OUT)               :: total_rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: only_nopaw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_core', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, ierr, &
                                                ikind, ithread, j, natom, &
                                                npme, nthread
    INTEGER(KIND=int_8)                      :: subpatch_pattern
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores
    LOGICAL                                  :: my_only_nopaw, paw_atom
    REAL(KIND=dp)                            :: alpha, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type)                     :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rhoc_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_rho

    CALL timeset(routineN,handle)
    NULLIFY(atomic_kind,cell,dft_control,pab,atomic_kind_set,particle_set,&
         atom_list,pw_env,rs_rho,auxbas_rs_desc,auxbas_pw_pool,cores)
    ALLOCATE (pab(1,1),STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "pab",dp_size)

    my_only_nopaw = .FALSE.
    IF(PRESENT(only_nopaw)) my_only_nopaw = only_nopaw

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    para_env=para_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env,auxbas_rs_desc=auxbas_rs_desc,auxbas_rs_grid=rs_rho,&
         auxbas_pw_pool=auxbas_pw_pool,error=error)
    cube_info=pw_env%cube_info(1)
    ! be careful in parallel nsmax is choosen with multigrid in mind!
    CALL rs_grid_retain(rs_rho, error=error)
    CALL rs_grid_zero(rs_rho)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           paw_atom=paw_atom,&
                           atom_list=atom_list,&
                           alpha_core_charge=alpha,&
                           ccore_charge=pab(1,1))

      IF(my_only_nopaw .AND. paw_atom ) CYCLE
      IF (alpha == 0.0_dp .OR. pab(1,1)== 0.0_dp) CYCLE

      nthread = 1
      ithread=0

      CALL reallocate ( cores, 1, natom )
      npme = 0
      cores = 0

      DO iatom = 1, natom
         atom_a = atom_list(iatom)
         ra(:) = pbc(particle_set(atom_a)%r,cell)
         IF(rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
            ! replicated realspace grid, split the atoms up between procs
            IF (MODULO(iatom,rs_rho%desc%group_size) == rs_rho % desc % my_pos ) THEN
               npme = npme + 1
               cores (npme) = iatom
            ENDIF
         ELSE
            npme = npme + 1
            cores (npme) = iatom
         ENDIF
      END DO

      IF(npme.GT.0)THEN
         DO j=1,npme

            iatom = cores(j)
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r,cell)
            subpatch_pattern=0
            CALL collocate_pgf_product_rspace(0,alpha,0,0,0.0_dp,0,ra,&
                  (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,-1.0_dp,pab,0,0,rs_rho,&
                  cell,cube_info,eps_rho_rspace,ga_gb_function=FUNC_AB,&
                  ithread=ithread,use_subpatch=.TRUE.,subpatch_pattern=subpatch_pattern,error=error)

         END DO
      ENDIF

    END DO

    IF (ASSOCIATED(cores)) THEN
      DEALLOCATE (cores,STAT=ierr)
      IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cores")
    END IF
    DEALLOCATE (pab,STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"pab")

    CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
         use_data=REALDATA3D,in_space=REALSPACE, error=error)

    CALL rs_pw_transfer(rs_rho,rhoc_r%pw,rs2pw,error=error)
    CALL rs_grid_release(rs_rho, error=error)

    total_rho = pw_integrate_function(rhoc_r%pw,isign=-1,error=error)

    CALL pw_transfer(rhoc_r%pw,rho_core%pw,error=error)

    CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw, error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_rho_core

! *****************************************************************************
!> \brief collocate a single gaussian on the grid
!> \par History
!>        12.2011 created
!> \author Dorothea Golze
! *****************************************************************************
  SUBROUTINE calculate_rho_single_gaussian(rho_gb,qs_env,iatom_in,error)

    TYPE(pw_p_type), INTENT(INOUT)           :: rho_gb
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom_in
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'calculate_rho_single_gaussian', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, ierr, &
                                                npme
    INTEGER(KIND=int_8)                      :: subpatch_pattern
    REAL(KIND=dp)                            :: eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rhoc_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_rho

    CALL timeset(routineN,handle)
    NULLIFY(cell,dft_control,pab,pw_env,rs_rho,auxbas_rs_desc,auxbas_pw_pool)
    ALLOCATE (pab(1,1),STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "pab",dp_size)

    CALL get_qs_env(qs_env=qs_env,&
                    cell=cell,&
                    dft_control=dft_control,&
                    para_env=para_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env,auxbas_rs_desc=auxbas_rs_desc,auxbas_rs_grid=rs_rho,&
         auxbas_pw_pool=auxbas_pw_pool,error=error)
    CALL rs_grid_retain(rs_rho, error=error)
    CALL rs_grid_zero(rs_rho)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    pab(1,1)=1.0_dp
    iatom=iatom_in

    npme = 0

    IF(rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
       IF (MODULO(iatom,rs_rho%desc%group_size) == rs_rho % desc % my_pos ) THEN
          npme = npme + 1
       ENDIF
    ELSE
       npme = npme + 1
    ENDIF

   IF(npme.GT.0)THEN
    atom_a = qs_env%qmmm_env_qm%image_charge_pot%image_mm_list(iatom)
    ra(:) = pbc(qs_env%qmmm_env_qm%image_charge_pot%particles_all(atom_a)%r,cell)
    subpatch_pattern=0
    CALL collocate_pgf_product_rspace(0,qs_env%qmmm_env_qm%image_charge_pot%eta,&
         0,0,0.0_dp,0,ra,(/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,1.0_dp,pab,0,0,rs_rho,&
         cell,pw_env%cube_info(1),eps_rho_rspace,ga_gb_function=FUNC_AB,&
         use_subpatch=.TRUE.,subpatch_pattern=subpatch_pattern,error=error)
   ENDIF

    DEALLOCATE (pab,STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"pab")

    CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
         use_data=REALDATA3D,in_space=REALSPACE, error=error)

    CALL rs_pw_transfer(rs_rho,rhoc_r%pw,rs2pw,error=error)
    CALL rs_grid_release(rs_rho, error=error)

    CALL pw_transfer(rhoc_r%pw,rho_gb%pw,error=error)

    CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw, error=error)

    CALL timestop(handle)

END SUBROUTINE calculate_rho_single_gaussian

! *****************************************************************************
!> \brief computes the image charge density on grid (including coeffcients)
!> \par History
!>        01.2012 created
!> \author Dorothea Golze
! *****************************************************************************
SUBROUTINE calculate_rho_metal(rho_metal,coeff,total_rho_metal,qs_env,error)

    TYPE(pw_p_type), INTENT(INOUT)           :: rho_metal
    REAL(KIND=dp), DIMENSION(:), POINTER     :: coeff
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: total_rho_metal
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_metal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, ierr, &
                                                j, natom, npme
    INTEGER(KIND=int_8)                      :: subpatch_pattern
    INTEGER, DIMENSION(:), POINTER           :: cores
    REAL(KIND=dp)                            :: eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rhoc_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_rho

  CALL timeset(routineN,handle)

  NULLIFY(cell,dft_control,pab,pw_env,rs_rho,auxbas_rs_desc,auxbas_pw_pool,&
          para_env,cores)

  ALLOCATE (pab(1,1),STAT=ierr)
  IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "pab",dp_size)

  CALL get_qs_env(qs_env=qs_env,&
                  cell=cell,&
                  dft_control=dft_control,&
                  para_env=para_env,pw_env=pw_env,error=error)
  CALL pw_env_get(pw_env,auxbas_rs_desc=auxbas_rs_desc,auxbas_rs_grid=rs_rho,&
       auxbas_pw_pool=auxbas_pw_pool,error=error)
  CALL rs_grid_retain(rs_rho, error=error)
  CALL rs_grid_zero(rs_rho)

  eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
  pab(1,1)=1.0_dp

  natom=SIZE(qs_env%qmmm_env_qm%image_charge_pot%image_mm_list)

  CALL reallocate ( cores, 1, natom )
   npme = 0
   cores = 0

  DO iatom=1, natom
   IF(rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
      IF (MODULO(iatom,rs_rho%desc%group_size) == rs_rho % desc % my_pos ) THEN
         npme = npme + 1
         cores (npme) = iatom
      ENDIF
   ELSE
      npme = npme + 1
      cores (npme) = iatom
   ENDIF
  ENDDO

  IF(npme.GT.0)THEN
   DO j=1,npme
    iatom = cores(j)
    atom_a = qs_env%qmmm_env_qm%image_charge_pot%image_mm_list(iatom)
    ra(:) = pbc(qs_env%qmmm_env_qm%image_charge_pot%particles_all(atom_a)%r,cell)
    subpatch_pattern=0
    CALL collocate_pgf_product_rspace(0,qs_env%qmmm_env_qm%image_charge_pot%eta,&
        0,0,0.0_dp,0,ra,(/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,coeff(iatom),pab,0,0,rs_rho,&
        cell,pw_env%cube_info(1),eps_rho_rspace,ga_gb_function=FUNC_AB,&
        use_subpatch=.TRUE.,subpatch_pattern=subpatch_pattern,error=error)
   ENDDO
  ENDIF

  DEALLOCATE (pab,cores,stat=ierr)
  IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"pab")

  CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
       use_data=REALDATA3D,in_space=REALSPACE, error=error)

  CALL rs_pw_transfer(rs_rho,rhoc_r%pw,rs2pw,error=error)
  CALL rs_grid_release(rs_rho, error=error)

  IF(PRESENT(total_rho_metal)) &
  !minus sign: account for the fact that rho_metal has opposite sign
  total_rho_metal = pw_integrate_function(rhoc_r%pw,isign=-1,error=error)

  CALL pw_transfer(rhoc_r%pw,rho_metal%pw,error=error)
  CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw, error=error)

  CALL timestop(handle)

END SUBROUTINE calculate_rho_metal
! *****************************************************************************
!> \brief collocate a single gaussian on the grid for periodic RESP fitting
!> \par History
!>        06.2012 created
!> \author Dorothea Golze
! *****************************************************************************
  SUBROUTINE calculate_rho_resp_single(rho_gb,qs_env,eta,iatom_in,error)
 
    TYPE(pw_p_type), INTENT(INOUT)           :: rho_gb
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(IN)                :: eta
    INTEGER, INTENT(IN)                      :: iatom_in
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_resp_single', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, ierr, npme
    INTEGER(KIND=int_8)                      :: subpatch_pattern
    REAL(KIND=dp)                            :: eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rhoc_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_rho

    CALL timeset(routineN,handle)
    NULLIFY(cell,dft_control,pab,pw_env,rs_rho,auxbas_rs_desc,auxbas_pw_pool,&
            particle_set)
    ALLOCATE (pab(1,1),STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "pab",dp_size)
 
    CALL get_qs_env(qs_env=qs_env,&
                    cell=cell,&
                    dft_control=dft_control,&
                    particle_set=particle_set,&
                    para_env=para_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env,auxbas_rs_desc=auxbas_rs_desc,auxbas_rs_grid=rs_rho,&
         auxbas_pw_pool=auxbas_pw_pool,error=error)
    CALL rs_grid_retain(rs_rho, error=error)
    CALL rs_grid_zero(rs_rho)
 
    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    pab(1,1)=1.0_dp
    iatom=iatom_in
 
    npme = 0
 
    IF(rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
       IF (MODULO(iatom,rs_rho%desc%group_size) == rs_rho % desc % my_pos ) THEN
          npme = npme + 1
       ENDIF
    ELSE
       npme = npme + 1
    ENDIF
 
   IF(npme.GT.0)THEN
    ra(:) = pbc(particle_set(iatom)%r,cell)
    subpatch_pattern=0
    CALL collocate_pgf_product_rspace(0,eta,0,0,0.0_dp,0,ra,&
         (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,1.0_dp,pab,0,0,rs_rho,&
         cell,pw_env%cube_info(1),eps_rho_rspace,ga_gb_function=FUNC_AB,&
         use_subpatch=.TRUE.,subpatch_pattern=subpatch_pattern,error=error)
   ENDIF
 
    DEALLOCATE (pab,STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"pab")
 
    CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
         use_data=REALDATA3D,in_space=REALSPACE, error=error)
 
    CALL rs_pw_transfer(rs_rho,rhoc_r%pw,rs2pw,error=error)
    CALL rs_grid_release(rs_rho, error=error)
 
    CALL pw_transfer(rhoc_r%pw,rho_gb%pw,error=error)
 
    CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw, error=error)

    CALL timestop(handle)
 
 END SUBROUTINE calculate_rho_resp_single
! *****************************************************************************
!> \brief computes the RESP charge density on grid based on the RESP charges
!> \par History
!>        01.2012 created
!> \author Dorothea Golze
! *****************************************************************************
SUBROUTINE calculate_rho_resp_all(rho_resp,coeff,natom,eta,qs_env,error)

    TYPE(pw_p_type), INTENT(INOUT)           :: rho_resp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: coeff
    INTEGER, INTENT(IN)                      :: natom
    REAL(KIND=dp), INTENT(IN)                :: eta
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_resp_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, ierr, j, npme
    INTEGER(KIND=int_8)                      :: subpatch_pattern
    INTEGER, DIMENSION(:), POINTER           :: cores
    REAL(KIND=dp)                            :: eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rhoc_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_rho

  CALL timeset(routineN,handle)

  NULLIFY(cell,dft_control,pab,pw_env,rs_rho,auxbas_rs_desc,auxbas_pw_pool,&
          para_env,particle_set)

  ALLOCATE (pab(1,1),STAT=ierr)
  IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "pab",dp_size)

  CALL get_qs_env(qs_env=qs_env,&
                  cell=cell,&
                  dft_control=dft_control,&
                  particle_set=particle_set,&
                  para_env=para_env,pw_env=pw_env,error=error)
  CALL pw_env_get(pw_env,auxbas_rs_desc=auxbas_rs_desc,auxbas_rs_grid=rs_rho,&
       auxbas_pw_pool=auxbas_pw_pool,error=error)
  CALL rs_grid_retain(rs_rho, error=error)
  CALL rs_grid_zero(rs_rho)

  eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
  pab(1,1)=1.0_dp

  CALL reallocate ( cores, 1, natom )
   npme = 0
   cores = 0

  DO iatom=1, natom
   IF(rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
      IF (MODULO(iatom,rs_rho%desc%group_size) == rs_rho % desc % my_pos ) THEN
         npme = npme + 1
         cores (npme) = iatom
      ENDIF
   ELSE
      npme = npme + 1
      cores (npme) = iatom
   ENDIF
  ENDDO

  IF(npme.GT.0)THEN
   DO j=1,npme
    iatom = cores(j)
    ra(:) = pbc(particle_set(iatom)%r,cell)
    subpatch_pattern=0
    CALL collocate_pgf_product_rspace(0,eta,&
        0,0,0.0_dp,0,ra,(/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,coeff(iatom),pab,0,0,rs_rho,&
        cell,pw_env%cube_info(1),eps_rho_rspace,ga_gb_function=FUNC_AB,&
        use_subpatch=.TRUE.,subpatch_pattern=subpatch_pattern,error=error)
   ENDDO
  ENDIF

  DEALLOCATE (pab,cores,stat=ierr)
  IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"pab")

  CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
       use_data=REALDATA3D,in_space=REALSPACE, error=error)

  CALL rs_pw_transfer(rs_rho,rhoc_r%pw,rs2pw,error=error)
  CALL rs_grid_release(rs_rho, error=error)

  CALL pw_transfer(rhoc_r%pw,rho_resp%pw,error=error)
  CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw, error=error)

  CALL timestop(handle)

END SUBROUTINE calculate_rho_resp_all

! *****************************************************************************
!> \brief computes the density corresponding to a given density matrix on the grid
!> \note
!>      both rho and rho_gspace contain the new rho
!>      (in real and g-space respectively)
!> \par History
!>      IAB (15-Feb-2010): Added OpenMP parallelisation to task loop
!>                         (c) The Numerical Algorithms Group (NAG) Ltd, 2010 on behalf of the HECToR project
! *****************************************************************************
  SUBROUTINE calculate_rho_elec(matrix_p, rho,rho_gspace, total_rho,&
       qs_env, soft_valid, compute_tau, compute_grad, basis_set_id, der_type, idir, &
       task_list_external, pw_env_external, error)

    TYPE(cp_dbcsr_type), POINTER             :: matrix_p
    TYPE(pw_p_type), INTENT(INOUT)           :: rho, rho_gspace
    REAL(KIND=dp), INTENT(OUT)               :: total_rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: soft_valid, compute_tau, &
                                                compute_grad
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id, der_type, idir
    TYPE(task_list_type), OPTIONAL, POINTER  :: task_list_external
    TYPE(pw_env_type), OPTIONAL, POINTER     :: pw_env_external
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_elec', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 2000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: bcol, brow, ga_gb_function, handle, i, iatom, iatom_old, &
      igrid_level, igrid_level_dummy, ikind, ikind_old, ipair, ipgf, iset, &
      iset_old, itask, ithread, jatom, jatom_old, jkind, jkind_old, jpgf, &
      jset, jset_old, lb, lb_l, lbr, lbw, maxco, maxpgf, maxset, maxsgf, &
      maxsgf_set, my_basis_set_id, my_der_type, my_idir, n, na1, na2, natoms, &
      nb1, nb2, nblock, ncoa, ncob, nkind, nn, nr, nrlevel, nseta, nsetb, &
      ntasks, nthread, nxy, nz, nzsize, segment, sgfa, sgfb, stat, ub
    INTEGER(kind=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_recv, atom_pair_send
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL :: atom_pair_changed, distributed_rs_grids, failure, found, &
      map_consistent, my_compute_grad, my_compute_tau, my_soft, use_subpatch
    REAL(KIND=dp)                            :: eps_rho_rspace, rab2, scale, &
                                                zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rab_inv, rb
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dist_ab, p_block, pab, &
                                                sphi_a, sphi_b, work, zeta, &
                                                zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pabt, workt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_type), POINTER             :: deltap
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(lgrid_type), POINTER                :: lgrid
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    TYPE(task_list_type), POINTER            :: task_list, task_list_soft

!$    INTEGER :: omp_get_thread_num, omp_get_max_threads

    failure=.FALSE.
    NULLIFY(atomic_kind, cell, dft_control, orb_basis_set, deltap, &
         atomic_kind_set, particle_set, rs_rho, pw_env, rs_descs, &
         para_env, dist_ab, la_max, la_min, &
         lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb, p_block, &
         sphi_a, sphi_b, zeta, zetb, first_sgfa, first_sgfb, tasks, pabt, &
         workt,lgrid)

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

    IF (PRESENT(compute_grad)) THEN
       my_compute_grad = compute_grad
    ELSE
       my_compute_grad = .FALSE.
    ENDIF


    CALL timeset(routineN,handle)
    my_der_type = 0
    my_idir = 0
    IF(PRESENT(der_type)) THEN
      my_der_type = der_type
      ga_gb_function = FUNC_DER(my_der_type)
    ELSE IF (my_compute_tau) THEN
       ga_gb_function = FUNC_DADB
    ELSE IF (my_compute_grad) THEN
        ga_gb_function = FUNC_DABpADB
    ELSE
       ga_gb_function = FUNC_AB
    ENDIF
    IF(PRESENT(idir)) my_idir = idir

    IF( PRESENT(basis_set_id) ) THEN
      my_basis_set_id = basis_set_id
    ELSE
      my_basis_set_id = use_orb_basis_set
    END IF


    SELECT CASE (my_basis_set_id)
    CASE (use_orb_basis_set)
      CALL get_qs_env(qs_env=qs_env,&
           atomic_kind_set=atomic_kind_set,&
           cell=cell,&
           dft_control=dft_control,&
           particle_set=particle_set,&
           para_env=para_env,&
           task_list=task_list,&
           task_list_soft=task_list_soft,&
           pw_env=pw_env,error=error)
    CASE (use_aux_fit_basis_set)
      CALL get_qs_env(qs_env=qs_env,&
           atomic_kind_set=atomic_kind_set,&
           cell=cell,&
           dft_control=dft_control,&
           particle_set=particle_set,&
           para_env=para_env,&
           task_list_aux_fit=task_list,&
           task_list_soft=task_list_soft,&
           pw_env=pw_env,error=error)
    END SELECT

    IF (PRESENT(pw_env_external)) pw_env => pw_env_external
    IF (PRESENT(task_list_external)) task_list=>task_list_external

    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info

    !   *** Allocate work storage ***
    nthread = 1
!$  nthread = omp_get_max_threads()
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxco=maxco,&
         maxsgf=maxsgf,&
         maxsgf_set=maxsgf_set,&
         basis_set_id=my_basis_set_id)
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
       SELECT CASE( my_basis_set_id)
       CASE (use_orb_basis_set)
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
              softb = my_soft, &
              orb_basis_set=orb_basis_set)
       CASE (use_aux_fit_basis_set)
           CALL get_atomic_kind(atomic_kind=atomic_kind,&
              softb = my_soft, &
              aux_fit_basis_set=orb_basis_set,&
              basis_set_id=my_basis_set_id)
       END SELECT

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            npgf=npgfa, nset=nseta )

       maxset=MAX(nseta,maxset)
       maxpgf=MAX(MAXVAL(npgfa),maxpgf)
    END DO

    ! get the task lists
    IF (my_soft) task_list=>task_list_soft
    CPPrecondition(ASSOCIATED(task_list),cp_failure_level,routineP,error,failure)
    tasks  =>task_list%tasks
    dist_ab=>task_list%dist_ab
    atom_pair_send=>task_list%atom_pair_send
    atom_pair_recv=>task_list%atom_pair_recv
    ntasks=task_list%ntasks

    ! *** set up the pw multi-grids
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_rho, lgrid=lgrid, error=error)
    ! *** set up the rs multi-grids
    CPPrecondition(ASSOCIATED(rs_rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rs_descs),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(lgrid),cp_failure_level,routineP,error,failure)
    distributed_rs_grids=.FALSE.

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid, error=error)
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
       IF ( rs_rho(igrid_level)%rs_grid%desc%distributed ) THEN
          distributed_rs_grids=.TRUE.
       ENDIF
    END DO

    IF ( nthread > 1 ) THEN
      IF (.NOT. ASSOCIATED(lgrid%r)) THEN
        CALL lgrid_allocate_grid(lgrid, nthread, error)
      ENDIF
    END IF

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    !   *** Initialize working density matrix ***
    ! distributed rs grids require a matrix that will be changed
    ! whereas this is not the case for replicated grids
    IF (distributed_rs_grids) THEN
       NULLIFY(deltap)
       !CALL replicate_matrix(matrix_p,deltap,target_name="DeltaP",error=error)
       ALLOCATE(deltap)
       CALL cp_dbcsr_init(deltap, error=error)
       CALL cp_dbcsr_copy(deltap,matrix_p,name="DeltaP",error=error)
       ! this matrix has no strict sparsity pattern in parallel
       !deltap%sparsity_id=-1
    ELSE
       deltap=>matrix_p
    ENDIF

    ! distribute the matrix
    IF (distributed_rs_grids) THEN
       CALL rs_distribute_matrix (rs_descs, deltap, atom_pair_send, atom_pair_recv, natoms, scatter=.TRUE., error=error)
    ENDIF

    ! map all tasks on the grids

!$omp parallel default(none), &
!$omp          shared(ntasks,tasks,natoms,maxset,maxpgf,particle_set,pabt,workt), &
!$omp          shared(my_basis_set_id,my_soft,deltap,maxco,dist_ab,ncoset,nthread), &
!$omp          shared(cell,cube_info,eps_rho_rspace,ga_gb_function, my_idir,map_consistent), &
!$omp          shared(rs_rho,lgrid,gridlevel_info,task_list), &
!$omp          private(igrid_level,iatom,jatom,iset,jset,ipgf,jpgf,ikind,jkind,pab,work), &
!$omp          private(iatom_old,jatom_old,iset_old,jset_old,ikind_old,jkind_old), &
!$omp          private(brow,bcol,atomic_kind,orb_basis_set,first_sgfa,la_max,la_min), &
!$omp          private(npgfa,nseta,nsgfa,sphi_a,zeta,first_sgfb,lb_max,lb_min,npgfb), &
!$omp          private(nsetb,nsgfb,sphi_b,zetb,p_block,found), &
!$omp          private(atom_pair_changed,ncoa,sgfa,ncob,sgfb,rab,rab2,ra,rb,zetp), &
!$omp          private(na1,na2,nb1,nb2,scale,use_subpatch,rab_inv,ithread,error,lb,ub,n), &
!$omp          private(lb_l,nn,segment,itask,i,nz,nxy,nzsize,nrlevel,nblock,lbw,lbr,nr,igrid_level_dummy)

    ithread = 0
!$  ithread = omp_get_thread_num()
    pab => pabt(:,:,ithread)
    work => workt(:,:,ithread)

    iatom_old = -1 ; jatom_old = -1 ; iset_old = -1 ; jset_old = -1
    ikind_old = -1 ; jkind_old = -1

    ! Loop over each gridlevel first, then loop and load balance over atom pairs
    ! We only need a single lgrid, which we sum back onto the rsgrid after each
    ! grid level is completed
    loop_gridlevels: DO igrid_level=1,gridlevel_info%ngrid_levels

    ! Only zero the region of the lgrid required for this grid level
    IF (nthread > 1) THEN
      lb = ithread*lgrid%ldim + 1
      ub = ithread*lgrid%ldim + rs_rho(igrid_level)%rs_grid%ngpts_local
      lgrid%r(lb:ub) = 0._dp
    END IF
!$omp do schedule(dynamic, MAX(1,task_list%npairs(igrid_level)/(nthread*50)))
    loop_pairs: DO ipair = 1, task_list%npairs(igrid_level)
    loop_tasks: DO itask = task_list%taskstart(ipair,igrid_level), task_list%taskstop(ipair,igrid_level)
       !decode the atom pair and basis info (igrid_level_dummy equals do loop variable by construction).
       CALL int2pair(tasks(3,itask),igrid_level_dummy,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
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
             SELECT CASE (my_basis_set_id)
               CASE (use_orb_basis_set)
                 CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                    softb = my_soft, &
                    orb_basis_set=orb_basis_set)
               CASE (use_aux_fit_basis_set)
                 CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                    softb = my_soft, &
                    aux_fit_basis_set=orb_basis_set,&
                    basis_set_id=my_basis_set_id)
             END SELECT
             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  first_sgf=first_sgfa,&
                  lmax=la_max,&
                  lmin=la_min,&
                  npgf=npgfa,&
                  nset=nseta,&
                  nsgf_set=nsgfa,&
                  sphi=sphi_a,&
                  zet=zeta)
          ENDIF

          IF (jkind .NE. jkind_old ) THEN
             SELECT CASE (my_basis_set_id)
             CASE (use_orb_basis_set)
               CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind,&
                    softb = my_soft, &
                    orb_basis_set=orb_basis_set)
             CASE (use_aux_fit_basis_set)
               CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind,&
                    softb = my_soft, &
                    aux_fit_basis_set=orb_basis_set,&
                    basis_set_id=my_basis_set_id)
             END SELECT
             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  first_sgf=first_sgfb,&
                  lmax=lb_max,&
                  lmin=lb_min,&
                  npgf=npgfb,&
                  nset=nsetb,&
                  nsgf_set=nsgfb,&
                  sphi=sphi_b,&
                  zet=zetb)
          ENDIF

          CALL cp_dbcsr_get_block_p(matrix=deltap,&
               row=brow,col=bcol,BLOCK=p_block,found=found)
          IF (.NOT.ASSOCIATED(p_block)) & 
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "p_block not associated in deltap")

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

       ! check whether we need to use fawzi's generalised collocation scheme
       IF(rs_rho(igrid_level)%rs_grid%desc%distributed)THEN
          !tasks(4,:) is 0 for replicated, 1 for distributed 2 for exceptional distributed tasks
          IF (tasks(4,itask) .EQ. 2 ) THEN
             use_subpatch = .TRUE.
          ELSE
             use_subpatch = .FALSE.
          ENDIF
       ELSE
          use_subpatch = .FALSE.
       ENDIF

       IF (nthread > 1) THEN
          IF (iatom <= jatom) THEN
             CALL collocate_pgf_product_rspace(&
                 la_max(iset),zeta(ipgf,iset),la_min(iset),&
                 lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                 ra,rab,rab2,scale,pab,na1-1,nb1-1,&
                 rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
                 eps_rho_rspace,&
                 ga_gb_function=ga_gb_function, &
                 idir=my_idir,&
                 lgrid=lgrid,ithread=ithread, &
                 map_consistent=map_consistent,use_subpatch=use_subpatch,&
                 subpatch_pattern=tasks(6,itask),error=error)
          ELSE
             rab_inv=-rab
             CALL collocate_pgf_product_rspace(&
                 lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                 la_max(iset),zeta(ipgf,iset),la_min(iset),&
                 rb,rab_inv,rab2,scale,pab,nb1-1,na1-1,&
                 rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
                 eps_rho_rspace,&
                 ga_gb_function=ga_gb_function, &
                 idir=my_idir,&
                 lgrid=lgrid,ithread=ithread, &
                 map_consistent=map_consistent,use_subpatch=use_subpatch,&
                 subpatch_pattern=tasks(6,itask),error=error)
          END IF
       ELSE
          IF (iatom <= jatom) THEN
             CALL collocate_pgf_product_rspace(&
                 la_max(iset),zeta(ipgf,iset),la_min(iset),&
                 lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                 ra,rab,rab2,scale,pab,na1-1,nb1-1,&
                 rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
                 eps_rho_rspace,&
                 ga_gb_function=ga_gb_function, &
                 idir=my_idir,&
                 map_consistent=map_consistent,use_subpatch=use_subpatch,&
                 subpatch_pattern=tasks(6,itask),error=error)
          ELSE
             rab_inv=-rab
             CALL collocate_pgf_product_rspace(&
                 lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                 la_max(iset),zeta(ipgf,iset),la_min(iset),&
                 rb,rab_inv,rab2,scale,pab,nb1-1,na1-1,&
                 rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
                 eps_rho_rspace,&
                 ga_gb_function=ga_gb_function, &
                 idir=my_idir,&
                 map_consistent=map_consistent,use_subpatch=use_subpatch,&
                 subpatch_pattern=tasks(6,itask),error=error)
          END IF
       END IF
    END DO loop_tasks
    END DO loop_pairs
!$omp end do

    ! Now sum the thread-local grids back into the rs_grid (in parallel, each thread writes to a section of the rs_grid at a time)
    IF (nthread > 1) THEN
        nz = (1 + rs_rho(igrid_level)%rs_grid%ub_local(3) &
                - rs_rho(igrid_level)%rs_grid%lb_local(3))
        nxy = (1 + rs_rho(igrid_level)%rs_grid%ub_local(1) &
                 - rs_rho(igrid_level)%rs_grid%lb_local(1)) &
             *(1 + rs_rho(igrid_level)%rs_grid%ub_local(2) &
                 - rs_rho(igrid_level)%rs_grid%lb_local(2))

        ! Work out the number of tree levels, and start the reduction

        nrlevel = CEILING(LOG10(1.0_dp*nthread)/LOG10(2.0_dp))

        DO nr = 1, nrlevel
           nblock = 2**nr
           nzsize = MIN(nblock, nthread - (ithread/nblock)*nblock)
           itask = MODULO(ithread, nblock)

           lb = (nz*itask)/nzsize
           ub = (nz*(itask+1))/nzsize - 1
           n = (ithread/nblock)*nblock
           lbw = 1 + n*lgrid%ldim + nxy*lb

           n = n + nblock/2
           IF (n < nthread) THEN
              lbr = 1 + n*lgrid%ldim + nxy*lb

              CALL daxpy(nxy*(1+ub-lb), 1.0_dp, lgrid%r(lbr), 1, lgrid%r(lbw), 1)
           END IF
!$omp barrier
        END DO

        ! Copy final result from first local grid to rs_grid
        lb = (nz*ithread)/nthread
        ub = (nz*(ithread+1))/nthread - 1
        nzsize = 1 + (ub - lb)

        lb = lb + rs_rho(igrid_level)%rs_grid%lb_local(3)
        ub = ub + rs_rho(igrid_level)%rs_grid%lb_local(3)
        lbr = 1 + nxy*(nz*ithread/nthread)

        CALL daxpy(nxy*nzsize, 1.0_dp, lgrid%r(lbr),1, &
                   rs_rho(igrid_level)%rs_grid%r(:,:,lb:ub), 1)
!$omp barrier
      END IF
    END DO loop_gridlevels
!$omp end parallel

    !   *** Release work storage ***

    IF (distributed_rs_grids) THEN
       CALL cp_dbcsr_deallocate_matrix ( deltap ,error=error)
    ENDIF

    DEALLOCATE (pabt,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"pabt")
    DEALLOCATE (workt,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"workt")

    CALL density_rs2pw(pw_env,rs_rho,rho,rho_gspace,error=error)

    total_rho = pw_integrate_function(rho%pw,isign=-1,error=error)
    CALL timestop(handle)
  END SUBROUTINE calculate_rho_elec

! *****************************************************************************
!> \brief computes the gradient of the density corresponding to a given
!>        density matrix on the grid
!> \note  this is an alternative to calculate the gradient through FFTs
! *****************************************************************************
  SUBROUTINE calculate_drho_elec(matrix_p,drho,drho_gspace,qs_env,&
             soft_valid,basis_set_id,error)

    TYPE(cp_dbcsr_type), POINTER             :: matrix_p
    TYPE(pw_p_type), DIMENSION(:), &
      INTENT(INOUT)                          :: drho, drho_gspace
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: soft_valid
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_drho_elec', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 2000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: bcol, brow, handle, i, iatom, iatom_old, idir, igrid_level, &
      ikind, ikind_old, ipgf, iset, iset_old, itask, ithread, jatom, &
      jatom_old, jkind, jkind_old, jpgf, jset, jset_old, maxco, maxpgf, &
      maxset, maxsgf, maxsgf_set, my_basis_set_id, na1, na2, natoms, nb1, &
      nb2, ncoa, ncob, nkind, nseta, nsetb, ntasks, nthread, sgfa, sgfb, stat
    INTEGER(kind=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_recv, atom_pair_send
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL :: atom_pair_changed, distributed_rs_grids, failure, found, &
      map_consistent, my_soft, use_subpatch
    REAL(KIND=dp)                            :: eps_rho_rspace, rab2, scale, &
                                                zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rab_inv, rb
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dist_ab, p_block, pab, &
                                                sphi_a, sphi_b, work, zeta, &
                                                zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pabt, workt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_type), POINTER             :: deltap
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_rho
    TYPE(task_list_type), POINTER            :: task_list, task_list_soft

    failure=.FALSE.
    NULLIFY(atomic_kind, cell, dft_control, orb_basis_set, deltap, &
         atomic_kind_set, sab_orb, particle_set, rs_rho, pw_env, rs_descs, &
         para_env, dist_ab, la_max, la_min, &
         lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb, p_block, &
         sphi_a, sphi_b, zeta, zetb, first_sgfa, first_sgfb, tasks, pabt, &
         workt)

    debug_count=debug_count+1

    ! by default, the full density is calculated
    my_soft=.FALSE.
    IF (PRESENT(soft_valid)) my_soft = soft_valid

    CALL timeset(routineN,handle)
    !
    IF( PRESENT(basis_set_id) ) THEN
      my_basis_set_id = basis_set_id
    ELSE
      my_basis_set_id = use_orb_basis_set
    END IF

    SELECT CASE (my_basis_set_id)
    CASE (use_orb_basis_set)
      CALL get_qs_env(qs_env=qs_env,&
           atomic_kind_set=atomic_kind_set,&
           cell=cell,&
           dft_control=dft_control,&
           particle_set=particle_set,&
           sab_orb=sab_orb,&
           para_env=para_env,&
           task_list=task_list,&
           task_list_soft=task_list_soft,&
           pw_env=pw_env,error=error)
    CASE (use_aux_fit_basis_set)
      CALL get_qs_env(qs_env=qs_env,&
           atomic_kind_set=atomic_kind_set,&
           cell=cell,&
           dft_control=dft_control,&
           particle_set=particle_set,&
           sab_aux_fit=sab_orb,&
           para_env=para_env,&
           task_list_aux_fit=task_list,&
           task_list_soft=task_list_soft,&
           pw_env=pw_env,error=error)
    END SELECT

    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info

    !   *** Allocate work storage ***
    nthread = 1
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxco=maxco,&
         maxsgf=maxsgf,&
         maxsgf_set=maxsgf_set,&
         basis_set_id=my_basis_set_id)
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
       SELECT CASE( my_basis_set_id)
       CASE (use_orb_basis_set)
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
              softb = my_soft, &
              orb_basis_set=orb_basis_set)
       CASE (use_aux_fit_basis_set)
           CALL get_atomic_kind(atomic_kind=atomic_kind,&
              softb = my_soft, &
              aux_fit_basis_set=orb_basis_set,&
              basis_set_id=my_basis_set_id)
       END SELECT

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            npgf=npgfa, nset=nseta )

       maxset=MAX(nseta,maxset)
       maxpgf=MAX(MAXVAL(npgfa),maxpgf)
    END DO

    ! get the task lists
    IF (my_soft) task_list=>task_list_soft
    CPPrecondition(ASSOCIATED(task_list),cp_failure_level,routineP,error,failure)
    tasks  =>task_list%tasks
    dist_ab=>task_list%dist_ab
    atom_pair_send=>task_list%atom_pair_send
    atom_pair_recv=>task_list%atom_pair_recv
    ntasks=task_list%ntasks

    ! *** set up the rs multi-grids
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_rho, error=error)
    DO igrid_level=1,gridlevel_info%ngrid_levels
      CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid, error=error)
      distributed_rs_grids=rs_rho(igrid_level)%rs_grid%desc%distributed
    END DO

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent = dft_control%qs_control%map_consistent
    !   *** Initialize working density matrix ***
    ! distributed rs grids require a matrix that will be changed
    ! whereas this is not the case for replicated grids
    IF (distributed_rs_grids) THEN
       NULLIFY(deltap)
       ALLOCATE(deltap)
       CALL cp_dbcsr_init(deltap, error=error)
       CALL cp_dbcsr_copy(deltap,matrix_p,name="DeltaP",error=error)
       !CALL replicate_matrix(matrix_p,deltap,target_name="DeltaP",error=error)
       ! this matrix has no strict sparsity pattern in parallel
       !deltap%sparsity_id=-1
    ELSE
       deltap=>matrix_p
    ENDIF

    ! distribute the matrix
    IF (distributed_rs_grids) THEN
       CALL rs_distribute_matrix (rs_descs, deltap, atom_pair_send, atom_pair_recv, &
                               natoms, scatter=.TRUE., error=error)
    ENDIF

    ! map all tasks on the grids

    ithread = 0
    pab => pabt(:,:,ithread)
    work => workt(:,:,ithread)

    loop_xyz: DO idir=1,3

      DO igrid_level=1,gridlevel_info%ngrid_levels
         CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
      END DO

      iatom_old = -1 ; jatom_old = -1 ; iset_old = -1 ; jset_old = -1
      ikind_old = -1 ; jkind_old = -1
      loop_tasks: DO itask=1,ntasks

        !decode the atom pair and basis info
        CALL int2pair(tasks(3,itask),igrid_level,iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
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
             SELECT CASE (my_basis_set_id)
               CASE (use_orb_basis_set)
                 CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                    softb = my_soft, &
                    orb_basis_set=orb_basis_set)
               CASE (use_aux_fit_basis_set)
                 CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                    softb = my_soft, &
                    aux_fit_basis_set=orb_basis_set,&
                    basis_set_id=my_basis_set_id)
             END SELECT
             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  first_sgf=first_sgfa,&
                  lmax=la_max,&
                  lmin=la_min,&
                  npgf=npgfa,&
                  nset=nseta,&
                  nsgf_set=nsgfa,&
                  sphi=sphi_a,&
                  zet=zeta)
          ENDIF

          IF (jkind .NE. jkind_old ) THEN
             SELECT CASE (my_basis_set_id)
             CASE (use_orb_basis_set)
               CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind,&
                    softb = my_soft, &
                    orb_basis_set=orb_basis_set)
             CASE (use_aux_fit_basis_set)
               CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind,&
                    softb = my_soft, &
                    aux_fit_basis_set=orb_basis_set,&
                    basis_set_id=my_basis_set_id)
             END SELECT
             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  first_sgf=first_sgfb,&
                  lmax=lb_max,&
                  lmin=lb_min,&
                  npgf=npgfb,&
                  nset=nsetb,&
                  nsgf_set=nsgfb,&
                  sphi=sphi_b,&
                  zet=zetb)
          ENDIF

          CALL cp_dbcsr_get_block_p(matrix=deltap,&
               row=brow,col=bcol,BLOCK=p_block,found=found)
          IF (.NOT.ASSOCIATED(p_block)) &
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "p_block not associated in deltap")

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

        ! check whether we need to use fawzi's generalised collocation scheme
        IF(rs_rho(igrid_level)%rs_grid%desc%distributed)THEN
          !tasks(4,:) is 0 for replicated, 1 for distributed 2 for exceptional distributed tasks
          IF (tasks(4,itask) .EQ. 2 ) THEN
             use_subpatch = .TRUE.
          ELSE
             use_subpatch = .FALSE.
          ENDIF
        ELSE
          use_subpatch = .FALSE.
        ENDIF
        IF (iatom <= jatom) THEN
          CALL collocate_pgf_product_rspace(&
               la_max(iset),zeta(ipgf,iset),la_min(iset),&
               lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
               ra,rab,rab2,scale,pab,na1-1,nb1-1,&
               rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
               eps_rho_rspace,&
               ga_gb_function=FUNC_DABpADB, &
               idir=idir, &
               map_consistent=map_consistent,use_subpatch=use_subpatch,subpatch_pattern=tasks(6,itask),error=error)
        ELSE
          rab_inv=-rab
          CALL collocate_pgf_product_rspace(&
               lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
               la_max(iset),zeta(ipgf,iset),la_min(iset),&
               rb,rab_inv,rab2,scale,pab,nb1-1,na1-1,&
               rs_rho(igrid_level)%rs_grid,cell,cube_info(igrid_level),&
               eps_rho_rspace,&
               ga_gb_function=FUNC_DABpADB, &
               idir=idir, &
               map_consistent=map_consistent,use_subpatch=use_subpatch,subpatch_pattern=tasks(6,itask),error=error)
        END IF

      END DO loop_tasks

      CALL density_rs2pw_basic(pw_env,rs_rho,drho(idir),drho_gspace(idir),error=error)

    END DO loop_xyz

    !   *** Release work storage ***
    IF (ASSOCIATED(rs_rho)) THEN
      DO i=1, SIZE(rs_rho)
        CALL rs_grid_release(rs_rho(i)%rs_grid, error=error)
      END DO
    END IF
    IF (distributed_rs_grids) THEN
      CALL cp_dbcsr_deallocate_matrix ( deltap ,error=error)
    ENDIF

    DEALLOCATE (pabt,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"pabt")
    DEALLOCATE (workt,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"workt")

    CALL timestop(handle)

  END SUBROUTINE calculate_drho_elec

! *****************************************************************************
!> \brief maps a given wavefunction on the grid
!> \note
!>      modified calculate_rho_elec, should write the wavefunction represented by
!>      it's presumably dominated by the FFT and the rs->pw and back routines
!>
!>    BUGS ???
!>      does it take correctly the periodic images of the basis functions into account
!>      i.e. is only correct if the basis functions are localised enough to be just in 1 cell ?
!> \par History
!>      08.2002 created [Joost VandeVondele]
!>      03.2006 made independent of qs_env [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE calculate_wavefunction(mo_vectors,ivector,rho,rho_gspace, &
                   atomic_kind_set,cell,dft_control,particle_set, &
                   pw_env,basis_set_id, external_vector, error)

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
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: external_vector
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_wavefunction', &
      routineP = moduleN//':'//routineN

    INTEGER :: dir, group, group_size, handle, i, iatom, igrid_level, ipgf, &
      iset, maxco, maxsgf_set, my_basis_set_id, my_pos, na1, na2, nao, natom, &
      ncoa, nseta, offset, sgfa, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: where_is_the_point
    INTEGER, DIMENSION(3)                    :: lb, location, tp, ub
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: failure, local, map_it_here
    REAL(KIND=dp)                            :: dab, eps_rho_rspace, rab2, &
                                                scale, zetp
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvector
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab, sphi_a, work, zeta
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: mgrid_gspace, mgrid_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(REALSPACE_GRID_P_TYPE), &
      DIMENSION(:), POINTER                  :: rs_rho

    IF( PRESENT( basis_set_id ) ) THEN
      my_basis_set_id = basis_set_id
    ELSE
      my_basis_set_id = use_orb_basis_set
    END IF

    failure=.FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(eigenvector,  orb_basis_set,&
         pab,work,la_max, la_min,&
         npgfa, nsgfa, &
         sphi_a, zeta, first_sgfa,&
         rs_rho,rs_descs,pw_pools,mgrid_rspace,mgrid_gspace)

    IF (PRESENT(external_vector)) THEN
       nao=SIZE(external_vector)
       ALLOCATE (eigenvector(nao),STAT=stat)
       IF (stat.NE.0) CALL stop_memory(routineN,moduleN,__LINE__,"eigenvector",dp_size*nao)
       eigenvector=external_vector
    ELSE
       CALL cp_fm_get_info(matrix=mo_vectors,nrow_global=nao, error=error)
       ALLOCATE (eigenvector(nao),STAT=stat)
       IF (stat.NE.0) CALL stop_memory(routineN,moduleN,__LINE__,"eigenvector",dp_size*nao)
       DO i=1,nao
          CALL cp_fm_get_element(mo_vectors,i,ivector,eigenvector(i),local)
       ENDDO
    ENDIF

    ! *** set up the pw multi-grids
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_rho, pw_pools=pw_pools, &
        cube_info=cube_info, gridlevel_info=gridlevel_info, error=error)

    CALL pw_pools_create_pws(pw_pools,mgrid_gspace,&
                use_data = COMPLEXDATA1D,&
                in_space = RECIPROCALSPACE, error=error)
    CALL pw_pools_create_pws(pw_pools,mgrid_rspace,&
                use_data = REALDATA3D,&
                in_space = REALSPACE, error=error)

    ! *** set up rs multi-grids
    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid, error=error)
       CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
    END DO


    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
!   *** Allocate work storage ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             natom=natom,&
                             maxsgf_set=maxsgf_set,&
                             basis_set_id=my_basis_set_id)

    ALLOCATE (pab(maxco,1),STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "pab",maxco*dp_size)
    ALLOCATE (work(maxco,1),STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "work",maxco*dp_size)

    offset=0
    group=mgrid_rspace(1)%pw%pw_grid%para%group
    my_pos=mgrid_rspace(1)%pw%pw_grid%para%my_pos
    group_size=mgrid_rspace(1)%pw%pw_grid%para%group_size
    ALLOCATE(where_is_the_point(0:group_size-1))

    DO iatom=1,natom
      SELECT CASE (my_basis_set_id)
      CASE (use_orb_basis_set)
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             orb_basis_set=orb_basis_set)
      CASE (use_aux_fit_basis_set)
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             aux_fit_basis_set=orb_basis_set)
      CASE (use_aux_basis_set)
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             aux_basis_set=orb_basis_set)
      CASE (use_ri_aux_basis_set)
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             ri_aux_basis_set=orb_basis_set)
      END SELECT
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

            IF (.NOT. ALL (rs_rho(igrid_level)%rs_grid%desc%perd == 1)) THEN
               DO dir = 1,3
                     ! bounds of local grid (i.e. removing the 'wings'), if periodic
                     tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),ra)*rs_rho(igrid_level)%rs_grid%desc%npts(dir))
                     tp(dir) = MODULO ( tp(dir), rs_rho(igrid_level)%rs_grid%desc%npts(dir) )
                     IF (rs_rho(igrid_level)%rs_grid%desc%perd(dir) .NE. 1) THEN
                        lb(dir) = rs_rho(igrid_level)%rs_grid%lb_local ( dir ) + rs_rho(igrid_level)%rs_grid%desc%border
                        ub(dir) = rs_rho(igrid_level)%rs_grid%ub_local ( dir ) - rs_rho(igrid_level)%rs_grid%desc%border
                     ELSE
                        lb(dir) = rs_rho(igrid_level)%rs_grid%lb_local ( dir )
                        ub(dir) = rs_rho(igrid_level)%rs_grid%ub_local ( dir )
                     ENDIF
                     ! distributed grid, only map if it is local to the grid
                     location(dir)=tp(dir)+rs_rho(igrid_level)%rs_grid%desc%lb(dir)
               ENDDO
               IF  (lb(1)<=location(1) .AND. location(1)<=ub(1) .AND. &
                    lb(2)<=location(2) .AND. location(2)<=ub(2) .AND. &
                    lb(3)<=location(3) .AND. location(3)<=ub(3)) THEN
                  map_it_here=.TRUE.
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
                 eps_rho_rspace,map_consistent=.TRUE.,ga_gb_function=FUNC_AB,error=error)

         END DO

         offset=offset+nsgfa(iset)

      END DO

    END DO

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid,&
            mgrid_rspace(igrid_level)%pw,rs2pw,error=error)
       CALL rs_grid_release(rs_rho(igrid_level)%rs_grid, error=error)
    ENDDO

    CALL pw_zero(rho_gspace%pw,error=error)
    DO igrid_level=1,gridlevel_info%ngrid_levels
      CALL pw_transfer(mgrid_rspace(igrid_level)%pw,&
           mgrid_gspace(igrid_level)%pw,error=error)
      CALL pw_axpy(mgrid_gspace(igrid_level)%pw,rho_gspace%pw,&
           error=error)
    END DO

    CALL pw_transfer(rho_gspace%pw,rho%pw,error=error)

    ! Release work storage
    DEALLOCATE (eigenvector,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"eigenvector")

    DEALLOCATE (pab,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"pab")

    DEALLOCATE (work,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"work")

    ! give back the pw multi-grids
    CALL pw_pools_give_back_pws(pw_pools,mgrid_gspace,error=error)
    CALL pw_pools_give_back_pws(pw_pools,mgrid_rspace,error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_wavefunction

! *****************************************************************************
!> \brief low level collocation of primitive gaussian functions
! *****************************************************************************
  SUBROUTINE collocate_pgf_product_rspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2,scale,pab,o1,o2,&
                                          rsgrid,cell,cube_info,&
                                          eps_rho_rspace,ga_gb_function,&
                                          lgrid,ithread,&
                                          map_consistent,&
                                          collocate_rho0,&
                                          rpgf0_s,idir,ir,rsgauge,rsbuf,&
                                          use_subpatch,subpatch_pattern,error)

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
      TYPE(realspace_grid_type), POINTER, OPTIONAL :: rsgauge,rsbuf
      LOGICAL, OPTIONAL                        :: use_subpatch
      INTEGER(KIND=int_8), OPTIONAL, INTENT(IN):: subpatch_pattern
      TYPE(cp_error_type), INTENT(INOUT)       :: error

      CHARACTER(len=*), PARAMETER :: routineN = 'collocate_pgf_product_rspace', &
        routineP = moduleN//':'//routineN

      INTEGER :: cmax, gridbounds(2,3), i, ico, icoef, ider1, ider2, ig, ithread_l, &
        jco, k, l, la_max_local, la_min_local, lb_max_local, lb_min_local, &
        length, lx, lx_max, lxa, lxb, lxy, lxy_max, lxyz, lxyz_max, lya, lyb, &
        lza, lzb, o1_local, o2_local, offset, start
      INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, ng, &
                                                  ub_cube
      INTEGER, DIMENSION(:), POINTER           :: ly_max, lz_max, sphere_bounds
      LOGICAL                                  :: failure, my_collocate_rho0, &
                                                  my_map_consistent,subpatch_collocate
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

      REAL(kind=dp), POINTER, DIMENSION(:,:,:) :: grid_tmp,gauge

    failure = .FALSE.

    subpatch_collocate = .FALSE.

    IF(PRESENT(use_subpatch)) THEN
       IF(use_subpatch)THEN
         subpatch_collocate = .TRUE.
         CPPrecondition(PRESENT(subpatch_pattern),cp_failure_level,routineP,error,failure)
       ENDIF
    ENDIF

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
    CASE(FUNC_DABpADB)
        CPPrecondition(PRESENT(idir),cp_failure_level,routineP,error,failure)
        la_max_local=la_max+1
        la_min_local=MAX(la_min-1,0)
        lb_max_local=lb_max+1
        lb_min_local=MAX(lb_min-1,0)
        ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
        ! is equivalent to mapping pab with
        !    pgf_a (nabla_{idir} pgf_b) + (nabla_{idir} pgf_a) pgf_b
        ! ( pgf_a ) (ddx pgf_b) + (ddx pgf_a)( pgf_b ) =
        !          pgf_a *(lbx pgf_{b-1x} + 2*zetb*pgf_{b+1x}) +
        !                   (lax pgf_{a-1x} + 2*zeta*pgf_{a+1x}) pgf_b

        ALLOCATE(pab_local(ncoset(la_max_local),ncoset(lb_max_local)))
        pab_local = 0.0_dp
        DO lxa=0,la_max
        DO lxb=0,lb_max
           DO lya=0,la_max-lxa
           DO lyb=0,lb_max-lxb
              DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
              DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
                 ! this element of pab results in 4 elements of pab_local
                 CALL prepare_dab_p_adb(pab_local,pab,idir,&
                      lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)
              END DO
              END DO
           END DO
           END DO
        END DO
        END DO
        o1_local=0
        o2_local=0
    CASE(FUNC_DX,FUNC_DY,FUNC_DZ)
        ider1 = ga_gb_function - 500
        la_max_local=la_max+1
        la_min_local=MAX(la_min-1,0)
        lb_max_local=lb_max+1
        lb_min_local=MAX(lb_min-1,0)
        ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
        ! is equivalent to mapping pab with
        !   d_{ider1} pgf_a d_{ider1} pgf_b
        ! dx pgf_a dx pgf_b =
        !        (lax pgf_{a-1x})*(lbx pgf_{b-1x}) - 2*zetb*lax*pgf_{a-1x}*pgf_{b+1x} -
        !         lbx pgf_{b-1x}*2*zeta*pgf_{a+1x}+ 4*zeta*zetab*pgf_{a+1x}pgf_{b+1x}

        ALLOCATE(pab_local(ncoset(la_max_local),ncoset(lb_max_local)))
        pab_local = 0.0_dp
        DO lxa=0,la_max
        DO lxb=0,lb_max
           DO lya=0,la_max-lxa
           DO lyb=0,lb_max-lxb
              DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
              DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
                 ! this element of pab results in 4 elements of pab_local
                 CALL prepare_dIadIb(pab_local,pab,ider1,&
                      lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)
              END DO
              END DO
           END DO
           END DO
        END DO
        END DO
        o1_local=0
        o2_local=0
    CASE(FUNC_DXDY,FUNC_DYDZ,FUNC_DZDX)
        ider1 = ga_gb_function - 600
        ider2 = ga_gb_function - 600 + 1
        IF(ider2>3) ider2 = ider1-2
        la_max_local=la_max+2
        la_min_local=MAX(la_min-2,0)
        lb_max_local=lb_max+2
        lb_min_local=MAX(lb_min-2,0)
        ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
        ! is equivalent to mapping pab with
        !   d_{ider1} pgf_a d_{ider1} pgf_b
        ALLOCATE(pab_local(ncoset(la_max_local),ncoset(lb_max_local)))
        pab_local = 0.0_dp
        DO lxa=0,la_max
        DO lxb=0,lb_max
           DO lya=0,la_max-lxa
           DO lyb=0,lb_max-lxb
              DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
              DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
                 ! this element of pab results in 16 elements of pab_local
                  CALL prepare_dijadijb (pab_local,pab,ider1,ider2,&
                      lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)
              END DO
              END DO
           END DO
           END DO
        END DO
        END DO
        o1_local=0
        o2_local=0
    CASE(FUNC_DXDX,FUNC_DYDY,FUNC_DZDZ)
        ider1 = ga_gb_function - 603
        la_max_local=la_max+2
        la_min_local=MAX(la_min-2,0)
        lb_max_local=lb_max+2
        lb_min_local=MAX(lb_min-2,0)
        ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
        ! is equivalent to mapping pab with
        !   dd_{ider1} pgf_a dd_{ider1} pgf_b

        ALLOCATE(pab_local(ncoset(la_max_local),ncoset(lb_max_local)))
        pab_local = 0.0_dp
        DO lxa=0,la_max
        DO lxb=0,lb_max
           DO lya=0,la_max-lxa
           DO lyb=0,lb_max-lxb
              DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
              DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
                 ! this element of pab results in 9 elements of pab_local
                  CALL prepare_diiadiib (pab_local,pab,ider1,&
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
    CASE(FUNC_ARB)
        CPPrecondition(PRESENT(ir),cp_failure_level,routineP,error,failure)
        la_max_local=la_max
        la_min_local=la_min
        lb_max_local=lb_max+1
        lb_min_local=lb_min
        ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
        ! is equivalent to mapping pab with
        ! pgf_a (r-Rb)_{ir} pgf_b = pgf_a * pgf_{b+1ir}
        ALLOCATE(pab_local(ncoset(la_max_local),ncoset(lb_max_local)))
        pab_local = 0.0_dp
        DO lxa=0,la_max
        DO lxb=0,lb_max
           DO lya=0,la_max-lxa
           DO lyb=0,lb_max-lxb
              DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
              DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
                 ! this element of pab results in 4 elements of pab_local
                 CALL prepare_arb(pab_local,pab,ir,lxa,lya,lza,lxb,lyb,lzb,o1,o2)
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

    ng(:) = rsgrid%desc%npts(:)
    grid => rsgrid%r(:,:,:)
    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)


    IF(PRESENT(ir).AND.PRESENT(rsgauge).AND.PRESENT(rsbuf)) THEN
       grid => rsbuf%r(:,:,:)
       grid_tmp => rsgrid%r(:,:,:)
       gauge => rsgauge%r(:,:,:)
    ENDIF

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

    IF (subpatch_collocate) THEN
        CALL collocate_general_subpatch()
    ELSE
       IF (rsgrid%desc%orthorhombic ) THEN
          CALL collocate_ortho()
          ! CALL collocate_general()
       ELSE
          CALL collocate_general_wings()
          !CALL collocate_general_opt()
       END IF
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
! *****************************************************************************
    SUBROUTINE collocate_ortho()


!   *** properties of the grid ***

    ! notice we're in the ortho case
    dr(1) = rsgrid%desc%dh(1,1)
    dr(2) = rsgrid%desc%dh(2,2)
    dr(3) = rsgrid%desc%dh(3,3)

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
    CALL compute_cube_center(cubecenter,rsgrid%desc,zeta,zetb,ra,rab)
    roffset(:)    = rp(:) - REAL(cubecenter(:),dp)*dr(:)
!   *** a mapping so that the ig corresponds to the right grid point
    DO i=1,3
      IF ( rsgrid % desc % perd ( i ) == 1 ) THEN
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
        offset=MODULO(cubecenter(i)+lb_cube(i)+rsgrid%desc%lb(i)-rsgrid%lb_local(i),ng(i))+1-lb_cube(i)
        ! check for out of bounds
        IF (ub_cube(i)+offset>UBOUND(grid,i).OR.lb_cube(i)+offset<LBOUND(grid,i)) THEN
           CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        ENDIF
        DO ig=lb_cube(i),ub_cube(i)
           map(ig,i) = ig+offset
        END DO
      END IF
    ENDDO
    ALLOCATE(pol_z(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_y(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_x(0:lp,-cmax:cmax))

    IF(PRESENT(ir).AND.PRESENT(rsgauge)) CALL collocate_ortho_set_to_0()

#include "prep.f90"

    IF ( PRESENT ( lgrid ) ) THEN
       ig = lgrid%ldim * ithread_l + 1
#include "call_collocate_omp.f90"
    ELSE
#include "call_collocate.f90"
    END IF

    IF(PRESENT(ir).AND.PRESENT(rsgauge)) CALL collocate_gauge_ortho()

    ! deallocation needed to pass around a pgi bug..
    DEALLOCATE(pol_z)
    DEALLOCATE(pol_y)
    DEALLOCATE(pol_x)
    DEALLOCATE(map)

    END SUBROUTINE collocate_ortho

    SUBROUTINE collocate_gauge_ortho()
    INTEGER                                  :: i, igmax, igmin, j, j2, jg, &
                                                jg2, jgmin, k, k2, kg, kg2, &
                                                kgmin, sci
    REAL(KIND=dp)                            :: point(3,4), res(4), x, y, y2, &
                                                z, z2

! notice we're in the ortho case

      dr(1) = rsgrid%desc%dh(1,1)
      dr(2) = rsgrid%desc%dh(2,2)
      dr(3) = rsgrid%desc%dh(3,3)
      !
      sci=1
      kgmin=sphere_bounds(sci)
      sci=sci+1
      DO kg=kgmin,0
         kg2=1-kg
         k=map(kg,3)
         k2=map(kg2,3)
         jgmin=sphere_bounds(sci)
         sci=sci+1
         z  = (REAL( kg,dp)+REAL(cubecenter(3),dp))*dr(3)
         z2 = (REAL(kg2,dp)+REAL(cubecenter(3),dp))*dr(3)
         DO jg=jgmin,0
            jg2=1-jg
            j=map(jg,2)
            j2=map(jg2,2)
            igmin=sphere_bounds(sci)
            sci=sci+1
            igmax=1-igmin
            y  = (REAL( jg,dp)+REAL(cubecenter(2),dp))*dr(2)
            y2 = (REAL(jg2,dp)+REAL(cubecenter(2),dp))*dr(2)
            DO ig=igmin,igmax
               i=map(ig,1)
               x = (REAL(ig,dp)+REAL(cubecenter(1),dp))*dr(1)
               point(1,1) = x;point(2,1) = y ;point(3,1) = z
               point(1,2) = x;point(2,2) = y2;point(3,2) = z
               point(1,3) = x;point(2,3) = y ;point(3,3) = z2
               point(1,4) = x;point(2,4) = y2;point(3,4) = z2
               !
               res(1) = ( point(ir,1) - rb(ir) ) - gauge(i, j, k)
               res(2) = ( point(ir,2) - rb(ir) ) - gauge(i,j2, k)
               res(3) = ( point(ir,3) - rb(ir) ) - gauge(i, j,k2)
               res(4) = ( point(ir,4) - rb(ir) ) - gauge(i,j2,k2)
               !
               grid_tmp(i, j, k) = grid_tmp(i, j, k) + grid(i, j, k) * res(1)
               grid_tmp(i,j2, k) = grid_tmp(i,j2, k) + grid(i,j2, k) * res(2)
               grid_tmp(i, j,k2) = grid_tmp(i, j,k2) + grid(i, j,k2) * res(3)
               grid_tmp(i,j2,k2) = grid_tmp(i,j2,k2) + grid(i,j2,k2) * res(4)
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE collocate_gauge_ortho

    SUBROUTINE collocate_ortho_set_to_0()
    INTEGER                                  :: i, igmax, igmin, j, j2, jg, &
                                                jg2, jgmin, k, k2, kg, kg2, &
                                                kgmin, sci

!

      dr(1) = rsgrid%desc%dh(1,1)
      dr(2) = rsgrid%desc%dh(2,2)
      dr(3) = rsgrid%desc%dh(3,3)
      !
      sci=1
      kgmin=sphere_bounds(sci)
      sci=sci+1
      DO kg=kgmin,0
         kg2=1-kg
         k=map(kg,3)
         k2=map(kg2,3)
         jgmin=sphere_bounds(sci)
         sci=sci+1
         DO jg=jgmin,0
            jg2=1-jg
            j=map(jg,2)
            j2=map(jg2,2)
            igmin=sphere_bounds(sci)
            sci=sci+1
            igmax=1-igmin
            DO ig=igmin,igmax
               i=map(ig,1)
               grid(i, j, k) = 0.0_dp
               grid(i,j2, k) = 0.0_dp
               grid(i, j,k2) = 0.0_dp
               grid(i,j2,k2) = 0.0_dp
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE collocate_ortho_set_to_0

!
!   this is a general 'optimized' routine to do the collocation
!
! *****************************************************************************
    SUBROUTINE collocate_general_opt()

    INTEGER :: i, i_index, il, ilx, ily, ilz, index_max(3), index_min(3), &
      ismax, ismin, j, j_index, jl, jlx, jly, jlz, k, k_index, kl, klx, kly, &
      klz, lpx, lpy, lpz, lx, ly, lz, offset(3)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: grid_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: coef_map
    REAL(KIND=dp) :: a, b, c, d, di, dip, dj, djp, dk, dkp, exp0i, exp1i, &
      exp2i, gp(3), hmatgrid(3,3), pointj(3), pointk(3), res, v(3)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: coef_ijk
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: hmatgridp

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
      hmatgrid=rsgrid%desc%dh

      ! center in grid coords
      gp=MATMUL(rsgrid%desc%dh_inv,rp)
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

      offset(:)=MODULO(index_min(:)+rsgrid%desc%lb(:)-rsgrid%lb_local(:),ng(:))+1

      ALLOCATE(grid_map(index_min(1):index_max(1)))
      DO i=index_min(1),index_max(1)
        grid_map(i)=MODULO(i,ng(1))+1
        IF (rsgrid % desc % perd ( 1 )==1) THEN
           grid_map(i)=MODULO(i,ng(1))+1
        ELSE
           grid_map(i)=i-index_min(1)+offset(1)
        ENDIF
      ENDDO

      ! go over the grid, but cycle if the point is not within the radius
      DO k=index_min(3),index_max(3)
        dk=k-gp(3)
        pointk=hmatgrid(:,3)*dk

        IF (rsgrid % desc % perd ( 3 )==1) THEN
           k_index=MODULO(k,ng(3))+1
        ELSE
           k_index=k-index_min(3)+offset(3)
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
          IF (rsgrid % desc % perd ( 2 )==1) THEN
             j_index=MODULO(j,ng(2))+1
          ELSE
             j_index=j-index_min(2)+offset(2)
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
              d=SQRT(d)
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
          exp2i=EXP((a*i+b)*i+c)
          exp1i=EXP(2*a*i+a+b)
          exp0i=EXP(2*a)

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

    SUBROUTINE collocate_general_subpatch
    INTEGER, DIMENSION(2, 3)                 :: local_b
    INTEGER, DIMENSION(3)                    :: local_s, periodic
    REAL(dp), &
      DIMENSION((lp+1)*(lp+2)*(lp+3)/6)      :: poly_d3

        periodic=1 ! cell%perd
        CALL poly_cp2k2d3(coef_xyz,lp,poly_d3,error)
        local_b(1,:)=rsgrid%lb_real-rsgrid%desc%lb
        local_b(2,:)=rsgrid%ub_real-rsgrid%desc%lb
        local_s=rsgrid%lb_real-rsgrid%lb_local
        IF (BTEST(subpatch_pattern,0)) local_b(1,1)=local_b(1,1)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,1)) local_b(2,1)=local_b(2,1)+rsgrid%desc%border
        IF (BTEST(subpatch_pattern,2)) local_b(1,2)=local_b(1,2)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,3)) local_b(2,2)=local_b(2,2)+rsgrid%desc%border
        IF (BTEST(subpatch_pattern,4)) local_b(1,3)=local_b(1,3)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,5)) local_b(2,3)=local_b(2,3)+rsgrid%desc%border
        IF (BTEST(subpatch_pattern,0)) local_s(1)=local_s(1)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,2)) local_s(2)=local_s(2)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,4)) local_s(3)=local_s(3)-rsgrid%desc%border
        IF ( PRESENT ( lgrid ) ) THEN
          CALL collocGauss(h=cell%hmat,h_inv=cell%h_inv,&
            grid=grid,poly=poly_d3,alphai=zetp,posi=rp,max_r2=radius*radius,&
            periodic=periodic,gdim=ng,local_bounds=local_b,local_shift=local_s,&
            lgrid=lgrid,error=error)
        ELSE
          CALL collocGauss(h=cell%hmat,h_inv=cell%h_inv,&
            grid=grid,poly=poly_d3,alphai=zetp,posi=rp,max_r2=radius*radius,&
            periodic=periodic,gdim=ng,local_bounds=local_b,local_shift=local_s,&
            error=error)
        END IF
        ! defaults: local_shift=(/0,0,0/),poly_shift=(/0.0_dp,0.0_dp,0.0_dp/),scale=1.0_dp,

    END SUBROUTINE

    SUBROUTINE collocate_general_wings
    INTEGER, DIMENSION(2, 3)                 :: local_b
    INTEGER, DIMENSION(3)                    :: periodic
    REAL(dp), &
      DIMENSION((lp+1)*(lp+2)*(lp+3)/6)      :: poly_d3
    REAL(dp), DIMENSION(3)                   :: local_shift, rShifted

        periodic=1 ! cell%perd
        CALL poly_cp2k2d3(coef_xyz,lp,poly_d3,error)
        local_b(1,:)=0
        local_b(2,:)=MIN(rsgrid%desc%npts-1,rsgrid%ub_local-rsgrid%lb_local)
        local_shift=REAL(rsgrid%desc%lb-rsgrid%lb_local,dp)/REAL(rsgrid%desc%npts,dp)
        rShifted(1)=rp(1)+cell%hmat(1,1)*local_shift(1)&
             +cell%hmat(1,2)*local_shift(2)&
             +cell%hmat(1,3)*local_shift(3)
        rShifted(2)=rp(2)+cell%hmat(2,1)*local_shift(1)&
             +cell%hmat(2,2)*local_shift(2)&
             +cell%hmat(2,3)*local_shift(3)
        rShifted(3)=rp(3)+cell%hmat(3,1)*local_shift(1)&
             +cell%hmat(3,2)*local_shift(2)&
             +cell%hmat(3,3)*local_shift(3)
        IF ( PRESENT ( lgrid ) ) THEN
          CALL collocGauss(h=cell%hmat,h_inv=cell%h_inv,&
            grid=grid,poly=poly_d3,alphai=zetp,posi=rShifted,max_r2=radius*radius,&
            periodic=periodic,gdim=ng,local_bounds=local_b,&
            lgrid=lgrid,error=error)
        ELSE
          CALL collocGauss(h=cell%hmat,h_inv=cell%h_inv,&
            grid=grid,poly=poly_d3,alphai=zetp,posi=rShifted,max_r2=radius*radius,&
            periodic=periodic,gdim=ng,local_bounds=local_b,&
            error=error)
        END IF
        ! defaults: local_shift=(/0,0,0/),poly_shift=(/0.0_dp,0.0_dp,0.0_dp/),scale=1.0_dp,

    END SUBROUTINE

!
!   this is a general 'reference' routine to do the collocation
!
! *****************************************************************************
    SUBROUTINE collocate_general()

    INTEGER                                  :: i, index_max(3), &
                                                index_min(3), ipoint(3), j, k
    REAL(KIND=dp)                            :: point(3)

! still hard-wired (see MODULO)

      CPPostcondition(ALL(rsgrid % desc % perd==1),cp_failure_level,routineP,error,failure)

      CALL return_cube_nonortho(cube_info,radius,index_min,index_max,rp)

      ! go over the grid, but cycle if the point is not within the radius
      DO k=index_min(3),index_max(3)
      DO j=index_min(2),index_max(2)
      DO i=index_min(1),index_max(1)
         ! point in real space
         point=MATMUL(cell%hmat,REAL((/i,j,k/),KIND=dp)/ng)
         ! skip if outside of the sphere
         IF (SUM((point-rp)**2)>radius**2) CYCLE
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

! *****************************************************************************
    FUNCTION primitive_value(point) RESULT(res)
    REAL(KIND=dp)                            :: point(3), res

    REAL(KIND=dp)                            :: dra(3), drap(3), drb(3), &
                                                drbp(3), myexp

        res=0.0_dp
        myexp=EXP(-zetp*SUM((point-rp)**2))*prefactor
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

! *****************************************************************************
!> \brief low level collcation of primitive gaussian functions in g-space
! *****************************************************************************
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
    INTEGER :: ax, ay, az, bx, by, bz, handle, i, ico, ig, ig2, jco, jg, kg, &
      la, lb, lb_cube_min, lb_grid, stat, ub_cube_max, ub_grid
    INTEGER, DIMENSION(3)                    :: lb_cube, ub_cube
    REAL(dp)                                 :: f, fa, fb, pij, prefactor, &
                                                rzetp, twozetp, zetp
    REAL(dp), DIMENSION(3)                   :: dg, expfactor, fap, fbp, rap, &
                                                rbp, rp
    REAL(dp), DIMENSION(:), POINTER          :: g

    CALL timeset(routineN,handle)

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

    DEALLOCATE (cubeaxis,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cubeaxis")
    DEALLOCATE (g,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"g")
    DEALLOCATE (rag,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rag")
    DEALLOCATE (rbg,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rbg")

    CALL timestop(handle)

  END SUBROUTINE collocate_pgf_product_gspace

END MODULE qs_collocate_density
