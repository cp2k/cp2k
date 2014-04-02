!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Build up the plane wave density by collocating the primitive Gaussian
!>      functions (pgf).
!> \par History
!>      Joost VandeVondele (02.2002)
!>            1) rewrote collocate_pgf for increased accuracy and speed
!>            2) collocate_core hack for PGI compiler
!>            3) added multiple grid feature
!>            4) new way to go over the grid
!>      Joost VandeVondele (05.2002)
!>            1) prelim. introduction of the real space grid type
!>      JGH [30.08.02] multigrid arrays independent from potential
!>      JGH [17.07.03] distributed real space code
!>      JGH [23.11.03] refactoring and new loop ordering
!>      JGH [04.12.03] OpneMP parallelization of main loops
!>      Joost VandeVondele (12.2003)
!>           1) modified to compute tau
!>      Joost removed incremental build feature
!>      Joost introduced map consistent
!>      Rewrote grid integration/collocation routines, [Joost VandeVondele,03.2007]
!> \author Matthias Krack (03.04.2001)
! *****************************************************************************
MODULE qs_integrate_potential
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE atprop_types,                    ONLY: atprop_array_init
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_col_block_sizes, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_distribution, cp_dbcsr_finalize, cp_dbcsr_get_block_p, &
       cp_dbcsr_get_data_size, cp_dbcsr_get_matrix_type, &
       cp_dbcsr_get_num_blocks, cp_dbcsr_init, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_work_create
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_add_block_node,&
                                             cp_dbcsr_deallocate_matrix
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_spline_utils,                 ONLY: pw_restrict_s3
  USE cube_utils,                      ONLY: cube_info_type
  USE dbcsr_dist_operations
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_has_threads
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE input_constants,                 ONLY: pw_interp,&
                                             spline3_pbc_interp,&
                                             use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE lri_environment_types,           ONLY: lri_kind_type
  USE mathconstants,                   ONLY: dfac,&
                                             pi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: coset,&
                                             nco,&
                                             ncoset,&
                                             nso,&
                                             nsoset
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_copy,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                             pw_pools_create_pws,&
                                             pw_pools_give_back_pws
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_integrate_potential_low,      ONLY: integrate_pgf_product_rspace
  USE realspace_grid_types,            ONLY: &
       pw2rs, realspace_grid_desc_p_type, realspace_grid_desc_type, &
       realspace_grid_p_type, realspace_grid_type, rs_grid_release, &
       rs_grid_retain, rs_grid_zero, rs_pw_transfer
  USE scptb_types,                     ONLY: get_scptb_parameter,&
                                             scp_vector_type,&
                                             scptb_parameter_type
  USE task_list_methods,               ONLY: int2pair,&
                                             rs_distribute_matrix
  USE task_list_types,                 ONLY: task_list_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  INTEGER :: debug_count=0

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_integrate_potential'

! *** Public subroutines ***

  PUBLIC :: integrate_v_rspace,&
            integrate_v_rspace_one_center,&
            integrate_v_core_rspace,&
            integrate_ppl_rspace,&
            integrate_scp_rspace,&
            integrate_pgf_product_rspace,&
            potential_pw2rs, integrate_rho_nlcc

CONTAINS

! *****************************************************************************
!> \brief integrates the SCP charge functions on the grid potential
! *****************************************************************************
  SUBROUTINE integrate_scp_rspace(scp_pot,qs_env,scpv,calculate_forces,error)

    TYPE(pw_p_type), INTENT(INOUT)           :: scp_pot
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(scp_vector_type), POINTER           :: scpv
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_scp_rspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, i, iatom, &
                                                ierr, ii, ikind, j, jj, l, &
                                                lmaxscp, natom_of_kind, ni, &
                                                nj, npme, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores
    LOGICAL                                  :: defined, failure, use_virial
    REAL(KIND=dp)                            :: alpha, dvol, eps_rho_rspace, &
                                                norm, pp
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra
    REAL(KIND=dp), DIMENSION(3, 3)           :: my_virial_a, my_virial_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab, pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_v
    TYPE(scptb_parameter_type), POINTER      :: scptb_kind
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(pw_env,auxbas_rs_desc,cores)

    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env=pw_env,auxbas_rs_desc=auxbas_rs_desc, &
                    auxbas_rs_grid=rs_v,error=error)
    CALL rs_grid_retain(rs_v,error=error)

    CALL rs_pw_transfer(rs_v,scp_pot%pw,pw2rs,error=error)

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         particle_set=particle_set,&
         para_env=para_env,pw_env=pw_env,&
         force=force,virial=virial,error=error)

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer).AND.calculate_forces

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    dvol = scp_pot%pw%pw_grid%dvol

    DO ikind=1,SIZE(atomic_kind_set)

       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            natom=natom_of_kind,atom_list=atom_list,&
                            scptb_parameter=scptb_kind)
       CALL get_scptb_parameter(scptb_kind,defined=defined,lmaxscp=lmaxscp,ag=alpha)
       IF (.NOT.defined) CYCLE

       ni = ncoset(lmaxscp)
       ALLOCATE(hab(ni,1),pab(ni,1),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       pab = 0._dp
       hab = 0._dp

       ALLOCATE(cores(natom_of_kind),STAT=ierr)
       CPPrecondition(ierr==0,cp_failure_level,routineP,error,failure)
       cores = 0
       npme = 0

       DO iatom = 1, natom_of_kind
          atom_a = atom_list(iatom)
          ra(:) = pbc(particle_set(atom_a)%r,cell)
          IF(rs_v%desc%parallel .AND. .NOT. rs_v%desc%distributed) THEN
              ! replicated realspace grid, split the atoms up between procs
              IF (MODULO(iatom,rs_v%desc%group_size) == rs_v % desc % my_pos ) THEN
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
         hab(:,1) = 0.0_dp
         IF (calculate_forces) THEN
            force_a(:) = 0.0_dp
            force_b(:) = 0.0_dp
            IF (use_virial) THEN
               my_virial_a = 0.0_dp
               my_virial_b = 0.0_dp
            END IF
            pab = 0._dp
            DO l=0,lmaxscp
               pp = (2._dp*l+3._dp)/2._dp
               norm = 2._dp**(l+2)/SQRT(pi)/dfac(2*l+1)
               norm = SQRT(0.25_dp*dfac(2*l+1)/pi) * norm*alpha**pp
               DO jj=1,nco(l)
                  j = ncoset(l-1) + jj
                  DO ii=1,nso(l)
                     i = nsoset(l-1) + ii
                     pab(j,1) = pab(j,1) + orbtramat(l)%s2c(ii,jj)*norm
                  END DO
               END DO
            END DO
         END IF

         CALL integrate_pgf_product_rspace(lmaxscp,alpha,0,&
              0,0.0_dp,0,ra,(/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,&
              rs_v,cell,pw_env%cube_info(1),hab,pab=pab,o1=0,o2=0,&
              eps_gvg_rspace=eps_rho_rspace,&
              calculate_forces=calculate_forces,force_a=force_a,&
              force_b=force_b,use_virial=use_virial,my_virial_a=my_virial_a,&
              my_virial_b=my_virial_b,use_subpatch=.TRUE.,subpatch_pattern=0_int_8,error=error)

         DO l=0,lmaxscp
            pp = (2._dp*l+3._dp)/2._dp
            norm = 2._dp**(l+2)/SQRT(pi)/dfac(2*l+1) * dvol
            norm = SQRT(0.25_dp*dfac(2*l+1)/pi) * norm*alpha**pp
            DO jj=1,nco(l)
               j = ncoset(l-1) + jj
               DO ii=1,nso(l)
                  i = nsoset(l-1) + ii
                  scpv%vector(ikind)%vmat(i,iatom) = scpv%vector(ikind)%vmat(i,iatom) + &
                      hab(j,1)*orbtramat(l)%c2s(ii,jj)*norm
               END DO
            END DO
         END DO

         IF (calculate_forces) THEN
            force(ikind)%rho_elec(:,iatom) = force(ikind)%rho_elec(:,iatom) + force_a(:)*dvol
            IF (use_virial) THEN
              virial%pv_virial = virial%pv_virial + my_virial_a*dvol
            END IF
         END IF
       END DO

       DEALLOCATE(cores,hab,pab,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    END DO

    CALL rs_grid_release(rs_v, error=error)

    CALL timestop(handle)

  END SUBROUTINE integrate_scp_rspace
! *****************************************************************************
!> \brief computes the forces/virial due to the local pseudopotential
! *****************************************************************************
  SUBROUTINE integrate_ppl_rspace(rho_rspace,qs_env,error)
    TYPE(pw_p_type), INTENT(INOUT)           :: rho_rspace
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_ppl_rspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, ikind, &
                                                j, lppl, n, natom_of_kind, &
                                                ni, npme, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores
    LOGICAL                                  :: failure, use_virial
    REAL(KIND=dp)                            :: alpha, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra
    REAL(KIND=dp), DIMENSION(3, 3)           :: my_virial_a, my_virial_b
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cexp_ppl
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab, pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_v
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(pw_env,auxbas_rs_desc,cores)

    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env=pw_env,auxbas_rs_desc=auxbas_rs_desc, &
                    auxbas_rs_grid=rs_v,error=error)
    CALL rs_grid_retain(rs_v,error=error)

    CALL rs_pw_transfer(rs_v,rho_rspace%pw,pw2rs,error=error)

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         particle_set=particle_set,&
         para_env=para_env,pw_env=pw_env,&
         force=force,virial=virial,error=error)

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

    DO ikind=1,SIZE(atomic_kind_set)

       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            natom=natom_of_kind,&
                            atom_list=atom_list,&
                            gth_potential=gth_potential)

       IF (.NOT.ASSOCIATED(gth_potential)) CYCLE
       CALL get_potential(potential=gth_potential,alpha_ppl=alpha,nexp_ppl=lppl,cexp_ppl=cexp_ppl)

       IF ( lppl <= 0 ) CYCLE

       ni = ncoset(2*lppl-2)
       ALLOCATE(hab(ni,1),pab(ni,1),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       pab = 0._dp

       CALL reallocate ( cores, 1, natom_of_kind )
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

       DO iatom = 1, natom_of_kind
          atom_a = atom_list(iatom)
          ra(:) = pbc(particle_set(atom_a)%r,cell)
          IF(rs_v%desc%parallel .AND. .NOT. rs_v%desc%distributed) THEN
              ! replicated realspace grid, split the atoms up between procs
              IF (MODULO(iatom,rs_v%desc%group_size) == rs_v % desc % my_pos ) THEN
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
         hab(:,1) = 0.0_dp
         force_a(:) = 0.0_dp
         force_b(:) = 0.0_dp
         IF (use_virial) THEN
            my_virial_a = 0.0_dp
            my_virial_b = 0.0_dp
         END IF
         ni = 2*lppl-2

         CALL integrate_pgf_product_rspace(ni,alpha,0,&
              0,0.0_dp,0,ra,(/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,&
              rs_v,cell,pw_env%cube_info(1),hab,pab=pab,o1=0,o2=0,&
              eps_gvg_rspace=eps_rho_rspace,&
              calculate_forces=.TRUE.,force_a=force_a,&
              force_b=force_b,use_virial=use_virial,my_virial_a=my_virial_a,&
              my_virial_b=my_virial_b,use_subpatch=.TRUE.,subpatch_pattern=0_int_8,error=error)

         force(ikind)%gth_ppl(:,iatom) =&
           force(ikind)%gth_ppl(:,iatom) + force_a(:)*rho_rspace%pw%pw_grid%dvol

         IF (use_virial) THEN
           virial%pv_virial = virial%pv_virial + my_virial_a*rho_rspace%pw%pw_grid%dvol
           CALL cp_unimplemented_error(fromWhere=routineP, &
                message="Virial not debuged for CORE_PPL", &
                error=error, error_level=cp_failure_level)
         END IF
       END DO

       DEALLOCATE(hab,pab,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    END DO

    CALL rs_grid_release(rs_v, error=error)

    DEALLOCATE(cores,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE integrate_ppl_rspace

! *****************************************************************************
!> \brief computes the forces/virial due to the nlcc pseudopotential
! *****************************************************************************
  SUBROUTINE integrate_rho_nlcc(rho_rspace,qs_env,error)
    TYPE(pw_p_type), INTENT(INOUT)           :: rho_rspace
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_rho_nlcc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, &
                                                iexp_nlcc, ikind, ithread, j, &
                                                n, natom, nc, nexp_nlcc, ni, &
                                                npme, nthread, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores, nct_nlcc
    LOGICAL                                  :: failure, nlcc, use_virial
    REAL(KIND=dp)                            :: alpha, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra
    REAL(KIND=dp), DIMENSION(3, 3)           :: my_virial_a, my_virial_b
    REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha_nlcc
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cval_nlcc, hab, pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_v
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(pw_env,auxbas_rs_desc,cores)

    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env=pw_env,auxbas_rs_desc=auxbas_rs_desc, &
                    auxbas_rs_grid=rs_v,error=error)
    CALL rs_grid_retain(rs_v,error=error)

    CALL rs_pw_transfer(rs_v,rho_rspace%pw,pw2rs,error=error)

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         particle_set=particle_set,&
         para_env=para_env,pw_env=pw_env,&
         force=force,virial=virial,error=error)

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

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

         nthread = 1
         ithread=0

         ALLOCATE(hab(ni,1),pab(ni,1),STAT=stat)
         CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
         pab = 0._dp

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
         IF(dft_control%nspins==2)pab=pab*0.5_dp

         DO iatom = 1, natom
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r,cell)
            IF(rs_v%desc%parallel .AND. .NOT. rs_v%desc%distributed) THEN
                ! replicated realspace grid, split the atoms up between procs
                IF (MODULO(iatom,rs_v%desc%group_size) == rs_v % desc % my_pos ) THEN
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
           hab(:,1) = 0.0_dp
           force_a(:) = 0.0_dp
           force_b(:) = 0.0_dp
           IF (use_virial) THEN
              my_virial_a = 0.0_dp
              my_virial_b = 0.0_dp
           END IF
           ni = 2*nc-2

           CALL integrate_pgf_product_rspace(ni,1/(2*alpha**2),0,&
                0,0.0_dp,0,ra,(/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,&
                rs_v,cell,pw_env%cube_info(1),hab,pab=pab,o1=0,o2=0,&
                eps_gvg_rspace=eps_rho_rspace,&
                calculate_forces=.TRUE.,force_a=force_a,&
                force_b=force_b,use_virial=use_virial,my_virial_a=my_virial_a,&
                my_virial_b=my_virial_b,use_subpatch=.TRUE.,subpatch_pattern=0_int_8,error=error)

           force(ikind)%gth_nlcc(:,iatom) =&
             force(ikind)%gth_nlcc(:,iatom) + force_a(:)*rho_rspace%pw%pw_grid%dvol

           IF (use_virial) THEN
             virial%pv_virial = virial%pv_virial + my_virial_a*rho_rspace%pw%pw_grid%dvol
           END IF
         END DO

         DEALLOCATE(hab,pab,STAT=stat)
         CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

      END DO

    END DO

    CALL rs_grid_release(rs_v, error=error)

    DEALLOCATE(cores,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE integrate_rho_nlcc

! *****************************************************************************
!> \brief computes the forces/virial due to the ionic cores with a potential on
!>      grid
! *****************************************************************************
  SUBROUTINE integrate_v_core_rspace(v_rspace,qs_env,error)
    TYPE(pw_p_type), INTENT(INOUT)           :: v_rspace
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_v_core_rspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, ikind, &
                                                j, natom, natom_of_kind, &
                                                npme, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores
    LOGICAL                                  :: atenergy, failure, paw_atom, &
                                                skip_fcore, use_virial
    REAL(KIND=dp)                            :: alpha_core_charge, &
                                                ccore_charge, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra
    REAL(KIND=dp), DIMENSION(3, 3)           :: my_virial_a, my_virial_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab, pab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc
    TYPE(realspace_grid_type), POINTER       :: rs_v
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    !If gapw, check for gpw kinds
    skip_fcore = .FALSE.
    IF(qs_env%dft_control%qs_control%gapw) THEN
      IF(.NOT. qs_env%dft_control%qs_control%gapw_control%nopaw_as_gpw) skip_fcore = .TRUE.
    END IF
    ! atomic energy contributions
    atenergy = .FALSE.
    IF(ASSOCIATED(qs_env%atprop)) THEN
      atenergy=qs_env%atprop%energy
    END IF

    failure = .FALSE.
    IF(.NOT. skip_fcore) THEN
        NULLIFY(pw_env,auxbas_rs_desc,hab,pab,cores)
        ALLOCATE (hab(1,1),STAT=stat)
        CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE (pab(1,1),STAT=stat)
        CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

        CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
        CALL pw_env_get(pw_env=pw_env,auxbas_rs_desc=auxbas_rs_desc, &
                        auxbas_rs_grid=rs_v,error=error)
        CALL rs_grid_retain(rs_v,error=error)

        CALL rs_pw_transfer(rs_v,v_rspace%pw,pw2rs,error=error)

        CALL get_qs_env(qs_env=qs_env,&
             atomic_kind_set=atomic_kind_set,&
             cell=cell,&
             dft_control=dft_control,&
             particle_set=particle_set,&
             para_env=para_env,pw_env=pw_env,&
             force=force,virial=virial,error=error)

        IF(atenergy) THEN
           natom = SIZE(particle_set)
           CALL atprop_array_init(qs_env%atprop%ateb,natom,error)
        END IF

        use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

        eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

        DO ikind=1,SIZE(atomic_kind_set)

           atomic_kind => atomic_kind_set(ikind)

           CALL get_atomic_kind(atomic_kind=atomic_kind,&
                natom=natom_of_kind,&
                paw_atom=paw_atom, &
                atom_list=atom_list,&
                alpha_core_charge=alpha_core_charge,&
                ccore_charge=ccore_charge)

           IF(paw_atom) THEN
                force(ikind)%rho_core(:,:) =  0.0_dp
                CYCLE
           END IF
           pab(1,1) = -ccore_charge

           IF (ASSOCIATED(force)) THEN
              force(ikind)%rho_core=0.0_dp
           ENDIF

           IF (alpha_core_charge == 0.0_dp .OR. pab(1,1)== 0.0_dp) CYCLE

           CALL reallocate ( cores, 1, natom_of_kind )
           npme = 0
           cores = 0

           DO iatom = 1, natom_of_kind
              atom_a = atom_list(iatom)
              ra(:) = pbc(particle_set(atom_a)%r,cell)
              IF(rs_v%desc%parallel .AND. .NOT. rs_v%desc%distributed) THEN
                  ! replicated realspace grid, split the atoms up between procs
                  IF (MODULO(iatom,rs_v%desc%group_size) == rs_v % desc % my_pos ) THEN
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
            hab(1,1) = 0.0_dp
            force_a(:) = 0.0_dp
            force_b(:) = 0.0_dp
            IF (use_virial) THEN
              my_virial_a = 0.0_dp
              my_virial_b = 0.0_dp
            END IF

            CALL integrate_pgf_product_rspace(0,alpha_core_charge,0,&
                 0,0.0_dp,0,ra,(/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,&
                 rs_v,cell,pw_env%cube_info(1),hab,pab=pab,o1=0,o2=0,&
                 eps_gvg_rspace=eps_rho_rspace,&
                 calculate_forces=.TRUE.,force_a=force_a,&
                 force_b=force_b,use_virial=use_virial,my_virial_a=my_virial_a,&
                 my_virial_b=my_virial_b,use_subpatch=.TRUE.,subpatch_pattern=0_int_8,error=error)

            IF (ASSOCIATED(force)) THEN
               force(ikind)%rho_core(:,iatom) =&
                 force(ikind)%rho_core(:,iatom) + force_a(:)
            ENDIF

            IF (use_virial) THEN
              virial%pv_virial = virial%pv_virial + my_virial_a
            END IF
            IF (atenergy) THEN
               qs_env%atprop%ateb(atom_a) = qs_env%atprop%ateb(atom_a) + 0.5_dp*hab(1,1)*pab(1,1)
            END IF

         END DO

        END DO

        CALL rs_grid_release(rs_v, error=error)

        DEALLOCATE (hab,pab,cores,STAT=stat)
        CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    END IF

    CALL timestop(handle)

  END SUBROUTINE integrate_v_core_rspace

! *****************************************************************************
!> \brief computes integrals of product of v_rspace times a one-center function  
!>        required for LRIGPW
!> \author Dorothea Golze
! *****************************************************************************
  SUBROUTINE integrate_v_rspace_one_center(v_rspace,qs_env,int_res,error)
    TYPE(pw_p_type), INTENT(INOUT)           :: v_rspace
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(lri_kind_type), DIMENSION(:), &
      POINTER                                :: int_res
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'integrate_v_rspace_one_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, handle, i, iatom, igrid_level, ikind, ipgf, iset, j, &
      maxco, maxsgf_set, na1, natom, natom_of_kind, ncoa, nkind, npme, nseta, &
      nsgfa, sgfa, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cores, la_max, &
                                                la_min, npgfa, nsgf_seta
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: distributed_grids, failure, &
                                                map_consistent, parallel, &
                                                use_virial
    REAL(KIND=dp)                            :: eps_gvg_rspace, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab, rpgfa, sphi_a, work, zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: lri_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_v
    TYPE(section_vals_type), POINTER         :: input, interp_section
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    NULLIFY(cell, cores,hab,input, interp_section,rs_descs,pw_env,para_env,&
            dft_control,rs_v,virial,work)

    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,&
                    input=input,error=error)

    CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_v, error=error)
    DO i=1,SIZE(rs_v)
      CALL rs_grid_retain(rs_v(i)%rs_grid,error=error)
      CALL rs_grid_zero(rs_v(i)%rs_grid)
    END DO

    gridlevel_info=>pw_env%gridlevel_info
    interp_section => section_vals_get_subs_vals(input,"DFT%MGRID%INTERPOLATOR",&
         error=error)
    
    CALL potential_pw2rs(rs_v,v_rspace,pw_env,interp_section,error)

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         nkind=nkind,&
         particle_set=particle_set,&
         para_env=para_env,pw_env=pw_env,&
         virial=virial,error=error)
    
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
    IF (use_virial) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Virial not implemented for LRI method", &
            error=error, error_level=cp_failure_level)
    END IF

    eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
    map_consistent=dft_control%qs_control%map_consistent
    IF (map_consistent) THEN
       eps_gvg_rspace = dft_control%qs_control%eps_rho_rspace 
    ELSE
       eps_gvg_rspace = dft_control%qs_control%eps_gvg_rspace
    ENDIF

    DO ikind=1,nkind

       atomic_kind => atomic_kind_set(ikind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            natom=natom_of_kind,&
            atom_list=atom_list,&
            lri_basis_set=lri_basis_set)
       CALL get_gto_basis_set(gto_basis_set=lri_basis_set,&
            first_sgf=first_sgfa,&
            lmax=la_max,&
            lmin=la_min,&
            maxco=maxco,&
            maxsgf_set=maxsgf_set,&
            npgf=npgfa,&
            nset=nseta,&
            nsgf=nsgfa,&  
            nsgf_set=nsgf_seta,&
            pgf_radius=rpgfa,&
            set_radius=set_radius_a,&
            sphi=sphi_a,&
            zet=zeta)
 
       CALL reallocate(hab,1,maxco,1,1)
       CALL reallocate (cores, 1, natom_of_kind)
       int_res(ikind)%v_int = 0.0_dp
       hab = 0._dp
       npme = 0
       cores = 0

       distributed_grids = .FALSE.
       parallel = .FALSE.
       DO igrid_level = 1, gridlevel_info%ngrid_levels
          IF (rs_descs(igrid_level)%rs_desc%distributed) distributed_grids = .TRUE.
          IF (rs_descs(igrid_level)%rs_desc%parallel) parallel = .TRUE.
       ENDDO

      DO iatom = 1, natom_of_kind
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

      DO j=1,npme

        iatom = cores(j)
        atom_a = atom_list(iatom)
        ra(:) = pbc(particle_set(atom_a)%r,cell)

        DO iset = 1,nseta
           sgfa = first_sgfa(1,iset)
           ncoa = npgfa(iset)*ncoset(la_max(iset))
           hab(:,1)=0._dp
           CALL reallocate(work,1,nsgf_seta(iset),1,1)
           work =0.0_dp
           
           DO ipgf = 1,npgfa(iset)
   
              na1 = (ipgf-1)*ncoset(la_max(iset))
              igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf,iset))

              CALL integrate_pgf_product_rspace(la_max=la_max(iset),&
                   zeta=zeta(ipgf,iset),la_min=la_min(iset),&
                   lb_max=0,zetb=0.0_dp,lb_min=0,&
                   ra=ra,rab=(/0.0_dp,0.0_dp,0.0_dp/),rab2=0.0_dp,&
                   rsgrid=rs_v(igrid_level)%rs_grid ,cell=cell,&
                   cube_info=pw_env%cube_info(igrid_level),&
                   hab=hab,o1=na1,o2=0,&
                   eps_gvg_rspace=eps_rho_rspace,calculate_forces=.FALSE.,&
                   map_consistent=map_consistent,&
                   error=error)

           ENDDO
           !contract hab here!!
           CALL dgemm("T","N",nsgf_seta(iset),1,ncoa, 1.0_dp,sphi_a(1,sgfa),&
                SIZE(sphi_a,1),hab(1,1),SIZE(hab,1), 0.0_dp,&
                work(1,1),SIZE(work,1))

           int_res(ikind)%v_int(iatom,sgfa:sgfa-1+nsgf_seta(iset))=work(1:nsgf_seta(iset),1)
        ENDDO

        DEALLOCATE (work, STAT=stat)
        CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
      END DO

      CALL mp_sum(int_res(ikind)%v_int,para_env%group)
     
      DEALLOCATE (cores,hab,STAT=stat)
      CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO

    DO i=1,SIZE(rs_v)
      CALL rs_grid_release(rs_v(i)%rs_grid, error=error)
    END DO
 
    CALL timestop(handle)

  END SUBROUTINE integrate_v_rspace_one_center

! *****************************************************************************
!> \brief computes matrix elements corresponding to a given potential
!> \note
!>     integrates a given potential (or other object on a real
!>     space grid) = v_rspace using a multi grid technique (mgrid_*)
!>     over the basis set producing a number for every element of h
!>     (should have the same sparsity structure of S)
!>     additional screening is available using the magnitude of the
!>     elements in p (? I'm not sure this is a very good idea)
!>     this argument is optional
!>     derivatives of these matrix elements with respect to the ionic
!>     coordinates can be computed as well
!> \par History
!>      IAB (29-Apr-2010): Added OpenMP parallelisation to task loop
!>                         (c) The Numerical Algorithms Group (NAG) Ltd, 2010 on behalf of the HECToR project
! *****************************************************************************
  SUBROUTINE integrate_v_rspace(v_rspace, p, h,qs_env,calculate_forces,compute_tau,gapw,&
       basis_set_id, pw_env_external, task_list_external, error)

    TYPE(pw_p_type)                          :: v_rspace
    TYPE(cp_dbcsr_p_type), INTENT(IN), &
      OPTIONAL                               :: p
    TYPE(cp_dbcsr_p_type), INTENT(INOUT)     :: h
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    LOGICAL, INTENT(IN), OPTIONAL            :: compute_tau, gapw
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    TYPE(pw_env_type), OPTIONAL, POINTER     :: pw_env_external
    TYPE(task_list_type), OPTIONAL, POINTER  :: task_list_external
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_v_rspace', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: add_tasks = 1000, &
                                                max_tasks = 3000
    REAL(kind=dp), PARAMETER                 :: mult_tasks = 2.0_dp

    INTEGER :: atom_a, atom_b, bcol, brow, handle, i, iatom, idir, &
      igrid_level, ikind, ikind_old, ilevel, ipair, ipgf, ipgf_new, iset, &
      iset_new, iset_old, itask, ithread, jatom, jkind, jkind_old, jpgf, &
      jpgf_new, jset, jset_new, jset_old, maxco, maxpgf, maxset, maxsgf_set, &
      my_basis_set_id, na1, na2, natom, nb1, nb2, ncoa, ncob, nkind, nseta, &
      nsetb, nthread, offs_dv, sgfa, sgfb, stat
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: atom_pair_recv, atom_pair_send
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: block_touched
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL :: atom_pair_changed, atom_pair_done, distributed_grids, failure, &
      found, h_duplicated, had_thread_dist, map_consistent, my_compute_tau, &
      my_gapw, new_set_pair_coming, p_duplicated, pab_required, scatter, &
      use_subpatch, use_virial
    REAL(KIND=dp)                            :: dab, eps_gvg_rspace, rab2, &
                                                zetp
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra, rab, &
                                                rab_inv, rb
    REAL(KIND=dp), DIMENSION(3, 3)           :: my_virial_a, my_virial_b
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dist_ab, h_block, hab, &
                                                p_block, pab, rpgfa, rpgfb, &
                                                sphi_a, sphi_b, work, zeta, &
                                                zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: habt, hadb, hdab, pabt, workt
    REAL(kind=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: hadbt, hdabt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_r_p_type), DIMENSION(3)       :: dv_block
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ddv
    TYPE(cp_dbcsr_type), POINTER             :: deltap, dh
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(dbcsr_distribution_obj), POINTER    :: dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(realspace_grid_desc_p_type), &
      DIMENSION(:), POINTER                  :: rs_descs
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_v
    TYPE(section_vals_type), POINTER         :: input, interp_section
    TYPE(task_list_type), POINTER            :: task_list, task_list_soft
    TYPE(virial_type), POINTER               :: virial

!$  INTEGER :: omp_get_max_threads, omp_get_thread_num

    failure=.FALSE.
    NULLIFY(pw_env, rs_descs, tasks, dist_ab)

    debug_count=debug_count+1

    offs_dv=0
    my_compute_tau = .FALSE.
    my_gapw = .FALSE.
    IF (PRESENT(compute_tau)) my_compute_tau = compute_tau
    IF (PRESENT(gapw)) my_gapw = gapw
    IF (PRESENT(basis_set_id)) THEN
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
         input=input,&
         task_list=task_list,&
         task_list_soft=task_list_soft,&
         force=force,pw_env=pw_env,&
         virial=virial,error=error)
    CASE (use_aux_fit_basis_set)
      CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         cell=cell,&
         dft_control=dft_control,&
         particle_set=particle_set,&
         para_env=para_env,&
         input=input,&
         task_list_aux_fit=task_list,&
         task_list_soft=task_list_soft,&
         force=force,pw_env=pw_env,&
         virial=virial,error=error)
    END SELECT

    ! maybe don't even pass qs_env ???
    IF (PRESENT(pw_env_external)) pw_env=>pw_env_external
    IF (PRESENT(task_list_external)) task_list=>task_list_external

    IF (my_compute_tau) THEN
       CALL timeset(routineN,handle)
    ELSE
       CALL timeset(routineN,handle)
    END IF

    ! get the task lists
    IF (my_gapw) task_list=>task_list_soft
    CPPrecondition(ASSOCIATED(task_list),cp_failure_level,routineP,error,failure)
    tasks  =>task_list%tasks
    dist_ab=>task_list%dist_ab
    atom_pair_send=>task_list%atom_pair_send
    atom_pair_recv=>task_list%atom_pair_recv

    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_v, error=error)
    DO i=1,SIZE(rs_v)
      CALL rs_grid_retain(rs_v(i)%rs_grid,error=error)
    END DO

    ! *** assign from pw_env
    gridlevel_info=>pw_env%gridlevel_info
    cube_info=>pw_env%cube_info

    interp_section => section_vals_get_subs_vals(input,"DFT%MGRID%INTERPOLATOR",&
         error=error)
    CALL potential_pw2rs(rs_v,v_rspace,pw_env,interp_section,error)

    !   *** having the potential on the rs_multigrids, just integrate ...
    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    IF (calculate_forces) THEN
       ALLOCATE (atom_of_kind(natom),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
            atom_of_kind=atom_of_kind)
    END IF

    map_consistent=dft_control%qs_control%map_consistent
    IF (map_consistent) THEN
       eps_gvg_rspace = dft_control%qs_control%eps_rho_rspace ! needs to be consistent with rho_rspace
    ELSE
       eps_gvg_rspace = dft_control%qs_control%eps_gvg_rspace
    ENDIF

    pab_required = PRESENT(p) .AND. (calculate_forces .OR. .NOT. map_consistent)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxco=maxco,&
         maxsgf_set=maxsgf_set,&
         basis_set_id=my_basis_set_id)

    distributed_grids = .FALSE.
    DO igrid_level = 1, gridlevel_info%ngrid_levels
       IF ( rs_v(igrid_level)%rs_grid%desc%distributed ) THEN
          distributed_grids = .TRUE.
       ENDIF
    ENDDO

    h_duplicated = .FALSE.
    dh => h%matrix
    IF (distributed_grids) THEN
       NULLIFY ( dh )
       h_duplicated = .TRUE.
       ALLOCATE(dh)
       CALL cp_dbcsr_init(dh, error=error)
       CALL cp_dbcsr_create(dh, 'LocalH', &
            cp_dbcsr_distribution (h%matrix),&
            cp_dbcsr_get_matrix_type (h%matrix), cp_dbcsr_row_block_sizes(h%matrix),&
            cp_dbcsr_col_block_sizes(h%matrix), cp_dbcsr_get_num_blocks(h%matrix),&
            cp_dbcsr_get_data_size(h%matrix),&
            error=error)

    END IF

    p_duplicated = .FALSE.
    IF ( pab_required ) THEN
       deltap => p%matrix
       IF (distributed_grids) THEN
          p_duplicated = .TRUE.
          NULLIFY ( deltap )
          ALLOCATE(deltap)
          CALL cp_dbcsr_init(deltap, error=error)
          CALL cp_dbcsr_copy(deltap,p%matrix,name="LocalP",error=error)
       END IF
    END IF

    nthread = 1
!$  nthread = omp_get_max_threads()

    !   *** Allocate work storage ***
    NULLIFY ( pabt, habt, workt )
    CALL reallocate(habt,1,maxco,1,maxco,0,nthread)
    CALL reallocate(workt,1,maxco,1,maxsgf_set,0,nthread)
    IF (pab_required) THEN
       CALL reallocate(pabt,1,maxco,1,maxco,0,nthread)
    ELSE
       IF (calculate_forces) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Need p for forces")
       END IF
    END IF

    NULLIFY(hdabt,hadbt,hdab,hadb)

    !   get maximum numbers
    natom = SIZE( particle_set )
    maxset=0
    maxpgf=0
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       SELECT CASE (my_basis_set_id)
       CASE (use_orb_basis_set)
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
              softb = my_gapw, &
              orb_basis_set=orb_basis_set)
       CASE (use_aux_fit_basis_set)
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
              softb = my_gapw, &
              aux_fit_basis_set=orb_basis_set,&
              basis_set_id=my_basis_set_id)
       END SELECT

       IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
            npgf=npgfa, nset=nseta )

       maxset=MAX(nseta,maxset)
       maxpgf=MAX(MAXVAL(npgfa),maxpgf)
    END DO

    IF (distributed_grids .AND. pab_required) THEN
        CALL rs_distribute_matrix (rs_descs, deltap, atom_pair_send, atom_pair_recv, natom, scatter=.TRUE., error=error)
    ENDIF

    IF (debug_this_module) &
      ALLOCATE(block_touched(natom,natom))

!$omp parallel default(none), &
!$omp shared(workt,habt,hdabt,hadbt,pabt,tasks,particle_set,natom,maxset), &
!$omp shared(maxpgf,my_basis_set_id,my_gapw,dh,ddv,deltap,use_virial), &
!$omp shared(pab_required,calculate_forces,ncoset,rs_v,cube_info,my_compute_tau), &
!$omp shared(map_consistent,eps_gvg_rspace,force,virial,cell,atom_of_kind,dist_ab), &
!$omp shared(gridlevel_info,task_list,failure,block_touched,nthread), &
!$omp private(ithread,work,hab,hdab,hadb,pab,iset_old,jset_old), &
!$omp private(ikind_old,jkind_old,iatom,jatom,iset,jset,ikind,jkind,ilevel,ipgf,jpgf), &
!$omp private(brow,bcol,orb_basis_set,first_sgfa,la_max,la_min,npgfa,nseta,nsgfa), &
!$omp private(rpgfa,set_radius_a,sphi_a,zeta,first_sgfb,lb_max,lb_min,npgfb), &
!$omp private(nsetb,nsgfb,rpgfb,set_radius_b,sphi_b,zetb,found,error,atom_a,atom_b), &
!$omp private(force_a,force_b,my_virial_a,my_virial_b,atom_pair_changed,h_block), &
!$omp private(dv_block,p_block,ncoa,sgfa,ncob,sgfb,rab,rab2,ra,rb,zetp,dab,igrid_level), &
!$omp private(na1,na2,nb1,nb2,use_subpatch,rab_inv,new_set_pair_coming,atom_pair_done), &
!$omp private(iset_new,jset_new,ipgf_new,jpgf_new,idir,dist), &
!$omp private(had_thread_dist,itask)

    ithread = 0
!$  ithread = omp_get_thread_num()
    work => workt(:,:,ithread)
    hab => habt(:,:,ithread)
    IF (pab_required) THEN
       pab => pabt(:,:,ithread)
    END IF

    iset_old = -1 ; jset_old = -1
    ikind_old = -1 ; jkind_old = -1


    ! Here we loop over gridlevels first, finalising the matrix after each grid level is
    ! completed.  On each grid level, we loop over atom pairs, which will only access
    ! a single block of each matrix, so with OpenMP, each matrix block is only touched
    ! by a single thread for each grid level
    loop_gridlevels: DO igrid_level = 1, gridlevel_info%ngrid_levels

    CALL cp_dbcsr_work_create(dh,work_mutable=.TRUE.,n=nthread,error=error)
!$  dist => dh%matrix%m%dist
!$  CALL cp_assert (dbcsr_distribution_has_threads(dist), cp_fatal_level,&
!$       cp_internal_error, routineN, "No thread distribution defined.",&
!$       error=error)
!$omp barrier

    IF (debug_this_module) THEN
!$omp single
      block_touched = -1
!$omp end single
!$omp flush
    END IF

!$omp do schedule (dynamic, MAX(1,task_list%npairs(igrid_level)/(nthread*50)))
    loop_pairs: DO ipair = 1, task_list%npairs(igrid_level)
    loop_tasks: DO itask = task_list%taskstart(ipair,igrid_level), task_list%taskstop(ipair,igrid_level)

       CALL int2pair(tasks(3,itask),ilevel,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)

       ! At the start of a block of tasks, get atom data (and kind data, if needed)
       IF (itask .EQ. task_list%taskstart(ipair,igrid_level) ) THEN

          ikind = particle_set(iatom)%atomic_kind%kind_number
          jkind = particle_set(jatom)%atomic_kind%kind_number

          ra(:) = pbc(particle_set(iatom)%r,cell)

          IF (iatom <= jatom) THEN
             brow = iatom
             bcol = jatom
          ELSE
             brow = jatom
             bcol = iatom
          END IF

          IF (ikind .NE. ikind_old ) THEN
             SELECT CASE (my_basis_set_id)
             CASE (use_orb_basis_set)
               CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                    softb = my_gapw, &
                    orb_basis_set=orb_basis_set)
             CASE (use_aux_fit_basis_set)
               CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                    softb = my_gapw, &
                    aux_fit_basis_set=orb_basis_set,&
                    basis_set_id = my_basis_set_id)
             END SELECT

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
             SELECT CASE (my_basis_set_id)
             CASE (use_orb_basis_set)
               CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind,&
                    softb = my_gapw, &
                    orb_basis_set=orb_basis_set)
             CASE (use_aux_fit_basis_set)
               CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind,&
                    softb = my_gapw, &
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
                  pgf_radius=rpgfb,&
                  set_radius=set_radius_b,&
                  sphi=sphi_b,&
                  zet=zetb)

          ENDIF

          IF (debug_this_module) THEN
!$omp critical (block_touched_critical)
            IF ((block_touched(brow,bcol).NE.ithread) .AND. (block_touched(brow,bcol).NE. -1) ) THEN
              CALL stop_program(routineN,moduleN,__LINE__,&
                                "Block has been modified by another thread")
            END IF
            block_touched(brow,bcol) = ithread
!$omp end critical (block_touched_critical)
          END IF

          NULLIFY(h_block)
          CALL cp_dbcsr_get_block_p(dh,brow,bcol,h_block,found)
          IF (.NOT.ASSOCIATED(h_block)) THEN
               CALL cp_dbcsr_add_block_node ( dh, brow, bcol, h_block ,error=error)
          END IF

          IF (pab_required) THEN
             CALL cp_dbcsr_get_block_p(matrix=deltap,&
                  row=brow,col=bcol,BLOCK=p_block,found=found)

             IF (.NOT.ASSOCIATED(p_block)) THEN
                CALL stop_program(routineN,moduleN,__LINE__,&
                                  "p_block not associated in deltap")
             END IF
          END IF

          IF (calculate_forces) THEN
             atom_a = atom_of_kind(iatom)
             atom_b = atom_of_kind(jatom)
             force_a(:) = 0.0_dp
             force_b(:) = 0.0_dp
          ENDIF
          IF (use_virial) THEN
             my_virial_a = 0.0_dp
             my_virial_b = 0.0_dp
          ENDIF

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
          IF (pab_required) THEN
             IF (iatom <= jatom) THEN
                CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                     1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                     p_block(sgfa,sgfb),SIZE(p_block,1),&
                     0.0_dp,work(1,1),SIZE(work,1))
                CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                     1.0_dp,work(1,1),SIZE(work,1),&
                     sphi_b(1,sgfb),SIZE(sphi_b,1),&
                     0.0_dp,pab(1,1),SIZE(pab,1))
             ELSE
                CALL dgemm("N","N",ncob,nsgfa(iset),nsgfb(jset),&
                     1.0_dp,sphi_b(1,sgfb),SIZE(sphi_b,1),&
                     p_block(sgfb,sgfa),SIZE(p_block,1),&
                     0.0_dp,work(1,1),SIZE(work,1))
                CALL dgemm("N","T",ncob,ncoa,nsgfa(iset),&
                     1.0_dp,work(1,1),SIZE(work,1),&
                     sphi_a(1,sgfa),SIZE(sphi_a,1),&
                     0.0_dp,pab(1,1),SIZE(pab,1))
             END IF
          END IF

          IF (iatom<=jatom) THEN
             hab(1:ncoa,1:ncob) = 0._dp
          ELSE
             hab(1:ncob,1:ncoa) = 0._dp
          ENDIF

          iset_old = iset
          jset_old = jset

       ENDIF

       rab(1) = dist_ab (1,itask)
       rab(2) = dist_ab (2,itask)
       rab(3) = dist_ab (3,itask)
       rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       rb(1) = ra(1) + rab(1)
       rb(2) = ra(2) + rab(2)
       rb(3) = ra(3) + rab(3)
       zetp = zeta(ipgf,iset) + zetb(jpgf,jset)
       dab=SQRT(rab2)

       na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
       na2 = ipgf*ncoset(la_max(iset))
       nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
       nb2 = jpgf*ncoset(lb_max(jset))

       ! check whether we need to use fawzi's generalised collocation scheme
       IF(rs_v(igrid_level)%rs_grid%desc%distributed)THEN
          !tasks(4,:) is 0 for replicated, 1 for distributed 2 for exceptional distributed tasks
          IF (tasks(4,itask) .EQ. 2 ) THEN
             use_subpatch = .TRUE.
          ELSE
             use_subpatch = .FALSE.
          ENDIF
       ELSE
          use_subpatch = .FALSE.
       ENDIF

       IF (pab_required) THEN
          IF (iatom <= jatom) THEN
             CALL integrate_pgf_product_rspace(&
                  la_max(iset),zeta(ipgf,iset),la_min(iset),&
                  lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                  ra,rab,rab2,rs_v(igrid_level)%rs_grid,cell,&
                  cube_info(igrid_level),&
                  hab,pab=pab,o1=na1-1,o2=nb1-1, &
                  eps_gvg_rspace=eps_gvg_rspace,&
                  calculate_forces=calculate_forces,&
                  force_a=force_a,force_b=force_b,&
                  compute_tau=my_compute_tau,map_consistent=map_consistent,&
                  use_virial=use_virial,my_virial_a=my_virial_a,&
                  my_virial_b=my_virial_b,use_subpatch=use_subpatch,subpatch_pattern=tasks(6,itask),error=error)
          ELSE
             rab_inv=-rab
             CALL integrate_pgf_product_rspace(&
                  lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                  la_max(iset),zeta(ipgf,iset),la_min(iset),&
                  rb,rab_inv,rab2,rs_v(igrid_level)%rs_grid,cell,&
                  cube_info(igrid_level),&
                  hab,pab=pab,o1=nb1-1,o2=na1-1, &
                  eps_gvg_rspace=eps_gvg_rspace,&
                  calculate_forces=calculate_forces,&
                  force_a=force_b,force_b=force_a,&
                  compute_tau=my_compute_tau,map_consistent=map_consistent,&
                  use_virial=use_virial,my_virial_a=my_virial_b,&
                  my_virial_b=my_virial_a,use_subpatch=use_subpatch,subpatch_pattern=tasks(6,itask),error=error)
          END IF
       ELSE
          IF (iatom <= jatom) THEN
             CALL integrate_pgf_product_rspace(&
                  la_max(iset),zeta(ipgf,iset),la_min(iset),&
                  lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                  ra,rab,rab2,rs_v(igrid_level)%rs_grid,cell,&
                  cube_info(igrid_level),&
                  hab,o1=na1-1,o2=nb1-1,&
                  eps_gvg_rspace=eps_gvg_rspace,&
                  calculate_forces=calculate_forces,&
                  force_a=force_a,force_b=force_b,&
                  compute_tau=my_compute_tau,&
                  map_consistent=map_consistent,use_subpatch=use_subpatch,subpatch_pattern=tasks(6,itask),error=error)
          ELSE
             rab_inv=-rab
             CALL integrate_pgf_product_rspace(&
                  lb_max(jset),zetb(jpgf,jset),lb_min(jset),&
                  la_max(iset),zeta(ipgf,iset),la_min(iset),&
                  rb,rab_inv,rab2,rs_v(igrid_level)%rs_grid,cell,&
                  cube_info(igrid_level),&
                  hab,o1=nb1-1,o2=na1-1,&
                  eps_gvg_rspace=eps_gvg_rspace,&
                  calculate_forces=calculate_forces,&
                  force_a=force_b,force_b=force_a, &
                  compute_tau=my_compute_tau,&
                  map_consistent=map_consistent,use_subpatch=use_subpatch,subpatch_pattern=tasks(6,itask),error=error)
          END IF
       END IF

       new_set_pair_coming=.FALSE.
       atom_pair_done = .FALSE.
       IF (itask < task_list%taskstop(ipair,igrid_level))  THEN
          CALL int2pair(tasks(3,itask+1),ilevel,iatom,jatom,iset_new,jset_new,ipgf_new,jpgf_new,natom,maxset,maxpgf)
          IF (iset_new .NE. iset .OR. jset_new .NE. jset) THEN
             new_set_pair_coming=.TRUE.
          ENDIF
       ELSE
          ! do not forget the last block
          new_set_pair_coming=.TRUE.
          atom_pair_done = .TRUE.
       ENDIF

       ! contract the block into h if we're done with the current set pair
       IF (new_set_pair_coming) THEN
          IF (iatom <= jatom) THEN
             CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                  1.0_dp,hab(1,1),SIZE(hab,1),&
                  sphi_b(1,sgfb),SIZE(sphi_b,1),&
                  0.0_dp,work(1,1),SIZE(work,1))
             CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                  1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                  work(1,1),SIZE(work,1),&
                  1.0_dp,h_block(sgfa,sgfb),SIZE(h_block,1))
          ELSE
             CALL dgemm("N","N",ncob,nsgfa(iset),ncoa,&
                  1.0_dp,hab(1,1),SIZE(hab,1),&
                  sphi_a(1,sgfa),SIZE(sphi_a,1),&
                  0.0_dp,work(1,1),SIZE(work,1))
             CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncob,&
                  1.0_dp,sphi_b(1,sgfb),SIZE(sphi_b,1),&
                  work(1,1),SIZE(work,1),&
                  1.0_dp,h_block(sgfb,sgfa),SIZE(h_block,1))
          END IF
       END IF

       IF (atom_pair_done) THEN
!$omp critical(force_critical)
          IF (calculate_forces) THEN
             force(ikind)%rho_elec(:,atom_a) =&
                  force(ikind)%rho_elec(:,atom_a) + 2.0_dp*force_a(:)
             IF (iatom /= jatom ) THEN
                force(jkind)%rho_elec(:,atom_b) =&
                     force(jkind)%rho_elec(:,atom_b) + 2.0_dp*force_b(:)
             END IF
          ENDIF
          IF (use_virial) THEN
             IF (use_virial .AND. calculate_forces) THEN
                virial%pv_virial = virial%pv_virial + 2.0_dp*my_virial_a
                IF (iatom /= jatom) THEN
                   virial%pv_virial = virial%pv_virial + 2.0_dp*my_virial_b
                END IF
             END IF
          END IF
!$omp end critical (force_critical)
       ENDIF
    END DO loop_tasks
    END DO loop_pairs
!$omp end do

    CALL cp_dbcsr_finalize(dh, error=error)

    END DO loop_gridlevels

!$omp end parallel

    IF (debug_this_module) &
      DEALLOCATE(block_touched)

    IF ( h_duplicated ) THEN
       ! Reconstruct H matrix if using distributed RS grids
       ! note send and recv direction reversed WRT collocate
       scatter = .FALSE.
       CALL rs_distribute_matrix (rs_descs, dh, atom_pair_recv, atom_pair_send,&
            natom, scatter, error, h%matrix)
       CALL cp_dbcsr_deallocate_matrix ( dh ,error=error)

    ELSE
       NULLIFY ( dh, ddv )
    END IF

    IF ( p_duplicated ) THEN
       CALL cp_dbcsr_deallocate_matrix ( deltap ,error=error)
    ELSE
       NULLIFY ( deltap )
    END IF

    !   *** Release work storage ***

    DEALLOCATE (habt,workt,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    IF ( pab_required ) THEN
       DEALLOCATE (pabt,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    IF (ASSOCIATED(rs_v)) THEN
      DO i=1,SIZE(rs_v)
        CALL rs_grid_release(rs_v(i)%rs_grid, error=error)
      END DO
    END IF

    IF (calculate_forces) THEN
       DEALLOCATE (atom_of_kind,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE integrate_v_rspace

! *****************************************************************************
!> \brief transfers a potential from a pw_grid to a vector of
!>      realspace multigrids
!> \param v_rspace INPUT : the potential on a planewave grid in Rspace
!> \param rs_v OUTPUT: the potential on the realspace multigrids
!> \note
!>      extracted from integrate_v_rspace
!>      should contain all parallel communication of integrate_v_rspace in the
!>      case of replicated grids.
!> \par History
!>      09.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE potential_pw2rs(rs_v,v_rspace,pw_env,interp_section,error)

    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_v
    TYPE(pw_p_type), INTENT(IN)              :: v_rspace
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(section_vals_type), POINTER         :: interp_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'potential_pw2rs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: auxbas_grid, handle, &
                                                igrid_level, interp_kind
    REAL(KIND=dp)                            :: scale
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: mgrid_gspace, mgrid_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools

    CALL timeset(routineN,handle)

    ! *** set up of the potential on the multigrids
    CALL pw_env_get(pw_env, pw_pools=pw_pools, gridlevel_info=gridlevel_info, &
           auxbas_grid = auxbas_grid, error=error)

    CALL pw_pools_create_pws(pw_pools,mgrid_rspace,&
                use_data = REALDATA3D,&
                in_space = REALSPACE, error=error)

    ! use either realspace or fft techniques to get the potential on the rs multigrids
    CALL section_vals_val_get(interp_section,"KIND",i_val=interp_kind,error=error)
    SELECT CASE(interp_kind)
    CASE (pw_interp)
       CALL pw_pools_create_pws(pw_pools,mgrid_gspace,&
                                 use_data = COMPLEXDATA1D,&
                                 in_space = RECIPROCALSPACE, error=error)
       CALL pw_transfer(v_rspace%pw,mgrid_gspace(auxbas_grid)%pw,error=error)
       DO igrid_level=1,gridlevel_info%ngrid_levels
         IF ( igrid_level /= auxbas_grid ) THEN
              CALL pw_copy(mgrid_gspace(auxbas_grid)%pw,mgrid_gspace(igrid_level)%pw,&
                   error=error)
              CALL pw_transfer(mgrid_gspace(igrid_level)%pw,mgrid_rspace(igrid_level)%pw,&
                   error=error)
         ELSE
              IF (mgrid_gspace(auxbas_grid)%pw%pw_grid%spherical) THEN
                  CALL pw_transfer(mgrid_gspace(auxbas_grid)%pw,mgrid_rspace(auxbas_grid)%pw,&
                       error=error)
              ELSE ! fft forward + backward should be identical
                  CALL pw_copy(v_rspace%pw,mgrid_rspace(auxbas_grid)%pw,error=error)
              ENDIF
         ENDIF
         ! *** Multiply by the grid volume element ratio ***
         IF ( igrid_level /= auxbas_grid ) THEN
            scale = mgrid_rspace(igrid_level)%pw%pw_grid%dvol/&
                    mgrid_rspace(auxbas_grid)%pw%pw_grid%dvol
            mgrid_rspace(igrid_level)%pw%cr3d = &
                                      scale*mgrid_rspace(igrid_level)%pw%cr3d
         END IF
       END DO
       CALL pw_pools_give_back_pws(pw_pools,mgrid_gspace,error=error)
    CASE(spline3_pbc_interp)
       CALL pw_copy(v_rspace%pw,mgrid_rspace(1)%pw,error=error)
       DO igrid_level=1,gridlevel_info%ngrid_levels-1
          CALL pw_zero(mgrid_rspace(igrid_level+1)%pw,error=error)
          CALL pw_restrict_s3(mgrid_rspace(igrid_level)%pw,&
               mgrid_rspace(igrid_level+1)%pw,pw_pools(igrid_level+1)%pool,&
               interp_section,error=error)
          ! *** Multiply by the grid volume element ratio
          mgrid_rspace(igrid_level+1) % pw % cr3d = &
                 mgrid_rspace(igrid_level+1) % pw % cr3d * 8._dp
       END DO
    CASE default
       CALL cp_unimplemented_error(routineN,"interpolation not supported "//&
            cp_to_string(interp_kind),error=error)
    END SELECT

    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL rs_pw_transfer(rs_v(igrid_level)%rs_grid,&
                           mgrid_rspace(igrid_level)%pw,pw2rs,error=error)
    ENDDO
    ! *** give back the pw multi-grids
    CALL pw_pools_give_back_pws(pw_pools,mgrid_rspace,error=error)

    CALL timestop(handle)

  END SUBROUTINE potential_pw2rs

END MODULE qs_integrate_potential
