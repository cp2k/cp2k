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
MODULE qs_integrate_potential_single
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE atprop_types,                    ONLY: atprop_array_init
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
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
  USE pw_types,                        ONLY: pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_integrate_potential_low,      ONLY: integrate_pgf_product_rspace
  USE realspace_grid_types,            ONLY: &
       pw2rs, realspace_grid_desc_p_type, realspace_grid_desc_type, &
       realspace_grid_p_type, realspace_grid_type, rs_grid_release, &
       rs_grid_retain, rs_grid_zero, rs_pw_transfer
  USE rs_pw_interface,                 ONLY: potential_pw2rs
  USE scptb_types,                     ONLY: get_scptb_parameter,&
                                             scp_vector_type,&
                                             scptb_parameter_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_integrate_potential_single'

! *** Public subroutines ***
! *** Don't include this routines directly, use the interface to
! *** qs_integrate_potential

  PUBLIC :: integrate_v_rspace_one_center,&
            integrate_v_core_rspace,&
            integrate_ppl_rspace,&
            integrate_scp_rspace,&
            integrate_rho_nlcc

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
  SUBROUTINE integrate_v_rspace_one_center(v_rspace,qs_env,int_res,&
                                           calculate_forces,error)
    TYPE(pw_p_type), INTENT(INOUT)           :: v_rspace
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(lri_kind_type), DIMENSION(:), &
      POINTER                                :: int_res
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'integrate_v_rspace_one_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, dir, group_size, handle, i, iatom, igrid_level, ikind, &
      ipgf, iset, m1, maxco, maxsgf_set, my_pos, na1, natom_of_kind, ncoa, &
      nkind, nseta, offset, sgfa, stat
    INTEGER, DIMENSION(3)                    :: lb, location, tp, ub
    INTEGER, DIMENSION(:), POINTER           :: atom_list, la_max, la_min, &
                                                npgfa, nsgf_seta
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: failure, map_consistent, &
                                                map_it_here, use_virial
    REAL(KIND=dp)                            :: eps_gvg_rspace, eps_rho_rspace
    REAL(KIND=dp), DIMENSION(3)              :: force_a, force_b, ra
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab, pab, rpgfa, sphi_a, &
                                                work_f, work_i, zeta
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
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    NULLIFY(atomic_kind, atomic_kind_set, atom_list, cell, dft_control,&
            first_sgfa, gridlevel_info, hab, la_max, la_min, lri_basis_set,&
            npgfa, nsgf_seta, pab, para_env, particle_set, pw_env, rpgfa, &
            rs_descs, rs_grid, rs_v, virial, set_radius_a, sphi_a,  work_f,&
            work_i, zeta)

    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)

    CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_v, error=error)
    DO i=1,SIZE(rs_v)
      CALL rs_grid_retain(rs_v(i)%rs_grid,error=error)
      CALL rs_grid_zero(rs_v(i)%rs_grid)
    END DO

    gridlevel_info=>pw_env%gridlevel_info
    
    CALL potential_pw2rs(rs_v,v_rspace,pw_env,error)

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

    offset = 0           
    my_pos=v_rspace%pw%pw_grid%para%my_pos           
    group_size=v_rspace%pw%pw_grid%para%group_size   

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
            nsgf_set=nsgf_seta,&
            pgf_radius=rpgfa,&
            set_radius=set_radius_a,&
            sphi=sphi_a,&
            zet=zeta)
 
       CALL reallocate(hab,1,maxco,1,1)
       CALL reallocate(pab,1,maxco,1,1)
       int_res(ikind)%v_int = 0.0_dp
       hab = 0._dp
       pab(:,1) = 0._dp

       DO iatom = 1,natom_of_kind

          atom_a = atom_list(iatom)
          ra(:) = pbc(particle_set(atom_a)%r,cell)
          force_a(:) = 0._dp
          force_b(:) = 0._dp

          DO iset = 1,nseta
             sgfa = first_sgfa(1,iset)
             ncoa = npgfa(iset)*ncoset(la_max(iset))
             hab(:,1) = 0._dp
             CALL reallocate(work_i,1,nsgf_seta(iset),1,1)
             work_i = 0.0_dp

             ! get fit coefficients for forces
             IF (calculate_forces) THEN
              m1=sgfa+nsgf_seta(iset)-1
              CALL reallocate(work_f,1,nsgf_seta(iset),1,1)
              work_f(1:nsgf_seta(iset),1) = int_res(ikind)%acoef(iatom,sgfa:m1)  
              CALL dgemm("N","N",ncoa,1,nsgf_seta(iset),1.0_dp,sphi_a(1,sgfa),&
                   SIZE(sphi_a,1),work_f(1,1),SIZE(work_f,1),0.0_dp,pab(1,1),&
                   SIZE(pab,1))
             ENDIF
             
             DO ipgf = 1,npgfa(iset)
                na1 = (ipgf-1)*ncoset(la_max(iset))
                igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf,iset))
                rs_grid => rs_v(igrid_level)%rs_grid          
  
                map_it_here=.FALSE.

                IF (.NOT. ALL (rs_grid%desc%perd == 1)) THEN
                   DO dir = 1,3
                      ! bounds of local grid (i.e. removing the 'wings'), if periodic
                      tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir,:),ra)*rs_grid%desc%npts(dir))
                      tp(dir) = MODULO ( tp(dir), rs_grid%desc%npts(dir) )
                      IF (rs_grid%desc%perd(dir) .NE. 1) THEN
                         lb(dir) = rs_grid%lb_local ( dir ) + rs_grid%desc%border
                         ub(dir) = rs_grid%ub_local ( dir ) - rs_grid%desc%border
                      ELSE
                         lb(dir) = rs_grid%lb_local ( dir )
                         ub(dir) = rs_grid%ub_local ( dir )
                      ENDIF
                      ! distributed grid, only map if it is local to the grid
                      location(dir)=tp(dir)+rs_grid%desc%lb(dir)
                   ENDDO
                   IF (lb(1)<=location(1) .AND. location(1)<=ub(1) .AND. &
                       lb(2)<=location(2) .AND. location(2)<=ub(2) .AND. &
                       lb(3)<=location(3) .AND. location(3)<=ub(3)) THEN
                      map_it_here=.TRUE.
                   ENDIF
                ELSE
                   ! not distributed, just a round-robin distribution over the full set of CPUs
                   IF (MODULO(offset,group_size)==my_pos) map_it_here=.TRUE.
                ENDIF

                IF (map_it_here) THEN
                   IF(.NOT.calculate_forces) THEN
                      CALL integrate_pgf_product_rspace(la_max=la_max(iset),&
                           zeta=zeta(ipgf,iset),la_min=la_min(iset),&
                           lb_max=0,zetb=0.0_dp,lb_min=0,&
                           ra=ra,rab=(/0.0_dp,0.0_dp,0.0_dp/),rab2=0.0_dp,&
                           rsgrid=rs_grid ,cell=cell,&
                           cube_info=pw_env%cube_info(igrid_level),&
                           hab=hab,o1=na1,o2=0, eps_gvg_rspace=eps_rho_rspace,&
                           calculate_forces=calculate_forces,&
                           map_consistent=map_consistent,&
                           error=error)
                   ELSE
                      CALL integrate_pgf_product_rspace(la_max=la_max(iset),&
                           zeta=zeta(ipgf,iset),la_min=la_min(iset),&
                           lb_max=0,zetb=0.0_dp,lb_min=0,&
                           ra=ra,rab=(/0.0_dp,0.0_dp,0.0_dp/),rab2=0.0_dp,&
                           rsgrid=rs_grid ,cell=cell,&
                           cube_info=pw_env%cube_info(igrid_level),&
                           hab=hab,pab=pab,o1=na1,o2=0, eps_gvg_rspace=eps_rho_rspace,&
                           calculate_forces=calculate_forces,force_a=force_a,&
                           force_b=force_b,map_consistent=map_consistent,&
                           error=error)
                   ENDIF
                ENDIF               
             ENDDO
             ! contract hab 
             CALL dgemm("T","N",nsgf_seta(iset),1,ncoa, 1.0_dp,sphi_a(1,sgfa),&
                  SIZE(sphi_a,1),hab(1,1),SIZE(hab,1), 0.0_dp,work_i(1,1),SIZE(work_i,1))

             int_res(ikind)%v_int(iatom,sgfa:sgfa-1+nsgf_seta(iset))=work_i(1:nsgf_seta(iset),1)

             offset=offset+nsgf_seta(iset) 

          ENDDO

          IF(calculate_forces) int_res(ikind)%v_dfdr(iatom,:)=force_a(:)

          DEALLOCATE (work_i, STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          IF(calculate_forces) THEN
           DEALLOCATE (work_f, STAT=stat)
           CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF

       END DO

      CALL mp_sum(int_res(ikind)%v_int,para_env%group)
     
      DEALLOCATE (hab,pab,STAT=stat)
      CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO

    DO i=1,SIZE(rs_v)
      CALL rs_grid_release(rs_v(i)%rs_grid, error=error)
    END DO
 
    CALL timestop(handle)

  END SUBROUTINE integrate_v_rspace_one_center

END MODULE qs_integrate_potential_single
