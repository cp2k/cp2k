!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH (11 May 2001) : cleaning up of support structures
!>      CJM & HAF (27 July 2001): fixed bug with handling of cutoff larger than
!>                                half the boxsize.
!>      07.02.2005: getting rid of scaled_to_real calls in force loop (MK)
!>      22.06.2013: OpenMP parallelisation of pair interaction loop (MK)
!> \author CJM
! *****************************************************************************
MODULE fist_nonbond_force

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE atprop_types,                    ONLY: atprop_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE mathconstants,                   ONLY: oorootpi,&
                                             sqrthalf
  USE mathlib,                         ONLY: matvec_3x3
  USE message_passing,                 ONLY: mp_sum
  USE pair_potential_coulomb,          ONLY: potential_coulomb
  USE pair_potential_types,            ONLY: nosh_nosh,&
                                             nosh_sh,&
                                             pair_potential_pp_type,&
                                             pair_potential_single_type,&
                                             sh_sh,&
                                             siepmann_type,&
                                             tersoff_type
  USE particle_types,                  ONLY: particle_type
  USE shell_potential_types,           ONLY: get_shell,&
                                             shell_kind_type
  USE splines_methods,                 ONLY: potential_s
  USE splines_types,                   ONLY: spline_data_p_type,&
                                             spline_factor_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_nonbond_force'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

  PUBLIC :: force_nonbond,&
            bonded_correct_gaussian

CONTAINS

! *****************************************************************************
!> \brief Calculates the force and the potential of the minimum image, and
!>      the pressure tensor
!> \param fist_nonbond_env ...
!> \param ewald_env ...
!> \param particle_set ...
!> \param cell ...
!> \param pot_nonbond ...
!> \param f_nonbond ...
!> \param pv_nonbond ...
!> \param fshell_nonbond ...
!> \param fcore_nonbond ...
!> \param atprop_env ...
!> \param atomic_kind_set ...
!> \param use_virial ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_nonbond(fist_nonbond_env, ewald_env, particle_set, cell, &
       pot_nonbond, f_nonbond, pv_nonbond, fshell_nonbond, fcore_nonbond, &
       atprop_env, atomic_kind_set, use_virial, error)

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(OUT)               :: pot_nonbond
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_nonbond, pv_nonbond
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: fshell_nonbond, fcore_nonbond
    TYPE(atprop_type), POINTER               :: atprop_env
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set(:)
    LOGICAL, INTENT(IN)                      :: use_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_nonbond', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, ewald_type, handle, i, iend, igrp, ikind, &
      ilist, ipair, istart, istat, j, kind_a, kind_b, nkind, npairs, shell_a, &
      shell_b, shell_type
    INTEGER, DIMENSION(:, :), POINTER        :: list
    LOGICAL                                  :: all_terms, do_multipoles, &
                                                failure, full_nl, &
                                                shell_present
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: is_shell_kind
    REAL(KIND=dp) :: alpha, beta, beta_a, beta_b, energy, etot, fac_ei, &
      fac_kind, fac_vdw, fscalar, mm_radius_a, mm_radius_b, qcore_a, qcore_b, &
      qeff_a, qeff_b, qshell_a, qshell_b, rab2, rab2_com, rab2_max
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: mm_radius, qcore, qeff, qshell
    REAL(KIND=dp), DIMENSION(3) :: cell_v, cvi, fatom_a, fatom_b, fcore_a, &
      fcore_b, fshell_a, fshell_b, rab, rab_cc, rab_com, rab_cs, rab_sc, &
      rab_ss
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv, pv_thread
    REAL(KIND=dp), DIMENSION(3, 4)           :: rab_list
    REAL(KIND=dp), DIMENSION(4)              :: rab2_list
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ij_kind_full_fac
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: ei_interaction_cutoffs
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_pp_type), POINTER    :: potparm, potparm14
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc, &
                                                rcore_last_update_pbc, &
                                                rshell_last_update_pbc
    TYPE(shell_kind_type), POINTER           :: shell_kind
    TYPE(spline_data_p_type), DIMENSION(:), &
      POINTER                                :: spline_data
    TYPE(spline_factor_type), POINTER        :: spl_f

    CALL timeset(routineN,handle)
    NULLIFY (logger)
    logger => cp_error_get_logger(error)
    NULLIFY(pot, rshell_last_update_pbc, spl_f, ij_kind_full_fac)
    CALL fist_nonbond_env_get(fist_nonbond_env, nonbonded=nonbonded, &
         potparm14=potparm14, potparm=potparm, r_last_update=r_last_update, &
         r_last_update_pbc=r_last_update_pbc,natom_types=nkind, &
         rshell_last_update_pbc=rshell_last_update_pbc, &
         rcore_last_update_pbc=rcore_last_update_pbc, &
         ij_kind_full_fac=ij_kind_full_fac, error=error)
    CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type, &
         do_multipoles=do_multipoles, &
         interaction_cutoffs=ei_interaction_cutoffs)

    ! Initializing the potential energy, pressure tensor and force
    pot_nonbond    = 0.0_dp
    f_nonbond(:,:) = 0.0_dp

    IF (use_virial) THEN
       pv_nonbond(:,:) = 0.0_dp
    END IF
    shell_present = .FALSE.
    IF (PRESENT(fshell_nonbond)) THEN
       CPPostcondition(PRESENT(fcore_nonbond),cp_failure_level,routineP,error,failure)
       fshell_nonbond = 0.0_dp
       fcore_nonbond = 0.0_dp
       shell_present = .TRUE.
    END IF
    ! Load atomic kind information
    ALLOCATE (mm_radius(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"mm_radius",dp_size*nkind)
    ALLOCATE (qeff(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"qeff",dp_size*nkind)
    ALLOCATE (qcore(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"qcore",dp_size*nkind)
    ALLOCATE (qshell(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"qshell",dp_size*nkind)
    ALLOCATE (is_shell_kind(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"is_shell_kind",int_size*nkind)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind,&
                            qeff=qeff(ikind),&
                            mm_radius=mm_radius(ikind),&
                            shell=shell_kind)
       is_shell_kind(ikind) = ASSOCIATED(shell_kind)
       IF (ASSOCIATED(shell_kind)) THEN
          CALL get_shell(shell=shell_kind,&
                         charge_core=qcore(ikind),&
                         charge_shell=qshell(ikind),&
                         error=error)
       ELSE
          qcore(ikind) = 0.0_dp
          qshell(ikind) = 0.0_dp
       END IF
    END DO
    ! Starting the force loop
    Lists: DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       npairs = neighbor_kind_pair%npairs
       IF (npairs == 0) CYCLE
       list  => neighbor_kind_pair%list
       cvi   =  neighbor_kind_pair%cell_vector
       CALL matvec_3x3(cell_v,cell%hmat,cvi)
       Kind_Group_Loop: DO igrp=1,neighbor_kind_pair%ngrp_kind
          istart = neighbor_kind_pair%grp_kind_start(igrp)
          iend   = neighbor_kind_pair%grp_kind_end(igrp)
          !$omp parallel default(none) &
          !$omp          private(ipair,atom_a,atom_b,kind_a,kind_b,fac_kind,pot) &
          !$omp          private(fac_ei,fac_vdw,atomic_kind,full_nl,qcore_a,qshell_a) &
          !$omp          private(qeff_a,qcore_b,qshell_b,qeff_b,mm_radius_a,mm_radius_b) &
          !$omp          private(shell_kind,beta,beta_a,beta_b,spl_f,spline_data) &
          !$omp          private(shell_type,all_terms,rab_cc,rab_cs,rab_sc,rab_ss) &
          !$omp          private(rab_list,rab2_list,rab_com,rab2_com,pv,pv_thread) &
          !$omp          private(rab,rab2,rab2_max,fscalar,energy,error,failure) &
          !$omp          private(shell_a,shell_b,etot,fatom_a,fatom_b) &
          !$omp          private(fcore_a,fcore_b,fshell_a,fshell_b,i,j) &
          !$omp          shared(shell_present) &
          !$omp          shared(istart,iend,list,particle_set,ij_kind_full_fac) &
          !$omp          shared(neighbor_kind_pair,atomic_kind_set,fist_nonbond_env) &
          !$omp          shared(potparm,potparm14,do_multipoles,r_last_update_pbc) &
          !$omp          shared(use_virial,ei_interaction_cutoffs,alpha,cell_v) &
          !$omp          shared(rcore_last_update_pbc,rshell_last_update_pbc) &
          !$omp          shared(f_nonbond,fcore_nonbond,fshell_nonbond,logger) &
          !$omp          shared(ewald_type,pot_nonbond,pv_nonbond,atprop_env) &
          !$omp          shared(is_shell_kind,mm_radius,qcore,qeff,qshell)
          IF (use_virial) pv_thread(:,:) = 0.0_dp
          !$omp do
          Pairs: DO ipair=istart,iend
             atom_a = list(1,ipair)
             atom_b = list(2,ipair)
             ! Get actual atomic kinds, since atom_a is not always of
             ! kind_a and atom_b of kind_b, ie. they might be swapped.
             kind_a = particle_set(atom_a)%atomic_kind%kind_number
             kind_b = particle_set(atom_b)%atomic_kind%kind_number

             fac_kind = ij_kind_full_fac(kind_a,kind_b)
             ! take the proper potential
             pot => potparm%pot(kind_a,kind_b)%pot
             IF (ipair <= neighbor_kind_pair%nscale) THEN
                IF (neighbor_kind_pair%is_onfo(ipair)) THEN
                   pot => potparm14%pot(kind_a,kind_b)%pot
                END IF
             END IF

             ! Determine the scaling factors
             fac_ei = fac_kind
             fac_vdw = fac_kind
             full_nl = ANY(pot%type == tersoff_type).OR.ANY(pot%type == siepmann_type)
             IF ((.NOT.full_nl).AND.(atom_a == atom_b)) THEN
                fac_ei = 0.5_dp*fac_ei
                fac_vdw = 0.5_dp*fac_vdw
             END IF
             ! decide which interactions to compute
             IF (do_multipoles) THEN
                fac_ei = 0.0_dp
             END IF
             IF (ipair <= neighbor_kind_pair%nscale) THEN
                fac_ei = fac_ei*neighbor_kind_pair%ei_scale(ipair)
                fac_vdw = fac_vdw*neighbor_kind_pair%vdw_scale(ipair)
             END IF

             IF (fac_ei > 0.0_dp) THEN
                ! Get the electrostatic parameters for the atoms a and b
                mm_radius_a = mm_radius(kind_a)
                mm_radius_b = mm_radius(kind_b)
                IF (ASSOCIATED(fist_nonbond_env%charges)) THEN
                   qeff_a = fist_nonbond_env%charges(atom_a)
                   qeff_b = fist_nonbond_env%charges(atom_b)
                ELSE
                   qeff_a = qeff(kind_a)
                   qeff_b = qeff(kind_b)
                END IF
                IF (is_shell_kind(kind_a)) THEN
                   qcore_a = qcore(kind_a)
                   qshell_a = qshell(kind_a)
                   IF ((qcore_a == 0.0_dp).AND.(qshell_a == 0.0_dp)) fac_ei = 0.0_dp
                ELSE
                   qcore_a = qeff_a
                   qshell_a = HUGE(0.0_dp)
                   IF (qeff_a == 0.0_dp) fac_ei = 0.0_dp
                END IF
                IF (is_shell_kind(kind_b)) THEN
                   qcore_b = qcore(kind_b)
                   qshell_b = qshell(kind_b)
                   IF ((qcore_b == 0.0_dp).AND.(qshell_b == 0.0_dp)) fac_ei = 0.0_dp
                ELSE
                   qcore_b = qeff_b
                   qshell_b = HUGE(0.0_dp)
                   IF (qeff_b == 0.0_dp) fac_ei = 0.0_dp
                END IF
                ! Derive beta parameters
                beta = 0.0_dp
                beta_a = 0.0_dp
                beta_b = 0.0_dp
                IF (mm_radius_a > 0) THEN
                   beta_a = sqrthalf/mm_radius_a
                END IF
                IF (mm_radius_b > 0) THEN
                   beta_b = sqrthalf/mm_radius_b
                END IF
                IF ((mm_radius_a > 0) .OR. (mm_radius_b > 0)) THEN
                   beta = sqrthalf/SQRT(mm_radius_a*mm_radius_a + mm_radius_b*mm_radius_b)
                END IF
             END IF

             ! In case we have only manybody potentials and no charges, this
             ! pair of atom types can be ignored here.
             IF (pot%no_pp.AND.(fac_ei == 0.0)) CYCLE

             ! Setup spline_data set
             spl_f => pot%spl_f
             spline_data => pot%pair_spline_data
             shell_type = pot%shell_type
             IF (shell_type /= nosh_nosh) THEN
                CPPrecondition(.NOT.do_multipoles,cp_failure_level,routineP,error,failure)
                CPPostcondition(shell_present,cp_failure_level,routineP,error,failure)
             END IF
             rab2_max = pot%rcutsq

             ! compute the relative vector(s) for this pair
             IF (shell_type /= nosh_nosh) THEN
                ! do shell
                all_terms = .TRUE.
                IF (shell_type == sh_sh) THEN
                   shell_a = particle_set(atom_a)%shell_index
                   shell_b = particle_set(atom_b)%shell_index
                   rab_cc =  rcore_last_update_pbc(shell_b)%r -  rcore_last_update_pbc(shell_a)%r
                   rab_cs = rshell_last_update_pbc(shell_b)%r -  rcore_last_update_pbc(shell_a)%r
                   rab_sc =  rcore_last_update_pbc(shell_b)%r - rshell_last_update_pbc(shell_a)%r
                   rab_ss = rshell_last_update_pbc(shell_b)%r - rshell_last_update_pbc(shell_a)%r
                   rab_list(1:3,1) = rab_cc(1:3) + cell_v(1:3)
                   rab_list(1:3,2) = rab_cs(1:3) + cell_v(1:3)
                   rab_list(1:3,3) = rab_sc(1:3) + cell_v(1:3)
                   rab_list(1:3,4) = rab_ss(1:3) + cell_v(1:3)
                ELSE IF ((shell_type == nosh_sh).AND.(particle_set(atom_a)%shell_index /= 0)) THEN
                   shell_a = particle_set(atom_a)%shell_index
                   shell_b = 0
                   rab_cc = r_last_update_pbc(atom_b)%r -  rcore_last_update_pbc(shell_a)%r
                   rab_sc = 0.0_dp
                   rab_cs = 0.0_dp
                   rab_ss = r_last_update_pbc(atom_b)%r - rshell_last_update_pbc(shell_a)%r
                   rab_list(1:3,1) = rab_cc(1:3) + cell_v(1:3)
                   rab_list(1:3,2) = 0.0_dp
                   rab_list(1:3,3) = 0.0_dp
                   rab_list(1:3,4) = rab_ss(1:3) + cell_v(1:3)
                ELSE IF ((shell_type == nosh_sh).AND.(particle_set(atom_b)%shell_index /= 0)) THEN
                   shell_b = particle_set(atom_b)%shell_index
                   shell_a = 0
                   rab_cc =  rcore_last_update_pbc(shell_b)%r - r_last_update_pbc(atom_a)%r
                   rab_sc = 0.0_dp
                   rab_cs = 0.0_dp
                   rab_ss = rshell_last_update_pbc(shell_b)%r - r_last_update_pbc(atom_a)%r
                   rab_list(1:3,1) = rab_cc(1:3) + cell_v(1:3)
                   rab_list(1:3,2) = 0.0_dp
                   rab_list(1:3,3) = 0.0_dp
                   rab_list(1:3,4) = rab_ss(1:3) + cell_v(1:3)
                END IF
                ! Compute the term only if all the pairs (cc,cs,sc,ss) are within the cut-off
                Check_terms: DO i=1,4
                   rab2_list(i) = rab_list(1,i)**2 + rab_list(2,i)**2 + rab_list(3,i)**2
                   IF (rab2_list(i) >= rab2_max) THEN
                      all_terms = .FALSE.
                      EXIT Check_terms
                   END IF
                END DO Check_terms
                rab_com = r_last_update_pbc(atom_b)%r - r_last_update_pbc(atom_a)%r
             ELSE
                ! not do shell
                rab_cc  = r_last_update_pbc(atom_b)%r - r_last_update_pbc(atom_a)%r
                rab_com = rab_cc
                shell_a = 0
                shell_b = 0
             END IF
             rab_com  = rab_com + cell_v
             rab2_com = rab_com(1)**2 + rab_com(2)**2 + rab_com(3)**2

             ! compute the interactions for the current pair
             etot = 0.0_dp
             fatom_a(:) = 0.0_dp
             fatom_b(:) = 0.0_dp
             fcore_a(:) = 0.0_dp
             fcore_b(:) = 0.0_dp
             fshell_a(:) = 0.0_dp
             fshell_b(:) = 0.0_dp
             IF (use_virial) pv(:,:) = 0.0_dp
             IF (shell_type /= nosh_nosh) THEN
                ! do shell
                IF ((rab2_com <= rab2_max).AND.all_terms) THEN
                   IF (fac_ei > 0) THEN
                      ! core-core or core-ion/ion-core: Coulomb only
                      rab = rab_list(:,1)
                      rab2 = rab2_list(1)
                      fscalar = 0.0_dp
                      IF (shell_a == 0) THEN
                         ! atom a is a plain ion and can have beta_a > 0
                         energy = potential_coulomb(rab2,fscalar,fac_ei*qeff_a*qcore_b,&
                                                    ewald_type,alpha,beta_a,&
                                                    ei_interaction_cutoffs(2,kind_a,kind_b))
                         CALL add_force_nonbond(fatom_a,fcore_b,pv,fscalar,rab,use_virial)
                      ELSE IF (shell_b == 0) THEN
                         ! atom b is a plain ion and can have beta_b > 0
                         energy = potential_coulomb(rab2,fscalar,fac_ei*qcore_a*qeff_b,&
                                                    ewald_type,alpha,beta_b,&
                                                    ei_interaction_cutoffs(2,kind_b,kind_a))
                         CALL add_force_nonbond(fcore_a,fatom_b,pv,fscalar,rab,use_virial)
                      ELSE
                         ! core-core interaction is always pure point charge
                         energy = potential_coulomb(rab2,fscalar,fac_ei*qcore_a*qcore_b,&
                                                    ewald_type,alpha,0.0_dp,&
                                                    ei_interaction_cutoffs(1,kind_a,kind_b))
                         CALL add_force_nonbond(fcore_a,fcore_b,pv,fscalar,rab,use_virial)
                      END IF
                      etot = etot + energy
                   END IF

                   IF (shell_type == sh_sh) THEN
                      ! shell-shell: VDW + Coulomb
                      rab  = rab_list(:,4)
                      rab2 = rab2_list(4)
                      fscalar = 0.0_dp
                      IF (fac_vdw > 0) THEN
                         energy = potential_s(spline_data,rab2,fscalar,spl_f,logger)
                         etot = etot + energy*fac_vdw
                         fscalar = fscalar*fac_vdw
                      END IF
                      IF (fac_ei > 0) THEN
                         ! note that potential_coulomb increments fscalar
                         energy = potential_coulomb(rab2,fscalar,fac_ei*qshell_a*qshell_b,&
                                                    ewald_type,alpha,beta,&
                                                    ei_interaction_cutoffs(3,kind_a,kind_b))
                         etot = etot + energy
                      END IF
                      CALL add_force_nonbond(fshell_a,fshell_b,pv,fscalar,rab,use_virial)

                      IF (fac_ei > 0) THEN
                         ! core-shell: Coulomb only
                         rab  = rab_list(:,2)
                         rab2 = rab2_list(2)
                         fscalar = 0.0_dp
                         ! swap kind_a and kind_b to get the right cutoff
                         energy = potential_coulomb(rab2,fscalar,fac_ei*qcore_a*qshell_b,&
                                                    ewald_type,alpha,beta_b,&
                                                    ei_interaction_cutoffs(2,kind_b,kind_a))
                         etot = etot + energy
                         CALL add_force_nonbond(fcore_a,fshell_b,pv,fscalar,rab,use_virial)

                         ! shell-core: Coulomb only
                         rab  = rab_list(:,3)
                         rab2 = rab2_list(3)
                         fscalar = 0.0_dp
                         energy = potential_coulomb(rab2,fscalar,fac_ei*qshell_a*qcore_b,&
                                                    ewald_type,alpha,beta_a,&
                                                    ei_interaction_cutoffs(2,kind_a,kind_b))
                         etot = etot + energy
                         CALL add_force_nonbond(fshell_a,fcore_b,pv,fscalar,rab,use_virial)
                      END IF
                   ELSE IF ((shell_type == nosh_sh).AND.(shell_a == 0)) THEN
                      ! ion-shell: VDW + Coulomb
                      rab = rab_list(:,4)
                      rab2 = rab2_list(4)
                      fscalar = 0.0_dp
                      IF (fac_vdw > 0) THEN
                         energy = potential_s(spline_data,rab2,fscalar,spl_f,logger)
                         etot = etot + energy*fac_vdw
                         fscalar = fscalar*fac_vdw
                      END IF
                      IF (fac_ei > 0) THEN
                         ! note that potential_coulomb increments fscalar
                         energy = potential_coulomb(rab2,fscalar,fac_ei*qeff_a*qshell_b,&
                                                    ewald_type,alpha,beta,&
                                                    ei_interaction_cutoffs(3,kind_a,kind_b))
                         etot = etot + energy
                      END IF
                      CALL add_force_nonbond(fatom_a,fshell_b,pv,fscalar,rab,use_virial)
                   ELSE IF ((shell_type == nosh_sh) .AND. (shell_b == 0)) THEN
                      ! shell-ion : VDW + Coulomb
                      rab = rab_list(:,4)
                      rab2 = rab2_list(4)
                      fscalar = 0.0_dp
                      IF (fac_vdw>0) THEN
                         energy = potential_s(spline_data,rab2,fscalar,spl_f,logger)
                         etot = etot + energy*fac_vdw
                         fscalar = fscalar*fac_vdw
                      END IF
                      IF (fac_ei > 0) THEN
                         ! note that potential_coulomb increments fscalar
                         energy = potential_coulomb(rab2,fscalar,fac_ei*qshell_a*qeff_b,&
                                                    ewald_type,alpha,beta,&
                                                    ei_interaction_cutoffs(3,kind_a,kind_b))
                         etot = etot + energy
                      END IF
                      CALL add_force_nonbond(fshell_a,fatom_b,pv,fscalar,rab,use_virial)
                   END IF
                END IF
             ELSE
                IF (rab2_com <= rab2_max) THEN
                   ! NO SHELL MODEL...
                   ! Ion-Ion: no shell model, VDW + coulomb
                   rab  = rab_com
                   rab2 = rab2_com
                   fscalar = 0.0_dp
                   IF (fac_vdw > 0) THEN
                      energy = potential_s(spline_data,rab2,fscalar,spl_f,logger)
                      etot = etot + energy*fac_vdw
                      fscalar = fscalar*fac_vdw
                   END IF
                   IF (fac_ei > 0) THEN
                      ! note that potential_coulomb increments fscalar
                      energy = potential_coulomb(rab2,fscalar,fac_ei*qeff_a*qeff_b,&
                                                 ewald_type,alpha,beta,&
                                                 ei_interaction_cutoffs(3,kind_a,kind_b))
                      etot = etot + energy
                   END IF
                   CALL add_force_nonbond(fatom_a,fatom_b,pv,fscalar,rab,use_virial)
                END IF
             END IF
             ! Nonbonded energy
             !$omp atomic
             pot_nonbond = pot_nonbond + etot
             IF (atprop_env%energy) THEN
                ! Update atomic energies
                !$omp atomic
                atprop_env%atener(atom_a) = atprop_env%atener(atom_a) + 0.5_dp*etot
                !$omp atomic
                atprop_env%atener(atom_b) = atprop_env%atener(atom_b) + 0.5_dp*etot
             END IF
             ! Nonbonded forces
             DO i=1,3
                !$omp atomic
                f_nonbond(i,atom_a) = f_nonbond(i,atom_a) + fatom_a(i)
                !$omp atomic
                f_nonbond(i,atom_b) = f_nonbond(i,atom_b) + fatom_b(i)
             END DO
             IF (shell_a > 0) THEN
                DO i=1,3
                   !$omp atomic
                   fcore_nonbond(i,shell_a) = fcore_nonbond(i,shell_a) + fcore_a(i)
                   !$omp atomic
                   fshell_nonbond(i,shell_a) = fshell_nonbond(i,shell_a) + fshell_a(i)
                END DO
             END IF
             IF (shell_b > 0) THEN
                DO i=1,3
                   !$omp atomic
                   fcore_nonbond(i,shell_b) = fcore_nonbond(i,shell_b) + fcore_b(i)
                   !$omp atomic
                   fshell_nonbond(i,shell_b) = fshell_nonbond(i,shell_b) + fshell_b(i)
                END DO
             END IF
             ! Add the contribution of the current pair to the total pressure tensor
             IF (use_virial) THEN
                DO i=1,3
                   DO j=1,3
                      pv_thread(j,i) = pv_thread(j,i) + pv(j,i)
                      ! Update atomic stress tensors
                      IF (atprop_env%stress) THEN
                         !$omp atomic
                         atprop_env%atstress(j,i,atom_a) = atprop_env%atstress(j,i,atom_a) +&
                                                           0.5_dp*pv(j,i)
                         !$omp atomic
                         atprop_env%atstress(j,i,atom_b) = atprop_env%atstress(j,i,atom_b) +&
                                                           0.5_dp*pv(j,i)
                      END IF
                   END DO
                END DO
             END IF
          END DO Pairs 
          !$omp end do
          IF (use_virial) THEN
             DO i=1,3
                DO j=1,3
                   !$omp atomic
                   pv_nonbond(j,i) = pv_nonbond(j,i) + pv_thread(j,i)
                END DO
             END DO
          END IF
          !$omp end parallel
       END DO Kind_Group_Loop
    END DO Lists

    DEALLOCATE (mm_radius,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"mm_radius")
    DEALLOCATE (qeff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"qeff")
    DEALLOCATE (qcore,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"qcore")
    DEALLOCATE (qshell,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"qshell")
    DEALLOCATE (is_shell_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"is_shell_kind")

    CALL timestop(handle)

  END SUBROUTINE force_nonbond

  ! *****************************************************************************
  !> \brief Adds a non-bonding contribution to the total force and optionally to
  !>        the virial.
  ! *****************************************************************************
! *****************************************************************************
!> \brief ...
!> \param f_nonbond_a ...
!> \param f_nonbond_b ...
!> \param pv ...
!> \param fscalar ...
!> \param rab ...
!> \param use_virial ...
! *****************************************************************************
  SUBROUTINE add_force_nonbond(f_nonbond_a,f_nonbond_b,pv,fscalar,rab,use_virial)

    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: f_nonbond_a, f_nonbond_b
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(INOUT)                          :: pv
    REAL(KIND=dp), INTENT(IN)                :: fscalar
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    LOGICAL, INTENT(IN)                      :: use_virial

    REAL(KIND=dp), DIMENSION(3)              :: fr

    fr(1) = fscalar*rab(1)
    fr(2) = fscalar*rab(2)
    fr(3) = fscalar*rab(3)
    f_nonbond_a(1) = f_nonbond_a(1) - fr(1)
    f_nonbond_a(2) = f_nonbond_a(2) - fr(2)
    f_nonbond_a(3) = f_nonbond_a(3) - fr(3)
    f_nonbond_b(1) = f_nonbond_b(1) + fr(1)
    f_nonbond_b(2) = f_nonbond_b(2) + fr(2)
    f_nonbond_b(3) = f_nonbond_b(3) + fr(3)
    IF (use_virial) THEN
       pv(1,1) = pv(1,1) + rab(1)*fr(1)
       pv(1,2) = pv(1,2) + rab(1)*fr(2)
       pv(1,3) = pv(1,3) + rab(1)*fr(3)
       pv(2,1) = pv(2,1) + rab(2)*fr(1)
       pv(2,2) = pv(2,2) + rab(2)*fr(2)
       pv(2,3) = pv(2,3) + rab(2)*fr(3)
       pv(3,1) = pv(3,1) + rab(3)*fr(1)
       pv(3,2) = pv(3,2) + rab(3)*fr(2)
       pv(3,3) = pv(3,3) + rab(3)*fr(3)
    END IF

  END SUBROUTINE

! *****************************************************************************
!> \brief corrects electrostatics for bonded terms
!> \param fist_nonbond_env ...
!> \param atomic_kind_set ...
!> \param local_particles ...
!> \param particle_set ...
!> \param ewald_env ...
!> \param v_bonded_corr ...
!> \param pv_bc ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param atprop_env ...
!> \param cell ...
!> \param use_virial ...
!> \param error ...
!> \par History
!>      Splitted routines to clean and to fix a bug with the tensor whose
!>      original definition was not correct for PBC.. [Teodoro Laino -06/2007]
! *****************************************************************************
  SUBROUTINE bonded_correct_gaussian(fist_nonbond_env, atomic_kind_set, &
             local_particles, particle_set, ewald_env, v_bonded_corr, pv_bc, &
             shell_particle_set, core_particle_set, atprop_env, cell, use_virial, error)

    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    REAL(KIND=dp), INTENT(OUT)               :: v_bonded_corr
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: pv_bc
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set( : ), &
                                                core_particle_set( : )
    TYPE(atprop_type), POINTER               :: atprop_env
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN)                      :: use_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'bonded_correct_gaussian', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, group, handle, iatom, iend, igrp, ilist, &
      ipair, istart, kind_a, kind_b, natoms_per_kind, nkind, npairs, shell_a, &
      shell_b
    INTEGER, DIMENSION(:, :), POINTER        :: list
    LOGICAL                                  :: a_is_shell, b_is_shell, &
                                                do_multipoles, full_nl, &
                                                shell_adiabatic
    REAL(KIND=dp)                            :: alpha, const, fac_cor, &
                                                fac_ei, qcore_a, qcore_b, &
                                                qeff_a, qeff_b, qshell_a, &
                                                qshell_b
    REAL(KIND=dp), DIMENSION(3)              :: rca, rcb, rsa, rsb
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ij_kind_full_fac
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pair_potential_pp_type), POINTER    :: potparm, potparm14
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(shell_kind_type), POINTER           :: shell_kind

    CALL timeset(routineN,handle)

    ! Initializing values
    IF (use_virial) pv_bc = 0.0_dp
    v_bonded_corr = 0.0_dp

    CALL fist_nonbond_env_get(fist_nonbond_env, nonbonded=nonbonded, &
         potparm14=potparm14, potparm=potparm, &
         ij_kind_full_fac=ij_kind_full_fac, error=error)
    CALL ewald_env_get(ewald_env, alpha=alpha, do_multipoles=do_multipoles, &
         group=group, error=error)
    ! Defining the constants
    const = 2.0_dp*alpha*oorootpi

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             shell_adiabatic=shell_adiabatic)


    Lists: DO ilist=1,nonbonded%nlists
       neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
       npairs = neighbor_kind_pair%nscale
       IF (npairs==0) CYCLE
       list => neighbor_kind_pair%list
       Kind_Group_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
          istart = neighbor_kind_pair%grp_kind_start(igrp)
          IF (istart > npairs) THEN
             EXIT
          END IF
          iend = MIN(npairs,neighbor_kind_pair%grp_kind_end(igrp))

          Pairs: DO ipair = istart, iend
             atom_a = list(1,ipair)
             atom_b = list(2,ipair)
             ! Get actual atomic kinds, since atom_a is not always of
             ! kind_a and atom_b of kind_b, ie. they might be swapped.
             kind_a = particle_set(atom_a)%atomic_kind%kind_number
             kind_b = particle_set(atom_b)%atomic_kind%kind_number

             ! take the proper potential, only for full_nl test
             pot => potparm%pot(kind_a,kind_b)%pot
             IF (ipair <= neighbor_kind_pair%nscale) THEN
                IF (neighbor_kind_pair%is_onfo(ipair)) THEN
                   pot => potparm14%pot(kind_a,kind_b)%pot
                END IF
             END IF

             ! Determine the scaling factors
             fac_ei = ij_kind_full_fac(kind_a, kind_b)
             full_nl = ANY(pot%type==tersoff_type).OR.ANY(pot%type==siepmann_type)
             IF ((.NOT.full_nl).AND.(atom_a==atom_b)) THEN
                fac_ei = fac_ei*0.5_dp
             END IF
             IF (ipair <= neighbor_kind_pair%nscale) THEN
                fac_ei = fac_ei*neighbor_kind_pair%ei_scale(ipair)
             END IF
             ! The amount of correction is related to the
             ! amount of scaling as follows:
             fac_cor = 1.0_dp - fac_ei
             IF (fac_cor <= 0.0_dp) CYCLE

             ! Parameters for kind a
             atomic_kind => atomic_kind_set(kind_a)
             CALL get_atomic_kind(atomic_kind, qeff=qeff_a, shell=shell_kind)
             IF (ASSOCIATED(fist_nonbond_env%charges)) qeff_a=fist_nonbond_env%charges(atom_a)
             a_is_shell = ASSOCIATED(shell_kind)
             IF (a_is_shell) THEN
                CALL get_shell(shell=shell_kind, charge_core=qcore_a, &
                               charge_shell=qshell_a, error=error)
                shell_a = particle_set(atom_a)%shell_index
                rca = core_particle_set(shell_a)%r
                rsa = shell_particle_set(shell_a)%r
             ELSE
                qcore_a = qeff_a
                qshell_a = HUGE(0.0_dp)
                shell_a = 0
                rca = particle_set(atom_a)%r
                rsa = 0.0_dp
             END IF

             ! Parameters for kind b
             atomic_kind => atomic_kind_set(kind_b)
             CALL get_atomic_kind(atomic_kind, qeff=qeff_b, shell=shell_kind)
             IF (ASSOCIATED(fist_nonbond_env%charges)) qeff_b=fist_nonbond_env%charges(atom_b)
             b_is_shell = ASSOCIATED(shell_kind)
             IF (b_is_shell) THEN
                CALL get_shell(shell=shell_kind, charge_core=qcore_b, &
                               charge_shell=qshell_b, error=error)
                shell_b = particle_set(atom_b)%shell_index
                rcb = core_particle_set(shell_b)%r
                rsb = shell_particle_set(shell_b)%r
             ELSE
                qcore_b = qeff_b
                qshell_b = HUGE(0.0_dp)
                shell_b = 0
                rcb = particle_set(atom_b)%r
                rsb = 0.0_dp
             END IF

             ! First part: take care of core/ion-core/ion correction
             IF (a_is_shell .AND. b_is_shell) THEN
                ! correct for core-core interaction
                CALL bonded_correct_gaussian_low(rca, rcb, cell, &
                     v_bonded_corr, core_particle_set, core_particle_set, &
                     shell_a, shell_b, .TRUE., alpha, qcore_a, qcore_b, &
                     const, fac_cor, pv_bc, atprop_env, use_virial)
             ELSE IF (a_is_shell) THEN
                ! correct for core-ion interaction
                CALL bonded_correct_gaussian_low(rca, rcb, cell, &
                     v_bonded_corr, core_particle_set, particle_set, &
                     shell_a, atom_b, .TRUE., alpha, qcore_a, qcore_b, &
                     const, fac_cor, pv_bc, atprop_env, use_virial)
             ELSE IF (b_is_shell) THEN
                ! correct for ion-core interaction
                CALL bonded_correct_gaussian_low(rca, rcb, cell, &
                     v_bonded_corr, particle_set, core_particle_set, &
                     atom_a, shell_b, .TRUE., alpha, qcore_a, qcore_b, &
                     const, fac_cor, pv_bc, atprop_env, use_virial)
             ELSE
                ! correct for ion-ion interaction
                CALL bonded_correct_gaussian_low(rca, rcb, cell, &
                     v_bonded_corr, particle_set, particle_set, &
                     atom_a, atom_b, .TRUE., alpha, qcore_a, qcore_b, &
                     const, fac_cor, pv_bc, atprop_env, use_virial)
             END IF

             ! Second part: take care of shell-shell, shell-core/ion and
             ! core/ion-shell corrections
             IF (a_is_shell .AND. b_is_shell) THEN
                ! correct for shell-shell interaction
                CALL bonded_correct_gaussian_low(rsa, rsa, cell, &
                     v_bonded_corr, shell_particle_set, shell_particle_set, &
                     shell_a, shell_b, shell_adiabatic, alpha, qshell_a, &
                     qshell_b, const, fac_cor, pv_bc, atprop_env, use_virial)
             END IF
             IF (a_is_shell) THEN
                IF (b_is_shell) THEN
                   ! correct for shell-core interaction
                   CALL bonded_correct_gaussian_low(rsa, rcb, cell, &
                        v_bonded_corr, shell_particle_set, core_particle_set, &
                        shell_a, shell_b, shell_adiabatic, alpha, qshell_a, qcore_b, &
                        const, fac_cor, pv_bc, atprop_env, use_virial)
                ELSE
                   ! correct for shell-ion interaction
                   CALL bonded_correct_gaussian_low(rsa, rcb, cell, &
                        v_bonded_corr, shell_particle_set, particle_set, &
                        shell_a, atom_b, shell_adiabatic, alpha, qshell_a, qcore_b, &
                        const, fac_cor, pv_bc,atprop_env, use_virial)
                END IF
             END IF
             IF (b_is_shell) THEN
                IF (a_is_shell) THEN
                   ! correct for core-shell interaction
                   CALL bonded_correct_gaussian_low(rca, rsb, cell, &
                        v_bonded_corr, core_particle_set, shell_particle_set, &
                        shell_a, shell_b, shell_adiabatic, alpha, qcore_a, qshell_b, &
                        const, fac_cor, pv_bc, atprop_env, use_virial)
                ELSE
                   ! correct for ion-shell interaction
                   CALL bonded_correct_gaussian_low(rca, rsb, cell, &
                        v_bonded_corr, particle_set, shell_particle_set, &
                        atom_a, shell_b, shell_adiabatic, alpha, qcore_a, qshell_b, &
                        const, fac_cor, pv_bc, atprop_env, use_virial)
                END IF
             END IF
          END DO Pairs
       END DO Kind_Group_Loop
    END DO Lists

    ! Always correct core-shell interaction within one atom.
    nkind = SIZE(atomic_kind_set)
    DO kind_a = 1, nkind
       ! parameters for kind a
       atomic_kind => atomic_kind_set(kind_a)
       CALL get_atomic_kind(atomic_kind, shell=shell_kind)
       IF (ASSOCIATED(shell_kind)) THEN
          CALL get_shell(shell=shell_kind, charge_core=qcore_a, &
                         charge_shell=qshell_a, error=error)

          natoms_per_kind = local_particles%n_el(kind_a)
          DO iatom = 1, natoms_per_kind

             ! Data for atom a
             atom_a = local_particles%list(kind_a)%array(iatom)
             shell_a = particle_set(atom_a)%shell_index
             rca = core_particle_set(shell_a)%r
             rsa = shell_particle_set(shell_a)%r

             CALL bonded_correct_gaussian_low_sh(rca, rsa, cell, &
                  v_bonded_corr, core_particle_set, shell_particle_set, &
                  shell_a, shell_adiabatic, alpha, qcore_a, qshell_a, &
                  const, pv_bc, atprop_env, use_virial)

          END DO
       END IF
    END DO

    CALL mp_sum(v_bonded_corr,group)

    CALL timestop(handle)

  END SUBROUTINE bonded_correct_gaussian

! *****************************************************************************
!> \brief ...
!> \param r1 ...
!> \param r2 ...
!> \param cell ...
!> \param v_bonded_corr ...
!> \param particle_set1 ...
!> \param particle_set2 ...
!> \param i ...
!> \param j ...
!> \param shell_adiabatic ...
!> \param alpha ...
!> \param q1 ...
!> \param q2 ...
!> \param const ...
!> \param fac_cor ...
!> \param pv_bc ...
!> \param atprop_env ...
!> \param use_virial ...
!> \par History
!>      Splitted routines to clean and to fix a bug with the tensor whose
!>      original definition was not correct for PBC..
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE bonded_correct_gaussian_low(r1, r2, cell, v_bonded_corr,&
       particle_set1, particle_set2, i, j, shell_adiabatic, alpha, q1, q2, &
       const, fac_cor, pv_bc, atprop_env, use_virial)
    REAL(KIND=dp), DIMENSION(3)              :: r1, r2
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: v_bonded_corr
    TYPE(particle_type), POINTER             :: particle_set1(:), &
                                                particle_set2(:)
    INTEGER, INTENT(IN)                      :: i, j
    LOGICAL, INTENT(IN)                      :: shell_adiabatic
    REAL(KIND=dp), INTENT(IN)                :: alpha, q1, q2, const, fac_cor
    REAL(KIND=dp), INTENT(INOUT)             :: pv_bc(3,3)
    TYPE(atprop_type), POINTER               :: atprop_env
    LOGICAL, INTENT(IN)                      :: use_virial

    CHARACTER(LEN=*), PARAMETER :: routineN = 'bonded_correct_gaussian_low', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER :: ac1 = 0.254829592_dp, ac2 = -0.284496736_dp, &
      ac3 = 1.421413741_dp, ac4 = -1.453152027_dp, ac5 = 1.061405429_dp, &
      pc = 0.3275911_dp

    INTEGER                                  :: iatom, jatom
    REAL(KIND=dp)                            :: arg, dij, e_arg_arg, errf, &
                                                fscalar, idij, rijsq, tc, vbc
    REAL(KIND=dp), DIMENSION(3)              :: fij_com, rij
    REAL(KIND=dp), DIMENSION(3, 3)           :: fbc

    rij = r1 - r2
    rij = pbc(rij,cell)
    rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
    idij = 1.0_dp/SQRT(rijsq)
    dij = rijsq*idij
    arg = alpha*dij
    e_arg_arg = EXP(-arg**2)
    tc = 1.0_dp/(1.0_dp + pc*arg)

    ! Defining errf=1-erfc
    errf = 1.0_dp - ((((ac5*tc + ac4)*tc + ac3)*tc + ac2)*tc + ac1)*tc*e_arg_arg

    ! Getting the potential
    vbc = -q1*q2*idij*errf*fac_cor
    v_bonded_corr = v_bonded_corr + vbc
    IF (atprop_env%energy) THEN
       iatom = particle_set1(i)%atom_index
       atprop_env%atener(iatom) = atprop_env%atener(iatom) + 0.5_dp*vbc
       jatom = particle_set2(j)%atom_index
       atprop_env%atener(jatom) = atprop_env%atener(jatom) + 0.5_dp*vbc
    END IF

    ! Subtracting the force from the total force
    fscalar = q1*q2*idij**2*(idij*errf - const*e_arg_arg)*fac_cor

    particle_set1(i)%f(1) = particle_set1(i)%f(1) - fscalar*rij(1)
    particle_set1(i)%f(2) = particle_set1(i)%f(2) - fscalar*rij(2)
    particle_set1(i)%f(3) = particle_set1(i)%f(3) - fscalar*rij(3)

    particle_set2(j)%f(1) = particle_set2(j)%f(1) + fscalar*rij(1)
    particle_set2(j)%f(2) = particle_set2(j)%f(2) + fscalar*rij(2)
    particle_set2(j)%f(3) = particle_set2(j)%f(3) + fscalar*rij(3)

    IF (use_virial.AND.shell_adiabatic) THEN
       fij_com = fscalar*rij
       fbc(1,1) = -fij_com(1)*rij(1)
       fbc(1,2) = -fij_com(1)*rij(2)
       fbc(1,3) = -fij_com(1)*rij(3)
       fbc(2,1) = -fij_com(2)*rij(1)
       fbc(2,2) = -fij_com(2)*rij(2)
       fbc(2,3) = -fij_com(2)*rij(3)
       fbc(3,1) = -fij_com(3)*rij(1)
       fbc(3,2) = -fij_com(3)*rij(2)
       fbc(3,3) = -fij_com(3)*rij(3)
       pv_bc(:,:) = pv_bc(:,:) + fbc(:,:)
       IF (atprop_env%stress) THEN
          ! Atomic stress tensors
          iatom = particle_set1(i)%atom_index
          atprop_env%atstress(:,:,iatom) = atprop_env%atstress(:,:,iatom) + 0.5_dp*fbc(:,:)
          jatom = particle_set2(j)%atom_index
          atprop_env%atstress(:,:,jatom) = atprop_env%atstress(:,:,jatom) + 0.5_dp*fbc(:,:)
       END IF
    END IF

  END SUBROUTINE bonded_correct_gaussian_low

! *****************************************************************************
!> \brief specific for shell models cleans the interaction core-shell on the same
!>      atom
!> \param r1 ...
!> \param r2 ...
!> \param cell ...
!> \param v_bonded_corr ...
!> \param core_particle_set ...
!> \param shell_particle_set ...
!> \param i ...
!> \param shell_adiabatic ...
!> \param alpha ...
!> \param q1 ...
!> \param q2 ...
!> \param const ...
!> \param pv_bc ...
!> \param atprop_env ...
!> \param use_virial ...
!> \par History
!>      Splitted routines to clean and to fix a bug with the tensor whose
!>      original definition was not correct for PBC..
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE bonded_correct_gaussian_low_sh(r1, r2, cell, v_bonded_corr, &
       core_particle_set, shell_particle_set, i, shell_adiabatic, alpha, q1, q2,&
       const, pv_bc, atprop_env, use_virial)
    REAL(KIND=dp), DIMENSION(3)              :: r1, r2
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: v_bonded_corr
    TYPE(particle_type), POINTER             :: core_particle_set(:), &
                                                shell_particle_set(:)
    INTEGER, INTENT(IN)                      :: i
    LOGICAL, INTENT(IN)                      :: shell_adiabatic
    REAL(KIND=dp), INTENT(IN)                :: alpha, q1, q2, const
    REAL(KIND=dp), INTENT(INOUT)             :: pv_bc(3,3)
    TYPE(atprop_type), POINTER               :: atprop_env
    LOGICAL, INTENT(IN)                      :: use_virial

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'bonded_correct_gaussian_low_sh', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER :: ac1 = 0.254829592_dp, ac2 = -0.284496736_dp, &
      ac3 = 1.421413741_dp, ac4 = -1.453152027_dp, ac5 = 1.061405429_dp, &
      pc = 0.3275911_dp

    INTEGER                                  :: iatom
    REAL(KIND=dp)                            :: arg, dij, e_arg_arg, efac, &
                                                errf, ffac, fscalar, idij, &
                                                rijsq, tc, tc2, tc4, vbc
    REAL(KIND=dp), DIMENSION(3)              :: fr, rij
    REAL(KIND=dp), DIMENSION(3, 3)           :: fbc

    rij = r1 - r2
    rij = pbc(rij,cell)
    rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
    dij = SQRT(rijsq)
    ! Two possible limiting cases according the value of dij
    arg = alpha*dij
    ! and this is a magic number.. it is related to the order expansion
    ! and to the value of the polynomial coefficients
    IF (arg > 0.355_dp) THEN
       idij = 1.0_dp/dij
       e_arg_arg = EXP(-arg*arg)
       tc = 1.0_dp/(1.0_dp + pc*arg)
       ! defining errf = 1 - erfc
       errf = 1.0_dp - ((((ac5*tc+ac4)*tc+ac3)*tc+ac2)*tc+ac1)*tc*e_arg_arg
       efac = idij*errf
       ffac = idij**2*(efac - const*e_arg_arg)
    ELSE
       tc  = arg*arg
       tc2 = tc*tc
       tc4 = tc2*tc2
       efac = const*(1.0_dp - tc/3.0_dp + tc2/10.0_dp - tc*tc2/42.0_dp + tc4/216.0_dp -&
              tc*tc4/1320.0_dp + tc2*tc4/9360.0_dp)
       ffac = const*alpha**2*(2.0_dp/3.0_dp - 2.0_dp*tc/5.0_dp + tc2/7.0_dp - tc*tc2/27.0_dp +&
              tc4/132.0_dp - tc*tc4/780.0_dp)
    END IF

    ! getting the potential
    vbc = -q1*q2*efac
    v_bonded_corr = v_bonded_corr + vbc
    IF (atprop_env%energy) THEN
       iatom = shell_particle_set(i)%atom_index
       atprop_env%atener(iatom) = atprop_env%atener(iatom) + vbc
    END IF

    ! subtracting the force from the total force
    fscalar = q1*q2*ffac
    fr(:) = fscalar*rij(:)

    core_particle_set(i)%f(1) = core_particle_set(i)%f(1) - fr(1)
    core_particle_set(i)%f(2) = core_particle_set(i)%f(2) - fr(2)
    core_particle_set(i)%f(3) = core_particle_set(i)%f(3) - fr(3)

    shell_particle_set(i)%f(1) = shell_particle_set(i)%f(1) + fr(1)
    shell_particle_set(i)%f(2) = shell_particle_set(i)%f(2) + fr(2)
    shell_particle_set(i)%f(3) = shell_particle_set(i)%f(3) + fr(3)

    IF (use_virial.AND.shell_adiabatic) THEN
       fbc(1,1) = -fr(1)*rij(1)
       fbc(1,2) = -fr(1)*rij(2)
       fbc(1,3) = -fr(1)*rij(3)
       fbc(2,1) = -fr(2)*rij(1)
       fbc(2,2) = -fr(2)*rij(2)
       fbc(2,3) = -fr(2)*rij(3)
       fbc(3,1) = -fr(3)*rij(1)
       fbc(3,2) = -fr(3)*rij(2)
       fbc(3,3) = -fr(3)*rij(3)
       pv_bc(:,:) = pv_bc(:,:) + fbc(:,:)
       IF (atprop_env%stress) THEN
          ! Atomic stress tensors
          iatom = shell_particle_set(i)%atom_index
          atprop_env%atstress(:,:,iatom) = atprop_env%atstress(:,:,iatom) + fbc(:,:)
       END IF
    END IF

  END SUBROUTINE bonded_correct_gaussian_low_sh

END MODULE fist_nonbond_force
