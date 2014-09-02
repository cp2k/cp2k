!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module that collects all Coulomb parts of the fock matrix construction
!>
!> \author Teodoro Laino (05.2009) [tlaino] - split and module reorganization
!> \par History
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : d-orbitals
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Speed-up
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Periodic SE
!>      Teodoro Laino (05.2009) [tlaino] - Stress Tensor
!>
! *****************************************************************************
MODULE scp_delta_fock_matrix
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             semi_empirical_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, &
       cp_dbcsr_deallocate_matrix, cp_dbcsr_deallocate_matrix_set, &
       cp_dbcsr_distribute, cp_dbcsr_get_block_diag, cp_dbcsr_get_block_p, &
       cp_dbcsr_init, cp_dbcsr_p_type, cp_dbcsr_replicate_all, cp_dbcsr_set, &
       cp_dbcsr_sum_replicated, cp_dbcsr_trace, cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                             ewald_pw_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_method_am1, do_method_mndo, do_method_mndod, do_method_pdg, &
       do_method_pm3, do_method_pm6, do_method_pnnl, do_method_rm1, &
       do_se_IS_slater
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type,&
                                             set_scp_env
  USE scp_fock_matrix_integrals,       ONLY: dscp2C,&
                                             dscp_1el,&
                                             scp2C,&
                                             scp2C_ew,&
                                             scp_1el,&
                                             scp_1el_ew
  USE scp_nddo_utils,                  ONLY: evaluate_scp_force
  USE semi_empirical_int_arrays,       ONLY: rij_threshold
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_type
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             se_int_control_type,&
                                             se_taper_type,&
                                             semi_empirical_p_type,&
                                             semi_empirical_type,&
                                             setup_se_int_control_type
  USE semi_empirical_utils,            ONLY: finalize_se_taper,&
                                             get_se_type,&
                                             initialize_se_taper
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_delta_fock_matrix'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module       = .FALSE.

  PUBLIC :: build_delta_fock_matrix

CONTAINS

! *****************************************************************************
!> \brief Construction of the Coulomb part of the Fock matrix
!> \param qs_env ...
!> \param ks_matrix ...
!> \param matrix_p ...
!> \param energy ...
!> \param calculate_forces ...
!> \param just_energy ...
!> \param store_int_env ...
!> \param error ...
!> \author JGH
! *****************************************************************************
  SUBROUTINE build_delta_fock_matrix (qs_env,ks_matrix,matrix_p,energy,calculate_forces,&
             just_energy,store_int_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_p
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_delta_fock_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iatom, ikind, inode, ispin, itype, &
      jatom, jkind, natom, natorb_a, natorb_a2, natorb_b, natorb_b2, nkind, &
      nspins, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: se_natorb
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: anag, check, defined, &
                                                failure, found, switch, &
                                                use_virial
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: se_defined
    REAL(KIND=dp)                            :: delta, dr1, ecore2, ecores, &
                                                ecoul_scp
    REAL(KIND=dp), DIMENSION(2025)           :: pa_a, pb_a, pscpa, pscpb
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    REAL(KIND=dp), DIMENSION(45, 45)         :: p_block_tot_a, p_block_tot_b, &
                                                pscp_block_tot_a, &
                                                pscp_block_tot_b
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ascp_diag
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: delta_ksa_block, &
      delta_ksa_tmp, delta_ksb_block, delta_ksb_tmp, k_scpa_block, &
      k_scpb_block, ksa_block_a, ksb_block_a, ktmp_a_block, ktmp_b_block, &
      pa_block_a, pb_block_a, pscpa_block, pscpb_block, ptmp_a_block, &
      ptmp_b_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: diagmat_ks, diagmat_p
    TYPE(cp_dbcsr_type), POINTER             :: delta_ks, delta_ks_tmp, &
                                                kmat_loc, ks_scp, pmat_loc, &
                                                pscp
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_se
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(scp_energy_type), POINTER           :: energy_scp
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(se_int_control_type)                :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_p_type), &
      DIMENSION(:), POINTER                  :: se_kind_list
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(virial_type), POINTER               :: virial

    failure=.FALSE.
    CALL timeset(routineN,handle)

    IF ( .NOT. failure ) THEN
       NULLIFY(dft_control,cell,force,particle_set,diagmat_ks,diagmat_p, pmat_loc, kmat_loc, &
               delta_ks, delta_ks_tmp, se_control,se_taper,pscp,ks_scp,scp_env,&
               energy_scp, ascp_diag)

       CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, se_taper=se_taper,&
            para_env=para_env, sab_se=sab_se, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
            particle_set=particle_set, virial=virial, error=error)

       ! Parameters
       CALL initialize_se_taper(se_taper,coulomb=.TRUE.,error=error)
       se_control  => dft_control%qs_control%se_control
       anag        =  se_control%analytical_gradients
       use_virial  = virial%pv_availability.AND.(.NOT.virial%pv_numer)

       ! Always do SLATER integral screening.
       ! Will need a smooth transition to the TT and ewald_gks.
       CALL setup_se_int_control_type(se_int_control, do_ewald_r3=se_control%do_ewald_r3,&
            do_ewald_gks=se_control%do_ewald_gks, integral_screening=do_se_IS_slater,&
            shortrange=(se_control%do_ewald.OR.se_control%do_ewald_gks),&
            max_multipole=se_control%max_multipole,pc_coulomb_int=.FALSE.)
       IF(se_control%do_ewald_gks) THEN
         CALL get_qs_env(qs_env=qs_env,ewald_env=ewald_env,ewald_pw=ewald_pw, error=error)
         CALL ewald_env_get (ewald_env, alpha=se_int_control%ewald_gks%alpha, error=error)
         CALL ewald_pw_get (ewald_pw, pw_big_pool=se_int_control%ewald_gks%pw_pool, &
                            dg=se_int_control%ewald_gks%dg)
       END IF

       nspins=dft_control%nspins
       CPPrecondition(ASSOCIATED(matrix_p),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)

       nkind = SIZE(atomic_kind_set)
       IF(calculate_forces) THEN
          CALL get_qs_env(qs_env=qs_env, force=force, error=error)
          natom = SIZE (particle_set)
          delta = se_control%delta
          ! Allocate atom index for kind
          ALLOCATE (atom_of_kind(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)
       END IF

       CALL cp_dbcsr_allocate_matrix_set(diagmat_ks,nspins,error=error)
       CALL cp_dbcsr_allocate_matrix_set(diagmat_p,nspins,error=error)
       ALLOCATE ( delta_ks )
       ALLOCATE ( delta_ks_tmp )
       ALLOCATE ( pmat_loc )
       ALLOCATE ( kmat_loc )

       DO ispin=1,nspins
          ! Allocate diagonal block matrices
          ALLOCATE(diagmat_p(ispin)%matrix,diagmat_ks(ispin)%matrix)!sm->dbcsr
          CALL cp_dbcsr_init(diagmat_p(ispin)%matrix, error=error)
          CALL cp_dbcsr_init(diagmat_ks(ispin)%matrix, error=error)
          CALL cp_dbcsr_get_block_diag(matrix_p(ispin)%matrix, diagmat_p(ispin)%matrix, error=error)
          CALL cp_dbcsr_get_block_diag(ks_matrix(ispin)%matrix, diagmat_ks(ispin)%matrix, error=error)
          CALL cp_dbcsr_set(diagmat_ks(ispin)%matrix, 0.0_dp, error=error)
          CALL cp_dbcsr_replicate_all(diagmat_p(ispin)%matrix, error=error)
          CALL cp_dbcsr_replicate_all(diagmat_ks(ispin)%matrix, error=error)
       END DO

       CALL cp_dbcsr_init(delta_ks_tmp, error=error)
       CALL cp_dbcsr_get_block_diag(ks_matrix(1)%matrix, delta_ks_tmp, error=error)
       CALL cp_dbcsr_set(delta_ks_tmp, 0.0_dp, error=error)
       CALL cp_dbcsr_replicate_all(delta_ks_tmp, error=error)

       CALL get_qs_env(qs_env=qs_env, scp_env=scp_env, error=error)
       CALL get_scp_env(scp_env=scp_env, pscp=pscp, ks_scp=ks_scp, energy=energy_scp, ascp_diag=ascp_diag, error=error)
       CALL cp_dbcsr_set ( ks_scp, 0.0_dp, error=error )

       CALL cp_dbcsr_init(delta_ks, error=error)
       CALL cp_dbcsr_init(pmat_loc, error=error)
       CALL cp_dbcsr_init(kmat_loc, error=error)

       CALL cp_dbcsr_copy(kmat_loc, ks_scp, error=error)
       CALL cp_dbcsr_copy(pmat_loc, ks_scp, error=error)
       CALL cp_dbcsr_copy(delta_ks, ks_scp, error=error)

       CALL cp_dbcsr_replicate_all(delta_ks, error=error)
       CALL cp_dbcsr_replicate_all(pmat_loc, error=error)
       CALL cp_dbcsr_replicate_all(kmat_loc, error=error)

       ecore2   = 0.0_dp
       itype    = get_se_type(dft_control%qs_control%method_id)
       ALLOCATE (se_defined(nkind),se_kind_list(nkind),se_natorb(nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO ikind=1,nkind
          CALL get_qs_kind(qs_kind_set(ikind), scp_se_parameter=se_kind_a)
          se_kind_list(ikind)%se_param => se_kind_a
          CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a)
          se_defined(ikind) = (defined .AND. natorb_a >= 1)
          se_natorb(ikind) = natorb_a
       END DO

       CALL neighbor_list_iterator_create(nl_iterator,sab_se)
       DO WHILE (neighbor_list_iterate(nl_iterator)==0)
          CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,inode=inode,r=rij)
          IF (.NOT.se_defined(ikind)) CYCLE
          IF (.NOT.se_defined(jkind)) CYCLE
          se_kind_a => se_kind_list(ikind)%se_param
          se_kind_b => se_kind_list(jkind)%se_param
          natorb_a = se_natorb(ikind)
          natorb_b = se_natorb(jkind)
          natorb_a2 = natorb_a**2
          natorb_b2 = natorb_b**2

          IF (inode==1) THEN
         ! density matrix
             CALL cp_dbcsr_get_block_p(matrix=diagmat_p(1)%matrix,&
                  row=iatom,col=iatom,BLOCK=pa_block_a,found=found)
             CPPostcondition(ASSOCIATED(pa_block_a),cp_failure_level,routineP,error,failure)
             CALL cp_dbcsr_get_block_p(matrix=pmat_loc,&
                  row=iatom,col=iatom,BLOCK=ptmp_a_block,found=found)
             CPPostcondition(ASSOCIATED(ptmp_a_block),cp_failure_level,routineP,error,failure)
          ! copy everything to structures of the correct size
             IF ( natorb_a == SIZE (pa_block_a,1 ) ) THEN
                ptmp_a_block = pa_block_a
             ELSE
                ptmp_a_block = 0.0_dp
                ptmp_a_block ( 1, 1 ) = pa_block_a ( 1, 1 )
             END IF
         ! KS matrix
             CALL cp_dbcsr_get_block_p(matrix=diagmat_ks(1)%matrix,&
                  row=iatom,col=iatom,BLOCK=ksa_block_a,found=found)
             CPPostcondition(ASSOCIATED(ksa_block_a),cp_failure_level,routineP,error,failure)
             CALL cp_dbcsr_get_block_p(matrix=kmat_loc,&
                  row=iatom,col=iatom,BLOCK=ktmp_a_block,found=found)
             CPPostcondition(ASSOCIATED(ktmp_a_block),cp_failure_level,routineP,error,failure)
          ! copy everything to structures of the correct size
             IF ( natorb_a == SIZE (ksa_block_a,1 ) ) THEN
                ktmp_a_block = ksa_block_a
             ELSE
                ktmp_a_block = 0.0_dp
                ktmp_a_block ( 1, 1 ) = ksa_block_a ( 1, 1 )
             END IF

             p_block_tot_a(1:natorb_a,1:natorb_a) = 2.0_dp * ptmp_a_block
             pa_a(1:natorb_a2)                    = RESHAPE(ptmp_a_block,(/natorb_a2/))

           ! pscp
             CALL cp_dbcsr_get_block_p(matrix=pscp,&
                  row=iatom,col=iatom,BLOCK=pscpa_block,found=found)
             CPPostcondition(ASSOCIATED(pscpa_block),cp_failure_level,routineP,error,failure)
             check = (SIZE(pscpa_block,1)==natorb_a).AND.(SIZE(pscpa_block,2)==natorb_a)
             CPPostcondition(check,cp_failure_level,routineP,error,failure)
             pscp_block_tot_a(1:natorb_a,1:natorb_a) = 2.0_dp*pscpa_block
             pscpa(1:natorb_a2)                    = RESHAPE(pscpa_block,(/natorb_a2/))

           ! kscp
             CALL cp_dbcsr_get_block_p(matrix=ks_scp,&
                  row=iatom,col=iatom,BLOCK=k_scpa_block,found=found)
             CPPostcondition(ASSOCIATED(k_scpa_block),cp_failure_level,routineP,error,failure)
             check = (SIZE(k_scpa_block,1)==natorb_a).AND.(SIZE(k_scpa_block,2)==natorb_a)
             CPPostcondition(check,cp_failure_level,routineP,error,failure)

           ! delta_ks
             CALL cp_dbcsr_get_block_p(matrix=delta_ks,&
                  row=iatom,col=iatom,BLOCK=delta_ksa_block,found=found)
             CPPostcondition(ASSOCIATED(delta_ksa_block),cp_failure_level,routineP,error,failure)
             check = (SIZE(delta_ksa_block,1)==natorb_a).AND.(SIZE(delta_ksa_block,2)==natorb_a)
             CPPostcondition(check,cp_failure_level,routineP,error,failure)

           ! delta_ks_tmp
             CALL cp_dbcsr_get_block_p(matrix=delta_ks_tmp,&
                  row=iatom,col=iatom,BLOCK=delta_ksa_tmp,found=found)
             CPPostcondition(ASSOCIATED(delta_ksa_tmp),cp_failure_level,routineP,error,failure)

          ENDIF

          dr1 = DOT_PRODUCT(rij,rij)
          IF ( dr1 > rij_threshold ) THEN
             ! Determine the order of the atoms, and in case switch them..
             IF (iatom <= jatom) THEN
                switch = .FALSE.
             ELSE
                switch = .TRUE.
             END IF
             ! density matrix
             CALL cp_dbcsr_get_block_p(matrix=diagmat_p(1)%matrix,&
                  row=jatom,col=jatom,BLOCK=pb_block_a,found=found)
             CPPostcondition(ASSOCIATED(pb_block_a),cp_failure_level,routineP,error,failure)
             CALL cp_dbcsr_get_block_p(matrix=pmat_loc,&
                  row=jatom,col=jatom,BLOCK=ptmp_b_block,found=found)
             CPPostcondition(ASSOCIATED(ptmp_b_block),cp_failure_level,routineP,error,failure)
          ! copy everything to structures of the correct size
             IF ( natorb_b == SIZE (pb_block_a,1 ) ) THEN
                ptmp_b_block = pb_block_a
             ELSE
                ptmp_b_block = 0.0_dp
                ptmp_b_block ( 1, 1 ) = pb_block_a ( 1, 1 )
             END IF
             ! KS matrix
             CALL cp_dbcsr_get_block_p(matrix=diagmat_ks(1)%matrix,&
                  row=jatom,col=jatom,BLOCK=ksb_block_a,found=found)
             CPPostcondition(ASSOCIATED(ksb_block_a),cp_failure_level,routineP,error,failure)
             CALL cp_dbcsr_get_block_p(matrix=kmat_loc,&
                  row=jatom,col=jatom,BLOCK=ktmp_b_block,found=found)
             CPPostcondition(ASSOCIATED(ktmp_b_block),cp_failure_level,routineP,error,failure)

          ! copy everything to structures of the correct size
             IF ( natorb_b == SIZE (ksb_block_a,1 ) ) THEN
                ktmp_b_block = ksb_block_a
             ELSE
                ktmp_b_block = 0.0_dp
                ktmp_b_block ( 1, 1 ) = ksb_block_a ( 1, 1 )
             END IF

             p_block_tot_b(1:natorb_b,1:natorb_b) = 2.0_dp * ptmp_b_block
             pb_a(1:natorb_b2)                    = RESHAPE(ptmp_b_block,(/natorb_b2/))

           ! pscp
             CALL cp_dbcsr_get_block_p(matrix=pscp,&
                  row=jatom,col=jatom,BLOCK=pscpb_block,found=found)
             CPPostcondition(ASSOCIATED(pscpb_block),cp_failure_level,routineP,error,failure)
             check = (SIZE(pscpb_block,1)==natorb_b).AND.(SIZE(pscpb_block,2)==natorb_b)
             CPPostcondition(check,cp_failure_level,routineP,error,failure)
             pscp_block_tot_b(1:natorb_b,1:natorb_b) = 2.0_dp*pscpb_block
             pscpb(1:natorb_b2)                    = RESHAPE(pscpb_block,(/natorb_b2/))

           ! ks_scp
             CALL cp_dbcsr_get_block_p(matrix=ks_scp,&
                  row=jatom,col=jatom,BLOCK=k_scpb_block,found=found)
             CPPostcondition(ASSOCIATED(k_scpb_block),cp_failure_level,routineP,error,failure)
             check = (SIZE(k_scpb_block,1)==natorb_b).AND.(SIZE(k_scpb_block,2)==natorb_b)
             CPPostcondition(check,cp_failure_level,routineP,error,failure)

           ! delta_ks
             CALL cp_dbcsr_get_block_p(matrix=delta_ks,&
                  row=jatom,col=jatom,BLOCK=delta_ksb_block,found=found)
             CPPostcondition(ASSOCIATED(delta_ksb_block),cp_failure_level,routineP,error,failure)
             check = (SIZE(delta_ksb_block,1)==natorb_b).AND.(SIZE(delta_ksb_block,2)==natorb_b)
             CPPostcondition(check,cp_failure_level,routineP,error,failure)

           ! delta_ks_tmp
             CALL cp_dbcsr_get_block_p(matrix=delta_ks_tmp,&
                  row=jatom,col=jatom,BLOCK=delta_ksb_tmp,found=found)
             CPPostcondition(ASSOCIATED(delta_ksb_tmp),cp_failure_level,routineP,error,failure)

             SELECT CASE (dft_control%qs_control%method_id)
             CASE (do_method_mndo, do_method_am1, do_method_pm3, do_method_pm6, do_method_pdg,&
                   do_method_rm1, do_method_mndod, do_method_pnnl)
                ! Two-centers One-electron terms  Pscp C Z
                CALL scp_1el (se_kind_a,se_kind_b,rij, k_scpa_block, k_scpb_block,&
                     pscpa, pscpb, ecore=ecore2, itype=itype, anag=anag, se_int_control=se_int_control,&
                     se_taper=se_taper, store_int_env=store_int_env, error=error)
                ! Coulomb Terms
                IF      ( nspins == 1 ) THEN
                   ! compute C Pscp, where delta_ks = C Pscp
                   CALL scp2C (se_kind_a, se_kind_b, rij, switch, pscp_block_tot_a, delta_ksa_block, pscp_block_tot_b,&
                        delta_ksb_block, factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                        store_int_env=store_int_env, error=error)
                   ! compute  C ( P + Z )
                   CALL scp2C (se_kind_a, se_kind_b, rij, switch, p_block_tot_a, k_scpa_block, p_block_tot_b,&
                        k_scpb_block, factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                        store_int_env=store_int_env, error=error)
                   ! copy everything to structures of the correct size
                   ! b block
                   IF ( natorb_b == SIZE (delta_ksb_tmp,1 ) ) THEN
                         delta_ksb_tmp =  delta_ksb_block
                   ELSE
                         delta_ksb_tmp ( 1, 1 ) =  delta_ksb_block ( 1, 1 )
                   END IF
                   ! a block
                   IF ( natorb_a == SIZE (delta_ksa_tmp,1 ) ) THEN
                         delta_ksa_tmp =  delta_ksa_block
                   ELSE
                         delta_ksa_tmp ( 1, 1 ) =  delta_ksa_block ( 1, 1 )
                   END IF
                END IF


                IF(calculate_forces) THEN
                   atom_a   = atom_of_kind(iatom)
                   atom_b   = atom_of_kind(jatom)

                   ! Derivatives of the Two-centre One-electron terms
                   force_ab = 0.0_dp

                   ! Derivatives of the Two-centre One-electron terms
                   CALL dscp_1el (se_kind_a,se_kind_b,rij, pscpa, pscpb, itype=itype, anag=anag,&
                        se_int_control=se_int_control, se_taper=se_taper, force=force_ab,&
                        delta=delta, error=error)


                   IF (use_virial) THEN
                      CALL virial_pair_force (virial%pv_virial, -1.0_dp, force_ab, rij, error)
                   END IF

                   ! Sum up force components
                   force(ikind)%all_potential(1,atom_a) = force(ikind)%all_potential(1,atom_a) - force_ab(1)
                   force(jkind)%all_potential(1,atom_b) = force(jkind)%all_potential(1,atom_b) + force_ab(1)

                   force(ikind)%all_potential(2,atom_a) = force(ikind)%all_potential(2,atom_a) - force_ab(2)
                   force(jkind)%all_potential(2,atom_b) = force(jkind)%all_potential(2,atom_b) + force_ab(2)

                   force(ikind)%all_potential(3,atom_a) = force(ikind)%all_potential(3,atom_a) - force_ab(3)
                   force(jkind)%all_potential(3,atom_b) = force(jkind)%all_potential(3,atom_b) + force_ab(3)

                   ! Derivatives of the Coulomb Terms
                   force_ab = 0._dp
                   IF      ( nspins == 1 ) THEN
                      CALL dscp2C(se_kind_a, se_kind_b, rij, switch, pscp_block_tot_a, pscp_block_tot_b, factor=0.25_dp,&
                           anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab, delta=delta,&
                           error=error)
                      CALL dscp2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, pscp_block_tot_b, factor=0.25_dp,&
                           anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab, delta=delta,&
                           error=error)
                      CALL dscp2C(se_kind_a, se_kind_b, rij, switch, pscp_block_tot_a, p_block_tot_b, factor=0.25_dp,&
                           anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab, delta=delta,&
                           error=error)
                   END IF
                   IF ( switch ) THEN
                      force_ab(1) = -force_ab(1)
                      force_ab(2) = -force_ab(2)
                      force_ab(3) = -force_ab(3)
                   END IF
                   IF (use_virial) THEN
                      CALL virial_pair_force ( virial%pv_virial, -1.0_dp, force_ab, rij, error)
                   END IF
                   ! Sum up force components
                   force(ikind)%rho_elec(1,atom_a) = force(ikind)%rho_elec(1,atom_a) - force_ab(1)
                   force(jkind)%rho_elec(1,atom_b) = force(jkind)%rho_elec(1,atom_b) + force_ab(1)

                   force(ikind)%rho_elec(2,atom_a) = force(ikind)%rho_elec(2,atom_a) - force_ab(2)
                   force(jkind)%rho_elec(2,atom_b) = force(jkind)%rho_elec(2,atom_b) + force_ab(2)

                   force(ikind)%rho_elec(3,atom_a) = force(ikind)%rho_elec(3,atom_a) - force_ab(3)
                   force(jkind)%rho_elec(3,atom_b) = force(jkind)%rho_elec(3,atom_b) + force_ab(3)
                END IF
             CASE DEFAULT
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
          ELSE
             IF ( se_int_control%do_ewald_gks ) THEN
                CPPostcondition(iatom==jatom,cp_failure_level,routineP,error,failure)
                ! Two-centers One-electron terms
                ecores=0._dp
                ! Pscp C Z
                CALL scp_1el_ew (se_kind_a,rij, k_scpa_block, pscpa, &
                     ecore=ecores, itype=itype, anag=anag, se_int_control=se_int_control,&
                     se_taper=se_taper, store_int_env=store_int_env, error=error)
                ecore2=ecore2+ecores
                ! Coulomb Terms
                IF      ( nspins == 1 ) THEN
                   ! compute C Pscp, where delta_ks = C Pscp
                   CALL scp2C_ew(se_kind_a, rij, pscp_block_tot_a, delta_ksa_block,&
                        factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                        store_int_env=store_int_env, error=error)
                   ! compute C ( P + Z )
                   CALL scp2C_ew(se_kind_a, rij, p_block_tot_a, k_scpa_block,&
                        factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                        store_int_env=store_int_env, error=error)
                   ! a block
                   IF ( natorb_a == SIZE (delta_ksa_tmp,1 ) ) THEN
                         delta_ksa_tmp =  delta_ksa_block
                   ELSE
                         delta_ksa_tmp ( 1, 1 ) =  delta_ksa_block ( 1, 1 )
                   END IF
                END IF
             END IF
          END IF
       END DO
       CALL neighbor_list_iterator_release(nl_iterator)

       DEALLOCATE(se_kind_list,se_defined,se_natorb,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ! sum up ks_scp, delta_ks, diagmat_ks on all procs
       CALL cp_dbcsr_sum_replicated(delta_ks, error=error)
       CALL cp_dbcsr_sum_replicated(delta_ks_tmp, error=error)
       CALL cp_dbcsr_sum_replicated( ks_scp, error=error )
       DO ispin = 1, nspins
         CALL cp_dbcsr_sum_replicated(diagmat_ks(ispin)%matrix, error=error)
       END DO
       ! Update the fock matrix.  Here K =  Pscp C
       CALL cp_dbcsr_add ( diagmat_ks(1)%matrix, delta_ks_tmp, 1.0_dp, 1.0_dp, error )
       ! add Kscp = Kscp (=C*(P+Z) ) + delta_ks (=C*Pscp).
       CALL cp_dbcsr_add ( ks_scp, delta_ks, 1.0_dp, 1.0_dp, error )
       ! ecoul_scp = p C (P + p + Z)
       CALL cp_dbcsr_trace ( ks_scp, pscp, trace=ecoul_scp, local_sum=.TRUE., error=error )
       CALL mp_sum(ecore2,para_env%group)
       energy_scp%e_scp_hartree = ecoul_scp + ecore2
       energy_scp%e_scp_hartree = 0.5_dp * energy_scp%e_scp_hartree
       CALL set_scp_env ( scp_env=scp_env, ks_scp = ks_scp, error = error )
       ! put diagmat_ks back into ks_matrix
       DO ispin = 1, nspins
          CALL cp_dbcsr_distribute(diagmat_ks(ispin)%matrix, error=error)
          CALL cp_dbcsr_add(ks_matrix(ispin)%matrix, diagmat_ks(ispin)%matrix,&
               1.0_dp,1.0_dp,error=error)
       END DO
!       CALL cp_dbcsr_print ( pscp, matlab_format = .TRUE., error=error )
       CALL cp_dbcsr_deallocate_matrix_set ( diagmat_p, error=error )
       CALL cp_dbcsr_deallocate_matrix_set ( diagmat_ks, error=error )
       CALL cp_dbcsr_deallocate_matrix ( delta_ks, error=error )
       CALL cp_dbcsr_deallocate_matrix ( delta_ks_tmp, error=error )
       CALL cp_dbcsr_deallocate_matrix ( pmat_loc, error=error )
       CALL cp_dbcsr_deallocate_matrix ( kmat_loc, error=error )

      ! evalulate kernel and forces on the SCP coeffs
       CALL evaluate_scp_force ( atomic_kind_set, qs_kind_set, pscp, ks_scp, ascp_diag, &
                                 just_energy, energy_scp%e_scp_kernel, error )

       IF (calculate_forces) THEN
          DEALLOCATE(atom_of_kind,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

! Add to total QS  energy
       energy% total =  energy_scp % e_scp_hartree + energy_scp % e_scp_kernel + energy_scp % e_scp_disp

    END IF
    CALL finalize_se_taper(se_taper,error=error)

    CALL timestop(handle)
  END SUBROUTINE build_delta_fock_matrix

END MODULE scp_delta_fock_matrix

