!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module that collects all Coulomb parts of the fock matrix construction
!>        
!> \author Teodoro Laino (05.2009) [tlaino] - split and module reorganization
!> \par History
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : d-orbitals
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Speed-up
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Periodic SE
!>      
! *****************************************************************************
MODULE se_fock_matrix_coulomb
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             semi_empirical_control_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                             ewald_pw_type
  USE ewalds_multipole,                ONLY: ewald_multipole_evaluate
  USE f77_blas
  USE fist_neighbor_list_control,      ONLY: list_control
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type
  USE input_constants,                 ONLY: &
       do_ewald_ewald, do_method_am1, do_method_mndo, do_method_mndod, &
       do_method_pdg, do_method_pm3, do_method_pm6, do_method_pnnl, &
       do_method_rm1, do_multipole_charge, do_multipole_dipole, &
       do_multipole_none, do_multipole_quadrupole, do_se_IS_slater
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type,&
                                             set_scp_env
  USE se_fock_matrix_integrals,        ONLY: &
       dfock2C, dfock2C_r3, dfock2_1el, dfock2_1el_r3, fock2C, fock2C_ew, &
       fock2C_r3, fock2_1el, fock2_1el_ew, fock2_1el_r3, &
       se_coulomb_ij_interaction
  USE semi_empirical_int_arrays,       ONLY: rij_threshold,&
                                             se_orbital_pointer
  USE semi_empirical_integrals,        ONLY: corecore_el,&
                                             dcorecore_el
  USE semi_empirical_mpole_methods,    ONLY: quadrupole_sph_to_cart
  USE semi_empirical_mpole_types,      ONLY: nddo_mpole_type,&
                                             semi_empirical_mpole_type
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_type
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             se_int_control_type,&
                                             se_taper_type,&
                                             semi_empirical_type,&
                                             setup_se_int_control_type
  USE semi_empirical_utils,            ONLY: finalize_se_taper,&
                                             get_se_type,&
                                             initialize_se_taper
  USE sparse_matrix_types,             ONLY: &
       allocate_matrix_set, cp_sm_scale_and_add, create_dbm_from_sm, &
       deallocate_dbm_from_sm, deallocate_matrix, get_block_node, mpsum_dbm, &
       real_matrix_p_type, real_matrix_type, transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_fock_matrix_coulomb'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module       = .FALSE.
  
  PUBLIC :: build_fock_matrix_coulomb, build_fock_matrix_coulomb_lr,&
            build_fock_matrix_coul_lr_r3

CONTAINS
  
! *****************************************************************************
!> \brief Construction of the Coulomb part of the Fock matrix
!> \author JGH
! *****************************************************************************
  SUBROUTINE build_fock_matrix_coulomb (qs_env,ks_matrix,matrix_p,energy,calculate_forces,&
       store_int_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_p
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_fock_matrix_coulomb', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iab, iatom, ikind, ilist, inode, &
      ispin, itype, jatom, jkind, natom, natorb_a, natorb_a2, natorb_b, &
      natorb_b2, nkind, nlist, nnode, nspins, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: anag, check, defined, &
                                                failure, scp_nddo, switch, &
                                                use_virial
    REAL(KIND=dp)                            :: delta, dr1, ecore2, ecores
    REAL(KIND=dp), DIMENSION(2025)           :: pa_a, pa_b, pb_a, pb_b
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    REAL(KIND=dp), DIMENSION(45, 45)         :: p_block_tot_a, p_block_tot_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: ksa_block_a, ksa_block_b, &
      ksb_block_a, ksb_block_b, pa_block_a, pa_block_b, pb_block_a, &
      pb_block_b, pscpa_block, pscpb_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_2c
    TYPE(neighbor_list_type), POINTER        :: sab_2c_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_2c_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: diagmat_ks, diagmat_p
    TYPE(real_matrix_type), POINTER          :: ks_scp, pscp
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(se_int_control_type)                :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(virial_type), POINTER               :: virial

    failure=.FALSE.
    CALL timeset(routineN,handle)
    IF ( .NOT. failure ) THEN
       NULLIFY(dft_control,cell,force,particle_set,diagmat_ks,diagmat_p,&
               se_control,se_taper,pscp,ks_scp,scp_env)
       
       CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, se_taper=se_taper,&
            para_env=para_env, sab_2c=sab_2c, atomic_kind_set=atomic_kind_set, &
            particle_set=particle_set, virial=virial, error=error)
      
       ! Parameters
       CALL initialize_se_taper(se_taper,coulomb=.TRUE.,error=error)
       se_control  => dft_control%qs_control%se_control
       anag        =  se_control%analytical_gradients
       scp_nddo    =  se_control%scp
       use_virial  = virial%pv_availability.AND.(.NOT.virial%pv_numer)

       CALL setup_se_int_control_type(se_int_control, do_ewald_r3=se_control%do_ewald_r3,&
            do_ewald_gks=se_control%do_ewald_gks, integral_screening=se_control%integral_screening,&
            shortrange=(se_control%do_ewald.OR.se_control%do_ewald_gks),&
            max_multipole=se_control%max_multipole,pc_coulomb_int=.FALSE.)
       IF(se_control%do_ewald_gks) THEN
         CALL get_qs_env(qs_env=qs_env,ewald_env=ewald_env,ewald_pw=ewald_pw,error=error)
         CALL ewald_env_get (ewald_env, alpha=se_int_control%ewald_gks%alpha, error=error)
         CALL ewald_pw_get (ewald_pw, pw_big_pool=se_int_control%ewald_gks%pw_pool, &
                            dg=se_int_control%ewald_gks%dg)
       END IF

       nspins=dft_control%nspins
       CPPrecondition(ASSOCIATED(matrix_p),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)
       
       CALL allocate_matrix_set(diagmat_ks,nspins,error=error)
       CALL allocate_matrix_set(diagmat_p,nspins,error=error)
       
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

       DO ispin=1,nspins
          ! Allocate diagonal block matrices
          CALL create_dbm_from_sm(diagmat_p(ispin)%matrix,matrix_p(ispin)%matrix,&
                                  "DIAGONAL BLOCK P MATRIX",.TRUE.,para_env,error=error)
          CALL create_dbm_from_sm(diagmat_ks(ispin)%matrix,ks_matrix(ispin)%matrix,&
                                  "DIAGONAL BLOCK KS MATRIX",.FALSE.,error=error)
       END DO

       IF (scp_nddo)  THEN 
          CALL get_qs_env(qs_env=qs_env, scp_env=scp_env, error=error)
          CALL get_scp_env(scp_env=scp_env, pscp=pscp, ks_scp=ks_scp, error=error)
       END IF
       
       ecore2   = 0.0_dp
       itype    = get_se_type(dft_control%qs_control%method_id)
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a)
 
          IF (.NOT.defined .OR. natorb_a < 1) CYCLE
          natorb_a2 = natorb_a**2
          DO jkind=1,nkind
             atomic_kind => atomic_kind_set(jkind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_b)
             CALL get_se_param(se_kind_b,defined=defined,natorb=natorb_b)

             IF (.NOT.defined .OR. natorb_b < 1) CYCLE
             natorb_b2 = natorb_b**2
             iab = ikind + nkind*(jkind - 1)
             IF (.NOT.ASSOCIATED(sab_2c(iab)%neighbor_list_set)) CYCLE
             CALL get_neighbor_list_set(neighbor_list_set=sab_2c(iab)%neighbor_list_set,nlist=nlist)
             NULLIFY(sab_2c_neighbor_list)

             DO ilist=1,nlist
                IF ( .NOT. ASSOCIATED(sab_2c_neighbor_list) ) THEN
                   sab_2c_neighbor_list => first_list(sab_2c(iab)%neighbor_list_set)
                ELSE
                   sab_2c_neighbor_list => next(sab_2c_neighbor_list)
                END IF
                CALL get_neighbor_list(neighbor_list=sab_2c_neighbor_list,atom=iatom, nnode=nnode)
                sab_2c_neighbor_node => first_node(sab_2c_neighbor_list)
 
                CALL get_block_node(matrix=diagmat_p(1)%matrix,&
                                    block_row=iatom,&
                                    block_col=iatom,&
                                    BLOCK=pa_block_a)
                CPPostcondition(ASSOCIATED(pa_block_a),cp_failure_level,routineP,error,failure)
                check = (SIZE(pa_block_a,1)==natorb_a).AND.(SIZE(pa_block_a,2)==natorb_a)
                CPPostcondition(check,cp_failure_level,routineP,error,failure)
                CALL get_block_node(matrix=diagmat_ks(1)%matrix,&
                                    block_row=iatom,&
                                    block_col=iatom,&
                                    BLOCK=ksa_block_a)
                CPPostcondition(ASSOCIATED(ksa_block_a),cp_failure_level,routineP,error,failure)
                p_block_tot_a(1:natorb_a,1:natorb_a) = 2.0_dp * pa_block_a
                pa_a(1:natorb_a2)                    = RESHAPE(pa_block_a,(/natorb_a2/)) 
                IF ( nspins == 2 ) THEN
                   CALL get_block_node(matrix=diagmat_p(2)%matrix,&
                                       block_row=iatom,&
                                       block_col=iatom,&
                                       BLOCK=pa_block_b)
                   CPPostcondition(ASSOCIATED(pa_block_b),cp_failure_level,routineP,error,failure)
                   check = (SIZE(pa_block_b,1)==natorb_a).AND.(SIZE(pa_block_b,2)==natorb_a)
                   CPPostcondition(check,cp_failure_level,routineP,error,failure)
                   CALL get_block_node(matrix=diagmat_ks(2)%matrix,&
                                       block_row=iatom,&
                                       block_col=iatom,&
                                       BLOCK=ksa_block_b)
                   CPPostcondition(ASSOCIATED(ksa_block_b),cp_failure_level,routineP,error,failure)
                   p_block_tot_a(1:natorb_a,1:natorb_a) = pa_block_a + pa_block_b
                   pa_b(1:natorb_a2)                    = RESHAPE(pa_block_b,(/natorb_a2/))
                END IF

                IF (scp_nddo) THEN
                   CALL get_block_node(matrix=pscp,&
                                       block_row=iatom,&
                                       block_col=iatom,&
                                       BLOCK=pscpa_block)
                   CPPostcondition(ASSOCIATED(pscpa_block),cp_failure_level,routineP,error,failure)
                   check = (SIZE(pscpa_block,1)==natorb_a).AND.(SIZE(pscpa_block,2)==natorb_a)
                   CPPostcondition(check,cp_failure_level,routineP,error,failure)
                   ! Add P_scp + P_block_tot_a to get P_tot
                   p_block_tot_a(1:natorb_a,1:natorb_a) = p_block_tot_a(1:natorb_a,1:natorb_a) + 2.0_dp * pscpa_block
                   ! Modify accordingly the density matrix for SCP
                        pa_a(1:natorb_a2)       = RESHAPE(pa_block_a+pscpa_block,(/natorb_a2/))
                   IF ( nspins == 2 ) &
                        pa_b(1:natorb_a2)       = RESHAPE(pa_block_b+pscpa_block,(/natorb_a2/))
                ENDIF
 
                DO inode=1,nnode
                   CALL get_neighbor_node(neighbor_node=sab_2c_neighbor_node,neighbor=jatom,r=rij)
 
                   dr1 = DOT_PRODUCT(rij,rij)
                   IF ( dr1 > rij_threshold ) THEN
                      ! Determine the order of the atoms, and in case switch them..
                      IF (iatom <= jatom) THEN
                         switch = .FALSE.
                      ELSE
                         switch = .TRUE.
                      END IF
                      ! Retrieve blocks for KS and P
                      CALL get_block_node(matrix=diagmat_p(1)%matrix,&
                                          block_row=jatom,&
                                          block_col=jatom,&
                                          BLOCK=pb_block_a)
                      CPPostcondition(ASSOCIATED(pb_block_a),cp_failure_level,routineP,error,failure)
                      check = (SIZE(pb_block_a,1)==natorb_b).AND.(SIZE(pb_block_a,2)==natorb_b)
                      CPPostcondition(check,cp_failure_level,routineP,error,failure)
                      CALL get_block_node(matrix=diagmat_ks(1)%matrix,&
                                          block_row=jatom,&
                                          block_col=jatom,&
                                          BLOCK=ksb_block_a)
                      CPPostcondition(ASSOCIATED(ksb_block_a),cp_failure_level,routineP,error,failure)
                      p_block_tot_b(1:natorb_b,1:natorb_b) = 2.0_dp * pb_block_a
                      pb_a(1:natorb_b2)                    = RESHAPE(pb_block_a,(/natorb_b2/)) 
                      ! Handle more than one configuration
                      IF ( nspins == 2 ) THEN
                         CALL get_block_node(matrix=diagmat_p(2)%matrix,&
                                             block_row=jatom,&
                                             block_col=jatom,&
                                             BLOCK=pb_block_b)
                         CPPostcondition(ASSOCIATED(pb_block_b),cp_failure_level,routineP,error,failure)
                         check = (SIZE(pb_block_b,1)==natorb_b).AND.(SIZE(pb_block_b,2)==natorb_b)
                         CPPostcondition(check,cp_failure_level,routineP,error,failure)
                         CALL get_block_node(matrix=diagmat_ks(2)%matrix,&
                                             block_row=jatom,&
                                             block_col=jatom,&
                                             BLOCK=ksb_block_b)
                         CPPostcondition(ASSOCIATED(ksb_block_b),cp_failure_level,routineP,error,failure)
                         check = (SIZE(pb_block_a,1)==SIZE(pb_block_b,1)).AND.(SIZE(pb_block_a,2)==SIZE(pb_block_b,2))
                         CPPostcondition(check,cp_failure_level,routineP,error,failure)
                         p_block_tot_b(1:natorb_b,1:natorb_b) = pb_block_a + pb_block_b
                         pb_b(1:natorb_b2)                    = RESHAPE(pb_block_b,(/natorb_b2/))
                      END IF

                      IF (scp_nddo) THEN
                         CALL get_block_node(matrix=pscp,&
                                             block_row=jatom,&
                                             block_col=jatom,&
                                             BLOCK=pscpb_block)
                         CPPostcondition(ASSOCIATED(pscpb_block),cp_failure_level,routineP,error,failure)
                         check = (SIZE(pscpb_block,1)==natorb_b).AND.(SIZE(pscpb_block,2)==natorb_b)
                         CPPostcondition(check,cp_failure_level,routineP,error,failure)
                         ! Add P_scp + P_block_tot_b to get P_tot
                         p_block_tot_b(1:natorb_b,1:natorb_b) = p_block_tot_b(1:natorb_b,1:natorb_b) + 2.0_dp * pscpb_block
                         ! Modify accordingly the density matrix for SCP
                         pb_a(1:natorb_b2) = RESHAPE(pb_block_a+pscpb_block,(/natorb_b2/))
                         IF (nspins == 2) pb_b(1:natorb_b2) = RESHAPE(pb_block_b+pscpb_block,(/natorb_b2/))
                      ENDIF 
 
                      SELECT CASE (dft_control%qs_control%method_id)
                      CASE (do_method_mndo, do_method_am1, do_method_pm3, do_method_pm6, do_method_pdg,&
                            do_method_rm1, do_method_mndod, do_method_pnnl)

                         ! Two-centers One-electron terms
                         IF      ( nspins == 1 ) THEN
                            CALL fock2_1el (se_kind_a,se_kind_b,rij, ksa_block_a, ksb_block_a,&
                                 pa_a, pb_a, ecore=ecore2, itype=itype, anag=anag, se_int_control=se_int_control,&
                                 se_taper=se_taper, store_int_env=store_int_env, error=error)
                         ELSE IF ( nspins == 2 ) THEN
                            CALL fock2_1el (se_kind_a,se_kind_b,rij, ksa_block_a, ksb_block_a,&
                                 pa_block_a, pb_block_a, ecore=ecore2, itype=itype, anag=anag,&
                                 se_int_control=se_int_control, se_taper=se_taper, store_int_env=store_int_env,&
                                 error=error)
                            CALL fock2_1el (se_kind_a,se_kind_b,rij, ksa_block_b, ksb_block_b,&
                                 pa_b, pb_b, ecore=ecore2, itype=itype, anag=anag, se_int_control=se_int_control,&
                                 se_taper=se_taper, store_int_env=store_int_env, error=error)
                         END IF
                         ! Coulomb Terms
                         IF      ( nspins == 1 ) THEN
                            CALL fock2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, ksa_block_a, p_block_tot_b,&
                                 ksb_block_a, factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                                 store_int_env=store_int_env, error=error)
                         ELSE IF ( nspins == 2 ) THEN
                            CALL fock2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, ksa_block_a, p_block_tot_b,&
                                 ksb_block_a, factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                                 store_int_env=store_int_env, error=error)
 
                            CALL fock2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, ksa_block_b, p_block_tot_b,&
                                 ksb_block_b, factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                                 store_int_env=store_int_env, error=error)
                         END IF
 
                         IF(calculate_forces) THEN
                            atom_a   = atom_of_kind(iatom)
                            atom_b   = atom_of_kind(jatom)
 
                            ! Derivatives of the Two-centre One-electron terms
                            force_ab = 0.0_dp
                            IF      ( nspins == 1 ) THEN
                               CALL dfock2_1el (se_kind_a,se_kind_b,rij, pa_a, pb_a, itype=itype, anag=anag,&
                                    se_int_control=se_int_control, se_taper=se_taper, force=force_ab,&
                                    delta=delta, error=error)
                            ELSE IF ( nspins == 2 ) THEN
                               CALL dfock2_1el (se_kind_a,se_kind_b,rij, pa_block_a, pb_block_a, itype=itype, anag=anag,&
                                    se_int_control=se_int_control, se_taper=se_taper, force=force_ab, delta=delta,&
                                    error=error)
                               CALL dfock2_1el (se_kind_a,se_kind_b,rij, pa_b, pb_b, itype=itype, anag=anag,&
                                    se_int_control=se_int_control, se_taper=se_taper, force=force_ab, delta=delta,&
                                    error=error)
                            END IF

                            IF (use_virial) THEN
                               CALL virial_pair_force ( virial%pv_virial, -1.0_dp, force_ab, rij, error)
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
                               CALL dfock2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, p_block_tot_b, factor=0.25_dp,&
                                    anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab, delta=delta,&
                                    error=error)
                            ELSE IF ( nspins == 2 ) THEN
                               CALL dfock2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, p_block_tot_b, factor=0.50_dp,&
                                    anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab, delta=delta,&
                                    error=error)
                               
                               CALL dfock2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, p_block_tot_b, factor=0.50_dp,&
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
                         IF      ( nspins == 1 ) THEN
                            CALL fock2_1el_ew (se_kind_a,rij, ksa_block_a, pa_a, &
                                 ecore=ecores, itype=itype, anag=anag, se_int_control=se_int_control,&
                                 se_taper=se_taper, store_int_env=store_int_env, error=error)
                         ELSE IF ( nspins == 2 ) THEN
                            CALL fock2_1el_ew (se_kind_a,rij, ksa_block_a, pa_block_a,&
                                 ecore=ecores, itype=itype, anag=anag, se_int_control=se_int_control, &
                                 se_taper=se_taper, store_int_env=store_int_env, error=error)
                            CALL fock2_1el_ew (se_kind_a,rij, ksa_block_b, pa_b,&
                                 ecore=ecores, itype=itype, anag=anag, se_int_control=se_int_control,&
                                 se_taper=se_taper, store_int_env=store_int_env, error=error)
                         END IF
                         ecore2=ecore2+ecores
                         ! Coulomb Terms
                         IF      ( nspins == 1 ) THEN
                            CALL fock2C_ew(se_kind_a, rij, p_block_tot_a, ksa_block_a,& 
                                 factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                                 store_int_env=store_int_env, error=error)
                         ELSE IF ( nspins == 2 ) THEN
                            CALL fock2C_ew(se_kind_a, rij, p_block_tot_a, ksa_block_a,&
                                 factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                                 store_int_env=store_int_env, error=error)
                            CALL fock2C_ew(se_kind_a, rij, p_block_tot_a, ksa_block_b,&
                                 factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                                 store_int_env=store_int_env, error=error)
                         END IF
                      END IF
                   END IF
                   sab_2c_neighbor_node => next(sab_2c_neighbor_node)
                END DO ! inode => jatom(atom B)
             END DO ! ilist => iatom(atom A)
          END DO ! jkind
       END DO ! ikind

       IF ( scp_nddo ) THEN
          ! Update the SCP coulomb fock matrix.  Should be K = (P + p )C + CZ
          CALL transfer_matrix ( diagmat_ks(1)%matrix, ks_scp, error )
          ! sum over all procs to get the replicated ks_scp
          CALL mpsum_dbm ( ks_scp, para_env, error )
          CALL set_scp_env ( scp_env=scp_env, ks_scp = ks_scp, error = error ) 
       ENDIF 
       
       DO ispin = 1, nspins
          CALL deallocate_dbm_from_sm(diagmat_p(ispin)%matrix,matrix_p(ispin)%matrix,.FALSE.,error=error)
          CALL deallocate_dbm_from_sm(diagmat_ks(ispin)%matrix,ks_matrix(ispin)%matrix,.TRUE.,para_env,error=error)
       END DO
       DEALLOCATE(diagmat_ks,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(diagmat_p,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       
       IF (calculate_forces) THEN 
          DEALLOCATE(atom_of_kind,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       
       ! Two-centers one-electron terms
       CALL mp_sum(ecore2,para_env%group)       
       energy%hartree = ecore2 - energy%core
    END IF
    CALL finalize_se_taper(se_taper,error=error)    
    CALL timestop(handle)
  END SUBROUTINE build_fock_matrix_coulomb

! *****************************************************************************
!> \brief  Long-Range part for SE Coulomb interactions
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   08.2008 [created]
! *****************************************************************************
  SUBROUTINE build_fock_matrix_coulomb_lr (qs_env, ks_matrix, matrix_p, energy,&
       calculate_forces, store_int_env, error)
   
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_p
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_fock_matrix_coulomb_lr', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, ewald_type, forces_g_size, handle, iatom, ikind, &
      ilist, indi, indj, ispin, itype, iw, jint, natoms, natorb_a, nkind, &
      nlocal_particles, node, nparticle_local, nspins, size_1c_int, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: anag, defined, failure, &
                                                scp_nddo, use_virial
    LOGICAL, DIMENSION(3)                    :: task
    REAL(KIND=dp)                            :: e_neut, e_self, energy_glob, &
                                                energy_local, enuc, fac, tmp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: forces_g, forces_r
    REAL(KIND=dp), DIMENSION(3)              :: force_a
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_glob, pv_local, qcart
    REAL(KIND=dp), DIMENSION(5)              :: qsph
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ksa_block_a, pa_block_a, &
                                                pscpa_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(fist_nonbond_env_type), POINTER     :: se_nonbond_env
    TYPE(nddo_mpole_type), POINTER           :: se_nddo_mpole
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: diagmat_ks, diagmat_p
    TYPE(real_matrix_type), POINTER          :: dbm, ks_scp, pscp
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER         :: se_section
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_mpole_type), POINTER :: mpole
    TYPE(semi_empirical_type), POINTER       :: se_kind_a
    TYPE(virial_type), POINTER               :: virial

    failure=.FALSE.
    CALL timeset(routineN,handle)
    IF ( .NOT. failure ) THEN
       NULLIFY(dft_control, cell, force,particle_set, diagmat_ks, diagmat_p,local_particles,&
            se_control, ewald_env, ewald_pw, se_nddo_mpole, se_nonbond_env, se_section, mpole,&
            logger,pscp,ks_scp,scp_env,dbm)

       logger => cp_error_get_logger(error)
       CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, para_env=para_env,&
            atomic_kind_set=atomic_kind_set, particle_set=particle_set, ewald_env=ewald_env,&
            local_particles=local_particles, ewald_pw=ewald_pw, se_nddo_mpole=se_nddo_mpole,&
            se_nonbond_env=se_nonbond_env, scp_env=scp_env, virial=virial, error=error)
       
       nlocal_particles = SUM(local_particles%n_el(:))
       natoms           = SIZE(particle_set)
       CALL ewald_env_get (ewald_env, ewald_type=ewald_type, error=error)
       SELECT CASE(ewald_type)
       CASE (do_ewald_ewald)
          forces_g_size= nlocal_particles
       CASE DEFAULT
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Periodic SE implemented only for standard EWALD sums.", &
               error=error, error_level=cp_failure_level)
       END SELECT

       ! Parameters
       se_section => section_vals_get_subs_vals(qs_env%input,"DFT%QS%SE",error=error)
       se_control => dft_control%qs_control%se_control
       scp_nddo   =  se_control%scp
       anag       =  se_control%analytical_gradients
       use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
       
       nspins=dft_control%nspins
       CPPrecondition(ASSOCIATED(matrix_p),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)
       
       CALL allocate_matrix_set(diagmat_ks,nspins,error=error)
       CALL allocate_matrix_set(diagmat_p,nspins,error=error)
       
       nkind = SIZE(atomic_kind_set)
       DO ispin=1,nspins
          ! Allocate diagonal block matrices
          CALL create_dbm_from_sm(diagmat_p(ispin)%matrix,matrix_p(ispin)%matrix,&
                                  "DIAGONAL BLOCK P MATRIX",.TRUE.,para_env,error=error)
          CALL create_dbm_from_sm(diagmat_ks(ispin)%matrix,ks_matrix(ispin)%matrix,&
                                  "DIAGONAL BLOCK KS MATRIX",.FALSE.,error=error)
       END DO

       ! Handle the case of SCP-NDDO
       IF (scp_nddo) THEN
          CALL create_dbm_from_sm(dbm,ks_matrix(1)%matrix,&
                                  "DUMMY BLOCK KS MATRIX",.FALSE.,error=error)
          CALL get_scp_env(scp_env=scp_env, pscp=pscp, ks_scp=ks_scp, error=error)
       END IF

       ! Check for implemented SE methods
       SELECT CASE (dft_control%qs_control%method_id)
       CASE(do_method_mndo, do_method_am1, do_method_pm3, do_method_pm6, do_method_pdg,&
            do_method_rm1, do_method_mndod, do_method_pnnl)
          itype    = get_se_type(dft_control%qs_control%method_id)
       CASE DEFAULT
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT

       ! Zero arrays and possibly build neighbor lists
       energy_local = 0.0_dp
       energy_glob  = 0.0_dp
       e_neut       = 0.0_dp
       e_self       = 0.0_dp
       task         = .FALSE.
       SELECT CASE(se_control%max_multipole)
       CASE (do_multipole_none)
          ! Do Nothing
       CASE(do_multipole_charge)
          task(1) = .TRUE.
       CASE(do_multipole_dipole)
          task    = .TRUE.
          task(3) = .FALSE.
       CASE(do_multipole_quadrupole)
          task    = .TRUE.
       CASE DEFAULT
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       
       ! Build-up neighbor lists for real-space part of Ewald multipoles
       CALL list_control ( atomic_kind_set, particle_set, local_particles, &
            cell, se_nonbond_env, para_env, se_section, error=error)

       enuc = 0.0_dp
       se_nddo_mpole%charge     = 0.0_dp
       se_nddo_mpole%dipole     = 0.0_dp
       se_nddo_mpole%quadrupole = 0.0_dp

       DO ispin = 1, nspins
          ! Compute the NDDO mpole expansion
          DO ikind=1,nkind
             atomic_kind => atomic_kind_set(ikind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_a)
             CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a)

             IF (.NOT.defined .OR. natorb_a < 1) CYCLE
             
             nparticle_local = local_particles%n_el(ikind)
             DO ilist=1, nparticle_local
                iatom = local_particles%list(ikind)%array(ilist)
                CALL get_block_node(matrix=diagmat_p(ispin)%matrix,&
                                    block_row=iatom,&
                                    block_col=iatom,&
                                    BLOCK=pa_block_a)
                CPPostcondition(ASSOCIATED(pa_block_a),cp_failure_level,routineP,error,failure)
                IF ( scp_nddo ) THEN
                   CALL get_block_node(matrix=pscp,&
                                       block_row=iatom,&
                                       block_col=iatom,&
                                       BLOCK=pscpa_block)
                   CPPostcondition(ASSOCIATED(pscpa_block),cp_failure_level,routineP,error,failure)
                   pa_block_a=pa_block_a + pscpa_block 
                END IF

                ! Nuclei
                IF (task(1).AND.ispin==1) se_nddo_mpole%charge(iatom)   = se_kind_a%zeff
                ! Electrons
                size_1c_int = SIZE(se_kind_a%w_mpole)
                DO jint = 1, size_1c_int
                   mpole => se_kind_a%w_mpole(jint)%mpole
                   indi  = se_orbital_pointer(mpole%indi)
                   indj  = se_orbital_pointer(mpole%indj)
                   fac   = 1.0_dp
                   IF (indi/=indj) fac = 2.0_dp

                   ! Charge
                   IF (mpole%task(1).AND.task(1)) THEN
                      se_nddo_mpole%charge(iatom)          = se_nddo_mpole%charge(iatom)   +&
                                                             fac*pa_block_a(indi,indj)*mpole%c
                   END IF

                   ! Dipole
                   IF (mpole%task(2).AND.task(2)) THEN
                      se_nddo_mpole%dipole(:,iatom)        = se_nddo_mpole%dipole(:,iatom) +&
                                                             fac*pa_block_a(indi,indj)*mpole%d(:)
                   END IF
                   
                   ! Quadrupole
                   IF (mpole%task(3).AND.task(3)) THEN
                      qsph = fac*mpole%qs * pa_block_a(indi,indj)
                      CALL quadrupole_sph_to_cart(qcart, qsph, error)
                      se_nddo_mpole%quadrupole(:,:,iatom)  = se_nddo_mpole%quadrupole(:,:,iatom) +&
                                                             qcart
                   END IF
                END DO
                ! Print some info about charge, dipole and quadrupole (debug purpose only)
                IF (debug_this_module) THEN
                   WRITE(*,'(I5,F12.6,5X,3F12.6,5X,9F12.6)')iatom, se_nddo_mpole%charge(iatom),&
                        se_nddo_mpole%dipole(:,iatom),  se_nddo_mpole%quadrupole(:,:,iatom)
                END IF
             END DO
          END DO
       END DO
       CALL mp_sum(se_nddo_mpole%charge, para_env%group)
       CALL mp_sum(se_nddo_mpole%dipole, para_env%group)
       CALL mp_sum(se_nddo_mpole%quadrupole, para_env%group)

       ! Ewald Multipoles Sum
       iw = cp_print_key_unit_nr(logger,se_section,"PRINT%EWALD_INFO",extension=".seLog",error=error)
       IF(calculate_forces) THEN
          CALL get_qs_env(qs_env=qs_env, force=force, error=error)
             
          ! Allocate atom index for kind
          ALLOCATE (atom_of_kind(natoms),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind)
             
          ! Allocate and zeroing arrays
          ALLOCATE ( forces_g(3, forces_g_size), STAT=stat )
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE ( forces_r(3, natoms), STAT=stat )
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          forces_g = 0.0_dp
          forces_r = 0.0_dp
          pv_glob  = 0.0_dp
          pv_local = 0.0_dp
             
          CALL ewald_multipole_evaluate(ewald_env, ewald_pw, se_nonbond_env, cell,&
                  particle_set, local_particles, energy_local, energy_glob, e_neut, e_self, task,&
                  do_forces=.TRUE., do_stress=use_virial, do_efield=.TRUE., charges=se_nddo_mpole%charge,&
                  dipoles=se_nddo_mpole%dipole, quadrupoles=se_nddo_mpole%quadrupole,&
                  forces_local=forces_g, forces_glob=forces_r, pv_glob=pv_glob, pv_local=pv_local,&
                  efield0=se_nddo_mpole%efield0,&
                  efield1=se_nddo_mpole%efield1, efield2=se_nddo_mpole%efield2, iw=iw, &
                  do_debug=.TRUE.,error=error)
          ! Only SR force have to be summed up.. the one in g-space are already fully local..
          CALL mp_sum(forces_r, para_env%group)
          IF ( use_virial ) THEN
             ! why do we not need pv_local?
             virial%pv_virial = virial%pv_virial + pv_glob 
          END IF
       ELSE
          CALL ewald_multipole_evaluate(ewald_env, ewald_pw, se_nonbond_env, cell,&
                  particle_set, local_particles, energy_local, energy_glob, e_neut, e_self, task,&
                  do_forces=.FALSE., do_stress=.FALSE., do_efield=.TRUE., charges=se_nddo_mpole%charge,&
                  dipoles=se_nddo_mpole%dipole, quadrupoles=se_nddo_mpole%quadrupole,&
                  efield0=se_nddo_mpole%efield0, efield1=se_nddo_mpole%efield1,&
                  efield2=se_nddo_mpole%efield2, iw=iw, do_debug=.TRUE.,error=error)
       END IF
       CALL cp_print_key_finished_output(iw,logger,se_section,"PRINT%EWALD_INFO",error=error)

       ! Apply correction only when the Integral Scheme is different from Slater
       IF (se_control%integral_screening/=do_se_IS_slater) THEN
          CALL build_fock_matrix_coul_lrc(qs_env, ks_matrix, matrix_p, energy, calculate_forces,&
               store_int_env, se_nddo_mpole, task, error)
       END IF

       ! Debug Statements
       IF (debug_this_module) THEN
          CALL mp_sum(energy_glob,para_env%group)
          WRITE(*,*)"TOTAL ENERGY AFTER EWALD:",energy_local+ energy_glob+ e_neut+ e_self,&
                  energy_local, energy_glob, e_neut, e_self
       END IF

       ! Modify the KS matrix and possibly compute derivatives
       node = 0 
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a)
             
          IF (.NOT.defined .OR. natorb_a < 1) CYCLE
             
          nparticle_local = local_particles%n_el(ikind)
          DO ilist=1, nparticle_local
             node  = node + 1
             iatom = local_particles%list(ikind)%array(ilist)
             DO ispin = 1, nspins
                CALL get_block_node(matrix=diagmat_ks(ispin)%matrix,&
                                    block_row=iatom,&
                                    block_col=iatom,&
                                    BLOCK=ksa_block_a)
                CPPostcondition(ASSOCIATED(ksa_block_a),cp_failure_level,routineP,error,failure)

                ! Modify Hamiltonian Matrix accordingly potential, field and electric field gradient
                size_1c_int = SIZE(se_kind_a%w_mpole)
                DO jint = 1, size_1c_int
                   tmp   =  0.0_dp
                   mpole => se_kind_a%w_mpole(jint)%mpole
                   indi  =  se_orbital_pointer(mpole%indi)
                   indj  =  se_orbital_pointer(mpole%indj)

                   ! Charge
                   IF (mpole%task(1).AND.task(1)) THEN
                      tmp   = tmp + mpole%c*se_nddo_mpole%efield0(iatom)
                   END IF
                   
                   ! Dipole
                   IF (mpole%task(2).AND.task(2)) THEN
                      tmp   = tmp - DOT_PRODUCT(mpole%d,se_nddo_mpole%efield1(:,iatom))
                   END IF
                   
                   ! Quadrupole
                   IF (mpole%task(3).AND.task(3)) THEN
                      tmp   = tmp - (1.0_dp/3.0_dp) * SUM(mpole%qc * RESHAPE(se_nddo_mpole%efield2(:,iatom),(/3,3/)))
                   END IF
                   
                   ksa_block_a(indi, indj) = ksa_block_a(indi, indj) + tmp
                   ksa_block_a(indj, indi) = ksa_block_a(indi, indj)
                END DO
             END DO

             ! Nuclear term and forces
             IF (task(1)) enuc = enuc + se_kind_a%zeff * se_nddo_mpole%efield0(iatom)
             IF(calculate_forces) THEN
                atom_a   = atom_of_kind(iatom)
                force_a  = forces_r(1:3,iatom) + forces_g(1:3,node)
                ! Derivatives of the periodic Coulomb Terms
                force(ikind)%all_potential(:,atom_a) = force(ikind)%all_potential(:,atom_a) - force_a(:)
             END IF                
          END DO
       END DO
       ! Sum nuclear energy contribution
       CALL mp_sum(enuc, para_env%group)
       energy%core_overlap = energy%core_overlap0 + 0.5_dp*enuc

       ! Debug Statements
       IF (debug_this_module) THEN
          WRITE(*,*)"ENUC: ",enuc*0.5_dp
       END IF

       ! Increment the SCP portion of the Fock matrix with the long-range coulomb portion of the Fock matrix
       IF ( scp_nddo ) THEN
          ! transfer to dummy dbm
          CALL transfer_matrix ( diagmat_ks ( 1 ) % matrix, dbm, error )
          ! replicate the dbm
          CALL mpsum_dbm ( dbm, para_env, error )
          ! add dbm (which is the long range part of (P+p)C + CZ )to the replicated ks_scp
          CALL cp_sm_scale_and_add ( ks_scp, 1.0_dp, dbm, 1.0_dp, error )  
          CALL set_scp_env ( scp_env=scp_env, ks_scp = ks_scp, error = error ) 
          CALL deallocate_matrix ( dbm, error=error )
       ENDIF
       
       ! Deallocate diagonal block matrix
       DO ispin = 1, nspins
          CALL deallocate_dbm_from_sm(diagmat_p(ispin)%matrix,matrix_p(ispin)%matrix,.FALSE.,error=error)
          CALL deallocate_dbm_from_sm(diagmat_ks(ispin)%matrix,ks_matrix(ispin)%matrix,.TRUE.,para_env,error=error)
       END DO
       DEALLOCATE(diagmat_ks,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(diagmat_p,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Set the Fock matrix contribution to SCP
       IF (calculate_forces) THEN 
          DEALLOCATE(atom_of_kind,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(forces_g,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(forces_r,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE build_fock_matrix_coulomb_lr

! *****************************************************************************
!> \brief When doing long-range SE calculation, this module computes the correction
!>        between the mismatch of point-like multipoles and multipoles represented
!>        with charges
!> \author Teodoro Laino [tlaino] - 05.2009
! *****************************************************************************
  SUBROUTINE build_fock_matrix_coul_lrc (qs_env,ks_matrix,matrix_p,energy,&
       calculate_forces, store_int_env,se_nddo_mpole,task,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_p
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(nddo_mpole_type), POINTER           :: se_nddo_mpole
    LOGICAL, DIMENSION(3), INTENT(IN)        :: task
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_fock_matrix_coul_lrc', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iab, iatom, ikind, ilist, inode, &
      ispin, itype, jatom, jkind, natom, natorb_a, natorb_a2, natorb_b, &
      natorb_b2, nkind, nlist, nnode, nspins, size1, size2, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: anag, check, defined, &
                                                failure, scp_nddo, switch, &
                                                use_virial
    REAL(KIND=dp)                            :: delta, dr1, ecore2, enuc
    REAL(KIND=dp), DIMENSION(2025)           :: pa_a, pa_b, pb_a, pb_b
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, force_ab0, rij
    REAL(KIND=dp), DIMENSION(45, 45)         :: p_block_tot_a, p_block_tot_b
    REAL(KIND=dp), DIMENSION(:), POINTER     :: efield0
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: efield1, efield2, ksa_block_a, &
      ksa_block_b, ksb_block_a, ksb_block_b, pa_block_a, pa_block_b, &
      pb_block_a, pb_block_b, pscpa_block, pscpb_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_lrc
    TYPE(neighbor_list_type), POINTER        :: sab_lrc_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_lrc_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: diagmat_ks, diagmat_p
    TYPE(real_matrix_type), POINTER          :: ks_scp, pscp
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(se_int_control_type)                :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(virial_type), POINTER               :: virial

    failure=.FALSE.
    ! Temporary Hack.. this routine requires debugging...
    RETURN
    ! 
    CALL timeset(routineN,handle)
    IF ( .NOT. failure ) THEN
       NULLIFY(dft_control,cell,force,particle_set,diagmat_ks,diagmat_p,&
               se_control,se_taper,pscp,ks_scp,scp_env,efield0,efield1,&
               efield2)
       
       CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, se_taper=se_taper,&
            para_env=para_env, sab_lrc=sab_lrc, atomic_kind_set=atomic_kind_set, &
            particle_set=particle_set, virial=virial, error=error)
      
       ! Parameters
       CALL initialize_se_taper(se_taper,lr_corr=.TRUE.,error=error)
       se_control  => dft_control%qs_control%se_control
       anag        =  se_control%analytical_gradients
       scp_nddo    =  se_control%scp
       use_virial  = virial%pv_availability.AND.(.NOT.virial%pv_numer)

       CALL setup_se_int_control_type(se_int_control, do_ewald_r3=se_control%do_ewald_r3,&
            do_ewald_gks=.FALSE., integral_screening=se_control%integral_screening,&
            shortrange=.FALSE., max_multipole=se_control%max_multipole,&
            pc_coulomb_int=.TRUE.)

       nspins=dft_control%nspins
       CPPrecondition(ASSOCIATED(matrix_p),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)
       
       CALL allocate_matrix_set(diagmat_ks,nspins,error=error)
       CALL allocate_matrix_set(diagmat_p,nspins,error=error)
       
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

       DO ispin=1,nspins
          ! Allocate diagonal block matrices
          CALL create_dbm_from_sm(diagmat_p(ispin)%matrix,matrix_p(ispin)%matrix,&
                                  "DIAGONAL BLOCK P MATRIX",.TRUE.,para_env,error=error)
          CALL create_dbm_from_sm(diagmat_ks(ispin)%matrix,ks_matrix(ispin)%matrix,&
                                  "DIAGONAL BLOCK KS MATRIX",.FALSE.,error=error)
       END DO

       ! SCP
       IF (scp_nddo)  THEN 
          CALL get_qs_env(qs_env=qs_env, scp_env=scp_env, error=error)
          CALL get_scp_env(scp_env=scp_env, pscp=pscp, ks_scp=ks_scp, error=error)
       END IF
       
       ! Allocate arrays for storing partial information on potential, field, field gradient
       size1 = SIZE(se_nddo_mpole%efield0)
       ALLOCATE (efield0(size1), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       efield0 = 0.0_dp
       size1 = SIZE(se_nddo_mpole%efield1,1)
       size2 = SIZE(se_nddo_mpole%efield1,2)
       ALLOCATE (efield1(size1,size2), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       efield1 = 0.0_dp
       size1 = SIZE(se_nddo_mpole%efield2,1)
       size2 = SIZE(se_nddo_mpole%efield2,2)
       ALLOCATE (efield2(size1,size2), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       efield2 = 0.0_dp

       ! Start of the loop for the correction of the pair interactions
       ecore2   = 0.0_dp
       itype    = get_se_type(dft_control%qs_control%method_id)
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a)
 
          IF (.NOT.defined .OR. natorb_a < 1) CYCLE
          natorb_a2 = natorb_a**2
          DO jkind=1,nkind
             atomic_kind => atomic_kind_set(jkind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_b)
             CALL get_se_param(se_kind_b,defined=defined,natorb=natorb_b)

             IF (.NOT.defined .OR. natorb_b < 1) CYCLE
             natorb_b2 = natorb_b**2
             iab = ikind + nkind*(jkind - 1)
             IF (.NOT.ASSOCIATED(sab_lrc(iab)%neighbor_list_set)) CYCLE
             CALL get_neighbor_list_set(neighbor_list_set=sab_lrc(iab)%neighbor_list_set,nlist=nlist)
             NULLIFY(sab_lrc_neighbor_list)

             DO ilist=1,nlist
                IF ( .NOT. ASSOCIATED(sab_lrc_neighbor_list) ) THEN
                   sab_lrc_neighbor_list => first_list(sab_lrc(iab)%neighbor_list_set)
                ELSE
                   sab_lrc_neighbor_list => next(sab_lrc_neighbor_list)
                END IF
                CALL get_neighbor_list(neighbor_list=sab_lrc_neighbor_list,atom=iatom, nnode=nnode)
                sab_lrc_neighbor_node => first_node(sab_lrc_neighbor_list)
 
                CALL get_block_node(matrix=diagmat_p(1)%matrix,&
                                    block_row=iatom,&
                                    block_col=iatom,&
                                    BLOCK=pa_block_a)
                CPPostcondition(ASSOCIATED(pa_block_a),cp_failure_level,routineP,error,failure)
                check = (SIZE(pa_block_a,1)==natorb_a).AND.(SIZE(pa_block_a,2)==natorb_a)
                CPPostcondition(check,cp_failure_level,routineP,error,failure)
                CALL get_block_node(matrix=diagmat_ks(1)%matrix,&
                                    block_row=iatom,&
                                    block_col=iatom,&
                                    BLOCK=ksa_block_a)
                CPPostcondition(ASSOCIATED(ksa_block_a),cp_failure_level,routineP,error,failure)
                p_block_tot_a(1:natorb_a,1:natorb_a) = 2.0_dp * pa_block_a
                pa_a(1:natorb_a2)                    = RESHAPE(pa_block_a,(/natorb_a2/)) 
                IF ( nspins == 2 ) THEN
                   CALL get_block_node(matrix=diagmat_p(2)%matrix,&
                                       block_row=iatom,&
                                       block_col=iatom,&
                                       BLOCK=pa_block_b)
                   CPPostcondition(ASSOCIATED(pa_block_b),cp_failure_level,routineP,error,failure)
                   check = (SIZE(pa_block_b,1)==natorb_a).AND.(SIZE(pa_block_b,2)==natorb_a)
                   CPPostcondition(check,cp_failure_level,routineP,error,failure)
                   CALL get_block_node(matrix=diagmat_ks(2)%matrix,&
                                       block_row=iatom,&
                                       block_col=iatom,&
                                       BLOCK=ksa_block_b)
                   CPPostcondition(ASSOCIATED(ksa_block_b),cp_failure_level,routineP,error,failure)
                   p_block_tot_a(1:natorb_a,1:natorb_a) = pa_block_a + pa_block_b
                   pa_b(1:natorb_a2)                    = RESHAPE(pa_block_b,(/natorb_a2/))
                END IF

                IF (scp_nddo) THEN
                   CALL get_block_node(matrix=pscp,&
                                       block_row=iatom,&
                                       block_col=iatom,&
                                       BLOCK=pscpa_block)
                   CPPostcondition(ASSOCIATED(pscpa_block),cp_failure_level,routineP,error,failure)
                   check = (SIZE(pscpa_block,1)==natorb_a).AND.(SIZE(pscpa_block,2)==natorb_a)
                   CPPostcondition(check,cp_failure_level,routineP,error,failure)
                   ! Add P_scp + P_block_tot_a to get P_tot
                   p_block_tot_a(1:natorb_a,1:natorb_a) = p_block_tot_a(1:natorb_a,1:natorb_a) + 2.0_dp * pscpa_block
                   ! Modify accordingly the density matrix for SCP
                        pa_a(1:natorb_a2)       = RESHAPE(pa_block_a+pscpa_block,(/natorb_a2/))
                   IF ( nspins == 2 ) &
                        pa_b(1:natorb_a2)       = RESHAPE(pa_block_b+pscpa_block,(/natorb_a2/))
                ENDIF
 
                DO inode=1,nnode
                   CALL get_neighbor_node(neighbor_node=sab_lrc_neighbor_node,neighbor=jatom,r=rij)
 
                   dr1 = DOT_PRODUCT(rij,rij)
                   IF ( dr1 > rij_threshold ) THEN
                      ! Determine the order of the atoms, and in case switch them..
                      IF (iatom <= jatom) THEN
                         switch = .FALSE.
                      ELSE
                         switch = .TRUE.
                      END IF

                      ! Pint-like interaction corrections
                      CALL se_coulomb_ij_interaction (iatom, jatom, task, do_forces=calculate_forces,&
                           do_efield=.TRUE., charges=se_nddo_mpole%charge, dipoles=se_nddo_mpole%dipole,&
                           quadrupoles=se_nddo_mpole%quadrupole, force_ab=force_ab0, efield0=efield0,&
                           efield1=efield1, efield2=efield2, rab2=dr1, rab=rij, error=error)
                      
                      ! Retrieve blocks for KS and P
                      CALL get_block_node(matrix=diagmat_p(1)%matrix,&
                                          block_row=jatom,&
                                          block_col=jatom,&
                                          BLOCK=pb_block_a)
                      CPPostcondition(ASSOCIATED(pb_block_a),cp_failure_level,routineP,error,failure)
                      check = (SIZE(pb_block_a,1)==natorb_b).AND.(SIZE(pb_block_a,2)==natorb_b)
                      CPPostcondition(check,cp_failure_level,routineP,error,failure)
                      CALL get_block_node(matrix=diagmat_ks(1)%matrix,&
                                          block_row=jatom,&
                                          block_col=jatom,&
                                          BLOCK=ksb_block_a)
                      CPPostcondition(ASSOCIATED(ksb_block_a),cp_failure_level,routineP,error,failure)
                      p_block_tot_b(1:natorb_b,1:natorb_b) = 2.0_dp * pb_block_a
                      pb_a(1:natorb_b2)                    = RESHAPE(pb_block_a,(/natorb_b2/)) 
                      ! Handle more than one configuration
                      IF ( nspins == 2 ) THEN
                         CALL get_block_node(matrix=diagmat_p(2)%matrix,&
                                             block_row=jatom,&
                                             block_col=jatom,&
                                             BLOCK=pb_block_b)
                         CPPostcondition(ASSOCIATED(pb_block_b),cp_failure_level,routineP,error,failure)
                         check = (SIZE(pb_block_b,1)==natorb_b).AND.(SIZE(pb_block_b,2)==natorb_b)
                         CPPostcondition(check,cp_failure_level,routineP,error,failure)
                         CALL get_block_node(matrix=diagmat_ks(2)%matrix,&
                                             block_row=jatom,&
                                             block_col=jatom,&
                                             BLOCK=ksb_block_b)
                         CPPostcondition(ASSOCIATED(ksb_block_b),cp_failure_level,routineP,error,failure)
                         check = (SIZE(pb_block_a,1)==SIZE(pb_block_b,1)).AND.(SIZE(pb_block_a,2)==SIZE(pb_block_b,2))
                         CPPostcondition(check,cp_failure_level,routineP,error,failure)
                         p_block_tot_b(1:natorb_b,1:natorb_b) = pb_block_a + pb_block_b
                         pb_b(1:natorb_b2)                    = RESHAPE(pb_block_b,(/natorb_b2/))
                      END IF

                      ! SCP
                      IF ( scp_nddo ) THEN
                         CALL get_block_node(matrix=pscp,&
                                             block_row=jatom,&
                                             block_col=jatom,&
                                             BLOCK=pscpb_block)
                         CPPostcondition(ASSOCIATED(pscpb_block),cp_failure_level,routineP,error,failure)
                         check = (SIZE(pscpb_block,1)==natorb_b).AND.(SIZE(pscpb_block,2)==natorb_b)
                         CPPostcondition(check,cp_failure_level,routineP,error,failure)
                         ! Add P_scp + P_block_tot_b to get P_tot
                         p_block_tot_b(1:natorb_b,1:natorb_b) = p_block_tot_b(1:natorb_b,1:natorb_b) + 2.0_dp * pscpb_block
                         ! Modify accordingly the density matrix for SCP
                              pb_a(1:natorb_b2) = RESHAPE(pb_block_a+pscpb_block,(/natorb_b2/))
                         IF ( nspins == 2 ) &
                              pb_b(1:natorb_b2) = RESHAPE(pb_block_b+pscpb_block,(/natorb_b2/))
                      ENDIF 
 
                      SELECT CASE (dft_control%qs_control%method_id)
                      CASE (do_method_mndo, do_method_am1, do_method_pm3, do_method_pm6, do_method_pdg,&
                            do_method_rm1, do_method_mndod, do_method_pnnl)
                         ! Evaluate nuclear contribution..
                         CALL corecore_el (se_kind_a,se_kind_b,rij,enuc=enuc,itype=itype,anag=anag,&
                              se_int_control=se_int_control, se_taper=se_taper, error=error)
                         ecore2 = ecore2 + enuc

                         ! Two-centers One-electron terms
                         IF      ( nspins == 1 ) THEN
                            CALL fock2_1el (se_kind_a,se_kind_b,rij, ksa_block_a, ksb_block_a,&
                                 pa_a, pb_a, ecore=ecore2, itype=itype, anag=anag, se_int_control=se_int_control,&
                                 se_taper=se_taper, store_int_env=store_int_env, error=error)
                         ELSE IF ( nspins == 2 ) THEN
                            CALL fock2_1el (se_kind_a,se_kind_b,rij, ksa_block_a, ksb_block_a,&
                                 pa_block_a, pb_block_a, ecore=ecore2, itype=itype, anag=anag,&
                                 se_int_control=se_int_control, se_taper=se_taper, store_int_env=store_int_env,&
                                 error=error)
                            CALL fock2_1el (se_kind_a,se_kind_b,rij, ksa_block_b, ksb_block_b,&
                                 pa_b, pb_b, ecore=ecore2, itype=itype, anag=anag, se_int_control=se_int_control,&
                                 se_taper=se_taper, store_int_env=store_int_env, error=error)
                         END IF
                         ! Coulomb Terms
                         IF      ( nspins == 1 ) THEN
                            CALL fock2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, ksa_block_a, p_block_tot_b,&
                                 ksb_block_a, factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                                 store_int_env=store_int_env, error=error)
                         ELSE IF ( nspins == 2 ) THEN
                            CALL fock2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, ksa_block_a, p_block_tot_b,&
                                 ksb_block_a, factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                                 store_int_env=store_int_env, error=error)
 
                            CALL fock2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, ksa_block_b, p_block_tot_b,&
                                 ksb_block_b, factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper,&
                                 store_int_env=store_int_env, error=error)
                         END IF
 
                         IF(calculate_forces) THEN
                            atom_a   = atom_of_kind(iatom)
                            atom_b   = atom_of_kind(jatom)

                            ! Evaluate nuclear contribution..
                            CALL dcorecore_el (se_kind_a,se_kind_b,rij,denuc=force_ab,itype=itype,delta=delta,&
                                 anag=anag,se_int_control=se_int_control,se_taper=se_taper,error=error)

                            ! Derivatives of the Two-centre One-electron terms
                            IF      ( nspins == 1 ) THEN
                               CALL dfock2_1el (se_kind_a,se_kind_b,rij, pa_a, pb_a, itype=itype, anag=anag,&
                                    se_int_control=se_int_control, se_taper=se_taper, force=force_ab,&
                                    delta=delta, error=error)
                            ELSE IF ( nspins == 2 ) THEN
                               CALL dfock2_1el (se_kind_a,se_kind_b,rij, pa_block_a, pb_block_a, itype=itype, anag=anag,&
                                    se_int_control=se_int_control, se_taper=se_taper, force=force_ab, delta=delta,&
                                    error=error)
                               CALL dfock2_1el (se_kind_a,se_kind_b,rij, pa_b, pb_b, itype=itype, anag=anag,&
                                    se_int_control=se_int_control, se_taper=se_taper, force=force_ab, delta=delta,&
                                    error=error)
                            END IF
                            force_ab = force_ab + force_ab0

                            IF (use_virial) THEN
                               CALL virial_pair_force ( virial%pv_virial, -1.0_dp, force_ab, rij, error)
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
                               CALL dfock2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, p_block_tot_b, factor=0.25_dp,&
                                    anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab, delta=delta,&
                                    error=error)
                            ELSE IF ( nspins == 2 ) THEN
                               CALL dfock2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, p_block_tot_b, factor=0.50_dp,&
                                    anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab, delta=delta,&
                                    error=error)
                               
                               CALL dfock2C(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, p_block_tot_b, factor=0.50_dp,&
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
                   END IF
                   sab_lrc_neighbor_node => next(sab_lrc_neighbor_node)
                END DO ! inode => jatom(atom B)
             END DO ! ilist => iatom(atom A)
          END DO ! jkind
       END DO ! ikind

       IF ( scp_nddo ) THEN
          ! Update the SCP coulomb fock matrix.  Should be K = (P + p )C + CZ
          CALL transfer_matrix ( diagmat_ks(1)%matrix, ks_scp, error )
          ! sum over all procs to get the replicated ks_scp
          CALL mpsum_dbm ( ks_scp, para_env, error )
          CALL set_scp_env ( scp_env=scp_env, ks_scp = ks_scp, error = error ) 
       ENDIF 
       
       DO ispin = 1, nspins
          CALL deallocate_dbm_from_sm(diagmat_p(ispin)%matrix,matrix_p(ispin)%matrix,.FALSE.,error=error)
          CALL deallocate_dbm_from_sm(diagmat_ks(ispin)%matrix,ks_matrix(ispin)%matrix,.TRUE.,para_env,error=error)
       END DO
       DEALLOCATE(diagmat_ks,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(diagmat_p,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       
       IF (calculate_forces) THEN 
          DEALLOCATE(atom_of_kind,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       ! collect information on potential, field, field gradient
       CALL mp_sum(efield0, para_env%group)
       CALL mp_sum(efield1, para_env%group)
       CALL mp_sum(efield2, para_env%group)
       se_nddo_mpole%efield0 = se_nddo_mpole%efield0 - efield0
       se_nddo_mpole%efield1 = se_nddo_mpole%efield1 - efield1
       se_nddo_mpole%efield2 = se_nddo_mpole%efield2 - efield2
       ! deallocate working arrays
       DEALLOCATE (efield0, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (efield1, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (efield2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       
       ! Corrections for two-centers one-electron terms + nuclear terms
       CALL mp_sum(ecore2,para_env%group)
       energy%hartree = energy%hartree + ecore2
    END IF
    CALL finalize_se_taper(se_taper,error=error)    
    CALL timestop(handle)
  END SUBROUTINE build_fock_matrix_coul_lrc

! *****************************************************************************
!> \brief Construction of the residual part (1/R^3) of the Coulomb long-range 
!>        term of the Fock matrix
!>        The 1/R^3 correction works in real-space strictly on the zero-cell,
!>        in order to avoid more parameters to be provided in the input..
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
  SUBROUTINE build_fock_matrix_coul_lr_r3(qs_env,ks_matrix,matrix_p,energy,calculate_forces,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_p
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_fock_matrix_coul_lr_r3', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, ewald_type, handle, iab, iatom, ikind, ilist, &
      inode, ispin, itype, jatom, jkind, natoms, natorb_a, natorb_a2, &
      natorb_b, natorb_b2, nkind, nlist, nlocal_particles, nnode, nspins, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: anag, check, defined, &
                                                failure, switch
    REAL(KIND=dp)                            :: dr1, ecore2, r2inv, r3inv, &
                                                rinv
    REAL(KIND=dp), DIMENSION(2025)           :: pa_a, pa_b, pb_a, pb_b
    REAL(KIND=dp), DIMENSION(3)              :: dr3inv, force_ab, rij
    REAL(KIND=dp), DIMENSION(45, 45)         :: p_block_tot_a, p_block_tot_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: ksa_block_a, ksa_block_b, &
      ksb_block_a, ksb_block_b, pa_block_a, pa_block_b, pb_block_a, pb_block_b
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(fist_nonbond_env_type), POINTER     :: se_nonbond_env
    TYPE(nddo_mpole_type), POINTER           :: se_nddo_mpole
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: diagmat_ks, diagmat_p
    TYPE(section_vals_type), POINTER         :: se_section
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b

    failure=.FALSE.
    CALL timeset(routineN,handle)
    IF ( .NOT. failure ) THEN
       NULLIFY(dft_control, cell, force, particle_set, diagmat_ks, &
            diagmat_p, local_particles, se_control, ewald_env, ewald_pw, &
            se_nddo_mpole, se_nonbond_env, se_section, sab_orb, logger)

       logger => cp_error_get_logger(error)
       CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, para_env=para_env,&
            atomic_kind_set=atomic_kind_set, particle_set=particle_set, ewald_env=ewald_env,&
            local_particles=local_particles, ewald_pw=ewald_pw, se_nddo_mpole=se_nddo_mpole,&
            se_nonbond_env=se_nonbond_env, sab_orb=sab_orb, error=error)
       
       nlocal_particles = SUM(local_particles%n_el(:))
       natoms           = SIZE(particle_set)
       CALL ewald_env_get (ewald_env, ewald_type=ewald_type, error=error)
       SELECT CASE(ewald_type)
       CASE (do_ewald_ewald)
          ! Do Nothing
       CASE DEFAULT
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Periodic SE implemented only for standard EWALD sums.", &
               error=error, error_level=cp_failure_level)
       END SELECT

       ! Parameters
       se_section => section_vals_get_subs_vals(qs_env%input,"DFT%QS%SE",error=error)
       se_control => dft_control%qs_control%se_control
       anag       =  se_control%analytical_gradients
       
       nspins=dft_control%nspins
       CPPrecondition(ASSOCIATED(matrix_p),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)
       
       CALL allocate_matrix_set(diagmat_ks,nspins,error=error)
       CALL allocate_matrix_set(diagmat_p,nspins,error=error)
       
       nkind = SIZE(atomic_kind_set)
       DO ispin=1,nspins
          ! Allocate diagonal block matrices
          CALL create_dbm_from_sm(diagmat_p(ispin)%matrix,matrix_p(ispin)%matrix,&
                                  "DIAGONAL BLOCK P MATRIX",.TRUE.,para_env,error=error)
          CALL create_dbm_from_sm(diagmat_ks(ispin)%matrix,ks_matrix(ispin)%matrix,&
                                  "DIAGONAL BLOCK KS MATRIX",.FALSE.,error=error)
       END DO

       ! Possibly compute forces
       IF(calculate_forces) THEN
          CALL get_qs_env(qs_env=qs_env,force=force,error=error)
          ALLOCATE (atom_of_kind(natoms),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)
       END IF
       itype = get_se_type(dft_control%qs_control%method_id)

       ecore2   = 0.0_dp
       ! Real space part of the 1/R^3 sum
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a)
          
          IF (.NOT.defined .OR. natorb_a < 1) CYCLE
          natorb_a2 = natorb_a**2
          DO jkind=1,nkind
             atomic_kind => atomic_kind_set(jkind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_b)
             CALL get_se_param(se_kind_b,defined=defined,natorb=natorb_b)

             IF (.NOT.defined .OR. natorb_b < 1) CYCLE
             natorb_b2 = natorb_b**2
             iab = ikind + nkind*(jkind - 1)
             IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE
             CALL get_neighbor_list_set(neighbor_list_set=sab_orb(iab)%neighbor_list_set,nlist=nlist)
             NULLIFY(sab_orb_neighbor_list)

             DO ilist=1,nlist
                IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
                   sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
                ELSE
                   sab_orb_neighbor_list => next(sab_orb_neighbor_list)
                END IF
                CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,atom=iatom, nnode=nnode)
                sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)
 
                CALL get_block_node(matrix=diagmat_p(1)%matrix,&
                                    block_row=iatom,&
                                    block_col=iatom,&
                                    BLOCK=pa_block_a)
                CPPostcondition(ASSOCIATED(pa_block_a),cp_failure_level,routineP,error,failure)
                check = (SIZE(pa_block_a,1)==natorb_a).AND.(SIZE(pa_block_a,2)==natorb_a)
                CPPostcondition(check,cp_failure_level,routineP,error,failure)
                CALL get_block_node(matrix=diagmat_ks(1)%matrix,&
                                    block_row=iatom,&
                                    block_col=iatom,&
                                    BLOCK=ksa_block_a)
                CPPostcondition(ASSOCIATED(ksa_block_a),cp_failure_level,routineP,error,failure)
                p_block_tot_a(1:natorb_a,1:natorb_a) = 2.0_dp * pa_block_a
                pa_a(1:natorb_a2)                    = RESHAPE(pa_block_a,(/natorb_a2/)) 
                IF ( nspins == 2 ) THEN
                   CALL get_block_node(matrix=diagmat_p(2)%matrix,&
                                       block_row=iatom,&
                                       block_col=iatom,&
                                       BLOCK=pa_block_b)
                   CPPostcondition(ASSOCIATED(pa_block_b),cp_failure_level,routineP,error,failure)
                   check = (SIZE(pa_block_b,1)==natorb_a).AND.(SIZE(pa_block_b,2)==natorb_a)
                   CPPostcondition(check,cp_failure_level,routineP,error,failure)
                   CALL get_block_node(matrix=diagmat_ks(2)%matrix,&
                                       block_row=iatom,&
                                       block_col=iatom,&
                                       BLOCK=ksa_block_b)
                   CPPostcondition(ASSOCIATED(ksa_block_b),cp_failure_level,routineP,error,failure)
                   p_block_tot_a(1:natorb_a,1:natorb_a) = pa_block_a + pa_block_b
                   pa_b(1:natorb_a2)                    = RESHAPE(pa_block_b,(/natorb_a2/))
                END IF
 
                DO inode=1,nnode
                   CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,neighbor=jatom,r=rij)
                   dr1 = DOT_PRODUCT(rij,rij)
                   IF ( dr1 > rij_threshold ) THEN
                      ! Determine the order of the atoms, and in case switch them..
                      IF (iatom <= jatom) THEN
                         switch = .FALSE.
                      ELSE
                         switch = .TRUE.
                      END IF
                      ! Retrieve blocks for KS and P
                      CALL get_block_node(matrix=diagmat_p(1)%matrix,&
                                          block_row=jatom,&
                                          block_col=jatom,&
                                          BLOCK=pb_block_a)
                      CPPostcondition(ASSOCIATED(pb_block_a),cp_failure_level,routineP,error,failure)
                      check = (SIZE(pb_block_a,1)==natorb_b).AND.(SIZE(pb_block_a,2)==natorb_b)
                      CPPostcondition(check,cp_failure_level,routineP,error,failure)
                      CALL get_block_node(matrix=diagmat_ks(1)%matrix,&
                                          block_row=jatom,&
                                          block_col=jatom,&
                                          BLOCK=ksb_block_a)
                      CPPostcondition(ASSOCIATED(ksb_block_a),cp_failure_level,routineP,error,failure)
                      p_block_tot_b(1:natorb_b,1:natorb_b) = 2.0_dp * pb_block_a
                      pb_a(1:natorb_b2)                    = RESHAPE(pb_block_a,(/natorb_b2/)) 
                      ! Handle more than one configuration
                      IF ( nspins == 2 ) THEN
                         CALL get_block_node(matrix=diagmat_p(2)%matrix,&
                                             block_row=jatom,&
                                             block_col=jatom,&
                                             BLOCK=pb_block_b)
                         CPPostcondition(ASSOCIATED(pb_block_b),cp_failure_level,routineP,error,failure)
                         check = (SIZE(pb_block_b,1)==natorb_b).AND.(SIZE(pb_block_b,2)==natorb_b)
                         CPPostcondition(check,cp_failure_level,routineP,error,failure)
                         CALL get_block_node(matrix=diagmat_ks(2)%matrix,&
                                             block_row=jatom,&
                                             block_col=jatom,&
                                             BLOCK=ksb_block_b)
                         CPPostcondition(ASSOCIATED(ksb_block_b),cp_failure_level,routineP,error,failure)
                         check = (SIZE(pb_block_a,1)==SIZE(pb_block_b,1)).AND.(SIZE(pb_block_a,2)==SIZE(pb_block_b,2))
                         CPPostcondition(check,cp_failure_level,routineP,error,failure)
                         p_block_tot_b(1:natorb_b,1:natorb_b) = pb_block_a + pb_block_b
                         pb_b(1:natorb_b2)                    = RESHAPE(pb_block_b,(/natorb_b2/))
                      END IF
 
                      SELECT CASE (dft_control%qs_control%method_id)
                      CASE (do_method_mndo, do_method_am1, do_method_pm3, do_method_pm6, do_method_pdg,&
                            do_method_rm1, do_method_mndod, do_method_pnnl)

                         ! Pre-compute some quantities..
                         r2inv = 1.0_dp/dr1
                         rinv  = SQRT(r2inv)
                         r3inv = rinv**3

                         ! Two-centers One-electron terms
                         IF      ( nspins == 1 ) THEN
                            CALL fock2_1el_r3(se_kind_a, se_kind_b, ksa_block_a, ksb_block_a, pa_a, pb_a,&
                                 ecore=ecore2, e1b=se_kind_a%expns3_int(jkind)%expns3%e1b,&
                                 e2a=se_kind_a%expns3_int(jkind)%expns3%e2a, rp=r3inv, error=error)
                         ELSE IF ( nspins == 2 ) THEN
                            CALL fock2_1el_r3(se_kind_a, se_kind_b, ksa_block_a, ksb_block_a, pa_block_a,&
                                 pb_block_a, ecore=ecore2, e1b=se_kind_a%expns3_int(jkind)%expns3%e1b,&
                                 e2a=se_kind_a%expns3_int(jkind)%expns3%e2a, rp=r3inv, error=error)

                            CALL fock2_1el_r3(se_kind_a, se_kind_b, ksa_block_b, ksb_block_b, pa_b, pb_b,&
                                 ecore=ecore2, e1b=se_kind_a%expns3_int(jkind)%expns3%e1b,&
                                 e2a=se_kind_a%expns3_int(jkind)%expns3%e2a, rp=r3inv, error=error)
                         END IF
                         ! Coulomb Terms
                         IF      ( nspins == 1 ) THEN
                            CALL fock2C_r3(se_kind_a, se_kind_b, switch, p_block_tot_a, ksa_block_a, p_block_tot_b,&
                                 ksb_block_a, factor=0.5_dp, w=se_kind_a%expns3_int(jkind)%expns3%w, rp=r3inv,&
                                 error=error)
                         ELSE IF ( nspins == 2 ) THEN
                            CALL fock2C_r3(se_kind_a, se_kind_b, switch, p_block_tot_a, ksa_block_a, p_block_tot_b,&
                                 ksb_block_a, factor=1.0_dp, w=se_kind_a%expns3_int(jkind)%expns3%w, rp=r3inv,&
                                 error=error)
 
                            CALL fock2C_r3(se_kind_a, se_kind_b, switch, p_block_tot_a, ksa_block_b, p_block_tot_b,&
                                 ksb_block_b, factor=1.0_dp, w=se_kind_a%expns3_int(jkind)%expns3%w, rp=r3inv,&
                                 error=error)
                         END IF

                         ! Compute forces if requested
                         IF(calculate_forces) THEN
                            dr3inv   = -3.0_dp*rij*r3inv*r2inv
                            atom_a   = atom_of_kind(iatom)
                            atom_b   = atom_of_kind(jatom)

                            force_ab = 0.0_dp
                            ! Derivatives of the One-centre One-electron terms
                            IF      ( nspins == 1 ) THEN
                               CALL dfock2_1el_r3(se_kind_a,se_kind_b, dr3inv, pa_a, pb_a, force_ab,&
                                    se_kind_a%expns3_int(jkind)%expns3%e1b, se_kind_a%expns3_int(jkind)%expns3%e2a,&
                                    error=error)
                            ELSE IF ( nspins == 2 ) THEN
                               CALL dfock2_1el_r3(se_kind_a,se_kind_b, dr3inv, pa_block_a, pb_block_a, force_ab,&
                                    se_kind_a%expns3_int(jkind)%expns3%e1b, se_kind_a%expns3_int(jkind)%expns3%e2a,&
                                    error=error)

                               CALL dfock2_1el_r3(se_kind_a,se_kind_b, dr3inv, pa_b, pb_b, force_ab,&
                                    se_kind_a%expns3_int(jkind)%expns3%e1b, se_kind_a%expns3_int(jkind)%expns3%e2a,&
                                    error=error)
                            END IF

                            ! Sum up force components
                            force(ikind)%all_potential(1,atom_a) = force(ikind)%all_potential(1,atom_a) - force_ab(1)
                            force(jkind)%all_potential(1,atom_b) = force(jkind)%all_potential(1,atom_b) + force_ab(1)

                            force(ikind)%all_potential(2,atom_a) = force(ikind)%all_potential(2,atom_a) - force_ab(2)
                            force(jkind)%all_potential(2,atom_b) = force(jkind)%all_potential(2,atom_b) + force_ab(2)

                            force(ikind)%all_potential(3,atom_a) = force(ikind)%all_potential(3,atom_a) - force_ab(3)
                            force(jkind)%all_potential(3,atom_b) = force(jkind)%all_potential(3,atom_b) + force_ab(3)
 
                            ! Derivatives of the Coulomb Terms
                            force_ab = 0.0_dp
                            IF      ( nspins == 1 ) THEN
                               CALL dfock2C_r3(se_kind_a, se_kind_b, switch, p_block_tot_a, p_block_tot_b, factor=0.25_dp,&
                                    w=se_kind_a%expns3_int(jkind)%expns3%w, drp=dr3inv, force=force_ab, error=error)
                            ELSE IF ( nspins == 2 ) THEN
                               CALL dfock2C_r3(se_kind_a, se_kind_b, switch, p_block_tot_a, p_block_tot_b, factor=0.50_dp,&
                                    w=se_kind_a%expns3_int(jkind)%expns3%w, drp=dr3inv, force=force_ab, error=error)
                               
                               CALL dfock2C_r3(se_kind_a, se_kind_b, switch, p_block_tot_a, p_block_tot_b, factor=0.50_dp,&
                                    w=se_kind_a%expns3_int(jkind)%expns3%w, drp=dr3inv, force=force_ab, error=error)
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
                   END IF
                   sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                END DO ! inode => jatom(atom B)
             END DO ! ilist => iatom(atom A)
          END DO ! jkind
       END DO ! ikind
       
       DO ispin = 1, nspins
          CALL deallocate_dbm_from_sm(diagmat_p(ispin)%matrix,matrix_p(ispin)%matrix,.FALSE.,error=error)
          CALL deallocate_dbm_from_sm(diagmat_ks(ispin)%matrix,ks_matrix(ispin)%matrix,.TRUE.,para_env,error=error)
       END DO
       DEALLOCATE(diagmat_ks,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(diagmat_p,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       
       IF (calculate_forces) THEN 
          DEALLOCATE(atom_of_kind,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       
       ! Two-centers one-electron terms
       CALL mp_sum(ecore2,para_env%group)       
       energy%hartree = energy%hartree + ecore2
    END IF
    CALL timestop(handle)
  END SUBROUTINE build_fock_matrix_coul_lr_r3
  
END MODULE se_fock_matrix_coulomb

