!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the Fock matrix for NDDO methods
!> \author JGH
!> \par History
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : d-orbitals
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Speed-up
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Periodic SE
! *****************************************************************************
MODULE nddo_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             semi_empirical_control_type
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
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
       do_multipole_none, do_multipole_quadrupole, do_se_IS_kdso, &
       do_se_IS_kdso_d
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type,&
                                             set_scp_env
  USE scp_nddo_utils,                  ONLY: evaluate_scp_force
  USE semi_empirical_int_arrays,       ONLY: rij_threshold,&
                                             se_orbital_pointer
  USE semi_empirical_integrals,        ONLY: corecore,&
                                             dcorecore,&
                                             drotint,&
                                             drotnuc,&
                                             rotint,&
                                             rotnuc
  USE semi_empirical_mpole_methods,    ONLY: quadrupole_sph_to_cart
  USE semi_empirical_mpole_types,      ONLY: nddo_mpole_type,&
                                             semi_empirical_mpole_type
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_finalize,&
                                             semi_empirical_si_initialize,&
                                             semi_empirical_si_type
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             se_int_control_type,&
                                             se_taper_type,&
                                             semi_empirical_type,&
                                             setup_se_int_control_type
  USE semi_empirical_utils,            ONLY: finalize_se_taper,&
                                             get_se_type,&
                                             initialize_se_taper
  USE sparse_matrix_types,             ONLY: &
       allocate_matrix_set, cp_sm_scale_and_add, cp_sm_sm_trace, &
       create_dbm_from_sm, deallocate_dbm_from_sm, deallocate_matrix, &
       get_block_node, mpsum_dbm, real_matrix_p_type, real_matrix_type, &
       set_matrix, transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nddo_methods'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
  LOGICAL, PARAMETER, PRIVATE          :: debug_coulomb_lr  = .FALSE.
  
  PUBLIC :: build_nddo_fock_matrix, se_core_core_interaction

CONTAINS

! *****************************************************************************
!> \brief Construction of the Fock matrix for NDDO methods
!> \author JGH
! *****************************************************************************
  SUBROUTINE build_nddo_fock_matrix (ks_env,qs_env,ks_matrix,rho,energy,&
       calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_nddo_fock_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, ncol_global, &
                                                nspins, output_unit
    LOGICAL                                  :: failure, scp_nddo
    REAL(dp), POINTER                        :: ascp_diag( : )
    REAL(KIND=dp)                            :: ecoul, ecoul_scp, qmmm_el
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p
    TYPE(real_matrix_type), POINTER          :: ks_scp, pscp
    TYPE(scp_energy_type), POINTER           :: energy_scp
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(matrix_h, dft_control, logger, scf_section, store_int_env, se_control)
    NULLIFY(scp_env, pscp, ks_scp, ascp_diag, energy_scp,atomic_kind_set)
    logger => cp_error_get_logger(error)
    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
    
    IF (.NOT.failure) THEN
       CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, matrix_h=matrix_h,&
            para_env=para_env, se_store_int_env=store_int_env, scp_env=scp_env,&
            atomic_kind_set=atomic_kind_set, error=error)

       se_control  => dft_control%qs_control%se_control
       scp_nddo    =  se_control%scp

       scf_section => section_vals_get_subs_vals(qs_env%input,"DFT%SCF",error=error)
       CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)
       nspins         = dft_control%nspins
       matrix_p       => rho%rho_ao
       energy%qmmm_el = 0.0_dp

       DO ispin=1,nspins
          ! Copy the core matrix into the fock matrix
          CALL transfer_matrix(matrix_h(1)%matrix,ks_matrix(ispin)%matrix,error=error)
       END DO

       ! Zero the ks_scp matrix.  Will get info from coulomb routines
       IF (scp_nddo) THEN
          CALL get_scp_env(scp_env, ks_scp=ks_scp, error=error)
          CALL set_matrix(ks_scp,0.0_dp)
       END IF

       ! Compute Exchange and Coulomb terms
       CALL semi_empirical_si_initialize(store_int_env, ks_env%s_mstruct_changed, error)
       CALL build_fock_matrix_exchange(qs_env,ks_matrix,matrix_p,calculate_forces,&
            store_int_env,error)
       CALL build_fock_matrix_coulomb(qs_env,ks_matrix,matrix_p,energy,calculate_forces,&
            store_int_env,error)
       CALL semi_empirical_si_finalize(store_int_env, ks_env%s_mstruct_changed, error)

       ! Long Range Electrostatic
       IF (se_control%do_ewald) THEN
          ! Debug statements
          IF (debug_coulomb_lr) THEN
             ! Zero structures only for debugging purpose
             energy%hartree = 0.0_dp
             DO ispin = 1, nspins
                CALL set_matrix(ks_matrix(ispin)%matrix,0.0_dp)
             END DO
          END IF

          ! Evaluate Coulomb Long-Range
          CALL build_fock_matrix_coulomb_lr(qs_env,ks_matrix,matrix_p,energy,calculate_forces,error)

          ! Possibly handle the slowly convergent term 1/R^3
          IF (se_control%do_ewald_r3) THEN
             CALL build_fock_matrix_coul_lr_r3(qs_env,ks_matrix,matrix_p,energy,calculate_forces,error)
          END IF
       END IF

       ! Compute the Hartree energy
       ! NOTE: If we are performing SCP-NDDO, ks_matrix contains coulomb piece from SCP.
       DO ispin=1,nspins
          CALL cp_sm_sm_trace(ks_matrix(ispin)%matrix,matrix_p(ispin)%matrix,trace=ecoul,&
               para_env=para_env,error=error)
          energy%hartree = energy%hartree + ecoul

          ! Debug statements
          IF (debug_coulomb_lr) THEN
             WRITE(*,*)ispin,"ECOUL ",ecoul
          END IF
       END DO

       ! Compute the SCP Hartree energy and polarization energy and update ks_scp
       IF (scp_nddo) THEN
          CALL get_scp_env ( scp_env, pscp=pscp, ks_scp=ks_scp, energy=energy_scp, &
                             ascp_diag= ascp_diag, error = error )
          ! Get energy
          CALL cp_sm_sm_trace(ks_scp,pscp,trace=ecoul_scp,error=error)
          energy_scp%e_scp_hartree = ecoul_scp*0.5_dp
          ! evalulate kernel and forces on the SCP coeffs
          CALL evaluate_scp_force ( atomic_kind_set, pscp, ks_scp, ascp_diag, para_env,&
                                    just_energy, energy_scp%e_scp_kernel, error )
          ! NOTE: If needed, at this point, forces on P_SCP are up-to-date and stored 
          !       in KS_SCP as a*ks_scp
       END IF

       ! Debug statements
       IF (debug_coulomb_lr) THEN
          WRITE(*,*)"TOTAL ENE",0.5_dp*energy%hartree
          STOP
       END IF

       ! QM/MM
       IF (qs_env%qmmm) THEN
          DO ispin=1,nspins
             ! If QM/MM sumup the 1el Hamiltonian
             CALL cp_sm_scale_and_add(matrix_a=ks_matrix(ispin)%matrix,&
                                      matrix_b=qs_env%ks_qmmm_env%matrix_h(1)%matrix,&
                                      error=error)
             ! Compute QM/MM Energy
             CALL cp_sm_sm_trace(qs_env%ks_qmmm_env%matrix_h(1)%matrix,&
                                 matrix_p(ispin)%matrix,trace=qmmm_el,&
                                 para_env=para_env,error=error)
             energy%qmmm_el = energy%qmmm_el + qmmm_el
          END DO
       END IF

       ! Collect all the energy terms
       energy%mulliken = 0.0_dp
       energy%exc      = 0.0_dp
       energy%total    = energy%core +&
                         energy%core_overlap +&
                         0.5_dp*energy%hartree + &
                         energy%qmmm_el +&
                         energy%mulliken

       ! Modify energy for SCP
       IF (scp_nddo) THEN 
          energy%total = energy%total + energy_scp%e_scp_hartree + energy_scp%e_scp_kernel
       END IF

       output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%DETAILED_ENERGY",&
            extension=".scfLog",error=error)
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
               "Core Hamiltonian energy:                       ",energy%core,&
               "Two-electron integral energy:                  ",energy%hartree
          IF (qs_env%qmmm) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "QM/MM Electrostatic energy:                    ",energy%qmmm_el
          END IF
          IF (scp_nddo) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "SCP electrostatic energy:                    ",energy_scp%e_scp_hartree
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "SCP kernel energy:                    ",energy_scp%e_scp_kernel
          END IF
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%DETAILED_ENERGY", error=error)

       ! Here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C (plus occupation numbers)
       IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
          CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array,error=error)
          DO ispin=1,SIZE(mo_derivs)
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff, &
                             occupation_numbers=occupation_numbers )
             CALL cp_fm_get_info(mo_coeff,ncol_global=ncol_global,error=error)
             CALL cp_sm_fm_multiply(ks_matrix(ispin)%matrix,mo_coeff,&
                                    mo_derivs(ispin)%matrix, ncol_global,error=error)
          ENDDO
       ENDIF
    END IF
    CALL timestop(handle)

  END SUBROUTINE build_nddo_fock_matrix

! *****************************************************************************
!> \brief Evaluates the core-core interactions for NDDO methods
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE se_core_core_interaction(qs_env, para_env, calculate_forces, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'se_core_core_interaction', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iab, iatom, ikind, ilist, inode, &
      itype, jatom, jkind, natom, nkind, nlist, nnode, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: anag, defined, failure
    REAL(KIND=dp)                            :: delta, dr1, dr3inv(3), enuc, &
                                                enuclear, r2inv, r3inv, rinv
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_2c
    TYPE(neighbor_list_type), POINTER        :: sab_2c_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_2c_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(se_int_control_type)                :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b

    failure  = .FALSE.
    enuclear = 0.0_dp
    NULLIFY(dft_control,cell,force,particle_set,se_control,se_taper,atomic_kind_set)

    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, se_taper=se_taper,&
         energy=energy, error=error)

    CALL initialize_se_taper(se_taper,coulomb=.TRUE.,error=error)
    ! Parameters
    se_control => dft_control%qs_control%se_control
    anag       =  se_control%analytical_gradients
    CALL setup_se_int_control_type(se_int_control, shortrange=se_control%do_ewald,&
         do_ewald_r3=se_control%do_ewald_r3,integral_screening=se_control%integral_screening,&
         max_multipole=se_control%max_multipole)
    IF(se_control%do_ewald_gks) THEN
      se_int_control%do_ewald_gks = .TRUE.
      se_int_control%shortrange = .TRUE.
      CALL get_qs_env(qs_env=qs_env,ewald_env=ewald_env,ewald_pw=ewald_pw,error=error)
      CALL ewald_env_get (ewald_env, alpha=se_int_control%ewald_gks%alpha, error=error)
      CALL ewald_pw_get (ewald_pw, pw_big_pool=se_int_control%ewald_gks%pw_pool, &
                         dg=se_int_control%ewald_gks%dg)
    END IF

    IF ( .NOT. failure ) THEN
       CALL get_qs_env(qs_env=qs_env,sab_2c=sab_2c,atomic_kind_set=atomic_kind_set,&
            error=error)

       nkind = SIZE(atomic_kind_set)
       ! Possibly compute forces
       IF(calculate_forces) THEN
          CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,force=force,error=error)
          natom = SIZE (particle_set)
          ALLOCATE (atom_of_kind(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          delta = se_control%delta
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)
       END IF

       itype    = get_se_type(dft_control%qs_control%method_id)

       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,defined=defined)

          IF (.NOT.defined) CYCLE

          DO jkind=1,nkind
             atomic_kind => atomic_kind_set(jkind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_b)
             CALL get_se_param(se_kind_b,defined=defined)

             IF (.NOT.defined) CYCLE

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
                CALL get_neighbor_list(neighbor_list=sab_2c_neighbor_list,atom=iatom,nnode=nnode)
                sab_2c_neighbor_node => first_node(sab_2c_neighbor_list)
                
                DO inode=1,nnode
                   CALL get_neighbor_node(neighbor_node=sab_2c_neighbor_node,neighbor=jatom,r=rij)
                   dr1 = DOT_PRODUCT(rij,rij)
                   IF ( dr1 > rij_threshold ) THEN
                      SELECT CASE (dft_control%qs_control%method_id)
                      CASE (do_method_mndo,do_method_am1,do_method_pm3,do_method_pm6,do_method_pdg,&
                           do_method_rm1,do_method_mndod, do_method_pnnl)

                         ! Core-Core energy term
                         CALL corecore (se_kind_a,se_kind_b,rij,enuc=enuc,itype=itype,anag=anag,&
                              se_int_control=se_int_control, se_taper=se_taper, error=error)
                         enuclear = enuclear + enuc
                         ! Residual integral (1/R^3) correction
                         IF (se_int_control%do_ewald_r3) THEN
                            r2inv = 1.0_dp/dr1
                            rinv  = SQRT(r2inv)
                            r3inv = rinv**3
                            
                            ! Core-Core term
                            enuclear = enuclear + se_kind_a%expns3_int(jkind)%expns3%core_core*r3inv
                         END IF
                         
                         ! Core-Core Derivatives
                         IF(calculate_forces) THEN
                            atom_a = atom_of_kind(iatom)
                            atom_b = atom_of_kind(jatom)
                            
                            CALL dcorecore (se_kind_a,se_kind_b,rij,denuc=force_ab,itype=itype,delta=delta,&
                                 anag=anag,se_int_control=se_int_control,se_taper=se_taper,error=error)

                            ! Residual integral (1/R^3) correction
                            IF (se_int_control%do_ewald_r3) THEN
                               dr3inv   = -3.0_dp*rij*r3inv*r2inv
                               ! Derivatives of core-core terms
                               force_ab = force_ab + se_kind_a%expns3_int(jkind)%expns3%core_core*dr3inv
                            END IF

                            ! Sum up force components
                            force(ikind)%all_potential(1,atom_a) = force(ikind)%all_potential(1,atom_a) - force_ab(1)
                            force(jkind)%all_potential(1,atom_b) = force(jkind)%all_potential(1,atom_b) + force_ab(1)

                            force(ikind)%all_potential(2,atom_a) = force(ikind)%all_potential(2,atom_a) - force_ab(2)
                            force(jkind)%all_potential(2,atom_b) = force(jkind)%all_potential(2,atom_b) + force_ab(2)

                            force(ikind)%all_potential(3,atom_a) = force(ikind)%all_potential(3,atom_a) - force_ab(3)
                            force(jkind)%all_potential(3,atom_b) = force(jkind)%all_potential(3,atom_b) + force_ab(3)
                         END IF
                      CASE DEFAULT
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                      END SELECT
                   END IF
                   sab_2c_neighbor_node => next(sab_2c_neighbor_node)
                END DO ! inode => jatom(atom B)               
             END DO ! ilist => iatom(atom A)
          END DO ! jkind
       END DO ! ikind
       
       IF (calculate_forces) THEN 
          DEALLOCATE(atom_of_kind,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       
       CALL mp_sum(enuclear,para_env%group)
       energy%core_overlap  = enuclear
       energy%core_overlap0 = enuclear
    END IF

    CALL finalize_se_taper(se_taper,error=error)
    CALL timestop(handle)
  END SUBROUTINE se_core_core_interaction
  
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
                                                failure, scp_nddo, switch
    REAL(KIND=dp)                            :: delta, dr1, ecore2
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

    failure=.FALSE.
    CALL timeset(routineN,handle)
    IF ( .NOT. failure ) THEN
       NULLIFY(dft_control,cell,force,particle_set,diagmat_ks,diagmat_p,&
               se_control,se_taper,pscp,ks_scp,scp_env)
       
       CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, se_taper=se_taper,&
            para_env=para_env, sab_2c=sab_2c, atomic_kind_set=atomic_kind_set, &
            particle_set=particle_set, error=error)
      

       ! Parameters
       CALL initialize_se_taper(se_taper,coulomb=.TRUE.,error=error)
       se_control  => dft_control%qs_control%se_control
       anag        =  se_control%analytical_gradients
       scp_nddo    =  se_control%scp

       CALL setup_se_int_control_type(se_int_control, shortrange=se_control%do_ewald,&
            do_ewald_r3=se_control%do_ewald_r3, integral_screening=se_control%integral_screening,&
            max_multipole=se_control%max_multipole)
       IF(se_control%do_ewald_gks) THEN
         se_int_control%do_ewald_gks = .TRUE.
         se_int_control%shortrange = .TRUE.
         CALL get_qs_env(qs_env=qs_env,ewald_env=ewald_env,ewald_pw=ewald_pw,error=error)
         CALL ewald_env_get (ewald_env, alpha=se_int_control%ewald_gks%alpha, error=error)
         CALL ewald_pw_get (ewald_pw, pw_big_pool=se_int_control%ewald_gks%pw_pool, &
                            dg=se_int_control%ewald_gks%dg)
       END IF

       IF (scp_nddo) CALL get_qs_env(qs_env=qs_env, scp_env=scp_env, error=error)
       
       nspins=dft_control%nspins
       CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
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

       IF (scp_nddo)  CALL get_scp_env(scp_env=scp_env, pscp=pscp, ks_scp=ks_scp, error=error)
       
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
                IF ( nspins >= 2 ) THEN
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
                   IF ( nspins >= 2 ) &
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
                      IF ( nspins >= 2 ) THEN
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
                         IF ( nspins >= 2 ) &
                              pb_b(1:natorb_b2) = RESHAPE(pb_block_b+pscpb_block,(/natorb_b2/))
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
 
                            ! Derivatives of the One-centre One-electron terms
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
       calculate_forces, error)
   
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_p
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_fock_matrix_coulomb_lr', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, ewald_type, forces_g_size, handle, iatom, ikind, &
      ilist, indi, indj, ispin, itype, iw, jint, natoms, natorb_a, nkind, &
      nlocal_particles, node, nparticle_local, nspins, size_1c_int, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: anag, defined, failure, &
                                                scp_nddo
    LOGICAL, DIMENSION(3)                    :: task
    REAL(KIND=dp)                            :: e_neut, e_self, energy_glob, &
                                                energy_local, enuc, fac, tmp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: forces_g, forces_r
    REAL(KIND=dp), DIMENSION(3)              :: force_a
    REAL(KIND=dp), DIMENSION(3, 3)           :: qcart
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

    failure=.FALSE.
    CALL timeset(routineN,handle)
    IF ( .NOT. failure ) THEN
       NULLIFY(dft_control, cell, force,particle_set, diagmat_ks, diagmat_p,local_particles,&
            se_control, ewald_env, ewald_pw, se_nddo_mpole, se_nonbond_env, se_section, mpole,&
            logger,pscp,ks_scp,scp_env)

       logger => cp_error_get_logger(error)
       CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, para_env=para_env,&
            atomic_kind_set=atomic_kind_set, particle_set=particle_set, ewald_env=ewald_env,&
            local_particles=local_particles, ewald_pw=ewald_pw, se_nddo_mpole=se_nddo_mpole,&
            se_nonbond_env=se_nonbond_env, scp_env=scp_env, error=error)
       
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
       
       nspins=dft_control%nspins
       CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
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
       DO ispin = 1, nspins
          se_nddo_mpole%charge     = 0.0_dp
          se_nddo_mpole%dipole     = 0.0_dp
          se_nddo_mpole%quadrupole = 0.0_dp

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
                IF (task(1)) se_nddo_mpole%charge(iatom)   = se_kind_a%zeff
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
             
             CALL ewald_multipole_evaluate(ewald_env, ewald_pw, se_nonbond_env, cell,&
                  particle_set, local_particles, energy_local, energy_glob, e_neut, e_self, task,&
                  do_forces=.TRUE., do_stress=.FALSE., do_efield=.TRUE., charges=se_nddo_mpole%charge,&
                  dipoles=se_nddo_mpole%dipole, quadrupoles=se_nddo_mpole%quadrupole,&
                  forces_local=forces_g, forces_glob=forces_r, efield0=se_nddo_mpole%efield0,&
                  efield1=se_nddo_mpole%efield1, efield2=se_nddo_mpole%efield2, iw=iw, &
                  do_debug=.TRUE.,error=error)
             ! Only SR force have to be summed up.. the one in g-space are already fully local..
             CALL mp_sum(forces_r, para_env%group)
          ELSE
             CALL ewald_multipole_evaluate(ewald_env, ewald_pw, se_nonbond_env, cell,&
                  particle_set, local_particles, energy_local, energy_glob, e_neut, e_self, task,&
                  do_forces=.FALSE., do_stress=.FALSE., do_efield=.TRUE., charges=se_nddo_mpole%charge,&
                  dipoles=se_nddo_mpole%dipole, quadrupoles=se_nddo_mpole%quadrupole,&
                  efield0=se_nddo_mpole%efield0, efield1=se_nddo_mpole%efield1,&
                  efield2=se_nddo_mpole%efield2, iw=iw, do_debug=.TRUE.,error=error)
          END IF
          CALL cp_print_key_finished_output(iw,logger,se_section,"PRINT%EWALD_INFO",error=error)

          ! Debug Statements
          IF (debug_coulomb_lr) THEN
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
                IF ( nspins >= 1 ) THEN
                   CALL get_block_node(matrix=diagmat_ks(ispin)%matrix,&
                                       block_row=iatom,&
                                       block_col=iatom,&
                                       BLOCK=ksa_block_a)
                   CPPostcondition(ASSOCIATED(ksa_block_a),cp_failure_level,routineP,error,failure)
                END IF

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
       END DO
       ! Sum nuclear energy contribution
       CALL mp_sum(enuc, para_env%group)
       energy%core_overlap = energy%core_overlap0 + 0.5_dp*enuc

       ! Debug Statements
       IF (debug_coulomb_lr) THEN
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
       CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
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
                IF ( nspins >= 2 ) THEN
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
                      IF ( nspins >= 2 ) THEN
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

! *****************************************************************************
!> \brief Construction of the Exchange part of the Fock matrix
!> \author JGH
! *****************************************************************************
  SUBROUTINE build_fock_matrix_exchange (qs_env, ks_matrix, matrix_p, calculate_forces,&
       store_int_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_p
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_fock_matrix_exchange', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iab, iatom, icol, ikind, ilist, inode, &
      integral_screening, irow, jatom, jkind, natom, natorb_a, natorb_b, &
      nkind, nlist, nnode, nspins, stat
    INTEGER, DIMENSION(2)                    :: size_p_block_a
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: anag, check, defined, &
                                                failure, switch
    REAL(KIND=dp)                            :: delta, dr
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    REAL(KIND=dp), DIMENSION(45, 45)         :: p_block_tot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ks_block_a, ks_block_b, &
                                                p_block_a, p_block_b
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(se_int_control_type)                :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(dft_control,cell,force,particle_set,se_control,se_taper)
    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,cell=cell,se_taper=se_taper,&
         para_env=para_env,error=error)

    CALL initialize_se_taper(se_taper,exchange=.TRUE.,error=error)
    se_control => dft_control%qs_control%se_control
    anag       =  se_control%analytical_gradients
    nspins=dft_control%nspins

    CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix_p),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)

    IF ( .NOT. failure ) THEN
       ! Identify proper integral screening (according user requests)
       integral_screening = se_control%integral_screening
       IF ((integral_screening==do_se_IS_kdso_d).AND.(.NOT.se_control%force_kdsod_EX)) THEN
          integral_screening = do_se_IS_kdso
       END IF
       CALL setup_se_int_control_type(se_int_control, shortrange=.FALSE., do_ewald_r3=.FALSE.,&
            integral_screening=integral_screening, max_multipole=do_multipole_none)

       CALL get_qs_env(qs_env=qs_env,sab_orb=sab_orb,atomic_kind_set=atomic_kind_set,error=error)

       nkind = SIZE(atomic_kind_set)
       IF(calculate_forces) THEN
          CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,force=force,error=error)
          natom = SIZE (particle_set)
          ALLOCATE (atom_of_kind(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          delta = se_control%delta
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)
       END IF

       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a)

          IF (.NOT.defined .OR. natorb_a < 1) CYCLE

          DO jkind=1,nkind
             atomic_kind => atomic_kind_set(jkind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_b)
             CALL get_se_param(se_kind_b,defined=defined,natorb=natorb_b)

             IF (.NOT.defined .OR. natorb_b < 1) CYCLE

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
                CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,atom=iatom,nnode=nnode)
                sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

                DO inode=1,nnode
                   CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,neighbor=jatom,r=rij)
                   IF (iatom <= jatom) THEN
                      irow = iatom
                      icol = jatom
                      switch = .FALSE.
                   ELSE
                      irow = jatom
                      icol = iatom
                      switch = .TRUE.
                   END IF
                   ! Retrieve blocks for KS and P
                   CALL get_block_node(matrix=ks_matrix(1)%matrix,&
                                       block_row=irow,&
                                       block_col=icol,&
                                       BLOCK=ks_block_a)
                   CPPostcondition(ASSOCIATED(ks_block_a),cp_failure_level,routineP,error,failure)
                   CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                       block_row=irow,&
                                       block_col=icol,&
                                       BLOCK=p_block_a)
                   CPPostcondition(ASSOCIATED(p_block_a),cp_failure_level,routineP,error,failure)
                   size_p_block_a(1) = SIZE(p_block_a,1)
                   size_p_block_a(2) = SIZE(p_block_a,2)
                   p_block_tot(1:size_p_block_a(1),1:size_p_block_a(2)) = 2.0_dp * p_block_a
                   
                   ! Handle more configurations
                   IF ( nspins >= 2 ) THEN
                      CALL get_block_node(matrix=ks_matrix(2)%matrix,&
                                          block_row=irow,&
                                          block_col=icol,&
                                          BLOCK=ks_block_b)
                      CPPostcondition(ASSOCIATED(ks_block_b),cp_failure_level,routineP,error,failure)
                      CALL get_block_node(matrix=matrix_p(2)%matrix,&
                                          block_row=irow,&
                                          block_col=icol,&
                                          BLOCK=p_block_b)
                      CPPostcondition(ASSOCIATED(p_block_b),cp_failure_level,routineP,error,failure)
                      check = (size_p_block_a(1)==SIZE(p_block_b,1)).AND.(size_p_block_a(2)==SIZE(p_block_b,2))
                      CPPostcondition(check,cp_failure_level,routineP,error,failure)
                      p_block_tot(1:SIZE(p_block_a,1),1:SIZE(p_block_a,2)) = p_block_a + p_block_b
                   END IF

                   dr = DOT_PRODUCT(rij,rij)
                   IF ( iatom == jatom .AND. dr < rij_threshold ) THEN
                      ! Once center - Two electron Terms
                      IF      ( nspins == 1 ) THEN
                         CALL fock1_2el(se_kind_a,p_block_tot,p_block_a,ks_block_a,factor=0.5_dp,error=error)
                      ELSE IF ( nspins == 2 ) THEN
                         CALL fock1_2el(se_kind_a,p_block_tot,p_block_a,ks_block_a,factor=1.0_dp,error=error)
                         CALL fock1_2el(se_kind_a,p_block_tot,p_block_b,ks_block_b,factor=1.0_dp,error=error)
                      END IF
                   ELSE
                      ! Exchange Terms
                      IF      ( nspins == 1 ) THEN
                         CALL fock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_tot, p_block_a, ks_block_a,&
                              factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, &
                              store_int_env=store_int_env, error=error)
                      ELSE IF ( nspins == 2 ) THEN
                         CALL fock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_tot, p_block_a, ks_block_a,&
                              factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, &
                              store_int_env=store_int_env, error=error)
                         
                         CALL fock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_tot, p_block_b, ks_block_b,&
                              factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, &
                              store_int_env=store_int_env, error=error)
                      END IF
                      IF(calculate_forces) THEN
                         atom_a = atom_of_kind(iatom)
                         atom_b = atom_of_kind(jatom)
                         force_ab = 0.0_dp
                         IF      ( nspins == 1 ) THEN
                            CALL dfock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_tot, p_block_a,&
                                 factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab,&
                                 delta=delta, error=error)
                         ELSE IF ( nspins == 2 ) THEN
                            CALL dfock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_tot, p_block_a,&
                                 factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab,&
                                 delta=delta, error=error)
                            
                            CALL dfock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_tot, p_block_b,&
                                 factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab,&
                                 delta=delta, error=error)
                         END IF
                         IF (switch) THEN
                            force_ab(1) = -force_ab(1)
                            force_ab(2) = -force_ab(2)
                            force_ab(3) = -force_ab(3)
                         END IF
                         force(ikind)%rho_elec(1,atom_a) = force(ikind)%rho_elec(1,atom_a) - force_ab(1)
                         force(jkind)%rho_elec(1,atom_b) = force(jkind)%rho_elec(1,atom_b) + force_ab(1)

                         force(ikind)%rho_elec(2,atom_a) = force(ikind)%rho_elec(2,atom_a) - force_ab(2)
                         force(jkind)%rho_elec(2,atom_b) = force(jkind)%rho_elec(2,atom_b) + force_ab(2)

                         force(ikind)%rho_elec(3,atom_a) = force(ikind)%rho_elec(3,atom_a) - force_ab(3)
                         force(jkind)%rho_elec(3,atom_b) = force(jkind)%rho_elec(3,atom_b) + force_ab(3)
                      END IF
                   END IF
                   sab_orb_neighbor_node => next(sab_orb_neighbor_node)
                END DO ! inode => jatom(atom B)
             END DO ! ilist => iatom(atom A)             
          END DO ! jkind
       END DO ! ikind
       
       IF (calculate_forces) THEN 
          DEALLOCATE(atom_of_kind,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    CALL finalize_se_taper(se_taper,error=error)
    CALL timestop(handle)
    
  END SUBROUTINE build_fock_matrix_exchange

! *****************************************************************************
!> \brief Construction of 2-center 1-electron Fock Matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock2_1el (sepi, sepj, rij, ksi_block, ksj_block, pi_block, pj_block,&
       ecore, itype, anag, se_int_control, se_taper, store_int_env, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(INOUT)                          :: ksi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(INOUT)                          :: ksj_block
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(IN)                             :: pi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(IN)                             :: pj_block
    REAL(KIND=dp), INTENT(INOUT)             :: ecore
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2_1el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2, j1, j1L
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(45)             :: e1b, e2a

    failure = .FALSE.
    ! Compute integrals
    CALL rotnuc (sepi, sepj, rij, e1b=e1b, e2a=e2a, itype=itype, anag=anag,&
         se_int_control=se_int_control, se_taper=se_taper, store_int_env=store_int_env,&
         error=error)
    !
    ! Add the electron-nuclear attraction term for atom sepi
    !
    i2 = 0
    DO i1L = 1, sepi%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          ksi_block(i1,j1) = ksi_block(i1,j1) + e1b(i2)
          ksi_block(j1,i1) = ksi_block(i1,j1)
          ecore = ecore + 2.0_dp * e1b(i2) * pi_block(i1,j1)
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       ksi_block(i1,j1) = ksi_block(i1,j1) + e1b(i2)
       ecore = ecore + e1b(i2) * pi_block(i1,j1)
    END DO
    !
    ! Add the electron-nuclear attraction term for atom sepj
    !
    i2 = 0
    DO i1L = 1, sepj%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          ksj_block(i1,j1) = ksj_block(i1,j1) + e2a(i2)
          ksj_block(j1,i1) = ksj_block(i1,j1)
          ecore = ecore + 2.0_dp * e2a(i2) * pj_block(i1,j1)
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       ksj_block(i1,j1) = ksj_block(i1,j1) + e2a(i2)
       ecore = ecore + e2a(i2) * pj_block(i1,j1)
    END DO

  END SUBROUTINE fock2_1el

! *****************************************************************************
!> \brief Derivatives of 2-center 1-electron Fock Matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2_1el (sepi, sepj, rij, pi_block, pj_block, itype, anag,&
       se_int_control, se_taper, force, delta, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(IN)                             :: pi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(IN)                             :: pj_block
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), INTENT(IN)                :: delta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2_1el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2, j1, j1L
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tmp
    REAL(KIND=dp), DIMENSION(3, 45)          :: de1b, de2a

    failure = .FALSE.
    ! Compute integrals
    CALL drotnuc (sepi, sepj, rij, de1b=de1b, de2a=de2a, itype=itype, anag=anag,&
         se_int_control=se_int_control, se_taper=se_taper, delta=delta, error=error)
    !
    ! Add the electron-nuclear attraction term for atom sepi
    !
    i2 = 0
    DO i1L = 1, sepi%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          tmp= 2.0_dp * pi_block(i1,j1)
          force(1) = force(1) + de1b(1,i2) * tmp 
          force(2) = force(2) + de1b(2,i2) * tmp
          force(3) = force(3) + de1b(3,i2) * tmp
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       force(1) = force(1) + de1b(1,i2) * pi_block(i1,j1)
       force(2) = force(2) + de1b(2,i2) * pi_block(i1,j1)
       force(3) = force(3) + de1b(3,i2) * pi_block(i1,j1)
    END DO
    !
    ! Add the electron-nuclear attraction term for atom sepj
    !
    i2 = 0
    DO i1L = 1, sepj%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          tmp= 2.0_dp * pj_block(i1,j1)
          force(1) = force(1) + de2a(1,i2) * tmp
          force(2) = force(2) + de2a(2,i2) * tmp
          force(3) = force(3) + de2a(3,i2) * tmp
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       force(1) = force(1) + de2a(1,i2) * pj_block(i1,j1)
       force(2) = force(2) + de2a(2,i2) * pj_block(i1,j1)
       force(3) = force(3) + de2a(3,i2) * pj_block(i1,j1)
    END DO
    
  END SUBROUTINE dfock2_1el

! *****************************************************************************
!> \brief Construction of 1-center 2-electron Fock Matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock1_2el(sep, p_tot, p_mat, f_mat, factor, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: p_tot
    REAL(KIND=dp), &
      DIMENSION(sep%natorb, sep%natorb), &
      INTENT(IN)                             :: p_mat
    REAL(KIND=dp), &
      DIMENSION(sep%natorb, sep%natorb), &
      INTENT(INOUT)                          :: f_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock1_2el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ijw, ikw, iL, im, j, jL, &
                                                jlw, jm, k, kL, klw, l, lL
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: sum

    failure = .FALSE.
    !   One-center coulomb and exchange terms for semiempirical_type sep
    !
    !  F(i,j)=F(i,j)+sum(k,l)((PA(k,l)+PB(k,l))*<i,j|k,l>
    !                        -(PA(k,l)        )*<i,k|j,l>), k,l on type sep.
    !
    DO iL = 1, sep%natorb
       i = se_orbital_pointer(iL)
       DO jL = 1, iL
          j = se_orbital_pointer(jL)

          !    `J' Address IJ in W
          ijw = (iL*(iL-1))/2 + jL
          sum = 0.0_dp
          DO kL = 1, sep%natorb
             k = se_orbital_pointer(kL)
             DO lL =  1, sep%natorb
                l = se_orbital_pointer(lL)

                !    `J' Address KL in W
                im = MAX(kL,lL)
                jm = MIN(kL,lL)
                klw = (im*(im-1))/2 + jm

                !    `K' Address IK in W
                im = MAX(kL,jL)
                jm = MIN(kL,jL)
                ikw = (im*(im-1))/2 + jm

                !    `K' Address JL in W
                im = MAX(lL,iL)
                jm = MIN(lL,iL)
                jlw = (im*(im-1))/2 + jm

                sum = sum + p_tot(k,l) * sep%w(ijw, klw) - p_mat(k,l) * sep%w(ikw, jlw)
             END DO
          END DO
          f_mat(i,j) = f_mat(i,j) + factor*sum
          f_mat(j,i) = f_mat(i,j)
       END DO
    END DO
  END SUBROUTINE fock1_2el

! *****************************************************************************
!> \brief Construction of 2-center Fock Matrix - Coulomb Terms
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock2C(sepi, sepj, rij, switch, pi_tot, fi_mat, pj_tot, fj_mat, &
       factor, anag, se_int_control, se_taper, store_int_env, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(INOUT)                          :: fi_mat
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pj_tot
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(INOUT)                          :: fj_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2C', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, aa, bb, irij(3)
    REAL(KIND=dp), DIMENSION(2025)           :: w

    failure = .FALSE.
    ! Evaluate integrals
    IF (.NOT.switch) THEN
       CALL rotint (sepi,sepj, rij,w,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,store_int_env=store_int_env, error=error)
    ELSE
       irij = -rij
       CALL rotint (sepj,sepi,irij,w,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,store_int_env=store_int_env, error=error)
    END IF
    kr = 0
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
    END IF
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb(2)
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr = kr + 1
                a = w(kr)*factor
                ! Coulomb
                IF (.NOT.switch) THEN
                   fi_mat(i,j) = fi_mat(i,j) + bb * a * pj_tot(k,l)
                   fj_mat(k,l) = fj_mat(k,l) + aa * a * pi_tot(i,j)
                   fi_mat(j,i) = fi_mat(i,j)
                   fj_mat(l,k) = fj_mat(k,l)
                ELSE
                   fj_mat(i,j) = fj_mat(i,j) + bb * a * pi_tot(k,l)
                   fi_mat(k,l) = fi_mat(k,l) + aa * a * pj_tot(i,j)
                   fj_mat(j,i) = fj_mat(i,j)
                   fi_mat(l,k) = fi_mat(k,l)
                END IF
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE fock2C

! *****************************************************************************
!> \brief Derivatives of 2-center Fock Matrix - Coulomb Terms
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2C(sepi, sepj, rij, switch, pi_tot, pj_tot, factor, anag,&
       se_int_control, se_taper, force, delta, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot, pj_tot
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), INTENT(IN)                :: delta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2C', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aa, bb, tmp
    REAL(KIND=dp), DIMENSION(3)              :: a, irij
    REAL(KIND=dp), DIMENSION(3, 2025)        :: dw

    failure = .FALSE.
    ! Evaluate integrals' derivatives
    IF (.NOT.switch) THEN
       CALL drotint (sepi,sepj, rij,dw,delta,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,error=error)
    ELSE
       irij = -rij
       CALL drotint (sepj,sepi,irij,dw,delta,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,error=error)
    END IF

    kr = 0
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
    END IF
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb(2)
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr = kr + 1
                a(1) = dw(1,kr)*factor
                a(2) = dw(2,kr)*factor
                a(3) = dw(3,kr)*factor
                ! Coulomb
                IF (.NOT.switch) THEN
                   tmp = bb * aa * pj_tot(k,l) * pi_tot(i,j)
                ELSE
                   tmp = bb * aa * pi_tot(k,l) * pj_tot(i,j)
                END IF
                force(1) = force(1) + a(1) * tmp
                force(2) = force(2) + a(2) * tmp
                force(3) = force(3) + a(3) * tmp
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE dfock2C

! *****************************************************************************
!> \brief Construction of 2-center Fock Matrix - General Driver
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock2E(sepi, sepj, rij, switch, isize, pi_tot, pi_mat, fi_mat, factor,&
       anag, se_int_control, se_taper, store_int_env, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    INTEGER, DIMENSION(2), INTENT(IN)        :: isize
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot
    REAL(KIND=dp), &
      DIMENSION(isize(1), isize(2)), &
      INTENT(IN)                             :: pi_mat
    REAL(KIND=dp), &
      DIMENSION(isize(1), isize(2)), &
      INTENT(INOUT)                          :: fi_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2E', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, aa, bb, irij(3)
    REAL(KIND=dp), DIMENSION(2025)           :: w

    failure = .FALSE.
    ! Evaluate integrals
    IF (.NOT.switch) THEN
       CALL rotint (sepi,sepj, rij,w,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,store_int_env=store_int_env, error=error)
    ELSE
       irij = -rij
       CALL rotint (sepj,sepi,irij,w,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,store_int_env=store_int_env, error=error)
    END IF
    kr = 0
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
    END IF
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb(2)
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr = kr + 1
                a = w(kr)*factor
                ! Exchange
                a = a * aa * bb * 0.25_dp
                fi_mat(i,k) = fi_mat(i,k) - a * pi_mat(j,l)
                fi_mat(i,l) = fi_mat(i,l) - a * pi_mat(j,k)
                fi_mat(j,k) = fi_mat(j,k) - a * pi_mat(i,l)
                fi_mat(j,l) = fi_mat(j,l) - a * pi_mat(i,k)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE fock2E

! *****************************************************************************
!> \brief Derivatives of 2-center Fock Matrix - General Driver
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2E(sepi, sepj, rij, switch, isize, pi_tot, pi_mat, factor, anag,&
       se_int_control, se_taper, force, delta, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    INTEGER, DIMENSION(2), INTENT(IN)        :: isize
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot
    REAL(KIND=dp), &
      DIMENSION(isize(1), isize(2)), &
      INTENT(IN)                             :: pi_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), INTENT(IN)                :: delta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2E', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aa, bb, tmp, tmp1, tmp2, &
                                                tmp3, tmp4
    REAL(KIND=dp), DIMENSION(3)              :: a, irij
    REAL(KIND=dp), DIMENSION(3, 2025)        :: dw

    failure = .FALSE.
    ! Evaluate integrals' derivatives
    IF (.NOT.switch) THEN
       CALL drotint (sepi,sepj, rij,dw,delta,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,error=error)
    ELSE
       irij = -rij
       CALL drotint (sepj,sepi,irij,dw,delta,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,error=error)
    END IF

    kr = 0
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
    END IF
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb(2)
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr   = kr + 1
                tmp  = factor * aa * bb * 0.25_dp
                a(1) = dw(1,kr)*tmp
                a(2) = dw(2,kr)*tmp
                a(3) = dw(3,kr)*tmp
                ! Exchange
                tmp1 = pi_mat(j,l) * pi_mat(i,k)
                tmp2 = pi_mat(j,k) * pi_mat(i,l)
                tmp3 = pi_mat(i,l) * pi_mat(j,k)
                tmp4 = pi_mat(i,k) * pi_mat(j,l)

                force(1) = force(1) - a(1) * tmp1
                force(1) = force(1) - a(1) * tmp2
                force(1) = force(1) - a(1) * tmp3
                force(1) = force(1) - a(1) * tmp4
                
                force(2) = force(2) - a(2) * tmp1
                force(2) = force(2) - a(2) * tmp2
                force(2) = force(2) - a(2) * tmp3
                force(2) = force(2) - a(2) * tmp4
                
                force(3) = force(3) - a(3) * tmp1
                force(3) = force(3) - a(3) * tmp2
                force(3) = force(3) - a(3) * tmp3
                force(3) = force(3) - a(3) * tmp4
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE dfock2E

! *****************************************************************************
!> \brief  Construction of 2-center 1-electron Fock Matrix for the residual 
!>         (1/R^3) integral part
!> \author Teodoro Laino [tlaino]
!> \date   12.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock2_1el_r3 (sepi, sepj, ksi_block, ksj_block, pi_block, pj_block,&
       e1b, e2a, ecore, rp, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(INOUT)                          :: ksi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(INOUT)                          :: ksj_block
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(IN)                             :: pi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(IN)                             :: pj_block
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: e1b, e2a
    REAL(KIND=dp), INTENT(INOUT)             :: ecore
    REAL(KIND=dp), INTENT(IN)                :: rp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2_1el_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2
    LOGICAL                                  :: failure

    failure = .FALSE.
    !
    ! Add the electron-nuclear residual attraction term for atom sepi
    !
    i2 = 0
    DO i1L = 1, sepi%natorb
       i2 = i2 + 1
       i1 = se_orbital_pointer(i1L)
       ksi_block(i1,i1) = ksi_block(i1,i1) + e1b(i2) * rp
       ecore = ecore + e1b(i2) * rp * pi_block(i1,i1)
    END DO
    !
    ! Add the electron-nuclear residual attraction term for atom sepj
    !
    i2 = 0
    DO i1L = 1, sepj%natorb
       i2 = i2 + 1
       i1 = se_orbital_pointer(i1L)
       ksj_block(i1,i1) = ksj_block(i1,i1) + e2a(i2) * rp
       ecore = ecore + e2a(i2) * rp * pj_block(i1,i1)
    END DO

  END SUBROUTINE fock2_1el_r3

! *****************************************************************************
!> \brief  Derivatives of 2-center 1-electron Fock Matrix residual (1/R^3) 
!>         integral part 
!> \author Teodoro Laino [tlaino]
!> \date   12.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2_1el_r3 (sepi, sepj, drp, pi_block, pj_block, force, e1b, e2a,&
       error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: drp
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(IN)                             :: pi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(IN)                             :: pj_block
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: e1b, e2a
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2_1el_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tmp

    failure = .FALSE.
    !
    ! Add the electron-nuclear residual attraction term for atom sepi
    !
    i2 = 0
    DO i1L = 1, sepi%natorb
       i1 = se_orbital_pointer(i1L)
       i2 = i2 + 1
       tmp = e1b(i2) * pi_block(i1,i1)
       force(1) = force(1) + tmp * drp(1)
       force(2) = force(2) + tmp * drp(2)
       force(3) = force(3) + tmp * drp(3)
    END DO
    !
    ! Add the electron-nuclear attraction term for atom sepj
    !
    i2 = 0
    DO i1L = 1, sepj%natorb
       i1 = se_orbital_pointer(i1L)
       i2 = i2 + 1
       tmp = e2a(i2) * pj_block(i1,i1)
       force(1) = force(1) + tmp * drp(1)
       force(2) = force(2) + tmp * drp(2)
       force(3) = force(3) + tmp * drp(3)
    END DO
    
  END SUBROUTINE dfock2_1el_r3

! *****************************************************************************
!> \brief  Construction of 2-center Fock Matrix - Coulomb Terms for the residual
!>         (1/R^3) integral part 
!> \author Teodoro Laino [tlaino]
!> \date   12.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock2C_r3(sepi, sepj, switch, pi_tot, fi_mat, pj_tot, fj_mat, &
       factor, w, rp, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(INOUT)                          :: fi_mat
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pj_tot
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(INOUT)                          :: fj_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    REAL(KIND=dp), DIMENSION(81), INTENT(IN) :: w
    REAL(KIND=dp), INTENT(IN)                :: rp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2C_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, ind, j, k, kL, kr, &
                                                natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, w_l(81)

    failure = .FALSE.
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
       ! Reshuffle the integral array (natural storage order is sepi/sepj)
       kr = 0
       DO i = 1, sepj%natorb
          DO j = 1, sepi%natorb
             kr  = kr + 1
             ind = (j-1)*sepj%natorb+i
             w_l(kr) = w(ind)
          END DO
       END DO
    ELSE
       w_l = w
    END IF

    ! Modify the Fock Matrix
    kr = 0
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       DO kL = 1, natorb(2)
          k = se_orbital_pointer(kL)
          kr = kr + 1
          a = w_l(kr) * factor * rp
          ! Coulomb
          IF (.NOT.switch) THEN
             fi_mat(i,i) = fi_mat(i,i) + a * pj_tot(k,k)
             fj_mat(k,k) = fj_mat(k,k) + a * pi_tot(i,i)
          ELSE
             fj_mat(i,i) = fj_mat(i,i) + a * pi_tot(k,k)
             fi_mat(k,k) = fi_mat(k,k) + a * pj_tot(i,i)
          END IF
       END DO
    END DO

  END SUBROUTINE fock2C_r3

! *****************************************************************************
!> \brief  Derivatives of 2-center Fock Matrix - Coulomb Terms for the residual
!>         (1/R^3) integral part 
!> \author Teodoro Laino [tlaino]
!> \date   12.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2C_r3(sepi, sepj, switch, pi_tot, pj_tot, factor,  w, drp,&
       force, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot, pj_tot
    REAL(KIND=dp), INTENT(IN)                :: factor
    REAL(KIND=dp), DIMENSION(81), INTENT(IN) :: w
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: drp
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2C_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, ind, j, k, kL, kr, &
                                                natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a(3), tmp, w_l(81)

    failure = .FALSE.
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
       ! Reshuffle the integral array (natural storage order is sepi/sepj)
       kr = 0
       DO i = 1, sepj%natorb
          DO j = 1, sepi%natorb
             kr  = kr + 1
             ind = (j-1)*sepj%natorb+i
             w_l(kr) = w(ind)
          END DO
       END DO
    ELSE
       w_l = w
    END IF

    ! Modify the Fock Matrix
    kr = 0
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       DO kL = 1, natorb(2)
          k = se_orbital_pointer(kL)
          kr = kr + 1
          tmp  = w_l(kr) * factor
          a(1) = tmp * drp(1)
          a(2) = tmp * drp(2)
          a(3) = tmp * drp(3)
          ! Coulomb
          IF (.NOT.switch) THEN
             tmp = pj_tot(k,k) * pi_tot(i,i)
          ELSE
             tmp = pi_tot(k,k) * pj_tot(i,i)
          END IF
          force(1) = force(1) + a(1) * tmp
          force(2) = force(2) + a(2) * tmp
          force(3) = force(3) + a(3) * tmp
       END DO
    END DO

  END SUBROUTINE dfock2C_r3
  
END MODULE nddo_methods

