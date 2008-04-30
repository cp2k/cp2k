!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the Fock matrix for NDDO methods
!> \author JGH
!> \par History
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : d-orbitals
! *****************************************************************************
MODULE nddo_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             semi_empirical_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE f77_blas
  USE input_constants,                 ONLY: do_method_am1,&
                                             do_method_mndo,&
                                             do_method_mndod,&
                                             do_method_pdg,&
                                             do_method_pm3,&
                                             do_method_pm6,&
                                             do_method_rm1
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
  USE semi_empirical_int_arrays,       ONLY: rij_threshold,&
                                             se_orbital_pointer
  USE semi_empirical_integrals,        ONLY: corecore,&
                                             dcorecore,&
                                             drotint,&
                                             drotnuc,&
                                             rotint,&
                                             rotnuc
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_type,&
                                             taper_type
  USE semi_empirical_utils,            ONLY: get_se_type
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             cp_sm_scale_and_add,&
                                             cp_sm_sm_trace,&
                                             create_dbm_from_sm,&
                                             deallocate_dbm_from_sm,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nddo_methods'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
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
                                                nspins, output_unit, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: qmmm_el
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers, &
                                                scaling_factor
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
    TYPE(section_vals_type), POINTER         :: scf_section

    CALL timeset(routineN,"I","",handle)
    NULLIFY(matrix_h, dft_control, logger, scf_section)
    logger => cp_error_get_logger(error)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF ( .NOT. failure ) THEN

       CALL get_qs_env(qs_env=qs_env,&
                      dft_control=dft_control,&
                      matrix_h=matrix_h,&
                      para_env=para_env,&
                      error=error)
       
       scf_section => section_vals_get_subs_vals(qs_env%input,"DFT%SCF",error=error)
       CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)
       nspins   = dft_control%nspins
       matrix_p => rho%rho_ao
       energy%qmmm_el = 0.0_dp
       DO ispin=1,nspins
          ! Copy the core matrix into the fock matrix
          CALL transfer_matrix(matrix_h(1)%matrix,ks_matrix(ispin)%matrix,error=error)
       END DO

       CALL build_nddo_fock_matrix_exchange(ks_env,qs_env,ks_matrix,rho,energy,&
            calculate_forces,just_energy,error)
       CALL build_nddo_fock_matrix_coulomb (ks_env,qs_env,ks_matrix,rho,energy,&
            calculate_forces,just_energy,error)
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

       energy%mulliken = 0.0_dp
       energy%exc      = 0.0_dp

       energy%total = energy%core +&
                      energy%core_overlap +&
                      0.5_dp*energy%hartree + &
                      energy%qmmm_el +&
                      energy%mulliken

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
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%DETAILED_ENERGY", error=error)
       ! here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C (plus occupation numbers)
       IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
          CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array,error=error)
          DO ispin=1,SIZE(mo_derivs)
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff, &
                             occupation_numbers=occupation_numbers )
             CALL cp_fm_get_info(mo_coeff,ncol_global=ncol_global,error=error)
             CALL cp_sm_fm_multiply(ks_matrix(ispin)%matrix,mo_coeff,&
                                    mo_derivs(ispin)%matrix, ncol_global,error=error)
             ALLOCATE(scaling_factor(SIZE(occupation_numbers)),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             scaling_factor=2.0_dp*occupation_numbers
             CALL cp_fm_column_scale(mo_derivs(ispin)%matrix,scaling_factor)
             DEALLOCATE(scaling_factor,stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDDO
       ENDIF

    END IF

    CALL timestop(0.0_dp,handle)

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
    LOGICAL                                  :: anag, defined, do_ewald, &
                                                failure
    REAL(KIND=dp)                            :: delta, dr1, enuc, enuclear
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_2c
    TYPE(neighbor_list_type), POINTER        :: sab_2c_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_2c_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(taper_type), POINTER                :: taper

    failure=.FALSE.
    NULLIFY(dft_control,cell,force,particle_set)
    NULLIFY(se_control,taper,atomic_kind_set)

    CALL timeset(routineN,"I","",handle)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, taper=taper,&
         energy=energy, error=error)

    ! Parameters
    se_control => dft_control%qs_control%se_control
    do_ewald   =  se_control%do_ewald
    anag       =  se_control%analytical_gradients

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

       enuclear = 0.0_dp
       itype    = get_se_type(dft_control%qs_control%method_id)

       IF (do_ewald) THEN
          !calculate 1/r contribution using Ewald type method
!deb      CALL nddo_ewald
       END IF

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
                          do_method_rm1,do_method_mndod)

                        ! Core-Core energy term
                        CALL corecore (se_kind_a,se_kind_b,rij,enuc=enuc,itype=itype,anag=anag,&
                             shortrange=do_ewald, taper=taper, error=error)
                        enuclear = enuclear + enuc

                        ! Core-Core Derivatives
                        IF(calculate_forces) THEN
                           atom_a = atom_of_kind(iatom)
                           atom_b = atom_of_kind(jatom)
                           
                           CALL dcorecore (se_kind_a,se_kind_b,rij,denuc=force_ab,itype=itype,delta=delta,&
                                anag=anag,shortrange=do_ewald,taper=taper,error=error)
                           force(ikind)%all_potential(:,atom_a) = force(ikind)%all_potential(:,atom_a) - force_ab(:)
                           force(jkind)%all_potential(:,atom_b) = force(jkind)%all_potential(:,atom_b) + force_ab(:)
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
      energy%core_overlap = enuclear      
   END IF
   
   CALL timestop(0.0_dp,handle)
 END SUBROUTINE se_core_core_interaction

! *****************************************************************************
!> \brief Construction of the Coulomb part of the Fock matrix
!> \author JGH
! *****************************************************************************
  SUBROUTINE build_nddo_fock_matrix_coulomb (ks_env,qs_env,ks_matrix,rho,energy,&
       calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'build_nddo_fock_matrix_coulomb', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iab, iatom, ikind, ilist, inode, &
      ispin, itype, jatom, jkind, natom, natorb_a, natorb_b, nkind, nlist, &
      nnode, nspins, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: anag, check, defined, &
                                                do_ewald, failure, switch
    REAL(KIND=dp)                            :: delta, dr1, ecore2, ecoul
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    REAL(KIND=dp), DIMENSION(45, 45)         :: p_block_tot_a, p_block_tot_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: ksa_block_a, ksa_block_b, &
      ksb_block_a, ksb_block_b, pa_block_a, pa_block_b, pb_block_a, pb_block_b
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_2c
    TYPE(neighbor_list_type), POINTER        :: sab_2c_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_2c_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: diagmat_ks, diagmat_p, &
                                                matrix_p
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(taper_type), POINTER                :: taper

    failure=.FALSE.
    NULLIFY(matrix_p,dft_control,cell,force,particle_set)
    NULLIFY(diagmat_ks,diagmat_p, se_control,taper)

    CALL timeset(routineN,"I","",handle)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
    
    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, taper=taper,&
         para_env=para_env, error=error)

    ! Parameters
    se_control => dft_control%qs_control%se_control
    do_ewald   =  se_control%do_ewald
    anag       =  se_control%analytical_gradients

    nspins=dft_control%nspins
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)
    
    IF ( .NOT. failure ) THEN
       CALL allocate_matrix_set(diagmat_ks,nspins,error=error)
       CALL allocate_matrix_set(diagmat_p,nspins,error=error)
       
       matrix_p => rho%rho_ao
       
       CALL get_qs_env(qs_env=qs_env,&
                       sab_2c=sab_2c,&
                       atomic_kind_set=atomic_kind_set,error=error)

       nkind = SIZE(atomic_kind_set)

       IF(calculate_forces) THEN
          CALL get_qs_env(qs_env=qs_env,&
                          particle_set=particle_set,&
                          force=force,error=error)
          natom = SIZE (particle_set)
          ALLOCATE (atom_of_kind(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          delta = se_control%delta
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)
       END IF

       DO ispin=1,nspins
          ! Allocate diagonal block matrices
          CALL create_dbm_from_sm(diagmat_p(ispin)%matrix,matrix_p(ispin)%matrix,&
                                  "DIAGONAL BLOCK P MATRIX",.TRUE.,para_env,error=error)
          CALL create_dbm_from_sm(diagmat_ks(ispin)%matrix,ks_matrix(ispin)%matrix,&
                                  "DIAGONAL BLOCK KS MATRIX",.FALSE.,error=error)
       END DO
       
       ecore2   = 0.0_dp
       itype    = get_se_type(dft_control%qs_control%method_id)

       IF (do_ewald) THEN
          !calculate 1/r contribution using Ewald type method
!deb      CALL nddo_ewald
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

               CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
               IF ( nspins >= 1 ) THEN
                  CALL get_block_node(matrix=diagmat_p(1)%matrix,&
                                      block_row=iatom,&
                                      block_col=iatom,&
                                      BLOCK=pa_block_a)
                  CPPostcondition(ASSOCIATED(pa_block_a),cp_failure_level,routineP,error,failure)
                  CALL get_block_node(matrix=diagmat_ks(1)%matrix,&
                                      block_row=iatom,&
                                      block_col=iatom,&
                                      BLOCK=ksa_block_a)
                  CPPostcondition(ASSOCIATED(ksa_block_a),cp_failure_level,routineP,error,failure)
                  p_block_tot_a(1:SIZE(pa_block_a,1),1:SIZE(pa_block_a,2)) = 2.0_dp * pa_block_a
               END IF
               IF ( nspins >= 2 ) THEN
                  CALL get_block_node(matrix=diagmat_p(2)%matrix,&
                                      block_row=iatom,&
                                      block_col=iatom,&
                                      BLOCK=pa_block_b)
                  CPPostcondition(ASSOCIATED(pa_block_b),cp_failure_level,routineP,error,failure)
                  CALL get_block_node(matrix=diagmat_ks(2)%matrix,&
                                      block_row=iatom,&
                                      block_col=iatom,&
                                      BLOCK=ksa_block_b)
                  CPPostcondition(ASSOCIATED(ksa_block_b),cp_failure_level,routineP,error,failure)
                  check = (SIZE(pa_block_a,1)==SIZE(pa_block_b,1)).AND.(SIZE(pa_block_a,2)==SIZE(pa_block_b,2))
                  CPPostcondition(check,cp_failure_level,routineP,error,failure)
                  p_block_tot_a(1:SIZE(pa_block_a,1),1:SIZE(pa_block_a,2)) = pa_block_a + pa_block_b
               END IF

               DO inode=1,nnode
                  CALL get_neighbor_node(neighbor_node=sab_2c_neighbor_node,neighbor=jatom,r=rij)

                  dr1 = DOT_PRODUCT(rij,rij)
                  IF ( dr1 > rij_threshold ) THEN

                     IF (iatom <= jatom) THEN
                        switch = .FALSE.
                     ELSE
                        switch = .TRUE.
                     END IF
                     CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                     IF ( nspins >= 1 ) THEN
                        CALL get_block_node(matrix=diagmat_p(1)%matrix,&
                                            block_row=jatom,&
                                            block_col=jatom,&
                                            BLOCK=pb_block_a)
                        CPPostcondition(ASSOCIATED(pb_block_a),cp_failure_level,routineP,error,failure)
                        CALL get_block_node(matrix=diagmat_ks(1)%matrix,&
                                            block_row=jatom,&
                                            block_col=jatom,&
                                            BLOCK=ksb_block_a)
                        CPPostcondition(ASSOCIATED(ksb_block_a),cp_failure_level,routineP,error,failure)
                        p_block_tot_b(1:SIZE(pb_block_a,1),1:SIZE(pb_block_a,2)) = 2.0_dp * pb_block_a
                     END IF
                     IF ( nspins >= 2 ) THEN
                        CALL get_block_node(matrix=diagmat_p(2)%matrix,&
                                            block_row=jatom,&
                                            block_col=jatom,&
                                            BLOCK=pb_block_b)
                        CPPostcondition(ASSOCIATED(pb_block_b),cp_failure_level,routineP,error,failure)
                        CALL get_block_node(matrix=diagmat_ks(2)%matrix,&
                                            block_row=jatom,&
                                            block_col=jatom,&
                                            BLOCK=ksb_block_b)
                        CPPostcondition(ASSOCIATED(ksb_block_b),cp_failure_level,routineP,error,failure)
                        check = (SIZE(pb_block_a,1)==SIZE(pb_block_b,1)).AND.(SIZE(pb_block_a,2)==SIZE(pb_block_b,2))
                        CPPostcondition(check,cp_failure_level,routineP,error,failure)
                        p_block_tot_b(1:SIZE(pb_block_a,1),1:SIZE(pb_block_a,2)) = pb_block_a + pb_block_b
                     END IF

                     SELECT CASE (dft_control%qs_control%method_id)
                     CASE (do_method_mndo, do_method_am1, do_method_pm3, do_method_pm6, do_method_pdg,&
                           do_method_rm1, do_method_mndod)

                        ! One-centre One-electron terms
                        IF      ( nspins == 1 ) THEN
                           CALL fock1_1el (se_kind_a,se_kind_b,rij, ksa_block_a, ksb_block_a,&
                                pa_block_a, pb_block_a, ecore=ecore2, itype=itype, anag=anag,&
                                shortrange=do_ewald, taper=taper, error=error)
                        ELSE IF ( nspins == 2 ) THEN
                           CALL fock1_1el (se_kind_a,se_kind_b,rij, ksa_block_a, ksb_block_a,&
                                pa_block_a, pb_block_a, ecore=ecore2, itype=itype, anag=anag,&
                                shortrange=do_ewald, taper=taper, error=error)

                           CALL fock1_1el (se_kind_a,se_kind_b,rij, ksa_block_b, ksb_block_b,&
                                pa_block_b, pb_block_b, ecore=ecore2, itype=itype, anag=anag,&
                                shortrange=do_ewald, taper=taper, error=error)
                        END IF
                        ! Coulomb Terms
                        IF      ( nspins == 1 ) THEN
                           CALL fock2(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, pa_block_a, ksa_block_a,&
                                p_block_tot_b, ksb_block_a, factor=0.5_dp, anag=anag, shortrange=do_ewald, taper=taper,&
                                error=error)
                        ELSE IF ( nspins == 2 ) THEN
                           CALL fock2(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, pa_block_a, ksa_block_a,&
                                p_block_tot_b, ksb_block_a, factor=1.0_dp, anag=anag, shortrange=do_ewald, taper=taper,&
                                error=error)

                           CALL fock2(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, pa_block_b, ksa_block_b,&
                                p_block_tot_b, ksb_block_b, factor=1.0_dp, anag=anag, shortrange=do_ewald, taper=taper,&
                                error=error)
                        END IF

                        IF(calculate_forces) THEN
                           atom_a   = atom_of_kind(iatom)
                           atom_b   = atom_of_kind(jatom)

                           ! Derivatives of the One-centre One-electron terms
                           force_ab = 0.0_dp
                           IF      ( nspins == 1 ) THEN
                              CALL dfock1_1el (se_kind_a,se_kind_b,rij, pa_block_a, pb_block_a, itype=itype, anag=anag,&
                                   shortrange=do_ewald, taper=taper, force=force_ab, delta=delta, error=error)
                           ELSE IF ( nspins == 2 ) THEN
                              CALL dfock1_1el (se_kind_a,se_kind_b,rij, pa_block_a, pb_block_a, itype=itype, anag=anag,&
                                   shortrange=do_ewald, taper=taper, force=force_ab, delta=delta, error=error)

                              CALL dfock1_1el (se_kind_a,se_kind_b,rij, pa_block_b, pb_block_b, itype=itype, anag=anag,&
                                   shortrange=do_ewald, taper=taper, force=force_ab, delta=delta, error=error)
                           END IF
                           force_ab = -force_ab
                           force(ikind)%all_potential(:,atom_a) = force(ikind)%all_potential(:,atom_a) - force_ab(:)
                           force(jkind)%all_potential(:,atom_b) = force(jkind)%all_potential(:,atom_b) + force_ab(:)

                           ! Derivatives of the Coulomb Terms
                           force_ab = 0._dp
                           IF      ( nspins == 1 ) THEN
                              CALL dfock2(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, pa_block_a, p_block_tot_b,&
                                   factor=0.25_dp, anag=anag, shortrange=do_ewald, taper=taper, force=force_ab, delta=delta,&
                                   error=error)
                           ELSE IF ( nspins == 2 ) THEN
                              CALL dfock2(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, pa_block_a, p_block_tot_b,&
                                   factor=0.5_dp, anag=anag, shortrange=do_ewald, taper=taper, force=force_ab, delta=delta,&
                                   error=error)
                              
                              CALL dfock2(se_kind_a, se_kind_b, rij, switch, p_block_tot_a, pa_block_b, p_block_tot_b,&
                                   factor=0.5_dp, anag=anag, shortrange=do_ewald, taper=taper, force=force_ab, delta=delta,&
                                   error=error)
                           END IF
                           IF ( switch ) force_ab = -force_ab
                           force(ikind)%rho_elec(:,atom_a) = force(ikind)%rho_elec(:,atom_a) - force_ab(:)
                           force(jkind)%rho_elec(:,atom_b) = force(jkind)%rho_elec(:,atom_b) + force_ab(:)
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
      
      DO ispin = 1, nspins
         CALL deallocate_dbm_from_sm(diagmat_p(ispin)%matrix,matrix_p(1)%matrix,.FALSE.,error=error)
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
      
      CALL mp_sum(ecore2,para_env%group)
      
      energy%hartree = ecore2 - energy%core
      DO ispin=1,nspins
         CALL cp_sm_sm_trace(ks_matrix(ispin)%matrix,matrix_p(ispin)%matrix,trace=ecoul,&
              para_env=para_env,error=error)
         energy%hartree=energy%hartree+ecoul
      END DO
   END IF
   
   CALL timestop(0.0_dp,handle)
 END SUBROUTINE build_nddo_fock_matrix_coulomb
 
! *****************************************************************************
!> \brief Construction of the Exchange part of the Fock matrix
!> \author JGH
! *****************************************************************************
  SUBROUTINE build_nddo_fock_matrix_exchange (ks_env,qs_env,ks_matrix,rho,energy,&
               calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'build_nddo_fock_matrix_exchange', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iab, iatom, icol, ikind, ilist, inode, &
      irow, jatom, jkind, natom, natorb_a, natorb_b, nkind, nlist, nnode, &
      nspins, stat
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
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(taper_type), POINTER                :: taper

    failure=.FALSE.
    NULLIFY(matrix_p,dft_control,cell,force,particle_set,se_control)
    NULLIFY(taper)
    CALL timeset(routineN,"I","",handle)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
         cell=cell,taper=taper,para_env=para_env,error=error)

    ! set values for tapering function
    se_control => dft_control%qs_control%se_control
    anag       =  se_control%analytical_gradients

    nspins=dft_control%nspins
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)

    IF ( .NOT. failure ) THEN

       matrix_p => rho%rho_ao

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
                   CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                   IF ( nspins >= 1 ) THEN
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
                      p_block_tot(1:SIZE(p_block_a,1),1:SIZE(p_block_a,2)) = 2.0_dp * p_block_a
                   END IF
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
                      check = (SIZE(p_block_a,1)==SIZE(p_block_b,1)).AND.(SIZE(p_block_a,2)==SIZE(p_block_b,2))
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
                         CALL fock2(se_kind_a, se_kind_b, rij, switch, p_block_tot, p_block_a, ks_block_a,&
                              factor=0.5_dp, anag=anag, shortrange=.FALSE., taper=taper, error=error)
                      ELSE IF ( nspins == 2 ) THEN
                         CALL fock2(se_kind_a, se_kind_b, rij, switch, p_block_tot, p_block_a, ks_block_a,&
                              factor=1.0_dp, anag=anag, shortrange=.FALSE., taper=taper, error=error)
                         
                         CALL fock2(se_kind_a, se_kind_b, rij, switch, p_block_tot, p_block_b, ks_block_b,&
                              factor=1.0_dp, anag=anag, shortrange=.FALSE., taper=taper, error=error)
                      END IF
                      IF(calculate_forces) THEN
                         atom_a = atom_of_kind(iatom)
                         atom_b = atom_of_kind(jatom)
                         force_ab = 0.0_dp
                         IF      ( nspins == 1 ) THEN
                            CALL dfock2(se_kind_a, se_kind_b, rij, switch, p_block_tot, p_block_a, factor=0.5_dp,&
                                 anag=anag, shortrange=.FALSE., taper=taper, force=force_ab, delta=delta, error=error)
                         ELSE IF ( nspins == 2 ) THEN
                            CALL dfock2(se_kind_a, se_kind_b, rij, switch, p_block_tot, p_block_a, factor=1.0_dp,&
                                 anag=anag, shortrange=.FALSE., taper=taper, force=force_ab, delta=delta, error=error)
                            
                            CALL dfock2(se_kind_a, se_kind_b, rij, switch, p_block_tot, p_block_b, factor=1.0_dp,&
                                 anag=anag, shortrange=.FALSE., taper=taper, force=force_ab, delta=delta, error=error)
                         END IF
                         IF (switch) force_ab = -force_ab
                         force(ikind)%rho_elec(:,atom_a) = force(ikind)%rho_elec(:,atom_a) - force_ab(:)
                         force(jkind)%rho_elec(:,atom_b) = force(jkind)%rho_elec(:,atom_b) + force_ab(:)
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
    CALL timestop(0.0_dp,handle)
    
  END SUBROUTINE build_nddo_fock_matrix_exchange

! *****************************************************************************
!> \brief Construction of 1-center 1-electron Fock Matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock1_1el (sepi, sepj, rij, ksi_block, ksj_block, pi_block,&
       pj_block, ecore, itype, anag, shortrange, taper, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: ksi_block, ksj_block
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: pi_block, pj_block
    REAL(KIND=dp), INTENT(INOUT)             :: ecore
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock1_1el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2, j1, j1L
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(45)             :: e1b, e2a

    failure = .FALSE.
    CALL rotnuc (sepi, sepj, rij, e1b=e1b, e2a=e2a, itype=itype,&
         anag=anag, shortrange=shortrange, taper=taper, error=error)
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

  END SUBROUTINE fock1_1el

! *****************************************************************************
!> \brief Derivatives of 1-center 1-electron Fock Matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dfock1_1el (sepi, sepj, rij, pi_block, pj_block, itype, anag,&
       shortrange, taper, force, delta, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: pi_block, pj_block
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), INTENT(IN)                :: delta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock1_1el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2, j1, j1L
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(45, 3)          :: de1b, de2a

    failure = .FALSE.
    CALL drotnuc (sepi, sepj, rij, de1b=de1b, de2a=de2a, itype=itype, anag=anag,&
         shortrange=shortrange, taper=taper, delta=delta, error=error)
    !
    ! Add the electron-nuclear attraction term for atom sepi
    !
    i2 = 0
    DO i1L = 1, sepi%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          force = force - 2.0_dp * de1b(i2,:) * pi_block(i1,j1)
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       force = force - de1b(i2,:) * pi_block(i1,j1)
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
          force = force - 2.0_dp * de2a(i2,:) * pj_block(i1,j1)
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       force = force - de2a(i2,:) * pj_block(i1,j1)
    END DO
    
  END SUBROUTINE dfock1_1el

! *****************************************************************************
!> \brief Construction of 1-center 2-electron Fock Matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock1_2el(sep, p_tot, p_mat, f_mat, factor, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: p_tot, p_mat
    REAL(KIND=dp), DIMENSION(:, :), &
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
!> \brief Construction of 2-center Fock Matrix - General Driver
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock2(sepi, sepj, rij, switch, pi_tot, pi_mat, fi_mat, pj_tot, fj_mat, &
       factor, anag, shortrange, taper, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: pi_tot, pi_mat
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: fi_mat
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: pj_tot
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: fj_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    LOGICAL                                  :: failure, task_coulomb
    REAL(KIND=dp)                            :: a, aa, bb
    REAL(KIND=dp), DIMENSION(2025)           :: w

    failure = .FALSE.
    CPPrecondition(PRESENT(pj_tot).EQV.PRESENT(fj_mat),cp_failure_level,routineP,error,failure)
    task_coulomb = PRESENT(pj_tot)
    
    ! Evaluate integrals
    IF (.NOT.switch) THEN
       CALL rotint (sepi,sepj, rij,w,anag=anag,shortrange=shortrange,taper=taper,error=error)
    ELSE
       CALL rotint (sepj,sepi,-rij,w,anag=anag,shortrange=shortrange,taper=taper,error=error)
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
                IF (task_coulomb) THEN
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
                ELSE
                   ! Exchange
                   a = a * aa * bb * 0.25_dp
                   fi_mat(i,k) = fi_mat(i,k) - a * pi_mat(j,l)
                   fi_mat(i,l) = fi_mat(i,l) - a * pi_mat(j,k)
                   fi_mat(j,k) = fi_mat(j,k) - a * pi_mat(i,l)
                   fi_mat(j,l) = fi_mat(j,l) - a * pi_mat(i,k)
                END IF
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE fock2

! *****************************************************************************
!> \brief Derivatives of 2-center Fock Matrix - General Driver
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2(sepi, sepj, rij, switch, pi_tot, pi_mat, pj_tot, &
       factor, anag, shortrange, taper, force, delta, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: pi_tot, pi_mat
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: pj_tot
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), INTENT(IN)                :: delta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    LOGICAL                                  :: failure, task_coulomb
    REAL(KIND=dp)                            :: aa, bb
    REAL(KIND=dp), DIMENSION(2025, 3)        :: dw
    REAL(KIND=dp), DIMENSION(3)              :: a

    failure = .FALSE.
    task_coulomb = PRESENT(pj_tot)
    
    ! Evaluate integrals
    IF (.NOT.switch) THEN
       CALL drotint (sepi,sepj, rij,dw,delta,anag=anag,shortrange=shortrange,taper=taper,error=error)
    ELSE
       CALL drotint (sepj,sepi,-rij,dw,delta,anag=anag,shortrange=shortrange,taper=taper,error=error)
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
                a = dw(kr,:)*factor
                IF (task_coulomb) THEN
                   ! Coulomb
                   IF (.NOT.switch) THEN
                      force = force + bb * aa * a * pj_tot(k,l) * pi_tot(i,j)
                   ELSE
                      force = force + bb * aa * a * pi_tot(k,l) * pj_tot(i,j)
                   END IF
                ELSE
                   ! Exchange
                   a = a * aa * bb * 0.25_dp
                   force = force - a * pi_mat(j,l) * pi_mat(i,k)
                   force = force - a * pi_mat(j,k) * pi_mat(i,l)
                   force = force - a * pi_mat(i,l) * pi_mat(j,k)
                   force = force - a * pi_mat(i,k) * pi_mat(j,l)
                END IF
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE dfock2

END MODULE nddo_methods

