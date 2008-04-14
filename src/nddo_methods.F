!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the Fock matrix for NDDO methods
!> \author JGH
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
  USE semi_empirical_integrals,        ONLY: corecore,&
                                             dcorecore,&
                                             drotint,&
                                             drotnuc,&
                                             rotint,&
                                             rotnuc
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_type,&
                                             taper_type
  USE semi_empirical_utils,            ONLY: contract,&
                                             expand,&
                                             get_se_type,&
                                             traceb
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

  PUBLIC :: build_nddo_fock_matrix, se_core_core_interaction

  INTERFACE fock2c
    MODULE PROCEDURE rfock2c,  ufock2c
  END INTERFACE

  INTERFACE dfock2c
    MODULE PROCEDURE rdfock2c, udfock2c
  END INTERFACE

  INTERFACE fock2e
    MODULE PROCEDURE rfock2e,  ufock2e
  END INTERFACE

  INTERFACE dfock2e
    MODULE PROCEDURE rdfock2e, udfock2e
  END INTERFACE

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

       energy%mulliken = 0._dp
       energy%exc      = 0._dp

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
                                                failure, switch
    REAL(KIND=dp)                            :: delta, dr1, enuc, enuclear
    REAL(KIND=dp), DIMENSION(3)              :: denuc, force_ab, rij
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_2c
    TYPE(neighbor_list_type), POINTER        :: sab_2c_neighbor_list, &
                                                sab_2c_neighbor_list_local
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
       CALL get_qs_env(qs_env=qs_env,&
                       sab_2c=sab_2c,&
                       atomic_kind_set=atomic_kind_set,&
                       error=error)

       nkind = SIZE(atomic_kind_set)
       ! Possibly compute forces
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

       enuclear = 0.0_dp
       itype    = get_se_type(dft_control%qs_control%method_id)

       IF (do_ewald) THEN
          !calculate 1/r contribution using Ewald type method
!deb      CALL nddo_ewald
       END IF

       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,defined=defined)

          IF (.NOT.defined) CYCLE

          DO jkind=1,nkind
             atomic_kind => atomic_kind_set(jkind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                  se_parameter=se_kind_b)
             CALL get_se_param(se_kind_b,defined=defined)

             IF (.NOT.defined) CYCLE

             iab = ikind + nkind*(jkind - 1)
             IF (.NOT.ASSOCIATED(sab_2c(iab)%neighbor_list_set)) CYCLE
             CALL get_neighbor_list_set(neighbor_list_set=sab_2c(iab)%neighbor_list_set,&
                                        nlist=nlist)
             NULLIFY(sab_2c_neighbor_list)

             DO ilist=1,nlist
                IF ( .NOT. ASSOCIATED(sab_2c_neighbor_list) ) THEN
                   sab_2c_neighbor_list => first_list(sab_2c(iab)%neighbor_list_set)
                ELSE
                   sab_2c_neighbor_list => next(sab_2c_neighbor_list)
               END IF
               sab_2c_neighbor_list_local => sab_2c_neighbor_list
               CALL get_neighbor_list(neighbor_list=sab_2c_neighbor_list_local,&
                                      atom=iatom,&
                                      nnode=nnode)

               sab_2c_neighbor_node => first_node(sab_2c_neighbor_list_local)

               DO inode=1,nnode
                  CALL get_neighbor_node(neighbor_node=sab_2c_neighbor_node,&
                                         neighbor=jatom,r=rij)
                  dr1 = SUM(rij(:)**2)
                  IF ( dr1 > 0.0_dp ) THEN
                     IF (iatom <= jatom) THEN
                        switch = .FALSE.
                     ELSE
                        switch = .TRUE.
                     END IF

                     SELECT CASE (dft_control%qs_control%method_id)
                     CASE (do_method_mndo,do_method_am1,do_method_pm3,do_method_pm6,do_method_pdg,&
                          do_method_rm1,do_method_mndod)
                        IF ( .NOT. switch ) THEN
                           CALL corecore (se_kind_a,se_kind_b,rij,enuc=enuc,itype=itype,anag=anag,&
                                shortrange=do_ewald, taper=taper, error=error)
                        ELSE
                           CALL corecore (se_kind_b,se_kind_a,-rij,enuc=enuc,itype=itype,anag=anag,&
                                shortrange=do_ewald, taper=taper, error=error)
                        END IF
                        enuclear = enuclear + enuc
                        ! one-centre one-electron terms
                        IF(calculate_forces) THEN
                           atom_a = atom_of_kind(iatom)
                           atom_b = atom_of_kind(jatom)
                           
                           IF ( .NOT. switch ) THEN
                              CALL dcorecore (se_kind_a,se_kind_b,rij,denuc=denuc,itype=itype,delta=delta,&
                                   anag=anag,shortrange=do_ewald,taper=taper,error=error)
                           ELSE
                              CALL dcorecore (se_kind_b,se_kind_a,-rij,denuc=denuc,itype=itype,delta=delta,&
                                   anag=anag,shortrange=do_ewald,taper=taper,error=error)
                           END IF
                           force_ab(1:3)=-denuc(1:3)
                           IF ( .NOT. switch ) force_ab = -force_ab
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

    INTEGER :: atom_a, atom_b, handle, iab, iatom, icor, ikind, ilist, inode, &
      ispin, itype, jatom, jkind, natom, natorb_a, natorb_b, nkind, nlist, &
      nnode, nspins, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: anag, defined, do_ewald, &
                                                failure, switch
    REAL(KIND=dp)                            :: delta, dr1, ecore2, ecoul
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    REAL(KIND=dp), DIMENSION(45)             :: e1b, e2a, pvec
    REAL(KIND=dp), DIMENSION(45, 3)          :: de1b, de2a
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
    TYPE(neighbor_list_type), POINTER        :: sab_2c_neighbor_list, &
                                                sab_2c_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_2c_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: diagmat_ks, diagmat_p, &
                                                matrix_h, matrix_p
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(taper_type), POINTER                :: taper

    failure=.FALSE.
    NULLIFY(matrix_h,matrix_p,dft_control,cell,force,particle_set)
    NULLIFY(diagmat_ks,diagmat_p, se_control,taper)

    CALL timeset(routineN,"I","",handle)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
    
    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, matrix_h=matrix_h,&
         cell=cell, taper=taper,  para_env=para_env, error=error)

    ! Parameters
    se_control => dft_control%qs_control%se_control
    do_ewald   =  se_control%do_ewald
    anag       =  se_control%analytical_gradients

    nspins=dft_control%nspins
    CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
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
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,&
                            defined=defined,&
                            natorb=natorb_a)

          IF (.NOT.defined .OR. natorb_a < 1) CYCLE

          DO jkind=1,nkind
             atomic_kind => atomic_kind_set(jkind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                  se_parameter=se_kind_b)
             CALL get_se_param(se_kind_b,&
                               defined=defined,&
                               natorb=natorb_b)

             IF (.NOT.defined .OR. natorb_b < 1) CYCLE

             iab = ikind + nkind*(jkind - 1)
             IF (.NOT.ASSOCIATED(sab_2c(iab)%neighbor_list_set)) CYCLE
             CALL get_neighbor_list_set(neighbor_list_set=&
                                        sab_2c(iab)%neighbor_list_set,&
                                        nlist=nlist)
             NULLIFY(sab_2c_neighbor_list)

             DO ilist=1,nlist
                IF ( .NOT. ASSOCIATED(sab_2c_neighbor_list) ) THEN
                   sab_2c_neighbor_list => first_list(sab_2c(iab)%neighbor_list_set)
                ELSE
                   sab_2c_neighbor_list => next(sab_2c_neighbor_list)
               END IF
               sab_2c_neighbor_list_local => sab_2c_neighbor_list
               CALL get_neighbor_list(neighbor_list=sab_2c_neighbor_list_local,&
                                      atom=iatom, nnode=nnode)

               sab_2c_neighbor_node => first_node(sab_2c_neighbor_list_local)

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
               END IF

               DO inode=1,nnode
                  CALL get_neighbor_node(neighbor_node=sab_2c_neighbor_node,&
                                         neighbor=jatom,r=rij)

                  dr1 = SUM(rij(:)**2)
                  IF ( dr1 > 0.0_dp ) THEN

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
                     END IF

                     SELECT CASE (dft_control%qs_control%method_id)
                     CASE (do_method_mndo, do_method_am1, do_method_pm3, do_method_pm6, do_method_pdg,&
                           do_method_rm1, do_method_mndod)
                        IF ( .NOT. switch ) THEN
                           CALL rotnuc (se_kind_a,se_kind_b,rij,e1b=e1b,e2a=e2a,itype=itype,&
                                anag=anag, shortrange=do_ewald, taper=taper, error=error)
                        ELSE
                           CALL rotnuc (se_kind_b,se_kind_a,-rij,e1b=e2a,e2a=e1b,itype=itype,&
                                anag=anag, shortrange=do_ewald, taper=taper, error=error)
                        END IF
                        ! one-centre one-electron terms
                        CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                        IF ( nspins >= 1 ) THEN
                           ! contribution to iatom
                           IF ( SIZE(ksa_block_a,1) == 1 ) THEN
                              ksa_block_a(1,1)=ksa_block_a(1,1)+e1b(1)
                              ecore2=ecore2+e1b(1)*pa_block_a(1,1)
                           ELSE
                              CALL expand ( ksa_block_a, e1b )
                              ecore2=ecore2+traceb(pa_block_a,e1b)
                           END IF
                           ! contribution to jatom
                           IF ( SIZE(ksb_block_a,1) == 1 ) THEN
                              ksb_block_a(1,1)=ksb_block_a(1,1)+e2a(1)
                              ecore2=ecore2+e2a(1)*pb_block_a(1,1)
                           ELSE
                              CALL expand ( ksb_block_a, e2a )
                              ecore2=ecore2+traceb(pb_block_a,e2a)
                           END IF
                        END IF
                        IF ( nspins >= 2 ) THEN
                           ! contribution to iatom
                           IF ( SIZE(ksa_block_b,1) == 1 ) THEN
                              ksa_block_b(1,1)=ksa_block_b(1,1)+e1b(1)
                              ecore2=ecore2+e1b(1)*pa_block_b(1,1)
                           ELSE
                              CALL expand ( ksa_block_b, e1b )
                              ecore2=ecore2+traceb(pa_block_b,e1b)
                           END IF
                           ! contribution to jatom
                           IF ( SIZE(ksb_block_b,1) == 1 ) THEN
                              ksb_block_b(1,1)=ksb_block_b(1,1)+e2a(1)
                              ecore2=ecore2+e2a(1)*pb_block_b(1,1)
                           ELSE
                              CALL expand ( ksb_block_b, e2a )
                              ecore2=ecore2+traceb(pb_block_b,e2a)
                           END IF
                        END IF
                        IF(calculate_forces) THEN
                           atom_a = atom_of_kind(iatom)
                           atom_b = atom_of_kind(jatom)
                           
                           IF ( .NOT. switch ) THEN
                              CALL drotnuc (se_kind_a,se_kind_b,rij,de1b=de1b,de2a=de2a,&
                                   itype=itype,delta=delta,anag=anag,shortrange=do_ewald,taper=taper,&
                                   error=error)
                           ELSE
                              CALL drotnuc (se_kind_b,se_kind_a,-rij,de1b=de2a,de2a=de1b,&
                                   itype=itype,delta=delta,anag=anag,shortrange=do_ewald,taper=taper,&
                                   error=error)
                           END IF
                           force_ab(1:3)=0.0_dp
                           CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                           IF ( nspins >= 1 ) THEN
                              IF ( SIZE(ksa_block_a,1) == 1 ) THEN
                                 force_ab(1:3)=force_ab(1:3)-de1b(1,1:3)*pa_block_a(1,1)
                              ELSE
                                 CALL contract ( pa_block_a, pvec )
                                 DO icor=1,3
                                    force_ab(icor)=force_ab(icor)-SUM(de1b(1:10,icor)*pvec(1:10))
                                 END DO
                              END IF
                              IF ( SIZE(ksb_block_a,1) == 1 ) THEN
                                 force_ab(1:3)=force_ab(1:3)-de2a(1,1:3)*pb_block_a(1,1)
                              ELSE
                                 CALL contract ( pb_block_a, pvec )
                                 DO icor=1,3
                                    force_ab(icor)=force_ab(icor)-SUM(de2a(1:10,icor)*pvec(1:10))
                                 END DO
                              END IF
                           END IF
                           IF ( nspins >= 2 ) THEN
                              IF ( SIZE(ksa_block_b,1) == 1 ) THEN
                                 force_ab(1:3)=force_ab(1:3)-de1b(1,1:3)*pa_block_b(1,1)
                              ELSE
                                 CALL contract ( pa_block_b, pvec )
                                 DO icor=1,3
                                    force_ab(icor)=force_ab(icor)-SUM(de1b(1:10,icor)*pvec(1:10))
                                 END DO
                              END IF
                              IF ( SIZE(ksb_block_b,1) == 1 ) THEN
                                 force_ab(1:3)=force_ab(1:3)-de2a(1,1:3)*pb_block_b(1,1)
                              ELSE
                                 CALL contract ( pb_block_b, pvec )
                                 DO icor=1,3
                                    force_ab(icor)=force_ab(icor)-SUM(de2a(1:10,icor)*pvec(1:10))
                                 END DO
                              END IF
                           END IF
                           IF ( .NOT. switch ) force_ab = -force_ab
                           force(ikind)%all_potential(:,atom_a) = force(ikind)%all_potential(:,atom_a) - force_ab(:)
                           force(jkind)%all_potential(:,atom_b) = force(jkind)%all_potential(:,atom_b) + force_ab(:)
                        END IF

                        CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                        IF ( nspins == 1 ) THEN
                           IF ( .NOT. switch ) THEN
                              CALL fock2c(se_kind_a,se_kind_b,rij,pa_block_a,pb_block_a,&
                                          ksa_block_a,ksb_block_a,anag,do_ewald,taper,error)
                           ELSE
                              CALL fock2c(se_kind_b,se_kind_a,-rij,pb_block_a,pa_block_a,&
                                          ksb_block_a,ksa_block_a,anag,do_ewald,taper,error)
                           ENDIF
                        ELSE IF ( nspins == 2 ) THEN
                           IF ( .NOT. switch ) THEN
                              CALL fock2c(se_kind_a,se_kind_b,rij,pa_block_a,pb_block_a,&
                                          ksa_block_a,ksb_block_a,pa_block_b,pb_block_b,&
                                          ksa_block_b,ksb_block_b,anag,do_ewald,taper,error)
                           ELSE
                              CALL fock2c(se_kind_b,se_kind_a,-rij,pb_block_a,pa_block_a,&
                                          ksb_block_a,ksa_block_a,pb_block_b,pa_block_b,&
                                          ksb_block_b,ksa_block_b,anag,do_ewald,taper,error)
                           ENDIF
                        END IF
                        IF(calculate_forces) THEN
                           force_ab = 0._dp
                           CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                           IF ( nspins == 1 ) THEN
                              IF ( .NOT. switch ) THEN
                                 CALL dfock2c(se_kind_a,se_kind_b,rij,pa_block_a,pb_block_a,&
                                              force_ab,delta,anag,do_ewald,taper,error)
                              ELSE
                                 CALL dfock2c(se_kind_b,se_kind_a,-rij,pb_block_a,pa_block_a,&
                                              force_ab,delta,anag,do_ewald,taper,error)
                              ENDIF
                           ELSE IF ( nspins == 2 ) THEN
                              IF ( .NOT. switch ) THEN
                                 CALL dfock2c(se_kind_a,se_kind_b,rij,pa_block_a,pb_block_a,&
                                              pa_block_b,pb_block_b,force_ab,delta,anag,do_ewald,&
                                              taper,error)
                              ELSE
                                 CALL dfock2c(se_kind_b,se_kind_a,-rij,pb_block_a,pa_block_a,&
                                              pb_block_b,pa_block_b,force_ab,delta,anag,do_ewald,&
                                              taper,error)
                              ENDIF
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
         CALL deallocate_dbm_from_sm(diagmat_p(ispin)%matrix,&
                                     matrix_p(1)%matrix,.FALSE.,error=error)
         CALL deallocate_dbm_from_sm(diagmat_ks(ispin)%matrix,&
                                     ks_matrix(ispin)%matrix,.TRUE.,para_env,error=error)
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
         CALL cp_sm_sm_trace(ks_matrix(ispin)%matrix,&
                             matrix_p(ispin)%matrix,trace=ecoul,&
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
    LOGICAL                                  :: anag, defined, failure
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
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(taper_type), POINTER                :: taper

    failure=.FALSE.
    NULLIFY(matrix_h,matrix_p,dft_control,cell,force,particle_set,se_control)
    NULLIFY(taper)
    CALL timeset(routineN,"I","",handle)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,matrix_h=matrix_h,&
         cell=cell,taper=taper,para_env=para_env,error=error)

    ! set values for tapering function
    se_control => dft_control%qs_control%se_control
    anag       =  se_control%analytical_gradients

    nspins=dft_control%nspins
    CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)

    IF ( .NOT. failure ) THEN

       matrix_p => rho%rho_ao

       CALL get_qs_env(qs_env=qs_env,&
                       sab_orb=sab_orb,&
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

       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a)

          IF (.NOT.defined .OR. natorb_a < 1) CYCLE

          DO jkind=1,nkind
             atomic_kind => atomic_kind_set(jkind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                  se_parameter=se_kind_b)
             CALL get_se_param(se_kind_b,&
                               defined=defined,&
                               natorb=natorb_b)

             IF (.NOT.defined .OR. natorb_b < 1) CYCLE

             iab = ikind + nkind*(jkind - 1)
             IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE
             CALL get_neighbor_list_set(neighbor_list_set=&
                                        sab_orb(iab)%neighbor_list_set,&
                                        nlist=nlist)
             NULLIFY(sab_orb_neighbor_list)

             DO ilist=1,nlist
                IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
                   sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
                ELSE
                   sab_orb_neighbor_list => next(sab_orb_neighbor_list)
                END IF
                sab_orb_neighbor_list_local => sab_orb_neighbor_list
                CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                       atom=iatom,nnode=nnode)

                sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

                DO inode=1,nnode
                   CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,neighbor=jatom,r=rij)
                   IF (iatom <= jatom) THEN
                      irow = iatom
                      icol = jatom
                   ELSE
                      irow = jatom
                      icol = iatom
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
                      CPPostcondition(SIZE(p_block_a,1)==SIZE(p_block_b,1),cp_failure_level,routineP,error,failure)
                      CPPostcondition(SIZE(p_block_a,2)==SIZE(p_block_b,2),cp_failure_level,routineP,error,failure)
                      p_block_tot(1:SIZE(p_block_a,1),1:SIZE(p_block_a,2)) = p_block_a + p_block_b
                   END IF

                   dr = SUM(rij(:)**2)
                   IF ( iatom == jatom .AND. dr < 0.001_dp ) THEN
                      CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                      IF ( nspins == 1 ) THEN
                         CALL fock1(se_kind_a,p_block_tot,p_block_a,ks_block_a,factor=0.5_dp,error=error)
                      ELSE IF ( nspins == 2 ) THEN
                         CALL fock1(se_kind_a,p_block_tot,p_block_a,ks_block_a,factor=1.0_dp,error=error)
                         CALL fock1(se_kind_a,p_block_tot,p_block_b,ks_block_b,factor=1.0_dp,error=error)
                      END IF
                   ELSE
                      CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                      IF ( nspins == 1 ) THEN
                         IF ( irow == iatom ) THEN
                            CALL fock2e(se_kind_a,se_kind_b,rij,p_block_a,ks_block_a,anag,taper,error)
                         ELSE
                            CALL fock2e(se_kind_b,se_kind_a,-rij,p_block_a,ks_block_a,anag,taper,error)
                         ENDIF
                      ELSE IF ( nspins == 2 ) THEN
                         IF ( irow == iatom ) THEN
                            CALL fock2e(se_kind_a,se_kind_b,rij,p_block_a,ks_block_a,&
                                        p_block_b,ks_block_b,anag,taper,error)
                         ELSE
                            CALL fock2e(se_kind_b,se_kind_a,-rij,p_block_a,ks_block_a,&
                                        p_block_b,ks_block_b,anag,taper,error)
                         ENDIF
                      END IF

                      IF(calculate_forces) THEN
                         force_ab(1:3) = 0._dp
                         CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                         IF ( nspins == 1 ) THEN
                            IF ( irow == iatom ) THEN
                               CALL dfock2e(se_kind_a,se_kind_b,rij,p_block_a,force_ab,delta,anag,taper,error)
                            ELSE
                               CALL dfock2e(se_kind_b,se_kind_a,-rij,p_block_a,force_ab,delta,anag,taper,error)
                            ENDIF
                         ELSE IF ( nspins == 2 ) THEN
                            IF ( irow == iatom ) THEN
                               CALL dfock2e(se_kind_a,se_kind_b,rij,p_block_a,p_block_b,force_ab,&
                                    delta,anag,taper,error)
                            ELSE
                               CALL dfock2e(se_kind_b,se_kind_a,-rij,p_block_a,p_block_b,force_ab,&
                                    delta,anag,taper,error)
                            ENDIF
                         END IF
                         atom_a = atom_of_kind(iatom)
                         atom_b = atom_of_kind(jatom)
                         IF ( irow /= iatom ) force_ab = -force_ab
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
!> \brief Construction of 1-center Fock Matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock1(sep, p_tot, p_mat, f_mat, factor, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: p_tot, p_mat
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ijw, ikw, im, ip, j, jlw, &
                                                jm, jp, k, klw, l
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: sum

    failure = .FALSE.
    !   One-center coulomb and exchange terms for semiempirical_type sep
    !
    !  F(i,j)=F(i,j)+sum(k,l)((PA(k,l)+PB(k,l))*<i,j|k,l>
    !                        -(PA(k,l)        )*<i,k|j,l>), k,l on type sep.
    !
    DO i = 1, sep%natorb
       DO j = 1, i
          !
          !    `J' Address IJ in W
          !
          ijw = (i*(i-1))/2 + j
          sum = 0.0_dp
          DO k = 1, sep%natorb
             DO l =  1, sep%natorb
                ip = MAX(k,l)
                jp = MIN(k,l)
                !
                !    `J' Address KL in W
                !
                im = MAX(k,l)
                jm = MIN(k,l)
                klw = (im*(im-1))/2 + jm
                !
                !    `K' Address IK in W
                !
                im = MAX(k,j)
                jm = MIN(k,j)
                ikw = (im*(im-1))/2 + jm
                !
                !    `K' Address JL in W
                !
                im = MAX(l,i)
                jm = MIN(l,i)
                jlw = (im*(im-1))/2 + jm

                sum = sum + p_tot(ip, jp) * sep%w(ijw, klw) - p_mat(ip,jp) * sep%w(ikw, jlw)
             END DO
          END DO
          f_mat(i,j) = f_mat(i,j) + factor*sum
          f_mat(j,i) = f_mat(i,j)
       END DO
    END DO
  END SUBROUTINE fock1

! *****************************************************************************
  SUBROUTINE rfock2c(sepa,sepb,rij,pa,pb,fa,fb,anag,shortrange,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pa, pb
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fa, fb
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rfock2c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(2025)           :: wint
    REAL(KIND=dp), DIMENSION(45)             :: ca, cb, ppa, ppb

    na = SIZE ( pa, 1 )
    nb = SIZE ( pb, 1 )

    CALL rotint (sepa,sepb,rij,wint,anag=anag,shortrange=shortrange,taper=taper,error=error)
    IF ( na==1 .AND. nb==1 ) THEN
      fa(1,1) = fa(1,1) + pb(1,1)*wint(1)
      fb(1,1) = fb(1,1) + pa(1,1)*wint(1)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      CALL contract( pb, ppb )
      fa(1,1) = fa(1,1) + SUM ( ppb(1:10)*wint(1:10) )
      ppa(1:10) = pa(1,1)*wint(1:10)
      CALL expand( fb, ppa )
    ELSE IF ( na==4 .AND. nb==1 ) THEN
      CALL contract( pa, ppa )
      fb(1,1) = fb(1,1) + SUM ( ppa(1:10)*wint(1:10) )
      ppb(1:10) = pb(1,1)*wint(1:10)
      CALL expand( fa, ppb )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL contract( pa, ppa )
      CALL contract( pb, ppb )
      DO i=1,10
        ii = (i-1)*10 + 1
        ca(i) = SUM ( wint(ii:ii+9)*ppb(1:10) )
        cb(i) = SUM ( wint(i:90+i:10)*ppa(1:10) )
      END DO
      CALL expand( fa, ca )
      CALL expand( fb, cb )
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE rfock2c

! *****************************************************************************
  SUBROUTINE ufock2c(sepa,sepb,rij,pa_a,pb_a,fa_a,fb_a,pa_b,pb_b,fa_b,fb_b,anag,&
       shortrange,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pa_a, pb_a
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fa_a, fb_a
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pa_b, pb_b
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fa_b, fb_b
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ufock2c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ess
    REAL(KIND=dp), DIMENSION(2025)           :: wint
    REAL(KIND=dp), DIMENSION(45)             :: ca, cb, ppa_a, ppa_b, ppb_a, &
                                                ppb_b

    na = SIZE ( pa_a, 1 )
    nb = SIZE ( pb_a, 1 )

    CALL rotint (sepa,sepb,rij,wint,anag=anag,shortrange=shortrange,taper=taper,error=error)
    IF ( na==1 .AND. nb==1 ) THEN
      fa_a(1,1) = fa_a(1,1) + (pb_a(1,1)+pb_b(1,1))*wint(1)
      fb_a(1,1) = fb_a(1,1) + (pa_a(1,1)+pa_b(1,1))*wint(1)
      fa_b(1,1) = fa_b(1,1) + (pb_a(1,1)+pb_b(1,1))*wint(1)
      fb_b(1,1) = fb_b(1,1) + (pa_a(1,1)+pa_b(1,1))*wint(1)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      CALL contract( pb_a, ppb_a )
      CALL contract( pb_b, ppb_b )
      ess = SUM ( (ppb_a(1:10)+ppb_b(1:10)) * wint(1:10) )
      fa_a(1,1) = fa_a(1,1) + ess
      fa_b(1,1) = fa_b(1,1) + ess
      ppa_a(1:10) = (pa_a(1,1)+pa_b(1,1))*wint(1:10)
      CALL expand( fb_a, ppa_a )
      CALL expand( fb_b, ppa_a )
    ELSE IF ( na==4 .AND. nb==1 ) THEN
      CALL contract( pa_a, ppa_a )
      CALL contract( pa_b, ppa_b )
      ess = SUM ( (ppa_a(1:10)+ppa_b(1:10)) * wint(1:10) )
      fb_a(1,1) = fb_a(1,1) + ess
      fb_b(1,1) = fb_b(1,1) + ess
      ppb_a(1:10) = (pb_a(1,1)+pb_b(1,1))*wint(1:10)
      CALL expand( fa_a, ppb_a )
      CALL expand( fa_b, ppb_a )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL contract( pa_a, ppa_a )
      CALL contract( pa_b, ppa_b )
      CALL contract( pb_a, ppb_a )
      CALL contract( pb_b, ppb_b )
      DO i=1,10
        ii = (i-1)*10 + 1
        ca(i) = SUM ( wint(ii:ii+9)*(ppb_a(1:10)+ppb_b(1:10)) )
        cb(i) = SUM ( wint(i:90+i:10)*(ppa_a(1:10)+ppa_b(1:10)) )
      END DO
      CALL expand( fa_a, ca )
      CALL expand( fa_b, ca )
      CALL expand( fb_a, cb )
      CALL expand( fb_b, cb )
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE ufock2c

! *****************************************************************************
  SUBROUTINE rdfock2c(sepa,sepb,rij,pa,pb,force,delta,anag,shortrange,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pa, pb
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: force
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rdfock2c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, j, na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(2025, 3)        :: dwint
    REAL(KIND=dp), DIMENSION(45)             :: ca, ppa, ppb

    na = SIZE ( pa, 1 )
    nb = SIZE ( pb, 1 )

    CALL drotint (sepa,sepb,rij,dwint,delta,anag=anag,shortrange=shortrange,taper=taper,error=error)
    IF ( na==1 .AND. nb==1 ) THEN
      force(:) = force(:) + pa(1,1)*pb(1,1)*dwint(1,:)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      CALL contract( pb, ppb )
      DO i=1,3
        force(i) = force(i) + pa(1,1) * SUM ( ppb(1:10)*dwint(1:10,i) )
      END DO
   ELSE IF ( na==4 .AND. nb==1 ) THEN
      CALL contract( pa, ppa )
      DO i=1,3
        force(i) = force(i) + pb(1,1) * SUM ( ppa(1:10)*dwint(1:10,i) )
      END DO
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL contract( pb, ppb )
      CALL contract( pa, ppa )
      DO j=1,3
        DO i=1,10
          ii = (i-1)*10 + 1
          ca(i) = SUM ( dwint(ii:ii+9,j)*ppb(1:10) )
        END DO
        force(j) = force(j) + SUM ( ca(1:10)*ppa(1:10) )
      END DO
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE rdfock2c

! *****************************************************************************
  SUBROUTINE udfock2c(sepa,sepb,rij,pa_a,pb_a,pa_b,pb_b,force,delta,anag,shortrange,&
       taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pa_a, pb_a, pa_b, pb_b
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: force
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'udfock2c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, j, na, nb
    LOGICAL                                  :: failure
    REAL(dp)                                 :: ess, pta, ptb
    REAL(KIND=dp), DIMENSION(2025, 3)        :: dwint
    REAL(KIND=dp), DIMENSION(45)             :: ca, ppa_a, ppa_b, ppb_a, ppb_b

    na = SIZE ( pa_a, 1 )
    nb = SIZE ( pb_a, 1 )

    CALL drotint (sepa,sepb,rij,dwint,delta,anag=anag,shortrange=shortrange,taper=taper,error=error)
    IF ( na==1 .AND. nb==1 ) THEN
      pta = pa_a(1,1)+pa_b(1,1)
      ptb = pb_a(1,1)+pb_b(1,1)
      force(:) = force(:) + pta*ptb*dwint(1,:)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      CALL contract( pb_a, ppb_a )
      CALL contract( pb_b, ppb_b )
      DO i=1,3
        ess = SUM ( (ppb_a(1:10)+ppb_b(1:10)) * dwint(1:10,i) )
        force(i) = force(i) + (pa_a(1,1)+pa_b(1,1)) * ess
      END DO
   ELSE IF ( na==4 .AND. nb==1 ) THEN
      CALL contract( pa_a, ppa_a )
      CALL contract( pa_b, ppa_b )
      DO i=1,3
        ess = SUM ( (ppa_a(1:10)+ppa_b(1:10)) * dwint(1:10,i) )
        force(i) = force(i) + (pb_a(1,1)+pb_b(1,1)) * ess
      END DO
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL contract( pa_a, ppa_a )
      CALL contract( pa_b, ppa_b )
      CALL contract( pb_a, ppb_a )
      CALL contract( pb_b, ppb_b )
      DO j=1,3
        DO i=1,10
          ii = (i-1)*10 + 1
          ca(i) = SUM ( dwint(ii:ii+9,j)*(ppb_a(1:10)+ppb_b(1:10)) )
        END DO
        force(j) = force(j) + SUM ( ca(1:10)*(ppa_a(1:10)+ppa_b(1:10)) )
      END DO
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE udfock2c

! *****************************************************************************
  SUBROUTINE rfock2e(sepa,sepb,rij,pab,fab,anag,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fab
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rfock2e', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(2025)           :: wint
    REAL(KIND=dp), DIMENSION(9, 9)           :: kab

    na = SIZE ( pab, 1 )
    nb = SIZE ( pab, 2 )

    CALL rotint (sepa,sepb,rij,wint,anag=anag,shortrange=.FALSE.,taper=taper,error=error)
    IF ( na==1 .AND. nb==1 ) THEN
      fab(1,1) = fab(1,1) - 0.5_dp * pab(1,1)*wint(1)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      fab(1,1) = fab(1,1) - 0.5_dp*( wint(1)*pab(1,1) +&
        wint(2)*pab(1,4) + wint(4)*pab(1,2) + wint(7)*pab(1,3) )
      fab(1,2) = fab(1,2) - 0.5_dp*( wint(4)*pab(1,1) +&
        wint(5)*pab(1,4) + wint(6)*pab(1,2) + wint(9)*pab(1,3) )
      fab(1,3) = fab(1,3) - 0.5_dp*( wint(7)*pab(1,1) +&
        wint(8)*pab(1,4) + wint(9)*pab(1,2) + wint(10)*pab(1,3) )
      fab(1,4) = fab(1,4) - 0.5_dp*( wint(2)*pab(1,1) +&
        wint(3)*pab(1,4) + wint(5)*pab(1,2) + wint(8)*pab(1,3) )
    ELSE IF ( na==4 .AND. nb==1 ) THEN
      fab(1,1) = fab(1,1) - 0.5_dp*( wint(1)*pab(1,1) +&
        wint(2)*pab(4,1) + wint(4)*pab(2,1) + wint(7)*pab(3,1) )
      fab(2,1) = fab(2,1) - 0.5_dp*( wint(4)*pab(1,1) +&
        wint(5)*pab(4,1) + wint(6)*pab(2,1) + wint(9)*pab(3,1) )
      fab(3,1) = fab(3,1) - 0.5_dp*( wint(7)*pab(1,1) +&
        wint(8)*pab(4,1) + wint(9)*pab(2,1) + wint(10)*pab(3,1) )
      fab(4,1) = fab(4,1) - 0.5_dp*( wint(2)*pab(1,1) +&
        wint(3)*pab(4,1) + wint(5)*pab(2,1) + wint(8)*pab(3,1) )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL exchange_sp ( kab, pab, wint )
      fab(1:4,1:4) = fab(1:4,1:4) - 0.5_dp*kab(1:4,1:4)
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE rfock2e

! *****************************************************************************
  SUBROUTINE ufock2e(sepa,sepb,rij,pab_a,fab_a,pab_b,fab_b,anag,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab_a
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fab_a
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab_b
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fab_b
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ufock2e', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(2025)           :: wint
    REAL(KIND=dp), DIMENSION(9, 9)           :: kab

    na = SIZE ( pab_a, 1 )
    nb = SIZE ( pab_a, 2 )

    CALL rotint (sepa,sepb,rij,wint,anag=anag,shortrange=.FALSE.,taper=taper,error=error)
    IF ( na==1 .AND. nb==1 ) THEN
      fab_a(1,1) = fab_a(1,1) - pab_a(1,1)*wint(1)
      fab_b(1,1) = fab_b(1,1) - pab_b(1,1)*wint(1)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      fab_a(1,1) = fab_a(1,1) - wint(1)*pab_a(1,1) -&
        wint(2)*pab_a(1,4) - wint(4)*pab_a(1,2) - wint(7)*pab_a(1,3)
      fab_a(1,2) = fab_a(1,2) - wint(4)*pab_a(1,1) -&
        wint(5)*pab_a(1,4) - wint(6)*pab_a(1,2) - wint(9)*pab_a(1,3)
      fab_a(1,3) = fab_a(1,3) - wint(7)*pab_a(1,1) -&
        wint(8)*pab_a(1,4) - wint(9)*pab_a(1,2) - wint(10)*pab_a(1,3)
      fab_a(1,4) = fab_a(1,4) - wint(2)*pab_a(1,1) -&
        wint(3)*pab_a(1,4) - wint(5)*pab_a(1,2) - wint(8)*pab_a(1,3)
      fab_b(1,1) = fab_b(1,1) - wint(1)*pab_b(1,1) -&
        wint(2)*pab_b(1,4) - wint(4)*pab_b(1,2) - wint(7)*pab_b(1,3)
      fab_b(1,2) = fab_b(1,2) - wint(4)*pab_b(1,1) -&
        wint(5)*pab_b(1,4) - wint(6)*pab_b(1,2) - wint(9)*pab_b(1,3)
      fab_b(1,3) = fab_b(1,3) - wint(7)*pab_b(1,1) -&
        wint(8)*pab_b(1,4) - wint(9)*pab_b(1,2) - wint(10)*pab_b(1,3)
      fab_b(1,4) = fab_b(1,4) - wint(2)*pab_b(1,1) -&
        wint(3)*pab_b(1,4) - wint(5)*pab_b(1,2) - wint(8)*pab_b(1,3)
    ELSE IF ( na==4 .AND. nb==1 ) THEN
      fab_a(1,1) = fab_a(1,1) - wint(1)*pab_a(1,1) -&
        wint(2)*pab_a(4,1) - wint(4)*pab_a(2,1) - wint(7)*pab_a(3,1)
      fab_a(2,1) = fab_a(2,1) - wint(4)*pab_a(1,1) -&
        wint(5)*pab_a(4,1) - wint(6)*pab_a(2,1) - wint(9)*pab_a(3,1)
      fab_a(3,1) = fab_a(3,1) - wint(7)*pab_a(1,1) -&
        wint(8)*pab_a(4,1) - wint(9)*pab_a(2,1) - wint(10)*pab_a(3,1)
      fab_a(4,1) = fab_a(4,1) - wint(2)*pab_a(1,1) -&
        wint(3)*pab_a(4,1) - wint(5)*pab_a(2,1) - wint(8)*pab_a(3,1)
      fab_b(1,1) = fab_b(1,1) - wint(1)*pab_b(1,1) -&
        wint(2)*pab_b(4,1) - wint(4)*pab_b(2,1) - wint(7)*pab_b(3,1)
      fab_b(2,1) = fab_b(2,1) - wint(4)*pab_b(1,1) -&
        wint(5)*pab_b(4,1) - wint(6)*pab_b(2,1) - wint(9)*pab_b(3,1)
      fab_b(3,1) = fab_b(3,1) - wint(7)*pab_b(1,1) -&
        wint(8)*pab_b(4,1) - wint(9)*pab_b(2,1) - wint(10)*pab_b(3,1)
      fab_b(4,1) = fab_b(4,1) - wint(2)*pab_b(1,1) -&
        wint(3)*pab_b(4,1) - wint(5)*pab_b(2,1) - wint(8)*pab_b(3,1)
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL exchange_sp ( kab, pab_a, wint )
      fab_a(1:4,1:4) = fab_a(1:4,1:4) - kab(1:4,1:4)
      CALL exchange_sp ( kab, pab_b, wint )
      fab_b(1:4,1:4) = fab_b(1:4,1:4) - kab(1:4,1:4)
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE ufock2e

! *****************************************************************************
  SUBROUTINE rdfock2e(sepa,sepb,rij,pab,force,delta,anag,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: force
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rdfock2e', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(2025, 3)        :: dwint
    REAL(KIND=dp), DIMENSION(45, 45)         :: kab

    na = SIZE ( pab, 1 )
    nb = SIZE ( pab, 2 )

    CALL drotint (sepa,sepb,rij,dwint,delta,anag=anag,shortrange=.FALSE.,taper=taper,error=error)
    IF ( na==1 .AND. nb==1 ) THEN
      force(:) = force(:) - 0.5_dp*pab(1,1)*pab(1,1)*dwint(1,:)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      force(:) = force(:) - 0.5_dp*pab(1,1)*( dwint(1,:)*pab(1,1) +&
        dwint(2,:)*pab(1,4) + dwint(4,:)*pab(1,2) + dwint(7,:)*pab(1,3) )
      force(:) = force(:) - 0.5_dp*pab(1,2)*( dwint(4,:)*pab(1,1) +&
        dwint(5,:)*pab(1,4) + dwint(6,:)*pab(1,2) + dwint(9,:)*pab(1,3) )
      force(:) = force(:) - 0.5_dp*pab(1,3)*( dwint(7,:)*pab(1,1) +&
        dwint(8,:)*pab(1,4) + dwint(9,:)*pab(1,2) + dwint(10,:)*pab(1,3) )
      force(:) = force(:) - 0.5_dp*pab(1,4)*( dwint(2,:)*pab(1,1) +&
        dwint(3,:)*pab(1,4) + dwint(5,:)*pab(1,2) + dwint(8,:)*pab(1,3) )
   ELSE IF ( na==4 .AND. nb==1 ) THEN
      force(:) = force(:) - 0.5_dp*pab(1,1)*( dwint(1,:)*pab(1,1) +&
        dwint(2,:)*pab(4,1) + dwint(4,:)*pab(2,1) + dwint(7,:)*pab(3,1) )
      force(:) = force(:) - 0.5_dp*pab(2,1)*( dwint(4,:)*pab(1,1) +&
        dwint(5,:)*pab(4,1) + dwint(6,:)*pab(2,1) + dwint(9,:)*pab(3,1) )
      force(:) = force(:) - 0.5_dp*pab(3,1)*( dwint(7,:)*pab(1,1) +&
        dwint(8,:)*pab(4,1) + dwint(9,:)*pab(2,1) + dwint(10,:)*pab(3,1) )
      force(:) = force(:) - 0.5_dp*pab(4,1)*( dwint(2,:)*pab(1,1) +&
        dwint(3,:)*pab(4,1) + dwint(5,:)*pab(2,1) + dwint(8,:)*pab(3,1) )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      DO j=1,3
        CALL exchange_sp ( kab, pab, dwint(:,j) )
        force(j) = force(j) - 0.5_dp * SUM(pab*kab(1:4,1:4))
      END DO
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE rdfock2e

! *****************************************************************************
  SUBROUTINE udfock2e(sepa,sepb,rij,pab_a,pab_b,force,delta,anag,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab_a, pab_b
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: force
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'udfock2e', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(2025, 3)        :: dwint
    REAL(KIND=dp), DIMENSION(45, 45)         :: kab

    na = SIZE ( pab_a, 1 )
    nb = SIZE ( pab_a, 2 )

    CALL drotint (sepa,sepb,rij,dwint,delta,anag=anag,shortrange=.FALSE.,taper=taper,error=error)
    IF ( na==1 .AND. nb==1 ) THEN
      force(:) = force(:) - pab_a(1,1)*pab_a(1,1)*dwint(1,:)
      force(:) = force(:) - pab_b(1,1)*pab_b(1,1)*dwint(1,:)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      force(:) = force(:) - pab_a(1,1)*( dwint(1,:)*pab_a(1,1) +&
        dwint(2,:)*pab_a(1,4) + dwint(4,:)*pab_a(1,2) + dwint(7,:)*pab_a(1,3) )
      force(:) = force(:) - pab_a(1,2)*( dwint(4,:)*pab_a(1,1) +&
        dwint(5,:)*pab_a(1,4) + dwint(6,:)*pab_a(1,2) + dwint(9,:)*pab_a(1,3) )
      force(:) = force(:) - pab_a(1,3)*( dwint(7,:)*pab_a(1,1) +&
        dwint(8,:)*pab_a(1,4) + dwint(9,:)*pab_a(1,2) + dwint(10,:)*pab_a(1,3) )
      force(:) = force(:) - pab_a(1,4)*( dwint(2,:)*pab_a(1,1) +&
        dwint(3,:)*pab_a(1,4) + dwint(5,:)*pab_a(1,2) + dwint(8,:)*pab_a(1,3) )
      force(:) = force(:) - pab_b(1,1)*( dwint(1,:)*pab_b(1,1) +&
        dwint(2,:)*pab_b(1,4) + dwint(4,:)*pab_b(1,2) + dwint(7,:)*pab_b(1,3) )
      force(:) = force(:) - pab_b(1,2)*( dwint(4,:)*pab_b(1,1) +&
        dwint(5,:)*pab_b(1,4) + dwint(6,:)*pab_b(1,2) + dwint(9,:)*pab_b(1,3) )
      force(:) = force(:) - pab_b(1,3)*( dwint(7,:)*pab_b(1,1) +&
        dwint(8,:)*pab_b(1,4) + dwint(9,:)*pab_b(1,2) + dwint(10,:)*pab_b(1,3) )
      force(:) = force(:) - pab_b(1,4)*( dwint(2,:)*pab_b(1,1) +&
        dwint(3,:)*pab_b(1,4) + dwint(5,:)*pab_b(1,2) + dwint(8,:)*pab_b(1,3) )
   ELSE IF ( na==4 .AND. nb==1 ) THEN
      force(:) = force(:) - pab_a(1,1)*( dwint(1,:)*pab_a(1,1) +&
        dwint(2,:)*pab_a(4,1) + dwint(4,:)*pab_a(2,1) + dwint(7,:)*pab_a(3,1) )
      force(:) = force(:) - pab_a(2,1)*( dwint(4,:)*pab_a(1,1) +&
        dwint(5,:)*pab_a(4,1) + dwint(6,:)*pab_a(2,1) + dwint(9,:)*pab_a(3,1) )
      force(:) = force(:) - pab_a(3,1)*( dwint(7,:)*pab_a(1,1) +&
        dwint(8,:)*pab_a(4,1) + dwint(9,:)*pab_a(2,1) + dwint(10,:)*pab_a(3,1) )
      force(:) = force(:) - pab_a(4,1)*( dwint(2,:)*pab_a(1,1) +&
        dwint(3,:)*pab_a(4,1) + dwint(5,:)*pab_a(2,1) + dwint(8,:)*pab_a(3,1) )
      force(:) = force(:) - pab_b(1,1)*( dwint(1,:)*pab_b(1,1) +&
        dwint(2,:)*pab_b(4,1) + dwint(4,:)*pab_b(2,1) + dwint(7,:)*pab_b(3,1) )
      force(:) = force(:) - pab_b(2,1)*( dwint(4,:)*pab_b(1,1) +&
        dwint(5,:)*pab_b(4,1) + dwint(6,:)*pab_b(2,1) + dwint(9,:)*pab_b(3,1) )
      force(:) = force(:) - pab_b(3,1)*( dwint(7,:)*pab_b(1,1) +&
        dwint(8,:)*pab_b(4,1) + dwint(9,:)*pab_b(2,1) + dwint(10,:)*pab_b(3,1) )
      force(:) = force(:) - pab_b(4,1)*( dwint(2,:)*pab_b(1,1) +&
        dwint(3,:)*pab_b(4,1) + dwint(5,:)*pab_b(2,1) + dwint(8,:)*pab_b(3,1) )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      DO j=1,3
        CALL exchange_sp ( kab, pab_a, dwint(:,j) )
        force(j) = force(j) - SUM(pab_a*kab(1:4,1:4))
        CALL exchange_sp ( kab, pab_b, dwint(:,j) )
        force(j) = force(j) - SUM(pab_b*kab(1:4,1:4))
      END DO
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE udfock2e

! *****************************************************************************
  SUBROUTINE exchange_sp ( kab, pab, w )
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: kab
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab
    REAL(dp), DIMENSION(:), INTENT(IN)       :: w

    REAL(dp), DIMENSION(16)                  :: pk

    pk( 1) = pab(1,1)
    pk( 2) = pab(1,4)
    pk( 3) = pab(1,2)
    pk( 4) = pab(1,3)
    pk( 5) = pab(4,1)
    pk( 6) = pab(4,4)
    pk( 7) = pab(4,2)
    pk( 8) = pab(4,3)
    pk( 9) = pab(2,1)
    pk(10) = pab(2,4)
    pk(11) = pab(2,2)
    pk(12) = pab(2,3)
    pk(13) = pab(3,1)
    pk(14) = pab(3,4)
    pk(15) = pab(3,2)
    pk(16) = pab(3,3)

    kab(1,1)= &
      +pk( 1)*w(  1)+pk( 2)*w(  2)+pk( 3)*w(  4)+pk( 4)*w(  7) &
      +pk( 5)*w( 11)+pk( 6)*w( 12)+pk( 7)*w( 14)+pk( 8)*w( 17) &
      +pk( 9)*w( 31)+pk(10)*w( 32)+pk(11)*w( 34)+pk(12)*w( 37) &
      +pk(13)*w( 61)+pk(14)*w( 62)+pk(15)*w( 64)+pk(16)*w( 67)
    kab(1,4)= &
      +pk( 1)*w(  2)+pk( 2)*w(  3)+pk( 3)*w(  5)+pk( 4)*w(  8) &
      +pk( 5)*w( 12)+pk( 6)*w( 13)+pk( 7)*w( 15)+pk( 8)*w( 18) &
      +pk( 9)*w( 32)+pk(10)*w( 33)+pk(11)*w( 35)+pk(12)*w( 38) &
      +pk(13)*w( 62)+pk(14)*w( 63)+pk(15)*w( 65)+pk(16)*w( 68)
    kab(1,2)= &
      +pk( 1)*w(  4)+pk( 2)*w(  5)+pk( 3)*w(  6)+pk( 4)*w(  9) &
      +pk( 5)*w( 14)+pk( 6)*w( 15)+pk( 7)*w( 16)+pk( 8)*w( 19) &
      +pk( 9)*w( 34)+pk(10)*w( 35)+pk(11)*w( 36)+pk(12)*w( 39) &
      +pk(13)*w( 64)+pk(14)*w( 65)+pk(15)*w( 66)+pk(16)*w( 69)
    kab(1,3)= &
      +pk( 1)*w(  7)+pk( 2)*w(  8)+pk( 3)*w(  9)+pk( 4)*w( 10) &
      +pk( 5)*w( 17)+pk( 6)*w( 18)+pk( 7)*w( 19)+pk( 8)*w( 20) &
      +pk( 9)*w( 37)+pk(10)*w( 38)+pk(11)*w( 39)+pk(12)*w( 40) &
      +pk(13)*w( 67)+pk(14)*w( 68)+pk(15)*w( 69)+pk(16)*w( 70)
    kab(4,1)= &
      +pk( 1)*w( 11)+pk( 2)*w( 12)+pk( 3)*w( 14)+pk( 4)*w( 17) &
      +pk( 5)*w( 21)+pk( 6)*w( 22)+pk( 7)*w( 24)+pk( 8)*w( 27) &
      +pk( 9)*w( 41)+pk(10)*w( 42)+pk(11)*w( 44)+pk(12)*w( 47) &
      +pk(13)*w( 71)+pk(14)*w( 72)+pk(15)*w( 74)+pk(16)*w( 77)
    kab(4,4)= &
      +pk( 1)*w( 12)+pk( 2)*w( 13)+pk( 3)*w( 15)+pk( 4)*w( 18) &
      +pk( 5)*w( 22)+pk( 6)*w( 23)+pk( 7)*w( 25)+pk( 8)*w( 28) &
      +pk( 9)*w( 42)+pk(10)*w( 43)+pk(11)*w( 45)+pk(12)*w( 48) &
      +pk(13)*w( 72)+pk(14)*w( 73)+pk(15)*w( 75)+pk(16)*w( 78)
    kab(4,2)= &
      +pk( 1)*w( 14)+pk( 2)*w( 15)+pk( 3)*w( 16)+pk( 4)*w( 19) &
      +pk( 5)*w( 24)+pk( 6)*w( 25)+pk( 7)*w( 26)+pk( 8)*w( 29) &
      +pk( 9)*w( 44)+pk(10)*w( 45)+pk(11)*w( 46)+pk(12)*w( 49) &
      +pk(13)*w( 74)+pk(14)*w( 75)+pk(15)*w( 76)+pk(16)*w( 79)
    kab(4,3)= &
      +pk( 1)*w( 17)+pk( 2)*w( 18)+pk( 3)*w( 19)+pk( 4)*w( 20) &
      +pk( 5)*w( 27)+pk( 6)*w( 28)+pk( 7)*w( 29)+pk( 8)*w( 30) &
      +pk( 9)*w( 47)+pk(10)*w( 48)+pk(11)*w( 49)+pk(12)*w( 50) &
      +pk(13)*w( 77)+pk(14)*w( 78)+pk(15)*w( 79)+pk(16)*w( 80)
    kab(2,1)= &
      +pk( 1)*w( 31)+pk( 2)*w( 32)+pk( 3)*w( 34)+pk( 4)*w( 37) &
      +pk( 5)*w( 41)+pk( 6)*w( 42)+pk( 7)*w( 44)+pk( 8)*w( 47) &
      +pk( 9)*w( 51)+pk(10)*w( 52)+pk(11)*w( 54)+pk(12)*w( 57) &
      +pk(13)*w( 81)+pk(14)*w( 82)+pk(15)*w( 84)+pk(16)*w( 87)
    kab(2,4)= &
      +pk( 1)*w( 32)+pk( 2)*w( 33)+pk( 3)*w( 35)+pk( 4)*w( 38) &
      +pk( 5)*w( 42)+pk( 6)*w( 43)+pk( 7)*w( 45)+pk( 8)*w( 48) &
      +pk( 9)*w( 52)+pk(10)*w( 53)+pk(11)*w( 55)+pk(12)*w( 58) &
      +pk(13)*w( 82)+pk(14)*w( 83)+pk(15)*w( 85)+pk(16)*w( 88)
    kab(2,2)= &
      +pk( 1)*w( 34)+pk( 2)*w( 35)+pk( 3)*w( 36)+pk( 4)*w( 39) &
      +pk( 5)*w( 44)+pk( 6)*w( 45)+pk( 7)*w( 46)+pk( 8)*w( 49) &
      +pk( 9)*w( 54)+pk(10)*w( 55)+pk(11)*w( 56)+pk(12)*w( 59) &
      +pk(13)*w( 84)+pk(14)*w( 85)+pk(15)*w( 86)+pk(16)*w( 89)
    kab(2,3)= &
      +pk( 1)*w( 37)+pk( 2)*w( 38)+pk( 3)*w( 39)+pk( 4)*w( 40) &
      +pk( 5)*w( 47)+pk( 6)*w( 48)+pk( 7)*w( 49)+pk( 8)*w( 50) &
      +pk( 9)*w( 57)+pk(10)*w( 58)+pk(11)*w( 59)+pk(12)*w( 60) &
      +pk(13)*w( 87)+pk(14)*w( 88)+pk(15)*w( 89)+pk(16)*w( 90)
    kab(3,1)= &
      +pk( 1)*w( 61)+pk( 2)*w( 62)+pk( 3)*w( 64)+pk( 4)*w( 67) &
      +pk( 5)*w( 71)+pk( 6)*w( 72)+pk( 7)*w( 74)+pk( 8)*w( 77) &
      +pk( 9)*w( 81)+pk(10)*w( 82)+pk(11)*w( 84)+pk(12)*w( 87) &
      +pk(13)*w( 91)+pk(14)*w( 92)+pk(15)*w( 94)+pk(16)*w( 97)
    kab(3,4)= &
      +pk( 1)*w( 62)+pk( 2)*w( 63)+pk( 3)*w( 65)+pk( 4)*w( 68) &
      +pk( 5)*w( 72)+pk( 6)*w( 73)+pk( 7)*w( 75)+pk( 8)*w( 78) &
      +pk( 9)*w( 82)+pk(10)*w( 83)+pk(11)*w( 85)+pk(12)*w( 88) &
      +pk(13)*w( 92)+pk(14)*w( 93)+pk(15)*w( 95)+pk(16)*w( 98)
    kab(3,2)= &
      +pk( 1)*w( 64)+pk( 2)*w( 65)+pk( 3)*w( 66)+pk( 4)*w( 69) &
      +pk( 5)*w( 74)+pk( 6)*w( 75)+pk( 7)*w( 76)+pk( 8)*w( 79) &
      +pk( 9)*w( 84)+pk(10)*w( 85)+pk(11)*w( 86)+pk(12)*w( 89) &
      +pk(13)*w( 94)+pk(14)*w( 95)+pk(15)*w( 96)+pk(16)*w( 99)
    kab(3,3)= &
      +pk( 1)*w( 67)+pk( 2)*w( 68)+pk( 3)*w( 69)+pk( 4)*w( 70) &
      +pk( 5)*w( 77)+pk( 6)*w( 78)+pk( 7)*w( 79)+pk( 8)*w( 80) &
      +pk( 9)*w( 87)+pk(10)*w( 88)+pk(11)*w( 89)+pk(12)*w( 90) &
      +pk(13)*w( 97)+pk(14)*w( 98)+pk(15)*w( 99)+pk(16)*w(100)

  END SUBROUTINE exchange_sp

END MODULE nddo_methods

