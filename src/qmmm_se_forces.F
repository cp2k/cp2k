!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****** cp2k/qmmm_se_forces [1.1] *
!!
!!   FUNCTION
!!     Calculation of the derivative of the QMMM Hamiltonian integral 
!!     matrix <a|\sum_i q_i|b> for semi-empirical methods
!!
!!   AUTHOR
!!     Teodoro Laino - 04.2007 [tlaino]
!!
!!***
!******************************************************************************

MODULE qmmm_se_forces
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: do_method_am1,&
                                             do_method_mndo,&
                                             do_method_pchg,&
                                             do_method_pdg,&
                                             do_method_pm3,&
                                             do_method_pm6,&
                                             do_method_rm1
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_type,&
                                             qmmm_pot_p_type,&
                                             qmmm_pot_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE semi_empirical_int_ana,          ONLY: rotnuc_ana
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE semi_empirical_utils,            ONLY: allocate_se_param,&
                                             contract,&
                                             deallocate_se_param,&
                                             get_se_param,&
                                             get_se_type,&
                                             se_param_set_default
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_se_forces'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.
  PUBLIC :: deriv_se_qmmm_matrix

CONTAINS

!!****f* qmmm_se_forces/deriv_se_qmmm_matrix [1.0] *
!!
!!   NAME
!!     deriv_se_qmmm_matrix
!!
!!   FUNCTION
!!     Constructs the derivative w.r.t. 1-el semi-empirical hamiltonian
!!     QMMM terms
!!
!!   NOTES
!!     -
!!   AUTHOR
!!     Teodoro Laino 04.2007 [created]
!!
!!*** **********************************************************************
  SUBROUTINE deriv_se_qmmm_matrix(qs_env,qmmm_env,mm_particles,mm_cell,para_env,&
       calc_force, Forces, Forces_added_charges, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(in), OPTIONAL            :: calc_force
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces, Forces_added_charges
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deriv_se_qmmm_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ikind, iqm, &
                                                ispin, itype, natom, &
                                                natorb_a, nkind, &
                                                number_qm_atoms, stat
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: defined, failure
    REAL(KIND=dp)                            :: enuclear
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces_QM, p_block_a
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_mm

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,"I","",handle)
    IF (calc_force) THEN
       NULLIFY (rho, atomic_kind_set)
       NULLIFY (se_kind_a, se_kind_mm)
       CALL get_qs_env(qs_env=qs_env,&
                       rho=rho,&
                       atomic_kind_set=atomic_kind_set,&
                       ks_qmmm_env=ks_qmmm_env_loc,&
                       dft_control=dft_control,error=error)
       number_qm_atoms = SIZE(qs_env%subsys%particles%els)
       SELECT CASE (dft_control%qs_control%method_id)
       CASE (do_method_rm1,do_method_am1,do_method_mndo,do_method_pdg,do_method_pm3,do_method_pm6)
          ! Go on with the calculation.. 
       CASE DEFAULT
          ! Otherwise stop..
          CALL stop_program(routineN,moduleN,__LINE__,"Method not available",para_env)
       END SELECT
       ! Create a fake semi-empirical type to handle the classical atom
       ALLOCATE(Forces_QM(3,number_qm_atoms),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL allocate_se_param(se_kind_mm)
       CALL se_param_set_default(se_kind_mm,0,do_method_pchg)
       itype    = get_se_type(se_kind_mm%typ)
       nkind    = SIZE(atomic_kind_set)
       enuclear = 0.0_dp
       Forces_QM= 0.0_dp
       matrix_p => rho%rho_ao
       DO ispin = 1, dft_control%nspins
          iqm = 0
          Kinds: DO ikind=1,nkind
             atomic_kind => atomic_kind_set(ikind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                  natom=natom,&
                                  se_parameter=se_kind_a,&
                                  atom_list=list)
             CALL get_se_param(se_kind_a,&
                               defined=defined,&
                               natorb=natorb_a)
             IF (.NOT.defined .OR. natorb_a < 1) CYCLE
             Atoms: DO i = 1, SIZE(list)
                iqm = iqm + 1
                iatom = list(i)
                ! Give back block
                NULLIFY(p_block_a)
                CALL get_block_node(matrix=matrix_p(ispin)%matrix,&
                                    block_row=iatom,&
                                    block_col=iatom,&
                                    BLOCK=p_block_a)

                IF (ASSOCIATED(p_block_a)) THEN
                   ! Expand derivative of geometrical factors
                   CALL deriv_se_qmmm_matrix_low(p_block_a,&
                                                 se_kind_a,&
                                                 se_kind_mm,&
                                                 qmmm_env%Potentials,&
                                                 mm_particles,&
                                                 qmmm_env%mm_atom_chrg,&
                                                 qmmm_env%mm_el_pot_radius,&
                                                 qmmm_env%mm_atom_index,&
                                                 qmmm_env%num_mm_atoms,&
                                                 mm_cell,&
                                                 qmmm_env%qm_atom_index(iatom),&
                                                 itype,&
                                                 para_env,&
                                                 Forces,&
                                                 Forces_QM(:,iqm),&
                                                 error)
                   ! Possibly added charges
                   IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
                      CALL deriv_se_qmmm_matrix_low(p_block_a,&
                                                    se_kind_a,&
                                                    se_kind_mm,&
                                                    qmmm_env%added_charges%potentials,&
                                                    qmmm_env%added_charges%added_particles,&
                                                    qmmm_env%added_charges%mm_atom_chrg,&
                                                    qmmm_env%added_charges%mm_el_pot_radius,&
                                                    qmmm_env%added_charges%mm_atom_index,&
                                                    qmmm_env%added_charges%num_mm_atoms,&
                                                    mm_cell,&
                                                    qmmm_env%qm_atom_index(iatom),&
                                                    itype,&
                                                    para_env,&
                                                    Forces_added_charges,&
                                                    Forces_QM(:,iqm),&
                                                    error)
                   END IF
                END IF
             END DO Atoms
          END DO Kinds
       END DO
       CPPrecondition(iqm==number_qm_atoms,cp_failure_level,routineP,error,failure)
       ! Transfer QM gradients to the QM particles..
       CALL mp_sum(Forces_QM, para_env%group)
       iqm = 0
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               se_parameter=se_kind_a,&
                               atom_list=list)
          CALL get_se_param(se_kind_a,&
                            defined=defined,&
                            natorb=natorb_a)
          IF (.NOT.defined .OR. natorb_a < 1) CYCLE
          DO i = 1, SIZE(list)
             iqm = iqm + 1
             iatom = qmmm_env%qm_atom_index(list(i))
             mm_particles(iatom)%f(:) = mm_particles(iatom)%f(:) + Forces_QM(:,iqm)
          END DO
       END DO
       ! MM forces will be handled directly from the QMMM module in the same way
       ! as for GPW/GAPW methods
       DEALLOCATE(Forces_QM,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL deallocate_se_param(se_kind_mm)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE deriv_se_qmmm_matrix

!!****f* qmmm_se_forces/deriv_se_qmmm_matrix_low [1.0] *
!!
!!   NAME
!!     deriv_se_qmmm_matrix_low
!!
!!   FUNCTION
!!     Low Level : Computes derivatives of the 1-el semi-empirical QMMM 
!!                 hamiltonian block w.r.t. MM and QM coordinates
!!
!!   NOTES
!!     -
!!   AUTHOR
!!     Teodoro Laino 04.2007 [created]
!!
!!*** **********************************************************************
  SUBROUTINE deriv_se_qmmm_matrix_low(p_block_a, se_kind_a, se_kind_mm,&
       potentials, mm_particles, mm_charges, mm_el_pot_radius, mm_atom_index, &
       num_mm_atoms, mm_cell, IndQM, itype, para_env, forces, forces_qm, error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block_a
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_mm
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: potentials
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges, mm_el_pot_radius
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    TYPE(cell_type), POINTER                 :: mm_cell
    INTEGER, INTENT(IN)                      :: IndQM, itype
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: forces
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: forces_qm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deriv_se_qmmm_matrix_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icor, Imm, Imp, &
                                                IndMM, Ipot
    REAL(KIND=dp)                            :: r, rt1, rt2, rt3
    REAL(KIND=dp), DIMENSION(10)             :: pvec
    REAL(KIND=dp), DIMENSION(10, 3)          :: de1b
    REAL(KIND=dp), DIMENSION(3)              :: denuc, force_ab, r_pbc, rij
    TYPE(qmmm_pot_type), POINTER             :: Pot

    CALL timeset(routineN,"I","",handle)
    ! Loop Over MM atoms - parallelization over MM atoms...
    ! Loop over Pot stores atoms with the same charge
    MainLoopPot: DO Ipot = 1, SIZE(Potentials)
       Pot    => Potentials(Ipot)%Pot
       ! Loop over atoms belonging to this type
       LoopMM: DO Imp = 1, SIZE(Pot%mm_atom_index)
          Imm = Pot%mm_atom_index(Imp)
          IndMM = mm_atom_index(Imm) 
          r_pbc=pbc(mm_particles(IndMM)%r-mm_particles(IndQM)%r, mm_cell)
          rt1= r_pbc(1) 
          rt2= r_pbc(2) 
          rt3= r_pbc(3) 
          rij = (/rt1,rt2,rt3/)
          se_kind_mm%zeff = mm_charges(Imm)
          ! Integrals derivatives involving QM - MM atoms
          CALL rotnuc_ana (se_kind_a,se_kind_mm,rij,itype,denuc=denuc,de1b=de1b,shortrange=.FALSE.)
          ! Nucler - Nuclear term
          force_ab(1:3)=-denuc(1:3)
          IF ( SIZE(p_block_a,1) == 1 ) THEN
             force_ab(1:3) = force_ab(1:3) - de1b(1,1:3)*p_block_a(1,1)
          ELSE
             CALL contract ( p_block_a, pvec )
             DO icor=1,3
                force_ab(icor) = force_ab(icor) - SUM(de1b(1:10,icor)*pvec(1:10))
             END DO
          END IF
          ! The array of QM forces are really the forces
          forces_qm(:)  = forces_qm(:)  - force_ab
          ! The one of MM atoms are instead gradients
          forces(:,Imm) = forces(:,Imm) - force_ab
       END DO LoopMM
    END DO MainLoopPot
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE deriv_se_qmmm_matrix_low

END MODULE qmmm_se_forces
