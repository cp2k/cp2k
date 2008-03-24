!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the QMMM Hamiltonian integral matrix <a|\sum_i q_i|b> for
!>      semi-empirical methods
!> \author Teodoro Laino 04.2007 [created]
! *****************************************************************************
MODULE qmmm_se_energy
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_method_am1, do_method_mndo, do_method_pchg, do_method_pdg, &
       do_method_pm3, do_method_pm6, do_method_rm1, do_qmmm_coulomb, &
       do_qmmm_gauss, do_qmmm_none, do_qmmm_swave
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_type,&
                                             qmmm_pot_p_type,&
                                             qmmm_pot_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
  USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
  USE qs_overlap,                      ONLY: build_overlap_matrix
  USE semi_empirical_int_ana,          ONLY: rotnuc_ana
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE semi_empirical_utils,            ONLY: allocate_se_param,&
                                             deallocate_se_param,&
                                             expand,&
                                             get_se_param,&
                                             get_se_type,&
                                             se_param_set_default
  USE sparse_matrix_output,            ONLY: write_sparse_matrix
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_se_energy'

  PUBLIC :: build_se_qmmm_matrix

CONTAINS

! *****************************************************************************
!> \brief Constructs the 1-el semi-empirical hamiltonian
!> \author Teodoro Laino 04.2007 [created]
! *****************************************************************************
  SUBROUTINE build_se_qmmm_matrix(qs_env,qmmm_env,mm_particles,mm_cell,para_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_se_qmmm_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ikind, &
                                                itype, iw, natom, natorb_a, &
                                                nkind
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: defined
    REAL(KIND=dp)                            :: enuclear
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block_a
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_s
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_mm

    CALL timeset(routineN,"I","",handle)
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    NULLIFY (matrix_h, matrix_s, atomic_kind_set)
    NULLIFY (se_kind_a, se_kind_mm)
    CALL build_qs_neighbor_lists(qs_env,para_env,force_env_section=qs_env%input,error=error)
    CALL build_overlap_matrix(qs_env=qs_env,para_env=para_env,error=error)

    CALL get_qs_env(qs_env=qs_env,&
                    matrix_s=matrix_s,&
                    atomic_kind_set=atomic_kind_set,&
                    ks_qmmm_env=ks_qmmm_env_loc,&
                    dft_control=dft_control,error=error)

    SELECT CASE (dft_control%qs_control%method_id)
    CASE (do_method_am1,do_method_rm1,do_method_mndo,do_method_pdg,do_method_pm3,do_method_pm6)
       ! Go on with the calculation.. 
    CASE DEFAULT
       ! Otherwise stop..
       CALL stop_program(routineN,moduleN,__LINE__,&
            "Method not available",para_env)
    END SELECT

    ! Allocate the core Hamiltonian matrix
    CALL allocate_matrix_set(ks_qmmm_env_loc%matrix_h,1,error=error)
    CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                    TARGET=ks_qmmm_env_loc%matrix_h(1)%matrix,&
                                    target_name="QMMM HAMILTONIAN MATRIX",error=error)
    SELECT CASE(qmmm_env%qmmm_coupl_type)
    CASE(do_qmmm_coulomb,do_qmmm_gauss,do_qmmm_swave)
       ! Create a fake semi-empirical type to handle the classical atom
       CALL allocate_se_param(se_kind_mm,error)
       CALL se_param_set_default(se_kind_mm,0,do_method_pchg)
       itype    = get_se_type(se_kind_mm%typ)
       nkind    = SIZE(atomic_kind_set)
       enuclear = 0.0_dp
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
             iatom = list(i)
             ! Give back block
             NULLIFY(h_block_a)
             CALL get_block_node(matrix=ks_qmmm_env_loc%matrix_h(1)%matrix,&
                                 block_row=iatom,&
                                 block_col=iatom,&
                                 BLOCK=h_block_a)

             IF (ASSOCIATED(h_block_a)) THEN
                h_block_a = 0.0_dp
                ! Real QM/MM computation
                CALL build_se_qmmm_matrix_low(h_block_a,&
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
                                              enuclear,&
                                              itype,&
                                              error)
                ! Possibly added charges
                IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
                   CALL build_se_qmmm_matrix_low(h_block_a,&
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
                                                 enuclear,&
                                                 itype,&
                                                 error)       
                END IF
             END IF
          END DO Atoms
       END DO Kinds
       CALL mp_sum(enuclear,para_env%group)
       qs_env%energy%qmmm_nu = enuclear
       CALL deallocate_se_param(se_kind_mm,error)
    CASE(do_qmmm_none)
       ! Zero Matrix 
       CALL set_matrix(ks_qmmm_env_loc%matrix_h(1)%matrix,0.0_dp)
    END SELECT
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"QMMM%PRINT%QMMM_MATRIX",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"QMMM%PRINT%QMMM_MATRIX",&
            extension=".Log",error=error)
       CALL write_sparse_matrix(ks_qmmm_env_loc%matrix_h(1)%matrix,4,6,qs_env,para_env,&
            scale=1.0_dp,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "QMMM%PRINT%QMMM_MATRIX", error=error)
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_se_qmmm_matrix

! *****************************************************************************
!> \brief Low Level : Constructs the 1-el semi-empirical hamiltonian block
!> \author Teodoro Laino 04.2007 [created]
! *****************************************************************************
  SUBROUTINE build_se_qmmm_matrix_low(h_block_a, se_kind_a, se_kind_mm, potentials,&
       mm_particles, mm_charges, mm_el_pot_radius, mm_atom_index, &
       num_mm_atoms, mm_cell, IndQM, enuclear, itype, error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block_a
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_mm
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: potentials
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges, mm_el_pot_radius
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    TYPE(cell_type), POINTER                 :: mm_cell
    INTEGER, INTENT(IN)                      :: IndQM
    REAL(KIND=dp), INTENT(INOUT)             :: enuclear
    INTEGER, INTENT(IN)                      :: itype
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_se_qmmm_matrix_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, Imm, Imp, IndMM, &
                                                Ipot
    REAL(KIND=dp)                            :: enuc, r, rt1, rt2, rt3
    REAL(KIND=dp), DIMENSION(10)             :: e1b
    REAL(KIND=dp), DIMENSION(3)              :: r_pbc, rij
    TYPE(qmmm_pot_type), POINTER             :: Pot

    CALL timeset(routineN,"I","",handle)
    ! Loop Over MM atoms
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
          CALL rotnuc_ana (se_kind_a,se_kind_mm,rij,itype,enuc,e1b,shortrange=.FALSE.,&
               error=error)
          enuclear = enuclear + enuc
          ! Contribution to the iatom block
          ! Computation of the QMMM core matrix
          IF ( SIZE(h_block_a,1) == 1 ) THEN
             h_block_a(1,1)=h_block_a(1,1)+e1b(1)
          ELSE
             CALL expand (h_block_a, e1b )
          END IF
       END DO LoopMM
    END DO MainLoopPot
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE build_se_qmmm_matrix_low

END MODULE qmmm_se_energy
