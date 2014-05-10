!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the derivative of the QMMM Hamiltonian integral
!>      matrix <a|\sum_i q_i|b> for semi-empirical methods
!> \author Teodoro Laino - 04.2007 [tlaino]
! *****************************************************************************
MODULE qmmm_se_forces
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_method_am1, do_method_mndo, do_method_mndod, do_method_pchg, &
       do_method_pdg, do_method_pm3, do_method_pm6, do_method_pnnl, &
       do_method_rm1, do_multipole_none
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_type,&
                                             qmmm_pot_p_type,&
                                             qmmm_pot_type
  USE qmmm_util,                       ONLY: spherical_cutoff_factor
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE semi_empirical_int_arrays,       ONLY: se_orbital_pointer
  USE semi_empirical_integrals,        ONLY: dcorecore,&
                                             drotnuc
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             se_int_control_type,&
                                             se_taper_type,&
                                             semi_empirical_create,&
                                             semi_empirical_release,&
                                             semi_empirical_type,&
                                             setup_se_int_control_type
  USE semi_empirical_utils,            ONLY: get_se_type,&
                                             se_param_set_default
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_se_forces'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
  PUBLIC :: deriv_se_qmmm_matrix

CONTAINS

! *****************************************************************************
!> \brief Constructs the derivative w.r.t. 1-el semi-empirical hamiltonian
!>      QMMM terms
!> \author Teodoro Laino 04.2007 [created]
! *****************************************************************************
  SUBROUTINE deriv_se_qmmm_matrix(qs_env,qmmm_env,particles_mm,mm_cell,para_env,&
       calc_force, Forces, Forces_added_charges, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
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
    LOGICAL                                  :: anag, defined, failure, found
    REAL(KIND=dp)                            :: delta, enuclear
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces_QM, p_block_a
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_qm
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(se_int_control_type)                :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_mm

    failure = .FALSE.
    CALL timeset(routineN,handle)
    IF (calc_force) THEN
       NULLIFY (rho, atomic_kind_set, se_taper)
       NULLIFY (se_kind_a, se_kind_mm, particles_qm)
       CALL get_qs_env(qs_env=qs_env,&
                       rho=rho,&
                       se_taper=se_taper,&
                       atomic_kind_set=atomic_kind_set,&
                       ks_qmmm_env=ks_qmmm_env_loc,&
                       dft_control=dft_control,&
                       particle_set=particles_qm,&
                       error=error)
       number_qm_atoms = SIZE(qs_env%subsys%particles%els)
       SELECT CASE (dft_control%qs_control%method_id)
       CASE (do_method_rm1,do_method_am1,do_method_mndo,do_method_pdg,&
             do_method_pm3,do_method_pm6,do_method_mndod,do_method_pnnl)
          ! Go on with the calculation..
       CASE DEFAULT
          ! Otherwise stop..
          CALL stop_program(routineN,moduleN,__LINE__,"Method not available",para_env)
       END SELECT
       anag = dft_control%qs_control%se_control%analytical_gradients
       delta= dft_control%qs_control%se_control%delta
       ! Setup SE integral control type
       CALL setup_se_int_control_type(se_int_control, shortrange=.FALSE., do_ewald_r3=.FALSE.,&
            do_ewald_gks=.FALSE., integral_screening=dft_control%qs_control%se_control%integral_screening,&
            max_multipole=do_multipole_none, pc_coulomb_int=.FALSE.)

       ! Create a fake semi-empirical type to handle the classical atom
       ALLOCATE(Forces_QM(3,number_qm_atoms),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL semi_empirical_create(se_kind_mm,error)
       CALL se_param_set_default(se_kind_mm,0,do_method_pchg,error)
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
                CALL cp_dbcsr_get_block_p(matrix=matrix_p(ispin)%matrix,&
                     row=iatom,col=iatom,BLOCK=p_block_a,found=found)

                IF (ASSOCIATED(p_block_a)) THEN
                   ! Expand derivative of geometrical factors
                   CALL deriv_se_qmmm_matrix_low(p_block_a,&
                                                 se_kind_a,&
                                                 se_kind_mm,&
                                                 qmmm_env%Potentials,&
                                                 particles_mm,&
                                                 qmmm_env%mm_atom_chrg,&
                                                 qmmm_env%mm_el_pot_radius,&
                                                 qmmm_env%mm_atom_index,&
                                                 qmmm_env%num_mm_atoms,&
                                                 mm_cell,&
                                                 iatom,&
                                                 itype,&
                                                 Forces,&
                                                 Forces_QM(:,iqm),&
                                                 se_taper,&
                                                 se_int_control,&
                                                 anag,&
                                                 delta,&
                                                 qmmm_env%spherical_cutoff,&
                                                 particles_qm,&
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
                                                    iatom,&
                                                    itype,&
                                                    Forces_added_charges,&
                                                    Forces_QM(:,iqm),&
                                                    se_taper,&
                                                    se_int_control,&
                                                    anag,&
                                                    delta,&
                                                    qmmm_env%spherical_cutoff,&
                                                    particles_qm,&
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
             particles_mm(iatom)%f(:) = particles_mm(iatom)%f(:) + Forces_QM(:,iqm)
          END DO
       END DO
       ! MM forces will be handled directly from the QMMM module in the same way
       ! as for GPW/GAPW methods
       DEALLOCATE(Forces_QM,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL semi_empirical_release(se_kind_mm,error)

    END IF
    CALL timestop(handle)
  END SUBROUTINE deriv_se_qmmm_matrix

! *****************************************************************************
!> \brief Low Level : Computes derivatives of the 1-el semi-empirical QMMM
!>                  hamiltonian block w.r.t. MM and QM coordinates
!> \author Teodoro Laino 04.2007 [created]
! *****************************************************************************
  SUBROUTINE deriv_se_qmmm_matrix_low(p_block_a, se_kind_a, se_kind_mm,&
       potentials, particles_mm, mm_charges, mm_el_pot_radius, mm_atom_index, &
       num_mm_atoms, mm_cell, IndQM, itype, forces, forces_qm, se_taper,&
       se_int_control, anag, delta, qmmm_spherical_cutoff, particles_qm, error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block_a
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_mm
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: potentials
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_mm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges, mm_el_pot_radius
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    TYPE(cell_type), POINTER                 :: mm_cell
    INTEGER, INTENT(IN)                      :: IndQM, itype
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: forces
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: forces_qm
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    LOGICAL, INTENT(IN)                      :: anag
    REAL(KIND=dp), INTENT(IN)                :: delta, &
                                                qmmm_spherical_cutoff(2)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_qm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deriv_se_qmmm_matrix_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i1, i1L, i2, Imm, &
                                                Imp, IndMM, Ipot, j1, j1L
    REAL(KIND=dp)                            :: rt1, rt2, rt3, sph_chrg_factor
    REAL(KIND=dp), DIMENSION(3)              :: denuc, force_ab, r_pbc, rij
    REAL(KIND=dp), DIMENSION(3, 45)          :: de1b
    TYPE(qmmm_pot_type), POINTER             :: Pot

    CALL timeset(routineN,handle)
    ! Loop Over MM atoms - parallelization over MM atoms...
    ! Loop over Pot stores atoms with the same charge
    MainLoopPot: DO Ipot = 1, SIZE(Potentials)
       Pot    => Potentials(Ipot)%Pot
       ! Loop over atoms belonging to this type
       LoopMM: DO Imp = 1, SIZE(Pot%mm_atom_index)
          Imm = Pot%mm_atom_index(Imp)
          IndMM = mm_atom_index(Imm)
          r_pbc=pbc(particles_mm(IndMM)%r-particles_qm(IndQM)%r, mm_cell)
          rt1= r_pbc(1)
          rt2= r_pbc(2)
          rt3= r_pbc(3)
          rij = (/rt1,rt2,rt3/)
          se_kind_mm%zeff = mm_charges(Imm)
          ! Computes the screening factor for the spherical cutoff
          IF (qmmm_spherical_cutoff(1)>0.0_dp) THEN
             CALL spherical_cutoff_factor(qmmm_spherical_cutoff, rij, sph_chrg_factor, error)
             se_kind_mm%zeff = se_kind_mm%zeff * sph_chrg_factor
          END IF
          IF (ABS(se_kind_mm%zeff)<=EPSILON(0.0_dp)) CYCLE
          ! Integrals derivatives involving QM - MM atoms
          CALL drotnuc(se_kind_a, se_kind_mm, rij, itype=itype, de1b=de1b,&
               se_int_control=se_int_control, anag=anag, delta=delta, &
               se_taper=se_taper, error=error)
          CALL dcorecore(se_kind_a, se_kind_mm, rij, itype=itype, denuc=denuc,&
               se_int_control=se_int_control, anag=anag, delta=delta, &
               se_taper=se_taper, error=error)
          ! Nucler - Nuclear term
          force_ab(1:3)=-denuc(1:3)
          ! Force contribution from the QMMM Hamiltonian
          i2 = 0
          DO i1L = 1, se_kind_a%natorb
             i1 = se_orbital_pointer(i1L)
             DO j1L = 1, i1L-1
                j1 = se_orbital_pointer(j1L)
                i2 = i2 + 1
                force_ab = force_ab - 2.0_dp * de1b(:,i2)*p_block_a(i1,j1)
             END DO
             j1 = se_orbital_pointer(j1L)
             i2 = i2 + 1
             force_ab = force_ab -  de1b(:,i2)*p_block_a(i1,j1)
          END DO
          ! The array of QM forces are really the forces
          forces_qm(:)  = forces_qm(:)  - force_ab
          ! The one of MM atoms are instead gradients
          forces(:,Imm) = forces(:,Imm) - force_ab
       END DO LoopMM
    END DO MainLoopPot
    CALL timestop(handle)
  END SUBROUTINE deriv_se_qmmm_matrix_low

END MODULE qmmm_se_forces
