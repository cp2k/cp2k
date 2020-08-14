!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief TB methods used with QMMM
!> \author JGH
! **************************************************************************************************
MODULE qmmm_tb_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type,&
                                              dftb_control_type,&
                                              xtb_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_get_block_p, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_set
   USE ewald_environment_types,         ONLY: ewald_env_create,&
                                              ewald_env_get,&
                                              ewald_env_release,&
                                              ewald_env_set,&
                                              ewald_environment_type,&
                                              read_ewald_section
   USE ewald_pw_types,                  ONLY: ewald_pw_create,&
                                              ewald_pw_release,&
                                              ewald_pw_type
   USE input_constants,                 ONLY: do_fist_pol_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE mulliken,                        ONLY: mulliken_charges
   USE particle_types,                  ONLY: allocate_particle_set,&
                                              deallocate_particle_set,&
                                              particle_type
   USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                              do_ewald_none,&
                                              do_ewald_pme,&
                                              do_ewald_spme
   USE qmmm_types_low,                  ONLY: qmmm_env_qm_type,&
                                              qmmm_pot_p_type,&
                                              qmmm_pot_type
   USE qmmm_util,                       ONLY: spherical_cutoff_factor
   USE qs_dftb_coulomb,                 ONLY: gamma_rab_sr
   USE qs_dftb_matrices,                ONLY: build_dftb_overlap
   USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
   USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE spme,                            ONLY: spme_forces,&
                                              spme_potential
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   ! small real number
   REAL(dp), PARAMETER                    :: rtiny = 1.e-10_dp
   ! eta(0) for mm atoms and non-scc qm atoms
   REAL(dp), PARAMETER                    :: eta_mm = 0.47_dp
   ! step size for qmmm finite difference
   REAL(dp), PARAMETER                    :: ddrmm = 0.0001_dp

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_tb_methods'

   PUBLIC :: build_tb_qmmm_matrix, build_tb_qmmm_matrix_zero, &
             build_tb_qmmm_matrix_pc, deriv_tb_qmmm_matrix, &
             deriv_tb_qmmm_matrix_pc

CONTAINS

! **************************************************************************************************
!> \brief Constructs the 1-el DFTB hamiltonian
!> \param qs_env ...
!> \param qmmm_env ...
!> \param particles_mm ...
!> \param mm_cell ...
!> \param para_env ...
!> \author JGH 10.2014 [created]
! **************************************************************************************************
   SUBROUTINE build_tb_qmmm_matrix(qs_env, qmmm_env, particles_mm, mm_cell, para_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qmmm_env_qm_type), POINTER                    :: qmmm_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_mm
      TYPE(cell_type), POINTER                           :: mm_cell
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_tb_qmmm_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, handle, i, iatom, ikind, jatom, &
                                                            natom, natorb, nkind
      INTEGER, DIMENSION(:), POINTER                     :: list
      LOGICAL                                            :: defined, do_dftb, do_xtb, found
      REAL(KIND=dp)                                      :: pc_ener, zeff
      REAL(KIND=dp), DIMENSION(0:3)                      :: eta_a
      REAL(KIND=dp), DIMENSION(:), POINTER               :: qpot
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hblock, sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_h, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dftb_control_type), POINTER                   :: dftb_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_qm
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_kind
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_ks_qmmm_env_type), POINTER                 :: ks_qmmm_env_loc
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particles_qm, &
                      qs_kind_set=qs_kind_set, &
                      rho=rho, &
                      natom=natom)
      dftb_control => dft_control%qs_control%dftb_control
      xtb_control => dft_control%qs_control%xtb_control

      IF (dft_control%qs_control%dftb) THEN
         do_dftb = .TRUE.
         do_xtb = .FALSE.
      ELSEIF (dft_control%qs_control%xtb) THEN
         do_dftb = .FALSE.
         do_xtb = .TRUE.
      ELSE
         CPABORT("TB method unknown")
      END IF

      CALL build_qs_neighbor_lists(qs_env, para_env, force_env_section=qs_env%input)

      NULLIFY (matrix_s)
      IF (do_dftb) THEN
         CALL build_dftb_overlap(qs_env, 0, matrix_s)
      ELSEIF (do_xtb) THEN
         CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, sab_orb=sab_nl)
         CALL build_overlap_matrix(ks_env, matrix_s, basis_type_a='ORB', basis_type_b='ORB', sab_nl=sab_nl)
      END IF

      ALLOCATE (qpot(natom))
      qpot = 0.0_dp
      pc_ener = 0.0_dp

      nkind = SIZE(atomic_kind_set)
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
         IF (do_dftb) THEN
            NULLIFY (dftb_kind)
            CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
            CALL get_dftb_atom_param(dftb_kind, zeff=zeff, &
                                     defined=defined, eta=eta_a, natorb=natorb)
            ! use mm charge smearing for non-scc cases
            IF (.NOT. dftb_control%self_consistent) eta_a(0) = eta_mm
            IF (.NOT. defined .OR. natorb < 1) CYCLE
         ELSEIF (do_xtb) THEN
            NULLIFY (xtb_kind)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, zeff=zeff)
            eta_a(0) = eta_mm
         END IF
         DO i = 1, SIZE(list)
            iatom = list(i)
            CALL build_mm_pot(qpot(iatom), 0, eta_a(0), qmmm_env%Potentials, particles_mm, &
                              qmmm_env%mm_atom_chrg, qmmm_env%mm_atom_index, mm_cell, iatom, &
                              qmmm_env%spherical_cutoff, particles_qm)
            ! Possibly added charges
            IF (qmmm_env%move_mm_charges .OR. qmmm_env%add_mm_charges) THEN
               CALL build_mm_pot(qpot(iatom), 0, eta_a(0), qmmm_env%added_charges%potentials, &
                                 qmmm_env%added_charges%added_particles, &
                                 qmmm_env%added_charges%mm_atom_chrg, &
                                 qmmm_env%added_charges%mm_atom_index, mm_cell, iatom, &
                                 qmmm_env%spherical_cutoff, &
                                 particles_qm)
            END IF
            pc_ener = pc_ener + qpot(iatom)*zeff
         END DO
      END DO

      ! Allocate the core Hamiltonian matrix
      CALL get_qs_env(qs_env=qs_env, ks_qmmm_env=ks_qmmm_env_loc)
      matrix_h => ks_qmmm_env_loc%matrix_h
      CALL dbcsr_allocate_matrix_set(matrix_h, 1)
      ALLOCATE (matrix_h(1)%matrix)
      CALL dbcsr_copy(matrix_h(1)%matrix, matrix_s(1)%matrix, &
                      name="QMMM HAMILTONIAN MATRIX")
      CALL dbcsr_set(matrix_h(1)%matrix, 0.0_dp)

      CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, iatom, jatom, sblock, blk)
         NULLIFY (hblock)
         CALL dbcsr_get_block_p(matrix=matrix_h(1)%matrix, &
                                row=iatom, col=jatom, block=hblock, found=found)
         CPASSERT(found)
         hblock = hblock - 0.5_dp*sblock*(qpot(iatom) + qpot(jatom))
      END DO
      CALL dbcsr_iterator_stop(iter)

      ks_qmmm_env_loc%matrix_h => matrix_h
      ks_qmmm_env_loc%pc_ener = pc_ener

      DEALLOCATE (qpot)

      CALL dbcsr_deallocate_matrix_set(matrix_s)

      CALL timestop(handle)

   END SUBROUTINE build_tb_qmmm_matrix

! **************************************************************************************************
!> \brief Constructs an empty 1-el DFTB hamiltonian
!> \param qs_env ...
!> \param para_env ...
!> \author JGH 10.2014 [created]
! **************************************************************************************************
   SUBROUTINE build_tb_qmmm_matrix_zero(qs_env, para_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_tb_qmmm_matrix_zero', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: do_dftb, do_xtb
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_h, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_ks_qmmm_env_type), POINTER                 :: ks_qmmm_env_loc

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)

      IF (dft_control%qs_control%dftb) THEN
         do_dftb = .TRUE.
         do_xtb = .FALSE.
      ELSEIF (dft_control%qs_control%xtb) THEN
         do_dftb = .FALSE.
         do_xtb = .TRUE.
      ELSE
         CPABORT("TB method unknown")
      END IF

      CALL build_qs_neighbor_lists(qs_env, para_env, force_env_section=qs_env%input)

      NULLIFY (matrix_s)
      IF (do_dftb) THEN
         CALL build_dftb_overlap(qs_env, 0, matrix_s)
      ELSEIF (do_xtb) THEN
         CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, sab_orb=sab_nl)
         CALL build_overlap_matrix(ks_env, matrix_s, basis_type_a='ORB', basis_type_b='ORB', sab_nl=sab_nl)
      END IF

      ! Allocate the core Hamiltonian matrix
      CALL get_qs_env(qs_env=qs_env, ks_qmmm_env=ks_qmmm_env_loc)
      matrix_h => ks_qmmm_env_loc%matrix_h
      CALL dbcsr_allocate_matrix_set(matrix_h, 1)
      ALLOCATE (matrix_h(1)%matrix)
      CALL dbcsr_copy(matrix_h(1)%matrix, matrix_s(1)%matrix, &
                      name="QMMM HAMILTONIAN MATRIX")
      CALL dbcsr_set(matrix_h(1)%matrix, 0.0_dp)
      ks_qmmm_env_loc%matrix_h => matrix_h
      ks_qmmm_env_loc%pc_ener = 0.0_dp

      CALL dbcsr_deallocate_matrix_set(matrix_s)

      CALL timestop(handle)

   END SUBROUTINE build_tb_qmmm_matrix_zero

! **************************************************************************************************
!> \brief Constructs the 1-el DFTB hamiltonian
!> \param qs_env ...
!> \param qmmm_env ...
!> \param particles_mm ...
!> \param mm_cell ...
!> \param para_env ...
!> \author JGH 10.2014 [created]
! **************************************************************************************************
   SUBROUTINE build_tb_qmmm_matrix_pc(qs_env, qmmm_env, particles_mm, mm_cell, para_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qmmm_env_qm_type), POINTER                    :: qmmm_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_mm
      TYPE(cell_type), POINTER                           :: mm_cell
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'build_tb_qmmm_matrix_pc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, do_ipol, ewald_type, handle, i, &
                                                            iatom, ikind, imm, imp, indmm, ipot, &
                                                            jatom, natom, natorb, nkind, nmm
      INTEGER, DIMENSION(:), POINTER                     :: list
      LOGICAL                                            :: defined, do_dftb, do_multipoles, do_xtb, &
                                                            found
      REAL(KIND=dp)                                      :: alpha, pc_ener, zeff
      REAL(KIND=dp), DIMENSION(0:3)                      :: eta_a
      REAL(KIND=dp), DIMENSION(2)                        :: rcutoff
      REAL(KIND=dp), DIMENSION(:), POINTER               :: charges_mm, qpot
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hblock, sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_h, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dftb_control_type), POINTER                   :: dftb_control
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(particle_type), DIMENSION(:), POINTER         :: atoms_mm, particles_qm
      TYPE(qmmm_pot_type), POINTER                       :: Pot
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_kind
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_ks_qmmm_env_type), POINTER                 :: ks_qmmm_env_loc
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: ewald_section, poisson_section, &
                                                            print_section
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particles_qm, &
                      qs_kind_set=qs_kind_set, &
                      rho=rho, &
                      natom=natom)
      dftb_control => dft_control%qs_control%dftb_control
      xtb_control => dft_control%qs_control%xtb_control

      IF (dft_control%qs_control%dftb) THEN
         do_dftb = .TRUE.
         do_xtb = .FALSE.
      ELSEIF (dft_control%qs_control%xtb) THEN
         do_dftb = .FALSE.
         do_xtb = .TRUE.
      ELSE
         CPABORT("TB method unknown")
      END IF

      CALL build_qs_neighbor_lists(qs_env, para_env, force_env_section=qs_env%input)

      NULLIFY (matrix_s)
      IF (do_dftb) THEN
         CALL build_dftb_overlap(qs_env, 0, matrix_s)
      ELSEIF (do_xtb) THEN
         CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, sab_orb=sab_nl)
         CALL build_overlap_matrix(ks_env, matrix_s, basis_type_a='ORB', basis_type_b='ORB', sab_nl=sab_nl)
      END IF

      ALLOCATE (qpot(natom))
      qpot = 0.0_dp
      pc_ener = 0.0_dp

      ! Create Ewald environments
      poisson_section => section_vals_get_subs_vals(qs_env%input, "MM%POISSON")
      CALL ewald_env_create(ewald_env, para_env)
      CALL ewald_env_set(ewald_env, poisson_section=poisson_section)
      ewald_section => section_vals_get_subs_vals(poisson_section, "EWALD")
      CALL read_ewald_section(ewald_env, ewald_section)
      print_section => section_vals_get_subs_vals(qs_env%input, "PRINT%GRID_INFORMATION")
      CALL ewald_pw_create(ewald_pw, ewald_env, mm_cell, mm_cell, print_section=print_section)

      CALL ewald_env_get(ewald_env, ewald_type=ewald_type, do_multipoles=do_multipoles, do_ipol=do_ipol)
      IF (do_multipoles) CPABORT("No multipole force fields allowed in TB QM/MM")
      IF (do_ipol /= do_fist_pol_none) CPABORT("No polarizable force fields allowed in TB QM/MM")

      SELECT CASE (ewald_type)
      CASE (do_ewald_pme)
         CPABORT("PME Ewald type not implemented for TB/QMMM")
      CASE (do_ewald_ewald, do_ewald_spme)
         DO ipot = 1, SIZE(qmmm_env%Potentials)
            Pot => qmmm_env%Potentials(ipot)%Pot
            nmm = SIZE(Pot%mm_atom_index)
            ! get a 'clean' mm particle set
            NULLIFY (atoms_mm)
            CALL allocate_particle_set(atoms_mm, nmm)
            ALLOCATE (charges_mm(nmm))
            DO Imp = 1, nmm
               Imm = Pot%mm_atom_index(Imp)
               IndMM = qmmm_env%mm_atom_index(Imm)
               atoms_mm(Imp)%r = particles_mm(IndMM)%r
               atoms_mm(Imp)%atomic_kind => particles_mm(IndMM)%atomic_kind
               charges_mm(Imp) = qmmm_env%mm_atom_chrg(Imm)
            END DO
            IF (ewald_type == do_ewald_ewald) THEN
               CPABORT("Ewald not implemented for TB/QMMM")
            ELSE IF (ewald_type == do_ewald_spme) THEN
               ! spme electrostatic potential
               CALL spme_potential(ewald_env, ewald_pw, mm_cell, atoms_mm, charges_mm, particles_qm, qpot)
            END IF
            CALL deallocate_particle_set(atoms_mm)
            DEALLOCATE (charges_mm)
         END DO
         IF (qmmm_env%move_mm_charges .OR. qmmm_env%add_mm_charges) THEN
            DO ipot = 1, SIZE(qmmm_env%added_charges%Potentials)
               Pot => qmmm_env%added_charges%Potentials(ipot)%Pot
               nmm = SIZE(Pot%mm_atom_index)
               ! get a 'clean' mm particle set
               NULLIFY (atoms_mm)
               CALL allocate_particle_set(atoms_mm, nmm)
               ALLOCATE (charges_mm(nmm))
               DO Imp = 1, nmm
                  Imm = Pot%mm_atom_index(Imp)
                  IndMM = qmmm_env%added_charges%mm_atom_index(Imm)
                  atoms_mm(Imp)%r = qmmm_env%added_charges%added_particles(IndMM)%r
                  atoms_mm(Imp)%atomic_kind => qmmm_env%added_charges%added_particles(IndMM)%atomic_kind
                  charges_mm(Imp) = qmmm_env%added_charges%mm_atom_chrg(Imm)
               END DO
               IF (ewald_type == do_ewald_ewald) THEN
                  CPABORT("Ewald not implemented for TB/QMMM")
               ELSE IF (ewald_type == do_ewald_spme) THEN
                  ! spme electrostatic potential
                  CALL spme_potential(ewald_env, ewald_pw, mm_cell, atoms_mm, charges_mm, particles_qm, qpot)
               END IF
               CALL deallocate_particle_set(atoms_mm)
               DEALLOCATE (charges_mm)
            END DO
         END IF
         CALL mp_sum(qpot, para_env%group)
         ! Add Ewald and TB short range corrections
         ! This is effectively using a minimum image convention!
         ! Set rcutoff to values compatible with alpha Ewald
         CALL ewald_env_get(ewald_env, rcut=rcutoff(1), alpha=alpha)
         rcutoff(2) = 0.025_dp*rcutoff(1)
         rcutoff(1) = 2.0_dp*rcutoff(1)
         nkind = SIZE(atomic_kind_set)
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
            IF (do_dftb) THEN
               NULLIFY (dftb_kind)
               CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
               CALL get_dftb_atom_param(dftb_kind, zeff=zeff, &
                                        defined=defined, eta=eta_a, natorb=natorb)
               ! use mm charge smearing for non-scc cases
               IF (.NOT. dftb_control%self_consistent) eta_a(0) = eta_mm
               IF (.NOT. defined .OR. natorb < 1) CYCLE
            ELSEIF (do_xtb) THEN
               NULLIFY (xtb_kind)
               CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
               CALL get_xtb_atom_param(xtb_kind, zeff=zeff)
               eta_a(0) = eta_mm
            END IF
            DO i = 1, SIZE(list)
               iatom = list(i)
               CALL build_mm_pot(qpot(iatom), 1, eta_a(0), qmmm_env%Potentials, particles_mm, &
                                 qmmm_env%mm_atom_chrg, qmmm_env%mm_atom_index, mm_cell, iatom, rcutoff, &
                                 particles_qm)
               CALL build_mm_pot(qpot(iatom), 2, alpha, qmmm_env%Potentials, particles_mm, &
                                 qmmm_env%mm_atom_chrg, qmmm_env%mm_atom_index, mm_cell, iatom, rcutoff, &
                                 particles_qm)
               ! Possibly added charges
               IF (qmmm_env%move_mm_charges .OR. qmmm_env%add_mm_charges) THEN
                  CALL build_mm_pot(qpot(iatom), 1, eta_a(0), qmmm_env%added_charges%potentials, &
                                    qmmm_env%added_charges%added_particles, &
                                    qmmm_env%added_charges%mm_atom_chrg, &
                                    qmmm_env%added_charges%mm_atom_index, mm_cell, iatom, rcutoff, &
                                    particles_qm)
                  CALL build_mm_pot(qpot(iatom), 2, alpha, qmmm_env%added_charges%potentials, &
                                    qmmm_env%added_charges%added_particles, &
                                    qmmm_env%added_charges%mm_atom_chrg, &
                                    qmmm_env%added_charges%mm_atom_index, mm_cell, iatom, rcutoff, &
                                    particles_qm)
               END IF
               pc_ener = pc_ener + qpot(iatom)*zeff
            END DO
         END DO
      CASE (do_ewald_none)
         ! Simply summing up charges with 1/R (DFTB corrected)
         nkind = SIZE(atomic_kind_set)
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
            IF (do_dftb) THEN
               NULLIFY (dftb_kind)
               CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
               CALL get_dftb_atom_param(dftb_kind, zeff=zeff, &
                                        defined=defined, eta=eta_a, natorb=natorb)
               ! use mm charge smearing for non-scc cases
               IF (.NOT. dftb_control%self_consistent) eta_a(0) = eta_mm
               IF (.NOT. defined .OR. natorb < 1) CYCLE
            ELSEIF (do_xtb) THEN
               NULLIFY (xtb_kind)
               CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
               CALL get_xtb_atom_param(xtb_kind, zeff=zeff)
               eta_a(0) = eta_mm
            END IF
            DO i = 1, SIZE(list)
               iatom = list(i)
               CALL build_mm_pot(qpot(iatom), 0, eta_a(0), qmmm_env%Potentials, particles_mm, &
                                 qmmm_env%mm_atom_chrg, qmmm_env%mm_atom_index, mm_cell, iatom, &
                                 qmmm_env%spherical_cutoff, particles_qm)
               ! Possibly added charges
               IF (qmmm_env%move_mm_charges .OR. qmmm_env%add_mm_charges) THEN
                  CALL build_mm_pot(qpot(iatom), 0, eta_a(0), qmmm_env%added_charges%potentials, &
                                    qmmm_env%added_charges%added_particles, &
                                    qmmm_env%added_charges%mm_atom_chrg, &
                                    qmmm_env%added_charges%mm_atom_index, mm_cell, iatom, &
                                    qmmm_env%spherical_cutoff, &
                                    particles_qm)
               END IF
               pc_ener = pc_ener + qpot(iatom)*zeff
            END DO
         END DO
      CASE DEFAULT
         CPABORT("Unknown Ewald type!")
      END SELECT

      ! Allocate the core Hamiltonian matrix
      CALL get_qs_env(qs_env=qs_env, ks_qmmm_env=ks_qmmm_env_loc)
      matrix_h => ks_qmmm_env_loc%matrix_h
      CALL dbcsr_allocate_matrix_set(matrix_h, 1)
      ALLOCATE (matrix_h(1)%matrix)
      CALL dbcsr_copy(matrix_h(1)%matrix, matrix_s(1)%matrix, &
                      name="QMMM HAMILTONIAN MATRIX")
      CALL dbcsr_set(matrix_h(1)%matrix, 0.0_dp)

      CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, iatom, jatom, sblock, blk)
         NULLIFY (hblock)
         CALL dbcsr_get_block_p(matrix=matrix_h(1)%matrix, &
                                row=iatom, col=jatom, block=hblock, found=found)
         CPASSERT(found)
         hblock = hblock - 0.5_dp*sblock*(qpot(iatom) + qpot(jatom))
      END DO
      CALL dbcsr_iterator_stop(iter)

      ks_qmmm_env_loc%matrix_h => matrix_h
      ks_qmmm_env_loc%pc_ener = pc_ener

      DEALLOCATE (qpot)

      ! Release Ewald environment
      CALL ewald_env_release(ewald_env)
      CALL ewald_pw_release(ewald_pw)

      CALL dbcsr_deallocate_matrix_set(matrix_s)

      CALL timestop(handle)

   END SUBROUTINE build_tb_qmmm_matrix_pc

! **************************************************************************************************
!> \brief Constructs the derivative w.r.t. 1-el DFTB hamiltonian QMMM terms
!> \param qs_env ...
!> \param qmmm_env ...
!> \param particles_mm ...
!> \param mm_cell ...
!> \param para_env ...
!> \param calc_force ...
!> \param Forces ...
!> \param Forces_added_charges ...
!> \author JGH 10.2014 [created]
! **************************************************************************************************
   SUBROUTINE deriv_tb_qmmm_matrix(qs_env, qmmm_env, particles_mm, mm_cell, para_env, &
                                   calc_force, Forces, Forces_added_charges)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qmmm_env_qm_type), POINTER                    :: qmmm_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_mm
      TYPE(cell_type), POINTER                           :: mm_cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(in), OPTIONAL                      :: calc_force
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: Forces, Forces_added_charges

      CHARACTER(len=*), PARAMETER :: routineN = 'deriv_tb_qmmm_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, blk, handle, i, iatom, ikind, &
                                                            iqm, jatom, natom, natorb, nkind, &
                                                            nspins, number_qm_atoms
      INTEGER, DIMENSION(:), POINTER                     :: list
      LOGICAL                                            :: defined, do_dftb, do_xtb, found
      REAL(KIND=dp)                                      :: fi, gmij, zeff
      REAL(KIND=dp), DIMENSION(0:3)                      :: eta_a
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mcharge, qpot
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: charges, dsblock, Forces_QM, pblock, &
                                                            sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dftb_control_type), POINTER                   :: dftb_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_qm
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_kind
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_ks_qmmm_env_type), POINTER                 :: ks_qmmm_env_loc
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)
      IF (calc_force) THEN
         NULLIFY (rho, atomic_kind_set, qs_kind_set, particles_qm)
         CALL get_qs_env(qs_env=qs_env, &
                         rho=rho, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set, &
                         ks_qmmm_env=ks_qmmm_env_loc, &
                         dft_control=dft_control, &
                         particle_set=particles_qm, &
                         natom=number_qm_atoms)
         dftb_control => dft_control%qs_control%dftb_control
         xtb_control => dft_control%qs_control%xtb_control

         IF (dft_control%qs_control%dftb) THEN
            do_dftb = .TRUE.
            do_xtb = .FALSE.
         ELSEIF (dft_control%qs_control%xtb) THEN
            do_dftb = .FALSE.
            do_xtb = .TRUE.
         ELSE
            CPABORT("TB method unknown")
         END IF

         NULLIFY (matrix_s)
         IF (do_dftb) THEN
            CALL build_dftb_overlap(qs_env, 1, matrix_s)
         ELSEIF (do_xtb) THEN
            CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, sab_orb=sab_nl)
            CALL build_overlap_matrix(ks_env, matrix_s, nderivative=1, &
                                      basis_type_a='ORB', basis_type_b='ORB', sab_nl=sab_nl)
         END IF

         CALL qs_rho_get(rho, rho_ao=matrix_p)

         nspins = dft_control%nspins
         nkind = SIZE(atomic_kind_set)
         ! Mulliken charges
         ALLOCATE (charges(number_qm_atoms, nspins))
         !
         CALL mulliken_charges(matrix_p, matrix_s(1)%matrix, para_env, charges)
         !
         ALLOCATE (mcharge(number_qm_atoms))
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
            IF (do_dftb) THEN
               CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
               CALL get_dftb_atom_param(dftb_kind, zeff=zeff)
            ELSEIF (do_xtb) THEN
               CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
               CALL get_xtb_atom_param(xtb_kind, zeff=zeff)
            END IF
            DO iatom = 1, natom
               atom_a = atomic_kind_set(ikind)%atom_list(iatom)
               mcharge(atom_a) = zeff - SUM(charges(atom_a, 1:nspins))
            END DO
         END DO
         DEALLOCATE (charges)

         ALLOCATE (qpot(number_qm_atoms))
         qpot = 0.0_dp
         ALLOCATE (Forces_QM(3, number_qm_atoms))
         Forces_QM = 0.0_dp

         ! calculate potential and forces from classical charges
         iqm = 0
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
            IF (do_dftb) THEN
               NULLIFY (dftb_kind)
               CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
               CALL get_dftb_atom_param(dftb_kind, &
                                        defined=defined, eta=eta_a, natorb=natorb)
               ! use mm charge smearing for non-scc cases
               IF (.NOT. dftb_control%self_consistent) eta_a(0) = eta_mm
               IF (.NOT. defined .OR. natorb < 1) CYCLE
            ELSEIF (do_xtb) THEN
               eta_a(0) = eta_mm
            END IF
            DO i = 1, SIZE(list)
               iatom = list(i)
               iqm = iqm + 1
               CALL build_mm_pot(qpot(iatom), 0, eta_a(0), qmmm_env%Potentials, particles_mm, &
                                 qmmm_env%mm_atom_chrg, qmmm_env%mm_atom_index, mm_cell, iatom, &
                                 qmmm_env%spherical_cutoff, particles_qm)
               CALL build_mm_dpot(mcharge(iatom), 0, eta_a(0), qmmm_env%Potentials, particles_mm, &
                                  qmmm_env%mm_atom_chrg, qmmm_env%mm_atom_index, &
                                  mm_cell, iatom, Forces, Forces_QM(:, iqm), &
                                  qmmm_env%spherical_cutoff, particles_qm)
               ! Possibly added charges
               IF (qmmm_env%move_mm_charges .OR. qmmm_env%add_mm_charges) THEN
                  CALL build_mm_pot(qpot(iatom), 0, eta_a(0), qmmm_env%added_charges%potentials, &
                                    qmmm_env%added_charges%added_particles, &
                                    qmmm_env%added_charges%mm_atom_chrg, &
                                    qmmm_env%added_charges%mm_atom_index, mm_cell, iatom, &
                                    qmmm_env%spherical_cutoff, &
                                    particles_qm)
                  CALL build_mm_dpot(mcharge(iatom), 0, eta_a(0), qmmm_env%added_charges%potentials, &
                                     qmmm_env%added_charges%added_particles, &
                                     qmmm_env%added_charges%mm_atom_chrg, &
                                     qmmm_env%added_charges%mm_atom_index, mm_cell, iatom, &
                                     Forces_added_charges, &
                                     Forces_QM(:, iqm), qmmm_env%spherical_cutoff, particles_qm)
               END IF
            END DO
         END DO

         ! Transfer QM gradients to the QM particles..
         iqm = 0
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
            IF (do_dftb) THEN
               NULLIFY (dftb_kind)
               CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
               CALL get_dftb_atom_param(dftb_kind, defined=defined, natorb=natorb)
               IF (.NOT. defined .OR. natorb < 1) CYCLE
            ELSEIF (do_xtb) THEN
               ! use all kinds
            END IF
            DO i = 1, SIZE(list)
               iqm = iqm + 1
               iatom = qmmm_env%qm_atom_index(list(i))
               particles_mm(iatom)%f(:) = particles_mm(iatom)%f(:) + Forces_QM(:, iqm)
            END DO
         END DO

         ! derivatives from qm charges
         Forces_QM = 0.0_dp
         IF (SIZE(matrix_p) == 2) THEN
            CALL dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         END IF
         !
         CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, iatom, jatom, sblock, blk)
            !
            IF (iatom == jatom) CYCLE
            !
            gmij = -0.5_dp*(qpot(iatom) + qpot(jatom))
            NULLIFY (pblock)
            CALL dbcsr_get_block_p(matrix=matrix_p(1)%matrix, &
                                   row=iatom, col=jatom, block=pblock, found=found)
            CPASSERT(found)
            DO i = 1, 3
               NULLIFY (dsblock)
               CALL dbcsr_get_block_p(matrix=matrix_s(1 + i)%matrix, &
                                      row=iatom, col=jatom, block=dsblock, found=found)
               CPASSERT(found)
               fi = -2.0_dp*gmij*SUM(pblock*dsblock)
               Forces_QM(i, iatom) = Forces_QM(i, iatom) + fi
               Forces_QM(i, jatom) = Forces_QM(i, jatom) - fi
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)
         !
         IF (SIZE(matrix_p) == 2) THEN
            CALL dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         END IF
         !
         ! Transfer QM gradients to the QM particles..
         CALL mp_sum(Forces_QM, para_env%group)
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
            DO i = 1, SIZE(list)
               iqm = list(i)
               iatom = qmmm_env%qm_atom_index(iqm)
               particles_mm(iatom)%f(:) = particles_mm(iatom)%f(:) + Forces_QM(:, iqm)
            END DO
         END DO
         !
         DEALLOCATE (mcharge)
         !
         ! MM forces will be handled directly from the QMMM module in the same way
         ! as for GPW/GAPW methods
         DEALLOCATE (Forces_QM)
         DEALLOCATE (qpot)

         CALL dbcsr_deallocate_matrix_set(matrix_s)

      END IF
      CALL timestop(handle)

   END SUBROUTINE deriv_tb_qmmm_matrix

! **************************************************************************************************
!> \brief Constructs the derivative w.r.t. 1-el DFTB hamiltonian QMMM terms
!> \param qs_env ...
!> \param qmmm_env ...
!> \param particles_mm ...
!> \param mm_cell ...
!> \param para_env ...
!> \param calc_force ...
!> \param Forces ...
!> \param Forces_added_charges ...
!> \author JGH 10.2014 [created]
! **************************************************************************************************
   SUBROUTINE deriv_tb_qmmm_matrix_pc(qs_env, qmmm_env, particles_mm, mm_cell, para_env, &
                                      calc_force, Forces, Forces_added_charges)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qmmm_env_qm_type), POINTER                    :: qmmm_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_mm
      TYPE(cell_type), POINTER                           :: mm_cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(in), OPTIONAL                      :: calc_force
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: Forces, Forces_added_charges

      CHARACTER(len=*), PARAMETER :: routineN = 'deriv_tb_qmmm_matrix_pc', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_a, blk, do_ipol, ewald_type, handle, i, iatom, ikind, imm, imp, indmm, ipot, &
         iqm, jatom, natom, natorb, nkind, nmm, nspins, number_qm_atoms
      INTEGER, DIMENSION(:), POINTER                     :: list
      LOGICAL                                            :: defined, do_dftb, do_multipoles, do_xtb, &
                                                            found
      REAL(KIND=dp)                                      :: alpha, fi, gmij, zeff
      REAL(KIND=dp), DIMENSION(0:3)                      :: eta_a
      REAL(KIND=dp), DIMENSION(2)                        :: rcutoff
      REAL(KIND=dp), DIMENSION(:), POINTER               :: charges_mm, mcharge, qpot
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: charges, dsblock, Forces_MM, Forces_QM, &
                                                            pblock, sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dftb_control_type), POINTER                   :: dftb_control
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(particle_type), DIMENSION(:), POINTER         :: atoms_mm, particles_qm
      TYPE(qmmm_pot_type), POINTER                       :: Pot
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_kind
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_ks_qmmm_env_type), POINTER                 :: ks_qmmm_env_loc
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: ewald_section, poisson_section, &
                                                            print_section
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)
      IF (calc_force) THEN
         NULLIFY (rho, atomic_kind_set, qs_kind_set, particles_qm)
         CALL get_qs_env(qs_env=qs_env, &
                         rho=rho, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set, &
                         ks_qmmm_env=ks_qmmm_env_loc, &
                         dft_control=dft_control, &
                         particle_set=particles_qm, &
                         natom=number_qm_atoms)
         dftb_control => dft_control%qs_control%dftb_control
         xtb_control => dft_control%qs_control%xtb_control

         IF (dft_control%qs_control%dftb) THEN
            do_dftb = .TRUE.
            do_xtb = .FALSE.
         ELSEIF (dft_control%qs_control%xtb) THEN
            do_dftb = .FALSE.
            do_xtb = .TRUE.
         ELSE
            CPABORT("TB method unknown")
         END IF

         NULLIFY (matrix_s)
         IF (do_dftb) THEN
            CALL build_dftb_overlap(qs_env, 1, matrix_s)
         ELSEIF (do_xtb) THEN
            CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, sab_orb=sab_nl)
            CALL build_overlap_matrix(ks_env, matrix_s, nderivative=1, &
                                      basis_type_a='ORB', basis_type_b='ORB', sab_nl=sab_nl)
         END IF
         CALL qs_rho_get(rho, rho_ao=matrix_p)

         nspins = dft_control%nspins
         nkind = SIZE(atomic_kind_set)
         ! Mulliken charges
         ALLOCATE (charges(number_qm_atoms, nspins))
         !
         CALL mulliken_charges(matrix_p, matrix_s(1)%matrix, para_env, charges)
         !
         ALLOCATE (mcharge(number_qm_atoms))
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
            IF (do_dftb) THEN
               CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
               CALL get_dftb_atom_param(dftb_kind, zeff=zeff)
            ELSEIF (do_xtb) THEN
               CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
               CALL get_xtb_atom_param(xtb_kind, zeff=zeff)
            END IF
            DO iatom = 1, natom
               atom_a = atomic_kind_set(ikind)%atom_list(iatom)
               mcharge(atom_a) = zeff - SUM(charges(atom_a, 1:nspins))
            END DO
         END DO
         DEALLOCATE (charges)

         ALLOCATE (qpot(number_qm_atoms))
         qpot = 0.0_dp
         ALLOCATE (Forces_QM(3, number_qm_atoms))
         Forces_QM = 0.0_dp

         ! Create Ewald environments
         poisson_section => section_vals_get_subs_vals(qs_env%input, "MM%POISSON")
         CALL ewald_env_create(ewald_env, para_env)
         CALL ewald_env_set(ewald_env, poisson_section=poisson_section)
         ewald_section => section_vals_get_subs_vals(poisson_section, "EWALD")
         CALL read_ewald_section(ewald_env, ewald_section)
         print_section => section_vals_get_subs_vals(qs_env%input, "PRINT%GRID_INFORMATION")
         CALL ewald_pw_create(ewald_pw, ewald_env, mm_cell, mm_cell, print_section=print_section)

         CALL ewald_env_get(ewald_env, ewald_type=ewald_type, do_multipoles=do_multipoles, do_ipol=do_ipol)
         IF (do_multipoles) CPABORT("No multipole force fields allowed in DFTB QM/MM")
         IF (do_ipol /= do_fist_pol_none) CPABORT("No polarizable force fields allowed in DFTB QM/MM")

         SELECT CASE (ewald_type)
         CASE (do_ewald_pme)
            CPABORT("PME Ewald type not implemented for DFTB/QMMM")
         CASE (do_ewald_ewald, do_ewald_spme)
            DO ipot = 1, SIZE(qmmm_env%Potentials)
               Pot => qmmm_env%Potentials(ipot)%Pot
               nmm = SIZE(Pot%mm_atom_index)
               ! get a 'clean' mm particle set
               NULLIFY (atoms_mm)
               CALL allocate_particle_set(atoms_mm, nmm)
               ALLOCATE (charges_mm(nmm))
               DO Imp = 1, nmm
                  Imm = Pot%mm_atom_index(Imp)
                  IndMM = qmmm_env%mm_atom_index(Imm)
                  atoms_mm(Imp)%r = particles_mm(IndMM)%r
                  atoms_mm(Imp)%atomic_kind => particles_mm(IndMM)%atomic_kind
                  charges_mm(Imp) = qmmm_env%mm_atom_chrg(Imm)
               END DO
               ! force array for mm atoms
               ALLOCATE (Forces_MM(3, nmm))
               Forces_MM = 0.0_dp
               IF (ewald_type == do_ewald_ewald) THEN
                  CPABORT("Ewald not implemented for DFTB/QMMM")
               ELSE IF (ewald_type == do_ewald_spme) THEN
                  ! spme electrostatic potential
                  CALL spme_potential(ewald_env, ewald_pw, mm_cell, atoms_mm, charges_mm, &
                                      particles_qm, qpot)
                  ! forces QM
                  CALL spme_forces(ewald_env, ewald_pw, mm_cell, atoms_mm, charges_mm, &
                                   particles_qm, mcharge, Forces_QM)
                  ! forces MM
                  CALL spme_forces(ewald_env, ewald_pw, mm_cell, particles_qm, mcharge, &
                                   atoms_mm, charges_mm, Forces_MM)
               END IF
               CALL deallocate_particle_set(atoms_mm)
               DEALLOCATE (charges_mm)
               ! transfer MM forces
               CALL mp_sum(Forces_MM, para_env%group)
               DO Imp = 1, nmm
                  Imm = Pot%mm_atom_index(Imp)
                  Forces(:, Imm) = Forces(:, Imm) - Forces_MM(:, Imp)
               END DO
               DEALLOCATE (Forces_MM)
            END DO

            IF (qmmm_env%move_mm_charges .OR. qmmm_env%add_mm_charges) THEN
               DO ipot = 1, SIZE(qmmm_env%added_charges%Potentials)
                  Pot => qmmm_env%added_charges%Potentials(ipot)%Pot
                  nmm = SIZE(Pot%mm_atom_index)
                  ! get a 'clean' mm particle set
                  NULLIFY (atoms_mm)
                  CALL allocate_particle_set(atoms_mm, nmm)
                  ALLOCATE (charges_mm(nmm))
                  DO Imp = 1, nmm
                     Imm = Pot%mm_atom_index(Imp)
                     IndMM = qmmm_env%added_charges%mm_atom_index(Imm)
                     atoms_mm(Imp)%r = qmmm_env%added_charges%added_particles(IndMM)%r
                     atoms_mm(Imp)%atomic_kind => qmmm_env%added_charges%added_particles(IndMM)%atomic_kind
                     charges_mm(Imp) = qmmm_env%added_charges%mm_atom_chrg(Imm)
                  END DO
                  ! force array for mm atoms
                  ALLOCATE (Forces_MM(3, nmm))
                  Forces_MM = 0.0_dp
                  IF (ewald_type == do_ewald_ewald) THEN
                     CPABORT("Ewald not implemented for DFTB/QMMM")
                  ELSE IF (ewald_type == do_ewald_spme) THEN
                     ! spme electrostatic potential
                     CALL spme_potential(ewald_env, ewald_pw, mm_cell, atoms_mm, &
                                         charges_mm, particles_qm, qpot)
                     ! forces QM
                     CALL spme_forces(ewald_env, ewald_pw, mm_cell, atoms_mm, charges_mm, &
                                      particles_qm, mcharge, Forces_QM)
                     ! forces MM
                     CALL spme_forces(ewald_env, ewald_pw, mm_cell, particles_qm, mcharge, &
                                      atoms_mm, charges_mm, Forces_MM)
                  END IF
                  CALL deallocate_particle_set(atoms_mm)
                  ! transfer MM forces
                  CALL mp_sum(Forces_MM, para_env%group)
                  DO Imp = 1, nmm
                     Imm = Pot%mm_atom_index(Imp)
                     Forces_added_charges(:, Imm) = Forces_added_charges(:, Imm) - Forces_MM(:, Imp)
                  END DO
                  DEALLOCATE (Forces_MM)
               END DO
            END IF
            CALL mp_sum(qpot, para_env%group)
            CALL mp_sum(Forces_QM, para_env%group)
            ! Add Ewald and DFTB short range corrections
            ! This is effectively using a minimum image convention!
            ! Set rcutoff to values compatible with alpha Ewald
            CALL ewald_env_get(ewald_env, rcut=rcutoff(1), alpha=alpha)
            rcutoff(2) = 0.025_dp*rcutoff(1)
            rcutoff(1) = 2.0_dp*rcutoff(1)
            nkind = SIZE(atomic_kind_set)
            iqm = 0
            DO ikind = 1, nkind
               CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
               IF (do_dftb) THEN
                  NULLIFY (dftb_kind)
                  CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
                  CALL get_dftb_atom_param(dftb_kind, &
                                           defined=defined, eta=eta_a, natorb=natorb)
                  ! use mm charge smearing for non-scc cases
                  IF (.NOT. dftb_control%self_consistent) eta_a(0) = eta_mm
                  IF (.NOT. defined .OR. natorb < 1) CYCLE
               ELSEIF (do_xtb) THEN
                  eta_a(0) = eta_mm
               END IF
               DO i = 1, SIZE(list)
                  iatom = list(i)
                  iqm = iqm + 1
                  CALL build_mm_pot(qpot(iatom), 1, eta_a(0), qmmm_env%Potentials, particles_mm, &
                                    qmmm_env%mm_atom_chrg, qmmm_env%mm_atom_index, &
                                    mm_cell, iatom, rcutoff, particles_qm)
                  CALL build_mm_dpot(mcharge(iatom), 1, eta_a(0), qmmm_env%Potentials, particles_mm, &
                                     qmmm_env%mm_atom_chrg, qmmm_env%mm_atom_index, &
                                     mm_cell, iatom, Forces, Forces_QM(:, iqm), &
                                     rcutoff, particles_qm)
                  CALL build_mm_pot(qpot(iatom), 2, alpha, qmmm_env%Potentials, particles_mm, &
                                    qmmm_env%mm_atom_chrg, qmmm_env%mm_atom_index, &
                                    mm_cell, iatom, rcutoff, particles_qm)
                  CALL build_mm_dpot(mcharge(iatom), 2, alpha, qmmm_env%Potentials, particles_mm, &
                                     qmmm_env%mm_atom_chrg, qmmm_env%mm_atom_index, &
                                     mm_cell, iatom, Forces, Forces_QM(:, iqm), &
                                     rcutoff, particles_qm)
                  ! Possibly added charges
                  IF (qmmm_env%move_mm_charges .OR. qmmm_env%add_mm_charges) THEN
                     CALL build_mm_pot(qpot(iatom), 1, eta_a(0), qmmm_env%added_charges%potentials, &
                                       qmmm_env%added_charges%added_particles, &
                                       qmmm_env%added_charges%mm_atom_chrg, &
                                       qmmm_env%added_charges%mm_atom_index, mm_cell, iatom, rcutoff, &
                                       particles_qm)
                     CALL build_mm_dpot( &
                        mcharge(iatom), 1, eta_a(0), qmmm_env%added_charges%potentials, &
                        qmmm_env%added_charges%added_particles, qmmm_env%added_charges%mm_atom_chrg, &
                        qmmm_env%added_charges%mm_atom_index, mm_cell, iatom, &
                        Forces_added_charges, Forces_QM(:, iqm), &
                        rcutoff, particles_qm)
                     CALL build_mm_pot(qpot(iatom), 2, alpha, qmmm_env%added_charges%potentials, &
                                       qmmm_env%added_charges%added_particles, &
                                       qmmm_env%added_charges%mm_atom_chrg, &
                                       qmmm_env%added_charges%mm_atom_index, mm_cell, iatom, rcutoff, &
                                       particles_qm)
                     CALL build_mm_dpot( &
                        mcharge(iatom), 2, alpha, qmmm_env%added_charges%potentials, &
                        qmmm_env%added_charges%added_particles, qmmm_env%added_charges%mm_atom_chrg, &
                        qmmm_env%added_charges%mm_atom_index, mm_cell, iatom, &
                        Forces_added_charges, Forces_QM(:, iqm), &
                        rcutoff, particles_qm)
                  END IF
               END DO
            END DO

         CASE (do_ewald_none)
            ! Simply summing up charges with 1/R (DFTB corrected)
            ! calculate potential and forces from classical charges
            iqm = 0
            DO ikind = 1, nkind
               CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
               IF (do_dftb) THEN
                  NULLIFY (dftb_kind)
                  CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
                  CALL get_dftb_atom_param(dftb_kind, &
                                           defined=defined, eta=eta_a, natorb=natorb)
                  ! use mm charge smearing for non-scc cases
                  IF (.NOT. dftb_control%self_consistent) eta_a(0) = eta_mm
                  IF (.NOT. defined .OR. natorb < 1) CYCLE
               ELSEIF (do_xtb) THEN
                  eta_a(0) = eta_mm
               END IF
               DO i = 1, SIZE(list)
                  iatom = list(i)
                  iqm = iqm + 1
                  CALL build_mm_pot(qpot(iatom), 0, eta_a(0), qmmm_env%Potentials, particles_mm, &
                                    qmmm_env%mm_atom_chrg, qmmm_env%mm_atom_index, mm_cell, iatom, &
                                    qmmm_env%spherical_cutoff, particles_qm)
                  CALL build_mm_dpot(mcharge(iatom), 0, eta_a(0), qmmm_env%Potentials, particles_mm, &
                                     qmmm_env%mm_atom_chrg, qmmm_env%mm_atom_index, &
                                     mm_cell, iatom, Forces, Forces_QM(:, iqm), &
                                     qmmm_env%spherical_cutoff, particles_qm)
                  ! Possibly added charges
                  IF (qmmm_env%move_mm_charges .OR. qmmm_env%add_mm_charges) THEN
                     CALL build_mm_pot(qpot(iatom), 0, eta_a(0), qmmm_env%added_charges%potentials, &
                                       qmmm_env%added_charges%added_particles, &
                                       qmmm_env%added_charges%mm_atom_chrg, &
                                       qmmm_env%added_charges%mm_atom_index, &
                                       mm_cell, iatom, qmmm_env%spherical_cutoff, &
                                       particles_qm)
                     CALL build_mm_dpot(mcharge(iatom), 0, eta_a(0), qmmm_env%added_charges%potentials, &
                                        qmmm_env%added_charges%added_particles, &
                                        qmmm_env%added_charges%mm_atom_chrg, &
                                        qmmm_env%added_charges%mm_atom_index, mm_cell, iatom, &
                                        Forces_added_charges, &
                                        Forces_QM(:, iqm), qmmm_env%spherical_cutoff, particles_qm)
                  END IF
               END DO
            END DO
         CASE DEFAULT
            CPABORT("Unknown Ewald type!")
         END SELECT

         ! Transfer QM gradients to the QM particles..
         iqm = 0
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
            IF (do_dftb) THEN
               NULLIFY (dftb_kind)
               CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
               CALL get_dftb_atom_param(dftb_kind, defined=defined, natorb=natorb)
               IF (.NOT. defined .OR. natorb < 1) CYCLE
            ELSEIF (do_xtb) THEN
               !
            END IF
            DO i = 1, SIZE(list)
               iqm = iqm + 1
               iatom = qmmm_env%qm_atom_index(list(i))
               particles_mm(iatom)%f(:) = particles_mm(iatom)%f(:) + Forces_QM(:, iqm)
            END DO
         END DO

         ! derivatives from qm charges
         Forces_QM = 0.0_dp
         IF (SIZE(matrix_p) == 2) THEN
            CALL dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         END IF
         !
         CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, iatom, jatom, sblock, blk)
            !
            IF (iatom == jatom) CYCLE
            !
            gmij = -0.5_dp*(qpot(iatom) + qpot(jatom))
            NULLIFY (pblock)
            CALL dbcsr_get_block_p(matrix=matrix_p(1)%matrix, &
                                   row=iatom, col=jatom, block=pblock, found=found)
            CPASSERT(found)
            DO i = 1, 3
               NULLIFY (dsblock)
               CALL dbcsr_get_block_p(matrix=matrix_s(1 + i)%matrix, &
                                      row=iatom, col=jatom, block=dsblock, found=found)
               CPASSERT(found)
               fi = -2.0_dp*gmij*SUM(pblock*dsblock)
               Forces_QM(i, iatom) = Forces_QM(i, iatom) + fi
               Forces_QM(i, jatom) = Forces_QM(i, jatom) - fi
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)
         !
         IF (SIZE(matrix_p) == 2) THEN
            CALL dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         END IF
         !
         ! Transfer QM gradients to the QM particles..
         CALL mp_sum(Forces_QM, para_env%group)
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list)
            DO i = 1, SIZE(list)
               iqm = list(i)
               iatom = qmmm_env%qm_atom_index(iqm)
               particles_mm(iatom)%f(:) = particles_mm(iatom)%f(:) + Forces_QM(:, iqm)
            END DO
         END DO
         !
         DEALLOCATE (mcharge)
         !
         ! MM forces will be handled directly from the QMMM module in the same way
         ! as for GPW/GAPW methods
         DEALLOCATE (Forces_QM)
         DEALLOCATE (qpot)

         ! Release Ewald environment
         CALL ewald_env_release(ewald_env)
         CALL ewald_pw_release(ewald_pw)

         CALL dbcsr_deallocate_matrix_set(matrix_s)

      END IF

      CALL timestop(handle)

   END SUBROUTINE deriv_tb_qmmm_matrix_pc

! **************************************************************************************************
!> \brief ...
!> \param qpot ...
!> \param pot_type ...
!> \param qm_alpha ...
!> \param potentials ...
!> \param particles_mm ...
!> \param mm_charges ...
!> \param mm_atom_index ...
!> \param mm_cell ...
!> \param IndQM ...
!> \param qmmm_spherical_cutoff ...
!> \param particles_qm ...
! **************************************************************************************************
   SUBROUTINE build_mm_pot(qpot, pot_type, qm_alpha, potentials, &
                           particles_mm, mm_charges, mm_atom_index, mm_cell, IndQM, &
                           qmmm_spherical_cutoff, particles_qm)

      REAL(KIND=dp), INTENT(INOUT)                       :: qpot
      INTEGER, INTENT(IN)                                :: pot_type
      REAL(KIND=dp), INTENT(IN)                          :: qm_alpha
      TYPE(qmmm_pot_p_type), DIMENSION(:), POINTER       :: potentials
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_mm
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mm_charges
      INTEGER, DIMENSION(:), POINTER                     :: mm_atom_index
      TYPE(cell_type), POINTER                           :: mm_cell
      INTEGER, INTENT(IN)                                :: IndQM
      REAL(KIND=dp), INTENT(IN)                          :: qmmm_spherical_cutoff(2)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_qm

      CHARACTER(len=*), PARAMETER                        :: routineN = 'build_mm_pot'
      REAL(KIND=dp), PARAMETER                           :: qsmall = 1.0e-15_dp

      INTEGER                                            :: handle, Imm, Imp, IndMM, Ipot
      REAL(KIND=dp)                                      :: dr, qeff, rt1, rt2, rt3, &
                                                            sph_chrg_factor, sr
      REAL(KIND=dp), DIMENSION(3)                        :: r_pbc, rij
      TYPE(qmmm_pot_type), POINTER                       :: Pot

      CALL timeset(routineN, handle)
      ! Loop Over MM atoms
      ! Loop over Pot stores atoms with the same charge
      MainLoopPot: DO Ipot = 1, SIZE(Potentials)
         Pot => Potentials(Ipot)%Pot
         ! Loop over atoms belonging to this type
         LoopMM: DO Imp = 1, SIZE(Pot%mm_atom_index)
            Imm = Pot%mm_atom_index(Imp)
            IndMM = mm_atom_index(Imm)
            r_pbc = pbc(particles_mm(IndMM)%r - particles_qm(IndQM)%r, mm_cell)
            rt1 = r_pbc(1)
            rt2 = r_pbc(2)
            rt3 = r_pbc(3)
            rij = (/rt1, rt2, rt3/)
            dr = SQRT(SUM(rij**2))
            qeff = mm_charges(Imm)
            ! Computes the screening factor for the spherical cutoff (if defined)
            IF (qmmm_spherical_cutoff(1) > 0.0_dp) THEN
               CALL spherical_cutoff_factor(qmmm_spherical_cutoff, rij, sph_chrg_factor)
               qeff = qeff*sph_chrg_factor
            END IF
            IF (ABS(qeff) <= qsmall) CYCLE
            IF (dr > rtiny) THEN
               IF (pot_type == 0) THEN
                  sr = gamma_rab_sr(dr, qm_alpha, eta_mm, 0.0_dp)
                  qpot = qpot + qeff*(1.0_dp/dr - sr)
               ELSE IF (pot_type == 1) THEN
                  sr = gamma_rab_sr(dr, qm_alpha, eta_mm, 0.0_dp)
                  qpot = qpot - qeff*sr
               ELSE IF (pot_type == 2) THEN
                  sr = erfc(qm_alpha*dr)/dr
                  qpot = qpot + qeff*sr
               ELSE
                  CPABORT("")
               END IF
            END IF
         END DO LoopMM
      END DO MainLoopPot
      CALL timestop(handle)
   END SUBROUTINE build_mm_pot

! **************************************************************************************************
!> \brief ...
!> \param qcharge ...
!> \param pot_type ...
!> \param qm_alpha ...
!> \param potentials ...
!> \param particles_mm ...
!> \param mm_charges ...
!> \param mm_atom_index ...
!> \param mm_cell ...
!> \param IndQM ...
!> \param forces ...
!> \param forces_qm ...
!> \param qmmm_spherical_cutoff ...
!> \param particles_qm ...
! **************************************************************************************************
   SUBROUTINE build_mm_dpot(qcharge, pot_type, qm_alpha, potentials, &
                            particles_mm, mm_charges, mm_atom_index, mm_cell, IndQM, &
                            forces, forces_qm, qmmm_spherical_cutoff, particles_qm)

      REAL(KIND=dp), INTENT(IN)                          :: qcharge
      INTEGER, INTENT(IN)                                :: pot_type
      REAL(KIND=dp), INTENT(IN)                          :: qm_alpha
      TYPE(qmmm_pot_p_type), DIMENSION(:), POINTER       :: potentials
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_mm
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mm_charges
      INTEGER, DIMENSION(:), POINTER                     :: mm_atom_index
      TYPE(cell_type), POINTER                           :: mm_cell
      INTEGER, INTENT(IN)                                :: IndQM
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: forces
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: forces_qm
      REAL(KIND=dp), INTENT(IN)                          :: qmmm_spherical_cutoff(2)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_qm

      CHARACTER(len=*), PARAMETER                        :: routineN = 'build_mm_dpot'
      REAL(KIND=dp), PARAMETER                           :: qsmall = 1.0e-15_dp

      INTEGER                                            :: handle, Imm, Imp, IndMM, Ipot
      REAL(KIND=dp)                                      :: dr, drm, drp, dsr, fsr, qeff, rt1, rt2, &
                                                            rt3, sph_chrg_factor
      REAL(KIND=dp), DIMENSION(3)                        :: force_ab, r_pbc, rij
      TYPE(qmmm_pot_type), POINTER                       :: Pot

      CALL timeset(routineN, handle)
      ! Loop Over MM atoms
      ! Loop over Pot stores atoms with the same charge
      MainLoopPot: DO Ipot = 1, SIZE(Potentials)
         Pot => Potentials(Ipot)%Pot
         ! Loop over atoms belonging to this type
         LoopMM: DO Imp = 1, SIZE(Pot%mm_atom_index)
            Imm = Pot%mm_atom_index(Imp)
            IndMM = mm_atom_index(Imm)
            r_pbc = pbc(particles_mm(IndMM)%r - particles_qm(IndQM)%r, mm_cell)
            rt1 = r_pbc(1)
            rt2 = r_pbc(2)
            rt3 = r_pbc(3)
            rij = (/rt1, rt2, rt3/)
            dr = SQRT(SUM(rij**2))
            qeff = mm_charges(Imm)
            ! Computes the screening factor for the spherical cutoff (if defined)
            ! We neglect derivative of cutoff function for gradients!!!
            IF (qmmm_spherical_cutoff(1) > 0.0_dp) THEN
               CALL spherical_cutoff_factor(qmmm_spherical_cutoff, rij, sph_chrg_factor)
               qeff = qeff*sph_chrg_factor
            END IF
            IF (ABS(qeff) <= qsmall) CYCLE
            IF (dr > rtiny) THEN
               drp = dr + ddrmm
               drm = dr - ddrmm
               IF (pot_type == 0) THEN
                  dsr = 0.5_dp*(gamma_rab_sr(drp, qm_alpha, eta_mm, 0.0_dp) - &
                                gamma_rab_sr(drm, qm_alpha, eta_mm, 0.0_dp))/ddrmm
                  fsr = qeff*qcharge*(-1.0_dp/(dr*dr) - dsr)
               ELSE IF (pot_type == 1) THEN
                  dsr = 0.5_dp*(gamma_rab_sr(drp, qm_alpha, eta_mm, 0.0_dp) - &
                                gamma_rab_sr(drm, qm_alpha, eta_mm, 0.0_dp))/ddrmm
                  fsr = -qeff*qcharge*dsr
               ELSE IF (pot_type == 2) THEN
                  dsr = 0.5_dp*(erfc(qm_alpha*drp)/drp - erfc(qm_alpha*drm)/drm)/ddrmm
                  fsr = qeff*qcharge*dsr
               ELSE
                  CPABORT("")
               END IF
               force_ab = -fsr*rij/dr
            ELSE
               force_ab = 0.0_dp
            END IF
            ! The array of QM forces are really the forces
            forces_qm(:) = forces_qm(:) - force_ab
            ! The one of MM atoms are instead gradients
            forces(:, Imm) = forces(:, Imm) - force_ab
         END DO LoopMM
      END DO MainLoopPot

      CALL timestop(handle)

   END SUBROUTINE build_mm_dpot

END MODULE qmmm_tb_methods

