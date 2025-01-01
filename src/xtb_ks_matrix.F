!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of KS matrix in xTB
!>        Reference: Stefan Grimme, Christoph Bannwarth, Philip Shushkov
!>                   JCTC 13, 1989-2009, (2017)
!>                   DOI: 10.1021/acs.jctc.7b00118
!> \author JGH
! **************************************************************************************************
MODULE xtb_ks_matrix
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_dot,&
                                              dbcsr_multiply,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE efield_tb_methods,               ONLY: efield_tb_matrix
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE mulliken,                        ONLY: ao_charges
   USE particle_types,                  ONLY: particle_type
   USE qs_charge_mixing,                ONLY: charge_mixing
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE xtb_coulomb,                     ONLY: build_xtb_coulomb
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_ks_matrix'

   PUBLIC :: build_xtb_ks_matrix

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param just_energy ...
!> \param ext_ks_matrix ...
! **************************************************************************************************
   SUBROUTINE build_xtb_ks_matrix(qs_env, calculate_forces, just_energy, ext_ks_matrix)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: ext_ks_matrix

      INTEGER                                            :: gfn_type
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      gfn_type = dft_control%qs_control%xtb_control%gfn_type

      SELECT CASE (gfn_type)
      CASE (0)
         CPASSERT(.NOT. PRESENT(ext_ks_matrix))
         CALL build_gfn0_xtb_ks_matrix(qs_env, calculate_forces, just_energy)
      CASE (1)
         CALL build_gfn1_xtb_ks_matrix(qs_env, calculate_forces, just_energy, ext_ks_matrix)
      CASE (2)
         CPABORT("gfn_type = 2 not yet available")
      CASE DEFAULT
         CPABORT("Unknown gfn_type")
      END SELECT

   END SUBROUTINE build_xtb_ks_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param just_energy ...
! **************************************************************************************************
   SUBROUTINE build_gfn0_xtb_ks_matrix(qs_env, calculate_forces, just_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy

      CHARACTER(len=*), PARAMETER :: routineN = 'build_gfn0_xtb_ks_matrix'

      INTEGER                                            :: handle, img, iounit, ispin, natom, nimg, &
                                                            nspins
      REAL(KIND=dp)                                      :: pc_ener, qmmm_el
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p1, mo_derivs
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix, matrix_h
      TYPE(dbcsr_type), POINTER                          :: mo_coeff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: scf_section

      CALL timeset(routineN, handle)

      MARK_USED(calculate_forces)

      NULLIFY (dft_control, logger, scf_section, ks_env, ks_matrix, rho, &
               energy)
      CPASSERT(ASSOCIATED(qs_env))

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_h_kp=matrix_h, &
                      para_env=para_env, &
                      ks_env=ks_env, &
                      matrix_ks_kp=ks_matrix, &
                      energy=energy)

      energy%hartree = 0.0_dp
      energy%qmmm_el = 0.0_dp

      scf_section => section_vals_get_subs_vals(qs_env%input, "DFT%SCF")
      nspins = dft_control%nspins
      nimg = dft_control%nimages
      CPASSERT(ASSOCIATED(matrix_h))
      CPASSERT(SIZE(ks_matrix) > 0)

      DO ispin = 1, nspins
         DO img = 1, nimg
            ! copy the core matrix into the fock matrix
            CALL dbcsr_copy(ks_matrix(ispin, img)%matrix, matrix_h(1, img)%matrix)
         END DO
      END DO

      IF (qs_env%qmmm) THEN
         CPABORT("gfn0 QMMM NYA")
         CALL get_qs_env(qs_env=qs_env, rho=rho, natom=natom)
         CPASSERT(SIZE(ks_matrix, 2) == 1)
         DO ispin = 1, nspins
            ! If QM/MM sumup the 1el Hamiltonian
            CALL dbcsr_add(ks_matrix(ispin, 1)%matrix, qs_env%ks_qmmm_env%matrix_h(1)%matrix, &
                           1.0_dp, 1.0_dp)
            CALL qs_rho_get(rho, rho_ao=matrix_p1)
            ! Compute QM/MM Energy
            CALL dbcsr_dot(qs_env%ks_qmmm_env%matrix_h(1)%matrix, &
                           matrix_p1(ispin)%matrix, qmmm_el)
            energy%qmmm_el = energy%qmmm_el + qmmm_el
         END DO
         pc_ener = qs_env%ks_qmmm_env%pc_ener
         energy%qmmm_el = energy%qmmm_el + pc_ener
      END IF

      energy%total = energy%core + energy%eeq + energy%efield + energy%qmmm_el + &
                     energy%repulsive + energy%dispersion + energy%kTS

      iounit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DETAILED_ENERGY", &
                                    extension=".scfLog")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,(T9,A,T60,F20.10))") &
            "Repulsive pair potential energy:               ", energy%repulsive, &
            "SRB Correction energy:                         ", energy%srb, &
            "Zeroth order Hamiltonian energy:               ", energy%core, &
            "Charge equilibration energy:                   ", energy%eeq, &
            "London dispersion energy:                      ", energy%dispersion
         IF (dft_control%qs_control%xtb_control%do_nonbonded) &
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
            "Correction for nonbonded interactions:         ", energy%xtb_nonbonded
         IF (ABS(energy%efield) > 1.e-12_dp) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "Electric field interaction energy:             ", energy%efield
         END IF
         IF (qs_env%qmmm) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "QM/MM Electrostatic energy:                    ", energy%qmmm_el
         END IF
      END IF
      CALL cp_print_key_finished_output(iounit, logger, scf_section, &
                                        "PRINT%DETAILED_ENERGY")
      ! here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C
      IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
         CPASSERT(SIZE(ks_matrix, 2) == 1)
         BLOCK
            TYPE(mo_set_type), DIMENSION(:), POINTER         :: mo_array
            CALL get_qs_env(qs_env, mo_derivs=mo_derivs, mos=mo_array)
            DO ispin = 1, SIZE(mo_derivs)
               CALL get_mo_set(mo_set=mo_array(ispin), mo_coeff_b=mo_coeff)
               IF (.NOT. mo_array(ispin)%use_mo_coeff_b) THEN
                  CPABORT("")
               END IF
               CALL dbcsr_multiply('n', 'n', 1.0_dp, ks_matrix(ispin, 1)%matrix, mo_coeff, &
                                   0.0_dp, mo_derivs(ispin)%matrix)
            END DO
         END BLOCK
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_gfn0_xtb_ks_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param just_energy ...
!> \param ext_ks_matrix ...
! **************************************************************************************************
   SUBROUTINE build_gfn1_xtb_ks_matrix(qs_env, calculate_forces, just_energy, ext_ks_matrix)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: ext_ks_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'build_gfn1_xtb_ks_matrix'

      INTEGER                                            :: atom_a, handle, iatom, ikind, img, &
                                                            iounit, is, ispin, na, natom, natorb, &
                                                            nimg, nkind, ns, nsgf, nspins
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(5)                              :: occ
      LOGICAL                                            :: do_efield, pass_check
      REAL(KIND=dp)                                      :: achrg, chmax, pc_ener, qmmm_el
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mcharge
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aocg, charges
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p1, mo_derivs, p_matrix
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix, matrix_h, matrix_p, matrix_s
      TYPE(dbcsr_type), POINTER                          :: mo_coeff, s_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: scf_section
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, logger, scf_section, matrix_p, particle_set, ks_env, &
               ks_matrix, rho, energy)
      CPASSERT(ASSOCIATED(qs_env))

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_h_kp=matrix_h, &
                      para_env=para_env, &
                      ks_env=ks_env, &
                      matrix_ks_kp=ks_matrix, &
                      rho=rho, &
                      energy=energy)

      IF (PRESENT(ext_ks_matrix)) THEN
         ! remap pointer to allow for non-kpoint external ks matrix
         ! ext_ks_matrix is used in linear response code
         ns = SIZE(ext_ks_matrix)
         ks_matrix(1:ns, 1:1) => ext_ks_matrix(1:ns)
      END IF

      energy%hartree = 0.0_dp
      energy%qmmm_el = 0.0_dp
      energy%efield = 0.0_dp

      scf_section => section_vals_get_subs_vals(qs_env%input, "DFT%SCF")
      nspins = dft_control%nspins
      nimg = dft_control%nimages
      CPASSERT(ASSOCIATED(matrix_h))
      CPASSERT(ASSOCIATED(rho))
      CPASSERT(SIZE(ks_matrix) > 0)

      DO ispin = 1, nspins
         DO img = 1, nimg
            ! copy the core matrix into the fock matrix
            CALL dbcsr_copy(ks_matrix(ispin, img)%matrix, matrix_h(1, img)%matrix)
         END DO
      END DO

      IF (dft_control%apply_period_efield .OR. dft_control%apply_efield .OR. &
          dft_control%apply_efield_field) THEN
         do_efield = .TRUE.
      ELSE
         do_efield = .FALSE.
      END IF

      IF (dft_control%qs_control%xtb_control%coulomb_interaction .OR. do_efield) THEN
         ! Mulliken charges
         CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, matrix_s_kp=matrix_s)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
         natom = SIZE(particle_set)
         ALLOCATE (mcharge(natom), charges(natom, 5))
         charges = 0.0_dp
         nkind = SIZE(atomic_kind_set)
         CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)
         ALLOCATE (aocg(nsgf, natom))
         aocg = 0.0_dp
         IF (nimg > 1) THEN
            CALL ao_charges(matrix_p, matrix_s, aocg, para_env)
         ELSE
            p_matrix => matrix_p(:, 1)
            s_matrix => matrix_s(1, 1)%matrix
            CALL ao_charges(p_matrix, s_matrix, aocg, para_env)
         END IF
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=na)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, occupation=occ)
            DO iatom = 1, na
               atom_a = atomic_kind_set(ikind)%atom_list(iatom)
               charges(atom_a, :) = REAL(occ(:), KIND=dp)
               DO is = 1, natorb
                  ns = lao(is) + 1
                  charges(atom_a, ns) = charges(atom_a, ns) - aocg(is, atom_a)
               END DO
            END DO
         END DO
         DEALLOCATE (aocg)
         ! charge mixing
         IF (dft_control%qs_control%do_ls_scf) THEN
            !
         ELSE
            CALL get_qs_env(qs_env=qs_env, scf_env=scf_env)
            CALL charge_mixing(scf_env%mixing_method, scf_env%mixing_store, &
                               charges, para_env, scf_env%iter_count)
         END IF

         DO iatom = 1, natom
            mcharge(iatom) = SUM(charges(iatom, :))
         END DO
      END IF

      IF (dft_control%qs_control%xtb_control%coulomb_interaction) THEN
         CALL build_xtb_coulomb(qs_env, ks_matrix, rho, charges, mcharge, energy, &
                                calculate_forces, just_energy)
      END IF

      IF (do_efield) THEN
         CALL efield_tb_matrix(qs_env, ks_matrix, rho, mcharge, energy, calculate_forces, just_energy)
      END IF

      IF (dft_control%qs_control%xtb_control%coulomb_interaction) THEN
         IF (dft_control%qs_control%xtb_control%check_atomic_charges) THEN
            pass_check = .TRUE.
            DO ikind = 1, nkind
               CALL get_atomic_kind(atomic_kind_set(ikind), natom=na)
               CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
               CALL get_xtb_atom_param(xtb_kind, chmax=chmax)
               DO iatom = 1, na
                  atom_a = atomic_kind_set(ikind)%atom_list(iatom)
                  achrg = mcharge(atom_a)
                  IF (ABS(achrg) > chmax) THEN
                     IF (iounit > 0) THEN
                        WRITE (iounit, "(A,A,I3,I6,A,F4.2,A,F6.2)") " Charge outside chemical range:", &
                           "  Kind Atom=", ikind, atom_a, "   Limit=", chmax, "  Charge=", achrg
                     END IF
                     pass_check = .FALSE.
                  END IF
               END DO
            END DO
            IF (.NOT. pass_check) THEN
               CALL cp_warn(__LOCATION__, "Atomic charges outside chemical range were detected."// &
                            " Switch-off CHECK_ATOMIC_CHARGES keyword in the &xTB section"// &
                            " if you want to force to continue the calculation.")
               CPABORT("xTB Charges")
            END IF
         END IF
      END IF

      IF (dft_control%qs_control%xtb_control%coulomb_interaction .OR. do_efield) THEN
         DEALLOCATE (mcharge, charges)
      END IF

      IF (qs_env%qmmm) THEN
         CPASSERT(SIZE(ks_matrix, 2) == 1)
         DO ispin = 1, nspins
            ! If QM/MM sumup the 1el Hamiltonian
            CALL dbcsr_add(ks_matrix(ispin, 1)%matrix, qs_env%ks_qmmm_env%matrix_h(1)%matrix, &
                           1.0_dp, 1.0_dp)
            CALL qs_rho_get(rho, rho_ao=matrix_p1)
            ! Compute QM/MM Energy
            CALL dbcsr_dot(qs_env%ks_qmmm_env%matrix_h(1)%matrix, &
                           matrix_p1(ispin)%matrix, qmmm_el)
            energy%qmmm_el = energy%qmmm_el + qmmm_el
         END DO
         pc_ener = qs_env%ks_qmmm_env%pc_ener
         energy%qmmm_el = energy%qmmm_el + pc_ener
      END IF

      energy%total = energy%core + energy%hartree + energy%efield + energy%qmmm_el + &
                     energy%repulsive + energy%dispersion + energy%dftb3 + energy%kTS

      iounit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DETAILED_ENERGY", &
                                    extension=".scfLog")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,(T9,A,T60,F20.10))") &
            "Repulsive pair potential energy:               ", energy%repulsive, &
            "Zeroth order Hamiltonian energy:               ", energy%core, &
            "Charge fluctuation energy:                     ", energy%hartree, &
            "London dispersion energy:                      ", energy%dispersion
         IF (dft_control%qs_control%xtb_control%xb_interaction) &
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
            "Correction for halogen bonding:                ", energy%xtb_xb_inter
         IF (dft_control%qs_control%xtb_control%do_nonbonded) &
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
            "Correction for nonbonded interactions:         ", energy%xtb_nonbonded
         IF (ABS(energy%efield) > 1.e-12_dp) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "Electric field interaction energy:             ", energy%efield
         END IF
         WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
            "DFTB3 3rd Order Energy Correction              ", energy%dftb3
         IF (qs_env%qmmm) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "QM/MM Electrostatic energy:                    ", energy%qmmm_el
         END IF
      END IF
      CALL cp_print_key_finished_output(iounit, logger, scf_section, &
                                        "PRINT%DETAILED_ENERGY")
      ! here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C
      IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
         CPASSERT(SIZE(ks_matrix, 2) == 1)
         BLOCK
            TYPE(mo_set_type), DIMENSION(:), POINTER         :: mo_array
            CALL get_qs_env(qs_env, mo_derivs=mo_derivs, mos=mo_array)
            DO ispin = 1, SIZE(mo_derivs)
               CALL get_mo_set(mo_set=mo_array(ispin), mo_coeff_b=mo_coeff)
               IF (.NOT. mo_array(ispin)%use_mo_coeff_b) THEN
                  CPABORT("")
               END IF
               CALL dbcsr_multiply('n', 'n', 1.0_dp, ks_matrix(ispin, 1)%matrix, mo_coeff, &
                                   0.0_dp, mo_derivs(ispin)%matrix)
            END DO
         END BLOCK
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_gfn1_xtb_ks_matrix

END MODULE xtb_ks_matrix

