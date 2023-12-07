!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE manybody_deepmd

   USE cp_log_handling, ONLY: cp_logger_get_default_io_unit
   USE atomic_kind_types, ONLY: atomic_kind_type
   USE bibliography, ONLY: Han2018, &
                           Zeng2023, &
                           cite_reference
   USE cell_types, ONLY: cell_type
   USE fist_nonbond_env_types, ONLY: fist_nonbond_env_get, &
                                     fist_nonbond_env_set, &
                                     fist_nonbond_env_type, &
                                     deepmd_data_type
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_para_env_type
   USE pair_potential_types, ONLY: pair_potential_pp_type, &
                                   pair_potential_single_type, &
                                   deepmd_pot_type, &
                                   deepmd_type
   USE particle_types, ONLY: particle_type
   USE physcon, ONLY: angstrom, &
                      evolt
#ifdef __DEEPMD
   USE ISO_C_BINDING, ONLY: C_PTR, C_LOC
   USE deepmd_wrapper, ONLY: dp_deep_pot, &
                             dp_deep_pot_compute
#endif

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC deepmd_energy_store_force_virial, deepmd_add_force_virial

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_deepmd'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param potparm ...
!> \param fist_nonbond_env ...
!> \param pot_deepmd ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE deepmd_energy_store_force_virial(particle_set, cell, atomic_kind_set, potparm, fist_nonbond_env, &
                                               pot_deepmd, para_env)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      REAL(kind=dp)                                      :: pot_deepmd
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deepmd_energy_store_force_virial'

#ifdef __DEEPMD
      INTEGER                                  :: i, iat, iat_use, ikind, &
                                                  jkind, n_atoms, n_atoms_use, &
                                                  output_unit
      LOGICAL, ALLOCATABLE                     :: use_atom(:)
      REAL(kind=dp)                            :: lattice(3, 3)
      TYPE(pair_potential_single_type), &
         POINTER                               :: pot
      TYPE(deepmd_data_type), POINTER          :: deepmd_data
      TYPE(deepmd_pot_type), POINTER           :: deepmd
      REAL(kind=dp), ALLOCATABLE               :: dpmd_force(:), dpmd_virial(:), &
                                                  dpmd_atomic_energy(:), dpmd_atomic_virial(:), &
                                                  dpmd_cell(:), dpmd_coord(:)
      INTEGER, ALLOCATABLE                     :: dpmd_atype(:), use_atom_type(:)
      INTEGER                                  :: dpmd_natoms

#endif
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
#ifndef __DEEPMD
      CPABORT("In order to use DeePMD-kit you need to download "// &
              "and install libdeepmd_c (check CP2K manual for details)")
      MARK_USED(particle_set)
      MARK_USED(cell)
      MARK_USED(atomic_kind_set)
      MARK_USED(potparm)
      MARK_USED(fist_nonbond_env)
      MARK_USED(pot_deepmd)
      MARK_USED(para_env)
#else
      n_atoms = SIZE(particle_set)
      ALLOCATE (use_atom(n_atoms))
      ALLOCATE (use_atom_type(n_atoms))
      use_atom = .FALSE.
      use_atom_type = 0

      NULLIFY (deepmd)

      DO ikind = 1, SIZE(atomic_kind_set)
         DO jkind = 1, SIZE(atomic_kind_set)
            pot => potparm%pot(ikind, jkind)%pot
            ! ensure that each atom is only used once
            IF (ikind /= jkind) CYCLE
            DO i = 1, SIZE(pot%type)
               IF (pot%type(i) /= deepmd_type) CYCLE
               DO iat = 1, n_atoms
                  IF (particle_set(iat)%atomic_kind%kind_number == ikind .OR. &
                      particle_set(iat)%atomic_kind%kind_number == jkind) THEN
                     use_atom(iat) = .TRUE.
                     use_atom_type(iat) = pot%set(i)%deepmd%atom_deepmd_type
                  END IF
               END DO ! iat
            END DO ! i
         END DO ! jkind
      END DO ! ikind

      n_atoms_use = COUNT(use_atom)
      dpmd_natoms = n_atoms_use
      ALLOCATE (dpmd_cell(9), dpmd_atype(dpmd_natoms), dpmd_coord(dpmd_natoms*3), &
                dpmd_force(dpmd_natoms*3), dpmd_virial(9), &
                dpmd_atomic_energy(dpmd_natoms), dpmd_atomic_virial(dpmd_natoms*9))

      iat_use = 0
      DO iat = 1, n_atoms
         IF (.NOT. use_atom(iat)) CYCLE
         iat_use = iat_use + 1
         dpmd_coord((iat_use - 1)*3 + 1:(iat_use - 1)*3 + 3) = particle_set(iat)%r*angstrom
         dpmd_atype(iat_use) = use_atom_type(iat)
      END DO
      IF (iat_use > 0) THEN
         CALL cite_reference(Han2018)
         CALL cite_reference(Zeng2023)
      END IF
      output_unit = cp_logger_get_default_io_unit()
      lattice = cell%hmat*angstrom

      ! change matrix to one d array
      DO i = 1, 3
         dpmd_cell((i - 1)*3 + 1:(i - 1)*3 + 3) = lattice(:, i)
      END DO

      ! get deepmd_data to save force, virial info
      CALL fist_nonbond_env_get(fist_nonbond_env, deepmd_data=deepmd_data)
      IF (.NOT. ASSOCIATED(deepmd_data)) THEN
         ALLOCATE (deepmd_data)
         CALL fist_nonbond_env_set(fist_nonbond_env, deepmd_data=deepmd_data)
         NULLIFY (deepmd_data%use_indices, deepmd_data%force)
         deepmd_data%model = dp_deep_pot(pot%set(1)%deepmd%deepmd_file_name)
      END IF
      IF (ASSOCIATED(deepmd_data%force)) THEN
         IF (SIZE(deepmd_data%force, 2) /= n_atoms_use) THEN
            DEALLOCATE (deepmd_data%force, deepmd_data%use_indices)
         END IF
      END IF
      IF (.NOT. ASSOCIATED(deepmd_data%force)) THEN
         ALLOCATE (deepmd_data%force(3, n_atoms_use))
         ALLOCATE (deepmd_data%use_indices(n_atoms_use))
      END IF

      CALL dp_deep_pot_compute(dp=deepmd_data%model, &
                               natom=dpmd_natoms, &
                               coord=dpmd_coord, &
                               atype=dpmd_atype, &
                               cell=dpmd_cell, &
                               energy=pot_deepmd, &
                               force=dpmd_force, &
                               virial=dpmd_virial, &
                               atomic_energy=dpmd_atomic_energy, &
                               atomic_virial=dpmd_atomic_virial)

      ! convert units
      pot_deepmd = pot_deepmd/evolt
      dpmd_force = dpmd_force/(evolt/angstrom)
      dpmd_virial = dpmd_virial/evolt

      ! account for double counting from multiple MPI processes
      IF (PRESENT(para_env)) THEN
         pot_deepmd = pot_deepmd/REAL(para_env%num_pe, dp)
         dpmd_force = dpmd_force/REAL(para_env%num_pe, dp)
         dpmd_virial = dpmd_virial/REAL(para_env%num_pe, dp)
      END IF

      ! save force, virial info
      iat_use = 0
      DO iat = 1, n_atoms
         IF (use_atom(iat)) THEN
            iat_use = iat_use + 1
            deepmd_data%use_indices(iat_use) = iat
            deepmd_data%force(1:3, iat_use) = dpmd_force((iat_use - 1)*3 + 1:(iat_use - 1)*3 + 3)
         END IF
      END DO
      DO i = 1, 3
         deepmd_data%virial(1:3, i) = dpmd_virial((i - 1)*3 + 1:(i - 1)*3 + 3)
      END DO
      DEALLOCATE (use_atom, use_atom_type, dpmd_coord, dpmd_force, &
                  dpmd_virial, dpmd_atomic_energy, dpmd_atomic_virial, &
                  dpmd_cell, dpmd_atype)
#endif
      CALL timestop(handle)
   END SUBROUTINE deepmd_energy_store_force_virial

! **************************************************************************************************
!> \brief ...
!> \param fist_nonbond_env ...
!> \param force ...
!> \param pv_nonbond ...
!> \param use_virial ...
! **************************************************************************************************
   SUBROUTINE deepmd_add_force_virial(fist_nonbond_env, force, pv_nonbond, use_virial)
      TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
      REAL(KIND=dp)                            :: force(:, :), pv_nonbond(3, 3)
      LOGICAL, OPTIONAL                        :: use_virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deepmd_add_force_virial'
      INTEGER                                  :: handle

#ifdef __DEEPMD
      INTEGER                                  :: iat, iat_use
      TYPE(deepmd_data_type), POINTER          :: deepmd_data
#endif

      CALL timeset(routineN, handle)

#ifndef __DEEPMD
      CPABORT("In order to use DeePMD-kit you need to download "// &
              "and install libdeepmd_c (check CP2K manual for details)")
      MARK_USED(fist_nonbond_env)
      MARK_USED(force)
      MARK_USED(pv_nonbond)
#else
      CALL fist_nonbond_env_get(fist_nonbond_env, deepmd_data=deepmd_data)

      IF (.NOT. ASSOCIATED(deepmd_data)) RETURN

      DO iat_use = 1, SIZE(deepmd_data%use_indices)
         iat = deepmd_data%use_indices(iat_use)
         CPASSERT(iat >= 1 .AND. iat <= SIZE(force, 2))
         force(1:3, iat) = force(1:3, iat) + deepmd_data%force(1:3, iat_use)
      END DO

      IF (use_virial) THEN
         pv_nonbond = pv_nonbond + deepmd_data%virial
      END IF

#endif
      CALL timestop(handle)
   END SUBROUTINE deepmd_add_force_virial

END MODULE manybody_deepmd
