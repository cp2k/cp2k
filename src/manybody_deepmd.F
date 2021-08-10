!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE manybody_deepmd

   USE cp_log_handling, ONLY: cp_logger_get_default_io_unit
   USE atomic_kind_types, ONLY: atomic_kind_type
   USE bibliography, ONLY: Deepmd2018, &
                           cite_reference
   USE cell_types, ONLY: cell_type
   USE cp_para_types, ONLY: cp_para_env_type
   USE fist_nonbond_env_types, ONLY: fist_nonbond_env_get, &
                                     fist_nonbond_env_set, &
                                     fist_nonbond_env_type, &
                                     deepmd_data_type
   USE kinds, ONLY: dp
   USE pair_potential_types, ONLY: pair_potential_pp_type, &
                                   pair_potential_single_type, &
                                   deepmd_pot_type, &
                                   deepmd_type
   USE particle_types, ONLY: particle_type
   USE physcon, ONLY: angstrom, &
                      evolt
#ifdef __DEEPMD
   USE deepmd_wrapper, ONLY: nnp, &
                             create_nnp, &
                             delete_nnp_c, &
                             compute_nnp
#endif

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   LOGICAL, SAVE, PRIVATE :: do_create(2) = .TRUE.
   INTEGER, SAVE, PRIVATE :: do_count = 0
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
      TYPE(particle_type), POINTER             :: particle_set(:)
      TYPE(cell_type), POINTER                 :: cell
      TYPE(atomic_kind_type), POINTER          :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER    :: potparm
      TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
      REAL(kind=dp), TARGET                            :: pot_deepmd
      TYPE(cp_para_env_type), OPTIONAL, &
         POINTER                                :: para_env

      CHARACTER(LEN=*), PARAMETER :: &
         routineN = 'deepmd_energy_store_force_virial', &
         routineP = moduleN//':'//routineN

#ifdef __DEEPMD
      INTEGER                                  :: i, iat, iat_use, ikind, &
                                                  jkind, n_atoms, n_atoms_use, &
                                                  output_unit
      LOGICAL                                  :: do_parallel
      LOGICAL, ALLOCATABLE                     :: use_atom(:)
      REAL(kind=dp)                            :: lattice(3, 3)
      TYPE(pair_potential_single_type), &
         POINTER                                :: pot
      TYPE(deepmd_data_type), POINTER            :: deepmd_data
      TYPE(deepmd_pot_type), POINTER             :: deepmd
      TYPE(nnp), SAVE                                :: deepmd_tmp_pot1, &
                                                        deepmd_tmp_pot2
      REAL(kind=dp), ALLOCATABLE, TARGET       :: dpmd_force(:), dpmd_virial(:), &
                                                  dpmd_atom_ener(:), dpmd_atom_virial(:), &
                                                  dpmd_box(:), dpmd_coord(:)
      INTEGER, ALLOCATABLE, TARGET            :: dpmd_atype(:), use_atom_type(:)
      INTEGER, TARGET                         :: vecsize
      REAL(kind=dp), POINTER                  :: dforce(:), dvirial(:), datom_ener(:), &
                                                 datom_virial(:), dbox(:), dcoord(:), dener
      INTEGER, POINTER                        :: datype(:)

#endif
#ifndef __DEEPMD
      MARK_USED(particle_set)
      MARK_USED(cell)
      MARK_USED(atomic_kind_set)
      MARK_USED(potparm)
      MARK_USED(fist_nonbond_env)
      MARK_USED(pot_deepmd)
      MARK_USED(para_env)
      CALL cp_abort(__LOCATION__, "In order to use DeePMD-kit you need to download "// &
                    "and install the DeePMD-kit (check CP2K manual for details)")
#else
      n_atoms = SIZE(particle_set)
      ALLOCATE (use_atom(n_atoms))
      ALLOCATE (use_atom_type(n_atoms))
      use_atom = .FALSE.
      use_atom_type = 100 ! just initial the array

      NULLIFY (deepmd)

      DO ikind = 1, SIZE(atomic_kind_set)
      DO jkind = 1, SIZE(atomic_kind_set)
         pot => potparm%pot(ikind, jkind)%pot
         DO i = 1, SIZE(pot%type)
            IF (pot%type(i) /= deepmd_type) CYCLE
            IF (.NOT. ASSOCIATED(deepmd)) deepmd => pot%set(i)%deepmd
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
!     print*, use_atom_type
      n_atoms_use = COUNT(use_atom)
      !      ALLOCATE (pos(3, n_atoms_use), force(3, n_atoms_use), elem_symbol(n_atoms_use))
      vecsize = n_atoms_use
      ALLOCATE (dpmd_box(9))
      ALLOCATE (dpmd_atype(vecsize))
      ALLOCATE (dpmd_coord(vecsize*3))
      ALLOCATE (dpmd_force(vecsize*3))
      ALLOCATE (dpmd_virial(9))
      ALLOCATE (dpmd_atom_ener(vecsize))
      ALLOCATE (dpmd_atom_virial(vecsize*9))

! initial the properties
      pot_deepmd = 0.0
      dpmd_force = 0.0
      dpmd_virial = 0.0
      dpmd_atom_ener = 0.0
      dpmd_atom_virial = 0.0
      iat_use = 0
      DO iat = 1, n_atoms
         IF (.NOT. use_atom(iat)) CYCLE
         iat_use = iat_use + 1
         dpmd_coord((iat_use - 1)*3 + 1:(iat_use - 1)*3 + 3) = particle_set(iat)%r*angstrom
         dpmd_atype(iat_use) = use_atom_type(iat)
      END DO
      IF (iat_use > 0) CALL cite_reference(Deepmd2018)
      output_unit = cp_logger_get_default_io_unit()
      lattice = cell%hmat*angstrom
! change matrix to one d array
      DO i = 1, 3
         dpmd_box((i - 1)*3 + 1:(i - 1)*3 + 3) = lattice(:, i)
      END DO
      do_parallel = .FALSE.
      IF (PRESENT(para_env)) THEN
         do_parallel = para_env%num_pe > 1
      END IF
! change to pointer type
      dener => pot_deepmd
      dforce => dpmd_force
      dvirial => dpmd_virial
      datom_ener => dpmd_atom_ener
      datom_virial => dpmd_atom_virial
      dcoord => dpmd_coord
      dbox => dpmd_box
      datype => dpmd_atype
!      IF (do_parallel) THEN
!         CALL quip_unified_wrapper( &
!            N=n_atoms_use, pos=pos, lattice=lattice, symbol=elem_symbol, &
!            quip_param_file=TRIM(quip%quip_file_name), &
!            quip_param_file_len=LEN_TRIM(quip%quip_file_name), &
!            init_args_str=TRIM(quip%init_args), &
!            init_args_str_len=LEN_TRIM(quip%init_args), &
!            calc_args_str=TRIM(quip%calc_args), &
!            calc_args_str_len=LEN_TRIM(quip%calc_args), &
!            energy=pot_quip, force=force, virial=virial, &
!            output_unit=output_unit, mpi_communicator=para_env%group)
!      ELSE
! initial the two different potential if needed
      IF (do_count == 0) THEN
         IF (do_create(1)) THEN
            deepmd_tmp_pot1 = create_nnp(TRIM(deepmd%deepmd_file_name))
            do_create(1) = .FALSE.
         END IF
         CALL compute_nnp( &
            pot=deepmd_tmp_pot1%ptr, &
            vecsize=vecsize, &
            dener=dener, &
            dforce=dforce, &
            dvirial=dvirial, &
            datom_ener=datom_ener, &
            datom_virial=datom_virial, &
            dcoord=dcoord, &
            datype=datype, &
            dbox=dbox)
         do_count = 1
      ELSE IF (do_count == 1) THEN
         IF (do_create(2)) THEN
            deepmd_tmp_pot2 = create_nnp(TRIM(deepmd%deepmd_file_name))
            do_create(2) = .FALSE.
         END IF
         CALL compute_nnp( &
            pot=deepmd_tmp_pot2%ptr, &
            vecsize=vecsize, &
            dener=dener, &
            dforce=dforce, &
            dvirial=dvirial, &
            datom_ener=datom_ener, &
            datom_virial=datom_virial, &
            dcoord=dcoord, &
            datype=datype, &
            dbox=dbox)
         do_count = 0
      END IF

!      ENDIF
      ! convert units
      pot_deepmd = pot_deepmd/evolt
      dpmd_force = dpmd_force/(evolt/angstrom)
      dpmd_virial = dpmd_virial/evolt
      ! account for double counting from multiple MPI processes
!      IF (PRESENT(para_env)) pot_quip = pot_quip/REAL(para_env%num_pe, dp)
!      IF (PRESENT(para_env)) force = force/REAL(para_env%num_pe, dp)
!      IF (PRESENT(para_env)) virial = virial/REAL(para_env%num_pe, dp)
      ! get quip_data to save force, virial info
      CALL fist_nonbond_env_get(fist_nonbond_env, deepmd_data=deepmd_data)
      IF (.NOT. ASSOCIATED(deepmd_data)) THEN
         ALLOCATE (deepmd_data)
         CALL fist_nonbond_env_set(fist_nonbond_env, deepmd_data=deepmd_data)
         NULLIFY (deepmd_data%use_indices, deepmd_data%force)
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
!      DEALLOCATE(dener, dforce, dvirial, datom_ener, datom_virial, dcoord, dbox, datype)
!      NULLIFY(dener)
!      NULLIFY(dforce)
!      NULLIFY(dvirial)
!      NULLIFY(datom_ener)
!      NULLIFY(datom_virial)
!      NULLIFY(dcoord)
!      NULLIFY(dbox)
!      NULLIFY(datype)
!       DEALLOCATE(use_atom, use_atom_type)
!      CALL delete_nnp_c(deepmd_tmp_pot%ptr)
      DEALLOCATE (use_atom, use_atom_type, dpmd_coord, dpmd_force, &
                  dpmd_virial, dpmd_atom_ener, dpmd_atom_virial, &
                  dpmd_box, dpmd_atype)
#endif
   END SUBROUTINE deepmd_energy_store_force_virial

! **************************************************************************************************
!> \brief ...
!> \param fist_nonbond_env ...
!> \param force ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE deepmd_add_force_virial(fist_nonbond_env, force, virial)
      TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
      REAL(KIND=dp)                            :: force(:, :), virial(3, 3)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deepmd_add_force_virial', &
                                     routineP = moduleN//':'//routineN

#ifdef __DEEPMD
      INTEGER                                  :: iat, iat_use
      TYPE(deepmd_data_type), POINTER            :: deepmd_data
#endif

#ifndef __DEEPMD
      MARK_USED(fist_nonbond_env)
      MARK_USED(force)
      MARK_USED(virial)
      RETURN
#else
      CALL fist_nonbond_env_get(fist_nonbond_env, deepmd_data=deepmd_data)
      IF (.NOT. ASSOCIATED(deepmd_data)) RETURN

      DO iat_use = 1, SIZE(deepmd_data%use_indices)
         iat = deepmd_data%use_indices(iat_use)
         CPASSERT(iat >= 1 .AND. iat <= SIZE(force, 2))
         force(1:3, iat) = force(1:3, iat) + deepmd_data%force(1:3, iat_use)
      END DO
      virial = virial + deepmd_data%virial
#endif
   END SUBROUTINE deepmd_add_force_virial

END MODULE manybody_deepmd
