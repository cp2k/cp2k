!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      HAF (16-Apr-2025) : Import into CP2K
!> \author HAF and yury-lysogorskiy and ralf-drautz
! **************************************************************************************************

MODULE manybody_ace

   USE ISO_C_BINDING,                   ONLY: C_ASSOCIATED
   USE ace_nlist,                       ONLY: ace_interface
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE bibliography,                    ONLY: Bochkarev2024,&
                                              Drautz2019,&
                                              Lysogorskiy2021,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE fist_nonbond_env_types,          ONLY: ace_data_type,&
                                              fist_nonbond_env_get,&
                                              fist_nonbond_env_set,&
                                              fist_nonbond_env_type
   USE kinds,                           ONLY: dp
   USE pair_potential_types,            ONLY: ace_type,&
                                              pair_potential_pp_type,&
                                              pair_potential_single_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC ace_energy_store_force_virial, ace_add_force_virial

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_ace'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param potparm ...
!> \param ace_data ...
! **************************************************************************************************
   SUBROUTINE init_ace_data(particle_set, atomic_kind_set, potparm, &
                            ace_data)

      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(ace_data_type), POINTER                       :: ace_data

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_ace_data'

      CHARACTER(2)                                       :: element_symbol
      INTEGER                                            :: ace_natom, handle, i, iat, iat_use, &
                                                            ikind, jkind, lkind, n_atoms
      INTEGER, ALLOCATABLE                               :: use_atom_type(:)
      INTEGER, DIMENSION(:), POINTER                     :: ak_alist
      LOGICAL, ALLOCATABLE                               :: use_atom(:)
      TYPE(pair_potential_single_type), POINTER          :: pot

      CALL timeset(routineN, handle)

      ! init ace_data
      IF (.NOT. ASSOCIATED(ace_data)) THEN
         ALLOCATE (ace_data)
      END IF

      n_atoms = SIZE(particle_set)
      ALLOCATE (use_atom(n_atoms))
      ALLOCATE (use_atom_type(n_atoms))
      use_atom = .FALSE.
      use_atom_type = 0

      DO ikind = 1, SIZE(atomic_kind_set)
         pot => potparm%pot(ikind, ikind)%pot
         DO i = 1, SIZE(pot%type)
            IF (pot%type(i) /= ace_type) CYCLE
            CALL get_atomic_kind(atomic_kind=atomic_kind_set(ikind), &
                                 element_symbol=element_symbol, &
                                 natom=lkind, atom_list=ak_alist)
            IF (lkind < 1) CYCLE
            ace_data%model = pot%set(i)%ace%model
            jkind = 0
            DO iat = 1, SIZE(ace_data%model%symbolc)
               IF (element_symbol == ace_data%model%symbolc(iat)) THEN
                  jkind = iat
                  EXIT
               END IF
            END DO
            CPASSERT(jkind > 0)
            DO iat = 1, lkind
               use_atom_type(ak_alist(iat)) = jkind
               use_atom(ak_alist(iat)) = .TRUE.
            END DO
         END DO ! i
      END DO ! ikind
      CPASSERT(C_ASSOCIATED(ace_data%model%c_ptr))

      ace_natom = COUNT(use_atom)

      IF (.NOT. ALLOCATED(ace_data%uctype)) THEN
         ALLOCATE (ace_data%uctype(1:ace_natom))
      END IF

      iat_use = 0
      DO iat = 1, n_atoms
         IF (.NOT. use_atom(iat)) CYCLE
         iat_use = iat_use + 1
         ace_data%uctype(iat_use) = use_atom_type(iat)
      END DO

      IF (iat_use > 0) THEN
         CALL cite_reference(Drautz2019)
         CALL cite_reference(Lysogorskiy2021)
         CALL cite_reference(Bochkarev2024)
      END IF

      IF (.NOT. ALLOCATED(ace_data%force)) THEN
         ALLOCATE (ace_data%force(3, ace_natom))
         ALLOCATE (ace_data%use_indices(ace_natom))
         ALLOCATE (ace_data%inverse_index_map(n_atoms))
      END IF
      CPASSERT(SIZE(ace_data%force, 2) == ace_natom)

      iat_use = 0
      ace_data%inverse_index_map(:) = 0
      DO iat = 1, n_atoms
         IF (use_atom(iat)) THEN
            iat_use = iat_use + 1
            ace_data%use_indices(iat_use) = iat
            ace_data%inverse_index_map(iat) = iat_use
         END IF
      END DO
      ace_data%natom = ace_natom
      DEALLOCATE (use_atom, use_atom_type)

      CALL timestop(handle)

   END SUBROUTINE init_ace_data

! **************************************************************************************************
!> \brief ...
!>     > \param particle_set ...
!> \param particle_set ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param potparm ...
!> \param fist_nonbond_env ...
!> \param pot_ace ...
! **************************************************************************************************
   SUBROUTINE ace_energy_store_force_virial(particle_set, cell, atomic_kind_set, potparm, &
                                            fist_nonbond_env, pot_ace)

      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      REAL(kind=dp), INTENT(OUT)                         :: pot_ace

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ace_energy_store_force_virial'

      INTEGER                                            :: handle
      REAL(kind=dp)                                      :: ace_virial(1:6)
      TYPE(ace_data_type), POINTER                       :: ace_data

      CALL timeset(routineN, handle)

      ! get ace_data to save force, virial info
      CALL fist_nonbond_env_get(fist_nonbond_env, ace_data=ace_data)
      IF (.NOT. ASSOCIATED(ace_data)) THEN
         ALLOCATE (ace_data)
         !initialize ace_data:
         CALL init_ace_data(particle_set, atomic_kind_set, potparm, ace_data)
         CALL fist_nonbond_env_set(fist_nonbond_env, ace_data=ace_data)
      END IF

      CALL ace_interface(ace_data%natom, ace_data%uctype, &
                         pot_ace, ace_data%force, ace_virial, &
                         fist_nonbond_env, cell, ace_data)

      ! convert units
      pot_ace = pot_ace/evolt
      ace_data%force = ace_data%force/(evolt/angstrom)
      ace_virial = ace_virial/evolt

      ! minus sign due to CP2K conventions
      ace_data%virial(1, 1) = -ace_virial(1)
      ace_data%virial(2, 2) = -ace_virial(2)
      ace_data%virial(3, 3) = -ace_virial(3)
      ace_data%virial(1, 2) = -ace_virial(4)
      ace_data%virial(2, 1) = -ace_virial(4)
      ace_data%virial(1, 3) = -ace_virial(5)
      ace_data%virial(3, 1) = -ace_virial(5)
      ace_data%virial(2, 3) = -ace_virial(6)
      ace_data%virial(3, 2) = -ace_virial(6)

      CALL timestop(handle)
   END SUBROUTINE ace_energy_store_force_virial

! **************************************************************************************************
!> \brief ...
!> \param fist_nonbond_env ...
!> \param force ...
!> \param pv_nonbond ...
!> \param use_virial ...
! **************************************************************************************************
   SUBROUTINE ace_add_force_virial(fist_nonbond_env, force, pv_nonbond, use_virial)
      TYPE(fist_nonbond_env_type), POINTER               :: fist_nonbond_env
      REAL(KIND=dp), INTENT(INOUT)                       :: force(:, :), pv_nonbond(3, 3)
      LOGICAL, OPTIONAL                                  :: use_virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ace_add_force_virial'

      INTEGER                                            :: handle, iat, iat_use
      TYPE(ace_data_type), POINTER                       :: ace_data

      CALL timeset(routineN, handle)

      CALL fist_nonbond_env_get(fist_nonbond_env, ace_data=ace_data)

      IF (.NOT. ASSOCIATED(ace_data)) RETURN

      DO iat_use = 1, SIZE(ace_data%use_indices)
         iat = ace_data%use_indices(iat_use)
         CPASSERT(iat >= 1 .AND. iat <= SIZE(force, 2))
         force(1:3, iat) = force(1:3, iat) + ace_data%force(1:3, iat_use)
      END DO

      IF (use_virial) THEN
         pv_nonbond = pv_nonbond + ace_data%virial
      END IF

      CALL timestop(handle)
   END SUBROUTINE ace_add_force_virial

END MODULE manybody_ace
