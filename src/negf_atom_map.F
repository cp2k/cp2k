!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Map atoms between various force environments.
!> \author Sergey Chulkov
! **************************************************************************************************

MODULE negf_atom_map
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              real_to_scaled
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_atom_map'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   PUBLIC :: negf_atom_map_type, negf_map_atomic_indices

! **************************************************************************************************
!> \brief Structure that maps the given atom in the sourse FORCE_EVAL section with another atom
!>        from the target FORCE_EVAL section.
! **************************************************************************************************
   TYPE negf_atom_map_type
      !> atomic index within the target FORCE_EVAL
      INTEGER                                            :: iatom
      !> cell replica
      INTEGER, DIMENSION(3)                              :: cell
   END TYPE negf_atom_map_type

   PRIVATE :: qs_kind_group_type, qs_kind_groups_create, qs_kind_groups_release

! **************************************************************************************************
!> \brief List of equivalent atoms.
! **************************************************************************************************
   TYPE qs_kind_group_type
      !> atomic symbol
      CHARACTER(len=2)                                   :: element_symbol
      !> number of atoms of this kind in 'atom_list'
      INTEGER                                            :: natoms
      !> number of spherical Gaussian functions per atom
      INTEGER                                            :: nsgf
      !> list of atomic indices
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_list
      !> atomic coordinates [3 x natoms]
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: r
   END TYPE qs_kind_group_type

CONTAINS
! **************************************************************************************************
!> \brief Map atoms in the cell 'subsys_device' listed in 'atom_list' with the corresponding
!>        atoms in the cell 'subsys_contact'.
!> \param atom_map        list of atoms in the cell 'subsys_contact' (initialised on exit)
!> \param atom_list       atomic indices of selected atoms in the cell 'subsys_device'
!> \param subsys_device   QuickStep subsystem of the device force environment
!> \param subsys_contact  QuickStep subsystem of the contact force environment
!> \param eps_geometry    accuracy in mapping atoms based on their Cartesian coordinates
!> \par History
!>   * 08.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_map_atomic_indices(atom_map, atom_list, subsys_device, subsys_contact, eps_geometry)
      TYPE(negf_atom_map_type), DIMENSION(:), &
         INTENT(out)                                     :: atom_map
      INTEGER, DIMENSION(:), INTENT(in)                  :: atom_list
      TYPE(qs_subsys_type), POINTER                      :: subsys_device, subsys_contact
      REAL(kind=dp), INTENT(in)                          :: eps_geometry

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_map_atomic_indices'

      CHARACTER(len=2)                                   :: element_device
      CHARACTER(len=default_string_length)               :: atom_str
      INTEGER                                            :: atom_index_device, handle, iatom, &
                                                            iatom_kind, ikind, ikind_contact, &
                                                            natoms, nkinds_contact, nsgf_device
      REAL(kind=dp), DIMENSION(3)                        :: coords, coords_error, coords_scaled
      TYPE(cell_type), POINTER                           :: cell_contact
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set_contact, particle_set_device
      TYPE(qs_kind_group_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: kind_groups_contact
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set_contact, qs_kind_set_device

      CALL timeset(routineN, handle)

      natoms = SIZE(atom_map, 1)
      CPASSERT(SIZE(atom_list) == natoms)

      CALL qs_subsys_get(subsys_device, particle_set=particle_set_device, qs_kind_set=qs_kind_set_device)
      CALL qs_subsys_get(subsys_contact, cell=cell_contact, particle_set=particle_set_contact, qs_kind_set=qs_kind_set_contact)

      CALL qs_kind_groups_create(kind_groups_contact, particle_set_contact, qs_kind_set_contact)
      nkinds_contact = SIZE(kind_groups_contact)

      DO iatom = 1, natoms
         atom_index_device = atom_list(iatom)
         CALL get_atomic_kind(particle_set_device(atom_index_device)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set_device(ikind), element_symbol=element_device, nsgf=nsgf_device)

         atom_map(iatom)%iatom = 0

         iterate_kind: DO ikind_contact = 1, nkinds_contact
            ! looking for an equivalent atomic kind (based on the element symbol and the number of atomic basis functions)
            IF (kind_groups_contact(ikind_contact)%element_symbol == element_device .AND. &
                kind_groups_contact(ikind_contact)%nsgf == nsgf_device) THEN

               ! loop over matching atoms
               DO iatom_kind = 1, kind_groups_contact(ikind_contact)%natoms
                  coords(1:3) = particle_set_device(atom_index_device)%r(1:3) - &
                                kind_groups_contact(ikind_contact)%r(1:3, iatom_kind)

                  CALL real_to_scaled(coords_scaled, coords, cell_contact)
                  coords_error = coords_scaled - REAL(NINT(coords_scaled), kind=dp)

                  IF (SQRT(DOT_PRODUCT(coords_error, coords_error)) < eps_geometry) THEN
                     atom_map(iatom)%iatom = kind_groups_contact(ikind_contact)%atom_list(iatom_kind)
                     atom_map(iatom)%cell = NINT(coords_scaled)
                     EXIT iterate_kind
                  END IF
               END DO
            END IF
         END DO iterate_kind

         IF (atom_map(iatom)%iatom == 0) THEN
            ! atom has not been found in the corresponding force_env
            WRITE (atom_str, '(A2,3(1X,F11.6))') element_device, particle_set_device(atom_index_device)%r

            CALL cp_abort(__LOCATION__, &
                          "Unable to map the atom ("//TRIM(atom_str)//") onto the atom from the corresponding FORCE_EVAL section")
         END IF
      END DO

      CALL qs_kind_groups_release(kind_groups_contact)

      CALL timestop(handle)
   END SUBROUTINE negf_map_atomic_indices

! **************************************************************************************************
!> \brief Group particles from 'particle_set' according to their atomic (QS) kind.
!> \param kind_groups   kind groups that will be created
!> \param particle_set  list of particles
!> \param qs_kind_set   list of QS kinds
!> \par History
!>   * 08.2017 created [Sergey Chulkov]
!> \note used within the subroutine negf_map_atomic_indices() in order to map atoms in
!>       a linear scalling fashion
! **************************************************************************************************
   SUBROUTINE qs_kind_groups_create(kind_groups, particle_set, qs_kind_set)
      TYPE(qs_kind_group_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(inout)                     :: kind_groups
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_kind_groups_create'

      INTEGER                                            :: handle, iatom, ikind, natoms, nkinds

      CALL timeset(routineN, handle)

      natoms = SIZE(particle_set)
      nkinds = 0

      DO iatom = 1, natoms
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         IF (nkinds < ikind) nkinds = ikind
      END DO

      ALLOCATE (kind_groups(nkinds))

      DO ikind = 1, nkinds
         kind_groups(ikind)%natoms = 0
         CALL get_qs_kind(qs_kind_set(ikind), element_symbol=kind_groups(ikind)%element_symbol, nsgf=kind_groups(ikind)%nsgf)
      END DO

      DO iatom = 1, natoms
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         kind_groups(ikind)%natoms = kind_groups(ikind)%natoms + 1
      END DO

      DO ikind = 1, nkinds
         ALLOCATE (kind_groups(ikind)%atom_list(kind_groups(ikind)%natoms))
         ALLOCATE (kind_groups(ikind)%r(3, kind_groups(ikind)%natoms))

         kind_groups(ikind)%natoms = 0
      END DO

      DO iatom = 1, natoms
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         kind_groups(ikind)%natoms = kind_groups(ikind)%natoms + 1

         kind_groups(ikind)%atom_list(kind_groups(ikind)%natoms) = iatom
         kind_groups(ikind)%r(1:3, kind_groups(ikind)%natoms) = particle_set(iatom)%r(1:3)
      END DO

      CALL timestop(handle)
   END SUBROUTINE qs_kind_groups_create

! **************************************************************************************************
!> \brief Release groups of particles.
!> \param kind_groups  kind groups to release
!> \par History
!>   * 08.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE qs_kind_groups_release(kind_groups)
      TYPE(qs_kind_group_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(inout)                     :: kind_groups

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_kind_groups_release'

      INTEGER                                            :: handle, ikind

      CALL timeset(routineN, handle)

      IF (ALLOCATED(kind_groups)) THEN
         DO ikind = SIZE(kind_groups), 1, -1
            IF (ALLOCATED(kind_groups(ikind)%atom_list)) DEALLOCATE (kind_groups(ikind)%atom_list)
            IF (ALLOCATED(kind_groups(ikind)%r)) DEALLOCATE (kind_groups(ikind)%r)
         END DO

         DEALLOCATE (kind_groups)
      END IF

      CALL timestop(handle)
   END SUBROUTINE qs_kind_groups_release
END MODULE negf_atom_map
