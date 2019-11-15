!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines that work on qs_subsys_type
!> \author Ole Schuett
! **************************************************************************************************
MODULE qs_subsys_methods
   USE atom_types,                      ONLY: lmat
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_methods,                    ONLY: read_cell,&
                                              write_cell
   USE cell_types,                      ONLY: cell_clone,&
                                              cell_create,&
                                              cell_release,&
                                              cell_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_subsys_methods,               ONLY: cp_subsys_create
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_release,&
                                              cp_subsys_set,&
                                              cp_subsys_type
   USE external_potential_types,        ONLY: all_potential_type,&
                                              get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                              molecule_kind_type,&
                                              set_molecule_kind
   USE qs_kind_types,                   ONLY: create_qs_kind_set,&
                                              get_qs_kind,&
                                              init_atom_electronic_state,&
                                              qs_kind_type
   USE qs_subsys_types,                 ONLY: qs_subsys_set,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_subsys_methods'

   PUBLIC :: qs_subsys_create

CONTAINS

! **************************************************************************************************
!> \brief Creates a qs_subsys. Optionally an existsing cp_subsys is used.
!> \param subsys ...
!> \param para_env ...
!> \param root_section ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \param cp_subsys ...
!> \param cell ...
!> \param cell_ref ...
!> \param elkind ...
! **************************************************************************************************
   SUBROUTINE qs_subsys_create(subsys, para_env, root_section, force_env_section, subsys_section, &
                               use_motion_section, cp_subsys, cell, cell_ref, elkind)
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), OPTIONAL, POINTER         :: root_section
      TYPE(section_vals_type), POINTER                   :: force_env_section, subsys_section
      LOGICAL, INTENT(IN)                                :: use_motion_section
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: cp_subsys
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell, cell_ref
      LOGICAL, INTENT(IN), OPTIONAL                      :: elkind

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_subsys_create', &
         routineP = moduleN//':'//routineN

      LOGICAL                                            :: use_ref_cell
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: my_cell, my_cell_ref
      TYPE(cp_subsys_type), POINTER                      :: my_cp_subsys
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: cell_section, kind_section

      NULLIFY (atomic_kind_set, qs_kind_set, cell_section, kind_section, my_cell, my_cell_ref, my_cp_subsys)

      IF (ASSOCIATED(subsys)) CPABORT("qs_subsys_create: subsys already associated")

      ! create cp_subsys
      IF (PRESENT(cp_subsys)) THEN
         my_cp_subsys => cp_subsys
      ELSE IF (PRESENT(root_section)) THEN
         CALL cp_subsys_create(my_cp_subsys, para_env, root_section=root_section, &
                               force_env_section=force_env_section, &
                               subsys_section=subsys_section, &
                               use_motion_section=use_motion_section, &
                               elkind=elkind)
      ELSE
         CPABORT("qs_subsys_create: cp_subsys or root_section needed")
      END IF

      ! create cp_subsys%cell
      !TODO: moved to cp_subsys_create(), needs further disentanglement of cell_ref.
      IF (PRESENT(cell)) THEN
         my_cell => cell
         IF (PRESENT(cell_ref)) THEN
            my_cell_ref => cell_ref
            use_ref_cell = .TRUE.
         ELSE
            CALL cell_create(my_cell_ref)
            CALL cell_clone(my_cell, my_cell_ref)
            use_ref_cell = .FALSE.
         END IF
      ELSE
         cell_section => section_vals_get_subs_vals(subsys_section, "CELL")
         CALL read_cell(my_cell, my_cell_ref, use_ref_cell=use_ref_cell, &
                        cell_section=cell_section, para_env=para_env)
      END IF
      CALL cp_subsys_set(my_cp_subsys, cell=my_cell)
      CALL write_cell(my_cell, subsys_section, cell_ref=my_cell_ref)

      ! setup qs_kinds
      CALL cp_subsys_get(my_cp_subsys, atomic_kind_set=atomic_kind_set)
      kind_section => section_vals_get_subs_vals(subsys_section, "KIND")
      CALL create_qs_kind_set(qs_kind_set, atomic_kind_set, kind_section, &
                              para_env, force_env_section)

      CALL num_ao_el_per_molecule(my_cp_subsys%molecule_kinds%els, &
                                  qs_kind_set)

      ! finally create qs_subsys
      ALLOCATE (subsys)
      CALL qs_subsys_set(subsys, &
                         cp_subsys=my_cp_subsys, &
                         cell_ref=my_cell_ref, &
                         use_ref_cell=use_ref_cell, &
                         qs_kind_set=qs_kind_set)

      IF (.NOT. PRESENT(cell)) CALL cell_release(my_cell)
      IF (.NOT. PRESENT(cell_ref)) CALL cell_release(my_cell_ref)
      IF (.NOT. PRESENT(cp_subsys)) CALL cp_subsys_release(my_cp_subsys)
   END SUBROUTINE qs_subsys_create

! **************************************************************************************************
!> \brief   Read a molecule kind data set from the input file.
!> \param molecule_kind_set ...
!> \param qs_kind_set ...
!> \date    22.11.2004
!> \par History
!>      Rustam Z. Khaliullin 10.2014 - charges and electrons of molecules
!>                                     are now in agreement with atomic guess
!> \author  MI
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE num_ao_el_per_molecule(molecule_kind_set, qs_kind_set)

      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CHARACTER(len=*), PARAMETER :: routineN = 'num_ao_el_per_molecule', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: arbitrary_spin, iatom, ikind, imol, &
                                                            n_ao, natom, nmol_kind, nsgf, nspins, &
                                                            z_molecule
      INTEGER, DIMENSION(0:lmat, 10)                     :: ne_core, ne_elem, ne_explicit
      INTEGER, DIMENSION(2)                              :: n_occ_alpha_and_beta
      REAL(KIND=dp)                                      :: charge_molecule, zeff, zeff_correction
      REAL(KIND=dp), DIMENSION(0:lmat, 10, 2)            :: edelta
      TYPE(all_potential_type), POINTER                  :: all_potential
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      IF (ASSOCIATED(molecule_kind_set)) THEN

         nspins = 2
         nmol_kind = SIZE(molecule_kind_set, 1)
         natom = 0

         !   *** Initialize the molecule kind data structure ***
         ARBITRARY_SPIN = 1
         DO imol = 1, nmol_kind

            molecule_kind => molecule_kind_set(imol)
            CALL get_molecule_kind(molecule_kind=molecule_kind, &
                                   natom=natom)
            !nelectron = 0
            n_ao = 0
            n_occ_alpha_and_beta(1:nspins) = 0
            z_molecule = 0

            DO iatom = 1, natom

               atomic_kind => molecule_kind%atom_list(iatom)%atomic_kind
               CALL get_atomic_kind(atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), &
                                basis_set=orb_basis_set, &
                                all_potential=all_potential, &
                                gth_potential=gth_potential, &
                                sgp_potential=sgp_potential)

               ! Obtain the electronic state of the atom
               ! The same state is used to calculate the ATOMIC GUESS
               ! It is great that we are consistent with ATOMIC_GUESS
               CALL init_atom_electronic_state(atomic_kind=atomic_kind, &
                                               qs_kind=qs_kind_set(ikind), &
                                               ncalc=ne_explicit, &
                                               ncore=ne_core, &
                                               nelem=ne_elem, &
                                               edelta=edelta)

               ! If &BS section is used ATOMIC_GUESS is calculated twice
               ! for two separate wfns with their own alpha-beta combinations
               ! This is done to break the spin symmetry of the initial wfn
               ! For now, only alpha part of &BS is used to count electrons on
               ! molecules
               ! Get the number of explicit electrons (i.e. with orbitals)
               ! For now, only the total number of electrons can be obtained
               ! from init_atom_electronic_state
               n_occ_alpha_and_beta(ARBITRARY_SPIN) = &
                  n_occ_alpha_and_beta(ARBITRARY_SPIN) + SUM(ne_explicit) + &
                  SUM(NINT(2*edelta(:, :, ARBITRARY_SPIN)))
               ! We need a way to specify the number of alpha and beta electrons
               ! on each molecule (i.e. multiplicity is not enough)
               !n_occ(ispin) = n_occ(ispin) + SUM(ne_explicit) + SUM(NINT(2*edelta(:, :, ispin)))

               IF (ASSOCIATED(all_potential)) THEN
                  CALL get_potential(potential=all_potential, zeff=zeff, &
                                     zeff_correction=zeff_correction)
               ELSE IF (ASSOCIATED(gth_potential)) THEN
                  CALL get_potential(potential=gth_potential, zeff=zeff, &
                                     zeff_correction=zeff_correction)
               ELSE IF (ASSOCIATED(sgp_potential)) THEN
                  CALL get_potential(potential=sgp_potential, zeff=zeff, &
                                     zeff_correction=zeff_correction)
               ELSE
                  zeff = 0.0_dp
                  zeff_correction = 0.0_dp
               END IF
               z_molecule = z_molecule + NINT(zeff - zeff_correction)

               ! this one does not work because nelem is not adjusted in the symmetry breaking code
               !CALL get_atomic_kind(atomic_kind,z=z)
               !z_molecule=z_molecule+z

               IF (ASSOCIATED(orb_basis_set)) THEN
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set, nsgf=nsgf)
               ELSE
                  nsgf = 0
               ENDIF
               n_ao = n_ao + nsgf

            END DO ! iatom

            ! At this point we have the number of electrons (alpha+beta) on the molecule
            !  as they are seen by the ATOMIC GUESS routines
            charge_molecule = REAL(z_molecule - n_occ_alpha_and_beta(ARBITRARY_SPIN), dp)
            CALL set_molecule_kind(molecule_kind=molecule_kind, &
                                   nelectron=n_occ_alpha_and_beta(ARBITRARY_SPIN), &
                                   charge=charge_molecule, &
                                   nsgf=n_ao)

         END DO ! imol
      END IF

   END SUBROUTINE num_ao_el_per_molecule

END MODULE qs_subsys_methods
