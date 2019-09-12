!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Set of routines handling the localization for molecular properties
! **************************************************************************************************
MODULE qs_loc_molecules
   USE cell_types,                      ONLY: pbc
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_should_output
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_max,&
                                              mp_minloc
   USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                              molecule_kind_type
   USE molecule_types,                  ONLY: molecule_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_loc_types,                    ONLY: qs_loc_env_new_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public ***
   PUBLIC :: wfc_to_molecule

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_loc_molecules'

CONTAINS

! **************************************************************************************************
!> \brief maps wfc's to molecules and also prints molecular dipoles
!> \param qs_env the qs_env in which the qs_env lives
!> \param qs_loc_env ...
!> \param loc_print_key ...
!> \param center ...
!> \param molecule_set ...
!> \param ispin ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE wfc_to_molecule(qs_env, qs_loc_env, loc_print_key, center, &
                              molecule_set, ispin, nspins)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_loc_env_new_type), INTENT(IN)              :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: loc_print_key
      REAL(KIND=dp), INTENT(IN)                          :: center(:, :)
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      INTEGER, INTENT(IN)                                :: ispin, nspins

      CHARACTER(len=*), PARAMETER :: routineN = 'wfc_to_molecule', &
         routineP = moduleN//':'//routineN

      INTEGER :: counter, first_atom, i, iatom, ikind, imol, imol_now, istate, k, local_location, &
         natom, natom_loc, natom_max, nkind, nmol, nstate, reference
      INTEGER, POINTER                                   :: wfc_to_atom_map(:)
      LOGICAL                                            :: do_berry
      REAL(KIND=dp)                                      :: dr(3), mydist(2), ria(3)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      REAL(KIND=dp), POINTER                             :: distance(:), r(:, :)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      logger => cp_get_default_logger()

      ! Molecular Dipoles availables only for nspin == 1
      IF (BTEST(cp_print_key_should_output(logger%iter_info, loc_print_key, &
                                           "MOLECULAR_DIPOLES"), cp_p_file)) THEN
         IF (nspins > 1) THEN
            CALL cp_abort(__LOCATION__, &
                          "Molecular Dipoles not implemented for SPIN multiplicity "// &
                          "larger than 1!")
         END IF

         ! Setup reference point and some warning..
         reference = section_get_ival(loc_print_key, keyword_name="MOLECULAR_DIPOLES%REFERENCE")
         CALL section_vals_val_get(loc_print_key, "MOLECULAR_DIPOLES%REF_POINT", r_vals=ref_point)
         CALL section_vals_val_get(loc_print_key, "MOLECULAR_DIPOLES%PERIODIC", l_val=do_berry)
      END IF
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      particle_set => qs_loc_env%particle_set
      para_env => qs_loc_env%para_env
      local_molecules => qs_loc_env%local_molecules
      nstate = SIZE(center, 2)
      ALLOCATE (wfc_to_atom_map(nstate))
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      nkind = SIZE(local_molecules%n_el)
      natom = 0
      natom_max = 0
      DO ikind = 1, nkind
         nmol = SIZE(local_molecules%list(ikind)%array)
         DO imol = 1, nmol
            i = local_molecules%list(ikind)%array(imol)
            molecule_kind => molecule_set(i)%molecule_kind
            CALL get_molecule_kind(molecule_kind=molecule_kind, natom=natom)
            natom_max = natom_max+natom
            IF (.NOT. ASSOCIATED(molecule_set(i)%lmi)) THEN
               ALLOCATE (molecule_set(i)%lmi(nspins))
               DO k = 1, nspins
                  NULLIFY (molecule_set(i)%lmi(k)%states)
               END DO
            ENDIF
            molecule_set(i)%lmi(ispin)%nstates = 0
            IF (ASSOCIATED(molecule_set(i)%lmi(ispin)%states)) THEN
               DEALLOCATE (molecule_set(i)%lmi(ispin)%states)
            END IF
         END DO
      END DO
      natom_loc = natom_max
      natom = natom_max

      CALL mp_max(natom_max, para_env%group)

      ALLOCATE (r(3, natom_max))

      ALLOCATE (distance(natom_max))

      !Zero all the stuff
      r(:, :) = 0.0_dp
      distance(:) = 1.E10_dp

      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      counter = 0
      nkind = SIZE(local_molecules%n_el)
      DO ikind = 1, nkind
         nmol = SIZE(local_molecules%list(ikind)%array)
         DO imol = 1, nmol
            i = local_molecules%list(ikind)%array(imol)
            molecule_kind => molecule_set(i)%molecule_kind
            first_atom = molecule_set(i)%first_atom
            CALL get_molecule_kind(molecule_kind=molecule_kind, natom=natom)

            DO iatom = 1, natom
               counter = counter+1
               r(:, counter) = particle_set(first_atom+iatom-1)%r(:)
            END DO
         END DO
      END DO

      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      DO istate = 1, nstate
         distance(:) = 1.E10_dp
         DO iatom = 1, natom_loc
            dr(1) = r(1, iatom)-center(1, istate)
            dr(2) = r(2, iatom)-center(2, istate)
            dr(3) = r(3, iatom)-center(3, istate)
            ria = pbc(dr, qs_loc_env%cell)
            distance(iatom) = SQRT(DOT_PRODUCT(ria, ria))
         END DO

         !combine distance() from all procs
         local_location = MAX(1, MINLOC(distance, DIM=1))

         mydist(1) = distance(local_location)
         mydist(2) = para_env%mepos

         CALL mp_minloc(mydist, para_env%group)

         IF (mydist(2) == para_env%mepos) THEN
            wfc_to_atom_map(istate) = local_location
         ELSE
            wfc_to_atom_map(istate) = 0
         END IF
      END DO
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      IF (natom_loc /= 0) THEN
         DO istate = 1, nstate
            iatom = wfc_to_atom_map(istate)
            IF (iatom /= 0) THEN
               counter = 0
               nkind = SIZE(local_molecules%n_el)
               DO ikind = 1, nkind
                  nmol = SIZE(local_molecules%list(ikind)%array)
                  DO imol = 1, nmol
                     imol_now = local_molecules%list(ikind)%array(imol)
                     molecule_kind => molecule_set(imol_now)%molecule_kind
                     CALL get_molecule_kind(molecule_kind=molecule_kind, natom=natom)
                     counter = counter+natom
                     IF (counter >= iatom) EXIT
                  END DO
                  IF (counter >= iatom) EXIT
               END DO
               i = molecule_set(imol_now)%lmi(ispin)%nstates
               i = i+1
               molecule_set(imol_now)%lmi(ispin)%nstates = i
               CALL reallocate(molecule_set(imol_now)%lmi(ispin)%states, 1, i)
               molecule_set(imol_now)%lmi(ispin)%states(i) = istate
            END IF
         END DO
      END IF

      DEALLOCATE (distance)
      DEALLOCATE (r)
      DEALLOCATE (wfc_to_atom_map)

   END SUBROUTINE wfc_to_molecule
   !------------------------------------------------------------------------------

END MODULE qs_loc_molecules

