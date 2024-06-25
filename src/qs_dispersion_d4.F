!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of dispersion using pair potentials
!> \author JGH
! **************************************************************************************************
MODULE qs_dispersion_d4

   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE dftd4,                           ONLY: d4_model,&
                                              damping_param,&
                                              get_dispersion,&
                                              get_rational_damping,&
                                              new,&
                                              new_d4_model,&
                                              realspace_cutoff,&
                                              structure_type
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_force_types,                  ONLY: qs_force_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_d4'

   PUBLIC :: calculate_dispersion_d4_pairpot

! **************************************************************************************************

CONTAINS

#if defined(__DFTD4)
   ! **************************************************************************************************
!> \brief ...
!> \param dispersion_env ...
!> \param particle_set ...
!> \param cell ...
!> \param energy ...
!> \param force ...
!> \param atom_of_kind ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE calculate_dispersion_d4_pairpot(dispersion_env, particle_set, cell, energy, force, atom_of_kind, virial)

      class(damping_param), ALLOCATABLE                  :: param
      TYPE(d4_model)                                     :: disp
      TYPE(structure_type)                               :: mol
      TYPE(realspace_cutoff)                             :: cutoff

      TYPE(qs_dispersion_type), POINTER, INTENT(IN)                       :: dispersion_env
      TYPE(particle_type), DIMENSION(:), POINTER, INTENT(IN)              :: particle_set
      TYPE(cell_type), POINTER, INTENT(IN)                                :: cell
      REAL(KIND=dp), INTENT(OUT)                                          :: energy
      TYPE(qs_force_type), DIMENSION(:), OPTIONAL, POINTER, INTENT(INOUT) :: force
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)                         :: atom_of_kind
      REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL, INTENT(OUT)               :: virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dispersion_d4_pairpot'

      INTEGER                                            :: iatom, natom, ind_atom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: el_num
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gradient, xyz
      REAL(KIND=dp), DIMENSION(3, 3)                     :: stress
      INTEGER                                            :: handle, ikind
      !REAL(KIND=dp), DIMENSION(3, 3)                    :: h
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL, DIMENSION(3)                              :: lperiod
      CALL timeset(routineN, handle)

      !get information about particles
      natom = SIZE(particle_set)
      ALLOCATE (xyz(3, natom))
      ALLOCATE (el_num(natom))
      DO iatom = 1, natom
         xyz(:, iatom) = particle_set(iatom)%r(:)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         el_num(iatom) = ikind
      END DO

      !get information about cell / lattice
      CALL get_cell(cell=cell, periodic=periodic)
      lperiod(1) = periodic(1) == 1
      lperiod(2) = periodic(2) == 1
      lperiod(3) = periodic(3) == 1

      !prepare for the call to the dispersion function
      CALL new(mol, el_num, xyz, lattice=cell%hmat, periodic=lperiod)
      CALL new_d4_model(disp, mol)
      CALL get_rational_damping(dispersion_env%ref_functional, param, s9=dispersion_env%s9)

      ! Coordination number cutoff
      cutoff%cn = dispersion_env%rc_cn
      ! Two-body interaction cutoff
      cutoff%disp2 = dispersion_env%rc_d4*2._dp
      ! Three-body interaction cutoff
      cutoff%disp3 = dispersion_env%rc_disp*2._dp

      !> Wrapper to handle the evaluation of dispersion energy and derivatives
      IF (PRESENT(force)) THEN
         IF (PRESENT(atom_of_kind)) THEN
            ALLOCATE (gradient(3, mol%nat))
            CALL get_dispersion(mol, disp, param, cutoff, energy, gradient, stress)
            IF (PRESENT(virial)) THEN
               virial = stress
            END IF
            DO iatom = 1, natom
               ikind = el_num(iatom)
               ind_atom = atom_of_kind(iatom)
               force(ikind)%dispersion(:, ind_atom) = force(ikind)%dispersion(:, ind_atom) + gradient(:, iatom)
            END DO
            DEALLOCATE (gradient)
         ELSE
            CPABORT("missing atom_of_kind")
         END IF
      ELSE
         IF (PRESENT(virial)) THEN
            CPABORT("missing force for virial term")
         END IF
         CALL get_dispersion(mol, disp, param, cutoff, energy)
      END IF

      DEALLOCATE (el_num, xyz)

      CALL timestop(handle)

   END SUBROUTINE calculate_dispersion_d4_pairpot

#else

   ! **************************************************************************************************
!> \brief ...
!> \param dispersion_env ...
!> \param particle_set ...
!> \param cell ...
!> \param energy ...
!> \param force ...
!> \param atom_of_kind ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE calculate_dispersion_d4_pairpot(dispersion_env, particle_set, cell, energy, force, atom_of_kind, virial)

      TYPE(qs_dispersion_type), INTENT(IN), POINTER      :: dispersion_env
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), INTENT(OUT)                         :: energy
      TYPE(qs_force_type), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL, POINTER                               :: force
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: atom_of_kind
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT), &
         OPTIONAL                                        :: virial

      CPABORT("Build without DFTD4")

   END SUBROUTINE calculate_dispersion_d4_pairpot

#endif

END MODULE qs_dispersion_d4
