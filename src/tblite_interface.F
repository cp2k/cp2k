!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief interface to tblite [Johann Pototschnig]
!> \History creation 09.2024
! **************************************************************************************************
MODULE tblite_interface

#if defined(__TBLITE)
   USE mctc_io,                         ONLY: structure_type, new
#endif

   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE cell_types,                      ONLY: cell_type,  get_cell
   USE particle_types,                  ONLY: particle_type
   USE kinds,                           ONLY: dp
   USE atomic_kind_types,               ONLY: get_atomic_kind

#include "./base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_interface'

   PUBLIC :: init_geometry


#if defined(__TBLITE)

    TYPE(structure_type)                               :: mol

#endif

CONTAINS

#if defined(__TBLITE)

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param energy ...
!> \param calculate_forces ...
!> \param atevdw ...
! **************************************************************************************************
   SUBROUTINE init_geometry(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env


      CHARACTER(LEN=*), PARAMETER :: routineN = 'tblite_init_geometry'


      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      INTEGER                                            :: iatom, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: el_num
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xyz
      INTEGER                                            :: handle, ikind
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL, DIMENSION(3)                              :: lperiod

      CALL timeset(routineN, handle)

      !get info from environment vaiarable
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell)

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


      DEALLOCATE (el_num, xyz)

      CALL timestop(handle)

   END SUBROUTINE init_geometry




#else


   SUBROUTINE init_geometry(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      MARK_USED(qs_env)

      CPABORT("Build without TBLITE")

   END SUBROUTINE init_geometry

#endif


END MODULE tblite_interface

