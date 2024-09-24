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
   USE mctc_env, ONLY: error_type
   USE mctc_io, ONLY: structure_type, new
   USE tblite_xtb_calculator, ONLY: xtb_calculator, new_xtb_calculator
   USE tblite_xtb_gfn1, ONLY: new_gfn1_calculator
   USE tblite_xtb_gfn2, ONLY: new_gfn2_calculator
   USE tblite_xtb_ipea1, ONLY: new_ipea1_calculator
   USE tblite_container, ONLY: container_cache
#endif

   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE cell_types, ONLY: cell_type, get_cell
   USE particle_types, ONLY: particle_type
   USE kinds, ONLY: dp
   USE atomic_kind_types, ONLY: get_atomic_kind

#include "./base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER                  :: gfn1xtb = 1, gfn2xtb = 2, ipea1xtb = 3

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_interface'

   PUBLIC :: tb_init_geometry, tb_type

   TYPE tb_type
      REAL(KIND=dp), ALLOCATABLE                         :: e_hal(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_rep(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_disp(:)
      REAL(KIND=dp), ALLOCATABLE                         :: grad(:, :)
      REAL(KIND=dp), ALLOCATABLE                         :: sigma(:, :)

#if defined(__TBLITE)
      TYPE(structure_type)                               :: mol
      TYPE(xtb_calculator)                               :: calc
      TYPE(container_cache)                              :: dcache
#endif

   END TYPE tb_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_init_geometry(qs_env, tb)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tb_type), POINTER                             :: tb

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tblite_init_geometry'

      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      INTEGER                                            :: iatom, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: el_num
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xyz
      INTEGER                                            :: handle, ikind
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL, DIMENSION(3)                              :: lperiod

#if defined(__TBLITE)

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
         IF (ikind < 1 .OR. ikind > 85) THEN
            CPABORT("only elements 1-85 are supported by tblite")
         END IF
         el_num(iatom) = ikind
      END DO

      !get information about cell / lattice
      CALL get_cell(cell=cell, periodic=periodic)
      lperiod(1) = periodic(1) == 1
      lperiod(2) = periodic(2) == 1
      lperiod(3) = periodic(3) == 1

      !prepare for the call to the dispersion function
      CALL new(tb%mol, el_num, xyz, lattice=cell%hmat, periodic=lperiod)

      !allocate quantities later required
      ALLOCATE (tb%e_hal(natom), tb%e_rep(natom), tb%e_disp(natom))

      DEALLOCATE (el_num, xyz)

      CALL timestop(handle)

#else
      CPABORT("Build without TBLITE")
#endif

   END SUBROUTINE tb_init_geometry

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param typ ...
! **************************************************************************************************
   SUBROUTINE tb_set_calculator(tb, typ)

      TYPE(tb_type), POINTER                             :: tb
      INTEGER                                            :: typ

#if defined(__TBLITE)
      TYPE(error_type), ALLOCATABLE                      :: error

      SELECT CASE (typ)
      CASE default
         CPABORT("Unknown xtb type")
      CASE (gfn1xtb)
         CALL new_gfn1_calculator(tb%calc, tb%mol, error)
      CASE (gfn2xtb)
         CALL new_gfn2_calculator(tb%calc, tb%mol, error)
      CASE (ipea1xtb)
         CALL new_ipea1_calculator(tb%calc, tb%mol, error)
      END SELECT

#else
      CPABORT("Build without TBLITE")
#endif

   END SUBROUTINE tb_set_calculator

! **************************************************************************************************
!> \brief ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_init_ham(tb)

      TYPE(tb_type), POINTER                             :: tb

      INTEGER                                            :: icase
#if defined(__TBLITE)
      TYPE(container_cache)                              :: hcache, rcache

      tb%e_hal = 0.0_dp
      tb%e_rep = 0.0_dp
      tb%e_hal = 0.0_dp
      icase = 0
      IF (ALLOCATED(tb%grad)) THEN
         tb%grad = 0.0_dp
         IF (ALLOCATED(tb%sigma)) THEN
            tb%sigma = 0.0_dp
            icase = 2
         ELSE
            icase = 1
         END IF
      END IF

      IF (ALLOCATED(tb%calc%halogen)) THEN
         CALL tb%calc%halogen%update(tb%mol, hcache)
         SELECT CASE (icase)
         CASE (0)
            CALL tb%calc%halogen%get_engrad(tb%mol, hcache, tb%e_hal)
         CASE (1)
            CALL tb%calc%halogen%get_engrad(tb%mol, hcache, tb%e_hal, tb%grad)
         CASE (2)
            CALL tb%calc%halogen%get_engrad(tb%mol, hcache, tb%e_hal, &
            & tb%grad, tb%sigma)
         END SELECT
      END IF

      IF (ALLOCATED(tb%calc%repulsion)) THEN
         CALL tb%calc%repulsion%update(tb%mol, rcache)
         SELECT CASE (icase)
         CASE (0)
            CALL tb%calc%repulsion%get_engrad(tb%mol, rcache, tb%e_rep)
         CASE (1)
            CALL tb%calc%repulsion%get_engrad(tb%mol, rcache, tb%e_rep, tb%grad)
         CASE (2)
            CALL tb%calc%repulsion%get_engrad(tb%mol, rcache, tb%e_rep, &
            & tb%grad, tb%sigma)
         END SELECT
      END IF

      IF (ALLOCATED(tb%calc%dispersion)) THEN
         CALL tb%calc%dispersion%update(tb%mol, tb%dcache)
         SELECT CASE (icase)
         CASE (0)
            CALL tb%calc%dispersion%get_engrad(tb%mol, tb%dcache, tb%e_disp)
         CASE (1)
            CALL tb%calc%dispersion%get_engrad(tb%mol, tb%dcache, tb%e_disp, tb%grad)
         CASE (2)
            CALL tb%calc%dispersion%get_engrad(tb%mol, tb%dcache, tb%e_disp, &
            & tb%grad, tb%sigma)
         END SELECT
      END IF

#else
      CPABORT("Build without TBLITE")
#endif

   END SUBROUTINE tb_init_ham

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param energies ...
! **************************************************************************************************
   SUBROUTINE get_energy(tb, energies)

      TYPE(tb_type), POINTER                             :: tb
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: energies

#if defined(__TBLITE)

      energies(:) = tb%e_hal + tb%e_rep + tb%e_disp

#else
      CPABORT("Build without TBLITE")
#endif

   END SUBROUTINE get_energy

END MODULE tblite_interface

