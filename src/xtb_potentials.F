!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief xTB (repulsive) pair potentials
!>        Reference: Stefan Grimme, Christoph Bannwarth, Philip Shushkov
!>                   JCTC 13, 1989-2009, (2017)
!>                   DOI: 10.1021/acs.jctc.7b00118
!> \author JGH
! **************************************************************************************************
MODULE xtb_potentials
   USE kinds,                           ONLY: dp
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   TYPE neighbor_atoms_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER     :: coord
      REAL(KIND=dp), DIMENSION(:), POINTER        :: rab
      INTEGER, DIMENSION(:), POINTER              :: katom
   END TYPE neighbor_atoms_type

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_potentials'

   PUBLIC :: xtb_pp_radius

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param ppradius ...
!> \param eps_pair ...
!> \param kfparam ...
! **************************************************************************************************
   SUBROUTINE xtb_pp_radius(qs_kind_set, ppradius, eps_pair, kfparam)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: ppradius
      REAL(KIND=dp), INTENT(IN)                          :: eps_pair, kfparam

      INTEGER                                            :: ikind, ir, jkind, nkind
      LOGICAL                                            :: defa, defb
      REAL(KIND=dp)                                      :: alphaa, alphab, erep, rab, rab0, rcova, &
                                                            rcovb, saa, zneffa, zneffb
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b

      ppradius = 0.0_dp
      nkind = SIZE(ppradius, 1)
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, rcov=rcova, alpha=alphaa, zneff=zneffa, defined=defa)
         IF (.NOT. defa) CYCLE
         DO jkind = ikind, nkind
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, rcov=rcovb, alpha=alphab, zneff=zneffb, defined=defb)
            IF (.NOT. defb) CYCLE
            rab = 0.0_dp
            DO ir = 1, 24
               rab = rab + 1.0_dp
               saa = SQRT(alphaa*alphab)
               erep = zneffa*zneffb/rab*EXP(-saa*rab**kfparam)
               IF (erep < eps_pair) EXIT
            END DO
            rab0 = rcova + rcovb
            rab = MAX(rab, rab0 + 2.0_dp)
!deb
            WRITE (6, *) "XTB RADIUS  ", ikind, jkind, rab
!deb
            ppradius(ikind, jkind) = rab
            ppradius(jkind, ikind) = ppradius(ikind, jkind)
         END DO
      END DO

   END SUBROUTINE xtb_pp_radius

END MODULE xtb_potentials

