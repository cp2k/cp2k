!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Define all structure types related to force field kinds
!> \par History
!>      10.2014 Moved kind types out of force_field_types.F [Ole Schuett]
!> \author Ole Schuett
! **************************************************************************************************
MODULE force_field_kind_types

   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_field_kind_types'

   INTEGER, PARAMETER, PUBLIC               :: do_ff_undef = 0, &
                                               do_ff_quartic = 1, &
                                               do_ff_g96 = 2, &
                                               do_ff_charmm = 3, &
                                               do_ff_harmonic = 4, &
                                               do_ff_g87 = 5, &
                                               do_ff_morse = 6, &
                                               do_ff_cubic = 7, &
                                               do_ff_mixed_bend_stretch = 8, &
                                               do_ff_amber = 9, &
                                               do_ff_mm2 = 10, &
                                               do_ff_mm3 = 11, &
                                               do_ff_mm4 = 12, &
                                               do_ff_fues = 13, &
                                               do_ff_legendre = 14, &
                                               do_ff_opls = 15

! *** Define the derived structure types ***

! **************************************************************************************************
   TYPE legendre_data_type
      INTEGER :: order = 0
      REAL(KIND=dp), DIMENSION(:), POINTER :: coeffs => NULL()
   END TYPE legendre_data_type

! **************************************************************************************************
   TYPE bond_kind_type
      INTEGER                        :: id_type = do_ff_undef
      REAL(KIND=dp)                :: k(3) = 0.0_dp, r0 = 0.0_dp, cs = 0.0_dp
      INTEGER                        :: kind_number = 0
   END TYPE bond_kind_type

! **************************************************************************************************
   TYPE bend_kind_type
      INTEGER                        :: id_type = do_ff_undef
      REAL(KIND=dp)                :: k = 0.0_dp, theta0 = 0.0_dp, cb = 0.0_dp
      REAL(KIND=dp)                :: r012 = 0.0_dp, r032 = 0.0_dp, kbs12 = 0.0_dp, kbs32 = 0.0_dp, kss = 0.0_dp
      TYPE(legendre_data_type)       :: legendre = legendre_data_type()
      INTEGER                        :: kind_number = 0
   END TYPE bend_kind_type

! **************************************************************************************************
   TYPE ub_kind_type
      INTEGER                        :: id_type = do_ff_undef
      REAL(KIND=dp)                :: k(3) = 0.0_dp, r0 = 0.0_dp
      INTEGER                        :: kind_number = 0
   END TYPE ub_kind_type

! **************************************************************************************************
   TYPE torsion_kind_type
      INTEGER                        :: id_type = do_ff_undef
      INTEGER                        :: nmul = 0
      INTEGER, POINTER               :: m(:) => NULL()
      REAL(KIND=dp), POINTER       :: k(:) => NULL(), phi0(:) => NULL()
      INTEGER                        :: kind_number = 0
   END TYPE torsion_kind_type

! **************************************************************************************************
   TYPE impr_kind_type
      INTEGER                        :: id_type = do_ff_undef
      REAL(KIND=dp)                :: k = 0.0_dp, phi0 = 0.0_dp
      INTEGER                        :: kind_number = 0
   END TYPE impr_kind_type

! **************************************************************************************************
   TYPE opbend_kind_type
      INTEGER                        :: id_type = do_ff_undef
      REAL(KIND=dp)                :: k = 0.0_dp, phi0 = 0.0_dp
      INTEGER                        :: kind_number = 0
   END TYPE opbend_kind_type

! *** Public subroutines ***

   PUBLIC :: allocate_bend_kind_set, &
             allocate_bond_kind_set, &
             allocate_ub_kind_set, &
             allocate_torsion_kind_set, &
             allocate_impr_kind_set, &
             allocate_opbend_kind_set, &
             deallocate_bend_kind_set, &
             deallocate_bond_kind_set, &
             torsion_kind_dealloc_ref, &
             impr_kind_dealloc_ref

! *** Public data types ***

   PUBLIC :: bend_kind_type, &
             bond_kind_type, &
             impr_kind_type, &
             torsion_kind_type, &
             opbend_kind_type, &
             ub_kind_type, &
             ub_kind_dealloc_ref, &
             legendre_data_type
CONTAINS

! **************************************************************************************************
!> \brief Allocate and initialize a bend kind set.
!> \param bend_kind_set ...
!> \param nkind ...
! **************************************************************************************************
   PURE SUBROUTINE allocate_bend_kind_set(bend_kind_set, nkind)

      TYPE(bend_kind_type), DIMENSION(:), POINTER        :: bend_kind_set
      INTEGER, INTENT(IN)                                :: nkind

      INTEGER                                            :: ikind

      NULLIFY (bend_kind_set)
      ALLOCATE (bend_kind_set(nkind))
      DO ikind = 1, nkind
         bend_kind_set(ikind)%kind_number = ikind
      END DO
   END SUBROUTINE allocate_bend_kind_set

! **************************************************************************************************
!> \brief Allocate and initialize a bond kind set.
!> \param bond_kind_set ...
!> \param nkind ...
! **************************************************************************************************
   PURE SUBROUTINE allocate_bond_kind_set(bond_kind_set, nkind)

      TYPE(bond_kind_type), DIMENSION(:), POINTER        :: bond_kind_set
      INTEGER, INTENT(IN)                                :: nkind

      INTEGER                                            :: ikind

      NULLIFY (bond_kind_set)
      ALLOCATE (bond_kind_set(nkind))
      DO ikind = 1, nkind
         bond_kind_set(ikind)%kind_number = ikind
      END DO
   END SUBROUTINE allocate_bond_kind_set

! **************************************************************************************************
!> \brief Allocate and initialize a torsion kind set.
!> \param torsion_kind_set ...
!> \param nkind ...
! **************************************************************************************************
   PURE SUBROUTINE allocate_torsion_kind_set(torsion_kind_set, nkind)

      TYPE(torsion_kind_type), DIMENSION(:), POINTER     :: torsion_kind_set
      INTEGER, INTENT(IN)                                :: nkind

      INTEGER                                            :: ikind

      NULLIFY (torsion_kind_set)
      ALLOCATE (torsion_kind_set(nkind))
      DO ikind = 1, nkind
         torsion_kind_set(ikind)%kind_number = ikind
      END DO
   END SUBROUTINE allocate_torsion_kind_set

! **************************************************************************************************
!> \brief Allocate and initialize a ub kind set.
!> \param ub_kind_set ...
!> \param nkind ...
! **************************************************************************************************
   PURE SUBROUTINE allocate_ub_kind_set(ub_kind_set, nkind)

      TYPE(ub_kind_type), DIMENSION(:), POINTER          :: ub_kind_set
      INTEGER, INTENT(IN)                                :: nkind

      INTEGER                                            :: ikind

      NULLIFY (ub_kind_set)
      ALLOCATE (ub_kind_set(nkind))
      DO ikind = 1, nkind
         ub_kind_set(ikind)%kind_number = ikind
      END DO
   END SUBROUTINE allocate_ub_kind_set

! **************************************************************************************************
!> \brief Allocate and initialize a impr kind set.
!> \param impr_kind_set ...
!> \param nkind ...
! **************************************************************************************************
   PURE SUBROUTINE allocate_impr_kind_set(impr_kind_set, nkind)

      TYPE(impr_kind_type), DIMENSION(:), POINTER        :: impr_kind_set
      INTEGER, INTENT(IN)                                :: nkind

      INTEGER                                            :: ikind

      NULLIFY (impr_kind_set)
      ALLOCATE (impr_kind_set(nkind))
      DO ikind = 1, nkind
         impr_kind_set(ikind)%kind_number = ikind
      END DO
   END SUBROUTINE allocate_impr_kind_set

! **************************************************************************************************
!> \brief Allocate and initialize a opbend kind set.
!> \param opbend_kind_set ...
!> \param nkind ...
! **************************************************************************************************
   PURE SUBROUTINE allocate_opbend_kind_set(opbend_kind_set, nkind)

      TYPE(opbend_kind_type), DIMENSION(:), POINTER      :: opbend_kind_set
      INTEGER, INTENT(IN)                                :: nkind

      INTEGER                                            :: ikind

      NULLIFY (opbend_kind_set)
      ALLOCATE (opbend_kind_set(nkind))
      DO ikind = 1, nkind
         opbend_kind_set(ikind)%kind_number = ikind
      END DO
   END SUBROUTINE allocate_opbend_kind_set

! **************************************************************************************************
!> \brief Deallocate a bend kind set.
!> \param bend_kind_set ...
! **************************************************************************************************
   PURE SUBROUTINE deallocate_bend_kind_set(bend_kind_set)

      TYPE(bend_kind_type), DIMENSION(:), POINTER        :: bend_kind_set

      INTEGER                                            :: i

      IF (.NOT. ASSOCIATED(bend_kind_set)) RETURN
      DO i = 1, SIZE(bend_kind_set)
         IF (ASSOCIATED(bend_kind_set(i)%legendre%coeffs)) THEN
            DEALLOCATE (bend_kind_set(i)%legendre%coeffs)
         END IF
      END DO
      DEALLOCATE (bend_kind_set)
   END SUBROUTINE deallocate_bend_kind_set

! **************************************************************************************************
!> \brief Deallocate a bond kind set.
!> \param bond_kind_set ...
! **************************************************************************************************
   PURE SUBROUTINE deallocate_bond_kind_set(bond_kind_set)

      TYPE(bond_kind_type), DIMENSION(:), POINTER        :: bond_kind_set

      DEALLOCATE (bond_kind_set)

   END SUBROUTINE deallocate_bond_kind_set

! **************************************************************************************************
!> \brief Deallocate a torsion kind element
!> \param torsion_kind ...
! **************************************************************************************************
   PURE SUBROUTINE torsion_kind_dealloc_ref(torsion_kind)

      TYPE(torsion_kind_type), INTENT(INOUT)             :: torsion_kind

      IF (ASSOCIATED(torsion_kind%k)) THEN
         DEALLOCATE (torsion_kind%k)
      END IF
      IF (ASSOCIATED(torsion_kind%m)) THEN
         DEALLOCATE (torsion_kind%m)
      END IF
      IF (ASSOCIATED(torsion_kind%phi0)) THEN
         DEALLOCATE (torsion_kind%phi0)
      END IF

   END SUBROUTINE torsion_kind_dealloc_ref

! **************************************************************************************************
!> \brief Deallocate a ub kind set.
!> \param ub_kind_set ...
! **************************************************************************************************
   PURE SUBROUTINE ub_kind_dealloc_ref(ub_kind_set)
      TYPE(ub_kind_type), DIMENSION(:), POINTER          :: ub_kind_set

      DEALLOCATE (ub_kind_set)

   END SUBROUTINE ub_kind_dealloc_ref

! **************************************************************************************************
!> \brief Deallocate a impr kind element
! **************************************************************************************************
   PURE SUBROUTINE impr_kind_dealloc_ref()

!
! Questa e' la migliore routine che mente umana abbia concepito! ;-)
! Translation to english: This is the best subroutine that humanity can imagine! ;-)
!

   END SUBROUTINE impr_kind_dealloc_ref

END MODULE force_field_kind_types
