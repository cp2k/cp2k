!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Define all structure types related to force field kinds
!> \par History
!>      10.2014 Moved kind types out of force_field_types.F [Ole Schuett]
!> \author Ole Schuett
! *****************************************************************************
MODULE force_field_kind_types
  
  USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_field_kind_types'

  INTEGER, PARAMETER, PUBLIC               :: do_ff_undef=0,&
                                              do_ff_quartic=1,&
                                              do_ff_g96=2,&
                                              do_ff_charmm=3,&
                                              do_ff_harmonic=4,&
                                              do_ff_g87=5,&
                                              do_ff_morse=6,&
                                              do_ff_cubic=7,&
                                              do_ff_mixed_bend_stretch=8,&
                                              do_ff_amber=9,&
                                              do_ff_mm2=10,&
                                              do_ff_mm3=11,&
                                              do_ff_mm4=12,&
                                              do_ff_fues=13,&
                                              do_ff_legendre=14,&
                                              do_ff_opls=15

! *** Define the derived structure types ***

! *****************************************************************************
  TYPE legendre_data_type
    INTEGER :: order
    REAL(KIND=dp), DIMENSION(:), POINTER :: coeffs
  END TYPE legendre_data_type

! *****************************************************************************
  TYPE bond_kind_type
     INTEGER                        :: id_type
     REAL(KIND = dp)                :: k(3),r0,cs
     INTEGER                        :: kind_number
  END TYPE bond_kind_type

! *****************************************************************************
  TYPE bend_kind_type
     INTEGER                        :: id_type
     REAL(KIND = dp)                :: k,theta0,cb
     REAL(KIND = dp)                :: r012,r032,kbs12,kbs32,kss
     TYPE(legendre_data_type)       :: legendre
     INTEGER                        :: kind_number
  END TYPE bend_kind_type

! *****************************************************************************
  TYPE ub_kind_type
     INTEGER                        :: id_type
     REAL(KIND = dp)                :: k(3),r0
     INTEGER                        :: kind_number
  END TYPE ub_kind_type

! *****************************************************************************
  TYPE torsion_kind_type
     INTEGER                        :: id_type
     INTEGER                        :: nmul
     INTEGER, POINTER               :: m(:)
     REAL(KIND = dp), POINTER       :: k(:),phi0(:)
     INTEGER                        :: kind_number
  END TYPE torsion_kind_type

! *****************************************************************************
  TYPE impr_kind_type
     INTEGER                        :: id_type
     REAL(KIND = dp)                :: k,phi0
     INTEGER                        :: kind_number
  END TYPE impr_kind_type

! *****************************************************************************
  TYPE opbend_kind_type
     INTEGER                        :: id_type
     REAL(KIND = dp)                :: k,phi0
     INTEGER                        :: kind_number
  END TYPE opbend_kind_type

! *** Public subroutines ***

  PUBLIC :: allocate_bend_kind_set,&
            allocate_bond_kind_set,&
            allocate_ub_kind_set,&
            allocate_torsion_kind_set,&
            allocate_impr_kind_set,&
            allocate_opbend_kind_set,&
            deallocate_bend_kind_set,&
            deallocate_bond_kind_set,&
            torsion_kind_dealloc_ref,&
            impr_kind_dealloc_ref

! *** Public data types ***

  PUBLIC :: bend_kind_type,&
            bond_kind_type,&
            impr_kind_type,&
            torsion_kind_type,&
            opbend_kind_type,&
            ub_kind_type,&
            ub_kind_dealloc_ref,&
            legendre_data_type
CONTAINS

! *****************************************************************************
!> \brief Allocate and initialize a bend kind set.
!> \param bend_kind_set ...
!> \param nkind ...
! *****************************************************************************
  SUBROUTINE allocate_bend_kind_set(bend_kind_set,nkind)

    TYPE(bend_kind_type), DIMENSION(:), &
      POINTER                                :: bend_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_bend_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind

    NULLIFY(bend_kind_set)
    ALLOCATE (bend_kind_set(nkind))
    DO ikind=1,nkind
       bend_kind_set(ikind)%id_type = do_ff_undef
       bend_kind_set(ikind)%k = 0.0_dp
       bend_kind_set(ikind)%theta0 = 0.0_dp
       bend_kind_set(ikind)%cb = 0.0_dp
       bend_kind_set(ikind)%r012 = 0.0_dp
       bend_kind_set(ikind)%r032 = 0.0_dp
       bend_kind_set(ikind)%kbs12 = 0.0_dp
       bend_kind_set(ikind)%kbs32 = 0.0_dp
       bend_kind_set(ikind)%kss = 0.0_dp
       bend_kind_set(ikind)%legendre%order = 0
       NULLIFY(bend_kind_set(ikind)%legendre%coeffs)
       bend_kind_set(ikind)%kind_number = ikind
    END DO
  END SUBROUTINE allocate_bend_kind_set

! *****************************************************************************
!> \brief Allocate and initialize a bond kind set.
!> \param bond_kind_set ...
!> \param nkind ...
! *****************************************************************************
  SUBROUTINE allocate_bond_kind_set(bond_kind_set,nkind)

    TYPE(bond_kind_type), DIMENSION(:), &
      POINTER                                :: bond_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_bond_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind

    NULLIFY(bond_kind_set)
    ALLOCATE (bond_kind_set(nkind))
    DO ikind=1,nkind
       bond_kind_set(ikind)%id_type = do_ff_undef
       bond_kind_set(ikind)%k(:) = 0.0_dp
       bond_kind_set(ikind)%r0 = 0.0_dp
       bond_kind_set(ikind)%cs = 0.0_dp
       bond_kind_set(ikind)%kind_number = ikind
    END DO
  END SUBROUTINE allocate_bond_kind_set

! *****************************************************************************
!> \brief Allocate and initialize a torsion kind set.
!> \param torsion_kind_set ...
!> \param nkind ...
! *****************************************************************************
  SUBROUTINE allocate_torsion_kind_set(torsion_kind_set,nkind)

    TYPE(torsion_kind_type), DIMENSION(:), &
      POINTER                                :: torsion_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_torsion_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind

    NULLIFY(torsion_kind_set)
    ALLOCATE (torsion_kind_set(nkind))

    DO ikind=1,nkind
       torsion_kind_set(ikind)%id_type = do_ff_undef
       torsion_kind_set(ikind)%nmul = 0
       NULLIFY(torsion_kind_set(ikind)%k)
       NULLIFY(torsion_kind_set(ikind)%m)
       NULLIFY(torsion_kind_set(ikind)%phi0)
       torsion_kind_set(ikind)%kind_number = ikind
    END DO
  END SUBROUTINE allocate_torsion_kind_set

! *****************************************************************************
!> \brief Allocate and initialize a ub kind set.
!> \param ub_kind_set ...
!> \param nkind ...
! *****************************************************************************
  SUBROUTINE allocate_ub_kind_set(ub_kind_set,nkind)

    TYPE(ub_kind_type), DIMENSION(:), &
      POINTER                                :: ub_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_ub_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind

    NULLIFY(ub_kind_set)
    ALLOCATE (ub_kind_set(nkind))
    DO ikind=1,nkind
      ub_kind_set(ikind)%id_type = do_ff_undef
      ub_kind_set(ikind)%k = 0.0_dp
      ub_kind_set(ikind)%r0 = 0.0_dp
      ub_kind_set(ikind)%kind_number = ikind
    END DO
  END SUBROUTINE allocate_ub_kind_set

! *****************************************************************************
!> \brief Allocate and initialize a impr kind set.
!> \param impr_kind_set ...
!> \param nkind ...
! *****************************************************************************
  SUBROUTINE allocate_impr_kind_set(impr_kind_set,nkind)

    TYPE(impr_kind_type), DIMENSION(:), &
      POINTER                                :: impr_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_impr_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind

    NULLIFY(impr_kind_set)
    ALLOCATE (impr_kind_set(nkind))
    DO ikind=1,nkind
      impr_kind_set(ikind)%id_type = do_ff_undef
      impr_kind_set(ikind)%k = 0.0_dp
      impr_kind_set(ikind)%phi0 = 0.0_dp
      impr_kind_set(ikind)%kind_number = ikind
    END DO
  END SUBROUTINE allocate_impr_kind_set

! *****************************************************************************
!> \brief Allocate and initialize a opbend kind set.
!> \param opbend_kind_set ...
!> \param nkind ...
! *****************************************************************************
  SUBROUTINE allocate_opbend_kind_set(opbend_kind_set,nkind)

    TYPE(opbend_kind_type), DIMENSION(:), &
      POINTER                                :: opbend_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_opbend_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind

    NULLIFY(opbend_kind_set)
    ALLOCATE (opbend_kind_set(nkind))
    DO ikind=1,nkind
      opbend_kind_set(ikind)%id_type = do_ff_undef
      opbend_kind_set(ikind)%k = 0.0_dp
      opbend_kind_set(ikind)%phi0 = 0.0_dp
      opbend_kind_set(ikind)%kind_number = ikind
    END DO
  END SUBROUTINE allocate_opbend_kind_set

! *****************************************************************************
!> \brief Deallocate a bend kind set.
!> \param bend_kind_set ...
! *****************************************************************************
  SUBROUTINE deallocate_bend_kind_set(bend_kind_set)

    TYPE(bend_kind_type), DIMENSION(:), &
      POINTER                                :: bend_kind_set

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_bend_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF (.NOT.ASSOCIATED(bend_kind_set)) RETURN
    DO i = 1, SIZE(bend_kind_set)
      IF (ASSOCIATED(bend_kind_set(i)%legendre%coeffs)) THEN
         DEALLOCATE(bend_kind_set(i)%legendre%coeffs)
         NULLIFY(bend_kind_set(i)%legendre%coeffs)
      END IF
    END DO
    DEALLOCATE (bend_kind_set)
  END SUBROUTINE deallocate_bend_kind_set

! *****************************************************************************
!> \brief Deallocate a bond kind set.
!> \param bond_kind_set ...
! *****************************************************************************
  SUBROUTINE deallocate_bond_kind_set(bond_kind_set)

    TYPE(bond_kind_type), DIMENSION(:), &
      POINTER                                :: bond_kind_set

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_bond_kind_set', &
      routineP = moduleN//':'//routineN

    DEALLOCATE (bond_kind_set)

  END SUBROUTINE deallocate_bond_kind_set

! *****************************************************************************
!> \brief Deallocate a torsion kind element
!> \param torsion_kind ...
! *****************************************************************************
  SUBROUTINE torsion_kind_dealloc_ref(torsion_kind)

    TYPE(torsion_kind_type), INTENT(INOUT)   :: torsion_kind

    CHARACTER(len=*), PARAMETER :: routineN = 'torsion_kind_dealloc_ref', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(torsion_kind%k)) THEN
       DEALLOCATE(torsion_kind%k)
    END IF
    IF(ASSOCIATED(torsion_kind%m)) THEN
       DEALLOCATE(torsion_kind%m)
    END IF
    IF(ASSOCIATED(torsion_kind%phi0)) THEN
       DEALLOCATE(torsion_kind%phi0)
    END IF

  END SUBROUTINE torsion_kind_dealloc_ref

! *****************************************************************************
!> \brief Deallocate a ub kind set.
!> \param ub_kind_set ...
! *****************************************************************************
  SUBROUTINE ub_kind_dealloc_ref(ub_kind_set)
    TYPE(ub_kind_type), DIMENSION(:), &
      POINTER                                :: ub_kind_set

    CHARACTER(len=*), PARAMETER :: routineN = 'ub_kind_dealloc_ref', &
      routineP = moduleN//':'//routineN

    DEALLOCATE(ub_kind_set)

  END SUBROUTINE ub_kind_dealloc_ref

! *****************************************************************************
!> \brief Deallocate a impr kind element
! *****************************************************************************
  SUBROUTINE impr_kind_dealloc_ref()


    CHARACTER(len=*), PARAMETER :: routineN = 'impr_kind_dealloc_ref', &
      routineP = moduleN//':'//routineN

!
! Questa e' la migliore routine che mente umana abbia concepito! ;-)
! Translation to english: This is the best subroutine that humanity can imagine! ;-)
!

  END SUBROUTINE impr_kind_dealloc_ref

END MODULE force_field_kind_types
