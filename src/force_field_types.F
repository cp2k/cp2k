!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/force_field_types [1.0] *
!!
!!   NAME
!!     force_field_types
!!
!!   FUNCTION
!!     Define all structures types related to force_fields
!!
!!MK   SHOULD BE MERGED WITH MODULE force_fields.
!!
!!   AUTHOR
!!     MK (29.08.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE force_field_types

  USE kinds,                           ONLY: dp,&
                                             int_size
  USE pair_potential_types,            ONLY: pair_potential_type
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "force_field_types"
  INTEGER, PARAMETER          :: max_name_length = 60

! *** Define the derived structure types ***

  TYPE bond_kind_type
    CHARACTER(LEN=max_name_length) :: name
    REAL(KIND = dp)                :: k,r0
    INTEGER                        :: kind_number
  END TYPE bond_kind_type

  TYPE bend_kind_type
    CHARACTER(LEN=max_name_length) :: name
    REAL(KIND = dp)                :: k,theta0
    INTEGER                        :: kind_number
  END TYPE bend_kind_type

  TYPE ub_kind_type
    CHARACTER(LEN=max_name_length) :: name
    REAL(KIND = dp)                :: k,r0
    INTEGER                        :: kind_number
  END TYPE ub_kind_type

  TYPE torsion_kind_type
    CHARACTER(LEN=max_name_length) :: name
    REAL(KIND = dp)                :: k,phi0
    REAL(KIND = dp), DIMENSION ( 0:3 ) :: a
    INTEGER                        :: kind_number
  END TYPE torsion_kind_type

  TYPE impr_kind_type
    CHARACTER(LEN=max_name_length) :: name
    REAL(KIND = dp)                :: k,phi0
    INTEGER                        :: kind_number
  END TYPE impr_kind_type

  TYPE onfo_kind_type
    CHARACTER(LEN=max_name_length) :: name
    TYPE ( pair_potential_type )   :: potparm
    INTEGER                        :: kind_number
  END TYPE onfo_kind_type

! *** Public subroutines ***

  PUBLIC :: allocate_bend_kind_set,&
            allocate_bond_kind_set,&
            allocate_ub_kind_set,&
            allocate_torsion_kind_set,&
            allocate_impr_kind_set,&
            allocate_onfo_kind_set

! *** Public data types ***

  PUBLIC :: bend_kind_type,&
            bond_kind_type,&
            impr_kind_type,&
            onfo_kind_type,&
            torsion_kind_type,&
            ub_kind_type

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_bend_kind_set(bend_kind_set,nkind)

!   Purpose: Allocate and initialize a bend kind set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(bend_kind_type), DIMENSION(:), &
      POINTER                                :: bend_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_bend_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(bend_kind_set)) THEN
      DEALLOCATE (bend_kind_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "bend_kind_set")
    END IF

    ALLOCATE (bend_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "bend_kind_set",nkind*int_size)

    DO ikind=1,nkind
      bend_kind_set(ikind)%name = ""
      bend_kind_set(ikind)%k = 0.0_dp
      bend_kind_set(ikind)%theta0 = 0.0_dp
      bend_kind_set(ikind)%kind_number = ikind
    END DO

  END SUBROUTINE allocate_bend_kind_set

! *****************************************************************************

  SUBROUTINE allocate_bond_kind_set(bond_kind_set,nkind)

!   Purpose: Allocate and initialize a bond kind set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(bond_kind_type), DIMENSION(:), &
      POINTER                                :: bond_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_bond_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(bond_kind_set)) THEN
      DEALLOCATE (bond_kind_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "bond_kind_set")
    END IF

    ALLOCATE (bond_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "bond_kind_set",nkind*int_size)

    DO ikind=1,nkind
      bond_kind_set(ikind)%name = ""
      bond_kind_set(ikind)%k = 0.0_dp
      bond_kind_set(ikind)%r0 = 0.0_dp
      bond_kind_set(ikind)%kind_number = ikind
    END DO

  END SUBROUTINE allocate_bond_kind_set

! *****************************************************************************

  SUBROUTINE allocate_torsion_kind_set(torsion_kind_set,nkind)

!   Purpose: Allocate and initialize a torsion kind set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(torsion_kind_type), DIMENSION(:), &
      POINTER                                :: torsion_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_torsion_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(torsion_kind_set)) THEN
      DEALLOCATE (torsion_kind_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "torsion_kind_set")
    END IF

    ALLOCATE (torsion_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "torsion_kind_set",nkind*int_size)

    DO ikind=1,nkind
      torsion_kind_set(ikind)%name = ""
      torsion_kind_set(ikind)%k = 0.0_dp
      torsion_kind_set(ikind)%phi0 = 0.0_dp
      torsion_kind_set(ikind)%kind_number = ikind
    END DO

  END SUBROUTINE allocate_torsion_kind_set

! *****************************************************************************

  SUBROUTINE allocate_ub_kind_set(ub_kind_set,nkind)

!   Purpose: Allocate and initialize a ub kind set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(ub_kind_type), DIMENSION(:), &
      POINTER                                :: ub_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_ub_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(ub_kind_set)) THEN
      DEALLOCATE (ub_kind_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "ub_kind_set")
    END IF

    ALLOCATE (ub_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "ub_kind_set",nkind*int_size)

    DO ikind=1,nkind
      ub_kind_set(ikind)%name = ""
      ub_kind_set(ikind)%k = 0.0_dp
      ub_kind_set(ikind)%r0 = 0.0_dp
      ub_kind_set(ikind)%kind_number = ikind
    END DO

  END SUBROUTINE allocate_ub_kind_set

! *****************************************************************************

  SUBROUTINE allocate_impr_kind_set(impr_kind_set,nkind)

!   Purpose: Allocate and initialize a impr kind set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(impr_kind_type), DIMENSION(:), &
      POINTER                                :: impr_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_impr_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(impr_kind_set)) THEN
      DEALLOCATE (impr_kind_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "impr_kind_set")
    END IF

    ALLOCATE (impr_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "impr_kind_set",nkind*int_size)

    DO ikind=1,nkind
      impr_kind_set(ikind)%name = ""
      impr_kind_set(ikind)%k = 0.0_dp
      impr_kind_set(ikind)%phi0 = 0.0_dp
      impr_kind_set(ikind)%kind_number = ikind
    END DO

  END SUBROUTINE allocate_impr_kind_set

! *****************************************************************************

  SUBROUTINE allocate_onfo_kind_set(onfo_kind_set,nkind)

!   Purpose: Allocate and initialize a onfo kind set.

!   History: - Creation (29.08.2003,MK)

!   ***************************************************************************

    TYPE(onfo_kind_type), DIMENSION(:), &
      POINTER                                :: onfo_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_onfo_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(onfo_kind_set)) THEN
      DEALLOCATE (onfo_kind_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "onfo_kind_set")
    END IF

    ALLOCATE (onfo_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "onfo_kind_set",nkind*int_size)

    DO ikind=1,nkind
      onfo_kind_set(ikind)%name = ""
      onfo_kind_set(ikind)%kind_number = ikind
    END DO

  END SUBROUTINE allocate_onfo_kind_set

! *****************************************************************************

END MODULE force_field_types
