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
                                             dbl,&
                                             default_string_length,&
                                             int_size
  USE pair_potential_types,            ONLY: pair_potential_type
  USE splines,                         ONLY: spline_environment_type
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "force_field_types"

! *** Define the derived structure types ***

  TYPE bond_kind_type
    CHARACTER(LEN=default_string_length) :: name
    REAL(KIND = dp)                :: k,r0
    INTEGER                        :: kind_number
  END TYPE bond_kind_type

  TYPE bend_kind_type
    CHARACTER(LEN=default_string_length) :: name
    REAL(KIND = dp)                :: k,theta0
    INTEGER                        :: kind_number
  END TYPE bend_kind_type

  TYPE ub_kind_type
    CHARACTER(LEN=default_string_length) :: name
    REAL(KIND = dp)                :: k,r0
    INTEGER                        :: kind_number
  END TYPE ub_kind_type

  TYPE torsion_kind_type
    CHARACTER(LEN=default_string_length) :: name
    INTEGER                        :: m
    REAL(KIND = dp)                :: k,cosphi0
    REAL(KIND = dp), DIMENSION ( 0:3 ) :: a
    INTEGER                        :: kind_number
  END TYPE torsion_kind_type

  TYPE impr_kind_type
    CHARACTER(LEN=default_string_length) :: name
    REAL(KIND = dp)                :: k,phi0
    INTEGER                        :: kind_number
  END TYPE impr_kind_type

  TYPE onfo_kind_type
    CHARACTER(LEN=default_string_length)                  :: name
    TYPE ( spline_environment_type ), POINTER       :: spline_env
    INTEGER                                         :: kind_number
  END TYPE onfo_kind_type

  TYPE input_info_type
    CHARACTER (LEN=default_string_length), POINTER :: charge_atm(:)
    REAL (dbl), POINTER                            :: charge(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_b(:)
    REAL (dbl), POINTER                            :: bond_k(:)
    REAL (dbl), POINTER                            :: bond_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_c(:)
    REAL (dbl), POINTER                            :: bend_k(:)
    REAL (dbl), POINTER                            :: bend_theta0(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_c(:)
    REAL (dbl), POINTER                            :: ub_k(:)
    REAL (dbl), POINTER                            :: ub_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_d(:)
    REAL (dbl), POINTER                            :: torsion_k(:)
    REAL (dbl), POINTER                            :: torsion_m(:)
    REAL (dbl), POINTER                            :: torsion_cosphi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_d(:)
    REAL (dbl), POINTER                            :: impr_k(:)
    REAL (dbl), POINTER                            :: impr_phi0(:)
    TYPE (pair_potential_type), POINTER            :: nonbonded(:)
    TYPE (pair_potential_type), POINTER            :: nonbonded14(:)
  END TYPE input_info_type

  TYPE charmm_info_type
    CHARACTER (LEN=default_string_length), POINTER :: charge_atm(:)
    REAL (dbl), POINTER                            :: charge(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_b(:)
    REAL (dbl), POINTER                            :: bond_k(:)
    REAL (dbl), POINTER                            :: bond_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_c(:)
    REAL (dbl), POINTER                            :: bend_k(:)
    REAL (dbl), POINTER                            :: bend_theta0(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_c(:)
    REAL (dbl), POINTER                            :: ub_k(:)
    REAL (dbl), POINTER                            :: ub_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_d(:)
    REAL (dbl), POINTER                            :: torsion_k(:)
    INTEGER, POINTER                               :: torsion_m(:)
    REAL (dbl), POINTER                            :: torsion_cosphi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_d(:)
    REAL (dbl), POINTER                            :: impr_k(:)
    REAL (dbl), POINTER                            :: impr_phi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: nonbond_a(:)
    REAL (dbl), POINTER                            :: nonbond_eps(:)
    REAL (dbl), POINTER                            :: nonbond_rmin2(:)
    REAL (dbl), POINTER                            :: nonbond_eps_14(:)
    REAL (dbl), POINTER                            :: nonbond_rmin2_14(:)
  END TYPE charmm_info_type

  TYPE amber_info_type
    INTEGER                                   :: nothing_here_yet
  END TYPE amber_info_type

  TYPE force_field_type
    TYPE (input_info_type),POINTER   :: inp_info
    TYPE (charmm_info_type),POINTER  :: chm_info
    TYPE (amber_info_type),POINTER   :: amb_info
    CHARACTER (LEN=default_string_length)     :: ff_file_name
    CHARACTER (LEN=default_string_length)     :: ff_type
  END TYPE force_field_type

! *** Public subroutines ***

  PUBLIC :: allocate_bend_kind_set,&
            allocate_bond_kind_set,&
            allocate_ub_kind_set,&
            allocate_torsion_kind_set,&
            allocate_impr_kind_set,&
            allocate_onfo_kind_set,&
            init_ff_type,&
            deallocate_ff_type

! *** Public data types ***

  PUBLIC :: bend_kind_type,&
            bond_kind_type,&
            impr_kind_type,&
            onfo_kind_type,&
            torsion_kind_type,&
            ub_kind_type,&
            amber_info_type,&
            charmm_info_type,&
            input_info_type,&
            force_field_type

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

    ALLOCATE (torsion_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "torsion_kind_set",nkind*int_size)

    DO ikind=1,nkind
      torsion_kind_set(ikind)%name = ""
      torsion_kind_set(ikind)%k = 0.0_dp
      torsion_kind_set(ikind)%m = 0
      torsion_kind_set(ikind)%cosphi0 = 0.0_dp
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

    ALLOCATE (onfo_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "onfo_kind_set",nkind*int_size)

    DO ikind=1,nkind
      onfo_kind_set(ikind)%name = ""
      onfo_kind_set(ikind)%kind_number = ikind
    END DO

  END SUBROUTINE allocate_onfo_kind_set

! *****************************************************************************






!******************************************************************************
!!****** force_field_types/init_force_field [1.0] *
!!
!!   NAME
!!     init_force_field
!!
!!   FUNCTION
!!     1. Just NULLIFY and zero all the stuff
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE init_ff_type (ff_type)
    TYPE(force_field_type), &
      INTENT(INOUT)                          :: ff_type

    INTEGER                                  :: istat

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Initialize 
  !-----------------------------------------------------------------------------
  ALLOCATE(ff_type%inp_info,STAT=istat)
  IF(istat/=0) CALL stop_memory ('init_ff_type','inp_info')
  ALLOCATE(ff_type%chm_info,STAT=istat)
  IF(istat/=0) CALL stop_memory ('init_ff_type','chm_info')
  ALLOCATE(ff_type%amb_info,STAT=istat)
  IF(istat/=0) CALL stop_memory ('init_ff_type','amb_info')


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Initialize and Nullify things in ff_type
  !-----------------------------------------------------------------------------
  ff_type%ff_file_name=''
  ff_type%ff_type='OFF'


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Initialize and Nullify things in ff_type%inp_info
  !-----------------------------------------------------------------------------
  NULLIFY(ff_type%inp_info%charge_atm)
  NULLIFY(ff_type%inp_info%charge)
  NULLIFY(ff_type%inp_info%bond_a)
  NULLIFY(ff_type%inp_info%bond_b)
  NULLIFY(ff_type%inp_info%bond_k)
  NULLIFY(ff_type%inp_info%bond_r0)
  NULLIFY(ff_type%inp_info%bend_a)
  NULLIFY(ff_type%inp_info%bend_b)
  NULLIFY(ff_type%inp_info%bend_c)
  NULLIFY(ff_type%inp_info%bend_k)
  NULLIFY(ff_type%inp_info%bend_theta0)
  NULLIFY(ff_type%inp_info%ub_a)
  NULLIFY(ff_type%inp_info%ub_b)
  NULLIFY(ff_type%inp_info%ub_c)
  NULLIFY(ff_type%inp_info%ub_k)
  NULLIFY(ff_type%inp_info%ub_r0)
  NULLIFY(ff_type%inp_info%torsion_a)
  NULLIFY(ff_type%inp_info%torsion_b)
  NULLIFY(ff_type%inp_info%torsion_c)
  NULLIFY(ff_type%inp_info%torsion_d)
  NULLIFY(ff_type%inp_info%torsion_k)
  NULLIFY(ff_type%inp_info%torsion_m)
  NULLIFY(ff_type%inp_info%torsion_cosphi0)
  NULLIFY(ff_type%inp_info%impr_a)
  NULLIFY(ff_type%inp_info%impr_b)
  NULLIFY(ff_type%inp_info%impr_c)
  NULLIFY(ff_type%inp_info%impr_d)
  NULLIFY(ff_type%inp_info%impr_k)
  NULLIFY(ff_type%inp_info%impr_phi0)
  NULLIFY(ff_type%inp_info%nonbonded)
  NULLIFY(ff_type%inp_info%nonbonded14)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. Initialize and Nullify things in ff_type%chm_info
  !-----------------------------------------------------------------------------
  NULLIFY(ff_type%chm_info%charge_atm)
  NULLIFY(ff_type%chm_info%charge)
  NULLIFY(ff_type%chm_info%bond_a)
  NULLIFY(ff_type%chm_info%bond_b)
  NULLIFY(ff_type%chm_info%bond_k)
  NULLIFY(ff_type%chm_info%bond_r0)
  NULLIFY(ff_type%chm_info%bend_a)
  NULLIFY(ff_type%chm_info%bend_b)
  NULLIFY(ff_type%chm_info%bend_c)
  NULLIFY(ff_type%chm_info%bend_k)
  NULLIFY(ff_type%chm_info%bend_theta0)
  NULLIFY(ff_type%chm_info%ub_a)
  NULLIFY(ff_type%chm_info%ub_b)
  NULLIFY(ff_type%chm_info%ub_c)
  NULLIFY(ff_type%chm_info%ub_k)
  NULLIFY(ff_type%chm_info%ub_r0)
  NULLIFY(ff_type%chm_info%torsion_a)
  NULLIFY(ff_type%chm_info%torsion_b)
  NULLIFY(ff_type%chm_info%torsion_c)
  NULLIFY(ff_type%chm_info%torsion_d)
  NULLIFY(ff_type%chm_info%torsion_k)
  NULLIFY(ff_type%chm_info%torsion_m)
  NULLIFY(ff_type%chm_info%torsion_cosphi0)
  NULLIFY(ff_type%chm_info%impr_a)
  NULLIFY(ff_type%chm_info%impr_b)
  NULLIFY(ff_type%chm_info%impr_c)
  NULLIFY(ff_type%chm_info%impr_d)
  NULLIFY(ff_type%chm_info%impr_k)
  NULLIFY(ff_type%chm_info%impr_phi0)
  NULLIFY(ff_type%chm_info%nonbond_a)
  NULLIFY(ff_type%chm_info%nonbond_eps)
  NULLIFY(ff_type%chm_info%nonbond_rmin2)
  NULLIFY(ff_type%chm_info%nonbond_eps_14)
  NULLIFY(ff_type%chm_info%nonbond_rmin2_14)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 8. Initialize and Nullify things in ff_type%amb_info
  !-----------------------------------------------------------------------------



!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE init_ff_type



!******************************************************************************
!!****** force_field_types/deallocate_ff_type [1.0] *
!!
!!   NAME
!!     deallocate_ff_type
!!
!!   FUNCTION
!!     1. Just DEALLOCATE all the stuff
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE deallocate_ff_type (ff_type)
    TYPE(force_field_type), &
      INTENT(INOUT)                          :: ff_type

    INTEGER                                  :: istat

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. DEALLOCATE things in ff_type%inp_info
  !-----------------------------------------------------------------------------
  IF(ASSOCIATED(ff_type%inp_info%charge_atm)) THEN
    DEALLOCATE(ff_type%inp_info%charge_atm,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','charge_atm')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%charge)) THEN
    DEALLOCATE(ff_type%inp_info%charge,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','charge')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%bond_a)) THEN
    DEALLOCATE(ff_type%inp_info%bond_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bond_a')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%bond_b)) THEN
    DEALLOCATE(ff_type%inp_info%bond_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bond_b')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%bond_k)) THEN
    DEALLOCATE(ff_type%inp_info%bond_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bond_k')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%bond_r0)) THEN
    DEALLOCATE(ff_type%inp_info%bond_r0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bond_r0')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%bend_a)) THEN
    DEALLOCATE(ff_type%inp_info%bend_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bend_a')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%bend_b)) THEN
    DEALLOCATE(ff_type%inp_info%bend_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bend_b')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%bend_c)) THEN
    DEALLOCATE(ff_type%inp_info%bend_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bend_c')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%bend_k)) THEN
    DEALLOCATE(ff_type%inp_info%bend_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bend_k')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%bend_theta0)) THEN
    DEALLOCATE(ff_type%inp_info%bend_theta0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bend_theta0')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%ub_a)) THEN
    DEALLOCATE(ff_type%inp_info%ub_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','ub_a')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%ub_b)) THEN
    DEALLOCATE(ff_type%inp_info%ub_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','ub_b')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%ub_c)) THEN
    DEALLOCATE(ff_type%inp_info%ub_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','ub_c')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%ub_k)) THEN
    DEALLOCATE(ff_type%inp_info%ub_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','ub_k')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%ub_r0)) THEN
    DEALLOCATE(ff_type%inp_info%ub_r0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','ub_r0')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%torsion_a)) THEN
    DEALLOCATE(ff_type%inp_info%torsion_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_a')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%torsion_b)) THEN
    DEALLOCATE(ff_type%inp_info%torsion_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_b')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%torsion_c)) THEN
    DEALLOCATE(ff_type%inp_info%torsion_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_c')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%torsion_d)) THEN
    DEALLOCATE(ff_type%inp_info%torsion_d,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_d')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%torsion_k)) THEN
    DEALLOCATE(ff_type%inp_info%torsion_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_k')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%torsion_m)) THEN
    DEALLOCATE(ff_type%inp_info%torsion_m,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_m')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%torsion_cosphi0)) THEN
    DEALLOCATE(ff_type%inp_info%torsion_cosphi0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_cosphi0')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%impr_a)) THEN
    DEALLOCATE(ff_type%inp_info%impr_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','impr_a')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%impr_b)) THEN
    DEALLOCATE(ff_type%inp_info%impr_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','impr_b')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%impr_c)) THEN
    DEALLOCATE(ff_type%inp_info%impr_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','impr_c')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%impr_d)) THEN
    DEALLOCATE(ff_type%inp_info%impr_d,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','impr_d')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%impr_k)) THEN
    DEALLOCATE(ff_type%inp_info%impr_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','impr_k')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%impr_phi0)) THEN
    DEALLOCATE(ff_type%inp_info%impr_phi0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','impr_phi0')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%nonbonded)) THEN
    DEALLOCATE(ff_type%inp_info%nonbonded,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','nonbonded')
  END IF
  IF(ASSOCIATED(ff_type%inp_info%nonbonded14)) THEN
    DEALLOCATE(ff_type%inp_info%nonbonded14,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','nonbonded14')
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. DEALLOCATE things in ff_type%chm_info
  !-----------------------------------------------------------------------------
  IF(ASSOCIATED(ff_type%chm_info%charge_atm)) THEN
    DEALLOCATE(ff_type%chm_info%charge_atm,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','charge_atm')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%charge)) THEN
    DEALLOCATE(ff_type%chm_info%charge,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','charge')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%bond_a)) THEN
    DEALLOCATE(ff_type%chm_info%bond_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bond_a')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%bond_b)) THEN
    DEALLOCATE(ff_type%chm_info%bond_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bond_b')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%bond_k)) THEN
    DEALLOCATE(ff_type%chm_info%bond_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bond_k')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%bond_r0)) THEN
    DEALLOCATE(ff_type%chm_info%bond_r0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bond_r0')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%bend_a)) THEN
    DEALLOCATE(ff_type%chm_info%bend_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bend_a')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%bend_b)) THEN
    DEALLOCATE(ff_type%chm_info%bend_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bend_b')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%bend_c)) THEN
    DEALLOCATE(ff_type%chm_info%bend_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bend_c')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%bend_k)) THEN
    DEALLOCATE(ff_type%chm_info%bend_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bend_k')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%bend_theta0)) THEN
    DEALLOCATE(ff_type%chm_info%bend_theta0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','bend_theta0')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%ub_a)) THEN
    DEALLOCATE(ff_type%chm_info%ub_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','ub_a')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%ub_b)) THEN
    DEALLOCATE(ff_type%chm_info%ub_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','ub_b')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%ub_c)) THEN
    DEALLOCATE(ff_type%chm_info%ub_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','ub_c')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%ub_k)) THEN
    DEALLOCATE(ff_type%chm_info%ub_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','ub_k')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%ub_r0)) THEN
    DEALLOCATE(ff_type%chm_info%ub_r0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','ub_r0')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%torsion_a)) THEN
    DEALLOCATE(ff_type%chm_info%torsion_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_a')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%torsion_b)) THEN
    DEALLOCATE(ff_type%chm_info%torsion_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_b')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%torsion_c)) THEN
    DEALLOCATE(ff_type%chm_info%torsion_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_c')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%torsion_d)) THEN
    DEALLOCATE(ff_type%chm_info%torsion_d,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_d')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%torsion_k)) THEN
    DEALLOCATE(ff_type%chm_info%torsion_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_k')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%torsion_m)) THEN
    DEALLOCATE(ff_type%chm_info%torsion_m,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_m')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%torsion_cosphi0)) THEN
    DEALLOCATE(ff_type%chm_info%torsion_cosphi0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','torsion_cosphi0')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%impr_a)) THEN
    DEALLOCATE(ff_type%chm_info%impr_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','impr_a')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%impr_b)) THEN
    DEALLOCATE(ff_type%chm_info%impr_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','impr_b')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%impr_c)) THEN
    DEALLOCATE(ff_type%chm_info%impr_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','impr_c')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%impr_d)) THEN
    DEALLOCATE(ff_type%chm_info%impr_d,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','impr_d')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%impr_k)) THEN
    DEALLOCATE(ff_type%chm_info%impr_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','impr_k')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%impr_phi0)) THEN
    DEALLOCATE(ff_type%chm_info%impr_phi0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','impr_phi0')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%nonbond_a)) THEN
    DEALLOCATE(ff_type%chm_info%nonbond_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','nonbond_a')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%nonbond_eps)) THEN
    DEALLOCATE(ff_type%chm_info%nonbond_eps,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','nonbond_eps')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%nonbond_rmin2)) THEN
    DEALLOCATE(ff_type%chm_info%nonbond_rmin2,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','nonbond_rmin2')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%nonbond_eps_14)) THEN
    DEALLOCATE(ff_type%chm_info%nonbond_eps_14,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','nonbond_eps_14')
  END IF
  IF(ASSOCIATED(ff_type%chm_info%nonbond_rmin2_14)) THEN
    DEALLOCATE(ff_type%chm_info%nonbond_rmin2_14,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','nonbond_rmin2_14')
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. DEALLOCATE things in ff_type%amb_info
  !-----------------------------------------------------------------------------



  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. DEALLOCATE things in ff_type
  !-----------------------------------------------------------------------------
  IF(ASSOCIATED(ff_type%inp_info)) THEN
    DEALLOCATE(ff_type%inp_info,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','inp_info')
  END IF
  IF(ASSOCIATED(ff_type%chm_info)) THEN
    DEALLOCATE(ff_type%chm_info,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','chm_info')
  END IF
  IF(ASSOCIATED(ff_type%amb_info)) THEN
    DEALLOCATE(ff_type%amb_info,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_ff_type','amb_info')
  END IF



!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE deallocate_ff_type




END MODULE force_field_types
