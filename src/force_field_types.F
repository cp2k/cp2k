!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp,&
                                             default_string_length,&
                                             int_size
  USE pair_potential_types,            ONLY: pair_potential_type
  USE splines,                         ONLY: spline_environment_type,&
                                             spline_env_release
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = "force_field_types"
  

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
    INTEGER                        :: nmul
    INTEGER, POINTER               :: m(:)
    REAL(KIND = dp), POINTER       :: k(:),cosphi0(:)
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
    REAL(KIND=dp), POINTER                            :: charge(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_b(:)
    REAL(KIND=dp), POINTER                            :: bond_k(:)
    REAL(KIND=dp), POINTER                            :: bond_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_c(:)
    REAL(KIND=dp), POINTER                            :: bend_k(:)
    REAL(KIND=dp), POINTER                            :: bend_theta0(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_c(:)
    REAL(KIND=dp), POINTER                            :: ub_k(:)
    REAL(KIND=dp), POINTER                            :: ub_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_d(:)
    REAL(KIND=dp), POINTER                            :: torsion_k(:)
    REAL(KIND=dp), POINTER                            :: torsion_m(:)
    REAL(KIND=dp), POINTER                            :: torsion_cosphi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_d(:)
    REAL(KIND=dp), POINTER                            :: impr_k(:)
    REAL(KIND=dp), POINTER                            :: impr_phi0(:)
    TYPE (pair_potential_type), POINTER            :: nonbonded(:)
    TYPE (pair_potential_type), POINTER            :: nonbonded14(:)
  END TYPE input_info_type

  TYPE charmm_info_type
    CHARACTER (LEN=default_string_length), POINTER :: charge_atm(:)
    REAL(KIND=dp), POINTER                            :: charge(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_b(:)
    REAL(KIND=dp), POINTER                            :: bond_k(:)
    REAL(KIND=dp), POINTER                            :: bond_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_c(:)
    REAL(KIND=dp), POINTER                            :: bend_k(:)
    REAL(KIND=dp), POINTER                            :: bend_theta0(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_c(:)
    REAL(KIND=dp), POINTER                            :: ub_k(:)
    REAL(KIND=dp), POINTER                            :: ub_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_d(:)
    REAL(KIND=dp), POINTER                            :: torsion_k(:)
    INTEGER, POINTER                               :: torsion_m(:)
    REAL(KIND=dp), POINTER                            :: torsion_cosphi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_d(:)
    REAL(KIND=dp), POINTER                            :: impr_k(:)
    REAL(KIND=dp), POINTER                            :: impr_phi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: nonbond_a(:)
    REAL(KIND=dp), POINTER                            :: nonbond_eps(:)
    REAL(KIND=dp), POINTER                            :: nonbond_rmin2(:)
    REAL(KIND=dp), POINTER                            :: nonbond_eps_14(:)
    REAL(KIND=dp), POINTER                            :: nonbond_rmin2_14(:)
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
    REAL(KIND=dp)                                :: rcut_nb
  END TYPE force_field_type

! *** Public subroutines ***

  PUBLIC :: allocate_bend_kind_set,&
            allocate_bond_kind_set,&
            allocate_ub_kind_set,&
            allocate_torsion_kind_set,&
            allocate_impr_kind_set,&
            allocate_onfo_kind_set,&
            init_ff_type,&
            deallocate_ff_type,&
            deallocate_bend_kind_set,&
            deallocate_bond_kind_set,&
            torsion_kind_dealloc_ref,&
            impr_kind_dealloc_ref

! *** Public data types ***

  PUBLIC :: bend_kind_type,&
            bond_kind_type,&
            impr_kind_type,&
            onfo_kind_type,&
            onfo_kind_dealloc_ref,&
            torsion_kind_type,&
            ub_kind_type,&
            deallocate_ub_kind_set,&
            amber_info_type,&
            charmm_info_type,&
            input_info_type,&
            force_field_type

! *****************************************************************************

CONTAINS

! *****************************************************************************

!******************************************************************************
!!****** force_field_types/allocate_bend_kind_set [1.0] *
!!
!!   NAME
!!     allocate_bend_kind_set
!!
!!   FUNCTION
!!     Allocate and initialize a bend kind set.
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE allocate_bend_kind_set(bend_kind_set,nkind)

    TYPE(bend_kind_type), DIMENSION(:), &
      POINTER                                :: bend_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_bend_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    ALLOCATE (bend_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
                                     "bend_kind_set",nkind*int_size)

    DO ikind=1,nkind
      bend_kind_set(ikind)%name = "UNSET"
      bend_kind_set(ikind)%k = 0.0_dp
      bend_kind_set(ikind)%theta0 = 0.0_dp
      bend_kind_set(ikind)%kind_number = ikind
    END DO

  END SUBROUTINE allocate_bend_kind_set

! *****************************************************************************

!******************************************************************************
!!****** force_field_types/allocate_bond_kind_set [1.0] *
!!
!!   NAME
!!     allocate_bond_kind_set
!!
!!   FUNCTION
!!     Allocate and initialize a bond kind set.
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE allocate_bond_kind_set(bond_kind_set,nkind)

    TYPE(bond_kind_type), DIMENSION(:), &
      POINTER                                :: bond_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_bond_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    ALLOCATE (bond_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
                                     "bond_kind_set",nkind*int_size)

    DO ikind=1,nkind
      bond_kind_set(ikind)%name = "UNSET"
      bond_kind_set(ikind)%k = 0.0_dp
      bond_kind_set(ikind)%r0 = 0.0_dp
      bond_kind_set(ikind)%kind_number = ikind
    END DO

  END SUBROUTINE allocate_bond_kind_set

! *****************************************************************************

!******************************************************************************
!!****** force_field_types/allocate_torsion_kind_set [1.0] *
!!
!!   NAME
!!     allocate_torsion_kind_set
!!
!!   FUNCTION
!!     Allocate and initialize a torsion kind set.
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE allocate_torsion_kind_set(torsion_kind_set,nkind)

    TYPE(torsion_kind_type), DIMENSION(:), &
      POINTER                                :: torsion_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_torsion_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    ALLOCATE (torsion_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
                                     "torsion_kind_set",nkind*int_size)

    DO ikind=1,nkind
      torsion_kind_set(ikind)%name = "UNSET"
      torsion_kind_set(ikind)%nmul = 0
      NULLIFY(torsion_kind_set(ikind)%k)
      NULLIFY(torsion_kind_set(ikind)%m)
      NULLIFY(torsion_kind_set(ikind)%cosphi0)
      torsion_kind_set(ikind)%kind_number = ikind
    END DO

  END SUBROUTINE allocate_torsion_kind_set

! *****************************************************************************

!******************************************************************************
!!****** force_field_types/allocate_ub_kind_set [1.0] *
!!
!!   NAME
!!     allocate_ub_kind_set
!!
!!   FUNCTION
!!     Allocate and initialize a ub kind set.
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE allocate_ub_kind_set(ub_kind_set,nkind)

    TYPE(ub_kind_type), DIMENSION(:), &
      POINTER                                :: ub_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_ub_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    ALLOCATE (ub_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
                                     "ub_kind_set",nkind*int_size)

    DO ikind=1,nkind
      ub_kind_set(ikind)%name = "UNSET"
      ub_kind_set(ikind)%k = 0.0_dp
      ub_kind_set(ikind)%r0 = 0.0_dp
      ub_kind_set(ikind)%kind_number = ikind
    END DO

  END SUBROUTINE allocate_ub_kind_set


!******************************************************************************
!!****** force_field_types/allocate_impr_kind_set [1.0] *
!!
!!   NAME
!!     allocate_impr_kind_set
!!
!!   FUNCTION
!!     Allocate and initialize a impr kind set.
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE allocate_impr_kind_set(impr_kind_set,nkind)

    TYPE(impr_kind_type), DIMENSION(:), &
      POINTER                                :: impr_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_impr_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(impr_kind_set)) THEN
       DEALLOCATE (impr_kind_set,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
            "impr_kind_set",nkind*int_size)
    END IF
    ALLOCATE (impr_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
                                     "impr_kind_set",nkind*int_size)

    DO ikind=1,nkind
      impr_kind_set(ikind)%name = "UNSET"
      impr_kind_set(ikind)%k = 0.0_dp
      impr_kind_set(ikind)%phi0 = 0.0_dp
      impr_kind_set(ikind)%kind_number = ikind
    END DO

  END SUBROUTINE allocate_impr_kind_set

! *****************************************************************************



!!****** force_field_types/allocate_onfo_kind_set [1.0] *
!!
!!   NAME
!!     allocate_onfo_kind_set
!!
!!   FUNCTION
!!     Allocate and initialize a onfo kind set.
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE allocate_onfo_kind_set(onfo_kind_set,nkind)

    TYPE(onfo_kind_type), DIMENSION(:), &
      POINTER                                :: onfo_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_onfo_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    ALLOCATE (onfo_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
                                     "onfo_kind_set",nkind*int_size)

    DO ikind=1,nkind
      onfo_kind_set(ikind)%name = ""
      onfo_kind_set(ikind)%kind_number = ikind
      NULLIFY(onfo_kind_set(ikind)%spline_env)
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
  ff_type%rcut_nb=-1.0_dp


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

! *****************************************************************************

!******************************************************************************
!!****** force_field_types/deallocate_bend_kind_set [1.0] *
!!
!!   NAME
!!     deallocate_bend_kind_set
!!
!!   FUNCTION
!!     Deallocate a bend kind set.
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE deallocate_bend_kind_set(bend_kind_set)

    TYPE(bend_kind_type), DIMENSION(:), &
      POINTER                                :: bend_kind_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_bend_kind_set"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    DEALLOCATE (bend_kind_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
                                     "bend_kind_set")

  END SUBROUTINE deallocate_bend_kind_set

! *****************************************************************************

!******************************************************************************
!!****** force_field_types/deallocate_bond_kind_set [1.0] *
!!
!!   NAME
!!     deallocate_bond_kind_set
!!
!!   FUNCTION
!!     Deallocate a bond kind set.
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE deallocate_bond_kind_set(bond_kind_set)

    TYPE(bond_kind_type), DIMENSION(:), &
      POINTER                                :: bond_kind_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_bond_kind_set"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    DEALLOCATE (bond_kind_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
                                     "bond_kind_set")

  END SUBROUTINE deallocate_bond_kind_set

! *****************************************************************************

!******************************************************************************
!!****** force_field_types/torsion_kind_dealloc_ref [1.0] *
!!
!!   NAME
!!     torsion_kind_dealloc_ref
!!
!!   FUNCTION
!!     Deallocate a torsion kind element
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE torsion_kind_dealloc_ref(torsion_kind,error)

    TYPE(torsion_kind_type), INTENT(INOUT)       :: torsion_kind
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    integer :: istat
    CHARACTER(len=*), PARAMETER :: routineN='torsion_kind_dealloc_ref',&
         routineP=moduleN//':'//routineN
    logical :: failure
    
    failure = .false.
!   ---------------------------------------------------------------------------

    DEALLOCATE(torsion_kind%k,STAT=istat)
    IF(istat/=0) CALL stop_memory(routineN,"k")
    DEALLOCATE(torsion_kind%m,STAT=istat)
    IF(istat/=0) CALL stop_memory(routineN,"m")
    DEALLOCATE(torsion_kind%cosphi0,STAT=istat)
    IF(istat/=0) CALL stop_memory(routineN,"cosphi0")

  END SUBROUTINE torsion_kind_dealloc_ref

! *****************************************************************************

!******************************************************************************
!!****** force_field_types/deallocate_ub_kind_set [1.0] *
!!
!!   NAME
!!     deallocate_ub_kind_set
!!
!!   FUNCTION
!!     Deallocate a ub kind set.
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  subroutine deallocate_ub_kind_set(ub_kind_set, error)
    TYPE(ub_kind_type), DIMENSION(:), &
         POINTER                                :: ub_kind_set
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN='deallocate_ub_kind_set',&
         routineP=moduleN//':'//routineN
    INTEGER                                  ::  stat
    LOGICAL :: failure

    failure = .false.
    deallocate(ub_kind_set, STAT=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)

  end subroutine deallocate_ub_kind_set
! *****************************************************************************

!******************************************************************************
!!****** force_field_types/impr_kind_dealloc_ref [1.0] *
!!
!!   NAME
!!     impr_kind_dealloc_ref
!!
!!   FUNCTION
!!     Deallocate a impr kind element
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE impr_kind_dealloc_ref(impr_kind, error)

    TYPE(impr_kind_type), INTENT(INOUT) :: impr_kind
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    integer :: istat
    CHARACTER(len=*), PARAMETER :: routineN='impr_kind_dealloc_ref',&
         routineP=moduleN//':'//routineN
    logical :: failure
    
    failure = .false.

    !
    ! Questa e' la migliore routine che mente umana abbia concepito! ;-)
    ! 

  END SUBROUTINE impr_kind_dealloc_ref

! *****************************************************************************

!!****** force_field_types/onfo_kind_dealloc_ref [1.0] *
!!
!!   NAME
!!     deallocate_onfo_kind_set
!!
!!   FUNCTION
!!     Deallocate a onfo kind set.
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE  onfo_kind_dealloc_ref(onfo_kind, error)
    implicit none
    TYPE(onfo_kind_type), INTENT(INOUT):: onfo_kind
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    integer :: stat
    CHARACTER(len=*), PARAMETER :: routineN='onfo_kind_dealloc_ref',&
         routineP=moduleN//':'//routineN
    logical :: failure
    
    failure = .false.
    CALL spline_env_release(onfo_kind%spline_env, error)
  END SUBROUTINE onfo_kind_dealloc_ref
! *****************************************************************************

END MODULE force_field_types
