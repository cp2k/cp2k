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
!!   AUTHOR
!!     MK (29.08.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE force_field_types

  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE pair_potential_types,            ONLY: pair_potential_p_release,&
                                             pair_potential_p_type
  USE splines,                         ONLY: spline_data_release,&
                                             spline_data_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_field_types'
  

! *** Define the derived structure types ***

  ! k(3) ... allow for IPBV bond types sum(i=2,4) k_i (r-r0)**i/i
  TYPE bond_kind_type
    CHARACTER(LEN=default_string_length) :: name
    REAL(KIND = dp)                :: k(3),r0   
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
    CHARACTER(LEN=default_string_length)            :: name
    TYPE ( spline_data_type ), POINTER             :: spline_data
    INTEGER                                         :: kind_number
  END TYPE onfo_kind_type

  TYPE input_info_type
    CHARACTER (LEN=default_string_length), POINTER :: charge_atm(:)
    REAL(KIND=dp), POINTER                         :: charge(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_kind(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_b(:)
    REAL(KIND=dp), POINTER                         :: bond_k(:,:)
    REAL(KIND=dp), POINTER                         :: bond_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_kind(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_c(:)
    REAL(KIND=dp), POINTER                         :: bend_k(:)
    REAL(KIND=dp), POINTER                         :: bend_theta0(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_kind(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_c(:)
    REAL(KIND=dp), POINTER                         :: ub_k(:)
    REAL(KIND=dp), POINTER                         :: ub_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_kind(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_d(:)
    REAL(KIND=dp), POINTER                         :: torsion_k(:)
    REAL(KIND=dp), POINTER                         :: torsion_m(:)
    REAL(KIND=dp), POINTER                         :: torsion_cosphi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_kind(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_d(:)
    REAL(KIND=dp), POINTER                         :: impr_k(:)
    REAL(KIND=dp), POINTER                         :: impr_phi0(:)
    TYPE (pair_potential_p_type), POINTER          :: nonbonded
    TYPE (pair_potential_p_type), POINTER          :: nonbonded14
  END TYPE input_info_type

  TYPE charmm_info_type
    CHARACTER (LEN=default_string_length), POINTER :: charge_atm(:)
    REAL(KIND=dp), POINTER                         :: charge(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bond_b(:)
    REAL(KIND=dp), POINTER                         :: bond_k(:)
    REAL(KIND=dp), POINTER                         :: bond_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: bend_c(:)
    REAL(KIND=dp), POINTER                         :: bend_k(:)
    REAL(KIND=dp), POINTER                         :: bend_theta0(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: ub_c(:)
    REAL(KIND=dp), POINTER                         :: ub_k(:)
    REAL(KIND=dp), POINTER                         :: ub_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: torsion_d(:)
    REAL(KIND=dp), POINTER                         :: torsion_k(:)
    INTEGER, POINTER                               :: torsion_m(:)
    REAL(KIND=dp), POINTER                         :: torsion_cosphi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: impr_d(:)
    REAL(KIND=dp), POINTER                         :: impr_k(:)
    REAL(KIND=dp), POINTER                         :: impr_phi0(:)
    CHARACTER (LEN=default_string_length), POINTER :: nonbond_a(:)
    REAL (KIND=dp), POINTER                        :: nonbond_eps(:)
    REAL (KIND=dp), POINTER                        :: nonbond_rmin2(:)
    CHARACTER (LEN=default_string_length), POINTER :: nonbond_a_14(:)
    REAL (KIND=dp), POINTER                        :: nonbond_eps_14(:)
    REAL (KIND=dp), POINTER                        :: nonbond_rmin2_14(:)
  END TYPE charmm_info_type

  TYPE force_field_type
    TYPE (input_info_type),POINTER   :: inp_info
    TYPE (charmm_info_type),POINTER  :: chm_info
    CHARACTER (LEN=default_string_length)     :: ff_file_name
    CHARACTER (LEN=default_string_length)     :: ff_type
    REAL(KIND=dp)                             :: rcut_nb
    REAL(KIND=dp)                             :: ei_scale14
    REAL(KIND=dp)                             :: vdw_scale14
    REAL(KIND=dp)                             :: eps_spline,&
                                                 max_energy,&
                                                 emax_spline,&
                                                 rlow_nb
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
            impr_kind_dealloc_ref,&
            init_chm_info,deallocate_chm_info,&
            init_inp_info,deallocate_inp_info

! *** Public data types ***

  PUBLIC :: bend_kind_type,&
            bond_kind_type,&
            impr_kind_type,&
            onfo_kind_type,&
            onfo_kind_dealloc_ref,&
            torsion_kind_type,&
            ub_kind_type,&
            ub_kind_dealloc_ref,&
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
  SUBROUTINE allocate_bend_kind_set(bend_kind_set,nkind,error)

    TYPE(bend_kind_type), DIMENSION(:), &
      POINTER                                :: bend_kind_set
    INTEGER, INTENT(IN)                      :: nkind
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_bend_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       NULLIFY(bend_kind_set)
       ALLOCATE (bend_kind_set(nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO ikind=1,nkind
          bend_kind_set(ikind)%name = "UNSET"
          bend_kind_set(ikind)%k = 0.0_dp
          bend_kind_set(ikind)%theta0 = 0.0_dp
          bend_kind_set(ikind)%kind_number = ikind
       END DO
    END IF
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
  SUBROUTINE allocate_bond_kind_set(bond_kind_set,nkind,error)

    TYPE(bond_kind_type), DIMENSION(:), &
      POINTER                                :: bond_kind_set
    INTEGER, INTENT(IN)                      :: nkind
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_bond_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    NULLIFY(bond_kind_set)
    IF (.NOT.failure) THEN
       ALLOCATE (bond_kind_set(nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO ikind=1,nkind
          bond_kind_set(ikind)%name = "UNSET"
          bond_kind_set(ikind)%k(:) = 0.0_dp
          bond_kind_set(ikind)%r0 = 0.0_dp
          bond_kind_set(ikind)%kind_number = ikind
       END DO
    END IF
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
  SUBROUTINE allocate_torsion_kind_set(torsion_kind_set,nkind,error)

    TYPE(torsion_kind_type), DIMENSION(:), &
      POINTER                                :: torsion_kind_set
    INTEGER, INTENT(IN)                      :: nkind
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_torsion_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    NULLIFY(torsion_kind_set)
    IF (.NOT.failure) THEN
       ALLOCATE (torsion_kind_set(nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       DO ikind=1,nkind
          torsion_kind_set(ikind)%name = "UNSET"
          torsion_kind_set(ikind)%nmul = 0
          torsion_kind_set(ikind)%a(:) = 0.0_dp
          NULLIFY(torsion_kind_set(ikind)%k)
          NULLIFY(torsion_kind_set(ikind)%m)
          NULLIFY(torsion_kind_set(ikind)%cosphi0)
          torsion_kind_set(ikind)%kind_number = ikind
       END DO
    END IF
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
  SUBROUTINE allocate_ub_kind_set(ub_kind_set,nkind,error)

    TYPE(ub_kind_type), DIMENSION(:), &
      POINTER                                :: ub_kind_set
    INTEGER, INTENT(IN)                      :: nkind
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_ub_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    NULLIFY(ub_kind_set)
    ALLOCATE (ub_kind_set(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
  SUBROUTINE allocate_impr_kind_set(impr_kind_set,nkind,error)

    TYPE(impr_kind_type), DIMENSION(:), &
      POINTER                                :: impr_kind_set
    INTEGER, INTENT(IN)                      :: nkind
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_impr_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    NULLIFY(impr_kind_set)
    ALLOCATE (impr_kind_set(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      impr_kind_set(ikind)%name = "UNSET"
      impr_kind_set(ikind)%k = 0.0_dp
      impr_kind_set(ikind)%phi0 = 0.0_dp
      impr_kind_set(ikind)%kind_number = ikind
    END DO
  END SUBROUTINE allocate_impr_kind_set

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
  SUBROUTINE allocate_onfo_kind_set(onfo_kind_set,nkind,error)

    TYPE(onfo_kind_type), DIMENSION(:), &
      POINTER                                :: onfo_kind_set
    INTEGER, INTENT(IN)                      :: nkind
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_onfo_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    NULLIFY(onfo_kind_set)
    ALLOCATE (onfo_kind_set(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      onfo_kind_set(ikind)%name = "UNSET"
      onfo_kind_set(ikind)%kind_number = ikind
      NULLIFY(onfo_kind_set(ikind)%spline_data)
    END DO
  END SUBROUTINE allocate_onfo_kind_set

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
  SUBROUTINE init_ff_type (ff_type,error)
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_ff_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Initialize 
    !-----------------------------------------------------------------------------
    NULLIFY(ff_type%inp_info,ff_type%chm_info)
    ALLOCATE(ff_type%inp_info,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ff_type%chm_info,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Initialize hardcoded ff_type variables
    !-----------------------------------------------------------------------------
    ff_type%rlow_nb=0.9_dp
    
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Initialize and Nullify things in ff_type%inp_info
    !-----------------------------------------------------------------------------
    CALL init_inp_info(ff_type%inp_info)
    
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. Initialize and Nullify things in ff_type%chm_info
    !-----------------------------------------------------------------------------
    CALL init_chm_info(ff_type%chm_info)
    
  END SUBROUTINE init_ff_type

!******************************************************************************
!!****** force_field_types/init_chm_info [1.0] *
!!
!!   NAME
!!     init_chm_info
!!
!!   FUNCTION
!!     Specific function to initialize the chm_info
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE init_chm_info(chm_info)
    TYPE(charmm_info_type), POINTER          :: chm_info

  NULLIFY(chm_info%charge_atm)
  NULLIFY(chm_info%charge)
  NULLIFY(chm_info%bond_a)
  NULLIFY(chm_info%bond_b)
  NULLIFY(chm_info%bond_k)
  NULLIFY(chm_info%bond_r0)
  NULLIFY(chm_info%bend_a)
  NULLIFY(chm_info%bend_b)
  NULLIFY(chm_info%bend_c)
  NULLIFY(chm_info%bend_k)
  NULLIFY(chm_info%bend_theta0)
  NULLIFY(chm_info%ub_a)
  NULLIFY(chm_info%ub_b)
  NULLIFY(chm_info%ub_c)
  NULLIFY(chm_info%ub_k)
  NULLIFY(chm_info%ub_r0)
  NULLIFY(chm_info%torsion_a)
  NULLIFY(chm_info%torsion_b)
  NULLIFY(chm_info%torsion_c)
  NULLIFY(chm_info%torsion_d)
  NULLIFY(chm_info%torsion_k)
  NULLIFY(chm_info%torsion_m)
  NULLIFY(chm_info%torsion_cosphi0)
  NULLIFY(chm_info%impr_a)
  NULLIFY(chm_info%impr_b)
  NULLIFY(chm_info%impr_c)
  NULLIFY(chm_info%impr_d)
  NULLIFY(chm_info%impr_k)
  NULLIFY(chm_info%impr_phi0)
  NULLIFY(chm_info%nonbond_a)
  NULLIFY(chm_info%nonbond_eps)
  NULLIFY(chm_info%nonbond_rmin2)
  NULLIFY(chm_info%nonbond_a_14)
  NULLIFY(chm_info%nonbond_eps_14)
  NULLIFY(chm_info%nonbond_rmin2_14)
  
END SUBROUTINE init_chm_info

!******************************************************************************
!!****** force_field_types/init_inp_info [1.0] *
!!
!!   NAME
!!     init_inp_info
!!
!!   FUNCTION
!!     Specific function to initialize the inp_info
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE init_inp_info(inp_info)
    TYPE(input_info_type), POINTER           :: inp_info

  NULLIFY(inp_info%charge_atm)
  NULLIFY(inp_info%charge)
  NULLIFY(inp_info%bond_kind)
  NULLIFY(inp_info%bond_a)
  NULLIFY(inp_info%bond_b)
  NULLIFY(inp_info%bond_k)
  NULLIFY(inp_info%bond_r0)
  NULLIFY(inp_info%bend_kind)
  NULLIFY(inp_info%bend_a)
  NULLIFY(inp_info%bend_b)
  NULLIFY(inp_info%bend_c)
  NULLIFY(inp_info%bend_k)
  NULLIFY(inp_info%bend_theta0)
  NULLIFY(inp_info%ub_kind)
  NULLIFY(inp_info%ub_a)
  NULLIFY(inp_info%ub_b)
  NULLIFY(inp_info%ub_c)
  NULLIFY(inp_info%ub_k)
  NULLIFY(inp_info%ub_r0)
  NULLIFY(inp_info%torsion_kind)
  NULLIFY(inp_info%torsion_a)
  NULLIFY(inp_info%torsion_b)
  NULLIFY(inp_info%torsion_c)
  NULLIFY(inp_info%torsion_d)
  NULLIFY(inp_info%torsion_k)
  NULLIFY(inp_info%torsion_m)
  NULLIFY(inp_info%torsion_cosphi0)
  NULLIFY(inp_info%impr_kind)
  NULLIFY(inp_info%impr_a)
  NULLIFY(inp_info%impr_b)
  NULLIFY(inp_info%impr_c)
  NULLIFY(inp_info%impr_d)
  NULLIFY(inp_info%impr_k)
  NULLIFY(inp_info%impr_phi0)
  NULLIFY(inp_info%nonbonded)
  NULLIFY(inp_info%nonbonded14)

END SUBROUTINE init_inp_info


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
SUBROUTINE deallocate_ff_type (ff_type,error)
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_ff_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. DEALLOCATE things in ff_type%inp_info
  !-----------------------------------------------------------------------------
    CALL deallocate_inp_info(ff_type%inp_info, error=error)

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. DEALLOCATE things in ff_type%chm_info
  !-----------------------------------------------------------------------------
    CALL deallocate_chm_info(ff_type%chm_info,error)

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. DEALLOCATE things in ff_type
  !-----------------------------------------------------------------------------
    IF(ASSOCIATED(ff_type%inp_info)) THEN
       DEALLOCATE(ff_type%inp_info,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(ff_type%chm_info)) THEN
       DEALLOCATE(ff_type%chm_info,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
    END IF

END SUBROUTINE deallocate_ff_type

!******************************************************************************
!!****** force_field_types/deallocate_chm_info [1.0] *
!!
!!   NAME
!!     deallocate_chm_info
!!
!!   FUNCTION
!!     Specific function to deallocate the chm_info
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE deallocate_chm_info(chm_info,error)
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_chm_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure = .FALSE.
  IF(ASSOCIATED(chm_info%charge_atm)) THEN
    DEALLOCATE(chm_info%charge_atm,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%charge)) THEN
    DEALLOCATE(chm_info%charge,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%bond_a)) THEN
    DEALLOCATE(chm_info%bond_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%bond_b)) THEN
    DEALLOCATE(chm_info%bond_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%bond_k)) THEN
    DEALLOCATE(chm_info%bond_k,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%bond_r0)) THEN
    DEALLOCATE(chm_info%bond_r0,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%bend_a)) THEN
    DEALLOCATE(chm_info%bend_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%bend_b)) THEN
    DEALLOCATE(chm_info%bend_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%bend_c)) THEN
    DEALLOCATE(chm_info%bend_c,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%bend_k)) THEN
    DEALLOCATE(chm_info%bend_k,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%bend_theta0)) THEN
    DEALLOCATE(chm_info%bend_theta0,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%ub_a)) THEN
    DEALLOCATE(chm_info%ub_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%ub_b)) THEN
    DEALLOCATE(chm_info%ub_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%ub_c)) THEN
    DEALLOCATE(chm_info%ub_c,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%ub_k)) THEN
    DEALLOCATE(chm_info%ub_k,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%ub_r0)) THEN
    DEALLOCATE(chm_info%ub_r0,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%torsion_a)) THEN
    DEALLOCATE(chm_info%torsion_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%torsion_b)) THEN
    DEALLOCATE(chm_info%torsion_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%torsion_c)) THEN
    DEALLOCATE(chm_info%torsion_c,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%torsion_d)) THEN
    DEALLOCATE(chm_info%torsion_d,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%torsion_k)) THEN
    DEALLOCATE(chm_info%torsion_k,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%torsion_m)) THEN
    DEALLOCATE(chm_info%torsion_m,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%torsion_cosphi0)) THEN
    DEALLOCATE(chm_info%torsion_cosphi0,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%impr_a)) THEN
    DEALLOCATE(chm_info%impr_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%impr_b)) THEN
    DEALLOCATE(chm_info%impr_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%impr_c)) THEN
    DEALLOCATE(chm_info%impr_c,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%impr_d)) THEN
    DEALLOCATE(chm_info%impr_d,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%impr_k)) THEN
    DEALLOCATE(chm_info%impr_k,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%impr_phi0)) THEN
    DEALLOCATE(chm_info%impr_phi0,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%nonbond_a)) THEN
    DEALLOCATE(chm_info%nonbond_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%nonbond_eps)) THEN
    DEALLOCATE(chm_info%nonbond_eps,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%nonbond_rmin2)) THEN
    DEALLOCATE(chm_info%nonbond_rmin2,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%nonbond_a_14)) THEN
    DEALLOCATE(chm_info%nonbond_a_14,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%nonbond_eps_14)) THEN
    DEALLOCATE(chm_info%nonbond_eps_14,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(chm_info%nonbond_rmin2_14)) THEN
    DEALLOCATE(chm_info%nonbond_rmin2_14,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF

END SUBROUTINE deallocate_chm_info

!******************************************************************************
!!****** force_field_types/deallocate_inp_info [1.0] *
!!
!!   NAME
!!     deallocate_inp_info
!!
!!   FUNCTION
!!     Specific function to deallocate the inp_info
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE deallocate_inp_info(inp_info, error)
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_inp_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure = .FALSE.
  IF(ASSOCIATED(inp_info%charge_atm)) THEN
    DEALLOCATE(inp_info%charge_atm,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%charge)) THEN
    DEALLOCATE(inp_info%charge,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%bond_kind)) THEN
    DEALLOCATE(inp_info%bond_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%bond_a)) THEN
    DEALLOCATE(inp_info%bond_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%bond_b)) THEN
    DEALLOCATE(inp_info%bond_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%bond_k)) THEN
    DEALLOCATE(inp_info%bond_k,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%bond_r0)) THEN
    DEALLOCATE(inp_info%bond_r0,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%bend_kind)) THEN
    DEALLOCATE(inp_info%bend_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%bend_a)) THEN
    DEALLOCATE(inp_info%bend_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%bend_b)) THEN
    DEALLOCATE(inp_info%bend_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%bend_c)) THEN
    DEALLOCATE(inp_info%bend_c,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%bend_k)) THEN
    DEALLOCATE(inp_info%bend_k,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%bend_theta0)) THEN
    DEALLOCATE(inp_info%bend_theta0,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%ub_kind)) THEN
    DEALLOCATE(inp_info%ub_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%ub_a)) THEN
    DEALLOCATE(inp_info%ub_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%ub_b)) THEN
    DEALLOCATE(inp_info%ub_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%ub_c)) THEN
    DEALLOCATE(inp_info%ub_c,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%ub_k)) THEN
    DEALLOCATE(inp_info%ub_k,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%ub_r0)) THEN
    DEALLOCATE(inp_info%ub_r0,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%torsion_kind)) THEN
    DEALLOCATE(inp_info%torsion_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%torsion_a)) THEN
    DEALLOCATE(inp_info%torsion_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%torsion_b)) THEN
    DEALLOCATE(inp_info%torsion_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%torsion_c)) THEN
    DEALLOCATE(inp_info%torsion_c,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%torsion_d)) THEN
    DEALLOCATE(inp_info%torsion_d,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%torsion_k)) THEN
    DEALLOCATE(inp_info%torsion_k,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%torsion_m)) THEN
    DEALLOCATE(inp_info%torsion_m,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%torsion_cosphi0)) THEN
    DEALLOCATE(inp_info%torsion_cosphi0,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%impr_kind)) THEN
    DEALLOCATE(inp_info%impr_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%impr_a)) THEN
    DEALLOCATE(inp_info%impr_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%impr_b)) THEN
    DEALLOCATE(inp_info%impr_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%impr_c)) THEN
    DEALLOCATE(inp_info%impr_c,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%impr_d)) THEN
    DEALLOCATE(inp_info%impr_d,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%impr_k)) THEN
    DEALLOCATE(inp_info%impr_k,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%impr_phi0)) THEN
    DEALLOCATE(inp_info%impr_phi0,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF(ASSOCIATED(inp_info%nonbonded)) THEN
     CALL pair_potential_p_release(inp_info%nonbonded,error=error)
  END IF
  IF(ASSOCIATED(inp_info%nonbonded14)) THEN
     CALL pair_potential_p_release(inp_info%nonbonded14,error=error)
  END IF

END SUBROUTINE deallocate_inp_info

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
  SUBROUTINE deallocate_bend_kind_set(bend_kind_set,error)

    TYPE(bend_kind_type), DIMENSION(:), &
      POINTER                                :: bend_kind_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_bend_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    DEALLOCATE (bend_kind_set,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
  SUBROUTINE deallocate_bond_kind_set(bond_kind_set,error)

    TYPE(bond_kind_type), DIMENSION(:), &
      POINTER                                :: bond_kind_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_bond_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    DEALLOCATE (bond_kind_set,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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

    TYPE(torsion_kind_type), INTENT(INOUT)   :: torsion_kind
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'torsion_kind_dealloc_ref', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF(ASSOCIATED(torsion_kind%k)) THEN
       DEALLOCATE(torsion_kind%k,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(torsion_kind%m)) THEN
       DEALLOCATE(torsion_kind%m,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(ASSOCIATED(torsion_kind%cosphi0)) THEN
       DEALLOCATE(torsion_kind%cosphi0,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    
  END SUBROUTINE torsion_kind_dealloc_ref

! *****************************************************************************

!******************************************************************************
!!****** force_field_types/ub_kind_dealloc_ref [1.0] *
!!
!!   NAME
!!     ub_kind_dealloc_ref
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
  SUBROUTINE ub_kind_dealloc_ref(ub_kind_set, error)
    TYPE(ub_kind_type), DIMENSION(:), &
      POINTER                                :: ub_kind_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ub_kind_dealloc_ref', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    DEALLOCATE(ub_kind_set, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE ub_kind_dealloc_ref
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

    TYPE(impr_kind_type), INTENT(INOUT)      :: impr_kind
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'impr_kind_dealloc_ref', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

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
  SUBROUTINE  onfo_kind_dealloc_ref(onfo_kind_set, error)
    TYPE(onfo_kind_type), DIMENSION(:), &
      POINTER                                :: onfo_kind_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'onfo_kind_dealloc_ref', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ionfo, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    DO ionfo=1,SIZE(onfo_kind_set)
      CALL spline_data_release(onfo_kind_set(ionfo)%spline_data, error)
    END DO
    DEALLOCATE(onfo_kind_set,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE onfo_kind_dealloc_ref
! *****************************************************************************

END MODULE force_field_types
