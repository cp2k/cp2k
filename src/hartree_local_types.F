!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/hartree_local_types [1.0] *
!!
!!   NAME
!!     hartree_local_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     
!!
!!   SOURCE
!******************************************************************************

MODULE hartree_local_types

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_size
  USE qs_rho_atom_types,               ONLY: rho_atom_coeff
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "hartree_local_types"

! *** Define the ecoul_1center_type ***

  TYPE ecoul_1center_type
    TYPE(rho_atom_coeff), POINTER   :: Vh1_h, Vh1_s
    TYPE(rho_atom_coeff), POINTER   :: aVh1b_h, aVh1b_s, aVh1b_0
    REAL(dp)                        :: ecoul_1_h, ecoul_1_s, &
                                       ecoul_1_z, ecoul_1_0
  END TYPE ecoul_1center_type

  TYPE hartree_local_type
    TYPE(ecoul_1center_type), &
            DIMENSION(:), POINTER   :: ecoul_1c
  END TYPE hartree_local_type

! *** Public subroutines ***

  PUBLIC :: allocate_ecoul_1center, get_hartree_local, &
            hartree_local_create,&
            hartree_local_release, set_hartree_local

! *** Public data types ***

  PUBLIC :: ecoul_1center_type, hartree_local_type

! *****************************************************************************

  CONTAINS

! *****************************************************************************

  SUBROUTINE  allocate_ecoul_1center(ecoul_1c,natom)

    TYPE(ecoul_1center_type), DIMENSION(:), &
                                 POINTER        :: ecoul_1c
    INTEGER, INTENT(IN)                         :: natom 

    CHARACTER(len=*), PARAMETER :: routine_name =  'allocate_ecoul_1center'

    INTEGER                                     :: iat, istat
!   ---------------------------------------------------------------------------

    IF(ASSOCIATED(ecoul_1c)) THEN
      CALL deallocate_ecoul_1center(ecoul_1c)
    END IF

    ALLOCATE(ecoul_1c(natom), STAT=istat) 
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                            "ecoul_1c",natom*int_size)

    DO iat = 1,natom
      ALLOCATE(ecoul_1c(iat)%Vh1_h)
      NULLIFY(ecoul_1c(iat)%Vh1_h%r_coef)
      ALLOCATE(ecoul_1c(iat)%Vh1_s)
      NULLIFY(ecoul_1c(iat)%Vh1_s%r_coef)
      ALLOCATE(ecoul_1c(iat)%aVh1b_h)
      NULLIFY(ecoul_1c(iat)%aVh1b_h%r_coef)
      ALLOCATE(ecoul_1c(iat)%aVh1b_s)
      NULLIFY(ecoul_1c(iat)%aVh1b_s%r_coef)
      ALLOCATE(ecoul_1c(iat)%aVh1b_0)
      NULLIFY(ecoul_1c(iat)%aVh1b_0%r_coef)
    END DO

  END SUBROUTINE allocate_ecoul_1center

! *****************************************************************************

  SUBROUTINE  deallocate_ecoul_1center(ecoul_1c)

    TYPE(ecoul_1center_type), DIMENSION(:), &
                                 POINTER        :: ecoul_1c

    CHARACTER(len=*), PARAMETER :: routine_name =  'deallocate_ecoul_1center'

    INTEGER                                     :: iat, istat, natom
!   ---------------------------------------------------------------------------

    natom = SIZE(ecoul_1c,1)

    DO iat= 1,natom
      DEALLOCATE(ecoul_1c(iat)%Vh1_h%r_coef,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "ecoul_1c%Vh1_h%r_coef")
      NULLIFY(ecoul_1c(iat)%Vh1_h)

      DEALLOCATE(ecoul_1c(iat)%Vh1_s%r_coef,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "ecoul_1c%Vh1_s%r_coef")
      NULLIFY(ecoul_1c(iat)%Vh1_s)

      DEALLOCATE(ecoul_1c(iat)%aVh1b_h%r_coef,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "ecoul_1c%aVh1b_h%r_coef")
      NULLIFY(ecoul_1c(iat)%aVh1b_h)

      DEALLOCATE(ecoul_1c(iat)%aVh1b_s%r_coef,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "ecoul_1c%aVh1b_s%r_coef")
      NULLIFY(ecoul_1c(iat)%aVh1b_s)

      DEALLOCATE(ecoul_1c(iat)%aVh1b_0%r_coef,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "ecoul_1c%aVh1b_0%r_coef")
      NULLIFY(ecoul_1c(iat)%aVh1b_0)

    END DO

    DEALLOCATE(ecoul_1c, STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "ecoul_1c")

  END SUBROUTINE deallocate_ecoul_1center

! *****************************************************************************

  SUBROUTINE get_hartree_local(hartree_local,ecoul_1c)

    TYPE(hartree_local_type), INTENT(IN), TARGET  :: hartree_local
    TYPE(ecoul_1center_type), DIMENSION(:), &
                            OPTIONAL, POINTER    :: ecoul_1c

    CHARACTER(len=*), PARAMETER :: routineN = 'get_hartree_local'

!   ---------------------------------------------------------------------------

    IF(PRESENT(ecoul_1c)) ecoul_1c => hartree_local%ecoul_1c
    
  END SUBROUTINE get_hartree_local

! *****************************************************************************

  SUBROUTINE hartree_local_create(hartree_local, error)

    TYPE(hartree_local_type), POINTER         :: hartree_local
    TYPE(cp_error_type), INTENT(inout), &
           OPTIONAL                           :: error

    CHARACTER(len=*), PARAMETER :: routineN = "hartree_local_create" 

    INTEGER                                    :: stat
    LOGICAL                                    :: failure

!   ---------------------------------------------------------------------------
    failure=.FALSE.

    ALLOCATE(hartree_local, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineN,error,failure)
  
    NULLIFY (hartree_local%ecoul_1c)

  END SUBROUTINE hartree_local_create

! *****************************************************************************

  SUBROUTINE hartree_local_release(hartree_local, error)

    TYPE(hartree_local_type), POINTER         :: hartree_local
    TYPE(cp_error_type), INTENT(inout), &
                OPTIONAL                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = "hartree_local_release" 

    INTEGER                                  :: stat

!   ---------------------------------------------------------------------------
 
    IF (ASSOCIATED(hartree_local)) THEN
      IF (ASSOCIATED(hartree_local%ecoul_1c)) THEN
        CALL deallocate_ecoul_1center(hartree_local%ecoul_1c)
      END IF
      DEALLOCATE(hartree_local,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineN,error)
    END IF

  END SUBROUTINE hartree_local_release

! *****************************************************************************

  SUBROUTINE set_hartree_local(hartree_local,ecoul_1c)

    TYPE(hartree_local_type), INTENT(INOUT), TARGET  :: hartree_local
    TYPE(ecoul_1center_type), DIMENSION(:), &
                            OPTIONAL, POINTER        :: ecoul_1c

    CHARACTER(len=*), PARAMETER :: routineN = 'set_hartree_local'

!   ---------------------------------------------------------------------------

    IF(PRESENT(ecoul_1c)) hartree_local%ecoul_1c => ecoul_1c
    
  END SUBROUTINE set_hartree_local

! *****************************************************************************

END MODULE hartree_local_types

! ****************************************************************************

