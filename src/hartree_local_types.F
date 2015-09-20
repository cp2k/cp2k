!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE hartree_local_types

  USE kinds,                           ONLY: dp
  USE qs_rho_atom_types,               ONLY: rho_atom_coeff
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hartree_local_types'

! *** Define the ecoul_1center_type ***

! *****************************************************************************
  TYPE ecoul_1center_type
    TYPE(rho_atom_coeff), POINTER   :: Vh1_h, Vh1_s
    REAL(dp)                        :: ecoul_1_h, ecoul_1_s, &
                                       ecoul_1_z, ecoul_1_0
  END TYPE ecoul_1center_type

! *****************************************************************************
  TYPE hartree_local_type
    TYPE(ecoul_1center_type), &
            DIMENSION(:), POINTER   :: ecoul_1c
  END TYPE hartree_local_type

! *** Public subroutines ***

  PUBLIC :: allocate_ecoul_1center,&
            get_hartree_local, hartree_local_create,&
            hartree_local_release, set_ecoul_1c,&
            set_hartree_local

! *** Public data types ***

  PUBLIC :: ecoul_1center_type, hartree_local_type

  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param ecoul_1c ...
!> \param natom ...
! *****************************************************************************
  SUBROUTINE  allocate_ecoul_1center(ecoul_1c,natom)

    TYPE(ecoul_1center_type), DIMENSION(:), &
      POINTER                                :: ecoul_1c
    INTEGER, INTENT(IN)                      :: natom

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_ecoul_1center', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, istat
    LOGICAL                                  :: failure

    IF(ASSOCIATED(ecoul_1c)) THEN
      CALL deallocate_ecoul_1center(ecoul_1c)
    END IF

    ALLOCATE(ecoul_1c(natom), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)

    DO iat = 1,natom
      ALLOCATE(ecoul_1c(iat)%Vh1_h,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
      NULLIFY(ecoul_1c(iat)%Vh1_h%r_coef)
      ALLOCATE(ecoul_1c(iat)%Vh1_s,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
      NULLIFY(ecoul_1c(iat)%Vh1_s%r_coef)
    END DO

  END SUBROUTINE allocate_ecoul_1center

! *****************************************************************************
!> \brief ...
!> \param ecoul_1c ...
! *****************************************************************************
  SUBROUTINE  deallocate_ecoul_1center(ecoul_1c)

    TYPE(ecoul_1center_type), DIMENSION(:), &
      POINTER                                :: ecoul_1c

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_ecoul_1center', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, istat, natom
    LOGICAL                                  :: failure

    natom = SIZE(ecoul_1c,1)

    DO iat= 1,natom
      IF (ASSOCIATED(ecoul_1c(iat)%Vh1_h%r_coef)) THEN
          DEALLOCATE(ecoul_1c(iat)%Vh1_h%r_coef,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)
      ENDIF
      DEALLOCATE(ecoul_1c(iat)%Vh1_h,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)

      IF (ASSOCIATED(ecoul_1c(iat)%Vh1_s%r_coef)) THEN
          DEALLOCATE(ecoul_1c(iat)%Vh1_s%r_coef,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)
      ENDIF
      DEALLOCATE(ecoul_1c(iat)%Vh1_s,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)

    END DO

    DEALLOCATE(ecoul_1c, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)

  END SUBROUTINE deallocate_ecoul_1center

! *****************************************************************************
!> \brief ...
!> \param ecoul_1c ...
!> \param iatom ...
!> \param Vh1_h ...
!> \param Vh1_s ...
! *****************************************************************************
  SUBROUTINE get_ecoul_1c(ecoul_1c,iatom,Vh1_h,Vh1_s)

    TYPE(ecoul_1center_type), DIMENSION(:), &
      POINTER                                :: ecoul_1c
    INTEGER, INTENT(IN), OPTIONAL            :: iatom
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: Vh1_h, Vh1_s

    IF(PRESENT(iatom)) THEN
      IF(PRESENT(Vh1_h)) Vh1_h => ecoul_1c(iatom)%Vh1_h%r_coef
      IF(PRESENT(Vh1_s)) Vh1_s => ecoul_1c(iatom)%Vh1_s%r_coef
    END IF

  END SUBROUTINE get_ecoul_1c

! *****************************************************************************
!> \brief ...
!> \param hartree_local ...
!> \param ecoul_1c ...
! *****************************************************************************
  SUBROUTINE get_hartree_local(hartree_local,ecoul_1c)

    TYPE(hartree_local_type), POINTER        :: hartree_local
    TYPE(ecoul_1center_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ecoul_1c

    CHARACTER(len=*), PARAMETER :: routineN = 'get_hartree_local', &
      routineP = moduleN//':'//routineN

    IF(PRESENT(ecoul_1c)) ecoul_1c => hartree_local%ecoul_1c

  END SUBROUTINE get_hartree_local

! *****************************************************************************
!> \brief ...
!> \param hartree_local ...
! *****************************************************************************
  SUBROUTINE hartree_local_create(hartree_local)

    TYPE(hartree_local_type), POINTER        :: hartree_local

    CHARACTER(len=*), PARAMETER :: routineN = 'hartree_local_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(hartree_local, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    NULLIFY (hartree_local%ecoul_1c)

  END SUBROUTINE hartree_local_create

! *****************************************************************************
!> \brief ...
!> \param hartree_local ...
! *****************************************************************************
  SUBROUTINE hartree_local_release(hartree_local)

    TYPE(hartree_local_type), POINTER        :: hartree_local

    CHARACTER(len=*), PARAMETER :: routineN = 'hartree_local_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

    IF (ASSOCIATED(hartree_local)) THEN
      IF (ASSOCIATED(hartree_local%ecoul_1c)) THEN
        CALL deallocate_ecoul_1center(hartree_local%ecoul_1c)
      END IF

      DEALLOCATE(hartree_local,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP)
    END IF

  END SUBROUTINE hartree_local_release

! *****************************************************************************
!> \brief ...
!> \param ecoul_1c ...
!> \param iatom ...
!> \param ecoul_1_h ...
!> \param ecoul_1_s ...
!> \param ecoul_1_z ...
!> \param ecoul_1_0 ...
! *****************************************************************************
  SUBROUTINE set_ecoul_1c(ecoul_1c,iatom,ecoul_1_h,ecoul_1_s,ecoul_1_z,ecoul_1_0)

    TYPE(ecoul_1center_type), DIMENSION(:), &
      POINTER                                :: ecoul_1c
    INTEGER, INTENT(IN), OPTIONAL            :: iatom
    REAL(dp), INTENT(IN), OPTIONAL           :: ecoul_1_h, ecoul_1_s, &
                                                ecoul_1_z, ecoul_1_0

    IF(PRESENT(iatom)) THEN
      IF(PRESENT(ecoul_1_h)) ecoul_1c(iatom)%ecoul_1_h=ecoul_1_h
      IF(PRESENT(ecoul_1_s)) ecoul_1c(iatom)%ecoul_1_s=ecoul_1_s
      IF(PRESENT(ecoul_1_0)) ecoul_1c(iatom)%ecoul_1_0=ecoul_1_0
      IF(PRESENT(ecoul_1_z)) ecoul_1c(iatom)%ecoul_1_z=ecoul_1_z
    END IF

    END SUBROUTINE set_ecoul_1c

! *****************************************************************************
!> \brief ...
!> \param hartree_local ...
!> \param ecoul_1c ...
! *****************************************************************************
  SUBROUTINE set_hartree_local(hartree_local,ecoul_1c)

    TYPE(hartree_local_type), POINTER        :: hartree_local
    TYPE(ecoul_1center_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ecoul_1c

    CHARACTER(len=*), PARAMETER :: routineN = 'set_hartree_local', &
      routineP = moduleN//':'//routineN

    IF(PRESENT(ecoul_1c)) hartree_local%ecoul_1c => ecoul_1c

  END SUBROUTINE set_hartree_local

END MODULE hartree_local_types

