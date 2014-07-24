!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE qs_rho_atom_types

  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho_atom_types'

  TYPE rho_atom_coeff
    REAL(dp), DIMENSION(:,:), POINTER :: r_coef
  END TYPE rho_atom_coeff

  TYPE rho_atom_type
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER  :: cpc_h, cpc_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER  :: rho_rad_h,&
                                           rho_rad_s,&
                                           vrho_rad_h,&
                                           vrho_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER  :: drho_rad_h,&
                                           drho_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:,:),&
                               POINTER  :: rho_rad_h_d,&
                                           rho_rad_s_d
    TYPE(rho_atom_coeff), DIMENSION(:,:),&
                               POINTER  :: trho_rad_h,&
                                           trho_rad_s

    INTEGER                             :: rhoa_of_atom
    REAL(dp)                            :: exc_h,&
                                           exc_s
    TYPE(rho_atom_coeff),DIMENSION(:),&
                              POINTER   :: ga_Vlocal_gb_h,&
                                           ga_Vlocal_gb_s
  END TYPE rho_atom_type

  TYPE rho_atom_p_type
     TYPE(rho_atom_type), POINTER       :: rho_atom
  END TYPE rho_atom_p_type

  PUBLIC :: allocate_rho_atom_set,&
            deallocate_rho_atom_set,&
            get_rho_atom

  PUBLIC :: rho_atom_coeff, rho_atom_type

 CONTAINS

! *****************************************************************************
!> \brief ...
!> \param rho_atom_set ...
!> \param natom ...
! *****************************************************************************
  SUBROUTINE allocate_rho_atom_set(rho_atom_set,natom)

    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    INTEGER, INTENT(IN)                      :: natom

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_rho_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, istat

    IF (ASSOCIATED(rho_atom_set)) THEN
      CALL deallocate_rho_atom_set(rho_atom_set)
    END IF

    ALLOCATE (rho_atom_set(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                            "allocate_rho_atom_set",natom*int_size)

    DO iat = 1,natom

      NULLIFY (rho_atom_set(iat)%cpc_h)
      NULLIFY (rho_atom_set(iat)%cpc_s)
      NULLIFY (rho_atom_set(iat)%rho_rad_h)
      NULLIFY (rho_atom_set(iat)%rho_rad_s)
      NULLIFY (rho_atom_set(iat)%drho_rad_h)
      NULLIFY (rho_atom_set(iat)%drho_rad_s)
      NULLIFY (rho_atom_set(iat)%rho_rad_h_d)
      NULLIFY (rho_atom_set(iat)%rho_rad_s_d)
      NULLIFY (rho_atom_set(iat)%vrho_rad_h)
      NULLIFY (rho_atom_set(iat)%vrho_rad_s)
      NULLIFY (rho_atom_set(iat)%trho_rad_h)
      NULLIFY (rho_atom_set(iat)%trho_rad_s)
      NULLIFY (rho_atom_set(iat)%ga_Vlocal_gb_h)
      NULLIFY (rho_atom_set(iat)%ga_Vlocal_gb_s)

    ENDDO

  END SUBROUTINE allocate_rho_atom_set

! *****************************************************************************
!> \brief ...
!> \param rho_atom_set ...
! *****************************************************************************
  SUBROUTINE deallocate_rho_atom_set(rho_atom_set)

    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_rho_atom_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iat, istat, j, n, natom

    IF (ASSOCIATED(rho_atom_set)) THEN

      natom = SIZE(rho_atom_set)

      DO iat = 1,natom
        IF(ASSOCIATED(rho_atom_set(iat)%cpc_h)) THEN
          IF(ASSOCIATED(rho_atom_set(iat)%cpc_h(1)%r_coef)) THEN
            n = SIZE(rho_atom_set(iat)%cpc_h,1)
            DO i = 1,n
              DEALLOCATE(rho_atom_set(iat)%cpc_h(i)%r_coef)
              DEALLOCATE(rho_atom_set(iat)%cpc_s(i)%r_coef)
            END DO
          END IF
          DEALLOCATE (rho_atom_set(iat)%cpc_h)
          DEALLOCATE (rho_atom_set(iat)%cpc_s)
        END IF
        IF(ASSOCIATED(rho_atom_set(iat)%ga_Vlocal_gb_h)) THEN
          IF(ASSOCIATED(rho_atom_set(iat)%ga_Vlocal_gb_h(1)%r_coef)) THEN
            n = SIZE(rho_atom_set(iat)%ga_Vlocal_gb_h,1)
            DO i = 1,n
              DEALLOCATE(rho_atom_set(iat)%ga_Vlocal_gb_h(i)%r_coef)
              DEALLOCATE(rho_atom_set(iat)%ga_Vlocal_gb_s(i)%r_coef)
            END DO
          END IF
          DEALLOCATE (rho_atom_set(iat)%ga_Vlocal_gb_h)
          DEALLOCATE (rho_atom_set(iat)%ga_Vlocal_gb_s)
        END IF

        IF(ASSOCIATED(rho_atom_set(iat)%drho_rad_h)) THEN
          IF(ASSOCIATED(rho_atom_set(iat)%drho_rad_h(1)%r_coef)) THEN
            n = SIZE(rho_atom_set(iat)%drho_rad_h,1)
            DO i = 1,n
              DEALLOCATE(rho_atom_set(iat)%drho_rad_h(i)%r_coef)
              DEALLOCATE(rho_atom_set(iat)%drho_rad_s(i)%r_coef)
              DO j = 1,3
                DEALLOCATE(rho_atom_set(iat)%rho_rad_h_d(j,i)%r_coef)
                DEALLOCATE(rho_atom_set(iat)%rho_rad_s_d(j,i)%r_coef)
              END DO
            ENDDO
          END IF
          DEALLOCATE (rho_atom_set(iat)%drho_rad_h)
          DEALLOCATE (rho_atom_set(iat)%drho_rad_s)
          DEALLOCATE (rho_atom_set(iat)%rho_rad_h_d)
          DEALLOCATE (rho_atom_set(iat)%rho_rad_s_d)
        ENDIF

        IF(ASSOCIATED(rho_atom_set(iat)%rho_rad_h)) THEN
          IF(ASSOCIATED(rho_atom_set(iat)%rho_rad_h(1)%r_coef)) THEN
            n = SIZE(rho_atom_set(iat)%rho_rad_h)
            DO i = 1,n
              DEALLOCATE(rho_atom_set(iat)%rho_rad_h(i)%r_coef)
            ENDDO
          ENDIF
          DEALLOCATE (rho_atom_set(iat)%rho_rad_h)
        ENDIF

        IF(ASSOCIATED(rho_atom_set(iat)%rho_rad_s)) THEN
          IF(ASSOCIATED(rho_atom_set(iat)%rho_rad_s(1)%r_coef)) THEN
            n = SIZE(rho_atom_set(iat)%rho_rad_s)
            DO i = 1,n
              DEALLOCATE(rho_atom_set(iat)%rho_rad_s(i)%r_coef)
            ENDDO
          ENDIF
          DEALLOCATE (rho_atom_set(iat)%rho_rad_s)
        ENDIF

        IF(ASSOCIATED(rho_atom_set(iat)%vrho_rad_h)) THEN
          IF(ASSOCIATED(rho_atom_set(iat)%vrho_rad_h(1)%r_coef)) THEN
            n = SIZE(rho_atom_set(iat)%vrho_rad_h)
            DO i = 1,n
              DEALLOCATE(rho_atom_set(iat)%vrho_rad_h(i)%r_coef)
            ENDDO
          END IF
          DEALLOCATE (rho_atom_set(iat)%vrho_rad_h)
        ENDIF

        IF(ASSOCIATED(rho_atom_set(iat)%vrho_rad_s)) THEN
          IF(ASSOCIATED(rho_atom_set(iat)%vrho_rad_s(1)%r_coef)) THEN
            n = SIZE(rho_atom_set(iat)%vrho_rad_s)
            DO i = 1,n
              DEALLOCATE(rho_atom_set(iat)%vrho_rad_s(i)%r_coef)
            ENDDO
          END IF
          DEALLOCATE (rho_atom_set(iat)%vrho_rad_s)
        ENDIF

        IF(ASSOCIATED(rho_atom_set(iat)%trho_rad_h)) THEN
          IF(ASSOCIATED(rho_atom_set(iat)%trho_rad_h(1,1)%r_coef)) THEN
            n = SIZE(rho_atom_set(iat)%trho_rad_h,2)
            DO i = 1,n
              DEALLOCATE(rho_atom_set(iat)%trho_rad_h(1,i)%r_coef)
              DEALLOCATE(rho_atom_set(iat)%trho_rad_h(2,i)%r_coef)
              DEALLOCATE(rho_atom_set(iat)%trho_rad_h(3,i)%r_coef)
            ENDDO
          ENDIF
          DEALLOCATE (rho_atom_set(iat)%trho_rad_h)
        ENDIF

        IF(ASSOCIATED(rho_atom_set(iat)%trho_rad_s)) THEN
          IF(ASSOCIATED(rho_atom_set(iat)%trho_rad_s(1,1)%r_coef)) THEN
            n = SIZE(rho_atom_set(iat)%trho_rad_s,2)
            DO i = 1,n
              DEALLOCATE(rho_atom_set(iat)%trho_rad_s(1,i)%r_coef)
              DEALLOCATE(rho_atom_set(iat)%trho_rad_s(2,i)%r_coef)
              DEALLOCATE(rho_atom_set(iat)%trho_rad_s(3,i)%r_coef)
            ENDDO
          ENDIF
          DEALLOCATE (rho_atom_set(iat)%trho_rad_s)
        ENDIF

      ENDDO

      DEALLOCATE (rho_atom_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "rho_atom_set")

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer rho_atom_set is not associated and "//&
                        "cannot be deallocated")

    END IF

  END SUBROUTINE deallocate_rho_atom_set

! *****************************************************************************
!> \brief ...
!> \param rho_atom ...
!> \param cpc_h ...
!> \param cpc_s ...
!> \param rho_rad_h ...
!> \param rho_rad_s ...
!> \param drho_rad_h ...
!> \param drho_rad_s ...
!> \param vrho_rad_h ...
!> \param vrho_rad_s ...
!> \param rho_rad_h_d ...
!> \param rho_rad_s_d ...
!> \param trho_rad_h ...
!> \param trho_rad_s ...
!> \param ga_Vlocal_gb_h ...
!> \param ga_Vlocal_gb_s ...
! *****************************************************************************
  SUBROUTINE get_rho_atom(rho_atom,cpc_h,cpc_s,rho_rad_h,rho_rad_s,&
                          drho_rad_h,drho_rad_s,vrho_rad_h,vrho_rad_s,&
                          rho_rad_h_d,rho_rad_s_d,trho_rad_h,trho_rad_s,&
                          ga_Vlocal_gb_h,ga_Vlocal_gb_s)

    TYPE(rho_atom_type), POINTER             :: rho_atom
    TYPE(rho_atom_coeff), DIMENSION(:), &
      OPTIONAL, POINTER                      :: cpc_h, cpc_s, rho_rad_h, &
                                                rho_rad_s, drho_rad_h, &
                                                drho_rad_s, vrho_rad_h, &
                                                vrho_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: rho_rad_h_d, rho_rad_s_d, &
                                                trho_rad_h, trho_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ga_Vlocal_gb_h, ga_Vlocal_gb_s

    CHARACTER(len=*), PARAMETER :: routineN = 'get_rho_atom', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(rho_atom)) THEN
      IF (PRESENT(cpc_h)) cpc_h => rho_atom%cpc_h
      IF (PRESENT(cpc_s)) cpc_s => rho_atom%cpc_s
      IF (PRESENT(rho_rad_h)) rho_rad_h => rho_atom%rho_rad_h
      IF (PRESENT(rho_rad_s)) rho_rad_s => rho_atom%rho_rad_s
      IF (PRESENT(drho_rad_h)) drho_rad_h => rho_atom%drho_rad_h
      IF (PRESENT(drho_rad_s)) drho_rad_s => rho_atom%drho_rad_s
      IF (PRESENT(rho_rad_h_d)) rho_rad_h_d => rho_atom%rho_rad_h_d
      IF (PRESENT(rho_rad_s_d)) rho_rad_s_d => rho_atom%rho_rad_s_d
      IF (PRESENT(vrho_rad_h)) vrho_rad_h => rho_atom%vrho_rad_h
      IF (PRESENT(vrho_rad_s)) vrho_rad_s => rho_atom%vrho_rad_s
      IF (PRESENT(trho_rad_h)) trho_rad_h => rho_atom%trho_rad_h
      IF (PRESENT(trho_rad_s)) trho_rad_s => rho_atom%trho_rad_s
      IF (PRESENT(ga_Vlocal_gb_h)) ga_Vlocal_gb_h => rho_atom%ga_Vlocal_gb_h
      IF (PRESENT(ga_Vlocal_gb_s)) ga_Vlocal_gb_s => rho_atom%ga_Vlocal_gb_s
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer rho_atom is not associated")
    END IF

  END SUBROUTINE  get_rho_atom

END MODULE qs_rho_atom_types
