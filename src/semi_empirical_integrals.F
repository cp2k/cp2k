!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Set of wrappers for semi-empirical analytical/numerical Integrals 
!>        routines
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   04.2008
! *****************************************************************************
MODULE semi_empirical_integrals

  USE f77_blas
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_ana,          ONLY: corecore_ana,&
                                             rotint_ana,&
                                             rotnuc_ana
  USE semi_empirical_int_num,          ONLY: corecore_num,&
                                             dcorecore_num,&
                                             drotint_num,&
                                             drotnuc_num,&
                                             rotint_num,&
                                             rotnuc_num
  USE semi_empirical_types,            ONLY: semi_empirical_type,&
                                             taper_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_integrals'
  PUBLIC :: rotint, rotnuc, drotint, drotnuc, corecore, dcorecore

CONTAINS

! *****************************************************************************
!> \brief wrapper for numerical/analytical routines
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   04.2008
! *****************************************************************************
  SUBROUTINE rotint (sepi,sepj,rij,w,anag,shortrange,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: w
    LOGICAL                                  :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotint', &
      routineP = moduleN//':'//routineN

    IF (anag) THEN
       CALL rotint_ana (sepi,sepj,rij,w,shortrange=shortrange,taper=taper,error=error)
    ELSE
       CALL rotint_num (sepi,sepj,rij,w,shortrange=shortrange,taper=taper,error=error)
    END IF
  END SUBROUTINE rotint

! *****************************************************************************
!> \brief wrapper for numerical/analytical routines
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   04.2008
! *****************************************************************************
  SUBROUTINE rotnuc (sepi,sepj,rij,e1b,e2a,itype,anag,shortrange,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: e1b, e2a
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotnuc', &
      routineP = moduleN//':'//routineN

    IF (anag) THEN
       CALL rotnuc_ana (sepi, sepj, rij, e1b=e1b, e2a=e2a, itype=itype,&
            shortrange=shortrange, taper=taper, error=error)
    ELSE
       CALL rotnuc_num (sepi, sepj, rij, e1b=e1b, e2a=e2a, itype=itype,&
            shortrange=shortrange, taper=taper, error=error)    
    END IF

  END SUBROUTINE rotnuc

! *****************************************************************************
!> \brief wrapper for numerical/analytical routines
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   04.2008
! *****************************************************************************
  SUBROUTINE corecore (sepi,sepj,rij,enuc,itype,anag,shortrange,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), INTENT(OUT)                    :: enuc
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'corecore', &
      routineP = moduleN//':'//routineN

    IF (anag) THEN
       CALL corecore_ana (sepi,sepj,rij,enuc=enuc,itype=itype,shortrange=shortrange,&
            taper=taper, error=error)
    ELSE
       CALL corecore_num (sepi,sepj,rij,enuc=enuc,itype=itype,shortrange=shortrange,&
            taper=taper, error=error)    
    END IF

  END SUBROUTINE corecore

! *****************************************************************************
!> \brief wrapper for numerical/analytical routines
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   04.2008
! *****************************************************************************
  SUBROUTINE drotint (sepi,sepj,rij,dw,delta,anag,shortrange,taper, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: dw
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'drotint', &
      routineP = moduleN//':'//routineN

    IF (anag) THEN
       CALL  rotint_ana(sepi,sepj,rij,dw=dw,   shortrange=shortrange,taper=taper,error=error)
    ELSE
       CALL drotint_num(sepi,sepj,rij,dw,delta,shortrange=shortrange,taper=taper,error=error)
    END IF

  END SUBROUTINE drotint

! *****************************************************************************
!> \brief wrapper for numerical/analytical routines
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   04.2008
! *****************************************************************************
  SUBROUTINE drotnuc (sepi,sepj,rij,de1b,de2a,itype,delta,anag,shortrange,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(OUT), &
      OPTIONAL                               :: de1b, de2a
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'drotnuc', &
      routineP = moduleN//':'//routineN

    IF (anag) THEN
       CALL rotnuc_ana  (sepi, sepj, rij, de1b=de1b, de2a=de2a, itype=itype,&
            shortrange=shortrange,taper=taper,error=error)
    ELSE
       CALL drotnuc_num (sepi, sepj, rij, de1b=de1b, de2a=de2a, itype=itype,&
            delta=delta,shortrange=shortrange,taper=taper,error=error)
    END IF

  END SUBROUTINE drotnuc

! *****************************************************************************
!> \brief wrapper for numerical/analytical routines
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   04.2008
! *****************************************************************************
  SUBROUTINE dcorecore (sepi,sepj,rij,denuc,itype,delta,anag,shortrange,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: denuc
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: anag, shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dcorecore', &
      routineP = moduleN//':'//routineN

    IF (anag) THEN
       CALL corecore_ana (sepi,sepj,rij,denuc=denuc,itype=itype,shortrange=shortrange,&
            taper=taper, error=error)
    ELSE
       CALL dcorecore_num (sepi,sepj,rij,denuc=denuc,delta=delta,itype=itype,&
            shortrange=shortrange,taper=taper, error=error)    
    END IF

  END SUBROUTINE dcorecore

END MODULE semi_empirical_integrals
