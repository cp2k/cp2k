!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****h* cp2k/l_utils *
!!
!!   NAME
!!     l_utils
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     JGH [08.07.2006] remove pol{xyz}, pzyx, ipzyx from type
!!                      allocate/deallocate arrays when needed
!!
!! ***
!****************************************************************************

MODULE l_utils
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'l_utils'

  PUBLIC :: l_info_type, init_l_info, return_l_info, destroy_l_info, &
            l_info_release, l_set_info_type, l_info_retain, finish_l_info

  TYPE l_set_info_type
        INTEGER                        :: lx_max,lxy_max,lxyz_max
        INTEGER, POINTER, DIMENSION(:) :: ly_max,lz_max

        INTEGER, POINTER, DIMENSION(:,:)          :: map
        REAL(KIND = dp), POINTER, DIMENSION(:,:)  :: alpha,dpy,dpz
        REAL(KIND = dp), POINTER, DIMENSION(:)    :: pzyx
  END TYPE

  TYPE l_info_type
        INTEGER :: lmax
        INTEGER :: cmax
        INTEGER :: nthread
        INTEGER :: ref_count
        TYPE(l_set_info_type), POINTER, DIMENSION(:,:,:) :: DATA
  END TYPE

!****************************************************************************

CONTAINS

!****************************************************************************

  SUBROUTINE return_l_info(l_info,la_min,la_max,lb_min,lb_max,ithread, &
                           lx_max,lxy_max,lxyz_max,ly_max,lz_max,map,&
                           polx,poly,polz,dpy,dpz,alpha,pzyx,cmax,ipzyx,error)

    TYPE(l_info_type)                        :: l_info
    INTEGER                                  :: la_min, la_max, lb_min, &
                                                lb_max, ithread, lx_max, &
                                                lxy_max, lxyz_max
    INTEGER, DIMENSION(:), POINTER           :: ly_max, lz_max
    INTEGER, DIMENSION(:, :), POINTER        :: map
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: polx, poly, polz, dpy, dpz, &
                                                alpha
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pzyx
    INTEGER                                  :: cmax
    INTEGER, DIMENSION(:, :, :, :, :, :), &
      OPTIONAL, POINTER                      :: ipzyx
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'return_l_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: coef_max, istat, la, lb, lxa, &
                                                lxb, lxyz, lya, lyb, lza, lzb
    LOGICAL                                  :: failure

!--------------------------------------------------------------------------!

    IF (la_max.gt.l_info%lmax .OR. lb_max.gt.l_info%lmax) THEN
      CALL stop_program("return_l_info","called with l too large")
    ENDIF
    IF (ithread.ge.l_info%nthread) THEN
      CALL stop_program("return_l_info","called with ithread too large")
    ENDIF
    la = (la_max*(la_max+1))/2 + la_min + 1
    lb = (lb_max*(lb_max+1))/2 + lb_min + 1
    lx_max=l_info%data(la,lb,ithread)%lx_max
    lxy_max=l_info%data(la,lb,ithread)%lxy_max
    lxyz_max=l_info%data(la,lb,ithread)%lxyz_max
    ly_max=>l_info%data(la,lb,ithread)%ly_max
    lz_max=>l_info%data(la,lb,ithread)%lz_max
    map=>l_info%data(la,lb,ithread)%map
    dpy=>l_info%data(la,lb,ithread)%dpy
    dpz=>l_info%data(la,lb,ithread)%dpz
    alpha=>l_info%data(la,lb,ithread)%alpha
    cmax=l_info%cmax

    coef_max=la_max+lb_max+1

    ALLOCATE ( polx(coef_max,-cmax:cmax), STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( poly(lxy_max,-cmax:cmax), STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( polz(lxyz_max,-cmax:cmax), STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( pzyx(lxyz_max), STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(ipzyx)) THEN
      ALLOCATE(ipzyx(0:la_max,0:la_max,0:la_max,&
                     0:lb_max,0:lb_max,0:lb_max), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      lxyz=0
      DO lxa=0,la_max
        DO lxb=0,lb_max
          DO lya=0,la_max-lxa
            DO lyb=0,lb_max-lxb
              DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
                DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
                  lxyz=lxyz+1
                  ipzyx(lxa,lya,lza,lxb,lyb,lzb)=lxyz
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    END IF

  END SUBROUTINE return_l_info

!****************************************************************************

  SUBROUTINE finish_l_info(polx,poly,polz,pzyx,ipzyx,error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: polx, poly, polz
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pzyx
    INTEGER, DIMENSION(:, :, :, :, :, :), &
      OPTIONAL, POINTER                      :: ipzyx
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'finish_l_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    IF ( ASSOCIATED(polx) ) THEN
       DEALLOCATE(polx, STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF ( ASSOCIATED(poly) ) THEN
       DEALLOCATE(poly, STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF ( ASSOCIATED(polz) ) THEN
       DEALLOCATE(polz, STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF ( ASSOCIATED(pzyx) ) THEN
       DEALLOCATE(pzyx, STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    IF (PRESENT(ipzyx)) THEN
      IF ( ASSOCIATED(ipzyx) ) THEN
         DEALLOCATE(ipzyx, STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
    END IF

  END SUBROUTINE finish_l_info

!****************************************************************************


! needs a nthread agrument because it provides pre-allocated data streactures
! for the collocate routines hence, in some way, these data structures need
! the 'thread-private' attribute really not so elegant...
! the nthread argument is optional, the right number of threads should be
! determined automatically

  SUBROUTINE init_l_info(l_info,lmax,cmax,nthread,lrho0,error)

    TYPE(l_info_type)                        :: l_info
    INTEGER                                  :: lmax, cmax
    INTEGER, OPTIONAL                        :: nthread
    LOGICAL, INTENT(in), OPTIONAL            :: lrho0
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: handle, ithread, la, la_max, &
                                                la_min, lb, lb_max, lb_min, &
                                                nthread_l
    LOGICAL                                  :: my_lrho0

!$  integer                                  :: omp_get_max_threads

    CHARACTER(len=*), PARAMETER :: routineN = 'init_l_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: istat

    !--------------------------------------------------------------------------!

    my_lrho0 = .FALSE.
    IF(PRESENT(lrho0)) my_lrho0 = lrho0

    CALL timeset("init_l_info","I","",handle)
    IF (PRESENT(nthread)) THEN
       nthread_l=nthread
    ELSE
       nthread_l=1
!$     nthread_l=omp_get_max_threads()
    ENDIF
    l_info%lmax=lmax
    l_info%cmax=cmax
    l_info%nthread=nthread_l
    l_info%ref_count=1

    IF(my_lrho0) THEN

      la = ((lmax+1)*(lmax+2))/2
      lb = 1
      ALLOCATE(l_info%data(la,lb,0:nthread_l-1), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO ithread=0,nthread_l-1
        DO lb_max=0,0
          DO lb_min=0,lb_max
            lb = (lb_max*(lb_max+1))/2 + lb_min + 1
            DO la_max=0,lmax
              DO la_min=0,la_max
                la = (la_max*(la_max+1))/2 + la_min + 1
                CALL init_l_set_info_type(la_min,la_max,lb_min,lb_max,cmax, &
                          l_info%data(la,lb,ithread),error=error)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ELSE

      la = ((lmax+1)*(lmax+2))/2
      lb = ((lmax+1)*(lmax+2))/2
      ALLOCATE(l_info%data(la,lb,0:nthread_l-1), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO ithread=0,nthread_l-1
        DO lb_max=0,lmax
          DO lb_min=0,lb_max
            lb = (lb_max*(lb_max+1))/2 + lb_min + 1
            DO la_max=0,lmax
              DO la_min=0,la_max
                la = (la_max*(la_max+1))/2 + la_min + 1
                CALL init_l_set_info_type(la_min,la_max,lb_min,lb_max,cmax, &
                          l_info%data(la,lb,ithread),error=error)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE init_l_info

!****************************************************************************

  SUBROUTINE destroy_l_info(l_info,error)

    TYPE(l_info_type)                        :: l_info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'destroy_l_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, ithread, la, la_max, &
                                                la_min, lb
    LOGICAL                                  :: failure

!--------------------------------------------------------------------------!

    DO ithread=0,l_info%nthread-1
      DO lb=1, SIZE(l_info%data,2)
        DO la_max=0,l_info%lmax
          DO la_min=0,la_max
             la = (la_max*(la_max+1))/2 + la_min + 1
             CALL destroy_l_set_info_type(l_info%data(la,lb,ithread),error=error)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(l_info%data,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE destroy_l_info

!****************************************************************************

  SUBROUTINE init_l_set_info_type(la_min,la_max,lb_min,lb_max,cmax,&
                                  l_set_info,error)

    INTEGER                                  :: la_min, la_max, lb_min, &
                                                lb_max, cmax
    TYPE(l_set_info_type)                    :: l_set_info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_l_set_info_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: coef_max, istat, lx, lxa, &
                                                lxb, lxy, lxyz, lya, lyb, &
                                                lza, lzb
    LOGICAL                                  :: failure

!--------------------------------------------------------------------------!

    coef_max=la_max+lb_max+1

    l_set_info%lxyz_max=0
    l_set_info%lxy_max=0
    l_set_info%lx_max=0
    DO lxa=0,la_max
      DO lxb=0,lb_max
        l_set_info%lx_max=l_set_info%lx_max+1
        DO lya=0,la_max-lxa
          DO lyb=0,lb_max-lxb
            l_set_info%lxy_max=l_set_info%lxy_max+1
            DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
              DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
                 l_set_info%lxyz_max=l_set_info%lxyz_max+1
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ALLOCATE(l_set_info%ly_max(l_set_info%lx_max), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(l_set_info%lz_max(l_set_info%lxy_max), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(l_set_info%map(-cmax:cmax,1:3), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(l_set_info%alpha(coef_max,l_set_info%lx_max), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(l_set_info%dpy(0:lb_max,0:la_max), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(l_set_info%dpz(0:lb_max,0:la_max), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    l_set_info%ly_max(:)=0
    l_set_info%lz_max(:)=0
    lx=0
    lxy=0
    lxyz=0
    DO lxa=0,la_max
      DO lxb=0,lb_max
        lx=lx+1
        DO lya=0,la_max-lxa
          DO lyb=0,lb_max-lxb
            lxy=lxy+1
            l_set_info%ly_max(lx)=l_set_info%ly_max(lx)+1
            DO lza=MAX(la_min-lxa-lya,0),la_max-lxa-lya
              DO lzb=MAX(lb_min-lxb-lyb,0),lb_max-lxb-lyb
                lxyz=lxyz+1
                l_set_info%lz_max(lxy)=l_set_info%lz_max(lxy)+1
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

   END SUBROUTINE init_l_set_info_type

!****************************************************************************

   SUBROUTINE destroy_l_set_info_type(l_set_info,error)

    TYPE(l_set_info_type)                    :: l_set_info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'destroy_l_set_info_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    DEALLOCATE(l_set_info%ly_max,l_set_info%lz_max,l_set_info%map,&
               l_set_info%dpy,l_set_info%dpz,l_set_info%alpha,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE destroy_l_set_info_type

!****************************************************************************

  SUBROUTINE l_info_retain(l_info,error)

    TYPE(l_info_type), POINTER               :: l_info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'l_info_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(l_info),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CPPreconditionNoFail(l_info%ref_count>0,cp_failure_level,routineP,error)
      l_info%ref_count=l_info%ref_count+1
    END IF

  END SUBROUTINE l_info_retain

!****************************************************************************

  SUBROUTINE l_info_release(l_info,error)

    TYPE(l_info_type), POINTER               :: l_info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'l_info_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(l_info)) THEN
       CPPreconditionNoFail(l_info%ref_count>0,cp_failure_level,routineP,error)
       l_info%ref_count=l_info%ref_count-1
       IF (l_info%ref_count==0) THEN
          CALL destroy_l_info(l_info,error=error)
       END IF
    END IF

  END SUBROUTINE l_info_release

!****************************************************************************

END MODULE

!****************************************************************************
