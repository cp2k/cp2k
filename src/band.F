
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/band [1.0] *
!!
!!   NAME
!!     band
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (based on CPMD code)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
!
    MODULE band
!
!------------------------------------------------------------------------------!
!
  USE brillouin,                       ONLY: kpoint_dealloc,&
                                             kpoint_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE fermi,                           ONLY: dealloc_fermi_dist,&
                                             fermi_distribution_type
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

      IMPLICIT NONE
!
      PRIVATE

      CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN="band"
!
      PUBLIC :: band_structure_type, init_band_structure, band_structure_info,&
           band_structure_preinit, band_structure_dealloc_ref
!
      TYPE band_structure_type
        TYPE (fermi_distribution_type), POINTER :: fd
        TYPE (kpoint_type), POINTER :: kpt
        REAL(KIND=dp), DIMENSION (:,:,:), POINTER :: ev
        REAL(KIND=dp), DIMENSION (:,:,:), POINTER :: oc
      END TYPE band_structure_type
!
!!*****
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!!****** band/band_structure_preinit [1.0] *
!!
!!   NAME
!!     band_structure_preinit
!!
!!   SYNOPSIS
!!     Subroutine band_structure_preinit(bs, error)
!!       Implicit None
!!       Type(band_structure_type), Intent (OUT):: bs
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine band_structure_preinit
!!
!!   FUNCTION
!!     initialization (to avoid leaks)
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
      SUBROUTINE band_structure_preinit(bs, error)
    TYPE(band_structure_type), INTENT(OUT)   :: bs
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

        NULLIFY(bs%fd, bs%kpt, bs%ev, bs%oc)
      END SUBROUTINE band_structure_preinit
!!*****
!------------------------------------------------------------------------------!
!!****** band/band_structure_dealloc_ref [1.0] *
!!
!!   NAME
!!     band_structure_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine band_structure_dealloc_ref(bs, error)
!!       Implicit None
!!       Type(band_structure_type), Intent (INOUT):: bs
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine band_structure_dealloc_ref
!!
!!   FUNCTION
!!     deallocates the memory allocated by this object
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
      SUBROUTINE band_structure_dealloc_ref(bs, error)
    TYPE(band_structure_type), INTENT(INOUT) :: bs
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "band_structure_dealloc_ref", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: stat

        CALL dealloc_fermi_dist(bs%fd)
        CALL kpoint_dealloc(bs%kpt)
        IF (ASSOCIATED(bs%ev)) THEN
           DEALLOCATE(bs%ev, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(bs%oc)) THEN
           DEALLOCATE(bs%oc, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
      END SUBROUTINE band_structure_dealloc_ref
!!*****
!------------------------------------------------------------------------------!
!!****** band/init_band_structure [1.0] *
!!
!!   NAME
!!     init_band_structure
!!
!!   SYNOPSIS
!!     Subroutine init_band_structure(bs, fd, kp)
!!       Implicit None
!!       Type(band_structure_type), Intent (INOUT):: bs
!!       Type(fermi_distribution_type), Intent (IN), Target:: fd
!!       Type(kpoint_type), Intent (IN), Target:: kp
!!     End Subroutine init_band_structure
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

      SUBROUTINE init_band_structure(bs,fd,kp)
    TYPE(band_structure_type), INTENT(INOUT) :: bs
    TYPE(fermi_distribution_type), &
      INTENT(IN), TARGET                     :: fd
    TYPE(kpoint_type), INTENT(IN), TARGET    :: kp

    INTEGER                                  :: isos, nk, ns, nspin

        bs%fd => fd
        bs%kpt => kp

        nspin = fd%spin_polarization
        ns = MAX(fd%nalpha,fd%nbeta)
        nk = kp%nkpt

        IF (ASSOCIATED(bs%ev)) THEN
           CPAssertNoErr(.FALSE.,cp_warning_level,"init_band_structure")
           NULLIFY (bs%ev) !FM memory leak?
        END IF
        IF (ASSOCIATED(bs%oc)) THEN
           CPAssertNoErr(.FALSE.,cp_warning_level,"init_band_structure")
           NULLIFY (bs%oc) !FM memory leak?
        END IF
        ALLOCATE (bs%ev(ns,nk,nspin),STAT=isos)
        IF (isos/=0) CALL stop_memory('init_band_structure', &
          'bs%ev',ns*nk*nspin)
        ALLOCATE (bs%oc(ns,nk,nspin),STAT=isos)
        IF (isos/=0) CALL stop_memory('init_band_structure', &
          'bs%oc',ns*nk*nspin)

      END SUBROUTINE init_band_structure
!!*****
!------------------------------------------------------------------------------!
!!****** band/band_structure_info [1.0] *
!!
!!   NAME
!!     band_structure_info
!!
!!   SYNOPSIS
!!     Subroutine band_structure_info(bs, punit)
!!       Implicit None
!!       Type(band_structure_type), Intent (IN):: bs
!!       Integer, Intent (IN):: punit
!!     End Subroutine band_structure_info
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
      SUBROUTINE band_structure_info(bs,punit)
    TYPE(band_structure_type), INTENT(IN)    :: bs
    INTEGER, INTENT(IN)                      :: punit

    INTEGER                                  :: i, ik, left, nk, nm
    LOGICAL                                  :: nokp

        WRITE (punit,'(/,1x,79("-"))')
        WRITE (punit,'(" -",77x,"-")')
        WRITE (punit,'(" -",25x,a,24x,"-")') ' B A N D  S T R U C T U R E '
        WRITE (punit,'(" -",77x,"-")')
        WRITE (punit,'(1x,79("-"))')
        WRITE (punit,'(" -",T23,a,f8.3,T80,"-")') &
          ' Chemical Potential [a.u.] = ', bs%fd%chemical_potential
        WRITE (punit,'(1x,79("-"))')
        nokp = bs%kpt%scheme=='GAMMA' .OR. bs%kpt%scheme=='NULL'
        IF (nokp) THEN
          nk = 1
        ELSE
          nk = bs%kpt%nkpt
        END IF
        DO ik = 1, nk
          IF (nokp) THEN
            WRITE (punit,'(1x,79("."))')
            WRITE (punit, &
              '(" .",a,i5,6x,a,f8.3,a,f8.3,a,f8.3,4x,a,f8.5,T80,".")') &
              ' K-point:', ik, ' x=', bs%kpt%xk(1,ik), ' y=', bs%kpt%xk(2,ik), &
              ' z=', bs%kpt%xk(3,ik), ' weight=', bs%kpt%weight(ik)
            WRITE (punit,'(1x,79("."))')
          END IF
          IF (bs%fd%spin_polarization==1) THEN
            WRITE (punit,'(A,A)') ' State    Occupation     Eigenvalue[au] ', &
              '  State    Occupation     Eigenvalue[au] '
            DO i = 1, bs%fd%nstate, 2
              left = bs%fd%nstate - i + 1
              IF (left>1) THEN
                WRITE (punit,'(I5,5x,F10.3,6x,F13.6,T42,I5,5x,F10.3,6x,F13.6)' &
                  ) i, bs%oc(i,ik,1), bs%ev(i,ik,1), i + 1, bs%oc(i+1,ik,1), &
                  bs%ev(i+1,ik,1)
              ELSE
                WRITE (punit,'(I5,5x,F10.3,6x,F13.6)') i, bs%oc(i,ik,1), &
                  bs%ev(i,ik,1)
              END IF
            END DO
          ELSEIF (bs%fd%spin_polarization==2) THEN
            WRITE (punit,'(T12,A,T54,A)') ' Alpha Electrons ', &
              ' Beta Electrons'
            WRITE (punit,'(A,A)') ' State    Occupation     Eigenvalue[au] ', &
              '  State    Occupation     Eigenvalue[au] '
            nm = MIN(bs%fd%nalpha,bs%fd%nbeta)
            DO i = 1, nm
              WRITE (punit,'(I5,5x,F10.3,6x,F13.6,T42,I5,5x,F10.3,6x,F13.6)') &
                i, bs%oc(i,ik,1), bs%ev(i,ik,1), i, bs%oc(i,ik,2), &
                bs%ev(i,ik,2)
            END DO
            DO i = nm + 1, bs%fd%nalpha
              WRITE (punit,'(I5,5x,F10.3,6x,F13.6)') i, bs%oc(i,ik,1), &
                bs%ev(i,ik,1)
            END DO
            DO i = nm + 1, bs%fd%nbeta
              WRITE (punit,'(T40,I5,5x,F10.3,6x,F13.6)') i, bs%oc(i,ik,2), &
                bs%ev(i,ik,2)
            END DO
          ELSE
            CALL stop_program('band_structure_info','not programmed')
          END IF
        END DO
        WRITE (punit,'(1x,79("-"),/)')

      END SUBROUTINE band_structure_info
!!*****
!------------------------------------------------------------------------------!
    END MODULE band
!------------------------------------------------------------------------------!
