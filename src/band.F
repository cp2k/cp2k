!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
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

      USE kinds, ONLY : dbl

      USE brillouin, ONLY : kpoint_type
      USE fermi, ONLY : fermi_distribution_type
      USE termination, ONLY : stop_memory, stop_program
!
      IMPLICIT NONE
!
      PRIVATE
!
      PUBLIC :: band_structure_type, init_band_structure, band_structure_info
!
      TYPE band_structure_type
        TYPE (fermi_distribution_type), POINTER :: fd
        TYPE (kpoint_type), POINTER :: kpt
        REAL (dbl), DIMENSION (:,:,:), POINTER :: ev
        REAL (dbl), DIMENSION (:,:,:), POINTER :: oc
      END TYPE band_structure_type
!
!!*****
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!!****** band/init_band_structure [1.0] *
!!
!!   NAME
!!     init_band_structure
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
        IMPLICIT NONE
        TYPE (band_structure_type), INTENT (INOUT) :: bs
        TYPE (fermi_distribution_type), INTENT (IN), TARGET :: fd
        TYPE (kpoint_type), INTENT (IN), TARGET :: kp

        INTEGER :: nspin, ns, nk, isos

        bs%fd => fd
        bs%kpt => kp

        nspin = fd%spin_polarization
        ns = max(fd%nalpha,fd%nbeta)
        nk = kp%nkpt

        IF (associated(bs%ev)) NULLIFY (bs%ev)
        IF (associated(bs%oc)) NULLIFY (bs%oc)
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
        IMPLICIT NONE
        TYPE (band_structure_type), INTENT (IN) :: bs
        INTEGER, INTENT (IN) :: punit
        INTEGER :: i, ik, left, nm, nk
        LOGICAL :: nokp

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
            nm = min(bs%fd%nalpha,bs%fd%nbeta)
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
