!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
! 
    MODULE band
! 
!------------------------------------------------------------------------------!
      USE kinds, ONLY : dbl
      USE stop_program, ONLY : stop_prg, stop_memory
      USE brillouin, ONLY : kpoint_type
      USE fermi, ONLY : fermi_distribution_type
!
      IMPLICIT NONE
! 
      PRIVATE
!
      PUBLIC :: band_structure_type, init_band_structure, band_structure_info, &
        adjust_mu, occupation
! 
      TYPE band_structure_type
        TYPE (fermi_distribution_type), POINTER :: fd
        TYPE (kpoint_type), POINTER :: kpt
        REAL (dbl), DIMENSION (:,:,:), POINTER :: ev
        REAL (dbl), DIMENSION (:,:,:), POINTER :: oc
      END TYPE band_structure_type
!
      REAL (dbl), PARAMETER :: precpar = 1.E-10_dbl, &
        delta = 1.E-6_dbl/27.212_dbl
      REAL (dbl), PARAMETER :: toll = 1.E-13_dbl
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
      SUBROUTINE init_band_structure(bs,fd,kp)
        IMPLICIT NONE
        TYPE (band_structure_type), INTENT (INOUT) :: bs
        TYPE (fermi_distribution_type), INTENT (IN), TARGET :: fd
        TYPE (kpoint_type), INTENT (IN), TARGET :: kp

        INTEGER :: nspin, ns, nk, isos

        bs%fd => fd
        bs%kpt => kp

        nspin = fd%spin_polarization + 1
        ns = max(fd%na,fd%nb)
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
!------------------------------------------------------------------------------!
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
          ' Chemical Potential [a.u.] = ', bs%fd%mu
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
          IF (bs%fd%spin_polarization==0) THEN
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
          ELSE
            WRITE (punit,'(T12,A,T54,A)') ' Alpha Electrons ', &
              ' Beta Electrons'
            WRITE (punit,'(A,A)') ' State    Occupation     Eigenvalue[au] ', &
              '  State    Occupation     Eigenvalue[au] '
            nm = min(bs%fd%na,bs%fd%nb)
            DO i = 1, nm
              WRITE (punit,'(I5,5x,F10.3,6x,F13.6,T42,I5,5x,F10.3,6x,F13.6)') &
                i, bs%oc(i,ik,1), bs%ev(i,ik,1), i, bs%oc(i,ik,2), &
                bs%ev(i,ik,2)
            END DO
            DO i = nm + 1, bs%fd%na
              WRITE (punit,'(I5,5x,F10.3,6x,F13.6)') i, bs%oc(i,ik,1), &
                bs%ev(i,ik,1)
            END DO
            DO i = nm + 1, bs%fd%nb
              WRITE (punit,'(T40,I5,5x,F10.3,6x,F13.6)') i, bs%oc(i,ik,2), &
                bs%ev(i,ik,2)
            END DO
          END IF
        END DO
        WRITE (punit,'(1x,79("-"),/)')

      END SUBROUTINE band_structure_info
!------------------------------------------------------------------------------!
      SUBROUTINE adjust_mu(bs)
!
! adjust the chemical potential with the bisection method
! eigenvalues do not have to be sorted
! routine adapted from CPMD
!
        IMPLICIT NONE
        TYPE (band_structure_type), INTENT (INOUT) :: bs
        REAL (dbl) :: amu, amu1, amu2, damu, sdeg
        REAL (dbl) :: rhint, rhint1, rhint2, precision
        INTEGER, PARAMETER :: it_max = 200
        INTEGER :: nk, it

!
        precision = precpar*log(float(bs%fd%nel+1))
!..min and max eigenvalue
        amu1 = minval(bs%ev)
        amu2 = maxval(bs%ev)
!
        DO
          rhint1 = rhoint(bs,amu1)
          IF (rhint1<0._dbl) EXIT
          amu1 = -2._dbl*abs(amu1)
        END DO
        rhint2 = rhoint(bs,amu2)
        IF (rhint2<0._dbl) CALL stop_prg('ADJUST_MU', &
          'failed to find chem. potential','number of states is too small')
        damu = amu2 - amu1
        amu = amu1
!
        DO it = 1, it_max
          amu = amu1 + 0.5_dbl*damu
          rhint = rhoint(bs,amu)
          IF (damu<delta) THEN
            IF (abs(rhint)<precision) THEN
              EXIT
            ELSE IF (abs(rhint1)<precision) THEN
              amu = amu1
              EXIT
            ELSE IF (abs(rhint1)<precision) THEN
              amu = amu2
              EXIT
            END IF
          END IF
          IF (rhint<0._dbl) THEN
            amu1 = amu
            rhint1 = rhint
          ELSE
            amu2 = amu
            rhint2 = rhint
          END IF
          damu = amu2 - amu1
        END DO
        IF (damu>delta) THEN
          CALL stop_prg('ADJUST_MU','failed to find chem. potential')
        END IF
        bs%fd%mu = amu
!
      END SUBROUTINE adjust_mu
!------------------------------------------------------------------------------!
      FUNCTION rhoint(bs,mu)
        IMPLICIT NONE
        TYPE (band_structure_type), INTENT (IN) :: bs
        REAL (dbl), INTENT (IN) :: mu
        REAL (dbl) :: rhoint, xlm, arg, argmin, argmax, sdeg, remain
        REAL (dbl) :: dramu, drwe, wbigwe
        INTEGER :: is, ik, i, n

        argmin = log(toll)
        argmax = -argmin
        IF (bs%fd%spin_polarization==0) THEN
          sdeg = 2._dbl
        ELSE
          sdeg = 1._dbl
        END IF
        dramu = dround(mu,delta)
        rhoint = 0._dbl
        wbigwe = 0._dbl
        DO is = 1, size(bs%ev(1,1,:))
          DO ik = 1, size(bs%ev(1,:,1))
            IF (is==1) n = bs%fd%na
            IF (is==2) n = bs%fd%nb
            DO i = 1, n
              arg = -bs%fd%betael*(bs%ev(i,ik,is)-mu)
              IF (arg>argmax) THEN
                xlm = 1._dbl
              ELSE IF (arg<argmin) THEN
                xlm = 0._dbl
              ELSE
                xlm = exp(arg)/(exp(arg)+1._dbl)
              END IF
              drwe = dround(bs%ev(i,ik,is),delta)
              IF (abs(drwe-dramu)<delta) THEN
                wbigwe = wbigwe + bs%kpt%weight(ik)
              ELSE
                rhoint = rhoint + bs%kpt%weight(ik)*xlm
              END IF
            END DO
          END DO
        END DO
        remain = float(bs%fd%nel) - rhoint
        IF (remain<0._dbl) remain = sdeg
        IF (wbigwe>0._dbl) THEN
          remain = remain/wbigwe
          remain = min(remain,sdeg)
          DO is = 1, size(bs%ev(1,1,:))
            DO ik = 1, size(bs%ev(1,:,1))
              IF (is==1) n = bs%fd%na
              IF (is==2) n = bs%fd%nb
              DO i = 1, n
                drwe = dround(bs%ev(i,ik,is),delta)
                IF (abs(drwe-dramu)<delta) rhoint = rhoint + &
                  bs%kpt%weight(ik)*remain
              END DO
            END DO
          END DO
        END IF
        rhoint = sdeg*rhoint - float(bs%fd%nel)
      END FUNCTION rhoint
!------------------------------------------------------------------------------!
      SUBROUTINE occupation(bs)
        IMPLICIT NONE
        TYPE (band_structure_type), INTENT (INOUT) :: bs
        REAL (dbl) :: xlm, arg, argmin, argmax, sdeg
        REAL (dbl) :: dramu, drwe, wbigwe, remain, rhotot
        INTEGER :: is, ik, i, n

        argmin = log(toll)
        argmax = -argmin
        IF (bs%fd%spin_polarization==0) THEN
          sdeg = 2._dbl
        ELSE
          sdeg = 1._dbl
        END IF
        dramu = dround(bs%fd%mu,delta)
        wbigwe = 0._dbl
        rhotot = 0._dbl
        DO is = 1, size(bs%oc(1,1,:))
          DO ik = 1, size(bs%oc(1,:,1))
            IF (is==1) n = bs%fd%na
            IF (is==2) n = bs%fd%nb
            DO i = 1, n
              arg = -bs%fd%betael*(bs%ev(i,ik,is)-bs%fd%mu)
              IF (arg>argmax) THEN
                xlm = 1._dbl
              ELSE IF (arg<argmin) THEN
                xlm = 0._dbl
              ELSE
                xlm = exp(arg)/(exp(arg)+1._dbl)
              END IF
              drwe = dround(bs%ev(i,ik,is),delta)
              IF (abs(drwe-dramu)<delta) THEN
                wbigwe = wbigwe + bs%kpt%weight(ik)
              ELSE
                rhotot = rhotot + bs%kpt%weight(ik)*xlm
                bs%oc(i,ik,is) = xlm*sdeg
              END IF
            END DO
          END DO
        END DO
        remain = float(bs%fd%nel) - rhotot
        IF (remain<0._dbl) remain = sdeg
        IF (wbigwe>0._dbl) THEN
          remain = remain/wbigwe
          remain = min(remain,sdeg)
          DO is = 1, size(bs%ev(1,1,:))
            DO ik = 1, size(bs%ev(1,:,1))
              IF (is==1) n = bs%fd%na
              IF (is==2) n = bs%fd%nb
              DO i = 1, n
                drwe = dround(bs%ev(i,ik,is),delta)
                IF (abs(drwe-dramu)<delta) bs%oc(i,ik,is) = remain
              END DO
            END DO
          END DO
        END IF
      END SUBROUTINE occupation
!------------------------------------------------------------------------------!
      FUNCTION dround(a,delta)
        IMPLICIT NONE
        REAL (dbl), INTENT (IN) :: a, delta
        REAL (dbl) :: dround

        dround = a - mod(a,delta)
      END FUNCTION dround
!------------------------------------------------------------------------------!
    END MODULE band
!------------------------------------------------------------------------------!
