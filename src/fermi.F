!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
!!****** cp2k/fermi [1.0] *
!!
!!   NAME
!!     fermi
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
!
    MODULE fermi
!
!------------------------------------------------------------------------------!

      USE kinds, ONLY : dbl

      USE termination, ONLY : stop_program
      USE physcon, ONLY : boltzmann, joule
!
      IMPLICIT NONE
!
      PRIVATE

      PUBLIC :: fermi_distribution_type, init_fermi_dist, fermi_info
!
      TYPE fermi_distribution_type
        REAL (dbl) :: electronic_temp
        REAL (dbl) :: betael
        REAL (dbl) :: mu
        INTEGER :: nel, nalpha, nbeta
        INTEGER :: spin_polarization
        INTEGER :: multiplicity
        INTEGER :: nstate, na, nb
      END TYPE fermi_distribution_type
!!*****
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!!****** fermi/init_fermi_dist [1.0] *
!!
!!   NAME
!!     init_fermi_dist
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
      SUBROUTINE init_fermi_dist(fdist,nel,etemp,lsd,mult,nstate)
        IMPLICIT NONE
        TYPE (fermi_distribution_type), INTENT (INOUT) :: fdist
        INTEGER, OPTIONAL, INTENT (IN) :: nel
        REAL (dbl), OPTIONAL, INTENT (IN) :: etemp
        LOGICAL, OPTIONAL, INTENT (IN) :: lsd
        INTEGER, OPTIONAL, INTENT (IN) :: mult
        INTEGER, OPTIONAL, INTENT (IN) :: nstate

        INTEGER :: m, isos

        fdist%nel = 0
        fdist%electronic_temp = 0._dbl
        fdist%spin_polarization = 0
        fdist%multiplicity = 1
        fdist%nstate = 0
        fdist%mu = 0._dbl
!..total number of electrons
        IF (present(nel)) THEN
          fdist%nel = nel
        END IF
!..electronic temperature
        IF (present(etemp)) THEN
          fdist%electronic_temp = etemp
        END IF
        IF (fdist%electronic_temp>0._dbl) THEN
          fdist%betael = joule/(boltzmann*fdist%electronic_temp)
        ELSE
          fdist%betael = 1.E33_dbl
        END IF
!..spin polarization
        IF (present(lsd)) THEN
          IF (lsd) fdist%spin_polarization = 1
        END IF
!..multiplicity
        IF (present(mult)) THEN
          IF (lsd) fdist%multiplicity = mult
        END IF
!..number of electronic states
        IF (present(nstate)) THEN
          fdist%nstate = nstate
        ELSE
          fdist%nstate = 0
        END IF
!..test the current setting and adjust defaults
!..number of electrons at least 1
        IF (fdist%nel<1) THEN
          CALL stop_program('INIT_FERMI_DIST','Number of electrons < 1')
        END IF
!..if lda only singlet states
        IF (fdist%spin_polarization==0 .AND. fdist%multiplicity/=1) THEN
          CALL stop_program('INIT_FERMI_DIST','LDA: multiplicity has to be 1')
        END IF
!..set number of alpha and beta electrons
        IF (fdist%spin_polarization==1) THEN
          m = fdist%multiplicity
          IF (mod(nel+m-1,2)/=0) THEN
            CALL stop_program('INIT_FERMI_DIST','NEL incons. with multiplicity')
          END IF
          fdist%nalpha = (nel+m-1)/2
          fdist%nbeta = nel - fdist%nalpha
        ELSE
          fdist%nalpha = nel
          fdist%nbeta = 0
        END IF
!..set number of states
        IF (fdist%nstate/=0) THEN
          IF (fdist%spin_polarization==1) THEN
            fdist%na = fdist%nalpha
            fdist%nb = fdist%nbeta
            m = fdist%nstate - (fdist%na+fdist%nb)
            IF (m<0) THEN
              CALL stop_program('INIT_FERMI_DIST','Not enough states')
            END IF
            fdist%na = fdist%na + (m+1)/2
            fdist%nb = fdist%nb + m/2
          ELSE
            fdist%na = fdist%nstate
            fdist%nb = 0
            IF (fdist%nstate<(fdist%nel+1)/2) THEN
              CALL stop_program('INIT_FERMI_DIST','Not enough states')
            END IF
          END IF
        ELSE
          IF (fdist%electronic_temp==0._dbl) THEN
            IF (fdist%spin_polarization==1) THEN
              fdist%na = fdist%nalpha
              fdist%nb = fdist%nbeta
              fdist%nstate = fdist%na + fdist%nb
            ELSE
              fdist%nstate = (fdist%nel+1)/2
              fdist%na = fdist%nstate
              fdist%nb = 0
            END IF
          ELSE
            CALL stop_program('INIT_FERMI_DIST','If etemp /= 0 =>'//&
              ' number of states has to be specified')
          END IF
        END IF
      END SUBROUTINE init_fermi_dist
!!*****
!------------------------------------------------------------------------------!
!!****** fermi/fermi_info [1.0] *
!!
!!   NAME
!!     fermi_info
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
      SUBROUTINE fermi_info(fdist,plevel,punit)
        IMPLICIT NONE
        TYPE (fermi_distribution_type), INTENT (IN) :: fdist
        INTEGER :: plevel, punit
        CHARACTER (len=10) :: mul(1:8)

        mul = '          '
        mul(1) = 'singlet'
        mul(2) = 'doublet'
        mul(3) = 'triplet'
        mul(4) = 'quartet'
        mul(5) = 'quintet'
        mul(6) = 'sextet'
        mul(7) = 'septet'
        mul(8) = 'octet'

        IF (plevel>0) THEN
          WRITE (punit,*)
          WRITE (punit,'(A,T71,I10)') ' FERMI| Total number of electrons ', &
            fdist%nel
          IF (fdist%spin_polarization==1) THEN
            WRITE (punit,'(A,T71,I10)') ' FERMI| Number of alpha electrons ', &
              fdist%nalpha
            WRITE (punit,'(A,T71,I10)') ' FERMI| Number of beta electrons ', &
              fdist%nbeta
          END IF
          WRITE (punit,'(A,T71,A)') ' FERMI| Multiplicity', &
            adjustr(mul(fdist%multiplicity))
          WRITE (punit,'(A,T71,F10.2)') ' FERMI| Electronic temperature [K]', &
            fdist%electronic_temp
          WRITE (punit,'(A,T71,I10)') &
            ' FERMI| Total number of electronic states ', fdist%nstate
          IF (fdist%spin_polarization==1) THEN
            WRITE (punit,'(A,T71,I10)') &
              ' FERMI| Number of alpha electronic states ', fdist%na
            WRITE (punit,'(A,T71,I10)') &
              ' FERMI| Number of beta electronic states ', fdist%nb
          END IF
        END IF
      END SUBROUTINE fermi_info
!!*****
!------------------------------------------------------------------------------!

    END MODULE fermi

!------------------------------------------------------------------------------!
