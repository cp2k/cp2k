!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
#include "cp_prep_globals.h"
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
!!     JGH (24.11.2001) : modifyed fermi_distribution_type
!!
!!   SOURCE
!******************************************************************************
!
    MODULE fermi
!
!------------------------------------------------------------------------------!
!
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dbl
  USE physcon,                         ONLY: boltzmann,&
                                             joule
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
      IMPLICIT NONE
!
      PRIVATE

      CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN="fermi"

      PUBLIC :: fermi_distribution_type, init_fermi_dist, fermi_info, &
           dealloc_fermi_dist
!
      TYPE fermi_distribution_type
        REAL (dbl) :: chemical_potential
        REAL (dbl) :: total_electrons, alpha_electrons, beta_electrons
        REAL (dbl), POINTER :: electrons(:)
        INTEGER :: spin_polarization
        INTEGER :: nstate, nalpha, nbeta
        INTEGER, POINTER :: states(:)
        LOGICAL :: broadening_perform
        CHARACTER ( LEN = 40 ) :: broadening_method
        REAL (dbl) :: electronic_temp
        REAL (dbl) :: betael
        LOGICAL :: fixedspinmoment
        REAL (dbl) :: fixedspinmoment_amount
        INTEGER :: multiplicity
        LOGICAL :: fixedoccupation
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
!!   SYNOPSIS
!!     Subroutine init_fermi_dist(fdist, total_electrons, electrons,&
!!         total_states, states, electronic_temp, broadening_method, spins,&
!!         spin_moment, fixed_occupation)
!!       Implicit None
!!       Type(fermi_distribution_type), Intent (INOUT):: fdist
!!       Real(Kind=dbl), Optional, Intent (IN):: total_electrons
!!       Real(Kind=dbl), Dimension(:), Optional, Intent (IN):: electrons
!!       Integer, Optional, Intent (IN):: total_states
!!       Integer, Dimension(:), Optional, Intent (IN):: states
!!       Real(Kind=dbl), Optional, Intent (IN):: electronic_temp
!!       Character(Len=*), Optional, Intent (IN):: broadening_method
!!       Integer, Optional, Intent (IN):: spins
!!       Real(Kind=dbl), Optional, Intent (IN):: spin_moment
!!       Logical, Optional, Intent (IN):: fixed_occupation
!!     End Subroutine init_fermi_dist
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

      SUBROUTINE init_fermi_dist(fdist,total_electrons,electrons,total_states,&
                 states,electronic_temp,broadening_method,spins,&
                 spin_moment,fixed_occupation)

        IMPLICIT NONE

!       Arguments
        TYPE (fermi_distribution_type), INTENT (INOUT) :: fdist
        REAL (dbl), OPTIONAL, INTENT (IN) :: total_electrons
        REAL (dbl), DIMENSION (:), OPTIONAL, INTENT (IN) :: electrons
        INTEGER, OPTIONAL, INTENT (IN) :: total_states
        INTEGER, DIMENSION (:), OPTIONAL, INTENT (IN) :: states
        REAL (dbl), OPTIONAL, INTENT (IN) :: electronic_temp
        CHARACTER (LEN=*), OPTIONAL, INTENT (IN) :: broadening_method
        INTEGER, OPTIONAL, INTENT (IN) :: spins
        REAL (dbl), OPTIONAL, INTENT (IN) :: spin_moment
        LOGICAL, OPTIONAL, INTENT (IN) :: fixed_occupation

!       Locals
        INTEGER :: m, isos
        REAL ( dbl ) :: e

!------------------------------------------------------------------------------!
!       set defaults
        fdist%chemical_potential = 0._dbl
        fdist%total_electrons = 0._dbl
        fdist%alpha_electrons = 0._dbl
        fdist%beta_electrons  = 0._dbl
        fdist%spin_polarization = 1
        fdist%nstate = 0
        fdist%nalpha = 0
        fdist%nbeta  = 0
        fdist%broadening_perform = .FALSE.
        fdist%broadening_method = "FERMI"
        fdist%electronic_temp = 0._dbl
        fdist%fixedspinmoment = .FALSE.
        fdist%fixedspinmoment_amount = 1.E33_dbl
        fdist%multiplicity = -1
        fdist%fixedoccupation = .FALSE.
        NULLIFY(fdist%electrons)
        NULLIFY(fdist%states)
!------------------------------------------------------------------------------!
!       total number of electrons
        IF (PRESENT(total_electrons)) THEN
          fdist%total_electrons = total_electrons
        END IF
!       electrons per set
        IF (PRESENT(electrons)) THEN
          m = SIZE(electrons)
          ALLOCATE (fdist%electrons(m),STAT=isos)
          IF (isos/=0) CALL stop_memory('init_fermi_dist','fermi%electrons',m)
          fdist%electrons = electrons
        END IF
!       number of electronic states
        IF (PRESENT(total_states)) THEN
          fdist%nstate = total_states
        END IF
!       states per set
        IF (PRESENT(states)) THEN
          m = SIZE(states)
          ALLOCATE (fdist%states(m),STAT=isos)
          IF (isos/=0) CALL stop_memory('init_fermi_dist','fermi%states',m)
          fdist%states = states
        END IF
!       electronic temperature
        IF (PRESENT(electronic_temp)) THEN
          fdist%electronic_temp = electronic_temp
        END IF
!       broadening method
        IF (PRESENT(broadening_method)) THEN
          fdist%broadening_method = broadening_method
        END IF
!       number of spins
        IF (PRESENT(spins)) THEN
          fdist%spin_polarization = spins
        END IF
!       fixed spin moment
        IF (PRESENT(spin_moment)) THEN
          fdist%fixedspinmoment_amount = spin_moment
        END IF
!       fixed occupation
        IF (PRESENT(fixed_occupation)) THEN
          fdist%fixedoccupation = fixed_occupation
        END IF
!------------------------------------------------------------------------------!
!       test the current setting and adjust defaults

!       electronic temperature
        IF (fdist%electronic_temp > 0._dbl) THEN
          fdist%betael = joule/(boltzmann*fdist%electronic_temp)
          fdist%broadening_perform = .TRUE.
        ELSE
          fdist%betael = 1.E33_dbl
          fdist%broadening_perform = .FALSE.
        END IF

!       test number of spin sets
        IF (fdist%spin_polarization < 1) THEN
          CALL stop_program('INIT_FERMI_DIST','Number of spins sets< 1')
        END IF
        IF (ASSOCIATED(fdist%electrons)) THEN
          m = SIZE(fdist%electrons)
          IF (m /= fdist%spin_polarization) THEN
            CALL stop_program('INIT_FERMI_DIST','Inconsistent number of spin sets')
          END IF
        ELSE
          m = fdist%spin_polarization
          ALLOCATE (fdist%electrons(m),STAT=isos)
          IF (isos/=0) CALL stop_memory('init_fermi_dist','fermi%electrons',m)
          fdist%electrons = 0._dbl
        END IF
        IF (ASSOCIATED(fdist%states)) THEN
          m = SIZE(fdist%states)
          IF (m /= fdist%spin_polarization) THEN
            CALL stop_program('INIT_FERMI_DIST','Inconsistent number of spin sets')
          END IF
        ELSE
          m = fdist%spin_polarization
          ALLOCATE (fdist%states(m),STAT=isos)
          IF (isos/=0) CALL stop_memory('init_fermi_dist','fermi%states',m)
          fdist%states = 0
        END IF

!       number of electrons must be positiv
        IF (fdist%total_electrons < 0._dbl) THEN
          CALL stop_program('INIT_FERMI_DIST','Number of electrons < 0')
        END IF
        e = SUM(fdist%electrons)
        IF ( e > 0._dbl ) THEN
          IF ( fdist%total_electrons > 0._dbl ) THEN
            IF ( ABS ( e - fdist%total_electrons ) > 1.E-12_dbl ) THEN
              CALL stop_program('INIT_FERMI_DIST','Inconsistent number of electrons')
            END IF
          END IF
        END IF

!       check number of states
        m = MAXVAL(fdist%states)
        IF ( m > 0 ) THEN
          IF ( fdist%nstate == 0 ) THEN
            fdist%nstate = m
          ELSEIF ( fdist%nstate < m ) THEN
            CALL stop_program('INIT_FERMI_DIST','Inconsistent number of states')
          END IF
        ELSE
          fdist%states = fdist%nstate
        END IF

!       multiplicity and spin moments
        IF (fdist%spin_polarization==1) THEN
          fdist%fixedspinmoment = .TRUE.
          IF (fdist%multiplicity > 1 ) THEN
            CALL stop_program('INIT_FERMI_DIST','LDA: multiplicity has to be 1')
          ELSE
            fdist%multiplicity = 1
          END IF
          IF (fdist%fixedspinmoment_amount /= 0._dbl .AND. &
              fdist%fixedspinmoment_amount < 1.E30_dbl ) THEN
            CALL stop_program('INIT_FERMI_DIST','LDA: spin moment has to be 0')
          ELSE
            fdist%fixedspinmoment_amount = 0._dbl
          END IF
        ELSE
          IF (fdist%fixedspinmoment_amount > 1.E30_dbl .AND. &
              fdist%multiplicity < 1) THEN
            fdist%fixedspinmoment = .FALSE.
          ELSE
            fdist%fixedspinmoment = .TRUE.
            IF (fdist%multiplicity > 0) THEN
              e = 0.5_dbl * (fdist%multiplicity - 1)
              IF ( fdist%fixedspinmoment_amount > 1.E30_dbl ) THEN
                fdist%fixedspinmoment_amount = e
              ELSEIF ( ABS(e-fdist%fixedspinmoment_amount) > 1.E-10 ) THEN
                CALL stop_program('INIT_FERMI_DIST','Inconsistent spin moments')
              END IF
            ELSEIF ( fdist%fixedspinmoment_amount < 1.E30_dbl ) THEN
              e = 2._dbl * fdist%fixedspinmoment_amount + 1._dbl
              fdist%multiplicity = 0
            END IF
          END IF
        END IF

!       number of alpha and beta electrons
        IF (fdist%spin_polarization==1) THEN
          IF ( fdist%electrons(1) == 0 ) THEN
             fdist%electrons(1) = fdist%total_electrons
          ELSEIF ( fdist%electrons(1) /= fdist%total_electrons ) THEN
             CALL stop_program('INIT_FERMI_DIST','Inconsistent number of electrons')
          END IF
          fdist%alpha_electrons = fdist%total_electrons
          fdist%beta_electrons = 0
        ELSEIF (fdist%spin_polarization==2) THEN
          e = SUM(fdist%electrons)
          IF ( e == 0._dbl ) THEN
            IF ( fdist%fixedspinmoment ) THEN
              m = fdist%multiplicity
              IF ( m == 0 ) THEN
                fdist%alpha_electrons = 0.5_dbl*(fdist%total_electrons - &
                                         fdist%fixedspinmoment_amount)
                fdist%beta_electrons = 0.5_dbl*(fdist%total_electrons + &
                                         fdist%fixedspinmoment_amount)
              ELSE
                fdist%alpha_electrons = (NINT(fdist%total_electrons)+m-1)/2
                fdist%beta_electrons = fdist%total_electrons - fdist%alpha_electrons
                IF (m /= (fdist%beta_electrons-fdist%alpha_electrons)) THEN
                  CALL stop_program('INIT_FERMI_DIST',&
                     'electrons inconsistent with multiplicity')
                END IF
              END IF
              fdist%electrons(1) = fdist%alpha_electrons
              fdist%electrons(2) = fdist%beta_electrons
            ELSE
              e = fdist%total_electrons/2._dbl
              fdist%alpha_electrons = e
              fdist%beta_electrons = e
              fdist%electrons(1) = e
              fdist%electrons(2) = e
            END IF
          ELSE
            fdist%alpha_electrons = fdist%electrons(1)
            fdist%beta_electrons = fdist%electrons(2)
            IF ( fdist%fixedspinmoment ) THEN
              e = fdist%beta_electrons - fdist%alpha_electrons
              IF ( ABS ( e - fdist%fixedspinmoment_amount ) > 1.E-12_dbl ) THEN
                CALL stop_program('INIT_FERMI_DIST','electrons and spin are inconsistent')
              END IF
            END IF
          END IF
        ELSE
          CALL stop_program('INIT_FERMI_DIST','not yet done')
        END IF

!..set number of states
        IF (fdist%nstate/=0) THEN
          IF (fdist%spin_polarization==1) THEN
            fdist%nalpha = fdist%nstate
            fdist%nbeta = 0
            IF (REAL(2*fdist%nstate,KIND=dbl) < fdist%total_electrons) THEN
              CALL stop_program('INIT_FERMI_DIST','Not enough states')
            END IF
          ELSEIF (fdist%spin_polarization==2) THEN
            fdist%nalpha = INT(fdist%alpha_electrons-0.1_dbl)+1
            fdist%nbeta = INT(fdist%beta_electrons-0.1_dbl)+1
            m = fdist%nstate - (fdist%nalpha+fdist%nbeta)
            IF (m<0) THEN
              CALL stop_program('INIT_FERMI_DIST','Not enough states')
            END IF
            fdist%nalpha = fdist%nalpha + m/2
            fdist%nbeta = fdist%nbeta + (m+1)/2
          ELSE
            CALL stop_program('INIT_FERMI_DIST','not yet done')
          END IF
        ELSE
          IF (fdist%electronic_temp==0._dbl) THEN
            IF (fdist%spin_polarization==1) THEN
              fdist%nstate = CEILING(0.5_dbl*fdist%total_electrons)
              fdist%nalpha = fdist%nstate
              fdist%nbeta = 0
            ELSEIF (fdist%spin_polarization==2) THEN
              fdist%nalpha = INT(fdist%alpha_electrons-0.1_dbl)+1
              fdist%nbeta = INT(fdist%beta_electrons-0.1_dbl)+1
              fdist%nstate = fdist%nalpha + fdist%nbeta
            ELSE
              CALL stop_program('INIT_FERMI_DIST','not yet done')
            END IF
          ELSE
            CALL stop_program('INIT_FERMI_DIST','If etemp /= 0 =>'//&
              ' number of states has to be specified')
          END IF
        END IF

!------------------------------------------------------------------------------!
      END SUBROUTINE init_fermi_dist
!------------------------------------------------------------------------------!

!!*****
!------------------------------------------------------------------------------!
!!****** fermi/fermi_info [1.0] *
!!
!!   NAME
!!     fermi_info
!!
!!   SYNOPSIS
!!     Subroutine fermi_info(fdist, plevel, punit)
!!       Implicit None
!!       Type(fermi_distribution_type), Intent (IN):: fdist
!!       Integer:: plevel, punit
!!     End Subroutine fermi_info
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

!       Arguments
        TYPE (fermi_distribution_type), INTENT (IN) :: fdist
        INTEGER :: plevel, punit

!       Local
        CHARACTER (len=10) :: mul(1:8)
        INTEGER :: i
!------------------------------------------------------------------------------!

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
          WRITE (punit,'(A,T71,F10.2)') ' FERMI| Total number of electrons ', &
            fdist%total_electrons
          IF (fdist%spin_polarization == 2) THEN
            WRITE (punit,'(A,T71,F10.2)') ' FERMI| Number of alpha electrons ', &
              fdist%alpha_electrons
            WRITE (punit,'(A,T71,F10.2)') ' FERMI| Number of beta electrons ', &
              fdist%beta_electrons
          ELSEIF (fdist%spin_polarization > 2) THEN
            DO i=1,fdist%spin_polarization
              WRITE (punit,'(A,T40,I3,T71,F10.2)') &
                   ' FERMI| Number of electrons in set', i, fdist%electrons(i)
            END DO
          END IF
          IF(fdist%fixedspinmoment) THEN
            IF(fdist%multiplicity > 0) THEN
              WRITE (punit,'(A,T71,A)') ' FERMI| Multiplicity', &
                ADJUSTR(mul(fdist%multiplicity))
            ELSE
              WRITE (punit,'(A,T71,F10.2)') ' FERMI| Spin moment', &
                fdist%fixedspinmoment_amount
            END IF
          END IF
          IF(fdist%broadening_perform) THEN
            WRITE (punit,'(A,T71,F10.2)') ' FERMI| Electronic temperature [K]', &
              fdist%electronic_temp
            WRITE (punit,'(A,T41,A)') ' FERMI| Broadening method', &
              ADJUSTR(fdist%broadening_method)
          END IF
          IF(fdist%fixedoccupation) THEN
            WRITE (punit,'(A,T71,I10)') &
              ' FERMI| Occupation of states is fixed '
          END IF
          WRITE (punit,'(A,T71,I10)') &
            ' FERMI| Total number of electronic states ', fdist%nstate
          IF (fdist%spin_polarization==2) THEN
            WRITE (punit,'(A,T71,I10)') &
              ' FERMI| Number of alpha electronic states ', fdist%nalpha
            WRITE (punit,'(A,T71,I10)') &
              ' FERMI| Number of beta electronic states ', fdist%nbeta
          ELSEIF (fdist%spin_polarization > 2) THEN
            DO i=1,fdist%spin_polarization
              WRITE (punit,'(A,T40,I3,T71,I10)') &
                   ' FERMI| Number of electronic states in set', i, fdist%states(i)
            END DO
          END IF
        END IF

!------------------------------------------------------------------------------!
      END SUBROUTINE fermi_info
!------------------------------------------------------------------------------!

!!*****

!!****f* fermi/dealloc_fermi_dist [1.0] *
!!
!!   NAME
!!     dealloc_fermi_dist
!!
!!   SYNOPSIS
!!     Subroutine dealloc_fermi_dist(fermi_dist, error)
!!       Type(fermi_distribution_type), Pointer:: fermi_dist
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine dealloc_fermi_dist
!!
!!   FUNCTION
!!     deallocates the fermi distribution
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - fermi_dist: the distribution to deallocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE dealloc_fermi_dist(fermi_dist,error)
  TYPE(fermi_distribution_type), POINTER :: fermi_dist
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='dealloc_fermi_dist',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  IF (ASSOCIATED(fermi_dist)) THEN
     IF (ASSOCIATED(fermi_dist%electrons)) THEN
        DEALLOCATE(fermi_dist%electrons, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
     IF (ASSOCIATED(fermi_dist%states)) THEN
        DEALLOCATE(fermi_dist%states, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
     DEALLOCATE(fermi_dist, stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF

END SUBROUTINE dealloc_fermi_dist
!***************************************************************************
!------------------------------------------------------------------------------!

    END MODULE fermi

!------------------------------------------------------------------------------!
