!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/xc_exchange_gga [1.0] *
!!
!!   NAME
!!     xc_exchange_gga
!!
!!   FUNCTION
!!     Calculate several different exchange energy functionals
!!     with a GGA form 
!!
!!   AUTHOR
!!     JGH (27.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!     fawzi (04.2004)  : adapted to the new xc interface
!!
!!*** *************************************************************************

MODULE xc_exchange_gga_new
! *****************************************************************************
  USE cp_array_r_utils,                ONLY: cp_3d_r_p_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_functionals_utilities,        ONLY: calc_wave_vector,&
                                             set_util
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
                                             xc_rho_set_get,&
                                             xc_rho_set_type

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: xgga_info, xgga_eval
! *** Global parameters ***

  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13, &
                          f53 = 5._dbl*f13

  INTEGER, PARAMETER, PUBLIC :: xgga_b88=1, xgga_pw86=2, xgga_pw91=3, xgga_pbe=4,&
       xgga_revpbe=5, xgga_opt=6

  REAL(dbl) :: cx, flda, flsd, sfac, t13
  REAL(dbl) :: fact, tact
  REAL(dbl) :: eps_rho
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="xc_exchange_gga_new"  

! *****************************************************************************

CONTAINS
!!****f* xc_exchange_gga/xgga_info [1.0] *
!!
!!   NAME
!!     xgga_info
!!
!!   SYNOPSIS
!!     Subroutine xgga_info(functional, lsd, reference, shortform, needs,&
!!         max_deriv, error)
!!       Integer, Intent (IN):: functional
!!       Logical, Intent (IN):: lsd
!!       Character(Len=*), Intent (OUT), Optional:: reference, shortform
!!       Type(xc_rho_cflags_type), Intent (INOUT), Optional:: needs
!!       Integer, Intent (OUT), Optional:: max_deriv
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine xgga_info
!!
!!   FUNCTION
!!     return various information on the xgga functionals
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - functional: integer selecting the xgga functional, it should be one of
!!       the constants defined in this module: xgga_b88, xgga_pw86,...
!!     - lsd: a logical that specifies if you are asking about the lsd or lda
!!       version of the functional
!!     - reference: string with the reference of the actual functional
!!     - shortform: string with the shortform of the functional name
!!     - needs: the components needed by this functional are set to
!!       true (does not set the unneeded components to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE xgga_info(functional,lsd,reference,shortform, needs, max_deriv, error)
    INTEGER, INTENT(in)                      :: functional
    LOGICAL, INTENT(in)                      :: lsd
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xgga_info', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( reference ) ) THEN
       SELECT CASE(functional)
       CASE(xgga_b88)
          reference = "A. Becke, Phys. Rev. A 38, 3098 (1988)"
       CASE(xgga_pw86)
          reference = "J.P. Perdew and Y. Wang, Phys. Rev. B, 33, 8800 (1986)"
       CASE(xgga_pw91)
          reference = "J.P. Perdew et al., Phys. Rev. B, 46, 6671 (1992)"
       CASE(xgga_pbe)
          reference = "J.P. Perdew, K. Burke, M Ernzerhof, Phys. Rev. Lett, 77, 3865 (1996)"
       CASE(xgga_revpbe)
          reference = "Y. Zang et al., PRL, 80, 890 (1998) (Revised PBEX)"
       CASE(xgga_opt)
          reference = "Wee-Meng Hoe, A.J. Cohen, N.C. Handy, CPL, 341, 319 (2001)"
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "Invalid functional requested ("//cp_to_string(functional)//")",&
               error)
       END SELECT
       IF (.not.lsd) THEN
          IF (LEN_TRIM(reference)+6<LEN(reference)) THEN
             reference(LEN_TRIM(reference):LEN_TRIM(reference)+6)=' {LDA}'
          END IF
       END IF
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       SELECT CASE(functional)
       CASE(xgga_b88)
          shortform = "Becke 1988 Exchange Functional"
       CASE(xgga_pw86)
          shortform = "Perdew-Wang 1986 Functional (exchange energy)"
       CASE(xgga_pw91)
          shortform = "Perdew-Wang 1991 Functional (exchange energy)"
       CASE(xgga_pbe)
          shortform = "PBE exchange energy functional"
       CASE(xgga_revpbe)
          shortform = "Revised PBEX by Zang et al."
       CASE(xgga_opt)
          shortform = "OPTX exchange energy functional"
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "Invalid functional requested ("//cp_to_string(functional)//")",&
               error)
       END SELECT
       IF (.not.lsd) THEN
          IF (LEN_TRIM(shortform)+6<LEN(shortform)) THEN
             shortform(LEN_TRIM(shortform):LEN_TRIM(shortform)+6)=' {LDA}'
          END IF
       END IF
    END IF
    IF (PRESENT(needs)) THEN
       IF (lsd) THEN
          needs%rho_spin=.TRUE.
          needs%rho_spin_1_3=.TRUE.
          needs%norm_drho_spin=.TRUE.
       ELSE
          needs%rho=.TRUE.
          needs%rho_1_3=.TRUE.
          needs%norm_drho=.TRUE.
       END IF
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE xgga_info
  !***************************************************************************

!!****f* xc_exchange_gga/xgga_eval [1.0] *
!!
!!   NAME
!!     xgga_eval
!!
!!   SYNOPSIS
!!     Subroutine xgga_eval(functional, lsd, rho_set, deriv_set, order,&
!!         error)
!!       Integer, Intent (IN):: functional
!!       Logical, Intent (IN):: lsd
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Type(xc_derivative_set_type), Pointer:: deriv_set
!!       Integer, Intent (IN):: order
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine xgga_eval
!!
!!   FUNCTION
!!     evaluates different exchange gga
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - functional: integer to select the functional that should be evaluated
!!     - lsd: if the lsd version of the functional should be used
!!     - rho_set: the density where you want to evaluate the functional
!!     - deriv_set: place where to store the functional derivatives (they are
!!       added to the derivatives)
!!     - grad_deriv: degree of the derivative that should be evalated,
!!       if positive all the derivatives up to the given degree are evaluated,
!!       if negative only the given degree is calculated
!!     - iparset: the parameter set that should be used (93,120,147,407)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE xgga_eval(functional,lsd,rho_set,deriv_set,order,error)
    INTEGER, INTENT(in)                      :: functional
    LOGICAL, INTENT(in)                      :: lsd
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: order
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xgga_eval', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=12), DIMENSION(2)          :: norm_drho_spin_name
    CHARACTER(len=6), DIMENSION(2)           :: rho_spin_name
    INTEGER                                  :: handle, ispin, m, npoints, &
                                                nspin, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(dbl), DIMENSION(:, :, :), POINTER :: e_0, e_ndrho, e_ndrho_ndrho, &
      e_ndrho_ndrho_ndrho, e_rho, e_rho_ndrho, e_rho_ndrho_ndrho, e_rho_rho, &
      e_rho_rho_ndrho, e_rho_rho_rho
    REAL(kind=dbl)                           :: drho_cutoff, rho_cutoff
    REAL(kind=dbl), ALLOCATABLE, &
      DIMENSION(:)                           :: s
    REAL(kind=dbl), ALLOCATABLE, &
      DIMENSION(:, :)                        :: fs
    TYPE(cp_3d_r_p_type), DIMENSION(2)       :: norm_drho, rho, rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(bo,e_0,e_ndrho, e_ndrho_ndrho, e_ndrho_ndrho_ndrho, e_rho_ndrho_ndrho, &
         e_rho_ndrho, e_rho_rho_ndrho, e_rho, e_rho_rho, e_rho_rho_rho)
    DO ispin=1,2
       NULLIFY(norm_drho(ispin)%array, rho(ispin)%array, rho_1_3(ispin)%array)
    END DO

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (lsd) THEN
          CALL xc_rho_set_get(rho_set,rhoa_1_3=rho_1_3(1)%array,&
               rhob_1_3=rho_1_3(2)%array,rhoa=rho(1)%array,&
               rhob=rho(2)%array,norm_drhoa=norm_drho(1)%array, &
               norm_drhob=norm_drho(2)%array,rho_cutoff=rho_cutoff,&
               drho_cutoff=drho_cutoff, local_bounds=bo, error=error)
          nspin=2
          rho_spin_name=(/"(rhoa)","(rhob)"/)
          norm_drho_spin_name=(/"(norm_drhoa)","(norm_drhob)"/)
       ELSE
          CALL xc_rho_set_get(rho_set,rho=rho(1)%array,rho_1_3=rho_1_3(1)%array,&
               norm_drho=norm_drho(1)%array,local_bounds=bo,rho_cutoff=rho_cutoff,&
               drho_cutoff=drho_cutoff,error=error)
          nspin=1
          rho_spin_name=(/"(rho) ","(---) "/)
          norm_drho_spin_name=(/"(norm_drho) ","(----_----) "/)
       END IF
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)
       m = ABS(order)
       CALL xgga_init(rho_cutoff)

       ALLOCATE ( s(npoints), STAT=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE ( fs(npoints,m+1), STAT=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       DO ispin=1,nspin
          IF (lsd) THEN
             fact = flsd
             tact = 1._dbl
             CALL calc_wave_vector ( "p", rho(ispin)%array, norm_drho(ispin)%array, s )
          ELSE
             fact = flda
             tact = t13
             CALL calc_wave_vector ( "u", rho(ispin)%array, &
                  norm_drho(ispin)%array, s )
          END IF

          SELECT CASE (functional)
          CASE (xgga_b88)
             CALL efactor_b88 ( s, fs, m )
          CASE (xgga_pw86)
             CALL efactor_pw86 ( s, fs, m )
          CASE (xgga_pw91)
             CALL efactor_pw91 ( s, fs, m )
          CASE (xgga_pbe)
             tact = t13
             CALL efactor_pbex ( s, fs, m, 1 )
             IF (lsd) tact = 1._dbl
          CASE (xgga_revpbe)
             tact = t13
             CALL efactor_pbex ( s, fs, m, 2 )
             IF (lsd) tact = 1._dbl
          CASE (xgga_opt)
             CALL efactor_optx ( s, fs, m )
          CASE DEFAULT
             CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT

          IF ( order>=0 ) THEN
             deriv => xc_dset_get_derivative(deriv_set,"",&
                  allocate_deriv=.TRUE., error=error)
             CALL xc_derivative_get(deriv, deriv_data=e_0,error=error)

             CALL x_p_0 ( rho(ispin)%array, rho_1_3(ispin)%array, fs, e_0,&
                  npoints, error)
          END IF
          IF ( order>=1 .OR. order==-1 ) THEN
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin),&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,norm_drho_spin_name(ispin),&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho,error=error)

             CALL x_p_1 ( rho(ispin)%array, norm_drho(ispin)%array, &
                  rho_1_3(ispin)%array, s, fs, e_rho, e_ndrho, npoints, error )
          END IF
          IF ( order>=2 .OR. order==-2 ) THEN
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin)//&
                  rho_spin_name(ispin),allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin)//&
                  norm_drho_spin_name(ispin),allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_ndrho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,norm_drho_spin_name(ispin)//&
                  norm_drho_spin_name(ispin), allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho,error=error)

             CALL x_p_2 ( rho(ispin)%array, norm_drho(ispin)%array, &
                  rho_1_3(ispin)%array, s, fs, e_rho_rho, e_rho_ndrho,&
                  e_ndrho_ndrho, npoints, error )
          END IF
          IF ( order>=3 .OR. order==-3 ) THEN
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin)//&
                  rho_spin_name(ispin)//rho_spin_name(ispin),&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin)//&
                  rho_spin_name(ispin)//norm_drho_spin_name(ispin),&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_ndrho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin)//&
                  norm_drho_spin_name(ispin)//norm_drho_spin_name(ispin), &
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_ndrho_ndrho,error=error)
             deriv => xc_dset_get_derivative(deriv_set,norm_drho_spin_name(ispin)//&
                  norm_drho_spin_name(ispin)//norm_drho_spin_name(ispin),&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho,error=error)

             CALL x_p_3 ( rho(ispin)%array, norm_drho(ispin)%array,&
                  rho_1_3(ispin)%array, s, fs, e_rho_rho_rho, &
                  e_rho_rho_ndrho, e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho,&
                  npoints, error )
          END IF
          IF ( order>3.OR.order<-3) THEN
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="derivatives bigger than 3 not implemented", &
                  error=error, error_level=cp_failure_level)
          END IF
       END DO

       DEALLOCATE ( s, STAT=stat )
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       DEALLOCATE ( fs, STAT=stat )
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    END IF
    CALL timestop(0.0_dbl,handle)
  END SUBROUTINE xgga_eval
  !***************************************************************************


! *****************************************************************************

  SUBROUTINE xgga_init ( cutoff )

    REAL(dbl), INTENT(IN)                    :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

    cx = -0.75_dbl*(3._dbl/pi)**f13
    t13 = 2._dbl**f13
    flda = cx
    flsd = cx * t13

    sfac = 1._dbl / (2._dbl*(3._dbl*pi*pi)**f13)

  END SUBROUTINE xgga_init

! *****************************************************************************

  SUBROUTINE x_p_0 ( rho, r13, fs, e_0, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, r13
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(*), INTENT(INOUT)   :: e_0
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip

!   ---------------------------------------------------------------------------

!$omp parallel do private(ip)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         e_0(ip) = e_0(ip) + fact*r13(ip)*rho(ip) * fs(ip,1)
      END IF

    END DO

  END SUBROUTINE x_p_0

! *****************************************************************************

  SUBROUTINE x_p_1 ( rho, grho, r13, s, fs, e_rho,e_ndrho,npoints,error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(*), INTENT(INOUT)   :: e_rho, e_ndrho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: a0, a1, sx, sy

!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,a0,a1,sx,sy)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*rho(ip)
         a1 = f43*fact*r13(ip)
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         e_rho(ip) = e_rho(ip) + a1*fs(ip,1) + a0*fs(ip,2)*sx
         e_ndrho(ip) = e_ndrho(ip) + a0*fs(ip,2)*sy

      END IF

    END DO

  END SUBROUTINE x_p_1

! *****************************************************************************

  SUBROUTINE x_p_2 ( rho, grho, r13, s, fs, e_rho_rho, e_rho_ndrho, &
       e_ndrho_ndrho, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(*), INTENT(INOUT)   :: e_rho_rho, e_rho_ndrho, &
                                                e_ndrho_ndrho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: a0, a1, a2, sx, sxx, sxy, sy

!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,a0,a1,a2,sx,sy,sxx,sxy)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*rho(ip)
         a1 = f43*fact*r13(ip)
         a2 = f13*f43*fact/(r13(ip)*r13(ip))
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         sxx= 28._dbl/9._dbl*s(ip)/(rho(ip)*rho(ip))
         sxy= -f43*sfac*tact/(r13(ip)*rho(ip)*rho(ip))
         e_rho_rho(ip) = e_rho_rho(ip) + a2*fs(ip,1) +2._dbl*a1*fs(ip,2)*sx + &
                     a0*fs(ip,3)*sx*sx + a0*fs(ip,2)*sxx
         e_rho_ndrho(ip) = e_rho_ndrho(ip) &
              + a1*fs(ip,2)*sy + a0*fs(ip,3)*sx*sy + a0*fs(ip,2)*sxy
         e_ndrho_ndrho(ip) = e_ndrho_ndrho(ip) + a0*fs(ip,3)*sy*sy

      END IF

    END DO

  END SUBROUTINE x_p_2

! *****************************************************************************

  SUBROUTINE x_p_3 ( rho, grho, r13, s, fs, e_rho_rho_rho, e_rho_rho_ndrho,&
       e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho, npoints,error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, grho, r13, s
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: fs
    REAL(dbl), DIMENSION(*), INTENT(INOUT)   :: e_rho_rho_rho, &
                                                e_rho_rho_ndrho, &
                                                e_rho_ndrho_ndrho, &
                                                e_ndrho_ndrho_ndrho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: a0, a1, a2, a3, sx, sxx, &
                                                sxxx, sxxy, sxy, sy

!   ---------------------------------------------------------------------------

!$omp parallel do private(ip,a0,a1,a2,a3,sx,sy,sxx,sxy,sxxx,sxxy)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*rho(ip)
         a1 = f43*fact*r13(ip)
         a2 = f13*f43*fact/(r13(ip)*r13(ip))
         a3 = -f23*f13*f43*fact/(r13(ip)*r13(ip)*rho(ip))
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         sxx= 28._dbl/9._dbl*s(ip)/(rho(ip)*rho(ip))
         sxy= -f43*sfac*tact/(r13(ip)*rho(ip)*rho(ip))
         sxxx= -280._dbl/27._dbl*s(ip)/(rho(ip)*rho(ip)*rho(ip))
         sxxy= 28._dbl/9._dbl*sfac*tact/(r13(ip)*rho(ip)*rho(ip)*rho(ip))
         e_rho_rho_rho(ip) = e_rho_rho_rho(ip)&
                     + a3*fs(ip,1) + 3._dbl*a2*fs(ip,2)*sx + &
                     3._dbl*a1*fs(ip,3)*sx*sx + 3._dbl*a1*fs(ip,2)*sxx + &
                     a0*fs(ip,4)*sx*sx*sx + 3._dbl*a0*fs(ip,3)*sx*sxx + &
                     a0*fs(ip,2)*sxxx
         e_rho_rho_ndrho(ip) = e_rho_rho_ndrho(ip)&
                     + a2*fs(ip,2)*sy + 2._dbl*a1*fs(ip,3)*sx*sy + &
                     2._dbl*a1*fs(ip,2)*sxy + a0*fs(ip,4)*sx*sx*sy + &
                     2._dbl*a0*fs(ip,3)*sx*sxy + a0*fs(ip,3)*sxx*sy + &
                     a0*fs(ip,2)*sxxy
         e_rho_ndrho_ndrho(ip) = e_rho_ndrho_ndrho(ip)&
                     + a1*fs(ip,3)*sy*sy + a0*fs(ip,4)*sx*sy*sy + &
                     2._dbl*a0*fs(ip,3)*sxy*sy
         e_ndrho_ndrho_ndrho(ip) = e_ndrho_ndrho_ndrho(ip)&
                     + a0*fs(ip,4)*sy*sy*sy

      END IF

    END DO

  END SUBROUTINE x_p_3

! *****************************************************************************
! Enhancement Factors
! *****************************************************************************
  SUBROUTINE efactor_b88 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl) :: as, asp, beta, bs, f0, p, q, sas, sbs, sbs3, t1, t10, t13, &
      t15, t16, t19, t2, t22, t24, t25, t32, t34, t36, t39, t4, t40, t41, &
      t44, t48, t49, t5, t6, t65, t8, t87, t9, x, ys

    beta = 0.0042_dbl
    f0 = 1._dbl/sfac
    p = -beta/flsd
    q = 6._dbl*beta
!$omp parallel do private(ip,x,bs,sbs,as,sas,ys,asp,sbs3) &
!$omp   private(t1,t2,t4,t5,t6,t8,t9,t10,t13,t15,t16,t19,t22,t24,t25,t32,t34) &
!$omp   private(t36,t39,t40,t41,t44,t48,t49,t65,t87)
    DO ip=1,SIZE(s)
      x = s(ip) * f0
      bs = beta*x
      sbs = SQRT(x*x+1._dbl)
      as = LOG(x+sbs)
      sas = x*as
      ys =  1._dbl/(1._dbl+q*sas)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + p*x*x*ys
      CASE (1)
         asp = as + x / sbs
         fs(ip,1) = 1._dbl + p*x*x*ys
         fs(ip,2) = ( 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys ) * f0
      CASE (2)
         asp = as + x / sbs
         sbs3 = 1._dbl/(sbs*sbs*sbs)
         fs(ip,1) = 1._dbl + p*x*x*ys
         fs(ip,2) = ( 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys ) * f0
         fs(ip,3) = -f0*f0*p*ys**3*sbs3 * ( q*x*x*x*x*(q*sas + 5._dbl &
                    - 2._dbl*q*sbs) + 2._dbl*( x*x*(q*q*sas & 
                    + 3._dbl*q - sbs) - sbs))
      CASE (3)
         asp = as + x / sbs
         sbs3 = 1._dbl/(sbs*sbs*sbs)
         fs(ip,1) = 1._dbl + p*x*x*ys
         fs(ip,2) = ( 2._dbl*p*x*ys - p*q*x*x*asp*ys*ys ) * f0
         fs(ip,3) = -f0*f0*p*ys**3*sbs3 * ( q*x*x*x*x*(q*sas + 5._dbl &
                    - 2._dbl*q*sbs) + 2._dbl*( x*x*(q*q*sas & 
                    + 3._dbl*q - sbs) - sbs))
      t1 = q*x
      t2 = x**2
      t4 = SQRT(1+t2)
      t5 = x+t4
      t6 = LOG(t5)
      t8 = 1+t1*t6
      t9 = t8**2
      t10 = 1/t9
      t13 = 1/t4
      t15 = 1+t13*x
      t16 = 1/t5
      t19 = q*t6+t1*t15*t16
      t22 = p*x
      t24 = 1/t9/t8
      t25 = t19**2
      t32 = t4**2
      t34 = 1/t32/t4
      t36 = -t34*t2+t13
      t39 = t15**2
      t40 = t5**2
      t41 = 1/t40
      t44 = 2*q*t15*t16+t1*t36*t16-t1*t39*t41
      t48 = p*t2
      t49 = t9**2
      t65 = t32**2
      t87 = -6*p*t10*t19+12*t22*t24*t25-6*t22*t10*t44-6*t48/t49*t25*t19+ &
      6*t48*t24*t19*t44-t48*t10*(3*q*t36*t16-3*q*t39*t41+3*t1*(1/t65/t4* &
      t2*x-t34*x)*t16-3*t1*t36*t41*t15+2*t1*t39*t15/t40/t5)

         fs(ip,4) = t87
         fs(ip,4) = f0*f0*f0*fs(ip,4)

      CASE DEFAULT
         CALL stop_program("efactor_b88","Illegal order")
      END SELECT
    END DO


  END SUBROUTINE efactor_b88
! *****************************************************************************
  SUBROUTINE efactor_pw86 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl)                                :: f15, p0, p1, p15, s2, s4, s6, &
                                                t1, t10, t12, t13, t14, t19, &
                                                t2, t25, t3, t8, t9

    t1 = 1.296_dbl
    t2 = 14._dbl
    t3 = 0.2_dbl
    f15 = 1._dbl/15._dbl
!$omp parallel do private(ip,s2,s4,s6,p0,p1,p15,f15)&
!$omp  private( t1, t2, t3, t8, t9, t10, t12, t13, t14, t19, t25)
    DO ip=1,SIZE(s)
      s2 = s(ip)*s(ip)
      s4 = s2*s2
      s6 = s2*s4
      SELECT CASE (m)
      CASE (0)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         fs(ip,1) = p0**f15
      CASE (1)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2._dbl*t1 + 4._dbl*t2*s2 + 6._dbl*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
      CASE (2)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2._dbl*t1 + 4._dbl*t2*s2 + 6._dbl*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
         t9 = p15**2; t10 = t9**2; t12 = t10**2; t13 = t12*t10*t9
         t25 = p1*p1
         fs(ip,3) = -14._dbl/225._dbl/t13/p0*t25 + &
                     1._dbl/t13*(2._dbl*t1+12*t2*s2+30._dbl*t3*s4)/15._dbl
      CASE (3)
         p0 = 1._dbl + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2._dbl*t1 + 4._dbl*t2*s2 + 6._dbl*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
         t9 = p15**2; t10 = t9**2; t12 = t10**2; t13 = t12*t10*t9
         t25 = p1*p1
         fs(ip,3) = -14._dbl/225._dbl/t13/p0*t25 + &
                     1._dbl/t13*(2._dbl*t1+12*t2*s2+30._dbl*t3*s4)/15._dbl
         t8 = p0**2; t9 = p0**f15; t14 = p0/t9; t19 = s2*s(ip)
         fs(ip,4) = 406._dbl/3375._dbl/t14/t8*p1*p1*p1-14._dbl/&
                    75._dbl/t14/p0*p1*(2*t1+12*t2*s2+30*t3*s4)+&
                    1/t14*(24*t2*s(ip)+120*t3*t19)*f15
      CASE DEFAULT
         CALL stop_program("efactor_ol1","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_pw86
! *****************************************************************************
  SUBROUTINE efactor_optx ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    REAL(dbl), PARAMETER                     :: a1 = 1.05151_dbl, &
                                                a2 = 1.43169_dbl, &
                                                gamma_bo = 0.006_dbl

    INTEGER                                  :: ip
    REAL(dbl)                                :: a, b, f0, x, y

    f0 = 1._dbl/sfac
    b = -a2/flsd
!$omp parallel do private(ip,x,a,y)
    DO ip=1,SIZE(s)
      x = s(ip) * f0
      a = gamma_bo*x*x
      y = 1._dbl / (1._dbl + a)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = a1 + b * a*a*y*y
      CASE (1)
         fs(ip,1) = a1 + b * a*a*y*y
         fs(ip,2) = 4._dbl*b*f0*a*gamma_bo*x * y*y*y
      CASE (2)
         fs(ip,1) = a1 + b * a*a*y*y
         fs(ip,2) = 4._dbl*b*f0*a*gamma_bo*x * y*y*y
         fs(ip,3) = -12._dbl * b*f0*f0 * gamma_bo*a*(a-1._dbl) * y*y*y*y
      CASE (3)
         fs(ip,1) = a1 + b * a*a*y*y
         fs(ip,2) = 4._dbl*b*f0*a*gamma_bo*x * y*y*y
         fs(ip,3) = -12._dbl * b*f0*f0 * gamma_bo*a*(a-1._dbl) * y*y*y*y
         fs(ip,4) = 24._dbl*b*f0*f0*f0* gamma_bo*gamma_bo*x * &
                    (1._dbl-5._dbl*a+2._dbl*a*a) * y*y*y*y*y
      CASE DEFAULT
         CALL stop_program("efactor_optx","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_optx
! *****************************************************************************
  SUBROUTINE efactor_pw91 ( s, fs, m )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dbl)                                :: a, a1, a2, a3, a4, b1, o, s1, &
                                                s2, s3, x

    o=1._dbl
    a1 = 0.19645_dbl
    a2 = 0.2743_dbl
    a3 = 0.1508_dbl
    a4 = 100._dbl
    a  = 0.8145161_dbl
    b1 = 0.004_dbl
    IF ( m >= 0 ) THEN
!$omp parallel do private(ip,x)
      DO ip=1,SIZE(s)
        x = s(ip)
        fs(ip,1) = (o+a1*x*LOG(a*x+SQRT(o+a**2*x**2))+a2-a3*EXP(-a4*x**2)*x**2)&
                   /(o+a1*x*LOG(a*x+SQRT(o+a**2*x**2))+b1*x**4)
      END DO
    END IF
    IF ( m >= 1 ) THEN
!$omp parallel do private(ip,x)
      DO ip=1,SIZE(s)
        x = s(ip)
        fs(ip,2) =&
           (a1*LOG(a*x+SQRT(o+a**2*x**2))+a1*x*a/SQRT(o+a**2*x**2)+2*a3*       &
           a4*x**3*EXP(-a4*x**2)-2*a3*EXP(-a4*x**2)*x)/(o+a1*x*LOG(a*x+SQRT(o  &
           +a**2*x**2))+b1*x**4)-(o+a1*x*LOG(a*x+SQRT(o+a**2*x**2))+a2-a3*exp  &
           (-a4*x**2)*x**2)/(o+a1*x*LOG(a*x+SQRT(o+a**2*x**2))+b1*x**4)**2*(a1 &
           *LOG(a*x+SQRT(o+a**2*x**2))+a1*x*a/SQRT(o+a**2*x**2)+4*b1*x**3)
      END DO
    END IF
    IF ( m >= 2 ) THEN
!$omp parallel do private(ip,x,s1)
      DO ip=1,SIZE(s)
        x = s(ip)
        s1 = (2*a1*a/SQRT(1+a**2*x**2)-a1*x**2*a**3/SQRT(1+a**2*x**2)**3+  &
             10*a3*a4*x**2*EXP(-a4*x**2)-4*a3*a4**2*x**4*EXP(-a4*x**2)     &
             -2*a3*EXP(-a4*x**2))/(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+      &
             b1*x**4)-2*(a1*LOG(a*x+SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2) &
             +2*a3*a4*x**3*EXP(-a4*x**2)-2*a3*EXP(-a4*x**2)*x)/(1+a1*x*    &
             LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**2*(a1*LOG(a*x+SQRT(1+a**2*x**2))&
             +a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)
        fs(ip,3) =&
             s1+2*(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-a3*EXP(-a4*x**2)*x**2)&
             /(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**3*(a1*LOG(a*x+     &
             SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)**2-       &
             (1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-a3*EXP(-a4*x**2)*x**2)     &
             /(1+a1*x*LOG(a*x +SQRT(1+a**2*x**2))+b1*x**4)**2*(2*a1*a/        &
             SQRT(1+a**2*x**2)-a1*x**2*a**3/SQRT(1+a**2*x**2)**3+12*b1*x**2)
      END DO
    END IF
    IF ( m >= 3 ) THEN
!$omp parallel do private(ip,x,s1,s2,s3)
      DO ip=1,SIZE(s)
        x = s(ip)
        s2 = (-4*a1*a**3/SQRT(1+a**2*x**2)**3*x+3*a1*x**3*a**5/         &
             SQRT(1+a**2*x**2)**5+24*a3*a4*x*EXP(-a4*x**2)-36*a3*a4**2* &
             x**3*EXP(-a4*x**2)+8*a3*a4**3*x**5*EXP(-a4*x**2))/(1+a1*x* &
             LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)
        s3 = -3*(2*a1*a/SQRT(1+a**2*x**2)-a1*x**2*a**3/SQRT(1+a**2*     &
             x**2)**3+10*a3*a4*x**2*EXP(-a4*x**2)-4*a3*a4**2*x**4*      &
             EXP(-a4*x**2)-2*a3*EXP(-a4*x**2))/(1+a1*x*LOG(a*x+         &
             SQRT(1+a**2*x**2))+b1*x**4)**2*(a1*LOG(a*x+SQRT(1+a**2*    &
             x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)+6*(a1*LOG(a*x+  &
             SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2)+2*a3*a4*x**3*  &
             EXP(-a4*x**2)-2*a3*EXP(-a4*x**2)*x)/(1+a1*x*LOG(a*x+       &
             SQRT(1+a**2*x**2))+b1*x**4)**3*(a1*LOG(a*x+SQRT(1+a**2*    &
             x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)**2
        s1 = s2+s3
        s2 = s1-3*(a1*LOG(a*x+SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2) &
             +2*a3*a4*x**3*EXP(-a4*x**2)-2*a3*EXP(-a4*x**2)*x)/(1+a1*x*   &
             LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**2*(2*a1*a/SQRT(1+a**2*  &
             x**2)-a1*x**2*a**3/SQRT(1+a**2*x**2)**3+12*b1*x**2)
        s3 = s2-6*(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-a3*EXP(-a4*x**2)* &
             x**2)/(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**4*(a1*    &
             LOG(a*x+SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*x**3)**3
        fs(ip,4) =&
             s3+6*(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-a3*EXP(-a4*x**2)  &
             *x**2)/(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+b1*x**4)**3*(a1*   &
             LOG(a*x+SQRT(1+a**2*x**2))+a1*x*a/SQRT(1+a**2*x**2)+4*b1*    &
             x**3)*(2*a1*a/SQRT(1+a**2*x**2)-a1*x**2*a**3/SQRT(1+a**2*    &
             x**2)**3+12*b1*x**2)-(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+a2-  &
             a3*EXP(-a4*x**2)*x**2)/(1+a1*x*LOG(a*x+SQRT(1+a**2*x**2))+   &
             b1*x**4)**2*(-4*a1*a**3/SQRT(1+a**2*x**2)**3*x+3*a1*x**3*    &
             a**5/SQRT(1+a**2*x**2)**5+24*b1*x)
      END DO
    END IF

  END SUBROUTINE efactor_pw91
! *****************************************************************************
  SUBROUTINE efactor_pbex ( s, fs, m, pset )
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: s
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fs
    INTEGER, INTENT(IN)                      :: m, pset

    REAL(dbl), PARAMETER                     :: kappa1 = 0.804_dbl, &
                                                kappa2 = 1.245_dbl, &
                                                mu = 0.2195149727645171_dbl

    INTEGER                                  :: ip
    REAL(dbl)                                :: f0, mk, x, x2, y

    IF( pset==1 ) mk=mu/kappa1
    IF( pset==2 ) mk=mu/kappa2

    f0 = 1._dbl/tact
!$omp parallel do private(ip,x,x2,y)
    DO ip=1,SIZE(s)
      x = s(ip)*f0
      x2 = x*x
      y = 1._dbl/(1._dbl+mk*x2)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1._dbl + mu*x2*y
      CASE (1)
         fs(ip,1) = 1._dbl + mu*x2*y
         fs(ip,2) = 2._dbl*mu*x*y*y*f0
      CASE (2)
         fs(ip,1) = 1._dbl + mu*x2*y
         fs(ip,2) = 2._dbl*mu*x*y*y*f0
         fs(ip,3) = -2._dbl*mu*(3._dbl*mk*x2-1._dbl)*y*y*y*f0*f0
      CASE (3)
         fs(ip,1) = 1._dbl + mu*x2*y
         fs(ip,2) = 2._dbl*mu*x*y*y*f0
         fs(ip,3) = -2._dbl*mu*(3._dbl*mk*x2-1._dbl)*y*y*y*f0*f0
         fs(ip,4) = 24._dbl*mu*mk*x*(mk*x2-1._dbl)*y*y*y*y*f0*f0*f0
      CASE DEFAULT
         CALL stop_program("efactor_pbex","Illegal order")
      END SELECT
    END DO

  END SUBROUTINE efactor_pbex

! *****************************************************************************

END MODULE xc_exchange_gga_new

! *****************************************************************************
