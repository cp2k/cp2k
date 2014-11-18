!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate several different exchange energy functionals
!>      with a GGA form
!> \par History
!>      JGH (26.02.2003) : OpenMP enabled
!>      fawzi (04.2004)  : adapted to the new xc interface
!> \author JGH (27.02.2002)
! *****************************************************************************
MODULE xc_exchange_gga

  USE cp_array_r_utils,                ONLY: cp_3d_r_p_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_functionals_utilities,        ONLY: calc_wave_vector,&
                                             set_util
  USE xc_input_constants,              ONLY: xgga_b88,&
                                             xgga_ev93,&
                                             xgga_opt,&
                                             xgga_pbe,&
                                             xgga_pw86,&
                                             xgga_pw91,&
                                             xgga_revpbe
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: xgga_info, xgga_eval

  REAL(KIND=dp), PARAMETER :: f13 = 1.0_dp/3.0_dp, &
                              f23 = 2.0_dp*f13, &
                              f43 = 4.0_dp*f13, &
                              f53 = 5.0_dp*f13

  REAL(KIND=dp) :: cx, flda, flsd, sfac, t13
  REAL(KIND=dp) :: fact, tact
  REAL(KIND=dp) :: eps_rho
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_exchange_gga'

CONTAINS

! *****************************************************************************
!> \brief return various information on the xgga functionals
!> \param functional integer selecting the xgga functional, it should be one of
!>        the constants defined in this module: xgga_b88, xgga_pw86,...
!> \param lsd a logical that specifies if you are asking about the lsd or lda
!>        version of the functional
!> \param reference string with the reference of the actual functional
!> \param shortform string with the shortform of the functional name
!> \param needs the components needed by this functional are set to
!>        true (does not set the unneeded components to false)
!> \param max_deriv ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE xgga_info(functional,lsd,reference,shortform, needs, max_deriv, error)
    INTEGER, INTENT(in)                      :: functional
    LOGICAL, INTENT(in)                      :: lsd
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

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
       CASE(xgga_ev93)
          reference = "E. Engel and S.H. Vosko, Phys. Rev. B, 47, 13164 (1993)"
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
       CASE(xgga_ev93)
          shortform = "Engel-Vosko exchange energy from virial relation"
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

! *****************************************************************************
!> \brief evaluates different exchange gga
!> \param functional integer to select the functional that should be evaluated
!> \param lsd if the lsd version of the functional should be used
!> \param rho_set the density where you want to evaluate the functional
!> \param deriv_set place where to store the functional derivatives (they are
!>        added to the derivatives)
!> \param order ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE xgga_eval(functional,lsd,rho_set,deriv_set,order,error)
    INTEGER, INTENT(in)                      :: functional
    LOGICAL, INTENT(in)                      :: lsd
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: order
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xgga_eval', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=12), DIMENSION(2)          :: norm_drho_spin_name
    CHARACTER(len=6), DIMENSION(2)           :: rho_spin_name
    INTEGER                                  :: handle, ispin, m, npoints, &
                                                nspin, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: drho_cutoff, rho_cutoff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: s
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: fs
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: e_0, e_ndrho, &
      e_ndrho_ndrho, e_ndrho_ndrho_ndrho, e_rho, e_rho_ndrho, &
      e_rho_ndrho_ndrho, e_rho_rho, e_rho_rho_ndrho, e_rho_rho_rho
    TYPE(cp_3d_r_p_type), DIMENSION(2)       :: norm_drho, rho, rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
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
             tact = 1.0_dp
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

             !! omp: note this is handled slightly differently to the
             !! other cases to prevent sprawling scope declarations
             !! in efactor_pw91()

             !$omp parallel default (none) shared(s, fs, m)
             CALL efactor_pw91 ( s, fs, m )
             !$omp end parallel

          CASE (xgga_pbe)
             tact = t13
             CALL efactor_pbex ( s, fs, m, 1 )
             IF (lsd) tact = 1.0_dp
          CASE (xgga_revpbe)
             tact = t13
             CALL efactor_pbex ( s, fs, m, 2 )
             IF (lsd) tact = 1.0_dp
          CASE (xgga_opt)
             CALL efactor_optx ( s, fs, m )
          CASE (xgga_ev93)
             tact = t13
             CALL efactor_ev93 ( s, fs, m )
             IF (lsd) tact = 1.0_dp
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
    CALL timestop(handle)
  END SUBROUTINE xgga_eval

! *****************************************************************************
!> \brief ...
!> \param cutoff ...
! *****************************************************************************
  SUBROUTINE xgga_init ( cutoff )

    REAL(KIND=dp), INTENT(IN)                :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

    cx = -0.75_dp*(3.0_dp/pi)**f13
    t13 = 2.0_dp**f13
    flda = cx
    flsd = cx * t13

    sfac = 1.0_dp / (2.0_dp*(3.0_dp*pi*pi)**f13)

  END SUBROUTINE xgga_init

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param r13 ...
!> \param fs ...
!> \param e_0 ...
!> \param npoints ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE x_p_0 ( rho, r13, fs, e_0, npoints, error )

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, r13
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: fs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_0
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip

    !$omp parallel do default (none) &
    !$omp             shared (npoints, rho, eps_rho, fact, r13, fs, e_0) &
    !$omp             private(ip)

    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         e_0(ip) = e_0(ip) + fact*r13(ip)*rho(ip) * fs(ip,1)
      END IF

    END DO

    !$omp end parallel do

  END SUBROUTINE x_p_0

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param grho ...
!> \param r13 ...
!> \param s ...
!> \param fs ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param npoints ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE x_p_1 ( rho, grho, r13, s, fs, e_rho,e_ndrho,npoints,error )

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, grho, r13, s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: fs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_rho, e_ndrho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: a0, a1, sx, sy

    !$omp parallel do default(none) &
    !$omp             shared(npoints, rho, eps_rho, fact, r13, tact, fs) &
    !$omp             shared(e_rho, e_ndrho, sfac, s) &
    !$omp             private(ip,a0,a1,sx,sy)

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

    !$omp end parallel do

  END SUBROUTINE x_p_1

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param grho ...
!> \param r13 ...
!> \param s ...
!> \param fs ...
!> \param e_rho_rho ...
!> \param e_rho_ndrho ...
!> \param e_ndrho_ndrho ...
!> \param npoints ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE x_p_2 ( rho, grho, r13, s, fs, e_rho_rho, e_rho_ndrho, &
       e_ndrho_ndrho, npoints, error )

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, grho, r13, s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: fs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_rho_rho, e_rho_ndrho, &
                                                e_ndrho_ndrho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: a0, a1, a2, sx, sxx, sxy, sy

    !$omp parallel do default(none) &
    !$omp             shared(npoints, rho, eps_rho, r13, fact, e_rho_rho) &
    !$omp             shared(e_rho_ndrho, e_ndrho_ndrho, fs, sfac, tact, s) &
    !$omp             private(ip,a0,a1,a2,sx,sy,sxx,sxy)

    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*rho(ip)
         a1 = f43*fact*r13(ip)
         a2 = f13*f43*fact/(r13(ip)*r13(ip))
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         sxx= 28.0_dp/9.0_dp*s(ip)/(rho(ip)*rho(ip))
         sxy= -f43*sfac*tact/(r13(ip)*rho(ip)*rho(ip))
         e_rho_rho(ip) = e_rho_rho(ip) + a2*fs(ip,1) +2.0_dp*a1*fs(ip,2)*sx + &
                     a0*fs(ip,3)*sx*sx + a0*fs(ip,2)*sxx
         e_rho_ndrho(ip) = e_rho_ndrho(ip) &
              + a1*fs(ip,2)*sy + a0*fs(ip,3)*sx*sy + a0*fs(ip,2)*sxy
         e_ndrho_ndrho(ip) = e_ndrho_ndrho(ip) + a0*fs(ip,3)*sy*sy

      END IF

    END DO

    !$omp end parallel do

  END SUBROUTINE x_p_2

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param grho ...
!> \param r13 ...
!> \param s ...
!> \param fs ...
!> \param e_rho_rho_rho ...
!> \param e_rho_rho_ndrho ...
!> \param e_rho_ndrho_ndrho ...
!> \param e_ndrho_ndrho_ndrho ...
!> \param npoints ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE x_p_3 ( rho, grho, r13, s, fs, e_rho_rho_rho, e_rho_rho_ndrho,&
       e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho, npoints,error )

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, grho, r13, s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: fs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_rho_rho_rho, &
                                                e_rho_rho_ndrho, &
                                                e_rho_ndrho_ndrho, &
                                                e_ndrho_ndrho_ndrho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: a0, a1, a2, a3, sx, sxx, &
                                                sxxx, sxxy, sxy, sy

    !$omp parallel do default (none) &
    !$omp             shared(npoints, rho, eps_rho, r13, fact, fs) &
    !$omp             shared(e_rho_rho_rho, e_rho_rho_ndrho) &
    !$omp             shared(e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho) &
    !$omp             shared(sfac, tact, s) &
    !$omp             private(ip,a0,a1,a2,a3,sx,sy,sxx,sxy,sxxx,sxxy)

    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN

         a0 = fact*r13(ip)*rho(ip)
         a1 = f43*fact*r13(ip)
         a2 = f13*f43*fact/(r13(ip)*r13(ip))
         a3 = -f23*f13*f43*fact/(r13(ip)*r13(ip)*rho(ip))
         sx = -f43*s(ip)/rho(ip)
         sy = sfac*tact/(r13(ip)*rho(ip))
         sxx= 28.0_dp/9.0_dp*s(ip)/(rho(ip)*rho(ip))
         sxy= -f43*sfac*tact/(r13(ip)*rho(ip)*rho(ip))
         sxxx= -280.0_dp/27.0_dp*s(ip)/(rho(ip)*rho(ip)*rho(ip))
         sxxy= 28.0_dp/9.0_dp*sfac*tact/(r13(ip)*rho(ip)*rho(ip)*rho(ip))
         e_rho_rho_rho(ip) = e_rho_rho_rho(ip)&
                     + a3*fs(ip,1) + 3.0_dp*a2*fs(ip,2)*sx + &
                     3.0_dp*a1*fs(ip,3)*sx*sx + 3.0_dp*a1*fs(ip,2)*sxx + &
                     a0*fs(ip,4)*sx*sx*sx + 3.0_dp*a0*fs(ip,3)*sx*sxx + &
                     a0*fs(ip,2)*sxxx
         e_rho_rho_ndrho(ip) = e_rho_rho_ndrho(ip)&
                     + a2*fs(ip,2)*sy + 2.0_dp*a1*fs(ip,3)*sx*sy + &
                     2.0_dp*a1*fs(ip,2)*sxy + a0*fs(ip,4)*sx*sx*sy + &
                     2.0_dp*a0*fs(ip,3)*sx*sxy + a0*fs(ip,3)*sxx*sy + &
                     a0*fs(ip,2)*sxxy
         e_rho_ndrho_ndrho(ip) = e_rho_ndrho_ndrho(ip)&
                     + a1*fs(ip,3)*sy*sy + a0*fs(ip,4)*sx*sy*sy + &
                     2.0_dp*a0*fs(ip,3)*sxy*sy
         e_ndrho_ndrho_ndrho(ip) = e_ndrho_ndrho_ndrho(ip)&
                     + a0*fs(ip,4)*sy*sy*sy

      END IF

    END DO

    !$omp end parallel do

  END SUBROUTINE x_p_3

! Enhancement Factors
! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE efactor_b88 ( s, fs, m )
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: fs
    INTEGER, INTENT(IN)                      :: m

    CHARACTER(len=*), PARAMETER :: routineN = 'efactor_b88', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(KIND=dp) :: as, asp, beta, bs, f0, p, q, sas, sbs, sbs3, t1, t10, &
      t13, t15, t16, t19, t2, t22, t24, t25, t32, t34, t36, t39, t4, t40, &
      t41, t44, t48, t49, t5, t6, t65, t8, t87, t9, x, ys

    beta = 0.0042_dp
    f0 = 1.0_dp/sfac
    p = -beta/flsd
    q = 6.0_dp*beta

    !$omp parallel do default(none) &
    !$omp             shared(s,fs,m,beta,f0,p,q) &
    !$omp             private(ip,x,bs,sbs,as,sas,ys,asp,sbs3) &
    !$omp             private(t1,t2,t4,t5,t6,t8,t9,t10,t13,t15,t16,t19,t22) &
    !$omp             private(t24,t25,t32,t34) &
    !$omp             private(t36,t39,t40,t41,t44,t48,t49,t65,t87)

    DO ip=1,SIZE(s)
      x = s(ip) * f0
      bs = beta*x
      sbs = SQRT(x*x+1.0_dp)
      as = LOG(x+sbs)
      sas = x*as
      ys =  1.0_dp/(1.0_dp+q*sas)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1.0_dp + p*x*x*ys
      CASE (1)
         asp = as + x / sbs
         fs(ip,1) = 1.0_dp + p*x*x*ys
         fs(ip,2) = ( 2.0_dp*p*x*ys - p*q*x*x*asp*ys*ys ) * f0
      CASE (2)
         asp = as + x / sbs
         sbs3 = 1.0_dp/(sbs*sbs*sbs)
         fs(ip,1) = 1.0_dp + p*x*x*ys
         fs(ip,2) = ( 2.0_dp*p*x*ys - p*q*x*x*asp*ys*ys ) * f0
         fs(ip,3) = -f0*f0*p*ys**3*sbs3 * ( q*x*x*x*x*(q*sas + 5.0_dp &
                    - 2.0_dp*q*sbs) + 2.0_dp*( x*x*(q*q*sas &
                    + 3.0_dp*q - sbs) - sbs))
      CASE (3)
         asp = as + x / sbs
         sbs3 = 1.0_dp/(sbs*sbs*sbs)
         fs(ip,1) = 1.0_dp + p*x*x*ys
         fs(ip,2) = ( 2.0_dp*p*x*ys - p*q*x*x*asp*ys*ys ) * f0
         fs(ip,3) = -f0*f0*p*ys**3*sbs3 * ( q*x*x*x*x*(q*sas + 5.0_dp &
                    - 2.0_dp*q*sbs) + 2.0_dp*( x*x*(q*q*sas &
                    + 3.0_dp*q - sbs) - sbs))
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
         CALL stop_program(routineN,moduleN,__LINE__,"Illegal order")
      END SELECT
    END DO

    !$omp end parallel do

  END SUBROUTINE efactor_b88
! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE efactor_pw86 ( s, fs, m )
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: fs
    INTEGER, INTENT(IN)                      :: m

    CHARACTER(len=*), PARAMETER :: routineN = 'efactor_pw86', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: f15, p0, p1, p15, s2, s4, s6, &
                                                t1, t10, t12, t13, t14, t19, &
                                                t2, t25, t3, t8, t9

    t1 = 1.296_dp
    t2 = 14.0_dp
    t3 = 0.2_dp
    f15 = 1.0_dp/15.0_dp

    !$omp parallel do default(none) &
    !$omp             shared(s, fs, m, t1, t2, t3, f15) &
    !$omp             private(ip,s2,s4,s6,p0,p1,p15)&
    !$omp             private(t8, t9, t10, t12, t13, t14, t19, t25)

    DO ip=1,SIZE(s)
      s2 = s(ip)*s(ip)
      s4 = s2*s2
      s6 = s2*s4
      SELECT CASE (m)
      CASE (0)
         p0 = 1.0_dp + t1*s2 + t2*s4 + t3*s6
         fs(ip,1) = p0**f15
      CASE (1)
         p0 = 1.0_dp + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2.0_dp*t1 + 4.0_dp*t2*s2 + 6.0_dp*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
      CASE (2)
         p0 = 1.0_dp + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2.0_dp*t1 + 4.0_dp*t2*s2 + 6.0_dp*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
         t9 = p15**2; t10 = t9**2; t12 = t10**2; t13 = t12*t10*t9
         t25 = p1*p1
         fs(ip,3) = -14.0_dp/225.0_dp/t13/p0*t25 + &
                     1.0_dp/t13*(2.0_dp*t1+12*t2*s2+30.0_dp*t3*s4)/15.0_dp
      CASE (3)
         p0 = 1.0_dp + t1*s2 + t2*s4 + t3*s6
         p1 = s(ip)*(2.0_dp*t1 + 4.0_dp*t2*s2 + 6.0_dp*t3*s4)
         p15 = p0**f15
         fs(ip,1) = p15
         fs(ip,2) = f15 * p1 * p15/p0
         t9 = p15**2; t10 = t9**2; t12 = t10**2; t13 = t12*t10*t9
         t25 = p1*p1
         fs(ip,3) = -14.0_dp/225.0_dp/t13/p0*t25 + &
                     1.0_dp/t13*(2.0_dp*t1+12*t2*s2+30.0_dp*t3*s4)/15.0_dp
         t8 = p0**2; t9 = p0**f15; t14 = p0/t9; t19 = s2*s(ip)
         fs(ip,4) = 406.0_dp/3375.0_dp/t14/t8*p1*p1*p1-14.0_dp/&
                    75.0_dp/t14/p0*p1*(2*t1+12*t2*s2+30*t3*s4)+&
                    1/t14*(24*t2*s(ip)+120*t3*t19)*f15
      CASE DEFAULT
         CALL stop_program(routineN,moduleN,__LINE__,"Illegal order")
      END SELECT
    END DO

    !$omp end parallel do

  END SUBROUTINE efactor_pw86
! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE efactor_ev93 ( s, fs, m )
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: fs
    INTEGER, INTENT(IN)                      :: m

    CHARACTER(len=*), PARAMETER :: routineN = 'efactor_ev93', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: a1, a2, a3, b1, b2, b3, d0, &
                                                d1, d2, d3, f0, f1, f2, n0, &
                                                n1, n2, n3, s2, s4, s6, &
                                                scale_s, ss

    a1 = 1.647127_dp
    a2 = 0.980118_dp
    a3 = 0.017399_dp
    b1 = 1.523671_dp
    b2 = 0.367229_dp
    b3 = 0.011282_dp
    scale_s=1._dp/tact

    !$omp parallel do default(none) &
    !$omp             shared(s, fs, m, a1, a2, a3, b1, b2, b3, scale_s) &
    !$omp             private(ip,ss,s2,s4,s6) &
    !$omp             private(n0,n1,n2,n3,d0,d1,d2,d3,f0,f1,f2)

    DO ip=1,SIZE(s)
!     ss = s(ip)
!
      ss = scale_s*s(ip)
      s2 = ss*ss
      s4 = s2*s2
      s6 = s2*s4
      SELECT CASE (m)
      CASE (0)
         n0 = 1._dp + a1*s2 + a2*s4 + a3*s6
         d0 = 1._dp + b1*s2 + b2*s4 + b3*s6
         fs(ip,1) = n0/d0
      CASE (1)
         n0 = 1._dp + a1*s2 + a2*s4 + a3*s6
         d0 = 1._dp + b1*s2 + b2*s4 + b3*s6
         n1 = ss*(2._dp*a1 + 4._dp*a2*s2 + 6._dp*a3*s4)
         d1 = ss*(2._dp*b1 + 4._dp*b2*s2 + 6._dp*b3*s4)
         f0 = n0/d0
         fs(ip,1) = f0
         fs(ip,2) = (n1-f0*d1)/d0 * scale_s
      CASE (2)
         n0 = 1._dp + a1*s2 + a2*s4 + a3*s6
         d0 = 1._dp + b1*s2 + b2*s4 + b3*s6
         n1 = ss*(2._dp*a1 + 4._dp*a2*s2 + 6._dp*a3*s4)
         d1 = ss*(2._dp*b1 + 4._dp*b2*s2 + 6._dp*b3*s4)
         n2 = 2._dp*a1 + 12._dp*a2*s2 + 30._dp*a3*s4
         d2 = 2._dp*b1 + 12._dp*b2*s2 + 30._dp*b3*s4
         f0 = n0/d0
         f1 = (n1-f0*d1)/d0
         fs(ip,1) = f0
         fs(ip,2) = f1 * scale_s
         fs(ip,3) = (n2-f0*d2-2._dp*f1*d1)/d0 * scale_s * scale_s
      CASE (3)
         n0 = 1._dp + a1*s2 + a2*s4 + a3*s6
         d0 = 1._dp + b1*s2 + b2*s4 + b3*s6
         n1 = ss*(2._dp*a1 + 4._dp*a2*s2 + 6._dp*a3*s4)
         d1 = ss*(2._dp*b1 + 4._dp*b2*s2 + 6._dp*b3*s4)
         n2 = 2._dp*a1 + 12._dp*a2*s2 + 30._dp*a3*s4
         d2 = 2._dp*b1 + 12._dp*b2*s2 + 30._dp*b3*s4
         n3 = ss*(24._dp*a2 + 120._dp*a3*s2)
         d3 = ss*(24._dp*b2 + 120._dp*b3*s2)
         f0 = n0/d0
         f1 = (n1-f0*d1)/d0
         f2 = (n2-f0*d2-2._dp*f1*d1)/d0
         fs(ip,1) = f0
         fs(ip,2) = f1 * scale_s
         fs(ip,3) = f2 * scale_s * scale_s
         fs(ip,4) = (n3-f0*d3-3._dp*f2*d1-3._dp*f1*d2)/d0 * scale_s * scale_s * scale_s
      CASE DEFAULT
         CALL stop_program(routineN,moduleN,__LINE__,"Illegal order")
      END SELECT
    END DO

    !$omp end parallel do

  END SUBROUTINE efactor_ev93
! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE efactor_optx ( s, fs, m )
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: fs
    INTEGER, INTENT(IN)                      :: m

    CHARACTER(len=*), PARAMETER :: routineN = 'efactor_optx', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: a1 = 1.05151_dp, &
                                                a2 = 1.43169_dp, &
                                                gamma_bo = 0.006_dp

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: a, b, f0, x, y

    f0 = 1.0_dp/sfac
    b = -a2/flsd

    !$omp parallel do default(none) &
    !$omp             shared(s, fs, m, f0, b) &
    !$omp             private(ip,x,a,y)

    DO ip=1,SIZE(s)
      x = s(ip) * f0
      a = gamma_bo*x*x
      y = 1.0_dp / (1.0_dp + a)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = a1 + b * a*a*y*y
      CASE (1)
         fs(ip,1) = a1 + b * a*a*y*y
         fs(ip,2) = 4.0_dp*b*f0*a*gamma_bo*x * y*y*y
      CASE (2)
         fs(ip,1) = a1 + b * a*a*y*y
         fs(ip,2) = 4.0_dp*b*f0*a*gamma_bo*x * y*y*y
         fs(ip,3) = -12.0_dp * b*f0*f0 * gamma_bo*a*(a-1.0_dp) * y*y*y*y
      CASE (3)
         fs(ip,1) = a1 + b * a*a*y*y
         fs(ip,2) = 4.0_dp*b*f0*a*gamma_bo*x * y*y*y
         fs(ip,3) = -12.0_dp * b*f0*f0 * gamma_bo*a*(a-1.0_dp) * y*y*y*y
         fs(ip,4) = 24.0_dp*b*f0*f0*f0* gamma_bo*gamma_bo*x * &
                    (1.0_dp-5.0_dp*a+2.0_dp*a*a) * y*y*y*y*y
      CASE DEFAULT
         CALL stop_program(routineN,moduleN,__LINE__,"Illegal order")
      END SELECT
    END DO

    !$omp end parallel do

  END SUBROUTINE efactor_optx

! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE efactor_pw91 ( s, fs, m )
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: fs
    INTEGER, INTENT(IN)                      :: m

    INTEGER                                  :: ip
    REAL(dp) :: t1, t10, t101, t106, t109, t111, t113, t119, t12, t123, t124, &
      t13, t14, t15, t16, t17, t18, t19, t2, t20, t21, t22, t23, t25, t26, &
      t27, t28, t29, t3, t30, t31, t33, t35, t37, t38, t39, t4, t40, t44, &
      t47, t48, t5, t50, t51, t53, t55, t56, t57, t58, t59, t6, t60, t64, &
      t65, t69, t7, t70, t71, t73, t77, t78, t8, t80, t82, t9, t90, t93, t94, &
      t96, t98
    REAL(KIND=dp)                            :: a1, a2, a3, a4, a5, b, o, x

    o=1.0_dp
    a1 = 0.19645_dp
    a2 = 0.2743_dp
    a3 = 0.1508_dp
    a4 = 100.0_dp
    b = 0.8145161_dp
    a5 = 0.004_dp
    IF ( m >= 0 ) THEN

      !$omp do

      DO ip=1,SIZE(s)
        x = s(ip)
        t3 = b**2
        t4 = x**2
        t7 = SQRT(o+t3*t4)
        t9 = LOG(b*x+t7)
        t10 = a1*x*t9
        t12 = EXP(-a4*t4)
        t17 = t4**2
        fs(ip,1) = (o+t10+(a2-a3*t12)*t4)/(o+t10+a5*t17)
      END DO

      !$omp end do

    END IF
    IF ( m >= 1 ) THEN

      !$omp do

      DO ip=1,SIZE(s)
        x = s(ip)
        t2 = b**2
        t3 = x**2
        t6 = SQRT(o+t2*t3)
        t7 = b*x+t6
        t8 = LOG(t7)
        t9 = a1*t8
        t10 = a1*x
        t17 = t10*(b+1/t6*t2*x)/t7
        t19 = t3*x
        t21 = EXP(-a4*t3)
        t26 = a2-a3*t21
        t30 = t10*t8
        t31 = t3**2
        t33 = o+t30+a5*t31
        t38 = t33**2
        fs(ip,2) =&
           (t9+t17+2._dp*a3*a4*t19*t21+2._dp*t26*x)/&
           t33-(o+t30+t26*t3)/t38*(t9+t17+4._dp*a5*t19)
      END DO

      !$omp end do

    END IF
    IF ( m >= 2 ) THEN

      !$omp do

      DO ip=1,SIZE(s)
        x = s(ip)
        t1 = b**2
        t2 = x**2
        t5 = SQRT(o+t1*t2)
        t7 = o/t5*t1
        t9 = b+t7*x
        t12 = b*x+t5
        t13 = o/t12
        t15 = 2._dp*a1*t9*t13
        t16 = a1*x
        t17 = t5**2
        t20 = t1**2
        t25 = t16*(-o/t17/t5*t20*t2+t7)*t13
        t26 = t9**2
        t27 = t12**2
        t30 = t16*t26/t27
        t31 = a3*a4
        t33 = EXP(-a4*t2)
        t37 = a4**2
        t39 = t2**2
        t44 = a3*t33
        t47 = LOG(t12)
        t48 = t16*t47
        t50 = o+t48+a5*t39
        t53 = a1*t47
        t55 = t16*t9*t13
        t56 = t2*x
        t60 = a2-t44
        t64 = t50**2
        t65 = o/t64
        t69 = t53+t55+4._dp*a5*t56
        t73 = o+t48+t60*t2
        t77 = t69**2
        fs(ip,3) = &
         (t15+t25-t30+10._dp*t31*t2*t33-4._dp*a3*t37*t39*t33+&
          2._dp*a2-2._dp*t44)/t50-2._dp*&
         (t53+t55+2._dp*t31*t56*t33+2._dp*t60*x)*&
         t65*t69+2._dp*t73/t64/t50*t77-t73*t65*(t15+t25-t30+12._dp*a5*t2)
      END DO

      !$omp end do

    END IF
    IF ( m >= 3 ) THEN

      !$omp do

      DO ip=1,SIZE(s)
        x = s(ip)
          t1 = b**2
          t2 = x**2
          t5 = SQRT(0.1e1_dp+t1*t2)
          t6 = t5**2
          t9 = t1**2
          t10 = 1/t6/t5*t9
          t13 = 1/t5*t1
          t14 = -t10*t2+t13
          t17 = b*x+t5
          t18 = 1/t17
          t20 = 3*a1*t14*t18
          t22 = b+t13*x
          t23 = t22**2
          t25 = t17**2
          t26 = 1/t25
          t28 = 3*a1*t23*t26
          t29 = a1*x
          t30 = t6**2
          t35 = t2*x
          t40 = 3*t29*(1/t30/t5*t1*t9*t35-t10*x)*t18
          t44 = 3*t29*t14*t26*t22
          t50 = 2*t29*t23*t22/t25/t17
          t51 = a3*a4
          t53 = EXP(-a4*t2)
          t57 = a4**2
          t58 = a3*t57
          t59 = t35*t53
          t64 = t2**2
          t70 = LOG(t17)
          t71 = t29*t70
          t73 = 0.1e1_dp+t71+a5*t64
          t78 = 2*a1*t22*t18
          t80 = t29*t14*t18
          t82 = t29*t23*t26
          t90 = a3*t53
          t93 = t73**2
          t94 = 1/t93
          t96 = a1*t70
          t98 = t29*t18*t22
          t101 = t96+t98+4*a5*t35
          t106 = a2-t90
          t109 = t96+t98+2*t51*t59+2*t106*x
          t111 = 1/t93/t73
          t113 = t101**2
          t119 = t78+t80-t82+12*a5*t2
          t123 = 0.1e1_dp+t71+t106*t2
          t124 = t93**2
          fs(ip,4) = &
            (t20-t28+t40-t44+t50+24*t51*x*t53-36._dp*t58*t59+8._dp*a3*t57*a4*t64*&
            x*t53)/t73-3._dp*(t78+t80-t82+10._dp*t51*t2*t53-&
            4._dp*t58*t64*t53+2._dp*a2-2._dp*t90)*t94*t101+&
            6._dp*t109*t111*t113-3._dp*t109*t94*t119-6*t123/t124*t113*t101+&
            6._dp*t123*t111*t101*t119-t123*t94*(t20-t28+t40-t44+t50+24._dp*a5*x)
      END DO

      !$omp end do

    END IF

  END SUBROUTINE efactor_pw91

! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param fs ...
!> \param m ...
!> \param pset ...
! *****************************************************************************
  SUBROUTINE efactor_pbex ( s, fs, m, pset )
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: s
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: fs
    INTEGER, INTENT(IN)                      :: m, pset

    CHARACTER(len=*), PARAMETER :: routineN = 'efactor_pbex', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: kappa1 = 0.804_dp, &
                                                kappa2 = 1.245_dp, &
                                                mu = 0.2195149727645171_dp

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: f0, mk, x, x2, y

    IF( pset==1 ) mk=mu/kappa1
    IF( pset==2 ) mk=mu/kappa2

    f0 = 1.0_dp/tact

    !$omp parallel do default(none) &
    !$omp             shared (s, fs, m, mk, f0) &
    !$omp             private(ip,x,x2,y)

    DO ip=1,SIZE(s)
      x = s(ip)*f0
      x2 = x*x
      y = 1.0_dp/(1.0_dp+mk*x2)
      SELECT CASE (m)
      CASE (0)
         fs(ip,1) = 1.0_dp + mu*x2*y
      CASE (1)
         fs(ip,1) = 1.0_dp + mu*x2*y
         fs(ip,2) = 2.0_dp*mu*x*y*y*f0
      CASE (2)
         fs(ip,1) = 1.0_dp + mu*x2*y
         fs(ip,2) = 2.0_dp*mu*x*y*y*f0
         fs(ip,3) = -2.0_dp*mu*(3.0_dp*mk*x2-1.0_dp)*y*y*y*f0*f0
      CASE (3)
         fs(ip,1) = 1.0_dp + mu*x2*y
         fs(ip,2) = 2.0_dp*mu*x*y*y*f0
         fs(ip,3) = -2.0_dp*mu*(3.0_dp*mk*x2-1.0_dp)*y*y*y*f0*f0
         fs(ip,4) = 24.0_dp*mu*mk*x*(mk*x2-1.0_dp)*y*y*y*y*f0*f0*f0
      CASE DEFAULT
         CALL stop_program(routineN,moduleN,__LINE__,"Illegal order")
      END SELECT
    END DO

    !$omp end parallel do

  END SUBROUTINE efactor_pbex

END MODULE xc_exchange_gga

