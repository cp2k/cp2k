!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/xc_xalpha [1.0] *
!!
!!   NAME
!!     xc_xalpha
!!
!!   FUNCTION
!!     Calculate the local exchange functional
!!
!!   NOTES
!!     Order of derivatives is: LDA 0; 1; 2; 3; 
!!                              LSD 0; a  b; aa bb; aaa bbb;
!!
!!   AUTHOR
!!     JGH (17.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!     fawzi (04.2004)  : adapted to the new xc interface
!!
!!*** *************************************************************************

MODULE xc_xalpha_new
! *****************************************************************************
  USE cp_array_r_utils,                ONLY: cp_3d_r_p_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_functionals_utilities,        ONLY: set_util
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp
  REAL(KIND=dp), PARAMETER :: f13 = 1.0_dp/3.0_dp, &
                          f23 = 2.0_dp*f13, &
                          f43 = 4.0_dp*f13

  PUBLIC :: xalpha_info, xalpha_lda_eval, xalpha_lsd_eval
  PUBLIC :: xalpha_init, xalpha_lda_0 ! use of theese should be avoided

  REAL(KIND=dp) :: xparam, flda, flsd
  REAL(KIND=dp) :: eps_rho
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='xc_xalpha_new'

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE xalpha_init ( cutoff, xalpha )

    REAL(KIND=dp), INTENT(IN)                    :: cutoff
    REAL(KIND=dp), INTENT(IN), OPTIONAL          :: xalpha

    eps_rho = cutoff
    CALL set_util ( cutoff )
    IF ( PRESENT ( xalpha ) ) THEN
      xparam = xalpha
    ELSE
      xparam = 2.0_dp / 3.0_dp
    END IF

    flda = -9.0_dp/8.0_dp * xparam * (3.0_dp/pi)**f13
    flsd = flda * 2.0_dp**f13

  END SUBROUTINE xalpha_init

! *****************************************************************************

  SUBROUTINE xalpha_info ( lsd, reference, shortform, needs, max_deriv, &
       xa_parameter,scaling,error)
    LOGICAL, INTENT(in)                      :: lsd
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    REAL(KIND=dp), INTENT(in), OPTIONAL     :: xa_parameter,scaling
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    REAL(KIND=dp)                           :: my_xparam, my_scaling

    my_xparam=2.0_dp/3.0_dp
    IF (PRESENT(xa_parameter)) my_xparam=xa_parameter
    my_scaling=1.0_dp
    if (present(scaling)) my_scaling=scaling

    IF ( PRESENT ( reference ) ) THEN
       IF (my_scaling/=1._dp) THEN
          WRITE (reference,'(A,F8.4,A,F8.4)') &
         "Dirac/Slater local exchange; parameter=",my_xparam," scaling=",my_scaling
else
       WRITE (reference,'(A,F8.4)') &
         "Dirac/Slater local exchange; parameter=",my_xparam
    END IF
       IF (.not.lsd) THEN
          IF (LEN_TRIM(reference)+6<LEN(reference)) THEN
             reference(LEN_TRIM(reference):LEN_TRIM(reference)+6)=' {LDA}'
          END IF
       END IF
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       IF (my_scaling/=1._dp) THEN
          WRITE (shortform,'(A,F8.4,F8.4)') "Dirac/Slater exchange",my_xparam,my_scaling
       ELSE
          WRITE (shortform,'(A,F8.4)') "Dirac/Slater exchange",my_xparam
       END IF
       IF (.NOT.lsd) THEN
          IF (LEN_TRIM(shortform)+6<LEN(shortform)) THEN
             shortform(LEN_TRIM(shortform):LEN_TRIM(shortform)+6)=' {LDA}'
          END IF
       END IF
    END IF
    IF (PRESENT(needs)) THEN
       IF (lsd) THEN
          needs%rho_spin=.TRUE.
          needs%rho_spin_1_3=.TRUE.
       ELSE
          needs%rho=.TRUE.
          needs%rho_1_3=.TRUE.
       END IF
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE xalpha_info

! *****************************************************************************

  SUBROUTINE xalpha_lda_eval(rho_set,deriv_set,order,xa_parameter,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: order
    REAL(KIND=dp), INTENT(in), OPTIONAL     :: xa_parameter
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xalpha_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                           :: epsilon_rho
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: e_0, e_rho, e_rho_rho, &
                                                e_rho_rho_rho, r13, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(bo)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rho_1_3=r13,rho=rho,&
            local_bounds=bo,rho_cutoff=epsilon_rho,&
            error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)
       CALL xalpha_init(epsilon_rho,xa_parameter)

       IF ( order>=0 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE., error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
          
          CALL xalpha_lda_0 ( npoints, rho, r13, e_0 )
       END IF
       IF ( order>=1.OR.order==-1 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)

          CALL xalpha_lda_1 ( npoints, rho, r13, e_rho )
       END IF
       IF ( order>=2.OR.order==-2 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)

          CALL xalpha_lda_2 ( npoints, rho, r13, e_rho_rho )
       END IF
       IF ( order>=3.OR.order==-3 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)

          CALL xalpha_lda_3 ( npoints, rho, r13, e_rho_rho_rho )
       END IF
       IF ( order>3.OR.order<-3) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 3 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE xalpha_lda_eval

! *****************************************************************************

  SUBROUTINE xalpha_lsd_eval (rho_set,deriv_set,order,xa_parameter,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: order
    REAL(KIND=dp), INTENT(in), OPTIONAL     :: xa_parameter
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xalpha_lsd_eval', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=6), DIMENSION(2) :: rho_spin_name = (/"(rhoa)","(rhob)"/)
    INTEGER                                  :: handle, i, ispin, npoints
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER   :: e_0, e_rho, e_rho_rho, &
                                                e_rho_rho_rho
    REAL(KIND=dp)                           :: epsilon_rho
    TYPE(cp_3d_r_p_type), DIMENSION(2)       :: rho, rho_1_3
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(deriv, bo)
    DO i=1,2
       NULLIFY(rho(i)%array, rho_1_3(i)%array)
    END DO

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rhoa_1_3=rho_1_3(1)%array,&
            rhob_1_3=rho_1_3(2)%array,rhoa=rho(1)%array,&
            rhob=rho(2)%array, rho_cutoff=epsilon_rho,&
            local_bounds=bo, error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)
       CALL xalpha_init(epsilon_rho,xa_parameter)

       DO ispin=1,2          
          IF ( order>=0 ) THEN
             deriv => xc_dset_get_derivative(deriv_set,"",&
                  allocate_deriv=.TRUE., error=error)
             CALL xc_derivative_get(deriv, deriv_data=e_0,error=error)

             CALL xalpha_lsd_0 ( npoints,rho(ispin)%array, rho_1_3(ispin)%array,&
                  e_0 )
          END IF
          IF ( order>=1.OR.order==-1 ) THEN
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin),&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)

             CALL xalpha_lsd_1( npoints, rho(ispin)%array, rho_1_3(ispin)%array,&
                  e_rho )
          END IF
          IF ( order>=2.OR.order==-2 ) THEN
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin)//&
                  rho_spin_name(ispin),allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)
             
             CALL xalpha_lsd_2( npoints, rho(ispin)%array, rho_1_3(ispin)%array,&
                  e_rho_rho )
          END IF
          IF ( order>=3 .OR. order==-3 ) THEN
             deriv => xc_dset_get_derivative(deriv_set,rho_spin_name(ispin)//&
                  rho_spin_name(ispin)//rho_spin_name(ispin),&
                  allocate_deriv=.TRUE.,error=error)
             CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)

             CALL xalpha_lsd_3( npoints, rho(ispin)%array, rho_1_3(ispin)%array,&
                  e_rho_rho_rho )
          END IF
          IF ( order>3.OR.order<-3) THEN
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="derivatives bigger than 3 not implemented", &
                  error=error, error_level=cp_failure_level)
          END IF
       END DO
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE xalpha_lsd_eval

! *****************************************************************************

  SUBROUTINE xalpha_lda_0(n, rho, r13, pot)

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)      :: rho, r13
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: pot

    INTEGER                                  :: ip

!   ---------------------------------------------------------------------------

!$omp parallel do private(ip)

    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + flda*r13(ip)*rho(ip)
      END IF
    END DO

  END SUBROUTINE xalpha_lda_0

! *****************************************************************************

  SUBROUTINE xalpha_lda_1(n, rho, r13, pot)
    
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)      :: rho, r13
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: pot

    INTEGER                                  :: ip
    REAL(KIND=dp)                                :: f

!   ---------------------------------------------------------------------------

    f = f43 * flda

!$omp parallel do private(ip)
    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + f * r13(ip)
      END IF
    END DO

  END SUBROUTINE xalpha_lda_1

! *****************************************************************************

  SUBROUTINE xalpha_lda_2(n, rho, r13, pot)

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)      :: rho, r13
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: pot

    INTEGER                                  :: ip
    REAL(KIND=dp)                                :: f

!   ---------------------------------------------------------------------------

    f = f13 * f43 * flda

!$omp parallel do private(ip)
    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + f * r13(ip) / rho(ip)
      END IF
    END DO

  END SUBROUTINE xalpha_lda_2

! *****************************************************************************

  SUBROUTINE xalpha_lda_3(n, rho, r13, pot)

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)      :: rho, r13
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: pot

    INTEGER                                  :: ip
    REAL(KIND=dp)                                :: f

!   ---------------------------------------------------------------------------

    f = -f23 * f13 * f43 * flda

!$omp parallel do private(ip)
    DO ip = 1, n
      IF ( rho(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + f * r13(ip) / ( rho(ip) * rho(ip) )
      END IF
    END DO

  END SUBROUTINE xalpha_lda_3

! *****************************************************************************

  SUBROUTINE xalpha_lsd_0 ( n, rhoa, r13a, pot )

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)      :: rhoa, r13a
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: pot

    INTEGER                                  :: ip

! number of points in array 
!   ---------------------------------------------------------------------------

!$omp parallel do private(ip)
    DO ip = 1, n

      IF ( rhoa(ip) > eps_rho ) THEN
         pot(ip) = pot(ip) + flsd * r13a(ip) * rhoa(ip)
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_0

! *****************************************************************************

  SUBROUTINE xalpha_lsd_1 ( n, rhoa, r13a, pota )

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)      :: rhoa, r13a
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: pota

    INTEGER                                  :: ip
    REAL(KIND=dp)                                :: f

! number of points in array 
!   ---------------------------------------------------------------------------

    f = f43 * flsd

!$omp parallel do private(ip)
    DO ip = 1, n

      IF ( rhoa(ip) > eps_rho ) THEN
         pota(ip) = pota(ip) + f * r13a(ip)
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_1

! *****************************************************************************

  SUBROUTINE xalpha_lsd_2 ( n, rhoa, r13a, potaa )

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)      :: rhoa, r13a
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: potaa

    INTEGER                                  :: ip
    REAL(KIND=dp)                                :: f

! number of points in array 
!   ---------------------------------------------------------------------------

    f = f13 * f43 * flsd

!$omp parallel do private(ip)
    DO ip = 1, n

      IF ( rhoa(ip) > eps_rho ) THEN
         potaa(ip) = potaa(ip) + f * r13a(ip)/rhoa(ip)
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_2

! *****************************************************************************

  SUBROUTINE xalpha_lsd_3 ( n, rhoa, r13a, potaaa )

    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)      :: rhoa, r13a
    REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)   :: potaaa

    INTEGER                                  :: ip
    REAL(KIND=dp)                                :: f

! number of points in array 
!   ---------------------------------------------------------------------------

    f = -f23 * f13 * f43 * flsd

!$omp parallel do private(ip)
    DO ip = 1, n

      IF ( rhoa(ip) > eps_rho ) THEN
         potaaa(ip) = potaaa(ip) + f * r13a(ip)/(rhoa(ip)*rhoa(ip))
      END IF

    END DO

  END SUBROUTINE xalpha_lsd_3

! *****************************************************************************

END MODULE xc_xalpha_new

! *****************************************************************************
