!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the Perdew-Wang correlation potential and
!>      energy density and ist derivatives with respect to
!>      the spin-up and spin-down densities up to 3rd order.
!> \par History
!>      18-MAR-2002, TCH, working version
!>      fawzi (04.2004)  : adapted to the new xc interface
!> \see functionals_utilities
! *****************************************************************************
MODULE xc_perdew_wang
  
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_functionals_utilities,        ONLY: calc_fx,&
                                             calc_rs,&
                                             calc_z,&
                                             set_util
  USE xc_input_constants,              ONLY: pw_dmc,&
                                             pw_orig,&
                                             pw_vmc
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: initialized = .FALSE.
  REAL(KIND=dp), DIMENSION(-1:1) :: A, a1, b1, b2, b3, b4
  REAL(KIND=dp), DIMENSION(-1:1) :: c0, c1, c2, c3
  REAL(KIND=dp), DIMENSION(-1:1) :: d0, d1
  REAL(KIND=dp), PARAMETER :: &
       epsilon = 5.E-13_dp, &
       fpp = 0.584822362263464620726223866376013788782_dp ! d^2f(0)/dz^2
  REAL(KIND=dp) :: eps_rho
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_perdew_wang'

  PUBLIC :: perdew_wang_info, perdew_wang_lda_eval, perdew_wang_lsd_eval

CONTAINS

! *****************************************************************************
!> \brief Return some info on the functionals.
!> \param method ...
!> \param lsd ...
!> \param reference full reference
!> \param shortform short reference
!> \param needs ...
!> \param max_deriv ...
!> \param scale ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE perdew_wang_info ( method, lsd, reference, shortform, needs, &
       max_deriv, scale, error)
    INTEGER, INTENT(in)                      :: method
    LOGICAL, INTENT(in)                      :: lsd
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    REAL(kind=dp), INTENT(in)                :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'perdew_wang_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=3)                         :: p_string

    SELECT CASE (method)
    CASE DEFAULT
       CPAssertNoFail(.FALSE.,cp_failure_level,routineP,error)
    CASE (pw_orig)
       p_string='PWO'
    CASE (pw_dmc)
       p_string='DMC'
    CASE (pw_vmc)
       p_string='VMC'
    END SELECT

    IF ( PRESENT ( reference ) ) THEN
       reference = "J. P. Perdew and Yue Wang," &
                   //" Phys. Rev. B 45, 13244 (1992)"&
                   //"["//TRIM(p_string)//"]"
       IF (scale/=1._dp) THEN
          WRITE(reference(LEN_TRIM(reference)+1:LEN(reference)),"('s=',f5.3)")&
               scale
       END IF
       IF (.not.lsd) THEN
          IF (LEN_TRIM(reference)+6<LEN(reference)) THEN
             reference(LEN_TRIM(reference)+1:LEN_TRIM(reference)+7)=' {LDA}'
          END IF
       END IF
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "J. P. Perdew et al., PRB 45, 13244 (1992)"&
                   //"["//TRIM(p_string)//"]"
       IF (scale/=1._dp) THEN
          WRITE(shortform(LEN_TRIM(shortform)+1:LEN(shortform)),"('s=',f5.3)")&
               scale
       END IF
       IF (.not.lsd) THEN
          IF (LEN_TRIM(shortform)+6<LEN(shortform)) THEN
             shortform(LEN_TRIM(shortform)+1:LEN_TRIM(shortform)+7)=' {LDA}'
          END IF
       END IF
    END IF
    IF (PRESENT(needs)) THEN
       IF (lsd) THEN
          needs%rho_spin=.TRUE.
       ELSE
          needs%rho=.TRUE.
       END IF
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE perdew_wang_info

! *****************************************************************************
!> \brief Initializes the functionals
!> \param method name of the method used for parameters
!> \param cutoff the cutoff density
! *****************************************************************************
  SUBROUTINE perdew_wang_init(method, cutoff)

    INTEGER, INTENT(IN)                      :: method
    REAL(KIND=dp), INTENT(IN)                :: cutoff

    CHARACTER(len=*), PARAMETER :: routineN = 'perdew_wang_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k

    CALL set_util(cutoff)

    eps_rho = cutoff

    initialized = .FALSE.

    ! values for -ac are the same for all methods
    A(-1)  = 0.016887_dp
    a1(-1) = 0.11125_dp
    b1(-1) = 10.357_dp
    b2(-1) = 3.6231_dp
    b3(-1) = 0.88026_dp
    b4(-1) = 0.49671_dp

    SELECT CASE (method)

    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unknown method")

    CASE (pw_orig)
       A(0)  = 0.031091_dp ; A(1)  = 0.015545_dp
       a1(0) = 0.21370_dp  ; a1(1) = 0.20548_dp
       b1(0) = 7.5957_dp   ; b1(1) = 14.1189_dp
       b2(0) = 3.5876_dp   ; b2(1) = 6.1977_dp
       b3(0) = 1.6382_dp   ; b3(1) = 3.3662_dp
       b4(0) = 0.49294_dp  ; b4(1) = 0.62517_dp

    CASE (pw_dmc)
       A(0)  = 0.031091_dp ;  A(1)  = 0.015545_dp
       a1(0) = 0.026481_dp ;  a1(1) = 0.022465_dp
       b1(0) = 7.5957_dp   ;  b1(1) = 14.1189_dp
       b2(0) = 3.5876_dp   ;  b2(1) = 6.1977_dp
       b3(0) = -0.46647_dp ;  b3(1) = -0.56043_dp
       b4(0) = 0.13354_dp  ;  b4(1) = 0.11313_dp

    CASE (pw_vmc)
       A(0)  = 0.031091_dp ; A(1)  = 0.015545_dp
       a1(0) = -0.002257_dp; a1(1) = -0.009797_dp
       b1(0) = 7.5957_dp   ; b1(1) = 14.1189_dp
       b2(0) = 3.5876_dp   ; b2(1) = 6.1977_dp
       b3(0) = -0.52669_dp ; b3(1) = -0.91381_dp
       b4(0) = 0.03755_dp  ; b4(1) = 0.01538_dp

    END SELECT

    DO k=-1, 1, 1
       c0(k) = A(k)
       c1(k) = -2.0_dp*c0(k)*LOG(2.0_dp*A(k)*b1(k))
       c2(k) = A(k)*a1(k)
       c3(k) = -2.0_dp*A(k)*(a1(k)*LOG(2.0_dp*A(k)*b1(k)) &
            - (b2(k)/b1(k))**2 + (b3(k)/b1(k)))
       d0(k) = a1(k)/b4(k)
       d1(k) = a1(k)*b3(k)/(b4(k)**2)
    END DO

    initialized = .TRUE.

  END SUBROUTINE perdew_wang_init

! *****************************************************************************
!> \brief Calculate the correlation energy and its derivatives
!>      wrt to rho (the electron density) up to 3rd order. This
!>      is the LDA version of the Perdew-Wang correlation energy
!>      If no order argument is given, then the routine calculates
!>      just the energy.
!> \param method ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order order of derivatives to calculate
!>      order must lie between -3 and 3. If it is negative then only
!>      that order will be calculated, otherwise all derivatives up to
!>      that order will be calculated.
!> \param scale ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE perdew_wang_lda_eval( method, rho_set, deriv_set, order, scale, error )

    INTEGER, INTENT(in)                      :: method
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: order
    REAL(kind=dp), INTENT(in)                :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'perdew_wang_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: npoints, stat, timer_handle
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rho_cutoff
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_rho, e_rho_rho, &
                                                e_rho_rho_rho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,timer_handle)
    failure=.FALSE.
    NULLIFY(bo,rho, e_0,e_rho,e_rho_rho,e_rho_rho_rho, dummy)
    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL xc_rho_set_get(rho_set,rho=rho,&
            local_bounds=bo,rho_cutoff=rho_cutoff,error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

       CALL perdew_wang_init(method, rho_cutoff)

       ! meaningful default for the arrays we don't need: let us make compiler
       ! and debugger happy...
       IF (cp_debug) THEN
          ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          dummy=> rho
       END IF

       e_0 => dummy
       e_rho => dummy
       e_rho_rho => dummy
       e_rho_rho_rho => dummy

       IF (order>=0) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE., error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
       END IF
       IF (order>=1.OR.order==-1) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)
       END IF
       IF (order>=2.OR.order==-2) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)
       END IF
       IF (order>=3.OR.order==-3) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)
       END IF
       IF (order>3.OR.order<-3) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 3 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       CALL perdew_wang_lda_calc(rho,e_0,e_rho, e_rho_rho, e_rho_rho_rho, &
            npoints, order,scale)

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF

    CALL timestop(timer_handle)

  END SUBROUTINE perdew_wang_lda_eval

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_rho_rho ...
!> \param e_rho_rho_rho ...
!> \param npoints ...
!> \param order ...
!> \param scale ...
! *****************************************************************************
SUBROUTINE perdew_wang_lda_calc(rho,e_0,e_rho, e_rho_rho, e_rho_rho_rho, npoints, order,scale)
  !FM low level calc routine
    REAL(KIND=dp), DIMENSION(*), INTENT(in)  :: rho
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(inout)                          :: e_0, e_rho, e_rho_rho, &
                                                e_rho_rho_rho
    INTEGER, INTENT(in)                      :: npoints, order
    REAL(kind=dp), INTENT(in)                :: scale

    INTEGER                                  :: abs_order, k
    REAL(KIND=dp), DIMENSION(0:3)            :: ed

  abs_order=ABS(order)

!$omp parallel do private (k, ed) default(none)&
!$omp shared(npoints,rho,eps_rho,abs_order,scale,e_0,e_rho,e_rho_rho,e_rho_rho_rho,order)
    DO k=1, npoints

       IF ( rho(k) > eps_rho ) THEN
!! order_ is positive as it must be in this case:
!! ec(:,2) needs ed(:,1) for example
         CALL pw_lda_ed_loc(rho(k), ed, abs_order)
         ed=scale*ed

         IF (order>=0) THEN
            e_0(k) = e_0(k) + rho(k)*ed(0)
         END IF
         IF (order>=1.OR.order==-1) THEN
            e_rho(k) = e_rho(k) + ed(0) + rho(k)*ed(1)
         END IF
         IF (order>=2.OR.order==-2) THEN
            e_rho_rho(k) = e_rho_rho(k) + 2.0_dp*ed(1) + rho(k)*ed(2)
         END IF
         IF (order>=3.OR.order==-3) THEN
            e_rho_rho_rho(k) = e_rho_rho_rho(k) + 3.0_dp*ed(2) + rho(k)*ed(3)
         END IF

      END IF

    END DO
!$omp end parallel do

  END SUBROUTINE perdew_wang_lda_calc

! *****************************************************************************
!> \brief Calculate the correlation energy and its derivatives
!>      wrt to rho (the electron density) up to 3rd order. This
!>      is the LSD version of the Perdew-Wang correlation energy
!>      If no order argument is given, then the routine calculates
!>      just the energy.
!> \param method ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order order of derivatives to calculate
!>      order must lie between -3 and 3. If it is negative then only
!>      that order will be calculated, otherwise all derivatives up to
!>      that order will be calculated.
!> \param scale ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE perdew_wang_lsd_eval( method, rho_set, deriv_set, order, scale, error )
    INTEGER, INTENT(in)                      :: method
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN), OPTIONAL            :: order
    REAL(kind=dp), INTENT(in)                :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'perdew_wang_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: npoints, stat, timer_handle
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rho_cutoff
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: a, b, dummy, e_0, ea, eaa, &
                                                eaaa, eaab, eab, eabb, eb, &
                                                ebb, ebbb
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,timer_handle)
    failure=.FALSE.
    NULLIFY(bo,a,b,e_0,ea,eb,eaa,eab,ebb,eaaa,eaab,eabb,ebbb)
    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL xc_rho_set_get(rho_set,rhoa=a,rhob=b,&
            local_bounds=bo,rho_cutoff=rho_cutoff,error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)

       CALL perdew_wang_init(method, rho_cutoff)

       ! meaningful default for the arrays we don't need: let us make compiler
       ! and debugger happy...
       IF (cp_debug) THEN
          ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE
          dummy=> a
       END IF

       e_0 => dummy
       ea => dummy; eb => dummy
       eaa => dummy; eab => dummy; ebb => dummy
       eaaa => dummy; eaab=> dummy; eabb => dummy; ebbb => dummy

       IF (order>=0) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE., error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
       END IF
       IF (order>=1.OR.order==-1) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=ea,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=eb,error=error)
       END IF
       IF (order>=2.OR.order==-2) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=eaa,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=eab,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=ebb,error=error)
       END IF
       IF (order>=3.OR.order==-3) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)(rhoa)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=eaaa,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=eaab,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=eabb,error=error)
          deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)(rhob)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=ebbb,error=error)
       END IF
       IF (order>3.OR.order<-3) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 3 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       CALL perdew_wang_lsd_calc(a,b,e_0,ea,eb,eaa,eab,ebb,eaaa,eaab,eabb,&
         ebbb, npoints, order,scale)

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF
    CALL timestop(timer_handle)

  END SUBROUTINE perdew_wang_lsd_eval

! *****************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param rhob ...
!> \param e_0 ...
!> \param ea ...
!> \param eb ...
!> \param eaa ...
!> \param eab ...
!> \param ebb ...
!> \param eaaa ...
!> \param eaab ...
!> \param eabb ...
!> \param ebbb ...
!> \param npoints ...
!> \param order ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE perdew_wang_lsd_calc(rhoa,rhob,e_0,ea,eb,eaa,eab,ebb,eaaa,eaab,eabb,&
       ebbb, npoints, order,scale)
      !FM low-level computation routine
    REAL(KIND=dp), DIMENSION(*), INTENT(in)  :: rhoa, rhob
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(inout)                          :: e_0, ea, eb, eaa, eab, ebb, &
                                                eaaa, eaab, eabb, ebbb
    INTEGER, INTENT(in)                      :: npoints, order
    REAL(kind=dp), INTENT(in)                :: scale

    INTEGER                                  :: abs_order, k
    REAL(KIND=dp)                            :: rho
    REAL(KIND=dp), DIMENSION(0:9)            :: ed

    abs_order=ABS(order)

!$omp parallel do private (k, rho, ed) default(none)&
!$omp shared(npoints,rhoa,rhob,eps_rho,abs_order,order,e_0,ea,eb,eaa,eab,ebb,eaaa,eaab,eabb,ebbb,scale)
    DO k=1, npoints

       rho = rhoa(k) + rhob(k)
       IF ( rho > eps_rho ) THEN

         CALL pw_lsd_ed_loc(rhoa(k), rhob(k), ed, abs_order)
         ed=ed*scale
         IF (order>=0) THEN
            e_0(k) = e_0(k) + rho*ed(0)
         END IF
         IF (order>=1.OR.order==-1) THEN
            ea(k)   = ea(k) + ed(0) + rho*ed(1)
            eb(k) = eb(k) + ed(0) + rho*ed(2)
         END IF
         IF (order>=2.OR.order==-2) THEN
            eaa(k) = eaa(k) + 2.0_dp*ed(1) + rho*ed(3)
            eab(k) = eab(k) + ed(1) + ed(2) + rho*ed(4)
            ebb(k) = ebb(k) + 2.0_dp*ed(2) + rho*ed(5)
         END IF
         IF (order>=3.OR.order==-3) THEN
            eaaa(k) = eaaa(k) + 3.0_dp*ed(3) + rho*ed(6)
            eaab(k) = eaab(k) + 2.0_dp*ed(4) + ed(3) + rho*ed(7)
            eabb(k) = eabb(k) + 2.0_dp*ed(4) + ed(5) + rho*ed(8)
            ebbb(k) = ebbb(k) + 3.0_dp*ed(5) + rho*ed(9)
         END IF

       END IF

    END DO

  END SUBROUTINE perdew_wang_lsd_calc

! *****************************************************************************
!> \brief ...
!> \param r ...
!> \param z ...
!> \param g ...
!> \param order ...
! *****************************************************************************
  SUBROUTINE calc_g(r, z, g, order)

!   Calculates g and its derivatives wrt r up to 3rd order, where:
!
!   g = .... for r < 1
!   g = .... for r > 100 and everywhere else
!   g = 2A(1+a1*r)ln(1+1/(2A(b1*r^1/2 + b2*r + b3*r^(3/2) + b4*r^2))).

    REAL(KIND=dp), INTENT(IN)                :: r
    INTEGER, INTENT(IN)                      :: z
    REAL(KIND=dp), DIMENSION(0:), &
      INTENT(OUT)                            :: g
    INTEGER, INTENT(IN)                      :: order

    REAL(KIND=dp)                            :: a1_, A_, b1_, b2_, b3_, b4_, &
                                                rr, rsr, sr, t11, t12, t14, &
                                                t15, t16, t20, t22, t3, t40, &
                                                t44, t45, t47, t48, t55, t56

    A_ = A(z); a1_ = a1(z)
    b1_ = b1(z); b2_ = b2(z); b3_ = b3(z); b4_ = b4(z)

    sr = SQRT(r)
    rsr = r*sr
    rr = r*r

    IF (r < 1.0_dp) THEN

       ! order 0 must always be calculated
       g(0) = c0(z)*LOG(r) - c1(z) + c2(z)*r*LOG(r) - c3(z)*r
       IF (order >= 1) g(1) = c0(z)/r + c2(z)*LOG(r) + c2(z) - c3(z)
       IF (order >= 2) g(2) = -c0(z)/rr + c2(z)/r
       IF (order >= 3) g(3) = 2.0_dp*c0(z)/(rr*r) - c2(z)/rr

    ELSE IF (r <= 100.0_dp) THEN

       t3 = 1.0_dp+a1_*r
       t11 = b1_*sr + b2_*r + b3_*rsr + b4_*rr
       t12 = t11**2
       t15 = 1.0_dp + 0.5_dp/A_/t11
       t16 = LOG(t15)
       t20 = 0.5_dp*b1_/sr + b2_ + 1.5_dp*b3_*sr + 2.0_dp*b4_*r

       ! order 0 must always be calculated
       g(0) = -2.0_dp*A_*t3*t16

       IF (order >= 1) THEN

          g(1) = -2.0_dp*A_*a1_*t16 + t3*t20/(t12*t15)

       END IF

       IF (order >= 2) THEN

          t40 = -0.25_dp*b1_/rsr + 0.75_dp*b3_/sr + 2.0_dp*b4_

          g(2) = 2.0_dp*a1_*t20/(t12*t15) &
               - 2.0_dp*(t20**2)*t3/(t12*t11*t15) &
               + t3*t40/(t12*t15) &
               + 0.5_dp*t3*(t20**2)/(A_*(t12**2)*(t15**2))

       END IF

       IF (order >= 3) THEN

          t14 = 1.0_dp/t12/t11
          t22 = t20**2
          t56 = t22*t20
          t47 = t15**2
          t48 = 1.0_dp/t47

          t44 = t12**2
          t45 = 1.0_dp/t44
          t55 = t3*t45

          g(3) = &
               - 6.0_dp*a1_*t14*t22/t15 &
               + 3.0_dp*a1_*t40/(t15*t12) &
               + 1.5_dp*a1_*t45*t22*t48/A_ &
               + 6.0_dp*t55*t56/t15 &
               - 6.0_dp*t3*t14*t20*t40/t15 &
               - 3.0_dp*t3*t56*t48/(A_*t44*t11) &
               + 0.375_dp*t3*(b1_/(rr*sr)-b3_/rsr)/(t12*t15) &
               + 1.5_dp*t55*t40*t48*t20/A_ &
               + 0.5_dp*t3*t56/((A_**2)*t44*t12*t47*t15)

       END IF

    ELSE

       ! order 0 must always be calculated
       g(0) = -d0(z)/r + d1(z)/rsr
       IF (order >= 1) g(1) = d0(z)/rr - 1.5_dp*d1(z)/(rsr*r)
       IF (order >= 2) g(2) = -2.0_dp*d0(z)/(rr*r) + 3.75_dp*d1(z)/(rsr*rr)
       IF (order >= 3) g(3) = 6.0_dp*d0(z)/(rr*rr) - 13.125_dp*d1(z)/(rsr*rr*r)

    END IF

  END SUBROUTINE calc_g

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param ed ...
!> \param order ...
! *****************************************************************************
SUBROUTINE pw_lda_ed_loc(rho, ed, order)

    REAL(KIND=dp), INTENT(IN)                :: rho
    REAL(KIND=dp), DIMENSION(0:), &
      INTENT(OUT)                            :: ed
    INTEGER, INTENT(IN)                      :: order

    INTEGER                                  :: m, order_
    LOGICAL, DIMENSION(0:3)                  :: calc
    REAL(KIND=dp), DIMENSION(0:3)            :: e0, r

    order_ = order
    ed = 0
    calc=.FALSE.

    IF (order_ >= 0) THEN
       calc(0:order_) = .TRUE.
    ELSE
       order_ = -1 * order_
       calc(order_) = .TRUE.
    END IF

    CALL calc_rs(rho, r(0))
    CALL calc_g(r(0), 0, e0, order_)

    IF (order_ >= 1) r(1) = (-1.0_dp/3.0_dp)*r(0)/rho
    IF (order_ >= 2) r(2) = (-4.0_dp/3.0_dp)*r(1)/rho
    IF (order_ >= 3) r(3) = (-7.0_dp/3.0_dp)*r(2)/rho

    m = 0
    IF (calc(0)) THEN
       ed(m) = e0(0)
       m = m + 1
    END IF
    IF (calc(1)) THEN
       ed(m) = e0(1)*r(1)
       m =  m + 1
    END IF
    IF (calc(2)) THEN
       ed(m) = e0(2)*r(1)**2 + e0(1)*r(2)
       m = m + 1
    END IF
    IF (calc(3)) THEN
       ed(m) = e0(3)*r(1)**3 + e0(2)*3.0_dp*r(1)*r(2) + e0(1)*r(3)
    END IF

  END SUBROUTINE  pw_lda_ed_loc

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param ed ...
!> \param order ...
! *****************************************************************************
  SUBROUTINE pw_lsd_ed_loc(a, b, ed, order)

    REAL(KIND=dp), INTENT(IN)                :: a, b
    REAL(KIND=dp), DIMENSION(0:), &
      INTENT(OUT)                            :: ed
    INTEGER, INTENT(IN)                      :: order

    INTEGER                                  :: m, order_
    LOGICAL, DIMENSION(0:3)                  :: calc
    REAL(KIND=dp)                            :: rho, tr, trr, trrr, trrz, &
                                                trz, trzz, tz, tzz, tzzz
    REAL(KIND=dp), DIMENSION(0:3)            :: ac, e0, e1, f, r
    REAL(KIND=dp), DIMENSION(0:3, 0:3)       :: z

    order_ = order
    calc = .FALSE.

    IF (order_ > 0) THEN
       calc(0:order_) = .TRUE.
    ELSE
       order_ = -1 * order_
       calc(order_) = .TRUE.
    END IF

    rho = a + b

    CALL calc_fx(a, b, f(0:order_), order_)
    CALL calc_rs(rho, r(0))
    CALL calc_g(r(0), -1, ac(0:order_), order_)
    CALL calc_g(r(0), 0, e0(0:order_), order_)
    CALL calc_g(r(0), 1, e1(0:order_), order_)
    CALL calc_z(a, b, z(0:order_,0:order_), order_)

!! calculate first partial derivatives
    IF (order_ >= 1) THEN
       r(1) = (-1.0_dp/3.0_dp)*r(0)/rho
       tr = e0(1) &
            + fpp*ac(1)*f(0) &
            - fpp*ac(1)*f(0)*z(0,0)**4 &
            + (e1(1)-e0(1))*f(0)*z(0,0)**4
       tz = fpp*ac(0)*f(1) &
            - fpp*ac(0)*f(1)*z(0,0)**4 &
            - fpp*ac(0)*f(0)*4.0_dp*z(0,0)**3 &
            + (e1(0)-e0(0))*f(1)*z(0,0)**4 &
            + (e1(0)-e0(0))*f(0)*4.0_dp*z(0,0)**3
    END IF

!! calculate second partial derivatives
    IF (order_ >= 2) THEN
       r(2) = (-4.0_dp/3.0_dp)*r(1)/rho
       trr = e0(2) &
            + fpp*ac(2)*f(0) &
            - fpp*ac(2)*f(0)*z(0,0)**4 &
            + (e1(2)-e0(2))*f(0)*z(0,0)**4
       trz = fpp*ac(1)*f(1) &
            - fpp*ac(1)*f(1)*z(0,0)**4 &
            - fpp*ac(1)*f(0)*4.0_dp*z(0,0)**3 &
            + (e1(1)-e0(1))*f(1)*z(0,0)**4 &
            + (e1(1)-e0(1))*f(0)*4.0_dp*z(0,0)**3
       tzz = fpp*ac(0)*f(2) &
            - fpp*ac(0)*f(2)*z(0,0)**4 &
            - fpp*ac(0)*f(1)*8.0_dp*z(0,0)**3 &
            - fpp*ac(0)*f(0)*12.0_dp*z(0,0)**2 &
            + (e1(0)-e0(0))*f(2)*z(0,0)**4 &
            + (e1(0)-e0(0))*f(1)*8.0_dp*z(0,0)**3 &
            + (e1(0)-e0(0))*f(0)*12.0_dp*z(0,0)**2
    END IF

!! calculate third derivatives
    IF (order_ >= 3) THEN

       r(3) = (-7.0_dp/3.0_dp)*r(2)/rho

       trrr = e0(3) &
            + fpp*ac(3)*f(0) &
            - fpp*ac(3)*f(0)*z(0,0)**4 &
            + (e1(3)-e0(3))*f(0)*z(0,0)**4

       trrz = fpp*ac(2)*f(1) &
            - fpp*ac(2)*f(1)*z(0,0)**4 &
            - fpp*ac(2)*f(0)*4.0_dp*z(0,0)**3 &
            + (e1(2)-e0(2))*f(1)*z(0,0)**4 &
            + (e1(2)-e0(2))*f(0)*4.0_dp*z(0,0)**3

       trzz = fpp*ac(1)*f(2) &
            - fpp*ac(1)*f(2)*z(0,0)**4 &
            - fpp*ac(1)*f(1)*8.0_dp*z(0,0)**3 &
            - fpp*ac(1)*f(0)*12.0_dp*z(0,0)**2 &
            + (e1(1)-e0(1))*f(2)*z(0,0)**4 &
            + (e1(1)-e0(1))*f(1)*8.0_dp*z(0,0)**3 &
            + (e1(1)-e0(1))*f(0)*12.0_dp*z(0,0)**2

       tzzz = fpp*ac(0)*f(3) &
            - fpp*ac(0)*f(3)*z(0,0)**4 &
            - fpp*ac(0)*f(2)*12.0_dp*z(0,0)**3 &
            - fpp*ac(0)*f(1)*36.0_dp*z(0,0)**2 &
            - fpp*ac(0)*f(0)*24.0_dp*z(0,0) &
            + (e1(0)-e0(0))*f(3)*z(0,0)**4 &
            + (e1(0)-e0(0))*f(2)*12.0_dp*z(0,0)**3 &
            + (e1(0)-e0(0))*f(1)*36.0_dp*z(0,0)**2 &
            + (e1(0)-e0(0))*f(0)*24.0_dp*z(0,0)
    END IF

    m=0
    IF (calc(0)) THEN
       ed(m) = e0(0) &
            + fpp*ac(0)*f(0)*(1.0_dp - z(0,0)**4) &
            + (e1(0)-e0(0))*f(0)*z(0,0)**4
       m = m + 1
    END IF
    IF (calc(1)) THEN
       ed(m)   = tr*r(1) + tz*z(1,0)
       ed(m+1) = tr*r(1) + tz*z(0,1)
       m = m + 2
    END IF
    IF (calc(2)) THEN
       ed(m)   = trr*r(1)**2 + 2.0_dp*trz*r(1)*z(1,0) &
            + tr*r(2) + tzz*z(1,0)**2 + tz*z(2,0)
       ed(m+1) = trr*r(1)**2 + trz*r(1)*(z(0,1)+z(1,0)) &
            + tr*r(2) + tzz*z(1,0)*z(0,1) + tz*z(1,1)
       ed(m+2) = trr*r(1)**2 + 2.0_dp*trz*r(1)*z(0,1) &
            + tr*r(2) + tzz*z(0,1)**2 + tz*z(0,2)
       m = m + 3
    END IF
    IF (calc(3)) THEN
       ed(m)   = &
            trrr*r(1)**3 + 3.0_dp*trrz*r(1)**2*z(1,0) &
            + 3.0_dp*trr*r(1)*r(2) + 3.0_dp*trz*r(2)*z(1,0) + tr*r(3) &
            + 3.0_dp*trzz*r(1)*z(1,0)**2 + tzzz*z(1,0)**3 &
            + 3.0_dp*trz*r(1)*z(2,0) &
            + 3.0_dp*tzz*z(1,0)*z(2,0) + tz*z(3,0)
       ed(m+1) = &
            trrr*r(1)**3 + trrz*r(1)**2*(2.0_dp*z(1,0)+z(0,1)) &
            + 2.0_dp*trzz*r(1)*z(1,0)*z(0,1) &
            + 2.0_dp*trz*(r(2)*z(1,0)+r(1)*z(1,1)) &
            + 3.0_dp*trr*r(2)*r(1) + trz*r(2)*z(0,1) + tr*r(3) &
            + trzz*r(1)*z(1,0)**2 + tzzz*z(1,0)**2*z(0,1) &
            + 2.0_dp*tzz*z(1,0)*z(1,1) &
            + trz*r(1)*z(2,0) + tzz*z(2,0)*z(0,1) + tz*z(2,1)
       ed(m+2) = &
            trrr*r(1)**3 + trrz*r(1)**2*(2.0_dp*z(0,1)+z(1,0)) &
            + 2.0_dp*trzz*r(1)*z(0,1)*z(1,0) &
            + 2.0_dp*trz*(r(2)*z(0,1)+r(1)*z(1,1)) &
            + 3.0_dp*trr*r(2)*r(1) + trz*r(2)*z(1,0) + tr*r(3) &
            + trzz*r(1)*z(0,1)**2 + tzzz*z(0,1)**2*z(1,0) &
            + 2.0_dp*tzz*z(0,1)*z(1,1) &
            + trz*r(1)*z(0,2) + tzz*z(0,2)*z(1,0) + tz*z(1,2)
       ed(m+3) = &
            trrr*r(1)**3 + 3.0_dp*trrz*r(1)**2*z(0,1) &
            + 3.0_dp*trr*r(1)*r(2) + 3.0_dp*trz*r(2)*z(0,1) + tr*r(3) &
            + 3.0_dp*trzz*r(1)*z(0,1)**2 + tzzz*z(0,1)**3 &
            + 3.0_dp*trz*r(1)*z(0,2) &
            + 3.0_dp*tzz*z(0,1)*z(0,2) + tz*z(0,3)
    END IF

  END SUBROUTINE pw_lsd_ed_loc

END MODULE xc_perdew_wang
