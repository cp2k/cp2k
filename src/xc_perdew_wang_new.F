!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****m* cp2k/xc_perdew_wang *
!!
!!   NAME
!!     xc_perdew_wang
!!
!!   FUNCTION
!!     Calculate the Perdew-Wang correlation potential and 
!!     energy density and ist derivatives with respect to
!!     the spin-up and spin-down densities up to 3rd order.
!!
!!   MODIFICATION HISTORY
!!     18-MAR-2002, TCH, working version
!!     fawzi (04.2004)  : adapted to the new xc interface
!!
!!   SEE ALSO
!!     functionals_utilities
!!
!!****

MODULE xc_perdew_wang_new
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_debug,&
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
  USE xc_functionals_utilities,        ONLY: calc_fx,&
                                             calc_rs,&
                                             calc_z,&
                                             set_util
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
                                             xc_rho_set_get,&
                                             xc_rho_set_type

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: initialized = .FALSE.
  REAL(DBL), DIMENSION(-1:1) :: A, a1, b1, b2, b3, b4
  REAL(DBL), DIMENSION(-1:1) :: c0, c1, c2, c3
  REAL(DBL), DIMENSION(-1:1) :: d0, d1
  REAL(DBL), PARAMETER :: &
       epsilon = 5.E-13_dbl, &
       fpp = 0.584822362263464620726223866376013788782_dbl ! d^2f(0)/dz^2
  REAL(DBL) :: eps_rho
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="xc_perdew_wang_new"


  PUBLIC :: perdew_wang_info, perdew_wang_lda_eval, perdew_wang_lsd_eval

CONTAINS

! *****************************************************************************

!!****f* cp2k/xc_perdew_wang/perdew_wang_info *
!!
!!   NAME
!!     perdew_wang_info
!!
!!   SYNOPSIS
!!     Subroutine perdew_wang_info(method, lsd, reference, shortform, needs,&
!!         max_deriv, error)
!!       Character(Len=*), Intent (IN):: method
!!       Logical, Intent (IN):: lsd
!!       Character(Len=*), Intent (OUT), Optional:: reference, shortform
!!       Type(xc_rho_cflags_type), Intent (INOUT), Optional:: needs
!!       Integer, Intent (OUT), Optional:: max_deriv
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine perdew_wang_info
!!
!!   FUNCTION
!!     Return some info on the functionals.
!!
!!   ARGUMENTS
!!     reference - CHARACTER(*), INTENT(OUT), OPTIONAL - full reference
!!     shortform - CHARACTER(*), INTENT(OUT), OPTIONAL - short reference
!!
!!****
  SUBROUTINE perdew_wang_info ( method, lsd, reference, shortform, needs, &
       max_deriv, error)
    CHARACTER(len=*), INTENT(in)             :: method
    LOGICAL, INTENT(in)                      :: lsd
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "perdew_wang_info", &
      routineP = moduleN//":"//routineN

    SELECT CASE (method)
    CASE DEFAULT
       CPAssertNoFail(.FALSE.,cp_failure_level,routineP,error)
    CASE ('PWO')
    CASE ('DMC')
    CASE ('VMC') 
    END SELECT

    IF ( PRESENT ( reference ) ) THEN
       reference = "J. P. Perdew and Yue Wang," &
                   //" Phys. Rev. B 45, 13244 (1992)"&
                   //"["//TRIM(method)//"]"
       IF (.not.lsd) THEN
          IF (LEN_TRIM(reference)+6<LEN(reference)) THEN
             reference(LEN_TRIM(reference):LEN_TRIM(reference)+6)=' {LDA}'
          END IF
       END IF
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "J. P. Perdew et al., PRB 45, 13244 (1992)"&
                   //"["//TRIM(method)//"]"
       IF (.not.lsd) THEN
          IF (LEN_TRIM(shortform)+6<LEN(shortform)) THEN
             shortform(LEN_TRIM(shortform):LEN_TRIM(shortform)+6)=' {LDA}'
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

!!****f* cp2k/xc_perdew_wang/perdew_wang_init *
!!
!!   NAME
!!     perdew_wang
!!
!!   FUNCTION
!!     Initializes the functionals
!!
!!   ARGUMENTS
!!     method - CHARACTER(3) - name of the method used for parameters
!!     cutoff - REAL(DBL) - the cutoff density
!!
!!****
  SUBROUTINE perdew_wang_init(method, cutoff)
    

    CHARACTER(len=*), INTENT(IN)                 :: method
    REAL(DBL), INTENT(IN)                    :: cutoff

    INTEGER                                  :: k

    CALL set_util(cutoff)

    eps_rho = cutoff
    
    initialized = .FALSE.

    ! values for -ac are the same for all methods
    A(-1)  = 0.016887_dbl
    a1(-1) = 0.11125_dbl
    b1(-1) = 10.357_dbl 
    b2(-1) = 3.6231_dbl
    b3(-1) = 0.88026_dbl
    b4(-1) = 0.49671_dbl

    SELECT CASE (method)

    CASE DEFAULT
       CALL stop_program("xc_perdew_wang/pw_init", "unknown method")

    CASE ('PWO')
       A(0)  = 0.031091_dbl ; A(1)  = 0.015545_dbl
       a1(0) = 0.21370_dbl  ; a1(1) = 0.20548_dbl
       b1(0) = 7.5957_dbl   ; b1(1) = 14.1189_dbl
       b2(0) = 3.5876_dbl   ; b2(1) = 6.1977_dbl
       b3(0) = 1.6382_dbl   ; b3(1) = 3.3662_dbl
       b4(0) = 0.49294_dbl  ; b4(1) = 0.62517_dbl

    CASE ('DMC')
       A(0)  = 0.031091_dbl ;  A(1)  = 0.015545_dbl
       a1(0) = 0.026481_dbl ;  a1(1) = 0.022465_dbl
       b1(0) = 7.5957_dbl   ;  b1(1) = 14.1189_dbl
       b2(0) = 3.5876_dbl   ;  b2(1) = 6.1977_dbl
       b3(0) = -0.46647_dbl ;  b3(1) = -0.56043_dbl
       b4(0) = 0.13354_dbl  ;  b4(1) = 0.11313_dbl
       
    CASE ('VMC') 
       A(0)  = 0.031091_dbl ; A(1)  = 0.015545_dbl
       a1(0) = -0.002257_dbl; a1(1) = -0.009797_dbl
       b1(0) = 7.5957_dbl   ; b1(1) = 14.1189_dbl
       b2(0) = 3.5876_dbl   ; b2(1) = 6.1977_dbl
       b3(0) = -0.52669_dbl ; b3(1) = -0.91381_dbl
       b4(0) = 0.03755_dbl  ; b4(1) = 0.01538_dbl
       
    END SELECT
          
    DO k=-1, 1, 1
       c0(k) = A(k)
       c1(k) = -2._dbl*c0(k)*LOG(2._dbl*A(k)*b1(k))
       c2(k) = A(k)*a1(k)
       c3(k) = -2._dbl*A(k)*(a1(k)*LOG(2._dbl*A(k)*b1(k)) &
            - (b2(k)/b1(k))**2 + (b3(k)/b1(k)))
       d0(k) = a1(k)/b4(k)
       d1(k) = a1(k)*b3(k)/(b4(k)**2)
    END DO

    initialized = .TRUE.

  END SUBROUTINE perdew_wang_init
  
! *****************************************************************************

!!****f* cp2k/xc_perdew_wang/perdew_wang_lda_eval *
!!
!!   NAME
!!     perdew_wang_lda_eval
!!
!!   SYNOPSIS
!!     Subroutine perdew_wang_lda_eval(method, rho_set, deriv_set, order,&
!!         error)
!!       Character(Len=*), Intent (IN):: method
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Type(xc_derivative_set_type), Pointer:: deriv_set
!!       Integer, Intent (IN):: order
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine perdew_wang_lda_eval
!!
!!   FUNCTION
!!     Calculate the correlation energy and its derivatives
!!     wrt to rho (the electron density) up to 3rd order. This
!!     is the LDA version of the Perdew-Wang correlation energy
!!     If no order argument is given, then the routine calculates
!!     just the energy.
!!
!!   ARGUMENTS
!!     order - INTEGER - order of derivatives to calculate
!!
!!     order must lie between -3 and 3. If it is negative then only
!!     that order will be calculated, otherwise all derivatives up to
!!     that order will be calculated.
!!
!!****
  SUBROUTINE perdew_wang_lda_eval( method, rho_set, deriv_set, order, error )
    
    CHARACTER(len=*), INTENT(in)             :: method
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: order
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "perdew_wang_lda_eval", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: npoints, stat, timer_handle
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(dbl)                                :: rho_cutoff
    REAL(dbl), DIMENSION(:, :, :), POINTER   :: dummy, e_0, e_rho, e_rho_rho, &
                                                e_rho_rho_rho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN, 'I', ' ', timer_handle)
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
            npoints, order)

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF

    CALL timestop(0.0_dbl, timer_handle)
    
  END SUBROUTINE perdew_wang_lda_eval

! *****************************************************************************

SUBROUTINE perdew_wang_lda_calc(rho,e_0,e_rho, e_rho_rho, e_rho_rho_rho, npoints, order)
  !FM low level calc routine
    REAL(dbl), DIMENSION(*), INTENT(in)      :: rho
    REAL(dbl), DIMENSION(*), INTENT(inout)   :: e_0, e_rho, e_rho_rho, &
                                                e_rho_rho_rho
    INTEGER, INTENT(in)                      :: npoints, order

    INTEGER                                  :: abs_order, k
    REAL(DBL), DIMENSION(0:3)                :: ed

  abs_order=ABS(order)

!$omp parallel do private (k, ed)
    DO k=1, npoints
       
       IF ( rho(k) > eps_rho ) THEN
!! order_ is positive as it must be in this case:
!! ec(:,2) needs ed(:,1) for example
         CALL pw_lda_ed_loc(rho(k), ed, abs_order)
  
         IF (order>=0) THEN
            e_0(k) = e_0(k) + rho(k)*ed(0)
         END IF
         IF (order>=1.OR.order==-1) THEN
            e_rho(k) = e_rho(k) + ed(0) + rho(k)*ed(1)
         END IF
         IF (order>=2.OR.order==-2) THEN
            e_rho_rho(k) = e_rho_rho(k) + 2._dbl*ed(1) + rho(k)*ed(2)
         END IF
         IF (order>=3.OR.order==-3) THEN
            e_rho_rho_rho(k) = e_rho_rho_rho(k) + 3._dbl*ed(2) + rho(k)*ed(3)
         END IF

      END IF
    
    END DO
!$omp end parallel do
    
  END SUBROUTINE perdew_wang_lda_calc
  
! *****************************************************************************

!!****f* cp2k/xc_perdew_wang/perdew_wang_lsd_eval *
!!
!!   NAME
!!     perdew_wang_lsd_eval
!!
!!   SYNOPSIS
!!     Subroutine perdew_wang_lsd_eval(method, rho_set, deriv_set, order,&
!!         error)
!!       Character(Len=*), Intent (IN):: method
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Type(xc_derivative_set_type), Pointer:: deriv_set
!!       Integer, Intent (IN), Optional:: order
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine perdew_wang_lsd_eval
!!
!!   FUNCTION
!!     Calculate the correlation energy and its derivatives
!!     wrt to rho (the electron density) up to 3rd order. This
!!     is the LSD version of the Perdew-Wang correlation energy
!!     If no order argument is given, then the routine calculates
!!     just the energy.
!!
!!   ARGUMENTS
!!     order - INTEGER, OPTIONAL - order of derivatives to calculate
!!
!!     order must lie between -3 and 3. If it is negative then only
!!     that order will be calculated, otherwise all derivatives up to
!!     that order will be calculated.
!!
!!****
  SUBROUTINE perdew_wang_lsd_eval( method, rho_set, deriv_set, order, error )
    CHARACTER(len=*), INTENT(in)             :: method
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN), OPTIONAL            :: order
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "perdew_wang_lsd_eval", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: npoints, stat, timer_handle
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(dbl)                                :: rho_cutoff
    REAL(dbl), DIMENSION(:, :, :), POINTER   :: a, b, dummy, e_0, ea, eaa, &
                                                eaaa, eaab, eab, eabb, eb, &
                                                ebb, ebbb
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN, 'I', ' ', timer_handle)
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
         ebbb, npoints, order)

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF
    CALL timestop(0.0_dbl, timer_handle)

  END SUBROUTINE perdew_wang_lsd_eval

  
  SUBROUTINE perdew_wang_lsd_calc(rhoa,rhob,e_0,ea,eb,eaa,eab,ebb,eaaa,eaab,eabb,&
       ebbb, npoints, order)
      !FM low-level computation routine
    REAL(dbl), DIMENSION(*), INTENT(in)      :: rhoa, rhob
    REAL(dbl), DIMENSION(*), INTENT(inout)   :: e_0, ea, eb, eaa, eab, ebb, &
                                                eaaa, eaab, eabb, ebbb
    INTEGER, INTENT(in)                      :: npoints, order

    INTEGER                                  :: abs_order, k
    REAL(DBL)                                :: rho
    REAL(DBL), DIMENSION(0:9)                :: ed

    abs_order=ABS(order)

!$omp parallel do private (k, rho, ed)
    DO k=1, npoints

       rho = rhoa(k) + rhob(k)
       IF ( rho > eps_rho ) THEN

         CALL pw_lsd_ed_loc(rhoa(k), rhob(k), ed, abs_order)

         IF (order>=0) THEN
            e_0(k) = e_0(k) + rho*ed(0)
         END IF
         IF (order>=1.OR.order==-1) THEN
            ea(k)   = ea(k) + ed(0) + rho*ed(1)
            eb(k) = eb(k) + ed(0) + rho*ed(2)
         END IF
         IF (order>=2.OR.order==-2) THEN
            eaa(k) = eaa(k) + 2._dbl*ed(1) + rho*ed(3)
            eab(k) = eab(k) + ed(1) + ed(2) + rho*ed(4)
            ebb(k) = ebb(k) + 2._dbl*ed(2) + rho*ed(5)
         END IF
         IF (order>=3.OR.order==-3) THEN
            eaaa(k) = eaaa(k) + 3._dbl*ed(3) + rho*ed(6)
            eaab(k) = eaab(k) + 2._dbl*ed(4) + ed(3) + rho*ed(7)
            eabb(k) = eabb(k) + 2._dbl*ed(4) + ed(5) + rho*ed(8)
            ebbb(k) = ebbb(k) + 3._dbl*ed(5) + rho*ed(9)
         END IF

       END IF
       
    END DO
    
  END SUBROUTINE perdew_wang_lsd_calc

! *****************************************************************************

  SUBROUTINE calc_g(r, z, g, order)

!   ---------------------------------------------------------------------------
!   Calculates g and its derivatives wrt r up to 3rd order, where:
!
!   g = .... for r < 1
!   g = .... for r > 100 and everywhere else
!   g = 2A(1+a1*r)ln(1+1/(2A(b1*r^1/2 + b2*r + b3*r^(3/2) + b4*r^2))).
!   ---------------------------------------------------------------------------


    REAL(DBL), INTENT(IN)                    :: r
    INTEGER, INTENT(IN)                      :: z
    REAL(DBL), DIMENSION(0:), INTENT(OUT)    :: g
    INTEGER, INTENT(IN)                      :: order

    REAL(DBL)                                :: a1_, A_, b1_, b2_, b3_, b4_, &
                                                rr, rsr, sr, t11, t12, t14, &
                                                t15, t16, t20, t22, t3, t40, &
                                                t44, t45, t47, t48, t55, t56

    A_ = A(z); a1_ = a1(z)
    b1_ = b1(z); b2_ = b2(z); b3_ = b3(z); b4_ = b4(z)
    
    sr = SQRT(r)
    rsr = r*sr
    rr = r*r

    IF (r < 1._dbl) THEN
       
       ! order 0 must always be calculated
       g(0) = c0(z)*LOG(r) - c1(z) + c2(z)*r*LOG(r) - c3(z)*r
       IF (order >= 1) g(1) = c0(z)/r + c2(z)*LOG(r) + c2(z) - c3(z)
       IF (order >= 2) g(2) = -c0(z)/rr + c2(z)/r
       IF (order >= 3) g(3) = 2._dbl*c0(z)/(rr*r) - c2(z)/rr

    ELSE IF (r <= 100._dbl) THEN

       t3 = 1._dbl+a1_*r
       t11 = b1_*sr + b2_*r + b3_*rsr + b4_*rr
       t12 = t11**2
       t15 = 1._dbl + 0.5_dbl/A_/t11
       t16 = LOG(t15)
       t20 = 0.5_dbl*b1_/sr + b2_ + 1.5_dbl*b3_*sr + 2._dbl*b4_*r

       ! order 0 must always be calculated
       g(0) = -2._dbl*A_*t3*t16

       IF (order >= 1) THEN

          g(1) = -2._dbl*A_*a1_*t16 + t3*t20/(t12*t15)

       END IF

       IF (order >= 2) THEN

          t40 = -0.25_dbl*b1_/rsr + 0.75_dbl*b3_/sr + 2._dbl*b4_

          g(2) = 2._dbl*a1_*t20/(t12*t15) &
               - 2._dbl*(t20**2)*t3/(t12*t11*t15) &
               + t3*t40/(t12*t15) &
               + 0.5_dbl*t3*(t20**2)/(A_*(t12**2)*(t15**2))

       END IF

       IF (order >= 3) THEN

          t14 = 1._dbl/t12/t11
          t22 = t20**2
          t56 = t22*t20
          t47 = t15**2
          t48 = 1._dbl/t47

          t44 = t12**2
          t45 = 1._dbl/t44
          t55 = t3*t45

          g(3) = &
               - 6._dbl*a1_*t14*t22/t15 &
               + 3._dbl*a1_*t40/(t15*t12) &
               + 1.5_dbl*a1_*t45*t22*t48/A_ &
               + 6._dbl*t55*t56/t15 &
               - 6._dbl*t3*t14*t20*t40/t15 &
               - 3._dbl*t3*t56*t48/(A_*t44*t11) &
               + 0.375_dbl*t3*(b1_/(rr*sr)-b3_/rsr)/(t12*t15) &
               + 1.5_dbl*t55*t40*t48*t20/A_ &
               + 0.5_dbl*t3*t56/((A_**2)*t44*t12*t47*t15)
          
       END IF

    ELSE

       ! order 0 must always be calculated
       g(0) = -d0(z)/r + d1(z)/rsr
       IF (order >= 1) g(1) = d0(z)/rr - 1.5_dbl*d1(z)/(rsr*r)
       IF (order >= 2) g(2) = -2._dbl*d0(z)/(rr*r) + 3.75_dbl*d1(z)/(rsr*rr)
       IF (order >= 3) g(3) = 6._dbl*d0(z)/(rr*rr) - 13.125_dbl*d1(z)/(rsr*rr*r)

    END IF

  END SUBROUTINE calc_g

! *****************************************************************************

SUBROUTINE pw_lda_ed_loc(rho, ed, order)


    REAL(DBL), INTENT(IN)                    :: rho
    REAL(DBL), DIMENSION(0:), INTENT(OUT)    :: ed
    INTEGER, INTENT(IN)                      :: order

    INTEGER                                  :: m, order_
    LOGICAL, DIMENSION(0:3)                  :: calc = .FALSE.
    REAL(DBL), DIMENSION(0:3)                :: e0, r

    order_ = order
    
    IF (order_ >= 0) THEN
       calc(0:order_) = .TRUE.
    ELSE
       order_ = -1 * order_
       calc(order_) = .TRUE.
    END IF

    CALL calc_rs(rho, r(0))
    CALL calc_g(r(0), 0, e0, order_)

    IF (order_ >= 1) r(1) = (-1._dbl/3._dbl)*r(0)/rho
    IF (order_ >= 2) r(2) = (-4._dbl/3._dbl)*r(1)/rho
    IF (order_ >= 3) r(3) = (-7._dbl/3._dbl)*r(2)/rho

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
       ed(m) = e0(3)*r(1)**3 + e0(2)*3._dbl*r(1)*r(2) + e0(1)*r(3)
    END IF
    
  END SUBROUTINE  pw_lda_ed_loc

! *****************************************************************************

  SUBROUTINE pw_lsd_ed_loc(a, b, ed, order)


    REAL(DBL), INTENT(IN)                    :: a, b
    REAL(DBL), DIMENSION(0:), INTENT(OUT)    :: ed
    INTEGER, INTENT(IN)                      :: order

    INTEGER                                  :: m, order_ = 0
    LOGICAL, DIMENSION(0:3)                  :: calc = .FALSE.
    REAL(DBL)                                :: rho, tr, trr, trrr, trrz, &
                                                trz, trzz, tz, tzz, tzzz
    REAL(DBL), DIMENSION(0:3)                :: ac, e0, e1, f, r
    REAL(DBL), DIMENSION(0:3, 0:3)           :: z

    order_ = order
    
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
       r(1) = (-1._dbl/3._dbl)*r(0)/rho
       tr = e0(1) &
            + fpp*ac(1)*f(0) &
            - fpp*ac(1)*f(0)*z(0,0)**4 &
            + (e1(1)-e0(1))*f(0)*z(0,0)**4
       tz = fpp*ac(0)*f(1) &
            - fpp*ac(0)*f(1)*z(0,0)**4 &
            - fpp*ac(0)*f(0)*4._dbl*z(0,0)**3 &
            + (e1(0)-e0(0))*f(1)*z(0,0)**4 &
            + (e1(0)-e0(0))*f(0)*4._dbl*z(0,0)**3
    END IF
    
!! calculate second partial derivatives
    IF (order_ >= 2) THEN
       r(2) = (-4._dbl/3._dbl)*r(1)/rho
       trr = e0(2) &
            + fpp*ac(2)*f(0) &
            - fpp*ac(2)*f(0)*z(0,0)**4 &
            + (e1(2)-e0(2))*f(0)*z(0,0)**4
       trz = fpp*ac(1)*f(1) &
            - fpp*ac(1)*f(1)*z(0,0)**4 &
            - fpp*ac(1)*f(0)*4._dbl*z(0,0)**3 &
            + (e1(1)-e0(1))*f(1)*z(0,0)**4 &
            + (e1(1)-e0(1))*f(0)*4._dbl*z(0,0)**3
       tzz = fpp*ac(0)*f(2) &
            - fpp*ac(0)*f(2)*z(0,0)**4 &
            - fpp*ac(0)*f(1)*8._dbl*z(0,0)**3 &
            - fpp*ac(0)*f(0)*12._dbl*z(0,0)**2 &
            + (e1(0)-e0(0))*f(2)*z(0,0)**4 &
            + (e1(0)-e0(0))*f(1)*8._dbl*z(0,0)**3 &
            + (e1(0)-e0(0))*f(0)*12._dbl*z(0,0)**2
    END IF
    
!! calculate third derivatives
    IF (order_ >= 3) THEN
       
       r(3) = (-7._dbl/3._dbl)*r(2)/rho
       
       trrr = e0(3) &
            + fpp*ac(3)*f(0) &
            - fpp*ac(3)*f(0)*z(0,0)**4 &
            + (e1(3)-e0(3))*f(0)*z(0,0)**4
       
       trrz = fpp*ac(2)*f(1) &
            - fpp*ac(2)*f(1)*z(0,0)**4 &
            - fpp*ac(2)*f(0)*4._dbl*z(0,0)**3 &
            + (e1(2)-e0(2))*f(1)*z(0,0)**4 &
            + (e1(2)-e0(2))*f(0)*4._dbl*z(0,0)**3
       
       trzz = fpp*ac(1)*f(2) &
            - fpp*ac(1)*f(2)*z(0,0)**4 &
            - fpp*ac(1)*f(1)*8._dbl*z(0,0)**3 &
            - fpp*ac(1)*f(0)*12._dbl*z(0,0)**2 &
            + (e1(1)-e0(1))*f(2)*z(0,0)**4 &
            + (e1(1)-e0(1))*f(1)*8._dbl*z(0,0)**3 &
            + (e1(1)-e0(1))*f(0)*12._dbl*z(0,0)**2
       
       tzzz = fpp*ac(0)*f(3) &
            - fpp*ac(0)*f(3)*z(0,0)**4 &
            - fpp*ac(0)*f(2)*12._dbl*z(0,0)**3 &
            - fpp*ac(0)*f(1)*36._dbl*z(0,0)**2 &
            - fpp*ac(0)*f(0)*24._dbl*z(0,0) &
            + (e1(0)-e0(0))*f(3)*z(0,0)**4 &
            + (e1(0)-e0(0))*f(2)*12._dbl*z(0,0)**3 &
            + (e1(0)-e0(0))*f(1)*36._dbl*z(0,0)**2 &
            + (e1(0)-e0(0))*f(0)*24._dbl*z(0,0)
    END IF

    m=0
    IF (calc(0)) THEN
       ed(m) = e0(0) &
            + fpp*ac(0)*f(0)*(1._dbl - z(0,0)**4) &
            + (e1(0)-e0(0))*f(0)*z(0,0)**4
       m = m + 1
    END IF
    IF (calc(1)) THEN
       ed(m)   = tr*r(1) + tz*z(1,0)
       ed(m+1) = tr*r(1) + tz*z(0,1)
       m = m + 2
    END IF
    IF (calc(2)) THEN
       ed(m)   = trr*r(1)**2 + 2._dbl*trz*r(1)*z(1,0) &
            + tr*r(2) + tzz*z(1,0)**2 + tz*z(2,0)
       ed(m+1) = trr*r(1)**2 + trz*r(1)*(z(0,1)+z(1,0)) &
            + tr*r(2) + tzz*z(1,0)*z(0,1) + tz*z(1,1)
       ed(m+2) = trr*r(1)**2 + 2._dbl*trz*r(1)*z(0,1) &
            + tr*r(2) + tzz*z(0,1)**2 + tz*z(0,2)
       m = m + 3
    END IF
    IF (calc(3)) THEN
       ed(m)   = &
            trrr*r(1)**3 + 3._dbl*trrz*r(1)**2*z(1,0) &
            + 3._dbl*trr*r(1)*r(2) + 3._dbl*trz*r(2)*z(1,0) + tr*r(3) &
            + 3._dbl*trzz*r(1)*z(1,0)**2 + tzzz*z(1,0)**3 &
            + 3._dbl*trz*r(1)*z(2,0) &
            + 3._dbl*tzz*z(1,0)*z(2,0) + tz*z(3,0)
       ed(m+1) = &
            trrr*r(1)**3 + trrz*r(1)**2*(2._dbl*z(1,0)+z(0,1)) &
            + 2._dbl*trzz*r(1)*z(1,0)*z(0,1) &
            + 2._dbl*trz*(r(2)*z(1,0)+r(1)*z(1,1)) &
            + 3._dbl*trr*r(2)*r(1) + trz*r(2)*z(0,1) + tr*r(3) &
            + trzz*r(1)*z(1,0)**2 + tzzz*z(1,0)**2*z(0,1) &
            + 2._dbl*tzz*z(1,0)*z(1,1) &
            + trz*r(1)*z(2,0) + tzz*z(2,0)*z(0,1) + tz*z(2,1)
       ed(m+2) = &
            trrr*r(1)**3 + trrz*r(1)**2*(2._dbl*z(0,1)+z(1,0)) &
            + 2._dbl*trzz*r(1)*z(0,1)*z(1,0) &
            + 2._dbl*trz*(r(2)*z(0,1)+r(1)*z(1,1)) &
            + 3._dbl*trr*r(2)*r(1) + trz*r(2)*z(1,0) + tr*r(3) &
            + trzz*r(1)*z(0,1)**2 + tzzz*z(0,1)**2*z(1,0) &
            + 2._dbl*tzz*z(0,1)*z(1,1) &
            + trz*r(1)*z(0,2) + tzz*z(0,2)*z(1,0) + tz*z(1,2)
       ed(m+3) = &
            trrr*r(1)**3 + 3._dbl*trrz*r(1)**2*z(0,1) &
            + 3._dbl*trr*r(1)*r(2) + 3._dbl*trz*r(2)*z(0,1) + tr*r(3) &
            + 3._dbl*trzz*r(1)*z(0,1)**2 + tzzz*z(0,1)**3 &
            + 3._dbl*trz*r(1)*z(0,2) &
            + 3._dbl*tzz*z(0,1)*z(0,2) + tz*z(0,3)
    END IF
       
  END SUBROUTINE pw_lsd_ed_loc

! *****************************************************************************

END MODULE xc_perdew_wang_new
