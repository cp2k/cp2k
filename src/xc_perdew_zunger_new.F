!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****m* cp2k/xc_perdew_zunger *
!!
!!   NAME
!!     xc_perdew_zunger
!!
!!   FUNCTION
!!     Calculate the Perdew-Zunger correlation potential and 
!!     energy density and ist derivatives with respect to
!!     the spin-up and spin-down densities up to 3rd order.
!!
!!   SYNOPSIS
!!     perdew_zunger ( method, cutoff )
!!     perdew_zunger ( reference, shortform )
!!     perdew_zunger ( rho, ec, order )
!!     perdew_zunger ( rhoa, rhob, ec, order )
!!     perdew_zunger_energy_density ( rho, ec, order )
!!     perdew_zunger_energy_density ( rhoa, rhob, ec, order )
!!
!!   MODIFICATION HISTORY
!!     18-MAR-2002, TCH, working version
!!     fawzi (04.2004)  : adapted to the new xc interface
!!
!!   SEE ALSO
!!     functionals_utilities
!!
!!****

MODULE xc_perdew_zunger_new
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
  USE mathconstants,                   ONLY: zero
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
  REAL(DBL), DIMENSION(0:1) :: A = (0,0), B = (0,0), C = (0,0), D = (0,0), &
                               b1 = (0,0), b2 = (0,0), ga = (0,0)

  REAL(DBL), PARAMETER :: epsilon = 5.E-13
  REAL(DBL) :: eps_rho

  PUBLIC :: pz_info, pz_lda_eval, pz_lsd_eval
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="xc_perdew_zunger_new"

CONTAINS

! *****************************************************************************

!!****f* cp2k/xc_perdew_zunger/pz_info *
!!
!!   NAME
!!     perdew_zunger
!!
!!   FUNCTION
!!     Return some info on the functionals.
!!
!!   ARGUMENTS
!!     reference - CHARACTER(*), INTENT(OUT), OPTIONAL - full reference
!!     shortform - CHARACTER(*), INTENT(OUT), OPTIONAL - short reference
!!
!!   MODIFICATION HISTORY
!!     18-MAR-2002, TCH, working version
!!
!!   SEE ALSO
!!     
!!
!!****
  SUBROUTINE pz_info ( method, lsd, reference, shortform, needs, max_deriv, error)
    CHARACTER(len=*), INTENT(in)             :: method
    LOGICAL, INTENT(in)                      :: lsd
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pz_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    SELECT CASE (method)
    CASE DEFAULT
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    CASE ('ORIG')
    CASE ('DMC')
    CASE ('VMC') 
    END SELECT

    IF ( PRESENT ( reference ) ) THEN
       reference = "J. P. Perdew and Alex Zunger," &
                   //" Phys. Rev. B 23, 5048 (1981)"&
                   //"["//TRIM(method)//"]"
       IF (.not.lsd) THEN
          IF (LEN_TRIM(reference)+6<LEN(reference)) THEN
             reference(LEN_TRIM(reference):LEN_TRIM(reference)+6)=' {LDA}'
          END IF
       END IF
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "J. P. Perdew et al., PRB 23, 5048 (1981)"&
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

  END SUBROUTINE pz_info

! *****************************************************************************

!!****f* cp2k/xc_perdew_zunger/pz_lda_eval *
!!
!!   NAME
!!     pz_lda_eval
!!
!!   SYNOPSIS
!!     Subroutine pz_lda_eval(method, rho_set, deriv_set, order, error)
!!       Character(Len=*), Intent (IN):: method
!!       Type(xc_rho_set_type), Pointer:: rho_set
!!       Type(xc_derivative_set_type), Pointer:: deriv_set
!!       Integer, Intent (IN):: order
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pz_lda_eval
!!
!!   FUNCTION
!!     Calculate the correlation energy and its derivatives
!!     wrt to rho (the electron density) up to 3rd order. This
!!     is the LDA version of the Perdew-Zunger correlation energy
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
  SUBROUTINE pz_lda_eval ( method, rho_set, deriv_set, order, error )
    
    CHARACTER(len=*), INTENT(in)             :: method
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: order
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "pz_lda_eval", &
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

       CALL pz_init(method, rho_cutoff)
       
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

       CALL pz_lda_calc(rho,e_0,e_rho, e_rho_rho, e_rho_rho_rho, npoints, order)

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF
    END IF

    CALL timestop(zero, timer_handle)
    
  END SUBROUTINE pz_lda_eval

! *****************************************************************************

SUBROUTINE pz_lda_calc(rho,e_0,e_rho, e_rho_rho, e_rho_rho_rho, npoints, order)
  !FM low level calc routine
    REAL(dbl), DIMENSION(*), INTENT(in)      :: rho
    REAL(dbl), DIMENSION(*), INTENT(inout)   :: e_0, e_rho, e_rho_rho, &
                                                e_rho_rho_rho
    INTEGER, INTENT(in)                      :: npoints, order

    INTEGER                                  :: k
    REAL(DBL), DIMENSION(0:3)                :: ed

!$omp parallel do private ( k, ed )
    DO k=1, npoints
       
       IF ( rho(k) > eps_rho ) THEN

         CALL pz_lda_ed_loc(rho(k), ed, ABS(order))
       
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

  END SUBROUTINE pz_lda_calc

! *****************************************************************************

!!****f* cp2k/xc_perdew_zunger/pz_lsd *
!!
!!   NAME
!!     perdew_zunger
!!
!!   FUNCTION
!!     Calculate the correlation energy and its derivatives
!!     wrt to rho (the electron density) up to 3rd order. This
!!     is the LSD version of the Perdew-Zunger correlation energy
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
  SUBROUTINE pz_lsd_eval ( method, rho_set, deriv_set, order, error )
    CHARACTER(len=*), INTENT(in)             :: method
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN), OPTIONAL            :: order
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "pz_lsd_eval", &
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

       CALL pz_init(method, rho_cutoff)
       
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
       
       CALL pz_lsd_calc(a,b,e_0,ea,eb,eaa,eab,ebb,eaaa,eaab,eabb,&
         ebbb, npoints, order)

       IF (cp_debug) THEN
          DEALLOCATE(dummy,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          NULLIFY(dummy)
       END IF

    END IF
    CALL timestop(zero, timer_handle)

  END SUBROUTINE pz_lsd_eval

  
    SUBROUTINE pz_lsd_calc(a,b,e_0,ea,eb,eaa,eab,ebb,eaaa,eaab,eabb,&
         ebbb, npoints, order)
      !FM low-level computation routine
    REAL(dbl), DIMENSION(*), INTENT(in)      :: a, b
    REAL(dbl), DIMENSION(*), INTENT(inout)   :: e_0, ea, eb, eaa, eab, ebb, &
                                                eaaa, eaab, eabb, ebbb
    INTEGER, INTENT(in)                      :: npoints, order

    INTEGER                                  :: k, order_
    REAL(DBL)                                :: rho
    REAL(DBL), DIMENSION(0:9)                :: ed

    order_=ABS(order)

!$omp parallel do private ( k, rho, ed ) shared(order,npoints,eps_rho)
    DO k=1, npoints

       rho = a(k) + b(k)

       IF ( rho > eps_rho ) THEN

          CALL pz_lsd_ed_loc(a(k), b(k), ed, order_)

          IF (order>=0) THEN
             e_0(k) = e_0(k) + rho*ed(0)
          END IF

          IF (order>=1.OR.order==-1) THEN
             ea(k) = ea(k) + ed(0) + rho*ed(1)
             ea(k) = ea(k) + ed(0) + rho*ed(2)
          END IF
       
          IF (order>=2.OR.order==-2) THEN
             eaa(k)   = eaa(k) + 2._dbl*ed(1) + rho*ed(3)
             eab(k) = eab(k) + ed(1) + ed(2) + rho*ed(4)
             ebb(k) = ebb(k) + 2._dbl*ed(2) + rho*ed(5)
          END IF
       
          IF (order>=3.OR.order==-3) THEN
             eaaa(k)   = eaaa(k) + 3._dbl*ed(3) + rho*ed(6)
             eaab(k) = eaab(k) + 2._dbl*ed(4) + ed(3) + rho*ed(7)
             eabb(k) = eabb(k) + 2._dbl*ed(4) + ed(5) + rho*ed(8)
             ebbb(k) = ebbb(k) + 3._dbl*ed(5) + rho*ed(9)
          END IF

       END IF

    END DO
!$omp end parallel do
       
  END SUBROUTINE pz_lsd_calc

! *****************************************************************************

!!****f* cp2k/xc_perdew_zunger/pz_init *
!!
!!   NAME
!!     perdew_zunger
!!
!!   FUNCTION
!!     Initializes the functionals
!!
!!   ARGUMENTS
!!     method - CHARACTER(3) - name of the method used for parameters
!!     cutoff - REAL(DBL) - the cutoff density
!!
!!   MODIFICATION HISTORY
!!     18-MAR-2002, TCH, working version
!!
!!   SEE ALSO
!!     functionals_utilities
!!
!!****
  SUBROUTINE pz_init ( method, cutoff )
    

    CHARACTER(*), INTENT(IN)                 :: method
    REAL(DBL), INTENT(IN)                    :: cutoff

    CALL set_util(cutoff)

    eps_rho = cutoff
    
    initialized = .FALSE.

    SELECT CASE (method)

    CASE DEFAULT
       CALL stop_program("xc_perdew_zunger/pz_init", "unknown method")
       
    CASE ('ORIG')
       ga(0)  = -0.1423_dbl;  ga(1)  = -0.0843_dbl
       b1(0) = 1.0529_dbl;    b1(1) = 1.3981_dbl
       b2(0) = 0.3334_dbl;    b2(1) = 0.2611_dbl
       A(0)  = 0.0311_dbl;    A(1)  = 0.01555_dbl
       B(0)  = -0.048_dbl;    B(1)  = -0.0269_dbl
       C(0)  = -0.0020_dbl;   C(1)  = -0.0007_dbl
       D(0)  = -0.0116_dbl;   D(1)  = -0.0048_dbl

    CASE ('DMC')
       ga(0)  = -0.103756_dbl; ga(1)  = -0.065951_dbl
       b1(0) = 0.56371_dbl;   b1(1) = 1.11846_dbl
       b2(0) = 0.27358_dbl;   b2(1) = 0.18797_dbl
       A(0)  = 0.031091_dbl;  A(1)  = 0.015545_dbl
       B(0)  = -0.046644_dbl; B(1)  = -0.025599_dbl
       C(0)  = -0.00419_dbl;  C(1)  = -0.00329_dbl
       D(0)  = -0.00983_dbl;  D(1)  = -0.00300_dbl
       
    CASE ('VMC') 
       ga(0)  = -0.093662_dbl; ga(1)  = -0.055331_dbl
       b1(0) = 0.49453_dbl;   b1(1) = 0.93766_dbl
       b2(0) = 0.25534_dbl;   b2(1) = 0.14829_dbl
       A(0)  = 0.031091_dbl;  A(1)  = 0.015545_dbl
       B(0)  = -0.046644_dbl; B(1)  = -0.025599_dbl
       C(0)  = -0.00884_dbl;  C(1)  = -0.00677_dbl
       D(0)  = -0.00688_dbl;  D(1)  = -0.00093_dbl
       
    END SELECT
          
    initialized = .TRUE.
    
  END SUBROUTINE pz_init

! *****************************************************************************

  SUBROUTINE calc_g(r, z, g, order)


    REAL(DBL), INTENT(IN)                    :: r
    INTEGER, INTENT(IN)                      :: z
    REAL(DBL), DIMENSION(0:), INTENT(OUT)    :: g
    INTEGER, INTENT(IN)                      :: order

    REAL(DBL)                                :: rr, rsr, sr

    IF (r >= 1.0_dbl) THEN

       sr = SQRT(r)
       ! order 0 must always be calculated
       g(0) = ga(z) / (1.0_dbl + b1(z)*sr + b2(z)*r)
       IF (order >= 1) THEN
          g(1) = -ga(z) * (b1(z)/(2.0_dbl*sr) + b2(z)) / &
               (1.0_dbl + b1(z)*sr + b2(z)*r)**2
       END IF
       IF (order >= 2) THEN
          rsr = r*sr
          g(2) = &
               2._dbl * ga(z) * (b1(z)/(2.0_dbl*sr) + b2(z))**2 &
               / (1._dbl + b1(z)*sr + b2(z)*r)**3 &
               + ga(z) * b1(z) &
               / (4.0_dbl * (1._dbl + b1(z)*sr + b2(z)*r)**2 * rsr)
       END IF
       IF (order >= 3) THEN
          g(3) = &
               -6.0_dbl * ga(z) * (b1(z)/(2.0_dbl*sr) + b2(z))**3 / &
               (1.0_dbl + b1(z)*sr + b2(z)*r)**4 &
               -(3._dbl/2._dbl) * ga(z) * (b1(z)/(2.0_dbl*sr) + b2(z)) * b1(z) / &
               ((1._dbl + b1(z)*sr + b2(z)*r)**3 * rsr) &
               -(3._dbl/8._dbl) * ga(z) * b1(z) / &
               ((1._dbl + b1(z)*sr + b2(z)*r)**2 * r*rsr)
       END IF

    ELSE

       ! order 0 must always be calculated
       g(0) = A(z)*LOG(r) + B(z) + C(z)*r*LOG(r) + D(z)*r
       IF (order >= 1) THEN
          g(1) = A(z)/r + C(z)*LOG(r) + C(z) + D(z)
       END IF
       IF (order >= 2) THEN
          rr = r*r
          g(2) = -A(z)/rr + C(z)/r
       END IF
       IF (order >= 3) THEN
          g(3) = 2.0_dbl*A(z)/(rr*r) - C(z)/rr
       END IF
       
    END IF

  END SUBROUTINE calc_g

! *****************************************************************************

  SUBROUTINE pz_lda_ed_loc ( rho, ed, order )
    
    
    REAL(DBL), INTENT(IN)                    :: rho
    REAL(DBL), DIMENSION(0:), INTENT(OUT)    :: ed
    INTEGER, INTENT(IN)                      :: order

    INTEGER                                  :: m, order_
    LOGICAL, DIMENSION(0:3)                  :: calc
    REAL(DBL), DIMENSION(0:3)                :: e0, r

    calc = .FALSE.

    order_ = order
    IF (order_ >= 0) THEN
       calc(0:order_) = .TRUE.
    ELSE
       order_ = -1 * order_
       calc(order_) = .TRUE.
    END IF

    CALL calc_rs(rho, r(0))
    CALL calc_g(r(0), 0, e0, order_)

    IF (order_ >= 1) r(1) = (-1._dbl/3._dbl) * r(0) / rho
    IF (order_ >= 2) r(2) = (-4._dbl/3._dbl) * r(1) / rho
    IF (order_ >= 3) r(3) = (-7._dbl/3._dbl) * r(2) / rho

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

  END SUBROUTINE pz_lda_ed_loc

! *****************************************************************************

  SUBROUTINE pz_lsd_ed_loc ( a, b, ed, order )
    
    
    REAL(DBL), INTENT(IN)                    :: a, b
    REAL(DBL), DIMENSION(0:), INTENT(OUT)    :: ed
    INTEGER, INTENT(IN), OPTIONAL            :: order

    INTEGER                                  :: m, order_
    LOGICAL, DIMENSION(0:3)                  :: calc
    REAL(DBL)                                :: rho, tr, trr, trrr, trrz, &
                                                trz, trzz, tz, tzz, tzzz
    REAL(DBL), DIMENSION(0:3)                :: e0, e1, f, r
    REAL(DBL), DIMENSION(0:3, 0:3)           :: z

    calc = .FALSE.

    order_ = 0
    IF (PRESENT(order)) order_ = order

    IF (order_ >= 0) THEN
       calc(0:order_) = .TRUE.
    ELSE
       order_ = -1 * order_
       calc(order_) = .TRUE.
    END IF

    rho = a + b

    CALL calc_fx(a, b, f(0:order_), order_)
    CALL calc_rs(rho, r(0))
    CALL calc_g(r(0), 0, e0(0:order_), order_)
    CALL calc_g(r(0), 1, e1(0:order_), order_)
    CALL calc_z(a, b, z(:,:), order_)

!! calculate first partial derivatives
    IF (order_ >= 1) THEN
       r(1) = (-1._dbl/3._dbl)*r(0)/rho
       tr = e0(1) + (e1(1)-e0(1))*f(0)
       tz = (e1(0)-e0(0))*f(1)
    END IF
    
!! calculate second partial derivatives
    IF (order_ >= 2) THEN
       r(2) = (-4._dbl/3._dbl)*r(1)/rho
       trr = e0(2) + (e1(2)-e0(2))*f(0)
       trz = (e1(1)-e0(1))*f(1)
       tzz = (e1(0)-e0(0))*f(2)
    END IF
    
!! calculate third derivatives
    IF (order_ >= 3) THEN
       r(3) = (-7._dbl/3._dbl)*r(2)/rho
       trrr = e0(3) + (e1(3)-e0(3))*f(0)
       trrz = (e1(2)-e0(2))*f(1)
       trzz = (e1(1)-e0(1))*f(2)
       tzzz = (e1(0)-e0(0))*f(3)
    END IF
    
    m=0
    IF (calc(0)) THEN
       ed(m) = e0(0) + (e1(0)-e0(0))*f(0)
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
       ed(m)   = trrr*r(1)**3 + 3._dbl*trrz*r(1)**2*z(1,0) &
            + 3._dbl*trr*r(1)*r(2) + 3._dbl*trz*r(2)*z(1,0) &
            + tr*r(3) + 3._dbl*trzz*r(1)*z(1,0)**2 &
            + tzzz*z(1,0)**3 + 3._dbl*trz*r(1)*z(2,0) &
            + 3._dbl*tzz*z(1,0)*z(2,0) + tz*z(3,0)
       ed(m+1) = trrr*r(1)**3 + trrz*r(1)**2*(2._dbl*z(1,0)+z(0,1)) &
            + 2._dbl*trzz*r(1)*z(1,0)*z(0,1) &
            + 2._dbl*trz*(r(2)*z(1,0)+r(1)*z(1,1)) &
            + 3._dbl*trr*r(2)*r(1) + trz*r(2)*z(0,1) + tr*r(3) &
            + trzz*r(1)*z(1,0)**2 + tzzz*z(1,0)**2*z(0,1) &
            + 2._dbl*tzz*z(1,0)*z(1,1) &
            + trz*r(1)*z(2,0) + tzz*z(2,0)*z(0,1) + tz*z(2,1)
       ed(m+2) = trrr*r(1)**3 + trrz*r(1)**2*(2._dbl*z(0,1)+z(1,0)) &
            + 2._dbl*trzz*r(1)*z(0,1)*z(1,0) &
            + 2._dbl*trz*(r(2)*z(0,1)+r(1)*z(1,1)) &
            + 3._dbl*trr*r(2)*r(1) + trz*r(2)*z(1,0) + tr*r(3) &
            + trzz*r(1)*z(0,1)**2 + tzzz*z(0,1)**2*z(1,0) &
            + 2._dbl*tzz*z(0,1)*z(1,1) &
            + trz*r(1)*z(0,2) + tzz*z(0,2)*z(1,0) + tz*z(1,2)
       ed(m+3) = trrr*r(1)**3 + 3._dbl*trrz*r(1)**2*z(0,1) &
            + 3._dbl*trr*r(1)*r(2) + 3._dbl*trz*r(2)*z(0,1) + tr*r(3) &
            + 3._dbl*trzz*r(1)*z(0,1)**2 + tzzz*z(0,1)**3 &
            + 3._dbl*trz*r(1)*z(0,2) &
            + 3._dbl*tzz*z(0,1)*z(0,2) + tz*z(0,3)
    END IF
    
  END SUBROUTINE pz_lsd_ed_loc

! *****************************************************************************

END MODULE xc_perdew_zunger_new
