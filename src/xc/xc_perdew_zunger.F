!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the Perdew-Zunger correlation potential and
!>      energy density and ist derivatives with respect to
!>      the spin-up and spin-down densities up to 3rd order.
!> \par History
!>      18-MAR-2002, TCH, working version
!>      fawzi (04.2004)  : adapted to the new xc interface
!> \see functionals_utilities
! *****************************************************************************
MODULE xc_perdew_zunger
  USE bibliography,                    ONLY: Ortiz1994,&
                                             Perdew1981,&
                                             cite_reference
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
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
  USE xc_input_constants,              ONLY: pz_dmc,&
                                             pz_orig,&
                                             pz_vmc
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL :: initialized = .FALSE.
  REAL(KIND=dp), DIMENSION(0:1) :: A = 0.0_dp, B = 0.0_dp, C = 0.0_dp, D = 0.0_dp, &
                               b1 = 0.0_dp, b2 = 0.0_dp, ga = 0.0_dp

  REAL(KIND=dp), PARAMETER :: epsilon = 5.E-13
  REAL(KIND=dp) :: eps_rho

  PUBLIC :: pz_info, pz_lda_eval, pz_lsd_eval
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_perdew_zunger'

CONTAINS

! *****************************************************************************
!> \brief Return some info on the functionals.
!> \param method ...
!> \param lsd ...
!> \param reference CHARACTER(*), INTENT(OUT), OPTIONAL - full reference
!> \param shortform CHARACTER(*), INTENT(OUT), OPTIONAL - short reference
!> \param needs ...
!> \param max_deriv ...
!> \param error ...
!> \par History
!>      18-MAR-2002, TCH, working version
! *****************************************************************************
  SUBROUTINE pz_info ( method, lsd, reference, shortform, needs, max_deriv, error)
    INTEGER, INTENT(in)                      :: method
    LOGICAL, INTENT(in)                      :: lsd
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pz_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=4)                         :: p_string
    LOGICAL                                  :: failure

    failure=.FALSE.
    SELECT CASE (method)
    CASE DEFAULT
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    CASE (pz_orig)
       p_string='ORIG'
    CASE (pz_dmc)
       p_string='DMC'
    CASE (pz_vmc)
       p_string='VMC'
    END SELECT

    IF ( PRESENT ( reference ) ) THEN
       reference = "J. P. Perdew and Alex Zunger," &
                   //" Phys. Rev. B 23, 5048 (1981)"&
                   //"["//TRIM(p_string)//"]"
       IF (.not.lsd) THEN
          IF (LEN_TRIM(reference)+6<LEN(reference)) THEN
             reference(LEN_TRIM(reference):LEN_TRIM(reference)+6)=' {LDA}'
          END IF
       END IF
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "J. P. Perdew et al., PRB 23, 5048 (1981)"&
                   //"["//TRIM(p_string)//"]"
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
!> \brief Calculate the correlation energy and its derivatives
!>      wrt to rho (the electron density) up to 3rd order. This
!>      is the LDA version of the Perdew-Zunger correlation energy
!>      If no order argument is given, then the routine calculates
!>      just the energy.
!> \param method ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order INTEGER, OPTIONAL - order of derivatives to calculate
!>        order must lie between -3 and 3. If it is negative then only
!>        that order will be calculated, otherwise all derivatives up to
!>        that order will be calculated.
!> \param pz_params input paramter (scaling)
!> \param error ...
!> \par History
!>     01.2007 added scaling [Manuel Guidon]
! *****************************************************************************
  SUBROUTINE pz_lda_eval ( method, rho_set, deriv_set, order, pz_params, error )

    INTEGER, INTENT(in)                      :: method
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: order
    TYPE(section_vals_type), POINTER         :: pz_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pz_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: npoints, stat, timer_handle
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rho_cutoff, sc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_rho, e_rho_rho, &
                                                e_rho_rho_rho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,timer_handle)
    failure=.FALSE.
    NULLIFY(bo,rho, e_0,e_rho,e_rho_rho,e_rho_rho_rho, dummy)

    CALL section_vals_val_get(pz_params,"scale_c",r_val=sc,error=error)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
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

    CALL pz_lda_calc(rho,e_0,e_rho, e_rho_rho, e_rho_rho_rho, npoints, order, sc)

    IF (cp_debug) THEN
       DEALLOCATE(dummy,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    ELSE
       NULLIFY(dummy)
    END IF

    CALL timestop(timer_handle)

  END SUBROUTINE pz_lda_eval

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param e_0 ...
!> \param e_rho ...
!> \param e_rho_rho ...
!> \param e_rho_rho_rho ...
!> \param npoints ...
!> \param order ...
!> \param sc ...
! *****************************************************************************
SUBROUTINE pz_lda_calc(rho,e_0,e_rho, e_rho_rho, e_rho_rho_rho, npoints, order, sc)
  !FM low level calc routine
    REAL(KIND=dp), DIMENSION(*), INTENT(in)  :: rho
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(inout)                          :: e_0, e_rho, e_rho_rho, &
                                                e_rho_rho_rho
    INTEGER, INTENT(in)                      :: npoints, order
    REAL(KIND=dp)                            :: sc

    INTEGER                                  :: k
    REAL(KIND=dp), DIMENSION(0:3)            :: ed

!$omp parallel do private ( k, ed ) default(none)&
!$omp shared(npoints,rho,eps_rho,order,e_0,e_rho,e_rho_rho,e_rho_rho_rho,sc)
    DO k=1, npoints

       IF ( rho(k) > eps_rho ) THEN

         CALL pz_lda_ed_loc(rho(k), ed, ABS(order), sc)

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

  END SUBROUTINE pz_lda_calc

! *****************************************************************************
!> \brief Calculate the correlation energy and its derivatives
!>      wrt to rho (the electron density) up to 3rd order. This
!>      is the LSD version of the Perdew-Zunger correlation energy
!>      If no order argument is given, then the routine calculates
!>      just the energy.
!> \param method ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order INTEGER, OPTIONAL - order of derivatives to calculate
!>        order must lie between -3 and 3. If it is negative then only
!>        that order will be calculated, otherwise all derivatives up to
!>        that order will be calculated.
!> \param pz_params input paramter (scaling)
!> \param error ...
!> \par History
!>      01.2007 added scaling [Manuel Guidon]
! *****************************************************************************
  SUBROUTINE pz_lsd_eval ( method, rho_set, deriv_set, order,pz_params, error )
    INTEGER, INTENT(in)                      :: method
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN), OPTIONAL            :: order
    TYPE(section_vals_type), POINTER         :: pz_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pz_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: npoints, stat, timer_handle
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rho_cutoff, sc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: a, b, dummy, e_0, ea, eaa, &
                                                eaaa, eaab, eab, eabb, ebb, &
                                                ebbb
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,timer_handle)
    failure=.FALSE.
    NULLIFY(bo,a,b,e_0,ea,eaa,eab,ebb,eaaa,eaab,eabb,ebbb)

    CALL section_vals_val_get(pz_params,"scale_c",r_val=sc,error=error)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
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
    ea => dummy; 
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

    CALL pz_lsd_calc(a,b,e_0,ea,eaa,eab,ebb,eaaa,eaab,eabb,&
      ebbb, npoints, order, sc)

    IF (cp_debug) THEN
       DEALLOCATE(dummy,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    ELSE
       NULLIFY(dummy)
    END IF

    CALL timestop(timer_handle)

  END SUBROUTINE pz_lsd_eval

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param e_0 ...
!> \param ea ...
!> \param eaa ...
!> \param eab ...
!> \param ebb ...
!> \param eaaa ...
!> \param eaab ...
!> \param eabb ...
!> \param ebbb ...
!> \param npoints ...
!> \param order ...
!> \param sc ...
! *****************************************************************************
    SUBROUTINE pz_lsd_calc(a,b,e_0,ea,eaa,eab,ebb,eaaa,eaab,eabb,&
         ebbb, npoints, order, sc)
      !FM low-level computation routine
    REAL(KIND=dp), DIMENSION(*), INTENT(in)  :: a, b
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(inout)                          :: e_0, ea, eaa, eab, ebb, eaaa, &
                                                eaab, eabb, ebbb
    INTEGER, INTENT(in)                      :: npoints, order
    REAL(KIND=dp), INTENT(IN)                :: sc

    INTEGER                                  :: k, order_
    REAL(KIND=dp)                            :: rho
    REAL(KIND=dp), DIMENSION(0:9)            :: ed

    order_=ABS(order)

!$omp parallel do private ( k, rho, ed ) default(none)&
!$omp shared(order_,order,npoints,eps_rho,a,b,sc,e_0,ea,eaa,eab,ebb,eaaa,eaab,eabb,ebbb)
    DO k=1, npoints

       rho = a(k) + b(k)

       IF ( rho > eps_rho ) THEN

          CALL pz_lsd_ed_loc(a(k), b(k), ed, order_, sc)

          IF (order>=0) THEN
             e_0(k) = e_0(k) + rho*ed(0)
          END IF

          IF (order>=1.OR.order==-1) THEN
             ea(k) = ea(k) + ed(0) + rho*ed(1)
             ea(k) = ea(k) + ed(0) + rho*ed(2)
          END IF

          IF (order>=2.OR.order==-2) THEN
             eaa(k)   = eaa(k) + 2.0_dp*ed(1) + rho*ed(3)
             eab(k) = eab(k) + ed(1) + ed(2) + rho*ed(4)
             ebb(k) = ebb(k) + 2.0_dp*ed(2) + rho*ed(5)
          END IF

          IF (order>=3.OR.order==-3) THEN
             eaaa(k)   = eaaa(k) + 3.0_dp*ed(3) + rho*ed(6)
             eaab(k) = eaab(k) + 2.0_dp*ed(4) + ed(3) + rho*ed(7)
             eabb(k) = eabb(k) + 2.0_dp*ed(4) + ed(5) + rho*ed(8)
             ebbb(k) = ebbb(k) + 3.0_dp*ed(5) + rho*ed(9)
          END IF

       END IF

    END DO
!$omp end parallel do

  END SUBROUTINE pz_lsd_calc

! *****************************************************************************
!> \brief Initializes the functionals
!> \param method CHARACTER(3) - name of the method used for parameters
!> \param cutoff REAL(KIND=dp) - the cutoff density
!> \par History
!>      18-MAR-2002, TCH, working version
!> \note see functionals_utilities
! *****************************************************************************
  SUBROUTINE pz_init ( method, cutoff )

    INTEGER, INTENT(IN)                      :: method
    REAL(KIND=dp), INTENT(IN)                :: cutoff

    CHARACTER(len=*), PARAMETER :: routineN = 'pz_init', &
      routineP = moduleN//':'//routineN

    CALL set_util(cutoff)

    eps_rho = cutoff

    initialized = .FALSE.

    SELECT CASE (method)

    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unknown method")

    CASE (pz_orig)
       CALL cite_reference(Perdew1981)
       ga(0)  = -0.1423_dp;  ga(1)  = -0.0843_dp
       b1(0) = 1.0529_dp;    b1(1) = 1.3981_dp
       b2(0) = 0.3334_dp;    b2(1) = 0.2611_dp
       A(0)  = 0.0311_dp;    A(1)  = 0.01555_dp
       B(0)  = -0.048_dp;    B(1)  = -0.0269_dp
       C(0)  = +0.0020_dp;   C(1)  = +0.0007_dp
       D(0)  = -0.0116_dp;   D(1)  = -0.0048_dp

    CASE (pz_dmc)
       CALL cite_reference(Ortiz1994)
       ga(0)  = -0.103756_dp; ga(1)  = -0.065951_dp
       b1(0) = 0.56371_dp;   b1(1) = 1.11846_dp
       b2(0) = 0.27358_dp;   b2(1) = 0.18797_dp
       A(0)  = 0.031091_dp;  A(1)  = 0.015545_dp
       B(0)  = -0.046644_dp; B(1)  = -0.025599_dp
       C(0)  = -0.00419_dp;  C(1)  = -0.00329_dp
       D(0)  = -0.00983_dp;  D(1)  = -0.00300_dp

    CASE (pz_vmc)
       CALL cite_reference(Ortiz1994)
       ga(0)  = -0.093662_dp; ga(1)  = -0.055331_dp
       b1(0) = 0.49453_dp;   b1(1) = 0.93766_dp
       b2(0) = 0.25534_dp;   b2(1) = 0.14829_dp
       A(0)  = 0.031091_dp;  A(1)  = 0.015545_dp
       B(0)  = -0.046644_dp; B(1)  = -0.025599_dp
       C(0)  = -0.00884_dp;  C(1)  = -0.00677_dp
       D(0)  = -0.00688_dp;  D(1)  = -0.00093_dp

    END SELECT

    initialized = .TRUE.

  END SUBROUTINE pz_init

! *****************************************************************************
!> \brief ...
!> \param r ...
!> \param z ...
!> \param g ...
!> \param order ...
! *****************************************************************************
  SUBROUTINE calc_g(r, z, g, order)

    REAL(KIND=dp), INTENT(IN)                :: r
    INTEGER, INTENT(IN)                      :: z
    REAL(KIND=dp), DIMENSION(0:), &
      INTENT(OUT)                            :: g
    INTEGER, INTENT(IN)                      :: order

    REAL(KIND=dp)                            :: rr, rsr, sr

    IF (r >= 1.0_dp) THEN

       sr = SQRT(r)
       ! order 0 must always be calculated
       g(0) = ga(z) / (1.0_dp + b1(z)*sr + b2(z)*r)
       IF (order >= 1) THEN
          g(1) = -ga(z) * (b1(z)/(2.0_dp*sr) + b2(z)) / &
               (1.0_dp + b1(z)*sr + b2(z)*r)**2
       END IF
       IF (order >= 2) THEN
          rsr = r*sr
          g(2) = &
               2.0_dp * ga(z) * (b1(z)/(2.0_dp*sr) + b2(z))**2 &
               / (1.0_dp + b1(z)*sr + b2(z)*r)**3 &
               + ga(z) * b1(z) &
               / (4.0_dp * (1.0_dp + b1(z)*sr + b2(z)*r)**2 * rsr)
       END IF
       IF (order >= 3) THEN
          g(3) = &
               -6.0_dp * ga(z) * (b1(z)/(2.0_dp*sr) + b2(z))**3 / &
               (1.0_dp + b1(z)*sr + b2(z)*r)**4 &
               -(3.0_dp/2.0_dp) * ga(z) * (b1(z)/(2.0_dp*sr) + b2(z)) * b1(z) / &
               ((1.0_dp + b1(z)*sr + b2(z)*r)**3 * rsr) &
               -(3.0_dp/8.0_dp) * ga(z) * b1(z) / &
               ((1.0_dp + b1(z)*sr + b2(z)*r)**2 * r*rsr)
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
          g(3) = 2.0_dp*A(z)/(rr*r) - C(z)/rr
       END IF

    END IF

  END SUBROUTINE calc_g

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param ed ...
!> \param order ...
!> \param sc ...
! *****************************************************************************
  SUBROUTINE pz_lda_ed_loc ( rho, ed, order, sc )

    REAL(KIND=dp), INTENT(IN)                :: rho
    REAL(KIND=dp), DIMENSION(0:), &
      INTENT(OUT)                            :: ed
    INTEGER, INTENT(IN)                      :: order
    REAL(dp), INTENT(IN)                     :: sc

    INTEGER                                  :: m, order_
    LOGICAL, DIMENSION(0:3)                  :: calc
    REAL(KIND=dp), DIMENSION(0:3)            :: e0, r

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

    IF (order_ >= 1) r(1) = (-1.0_dp/3.0_dp) * r(0) / rho
    IF (order_ >= 2) r(2) = (-4.0_dp/3.0_dp) * r(1) / rho
    IF (order_ >= 3) r(3) = (-7.0_dp/3.0_dp) * r(2) / rho

    m = 0
    IF (calc(0)) THEN
       ed(m) = sc * e0(0)
       m = m + 1
    END IF
    IF (calc(1)) THEN
       ed(m) = sc * e0(1)*r(1)
       m =  m + 1
    END IF
    IF (calc(2)) THEN
       ed(m) = sc * e0(2)*r(1)**2 + sc * e0(1)*r(2)
       m = m + 1
    END IF
    IF (calc(3)) THEN
       ed(m) = sc * e0(3)*r(1)**3 + sc * e0(2)*3.0_dp*r(1)*r(2) + sc * e0(1)*r(3)
    END IF

  END SUBROUTINE pz_lda_ed_loc

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param ed ...
!> \param order ...
!> \param sc ...
! *****************************************************************************
  SUBROUTINE pz_lsd_ed_loc ( a, b, ed, order, sc )

    REAL(KIND=dp), INTENT(IN)                :: a, b
    REAL(KIND=dp), DIMENSION(0:), &
      INTENT(OUT)                            :: ed
    INTEGER, INTENT(IN), OPTIONAL            :: order
    REAL(KIND=dp), INTENT(IN)                :: sc

    INTEGER                                  :: m, order_
    LOGICAL, DIMENSION(0:3)                  :: calc
    REAL(KIND=dp)                            :: rho, tr, trr, trrr, trrz, &
                                                trz, trzz, tz, tzz, tzzz
    REAL(KIND=dp), DIMENSION(0:3)            :: e0, e1, f, r
    REAL(KIND=dp), DIMENSION(0:3, 0:3)       :: z

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
       r(1) = (-1.0_dp/3.0_dp)*r(0)/rho
       tr = e0(1) + (e1(1)-e0(1))*f(0)
       tz = (e1(0)-e0(0))*f(1)
    END IF

!! calculate second partial derivatives
    IF (order_ >= 2) THEN
       r(2) = (-4.0_dp/3.0_dp)*r(1)/rho
       trr = e0(2) + (e1(2)-e0(2))*f(0)
       trz = (e1(1)-e0(1))*f(1)
       tzz = (e1(0)-e0(0))*f(2)
    END IF

!! calculate third derivatives
    IF (order_ >= 3) THEN
       r(3) = (-7.0_dp/3.0_dp)*r(2)/rho
       trrr = e0(3) + (e1(3)-e0(3))*f(0)
       trrz = (e1(2)-e0(2))*f(1)
       trzz = (e1(1)-e0(1))*f(2)
       tzzz = (e1(0)-e0(0))*f(3)
    END IF

    m=0
    IF (calc(0)) THEN
       ed(m) = e0(0) + (e1(0)-e0(0))*f(0)
       ed(m) = ed(m) * sc
       m = m + 1
    END IF

    IF (calc(1)) THEN
       ed(m)   = tr*r(1) + tz*z(1,0)
       ed(m)   = ed(m) * sc
       ed(m+1) = tr*r(1) + tz*z(0,1)
       ed(m+1) = ed(m+1) * sc
       m = m + 2
    END IF

    IF (calc(2)) THEN
       ed(m)   = trr*r(1)**2 + 2.0_dp*trz*r(1)*z(1,0) &
            + tr*r(2) + tzz*z(1,0)**2 + tz*z(2,0)
       ed(m)   = ed(m) * sc
       ed(m+1) = trr*r(1)**2 + trz*r(1)*(z(0,1)+z(1,0)) &
            + tr*r(2) + tzz*z(1,0)*z(0,1) + tz*z(1,1)
       ed(m+1) = ed(m+1) * sc
       ed(m+2) = trr*r(1)**2 + 2.0_dp*trz*r(1)*z(0,1) &
            + tr*r(2) + tzz*z(0,1)**2 + tz*z(0,2)
       ed(m+2) = ed(m+2) * sc
       m = m + 3
    END IF

    IF (calc(3)) THEN
       ed(m)   = trrr*r(1)**3 + 3.0_dp*trrz*r(1)**2*z(1,0) &
            + 3.0_dp*trr*r(1)*r(2) + 3.0_dp*trz*r(2)*z(1,0) &
            + tr*r(3) + 3.0_dp*trzz*r(1)*z(1,0)**2 &
            + tzzz*z(1,0)**3 + 3.0_dp*trz*r(1)*z(2,0) &
            + 3.0_dp*tzz*z(1,0)*z(2,0) + tz*z(3,0)
       ed(m) = ed(m) * sc
       ed(m+1) = trrr*r(1)**3 + trrz*r(1)**2*(2.0_dp*z(1,0)+z(0,1)) &
            + 2.0_dp*trzz*r(1)*z(1,0)*z(0,1) &
            + 2.0_dp*trz*(r(2)*z(1,0)+r(1)*z(1,1)) &
            + 3.0_dp*trr*r(2)*r(1) + trz*r(2)*z(0,1) + tr*r(3) &
            + trzz*r(1)*z(1,0)**2 + tzzz*z(1,0)**2*z(0,1) &
            + 2.0_dp*tzz*z(1,0)*z(1,1) &
            + trz*r(1)*z(2,0) + tzz*z(2,0)*z(0,1) + tz*z(2,1)
       ed(m+1) = ed(m+1) * sc
       ed(m+2) = trrr*r(1)**3 + trrz*r(1)**2*(2.0_dp*z(0,1)+z(1,0)) &
            + 2.0_dp*trzz*r(1)*z(0,1)*z(1,0) &
            + 2.0_dp*trz*(r(2)*z(0,1)+r(1)*z(1,1)) &
            + 3.0_dp*trr*r(2)*r(1) + trz*r(2)*z(1,0) + tr*r(3) &
            + trzz*r(1)*z(0,1)**2 + tzzz*z(0,1)**2*z(1,0) &
            + 2.0_dp*tzz*z(0,1)*z(1,1) &
            + trz*r(1)*z(0,2) + tzz*z(0,2)*z(1,0) + tz*z(1,2)
       ed(m+2) = ed(m+2) * sc
       ed(m+3) = trrr*r(1)**3 + 3.0_dp*trrz*r(1)**2*z(0,1) &
            + 3.0_dp*trr*r(1)*r(2) + 3.0_dp*trz*r(2)*z(0,1) + tr*r(3) &
            + 3.0_dp*trzz*r(1)*z(0,1)**2 + tzzz*z(0,1)**3 &
            + 3.0_dp*trz*r(1)*z(0,2) &
            + 3.0_dp*tzz*z(0,1)*z(0,2) + tz*z(0,3)
       ed(m+3) = ed(m+3) * sc
    END IF

  END SUBROUTINE pz_lsd_ed_loc

END MODULE xc_perdew_zunger
