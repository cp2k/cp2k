!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the LDA functional according to Vosk, Wilk and Nusair
!>      Literature: S. H. Vosko, L. Wilk and M. Nusair,
!>                  Can. J. Phys. 58, 1200 (1980)
!> \note
!>      Order of derivatives is: LDA 0; 1; 2; 3;
!> \par History
!>      JGH (26.02.2003) : OpenMP enabled
!>      fawzi (04.2004)  : adapted to the new xc interface
!>      MG 01.2007       : added scaling
!>      MG 01.2007       : bug fix in vwn_lda_xx
!> \author JGH (26.02.2002)
! *****************************************************************************
MODULE xc_vwn
  USE bibliography,                    ONLY: Vosko1980,&
                                             cite_reference
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_functionals_utilities,        ONLY: calc_srs_pw,&
                                             set_util
  USE xc_input_constants,              ONLY: do_vwn3,&
                                             do_vwn5
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(KIND=dp), PARAMETER :: a   = 0.0310907_dp,&
                              af = 0.01554535_dp,&
                              aa = -0.0168868639404_dp !-1/(6pi^2)

  PUBLIC :: vwn_lda_info, vwn_lda_eval, vwn_lsd_eval

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_vwn'

  REAL(KIND=dp) :: eps_rho, b, c, x0, bf, cf, x0f, ba, ca, x0a

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param cutoff ...
!> \param vwn_params ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE vwn_init ( cutoff, vwn_params, error )

    REAL(KIND=dp), INTENT(IN)                :: cutoff
    TYPE(section_vals_type), POINTER         :: vwn_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'vwn_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: TYPE
    LOGICAL                                  :: failure

    CALL section_vals_val_get(vwn_params,"functional_type",i_val=TYPE,error=error)
    eps_rho = cutoff
    CALL set_util ( cutoff )

    CALL cite_reference(Vosko1980)

    SELECT CASE(TYPE)
       CASE(do_vwn5)
          b    = 3.72744_dp
          c    = 12.9352_dp
          x0   = -0.10498_dp
          bf  = 7.06042_dp
          cf  = 18.0578_dp
          x0f = -0.32500_dp
          ba  = 1.13107_dp
          ca  = 13.0045_dp
          x0a = -0.00475840_dp
       CASE(do_vwn3)
          b  = 13.0720_dp
          c  = 42.7198_dp
          x0 = -0.409286_dp
          bf  = 20.1231_dp
          cf  = 101.578_dp
          x0f = -0.743294_dp
          ba  = 1.13107_dp
          ca  = 13.0045_dp
          x0a = -0.00475840_dp
       CASE DEFAULT
         failure = .FALSE.
         CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " Only functionals VWN3 and VWN5 are supported" ,error,failure)

    END SELECT

  END SUBROUTINE vwn_init

! *****************************************************************************
!> \brief ...
!> \param reference ...
!> \param shortform ...
!> \param needs ...
!> \param max_deriv ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE vwn_lda_info ( reference, shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout)       :: error

    IF ( PRESENT ( reference ) ) THEN
       reference = "S. H. Vosko, L. Wilk and M. Nusair,"// &
                   " Can. J. Phys. 58, 1200 (1980) {LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
       shortform = "Vosko-Wilk-Nusair Functional {LDA}"
    END IF
    IF (PRESENT(needs)) THEN
       needs%rho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE vwn_lda_info

! *****************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order ...
!> \param vwn_params ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE vwn_lda_eval(rho_set,deriv_set,order,vwn_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: order
    TYPE(section_vals_type), POINTER         :: vwn_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vwn_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: epsilon_rho, sc
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: x
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: e_0, e_rho, e_rho_rho, &
                                                e_rho_rho_rho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(bo)

    CALL section_vals_val_get(vwn_params,"scale_c",r_val=sc,error=error)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CALL xc_rho_set_get(rho_set,rho=rho,&
         local_bounds=bo,rho_cutoff=epsilon_rho,&
         error=error)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)
    CALL vwn_init(epsilon_rho, vwn_params, error)

    ALLOCATE ( x(npoints), STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL calc_srs_pw ( rho, x, npoints )

    IF (order>=1) THEN
       deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
       deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)

       CALL vwn_lda_01 ( rho, x, e_0,e_rho,npoints, sc, error )
    ELSEIF ( order>=0 ) THEN
       deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE., error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)

       CALL vwn_lda_0 ( rho, x, e_0, npoints, sc, error )
    ELSEIF ( order==-1 ) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)

       CALL vwn_lda_1 ( rho, x, e_rho, npoints, sc, error )
    END IF
    IF ( order>=2.OR.order==-2 ) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)

       CALL vwn_lda_2 ( rho, x, e_rho_rho, npoints, sc, error )
    END IF
    IF ( order>=3.OR.order==-3 ) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)

       CALL vwn_lda_3 ( rho, x, e_rho_rho_rho, npoints, sc,&
            error)
    END IF
    IF ( order>3.OR.order<-3) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="derivatives bigger than 3 not implemented", &
            error=error, error_level=cp_failure_level)
    END IF

    DEALLOCATE(x,stat=stat)
    CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    CALL timestop(handle)

  END SUBROUTINE vwn_lda_eval

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param x ...
!> \param e_0 ...
!> \param npoints ...
!> \param sc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE vwn_lda_0 ( rho, x, e_0, npoints, sc, error )

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, x
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_0
    INTEGER, INTENT(in)                      :: npoints
    REAL(KIND=dp)                            :: sc
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: at, dpx, ln1, ln2, px, px0, &
                                                q, xb

    q = SQRT ( 4.0_dp*c - b*b )
    xb = 2.0_dp*x0+b
    px0 = x0*x0+b*x0+c

    !$omp parallel do default(none) &
    !$omp             shared (npoints, rho, eps_rho, x, c, b, e_0, sc) &
    !$omp             shared (q, x0, px0, xb) &
    !$omp             private(ip,px,dpx,at,ln1,ln2)

    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2.0_dp*x(ip)+b
         at = 2.0_dp/q*ATAN(q/dpx)
         ln1 = LOG(x(ip)*x(ip)/px)
         ln2 = LOG((x(ip)-x0)**2/px)
         e_0(ip) = e_0(ip) + a * ( ln1 + b*at - b*x0/px0* ( ln2 + xb * at ) ) * rho(ip) * sc
      END IF

    END DO

    !$omp end parallel do

  END SUBROUTINE vwn_lda_0

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param x ...
!> \param e_rho ...
!> \param npoints ...
!> \param sc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE vwn_lda_1 ( rho, x, e_rho, npoints, sc, error )

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, x
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_rho
    INTEGER, INTENT(in)                      :: npoints
    REAL(KIND=dp)                            :: sc
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: at, dat, dex, dln1, dln2, &
                                                dpx, ex, ln1, ln2, pa, px, &
                                                px0, q, xb

    q = SQRT ( 4.0_dp*c - b*b )
    xb = 2.0_dp*x0+b
    px0 = x0*x0+b*x0+c

    !$omp parallel do default(none) &
    !$omp             shared(npoints, rho, eps_rho, x, b, sc, e_rho) &
    !$omp             shared(c, x0, q, xb, px0) &
    !$omp             private(ip,px,dpx,at,pa,dat,ln1,dln1,ln2,dln2,ex,dex)

    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2.0_dp*x(ip)+b
         at = 2.0_dp/q*ATAN(q/dpx)
         pa = 4.0_dp*x(ip)*x(ip)+4.0_dp*b*x(ip)+b*b+q*q
         dat= -4.0_dp/pa
         ln1 = LOG(x(ip)*x(ip)/px)
         dln1 = (b*x(ip)+2.0_dp*c)/(x(ip)*px)
         ln2 = LOG((x(ip)-x0)**2/px)
         dln2 = (b*x(ip)+2.0_dp*c+2.0_dp*x0*x(ip)+x0*b)/((x(ip)-x0)*px)
         ex = a * ( ln1 + b*at - b*x0/px0* ( ln2 + xb * at ) )
         dex = a * ( dln1 + b*dat  - b*x0/px0* ( dln2 + xb * dat ))
         e_rho(ip) = e_rho(ip) +( ex - x(ip)*dex/6.0_dp ) * sc
      END IF

    END DO

    !$omp end parallel do

  END SUBROUTINE vwn_lda_1

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param x ...
!> \param e_0 ...
!> \param e_rho ...
!> \param npoints ...
!> \param sc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE vwn_lda_01 ( rho, x, e_0, e_rho, npoints, sc, error )

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, x
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_0, e_rho
    INTEGER, INTENT(in)                      :: npoints
    REAL(KIND=dp)                            :: sc
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: at, dat, dex, dln1, dln2, &
                                                dpx, ex, ln1, ln2, pa, px, &
                                                px0, q, xb

    q = SQRT ( 4.0_dp*c - b*b )
    xb = 2.0_dp*x0+b
    px0 = x0*x0+b*x0+c

    !$omp parallel do default(none) &
    !$omp             shared(npoints, rho, eps_rho, x, b, c, e_0, sc) &
    !$omp             shared(x0, q, xb, px0, e_rho) &
    !$omp             private(ip,px,dpx,at,pa,dat,ln1,dln1,ln2,dln2,ex,dex)

    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2.0_dp*x(ip)+b
         at = 2.0_dp/q*ATAN(q/dpx)
         pa = 4.0_dp*x(ip)*x(ip)+4.0_dp*b*x(ip)+b*b+q*q
         dat= -4.0_dp/pa
         ln1 = LOG(x(ip)*x(ip)/px)
         dln1 = (b*x(ip)+2.0_dp*c)/(x(ip)*px)
         ln2 = LOG((x(ip)-x0)**2/px)
         dln2 = (x(ip)*(b+2.0_dp*x0)+2.0_dp*c+x0*b)/((x(ip)-x0)*px)
         ex = a * ( ln1 + b*at - b*x0/px0* ( ln2 + xb * at ) )
         dex = a * ( dln1 + b*dat  - b*x0/px0* ( dln2 + xb * dat ) )
         e_0(ip) = e_0(ip) + ex * rho(ip) * sc
         e_rho(ip) = e_rho(ip) + ( ex - x(ip)*dex/6.0_dp) * sc
      END IF

    END DO

    !$omp end parallel do

  END SUBROUTINE vwn_lda_01

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param x ...
!> \param e_rho_rho ...
!> \param npoints ...
!> \param sc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE vwn_lda_2 ( rho, x, e_rho_rho, npoints, sc, error )

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, x
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_rho_rho
    INTEGER, INTENT(in)                      :: npoints
    REAL(KIND=dp)                            :: sc
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: at, d2at, d2ex, d2ln1, d2ln2, &
                                                dat, dex, dln1, dln2, dpx, &
                                                fp, ln1, ln2, pa, px, px0, q, &
                                                xb

    q = SQRT ( 4.0_dp*c - b*b )
    xb = 2.0_dp*x0+b
    px0 = x0*x0+b*x0+c
    fp = -b*x0/px0

    !$omp parallel do default(none) &
    !$omp             private(ip,px,dpx,at,pa,dat,d2at,ln1,dln1,d2ln1) &
    !$omp             private(ln2,dln2,d2ln2,dex,d2ex) &
    !$omp             shared(npoints, rho, fp, eps_rho, x, b, c, q, x0) &
    !$omp             shared(xb, e_rho_rho, sc)

    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2.0_dp*x(ip)+b
         at = 2.0_dp/q*ATAN(q/dpx)
         pa = 4.0_dp*x(ip)*x(ip)+4.0_dp*b*x(ip)+b*b+q*q
         dat= -4.0_dp/pa
         d2at= 16.0_dp*dpx/(pa*pa)
         ln1 = LOG(x(ip)*x(ip)/px)
         dln1 = (b*x(ip)+2.0_dp*c)/(x(ip)*px)
         d2ln1 = b/(x(ip)*px) - (b*x(ip)+2.0_dp*c)/(x(ip)*px)**2 * (px+x(ip)*dpx)
         ln2 = LOG((x(ip)-x0)**2/px)
         dln2 = (x(ip)*xb+2.0_dp*c+x0*b)/((x(ip)-x0)*px)
         d2ln2 = xb/((x(ip)-x0)*px)-(x(ip)*xb+2.0_dp*c+x0*b)/((x(ip)-x0)*px)**2 &
                 * (px+(x(ip)-x0)*dpx)
         dex = a * ( dln1 + b*dat + fp * ( dln2 + xb * dat ) )
         d2ex= a * ( d2ln1 + b*d2at  + fp * ( d2ln2 + xb * d2at ) )
         e_rho_rho(ip) = e_rho_rho(ip)&
              + ( x(ip)/(36.0_dp*rho(ip))*(x(ip)*d2ex-5.0_dp*dex) ) * sc
      END IF

    END DO
    !$omp end parallel do

  END SUBROUTINE vwn_lda_2

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param x ...
!> \param e_rho_rho_rho ...
!> \param npoints ...
!> \param sc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE vwn_lda_3 ( rho, x, e_rho_rho_rho, npoints, sc, error )

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, x
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_rho_rho_rho
    INTEGER, INTENT(in)                      :: npoints
    REAL(KIND=dp)                            :: sc
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp) :: at, ax, bx, cx, d2at, d2bx, d2dx, d2ex, d2ln1, d2ln2, &
      d3at, d3ex, d3ln1, d3ln2, dat, dbx, ddx, dex, dln1, dln2, dpx, dx, fp, &
      ln1, ln2, pa, px, px0, q, xb

    q = SQRT ( 4.0_dp*c - b*b )
    xb = 2.0_dp*x0+b
    px0 = x0*x0+b*x0+c
    fp = -b*x0/px0

    !$omp parallel do default(none) &
    !$omp             shared(npoints, rho, eps_rho, b, c, x, x0, sc) &
    !$omp             shared(q, xb, px0, fp, e_rho_rho_rho) &
    !$omp             private(ip,px,dpx,at,pa,dat,d2at,d3at,ln1,ax,bx) &
    !$omp             private(dbx,d2bx,dln1,d2ln1,d3ln1,ln2,cx,dx,ddx,d2dx) &
    !$omp             private(dln2,d2ln2,d3ln2,dex,d2ex,d3ex)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2.0_dp*x(ip)+b
         at = 2.0_dp/q*ATAN(q/dpx)
         pa = 4.0_dp*x(ip)*x(ip)+4.0_dp*b*x(ip)+b*b+q*q
         dat= -4.0_dp/pa
         d2at= 16.0_dp*dpx/(pa*pa)
         d3at= 32.0_dp/(pa*pa)*(1.0_dp-4.0_dp*dpx*dpx/pa)
         ln1 = LOG(x(ip)*x(ip)/px)
         ax = b*x(ip)+2.0_dp*c
         bx = x(ip)*px
         dbx = px + x(ip)*dpx
         d2bx = 2.0_dp*(dpx+x(ip))
         dln1 = ax/bx
         d2ln1 = (b*bx-ax*dbx)/(bx*bx)
         d3ln1 = -ax*d2bx/(bx*bx) - 2.0_dp*d2ln1*dbx/bx
         ln2 = LOG((x(ip)-x0)**2/px)
         cx = x(ip)*xb+2.0_dp*c+x0*b
         dx = (x(ip)-x0)*px
         ddx = px + (x(ip)-x0)*dpx
         d2dx = 2.0_dp*(dpx + (x(ip)-x0))
         dln2 = cx/dx
         d2ln2 = (xb*dx-cx*ddx)/(dx*dx)
         d3ln2 = -cx*d2dx/(dx*dx) - 2.0_dp*d2ln2*ddx/dx
         dex = a * ( dln1 + b*dat + fp * ( dln2 + xb * dat ) )
         d2ex= a * ( d2ln1 + b*d2at + fp * ( d2ln2 + xb * d2at ) )
         d3ex= a * ( d3ln1 + b*d3at + fp * ( d3ln2 + xb * d3at ) )
         e_rho_rho_rho(ip) = e_rho_rho_rho(ip)&
              - ( 7.0_dp*x(ip)/(216.0_dp*rho(ip)*rho(ip))*(x(ip)*d2ex-5.0_dp*dex) + &
              x(ip)*x(ip)/(216.0_dp*rho(ip)*rho(ip))*(x(ip)*d3ex-4.0_dp*d2ex) ) * sc
      END IF

    END DO
    !$omp end parallel do

  END SUBROUTINE vwn_lda_3

! *****************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order ...
!> \param vwn_params ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE vwn_lsd_eval(rho_set,deriv_set,order,vwn_params,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: order
    TYPE(section_vals_type), POINTER         :: vwn_params
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vwn_lsd_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat, TYPE
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: epsilon_rho, sc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dummy, e_0, e_a, e_aa, e_aaa, &
                                                e_aab, e_ab, e_abb, e_b, &
                                                e_bb, e_bbb, rhoa, rhob
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(bo)

    CALL section_vals_val_get(vwn_params,"scale_c",r_val=sc,error=error)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CALL xc_rho_set_get(rho_set,rhoa=rhoa,rhob=rhob,&
         local_bounds=bo,rho_cutoff=epsilon_rho,&
         error=error)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)
    CALL vwn_init(epsilon_rho, vwn_params, error)

    IF (cp_debug) THEN
      ALLOCATE(dummy(bo(1,1):bo(2,1),bo(1,2):bo(2,2),bo(1,3):bo(2,3)),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      dummy=> rhoa
    END IF

    e_0 => dummy
    e_a => dummy
    e_b => dummy
    e_aa => dummy
    e_bb => dummy
    e_ab => dummy
    e_aaa => dummy
    e_bbb => dummy
    e_aab => dummy
    e_abb => dummy

    IF (order>=0) THEN
      deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE., error=error)
      CALL xc_derivative_get(deriv, deriv_data=e_0,error=error)
    END IF
    IF (order>=1.OR.order==-1) THEN
      deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",&
            allocate_deriv=.TRUE.,error=error)
      CALL xc_derivative_get(deriv,deriv_data=e_a,error=error)

      deriv => xc_dset_get_derivative(deriv_set,"(rhob)",&
            allocate_deriv=.TRUE.,error=error)
       CALL xc_derivative_get(deriv,deriv_data=e_b,error=error)
    END IF
    IF (order>=2.OR.order==-2) THEN
      deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)",&
            allocate_deriv=.TRUE.,error=error)
      CALL xc_derivative_get(deriv,deriv_data=e_aa,error=error)

      deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)",&
            allocate_deriv=.TRUE.,error=error)
      CALL xc_derivative_get(deriv,deriv_data=e_bb,error=error)

      deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)",&
            allocate_deriv=.TRUE.,error=error)
      CALL xc_derivative_get(deriv,deriv_data=e_ab,error=error)
    END IF
    IF (order>=3.OR.order==-3) THEN
      deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)(rhoa)",&
            allocate_deriv=.TRUE.,error=error)
      CALL xc_derivative_get(deriv,deriv_data=e_aaa,error=error)

      deriv => xc_dset_get_derivative(deriv_set,"(rhob)(rhob)(rhob)",&
            allocate_deriv=.TRUE.,error=error)
      CALL xc_derivative_get(deriv,deriv_data=e_bbb,error=error)

      deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhoa)(rhob)",&
            allocate_deriv=.TRUE.,error=error)
      CALL xc_derivative_get(deriv,deriv_data=e_aab,error=error)

      deriv => xc_dset_get_derivative(deriv_set,"(rhoa)(rhob)(rhob)",&
            allocate_deriv=.TRUE.,error=error)
      CALL xc_derivative_get(deriv,deriv_data=e_abb,error=error)

    END IF
    IF ( order>3.OR.order<-3) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
           message="derivatives bigger than 3 not implemented", &
           error=error, error_level=cp_failure_level)
    END IF

    CALL section_vals_val_get(vwn_params,"functional_type",i_val=TYPE,error=error)
    SELECT CASE(TYPE)
      CASE(do_vwn3)

        !$omp parallel default(none) &
        !$omp          shared(rhoa, rhob, e_0, e_a, e_b, e_aa, e_bb, e_ab) &
        !$omp          shared(e_aaa, e_bbb, e_aab, e_abb, order, npoints) &
        !$omp          shared(sc, error)

        CALL vwn3_lsd_calc(&
               rhoa=rhoa, rhob=rhob, e_0=e_0, &
               e_a=e_a, e_b=e_b,&
               e_aa=e_aa, e_bb=e_bb,e_ab=e_ab,&
               e_aaa=e_aaa, e_bbb=e_bbb, e_aab=e_aab, e_abb=e_abb,&
               order=order, npoints=npoints, &
               sc=sc, error=error)

        !$omp end parallel

      CASE(do_vwn5)

        !$omp parallel default(none) &
        !$omp          shared(rhoa, rhob, e_0, e_a, e_b, e_aa, e_bb, e_ab) &
        !$omp          shared(e_aaa, e_bbb, e_aab, e_abb, order, npoints) &
        !$omp          shared(sc, error)

        CALL vwn5_lsd_calc(&
               rhoa=rhoa, rhob=rhob, e_0=e_0, &
               e_a=e_a, e_b=e_b,&
               e_aa=e_aa, e_bb=e_bb,e_ab=e_ab,&
               e_aaa=e_aaa, e_bbb=e_bbb, e_aab=e_aab, e_abb=e_abb,&
               order=order, npoints=npoints, &
               sc=sc, error=error)

        !$omp end parallel

      CASE DEFAULT
        failure = .FALSE.
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                      " Only functionals VWN3 and VWN5 are supported" ,error,failure)
    END SELECT

    IF (cp_debug) THEN
      DEALLOCATE(dummy,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    ELSE
      NULLIFY(dummy)
    END IF

    CALL timestop(handle)

  END SUBROUTINE vwn_lsd_eval

! *****************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param rhob ...
!> \param e_0 ...
!> \param e_a ...
!> \param e_b ...
!> \param e_aa ...
!> \param e_bb ...
!> \param e_ab ...
!> \param e_aaa ...
!> \param e_bbb ...
!> \param e_aab ...
!> \param e_abb ...
!> \param order ...
!> \param npoints ...
!> \param sc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE vwn3_lsd_calc ( rhoa, rhob, e_0, e_a, e_b, e_aa, e_bb, e_ab,&
                            e_aaa, e_bbb, e_aab, e_abb, order, npoints, sc, error )

    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: rhoa, rhob, e_0, e_a, e_b, &
                                                e_aa, e_bb, e_ab, e_aaa, &
                                                e_bbb, e_aab, e_abb
    INTEGER, INTENT(in)                      :: order, npoints
    REAL(KIND=dp)                            :: sc
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp) :: ap, bp, cp, myrho, myrhoa, myrhob, Qf, Qp, t1, t10, &
      t100, t101, t1019, t102, t1031, t105, t108, t109, t11, t110, t113, &
      t114, t115, t116, t117, t119, t120, t121, t124, t125, t126, t127, t130, &
      t132, t133, t134, t136, t137, t14, t144, t145, t146, t147, t15, t150, &
      t151, t154, t155, t156, t159, t16, t160, t161, t162, t165, t168, t169, &
      t172, t173, t174, t175, t176, t178, t179, t180, t183, t184, t187, t189, &
      t19, t190, t191, t193, t194, t2, t202, t203, t206, t209, t21, t212, &
      t213, t216, t217, t218, t219, t220, t221, t222, t223, t225, t226, t230, &
      t231, t235, t236, t237, t24, t240
    REAL(KIND=dp) :: t241, t242, t244, t245, t246, t25, t251, t254, t255, &
      t258, t259, t26, t260, t261, t267, t268, t269, t27, t270, t273, t276, &
      t279, t281, t282, t283, t284, t286, t287, t29, t291, t292, t295, t296, &
      t299, t3, t302, t306, t307, t309, t310, t311, t315, t316, t319, t32, &
      t324, t325, t332, t333, t334, t337, t339, t342, t343, t346, t349, t350, &
      t351, t352, t353, t354, t356, t357, t36, t363, t364, t365, t368, t373, &
      t376, t377, t38, t380, t381, t382, t388, t389, t390, t391, t394, t397, &
      t4, t400, t402, t403, t404, t405, t407, t408, t412, t413, t42, t420, &
      t424, t425, t427, t428, t429, t43
    REAL(KIND=dp) :: t433, t434, t437, t44, t442, t443, t45, t451, t452, &
      t455, t457, t458, t46, t461, t464, t467, t470, t471, t472, t475, t479, &
      t48, t482, t485, t488, t489, t49, t490, t494, t497, t499, t5, t500, &
      t501, t504, t505, t508, t509, t51, t510, t513, t516, t519, t522, t525, &
      t526, t528, t531, t536, t538, t54, t540, t541, t549, t55, t559, t563, &
      t565, t567, t570, t573, t58, t583, t589, t59, t595, t596, t6, t601, &
      t603, t613, t62, t625, t64, t665, t67, t68, t69, t690, t693, t696, t7, &
      t705, t71, t710, t719, t720, t721, t727, t729, t732, t74, t743, t745, &
      t748, t752, t755, t758, t769, t78, t792
    REAL(KIND=dp) :: t793, t798, t80, t800, t804, t807, t808, t810, t814, &
      t820, t823, t835, t836, t838, t841, t85, t851, t853, t86, t860, t863, &
      t872, t879, t88, t89, t90, t91, t92, t947, t95, t950, t953, t956, t96, &
      t967, t97, t98, t984, t986, t990, x0p

    cp  = c
    ap  = a
    bp  = b
    x0p = x0
    Qp = SQRT(4.0_dp*cp-bp*bp)
    Qf = SQRT(4.0_dp*cf-bf*bf)

    !$omp do
    DO ip = 1, npoints
      myrhoa = MAX(rhoa(ip),0.0_dp)
      myrhob = MAX(rhob(ip),0.0_dp)
      myrho  = myrhoa + myrhob
      IF (myrho > eps_rho ) THEN
        IF (order>=0) THEN
          t1 = 0.1e1_dp / 0.3141592654e1_dp
          t2 = myrho
          t3 = 0.1e1_dp / t2
          t4 = t1 * t3
          t5 = t4 ** 0.3333333332e0_dp
          t6 = 0.9085602964e0_dp * t5
          t7 = t4 ** 0.1666666666e0_dp
          t10 = t6 + 0.9531842930e0_dp * bp * t7 + cp
          t11 = 0.1e1_dp / t10
          t14 = LOG(0.9085602964e0_dp * t5 * t11)
          t15 = 0.1906368586e1_dp * t7
          t16 = t15 + bp
          t19 = ATAN(Qp / t16)
          t21 = 0.1e1_dp / Qp
          t24 = bp * x0p
          t25 = 0.9531842930e0_dp * t7
          t26 = t25 - x0p
          t27 = t26 ** 0.20e1_dp
          t29 = LOG(t27 * t11)
          t32 = 0.20e1_dp * bp + 0.40e1_dp * x0p
          t36 = x0p ** 0.20e1_dp
          t38 = 0.1e1_dp / (t36 + t24 + cp)
          t42 = ap * (t14 + 0.20e1_dp * bp * t19 * t21 - t24 * (t29 + t32 * t19  &
                * t21) * t38)
          t43 = myrhoa - myrhob
          t44 = t43 * t3
          t45 = 0.10e1_dp + t44
          t46 = t45 ** 0.1333333333e1_dp
          t48 = 0.10e1_dp - t44
          t49 = t48 ** 0.1333333333e1_dp
          t51 = 0.1923661050e1_dp * t46 + 0.1923661050e1_dp * t49 - 0.3847322100e1_dp
          t54 = t6 + 0.9531842930e0_dp * bf * t7 + cf
          t55 = 0.1e1_dp / t54
          t58 = LOG(0.9085602964e0_dp * t5 * t55)
          t59 = t15 + bf
          t62 = ATAN(Qf / t59)
          t64 = 0.1e1_dp / Qf
          t67 = bf * x0f
          t68 = t25 - x0f
          t69 = t68 ** 0.20e1_dp
          t71 = LOG(t69 * t55)
          t74 = 0.20e1_dp * bf + 0.40e1_dp * x0f
          t78 = x0f ** 0.20e1_dp
          t80 = 0.1e1_dp / (t78 + t67 + cf)
          t85 = af * (t58 + 0.20e1_dp * bf * t62 * t64 - t67 * (t71 + t74 * t62  &
                * t64) * t80) - t42
          t86 = t51 * t85

          e_0(ip) = e_0(ip) + ( (t42 + t86) * t2 ) * sc

        END IF
        IF (order>=1.OR.order==-1) THEN
          t88 = t4 ** (-0.6666666668e0_dp)
          t89 = t88 * t11
          t90 = t2 ** 2
          t91 = 0.1e1_dp / t90
          t92 = t1 * t91
          t95 = t10 ** 2
          t96 = 0.1e1_dp / t95
          t97 = t5 * t96
          t98 = t88 * t1
          t100 = 0.3028534320e0_dp * t98 * t91
          t101 = t4 ** (-0.8333333334e0_dp)
          t102 = bp * t101
          t105 = -t100 - 0.1588640488e0_dp * t102 * t92
          t108 = -0.3028534320e0_dp * t89 * t92 - 0.9085602964e0_dp * t97 * t105
          t109 = t4 ** (-0.3333333332e0_dp)
          t110 = t108 * t109
          t113 = t16 ** 2
          t114 = 0.1e1_dp / t113
          t115 = bp * t114
          t116 = t115 * t101
          t117 = Qp ** 2
          t119 = 0.1e1_dp + t117 * t114
          t120 = 0.1e1_dp / t119
          t121 = t92 * t120
          t124 = t26 ** 0.10e1_dp
          t125 = t124 * t11
          t126 = t101 * t1
          t127 = t126 * t91
          t130 = t27 * t96
          t132 = -0.3177280976e0_dp * t125 * t127 - t130 * t105
          t133 = t26 ** (-0.20e1_dp)
          t134 = t132 * t133
          t136 = t32 * t114
          t137 = t136 * t101
          t144 = ap * (0.1100642416e1_dp * t110 * t10 + 0.6354561950e0_dp * t116 *  &
                 t121 - t24 * (t134 * t10 + 0.3177280975e0_dp * t137 * t121) * t38)
          t145 = t45 ** 0.333333333e0_dp
          t146 = t43 * t91
          t147 = t3 - t146
          t150 = t48 ** 0.333333333e0_dp
          t151 = -t147
          t154 = 0.2564881399e1_dp * t145 * t147 + 0.2564881399e1_dp * t150 * t151
          t155 = t154 * t85
          t156 = t88 * t55
          t159 = t54 ** 2
          t160 = 0.1e1_dp / t159
          t161 = t5 * t160
          t162 = bf * t101
          t165 = -t100 - 0.1588640488e0_dp * t162 * t92
          t168 = -0.3028534320e0_dp * t156 * t92 - 0.9085602964e0_dp * t161 * t165
          t169 = t168 * t109
          t172 = t59 ** 2
          t173 = 0.1e1_dp / t172
          t174 = bf * t173
          t175 = t174 * t101
          t176 = Qf ** 2
          t178 = 0.1e1_dp + t176 * t173
          t179 = 0.1e1_dp / t178
          t180 = t92 * t179
          t183 = t68 ** 0.10e1_dp
          t184 = t183 * t55
          t187 = t69 * t160
          t189 = -0.3177280976e0_dp * t184 * t127 - t187 * t165
          t190 = t68 ** (-0.20e1_dp)
          t191 = t189 * t190
          t193 = t74 * t173
          t194 = t193 * t101
          t202 = af * (0.1100642416e1_dp * t169 * t54 + 0.6354561950e0_dp * t175 *  &
                 t180 - t67 * (t191 * t54 + 0.3177280975e0_dp * t194 * t180) * t80) -  &
                 t144
          t203 = t51 * t202

          e_a(ip) = e_a(ip) + ( (t144 + t155 + t203) * t2 + t42 + t86 ) * sc

          t206 = -t3 - t146
          t209 = -t206
          t212 = 0.2564881399e1_dp * t145 * t206 + 0.2564881399e1_dp * t150 * t209
          t213 = t212 * t85

          e_b(ip) = e_b(ip) + ( (t144 + t213 + t203) * t2 + t42 + t86 ) * sc

        END IF
        IF (order>=2.OR.order==-2) THEN
          t216 = t4 ** (-0.1666666667e1_dp)
          t217 = t216 * t11
          t218 = 0.3141592654e1_dp ** 2
          t219 = 0.1e1_dp / t218
          t220 = t90 ** 2
          t221 = 0.1e1_dp / t220
          t222 = t219 * t221
          t223 = t217 * t222
          t225 = t88 * t96
          t226 = t92 * t105
          t230 = 0.1e1_dp / t90 / t2
          t231 = t1 * t230
          t235 = 0.1e1_dp / t95 / t10
          t236 = t5 * t235
          t237 = t105 ** 2
          t240 = t216 * t219
          t241 = t240 * t221
          t242 = 0.2019022880e0_dp * t241
          t244 = 0.6057068640e0_dp * t98 * t230
          t245 = t4 ** (-0.1833333333e1_dp)
          t246 = bp * t245
          t251 = -t242 + t244 - 0.1323867073e0_dp * t246 * t222 + 0.3177280976e0_dp  &
                 * t102 * t231
          t254 = -0.2019022880e0_dp * t223 + 0.6057068640e0_dp * t225 * t226 + 0.6057068640e0_dp  &
                 * t89 * t231 + 0.1817120593e1_dp * t236 * t237 - 0.9085602964e0_dp  &
                 * t97 * t251
          t255 = t254 * t109
          t258 = t4 ** (-0.1333333333e1_dp)
          t259 = t108 * t258
          t260 = t10 * t1
          t261 = t260 * t91
          t267 = 0.1e1_dp / t113 / t16
          t268 = bp * t267
          t269 = t268 * t216
          t270 = t222 * t120
          t273 = t115 * t245
          t276 = t231 * t120
          t279 = t113 ** 2
          t281 = 0.1e1_dp / t279 / t16
          t282 = bp * t281
          t283 = t282 * t216
          t284 = t119 ** 2
          t286 = 0.1e1_dp / t284 * t117
          t287 = t222 * t286
          t291 = t124 * t96
          t292 = t291 * t101
          t295 = t245 * t219
          t296 = t295 * t221
          t299 = t126 * t230
          t302 = t27 * t235
          t306 = 0.5047557200e-1_dp * t223 + 0.6354561952e0_dp * t292 * t226 - 0.2647734147e0_dp  &
                 * t125 * t296 + 0.6354561952e0_dp * t125 * t299 + 0.2e1_dp *  &
                 t302 * t237 - t130 * t251
          t307 = t306 * t133
          t309 = t26 ** (-0.30e1_dp)
          t310 = t132 * t309
          t311 = t310 * t10
          t315 = t32 * t267
          t316 = t315 * t216
          t319 = t136 * t245
          t324 = t32 * t281
          t325 = t324 * t216
          t332 = ap * (0.1100642416e1_dp * t255 * t10 + 0.3668808052e0_dp * t259 *  &
                 t261 + 0.1100642416e1_dp * t110 * t105 + 0.4038045758e0_dp * t269 * t270  &
                 + 0.5295468292e0_dp * t273 * t270 - 0.1270912390e1_dp * t116 * t276 - 0.4038045758e0_dp  &
                 * t283 * t287 - t24 * (t307 * t10 + 0.3177280976e0_dp *  &
                 t311 * t127 + t134 * t105 + 0.2019022879e0_dp * t316 * t270 + 0.2647734146e0_dp  &
                 * t319 * t270 - 0.6354561950e0_dp * t137 * t276 - 0.2019022879e0_dp  &
                 * t325 * t287) * t38)
          t333 = t45 ** (-0.666666667e0_dp)
          t334 = t147 ** 2
          t337 = t43 * t230
          t339 = -0.2e1_dp * t91 + 0.2e1_dp * t337
          t342 = t48 ** (-0.666666667e0_dp)
          t343 = t151 ** 2
          t346 = -t339
          t349 = 0.8549604655e0_dp * t333 * t334 + 0.2564881399e1_dp * t145 * t339  &
                 + 0.8549604655e0_dp * t342 * t343 + 0.2564881399e1_dp * t150 * t346
          t350 = t349 * t85
          t351 = t154 * t202
          t352 = 0.2e1_dp * t351
          t353 = t216 * t55
          t354 = t353 * t222
          t356 = t88 * t160
          t357 = t92 * t165
          t363 = 0.1e1_dp / t159 / t54
          t364 = t5 * t363
          t365 = t165 ** 2
          t368 = bf * t245
          t373 = -t242 + t244 - 0.1323867073e0_dp * t368 * t222 + 0.3177280976e0_dp  &
                 * t162 * t231
          t376 = -0.2019022880e0_dp * t354 + 0.6057068640e0_dp * t356 * t357 + 0.6057068640e0_dp  &
                 * t156 * t231 + 0.1817120593e1_dp * t364 * t365 - 0.9085602964e0_dp  &
                 * t161 * t373
          t377 = t376 * t109
          t380 = t168 * t258
          t381 = t54 * t1
          t382 = t381 * t91
          t388 = 0.1e1_dp / t172 / t59
          t389 = bf * t388
          t390 = t389 * t216
          t391 = t222 * t179
          t394 = t174 * t245
          t397 = t231 * t179
          t400 = t172 ** 2
          t402 = 0.1e1_dp / t400 / t59
          t403 = bf * t402
          t404 = t403 * t216
          t405 = t178 ** 2
          t407 = 0.1e1_dp / t405 * t176
          t408 = t222 * t407
          t412 = t183 * t160
          t413 = t412 * t101
          t420 = t69 * t363
          t424 = 0.5047557200e-1_dp * t354 + 0.6354561952e0_dp * t413 * t357 - 0.2647734147e0_dp  &
                 * t184 * t296 + 0.6354561952e0_dp * t184 * t299 + 0.2e1_dp *  &
                 t420 * t365 - t187 * t373
          t425 = t424 * t190
          t427 = t68 ** (-0.30e1_dp)
          t428 = t189 * t427
          t429 = t428 * t54
          t433 = t74 * t388
          t434 = t433 * t216
          t437 = t193 * t245
          t442 = t74 * t402
          t443 = t442 * t216
          t451 = af * (0.1100642416e1_dp * t377 * t54 + 0.3668808052e0_dp * t380 *  &
                 t382 + 0.1100642416e1_dp * t169 * t165 + 0.4038045758e0_dp * t390 * t391  &
                 + 0.5295468292e0_dp * t394 * t391 - 0.1270912390e1_dp * t175 * t397 - 0.4038045758e0_dp  &
                 * t404 * t408 - t67 * (t425 * t54 + 0.3177280976e0_dp *  &
                 t429 * t127 + t191 * t165 + 0.2019022879e0_dp * t434 * t391 + 0.2647734146e0_dp  &
                 * t437 * t391 - 0.6354561950e0_dp * t194 * t397 - 0.2019022879e0_dp  &
                 * t443 * t408) * t80) - t332
          t452 = t51 * t451
          t455 = 0.2e1_dp * t144
          t457 = 0.2e1_dp * t203

          e_aa(ip) = e_aa(ip) + ( (t332 + t350 + t352 + t452) * t2 + t455 + 0.2e1_dp * t155 + t457 ) * sc

          t458 = t333 * t147
          t461 = t145 * t43
          t464 = t342 * t151
          t467 = t150 * t43
          t470 = 0.8549604655e0_dp * t458 * t206 + 0.5129762798e1_dp * t461 * t230  &
                 + 0.8549604655e0_dp * t464 * t209 - 0.5129762798e1_dp * t467 * t230
          t471 = t470 * t85
          t472 = t212 * t202

          e_ab(ip) = e_ab(ip) + ( (t332 + t471 + t351 + t472 + t452) * t2 + t455 + t155 + t457  &
                     + t213 ) * sc
          t475 = t206 ** 2
          t479 = 0.2e1_dp * t91 + 0.2e1_dp * t337
          t482 = t209 ** 2
          t485 = -t479
          t488 = 0.8549604655e0_dp * t333 * t475 + 0.2564881399e1_dp * t145 * t479  &
                 + 0.8549604655e0_dp * t342 * t482 + 0.2564881399e1_dp * t150 * t485
          t489 = t488 * t85
          t490 = 0.2e1_dp * t472

          e_bb(ip) = e_bb(ip) + ( (t332 + t489 + t490 + t452) * t2 + t455 + 0.2e1_dp * t213 + t457 ) * sc

        END IF
        IF (order>=3.OR.order==-3) THEN
          t494 = t4 ** (-0.2666666667e1_dp)
          t497 = 0.1e1_dp / t218 / 0.3141592654e1_dp
          t499 = 0.1e1_dp / t220 / t90
          t500 = t497 * t499
          t501 = t494 * t11 * t500
          t504 = t222 * t105
          t505 = t216 * t96 * t504
          t508 = 0.1e1_dp / t220 / t2
          t509 = t219 * t508
          t510 = t217 * t509
          t513 = t92 * t237
          t516 = t231 * t105
          t519 = t92 * t251
          t522 = t1 * t221
          t525 = t95 ** 2
          t526 = 0.1e1_dp / t525
          t528 = t237 * t105
          t531 = t105 * t251
          t536 = 0.3365038134e0_dp * t494 * t497 * t499
          t538 = 0.1211413728e1_dp * t240 * t508
          t540 = 0.1817120592e1_dp * t98 * t221
          t541 = t4 ** (-0.2833333333e1_dp)
          t549 = -t536 + t538 - t540 - 0.2427089633e0_dp * bp * t541 * t500 + 0.7943202439e0_dp  &
                 * t246 * t509 - 0.9531842928e0_dp * t102 * t522
          t559 = t500 * t120
          t563 = 0.1e1_dp / t279 / t113
          t565 = t4 ** (-0.2500000000e1_dp)
          t567 = t500 * t286
          t570 = t509 * t120
          t573 = t4 ** (-0.2666666666e1_dp)
          t583 = t522 * t120
          t589 = t509 * t286
          t595 = t279 ** 2
          t596 = 0.1e1_dp / t595
          t601 = t117 ** 2
          t603 = t500 / t284 / t119 * t601
          t613 = t4 ** (-0.2333333333e1_dp)
          t625 = 0.1e1_dp / t279
          t665 = t26 ** (-0.40e1_dp)
          t690 = t541 * t497 * t499
          t693 = t295 * t508
          t696 = t126 * t221
          t705 = 0.8412595335e-1_dp * t501 - 0.1514267160e0_dp * t505 - 0.3028534320e0_dp  &
                 * t510 - 0.1906368585e1_dp * t124 * t235 * t101 * t513 + 0.7943202441e0_dp  &
                 * t291 * t245 * t504 - 0.1906368585e1_dp * t292 * t516 + 0.9531842928e0_dp  &
                 * t292 * t519 + 0.4206297667e-1_dp * t11 * t573 * t500 - 0.4854179269e0_dp  &
                 * t125 * t690 + 0.1588640488e1_dp * t125 * t693 - 0.1906368586e1_dp  &
                 * t125 * t696 - 0.6e1_dp * t27 * t526 * t528 + 0.6e1_dp * t302  &
                 * t531 - t130 * t549
          t710 = 0.6354561952e0_dp * t306 * t309 * t10 * t127 - 0.1211413727e1_dp *  &
                 t316 * t570 + 0.1924500894e0_dp * t32 * t625 * t565 * t559 + 0.3365038132e0_dp  &
                 * t315 * t494 * t559 - 0.4490502088e0_dp * t32 * t563 * t565  &
                 * t567 - 0.1588640487e1_dp * t319 * t570 + 0.1682519066e0_dp * t315 * t573  &
                 * t559 + 0.4854179267e0_dp * t136 * t541 * t559 - 0.1682519066e0_dp *  &
                 t324 * t573 * t567 + 0.1906368585e1_dp * t137 * t583 + 0.1211413727e1_dp  &
                 * t325 * t589 - 0.3365038132e0_dp * t324 * t494 * t567 + 0.2566001193e0_dp  &
                 * t32 * t596 * t565 * t603 + t134 * t251 + 0.6354561952e0_dp * t310  &
                 * t105 * t127 - 0.6354561952e0_dp * t311 * t299 + 0.1514267160e0_dp  &
                 * t132 * t665 * t10 * t241 + 0.2647734147e0_dp * t311 * t296 + t705 *  &
                 t133 * t10 + 0.2e1_dp * t307 * t105
          t719 = 0.1100642416e1_dp * (-0.3365038134e0_dp * t501 + 0.6057068641e0_dp *  &
                 t505 + 0.1211413728e1_dp * t510 - 0.1817120592e1_dp * t88 * t235 * t513  &
                 - 0.1817120592e1_dp * t225 * t516 + 0.9085602960e0_dp * t225 * t519 - 0.1817120592e1_dp  &
                 * t89 * t522 - 0.5451361779e1_dp * t5 * t526 * t528 + 0.5451361779e1_dp  &
                 * t236 * t531 - 0.9085602964e0_dp * t97 * t549) * t109 *  &
                 t10 + 0.2201284832e1_dp * t255 * t105 + 0.6730076265e0_dp * t268 * t494  &
                 * t559 - 0.8981004177e0_dp * bp * t563 * t565 * t567 - 0.3177280975e1_dp  &
                 * t273 * t570 + 0.3365038132e0_dp * t268 * t573 * t559 + 0.9708358534e0_dp  &
                 * t115 * t541 * t559 - 0.3365038132e0_dp * t282 * t573 * t567 +  &
                 0.3812737170e1_dp * t116 * t583 + 0.7337616104e0_dp * t254 * t258 * t261  &
                 + 0.2422827454e1_dp * t283 * t589 - 0.6730076265e0_dp * t282 * t494 * t567  &
                 + 0.5132002385e0_dp * bp * t596 * t565 * t603 + 0.7337616104e0_dp *  &
                 t259 * t226 + 0.1100642416e1_dp * t110 * t251 - 0.7337616104e0_dp * t259  &
                 * t260 * t230 + 0.4891744068e0_dp * t108 * t613 * t10 * t219 * t221  &
                 - t24 * t710 * t38 - 0.2422827454e1_dp * t269 * t570 + 0.3849001789e0_dp  &
                 * bp * t625 * t565 * t559
          t720 = ap * t719
          t721 = t45 ** (-0.1666666667e1_dp)
          t727 = t43 * t221
          t729 = 0.6e1_dp * t230 - 0.6e1_dp * t727
          t732 = t48 ** (-0.1666666667e1_dp)
          t743 = t349 * t202
          t745 = t154 * t451
          t748 = t500 * t179
          t752 = t500 * t407
          t755 = t522 * t179
          t758 = t509 * t407
          t769 = t68 ** (-0.40e1_dp)
          t792 = t400 ** 2
          t793 = 0.1e1_dp / t792
          t798 = t176 ** 2
          t800 = t500 / t405 / t178 * t798
          t804 = t494 * t55 * t500
          t807 = t222 * t165
          t808 = t216 * t160 * t807
          t810 = t353 * t509
          t814 = t92 * t365
          t820 = t231 * t165
          t823 = t92 * t373
          t835 = t159 ** 2
          t836 = 0.1e1_dp / t835
          t838 = t365 * t165
          t841 = t165 * t373
          t851 = -t536 + t538 - t540 - 0.2427089633e0_dp * bf * t541 * t500 + 0.7943202439e0_dp  &
                 * t368 * t509 - 0.9531842928e0_dp * t162 * t522
          t853 = 0.8412595335e-1_dp * t804 - 0.1514267160e0_dp * t808 - 0.3028534320e0_dp  &
                 * t810 - 0.1906368585e1_dp * t183 * t363 * t101 * t814 + 0.7943202441e0_dp  &
                 * t412 * t245 * t807 - 0.1906368585e1_dp * t413 * t820 + 0.9531842928e0_dp  &
                 * t413 * t823 + 0.4206297667e-1_dp * t55 * t573 * t500 - 0.4854179269e0_dp  &
                 * t184 * t690 + 0.1588640488e1_dp * t184 * t693 - 0.1906368586e1_dp  &
                 * t184 * t696 - 0.6e1_dp * t69 * t836 * t838 + 0.6e1_dp * t420  &
                 * t841 - t187 * t851
          t860 = t509 * t179
          t863 = 0.1e1_dp / t400
          t872 = 0.1e1_dp / t400 / t172
          t879 = 0.2e1_dp * t425 * t165 + t191 * t373 + 0.6354561952e0_dp * t428 *  &
                 t165 * t127 - 0.6354561952e0_dp * t429 * t299 + 0.1514267160e0_dp * t189  &
                 * t769 * t54 * t241 + 0.2647734147e0_dp * t429 * t296 + 0.1682519066e0_dp  &
                 * t433 * t573 * t748 + 0.4854179267e0_dp * t193 * t541 * t748 - 0.1682519066e0_dp  &
                 * t442 * t573 * t752 + 0.1906368585e1_dp * t194 * t755 +  &
                 0.1211413727e1_dp * t443 * t758 - 0.3365038132e0_dp * t442 * t494 * t752  &
                 + 0.2566001193e0_dp * t74 * t793 * t565 * t800 + t853 * t190 * t54  &
                 + 0.6354561952e0_dp * t424 * t427 * t54 * t127 - 0.1211413727e1_dp * t434  &
                 * t860 + 0.1924500894e0_dp * t74 * t863 * t565 * t748 + 0.3365038132e0_dp  &
                 * t433 * t494 * t748 - 0.4490502088e0_dp * t74 * t872 * t565 * t752  &
                 - 0.1588640487e1_dp * t437 * t860
          t947 = 0.9708358534e0_dp * t174 * t541 * t748 - 0.3365038132e0_dp * t403  &
                 * t573 * t752 + 0.3812737170e1_dp * t175 * t755 + 0.2422827454e1_dp * t404  &
                 * t758 - t67 * t879 * t80 - 0.6730076265e0_dp * t403 * t494 * t752  &
                 + 0.5132002385e0_dp * bf * t793 * t565 * t800 + 0.2201284832e1_dp * t377  &
                 * t165 + 0.1100642416e1_dp * (-0.3365038134e0_dp * t804 + 0.6057068641e0_dp  &
                 * t808 + 0.1211413728e1_dp * t810 - 0.1817120592e1_dp * t88 * t363 *  &
                 t814 - 0.1817120592e1_dp * t356 * t820 + 0.9085602960e0_dp * t356 * t823  &
                 - 0.1817120592e1_dp * t156 * t522 - 0.5451361779e1_dp * t5 * t836 * t838  &
                 + 0.5451361779e1_dp * t364 * t841 - 0.9085602964e0_dp * t161 * t851) *  &
                 t109 * t54 + 0.7337616104e0_dp * t376 * t258 * t382 + 0.1100642416e1_dp  &
                 * t169 * t373 + 0.7337616104e0_dp * t380 * t357 - 0.7337616104e0_dp * t380  &
                 * t381 * t230 + 0.4891744068e0_dp * t168 * t613 * t54 * t219 * t221  &
                 - 0.2422827454e1_dp * t390 * t860 + 0.3849001789e0_dp * bf * t863 * t565  &
                 * t748 + 0.6730076265e0_dp * t389 * t494 * t748 - 0.8981004177e0_dp  &
                 * bf * t872 * t565 * t752 - 0.3177280975e1_dp * t394 * t860 + 0.3365038132e0_dp  &
                 * t389 * t573 * t748
          t950 = t51 * (af * t947 - t720)
          t953 = 0.3e1_dp * t332
          t956 = 0.3e1_dp * t452

          e_aaa(ip) = e_aaa(ip) + ( (t720 + (-0.5699736440e0_dp * t721 * t334 * t147 + 0.2564881396e1_dp  &
                      * t458 * t339 + 0.2564881399e1_dp * t145 * t729 - 0.5699736440e0_dp *  &
                      t732 * t343 * t151 + 0.2564881396e1_dp * t464 * t346 - 0.2564881399e1_dp  &
                      * t150 * t729) * t85 + 0.3e1_dp * t743 + 0.3e1_dp * t745 + t950) * t2  &
                      + t953 + 0.3e1_dp * t350 + 0.6e1_dp * t351 + t956 ) * sc

          t967 = 0.2e1_dp * t230 - 0.6e1_dp * t727
          t984 = 0.2e1_dp * t470 * t202
          t986 = t212 * t451
          t990 = 0.2e1_dp * t471

          e_aab(ip) = e_aab(ip) + ( (t720 + (-0.5699736440e0_dp * t721 * t334 * t206 + 0.3419841862e1_dp  &
                      * t458 * t337 + 0.8549604655e0_dp * t333 * t339 * t206 + 0.2564881399e1_dp  &
                      * t145 * t967 - 0.5699736440e0_dp * t732 * t343 * t209 - 0.3419841862e1_dp  &
                      * t464 * t337 + 0.8549604655e0_dp * t342 * t346 * t209 - 0.2564881399e1_dp  &
                      * t150 * t967) * t85 + t743 + t984 + 0.2e1_dp * t745 + t986  &
                      + t950) * t2 + t953 + t350 + 0.4e1_dp * t351 + t956 + t990 + t490 ) * sc

          t1019 = t488 * t202

          e_abb(ip) = e_abb(ip) + ( (t720 + (-0.5699736440e0_dp * t721 * t147 * t475 + 0.3419841862e1_dp  &
                      * t333 * t43 * t230 * t206 + 0.8549604655e0_dp * t458 * t479 - 0.5129762798e1_dp  &
                      * t145 * t230 - 0.1538928839e2_dp * t461 * t221 - 0.5699736440e0_dp  &
                      * t732 * t151 * t482 - 0.3419841862e1_dp * t342 * t43 * t230  &
                      * t209 + 0.8549604655e0_dp * t464 * t485 + 0.5129762798e1_dp * t150 * t230  &
                      + 0.1538928839e2_dp * t467 * t221) * t85 + t984 + t745 + t1019 + 0.2e1_dp  &
                      * t986 + t950) * t2 + t953 + t990 + t352 + 0.4e1_dp * t472 + t956  &
                      + t489 ) * sc

          t1031 = -0.6e1_dp * t230 - 0.6e1_dp * t727

          e_bbb(ip) = e_bbb(ip) + ( (t720 + (-0.5699736440e0_dp * t721 * t475 * t206 + 0.2564881396e1_dp  &
                      * t333 * t206 * t479 + 0.2564881399e1_dp * t145 * t1031 - 0.5699736440e0_dp  &
                      * t732 * t482 * t209 + 0.2564881396e1_dp * t342 * t209 * t485  &
                      - 0.2564881399e1_dp * t150 * t1031) * t85 + 0.3e1_dp * t1019 + 0.3e1_dp * t986  &
                      + t950) * t2 + t953 + 0.3e1_dp * t489 + 0.6e1_dp * t472 + t956 ) * sc

        END IF
      END IF
    END DO

    !$omp end do

  END SUBROUTINE vwn3_lsd_calc

! *****************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param rhob ...
!> \param e_0 ...
!> \param e_a ...
!> \param e_b ...
!> \param e_aa ...
!> \param e_bb ...
!> \param e_ab ...
!> \param e_aaa ...
!> \param e_bbb ...
!> \param e_aab ...
!> \param e_abb ...
!> \param order ...
!> \param npoints ...
!> \param sc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE vwn5_lsd_calc ( rhoa, rhob, e_0, e_a, e_b, e_aa, e_bb, e_ab,&
                            e_aaa, e_bbb, e_aab, e_abb, order, npoints, sc, error )

    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: rhoa, rhob, e_0, e_a, e_b, &
                                                e_aa, e_bb, e_ab, e_aaa, &
                                                e_bbb, e_aab, e_abb
    INTEGER, INTENT(in)                      :: order, npoints
    REAL(KIND=dp)                            :: sc
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ip
    REAL(KIND=dp) :: ap, bp, cp, d2f0, myrho, myrhoa, myrhob, Qa, Qf, Qp, t1, &
      t100, t101, t1010, t1013, t1019, t1020, t1025, t1027, t1030, t104, &
      t1043, t106, t1084, t1085, t1086, t1088, t1089, t109, t1094, t1096, &
      t1099, t11, t110, t1108, t1109, t111, t1110, t1111, t1113, t1114, &
      t1116, t1118, t1119, t1120, t1125, t113, t1137, t1142, t1145, t1148, &
      t1149, t1151, t1154, t1157, t116, t1160, t1165, t1166, t1168, t1171, &
      t1181, t12, t120, t1205, t1208, t1211, t1218, t122, t1221, t1235, &
      t1238, t1244, t1245, t125, t1250, t1252, t126, t127, t1284, t129, &
      t1299, t13, t130, t131, t132, t133, t1341, t1344
    REAL(KIND=dp) :: t1346, t1358, t136, t1360, t1361, t1363, t1365, t1368, &
      t137, t1371, t1372, t1374, t1377, t1380, t1383, t1388, t1389, t1391, &
      t1394, t140, t1404, t141, t142, t143, t144, t1455, t1469, t147, t1477, &
      t148, t1480, t1483, t149, t1490, t1493, t15, t150, t1507, t151, t1510, &
      t1516, t1517, t1522, t1524, t1527, t154, t155, t156, t1567, t157, &
      t1570, t1576, t1578, t1580, t1582, t1588, t159, t1590, t1593, t1599, &
      t16, t160, t1602, t1603, t1604, t1605, t1606, t1609, t161, t1610, &
      t1611, t1613, t1615, t1620, t1622, t1626, t1639, t164, t1653, t1654, &
      t1657, t1659, t1664, t1666, t1668, t167
    REAL(KIND=dp) :: t1671, t1672, t1675, t168, t169, t1690, t1694, t1696, &
      t1697, t1698, t17, t1701, t172, t1726, t173, t1730, t174, t175, t1750, &
      t1754, t176, t1764, t1765, t1768, t1775, t1776, t1778, t178, t1782, &
      t1783, t1785, t1786, t1789, t179, t1795, t1797, t18, t180, t1804, &
      t1826, t1829, t183, t1832, t1836, t1838, t1839, t184, t185, t1850, &
      t1853, t1858, t186, t1868, t1870, t1872, t1884, t1886, t189, t1891, &
      t1893, t1894, t1897, t19, t1901, t1904, t191, t192, t193, t1939, t1945, &
      t1949, t195, t196, t1975, t1979, t1986, t1998, t1999, t2, t20, t2010, &
      t2011, t2014, t2019, t202, t2025, t203, t204
    REAL(KIND=dp) :: t2048, t205, t206, t207, t208, t211, t212, t213, t214, &
      t217, t220, t221, t224, t225, t226, t227, t228, t23, t230, t231, t232, &
      t235, t236, t239, t24, t241, t242, t243, t245, t246, t252, t253, t254, &
      t255, t256, t257, t258, t259, t260, t263, t264, t265, t266, t269, t27, &
      t272, t273, t276, t277, t278, t279, t28, t280, t282, t283, t284, t287, &
      t288, t29, t291, t293, t294, t295, t297, t298, t3, t304, t305, t306, &
      t309, t312, t315, t316, t317, t32, t320, t321, t322, t323, t324, t325, &
      t326, t327, t328, t329, t332, t333, t337, t338, t34, t340, t343, t344, &
      t347, t350, t351, t352, t353
    REAL(KIND=dp) :: t355, t356, t357, t359, t362, t363, t364, t365, t366, &
      t367, t368, t369, t37, t370, t371, t372, t373, t375, t376, t379, t38, &
      t383, t384, t385, t388, t389, t39, t390, t392, t393, t394, t399, t4, &
      t40, t402, t403, t406, t407, t408, t409, t415, t416, t417, t418, t42, &
      t421, t424, t427, t429, t430, t431, t432, t434, t435, t439, t440, t443, &
      t444, t447, t45, t450, t454, t455, t457, t458, t459, t463, t464, t467, &
      t472, t473, t479, t480, t481, t482, t483, t484, t485, t486, t487, t488, &
      t489, t49, t490, t491, t492, t493, t494, t496, t497, t5, t503, t504, &
      t505, t508, t51, t513, t516, t517
    REAL(KIND=dp) :: t520, t521, t522, t528, t529, t530, t531, t534, t537, &
      t54, t540, t542, t543, t544, t545, t547, t548, t55, t552, t553, t56, &
      t560, t564, t565, t567, t568, t569, t57, t573, t574, t577, t582, t583, &
      t591, t592, t593, t594, t595, t596, t597, t598, t599, t6, t60, t600, &
      t601, t602, t603, t604, t605, t606, t607, t608, t61, t610, t611, t617, &
      t618, t619, t622, t627, t630, t631, t634, t635, t636, t64, t642, t643, &
      t644, t645, t648, t65, t651, t654, t656, t657, t658, t659, t661, t662, &
      t666, t667, t674, t678, t679, t68, t681, t682, t683, t687, t688, t691, &
      t696, t697, t70, t704, t705, t706
    REAL(KIND=dp) :: t709, t712, t715, t716, t719, t722, t725, t728, t729, &
      t73, t730, t732, t733, t735, t74, t741, t742, t743, t744, t745, t746, &
      t748, t75, t750, t751, t752, t753, t755, t756, t757, t758, t759, t762, &
      t763, t764, t766, t767, t769, t77, t771, t772, t773, t774, t777, t778, &
      t779, t780, t782, t783, t784, t786, t789, t793, t796, t799, t8, t80, &
      t802, t803, t804, t806, t808, t811, t812, t813, t814, t815, t816, t817, &
      t818, t819, t820, t821, t822, t823, t824, t825, t826, t827, t828, t829, &
      t830, t831, t832, t833, t834, t835, t84, t842, t851, t857, t859, t86, &
      t862, t864, t865, t866, t869
    REAL(KIND=dp) :: t870, t873, t874, t875, t878, t881, t884, t887, t89, &
      t890, t891, t893, t896, t9, t90, t901, t903, t905, t906, t91, t914, &
      t92, t93, t937, t942, t945, t948, t957, t96, t97, t972, t979, t982, &
      t984, t986, t993, t996, x0p

    cp  = c
    ap  = a
    bp  = b
    x0p = x0
    Qp = SQRT(4.0_dp*cp-bp*bp)
    Qa = SQRT(4.0_dp*ca-ba*ba)
    Qf = SQRT(4.0_dp*cf-bf*bf)
    d2f0 = 4.0_dp/(9.0_dp*(2.0_dp**(1.0_dp/3.0_dp)-1.0_dp))

    !$omp do

    DO ip = 1, npoints
      myrhoa = MAX(rhoa(ip),0.0_dp)
      myrhob = MAX(rhob(ip),0.0_dp)
      myrho = myrhoa + myrhob
      IF ( myrho > eps_rho ) THEN
        IF (order>=0) THEN
          t1 = myrhoa - myrhob
          t2 = myrho
          t3 = 0.1e1_dp / t2
          t4 = t1 * t3
          t5 = 0.10e1_dp + t4
          t6 = t5 ** 0.1333333333e1_dp
          t8 = 0.10e1_dp - t4
          t9 = t8 ** 0.1333333333e1_dp
          t11 = 0.1923661050e1_dp * t6 + 0.1923661050e1_dp * t9 - 0.3847322100e1_dp
          t12 = t4 ** 0.40e1_dp
          t13 = t11 * t12
          t15 = (0.10e1_dp - t13) * ap
          t16 = 0.1e1_dp / 0.3141592654e1_dp
          t17 = t16 * t3
          t18 = t17 ** 0.3333333332e0_dp
          t19 = 0.9085602964e0_dp * t18
          t20 = t17 ** 0.1666666666e0_dp
          t23 = t19 + 0.9531842930e0_dp * bp * t20 + cp
          t24 = 0.1e1_dp / t23
          t27 = LOG(0.9085602964e0_dp * t18 * t24)
          t28 = 0.1906368586e1_dp * t20
          t29 = t28 + bp
          t32 = ATAN(Qp / t29)
          t34 = 0.1e1_dp / Qp
          t37 = bp * x0p
          t38 = 0.9531842930e0_dp * t20
          t39 = t38 - x0p
          t40 = t39 ** 0.20e1_dp
          t42 = LOG(t40 * t24)
          t45 = 0.20e1_dp * bp + 0.40e1_dp * x0p
          t49 = x0p ** 0.20e1_dp
          t51 = 0.1e1_dp / (t49 + t37 + cp)
          t54 = t27 + 0.20e1_dp * bp * t32 * t34 - t37 * (t42 + t45 * t32 * t34)  &
                * t51
          t55 = t15 * t54
          t56 = 0.10e1_dp - t12
          t57 = t11 * t56
          t60 = t19 + 0.9531842930e0_dp * ba * t20 + ca
          t61 = 0.1e1_dp / t60
          t64 = LOG(0.9085602964e0_dp * t18 * t61)
          t65 = t28 + ba
          t68 = ATAN(Qa / t65)
          t70 = 0.1e1_dp / Qa
          t73 = ba * x0a
          t74 = t38 - x0a
          t75 = t74 ** 0.20e1_dp
          t77 = LOG(t75 * t61)
          t80 = 0.20e1_dp * ba + 0.40e1_dp * x0a
          t84 = x0a ** 0.20e1_dp
          t86 = 0.1e1_dp / (t84 + t73 + ca)
          t89 = t64 + 0.20e1_dp * ba * t68 * t70 - t73 * (t77 + t80 * t68 * t70)  &
                * t86
          t90 = aa * t89
          t91 = 0.1e1_dp / d2f0
          t92 = t90 * t91
          t93 = t57 * t92
          t96 = t19 + 0.9531842930e0_dp * bf * t20 + cf
          t97 = 0.1e1_dp / t96
          t100 = LOG(0.9085602964e0_dp * t18 * t97)
          t101 = t28 + bf
          t104 = ATAN(Qf / t101)
          t106 = 0.1e1_dp / Qf
          t109 = bf * x0f
          t110 = t38 - x0f
          t111 = t110 ** 0.20e1_dp
          t113 = LOG(t111 * t97)
          t116 = 0.20e1_dp * bf + 0.40e1_dp * x0f
          t120 = x0f ** 0.20e1_dp
          t122 = 0.1e1_dp / (t120 + t109 + cf)
          t125 = t100 + 0.20e1_dp * bf * t104 * t106 - t109 * (t113 + t116 * t104  &
                 * t106) * t122
          t126 = af * t125
          t127 = t13 * t126

          e_0(ip) = e_0(ip) + ( (t55 + t93 + t127) * t2 ) * sc

        END IF
        IF (order>=1.OR.order==-1) THEN
          t129 = t5 ** 0.333333333e0_dp
          t130 = t2 ** 2
          t131 = 0.1e1_dp / t130
          t132 = t1 * t131
          t133 = t3 - t132
          t136 = t8 ** 0.333333333e0_dp
          t137 = -t133
          t140 = 0.2564881399e1_dp * t129 * t133 + 0.2564881399e1_dp * t136 * t137
          t141 = t140 * t12
          t142 = t4 ** 0.30e1_dp
          t143 = t11 * t142
          t144 = t143 * t133
          t147 = (-t141 - 0.40e1_dp * t144) * ap
          t148 = t147 * t54
          t149 = t17 ** (-0.6666666668e0_dp)
          t150 = t149 * t24
          t151 = t16 * t131
          t154 = t23 ** 2
          t155 = 0.1e1_dp / t154
          t156 = t18 * t155
          t157 = t149 * t16
          t159 = 0.3028534320e0_dp * t157 * t131
          t160 = t17 ** (-0.8333333334e0_dp)
          t161 = bp * t160
          t164 = -t159 - 0.1588640488e0_dp * t161 * t151
          t167 = -0.3028534320e0_dp * t150 * t151 - 0.9085602964e0_dp * t156 * t164
          t168 = t17 ** (-0.3333333332e0_dp)
          t169 = t167 * t168
          t172 = t29 ** 2
          t173 = 0.1e1_dp / t172
          t174 = bp * t173
          t175 = t174 * t160
          t176 = Qp ** 2
          t178 = 0.1e1_dp + t176 * t173
          t179 = 0.1e1_dp / t178
          t180 = t151 * t179
          t183 = t39 ** 0.10e1_dp
          t184 = t183 * t24
          t185 = t160 * t16
          t186 = t185 * t131
          t189 = t40 * t155
          t191 = -0.3177280976e0_dp * t184 * t186 - t189 * t164
          t192 = t39 ** (-0.20e1_dp)
          t193 = t191 * t192
          t195 = t45 * t173
          t196 = t195 * t160
          t202 = 0.1100642416e1_dp * t169 * t23 + 0.6354561950e0_dp * t175 * t180 -  &
                 t37 * (t193 * t23 + 0.3177280975e0_dp * t196 * t180) * t51
          t203 = t15 * t202
          t204 = t140 * t56
          t205 = t204 * t92
          t206 = t144 * t92
          t207 = 0.40e1_dp * t206
          t208 = t149 * t61
          t211 = t60 ** 2
          t212 = 0.1e1_dp / t211
          t213 = t18 * t212
          t214 = ba * t160
          t217 = -t159 - 0.1588640488e0_dp * t214 * t151
          t220 = -0.3028534320e0_dp * t208 * t151 - 0.9085602964e0_dp * t213 * t217
          t221 = t220 * t168
          t224 = t65 ** 2
          t225 = 0.1e1_dp / t224
          t226 = ba * t225
          t227 = t226 * t160
          t228 = Qa ** 2
          t230 = 0.1e1_dp + t228 * t225
          t231 = 0.1e1_dp / t230
          t232 = t151 * t231
          t235 = t74 ** 0.10e1_dp
          t236 = t235 * t61
          t239 = t75 * t212
          t241 = -0.3177280976e0_dp * t236 * t186 - t239 * t217
          t242 = t74 ** (-0.20e1_dp)
          t243 = t241 * t242
          t245 = t80 * t225
          t246 = t245 * t160
          t252 = 0.1100642416e1_dp * t221 * t60 + 0.6354561950e0_dp * t227 * t232 -  &
                 t73 * (t243 * t60 + 0.3177280975e0_dp * t246 * t232) * t86
          t253 = aa * t252
          t254 = t253 * t91
          t255 = t57 * t254
          t256 = t141 * t126
          t257 = t126 * t133
          t258 = t143 * t257
          t259 = 0.40e1_dp * t258
          t260 = t149 * t97
          t263 = t96 ** 2
          t264 = 0.1e1_dp / t263
          t265 = t18 * t264
          t266 = bf * t160
          t269 = -t159 - 0.1588640488e0_dp * t266 * t151
          t272 = -0.3028534320e0_dp * t260 * t151 - 0.9085602964e0_dp * t265 * t269
          t273 = t272 * t168
          t276 = t101 ** 2
          t277 = 0.1e1_dp / t276
          t278 = bf * t277
          t279 = t278 * t160
          t280 = Qf ** 2
          t282 = 0.1e1_dp + t280 * t277
          t283 = 0.1e1_dp / t282
          t284 = t151 * t283
          t287 = t110 ** 0.10e1_dp
          t288 = t287 * t97
          t291 = t111 * t264
          t293 = -0.3177280976e0_dp * t288 * t186 - t291 * t269
          t294 = t110 ** (-0.20e1_dp)
          t295 = t293 * t294
          t297 = t116 * t277
          t298 = t297 * t160
          t304 = 0.1100642416e1_dp * t273 * t96 + 0.6354561950e0_dp * t279 * t284 -  &
                 t109 * (t295 * t96 + 0.3177280975e0_dp * t298 * t284) * t122
          t305 = af * t304
          t306 = t13 * t305

          e_a(ip) = e_a(ip) + ( (t148 + t203 + t205 - t207 + t255 + t256 + t259 + t306) * t2 +  &
                    t55 + t93 + t127 ) * sc

          t309 = -t3 - t132
          t312 = -t309
          t315 = 0.2564881399e1_dp * t129 * t309 + 0.2564881399e1_dp * t136 * t312
          t316 = t315 * t12
          t317 = t143 * t309
          t320 = (-t316 - 0.40e1_dp * t317) * ap
          t321 = t320 * t54
          t322 = t315 * t56
          t323 = t322 * t92
          t324 = t317 * t92
          t325 = 0.40e1_dp * t324
          t326 = t316 * t126
          t327 = t126 * t309
          t328 = t143 * t327
          t329 = 0.40e1_dp * t328

          e_b(ip) = e_b(ip) + ( (t321 + t203 + t323 - t325 + t255 + t326 + t329 + t306) * t2 +  &
                   t55 + t93 + t127 ) * sc

        END IF
        IF (order>=2.OR.order==-2) THEN
          t332 = t5 ** (-0.666666667e0_dp)
          t333 = t133 ** 2
          t337 = 0.1e1_dp / t130 / t2
          t338 = t1 * t337
          t340 = -0.2e1_dp * t131 + 0.2e1_dp * t338
          t343 = t8 ** (-0.666666667e0_dp)
          t344 = t137 ** 2
          t347 = -t340
          t350 = 0.8549604655e0_dp * t332 * t333 + 0.2564881399e1_dp * t129 * t340  &
                 + 0.8549604655e0_dp * t343 * t344 + 0.2564881399e1_dp * t136 * t347
          t351 = t350 * t12
          t352 = t140 * t142
          t353 = t352 * t133
          t355 = t4 ** 0.20e1_dp
          t356 = t11 * t355
          t357 = t356 * t333
          t359 = t143 * t340
          t362 = (-t351 - 0.80e1_dp * t353 - 0.1200e2_dp * t357 - 0.40e1_dp * t359) *  &
                 ap
          t363 = t362 * t54
          t364 = t147 * t202
          t365 = 0.2e1_dp * t364
          t366 = t17 ** (-0.1666666667e1_dp)
          t367 = t366 * t24
          t368 = 0.3141592654e1_dp ** 2
          t369 = 0.1e1_dp / t368
          t370 = t130 ** 2
          t371 = 0.1e1_dp / t370
          t372 = t369 * t371
          t373 = t367 * t372
          t375 = t149 * t155
          t376 = t151 * t164
          t379 = t16 * t337
          t383 = 0.1e1_dp / t154 / t23
          t384 = t18 * t383
          t385 = t164 ** 2
          t388 = t366 * t369
          t389 = t388 * t371
          t390 = 0.2019022880e0_dp * t389
          t392 = 0.6057068640e0_dp * t157 * t337
          t393 = t17 ** (-0.1833333333e1_dp)
          t394 = bp * t393
          t399 = -t390 + t392 - 0.1323867073e0_dp * t394 * t372 + 0.3177280976e0_dp  &
                 * t161 * t379
          t402 = -0.2019022880e0_dp * t373 + 0.6057068640e0_dp * t375 * t376 + 0.6057068640e0_dp  &
                 * t150 * t379 + 0.1817120593e1_dp * t384 * t385 - 0.9085602964e0_dp  &
                 * t156 * t399
          t403 = t402 * t168
          t406 = t17 ** (-0.1333333333e1_dp)
          t407 = t167 * t406
          t408 = t23 * t16
          t409 = t408 * t131
          t415 = 0.1e1_dp / t172 / t29
          t416 = bp * t415
          t417 = t416 * t366
          t418 = t372 * t179
          t421 = t174 * t393
          t424 = t379 * t179
          t427 = t172 ** 2
          t429 = 0.1e1_dp / t427 / t29
          t430 = bp * t429
          t431 = t430 * t366
          t432 = t178 ** 2
          t434 = 0.1e1_dp / t432 * t176
          t435 = t372 * t434
          t439 = t183 * t155
          t440 = t439 * t160
          t443 = t393 * t369
          t444 = t443 * t371
          t447 = t185 * t337
          t450 = t40 * t383
          t454 = 0.5047557200e-1_dp * t373 + 0.6354561952e0_dp * t440 * t376 - 0.2647734147e0_dp  &
                 * t184 * t444 + 0.6354561952e0_dp * t184 * t447 + 0.2e1_dp *  &
                 t450 * t385 - t189 * t399
          t455 = t454 * t192
          t457 = t39 ** (-0.30e1_dp)
          t458 = t191 * t457
          t459 = t458 * t23
          t463 = t45 * t415
          t464 = t463 * t366
          t467 = t195 * t393
          t472 = t45 * t429
          t473 = t472 * t366
          t479 = 0.1100642416e1_dp * t403 * t23 + 0.3668808052e0_dp * t407 * t409 +  &
                 0.1100642416e1_dp * t169 * t164 + 0.4038045758e0_dp * t417 * t418 + 0.5295468292e0_dp  &
                 * t421 * t418 - 0.1270912390e1_dp * t175 * t424 - 0.4038045758e0_dp  &
                 * t431 * t435 - t37 * (t455 * t23 + 0.3177280976e0_dp * t459  &
                 * t186 + t193 * t164 + 0.2019022879e0_dp * t464 * t418 + 0.2647734146e0_dp  &
                 * t467 * t418 - 0.6354561950e0_dp * t196 * t424 - 0.2019022879e0_dp *  &
                 t473 * t435) * t51
          t480 = t15 * t479
          t481 = t350 * t56
          t482 = t481 * t92
          t483 = t353 * t92
          t484 = 0.80e1_dp * t483
          t485 = t204 * t254
          t486 = 0.2e1_dp * t485
          t487 = t357 * t92
          t488 = 0.1200e2_dp * t487
          t489 = t359 * t92
          t490 = 0.40e1_dp * t489
          t491 = t144 * t254
          t492 = 0.80e1_dp * t491
          t493 = t366 * t61
          t494 = t493 * t372
          t496 = t149 * t212
          t497 = t151 * t217
          t503 = 0.1e1_dp / t211 / t60
          t504 = t18 * t503
          t505 = t217 ** 2
          t508 = ba * t393
          t513 = -t390 + t392 - 0.1323867073e0_dp * t508 * t372 + 0.3177280976e0_dp  &
                 * t214 * t379
          t516 = -0.2019022880e0_dp * t494 + 0.6057068640e0_dp * t496 * t497 + 0.6057068640e0_dp  &
                 * t208 * t379 + 0.1817120593e1_dp * t504 * t505 - 0.9085602964e0_dp  &
                 * t213 * t513
          t517 = t516 * t168
          t520 = t220 * t406
          t521 = t60 * t16
          t522 = t521 * t131
          t528 = 0.1e1_dp / t224 / t65
          t529 = ba * t528
          t530 = t529 * t366
          t531 = t372 * t231
          t534 = t226 * t393
          t537 = t379 * t231
          t540 = t224 ** 2
          t542 = 0.1e1_dp / t540 / t65
          t543 = ba * t542
          t544 = t543 * t366
          t545 = t230 ** 2
          t547 = 0.1e1_dp / t545 * t228
          t548 = t372 * t547
          t552 = t235 * t212
          t553 = t552 * t160
          t560 = t75 * t503
          t564 = 0.5047557200e-1_dp * t494 + 0.6354561952e0_dp * t553 * t497 - 0.2647734147e0_dp  &
                 * t236 * t444 + 0.6354561952e0_dp * t236 * t447 + 0.2e1_dp *  &
                 t560 * t505 - t239 * t513
          t565 = t564 * t242
          t567 = t74 ** (-0.30e1_dp)
          t568 = t241 * t567
          t569 = t568 * t60
          t573 = t80 * t528
          t574 = t573 * t366
          t577 = t245 * t393
          t582 = t80 * t542
          t583 = t582 * t366
          t591 = aa * (0.1100642416e1_dp * t517 * t60 + 0.3668808052e0_dp * t520 *  &
                 t522 + 0.1100642416e1_dp * t221 * t217 + 0.4038045758e0_dp * t530 * t531  &
                 + 0.5295468292e0_dp * t534 * t531 - 0.1270912390e1_dp * t227 * t537 - 0.4038045758e0_dp  &
                 * t544 * t548 - t73 * (t565 * t60 + 0.3177280976e0_dp *  &
                 t569 * t186 + t243 * t217 + 0.2019022879e0_dp * t574 * t531 + 0.2647734146e0_dp  &
                 * t577 * t531 - 0.6354561950e0_dp * t246 * t537 - 0.2019022879e0_dp  &
                 * t583 * t548) * t86) * t91
          t592 = t57 * t591
          t593 = t351 * t126
          t594 = t352 * t257
          t595 = 0.80e1_dp * t594
          t596 = t141 * t305
          t597 = 0.2e1_dp * t596
          t598 = t126 * t333
          t599 = t356 * t598
          t600 = 0.1200e2_dp * t599
          t601 = t305 * t133
          t602 = t143 * t601
          t603 = 0.80e1_dp * t602
          t604 = t126 * t340
          t605 = t143 * t604
          t606 = 0.40e1_dp * t605
          t607 = t366 * t97
          t608 = t607 * t372
          t610 = t149 * t264
          t611 = t151 * t269
          t617 = 0.1e1_dp / t263 / t96
          t618 = t18 * t617
          t619 = t269 ** 2
          t622 = bf * t393
          t627 = -t390 + t392 - 0.1323867073e0_dp * t622 * t372 + 0.3177280976e0_dp  &
                 * t266 * t379
          t630 = -0.2019022880e0_dp * t608 + 0.6057068640e0_dp * t610 * t611 + 0.6057068640e0_dp  &
                 * t260 * t379 + 0.1817120593e1_dp * t618 * t619 - 0.9085602964e0_dp  &
                 * t265 * t627
          t631 = t630 * t168
          t634 = t272 * t406
          t635 = t96 * t16
          t636 = t635 * t131
          t642 = 0.1e1_dp / t276 / t101
          t643 = bf * t642
          t644 = t643 * t366
          t645 = t372 * t283
          t648 = t278 * t393
          t651 = t379 * t283
          t654 = t276 ** 2
          t656 = 0.1e1_dp / t654 / t101
          t657 = bf * t656
          t658 = t657 * t366
          t659 = t282 ** 2
          t661 = 0.1e1_dp / t659 * t280
          t662 = t372 * t661
          t666 = t287 * t264
          t667 = t666 * t160
          t674 = t111 * t617
          t678 = 0.5047557200e-1_dp * t608 + 0.6354561952e0_dp * t667 * t611 - 0.2647734147e0_dp  &
                 * t288 * t444 + 0.6354561952e0_dp * t288 * t447 + 0.2e1_dp *  &
                 t674 * t619 - t291 * t627
          t679 = t678 * t294
          t681 = t110 ** (-0.30e1_dp)
          t682 = t293 * t681
          t683 = t682 * t96
          t687 = t116 * t642
          t688 = t687 * t366
          t691 = t297 * t393
          t696 = t116 * t656
          t697 = t696 * t366
          t704 = af * (0.1100642416e1_dp * t631 * t96 + 0.3668808052e0_dp * t634 *  &
                 t636 + 0.1100642416e1_dp * t273 * t269 + 0.4038045758e0_dp * t644 * t645  &
                 + 0.5295468292e0_dp * t648 * t645 - 0.1270912390e1_dp * t279 * t651 - 0.4038045758e0_dp  &
                 * t658 * t662 - t109 * (t679 * t96 + 0.3177280976e0_dp  &
                 * t683 * t186 + t295 * t269 + 0.2019022879e0_dp * t688 * t645 + 0.2647734146e0_dp  &
                 * t691 * t645 - 0.6354561950e0_dp * t298 * t651 - 0.2019022879e0_dp  &
                 * t697 * t662) * t122)
          t705 = t13 * t704
          t706 = t363 + t365 + t480 + t482 - t484 + t486 - t488 - t490 - t492  &
                 + t592 + t593 + t595 + t597 + t600 + t603 + t606 + t705
          t709 = 0.2e1_dp * t203
          t712 = 0.2e1_dp * t255
          t715 = 0.2e1_dp * t306

          e_aa(ip) = e_aa(ip) + ( t706 * t2 + 0.2e1_dp * t148 + t709 + 0.2e1_dp * t205 - 0.80e1_dp * t206  &
                     + t712 + 0.2e1_dp * t256 + 0.80e1_dp * t258 + t715 ) * sc

          t716 = t332 * t133
          t719 = t129 * t1
          t722 = t343 * t137
          t725 = t136 * t1
          t728 = 0.8549604655e0_dp * t716 * t309 + 0.5129762798e1_dp * t719 * t337  &
                 + 0.8549604655e0_dp * t722 * t312 - 0.5129762798e1_dp * t725 * t337
          t729 = t728 * t12
          t730 = t352 * t309
          t732 = t315 * t142
          t733 = t732 * t133
          t735 = t133 * t309
          t741 = (-t729 - 0.40e1_dp * t730 - 0.40e1_dp * t733 - 0.1200e2_dp * t356 * t735  &
                 - 0.80e1_dp * t143 * t338) * ap
          t742 = t741 * t54
          t743 = t320 * t202
          t744 = t728 * t56
          t745 = t744 * t92
          t746 = t730 * t92
          t748 = t733 * t92
          t750 = t356 * t133
          t751 = t91 * t309
          t752 = t90 * t751
          t753 = t750 * t752
          t755 = t143 * t1
          t756 = t337 * aa
          t757 = t89 * t91
          t758 = t756 * t757
          t759 = t755 * t758
          t762 = t322 * t254
          t763 = t742 + t364 + t743 + t480 + t745 - 0.40e1_dp * t746 + t485 - 0.40e1_dp  &
                 * t748 - 0.1200e2_dp * t753 - 0.80e1_dp * t759 - 0.40e1_dp * t491 + t762
          t764 = t317 * t254
          t766 = t729 * t126
          t767 = t352 * t327
          t769 = t732 * t257
          t771 = t356 * af
          t772 = t125 * t133
          t773 = t772 * t309
          t774 = t771 * t773
          t777 = t143 * af
          t778 = t125 * t1
          t779 = t778 * t337
          t780 = t777 * t779
          t782 = t316 * t305
          t783 = t305 * t309
          t784 = t143 * t783
          t786 = -0.40e1_dp * t764 + t592 + t766 + 0.40e1_dp * t767 + t596 + 0.40e1_dp  &
                 * t769 + 0.1200e2_dp * t774 + 0.40e1_dp * t602 + 0.80e1_dp * t780 + t782 +  &
                 0.40e1_dp * t784 + t705

          e_ab(ip) = e_ab(ip) + ( (t763 + t786) * t2 + t148 + t709 + t205 - t207 + t712 + t256  &
                     + t259 + t715 + t321 + t323 - t325 + t326 + t329 ) * sc
          t789 = t309 ** 2
          t793 = 0.2e1_dp * t131 + 0.2e1_dp * t338
          t796 = t312 ** 2
          t799 = -t793
          t802 = 0.8549604655e0_dp * t332 * t789 + 0.2564881399e1_dp * t129 * t793  &
                 + 0.8549604655e0_dp * t343 * t796 + 0.2564881399e1_dp * t136 * t799
          t803 = t802 * t12
          t804 = t732 * t309
          t806 = t356 * t789
          t808 = t143 * t793
          t811 = (-t803 - 0.80e1_dp * t804 - 0.1200e2_dp * t806 - 0.40e1_dp * t808) *  &
                 ap
          t812 = t811 * t54
          t813 = 0.2e1_dp * t743
          t814 = t802 * t56
          t815 = t814 * t92
          t816 = t804 * t92
          t817 = 0.80e1_dp * t816
          t818 = 0.2e1_dp * t762
          t819 = t806 * t92
          t820 = 0.1200e2_dp * t819
          t821 = t808 * t92
          t822 = 0.40e1_dp * t821
          t823 = 0.80e1_dp * t764
          t824 = t803 * t126
          t825 = t732 * t327
          t826 = 0.80e1_dp * t825
          t827 = 0.2e1_dp * t782
          t828 = t126 * t789
          t829 = t356 * t828
          t830 = 0.1200e2_dp * t829
          t831 = 0.80e1_dp * t784
          t832 = t126 * t793
          t833 = t143 * t832
          t834 = 0.40e1_dp * t833
          t835 = t812 + t813 + t480 + t815 - t817 + t818 - t820 - t822 - t823  &
                 + t592 + t824 + t826 + t827 + t830 + t831 + t834 + t705

          e_bb(ip) = e_bb(ip) + ( t835 * t2 + 0.2e1_dp * t321 + t709 + 0.2e1_dp * t323 - 0.80e1_dp * t324  &
                     + t712 + 0.2e1_dp * t326 + 0.80e1_dp * t328 + t715 ) * sc

        END IF
        IF (order>=3.OR.order==-3) THEN
          t842 = 0.3e1_dp * t480
          t851 = 0.3e1_dp * t592
          t857 = 0.3e1_dp * t705
          t859 = t17 ** (-0.2666666667e1_dp)
          t862 = 0.1e1_dp / t368 / 0.3141592654e1_dp
          t864 = 0.1e1_dp / t370 / t130
          t865 = t862 * t864
          t866 = t859 * t24 * t865
          t869 = t372 * t164
          t870 = t366 * t155 * t869
          t873 = 0.1e1_dp / t370 / t2
          t874 = t369 * t873
          t875 = t367 * t874
          t878 = t151 * t385
          t881 = t379 * t164
          t884 = t151 * t399
          t887 = t16 * t371
          t890 = t154 ** 2
          t891 = 0.1e1_dp / t890
          t893 = t385 * t164
          t896 = t164 * t399
          t901 = 0.3365038134e0_dp * t859 * t862 * t864
          t903 = 0.1211413728e1_dp * t388 * t873
          t905 = 0.1817120592e1_dp * t157 * t371
          t906 = t17 ** (-0.2833333333e1_dp)
          t914 = -t901 + t903 - t905 - 0.2427089633e0_dp * bp * t906 * t865 + 0.7943202439e0_dp  &
                 * t394 * t874 - 0.9531842928e0_dp * t161 * t887
          t937 = t17 ** (-0.2666666666e1_dp)
          t942 = t906 * t862 * t864
          t945 = t443 * t873
          t948 = t185 * t371
          t957 = 0.8412595335e-1_dp * t866 - 0.1514267160e0_dp * t870 - 0.3028534320e0_dp  &
                 * t875 - 0.1906368585e1_dp * t183 * t383 * t160 * t878 + 0.7943202441e0_dp  &
                 * t439 * t393 * t869 - 0.1906368585e1_dp * t440 * t881 + 0.9531842928e0_dp  &
                 * t440 * t884 + 0.4206297667e-1_dp * t24 * t937 * t865 - 0.4854179269e0_dp  &
                 * t184 * t942 + 0.1588640488e1_dp * t184 * t945 - 0.1906368586e1_dp  &
                 * t184 * t948 - 0.6e1_dp * t40 * t891 * t893 + 0.6e1_dp * t450  &
                 * t896 - t189 * t914
          t972 = t39 ** (-0.40e1_dp)
          t979 = t874 * t179
          t982 = 0.1e1_dp / t427
          t984 = t17 ** (-0.2500000000e1_dp)
          t986 = t865 * t179
          t993 = 0.1e1_dp / t427 / t172
          t996 = t865 * t434
          t1010 = t887 * t179
          t1013 = t874 * t434
          t1019 = t427 ** 2
          t1020 = 0.1e1_dp / t1019
          t1025 = t176 ** 2
          t1027 = t865 / t432 / t178 * t1025
          t1030 = t957 * t192 * t23 + 0.2e1_dp * t455 * t164 + 0.6354561952e0_dp *  &
                  t454 * t457 * t23 * t186 + t193 * t399 + 0.6354561952e0_dp * t458 * t164  &
                  * t186 - 0.6354561952e0_dp * t459 * t447 + 0.1514267160e0_dp * t191  &
                  * t972 * t23 * t389 + 0.2647734147e0_dp * t459 * t444 - 0.1211413727e1_dp  &
                  * t464 * t979 + 0.1924500894e0_dp * t45 * t982 * t984 * t986 + 0.3365038132e0_dp  &
                  * t463 * t859 * t986 - 0.4490502088e0_dp * t45 * t993 * t984  &
                  * t996 - 0.1588640487e1_dp * t467 * t979 + 0.1682519066e0_dp * t463 *  &
                  t937 * t986 + 0.4854179267e0_dp * t195 * t906 * t986 - 0.1682519066e0_dp  &
                  * t472 * t937 * t996 + 0.1906368585e1_dp * t196 * t1010 + 0.1211413727e1_dp  &
                  * t473 * t1013 - 0.3365038132e0_dp * t472 * t859 * t996 + 0.2566001193e0_dp  &
                  * t45 * t1020 * t984 * t1027
          t1043 = t17 ** (-0.2333333333e1_dp)
          t1084 = 0.1100642416e1_dp * (-0.3365038134e0_dp * t866 + 0.6057068641e0_dp *  &
                  t870 + 0.1211413728e1_dp * t875 - 0.1817120592e1_dp * t149 * t383 * t878  &
                  - 0.1817120592e1_dp * t375 * t881 + 0.9085602960e0_dp * t375 * t884 -  &
                  0.1817120592e1_dp * t150 * t887 - 0.5451361779e1_dp * t18 * t891 * t893  &
                  + 0.5451361779e1_dp * t384 * t896 - 0.9085602964e0_dp * t156 * t914) * t168  &
                  * t23 + 0.2201284832e1_dp * t403 * t164 - t37 * t1030 * t51 + 0.7337616104e0_dp  &
                  * t402 * t406 * t409 + 0.1100642416e1_dp * t169 * t399 +  &
                  0.7337616104e0_dp * t407 * t376 - 0.7337616104e0_dp * t407 * t408 * t337  &
                  + 0.4891744068e0_dp * t167 * t1043 * t23 * t369 * t371 - 0.2422827454e1_dp  &
                  * t417 * t979 + 0.3849001789e0_dp * bp * t982 * t984 * t986 + 0.6730076265e0_dp  &
                  * t416 * t859 * t986 - 0.8981004177e0_dp * bp * t993 * t984  &
                  * t996 - 0.3177280975e1_dp * t421 * t979 + 0.3365038132e0_dp * t416 *  &
                  t937 * t986 + 0.9708358534e0_dp * t174 * t906 * t986 - 0.3365038132e0_dp  &
                  * t430 * t937 * t996 + 0.3812737170e1_dp * t175 * t1010 + 0.2422827454e1_dp  &
                  * t431 * t1013 - 0.6730076265e0_dp * t430 * t859 * t996 + 0.5132002385e0_dp  &
                  * bp * t1020 * t984 * t1027
          t1085 = t15 * t1084
          t1086 = t147 * t479
          t1088 = t5 ** (-0.1666666667e1_dp)
          t1089 = t333 * t133
          t1094 = t1 * t371
          t1096 = 0.6e1_dp * t337 - 0.6e1_dp * t1094
          t1099 = t8 ** (-0.1666666667e1_dp)
          t1108 = -0.5699736440e0_dp * t1088 * t1089 + 0.2564881396e1_dp * t716 * t340  &
                  + 0.2564881399e1_dp * t129 * t1096 - 0.5699736440e0_dp * t1099 * t344  &
                  * t137 + 0.2564881396e1_dp * t722 * t347 - 0.2564881399e1_dp * t136 *  &
                  t1096
          t1109 = t1108 * t12
          t1110 = t350 * t142
          t1111 = t1110 * t133
          t1113 = t140 * t355
          t1114 = t1113 * t333
          t1116 = t352 * t340
          t1118 = t4 ** 0.10e1_dp
          t1119 = t11 * t1118
          t1120 = t1119 * t1089
          t1125 = t143 * t1096
          t1137 = t351 * t305
          t1142 = t352 * t601
          t1145 = t859 * t97 * t865
          t1148 = t372 * t269
          t1149 = t366 * t264 * t1148
          t1151 = t607 * t874
          t1154 = t151 * t619
          t1157 = t379 * t269
          t1160 = t151 * t627
          t1165 = t263 ** 2
          t1166 = 0.1e1_dp / t1165
          t1168 = t619 * t269
          t1171 = t269 * t627
          t1181 = -t901 + t903 - t905 - 0.2427089633e0_dp * bf * t906 * t865 + 0.7943202439e0_dp  &
                  * t622 * t874 - 0.9531842928e0_dp * t266 * t887
          t1205 = t874 * t283
          t1208 = 0.1e1_dp / t654
          t1211 = t865 * t283
          t1218 = 0.1e1_dp / t654 / t276
          t1221 = t865 * t661
          t1235 = t887 * t283
          t1238 = t874 * t661
          t1244 = t654 ** 2
          t1245 = 0.1e1_dp / t1244
          t1250 = t280 ** 2
          t1252 = t865 / t659 / t282 * t1250
          t1284 = 0.8412595335e-1_dp * t1145 - 0.1514267160e0_dp * t1149 - 0.3028534320e0_dp  &
                  * t1151 - 0.1906368585e1_dp * t287 * t617 * t160 * t1154 + 0.7943202441e0_dp  &
                  * t666 * t393 * t1148 - 0.1906368585e1_dp * t667 * t1157  &
                  + 0.9531842928e0_dp * t667 * t1160 + 0.4206297667e-1_dp * t97 * t937 * t865  &
                  - 0.4854179269e0_dp * t288 * t942 + 0.1588640488e1_dp * t288 * t945  &
                  - 0.1906368586e1_dp * t288 * t948 - 0.6e1_dp * t111 * t1166 * t1168 + 0.6e1_dp  &
                  * t674 * t1171 - t291 * t1181
          t1299 = t110 ** (-0.40e1_dp)
          t1341 = t1284 * t294 * t96 + 0.2e1_dp * t679 * t269 + 0.6354561952e0_dp *  &
                  t678 * t681 * t96 * t186 + t295 * t627 + 0.6354561952e0_dp * t682 *  &
                  t269 * t186 - 0.6354561952e0_dp * t683 * t447 + 0.1514267160e0_dp * t293  &
                  * t1299 * t96 * t389 + 0.2647734147e0_dp * t683 * t444 - 0.1211413727e1_dp  &
                  * t688 * t1205 + 0.1924500894e0_dp * t116 * t1208 * t984 * t1211  &
                  + 0.3365038132e0_dp * t687 * t859 * t1211 - 0.4490502088e0_dp * t116 * t1218  &
                  * t984 * t1221 - 0.1588640487e1_dp * t691 * t1205 + 0.1682519066e0_dp  &
                  * t687 * t937 * t1211 + 0.4854179267e0_dp * t297 * t906 * t1211 -  &
                  0.1682519066e0_dp * t696 * t937 * t1221 + 0.1906368585e1_dp * t298 * t1235  &
                  + 0.1211413727e1_dp * t697 * t1238 - 0.3365038132e0_dp * t696 * t859  &
                  * t1221 + 0.2566001193e0_dp * t116 * t1245 * t984 * t1252
          t1344 = 0.1100642416e1_dp * (-0.3365038134e0_dp * t1145 + 0.6057068641e0_dp  &
                  * t1149 + 0.1211413728e1_dp * t1151 - 0.1817120592e1_dp * t149 * t617 *  &
                  t1154 - 0.1817120592e1_dp * t610 * t1157 + 0.9085602960e0_dp * t610 * t1160  &
                  - 0.1817120592e1_dp * t260 * t887 - 0.5451361779e1_dp * t18 * t1166  &
                  * t1168 + 0.5451361779e1_dp * t618 * t1171 - 0.9085602964e0_dp * t265 *  &
                  t1181) * t168 * t96 + 0.2201284832e1_dp * t631 * t269 + 0.7337616104e0_dp  &
                  * t630 * t406 * t636 + 0.1100642416e1_dp * t273 * t627 + 0.7337616104e0_dp  &
                  * t634 * t611 - 0.7337616104e0_dp * t634 * t635 * t337 + 0.4891744068e0_dp  &
                  * t272 * t1043 * t96 * t369 * t371 - 0.2422827454e1_dp * t644  &
                  * t1205 + 0.3849001789e0_dp * bf * t1208 * t984 * t1211 + 0.6730076265e0_dp  &
                  * t643 * t859 * t1211 - 0.8981004177e0_dp * bf * t1218 * t984 *  &
                  t1221 - 0.3177280975e1_dp * t648 * t1205 + 0.3365038132e0_dp * t643 * t937  &
                  * t1211 + 0.9708358534e0_dp * t278 * t906 * t1211 - 0.3365038132e0_dp  &
                  * t657 * t937 * t1221 + 0.3812737170e1_dp * t279 * t1235 + 0.2422827454e1_dp  &
                  * t658 * t1238 - 0.6730076265e0_dp * t657 * t859 * t1221 + 0.5132002385e0_dp  &
                  * bf * t1245 * t984 * t1252 - t109 * t1341 * t122
          t1346 = t13 * af * t1344
          t1358 = t356 * t305 * t333
          t1360 = t1085 + 0.3e1_dp * t1086 + (-t1109 - 0.120e2_dp * t1111 - 0.3600e2_dp  &
                  * t1114 - 0.120e2_dp * t1116 - 0.24000e2_dp * t1120 - 0.3600e2_dp * t356  &
                  * t133 * t340 - 0.40e1_dp * t1125) * ap * t54 + t1109 * t126 - 0.24000e2_dp  &
                  * t1120 * t92 + 0.120e2_dp * t1110 * t257 - 0.120e2_dp * t1116 * t92  &
                  + 0.3e1_dp * t1137 + 0.3600e2_dp * t771 * t772 * t340 + 0.240e2_dp * t1142  &
                  + t1346 - 0.120e2_dp * t1111 * t92 - 0.3600e2_dp * t1114 * t92 - 0.3600e2_dp  &
                  * t750 * t90 * t91 * t340 - 0.40e1_dp * t1125 * t92 + 0.3600e2_dp * t1358
          t1361 = t362 * t202
          t1363 = t353 * t254
          t1365 = t144 * t591
          t1368 = t859 * t61 * t865
          t1371 = t372 * t217
          t1372 = t366 * t212 * t1371
          t1374 = t493 * t874
          t1377 = t151 * t505
          t1380 = t379 * t217
          t1383 = t151 * t513
          t1388 = t211 ** 2
          t1389 = 0.1e1_dp / t1388
          t1391 = t505 * t217
          t1394 = t217 * t513
          t1404 = -t901 + t903 - t905 - 0.2427089633e0_dp * ba * t906 * t865 + 0.7943202439e0_dp  &
                  * t508 * t874 - 0.9531842928e0_dp * t214 * t887
                  t1455 = 0.8412595335e-1_dp * t1368 - 0.1514267160e0_dp * t1372 - 0.3028534320e0_dp  &
                  * t1374 - 0.1906368585e1_dp * t235 * t503 * t160 * t1377 + 0.7943202441e0_dp  &
                  * t552 * t393 * t1371 - 0.1906368585e1_dp * t553 * t1380  &
                  + 0.9531842928e0_dp * t553 * t1383 + 0.4206297667e-1_dp * t61 * t937 * t865  &
                  - 0.4854179269e0_dp * t236 * t942 + 0.1588640488e1_dp * t236 * t945  &
                  - 0.1906368586e1_dp * t236 * t948 - 0.6e1_dp * t75 * t1389 * t1391 + 0.6e1_dp  &
                  * t560 * t1394 - t239 * t1404
          t1469 = t74 ** (-0.40e1_dp)
          t1477 = t874 * t231
          t1480 = 0.1e1_dp / t540
          t1483 = t865 * t231
          t1490 = 0.1e1_dp / t540 / t224
          t1493 = t865 * t547
          t1507 = t887 * t231
          t1510 = t874 * t547
          t1516 = t540 ** 2
          t1517 = 0.1e1_dp / t1516
          t1522 = t228 ** 2
          t1524 = t865 / t545 / t230 * t1522
          t1527 = t1455 * t242 * t60 + 0.2e1_dp * t565 * t217 + 0.6354561952e0_dp *  &
                  t564 * t567 * t60 * t186 + 0.6354561952e0_dp * t568 * t217 * t186 -  &
                  0.6354561952e0_dp * t569 * t447 + 0.1514267160e0_dp * t241 * t1469 * t60  &
                  * t389 + 0.2647734147e0_dp * t569 * t444 + t243 * t513 - 0.1211413727e1_dp  &
                  * t574 * t1477 + 0.1924500894e0_dp * t80 * t1480 * t984 * t1483 +  &
                  0.3365038132e0_dp * t573 * t859 * t1483 - 0.4490502088e0_dp * t80 * t1490  &
                  * t984 * t1493 - 0.1588640487e1_dp * t577 * t1477 + 0.1682519066e0_dp  &
                  * t573 * t937 * t1483 + 0.4854179267e0_dp * t245 * t906 * t1483 - 0.1682519066e0_dp  &
                  * t582 * t937 * t1493 + 0.1906368585e1_dp * t246 * t1507  &
                  + 0.1211413727e1_dp * t583 * t1510 - 0.3365038132e0_dp * t582 * t859 *  &
                  t1493 + 0.2566001193e0_dp * t80 * t1517 * t984 * t1524
          t1567 = 0.1100642416e1_dp * (-0.3365038134e0_dp * t1368 + 0.6057068641e0_dp  &
                  * t1372 + 0.1211413728e1_dp * t1374 - 0.1817120592e1_dp * t149 * t503 *  &
                  t1377 - 0.1817120592e1_dp * t496 * t1380 + 0.9085602960e0_dp * t496 * t1383  &
                  - 0.1817120592e1_dp * t208 * t887 - 0.5451361779e1_dp * t18 * t1389  &
                  * t1391 + 0.5451361779e1_dp * t504 * t1394 - 0.9085602964e0_dp * t213 *  &
                  t1404) * t168 * t60 + 0.2201284832e1_dp * t517 * t217 + 0.7337616104e0_dp  &
                  * t516 * t406 * t522 + 0.7337616104e0_dp * t520 * t497 - 0.7337616104e0_dp  &
                  * t520 * t521 * t337 + 0.4891744068e0_dp * t220 * t1043 * t60 *  &
                  t369 * t371 - t73 * t1527 * t86 + 0.1100642416e1_dp * t221 * t513 - 0.2422827454e1_dp  &
                  * t530 * t1477 + 0.3849001789e0_dp * ba * t1480 * t984  &
                  * t1483 + 0.6730076265e0_dp * t529 * t859 * t1483 - 0.8981004177e0_dp *  &
                  ba * t1490 * t984 * t1493 - 0.3177280975e1_dp * t534 * t1477 + 0.3365038132e0_dp  &
                  * t529 * t937 * t1483 + 0.9708358534e0_dp * t226 * t906 * t1483  &
                  - 0.3365038132e0_dp * t543 * t937 * t1493 + 0.3812737170e1_dp * t227  &
                  * t1507 + 0.2422827454e1_dp * t544 * t1510 - 0.6730076265e0_dp * t543 *  &
                  t859 * t1493 + 0.5132002385e0_dp * ba * t1517 * t984 * t1524
          t1570 = t57 * aa * t1567 * t91
          t1576 = t481 * t254
          t1578 = t357 * t254
          t1580 = t204 * t591
          t1582 = t359 * t254
          t1588 = t143 * t305 * t340
          t1590 = t141 * t704
          t1593 = t143 * t704 * t133
          t1599 = 0.3e1_dp * t1361 - 0.240e2_dp * t1363 - 0.120e2_dp * t1365 + t1570 +  &
                  0.40e1_dp * t143 * t126 * t1096 + t1108 * t56 * t92 + 0.3e1_dp * t1576  &
                  - 0.3600e2_dp * t1578 + 0.3e1_dp * t1580 - 0.120e2_dp * t1582 + 0.24000e2_dp *  &
                  t1119 * t126 * t1089 + 0.120e2_dp * t1588 + 0.3e1_dp * t1590 + 0.120e2_dp  &
                  * t1593 + 0.120e2_dp * t352 * t604 + 0.3600e2_dp * t1113 * t598

          e_aaa(ip) = e_aaa(ip) + ( t842 + 0.6e1_dp * t364 + 0.3e1_dp * t363 + 0.3e1_dp * t482 + 0.6e1_dp *  &
                      t485 - 0.240e2_dp * t483 - 0.3600e2_dp * t487 - 0.120e2_dp * t489 - 0.240e2_dp  &
                      * t491 + t851 + 0.3e1_dp * t593 + 0.6e1_dp * t596 + 0.240e2_dp * t594 + 0.3600e2_dp  &
                      * t599 + 0.240e2_dp * t602 + t857 + 0.120e2_dp * t605 + (t1360 +  &
                      t1599) * t2 ) * sc

          t1602 = 0.2400e2_dp * t753
          t1603 = 0.2e1_dp * t766
          t1604 = 0.80e1_dp * t746
          t1605 = 0.2e1_dp * t745
          t1606 = 0.80e1_dp * t748
          t1609 = 0.160e2_dp * t759
          t1610 = -t1602 + t1603 - t1604 + t1605 + t818 - t490 - t1606 + t827  &
                  + t363 - t488 + t831 + t813 - t823 - 0.160e2_dp * t491 + 0.4e1_dp * t364  &
                  + t842 - t1609
          t1611 = 0.80e1_dp * t767
          t1613 = t322 * t591
          t1615 = 0.80e1_dp * t732 * t601
          t1620 = 0.80e1_dp * t733 * t254
          t1622 = 0.80e1_dp * t352 * t783
          t1626 = t732 * t340
          t1639 = 0.2e1_dp * t337 - 0.6e1_dp * t1094
          t1653 = -0.5699736440e0_dp * t1088 * t333 * t309 + 0.3419841862e1_dp * t716  &
                  * t338 + 0.8549604655e0_dp * t332 * t340 * t309 + 0.2564881399e1_dp *  &
                  t129 * t1639 - 0.5699736440e0_dp * t1099 * t344 * t312 - 0.3419841862e1_dp  &
                  * t722 * t338 + 0.8549604655e0_dp * t343 * t347 * t312 - 0.2564881399e1_dp  &
                  * t136 * t1639
          t1654 = t1653 * t12
          t1657 = t143 * t704 * t309
          t1659 = t1085 + 0.2e1_dp * t1086 + t1613 + t1615 - 0.160e2_dp * t352 * t1  &
                  * t758 - t1620 + t1622 + 0.40e1_dp * t1110 * t327 + t1137 + 0.80e1_dp *  &
                  t1142 - 0.40e1_dp * t1626 * t92 + t1346 + t1654 * t126 + 0.40e1_dp * t1657
          t1664 = 0.160e2_dp * t777 * t304 * t1 * t337
          t1666 = 0.80e1_dp * t730 * t254
          t1668 = t143 * t1639
          t1671 = t728 * t142
          t1672 = t1671 * t133
          t1675 = t1110 * t309
          t1690 = 0.2400e2_dp * t771 * t304 * t133 * t309
          t1694 = 0.1200e2_dp * t1358 + t1664 + t1361 - t1666 - 0.80e1_dp * t1363 -  &
                  0.40e1_dp * t1668 * t92 - 0.80e1_dp * t1672 * t92 - 0.40e1_dp * t1675 * t92  &
                  - 0.2400e2_dp * t1113 * t133 * t752 - 0.80e1_dp * t1365 + t1570 - 0.4800e2_dp  &
                  * t356 * t133 * aa * t757 * t338 + t1690 + 0.40e1_dp * t143 * t126  &
                  * t1639
          t1696 = t316 * t704
          t1697 = t315 * t355
          t1698 = t1697 * t333
          t1701 = t1119 * af
          t1726 = t1696 - 0.1200e2_dp * t1698 * t92 + 0.24000e2_dp * t1701 * t125 *  &
                  t333 * t309 + t1653 * t56 * t92 + 0.2400e2_dp * t1113 * af * t773 +  &
                  0.4800e2_dp * t771 * t772 * t338 + 0.160e2_dp * t352 * af * t779 + 0.40e1_dp  &
                  * t732 * t604 + t1576 - 0.1200e2_dp * t1578 + 0.1200e2_dp * t1697 * t598  &
                  + 0.2e1_dp * t1580 - 0.40e1_dp * t1582 + 0.80e1_dp * t1671 * t257
          t1730 = 0.160e2_dp * t755 * t756 * t252 * t91
          t1750 = -t1654 - 0.40e1_dp * t1675 - 0.80e1_dp * t1672 - 0.2400e2_dp * t1113  &
                  * t735 - 0.160e2_dp * t352 * t338 - 0.1200e2_dp * t1698 - 0.24000e2_dp * t1119  &
                  * t333 * t309 - 0.4800e2_dp * t356 * t133 * t1 * t337 - 0.40e1_dp *  &
                  t1626 - 0.1200e2_dp * t356 * t340 * t309 - 0.40e1_dp * t1668
          t1754 = 0.2e1_dp * t744 * t254
          t1764 = 0.2e1_dp * t741 * t202
          t1765 = t320 * t479
          t1768 = 0.2400e2_dp * t750 * t253 * t751
          t1775 = 0.2e1_dp * t729 * t305
          t1776 = t317 * t591
          t1778 = -t1730 + t1750 * ap * t54 + t1754 - 0.24000e2_dp * t1119 * t333  &
                  * t752 + 0.40e1_dp * t1588 - 0.1200e2_dp * t356 * t340 * t752 + 0.2e1_dp  &
                  * t1590 + t1764 + t1765 - t1768 + 0.80e1_dp * t1593 + 0.1200e2_dp * t771  &
                  * t125 * t340 * t309 + t1775 - 0.40e1_dp * t1776
          t1782 = 0.2e1_dp * t742
          t1783 = 0.160e2_dp * t780
          t1785 = 0.2400e2_dp * t774
          t1786 = 0.80e1_dp * t769
          t1789 = t606 + t1611 + t600 + t482 + t851 + t595 + (t1659 + t1694 +  &
                  t1726 + t1778) * t2 + t1782 + t1783 + t593 + 0.4e1_dp * t596 + t857  &
                  - t484 + t1785 + t1786 + 0.4e1_dp * t485 + 0.160e2_dp * t602

          e_aab(ip) = e_aab(ip) + ( t1610 + t1789 ) * sc

          t1795 = -t1602 + t826 + t812 + t824 + t834 + t1603 - t1604 + t1605  &
                  + 0.4e1_dp * t762 - t1606 + 0.4e1_dp * t782 + t830 + 0.160e2_dp * t784 + 0.4e1_dp  &
                  * t743 - 0.160e2_dp * t764 - t492 + t365
          t1797 = t143 * t305 * t793
          t1804 = t337 * t309
          t1826 = -0.5699736440e0_dp * t1088 * t133 * t789 + 0.3419841862e1_dp * t332  &
                  * t1 * t1804 + 0.8549604655e0_dp * t716 * t793 - 0.5129762798e1_dp *  &
                  t129 * t337 - 0.1538928839e2_dp * t719 * t371 - 0.5699736440e0_dp * t1099  &
                  * t137 * t796 - 0.3419841862e1_dp * t343 * t1 * t337 * t312 + 0.8549604655e0_dp  &
                  * t722 * t799 + 0.5129762798e1_dp * t136 * t337 + 0.1538928839e2_dp  &
                  * t725 * t371
          t1829 = t352 * t793
          t1832 = t814 * t254
          t1836 = t732 * t783
          t1838 = t811 * t202
          t1839 = t1085 + t1086 + 0.40e1_dp * t1797 + 0.2e1_dp * t1613 + t1615 + t1826  &
                  * t56 * t92 - 0.40e1_dp * t1829 * t92 - t1620 + t1832 + t1622 + 0.24000e2_dp  &
                  * t1701 * t772 * t789 + 0.80e1_dp * t1836 + t1346 + t1838
          t1850 = t806 * t254
          t1853 = t356 * t305 * t789
          t1858 = t802 * t142
          t1868 = -0.80e1_dp * t143 * t126 * t337 + 0.240e2_dp * t755 * t371 * aa *  &
          t757 - 0.2400e2_dp * t1697 * t133 * t752 - 0.1200e2_dp * t1850 + 0.1200e2_dp  &
                  * t1853 + 0.80e1_dp * t1657 + t1664 - t1666 - 0.40e1_dp * t1365 + t1570  &
                  + t1690 + 0.2e1_dp * t1696 + 0.40e1_dp * t1858 * t257 - 0.24000e2_dp * t1119  &
                  * t133 * t90 * t91 * t789 + 0.80e1_dp * t1671 * t327
          t1870 = t804 * t254
          t1872 = t143 * t337
          t1884 = t1826 * t12
          t1886 = t1671 * t309
          t1891 = t808 * t254
          t1893 = t803 * t305
          t1894 = t1113 * t789
          t1897 = t1858 * t133
          t1901 = t90 * t91 * t793
          t1904 = -0.80e1_dp * t1870 + 0.80e1_dp * t1872 * t92 - 0.160e2_dp * t732 * t1  &
                  * t758 + 0.1200e2_dp * t771 * t772 * t793 + 0.2400e2_dp * t1697 * af *  &
                  t773 + t1884 * t126 - 0.80e1_dp * t1886 * t92 + 0.40e1_dp * t352 * t832  &
                  - 0.40e1_dp * t1891 + t1893 + t1580 - 0.1200e2_dp * t1894 * t92 - 0.40e1_dp  &
                  * t1897 * t92 - 0.1200e2_dp * t750 * t1901
          t1939 = -t1884 - 0.80e1_dp * t1886 - 0.1200e2_dp * t1894 - 0.40e1_dp * t1829  &
                  - 0.40e1_dp * t1897 - 0.2400e2_dp * t1697 * t735 - 0.160e2_dp * t732 * t338  &
                  - 0.24000e2_dp * t1119 * t133 * t789 - 0.4800e2_dp * t356 * t338 * t309  &
                  - 0.1200e2_dp * t356 * t133 * t793 + 0.80e1_dp * t1872 + 0.240e2_dp * t143  &
                  * t1094
          t1945 = -t1730 - 0.240e2_dp * t777 * t778 * t371 + t1754 - 0.4800e2_dp *  &
                  t356 * t338 * t752 + 0.160e2_dp * t732 * af * t779 + t1590 + t1764 +  &
                  0.2e1_dp * t1765 - t1768 + 0.40e1_dp * t1593 + 0.4800e2_dp * t771 * t778 *  &
                  t1804 + t1939 * ap * t54 + 0.1200e2_dp * t1113 * t828 + t1775 - 0.80e1_dp  &
                  * t1776
          t1949 = t842 - t1609 + t815 + t1611 + t851 + t1782 + t1783 + t597 +  &
                  t857 + (t1839 + t1868 + t1904 + t1945) * t2 - t817 + t1785 + t1786  &
                  - t820 + t486 + t603 - t822

          e_abb(ip) = e_abb(ip) + ( t1795 + t1949 ) * sc

          t1975 = t1697 * t789
          t1979 = t789 * t309
          t1986 = -0.6e1_dp * t337 - 0.6e1_dp * t1094
          t1998 = -0.5699736440e0_dp * t1088 * t1979 + 0.2564881396e1_dp * t332 * t309  &
                  * t793 + 0.2564881399e1_dp * t129 * t1986 - 0.5699736440e0_dp * t1099  &
                  * t796 * t312 + 0.2564881396e1_dp * t343 * t312 * t799 - 0.2564881399e1_dp  &
                  * t136 * t1986
          t1999 = t1998 * t12
          t2010 = t1085 + 0.120e2_dp * t1797 + 0.3e1_dp * t1613 - 0.3600e2_dp * t356 *  &
                  t309 * t1901 + 0.3600e2_dp * t771 * t125 * t309 * t793 + 0.3e1_dp * t1832  &
                  + 0.240e2_dp * t1836 - 0.3600e2_dp * t1975 * t92 + t1346 + 0.3e1_dp * t1838  &
                  + t1999 * t126 + 0.40e1_dp * t143 * t126 * t1986 - 0.3600e2_dp * t1850  &
                  + 0.3600e2_dp * t1853 + 0.120e2_dp * t1657 + 0.24000e2_dp * t1119 * t126  &
                  * t1979
          t2011 = t1858 * t309
          t2014 = t732 * t793
          t2019 = t143 * t1986
          t2025 = t1119 * t1979
          t2048 = -0.120e2_dp * t2011 * t92 + t1570 - 0.120e2_dp * t2014 * t92 + 0.3e1_dp  &
                  * t1696 - 0.240e2_dp * t1870 - 0.40e1_dp * t2019 * t92 + (-t1999 - 0.120e2_dp  &
                  * t2011 - 0.3600e2_dp * t1975 - 0.120e2_dp * t2014 - 0.24000e2_dp *  &
                  t2025 - 0.3600e2_dp * t356 * t309 * t793 - 0.40e1_dp * t2019) * ap * t54  &
                  - 0.24000e2_dp * t2025 * t92 - 0.120e2_dp * t1891 + 0.3e1_dp * t1893 + 0.3600e2_dp  &
                  * t1697 * t828 + 0.120e2_dp * t732 * t832 + t1998 * t56 * t92 +  &
                  0.3e1_dp * t1765 + 0.120e2_dp * t1858 * t327 - 0.120e2_dp * t1776

          e_bbb(ip) = e_bbb(ip) + ( 0.3e1_dp * t812 + 0.6e1_dp * t743 + t842 + t851 + t857 + 0.6e1_dp * t762  &
                      - 0.240e2_dp * t764 + 0.6e1_dp * t782 + 0.240e2_dp * t784 + 0.3e1_dp * t815  &
                      - 0.240e2_dp * t816 - 0.3600e2_dp * t819 - 0.120e2_dp * t821 + 0.3e1_dp * t824  &
                      + 0.240e2_dp * t825 + 0.3600e2_dp * t829 + 0.120e2_dp * t833 + (t2010 +  &
                      t2048) * t2 ) * sc

        END IF
      END IF
    END DO

    !$omp end do

  END SUBROUTINE vwn5_lsd_calc

END MODULE xc_vwn

