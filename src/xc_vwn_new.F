!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/xc_vwn [1.0] *
!!
!!   NAME
!!     xc_vwn
!!
!!   FUNCTION
!!     Calculate the LDA functional according to Vosk, Wilk and Nusair
!!     Literature: S. H. Vosko, L. Wilk and M. Nusair,
!!                 Can. J. Phys. 58, 1200 (1980)
!!
!!   NOTES
!!     Order of derivatives is: LDA 0; 1; 2; 3;
!!
!!   AUTHOR
!!     JGH (26.02.2002)
!!
!!   MODIFICATION HISTORY
!!     JGH (26.02.2003) : OpenMP enabled
!!     fawzi (04.2004)  : adapted to the new xc interface
!!
!!*** *************************************************************************

MODULE xc_vwn_new
! *****************************************************************************
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
  USE kinds,                           ONLY: dbl
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_functionals_utilities,        ONLY: calc_srs_pw,&
                                             set_util
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
                                             xc_rho_set_get,&
                                             xc_rho_set_type

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f43 = 4._dbl*f13

  REAL(dbl), PARAMETER :: a  = 0.0310907_dbl,&
                          b  = 3.72744_dbl,&
                          c  = 12.9352_dbl,&
                          x0 = -0.10498_dbl

  PUBLIC :: vwn_lda_info, vwn_lda_eval

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='xc_vwn_new'

  REAL(dbl) :: eps_rho

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE vwn_init ( cutoff )

    REAL(dbl), INTENT(IN)                    :: cutoff

    eps_rho = cutoff
    CALL set_util ( cutoff )

  END SUBROUTINE vwn_init

! *****************************************************************************

  SUBROUTINE vwn_lda_info ( reference, shortform, needs, max_deriv, error)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

  SUBROUTINE vwn_lda_eval(rho_set,deriv_set,order,error)
    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(in)                      :: order
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vwn_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: x
    REAL(kind=dbl)                           :: epsilon_rho
    REAL(kind=dbl), DIMENSION(:, :, :), &
      POINTER                                :: e_0, e_rho, e_rho_rho, &
                                                e_rho_rho_rho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(bo)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(deriv_set%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL xc_rho_set_get(rho_set,rho=rho,&
            local_bounds=bo,rho_cutoff=epsilon_rho,&
            error=error)
       npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)
       CALL vwn_init(epsilon_rho)

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

          CALL vwn_lda_01 ( rho, x, e_0,e_rho,npoints,error )
       ELSEIF ( order>=0 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"",&
               allocate_deriv=.TRUE., error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_0,error=error)
          
          CALL vwn_lda_0 ( rho, x, e_0,npoints,error )
       ELSEIF ( order==-1 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho,error=error)

          CALL vwn_lda_1 ( rho, x, e_rho, npoints, error )
       END IF
       IF ( order>=2.OR.order==-2 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho,error=error)

          CALL vwn_lda_2 ( rho, x, e_rho_rho, npoints, error )
       END IF
       IF ( order>=3.OR.order==-3 ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
               allocate_deriv=.TRUE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho,error=error)

          CALL vwn_lda_3 ( rho, x, e_rho_rho_rho, npoints,&
               error)
       END IF
       IF ( order>3.OR.order<-3) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="derivatives bigger than 3 not implemented", &
               error=error, error_level=cp_failure_level)
       END IF

       DEALLOCATE(x,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE vwn_lda_eval

! *****************************************************************************

  SUBROUTINE vwn_lda_0 ( rho, x, e_0, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, x
    REAL(dbl), DIMENSION(*), INTENT(INOUT)   :: e_0
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: at, dpx, ln1, ln2, px, q, xb

!   ---------------------------------------------------------------------------

    q = SQRT ( 4._dbl*c - b*b )
    xb = 2._dbl*x0+b

!$omp parallel do private(ip,px,dpx,at,ln1,ln2)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2._dbl*x(ip)+b
         at = 2._dbl/q*ATAN(q/dpx)
         ln1 = LOG(x(ip)*x(ip)/px)
         ln2 = LOG((x(ip)-x0)**2/px)
         e_0(ip) = e_0(ip) + a * ( ln1 + b*at - b*x0/px* ( ln2 + xb * at ) ) * rho(ip)
      END IF

    END DO

  END SUBROUTINE vwn_lda_0

! *****************************************************************************

  SUBROUTINE vwn_lda_1 ( rho, x, e_rho, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, x
    REAL(dbl), DIMENSION(*), INTENT(INOUT)   :: e_rho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: at, dat, dex, dln1, dln2, &
                                                dpx, ex, ln1, ln2, pa, px, q, &
                                                xb

!   ---------------------------------------------------------------------------

    q = SQRT ( 4._dbl*c - b*b )
    xb = 2._dbl*x0+b
!$omp parallel do private(ip,px,dpx,at,pa,dat,ln1,dln1,ln2,dln2,ex,dex)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2._dbl*x(ip)+b
         at = 2._dbl/q*ATAN(q/dpx)
         pa = 4._dbl*x(ip)*x(ip)+4._dbl*b*x(ip)+b*b+q*q
         dat= -4._dbl/pa
         ln1 = LOG(x(ip)*x(ip)/px)
         dln1 = (b*x(ip)+2._dbl*c)/(x(ip)*px)
         ln2 = LOG((x(ip)-x0)**2/px)
         dln2 = (b*x(ip)+2._dbl*c+2._dbl*x0*x(ip)+x0*b)/((x(ip)-x0)*px)
         ex = a * ( ln1 + b*at - b*x0/px* ( ln2 + xb * at ) )
         dex = a * ( dln1 + b*dat + b*x0*dpx/(px*px) * ( ln2 + xb * at ) - &
               b*x0/px* ( dln2 + xb * dat ) )
         e_rho(ip) = e_rho(ip) + ex - x(ip)*dex/6._dbl
      END IF

    END DO

  END SUBROUTINE vwn_lda_1

! *****************************************************************************

  SUBROUTINE vwn_lda_01 ( rho, x, e_0, e_rho, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, x
    REAL(dbl), DIMENSION(*), INTENT(INOUT)   :: e_0, e_rho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: at, dat, dex, dln1, dln2, &
                                                dpx, ex, ln1, ln2, pa, px, q, &
                                                xb

!   ---------------------------------------------------------------------------

    q = SQRT ( 4._dbl*c - b*b )
    xb = 2._dbl*x0+b
!$omp parallel do private(ip,px,dpx,at,pa,dat,ln1,dln1,ln2,dln2,ex,dex)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2._dbl*x(ip)+b
         at = 2._dbl/q*ATAN(q/dpx)
         pa = 4._dbl*x(ip)*x(ip)+4._dbl*b*x(ip)+b*b+q*q
         dat= -4._dbl/pa
         ln1 = LOG(x(ip)*x(ip)/px)
         dln1 = (b*x(ip)+2._dbl*c)/(x(ip)*px)
         ln2 = LOG((x(ip)-x0)**2/px)
         dln2 = (x(ip)*(b+2._dbl*x0)+2._dbl*c+x0*b)/((x(ip)-x0)*px)
         ex = a * ( ln1 + b*at - b*x0/px* ( ln2 + xb * at ) )
         dex = a * ( dln1 + b*dat + b*x0*dpx/(px*px) * ( ln2 + xb * at ) - &
               b*x0/px* ( dln2 + xb * dat ) )
         e_0(ip) = e_0(ip) + ex * rho(ip)
         e_rho(ip) = e_rho(ip) + ex - x(ip)*dex/6._dbl
      END IF

    END DO


  END SUBROUTINE vwn_lda_01

! *****************************************************************************

  SUBROUTINE vwn_lda_2 ( rho, x, e_rho_rho, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, x
    REAL(dbl), DIMENSION(*), INTENT(INOUT)   :: e_rho_rho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip
    REAL(dbl)                                :: at, d2at, d2ex, d2fp, d2ln1, &
                                                d2ln2, dat, dex, dfp, dln1, &
                                                dln2, dpx, fp, ln1, ln2, pa, &
                                                px, q, xb

!   ---------------------------------------------------------------------------

    q = SQRT ( 4._dbl*c - b*b )
    xb = 2._dbl*x0+b

!$omp parallel do private(ip,px,dpx,at,pa,dat,d2at,ln1,dln1,d2ln1) &
!$omp             private(ln2,dln2,d2ln2,fp,dfp,d2fp,dex,d2ex)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2._dbl*x(ip)+b
         at = 2._dbl/q*ATAN(q/dpx)
         pa = 4._dbl*x(ip)*x(ip)+4._dbl*b*x(ip)+b*b+q*q
         dat= -4._dbl/pa
         d2at= 16._dbl*dpx/(pa*pa)
         ln1 = LOG(x(ip)*x(ip)/px)
         dln1 = (b*x(ip)+2._dbl*c)/(x(ip)*px)
         d2ln1 = b/(x(ip)*px) - (b*x(ip)+2._dbl*c)/(x(ip)*px)**2 * (px+x(ip)*dpx)
         ln2 = LOG((x(ip)-x0)**2/px)
         dln2 = (x(ip)*xb+2._dbl*c+x0*b)/((x(ip)-x0)*px)
         d2ln2 = xb/((x(ip)-x0)*px)-(x(ip)*xb+2._dbl*c+x0*b)/((x(ip)-x0)*px)**2 &
                 * (px+(x(ip)-x0)*dpx)
         fp = -b*x0/px
         dfp = b*x0*dpx/(px*px)
         d2fp = b*x0*2._dbl*(1._dbl-dpx*dpx/px)/(px*px)
         dex = a * ( dln1 + b*dat + dfp * ( ln2 + xb * at ) + fp * ( dln2 + xb * dat ) )
         d2ex= a * ( d2ln1 + b*d2at + d2fp * ( ln2 + xb * at ) + &
               2._dbl*dfp * ( dln2 + xb * dat ) + fp * ( d2ln2 + xb * d2at ) )
         e_rho_rho(ip) = e_rho_rho(ip)&
              + x(ip)/(36._dbl*rho(ip))*(x(ip)*d2ex-5._dbl*dex)
      END IF

    END DO

  END SUBROUTINE vwn_lda_2

! *****************************************************************************

  SUBROUTINE vwn_lda_3 ( rho, x, e_rho_rho_rho, npoints, error )

    REAL(dbl), DIMENSION(*), INTENT(IN)      :: rho, x
    REAL(dbl), DIMENSION(*), INTENT(INOUT)   :: e_rho_rho_rho
    INTEGER, INTENT(in)                      :: npoints
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ip
    REAL(dbl) :: at, ax, bx, cx, d2at, d2bx, d2dx, d2ex, d2fp, d2ln1, d2ln2, &
      d3at, d3ex, d3fp, d3ln1, d3ln2, dat, dbx, ddx, dex, dfp, dln1, dln2, &
      dpx, dx, fp, ln1, ln2, pa, px, q, xb

!   ---------------------------------------------------------------------------

    q = SQRT ( 4._dbl*c - b*b )
    xb = 2._dbl*x0+b

!$omp parallel do private(ip,px,dpx,at,pa,dat,d2at,d3at,ln1,ax,bx) &
!$omp             private(dbx,d2bx,dln1,d2ln1,d3ln1,ln2,cx,dx,ddx,d2dx,dln2,d2ln2,d3ln2) &
!$omp             private(fp,dfp,d2fp,d3fp,dex,d2ex,d3ex)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         px = x(ip)*x(ip)+b*x(ip)+c
         dpx = 2._dbl*x(ip)+b
         at = 2._dbl/q*ATAN(q/dpx)
         pa = 4._dbl*x(ip)*x(ip)+4._dbl*b*x(ip)+b*b+q*q
         dat= -4._dbl/pa
         d2at= 16._dbl*dpx/(pa*pa)
         d3at= 32._dbl/(pa*pa)*(1._dbl-4._dbl*dpx*dpx/pa)
         ln1 = LOG(x(ip)*x(ip)/px)
         ax = b*x(ip)+2._dbl*c
         bx = x(ip)*px
         dbx = px + x(ip)*dpx
         d2bx = 2._dbl*(dpx+x(ip))
         dln1 = ax/bx
         d2ln1 = (b*bx-ax*dbx)/(bx*bx)
         d3ln1 = -ax*d2bx/(bx*bx) - 2._dbl*d2ln1*dbx/bx
         ln2 = LOG((x(ip)-x0)**2/px)
         cx = x(ip)*xb+2._dbl*c+x0*b
         dx = (x(ip)-x0)*px
         ddx = px + (x(ip)-x0)*dpx
         d2dx = 2._dbl*(dpx + (x(ip)-x0))
         dln2 = cx/dx
         d2ln2 = (xb*dx-cx*ddx)/(dx*dx)
         d3ln2 = -cx*d2dx/(dx*dx) - 2._dbl*d2ln2*ddx/dx
         fp = -b*x0/px
         dfp = b*x0*dpx/(px*px)
         d2fp = b*x0*2._dbl*(1._dbl-dpx*dpx/px)/(px*px)
         d3fp = b*x0*6._dbl*dpx*(dpx*dpx/px-2._dbl)/(px*px*px)
         dex = a * ( dln1 + b*dat + dfp * ( ln2 + xb * at ) + fp * ( dln2 + xb * dat ) )
         d2ex= a * ( d2ln1 + b*d2at + d2fp * ( ln2 + xb * at ) + &
               2._dbl*dfp * ( dln2 + xb * dat ) + fp * ( d2ln2 + xb * d2at ) )
         d3ex= a * ( d3ln1 + b*d3at + d3fp * ( ln2 + xb * at ) + &
               3._dbl*d2fp * ( dln2 + xb * dat ) + 3._dbl*dfp * ( d2ln2 + xb * d2at ) + &
               fp * ( d3ln2 + xb * d3at ) )
         e_rho_rho_rho(ip) = e_rho_rho_rho(ip)&
              - 7._dbl*x(ip)/(216._dbl*rho(ip)*rho(ip))*(x(ip)*d2ex-5._dbl*dex) - &
              x(ip)*x(ip)/(216._dbl*rho(ip)*rho(ip))*(x(ip)*d3ex-4._dbl*d2ex)
      END IF

    END DO

  END SUBROUTINE vwn_lda_3

! *****************************************************************************

END MODULE xc_vwn_new

! *****************************************************************************
