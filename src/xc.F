!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001 - 2002  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****s* cp2k/xc [1.0] *
!!
!!   NAME
!!     xc
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations
!!
!!   AUTHOR
!!     JGH (13-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE xc
  USE coefficient_types,               ONLY: coeff_add,&
                                             coeff_allocate,&
                                             coeff_copy,&
                                             coeff_deallocate,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error,&
                                             cp_unimplemented_error_nr
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_output_handling,              ONLY: cp_finished_output,&
                                             cp_unitnr_for_output,&
                                             cp_would_output
  USE dft_types,                       ONLY: dft_control_type
  USE kinds,                           ONLY: dbl,&
                                             wp => dp
  USE message_passing,                 ONLY: mp_sum
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_compare
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             NOSPACE,&
                                             REALDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_allocate,&
                                             pw_copy,&
                                             pw_deallocate,&
                                             pw_derive,&
                                             pw_p_type,&
                                             pw_sumup,&
                                             pw_transfer,&
                                             pw_type,&
                                             pw_zero
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_functionals,                  ONLY: xc_get_pot_size,&
                                             xc_lda,&
                                             xc_set
  PRIVATE
  PUBLIC :: xc_calculate

  INTERFACE xc_calculate
    MODULE PROCEDURE xc_calculate_pw_lda
  END INTERFACE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='xc'

!!*****
!******************************************************************************
  
CONTAINS

!******************************************************************************
!!****** xc/xc_calculate_pw [1.0] *
!!
!!   NAME
!!     xc_calculate_pw
!!
!!   FUNCTION
!!     Exchange and Correlation functional calculations
!!
!!   AUTHOR
!!     JGH (13-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (13-Jun-2002): adaptation to new functionals
!!     Fawzi (11.2002): drho_g(1:3)->drho_g
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE xc_calculate_pw_lda ( rho_r, rho_g, vxc, exc, dft_control, error )
  
  IMPLICIT NONE

! Arguments
  TYPE(pw_type), INTENT(IN) :: rho_r
  TYPE(pw_type), INTENT(IN) :: rho_g
  TYPE(pw_type), INTENT(INOUT) :: vxc
  REAL(dbl), INTENT(OUT) :: exc
  TYPE(dft_control_type), INTENT(IN) :: dft_control
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

! Local
  TYPE(pw_type), DIMENSION(3) :: drho_r
  TYPE(pw_type) :: vxc_g, drho_g
  REAL ( dbl ) , DIMENSION (:), ALLOCATABLE :: rgrid, drgrid
  REAL ( dbl ) , DIMENSION (:,:), ALLOCATABLE :: pot
  INTEGER :: stat, ntot, n(3), i, j, k, ii, order, ik, npot
  INTEGER, DIMENSION (3,3) :: nd 
  INTEGER, DIMENSION (2,3) :: bo
  LOGICAL :: failure, condition, spezgrid
  CHARACTER(len=*), PARAMETER :: routineN='xc_calculate_pw_lda',&
         routineP=moduleN//':'//routineN
!------------------------------------------------------------------------------
  failure=.FALSE.
  nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))
  ! check if the correct density is provided
  condition=(rho_r%in_use==REALDATA3D .AND. rho_r%in_space==REALSPACE)
  CPPrecondition(condition,cp_failure_level,routineP,error,failure)
  ! CPPrecondition(dft_control%nspins==1,cp_failure_level,routineP,error,failure)
  IF ( dft_control%gradient_functional ) THEN
    condition=(rho_g%in_use==COMPLEXDATA1D .AND. rho_g%in_space==RECIPROCALSPACE)
    CPPrecondition(condition,cp_failure_level,routineP,error,failure)
  END IF
  IF ( failure ) CALL stop_program ( routineP, "Assertions" )
  ! is a special grid needed
  spezgrid = .NOT. dft_control%integration_grid%pw_default
  ! transform to higher cutoff
  IF ( spezgrid ) THEN
     CALL stop_program ( routineP, "NYI" )
  ELSE
     ! allocate and copy to grids
     n(1:3) = rho_r%pw_grid%npts_local (1:3)
     ntot = PRODUCT ( n )
     ALLOCATE ( rgrid ( ntot ), STAT = stat )
     IF ( stat /= 0 ) CALL stop_memory ( routineP, "rgrid", ntot )
     bo = rho_r%pw_grid%bounds_local  
     ii = 0
     DO k = bo(1,3), bo(2,3)
        DO j = bo(1,2), bo(2,2)
           DO i = bo(1,1), bo(2,1)
              ii = ii + 1
              rgrid(ii) = rho_r%cr3d(i,j,k)
           END DO
        END DO
     END DO
     IF ( dft_control%gradient_functional ) THEN
        ALLOCATE ( drgrid ( ntot ), STAT = stat )
        IF ( stat /= 0 ) CALL stop_memory ( routineP, "drgrid", ntot )
        ! calculate gradients if needed
        CALL pw_allocate ( drho_g, rho_g%pw_grid, COMPLEXDATA1D )
        drho_g%in_space=RECIPROCALSPACE
        DO ik = 1, 3
           CALL pw_copy ( rho_g, drho_g )
           CALL pw_derive ( drho_g, nd(:,ik) )
           CALL pw_allocate ( drho_r(ik), rho_g%pw_grid, REALDATA3D )
           drho_r(ik)%in_space=REALSPACE
           CALL pw_transfer ( drho_g, drho_r(ik) )
        END DO
        ii = 0
        DO k = bo(1,3), bo(2,3)
           DO j = bo(1,2), bo(2,2)
              DO i = bo(1,1), bo(2,1)
                 ii = ii + 1
                 drgrid(ii) = SQRT(drho_r(1)%cr3d(i,j,k)**2 + &
                                   drho_r(2)%cr3d(i,j,k)**2 + &
                                   drho_r(3)%cr3d(i,j,k)**2)
              END DO
           END DO
        END DO
     END IF
  END IF
  order = 1
  ! allocate potential
  npot = xc_get_pot_size ( .FALSE., dft_control%gradient_functional, .FALSE., order )
  ALLOCATE ( pot ( ntot, 0:npot-1 ), STAT = stat )
  IF ( stat /= 0 ) CALL stop_memory ( routineP, "pot", ntot*npot )
  pot = 0._dbl
! calculate functionals
  CALL xc_set ( dft_control%density_cut, dft_control%gradient_cut )
  IF (dft_control%gradient_functional) THEN
     CALL xc_lda ( functional=dft_control%functional, &
          ftags=dft_control%ftags, &
          rho=rgrid, drho=drgrid, pot=pot, order=order, error=error )
  ELSE
     CALL xc_lda ( functional=dft_control%functional, &
          ftags=dft_control%ftags, &
          rho=rgrid, pot=pot, order=order, error=error )
  END IF
  exc = SUM ( pot(:,0) )
  IF ( rho_r%pw_grid%para%mode == 1 ) THEN
     CALL mp_sum ( exc, rho_r%pw_grid%para%group )
  END IF
  ! copy back to structure
  IF ( spezgrid ) THEN
     ! transform to lower cutoff
     CALL stop_program ( routineP, "NYI" )
  ELSE
     ii = 0
     DO k = bo(1,3), bo(2,3)
        DO j = bo(1,2), bo(2,2)
           DO i = bo(1,1), bo(2,1)
              ii = ii + 1
              vxc%cr3d(i,j,k) = pot(ii,1)
           END DO
        END DO
     END DO
     IF ( dft_control%gradient_functional ) THEN
        CALL pw_allocate ( vxc_g, rho_g%pw_grid, COMPLEXDATA1D )
        vxc_g%in_space=RECIPROCALSPACE
        CALL pw_transfer ( vxc, vxc_g )
        DO i = 1, ntot
          IF ( drgrid(i) > dft_control%gradient_cut ) THEN
            pot(i,2) = -pot(i,2) / drgrid(i)
          ELSE
            pot(i,2) = 0._dbl
          END IF
        END DO
        DO ik = 1, 3
          ii = 0
          DO k = bo(1,3), bo(2,3)
             DO j = bo(1,2), bo(2,2)
                DO i = bo(1,1), bo(2,1)
                   ii = ii + 1
                   drho_r(ik)%cr3d(i,j,k) = pot(ii,2)*drho_r(ik)%cr3d(i,j,k)
                END DO
             END DO
          END DO
          CALL pw_transfer ( drho_r(ik), drho_g )
          CALL pw_derive ( drho_g, nd(:,ik) )
          CALL pw_sumup ( drho_g, vxc_g )
        END DO
        CALL pw_transfer ( vxc_g, vxc )
        CALL pw_deallocate ( vxc_g )
     END IF
     DEALLOCATE ( rgrid, STAT = stat )
     IF ( stat /= 0 ) CALL stop_memory ( routineP, "rgrid" )
     IF ( dft_control%gradient_functional ) THEN
        DEALLOCATE ( drgrid, STAT = stat )
        IF ( stat /= 0 ) CALL stop_memory ( routineP, "drgrid" )
     END IF
  END IF
  DEALLOCATE ( pot, STAT = stat )
  IF ( stat /= 0 ) CALL stop_memory ( routineP, "pot" )
  IF ( dft_control%gradient_functional ) THEN
     CALL pw_deallocate ( drho_g )
     DO ik = 1, 3
        CALL pw_deallocate ( drho_r(ik) )
     END DO
  END IF

END SUBROUTINE xc_calculate_pw_lda

!!*****
!FM ! remove??
!FMsubroutine b88fun(ex,v1,v2,rho,drho2)
!FM  implicit none
!FM  real(dbl) :: rho,drho2,ex,v1,v2
!FM  real(dbl) :: b1,alpha,f1,f43,rs,ob3
!FM  real(dbl) :: two13,aa,a,br1,br2,br4,xs,xs2,sa2b8,shm1,dd,dd2,ee,sx,v1x,v2x
!FM
!FM      IF(RHO.LE.1.e-10_dbl) THEN
!FM        ex = 0.0D0
!FM        v1 = 0.0D0
!FM        v2 = 0.0D0
!FM      ELSE
!FM        b1    = 0.0042_dbl
!FM        alpha = 2._dbl/3._dbl
!FM        f1    = -1.10783814957303361_dbl
!FM        f43   = 4._dbl/3._dbl
!FM
!FM        RS = RHO**(1._dbl/3._dbl)
!FM        ex = F1*ALPHA*RS
!FM        v1 = F43*F1*ALPHA*RS
!FM        v2 = 0._dbl
!FM
!FM        TWO13 = 2._dbl**(1._dbl/3._dbl)
!FM        AA    = drho2
!FM        A     = SQRT(AA)
!FM        BR1   = RS
!FM        BR2   = BR1*BR1
!FM        BR4   = BR2*BR2
!FM        XS    = TWO13*A/BR4
!FM        XS2   = XS*XS
!FM        SA2B8 = SQRT(1._dbl+XS2)
!FM        SHM1  = LOG(XS+SA2B8)
!FM        DD    = 1._dbl + 6._dbl*B1*XS*SHM1
!FM        DD2   = DD*DD
!FM        EE    = 6._dbl*B1*XS2/SA2B8 - 1._dbl
!FM        SX    = TWO13*AA/BR4*(-B1/DD)
!FM        V1X   = -(4._dbl/3._dbl)/TWO13*XS2*B1*BR1*EE/DD2
!FM        V2X   = TWO13*B1*(EE-DD)/(BR4*DD2)
!FM
!FM        ex    = rho*ex+sx
!FM        v1    = v1 + v1x
!FM        v2    = v2x
!FM      ENDIF
!FM
!FMend subroutine b88fun
!FM!******************************************************************************

END MODULE xc

!******************************************************************************
