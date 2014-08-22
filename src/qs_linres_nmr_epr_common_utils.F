!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief given the response wavefunctions obtained by the application
!>      of the (rxp), p, and ((dk-dl)xp) operators,
!>      here the current density vector (jx, jy, jz)
!>      is computed for the 3 directions of the magnetic field (Bx, By, Bz)
!> \par History
!>      created 02-2006 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_linres_nmr_epr_common_utils
  USE cell_types,                      ONLY: cell_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: gaussi
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_transfer
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: mult_G_ov_G2_grid

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_nmr_epr_common_utils'

CONTAINS

! *****************************************************************************
!> \brief Given the current density on the PW grid in reciprcal space
!>       (obtained by FFT), calculate the integral
!>         \int_{r}[ ((r-r') x j(r))/|r-r'|^3 ] = Bind(r') 
!>       which in reciprcal space reads  (for G/=0) 
!>          i G/|G|^2 x J(G)
!> \param cell ...
!> \param pw_pool ...
!> \param rho_gspace ...
!> \param funcG_times_rho ...
!> \param idir ...
!> \param my_chi ...
!> \param error ...
!> \author MI
!> \note
!>      The G=0 component is not comnputed here, but can be evaluated
!>      through the susceptibility and added to the shift in a second time
!>
!>      This method would not work for a non periodic system
!>      It should be generalized like the calculation of Hartree
! *****************************************************************************
  SUBROUTINE  mult_G_ov_G2_grid(cell,pw_pool,rho_gspace,funcG_times_rho,idir,my_chi,error)

    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_p_type), POINTER                 :: rho_gspace
    TYPE(pw_p_type)                          :: funcG_times_rho
    INTEGER, INTENT(IN)                      :: idir
    REAL(dp), INTENT(IN)                     :: my_chi
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: handle, ig, ng
    LOGICAL                                  :: failure
    REAL(dp)                                 :: g2
    TYPE(pw_grid_type), POINTER              :: grid
    CHARACTER(len=*), PARAMETER :: routineN = 'mult_G_ov_G2_grid', &
      routineP = moduleN//':'//routineN

    TYPE(pw_type), POINTER                   :: frho, influence_fn

    failure = .FALSE.
    CALL timeset(routineN,handle)

    CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)

    CALL pw_pool_create_pw ( pw_pool, influence_fn,&
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE ,error=error)

    grid => influence_fn % pw_grid
    DO ig = grid % first_gne0, grid % ngpts_cut_local
       g2 = grid % gsq ( ig )
       influence_fn%cc(ig) = gaussi * grid % g(idir,ig)/g2
    END DO  ! ig
    IF ( grid % have_g0 )  influence_fn%cc ( 1 ) = 0.0_dp

    frho => funcG_times_rho%pw
    CALL pw_transfer (rho_gspace%pw,frho,error=error)

    ng = SIZE(grid % gsq)
    frho%cc(1:ng) = frho%cc(1:ng)*influence_fn % cc  ( 1 : ng )
    IF ( grid % have_g0 ) frho%cc(1) = my_chi

    CALL pw_pool_give_back_pw(pw_pool,influence_fn,&
         accept_non_compatible=.TRUE.,error=error)

    CALL timestop(handle)

  END SUBROUTINE  mult_G_ov_G2_grid

END MODULE qs_linres_nmr_epr_common_utils
