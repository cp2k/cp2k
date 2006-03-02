
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****m* cp2k/qs_linres_nmr_shift *
!!
!!   NAME 
!!     qs_linres_nmr_current
!!
!!   FUNCTION
!!     from the response current density calculates the shift tensor
!!     and the susceptibility
!!     
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     created 02-2006 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_nmr_shift

  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE coefficient_types,               ONLY: coeff_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: Eval_Interp_Spl3_pbc
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE,&
                                             pw_transfer,&
                                             pw_type
  USE qs_linres_types,                 ONLY: nmr_env_type
#include "cp_common_uses.h"

  IMPLICIT NONE


  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: interpolate_shift_pwgrid, mult_G_ov_G2_grid

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_nmr_shift'

!!***
! **************

CONTAINS

!!****f* qs_linres_nmr_shift%interpolate_shift_pwgrid
!!
!!  NAME 
!!      nmr_interpolate_shift_pwgrid
!!
!!  FUNCTION
!!      interpolate the shift calculated on the PW grid in order to ger 
!!      the value on arbitrary points in real space
!!
!!  ARGUMENTS
!!      nme_env : to get the shift tensor and the list of additional points 
!!      particle_set : for the atomic position
!!      cell : to take into account the pbs, and to have the volume
!!      shift_pw_rspace : specific component of the shift tensor on the pw grid
!!      i_B : component of the magnetic field for which the shift is calculated (row)
!!      idir : component of the vector \int_{r}[ ((r-r') x j(r))/|r-r'|^3 ] = Bind(r')
!!
!!  NOTES
!!
!!  AUTHOR
!!    MI
!!
!!*** **********************************************************************

  SUBROUTINE interpolate_shift_pwgrid(nmr_env,particle_set,cell,shift_pw_rspace,&
             i_B,idir,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: shift_pw_rspace
    INTEGER, INTENT(IN)                      :: i_B, idir
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'interpolate_shift_pwgrid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, natom
    LOGICAL                                  :: failure
    REAL(dp)                                 :: R_iatom(3), ra(3), shift_val

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(shift_pw_rspace%pw),cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN
      natom = SIZE(particle_set,1)

      DO iatom = 1,natom
        R_iatom = pbc(particle_set(iatom)%r,cell)
        shift_val = Eval_Interp_Spl3_pbc(R_iatom,shift_pw_rspace%pw,error)
        nmr_env%chemical_shift(i_B,idir,iatom)= nmr_env%chemical_shift(i_B,idir,iatom)+&
            nmr_env%shift_factor*shift_val
      END DO 

      IF(nmr_env%do_nics) THEN

        DO iatom = 1,nmr_env%n_nics
          ra(1:3) = nmr_env%r_nics(1:3,iatom)
          R_iatom = pbc(ra,cell)
          shift_val = Eval_Interp_Spl3_pbc(R_iatom,shift_pw_rspace%pw,error)
          nmr_env%chemical_shift_nics(i_B,idir,iatom) = &
             nmr_env%chemical_shift_nics(i_B,idir,iatom) + &
             nmr_env%shift_factor*shift_val
        END DO
      END IF

    END IF

  END SUBROUTINE interpolate_shift_pwgrid

!!****f* qs_linres_nmr_shift%mult_G_ov_G2_grid
!!
!!  NAME 
!!      mult_G_ov_G2_grid
!!
!!  FUNCTION
!!      Given the current density on the PW grid in reciprcal space
!!      (obtained by FFT), calculate the integral  
!!        \int_{r}[ ((r-r') x j(r))/|r-r'|^3 ] = Bind(r')
!!      which in reciprcal space reads  (for G/=0)
!!         G/|G|^2 x J(G)    
!!
!!  ARGUMENTS
!!
!!  NOTES
!!     The G=0 component is not comnputed here, but can be evaluated
!!     through the susceptibility and added to the shift in a second time
!!
!!  AUTHOR
!!    MI
!!
!!*** **********************************************************************


  SUBROUTINE  mult_G_ov_G2_grid(cell,pw_pool,rho_gspace,funcG_times_rho,idir,error)


    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(coeff_type), POINTER                :: rho_gspace, funcG_times_rho
    INTEGER, INTENT(IN)                      :: idir
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mult_G_ov_G2_grid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ig, ng
    LOGICAL                                  :: failure
    REAL(dp)                                 :: g2
    TYPE(pw_grid_type), POINTER              :: grid
    TYPE(pw_type), POINTER                   :: frho, influence_fn

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)

    CALL pw_pool_create_pw ( pw_pool, influence_fn,&
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE )

    grid => influence_fn % pw_grid

    DO ig = grid % first_gne0, grid % ngpts_cut_local
       g2 = grid % gsq ( ig )
       influence_fn%cc(ig) = grid % g(idir,ig)/g2
    END DO  ! ig
    IF ( grid % have_g0 )  influence_fn%cc ( 1 ) = 0.0_dp

    frho => funcG_times_rho%pw
    CALL pw_transfer (rho_gspace%pw,frho)
    
    ng = SIZE(grid % gsq)
    frho%cc(1:ng) = frho%cc(1:ng)*influence_fn % cc  ( 1 : ng )

    CALL pw_pool_give_back_pw(pw_pool,influence_fn,&
         accept_non_compatible=.TRUE.,error=error)
 
  END SUBROUTINE  mult_G_ov_G2_grid

END MODULE qs_linres_nmr_shift



