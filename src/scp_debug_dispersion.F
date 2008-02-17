!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Perform the debug calculation of the module scp_dispersion 
!> \author CJM
!******************************************************************************

MODULE scp_debug_dispersion
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE scp_dispersion,                  ONLY: scp_qs_dispersion
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_debug_dispersion'
  PUBLIC :: scp_qs_disp_debug


CONTAINS

! *****************************************************************************
!> \brief debug calculation of the scp_qs_dispersion 
!> \author CJM
! *****************************************************************************


SUBROUTINE scp_qs_disp_debug ( qs_env,  error )

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE scp_qs_debug (MODULE scp_debug)", &
      routineN = 'scp_qs_disp_debug', routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, idim, natoms
    REAL(dp)                                 :: em, ep
    REAL(dp), ALLOCATABLE                    :: f( :, : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(scp_energy_type), POINTER           :: energy
    TYPE(scp_environment_type), POINTER      :: scp_env

  NULLIFY(particle_set, scp_env, energy)

  WRITE ( *, * ) '******************************START DEBUG***********************************'
!-------------
   CALL get_qs_env(qs_env=qs_env, &
                   scp_env = scp_env, &
                   particle_set=particle_set,&
                   error=error)
   CALL get_scp_env ( scp_env, energy = energy, error = error )
   natoms = SIZE ( particle_set )
   ALLOCATE ( f ( natoms, 3 ) )
   CALL scp_qs_dispersion ( qs_env, .FALSE., error )
   WRITE ( *, * ) 'DISPERSION_ENERGY', energy % e_scp_disp
   DO iatom = 1, natoms
     DO idim = 1, 3
       particle_set ( iatom ) % r ( idim ) = particle_set ( iatom ) % r ( idim ) + 1.e-5_dp
       CALL scp_qs_dispersion ( qs_env, .FALSE., error )
       ep = energy % e_scp_disp 
       particle_set ( iatom ) % r ( idim ) = particle_set ( iatom ) % r ( idim ) - 2.e-5_dp
       CALL scp_qs_dispersion ( qs_env, .FALSE., error )
       em = energy % e_scp_disp 
       particle_set ( iatom ) % r ( idim ) = particle_set ( iatom ) % r ( idim ) + 1.e-5_dp
       f ( iatom, idim ) = ( em - ep )/2._dp/1.e-5_dp
     END DO 
     WRITE ( *, * ) 'FDEBUG_DISPERSION', iatom, f ( iatom, : )
   END DO
   DEALLOCATE ( f )
  WRITE ( *, * ) '**********************************END DEBUG*********************************'

END SUBROUTINE scp_qs_disp_debug

!******************************************************************************
END MODULE scp_debug_dispersion

