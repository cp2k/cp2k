!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/virial_methods [1.0] *
!!
!!   NAME
!!     virial_methods
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE virial_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC:: virial_evaluate

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'virial_methods'

!!*****
!******************************************************************************

CONTAINS

!!*****
!******************************************************************************
!!****** virial_methods/virial_evaluate [1.0] *
!!
!!   NAME
!!     virial_evaluate
!!
!!   SYNOPSIS
!!     Subroutine virial_evaluate( atomic_kind_set, particle_set,
!!                                 local_particle, virial )
!!       Implicit None
!!       Type(atomic_kind_type), Dimension(:), Intent (IN):: atomic_kind_set
!!       Type(particle_type), Dimension(:), Intent (IN):: particle_set
!!       Type(distribution_1d_type), Intent (IN):: local_particles
!!       Type(virial_type), Intent (INOUT):: virial
!!     End Subroutine virial_evaluate
!!
!!   FUNCTION
!!     Computes the kinetic part of the pressure tensor and updates
!!     the full VIRIAL (PV)
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE virial_evaluate ( atomic_kind_set, particle_set, &
                             local_particles, virial, igroup  )


    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(virial_type), INTENT(INOUT)         :: virial
    INTEGER, INTENT(IN)                      :: igroup

    INTEGER :: handle, i, iparticle, iparticle_kind, iparticle_local, j, &
      nnodes, nparticle_kind, nparticle_local
    REAL(KIND=dp)                            :: mass, mfl
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!
!------------------------------------------------------------------------------

   IF ( virial % pv_availability ) THEN
     CALL timeset ( 'PRESSURE', 'E', 'Mflops', handle )
     mfl = 0.0_dp
     nparticle_kind = SIZE ( atomic_kind_set )
     virial % pv_kinetic = 0.0_dp
     DO i = 1, 3
        DO j = 1, 3
           nnodes = 0
           DO iparticle_kind=1,nparticle_kind
              atomic_kind => atomic_kind_set(iparticle_kind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
              nparticle_local = local_particles%n_el(iparticle_kind)
              DO iparticle_local=1,nparticle_local
                nnodes = nnodes + 1
                iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                virial % pv_kinetic ( i, j ) = virial % pv_kinetic ( i, j ) + &
                                  mass * particle_set ( iparticle ) % v ( i )* &
                                         particle_set ( iparticle ) % v ( j )
              END DO
           END DO
        END DO
     END DO
     mfl = REAL( 9 * nnodes,KIND=dp) * 2.0_dp * 1.e-6_dp

     CALL mp_sum(virial%pv_kinetic,igroup)

! total virial
     virial%pv_total = virial%pv_virial + virial%pv_kinetic + virial%pv_constraint

     CALL timestop ( mfl, handle )
   ENDIF

END SUBROUTINE virial_evaluate
!!*****
END MODULE virial_methods
