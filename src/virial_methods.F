!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      JGH [04042007] code refactoring
! **************************************************************************************************
MODULE virial_methods

   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: det_3x3
   USE message_passing,                 ONLY: mp_sum
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC:: virial_evaluate, virial_pair_force, virial_update, write_stress_components

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'virial_methods'

CONTAINS
! **************************************************************************************************
!> \brief Updates the virial given the virial and subsys
!> \param virial ...
!> \param subsys ...
!> \param para_env ...
!> \par History
!>      none
!> \author Teodoro Laino [tlaino] - 03.2008 - Zurich University
! **************************************************************************************************
   SUBROUTINE virial_update(virial, subsys, para_env)
      TYPE(virial_type), INTENT(INOUT)                   :: virial
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_update', routineP = moduleN//':'//routineN

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(particle_list_type), POINTER                  :: particles

      CALL cp_subsys_get(subsys, local_particles=local_particles, atomic_kinds=atomic_kinds, &
                         particles=particles)

      CALL virial_evaluate(atomic_kinds%els, particles%els, local_particles, &
                           virial, para_env%group)

   END SUBROUTINE virial_update

! **************************************************************************************************
!> \brief Computes the kinetic part of the pressure tensor and updates
!>      the full VIRIAL (PV)
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param virial ...
!> \param igroup ...
!> \par History
!>      none
!> \author CJM
! **************************************************************************************************
   SUBROUTINE virial_evaluate(atomic_kind_set, particle_set, local_particles, &
                              virial, igroup)

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(virial_type), INTENT(INOUT)                   :: virial
      INTEGER, INTENT(IN)                                :: igroup

      CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_evaluate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iparticle, iparticle_kind, &
                                                            iparticle_local, j, nparticle_kind, &
                                                            nparticle_local
      REAL(KIND=dp)                                      :: mass
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind

      IF (virial%pv_availability) THEN
         CALL timeset(routineN, handle)
         NULLIFY (atomic_kind)
         nparticle_kind = SIZE(atomic_kind_set)
         virial%pv_kinetic = 0.0_dp
         DO i = 1, 3
            DO j = 1, i
               DO iparticle_kind = 1, nparticle_kind
                  atomic_kind => atomic_kind_set(iparticle_kind)
                  CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
                  nparticle_local = local_particles%n_el(iparticle_kind)
                  DO iparticle_local = 1, nparticle_local
                     iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                     virial%pv_kinetic(i, j) = virial%pv_kinetic(i, j) + &
                                               mass*particle_set(iparticle)%v(i)*particle_set(iparticle)%v(j)
                  END DO
               END DO
               virial%pv_kinetic(j, i) = virial%pv_kinetic(i, j)
            END DO
         END DO

         CALL mp_sum(virial%pv_kinetic, igroup)

         ! total virial
         virial%pv_total = virial%pv_virial + virial%pv_kinetic + virial%pv_constraint

         CALL timestop(handle)
      ENDIF

   END SUBROUTINE virial_evaluate

! **************************************************************************************************
!> \brief Computes the contribution to the stress tensor from two-body
!>      pair-wise forces
!> \param pv_virial ...
!> \param f0 ...
!> \param force ...
!> \param rab ...
!> \par History
!>      none
!> \author JGH
! **************************************************************************************************
   PURE SUBROUTINE virial_pair_force(pv_virial, f0, force, rab)

      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: pv_virial
      REAL(KIND=dp), INTENT(IN)                          :: f0
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: force, rab

      CHARACTER(LEN=*), PARAMETER :: routineN = 'virial_pair_force', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j

      DO i = 1, 3
         DO j = 1, 3
            pv_virial(i, j) = pv_virial(i, j) + f0*force(i)*rab(j)
         END DO
      END DO

   END SUBROUTINE virial_pair_force

! **************************************************************************************************
!> \brief ...
!> \param virial ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE write_stress_components(virial, iounit)
      TYPE(virial_type), INTENT(IN)                      :: virial
      INTEGER, INTENT(IN)                                :: iounit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_stress_components', &
         routineP = moduleN//':'//routineN

      REAL(KIND=dp), DIMENSION(3, 3)                     :: pv

      WRITE (UNIT=iounit, FMT="(/,T3,A)") REPEAT("=", 78)
      WRITE (UNIT=iounit, FMT="(T4,A)") "Stress Tensor Components (GPW/GAPW)"
      WRITE (UNIT=iounit, FMT="(T42,A,T64,A)") "       1/3 Trace", "     Determinant"
      pv = virial%pv_ekin
      WRITE (UNIT=iounit, FMT="(T4,A,T42,F16.8,T64,F16.8)") "Kinetic Energy Stress", &
         (pv(1, 1) + pv(2, 2) + pv(3, 3))/3.0_dp, det_3x3(pv)
      pv = virial%pv_overlap
      WRITE (UNIT=iounit, FMT="(T4,A,T42,F16.8,T64,F16.8)") "Basis Overlap Stress", &
         (pv(1, 1) + pv(2, 2) + pv(3, 3))/3.0_dp, det_3x3(pv)
      pv = virial%pv_hartree + virial%pv_exc
      WRITE (UNIT=iounit, FMT="(T4,A,T42,F16.8,T64,F16.8)") "ES + XC Stress", &
         (pv(1, 1) + pv(2, 2) + pv(3, 3))/3.0_dp, det_3x3(pv)
      pv = virial%pv_vdw
      WRITE (UNIT=iounit, FMT="(T4,A,T42,F16.8,T64,F16.8)") "vdW correction (ff) Stress", &
         (pv(1, 1) + pv(2, 2) + pv(3, 3))/3.0_dp, det_3x3(pv)
      pv = virial%pv_ppl
      WRITE (UNIT=iounit, FMT="(T4,A,T42,F16.8,T64,F16.8)") "Local Pseudopotential/Core Stress", &
         (pv(1, 1) + pv(2, 2) + pv(3, 3))/3.0_dp, det_3x3(pv)
      pv = virial%pv_ppnl
      WRITE (UNIT=iounit, FMT="(T4,A,T42,F16.8,T64,F16.8)") "Nonlocal Pseudopotential Stress", &
         (pv(1, 1) + pv(2, 2) + pv(3, 3))/3.0_dp, det_3x3(pv)
      pv = -virial%pv_fock_4c
      WRITE (UNIT=iounit, FMT="(T4,A,T42,F16.8,T64,F16.8)") "Exact Exchange Stress", &
         (pv(1, 1) + pv(2, 2) + pv(3, 3))/3.0_dp, det_3x3(pv)
      !
      pv = virial%pv_ekin + virial%pv_overlap + virial%pv_hartree + virial%pv_exc + virial%pv_vdw + &
           virial%pv_ppl + virial%pv_ppnl - virial%pv_fock_4c
      WRITE (UNIT=iounit, FMT="(T4,A,T42,F16.8,T64,F16.8)") "Sum of Parts Stress", &
         (pv(1, 1) + pv(2, 2) + pv(3, 3))/3.0_dp, det_3x3(pv)
      !
      pv = virial%pv_virial
      WRITE (UNIT=iounit, FMT="(T4,A,T42,F16.8,T64,F16.8)") "Total Stress", &
         (pv(1, 1) + pv(2, 2) + pv(3, 3))/3.0_dp, det_3x3(pv)
      WRITE (UNIT=iounit, FMT="(T3,A)") REPEAT("=", 78)

   END SUBROUTINE write_stress_components

END MODULE virial_methods

