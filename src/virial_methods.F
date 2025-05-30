!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      JGH [04042007] code refactoring
! **************************************************************************************************
MODULE virial_methods

   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathlib,                         ONLY: det_3x3,&
                                              jacobi
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC:: one_third_sum_diag, virial_evaluate, virial_pair_force, virial_update, &
            write_stress_tensor, write_stress_tensor_components

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
      TYPE(mp_para_env_type), POINTER                    :: para_env

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(particle_list_type), POINTER                  :: particles

      CALL cp_subsys_get(subsys, local_particles=local_particles, atomic_kinds=atomic_kinds, &
                         particles=particles)

      CALL virial_evaluate(atomic_kinds%els, particles%els, local_particles, &
                           virial, para_env)

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

      CLASS(mp_comm_type), INTENT(IN)                     :: igroup

      CHARACTER(LEN=*), PARAMETER                        :: routineN = "virial_evaluate"

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

         CALL igroup%sum(virial%pv_kinetic)

         ! total virial
         virial%pv_total = virial%pv_virial + virial%pv_kinetic + virial%pv_constraint

         CALL timestop(handle)
      END IF

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
!> \param iw ...
!> \param cell ...
!> \param unit_string ...
!> \par History
!>      - Revised virial components (14.10.2020, MK)
!> \author JGH
! **************************************************************************************************
   SUBROUTINE write_stress_tensor_components(virial, iw, cell, unit_string)

      TYPE(virial_type), INTENT(IN)                      :: virial
      INTEGER, INTENT(IN)                                :: iw
      TYPE(cell_type), POINTER                           :: cell
      CHARACTER(LEN=default_string_length), INTENT(IN)   :: unit_string

      CHARACTER(LEN=*), PARAMETER :: fmt = "(T2,A,T41,2(1X,ES19.11))"

      REAL(KIND=dp)                                      :: fconv
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pv

      IF (iw > 0) THEN
         CPASSERT(ASSOCIATED(cell))
         ! Conversion factor
         fconv = cp_unit_from_cp2k(1.0_dp/cell%deth, TRIM(ADJUSTL(unit_string)))
         WRITE (UNIT=iw, FMT="(/,T2,A)") &
            "STRESS| Stress tensor components (GPW/GAPW) ["//TRIM(ADJUSTL(unit_string))//"]"
         WRITE (UNIT=iw, FMT="(T2,A,T52,A,T70,A)") &
            "STRESS|", "1/3 Trace", "Determinant"
         pv = fconv*virial%pv_overlap
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Overlap", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_ekinetic
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Kinetic energy", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_ppl
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Local pseudopotential/core", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_ppnl
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Nonlocal pseudopotential", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_ecore_overlap
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Core charge overlap", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_ehartree
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Hartree", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_exc
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Exchange-correlation", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_exx
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Exact exchange (EXX)", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_vdw
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| vdW correction", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_mp2
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Moller-Plesset (MP2)", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_nlcc
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Nonlinear core correction (NLCC)", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_gapw
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Local atomic parts (GAPW)", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_lrigpw
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Resolution-of-the-identity (LRI)", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*(virial%pv_overlap + virial%pv_ekinetic + virial%pv_ppl + virial%pv_ppnl + &
                     virial%pv_ecore_overlap + virial%pv_ehartree + virial%pv_exc + &
                     virial%pv_exx + virial%pv_vdw + virial%pv_mp2 + virial%pv_nlcc + &
                     virial%pv_gapw + virial%pv_lrigpw)
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Sum of components", one_third_sum_diag(pv), det_3x3(pv)
         pv = fconv*virial%pv_virial
         WRITE (UNIT=iw, FMT=fmt) &
            "STRESS| Total", one_third_sum_diag(pv), det_3x3(pv)
      END IF

   END SUBROUTINE write_stress_tensor_components

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION one_third_sum_diag(a) RESULT(p)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: a
      REAL(KIND=dp)                                      :: p

      p = (a(1, 1) + a(2, 2) + a(3, 3))/3.0_dp
   END FUNCTION one_third_sum_diag

! **************************************************************************************************
!> \brief Print stress tensor to output file
!>
!> \param pv_virial ...
!> \param iw ...
!> \param cell ...
!> \param unit_string ...
!> \param numerical ...
!> \author MK (26.08.2010)
! **************************************************************************************************
   SUBROUTINE write_stress_tensor(pv_virial, iw, cell, unit_string, numerical)

      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: pv_virial
      INTEGER, INTENT(IN)                                :: iw
      TYPE(cell_type), POINTER                           :: cell
      CHARACTER(LEN=default_string_length), INTENT(IN)   :: unit_string
      LOGICAL, INTENT(IN)                                :: numerical

      REAL(KIND=dp)                                      :: fconv
      REAL(KIND=dp), DIMENSION(3)                        :: eigval
      REAL(KIND=dp), DIMENSION(3, 3)                     :: eigvec, stress_tensor

      IF (iw > 0) THEN
         CPASSERT(ASSOCIATED(cell))
         ! Conversion factor
         fconv = cp_unit_from_cp2k(1.0_dp/cell%deth, TRIM(ADJUSTL(unit_string)))
         stress_tensor(:, :) = fconv*pv_virial(:, :)
         IF (numerical) THEN
            WRITE (UNIT=iw, FMT="(/,T2,A)") &
               "STRESS| Numerical stress tensor ["//TRIM(ADJUSTL(unit_string))//"]"
         ELSE
            WRITE (UNIT=iw, FMT="(/,T2,A)") &
               "STRESS| Analytical stress tensor ["//TRIM(ADJUSTL(unit_string))//"]"
         END IF
         WRITE (UNIT=iw, FMT="(T2,A,T14,3(19X,A1))") &
            "STRESS|", "x", "y", "z"
         WRITE (UNIT=iw, FMT="(T2,A,T21,3(1X,ES19.11))") &
            "STRESS|      x", stress_tensor(1, 1:3)
         WRITE (UNIT=iw, FMT="(T2,A,T21,3(1X,ES19.11))") &
            "STRESS|      y", stress_tensor(2, 1:3)
         WRITE (UNIT=iw, FMT="(T2,A,T21,3(1X,ES19.11))") &
            "STRESS|      z", stress_tensor(3, 1:3)
         WRITE (UNIT=iw, FMT="(T2,A,T61,ES20.11)") &
            "STRESS| 1/3 Trace", (stress_tensor(1, 1) + &
                                  stress_tensor(2, 2) + &
                                  stress_tensor(3, 3))/3.0_dp
         WRITE (UNIT=iw, FMT="(T2,A,T61,ES20.11)") &
            "STRESS| Determinant", det_3x3(stress_tensor(1:3, 1), &
                                           stress_tensor(1:3, 2), &
                                           stress_tensor(1:3, 3))
         eigval(:) = 0.0_dp
         eigvec(:, :) = 0.0_dp
         CALL jacobi(stress_tensor, eigval, eigvec)
         IF (numerical) THEN
            WRITE (UNIT=iw, FMT="(/,T2,A)") &
               "STRESS| Eigenvectors and eigenvalues of the numerical stress tensor ["// &
               TRIM(ADJUSTL(unit_string))//"]"
         ELSE
            WRITE (UNIT=iw, FMT="(/,T2,A)") &
               "STRESS| Eigenvectors and eigenvalues of the analytical stress tensor ["// &
               TRIM(ADJUSTL(unit_string))//"]"
         END IF
         WRITE (UNIT=iw, FMT="(T2,A,T14,3(1X,I19))") &
            "STRESS|", 1, 2, 3
         WRITE (UNIT=iw, FMT="(T2,A,T21,3(1X,ES19.11))") &
            "STRESS| Eigenvalues", eigval(1:3)
         WRITE (UNIT=iw, FMT="(T2,A,T21,3(1X,F19.12))") &
            "STRESS|      x", eigvec(1, 1:3)
         WRITE (UNIT=iw, FMT="(T2,A,T21,3(1X,F19.12))") &
            "STRESS|      y", eigvec(2, 1:3)
         WRITE (UNIT=iw, FMT="(T2,A,T21,3(1X,F19.12))") &
            "STRESS|      z", eigvec(3, 1:3)
      END IF

   END SUBROUTINE write_stress_tensor

END MODULE virial_methods
