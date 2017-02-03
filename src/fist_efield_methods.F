!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!> \author JGH
! **************************************************************************************************
MODULE fist_efield_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_result_methods,               ONLY: cp_results_erase,&
                                              put_results
   USE cp_result_types,                 ONLY: cp_result_type
   USE fist_efield_types,               ONLY: fist_efield_type
   USE fist_environment_types,          ONLY: fist_env_get,&
                                              fist_environment_type
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: twopi
   USE moments_utils,                   ONLY: get_reference_point
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: debye
#include "./base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_efield_methods'

   PRIVATE

   PUBLIC :: fist_dipole, fist_efield_energy_force

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qenergy ...
!> \param qforce ...
!> \param qpv ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param cell ...
!> \param efield ...
!> \param use_virial ...
!> \param iunit ...
!> \param charges ...
! **************************************************************************************************
   SUBROUTINE fist_efield_energy_force(qenergy, qforce, qpv, atomic_kind_set, particle_set, cell, &
                                       efield, use_virial, iunit, charges)
      REAL(KIND=dp), INTENT(OUT)                         :: qenergy
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: qforce
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: qpv
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(fist_efield_type), POINTER                    :: efield
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_virial
      INTEGER, INTENT(IN), OPTIONAL                      :: iunit
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: charges

      COMPLEX(KIND=dp)                                   :: zeta
      COMPLEX(KIND=dp), DIMENSION(3)                     :: ggamma
      INTEGER                                            :: i, ii, iparticle_kind, iw, j
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: use_charges, virial
      REAL(KIND=dp)                                      :: q, theta
      REAL(KIND=dp), DIMENSION(3)                        :: ci, dfilter, di, dipole, fieldpol, fq, &
                                                            gvec, ria, tmp
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind

      qenergy = 0.0_dp
      qforce = 0.0_dp
      qpv = 0.0_dp

      use_charges = .FALSE.
      IF (PRESENT(charges)) THEN
         IF (ASSOCIATED(charges)) use_charges = .TRUE.
      END IF

      IF (PRESENT(iunit)) THEN
         iw = iunit
      ELSE
         iw = -1
      END IF

      IF (PRESENT(use_virial)) THEN
         virial = use_virial
      ELSE
         virial = .FALSE.
      END IF

      fieldpol = efield%polarisation
      fieldpol = fieldpol/SQRT(DOT_PRODUCT(fieldpol, fieldpol))
      fieldpol = -fieldpol*efield%strength

      dfilter = efield%dfilter

      dipole = 0.0_dp
      ggamma = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      DO iparticle_kind = 1, SIZE(atomic_kind_set)
         atomic_kind => atomic_kind_set(iparticle_kind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=q, atom_list=atom_list)
         ! TODO parallelization over atoms (local_particles)
         DO i = 1, SIZE(atom_list)
            ii = atom_list(i)
            ria = particle_set(ii)%r(:)
            ria = pbc(ria, cell)
            IF (use_charges) q = charges(ii)
            DO j = 1, 3
               gvec = twopi*cell%h_inv(j, :)
               theta = SUM(ria(:)*gvec(:))
               zeta = CMPLX(COS(theta), SIN(theta), KIND=dp)**(-q)
               ggamma(j) = ggamma(j)*zeta
            END DO
            qforce(1:3, ii) = q
         ENDDO
      END DO

      IF (ALL(REAL(ggamma, KIND=dp) /= 0.0_dp)) THEN
         tmp = AIMAG(ggamma)/REAL(ggamma, KIND=dp)
         ci = ATAN(tmp)
         dipole = MATMUL(cell%hmat, ci)/twopi
      END IF

      IF (efield%displacement) THEN
         ! E = (omega/8Pi)(D - 4Pi*P)^2
         di = dipole/cell%deth
         DO i = 1, 3
            theta = fieldpol(i)-2._dp*twopi*di(i)
            qenergy = qenergy+dfilter(i)*theta**2
            fq(i) = -dfilter(i)*theta
         END DO
         qenergy = 0.25_dp*cell%deth/twopi*qenergy
         DO i = 1, SIZE(qforce, 2)
            qforce(1:3, i) = fq(i)*qforce(1:3, i)
         END DO
      ELSE
         ! E = -omega*E*P
         qenergy = SUM(fieldpol*dipole)
         DO i = 1, SIZE(qforce, 2)
            qforce(1:3, i) = fieldpol(1:3)*qforce(1:3, i)
         END DO
      END IF

      IF (virial) THEN
         DO iparticle_kind = 1, SIZE(atomic_kind_set)
            atomic_kind => atomic_kind_set(iparticle_kind)
            CALL get_atomic_kind(atomic_kind=atomic_kind, atom_list=atom_list)
            DO i = 1, SIZE(atom_list)
               ii = atom_list(i)
               ria = particle_set(ii)%r(:)
               ria = pbc(ria, cell)
               DO j = 1, 3
                  qpv(j, 1:3) = qpv(j, 1:3)+qforce(j, ii)*ria(1:3)
               END DO
            END DO
         END DO
         ! Stress tensor for constant D needs further investigation
         IF (efield%displacement) THEN
            CPABORT("Stress Tensor for constant D simulation is not working")
         END IF
      END IF

   END SUBROUTINE fist_efield_energy_force
! **************************************************************************************************
!> \brief Evaluates the Dipole of a classical charge distribution(point-like)
!>      possibly using the berry phase formalism
!> \param fist_env ...
!> \param print_section ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param cell ...
!> \param unit_nr ...
!> \param charges ...
!> \par History
!>      [01.2006] created
!>      [12.2007] tlaino - University of Zurich - debug and extended
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE fist_dipole(fist_env, print_section, atomic_kind_set, particle_set, &
                          cell, unit_nr, charges)
      TYPE(fist_environment_type), POINTER               :: fist_env
      TYPE(section_vals_type), POINTER                   :: print_section
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: charges

      CHARACTER(LEN=default_string_length)               :: description, dipole_type
      COMPLEX(KIND=dp)                                   :: dzeta, dzphase(3), zeta, zphase(3)
      COMPLEX(KIND=dp), DIMENSION(3)                     :: dggamma, ggamma
      INTEGER                                            :: i, iparticle_kind, j, reference
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: do_berry, use_charges
      REAL(KIND=dp) :: charge_tot, ci(3), dci(3), dipole(3), dipole_deriv(3), drcc(3), dria(3), &
         dtheta, gvec(3), q, rcc(3), ria(3), theta, tmp(3), via(3)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cp_result_type), POINTER                      :: results

      NULLIFY (atomic_kind)
      ! Reference point
      reference = section_get_ival(print_section, keyword_name="DIPOLE%REFERENCE")
      NULLIFY (ref_point)
      description = '[DIPOLE]'
      CALL section_vals_val_get(print_section, "DIPOLE%REF_POINT", r_vals=ref_point)
      CALL section_vals_val_get(print_section, "DIPOLE%PERIODIC", l_val=do_berry)
      use_charges = .FALSE.
      IF (PRESENT(charges)) THEN
         IF (ASSOCIATED(charges)) use_charges = .TRUE.
      END IF

      CALL get_reference_point(rcc, drcc, fist_env=fist_env, reference=reference, ref_point=ref_point)

      ! Dipole deriv will be the derivative of the Dipole(dM/dt=\sum e_j v_j)
      dipole_deriv = 0.0_dp
      dipole = 0.0_dp
      IF (do_berry) THEN
         dipole_type = "[BERRY PHASE]"
         rcc = pbc(rcc, cell)
         charge_tot = 0._dp
         IF (use_charges) THEN
            charge_tot = SUM(charges)
         ELSE
            DO i = 1, SIZE(particle_set)
               atomic_kind => particle_set(i)%atomic_kind
               CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=q)
               charge_tot = charge_tot+q
            END DO
         END IF
         ria = twopi*MATMUL(cell%h_inv, rcc)
         zphase = CMPLX(COS(ria), SIN(ria), dp)**charge_tot

         dria = twopi*MATMUL(cell%h_inv, drcc)
         dzphase = charge_tot*CMPLX(-SIN(ria), COS(ria), dp)**(charge_tot-1.0_dp)*dria

         ggamma = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
         dggamma = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO iparticle_kind = 1, SIZE(atomic_kind_set)
            atomic_kind => atomic_kind_set(iparticle_kind)
            CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=q, atom_list=atom_list)

            DO i = 1, SIZE(atom_list)
               ria = particle_set(atom_list(i))%r(:)
               ria = pbc(ria, cell)
               via = particle_set(atom_list(i))%v(:)
               IF (use_charges) q = charges(atom_list(i))
               DO j = 1, 3
                  gvec = twopi*cell%h_inv(j, :)
                  theta = SUM(ria(:)*gvec(:))
                  dtheta = SUM(via(:)*gvec(:))
                  zeta = CMPLX(COS(theta), SIN(theta), KIND=dp)**(-q)
                  dzeta = -q*CMPLX(-SIN(theta), COS(theta), KIND=dp)**(-q-1.0_dp)*dtheta
                  dggamma(j) = dggamma(j)*zeta+ggamma(j)*dzeta
                  ggamma(j) = ggamma(j)*zeta
               END DO
            ENDDO
         END DO
         dggamma = dggamma*zphase+ggamma*dzphase
         ggamma = ggamma*zphase
         IF (ALL(REAL(ggamma, KIND=dp) /= 0.0_dp)) THEN
            tmp = AIMAG(ggamma)/REAL(ggamma, KIND=dp)
            ci = ATAN(tmp)
            dci = (1.0_dp/(1.0_dp+tmp**2))* &
                  (AIMAG(dggamma)*REAL(ggamma, KIND=dp)-AIMAG(ggamma)*REAL(dggamma, KIND=dp))/(REAL(ggamma, KIND=dp))**2

            dipole = MATMUL(cell%hmat, ci)/twopi
            dipole_deriv = MATMUL(cell%hmat, dci)/twopi
         END IF
         CALL fist_env_get(fist_env=fist_env, results=results)
         CALL cp_results_erase(results, description)
         CALL put_results(results, description, dipole)
      ELSE
         dipole_type = "[Non Periodic]"
         DO i = 1, SIZE(particle_set)
            atomic_kind => particle_set(i)%atomic_kind
            ria = particle_set(i)%r(:) ! no pbc(particle_set(i)%r(:),cell) so that the total dipole
            ! is the sum of the molecular dipoles
            CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=q)
            IF (use_charges) q = charges(i)
            dipole = dipole-q*(ria-rcc)
            dipole_deriv(:) = dipole_deriv(:)-q*(particle_set(i)%v(:)-drcc)
         END DO
         CALL fist_env_get(fist_env=fist_env, results=results)
         CALL cp_results_erase(results, description)
         CALL put_results(results, description, dipole)
      END IF
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(1X,A,T48,3F11.6)') "MM DIPOLE "//TRIM(dipole_type)//"(A.U.)|", dipole
         WRITE (unit_nr, '(1X,A,T48,3F11.6)') "MM DIPOLE "//TRIM(dipole_type)//"(Debye)|", dipole*debye
         WRITE (unit_nr, '(1X,A,T48,3F11.6)') "MM DIPOLE "//TRIM(dipole_type)//" DERIVATIVE(A.U.)|", dipole_deriv
      END IF

   END SUBROUTINE fist_dipole

END MODULE fist_efield_methods
