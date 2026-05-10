!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Local FCI solver interface for active-space embedding.
! **************************************************************************************************
MODULE qs_active_space_fci
#if defined(__LIBFCI)
   USE ISO_C_BINDING, ONLY: C_CHAR, &
                            C_DOUBLE, &
                            C_INT, &
                            C_NULL_CHAR
#endif
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_para_env_type
   USE qs_active_space_types, ONLY: active_space_type
#if defined(__LIBFCI)
   USE qs_active_space_utils, ONLY: eri_to_array, &
                                    subspace_matrix_to_array
#endif
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_active_space_fci'

   PUBLIC :: solve_active_space_fci

#if defined(__LIBFCI)
   INTERFACE
      FUNCTION libfci_solve(nspins, norb, nelec, ms2, h1_a, h1_b, eri_aa, eri_ab, eri_bb, &
                            max_iter, threshold, pspace, project_spin, verbose, energy, &
                            rdm1_a, rdm1_b, err_msg, err_msg_len) RESULT(status) &
         BIND(C, name="libfci_solve")
         IMPORT :: C_CHAR, C_DOUBLE, C_INT
         INTEGER(C_INT), VALUE                         :: nspins, norb, nelec, ms2
         REAL(C_DOUBLE), DIMENSION(*), INTENT(IN)      :: h1_a, h1_b, eri_aa, eri_ab, eri_bb
         INTEGER(C_INT), VALUE                         :: max_iter, pspace, project_spin, verbose
         REAL(C_DOUBLE), VALUE                         :: threshold
         REAL(C_DOUBLE), INTENT(OUT)                   :: energy
         REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT)     :: rdm1_a, rdm1_b
         CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(OUT) :: err_msg
         INTEGER(C_INT), VALUE                         :: err_msg_len
         INTEGER(C_INT)                                :: status
      END FUNCTION libfci_solve
   END INTERFACE
#endif

CONTAINS

! **************************************************************************************************
!> \brief Solve the current active-space Hamiltonian with the local FCI kernel.
!> \param active_space_env active-space environment
!> \param para_env parallel environment
!> \param p_act_mo_a alpha or spin-summed active-space 1-RDM
!> \param p_act_mo_b beta active-space 1-RDM for unrestricted calculations
! **************************************************************************************************
   SUBROUTINE solve_active_space_fci(active_space_env, para_env, p_act_mo_a, p_act_mo_b)
      TYPE(active_space_type), INTENT(INOUT), POINTER    :: active_space_env
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: p_act_mo_a
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT), OPTIONAL                           :: p_act_mo_b

#if defined(__LIBFCI)
      CHARACTER(KIND=C_CHAR), DIMENSION(512)             :: c_error
      CHARACTER(LEN=512)                                 :: error_text
      INTEGER                                            :: max_iter, ms2, n2, n4, nmo_active, &
                                                            nspins, status
      LOGICAL                                            :: ionode
      REAL(KIND=dp)                                      :: energy_active, threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eri_aa, eri_ab, eri_bb, fock_a, fock_b, &
                                                            p_beta

      nmo_active = active_space_env%nmo_active
      nspins = active_space_env%nspins
      n2 = nmo_active*nmo_active
      n4 = n2*n2
      ionode = para_env%is_source()

      ALLOCATE (fock_a(n2), eri_aa(n4), p_act_mo_a(n2))
      ALLOCATE (fock_b(MAX(1, n2)), eri_ab(MAX(1, n4)), eri_bb(MAX(1, n4)), p_beta(MAX(1, n2)))
      fock_b(:) = 0.0_dp
      eri_ab(:) = 0.0_dp
      eri_bb(:) = 0.0_dp
      p_beta(:) = 0.0_dp

      ASSOCIATE (act_indices => active_space_env%active_orbitals(:, 1))
         CALL subspace_matrix_to_array(active_space_env%fock_sub(1), fock_a, act_indices, act_indices)
      END ASSOCIATE
      CALL eri_to_array(active_space_env%eri, eri_aa, active_space_env%active_orbitals, 1, 1)

      IF (nspins == 2) THEN
         IF (.NOT. PRESENT(p_act_mo_b)) THEN
            CPABORT("Missing beta output buffer for local active-space FCI solver.")
         END IF
         ALLOCATE (p_act_mo_b(n2))
         ASSOCIATE (act_indices => active_space_env%active_orbitals(:, 2))
            CALL subspace_matrix_to_array(active_space_env%fock_sub(2), fock_b(1:n2), act_indices, act_indices)
         END ASSOCIATE
         CALL eri_to_array(active_space_env%eri, eri_ab(1:n4), active_space_env%active_orbitals, 1, 2)
         CALL eri_to_array(active_space_env%eri, eri_bb(1:n4), active_space_env%active_orbitals, 2, 2)
      END IF

      status = 0
      energy_active = 0.0_dp
      error_text = ""
      c_error(:) = C_NULL_CHAR
      ms2 = active_space_env%multiplicity - 1
      max_iter = 2048
      threshold = 1.0E-8_dp

      IF (ionode) THEN
         status = libfci_solve(INT(nspins, C_INT), INT(nmo_active, C_INT), &
                               INT(active_space_env%nelec_active, C_INT), INT(ms2, C_INT), &
                               fock_a, fock_b, eri_aa, eri_ab, eri_bb, INT(max_iter, C_INT), &
                               REAL(threshold, C_DOUBLE), INT(200, C_INT), INT(1, C_INT), INT(0, C_INT), &
                               energy_active, p_act_mo_a, p_beta, c_error, INT(SIZE(c_error), C_INT))
         error_text = c_string(c_error)
      END IF

      CALL para_env%bcast(status, para_env%source)
      CALL para_env%bcast(error_text, para_env%source)
      IF (status /= 0) THEN
         IF (LEN_TRIM(error_text) > 0) THEN
            CPABORT("Local active-space FCI solver failed: "//TRIM(error_text))
         ELSE
            CPABORT("Local active-space FCI solver failed.")
         END IF
      END IF

      CALL para_env%bcast(energy_active, para_env%source)
      CALL para_env%bcast(p_act_mo_a, para_env%source)
      active_space_env%energy_active = energy_active
      IF (nspins == 2) THEN
         p_act_mo_b(:) = p_beta(1:n2)
         CALL para_env%bcast(p_act_mo_b, para_env%source)
      END IF

      DEALLOCATE (fock_a, fock_b, eri_aa, eri_ab, eri_bb, p_beta)
#else
      MARK_USED(active_space_env)
      MARK_USED(para_env)
      MARK_USED(p_act_mo_a)
      IF (PRESENT(p_act_mo_b)) THEN
         MARK_USED(p_act_mo_b)
      END IF
      CPABORT("AS_SOLVER FCI requires CP2K to be built with LibFCI support.")
#endif

   END SUBROUTINE solve_active_space_fci

#if defined(__LIBFCI)
! **************************************************************************************************
!> \brief Convert a C null-terminated string buffer to a Fortran string.
!> \param buffer ...
!> \return ...
! **************************************************************************************************
   FUNCTION c_string(buffer) RESULT(text)
      CHARACTER(KIND=C_CHAR), DIMENSION(:), INTENT(IN)   :: buffer
      CHARACTER(LEN=512)                                 :: text

      INTEGER                                            :: i

      text = ""
      DO i = 1, MIN(SIZE(buffer), LEN(text))
         IF (buffer(i) == C_NULL_CHAR) EXIT
         text(i:i) = buffer(i)
      END DO
   END FUNCTION c_string
#endif

END MODULE qs_active_space_fci
