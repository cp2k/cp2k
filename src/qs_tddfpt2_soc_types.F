!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE qs_tddfpt2_soc_types
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_release
   USE kinds,                           ONLY: dp
   USE libint_wrapper,                  ONLY: cp_libint_static_cleanup
   USE qs_grid_atom,                    ONLY: deallocate_grid_atom,&
                                              grid_atom_type
   USE qs_harmonics_atom,               ONLY: deallocate_harmonics_atom,&
                                              harmonics_atom_type

   IMPLICIT NONE

   PUBLIC  :: soc_env_type, soc_env_create, &
              soc_env_release, soc_atom_env_type, &
              soc_atom_create, soc_atom_release

   PRIVATE

   !! Helper copied from xas_tdp_types
   TYPE grid_atom_p_type
      TYPE(grid_atom_type), POINTER                   :: grid_atom => Null()
   END TYPE grid_atom_p_type

   TYPE harmonics_atom_p_type
      TYPE(harmonics_atom_type), POINTER              :: harmonics_atom => Null()
   END TYPE harmonics_atom_p_type

!*************************************************************************************************
! \brief: This structure contains the static matrices for the soc-correction
! \param dbcsr_soc: the ZORA-operator within the ao-basis
! \param dipmat: the dipole-operator within the ao-basis
! \paramn evals_a: spin-conservin/singlet excitation energies
! \param evals_b: spin-filp/triplet exciation energies
! \param a_coeff: spin-conservin/singlet excitation vector
! \param b_coeff: spin-filp/triplet exciation vectors
! \param soc_evals: SOC-Corrected eigenvalues
! \param soc_osc: ozillatorstrength of soc-corrected excitations
!*************************************************************************************************
   TYPE soc_env_type
      !! a :: singlet or spin-conserving b :: triplet or spin flip
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: orb_soc => Null()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: dipmat => Null()
      REAL(dp), POINTER, DIMENSION(:)                 :: evals_a => Null(), &
                                                         evals_b => Null()
      TYPE(cp_fm_type)                                :: a_coeff = cp_fm_type(), &
                                                         b_coeff = cp_fm_type()
      REAL(dp), ALLOCATABLE, DIMENSION(:)             :: soc_evals, &
                                                         soc_osc
   END TYPE

! ************************************************************************************************
! \bief: an environment type analog to the xas_atom_env in xas_tdp_types:343
!        All unused parameters have been droped
! \param nspins: number of spins
! \param grid_atom_set:
! \param harmonics_atom_set:
! \param orb_sphi_so contains the coefficient for direct contraction from so to sgf, for the orb basis
!*************************************************************************************************
   TYPE soc_atom_env_type
      INTEGER                                            :: nspins = -1
      TYPE(grid_atom_p_type), DIMENSION(:), POINTER      :: grid_atom_set => Null()
      TYPE(harmonics_atom_p_type), DIMENSION(:), POINTER :: harmonics_atom_set => Null()
      TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER        :: orb_sphi_so => Null()
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: soc_pp => Null()
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param soc_env ...
! **************************************************************************************************
   SUBROUTINE soc_env_create(soc_env)
      TYPE(soc_env_type), TARGET                         :: soc_env

      NULLIFY (soc_env%orb_soc)
      NULLIFY (soc_env%evals_a)
      NULLIFY (soc_env%evals_b)
      NULLIFY (soc_env%dipmat)

   END SUBROUTINE soc_env_create

! **************************************************************************************************
!> \brief ...
!> \param soc_env ...
! **************************************************************************************************
   SUBROUTINE soc_env_release(soc_env)
      TYPE(soc_env_type), TARGET                         :: soc_env

      INTEGER                                            :: i

      IF (ASSOCIATED(soc_env%orb_soc)) THEN
         DO i = 1, SIZE(soc_env%orb_soc)
            CALL dbcsr_release(soc_env%orb_soc(i)%matrix)
            DEALLOCATE (soc_env%orb_soc(i)%matrix)
         END DO
         DEALLOCATE (soc_env%orb_soc)
      END IF

      CALL cp_fm_release(soc_env%a_coeff)
      CALL cp_fm_release(soc_env%b_coeff)
      IF (ASSOCIATED(soc_env%evals_a)) NULLIFY (soc_env%evals_a)
      IF (ASSOCIATED(soc_env%evals_b)) NULLIFY (soc_env%evals_b)
      IF (ASSOCIATED(soc_env%dipmat)) THEN
         DO i = 1, SIZE(soc_env%dipmat)
            CALL dbcsr_release(soc_env%dipmat(i)%matrix)
            DEALLOCATE (soc_env%dipmat(i)%matrix)
         END DO
         DEALLOCATE (soc_env%dipmat)
      END IF
      IF (ALLOCATED(soc_env%soc_evals)) DEALLOCATE (soc_env%soc_evals)
      IF (ALLOCATED(soc_env%soc_osc)) DEALLOCATE (soc_env%soc_osc)

   END SUBROUTINE soc_env_release

! ************************************************************************************************
!*************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param soc_atom_env ...
! **************************************************************************************************
   SUBROUTINE soc_atom_create(soc_atom_env)
      TYPE(soc_atom_env_type), POINTER                   :: soc_atom_env

      ALLOCATE (soc_atom_env)

      soc_atom_env%nspins = 1
      NULLIFY (soc_atom_env%grid_atom_set)
      NULLIFY (soc_atom_env%harmonics_atom_set)
      NULLIFY (soc_atom_env%orb_sphi_so)
      NULLIFY (soc_atom_env%soc_pp)

   END SUBROUTINE soc_atom_create

! **************************************************************************************************
!> \brief ...
!> \param soc_atom_env ...
! **************************************************************************************************
   SUBROUTINE soc_atom_release(soc_atom_env)
      TYPE(soc_atom_env_type), POINTER                   :: soc_atom_env

      INTEGER                                            :: i

      IF (ASSOCIATED(soc_atom_env%grid_atom_set)) THEN
         DO i = 1, SIZE(soc_atom_env%grid_atom_set)
            IF (ASSOCIATED(soc_atom_env%grid_atom_set(i)%grid_atom)) THEN
               CALL deallocate_grid_atom(soc_atom_env%grid_atom_set(i)%grid_atom)
            END IF
         END DO
         DEALLOCATE (soc_atom_env%grid_atom_set)
      END IF

      IF (ASSOCIATED(soc_atom_env%harmonics_atom_set)) THEN
         DO i = 1, SIZE(soc_atom_env%harmonics_atom_set)
            IF (ASSOCIATED(soc_atom_env%harmonics_atom_set(i)%harmonics_atom)) THEN
               CALL deallocate_harmonics_atom(soc_atom_env%harmonics_atom_set(i)%harmonics_atom)
            END IF
         END DO
         DEALLOCATE (soc_atom_env%harmonics_atom_set)
      END IF

      IF (ASSOCIATED(soc_atom_env%orb_sphi_so)) THEN
         DO i = 1, SIZE(soc_atom_env%orb_sphi_so)
            IF (ASSOCIATED(soc_atom_env%orb_sphi_so(i)%array)) THEN
               DEALLOCATE (soc_atom_env%orb_sphi_so(i)%array)
            END IF
         END DO
         DEALLOCATE (soc_atom_env%orb_sphi_so)
      END IF

      IF (ASSOCIATED(soc_atom_env%soc_pp)) CALL dbcsr_deallocate_matrix_set(soc_atom_env%soc_pp)

      !Clean-up libint
      CALL cp_libint_static_cleanup()

      DEALLOCATE (soc_atom_env)

   END SUBROUTINE soc_atom_release

END MODULE qs_tddfpt2_soc_types
