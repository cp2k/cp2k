!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Owned UPF dataset and radial conventions shared by ATOM and Quickstep.
!>        Dataset availability does not select an orbital representation or calculation method.
! **************************************************************************************************
MODULE upf_types
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! use same value as in atom_types!
   INTEGER, PARAMETER                                :: lmat = 3

   TYPE upf_partial_waves_type
      ! Original r*phi(r), one AE/PS pair per beta projector; not the atomic guesses.
      REAL(dp), ALLOCATABLE :: ae_wfc(:, :), ps_wfc(:, :)
      INTEGER, ALLOCATABLE :: angular_momenta(:)
   END TYPE upf_partial_waves_type

   TYPE upf_gipaw_orbitals_type
      ! Signed u=r*R(r), without a frozen-core assignment or assumed occupations.
      INTEGER :: data_format = -1
      REAL(dp), ALLOCATABLE :: wfc(:, :)
      INTEGER, ALLOCATABLE :: indices(:), principal_quantum_numbers(:), angular_momenta(:)
   END TYPE upf_gipaw_orbitals_type

   TYPE, EXTENDS(upf_partial_waves_type) :: upf_paw_type
      INTEGER :: data_format = -1, cutoff_index = -1, lmax_aug = -1
      CHARACTER(LEN=default_string_length) :: augmentation_shape = ""
      REAL(dp) :: cutoff_radius = -1.0_dp, augmentation_epsilon = 1.e-12_dp
      REAL(dp) :: core_energy = 0.0_dp
      LOGICAL :: has_core_energy = .FALSE.
      REAL(dp), ALLOCATABLE :: occupations(:), multipoles(:, :, :)
      ! AE core density per volume; AE local potential and core energy in Hartree.
      REAL(dp), ALLOCATABLE :: rho_core(:), vlocal(:)
      TYPE(upf_gipaw_orbitals_type), ALLOCATABLE :: gipaw_core_orbitals
   END TYPE upf_paw_type

   TYPE atom_upfpot_type
      CHARACTER(LEN=2)                               :: symbol = ""
      CHARACTER(LEN=default_string_length)           :: pname = ""
      INTEGER, DIMENSION(0:lmat)                     :: econf = 0
      REAL(dp)                                       :: zion = 0.0_dp
      CHARACTER(LEN=default_string_length)           :: version = ""
      CHARACTER(LEN=default_path_length)             :: filename = ""
      ! <INFO>
      INTEGER                                        :: maxinfo = 100
      CHARACTER(LEN=default_string_length), DIMENSION(100) &
         :: info = ""
      ! <HEADER>
      CHARACTER(LEN=default_string_length)           :: generated = ""
      CHARACTER(LEN=default_string_length)           :: author = ""
      CHARACTER(LEN=default_string_length)           :: date = ""
      CHARACTER(LEN=default_string_length)           :: comment = ""
      CHARACTER(LEN=4)                               :: pseudo_type = ""
      CHARACTER(LEN=15)                              :: relativistic = ""
      CHARACTER(LEN=default_string_length)           :: functional = ""
      LOGICAL                                        :: is_ultrasoft = .FALSE.
      LOGICAL                                        :: is_paw = .FALSE.
      LOGICAL                                        :: is_coulomb = .FALSE.
      LOGICAL                                        :: has_so = .FALSE.
      LOGICAL                                        :: has_wfc = .FALSE.
      LOGICAL                                        :: has_gipaw = .FALSE.
      LOGICAL                                        :: paw_as_gipaw = .FALSE.
      LOGICAL                                        :: core_correction = .FALSE.
      LOGICAL                                        :: with_metagga_info = .FALSE.
      REAL(dp)                                       :: total_psenergy = 0.0_dp
      REAL(dp)                                       :: wfc_cutoff = 0.0_dp
      REAL(dp)                                       :: rho_cutoff = 0.0_dp
      INTEGER                                        :: l_max = -100
      INTEGER                                        :: l_max_rho = -1
      INTEGER                                        :: l_local = -1
      INTEGER                                        :: mesh_size = -1
      INTEGER                                        :: number_of_wfc = -1
      INTEGER                                        :: number_of_proj = -1
      ! <MESH>
      REAL(dp)                                       :: dx = 0.0_dp
      REAL(dp)                                       :: xmin = 0.0_dp
      REAL(dp)                                       :: rmax = 0.0_dp
      REAL(dp)                                       :: zmesh = 0.0_dp
      REAL(dp), DIMENSION(:), ALLOCATABLE            :: r, rab
      ! <NLCC>
      REAL(dp), DIMENSION(:), ALLOCATABLE            :: rho_nlcc
      ! Preserve supplied samples. Only an explicit canonical convention permits XC use.
      REAL(dp), ALLOCATABLE :: tau_nlcc(:), tau_atom(:)
      CHARACTER(LEN=default_string_length) :: tau_nlcc_convention = "", tau_atom_convention = ""
      ! <LOCAL>
      REAL(dp), DIMENSION(:), ALLOCATABLE            :: vlocal
      ! <NONLOCAL>: keep the tabulated projector amplitude; D includes Ry -> Hartree.
      ! beta stores beta(r), with the radial prefactor in PP_BETA removed.
      REAL(dp), DIMENSION(:, :), ALLOCATABLE         :: dion
      REAL(dp), DIMENSION(:, :), ALLOCATABLE         :: beta
      INTEGER, DIMENSION(:), ALLOCATABLE             :: lbeta
      ! <AUGMENTATION>: PP_Q and r**2*Q_ij^L(r), in the original beta convention.
      ! The pair index is j*(j-1)/2+i for i<=j; the last dimension starts at L=0.
      REAL(dp), DIMENSION(:, :), ALLOCATABLE         :: qij
      REAL(dp), DIMENSION(:, :, :), ALLOCATABLE      :: qijl
      TYPE(upf_partial_waves_type), ALLOCATABLE :: partial_waves
      TYPE(upf_paw_type), ALLOCATABLE :: paw
      ! <SEMILOCAL>
      REAL(dp), DIMENSION(:, :), ALLOCATABLE         :: vsemi
   END TYPE atom_upfpot_type

   PUBLIC :: atom_upfpot_type, upf_partial_waves_type, upf_radial_weights
   PUBLIC :: upf_paw_type, upf_gipaw_orbitals_type
   PUBLIC :: upf_release, upf_has_complete_partial_waves
   PUBLIC :: upf_nlcc_kinetic, upf_uses_model_core_tau
CONTAINS

! **************************************************************************************************
!> \brief Select model-core tau for USPP and for NC with an explicitly supplied table.
!> \param pot Original UPF metadata and optional model-core table.
!> \return Whether this dataset uses the model-core convention instead of PAW reconstruction.
! **************************************************************************************************
   PURE FUNCTION upf_uses_model_core_tau(pot) RESULT(use_tau)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      LOGICAL                                            :: use_tau

      use_tau = .NOT. pot%is_paw .AND. (pot%is_ultrasoft .OR. ALLOCATED(pot%tau_nlcc))
   END FUNCTION upf_uses_model_core_tau

! **************************************************************************************************
!> \brief Obtain the independently supplied model-core kinetic density in Hartree per cubic bohr.
!>        Charge augmentation and NLCC charge do not determine this field. Unlabelled samples
!>        remain available as raw data, without an assumed generator normalization.
!> \param pot Original UPF dataset.
!> \param kinetic Spin-summed canonical core tau, zero for an absent partial core.
!> \param available Whether the canonical field is defined for this dataset.
!> \param reason Explanation when the field cannot be used.
! **************************************************************************************************
   SUBROUTINE upf_nlcc_kinetic(pot, kinetic, available, reason)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: kinetic(:)
      LOGICAL, INTENT(OUT)                               :: available
      CHARACTER(LEN=*), INTENT(OUT)                      :: reason

      LOGICAL                                            :: has_core

      CPASSERT(ALLOCATED(pot%r))
      ALLOCATE (kinetic(SIZE(pot%r)), SOURCE=0.0_dp)
      available = .FALSE.
      reason = "UPF model-core tau is missing. Supply canonical PP_TAUMOD or prepare an explicit model with add_upf_core_tau.py"
      has_core = .FALSE.
      IF (ALLOCATED(pot%rho_nlcc)) has_core = ANY(pot%rho_nlcc /= 0.0_dp)
      IF (.NOT. has_core) THEN
         IF (ALLOCATED(pot%tau_nlcc)) THEN
            IF (ANY(pot%tau_nlcc /= 0.0_dp)) THEN
               reason = "UPF has a nonzero model-core kinetic density without a nonzero NLCC density"
               RETURN
            END IF
         END IF
         available = .TRUE.
         reason = ""
         RETURN
      END IF
      IF (.NOT. pot%core_correction) THEN
         reason = "UPF has a nonzero NLCC density with core_correction disabled"
         RETURN
      END IF
      IF (.NOT. ALLOCATED(pot%tau_nlcc)) RETURN
      IF (TRIM(pot%tau_nlcc_convention) /= "HARTREE_PER_VOLUME") THEN
         reason = "UPF PP_TAUMOD requires an explicit hartree_per_volume convention before use in XC"
         RETURN
      END IF
      CPASSERT(SIZE(pot%tau_nlcc) == SIZE(pot%r))
      CPASSERT(ALL(ieee_is_finite(pot%tau_nlcc)))
      kinetic(:) = pot%tau_nlcc
      available = .TRUE.
      reason = ""
   END SUBROUTINE upf_nlcc_kinetic

! **************************************************************************************************
!> \brief Composite Simpson weights for dr on the supplied UPF mesh index.
!> \param pot Original radial coordinates and dr/d(index).
!> \return Radial weights, with a last trapezoidal interval for an even mesh size.
! **************************************************************************************************
   FUNCTION upf_radial_weights(pot) RESULT(weights)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      REAL(dp)                                           :: weights(SIZE(pot%r))

      INTEGER                                            :: n, ns

      n = SIZE(pot%r)
      CPASSERT(n >= 3 .AND. SIZE(pot%rab) == n)
      CPASSERT(ALL(pot%r(2:) > pot%r(:n - 1)) .AND. ALL(pot%rab > 0.0_dp))
      ns = n - 1 + MODULO(n, 2)
      weights = 0.0_dp
      weights(1:ns:2) = 2.0_dp/3.0_dp
      weights(2:ns - 1:2) = 4.0_dp/3.0_dp
      weights(1) = 1.0_dp/3.0_dp
      weights(ns) = 1.0_dp/3.0_dp
      IF (ns < n) weights(n - 1:n) = weights(n - 1:n) + 0.5_dp
      weights(:) = weights*pot%rab
   END FUNCTION upf_radial_weights

! **************************************************************************************************
!> \brief Check for a complete projector-indexed pair of original AE and PS waves.
!> \param pot Original UPF dataset; atomic guess waves do not satisfy this query.
!> \return Whether every beta projector has a compatible full-wave pair.
! **************************************************************************************************
   PURE FUNCTION upf_has_complete_partial_waves(pot) RESULT(complete)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      LOGICAL                                            :: complete

      INTEGER                                            :: n, nb

      complete = .FALSE.
      IF (.NOT. ALLOCATED(pot%lbeta)) RETURN
      IF (.NOT. ALLOCATED(pot%r) .OR. .NOT. ALLOCATED(pot%rab)) RETURN
      n = pot%mesh_size
      nb = pot%number_of_proj
      IF (n < 3 .OR. nb < 1 .OR. SIZE(pot%lbeta) /= nb) RETURN
      IF (SIZE(pot%r) /= n .OR. SIZE(pot%rab) /= n) RETURN
      IF (pot%is_paw) THEN
         IF (.NOT. ALLOCATED(pot%paw)) RETURN
         complete = compatible_waves(pot%paw)
      ELSE
         IF (.NOT. ALLOCATED(pot%partial_waves)) RETURN
         complete = compatible_waves(pot%partial_waves)
      END IF
   CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param waves ...
!> \return ...
! **************************************************************************************************
      PURE FUNCTION compatible_waves(waves) RESULT(compatible)
         CLASS(upf_partial_waves_type), INTENT(IN) :: waves
         LOGICAL :: compatible

         compatible = .FALSE.
         IF (.NOT. ALLOCATED(waves%ae_wfc) .OR. .NOT. ALLOCATED(waves%ps_wfc)) RETURN
         IF (.NOT. ALLOCATED(waves%angular_momenta)) RETURN
         IF (ANY(SHAPE(waves%ae_wfc) /= [n, nb])) RETURN
         IF (ANY(SHAPE(waves%ps_wfc) /= [n, nb])) RETURN
         IF (SIZE(waves%angular_momenta) /= nb) RETURN
         compatible = ALL(waves%angular_momenta == pot%lbeta)
      END FUNCTION compatible_waves
   END FUNCTION upf_has_complete_partial_waves

! **************************************************************************************************
!> \brief Release every owned field and reset metadata before reusing a dataset.
!> \param pot Owned UPF dataset.
! **************************************************************************************************
   PURE SUBROUTINE upf_release(pot)
      TYPE(atom_upfpot_type), INTENT(INOUT)              :: pot

      pot = atom_upfpot_type()
   END SUBROUTINE upf_release
END MODULE upf_types
