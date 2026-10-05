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
      ! <SEMILOCAL>
      REAL(dp), DIMENSION(:, :), ALLOCATABLE         :: vsemi
   END TYPE atom_upfpot_type

   PUBLIC :: atom_upfpot_type, upf_partial_waves_type, upf_radial_weights
   PUBLIC :: upf_release, upf_has_complete_partial_waves
CONTAINS

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
      IF (.NOT. ALLOCATED(pot%partial_waves)) RETURN
      IF (.NOT. ALLOCATED(pot%lbeta)) RETURN
      IF (.NOT. ALLOCATED(pot%r) .OR. .NOT. ALLOCATED(pot%rab)) RETURN
      n = pot%mesh_size
      nb = pot%number_of_proj
      IF (n < 3 .OR. nb < 1 .OR. SIZE(pot%lbeta) /= nb) RETURN
      IF (SIZE(pot%r) /= n .OR. SIZE(pot%rab) /= n) RETURN
      ASSOCIATE (waves => pot%partial_waves)
         IF (.NOT. ALLOCATED(waves%ae_wfc) .OR. .NOT. ALLOCATED(waves%ps_wfc)) RETURN
         IF (.NOT. ALLOCATED(waves%angular_momenta)) RETURN
         IF (ANY(SHAPE(waves%ae_wfc) /= [n, nb])) RETURN
         IF (ANY(SHAPE(waves%ps_wfc) /= [n, nb])) RETURN
         IF (SIZE(waves%angular_momenta) /= nb) RETURN
         complete = ALL(waves%angular_momenta == pot%lbeta)
      END ASSOCIATE
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
