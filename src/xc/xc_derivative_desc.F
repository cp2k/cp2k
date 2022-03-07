!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Module with functions to handle derivative descriptors.
!>      derivative description are strings have the following form
!>      "rhorhorhoa" which means that it is a forth order
!>      derivative, twice with respect to rho, once with respect to rhoa
!>      and once with respect to drhoa.
!>      Possible derivatives are:
!>      - rho: total density
!>      - norm_drho: norm of the gradient of the total density
!>      - rhoa, rhob: alpha and beta spin density (with LSD)
!>      - norm_drhoa, norm_drhob: norm of the gradient of the alpha and beta
!>        spin density
!>      - drhoa_drhob: scalar product between the gradient of rhoa and the
!>        grandient of rhob
!>      - tau: the local kinetic part
!>      - taua, taub: the kinetic part of the different spins
!> \note
!>       add drhox, drhoy, drhoz, drhoax,...?
!> \author thomas & fawzi
! **************************************************************************************************
MODULE xc_derivative_desc

   USE util,                            ONLY: sort
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: &
      deriv_rho = 1, &
      deriv_rhoa = 2, &
      deriv_rhob = 3, &
      deriv_norm_drho = 4, &
      deriv_norm_drhoa = 5, &
      deriv_norm_drhob = 6, &
      deriv_tau = 7, &
      deriv_tau_a = 8, &
      deriv_tau_b = 9, &
      deriv_laplace_rho = 10, &
      deriv_laplace_rhoa = 11, &
      deriv_laplace_rhob = 12

   INTEGER, PARAMETER :: MAX_LABEL_LENGTH = 12

   LOGICAL, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_derivative_desc'

   PUBLIC :: desc_to_id, id_to_desc, create_split_desc, standardize_desc

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param desc ...
!> \return ...
! **************************************************************************************************
   FUNCTION desc_to_id(desc) RESULT(id)
      CHARACTER(LEN=*), INTENT(IN)                       :: desc
      INTEGER                                            :: id

      SELECT CASE (TRIM(desc))
      CASE ("rho")
         id = deriv_rho
      CASE ("rhoa")
         id = deriv_rhoa
      CASE ("rhob")
         id = deriv_rhob
      CASE ("norm_drho")
         id = deriv_norm_drho
      CASE ("norm_drhoa")
         id = deriv_norm_drhoa
      CASE ("norm_drhob")
         id = deriv_norm_drhob
      CASE ("tau")
         id = deriv_tau
      CASE ("tau_a")
         id = deriv_tau_a
      CASE ("tau_b")
         id = deriv_tau_b
      CASE ("laplace_rho")
         id = deriv_laplace_rho
      CASE ("laplace_rhoa")
         id = deriv_laplace_rhoa
      CASE ("laplace_rhob")
         id = deriv_laplace_rhob
      CASE DEFAULT
         CPABORT("Unknown derivative variable: "//desc)
      END SELECT

   END FUNCTION desc_to_id

! **************************************************************************************************
!> \brief ...
!> \param id ...
!> \return ...
! **************************************************************************************************
   FUNCTION id_to_desc(id) RESULT(desc)
      INTEGER, INTENT(IN)                                :: id
      CHARACTER(LEN=MAX_LABEL_LENGTH)                    :: desc

      SELECT CASE (id)
      CASE (deriv_rho)
         desc = "rho"
      CASE (deriv_rhoa)
         desc = "rhoa"
      CASE (deriv_rhob)
         desc = "rhob"
      CASE (deriv_norm_drho)
         desc = "norm_drho"
      CASE (deriv_norm_drhoa)
         desc = "norm_drhoa"
      CASE (deriv_norm_drhob)
         desc = "norm_drhob"
      CASE (deriv_tau)
         desc = "tau"
      CASE (deriv_tau_a)
         desc = "tau_a"
      CASE (deriv_tau_b)
         desc = "tau_b"
      CASE (deriv_laplace_rho)
         desc = "laplace_rho"
      CASE (deriv_laplace_rhoa)
         desc = "laplace_rhoa"
      CASE (deriv_laplace_rhob)
         desc = "laplace_rhob"
      CASE DEFAULT
         CPABORT("Unknown derivative id!")
      END SELECT

   END FUNCTION id_to_desc

! **************************************************************************************************
!> \brief ...
!> \param desc ...
!> \param split_desc ...
! **************************************************************************************************
   SUBROUTINE create_split_desc(desc, split_desc)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: desc
      INTEGER, DIMENSION(:), POINTER                     :: split_desc

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: indices

      ALLOCATE (split_desc(SIZE(desc)))
      IF (SIZE(desc) > 0) THEN
         ALLOCATE (indices(SIZE(desc)))
         split_desc = desc
         CALL sort(split_desc, SIZE(desc), indices)
         DEALLOCATE (indices)
      END IF

   END SUBROUTINE create_split_desc

! **************************************************************************************************
!> \brief ...
!> \param desc ...
!> \param split_desc ...
! **************************************************************************************************
   SUBROUTINE standardize_desc(desc, split_desc)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: desc
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: split_desc

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: indices

      ALLOCATE (split_desc(SIZE(desc)))
      IF (SIZE(desc) > 0) THEN
         ALLOCATE (indices(SIZE(desc)))
         split_desc(:) = desc
         CALL sort(split_desc, SIZE(desc), indices)
         DEALLOCATE (indices)
      END IF

   END SUBROUTINE standardize_desc

END MODULE xc_derivative_desc
