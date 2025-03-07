!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Definition and initialisation of the ps_wavelet data type.
!> \author Florian Schiffmann (09.2007,fschiff)
! **************************************************************************************************
MODULE ps_wavelet_types

   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_wavelet_types'

   INTEGER, PARAMETER, PUBLIC               :: WAVELET3D = 1400, &
                                               WAVELET2D = 1401, &
                                               WAVELET1D = 1402, &
                                               WAVELET0D = 1403

   PUBLIC :: ps_wavelet_type, &
             ps_wavelet_release

! **************************************************************************************************
!> \par History
!>      09.2007 created [Florian Schiffmann]
!> \author fschiff
! **************************************************************************************************
   TYPE ps_wavelet_type
      CHARACTER(LEN=1)                                  :: geocode = " "
      CHARACTER(LEN=1)                                  :: datacode = " "
      INTEGER                                           :: itype_scf = 0
      INTEGER                                           :: method = WAVELET3D, special_dimension = 0
      REAL(kind=dp), POINTER, DIMENSION(:)             :: karray => NULL()
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER    :: rho_z_sliced => NULL()
      INTEGER, DIMENSION(3)                              :: PS_grid = 0
   END TYPE ps_wavelet_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param wavelet ...
! **************************************************************************************************
   SUBROUTINE ps_wavelet_release(wavelet)

      TYPE(ps_wavelet_type), POINTER                     :: wavelet

      IF (ASSOCIATED(wavelet)) THEN
         IF (ASSOCIATED(wavelet%karray)) &
            DEALLOCATE (wavelet%karray)
         IF (ASSOCIATED(wavelet%rho_z_sliced)) &
            DEALLOCATE (wavelet%rho_z_sliced)
         DEALLOCATE (wavelet)
      END IF
   END SUBROUTINE ps_wavelet_release

END MODULE ps_wavelet_types
