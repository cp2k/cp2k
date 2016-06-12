!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for encoding and decoding the output presented to the learning algorithm
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_ml_codec
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_ml_codec'

   PUBLIC :: pao_ml_encode_xblock, pao_ml_decode_xblock

CONTAINS

! **************************************************************************************************
!> \brief Encodes block_X
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param xblock ...
!> \param output ...
! **************************************************************************************************
   SUBROUTINE pao_ml_encode_xblock(particle_set, qs_kind_set, cell, iatom, xblock, output)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: xblock
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: output

      INTEGER                                            :: nparams

      MARK_USED(particle_set)
      MARK_USED(qs_kind_set)
      MARK_USED(cell)
      MARK_USED(iatom)

      CPASSERT(.NOT. ALLOCATED(output))
      nparams = SIZE(xblock)
      ALLOCATE (output(nparams))
      output(:) = xblock(:, 1)

   END SUBROUTINE pao_ml_encode_xblock

! **************************************************************************************************
!> \brief Decodes block_X
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param output ...
!> \param xblock ...
! **************************************************************************************************
   SUBROUTINE pao_ml_decode_xblock(particle_set, qs_kind_set, cell, iatom, output, xblock)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: output
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: xblock

      MARK_USED(particle_set)
      MARK_USED(qs_kind_set)
      MARK_USED(cell)
      MARK_USED(iatom)

      xblock(:, 1) = output(:)

   END SUBROUTINE pao_ml_decode_xblock

END MODULE pao_ml_codec
