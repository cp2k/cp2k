!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for encoding and decoding the output presented to the learning algorithm
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_ml_codec
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE pao_potentials,                  ONLY: pao_build_V_core,&
                                              pao_build_V_neighbors
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_ml_codec'

   PUBLIC :: pao_ml_encode_xblock, pao_ml_decode_xblock

CONTAINS

! **************************************************************************************************
!> \brief Uses U0 to encode the sign of rotation angles in a contiguous fashion.
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param u0block ...
!> \param xblock ...
!> \param output ...
! **************************************************************************************************
   SUBROUTINE pao_ml_encode_xblock(particle_set, qs_kind_set, cell, iatom, u0block, xblock, output)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: u0block, xblock
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: output

      INTEGER, PARAMETER                                 :: L = 3

      INTEGER                                            :: a, b, i, ikind, M, N, nparams
      REAL(dp)                                           :: Vab_core, Vab_neighbors
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: V_core, V_neighbors

      CPASSERT(.NOT. ALLOCATED(output))

      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), pao_basis_size=M)
      N = SIZE(u0block, 1) ! primary basis-size

      ! build potentials
      CALL pao_build_V_neighbors(particle_set, qs_kind_set, cell, iatom, V_neighbors)
      CALL pao_build_V_core(particle_set, qs_kind_set, iatom, V_core)

      nparams = SIZE(xblock)
      ALLOCATE (output(L*nparams))

      DO i = 1, nparams
         a = MOD(i-1, M)+1
         b = M+(i-1)/M+1

         Vab_neighbors = SUM(MATMUL(MATMUL(TRANSPOSE(u0block(:, a:a)), V_neighbors), u0block(:, b:b)))
         Vab_core = SUM(MATMUL(MATMUL(TRANSPOSE(u0block(:, a:a)), V_core), u0block(:, b:b)))

         output(L*(i-1)+1) = ABS(xblock(i, 1))
         output(L*(i-1)+2) = SIGN(1.0_dp, xblock(i, 1))*Vab_core
         output(L*(i-1)+3) = SIGN(1.0_dp, xblock(i, 1))*Vab_neighbors
      ENDDO

      DEALLOCATE (V_neighbors, V_core)

   END SUBROUTINE pao_ml_encode_xblock

! **************************************************************************************************
!> \brief Uses U0 to deencode the sign of rotation angles
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param iatom ...
!> \param u0block ...
!> \param output ...
!> \param xblock ...
! **************************************************************************************************
   SUBROUTINE pao_ml_decode_xblock(particle_set, qs_kind_set, cell, iatom, u0block, output, xblock)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: iatom
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: u0block
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: output
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: xblock

      INTEGER, PARAMETER                                 :: L = 3

      INTEGER                                            :: a, b, i, ikind, M, N, nparams
      REAL(dp)                                           :: sign_quenched, sign_vote, Vab_core, &
                                                            Vab_neighbors
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: V_core, V_neighbors

      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), pao_basis_size=M)
      N = SIZE(u0block, 1) ! primary basis-size

      ! build potentials
      CALL pao_build_V_neighbors(particle_set, qs_kind_set, cell, iatom, V_neighbors)
      CALL pao_build_V_core(particle_set, qs_kind_set, iatom, V_core)

      nparams = SIZE(xblock)

      xblock(:, :) = 0.0_dp
      DO i = 1, nparams
         a = MOD(i-1, M)+1
         b = M+(i-1)/M+1

         Vab_neighbors = SUM(MATMUL(MATMUL(TRANSPOSE(u0block(:, a:a)), V_neighbors), u0block(:, b:b)))
         Vab_core = SUM(MATMUL(MATMUL(TRANSPOSE(u0block(:, a:a)), V_core), u0block(:, b:b)))

         sign_vote = Vab_core*output(L*(i-1)+2)+Vab_neighbors*output(L*(i-1)+3)
         sign_quenched = TANH(1e5_dp*sign_vote) !TODO make scaling factor input param
         xblock(i, 1) = sign_quenched*output(L*(i-1)+1)

      ENDDO

      DEALLOCATE (V_neighbors, V_core)
   END SUBROUTINE pao_ml_decode_xblock

END MODULE pao_ml_codec
