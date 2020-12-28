!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contraction of integrals over primitive Cartesian Gaussians based on the contraction
!>        matrix sphi which is part of the gto_basis_set_type
!> \par History
!>      -added libxsmm_abc_contract routine, A. Bussy (04.2020)
!> \author Dorothea Golze (05.2016)
! **************************************************************************************************
MODULE ai_contraction_sphi

   USE kinds, ONLY: dp
#if(__LIBXSMM)
   USE libxsmm, ONLY: LIBXSMM_PREFETCH_NONE, &
                      libxsmm_blasint_kind, &
                      libxsmm_dgemm, &
                      libxsmm_dispatch, &
                      libxsmm_dmmcall, &
                      libxsmm_dmmfunction
#endif

#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_contraction_sphi'

   PUBLIC :: ab_contract, abc_contract, abcd_contract, libxsmm_abc_contract

CONTAINS

! **************************************************************************************************
!> \brief contract overlap integrals (a,b) and transfer to spherical Gaussians
!> \param abint contracted, normalized integrals of spherical Gaussians
!> \param sab uncontracted, unnormalized integrals of primitive Cartesian Gaussians
!> \param sphi_a contraction matrix for center a
!> \param sphi_b contraction matrix for center b
!> \param ncoa number of cartesian orbitals on a
!> \param ncob number of cartesian orbitals on b
!> \param nsgfa number of spherical Gaussian functions on a
!> \param nsgfb number of spherical Gaussian functions on b
! **************************************************************************************************
   SUBROUTINE ab_contract(abint, sab, sphi_a, sphi_b, ncoa, ncob, nsgfa, nsgfb)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: abint
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sab, sphi_a, sphi_b
      INTEGER, INTENT(IN)                                :: ncoa, ncob, nsgfa, nsgfb

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ab_contract'

      INTEGER                                            :: handle, m1, m2, msphia, msphib, nn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cpp

      CALL timeset(routineN, handle)

      msphia = SIZE(sphi_a, 1)
      msphib = SIZE(sphi_b, 1)

      m1 = SIZE(sab, 1)
      m2 = SIZE(sab, 2)

      nn = SIZE(abint, 1)

      ALLOCATE (cpp(nsgfa, m2))

      CALL dgemm("T", "N", nsgfa, m2, ncoa, 1._dp, sphi_a, msphia, sab, m1, 0.0_dp, cpp, nsgfa)
      CALL dgemm("N", "N", nsgfa, nsgfb, ncob, 1._dp, cpp, nsgfa, sphi_b, msphib, 0.0_dp, &
                 abint, nn)

      DEALLOCATE (cpp)

      CALL timestop(handle)

   END SUBROUTINE ab_contract

! **************************************************************************************************
!> \brief contract three-center overlap integrals (a,b,c) and transfer
!>        to spherical Gaussians
!> \param abcint contracted, normalized integrals of spherical Gaussians
!> \param sabc uncontracted, unnormalized integrals of primitive Cartesian Gaussians
!> \param sphi_a contraction matrix for center a
!> \param sphi_b contraction matrix for center b
!> \param sphi_c contraction matrix for center c
!> \param ncoa number of cartesian orbitals on a
!> \param ncob number of cartesian orbitals on b
!> \param ncoc number of cartesian orbitals on c
!> \param nsgfa number of spherical Gaussian functions on a
!> \param nsgfb number of spherical Gaussian functions on b
!> \param nsgfc number of spherical Gaussian functions on c
! **************************************************************************************************
   SUBROUTINE abc_contract(abcint, sabc, sphi_a, sphi_b, sphi_c, ncoa, ncob, ncoc, &
                           nsgfa, nsgfb, nsgfc)

      REAL(KIND=dp), DIMENSION(:, :, :)                  :: abcint, sabc
      REAL(KIND=dp), DIMENSION(:, :)                     :: sphi_a, sphi_b, sphi_c
      INTEGER, INTENT(IN)                                :: ncoa, ncob, ncoc, nsgfa, nsgfb, nsgfc

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'abc_contract'

      INTEGER                                            :: handle, i, m1, m2, m3, msphia, msphib, &
                                                            msphic
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: cpc, cpp

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(abcint, 1) == nsgfa)
      CPASSERT(SIZE(abcint, 2) == nsgfb)

      msphia = SIZE(sphi_a, 1)
      msphib = SIZE(sphi_b, 1)
      msphic = SIZE(sphi_c, 1)

      m1 = SIZE(sabc, 1)
      m2 = SIZE(sabc, 2)
      m3 = SIZE(sabc, 3)

      ALLOCATE (cpp(nsgfa, m2, m3), cpc(nsgfa, nsgfb, m3))

      CALL dgemm("T", "N", nsgfa, m2*m3, ncoa, 1._dp, sphi_a, msphia, sabc, m1, 0.0_dp, cpp, nsgfa)
      DO i = 1, m3
         CALL dgemm("N", "N", nsgfa, nsgfb, ncob, 1._dp, cpp(:, :, i), nsgfa, sphi_b, msphib, &
                    0.0_dp, cpc(:, :, i), nsgfa)
      END DO
      CALL dgemm("N", "N", nsgfa*nsgfb, nsgfc, ncoc, 1._dp, cpc, nsgfa*nsgfb, sphi_c, msphic, 0.0_dp, &
                 abcint, nsgfa*nsgfb)

      DEALLOCATE (cpp, cpc)

      CALL timestop(handle)

   END SUBROUTINE abc_contract

! **************************************************************************************************
!> \brief contract four-center overlap integrals (a,b,c,d) and transfer
!>        to spherical Gaussians
!> \param abcdint contracted, normalized integrals of spherical Gaussians
!> \param sabcd uncontracted, unnormalized integrals of primitive Cartesian Gaussians
!> \param sphi_a contraction matrix for center a
!> \param sphi_b contraction matrix for center b
!> \param sphi_c contraction matrix for center c
!> \param sphi_d contraction matrix for center d
!> \param ncoa number of cartesian orbitals on a
!> \param ncob number of cartesian orbitals on b
!> \param ncoc number of cartesian orbitals on c
!> \param ncod number of cartesian orbitals on d
!> \param nsgfa number of spherical Gaussian functions on a
!> \param nsgfb number of spherical Gaussian functions on b
!> \param nsgfc number of spherical Gaussian functions on c
!> \param nsgfd number of spherical Gaussian functions on d
! **************************************************************************************************
   SUBROUTINE abcd_contract(abcdint, sabcd, sphi_a, sphi_b, sphi_c, sphi_d, ncoa, ncob, &
                            ncoc, ncod, nsgfa, nsgfb, nsgfc, nsgfd)

      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT)                                   :: abcdint
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: sabcd
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sphi_a, sphi_b, sphi_c, sphi_d
      INTEGER, INTENT(IN)                                :: ncoa, ncob, ncoc, ncod, nsgfa, nsgfb, &
                                                            nsgfc, nsgfd

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'abcd_contract'

      INTEGER                                            :: handle, isgfc, isgfd, m1, m2, m3, m4, &
                                                            msphia, msphib, msphic, msphid
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: temp_cccc, work_cpcc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: temp_cpcc, work_cppc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: cpcc, cppc, cppp

      CALL timeset(routineN, handle)

      msphia = SIZE(sphi_a, 1)
      msphib = SIZE(sphi_b, 1)
      msphic = SIZE(sphi_c, 1)
      msphid = SIZE(sphi_d, 1)

      m1 = SIZE(sabcd, 1)
      m2 = SIZE(sabcd, 2)
      m3 = SIZE(sabcd, 3)
      m4 = SIZE(sabcd, 4)

      ALLOCATE (cppp(nsgfa, m2, m3, m4), cppc(nsgfa, m2, m3, nsgfd), &
                cpcc(nsgfa, m2, nsgfc, nsgfd))

      ALLOCATE (work_cppc(nsgfa, m2, m3), temp_cpcc(nsgfa, m2, nsgfc))
      work_cppc = 0._dp
      temp_cpcc = 0._dp

      ALLOCATE (work_cpcc(nsgfa, m2), temp_cccc(nsgfa, nsgfb))
      work_cpcc = 0._dp
      temp_cccc = 0._dp

      CALL dgemm("T", "N", nsgfa, m2*m3*m4, ncoa, 1._dp, sphi_a, msphia, sabcd, m1, &
                 0.0_dp, cppp, nsgfa)
      CALL dgemm("N", "N", nsgfa*m2*m3, nsgfd, ncod, 1._dp, cppp, nsgfa*m2*m3, &
                 sphi_d, msphid, 0.0_dp, cppc, nsgfa*m2*m3)

      DO isgfd = 1, nsgfd
         work_cppc(:, :, :) = cppc(:, :, :, isgfd)
         CALL dgemm("N", "N", nsgfa*m2, nsgfc, ncoc, 1._dp, work_cppc, nsgfa*m2, &
                    sphi_c, msphic, 0.0_dp, temp_cpcc, nsgfa*m2)
         cpcc(:, :, :, isgfd) = temp_cpcc(:, :, :)
         DO isgfc = 1, nsgfc
            work_cpcc(:, :) = cpcc(:, :, isgfc, isgfd)
            CALL dgemm("N", "N", nsgfa, nsgfb, ncob, 1._dp, work_cpcc, nsgfa, sphi_b, &
                       msphib, 0.0_dp, temp_cccc, nsgfa)
            abcdint(:, :, isgfc, isgfd) = temp_cccc(:, :)
         END DO
      END DO

      DEALLOCATE (cpcc, cppc, cppp)
      DEALLOCATE (work_cpcc, work_cppc, temp_cpcc, temp_cccc)

      CALL timestop(handle)

   END SUBROUTINE abcd_contract

! **************************************************************************************************
!> \brief 3-center contraction routine from primitive cartesain Gaussians to spherical Gaussian
!>        functions. Exploits LIBXSMM for performance, falls back to BLAS if LIBXSMM not available.
!>        Requires pre-allocation of work buffers and pre-transposition of the sphi_a array. Requires
!>        the LIBXSMM library to be initialized somewhere before this routine is called.
!> \param abcint contracted integrals
!> \param sabc uncontracted integrals
!> \param tsphi_a assumed to have dimensions nsgfa x ncoa
!> \param sphi_b assumed to have dimensions ncob x nsgfb
!> \param sphi_c assumed to have dimensions ncoc x nsgfc
!> \param ncoa ...
!> \param ncob ...
!> \param ncoc ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param cpp_buffer ...
!> \param ccp_buffer ...
!> \note tested from version 1.9.0 of libxsmm
! **************************************************************************************************
   SUBROUTINE libxsmm_abc_contract(abcint, sabc, tsphi_a, sphi_b, sphi_c, ncoa, ncob, ncoc, &
                                   nsgfa, nsgfb, nsgfc, cpp_buffer, ccp_buffer)

      REAL(dp), DIMENSION(*)                             :: abcint, sabc
      REAL(dp), DIMENSION(:, :)                          :: tsphi_a, sphi_b, sphi_c
      INTEGER, INTENT(IN)                                :: ncoa, ncob, ncoc, nsgfa, nsgfb, nsgfc
      REAL(dp), DIMENSION(nsgfa*ncob)                    :: cpp_buffer
      REAL(dp), DIMENSION(nsgfa*nsgfb*ncoc)              :: ccp_buffer

      CHARACTER(LEN=*), PARAMETER :: routineN = 'libxsmm_abc_contract'

      INTEGER                                            :: handle, i
#if(__LIBXSMM)
      INTEGER(libxsmm_blasint_kind)                      :: m, n, k
      TYPE(libxsmm_dmmfunction)                          :: xmm1, xmm2
#endif

      CALL timeset(routineN, handle)

#if(__LIBXSMM)
      !We make use of libxsmm code dispatching feature and call the same kernel multiple times
      !We loop over the last index of the matrix and call libxsmm each time

      !pre-fetch the kernels
      m = nsgfa; n = ncob; k = ncoa
      CALL libxsmm_dispatch(xmm1, m, n, k, beta=0.0_dp, prefetch=LIBXSMM_PREFETCH_NONE)
      m = nsgfa; n = nsgfb; k = ncob
      CALL libxsmm_dispatch(xmm2, m, n, k, beta=0.0_dp, prefetch=LIBXSMM_PREFETCH_NONE)

      !contractions over a and b
      DO i = 1, ncoc
         CALL libxsmm_dmmcall(xmm1, tsphi_a, sabc((i - 1)*ncoa*ncob + 1), cpp_buffer)
         CALL libxsmm_dmmcall(xmm2, cpp_buffer, sphi_b, ccp_buffer((i - 1)*nsgfa*nsgfb + 1))
      END DO

#else

      !we follow the same flow as for libxsmm above, but use BLAS dgemm

      !contractions over a and b
      DO i = 1, ncoc
         CALL dgemm("N", "N", nsgfa, ncob, ncoa, 1.0_dp, tsphi_a, nsgfa, sabc((i - 1)*ncoa*ncob + 1), &
                    ncoa, 0.0_dp, cpp_buffer, nsgfa)
         CALL dgemm("N", "N", nsgfa, nsgfb, ncob, 1.0_dp, cpp_buffer, nsgfa, sphi_b, ncob, 0.0_dp, &
                    ccp_buffer((i - 1)*nsgfa*nsgfb + 1), nsgfa)
      END DO

#endif

      !last contraciton, over c, as a larger MM
      CALL dgemm("N", "N", nsgfa*nsgfb, nsgfc, ncoc, 1.0_dp, ccp_buffer, nsgfa*nsgfb, sphi_c, ncoc, &
                 0.0_dp, abcint, nsgfa*nsgfb)

      CALL timestop(handle)

   END SUBROUTINE libxsmm_abc_contract

END MODULE ai_contraction_sphi
