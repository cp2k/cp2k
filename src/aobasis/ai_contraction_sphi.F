!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contraction of integrals over primitive Cartesian Gaussians based on the contraction
!>        matrix sphi which is part of the gto_basis_set_type
!> \par History
!>      -added abc_contract_xsmm routine, A. Bussy (04.2020)
!> \author Dorothea Golze (05.2016)
! **************************************************************************************************
MODULE ai_contraction_sphi

   USE kinds, ONLY: dp, int_8
#if defined(__LIBXSMM)
   USE libxsmm, ONLY: LIBXSMM_PREFETCH_NONE, &
                      libxsmm_blasint_kind, &
                      libxsmm_dgemm, &
                      libxsmm_dispatch, &
                      libxsmm_available, &
                      libxsmm_dmmcall, &
                      libxsmm_dmmfunction
#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_contraction_sphi'

   PUBLIC :: ab_contract, abc_contract, abcd_contract, abc_contract_xsmm

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

      INTEGER                                            :: handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cpp

      CALL timeset(routineN, handle)

      CPASSERT(ncob <= SIZE(sab, 2))

      IF ((nsgfa*ncob*(ncoa + nsgfb)) <= (nsgfb*ncoa*(ncob + nsgfa))) THEN ! (sphi_a x sab) x sphi_b
         ALLOCATE (cpp(nsgfa, ncob))
         ! [nsgfa,ncoa] x [ncoa,ncob] -> [nsgfa,ncob]
         CALL dgemm("T", "N", nsgfa, ncob, ncoa, 1._dp, sphi_a, SIZE(sphi_a, 1), sab, SIZE(sab, 1), 0.0_dp, cpp, nsgfa)
         ! [nsgfa,ncob] x [ncob,nsgfb] -> [nsgfa,nsgfb]
         CALL dgemm("N", "N", nsgfa, nsgfb, ncob, 1._dp, cpp, nsgfa, sphi_b, SIZE(sphi_b, 1), 0.0_dp, &
                    abint, SIZE(abint, 1))
      ELSE ! sphi_a x (sab x sphi_b)
         ALLOCATE (cpp(ncoa, nsgfb))
         ! [ncoa,ncob] x [ncob,nsgfb] -> [ncoa,nsgfb]
         CALL dgemm("N", "N", ncoa, nsgfb, ncob, 1._dp, sab, SIZE(sab, 1), sphi_b, SIZE(sphi_b, 1), 0.0_dp, cpp, ncoa)
         ! [nsgfa,ncoa] x [ncoa,nsgfb] -> [nsgfa,nsgfb]
         CALL dgemm("T", "N", nsgfa, nsgfb, ncoa, 1._dp, sphi_a, SIZE(sphi_a, 1), cpp, ncoa, 0.0_dp, &
                    abint, SIZE(abint, 1))
      END IF

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
                                                            msphic, mx
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: tmp

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(abcint, 1) == nsgfa)
      CPASSERT(SIZE(abcint, 2) == nsgfb)

      msphia = SIZE(sphi_a, 1)
      msphib = SIZE(sphi_b, 1)
      msphic = SIZE(sphi_c, 1)

      m1 = SIZE(sabc, 1)
      m2 = SIZE(sabc, 2)
      m3 = SIZE(sabc, 3)
      mx = MAX(m2, nsgfb)

      ! ALLOCATE (cpp(nsgfa, m2, m3), cpc(nsgfa, nsgfb, m3))
      ALLOCATE (tmp(nsgfa, mx, m3 + 1))

      CALL dgemm("T", "N", nsgfa, m2*m3, ncoa, 1._dp, sphi_a, msphia, sabc, m1, 0.0_dp, tmp(:, :, 2), nsgfa)
      DO i = 1, m3
         CALL dgemm("N", "N", nsgfa, nsgfb, ncob, 1._dp, tmp(:, :, i + 1), nsgfa, sphi_b, msphib, &
                    0.0_dp, tmp(:, :, i), nsgfa)
      END DO
      CALL dgemm("N", "N", nsgfa*nsgfb, nsgfc, ncoc, 1._dp, tmp, nsgfa*mx, sphi_c, msphic, 0.0_dp, &
                 abcint, nsgfa*nsgfb)

      DEALLOCATE (tmp)

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
!> \brief 3-center contraction routine from primitive cartesian Gaussians to spherical Gaussian
!>        functions; can use LIBXSMM (falls back to BLAS otherwise).
!>        Requires pre-transposition of the sphi_a array. The call-side shall DEALLOCATE buffers
!>        end of scope or after last use. This function ALLOCATEs or grows the work buffers
!>        as necessary. LIBXSMM may be initialized upfront (elsewhere).
!> \param abcint contracted integrals
!> \param sabc uncontracted integrals
!> \param sphi_a assumed to have dimensions nsgfa x ncoa
!> \param sphi_b assumed to have dimensions ncob x nsgfb
!> \param sphi_c assumed to have dimensions ncoc x nsgfc
!> \param ncoa ...
!> \param ncob ...
!> \param ncoc ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param cpp_buffer Buffer used for intermediate results (automatically allocated).
!> \param ccp_buffer Buffer used for intermediate results (automatically allocated).
!> \param prefac Prefactor which is finally multiplied into abcint (default: 1.0).
!> \param pstfac Factor used to consider initial abcint (default: 0.0).
! **************************************************************************************************
   SUBROUTINE abc_contract_xsmm(abcint, sabc, sphi_a, sphi_b, sphi_c, ncoa, ncob, ncoc, &
                                nsgfa, nsgfb, nsgfc, cpp_buffer, ccp_buffer, prefac, pstfac)

      REAL(KIND=dp), DIMENSION(:, :, :)           :: abcint
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)     :: sabc
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)  :: sphi_a, sphi_b, sphi_c
      INTEGER, INTENT(IN)                         :: ncoa, ncob, ncoc, nsgfa, nsgfb, nsgfc
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: cpp_buffer, ccp_buffer
      REAL(KIND=dp), INTENT(IN), OPTIONAL         :: prefac, pstfac

      CHARACTER(LEN=*), PARAMETER :: routineN = 'abc_contract_xsmm'

      REAL(KIND=dp)                                           :: alpha, beta
      INTEGER(KIND=int_8)                                     :: cpp_size, ccp_size
      INTEGER                                                 :: handle, i
      LOGICAL                                                 :: ab_first
#if defined(__LIBXSMM)
      TYPE(libxsmm_dmmfunction)                               :: xmm1, xmm2
#endif

      CALL timeset(routineN, handle)

      alpha = 1.0_dp
      IF (PRESENT(prefac)) alpha = prefac

      beta = 0.0_dp
      IF (PRESENT(pstfac)) beta = pstfac

      ! M*N*K FLOPS are used to decide if contracting (AB)C vs A(BC)
      IF ((nsgfa*ncob*(ncoa + nsgfb)) <= (ncoa*nsgfb*(ncob + nsgfa))) THEN
         cpp_size = nsgfa*ncob
         ab_first = .TRUE.
      ELSE
         cpp_size = ncoa*nsgfb
         ab_first = .FALSE.
      END IF

      ccp_size = nsgfa*nsgfb*ncoc
      IF (.NOT. ALLOCATED(ccp_buffer)) THEN
         ALLOCATE (ccp_buffer(ccp_size))
      ELSE IF (SIZE(ccp_buffer) < ccp_size) THEN
         DEALLOCATE (ccp_buffer)
         ALLOCATE (ccp_buffer(ccp_size))
      END IF

      IF (.NOT. ALLOCATED(cpp_buffer)) THEN
         ALLOCATE (cpp_buffer(cpp_size))
      ELSE IF (SIZE(cpp_buffer) < cpp_size) THEN
         DEALLOCATE (cpp_buffer)
         ALLOCATE (cpp_buffer(cpp_size))
      END IF

      ! loop over the last index of the matrix and call LIBXSMM/BLAS to contract over a and b
#if defined(__LIBXSMM)
      IF (ab_first) THEN ! (AB)C: dispatch kernels
         CALL libxsmm_dispatch(xmm1, nsgfa, ncob, ncoa, beta=0.0_dp, prefetch=LIBXSMM_PREFETCH_NONE)
         CALL libxsmm_dispatch(xmm2, nsgfa, nsgfb, ncob, beta=0.0_dp, prefetch=LIBXSMM_PREFETCH_NONE)
      ELSE ! A(BC): dispatch kernels
         CALL libxsmm_dispatch(xmm1, ncoa, nsgfb, ncob, beta=0.0_dp, prefetch=LIBXSMM_PREFETCH_NONE)
         CALL libxsmm_dispatch(xmm2, nsgfa, nsgfb, ncoa, beta=0.0_dp, prefetch=LIBXSMM_PREFETCH_NONE)
      END IF

      IF (libxsmm_available(xmm1) .AND. libxsmm_available(xmm2)) THEN
         IF (ab_first) THEN ! (AB)C
            DO i = 0, ncoc - 1 ! contractions over a and b
               ! [nsgfa,ncoa] x [ncoa,ncob] -> [nsgfa,ncob]
               CALL libxsmm_dmmcall(xmm1, sphi_a, sabc(i*ncoa*ncob + 1), cpp_buffer)
               ! [nsgfa,ncob] x [ncob,nsgfb] -> [nsgfa,nsgfb]
               CALL libxsmm_dmmcall(xmm2, cpp_buffer, sphi_b, ccp_buffer(i*nsgfa*nsgfb + 1))
            END DO
         ELSE ! A(BC)
            DO i = 0, ncoc - 1 ! contractions over a and b
               ! [ncoa,ncob] x [ncob,nsgfb] -> [ncoa,nsgfb]
               CALL libxsmm_dmmcall(xmm1, sabc(i*ncoa*ncob + 1), sphi_b, cpp_buffer)
               ! [nsgfa,ncoa] x [ncoa,nsgfb] -> [nsgfa,nsgfb]
               CALL libxsmm_dmmcall(xmm2, sphi_a, cpp_buffer, ccp_buffer(i*nsgfa*nsgfb + 1))
            END DO
         END IF
      ELSE
#endif
         IF (ab_first) THEN ! (AB)C
            DO i = 0, ncoc - 1 ! contractions over a and b
               CALL dgemm("N", "N", nsgfa, ncob, ncoa, 1.0_dp, sphi_a, nsgfa, sabc(i*ncoa*ncob + 1), &
                          ncoa, 0.0_dp, cpp_buffer, nsgfa) ! [nsgfa,ncoa] x [ncoa,ncob] -> [nsgfa,ncob]
               CALL dgemm("N", "N", nsgfa, nsgfb, ncob, 1.0_dp, cpp_buffer, nsgfa, sphi_b, ncob, 0.0_dp, &
                          ccp_buffer(i*nsgfa*nsgfb + 1), nsgfa) ! [nsgfa,ncob] x [ncob,nsgfb] -> [nsgfa,nsgfb]
            END DO
         ELSE ! A(BC)
            DO i = 0, ncoc - 1 ! contractions over a and b
               CALL dgemm("N", "N", ncoa, nsgfb, ncob, 1.0_dp, sabc(i*ncoa*ncob + 1), ncoa, sphi_b, &
                          ncob, 0.0_dp, cpp_buffer, ncoa) ! [ncoa,ncob] x [ncob,nsgfb] -> [ncoa,nsgfb]
               CALL dgemm("N", "N", nsgfa, nsgfb, ncoa, 1.0_dp, sphi_a, nsgfa, cpp_buffer, ncoa, 0.0_dp, &
                          ccp_buffer(i*nsgfa*nsgfb + 1), nsgfa) ! [nsgfa,ncoa] x [ncoa,nsgfb] -> [nsgfa,nsgfb]
            END DO
         END IF
#if defined(__LIBXSMM)
      END IF
#endif
      ! contractions over c: [nsgfa*nsgfb,ncoc] x [ncoc,nsgfc] -> [sgfa*nsgfb,nsgfc]
      CALL dgemm("N", "N", nsgfa*nsgfb, nsgfc, ncoc, alpha, ccp_buffer, nsgfa*nsgfb, &
                 sphi_c, ncoc, beta, abcint, nsgfa*nsgfb)

      CALL timestop(handle)

   END SUBROUTINE abc_contract_xsmm

END MODULE ai_contraction_sphi
