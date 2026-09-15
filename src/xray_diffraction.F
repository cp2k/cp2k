!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par Literature
!>      M. Krack, A. Gambirasio, and M. Parrinello,
!>      "Ab-initio x-ray scattering of liquid water",
!>      J. Chem. Phys. 117, 9409 (2002)
!> \author Matthias Krack
!> \date   30.11.2005
! **************************************************************************************************
MODULE xray_diffraction
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Krack2002,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE gapw_gspace_reconstruction,      ONLY: calculate_rhotot_elec_gspace
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathconstants,                   ONLY: twopi
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grids,                        ONLY: get_pw_grid_info
   USE pw_methods,                      ONLY: pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xray_diffraction'

   PUBLIC :: calculate_rhotot_elec_gspace, &
             xray_diffraction_spectrum

CONTAINS

! **************************************************************************************************
!> \brief Calculate the coherent X-ray diffraction spectrum using the total
!>        electronic density in reciprocal space (g-space).
!> \param qs_env ...
!> \param unit_number ...
!> \param q_max ...
!> \date   30.11.2005
!> \author Matthias Krack
! **************************************************************************************************
   SUBROUTINE xray_diffraction_spectrum(qs_env, unit_number, q_max)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_number
      REAL(KIND=dp), INTENT(IN)                          :: q_max

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xray_diffraction_spectrum'
      INTEGER, PARAMETER                                 :: nblock = 100

      INTEGER                                            :: handle, i, ig, ig_shell, ipe, ishell, &
                                                            jg, ng, npe, nshell, nshell_gather
      INTEGER(KIND=int_8)                                :: ngpts
      INTEGER, DIMENSION(3)                              :: npts
      INTEGER, DIMENSION(:), POINTER                     :: aux_index, ng_shell, ng_shell_gather, &
                                                            nshell_pe, offset_pe
      REAL(KIND=dp)                                      :: cutoff, f, f2, q, rho_hard, rho_soft, &
                                                            rho_total
      REAL(KIND=dp), DIMENSION(3)                        :: dg, dr
      REAL(KIND=dp), DIMENSION(:), POINTER :: f2sum, f2sum_gather, f4sum, f4sum_gather, fmax, &
         fmax_gather, fmin, fmin_gather, fsum, fsum_gather, gsq, q_shell, q_shell_gather
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type)                               :: rhotot_elec_gspace
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set

      CPASSERT(ASSOCIATED(qs_env))

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set)
      NULLIFY (aux_index)
      NULLIFY (auxbas_pw_pool)
      NULLIFY (dft_control)
      NULLIFY (f2sum)
      NULLIFY (f2sum_gather)
      NULLIFY (f4sum)
      NULLIFY (f4sum_gather)
      NULLIFY (fmax)
      NULLIFY (fmax_gather)
      NULLIFY (fmin)
      NULLIFY (fmin_gather)
      NULLIFY (fsum)
      NULLIFY (fsum_gather)
      NULLIFY (gsq)
      NULLIFY (ng_shell)
      NULLIFY (ng_shell_gather)
      NULLIFY (nshell_pe)
      NULLIFY (offset_pe)
      NULLIFY (para_env)
      NULLIFY (particle_set)
      NULLIFY (pw_env)
      NULLIFY (q_shell)
      NULLIFY (q_shell_gather)
      NULLIFY (rho)
      NULLIFY (rho_atom_set)

      CALL cite_reference(Krack2002)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      particle_set=particle_set, &
                      pw_env=pw_env, &
                      rho=rho, &
                      rho_atom_set=rho_atom_set)

      CALL pw_env_get(pw_env=pw_env, &
                      auxbas_pw_pool=auxbas_pw_pool)

      npe = para_env%num_pe

      ! Plane waves grid to assemble the total electronic density

      CALL auxbas_pw_pool%create_pw(pw=rhotot_elec_gspace)
      CALL pw_zero(rhotot_elec_gspace)

      CALL get_pw_grid_info(pw_grid=rhotot_elec_gspace%pw_grid, &
                            dr=dr, &
                            npts=npts, &
                            cutoff=cutoff, &
                            ngpts=ngpts, &
                            gsquare=gsq)

      dg(:) = twopi/(npts(:)*dr(:))

      ! Build the total electronic density in reciprocal space

      CALL calculate_rhotot_elec_gspace(qs_env=qs_env, &
                                        auxbas_pw_pool=auxbas_pw_pool, &
                                        rhotot_elec_gspace=rhotot_elec_gspace, &
                                        q_max=q_max, &
                                        rho_hard=rho_hard, &
                                        rho_soft=rho_soft)

      rho_total = rho_hard + rho_soft

      ! Calculate the coherent X-ray spectrum

      ! Now we have to gather the data from all processes, since each
      ! process has only worked his sub-grid

      ! Scan the g-vector shells

      CALL reallocate(q_shell, 1, nblock)
      CALL reallocate(ng_shell, 1, nblock)

      ng = SIZE(gsq)

      jg = 1
      nshell = 1
      q_shell(1) = SQRT(gsq(1))
      ng_shell(1) = 1

      DO ig = 2, ng
         CPASSERT(gsq(ig) >= gsq(jg))
         IF (ABS(gsq(ig) - gsq(jg)) > 1.0E-12_dp) THEN
            nshell = nshell + 1
            IF (nshell > SIZE(q_shell)) THEN
               CALL reallocate(q_shell, 1, SIZE(q_shell) + nblock)
               CALL reallocate(ng_shell, 1, SIZE(ng_shell) + nblock)
            END IF
            q = SQRT(gsq(ig))
            IF (q > q_max) THEN
               nshell = nshell - 1
               EXIT
            END IF
            q_shell(nshell) = q
            ng_shell(nshell) = 1
            jg = ig
         ELSE
            ng_shell(nshell) = ng_shell(nshell) + 1
         END IF
      END DO

      CALL reallocate(q_shell, 1, nshell)
      CALL reallocate(ng_shell, 1, nshell)
      CALL reallocate(fmin, 1, nshell)
      CALL reallocate(fmax, 1, nshell)
      CALL reallocate(fsum, 1, nshell)
      CALL reallocate(f2sum, 1, nshell)
      CALL reallocate(f4sum, 1, nshell)

      ig = 0
      DO ishell = 1, nshell
         fmin(ishell) = HUGE(0.0_dp)
         fmax(ishell) = 0.0_dp
         fsum(ishell) = 0.0_dp
         f2sum(ishell) = 0.0_dp
         f4sum(ishell) = 0.0_dp
         DO ig_shell = 1, ng_shell(ishell)
            f = ABS(rhotot_elec_gspace%array(ig + ig_shell))
            fmin(ishell) = MIN(fmin(ishell), f)
            fmax(ishell) = MAX(fmax(ishell), f)
            fsum(ishell) = fsum(ishell) + f
            f2 = f*f
            f2sum(ishell) = f2sum(ishell) + f2
            f4sum(ishell) = f4sum(ishell) + f2*f2
         END DO
         ig = ig + ng_shell(ishell)
      END DO

      CALL reallocate(nshell_pe, 0, npe - 1)
      CALL reallocate(offset_pe, 0, npe - 1)

      ! Root (source) process gathers the number of shell of each process

      CALL para_env%gather(nshell, nshell_pe)

      ! Only the root process which has to print the full spectrum has to
      ! allocate here the receive buffers with their real sizes

      IF (unit_number > 0) THEN
         nshell_gather = SUM(nshell_pe)
         offset_pe(0) = 0
         DO ipe = 1, npe - 1
            offset_pe(ipe) = offset_pe(ipe - 1) + nshell_pe(ipe - 1)
         END DO
      ELSE
         nshell_gather = 1 ! dummy value for the non-root processes
      END IF

      CALL reallocate(q_shell_gather, 1, nshell_gather)
      CALL reallocate(ng_shell_gather, 1, nshell_gather)
      CALL reallocate(fmin_gather, 1, nshell_gather)
      CALL reallocate(fmax_gather, 1, nshell_gather)
      CALL reallocate(fsum_gather, 1, nshell_gather)
      CALL reallocate(f2sum_gather, 1, nshell_gather)
      CALL reallocate(f4sum_gather, 1, nshell_gather)

      CALL para_env%gatherv(q_shell, q_shell_gather, nshell_pe, offset_pe)
      CALL para_env%gatherv(ng_shell, ng_shell_gather, nshell_pe, offset_pe)
      CALL para_env%gatherv(fmax, fmax_gather, nshell_pe, offset_pe)
      CALL para_env%gatherv(fmin, fmin_gather, nshell_pe, offset_pe)
      CALL para_env%gatherv(fsum, fsum_gather, nshell_pe, offset_pe)
      CALL para_env%gatherv(f2sum, f2sum_gather, nshell_pe, offset_pe)
      CALL para_env%gatherv(f4sum, f4sum_gather, nshell_pe, offset_pe)

      IF (ASSOCIATED(offset_pe)) THEN
         DEALLOCATE (offset_pe)
      END IF

      IF (ASSOCIATED(nshell_pe)) THEN
         DEALLOCATE (nshell_pe)
      END IF

      ! Print X-ray diffraction spectrum (I/O node only)

      IF (unit_number > 0) THEN

         CALL reallocate(aux_index, 1, nshell_gather)

         ! Sort the gathered shells

         CALL sort(q_shell_gather, nshell_gather, aux_index)

         ! Allocate final arrays of sufficient size, i.e. nshell_gather
         ! is always greater or equal the final nshell value

         CALL reallocate(q_shell, 1, nshell_gather)
         CALL reallocate(ng_shell, 1, nshell_gather)
         CALL reallocate(fmin, 1, nshell_gather)
         CALL reallocate(fmax, 1, nshell_gather)
         CALL reallocate(fsum, 1, nshell_gather)
         CALL reallocate(f2sum, 1, nshell_gather)
         CALL reallocate(f4sum, 1, nshell_gather)

         jg = 1
         nshell = 1
         q_shell(1) = q_shell_gather(1)
         i = aux_index(1)
         ng_shell(1) = ng_shell_gather(i)
         fmin(1) = fmin_gather(i)
         fmax(1) = fmax_gather(i)
         fsum(1) = fsum_gather(i)
         f2sum(1) = f2sum_gather(i)
         f4sum(1) = f4sum_gather(i)

         DO ig = 2, nshell_gather
            i = aux_index(ig)
            IF (ABS(q_shell_gather(ig) - q_shell_gather(jg)) > 1.0E-12_dp) THEN
               nshell = nshell + 1
               q_shell(nshell) = q_shell_gather(ig)
               ng_shell(nshell) = ng_shell_gather(i)
               fmin(nshell) = fmin_gather(i)
               fmax(nshell) = fmax_gather(i)
               fsum(nshell) = fsum_gather(i)
               f2sum(nshell) = f2sum_gather(i)
               f4sum(nshell) = f4sum_gather(i)
               jg = ig
            ELSE
               ng_shell(nshell) = ng_shell(nshell) + ng_shell_gather(i)
               fmin(nshell) = MIN(fmin(nshell), fmin_gather(i))
               fmax(nshell) = MAX(fmax(nshell), fmax_gather(i))
               fsum(nshell) = fsum(nshell) + fsum_gather(i)
               f2sum(nshell) = f2sum(nshell) + f2sum_gather(i)
               f4sum(nshell) = f4sum(nshell) + f4sum_gather(i)
            END IF
         END DO

         ! The auxiliary index array is no longer needed now

         IF (ASSOCIATED(aux_index)) THEN
            DEALLOCATE (aux_index)
         END IF

         ! Allocate the final arrays for printing with their real size

         CALL reallocate(q_shell, 1, nshell)
         CALL reallocate(ng_shell, 1, nshell)
         CALL reallocate(fmin, 1, nshell)
         CALL reallocate(fmax, 1, nshell)
         CALL reallocate(fsum, 1, nshell)
         CALL reallocate(f2sum, 1, nshell)
         CALL reallocate(f4sum, 1, nshell)

         ! Write the X-ray diffraction spectrum to the specified file

         WRITE (UNIT=unit_number, FMT="(A)") &
            "#", &
            "# Coherent X-ray diffraction spectrum", &
            "#"
         WRITE (UNIT=unit_number, FMT="(A,1X,F20.10)") &
            "# Soft electronic charge (G-space) :", rho_soft, &
            "# Hard electronic charge (G-space) :", rho_hard, &
            "# Total electronic charge (G-space):", rho_total, &
            "# Density cutoff [Rydberg]         :", 2.0_dp*cutoff, &
            "# q(min) [1/Angstrom]              :", q_shell(2)/angstrom, &
            "# q(max) [1/Angstrom]              :", q_shell(nshell)/angstrom, &
            "# q(max) [1/Angstrom] (requested)  :", q_max/angstrom
         WRITE (UNIT=unit_number, FMT="(A,2X,I8)") &
            "# Number of g-vectors (grid points):", ngpts, &
            "# Number of g-vector shells        :", nshell
         WRITE (UNIT=unit_number, FMT="(A,3(1X,I6))") &
            "# Grid size (a,b,c)                :", npts(1:3)
         WRITE (UNIT=unit_number, FMT="(A,3F7.3)") &
            "# dg [1/Angstrom]                  :", dg(1:3)/angstrom, &
            "# dr [Angstrom]                    :", dr(1:3)*angstrom
         WRITE (UNIT=unit_number, FMT="(A)") &
            "#", &
            "# shell  points         q [1/A]      <|F(q)|^2>     Min(|F(q)|)"// &
            "     Max(|F(q)|)      <|F(q)|>^2      <|F(q)|^4>"

         DO ishell = 1, nshell
            WRITE (UNIT=unit_number, FMT="(T2,I6,2X,I6,5(1X,F15.6),1X,ES15.6)") &
               ishell, &
               ng_shell(ishell), &
               q_shell(ishell)/angstrom, &
               f2sum(ishell)/REAL(ng_shell(ishell), KIND=dp), &
               fmin(ishell), &
               fmax(ishell), &
               (fsum(ishell)/REAL(ng_shell(ishell), KIND=dp))**2, &
               f4sum(ishell)/REAL(ng_shell(ishell), KIND=dp)
         END DO

      END IF

      ! Release work storage

      IF (ASSOCIATED(fmin)) THEN
         DEALLOCATE (fmin)
      END IF

      IF (ASSOCIATED(fmax)) THEN
         DEALLOCATE (fmax)
      END IF

      IF (ASSOCIATED(fsum)) THEN
         DEALLOCATE (fsum)
      END IF

      IF (ASSOCIATED(f2sum)) THEN
         DEALLOCATE (f2sum)
      END IF

      IF (ASSOCIATED(f4sum)) THEN
         DEALLOCATE (f4sum)
      END IF

      IF (ASSOCIATED(ng_shell)) THEN
         DEALLOCATE (ng_shell)
      END IF

      IF (ASSOCIATED(q_shell)) THEN
         DEALLOCATE (q_shell)
      END IF

      IF (ASSOCIATED(fmin_gather)) THEN
         DEALLOCATE (fmin_gather)
      END IF

      IF (ASSOCIATED(fmax_gather)) THEN
         DEALLOCATE (fmax_gather)
      END IF

      IF (ASSOCIATED(fsum_gather)) THEN
         DEALLOCATE (fsum_gather)
      END IF

      IF (ASSOCIATED(f2sum_gather)) THEN
         DEALLOCATE (f2sum_gather)
      END IF

      IF (ASSOCIATED(f4sum_gather)) THEN
         DEALLOCATE (f4sum_gather)
      END IF

      IF (ASSOCIATED(ng_shell_gather)) THEN
         DEALLOCATE (ng_shell_gather)
      END IF

      IF (ASSOCIATED(q_shell_gather)) THEN
         DEALLOCATE (q_shell_gather)
      END IF

      CALL auxbas_pw_pool%give_back_pw(rhotot_elec_gspace)

      CALL timestop(handle)

   END SUBROUTINE xray_diffraction_spectrum

END MODULE xray_diffraction
