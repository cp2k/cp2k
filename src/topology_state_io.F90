!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Versioned Gaussian-basis snapshots for cross-geometry topology overlaps.
!>        Atomic units, full contraction coefficients, CP2K periodic atom images.
! **************************************************************************************************
MODULE topology_state_io
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_files,                        ONLY: open_file
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_state_io'
   PUBLIC :: topology_state_begin, topology_state_point
CONTAINS

! **************************************************************************************************
!> \brief Write basis metadata; all ranks call, only the source writes.
!> \param qs_env electronic structure environment
!> \param filename snapshot file
!> \param nao number of AO functions
!> \param nband retained states
!> \param npoint k-points
!> \param nspin spinor components, one for a selected collinear channel
!> \param ntotal all computed eigenvalues
!> \param spin_channel collinear channel label
!> \param bands retained band indices
!> \param unit source output unit, -1 on other ranks
! **************************************************************************************************
   SUBROUTINE topology_state_begin(qs_env, filename, nao, nband, npoint, nspin, ntotal, spin_channel, bands, unit)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(len=*), INTENT(IN)                       :: filename
      INTEGER, INTENT(IN)                                :: nao, nband, npoint, nspin, ntotal, &
                                                            spin_channel, bands(:)
      INTEGER, INTENT(OUT)                               :: unit

      INTEGER                                            :: count_ao, first, ia, ic, ik, ip, iset, &
                                                            j, ncart, ns, row
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds_set(:)

      CALL get_qs_env(qs_env, cell=cell, particle_set=particles, qs_kind_set=kinds_set, para_env=para_env)
      unit = -1
      IF (.NOT. para_env%is_source()) RETURN
      CALL open_file(filename, unit_number=unit, file_status="UNKNOWN", file_action="WRITE")
      WRITE (unit, '(A)') 'CP2K_TOPOLOGY_STATE 1'
      WRITE (unit, '(7I12)') SIZE(particles), nao, nband, npoint, nspin, ntotal, spin_channel
      WRITE (unit, '(*(I12))') bands
      DO j = 1, 3
         WRITE (unit, '(3ES26.17)') cell%hmat(:, j)
      END DO
      WRITE (unit, '(3I12)') cell%perd
      count_ao = 0
      DO ia = 1, SIZE(particles)
         CALL get_atomic_kind(particles(ia)%atomic_kind, kind_number=ik)
         CALL get_qs_kind(kinds_set(ik), basis_set=basis)
         IF (.NOT. ASSOCIATED(basis)) THEN
            CPABORT("STATE_EXPORT requires a Gaussian orbital basis on every atom.")
         END IF
         WRITE (unit, '(4I12)') ia, ik, basis%nset, basis%nsgf
         WRITE (unit, '(3ES26.17)') pbc(particles(ia)%r, cell)
         DO iset = 1, basis%nset
            first = basis%first_sgf(1, iset)
            ns = basis%nsgf_set(iset)
            ncart = ncoset(basis%lmax(iset))
            WRITE (unit, '(5I12,ES26.17)') first, ns, basis%npgf(iset), ncart, basis%lmin(iset), &
               basis%set_radius(iset)
            DO ip = 1, basis%npgf(iset)
               WRITE (unit, '(2ES26.17)') basis%zet(ip, iset), basis%pgf_radius(ip, iset)
               DO ic = 1, ncart
                  row = (ip - 1)*ncart + ic
                  WRITE (unit, '(3I5,*(ES26.17))') indco(:, ic), basis%sphi(row, first:first + ns - 1)
               END DO
            END DO
         END DO
         count_ao = count_ao + basis%nsgf
      END DO
      IF (count_ao /= nao) CPABORT("STATE_EXPORT: inconsistent atom/AO dimensions.")
   END SUBROUTINE topology_state_begin

! **************************************************************************************************
!> \brief Write one point with eigenvalues in Hartree and physical AO coefficients.
!> \param unit snapshot output unit, -1 skips output
!> \param point one-based point index
!> \param k fractional reciprocal coordinates
!> \param eigenvalues all computed eigenvalues, including excluded states
!> \param coefficients selected states; spinor components are stacked by AO
! **************************************************************************************************
   SUBROUTINE topology_state_point(unit, point, k, eigenvalues, coefficients)
      INTEGER, INTENT(IN)                                :: unit, point
      REAL(KIND=dp), INTENT(IN)                          :: k(3), eigenvalues(:)
      COMPLEX(KIND=dp), INTENT(IN)                       :: coefficients(:, :)

      INTEGER                                            :: i, j

      IF (unit < 0) RETURN
      WRITE (unit, '(I12,3ES26.17)') point, k
      WRITE (unit, '(*(ES26.17))') eigenvalues
      DO j = 1, SIZE(coefficients, 2)
         DO i = 1, SIZE(coefficients, 1)
            WRITE (unit, '(2ES26.17)') REAL(coefficients(i, j), dp), AIMAG(coefficients(i, j))
         END DO
      END DO
   END SUBROUTINE topology_state_point
END MODULE topology_state_io
