!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the full diagonalization of GW + Bethe-Salpeter for computing
!> electronic excitations
!> \par History
!>      10.2023 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_print

   USE bse_util,                        ONLY: compute_absorption_spectrum,&
                                              filter_eigvec_contrib
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE mp2_types,                       ONLY: mp2_type
   USE physcon,                         ONLY: evolt
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_print'

   PUBLIC :: print_BSE_start_flag, fm_write_thresh, print_excitation_energies, &
             print_output_header, print_transition_amplitudes, print_optical_properties

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param bse_tda ...
!> \param bse_abba ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_BSE_start_flag(bse_tda, bse_abba, unit_nr)

      LOGICAL, INTENT(IN)                                :: bse_tda, bse_abba
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_BSE_start_flag'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) ' '
         WRITE (unit_nr, '(T2,A79)') '*******************************************************************************'
         WRITE (unit_nr, '(T2,A79)') '**                                                                           **'
         WRITE (unit_nr, '(T2,A79)') '**           Bethe Salpeter equation (BSE) for excitation energies           **'
         IF (bse_tda .AND. bse_abba) THEN
            WRITE (unit_nr, '(T2,A79)') '**          solved with and without Tamm-Dancoff approximation (TDA)         **'
         ELSE IF (bse_tda) THEN
            WRITE (unit_nr, '(T2,A79)') '**                solved with Tamm-Dancoff approximation (TDA)               **'
         ELSE
            WRITE (unit_nr, '(T2,A79)') '**               solved without Tamm-Dancoff approximation (TDA)             **'
         END IF

         WRITE (unit_nr, '(T2,A79)') '**                                                                           **'
         WRITE (unit_nr, '(T2,A79)') '*******************************************************************************'
         WRITE (unit_nr, *) ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param homo ...
!> \param virtual ...
!> \param homo_irred ...
!> \param flag_TDA ...
!> \param multiplet ...
!> \param alpha ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_output_header(homo, virtual, homo_irred, flag_TDA, &
                                  multiplet, alpha, unit_nr)

      INTEGER, INTENT(IN)                                :: homo, virtual, homo_irred
      LOGICAL, INTENT(IN)                                :: flag_TDA
      CHARACTER(LEN=10), INTENT(IN)                      :: multiplet
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_output_header'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         IF (flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '**************************************************************************'
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '*   Bethe Salpeter equation (BSE) with Tamm Dancoff approximation (TDA)  *'
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '**************************************************************************'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A48,A23)') 'BSE|', 'The excitations are calculated by diagonalizing ', &
               'the BSE within the TDA:'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T29,A16)') 'BSE|', 'A X^n = Ω^n X^n'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A23)') 'BSE|', 'i.e. in index notation:'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A41)') 'BSE|', 'sum_jb ( A_ia,jb   X_jb^n ) = Ω^n X_ia^n'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A30)') 'BSE|', 'prelim Ref.: Eq. (36) with B=0'
            WRITE (unit_nr, '(T2,A4,T7,A71)') 'BSE|', 'in PRB 92,045209 (2015); http://dx.doi.org/10.1103/PhysRevB.92.045209 .'
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '**************************************************************************'
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '*      Full ("ABBA") Bethe Salpeter equation (BSE) (i.e. without TDA)    *'
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', '**************************************************************************'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A48,A24)') 'BSE|', 'The excitations are calculated by diagonalizing ', &
               'the BSE without the TDA:'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T22,A30)') 'BSE|', '|A B| |X^n|       |1  0| |X^n|'
            WRITE (unit_nr, '(T2,A4,T22,A31)') 'BSE|', '|B A| |Y^n| = Ω^n |0 -1| |Y^n|'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A23)') 'BSE|', 'i.e. in index notation:'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A62)') 'BSE|', '  sum_jb ( A_ia,jb   X_jb^n + B_ia,jb   Y_jb^n ) = Ω^n X_ia^n'
            WRITE (unit_nr, '(T2,A4,T7,A62)') 'BSE|', '- sum_jb ( B_ia,jb   X_jb^n + A_ia,jb   Y_jb^n ) = Ω^n Y_ia^n'
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A4,T18,A42,T70,A1,I4,A1,I4,A1)') 'BSE|', 'i,j:', &
            'occupied molecular orbitals, i.e. state in', '[', homo_irred - homo + 1, ',', homo_irred, ']'
         WRITE (unit_nr, '(T2,A4,T7,A4,T18,A44,T70,A1,I4,A1,I4,A1)') 'BSE|', 'a,b:', &
            'unoccupied molecular orbitals, i.e. state in', '[', homo_irred + 1, ',', homo_irred + virtual, ']'
         WRITE (unit_nr, '(T2,A4,T7,A2,T18,A16)') 'BSE|', 'n:', 'Excitation index'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A58)') 'BSE|', 'A_ia,jb = (ε_a-ε_i) δ_ij δ_ab + α * v_ia,jb - W_ij,ab'
         IF (.NOT. flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4,T7,A32)') 'BSE|', 'B_ia,jb = α * v_ia,jb - W_ib,aj'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A35)') 'BSE|', 'prelim Ref.: Eqs. (24-27),(30),(35)'
            WRITE (unit_nr, '(T2,A4,T7,A71)') 'BSE|', 'in PRB 92,045209 (2015); http://dx.doi.org/10.1103/PhysRevB.92.045209 .'
         END IF
         IF (.NOT. flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A74)') 'BSE|', 'The BSE is solved for Ω^n and X_ia^n as a hermitian problem, e.g. Eq.(42)'
            WRITE (unit_nr, '(T2,A4,T7,A71)') 'BSE|', 'in PRB 92,045209 (2015); http://dx.doi.org/10.1103/PhysRevB.92.045209 .'
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A7,T19,A23)') 'BSE|', 'ε_...:', 'GW quasiparticle energy'
         WRITE (unit_nr, '(T2,A4,T7,A7,T19,A15)') 'BSE|', 'δ_...:', 'Kronecker delta'
         WRITE (unit_nr, '(T2,A4,T7,A3,T19,A21)') 'BSE|', 'α:', 'spin-dependent factor (Singlet/Triplet)'
         WRITE (unit_nr, '(T2,A4,T7,A6,T18,A34)') 'BSE|', 'v_...:', 'Electron-hole exchange interaction'
         WRITE (unit_nr, '(T2,A4,T7,A6,T18,A27)') 'BSE|', 'W_...:', 'Screened direct interaction'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A47,A7,A9,F3.1)') 'BSE|', &
            'The spin-dependent factor is for the requested ', multiplet, " is α = ", alpha
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param Exc_ens ...
!> \param homo ...
!> \param virtual ...
!> \param flag_TDA ...
!> \param multiplet ...
!> \param info_approximation ...
!> \param mp2_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_excitation_energies(Exc_ens, homo, virtual, flag_TDA, multiplet, &
                                        info_approximation, mp2_env, unit_nr)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens
      INTEGER, INTENT(IN)                                :: homo, virtual
      LOGICAL, INTENT(IN)                                :: flag_TDA
      CHARACTER(LEN=10), INTENT(IN)                      :: multiplet, info_approximation
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_excitation_energies'

      INTEGER                                            :: handle, i_exc

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         IF (flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4,T7,A56)') 'BSE|', 'Excitation energies from solving the BSE within the TDA:'
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A57)') 'BSE|', 'Excitation energies from solving the BSE without the TDA:'
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T11,A12,T26,A11,T44,A8,T55,A27)') 'BSE|', &
            'Excitation n', "Spin Config", 'TDA/ABBA', 'Excitation energy Ω^n (eV)'
      END IF
      !prints actual energies values
      IF (unit_nr > 0) THEN
         DO i_exc = 1, MIN(homo*virtual, mp2_env%bse%num_print_exc)
            WRITE (unit_nr, '(T2,A4,T7,I16,T30,A7,T46,A6,T59,F22.4)') &
               'BSE|', i_exc, multiplet, info_approximation, Exc_ens(i_exc)*evolt
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_excitation_energies

! **************************************************************************************************
!> \brief ...
!> \param fm_eigvec_X ...
!> \param homo ...
!> \param virtual ...
!> \param homo_irred ...
!> \param info_approximation ...
!> \param mp2_env ...
!> \param unit_nr ...
!> \param fm_eigvec_Y ...
! **************************************************************************************************
   SUBROUTINE print_transition_amplitudes(fm_eigvec_X, homo, virtual, homo_irred, &
                                          info_approximation, mp2_env, unit_nr, fm_eigvec_Y)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_eigvec_X
      INTEGER, INTENT(IN)                                :: homo, virtual, homo_irred
      CHARACTER(LEN=10), INTENT(IN)                      :: info_approximation
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: fm_eigvec_Y

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_transition_amplitudes'

      INTEGER                                            :: handle, i_exc

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A61)') &
            'BSE|', "Single-particle transitions are built up by (de-)excitations,"
         WRITE (unit_nr, '(T2,A4,T7,A18)') &
            'BSE|', "which we denote by"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T20,A2,T30,A40)') &
            'BSE|', "=>", "for excitations, i.e. entries of X_ia^n,"
         WRITE (unit_nr, '(T2,A4,T20,A2,T30,A42)') &
            'BSE|', "<=", "for deexcitations, i.e. entries of Y_ia^n."
         WRITE (unit_nr, '(T2,A4)') &
            'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A73)') &
            'BSE|', "The following single-particle transitions have significant contributions,"
         WRITE (unit_nr, '(T2,A4,T7,A16,F5.3,A15,F5.3,A16)') &
            'BSE|', "i.e. |X_ia^n| > ", mp2_env%bse%eps_x, " or |Y_ia^n| > ", &
            mp2_env%bse%eps_x, ", respectively :"

         WRITE (unit_nr, '(T2,A4,T15,A27,I5,A13,I5,A3)') 'BSE|', '-- Quick reminder: HOMO i =', &
            homo_irred, ' and LUMO a =', homo_irred + 1, " --"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A12,T30,A1,T32,A5,T42,A1,T49,A8,T64,A17)') &
            "BSE|", "Excitation n", "i", "=>/<=", "a", 'TDA/ABBA', "|X_ia^n|/|Y_ia^n|"
      END IF
      DO i_exc = 1, MIN(homo*virtual, mp2_env%bse%num_print_exc)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
         END IF
         !Iterate through eigenvector and print values above threshold
         CALL print_transition_amplitudes_core(fm_eigvec_X, "=>", info_approximation, &
                                               i_exc, virtual, homo, homo_irred, &
                                               unit_nr, mp2_env)
         IF (PRESENT(fm_eigvec_Y)) THEN
            CALL print_transition_amplitudes_core(fm_eigvec_Y, "<=", info_approximation, &
                                                  i_exc, virtual, homo, homo_irred, &
                                                  unit_nr, mp2_env)
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE print_transition_amplitudes

! **************************************************************************************************
!> \brief ...
!> \param Exc_ens ...
!> \param oscill_str ...
!> \param trans_mom_bse ...
!> \param polarizability_residues ...
!> \param homo ...
!> \param virtual ...
!> \param flag_TDA ...
!> \param info_approximation ...
!> \param mp2_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_optical_properties(Exc_ens, oscill_str, trans_mom_bse, polarizability_residues, &
                                       homo, virtual, flag_TDA, &
                                       info_approximation, mp2_env, unit_nr)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens, oscill_str
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: trans_mom_bse, polarizability_residues
      INTEGER, INTENT(IN)                                :: homo, virtual
      LOGICAL, INTENT(IN)                                :: flag_TDA
      CHARACTER(LEN=10), INTENT(IN)                      :: info_approximation
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_optical_properties'

      INTEGER                                            :: handle, i_exc

      CALL timeset(routineN, handle)

      ! Discriminate between singlet and triplet, since triplet state can't couple to light
      ! and therefore calculations of dipoles etc are not necessary
      IF (mp2_env%bse%bse_spin_config == 0) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A60)') &
               'BSE|', "Transition moments d_r^n (with r∈(x,y,z), in atomic units)"
            WRITE (unit_nr, '(T2,A4,T7,A67)') &
               'BSE|', "and oscillator strength f^n of excitation level n are obtained from"
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            IF (flag_TDA) THEN
               WRITE (unit_nr, '(T2,A4,T10,A52)') &
                  'BSE|', "d_r^n = sqrt(2) sum_ia < ψ_i | r | ψ_a >  X_ia^n"
            ELSE
               WRITE (unit_nr, '(T2,A4,T10,A64)') &
                  'BSE|', "d_r^n = sqrt(2) sum_ia < ψ_i | r | ψ_a > ( X_ia^n + Y_ia^n )"
            END IF
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T14,A44)') &
               'BSE|', "f^n = 2/3 * Ω^n sum_r∈(x,y,z) ( d_r^n )^2"
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A19)') &
               'BSE|', "where we introduced"
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A5,T15,A28)') &
               'BSE|', "ψ_i:", "occupied molecular orbitals,"
            WRITE (unit_nr, '(T2,A4,T7,A5,T15,A28)') &
               'BSE|', "ψ_a:", "empty molecular orbitals and"
            WRITE (unit_nr, '(T2,A4,T9,A2,T14,A18)') &
               'BSE|', "r:", "position operator."
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A28)') &
               'BSE|', "prelim Ref.: Eqs. (23), (24)"
            WRITE (unit_nr, '(T2,A4,T7,A71)') &
               'BSE|', "in J. Chem. Phys. 152, 044105 (2020); https://doi.org/10.1063/1.5123290"
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            IF (flag_TDA) THEN
               WRITE (unit_nr, '(T2,A4,T7,A55)') 'BSE|', &
                  'Optical properties from solving the BSE within the TDA:'
            ELSE
               WRITE (unit_nr, '(T2,A4,T7,A56)') 'BSE|', &
                  'Optical properties from solving the BSE without the TDA:'
            END IF
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T8,A12,T22,A8,T38,A5,T48,A5,T58,A5,T65,A16)') 'BSE|', &
               'Excitation n', "TDA/ABBA", "d_x^n", "d_y^n", "d_z^n", 'Osc. strength'
            DO i_exc = 1, MIN(homo*virtual, mp2_env%bse%num_print_exc)
               WRITE (unit_nr, '(T2,A4,T8,I12,T24,A6,T35,F8.3,T45,F8.3,T55,F8.3,T65,F16.3)') &
                  'BSE|', i_exc, info_approximation, trans_mom_bse(1, 1, i_exc), trans_mom_bse(2, 1, i_exc), &
                  trans_mom_bse(3, 1, i_exc), oscill_str(i_exc)
            END DO
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
         END IF

         ! Compute and print the absorption spectrum to external file
         IF (mp2_env%bse%bse_print_spectrum) THEN
            CALL compute_absorption_spectrum(oscill_str, polarizability_residues, Exc_ens, &
                                             info_approximation, unit_nr, mp2_env)
         END IF

      ELSE
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            CALL cp_warn(__LOCATION__, &
                         "Requested triplet excitation cannot couple to light. "// &
                         "Skipping calculation of transition moments, "// &
                         "oscillator strengths, and spectrum.")
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_optical_properties

! **************************************************************************************************
!> \brief ...
!> \param fm_eigvec ...
!> \param direction_excitation ...
!> \param info_approximation ...
!> \param i_exc ...
!> \param virtual ...
!> \param homo ...
!> \param homo_irred ...
!> \param unit_nr ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE print_transition_amplitudes_core(fm_eigvec, direction_excitation, info_approximation, &
                                               i_exc, virtual, homo, homo_irred, &
                                               unit_nr, mp2_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_eigvec
      CHARACTER(LEN=2), INTENT(IN)                       :: direction_excitation
      CHARACTER(LEN=10), INTENT(IN)                      :: info_approximation
      INTEGER                                            :: i_exc, virtual, homo, homo_irred
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_transition_amplitudes_core'

      INTEGER                                            :: handle, k, num_entries
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: idx_homo, idx_virt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigvec_entries

      CALL timeset(routineN, handle)

      CALL filter_eigvec_contrib(fm_eigvec, idx_homo, idx_virt, eigvec_entries, &
                                 i_exc, virtual, num_entries, mp2_env)
      ! direction_excitation can be either => (means excitation; from fm_eigvec_X)
      ! or <= (means deexcitation; from fm_eigvec_Y)
      IF (unit_nr > 0) THEN
         DO k = 1, num_entries
            WRITE (unit_nr, '(T2,A4,T14,I5,T26,I5,T35,A2,T38,I5,T51,A6,T65,F16.4)') &
               "BSE|", i_exc, homo_irred - homo + idx_homo(k), direction_excitation, &
               homo_irred + idx_virt(k), info_approximation, ABS(eigvec_entries(k))
         END DO
      END IF
      DEALLOCATE (idx_homo)
      DEALLOCATE (idx_virt)
      DEALLOCATE (eigvec_entries)
      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Debug function to write elements of a full matrix to file, if they are larger than a given threshold
!> \param fm ...
!> \param thresh ...
!> \param header ...
!> \param unit_nr ...
!> \param abs_vals ...
! **************************************************************************************************
   SUBROUTINE fm_write_thresh(fm, thresh, header, unit_nr, abs_vals)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm
      REAL(KIND=dp), INTENT(IN)                          :: thresh
      CHARACTER(LEN=*), INTENT(IN)                       :: header
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, OPTIONAL                                  :: abs_vals

      CHARACTER(LEN=*), PARAMETER :: my_footer = " | ENDING WRITING OF MATRIX", &
         routineN = 'fm_write_thresh'

      INTEGER                                            :: handle, i, j, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_abs_vals

      CALL timeset(routineN, handle)

      IF (PRESENT(abs_vals)) THEN
         my_abs_vals = abs_vals
      ELSE
         my_abs_vals = .FALSE.
      END IF

      CALL cp_fm_get_info(matrix=fm, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) header
      END IF
      IF (my_abs_vals) THEN
         DO i = 1, nrow_local
            DO j = 1, ncol_local
               IF (ABS(fm%local_data(i, j)) > thresh) THEN
                  WRITE (unit_nr, "(A7,T10,I5,T20,I5,T30,F13.5)") header, row_indices(i), col_indices(j), &
                     ABS(fm%local_data(i, j))
               END IF
            END DO
         END DO
      ELSE
         DO i = 1, nrow_local
            DO j = 1, ncol_local
               IF (ABS(fm%local_data(i, j)) > thresh) THEN
                  WRITE (unit_nr, "(A7,T10,I5,T20,I5,T30,F13.5)") header, row_indices(i), col_indices(j), &
                     fm%local_data(i, j)
               END IF
            END DO
         END DO
      END IF
      CALL fm%matrix_struct%para_env%sync()
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) my_footer
      END IF

      CALL timestop(handle)

   END SUBROUTINE

END MODULE bse_print
