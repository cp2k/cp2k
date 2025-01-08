!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for printing information in context of the BSE calculation
!> \par History
!>      10.2024 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_print

   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE bse_properties,                  ONLY: compute_absorption_spectrum,&
                                              exciton_descr_type
   USE bse_util,                        ONLY: filter_eigvec_contrib
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE input_constants,                 ONLY: bse_screening_alpha,&
                                              bse_screening_rpa,&
                                              bse_screening_tdhf,&
                                              bse_screening_w0
   USE kinds,                           ONLY: dp
   USE mp2_types,                       ONLY: mp2_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_print'

   PUBLIC :: print_BSE_start_flag, fm_write_thresh, print_excitation_energies, &
             print_output_header, print_transition_amplitudes, print_optical_properties, &
             print_exciton_descriptors

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
!> \param mp2_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_output_header(homo, virtual, homo_irred, flag_TDA, &
                                  multiplet, alpha, mp2_env, unit_nr)

      INTEGER, INTENT(IN)                                :: homo, virtual, homo_irred
      LOGICAL, INTENT(IN)                                :: flag_TDA
      CHARACTER(LEN=10), INTENT(IN)                      :: multiplet
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
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
         IF (mp2_env%bse%screening_method == bse_screening_w0) THEN
            WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'A_ia,jb = (ε_a-ε_i) δ_ij δ_ab + α * v_ia,jb - W_ij,ab'
         ELSE IF (mp2_env%bse%screening_method == bse_screening_rpa) THEN
            WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'A_ia,jb = (ε_a-ε_i) δ_ij δ_ab + α * v_ia,jb'
         END IF
         IF (.NOT. flag_TDA) THEN
            IF (mp2_env%bse%screening_method == bse_screening_w0) THEN
               WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'B_ia,jb = α * v_ia,jb - W_ib,aj'
            ELSE IF (mp2_env%bse%screening_method == bse_screening_rpa) THEN
               WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', 'B_ia,jb = α * v_ia,jb'
            END IF
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
         WRITE (unit_nr, '(T2,A4,T7,A7,T31,A23)') 'BSE|', 'ε_...:', 'GW quasiparticle energy'
         WRITE (unit_nr, '(T2,A4,T7,A7,T31,A15)') 'BSE|', 'δ_...:', 'Kronecker delta'
         WRITE (unit_nr, '(T2,A4,T7,A3,T31,A21)') 'BSE|', 'α:', 'spin-dependent factor (Singlet/Triplet)'
         WRITE (unit_nr, '(T2,A4,T7,A6,T30,A34)') 'BSE|', 'v_...:', 'Electron-hole exchange interaction'
         IF (mp2_env%bse%screening_method == bse_screening_w0) THEN
            WRITE (unit_nr, '(T2,A4,T7,A,T31,A)') 'BSE|', 'W_... = 1/ϵ v_...:', &
               'Direct interaction screened by '
            WRITE (unit_nr, '(T2,A4,T30,A)') 'BSE|', &
               'dielectric function ϵ(ω=0)'
         ELSE IF (mp2_env%bse%screening_method == bse_screening_tdhf) THEN
            WRITE (unit_nr, '(T2,A4,T7,A,T30,A)') 'BSE|', 'W_... = v_...:', 'Direct interaction without screening'
         ELSE IF (mp2_env%bse%screening_method == bse_screening_alpha) THEN
            WRITE (unit_nr, '(T2,A4,T7,A,T31,A,F5.2)') 'BSE|', 'W_... = γ v_...:', &
               'Direct interaction with artificial screening γ=', mp2_env%bse%screening_factor
         END IF
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
!> \param homo_irred ...
!> \param flag_TDA ...
!> \param info_approximation ...
!> \param mp2_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_optical_properties(Exc_ens, oscill_str, trans_mom_bse, polarizability_residues, &
                                       homo, virtual, homo_irred, flag_TDA, &
                                       info_approximation, mp2_env, unit_nr)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens, oscill_str
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: trans_mom_bse, polarizability_residues
      INTEGER, INTENT(IN)                                :: homo, virtual, homo_irred
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
               WRITE (unit_nr, '(T2,A4,T10,A)') &
                  'BSE|', "d_r^n = sqrt(2) sum_ia < ψ_i | r | ψ_a >  X_ia^n"
            ELSE
               WRITE (unit_nr, '(T2,A4,T10,A)') &
                  'BSE|', "d_r^n = sum_ia sqrt(2) < ψ_i | r | ψ_a > ( X_ia^n + Y_ia^n )"
            END IF
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T14,A)') &
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
            WRITE (unit_nr, '(T2,A4,T8,A12,T22,A8,T38,A5,T48,A5,T58,A5,T64,A17)') 'BSE|', &
               'Excitation n', "TDA/ABBA", "d_x^n", "d_y^n", "d_z^n", 'Osc. strength f^n'
            DO i_exc = 1, MIN(homo*virtual, mp2_env%bse%num_print_exc)
               WRITE (unit_nr, '(T2,A4,T8,I12,T24,A6,T35,F8.3,T45,F8.3,T55,F8.3,T65,F16.3)') &
                  'BSE|', i_exc, info_approximation, trans_mom_bse(1, 1, i_exc), trans_mom_bse(2, 1, i_exc), &
                  trans_mom_bse(3, 1, i_exc), oscill_str(i_exc)
            END DO
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', &
               'Check for Thomas-Reiche-Kuhn sum rule'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T35,A)') 'BSE|', &
               'N_e = Σ_n f^n'
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            WRITE (unit_nr, '(T2,A4,T7,A24,T65,I16)') 'BSE|', &
               'Number of electrons N_e:', homo_irred*2
            WRITE (unit_nr, '(T2,A4,T7,A,T66,F16.3)') 'BSE|', &
               'Sum over oscillator strengths Σ_n f^n :', SUM(oscill_str)
            WRITE (unit_nr, '(T2,A4)') 'BSE|'
            IF (mp2_env%bse%bse_cutoff_occ > 0 .OR. mp2_env%bse%bse_cutoff_empty > 0) THEN
               CALL cp_warn(__LOCATION__, &
                            "Accuracy of TRK sum rule might suffer from cutoffs.")
            END IF
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
!> \brief                              Prints exciton descriptors, cf. Mewes et al., JCTC 14, 710-725 (2018)
!> \param exc_descr                    Exciton descriptors with size of num_print_exc_descr
!> \param ref_point_multipole          Reference point for computation of multipole moments, e.g. center of mass
!> \param unit_nr ...
!> \param num_print_exc_descr          Number of excitation levels for which descriptors are printed
!> \param print_checkvalue             Flag, which determines if values for regtests should be printed
!> \param print_directional_exc_descr  Flag, which activates printing of directional descriptors
!> \param prefix_output                String, which is put in front of prints, i.e. "BSE|" or "" for TDDFPT
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE print_exciton_descriptors(exc_descr, ref_point_multipole, unit_nr, &
                                        num_print_exc_descr, print_checkvalue, print_directional_exc_descr, &
                                        prefix_output, qs_env)

      TYPE(exciton_descr_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: exc_descr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: ref_point_multipole
      INTEGER, INTENT(IN)                                :: unit_nr, num_print_exc_descr
      LOGICAL, INTENT(IN)                                :: print_checkvalue, &
                                                            print_directional_exc_descr
      CHARACTER(LEN=4), INTENT(IN)                       :: prefix_output
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_exciton_descriptors'

      CHARACTER(LEN=1), DIMENSION(3)                     :: array_direction_str
      CHARACTER(LEN=5)                                   :: method_name
      INTEGER                                            :: handle, i_dir, i_exc
      REAL(KIND=dp)                                      :: d_eh_dir, d_exc_dir, sigma_e_dir, &
                                                            sigma_h_dir
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      IF (prefix_output == 'BSE|') THEN
         method_name = 'BSE'
      ELSE
         method_name = 'TDDFT'
      END IF

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, particle_set=particle_set)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A)') prefix_output, &
            'Exciton descriptors for excitation level n are given by'
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
            'd_eh   = | <r_h - r_e>_exc |'
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
            'σ_e    = sqrt( <r_e^2>_exc - <r_e>_exc^2 )'
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
            'σ_h    = sqrt( <r_h^2>_exc - <r_h>_exc^2 )'
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
            'COV_eh = <r_e r_h>_exc - <r_e>_exc <r_h>_exc'
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
            'd_exc  = sqrt( | < |r_h - r_e|^2 >_exc )'
         WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
            '       = sqrt( d_eh^2 + σ_e^2 + σ_h^2 - 2 * COV_eh )'
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
            'R_eh   = COV_eh / (σ_e * σ_h)'
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T7,A)') prefix_output, &
            'where the expectation values <.>_exc are taken with respect to the '
         WRITE (unit_nr, '(T2,A4,T7,A)') prefix_output, &
            'exciton wavefunction  of excitation n:'
         WRITE (unit_nr, '(T2,A4)') prefix_output

         IF (exc_descr(1)%flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4,T20,A)') prefix_output, &
               '𝚿_n(r_e,r_h) = Σ_{i,a} X_ia^n ψ_i(r_h) ψ_a(r_e)  ,'
         ELSE
            WRITE (unit_nr, '(T2,A4,T20,A)') prefix_output, &
               '𝚿_n(r_e,r_h) = Σ_{i,a} X_ia^n ψ_i(r_h) ψ_a(r_e)'
            WRITE (unit_nr, '(T2,A4,T40,A)') prefix_output, &
               '+ Y_ia^n ψ_a(r_h) ψ_i(r_e) ,'
         END IF
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T7,A)') prefix_output, &
            'i.e.'
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T20,A)') prefix_output, &
            '< O >_exc = < 𝚿_n | O | 𝚿_n > / < 𝚿_n | 𝚿_n >  ,'
         WRITE (unit_nr, '(T2,A4)') prefix_output
         IF (exc_descr(1)%flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4,T7,A)') prefix_output, &
               'where c_n = < 𝚿_n | 𝚿_n > = 1 within TDA.'
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A)') prefix_output, &
               'where c_n = < 𝚿_n | 𝚿_n > deviates from 1 without TDA.'
         END IF
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T7,A)') prefix_output, &
            'Here, we introduced'
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T7,A5,T15,A)') &
            prefix_output, "ψ_i:", "occupied molecular orbitals,"
         WRITE (unit_nr, '(T2,A4,T7,A5,T15,A)') &
            prefix_output, "ψ_a:", "empty molecular orbitals and"
         WRITE (unit_nr, '(T2,A4,T9,A2,T14,A)') &
            prefix_output, "r:", "position operator."
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T7,A)') prefix_output, &
            'prelim Ref.: Eqs. (15)-(22)'
         WRITE (unit_nr, '(T2,A4,T7,A,A)') prefix_output, &
            'JCTC 2018, 14, 710-725; ', &
            'http://doi.org/10.1021/acs.jctc.7b01145'
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4)') prefix_output
         IF (exc_descr(1)%flag_TDA) THEN
            WRITE (unit_nr, '(T2,A4,T7,A,A,A)') prefix_output, &
               'Exciton descriptors from solving the ', method_name, ' within the TDA:'
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A,A,A)') prefix_output, &
               'Exciton descriptors from solving the ', method_name, ' without the TDA:'
         END IF
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T10,A1,6X,A3,1X,4X,A10,5X,A10,5X,A10,3X,A11,8X,A4)') prefix_output, &
            'n', 'c_n', 'd_eh [Å]', 'σ_e [Å]', 'σ_h [Å]', 'd_exc [Å]', 'R_eh'
         DO i_exc = 1, num_print_exc_descr
            WRITE (unit_nr, '(T2,A4,T7,I4,4X,F5.3,1X,5(2X,F10.4))') &
               prefix_output, i_exc, exc_descr(i_exc)%norm_XpY, &
               exc_descr(i_exc)%diff_r_abs*angstrom, &
               exc_descr(i_exc)%sigma_e*angstrom, exc_descr(i_exc)%sigma_h*angstrom, &
               exc_descr(i_exc)%diff_r_sqr*angstrom, exc_descr(i_exc)%corr_e_h*angstrom
         END DO
         WRITE (unit_nr, '(T2,A4)') prefix_output
         ! For debug runs, print first d_exc separately to allow the regtests to read in
         IF (print_checkvalue) THEN
            WRITE (unit_nr, '(T2)')
            WRITE (unit_nr, '(T2,A28,T65,F16.4)') 'Checksum exciton descriptors', &
               exc_descr(1)%diff_r_sqr*angstrom
            WRITE (unit_nr, '(T2)')
         END IF
         WRITE (unit_nr, '(T2,A4)') prefix_output
         ! Print exciton descriptor resolved per direction
         IF (print_directional_exc_descr) THEN
            WRITE (unit_nr, '(T2,A4,T7,A)') prefix_output, &
               'We can restrict the exciton descriptors to a specific direction,'
            WRITE (unit_nr, '(T2,A4,T7,A)') prefix_output, &
               'e.g. the x-components are:'
            WRITE (unit_nr, '(T2,A4)') prefix_output
            WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
               'd_eh^x   = | <x_h - x_e>_exc |'
            WRITE (unit_nr, '(T2,A4)') prefix_output
            WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
               'σ_e^x    = sqrt( <x_e^2>_exc - <x_e>_exc^2 )'
            WRITE (unit_nr, '(T2,A4)') prefix_output
            WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
               'σ_h^x    = sqrt( <x_h^2>_exc - <x_h>_exc^2 )'
            WRITE (unit_nr, '(T2,A4)') prefix_output
            WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
               'COV_eh^x = <x_e x_h>_exc - <x_e>_exc <x_h>_exc'
            WRITE (unit_nr, '(T2,A4)') prefix_output
            WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
               'd_exc^x  = sqrt( | < |x_h - x_e|^2 >_exc )'
            WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
               '         = sqrt( (d_eh^x)^2 + (σ_e^x)^2'
            WRITE (unit_nr, '(T2,A4,T15,A)') prefix_output, &
               '                 + (σ_h^x)^2 - 2 * (COV_eh^x) )'
            WRITE (unit_nr, '(T2,A4)') prefix_output
            WRITE (unit_nr, '(T2,A4)') prefix_output
            IF (exc_descr(1)%flag_TDA) THEN
               WRITE (unit_nr, '(T2,A4,T7,A,A,A)') prefix_output, &
                  'Exciton descriptors per direction from solving the ', method_name, ' within the TDA:'
            ELSE
               WRITE (unit_nr, '(T2,A4,T7,A,A,A)') prefix_output, &
                  'Exciton descriptors per direction from solving the ', method_name, ' without the TDA:'
            END IF
            WRITE (unit_nr, '(T2,A4)') prefix_output
            WRITE (unit_nr, '(T2,A4,T12,A1,2X,A9,5X,A12,5X,A12,5X,A12,3X,A13)') prefix_output, &
               'n', 'r = x/y/z', 'd_eh^r [Å]', 'σ_e^r [Å]', 'σ_h^r [Å]', 'd_exc^r [Å]'
            DO i_exc = 1, num_print_exc_descr
               DO i_dir = 1, 3
                  array_direction_str = (/"x", "y", "z"/)
                  d_eh_dir = ABS(exc_descr(i_exc)%r_h(i_dir) - exc_descr(i_exc)%r_e(i_dir))
                  sigma_e_dir = SQRT(exc_descr(i_exc)%r_e_sq(i_dir) - exc_descr(i_exc)%r_e(i_dir)**2)
                  sigma_h_dir = SQRT(exc_descr(i_exc)%r_h_sq(i_dir) - exc_descr(i_exc)%r_h(i_dir)**2)
                  d_exc_dir = SQRT(d_eh_dir**2 + sigma_e_dir**2 + sigma_h_dir**2 &
                                   - 2*exc_descr(i_exc)%cov_e_h(i_dir))
                  WRITE (unit_nr, '(T2,A4,T9,I4,10X,A1,1X,4(4X,F10.4))') &
                     prefix_output, i_exc, array_direction_str(i_dir), &
                     d_eh_dir*angstrom, &
                     sigma_e_dir*angstrom, sigma_h_dir*angstrom, &
                     d_exc_dir*angstrom
               END DO
               WRITE (unit_nr, '(T2,A4)') prefix_output
            END DO
            WRITE (unit_nr, '(T2,A4)') prefix_output
         END IF
         ! Print the reference atomic geometry for the exciton descriptors
         WRITE (unit_nr, '(T2,A4,T7,A)') prefix_output, &
            'With the center of charge as reference point r_0,'
         WRITE (unit_nr, '(T2,A4,T15,A7,F10.4,A2,F10.4,A2,F10.4,A1)') prefix_output, &
            'r_0 = (', ref_point_multipole(1)*angstrom, ', ', ref_point_multipole(2)*angstrom, ', ', &
            ref_point_multipole(3)*angstrom, ')'
         WRITE (unit_nr, '(T2,A4,T7,A,A,A)') prefix_output, &
            'we further obtain r_e and r_h from solving the ', method_name, ' within the TDA'
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T8,A12,1X,13X,A9,13X,A9,13X,A9)') prefix_output, &
            'Excitation n', 'x_e [Å]', 'y_e [Å]', 'z_e [Å]'
         DO i_exc = 1, num_print_exc_descr
            WRITE (unit_nr, '(T2,A4,T8,I12,1X,3(5X,F15.4))') &
               prefix_output, i_exc, &
               exc_descr(i_exc)%r_e_shift(:)*angstrom
         END DO
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T8,A12,1X,13X,A9,13X,A9,13X,A9)') prefix_output, &
            'Excitation n', 'x_h [Å]', 'y_h [Å]', 'z_h [Å]'
         DO i_exc = 1, num_print_exc_descr
            WRITE (unit_nr, '(T2,A4,T8,I12,1X,3(5X,F15.4))') &
               prefix_output, i_exc, &
               exc_descr(i_exc)%r_h_shift(:)*angstrom
         END DO
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T7,A)') prefix_output, &
            'The reference atomic geometry for these values is given by'
      END IF
      CALL write_qs_particle_coordinates_bse(particle_set, unit_nr, prefix_output)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') prefix_output
      END IF
      CALL timestop(handle)

   END SUBROUTINE print_exciton_descriptors

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

! **************************************************************************************************
!> \brief Write the atomic coordinates to the output unit.
!> \param particle_set ...
!>       \note Adapted from particle_methods.F [MG]
!> \param unit_nr ...
!> \param prefix_output ...
! **************************************************************************************************
   SUBROUTINE write_qs_particle_coordinates_bse(particle_set, unit_nr, prefix_output)

      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(LEN=4), INTENT(IN)                       :: prefix_output

      CHARACTER(len=*), PARAMETER :: routineN = 'write_qs_particle_coordinates_bse'

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: handle, iatom, natom

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') prefix_output
         WRITE (unit_nr, '(T2,A4,T13,A7,16X,A7,15X,A7,15X,A7)') prefix_output, &
            'Element', 'x [Å]', 'y [Å]', 'z [Å]'
         natom = SIZE(particle_set)
         DO iatom = 1, natom
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                 element_symbol=element_symbol)
            WRITE (unit_nr, '(T2,A4,T8,A12,1X,3(5X,F15.4))') &
               prefix_output, element_symbol, particle_set(iatom)%r(1:3)*angstrom
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_qs_particle_coordinates_bse

END MODULE bse_print
