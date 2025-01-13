!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for computing excitonic properties, e.g. exciton diameter, from the BSE
!> \par History
!>      10.2024 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_properties
   USE bse_util,                        ONLY: fm_general_add_bse,&
                                              print_bse_nto_cubes,&
                                              reshuffle_eigvec,&
                                              trace_exciton_descr
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_trace
   USE cp_fm_diag,                      ONLY: cp_fm_svd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE mp2_types,                       ONLY: mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: c_light_au,&
                                              evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              init_mo_set,&
                                              mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_properties'

   PUBLIC :: exciton_descr_type

   PUBLIC :: get_exciton_descriptors, get_oscillator_strengths, compute_absorption_spectrum, &
             calculate_NTOs

! TYPE definitions for exciton wavefunction descriptors

   TYPE exciton_descr_type
      REAL(KIND=dp), DIMENSION(3)                        :: r_e = 0.0_dp, &
                                                            r_h = 0.0_dp, &
                                                            r_e_sq = 0.0_dp, &
                                                            r_h_sq = 0.0_dp, &
                                                            r_e_h = 0.0_dp, &
                                                            r_e_shift = 0.0_dp, &
                                                            r_h_shift = 0.0_dp, &
                                                            cov_e_h = 0.0_dp
      REAL(KIND=dp)                                      :: sigma_e = 0.0_dp, &
                                                            sigma_h = 0.0_dp, &
                                                            cov_e_h_sum = 0.0_dp, &
                                                            corr_e_h = 0.0_dp, &
                                                            diff_r_abs = 0.0_dp, &
                                                            diff_r_sqr = 0.0_dp, &
                                                            norm_XpY = 0.0_dp
      LOGICAL                                           :: flag_TDA = .FALSE.
   END TYPE exciton_descr_type

CONTAINS

! **************************************************************************************************
!> \brief Compute and return BSE dipoles d_r^n = sqrt(2) sum_ia < ψ_i | r | ψ_a > ( X_ia^n + Y_ia^n )
!>    and oscillator strengths f^n = 2/3 * Ω^n sum_r∈(x,y,z) ( d_r^n )^2
!>    Prelim Ref.: Eqs. (23), (24)
!>    in J. Chem. Phys. 152, 044105 (2020); https://doi.org/10.1063/1.5123290
!> \param fm_eigvec_X ...
!> \param Exc_ens ...
!> \param fm_dipole_ai_trunc ...
!> \param trans_mom_bse BSE dipole vectors in real space per excitation level
!> \param oscill_str Oscillator strength per excitation level
!> \param polarizability_residues Residues of polarizability ("tensorial oscillator strength")
!>          per excitation level
!> \param mp2_env ...
!> \param homo_red ...
!> \param virtual_red ...
!> \param unit_nr ...
!> \param fm_eigvec_Y ...
! **************************************************************************************************
   SUBROUTINE get_oscillator_strengths(fm_eigvec_X, Exc_ens, fm_dipole_ai_trunc, &
                                       trans_mom_bse, oscill_str, polarizability_residues, &
                                       mp2_env, homo_red, virtual_red, unit_nr, &
                                       fm_eigvec_Y)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_eigvec_X
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: Exc_ens
      TYPE(cp_fm_type), DIMENSION(3)                     :: fm_dipole_ai_trunc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: trans_mom_bse
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: oscill_str
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: polarizability_residues
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: homo_red, virtual_red, unit_nr
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: fm_eigvec_Y

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_oscillator_strengths'

      INTEGER                                            :: handle, idir, jdir, n
      TYPE(cp_fm_struct_type), POINTER :: fm_struct_dipole_MO_trunc_reordered, &
         fm_struct_trans_mom_bse
      TYPE(cp_fm_type)                                   :: fm_eigvec_XYsum
      TYPE(cp_fm_type), DIMENSION(3)                     :: fm_dipole_MO_trunc_reordered, &
                                                            fm_dipole_per_dir, fm_trans_mom_bse

      CALL timeset(routineN, handle)

      CALL cp_fm_struct_create(fm_struct_dipole_MO_trunc_reordered, fm_eigvec_X%matrix_struct%para_env, &
                               fm_eigvec_X%matrix_struct%context, 1, homo_red*virtual_red)
      CALL cp_fm_struct_create(fm_struct_trans_mom_bse, fm_eigvec_X%matrix_struct%para_env, &
                               fm_eigvec_X%matrix_struct%context, 1, homo_red*virtual_red)

      ! Include excitonic amplitudes in dipoles, i.e. obtain "BSE dipoles":
      ! \vec{D}_n = sqrt(2) * sum_{i,a} \vec{D}_ai (X_{ai}^{(n)} + Y_{ai}^{(n)})

      ! Reorder dipoles in order to execute the sum over i and a by parallel gemm
      DO idir = 1, 3
         CALL cp_fm_create(fm_dipole_MO_trunc_reordered(idir), matrix_struct=fm_struct_dipole_MO_trunc_reordered, &
                           name="dipoles_mo_reordered")
         CALL cp_fm_set_all(fm_dipole_MO_trunc_reordered(idir), 0.0_dp)
         CALL fm_general_add_bse(fm_dipole_MO_trunc_reordered(idir), fm_dipole_ai_trunc(idir), 1.0_dp, &
                                 1, 1, &
                                 1, virtual_red, &
                                 unit_nr, (/2, 4, 3, 1/), mp2_env)
         CALL cp_fm_release(fm_dipole_per_dir(idir))
      END DO

      DO idir = 1, 3
         CALL cp_fm_create(fm_trans_mom_bse(idir), matrix_struct=fm_struct_trans_mom_bse, &
                           name="excitonic_dipoles")
         CALL cp_fm_set_all(fm_trans_mom_bse(idir), 0.0_dp)
      END DO

      ! If TDA is invoked, Y is not present as it is simply 0
      CALL cp_fm_create(fm_eigvec_XYsum, matrix_struct=fm_eigvec_X%matrix_struct, name="excit_amplitude_sum")
      CALL cp_fm_set_all(fm_eigvec_XYsum, 0.0_dp)
      CALL cp_fm_to_fm(fm_eigvec_X, fm_eigvec_XYsum)
      IF (PRESENT(fm_eigvec_Y)) THEN
         CALL cp_fm_scale_and_add(1.0_dp, fm_eigvec_XYsum, 1.0_dp, fm_eigvec_Y)
      END IF
      DO idir = 1, 3
         CALL parallel_gemm('N', 'N', 1, homo_red*virtual_red, homo_red*virtual_red, SQRT(2.0_dp), &
                            fm_dipole_MO_trunc_reordered(idir), fm_eigvec_XYsum, 0.0_dp, fm_trans_mom_bse(idir))
      END DO

      ! Get oscillator strengths themselves
      ALLOCATE (oscill_str(homo_red*virtual_red))
      ALLOCATE (trans_mom_bse(3, 1, homo_red*virtual_red))
      ALLOCATE (polarizability_residues(3, 3, homo_red*virtual_red))
      trans_mom_bse(:, :, :) = 0.0_dp

      ! Sum over all directions
      DO idir = 1, 3
         CALL cp_fm_get_submatrix(fm_trans_mom_bse(idir), trans_mom_bse(idir, :, :))
      END DO

      DO n = 1, homo_red*virtual_red
         DO idir = 1, 3
            DO jdir = 1, 3
               polarizability_residues(idir, jdir, n) = 2.0_dp*Exc_ens(n)*trans_mom_bse(idir, 1, n)*trans_mom_bse(jdir, 1, n)
            END DO
         END DO
         oscill_str(n) = 2.0_dp/3.0_dp*Exc_ens(n)*SUM(ABS(trans_mom_bse(:, :, n))**2)
      END DO

      CALL cp_fm_struct_release(fm_struct_dipole_MO_trunc_reordered)
      CALL cp_fm_struct_release(fm_struct_trans_mom_bse)
      DO idir = 1, 3
         CALL cp_fm_release(fm_dipole_MO_trunc_reordered(idir))
         CALL cp_fm_release(fm_trans_mom_bse(idir))
         CALL cp_fm_release(fm_dipole_ai_trunc(idir))
      END DO
      CALL cp_fm_release(fm_eigvec_XYsum)

      CALL timestop(handle)

   END SUBROUTINE get_oscillator_strengths

! **************************************************************************************************
!> \brief Computes and returns absorption spectrum for the frequency range and broadening
!>    provided by the user.
!>    Prelim Ref.: C. Ullrich, Time-Dependent Density-Functional Theory: Concepts and Applications
!>    (Oxford University Press, Oxford, 2012), Eq. 7.51
!> \param oscill_str ...
!> \param polarizability_residues ...
!> \param Exc_ens ...
!> \param info_approximation ...
!> \param unit_nr ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE compute_absorption_spectrum(oscill_str, polarizability_residues, Exc_ens, &
                                          info_approximation, unit_nr, mp2_env)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: oscill_str
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: polarizability_residues
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: Exc_ens
      CHARACTER(LEN=10)                                  :: info_approximation
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_absorption_spectrum'

      CHARACTER(LEN=10)                                  :: eta_str, width_eta_format_str
      CHARACTER(LEN=40)                                  :: file_name_crosssection, &
                                                            file_name_spectrum
      INTEGER                                            :: handle, i, idir, j, jdir, k, num_steps, &
                                                            unit_nr_file, width_eta
      REAL(KIND=dp)                                      :: eta, freq_end, freq_start, freq_step, &
                                                            omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: abs_cross_section, abs_spectrum
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eta_list
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      freq_step = mp2_env%bse%bse_spectrum_freq_step_size
      freq_start = mp2_env%bse%bse_spectrum_freq_start
      freq_end = mp2_env%bse%bse_spectrum_freq_end
      eta_list => mp2_env%bse%bse_eta_spectrum_list

      ! Calculate number of steps to fit given frequency range
      num_steps = NINT((freq_end - freq_start)/freq_step) + 1

      DO k = 1, SIZE(eta_list)
         eta = eta_list(k)

         ! Some magic to get a nice formatting of the eta value in filenames
         width_eta = MAX(1, INT(LOG10(eta)) + 1) + 4
         WRITE (width_eta_format_str, "(A2,I0,A3)") '(F', width_eta, '.3)'
         WRITE (eta_str, width_eta_format_str) eta*evolt
         ! Filename itself
         file_name_spectrum = 'BSE'//TRIM(ADJUSTL(info_approximation))//'eta='//TRIM(eta_str)//'.spectrum'
         file_name_crosssection = 'BSE'//TRIM(ADJUSTL(info_approximation))//'eta='//TRIM(eta_str)//'.crosssection'

         ! First column is frequency in eV, second column is imaginary part of the trace of the polarizability
         ! The following 9 columns are the entries of the polarizability tensor
         ALLOCATE (abs_spectrum(num_steps, 11))
         abs_spectrum(:, :) = 0.0_dp
         ! Also calculate and print the photoabsorption cross section tensor
         ! σ_{µ µ'}(ω) = 4πω Im{α_{µ µ'}(ω)} / c
         ALLOCATE (abs_cross_section(num_steps, 11))
         abs_cross_section(:, :) = 0.0_dp

         ! Calculate the imaginary part of the mean dipole polarizability α_{avg}(ω)
         ! which is given by (cf. C. Ullrichs Book on TDDFT, Eq. 7.51)
         ! We introduce an additional - due to his convention for charge vs particle density, see also:
         ! Computer Physics Communications, 208:149–161, November 2016
         ! https://doi.org/10.1016/j.cpc.2016.06.019
         ! α_{avg}(ω) = - \sum_{n=1}^{N_exc} \frac{f_n}{(ω+iη)² - (Ω^n)²}
         ! and then the imaginary part is (in the limit η -> 0)
         ! Im[α_{avg}(ω)] = - \sum_{n=1}^{N_exc} f_n * η / ((ω - Ω^n)² + η²)
         ! where f_n are the oscillator strengths and E_exc the excitation energies
         ! For the full polarizability tensor, we have
         ! α_{µ µ'}(ω) = - \sum_n [2 Ω^n d^n_µ d^n_µ'] / [(ω+iη)^2- (Ω^n)^2]
         !             = - \sum_n "polarizability_residues" / [(ω+iη)^2- (Ω^n)^2]
         DO i = 1, num_steps
            omega = freq_start + (i - 1)*freq_step
            abs_spectrum(i, 1) = omega
            DO j = 1, SIZE(oscill_str)
               abs_spectrum(i, 2) = abs_spectrum(i, 2) - oscill_str(j)* &
                                    AIMAG(1/((omega + CMPLX(0.0, eta, kind=dp))**2 - Exc_ens(j)**2))
               DO idir = 1, 3
                  DO jdir = 1, 3
                     ! Factor 2 from formula for tensor is already in the polarizability_residues
                     !  to follow the same convention as the oscillator strengths
                     abs_spectrum(i, 2 + (idir - 1)*3 + jdir) = abs_spectrum(i, 2 + (idir - 1)*3 + jdir) &
                                                                - polarizability_residues(idir, jdir, j)* &
                                                                AIMAG(1/((omega + CMPLX(0.0, eta, kind=dp))**2 - Exc_ens(j)**2))
                  END DO
               END DO
            END DO
         END DO

         ! Extract cross section σ from polarizability tensor
         DO i = 1, num_steps
            omega = abs_spectrum(i, 1)
            abs_cross_section(i, 1) = omega
            abs_cross_section(i, 2:) = 4.0_dp*pi*abs_spectrum(i, 2:)*omega/c_light_au
         END DO

         !For debug runs: Export an entry of the two tensors to allow regtests on spectra
         IF (mp2_env%bse%bse_debug_print) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A10,T13,A,T65,F16.4)') 'BSE|DEBUG|', &
                  'Averaged dynamical dipole polarizability at 8.2 eV:', &
                  abs_spectrum(83, 2)
               WRITE (unit_nr, '(T2,A10,T13,A,T65,F16.4)') 'BSE|DEBUG|', &
                  'Averaged photoabsorption cross section at 8.2 eV:', &
                  abs_cross_section(83, 2)
            END IF
         END IF

         ! Print it to file
         logger => cp_get_default_logger()
         IF (logger%para_env%is_source()) THEN
            unit_nr_file = cp_logger_get_default_unit_nr()
         ELSE
            unit_nr_file = -1
         END IF

         IF (unit_nr_file > 0) THEN
            CALL open_file(file_name_crosssection, unit_number=unit_nr_file, &
                           file_status="UNKNOWN", file_action="WRITE")
            WRITE (unit_nr_file, '(A,A6)') "# Photoabsorption cross section  σ_{µ µ'}(ω) =  -4πω/c * Im[ \sum_n "// &
               "[2 Ω^n d_µ^n d_µ'^n] / [(ω+iη)²- (Ω^n)²] ] from Bethe Salpeter equation for method ", &
               TRIM(ADJUSTL(info_approximation))
            WRITE (unit_nr_file, '(A20,1X,10(2X,A20,1X))') "#     Frequency (eV)", "σ_{avg}(ω)", "σ_xx(ω)", &
               "σ_xy(ω)", "σ_xz(ω)", "σ_yx(ω)", "σ_yy(ω)", "σ_yz(ω)", "σ_zx(ω)", &
               "σ_zy(ω)", "σ_zz(ω)"
            DO i = 1, num_steps
               WRITE (unit_nr_file, '(11(F20.8,1X))') abs_cross_section(i, 1)*evolt, abs_cross_section(i, 2:11)
            END DO
            CALL close_file(unit_nr_file)
         END IF
         DEALLOCATE (abs_cross_section)

         IF (unit_nr_file > 0) THEN
            CALL open_file(file_name_spectrum, unit_number=unit_nr_file, &
                           file_status="UNKNOWN", file_action="WRITE")
            WRITE (unit_nr_file, '(A,A6)') "# Imaginary part of polarizability α_{µ µ'}(ω) = -\sum_n "// &
               "[2 Ω^n d_µ^n d_µ'^n] / [(ω+iη)²- (Ω^n)²] from Bethe Salpeter equation for method ", &
               TRIM(ADJUSTL(info_approximation))
            WRITE (unit_nr_file, '(A20,1X,10(2X,A20,1X))') "#     Frequency (eV)", "Im{α_{avg}(ω)}", "Im{α_xx(ω)}", &
               "Im{α_xy(ω)}", "Im{α_xz(ω)}", "Im{α_yx(ω)}", "Im{α_yy(ω)}", "Im{α_yz(ω)}", "Im{α_zx(ω)}", &
               "Im{α_zy(ω)}", "Im{α_zz(ω)}"
            DO i = 1, num_steps
               WRITE (unit_nr_file, '(11(F20.8,1X))') abs_spectrum(i, 1)*evolt, abs_spectrum(i, 2:11)
            END DO
            CALL close_file(unit_nr_file)
         END IF
         DEALLOCATE (abs_spectrum)
      END DO

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A,A)') &
            'BSE|', "Printed optical absorption spectrum to local files, e.g. ", file_name_spectrum
         WRITE (unit_nr, '(T2,A4,T7,A,A)') &
            'BSE|', "as well as photoabsorption cross section to, e.g. ", file_name_crosssection
         WRITE (unit_nr, '(T2,A4,T7,A52)') &
            'BSE|', "using the Eq. (7.51) from C. Ullrichs Book on TDDFT:"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T10,A73)') &
            'BSE|', "Im{α_{avg}(ω)} = -Im{\sum_{n=1}^{N_exc} \frac{f_n}{(ω+iη)² - (Ω^n)²}}"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A)') &
            'BSE|', "or for the full polarizability tensor:"
         WRITE (unit_nr, '(T2,A4,T10,A)') &
            'BSE|', "α_{µ µ'}(ω) =  -\sum_n [2 Ω^n d_µ^n d_µ'^n] / [(ω+iη)²- (Ω^n)²]"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A)') &
            'BSE|', "as well as Eq. (7.48):"
         WRITE (unit_nr, '(T2,A4,T10,A)') &
            'BSE|', "σ_{µ µ'}(ω) = 4πω Im{α_{µ µ'}(ω)} / c"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A)') &
            'BSE|', "with transition moments d_µ^n, oscillator strengths f_n,"
         WRITE (unit_nr, '(T2,A4,T7,A)') &
            'BSE|', "excitation energies Ω^n and the speed of light c."
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A)') &
            'BSE|', "Please note that we adopt an additional minus sign for both quantities,"
         WRITE (unit_nr, '(T2,A4,T7,A)') &
            'BSE|', "due to the convention for charge vs particle density as done in MolGW:"
         WRITE (unit_nr, '(T2,A4,T7,A)') &
            'BSE|', "https://doi.org/10.1016/j.cpc.2016.06.019."
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param fm_X ...
!> \param fm_Y ...
!> \param mo_coeff ...
!> \param homo ...
!> \param virtual ...
!> \param info_approximation ...
!> \param oscill_str ...
!> \param qs_env ...
!> \param unit_nr ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE calculate_NTOs(fm_X, fm_Y, &
                             mo_coeff, homo, virtual, &
                             info_approximation, &
                             oscill_str, &
                             qs_env, unit_nr, mp2_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_X
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: fm_Y
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      INTEGER, INTENT(IN)                                :: homo, virtual
      CHARACTER(LEN=10)                                  :: info_approximation
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: oscill_str
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calculate_NTOs'

      CHARACTER(LEN=20), DIMENSION(2)                    :: nto_name
      INTEGER                                            :: handle, homo_irred, i, i_nto, info_svd, &
                                                            j, n_exc, n_nto, nao_full, nao_trunc
      INTEGER, DIMENSION(:), POINTER                     :: stride
      LOGICAL                                            :: append_cube, cube_file
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigval_svd_squ
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigval_svd
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_m, fm_struct_mo_coeff, &
                                                            fm_struct_nto_holes, &
                                                            fm_struct_nto_particles, &
                                                            fm_struct_nto_set
      TYPE(cp_fm_type) :: fm_eigvl, fm_eigvr_t, fm_m, fm_mo_coeff, fm_nto_coeff_holes, &
         fm_nto_coeff_particles, fm_nto_set, fm_X_ia, fm_Y_ai
      TYPE(mo_set_type), ALLOCATABLE, DIMENSION(:)       :: nto_set
      TYPE(section_vals_type), POINTER                   :: bse_section, input, nto_section

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env=qs_env, &
                      input=input)
      bse_section => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%GW%BSE")

      nao_full = qs_env%mos(1)%nao
      nao_trunc = homo + virtual
      ! This is not influenced by the BSE cutoff
      homo_irred = qs_env%mos(1)%homo
      ! M will have a block structure and is quadratic in homo+virtual, i.e.
      !                          occ   virt
      !                       |   0    X_i,a |  occ  = homo
      !     M        =        | Y_a,i    0   |  virt = virtual
      !
      ! X and Y are here not the eigenvectors X_ia,n - instead we fix n and reshape the combined ia index
      ! Notice the index structure of the lower block, i.e. X is transposed
      CALL cp_fm_struct_create(fm_struct_m, &
                               fm_X%matrix_struct%para_env, fm_X%matrix_struct%context, &
                               nao_trunc, nao_trunc)
      CALL cp_fm_struct_create(fm_struct_mo_coeff, &
                               fm_X%matrix_struct%para_env, fm_X%matrix_struct%context, &
                               nao_full, nao_trunc)
      CALL cp_fm_struct_create(fm_struct_nto_holes, &
                               fm_X%matrix_struct%para_env, fm_X%matrix_struct%context, &
                               nao_full, nao_trunc)
      CALL cp_fm_struct_create(fm_struct_nto_particles, &
                               fm_X%matrix_struct%para_env, fm_X%matrix_struct%context, &
                               nao_full, nao_trunc)

      CALL cp_fm_create(fm_mo_coeff, matrix_struct=fm_struct_mo_coeff, &
                        name="mo_coeff")
      ! Here, we take care of possible cutoffs
      ! Simply truncating the matrix causes problems with the print function
      ! Therefore, we keep the dimension, but set the coefficients of truncated indices to 0
      CALL cp_fm_to_fm_submat_general(mo_coeff(1), fm_mo_coeff, &
                                      nao_full, nao_trunc, &
                                      1, homo_irred - homo + 1, &
                                      1, 1, &
                                      mo_coeff(1)%matrix_struct%context)

      ! Print some information about the NTOs
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', &
            'The Natural Transition Orbital (NTO) pairs φ_I(r_e) and χ_I(r_h) for a fixed'
         WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', &
            'excitation index n are obtained by singular value decomposition of T'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T15,A)') 'BSE|', &
            '        = (0   X)'
         IF (PRESENT(fm_Y)) THEN
            WRITE (unit_nr, '(T2,A4,T15,A)') 'BSE|', &
               'T       = (Y^T 0)'
         ELSE
            WRITE (unit_nr, '(T2,A4,T15,A)') 'BSE|', &
               'T       = (0   0)'
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T15,A)') 'BSE|', &
            'T        = U Λ V^T'
         WRITE (unit_nr, '(T2,A4,T15,A)') 'BSE|', &
            'φ_I(r_e) = \sum_p V_pI ψ_p(r_e)'
         WRITE (unit_nr, '(T2,A4,T15,A)') 'BSE|', &
            'χ_I(r_h) = \sum_p U_pI ψ_p(r_e)'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', &
            'where we have introduced'
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A,T20,A)') &
            'BSE|', "ψ_p:", "occupied and virtual molecular orbitals,"
         WRITE (unit_nr, '(T2,A4,T7,A,T20,A)') &
            'BSE|', "φ_I(r_e):", "NTO state for the electron,"
         WRITE (unit_nr, '(T2,A4,T7,A,T20,A)') &
            'BSE|', "χ_I(r_h):", "NTO state for the hole,"
         WRITE (unit_nr, '(T2,A4,T7,A,T20,A)') &
            'BSE|', "Λ:", "diagonal matrix of NTO weights λ_I,"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', &
            "The NTOs are calculated with the following settings:"
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'Number of excitations, for which NTOs are computed', &
            mp2_env%bse%num_print_exc_ntos
         IF (mp2_env%bse%eps_nto_osc_str > 0.0_dp) THEN
            WRITE (unit_nr, '(T2,A4,T7,A,T71,F10.3)') 'BSE|', 'Threshold for oscillator strength f^n', &
               mp2_env%bse%eps_nto_osc_str
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A,T71,A10)') 'BSE|', 'Threshold for oscillator strength f^n', &
               ADJUSTL("---")
         END IF
         WRITE (unit_nr, '(T2,A4,T7,A,T72,F10.3)') 'BSE|', 'Threshold for NTO weights (λ_I)^2', &
            mp2_env%bse%eps_nto_eigval
      END IF

      ! Write the header of NTO info table
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         IF (.NOT. PRESENT(fm_Y)) THEN
            WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', &
               'NTOs from solving the BSE within the TDA:'
         ELSE
            WRITE (unit_nr, '(T2,A4,T7,A)') 'BSE|', &
               'NTOs from solving the BSE without the TDA:'
         END IF
         WRITE (unit_nr, '(T2,A4)') 'BSE|'
         WRITE (unit_nr, '(T2,A4,T8,A12,T22,A8,T33,A14,T62,A)') 'BSE|', &
            'Excitation n', "TDA/ABBA", "Index of NTO I", 'NTO weights (λ_I)^2'
      END IF

      DO j = 1, mp2_env%bse%num_print_exc_ntos
         n_exc = mp2_env%bse%bse_nto_state_list_final(j)
         ! Takes care of unallocated oscill_str array in case of Triplet
         IF (mp2_env%bse%eps_nto_osc_str > 0.0_dp) THEN
            ! Check actual values
            IF (oscill_str(n_exc) < mp2_env%bse%eps_nto_osc_str) THEN
               ! Print skipped levels to table
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, '(T2,A4)') 'BSE|'
                  WRITE (unit_nr, '(T2,A4,T8,I12,T24,A6,T42,A39)') 'BSE|', &
                     n_exc, info_approximation, "Skipped (Oscillator strength too small)"
               END IF
               CYCLE
            END IF
         END IF

         CALL cp_fm_create(fm_m, matrix_struct=fm_struct_m, &
                           name="single_part_trans_dm")
         CALL cp_fm_set_all(fm_m, 0.0_dp)

         CALL cp_fm_create(fm_nto_coeff_holes, matrix_struct=fm_struct_nto_holes, &
                           name="nto_coeffs_holes")
         CALL cp_fm_set_all(fm_nto_coeff_holes, 0.0_dp)

         CALL cp_fm_create(fm_nto_coeff_particles, matrix_struct=fm_struct_nto_particles, &
                           name="nto_coeffs_particles")
         CALL cp_fm_set_all(fm_nto_coeff_particles, 0.0_dp)

         ! Reshuffle from X_ia,n_exc to X_i,a
         CALL reshuffle_eigvec(fm_X, fm_X_ia, homo, virtual, n_exc, &
                               .FALSE., unit_nr, mp2_env)

         ! Copy X to upper block in M, i.e. starting from column homo+1
         CALL cp_fm_to_fm_submat(fm_X_ia, fm_m, &
                                 homo, virtual, &
                                 1, 1, &
                                 1, homo + 1)
         CALL cp_fm_release(fm_X_ia)
         ! Copy Y if present
         IF (PRESENT(fm_Y)) THEN
            ! Reshuffle from Y_ia,n_exc to Y_a,i
            CALL reshuffle_eigvec(fm_Y, fm_Y_ai, homo, virtual, n_exc, &
                                  .TRUE., unit_nr, mp2_env)

            ! Copy Y^T to lower block in M, i.e. starting from row homo+1
            CALL cp_fm_to_fm_submat(fm_Y_ai, fm_m, &
                                    virtual, homo, &
                                    1, 1, &
                                    homo + 1, 1)

            CALL cp_fm_release(fm_Y_ai)

         END IF

         ! Now we compute the SVD of M_{occ+virt,occ+virt}, which yields
         ! M = U * Lambda * V^T
         ! Initialize matrices and arrays to store left/right eigenvectors and singular values
         CALL cp_fm_create(matrix=fm_eigvl, &
                           matrix_struct=fm_m%matrix_struct, &
                           name="LEFT_SINGULAR_MATRIX")
         CALL cp_fm_set_all(fm_eigvl, alpha=0.0_dp)
         CALL cp_fm_create(matrix=fm_eigvr_t, &
                           matrix_struct=fm_m%matrix_struct, &
                           name="RIGHT_SINGULAR_MATRIX")
         CALL cp_fm_set_all(fm_eigvr_t, alpha=0.0_dp)

         ALLOCATE (eigval_svd(nao_trunc))
         eigval_svd(:) = 0.0_dp
         info_svd = 0
         CALL cp_fm_svd(fm_m, fm_eigvl, fm_eigvr_t, eigval_svd, info_svd)
         IF (info_svd /= 0) THEN
            IF (unit_nr > 0) THEN
               CALL cp_warn(__LOCATION__, &
                            "SVD for computation of NTOs not successful. "// &
                            "Skipping print of NTOs.")
               IF (info_svd > 0) THEN
                  CALL cp_warn(__LOCATION__, &
                               "PDGESVD detected heterogeneity. "// &
                               "Decreasing number of MPI ranks might solve this issue.")
               END IF
            END IF
            ! Release matrices to avoid memory leaks
            CALL cp_fm_release(fm_m)
            CALL cp_fm_release(fm_nto_coeff_holes)
            CALL cp_fm_release(fm_nto_coeff_particles)
         ELSE
            ! Rescale singular values as done in Martin2003 (10.1063/1.1558471)
            ALLOCATE (eigval_svd_squ(nao_trunc))
            eigval_svd_squ(:) = eigval_svd(:)**2
            ! Sanity check for TDA: In case of TDA, the sum should be \sum_ia |X_ia|^2 = 1
            IF (.NOT. PRESENT(fm_Y)) THEN
               IF (ABS(SUM(eigval_svd_squ) - 1) >= 1e-5) THEN
                  CPWARN("Sum of NTO coefficients deviates from 1!")
               END IF
            END IF

            ! Create NTO coefficients for later print to grid via TDDFT routine
            ! Apply U = fm_eigvl to MO coeffs, which yields hole states
            CALL parallel_gemm("N", "N", nao_full, nao_trunc, nao_trunc, 1.0_dp, fm_mo_coeff, fm_eigvl, 0.0_dp, &
                               fm_nto_coeff_holes)

            ! Apply V^T = fm_eigvr_t to MO coeffs, which yields particle states
            CALL parallel_gemm("N", "T", nao_full, nao_trunc, nao_trunc, 1.0_dp, fm_mo_coeff, fm_eigvr_t, 0.0_dp, &
                               fm_nto_coeff_particles)

            !Release intermediary work matrices
            CALL cp_fm_release(fm_m)
            CALL cp_fm_release(fm_eigvl)
            CALL cp_fm_release(fm_eigvr_t)

            ! Transfer NTO coefficients to sets
            nto_name(1) = 'Hole_coord'
            nto_name(2) = 'Particle_coord'
            ALLOCATE (nto_set(2))
            ! Extract number of significant NTOs
            n_nto = 0
            DO i_nto = 1, nao_trunc
               IF (eigval_svd_squ(i_nto) > mp2_env%bse%eps_nto_eigval) THEN
                  n_nto = n_nto + 1
               ELSE
                  ! Since svd orders in descending order, we can exit the loop if smaller
                  EXIT
               END IF
            END DO

            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A4)') 'BSE|'
               DO i_nto = 1, n_nto
                  WRITE (unit_nr, '(T2,A4,T8,I12,T24,A6,T41,I6,T71,F10.5)') 'BSE|', &
                     n_exc, info_approximation, i_nto, eigval_svd_squ(i_nto)
               END DO
            END IF

            CALL cp_fm_struct_create(fm_struct_nto_set, template_fmstruct=fm_struct_nto_holes, &
                                     ncol_global=n_nto)
            CALL cp_fm_create(fm_nto_set, fm_struct_nto_set)
            DO i = 1, 2
               CALL allocate_mo_set(nto_set(i), nao_trunc, n_nto, 0, 0.0_dp, 2.0_dp, 0.0_dp)
               CALL init_mo_set(nto_set(i), fm_ref=fm_nto_set, name=nto_name(i))
            END DO
            CALL cp_fm_release(fm_nto_set)
            CALL cp_fm_struct_release(fm_struct_nto_set)

            ! Fill NTO sets
            CALL cp_fm_to_fm(fm_nto_coeff_holes, nto_set(1)%mo_coeff, ncol=n_nto)
            CALL cp_fm_to_fm(fm_nto_coeff_particles, nto_set(2)%mo_coeff, ncol=n_nto)

            ! Cube files
            nto_section => section_vals_get_subs_vals(bse_section, "NTO_ANALYSIS")
            CALL section_vals_val_get(nto_section, "CUBE_FILES", l_val=cube_file)
            CALL section_vals_val_get(nto_section, "STRIDE", i_vals=stride)
            CALL section_vals_val_get(nto_section, "APPEND", l_val=append_cube)
            IF (cube_file) THEN
               CALL print_bse_nto_cubes(qs_env, nto_set, n_exc, info_approximation, &
                                        stride, append_cube, nto_section)
            END IF

            CALL cp_fm_release(fm_nto_coeff_holes)
            CALL cp_fm_release(fm_nto_coeff_particles)
            DEALLOCATE (eigval_svd)
            DEALLOCATE (eigval_svd_squ)
            DO i = 1, 2
               CALL deallocate_mo_set(nto_set(i))
            END DO
            DEALLOCATE (nto_set)
         END IF
      END DO

      CALL cp_fm_release(fm_mo_coeff)
      CALL cp_fm_struct_release(fm_struct_m)
      CALL cp_fm_struct_release(fm_struct_nto_holes)
      CALL cp_fm_struct_release(fm_struct_nto_particles)
      CALL cp_fm_struct_release(fm_struct_mo_coeff)

      CALL timestop(handle)

   END SUBROUTINE calculate_NTOs

! **************************************************************************************************
!> \brief ...
!> \param exc_descr Allocated and initialized on exit
!> \param fm_X_ia ...
!> \param fm_multipole_ij_trunc ...
!> \param fm_multipole_ab_trunc ...
!> \param fm_multipole_ai_trunc ...
!> \param i_exc ...
!> \param homo ...
!> \param virtual ...
!> \param fm_Y_ia ...
! **************************************************************************************************
   SUBROUTINE get_exciton_descriptors(exc_descr, fm_X_ia, &
                                      fm_multipole_ij_trunc, fm_multipole_ab_trunc, &
                                      fm_multipole_ai_trunc, &
                                      i_exc, homo, virtual, &
                                      fm_Y_ia)

      TYPE(exciton_descr_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: exc_descr
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_X_ia
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: fm_multipole_ij_trunc, &
                                                            fm_multipole_ab_trunc, &
                                                            fm_multipole_ai_trunc
      INTEGER, INTENT(IN)                                :: i_exc, homo, virtual
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: fm_Y_ia

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_exciton_descriptors'

      INTEGER                                            :: handle, i_dir
      INTEGER, DIMENSION(3)                              :: mask_quadrupole
      LOGICAL                                            :: flag_TDA
      REAL(KIND=dp)                                      :: norm_X, norm_XpY, norm_Y
      REAL(KIND=dp), DIMENSION(3)                        :: r_e_h_XX, r_e_h_XY, r_e_h_YX, r_e_h_YY, &
                                                            r_e_sq_X, r_e_sq_Y, r_e_X, r_e_Y, &
                                                            r_h_sq_X, r_h_sq_Y, r_h_X, r_h_Y
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_ab, fm_struct_ia
      TYPE(cp_fm_type)                                   :: fm_work_ba, fm_work_ia, fm_work_ia_2

      CALL timeset(routineN, handle)
      IF (PRESENT(fm_Y_ia)) THEN
         flag_TDA = .FALSE.
      ELSE
         flag_TDA = .TRUE.
      END IF

      ! translates 1,2,3 to diagonal entries of quadrupoles xx, yy, zz
      ! Ordering in quadrupole moments is x, y, z, xx, xy, xz, yy, yz, zz
      mask_quadrupole = (/4, 7, 9/)

      CALL cp_fm_struct_create(fm_struct_ia, &
                               context=fm_X_ia%matrix_struct%context, nrow_global=homo, ncol_global=virtual)
      CALL cp_fm_struct_create(fm_struct_ab, &
                               context=fm_X_ia%matrix_struct%context, nrow_global=virtual, ncol_global=virtual)

      r_e_X(:) = 0.0_dp
      r_e_Y(:) = 0.0_dp
      r_h_X(:) = 0.0_dp
      r_h_Y(:) = 0.0_dp
      r_e_sq_X(:) = 0.0_dp
      r_h_sq_X(:) = 0.0_dp
      r_e_sq_Y(:) = 0.0_dp
      r_h_sq_Y(:) = 0.0_dp
      r_e_h_XX(:) = 0.0_dp
      r_e_h_XY(:) = 0.0_dp
      r_e_h_YX(:) = 0.0_dp
      r_e_h_YY(:) = 0.0_dp

      norm_X = 0.0_dp
      norm_Y = 0.0_dp
      norm_XpY = 0.0_dp

      ! Initialize values of exciton descriptors
      exc_descr(i_exc)%r_e(:) = 0.0_dp
      exc_descr(i_exc)%r_h(:) = 0.0_dp
      exc_descr(i_exc)%r_e_sq(:) = 0.0_dp
      exc_descr(i_exc)%r_h_sq(:) = 0.0_dp
      exc_descr(i_exc)%r_e_h(:) = 0.0_dp

      exc_descr(i_exc)%flag_TDA = flag_TDA
      exc_descr(i_exc)%norm_XpY = 0.0_dp

      ! Norm of X
      CALL cp_fm_trace(fm_X_ia, fm_X_ia, norm_X)
      norm_XpY = norm_X
      ! Norm of Y
      IF (.NOT. flag_TDA) THEN
         CALL cp_fm_trace(fm_Y_ia, fm_Y_ia, norm_Y)
         norm_XpY = norm_XpY + norm_Y
      END IF

      exc_descr(i_exc)%norm_XpY = norm_XpY

      ! <r_h>_X = Tr[ X^T µ_ij X + Y µ_ab Y^T ] = X_ai µ_ij X_ja + Y_ia  µ_ab Y_bi
      DO i_dir = 1, 3
         ! <r_h>_X = X_ai µ_ij X_ja + ...
         CALL trace_exciton_descr(fm_X_ia, fm_multipole_ij_trunc(i_dir), fm_X_ia, r_h_X(i_dir))
         r_h_X(i_dir) = r_h_X(i_dir)/norm_XpY
         IF (.NOT. flag_TDA) THEN
            ! <r_h>_X = ... + Y_ia  µ_ab Y_bi
            CALL trace_exciton_descr(fm_Y_ia, fm_Y_ia, fm_multipole_ab_trunc(i_dir), r_h_Y(i_dir))
            r_h_Y(i_dir) = r_h_Y(i_dir)/norm_XpY
         END IF
      END DO
      exc_descr(i_exc)%r_h(:) = r_h_X(:) + r_h_Y(:)

      ! <r_e>_X = Tr[ X µ_ab X^T + Y^T µ_ij Y ] = X_ia µ_ab X_bi + Y_ai µ_ij Y_ja
      DO i_dir = 1, 3
         ! <r_e>_X = work_ib X_bi + ... = X_ib^T work_ib + ...
         CALL trace_exciton_descr(fm_X_ia, fm_X_ia, fm_multipole_ab_trunc(i_dir), r_e_X(i_dir))
         r_e_X(i_dir) = r_e_X(i_dir)/norm_XpY
         IF (.NOT. flag_TDA) THEN
            ! <r_e>_X = ... + Y_ai µ_ij Y_ja
            CALL trace_exciton_descr(fm_Y_ia, fm_multipole_ij_trunc(i_dir), fm_Y_ia, r_e_Y(i_dir))
            r_e_Y(i_dir) = r_e_Y(i_dir)/norm_XpY
         END IF
      END DO
      exc_descr(i_exc)%r_e(:) = r_e_X(:) + r_e_Y(:)

      ! <r_h^2>_X = Tr[ X^T M_ij X + Y M_ab Y^T ] = X_ai M_ij X_ja + Y_ia  M_ab Y_bi
      DO i_dir = 1, 3
         ! <r_h^2>_X = X_ai M_ij X_ja + ...
         CALL trace_exciton_descr(fm_X_ia, fm_multipole_ij_trunc(mask_quadrupole(i_dir)), &
                                  fm_X_ia, r_h_sq_X(i_dir))
         r_h_sq_X(i_dir) = r_h_sq_X(i_dir)/norm_XpY
         IF (.NOT. flag_TDA) THEN
            ! <r_h^2>_X = ... + Y_ia  M_ab Y_bi
            CALL trace_exciton_descr(fm_Y_ia, fm_Y_ia, &
                                     fm_multipole_ab_trunc(mask_quadrupole(i_dir)), r_h_sq_Y(i_dir))
            r_h_sq_Y(i_dir) = r_h_sq_Y(i_dir)/norm_XpY
         END IF
      END DO
      exc_descr(i_exc)%r_h_sq(:) = r_h_sq_X(:) + r_h_sq_Y(:)

      ! <r_e^2>_X = Tr[ X M_ab X^T + Y^T M_ij Y ] = X_ia M_ab X_bi + Y_ai M_ij Y_ja
      DO i_dir = 1, 3
         ! <r_e^2>_X = work_ib X_bi + ... = X_ib^T work_ib + ...
         CALL trace_exciton_descr(fm_X_ia, fm_X_ia, &
                                  fm_multipole_ab_trunc(mask_quadrupole(i_dir)), r_e_sq_X(i_dir))
         r_e_sq_X(i_dir) = r_e_sq_X(i_dir)/norm_XpY
         IF (.NOT. flag_TDA) THEN
            ! <r_e^2>_X = ... + Y_ai M_ij Y_ja
            CALL trace_exciton_descr(fm_Y_ia, fm_multipole_ij_trunc(mask_quadrupole(i_dir)), &
                                     fm_Y_ia, r_e_sq_Y(i_dir))
            r_e_sq_Y(i_dir) = r_e_sq_Y(i_dir)/norm_XpY
         END IF
      END DO
      exc_descr(i_exc)%r_e_sq(:) = r_e_sq_X(:) + r_e_sq_Y(:)

      ! <r_h r_e>_X
      !   = Tr[ X^T µ_ij X µ_ab  +  Y^T µ_ij Y µ_ab + X µ_ai Y µ_ai + Y µ_ai X µ_ai]
      !   = X_bj µ_ji X_ia µ_ab + Y_bj µ_ji Y_ia µ_ab + X_ia µ_aj Y_jb µ_bi + Y_ia µ_aj X_jb µ_bi
      CALL cp_fm_create(fm_work_ia, fm_struct_ia)
      CALL cp_fm_create(fm_work_ia_2, fm_struct_ia)
      CALL cp_fm_create(fm_work_ba, fm_struct_ab)
      DO i_dir = 1, 3
         ! First term - X^T µ_ij X µ_ab
         CALL cp_fm_set_all(fm_work_ia, 0.0_dp)
         CALL cp_fm_set_all(fm_work_ia_2, 0.0_dp)
         ! work_ib = X_ia µ_ab
         CALL parallel_gemm("N", "N", homo, virtual, virtual, 1.0_dp, &
                            fm_X_ia, fm_multipole_ab_trunc(i_dir), 0.0_dp, fm_work_ia)
         ! work_ja_2 = µ_ji work_ia
         CALL parallel_gemm("N", "N", homo, virtual, homo, 1.0_dp, &
                            fm_multipole_ij_trunc(i_dir), fm_work_ia, 0.0_dp, fm_work_ia_2)
         ! <r_h r_e>_X = work_ia_2 X_bj + ... = X^T work_ia_2 + ...
         CALL cp_fm_trace(fm_X_ia, fm_work_ia_2, r_e_h_XX(i_dir))
         r_e_h_XX(i_dir) = r_e_h_XX(i_dir)/norm_XpY
         IF (.NOT. flag_TDA) THEN
            ! Second term -  Y^T µ_ij Y µ_ab
            CALL cp_fm_set_all(fm_work_ia, 0.0_dp)
            CALL cp_fm_set_all(fm_work_ia_2, 0.0_dp)
            ! work_ib = Y_ia µ_ab
            CALL parallel_gemm("N", "N", homo, virtual, virtual, 1.0_dp, &
                               fm_Y_ia, fm_multipole_ab_trunc(i_dir), 0.0_dp, fm_work_ia)
            ! work_ja_2 = µ_ji work_ia
            CALL parallel_gemm("N", "N", homo, virtual, homo, 1.0_dp, &
                               fm_multipole_ij_trunc(i_dir), fm_work_ia, 0.0_dp, fm_work_ia_2)
            ! <r_h r_e>_X = work_ia_2 Y_bj + ... = Y^T work_ia_2 + ...
            CALL cp_fm_trace(fm_Y_ia, fm_work_ia_2, r_e_h_YY(i_dir))
            r_e_h_YY(i_dir) = r_e_h_YY(i_dir)/norm_XpY

            ! Third term - X µ_ai Y µ_ai = X_ia µ_aj Y_jb µ_bi
            !     Reshuffle for usage of trace (where first argument is transposed)
            !     = µ_aj Y_jb µ_bi X_ia =
            !        \__________/
            !         fm_work_ai
            !     fm_work_ai = µ_aj Y_jb µ_bi
            !     fm_work_ia = µ_ib Y_bj µ_ja
            !                        \_____/
            !                       fm_work_ba
            CALL cp_fm_set_all(fm_work_ba, 0.0_dp)
            CALL cp_fm_set_all(fm_work_ia, 0.0_dp)
            ! fm_work_ba = Y_bj µ_ja
            CALL parallel_gemm("T", "T", virtual, virtual, homo, 1.0_dp, &
                               fm_Y_ia, fm_multipole_ai_trunc(i_dir), 0.0_dp, fm_work_ba)
            ! fm_work_ia = µ_ib fm_work_ba
            CALL parallel_gemm("T", "N", homo, virtual, virtual, 1.0_dp, &
                               fm_multipole_ai_trunc(i_dir), fm_work_ba, 0.0_dp, fm_work_ia)
            ! <r_h r_e>_X = ... + X_ia µ_aj Y_jb µ_bi
            CALL cp_fm_trace(fm_work_ia, fm_X_ia, r_e_h_XY(i_dir))
            r_e_h_XY(i_dir) = r_e_h_XY(i_dir)/norm_XpY

            ! Fourth term - Y µ_ai X µ_ai = Y_ia µ_aj X_jb µ_bi
            !     Reshuffle for usage of trace (where first argument is transposed)
            !     = µ_aj X_jb µ_bi Y_ia =
            !        \__________/
            !         fm_work_ai
            !     fm_work_ai = µ_aj X_jb µ_bi
            !     fm_work_ia = µ_ib X_bj µ_ja
            !                        \_____/
            !                       fm_work_ba
            CALL cp_fm_set_all(fm_work_ba, 0.0_dp)
            CALL cp_fm_set_all(fm_work_ia, 0.0_dp)
            ! fm_work_ba = Y_bj µ_ja
            CALL parallel_gemm("T", "T", virtual, virtual, homo, 1.0_dp, &
                               fm_X_ia, fm_multipole_ai_trunc(i_dir), 0.0_dp, fm_work_ba)
            ! fm_work_ia = µ_ib fm_work_ba
            CALL parallel_gemm("T", "N", homo, virtual, virtual, 1.0_dp, &
                               fm_multipole_ai_trunc(i_dir), fm_work_ba, 0.0_dp, fm_work_ia)
            ! <r_h r_e>_X = ... + X_ia µ_aj Y_jb µ_bi
            CALL cp_fm_trace(fm_work_ia, fm_Y_ia, r_e_h_YX(i_dir))
            r_e_h_YX(i_dir) = r_e_h_YX(i_dir)/norm_XpY
         END IF
      END DO
      exc_descr(i_exc)%r_e_h(:) = r_e_h_XX(:) + r_e_h_XY(:) + r_e_h_YX(:) + r_e_h_YY(:)

      CALL cp_fm_release(fm_work_ia)
      CALL cp_fm_release(fm_work_ia_2)
      CALL cp_fm_release(fm_work_ba)

      ! diff_r_abs = |<r_h>_X - <r_e>_X|
      exc_descr(i_exc)%diff_r_abs = SQRT(SUM((exc_descr(i_exc)%r_h(:) - exc_descr(i_exc)%r_e(:))**2))

      ! σ_e = sqrt( <r_e^2>_X - <r_e>_X^2 )
      exc_descr(i_exc)%sigma_e = SQRT(SUM(exc_descr(i_exc)%r_e_sq(:)) - SUM(exc_descr(i_exc)%r_e(:)**2))

      ! σ_h = sqrt( <r_h^2>_X - <r_h>_X^2 )
      exc_descr(i_exc)%sigma_h = SQRT(SUM(exc_descr(i_exc)%r_h_sq(:)) - SUM(exc_descr(i_exc)%r_h(:)**2))

      ! COV(r_e, r_h) = < r_e r_h >_X - < r_e >_X < r_h >_X
      exc_descr(i_exc)%cov_e_h_sum = 0.0_dp
      exc_descr(i_exc)%cov_e_h(:) = 0.0_dp
      DO i_dir = 1, 3
         exc_descr(i_exc)%cov_e_h(i_dir) = exc_descr(i_exc)%r_e_h(i_dir) &
                                           - exc_descr(i_exc)%r_e(i_dir)*exc_descr(i_exc)%r_h(i_dir)
         exc_descr(i_exc)%cov_e_h_sum = exc_descr(i_exc)%cov_e_h_sum + &
                                        exc_descr(i_exc)%r_e_h(i_dir) - exc_descr(i_exc)%r_e(i_dir)*exc_descr(i_exc)%r_h(i_dir)
      END DO

      ! root-mean-square e-h separation
      exc_descr(i_exc)%diff_r_sqr = SQRT(exc_descr(i_exc)%diff_r_abs**2 + &
                                         exc_descr(i_exc)%sigma_e**2 + exc_descr(i_exc)%sigma_h**2 &
                                         - 2*exc_descr(i_exc)%cov_e_h_sum)

      ! e-h-correlation coefficient R_eh = COV(r_e, r_h) / ( σ_e σ_h )
      exc_descr(i_exc)%corr_e_h = exc_descr(i_exc)%cov_e_h_sum/(exc_descr(i_exc)%sigma_e*exc_descr(i_exc)%sigma_h)

      ! Expectation values of r_e and r_h
      exc_descr(i_exc)%r_e_shift(:) = exc_descr(i_exc)%r_e(:)
      exc_descr(i_exc)%r_h_shift(:) = exc_descr(i_exc)%r_h(:)

      CALL cp_fm_struct_release(fm_struct_ia)
      CALL cp_fm_struct_release(fm_struct_ab)

      CALL timestop(handle)

   END SUBROUTINE get_exciton_descriptors

END MODULE bse_properties
