!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Finite-volume Kubo-Greenwood transport from converged Quickstep matrices.
! **************************************************************************************************
MODULE qs_kubo_transport
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: a_bohr,&
                                              angstrom,&
                                              e_charge,&
                                              evolt,&
                                              h_bar,&
                                              kelvin
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kubo_transport'

   PUBLIC :: qs_scf_post_kubo_transport

CONTAINS

! **************************************************************************************************
!> \brief Compute and print the finite-volume Kubo-Greenwood conductivity tensor.
!> \param qs_env Quickstep environment
! **************************************************************************************************
   SUBROUTINE qs_scf_post_kubo_transport(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_kubo_transport'

      INTEGER                                            :: handle, imu, ispin, nao, natom, &
                                                            nelectron_total, neutral_grid, nmu, &
                                                            nspins, output_unit
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_atom
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes
      LOGICAL                                            :: do_kpoints, explicit, neutral_mu_explicit
      REAL(KIND=dp)                                      :: density_factor, dissipation, emax, emin, &
                                                            mu, mu0, mu_step, ne, nh, sigma_iso, &
                                                            temperature, volume, volume_m3
      REAL(KIND=dp), DIMENSION(:), POINTER               :: energy_range
      REAL(KIND=dp), DIMENSION(3, 3)                     :: sigma
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: maxocc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eig, s_dense
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: dHH
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: kubo_section

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, cell, energy_range, kubo_section, logger, matrix_ks, matrix_s, mos, para_env, &
               particle_set, row_blk_sizes)

      CALL get_qs_env(qs_env, input=kubo_section)
      kubo_section => section_vals_get_subs_vals(kubo_section, "PROPERTIES%KUBO_TRANSPORT")
      CALL section_vals_get(kubo_section, explicit=explicit)
      IF (.NOT. explicit) THEN
         CALL timestop(handle)
         RETURN
      END IF

      CALL get_qs_env(qs_env, blacs_env=blacs_env, cell=cell, do_kpoints=do_kpoints, &
                      matrix_ks=matrix_ks, matrix_s=matrix_s, mos=mos, natom=natom, &
                      nelectron_total=nelectron_total, para_env=para_env, particle_set=particle_set)

      IF (do_kpoints) CPABORT("KUBO_TRANSPORT currently expects a finite Gamma-point supercell.")
      CPASSERT(ASSOCIATED(matrix_s))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(particle_set))

      CALL section_vals_val_get(kubo_section, "TEMPERATURE", r_val=temperature)
      CALL section_vals_val_get(kubo_section, "DISSIPATION", r_val=dissipation)
      CALL section_vals_val_get(kubo_section, "ENERGY_RANGE", r_vals=energy_range)
      CALL section_vals_val_get(kubo_section, "NEUTRAL_MU", explicit=neutral_mu_explicit)
      IF (neutral_mu_explicit) CALL section_vals_val_get(kubo_section, "NEUTRAL_MU", r_val=mu0)
      CALL section_vals_val_get(kubo_section, "N_MU", i_val=nmu)
      CALL section_vals_val_get(kubo_section, "NEUTRAL_GRID", i_val=neutral_grid)

      IF (temperature <= 0.0_dp) CPABORT("KUBO_TRANSPORT TEMPERATURE must be positive.")
      IF (dissipation <= 0.0_dp) CPABORT("KUBO_TRANSPORT DISSIPATION must be positive.")
      IF (nmu < 1) CPABORT("KUBO_TRANSPORT N_MU must be at least one.")
      IF (neutral_grid < 10) CPABORT("KUBO_TRANSPORT NEUTRAL_GRID must be at least ten.")

      nspins = SIZE(matrix_ks)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao, row_blk_size=row_blk_sizes)
      CALL build_ao_atom_map(row_blk_sizes, natom, ao_atom)
      CALL dbcsr_to_dense(matrix_s(1)%matrix, blacs_env, para_env, s_dense)

      ALLOCATE (eig(nao, nspins), dHH(nao, nao, 3, nspins), maxocc(nspins))

      DO ispin = 1, nspins
         maxocc(ispin) = mos(ispin)%maxocc
         CALL setup_spin_transport(matrix_ks(ispin)%matrix, s_dense, blacs_env, para_env, &
                                   particle_set, cell, ao_atom, eig(:, ispin), dHH(:, :, :, ispin))
      END DO

      emin = MINVAL(eig)
      emax = MAXVAL(eig)
      IF (.NOT. neutral_mu_explicit) THEN
         CALL find_neutral_mu(eig, maxocc, temperature, REAL(nelectron_total, KIND=dp), neutral_grid, &
                              emin, emax, mu0)
      END IF

      IF (energy_range(1) < energy_range(2)) THEN
         emin = energy_range(1)
         emax = energy_range(2)
      END IF

      volume = cell%deth
      volume_m3 = volume*a_bohr**3
      density_factor = 1.0_dp/volume_m3

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(/,T2,A)') "Kubo-Greenwood finite-volume transport"
         WRITE (output_unit, '(T3,A,T38,I12)') "Number of atoms:", natom
         WRITE (output_unit, '(T3,A,T38,I12)') "Number of AO functions:", nao
         WRITE (output_unit, '(T3,A,T38,F12.4)') "Temperature [K]:", temperature*kelvin
         WRITE (output_unit, '(T3,A,T38,F12.4)') "Dissipation [K]:", dissipation*kelvin
         WRITE (output_unit, '(T3,A,T38,F12.6)') "Dissipation [eV]:", dissipation*evolt
         WRITE (output_unit, '(T3,A,T38,F12.6)') "Neutral chemical potential [eV]:", mu0*evolt
         WRITE (output_unit, '(/,T3,A)') "mu[Ha]        mu-mu0[eV]        Ne[m^-3]        Nh[m^-3]        "// &
            "sigma_iso[S/cm]   sigma_xx[S/m]     sigma_xy[S/m]     sigma_xz[S/m]     sigma_yx[S/m]     "// &
            "sigma_yy[S/m]     sigma_yz[S/m]     sigma_zx[S/m]     sigma_zy[S/m]     sigma_zz[S/m]"
      END IF

      IF (nmu == 1) THEN
         mu_step = 0.0_dp
      ELSE
         mu_step = (emax - emin)/REAL(nmu - 1, KIND=dp)
      END IF

      DO imu = 1, nmu
         mu = emin + REAL(imu - 1, KIND=dp)*mu_step
         CALL electron_hole_density(eig, maxocc, mu, mu0, temperature, ne, nh)
         CALL conductivity_at_mu(eig, dHH, maxocc, mu, temperature, dissipation, volume, sigma)
         IF (output_unit > 0) THEN
            sigma_iso = (sigma(1, 1) + sigma(2, 2) + sigma(3, 3))/3.0_dp*1.0E8_dp
            WRITE (output_unit, '(T3,F12.6,F14.6,2ES16.7,10ES17.8)') mu, (mu - mu0)*evolt, &
               ne*density_factor, nh*density_factor, &
               sigma_iso, &
               sigma(1, 1)*1.0E10_dp, sigma(1, 2)*1.0E10_dp, sigma(1, 3)*1.0E10_dp, &
               sigma(2, 1)*1.0E10_dp, sigma(2, 2)*1.0E10_dp, sigma(2, 3)*1.0E10_dp, &
               sigma(3, 1)*1.0E10_dp, sigma(3, 2)*1.0E10_dp, sigma(3, 3)*1.0E10_dp
         END IF
      END DO

      DEALLOCATE (ao_atom, dHH, eig, maxocc, s_dense)

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_kubo_transport

! **************************************************************************************************
!> \brief Build a global AO -> atom map from DBCSR atomic block sizes.
! **************************************************************************************************
   SUBROUTINE build_ao_atom_map(row_blk_sizes, natom, ao_atom)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: row_blk_sizes
      INTEGER, INTENT(IN)                                :: natom
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: ao_atom

      INTEGER                                            :: iatom, iao, nao, u

      nao = SUM(row_blk_sizes(1:natom))
      ALLOCATE (ao_atom(nao))
      iao = 0
      DO iatom = 1, natom
         DO u = 1, row_blk_sizes(iatom)
            iao = iao + 1
            ao_atom(iao) = iatom
         END DO
      END DO

   END SUBROUTINE build_ao_atom_map

! **************************************************************************************************
!> \brief Copy a DBCSR matrix to a replicated dense matrix.
! **************************************************************************************************
   SUBROUTINE dbcsr_to_dense(matrix, blacs_env, para_env, dense)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: dense

      INTEGER                                            :: ncol, nrow
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm

      NULLIFY (fm_struct)

      CALL dbcsr_get_info(matrix, nfullrows_total=nrow, nfullcols_total=ncol)
      ALLOCATE (dense(nrow, ncol))
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nrow, ncol_global=ncol)
      CALL cp_fm_create(fm, fm_struct)
      CALL copy_dbcsr_to_fm(matrix, fm)
      CALL cp_fm_get_submatrix(fm, dense, n_rows=nrow, n_cols=ncol)
      CALL cp_fm_release(fm)
      CALL cp_fm_struct_release(fm_struct)

   END SUBROUTINE dbcsr_to_dense

! **************************************************************************************************
!> \brief Prepare eigenvalues and the transformed [X,H] matrices for one spin channel.
! **************************************************************************************************
   SUBROUTINE setup_spin_transport(ks_matrix, s_dense, blacs_env, para_env, particle_set, cell, &
                                   ao_atom, eig, dHH)
      TYPE(dbcsr_type), INTENT(IN)                       :: ks_matrix
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: s_dense
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_atom
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eig
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: dHH

      INTEGER                                            :: idir, nao
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gmat, h_dense, horth, op, op_orth, uvec, work

      nao = SIZE(s_dense, 1)
      ALLOCATE (gmat(nao, nao), horth(nao, nao), op(nao, nao), op_orth(nao, nao), uvec(nao, nao), &
                work(nao, nao))

      CALL dbcsr_to_dense(ks_matrix, blacs_env, para_env, h_dense)
      CALL inverse_sqrt_symmetric(s_dense, gmat)

      work = MATMUL(h_dense, gmat)
      horth = MATMUL(gmat, work)
      CALL symmetrize(horth)

      uvec = horth
      CALL diagonalize_symmetric(uvec, eig)

      DO idir = 1, 3
         CALL build_commutator_kernel(h_dense, particle_set, cell, ao_atom, idir, op)
         work = MATMUL(op, gmat)
         op_orth = MATMUL(gmat, work)
         CALL transform_to_eigenbasis(op_orth, uvec, dHH(:, :, idir), work)
      END DO

      DEALLOCATE (gmat, h_dense, horth, op, op_orth, uvec, work)

   END SUBROUTINE setup_spin_transport

! **************************************************************************************************
!> \brief Symmetric inverse square root via dense diagonalization.
! **************************************************************************************************
   SUBROUTINE inverse_sqrt_symmetric(matrix, inverse_sqrt)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: matrix
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: inverse_sqrt

      INTEGER                                            :: i, info, lwork, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eig, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: scaled, vectors

      n = SIZE(matrix, 1)
      lwork = MAX(1, 8*n)
      ALLOCATE (eig(n), scaled(n, n), vectors(n, n), work(lwork))

      vectors = matrix
      CALL dsyev("V", "L", n, vectors, n, eig, work, lwork, info)
      IF (info /= 0) CPABORT("Diagonalization of the overlap matrix failed.")
      IF (MINVAL(eig) <= 0.0_dp) CPABORT("Overlap matrix is not positive definite.")

      scaled = vectors
      DO i = 1, n
         scaled(:, i) = scaled(:, i)/SQRT(eig(i))
      END DO
      CALL dgemm("N", "T", n, n, n, 1.0_dp, scaled, n, vectors, n, 0.0_dp, inverse_sqrt, n)

      DEALLOCATE (eig, scaled, vectors, work)

   END SUBROUTINE inverse_sqrt_symmetric

! **************************************************************************************************
!> \brief Dense symmetric diagonalization.
! **************************************************************************************************
   SUBROUTINE diagonalize_symmetric(matrix, eig)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eig

      INTEGER                                            :: info, lwork, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: work

      n = SIZE(matrix, 1)
      lwork = MAX(1, 8*n)
      ALLOCATE (work(lwork))
      CALL dsyev("V", "L", n, matrix, n, eig, work, lwork, info)
      IF (info /= 0) CPABORT("Diagonalization of the orthogonalized KS matrix failed.")
      DEALLOCATE (work)

   END SUBROUTINE diagonalize_symmetric

! **************************************************************************************************
!> \brief Enforce exact symmetry after dense products.
! **************************************************************************************************
   SUBROUTINE symmetrize(matrix)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: matrix

      INTEGER                                            :: i, j, n
      REAL(KIND=dp)                                      :: value

      n = SIZE(matrix, 1)
      DO j = 1, n
         DO i = j + 1, n
            value = 0.5_dp*(matrix(i, j) + matrix(j, i))
            matrix(i, j) = value
            matrix(j, i) = value
         END DO
      END DO

   END SUBROUTINE symmetrize

! **************************************************************************************************
!> \brief Build element-wise kernel (r_col-r_row)_idir * matrix(row,col).
! **************************************************************************************************
   SUBROUTINE build_commutator_kernel(matrix, particle_set, cell, ao_atom, idir, op)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: matrix
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_atom
      INTEGER, INTENT(IN)                                :: idir
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: op

      INTEGER                                            :: i, j, n
      REAL(KIND=dp), DIMENSION(3)                        :: dr

      n = SIZE(matrix, 1)
      DO j = 1, n
         DO i = 1, n
            dr(:) = pbc(particle_set(ao_atom(j))%r(:) - particle_set(ao_atom(i))%r(:), cell)
            op(i, j) = dr(idir)*matrix(i, j)
         END DO
      END DO

   END SUBROUTINE build_commutator_kernel

! **************************************************************************************************
!> \brief Transform an AO matrix to the eigenvector basis.
! **************************************************************************************************
   SUBROUTINE transform_to_eigenbasis(op, uvec, transformed, work)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: op, uvec
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: transformed
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: work

      INTEGER                                            :: n

      n = SIZE(op, 1)
      CALL dgemm("N", "N", n, n, n, 1.0_dp, op, n, uvec, n, 0.0_dp, work, n)
      CALL dgemm("T", "N", n, n, n, 1.0_dp, uvec, n, work, n, 0.0_dp, transformed, n)

   END SUBROUTINE transform_to_eigenbasis

! **************************************************************************************************
!> \brief Fermi occupation with overflow protection.
! **************************************************************************************************
   PURE FUNCTION fermi_occ(eig, mu, kT, maxocc) RESULT(occ)
      REAL(KIND=dp), INTENT(IN)                          :: eig, mu, kT, maxocc
      REAL(KIND=dp)                                      :: occ

      REAL(KIND=dp)                                      :: x

      x = (eig - mu)/kT
      IF (x > 80.0_dp) THEN
         occ = 0.0_dp
      ELSEIF (x < -80.0_dp) THEN
         occ = maxocc
      ELSE
         occ = maxocc/(1.0_dp + EXP(x))
      END IF

   END FUNCTION fermi_occ

! **************************************************************************************************
!> \brief Fermi vacancy below the neutral point.
! **************************************************************************************************
   PURE FUNCTION fermi_hole(eig, mu, kT, maxocc) RESULT(hole)
      REAL(KIND=dp), INTENT(IN)                          :: eig, mu, kT, maxocc
      REAL(KIND=dp)                                      :: hole

      hole = maxocc - fermi_occ(eig, mu, kT, maxocc)

   END FUNCTION fermi_hole

! **************************************************************************************************
!> \brief Total electron count at a chemical potential.
! **************************************************************************************************
   PURE FUNCTION electron_count(eig, maxocc, mu, kT) RESULT(count)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eig
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: maxocc
      REAL(KIND=dp), INTENT(IN)                          :: mu, kT
      REAL(KIND=dp)                                      :: count

      INTEGER                                            :: ispin, n

      count = 0.0_dp
      DO ispin = 1, SIZE(eig, 2)
         DO n = 1, SIZE(eig, 1)
            count = count + fermi_occ(eig(n, ispin), mu, kT, maxocc(ispin))
         END DO
      END DO

   END FUNCTION electron_count

! **************************************************************************************************
!> \brief Determine the neutral chemical potential using the legacy two-step definition.
! **************************************************************************************************
   SUBROUTINE find_neutral_mu(eig, maxocc, kT, neutral_electrons, neutral_grid, emin, emax, mu0)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eig
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: maxocc
      REAL(KIND=dp), INTENT(IN)                          :: kT, neutral_electrons, emin, emax
      INTEGER, INTENT(IN)                                :: neutral_grid
      REAL(KIND=dp), INTENT(OUT)                         :: mu0

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: best, count, mu, mup, previous_count, ne, nh

      mup = 0.5_dp*(emin + emax)
      previous_count = electron_count(eig, maxocc, emin, kT)
      DO i = 1, neutral_grid
         mu = emin + (emax - emin)*REAL(i, KIND=dp)/REAL(neutral_grid, KIND=dp)
         count = electron_count(eig, maxocc, mu, kT)
         IF ((count - neutral_electrons)*(previous_count - neutral_electrons) <= 0.0_dp) mup = mu
         previous_count = count
      END DO

      mu0 = mup
      best = HUGE(1.0_dp)
      DO i = 1, neutral_grid
         mu = emin + (emax - emin)*REAL(i, KIND=dp)/REAL(neutral_grid, KIND=dp)
         CALL electron_hole_density(eig, maxocc, mu, mup, kT, ne, nh)
         IF (ABS(nh - ne) <= best) THEN
            mu0 = mu
            best = ABS(nh - ne)
         END IF
      END DO

   END SUBROUTINE find_neutral_mu

! **************************************************************************************************
!> \brief Electron and hole counts relative to a neutral point.
! **************************************************************************************************
   SUBROUTINE electron_hole_density(eig, maxocc, mu, mu0, kT, ne, nh)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eig
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: maxocc
      REAL(KIND=dp), INTENT(IN)                          :: mu, mu0, kT
      REAL(KIND=dp), INTENT(OUT)                         :: ne, nh

      INTEGER                                            :: ispin, n

      ne = 0.0_dp
      nh = 0.0_dp
      DO ispin = 1, SIZE(eig, 2)
         DO n = 1, SIZE(eig, 1)
            IF (eig(n, ispin) <= mu0) nh = nh + fermi_hole(eig(n, ispin), mu, kT, maxocc(ispin))
            IF (eig(n, ispin) >= mu0) ne = ne + fermi_occ(eig(n, ispin), mu, kT, maxocc(ispin))
         END DO
      END DO

   END SUBROUTINE electron_hole_density

! **************************************************************************************************
!> \brief Conductivity tensor at one chemical potential.
! **************************************************************************************************
   SUBROUTINE conductivity_at_mu(eig, dHH, maxocc, mu, kT, dissipation, volume, sigma)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eig
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: dHH
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: maxocc
      REAL(KIND=dp), INTENT(IN)                          :: mu, kT, dissipation, volume
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: sigma

      INTEGER                                            :: idir, ispin, jdir, m, n, nao
      REAL(KIND=dp)                                      :: delta, eps, factor, g0
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: occ

      nao = SIZE(eig, 1)
      eps = 1.0E-12_dp
      g0 = e_charge**2/(pi*h_bar)
      sigma = 0.0_dp

      ALLOCATE (occ(nao))

      DO ispin = 1, SIZE(eig, 2)
         DO n = 1, nao
            occ(n) = fermi_occ(eig(n, ispin), mu, kT, maxocc(ispin))
         END DO

         DO idir = 1, 3
            DO jdir = 1, 3
               DO n = 1, nao
                  DO m = 1, nao
                     delta = eig(n, ispin) - eig(m, ispin)
                     factor = (occ(n) - occ(m))/(delta + eps)*dissipation/(dissipation**2 + delta**2)
                     sigma(jdir, idir) = sigma(jdir, idir) + &
                        dHH(m, n, jdir, ispin)*dHH(n, m, idir, ispin)*factor
                  END DO
               END DO
            END DO
         END DO
      END DO

      sigma = pi*g0*sigma/(volume*angstrom)

      DEALLOCATE (occ)

   END SUBROUTINE conductivity_at_mu

END MODULE qs_kubo_transport
