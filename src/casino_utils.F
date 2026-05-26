!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Writer for CASINO gwfn.data files.
!> \par History
!>      05.2026 created [Codex]
! **************************************************************************************************
MODULE casino_utils

   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp2k_info,                       ONLY: cp2k_version
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: kpoint_env_initialize,&
                                              kpoint_init_cell_index,&
                                              kpoint_initialize,&
                                              kpoint_initialize_mo_set,&
                                              kpoint_initialize_mos
   USE kpoint_types,                    ONLY: get_kpoint_env,&
                                              get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_env_p_type,&
                                              kpoint_release,&
                                              kpoint_type
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: nso
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_scf_diagonalization,          ONLY: do_general_diag_kp
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_wannier90,                    ONLY: prepare_wannier90_scf_mos
   USE scf_control_types,               ONLY: scf_control_type
   USE string_utilities,                ONLY: lowercase
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'casino_utils'
   INTEGER, PARAMETER, PRIVATE          :: max_casino_l = 4

   PUBLIC :: write_casino

CONTAINS

! **************************************************************************************************
!> \brief Write a CASINO gwfn.data file from the converged GPW/GAPW wavefunction.
!> \param qs_env the QS environment
!> \param casino_section the DFT%PRINT%CASINO input section
! **************************************************************************************************
   SUBROUTINE write_casino(qs_env, casino_section)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(section_vals_type), INTENT(IN), POINTER       :: casino_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'write_casino'

      CHARACTER(len=default_path_length)                 :: filename
      INTEGER :: ao_num, col_offset, handle, iao, iatom, ikind, ikp, ikp_loc, ikp_out, imo, ipgf, &
         iset, ishell, ishell_loc, ispin, iw, k, l, mo_num, nao_shell, natoms, nel_tot, &
         ngth_pseudo, nkp, nkp_mo, nkp_out, nmo, npseudo_atoms, nreal_k, nset, nsgf, nsgp_pseudo, &
         nspins, output_unit, periodicity, prim_num, shell_num, zatom
      INTEGER, ALLOCATABLE, DIMENSION(:) :: agauge, ao_to_atom, atomic_number, cp2k_to_casino_ao, &
         first_shell, kp_order, prim_per_shell, shell_ang_mom, shell_type
      INTEGER, DIMENSION(2)                              :: kp_range, nmo_spin
      INTEGER, DIMENSION(:), POINTER                     :: npgf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: l_shell_set
      LOGICAL                                            :: casino_kpoints_created, do_kpoints, &
                                                            ionode, periodic, use_real_wfn, &
                                                            write_pseudos
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: kp_real
      REAL(KIND=dp)                                      :: cval, e_nn, eps_kpoint_real, kdotg, &
                                                            pseudo_tol, sval, zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coefficients, exponents, mo_scale, &
                                                            valence_charge
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: coord, kvec, mo_energy, mos_sgf, &
                                                            mos_sgf_im, shell_position
      REAL(KIND=dp), DIMENSION(3)                        :: scoord
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp, zetas
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_dummy, fm_mo_coeff, fm_mo_coeff_im
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(kpoint_env_p_type), DIMENSION(:), POINTER     :: kp_env
      TYPE(kpoint_type), POINTER                         :: casino_kpoints, kpoints
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mo_set_type), DIMENSION(:, :), POINTER        :: mos_kp
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_inter_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kind_set
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CALL timeset(routineN, handle)

      NULLIFY (basis_set, blacs_env, casino_kpoints, cell, dft_control, fm_struct, gcc, &
               gth_potential, kind_set, kp_env, kpoints, l_shell_set, logger, mos, mos_kp, npgf, &
               nshell, para_env, para_env_inter_kp, particle_set, sgp_potential, xkp, zetas)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      CPASSERT(ASSOCIATED(qs_env))

      CALL section_vals_val_get(casino_section, "FILENAME", c_val=filename)
      IF (LEN_TRIM(filename) == 0) filename = "gwfn.data"
      CALL section_vals_val_get(casino_section, "EPS_KPOINT_REAL", r_val=eps_kpoint_real)
      CALL section_vals_val_get(casino_section, "WRITE_PSEUDOPOTENTIALS", l_val=write_pseudos)

      CALL get_qs_env(qs_env, para_env=para_env)
      ionode = para_env%is_source()

      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set, qs_kind_set=kind_set, &
                      natom=natoms, dft_control=dft_control, nelectron_total=nel_tot, &
                      do_kpoints=do_kpoints, kpoints=kpoints, blacs_env=blacs_env)
      casino_kpoints => kpoints
      casino_kpoints_created = .FALSE.
      CALL prepare_casino_kpoint_grid(qs_env, casino_section, do_kpoints, kpoints, &
                                      casino_kpoints, casino_kpoints_created)
      nspins = dft_control%nspins
      IF (nspins > 2) CPABORT("CASINO gwfn.data supports at most two spin channels.")

      periodicity = COUNT(cell%perd /= 0)
      periodic = periodicity > 0
      pseudo_tol = 1.0E-8_dp

      ALLOCATE (coord(3, natoms), atomic_number(natoms), valence_charge(natoms))
      npseudo_atoms = 0
      ngth_pseudo = 0
      nsgp_pseudo = 0
      DO iatom = 1, natoms
         coord(:, iatom) = particle_set(iatom)%r(1:3)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(kind_set(ikind), zatom=zatom, zeff=zeff, &
                          gth_potential=gth_potential, sgp_potential=sgp_potential)
         IF (ABS(zeff) < pseudo_tol) zeff = REAL(zatom, KIND=dp)
         atomic_number(iatom) = zatom
         IF (ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential)) THEN
            atomic_number(iatom) = zatom + 200
            npseudo_atoms = npseudo_atoms + 1
            IF (ASSOCIATED(gth_potential)) ngth_pseudo = ngth_pseudo + 1
            IF (ASSOCIATED(sgp_potential)) nsgp_pseudo = nsgp_pseudo + 1
         ELSE IF (ABS(zeff - REAL(zatom, KIND=dp)) > pseudo_tol) THEN
            atomic_number(iatom) = zatom + 200
            npseudo_atoms = npseudo_atoms + 1
         END IF
         valence_charge(iatom) = zeff
      END DO

      IF (ionode .AND. write_pseudos .AND. nsgp_pseudo > 0) THEN
         CALL write_casino_sgp_pseudopotentials(kind_set, particle_set, natoms, filename, output_unit)
      END IF

      IF (periodic) THEN
         CALL periodic_nuclear_repulsion_energy(cell, periodicity, coord, valence_charge, e_nn)
      ELSE
         CALL nuclear_repulsion_energy(particle_set, kind_set, e_nn)
      END IF
      e_nn = e_nn/REAL(natoms, KIND=dp)

      IF (do_kpoints) THEN
         CALL get_kpoint_info(casino_kpoints, nkp=nkp, xkp=xkp, use_real_wfn=use_real_wfn)
      ELSE
         nkp = 1
         use_real_wfn = .TRUE.
      END IF
      nkp_mo = MERGE(nkp, 1, do_kpoints)

      ALLOCATE (kp_order(nkp_mo), kp_real(nkp_mo), kvec(3, nkp_mo))
      CALL build_kpoint_order(cell, periodic, do_kpoints, nkp, xkp, eps_kpoint_real, &
                              kp_order, kp_real, nkp_out, nreal_k, kvec)
      IF (do_kpoints .AND. use_real_wfn .AND. ANY(.NOT. kp_real(1:nkp_out))) THEN
         CPABORT("CASINO complex k-points require CP2K complex k-point wavefunctions.")
      END IF

      CALL get_qs_kind_set(kind_set, nshell=shell_num, npgf_seg=prim_num, nsgf=nsgf)
      ao_num = nsgf

      ALLOCATE (shell_type(shell_num), prim_per_shell(shell_num), first_shell(natoms + 1), &
                shell_ang_mom(shell_num), shell_position(3, shell_num), &
                exponents(prim_num), coefficients(prim_num), ao_to_atom(ao_num), &
                cp2k_to_casino_ao(ao_num), mo_scale(ao_num))

      ishell = 0
      ipgf = 0
      iao = 0
      DO iatom = 1, natoms
         first_shell(iatom) = ishell + 1
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(kind_set(ikind), basis_set=basis_set, basis_type="ORB")
         CALL get_gto_basis_set(basis_set, nset=nset, nshell=nshell, npgf=npgf, &
                                zet=zetas, gcc=gcc, l=l_shell_set)
         DO iset = 1, nset
            DO ishell_loc = 1, nshell(iset)
               ishell = ishell + 1
               l = l_shell_set(ishell_loc, iset)
               IF (l > max_casino_l) THEN
                  CPABORT("CASINO writer currently supports harmonic Gaussian shells up to g.")
               END IF
               shell_ang_mom(ishell) = l
               shell_type(ishell) = casino_shell_type(l)
               prim_per_shell(ishell) = npgf(iset)
               shell_position(:, ishell) = particle_set(iatom)%r(1:3)
               CALL casino_shell_coefficients(l, npgf(iset), zetas(1:npgf(iset), iset), &
                                              gcc(1:npgf(iset), ishell_loc, iset), &
                                              exponents(ipgf + 1:ipgf + npgf(iset)), &
                                              coefficients(ipgf + 1:ipgf + npgf(iset)))
               nao_shell = nso(l)
               DO k = 1, nao_shell
                  cp2k_to_casino_ao(iao + k) = iao + casino_cp2k_index(l, k)
                  mo_scale(iao + k) = casino_mo_scale(l, k)
                  ao_to_atom(iao + k) = iatom
               END DO
               ipgf = ipgf + npgf(iset)
               iao = iao + nao_shell
            END DO
         END DO
      END DO
      first_shell(natoms + 1) = shell_num + 1
      CPASSERT(ishell == shell_num)
      CPASSERT(ipgf == prim_num)
      CPASSERT(iao == ao_num)

      ALLOCATE (mo_energy(ao_num, nkp_mo*nspins))
      mo_energy(:, :) = 0.0_dp
      nmo_spin(:) = 0

      IF (do_kpoints) THEN
         CALL get_kpoint_info(casino_kpoints, kp_env=kp_env, kp_range=kp_range, nkp=nkp)
         CALL get_kpoint_env(kp_env(1)%kpoint_env, mos=mos_kp)
         DO ispin = 1, nspins
            CALL get_mo_set(mos_kp(1, ispin), nmo=nmo)
            IF (nmo < ao_num) THEN
               CPABORT("CASINO gwfn.data requires a complete MO set. Increase ADDED_MOS.")
            END IF
            nmo_spin(ispin) = nmo
         END DO
         mo_num = nkp*SUM(nmo_spin)
         CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=nsgf, ncol_global=mo_num)
         CALL cp_fm_create(fm_mo_coeff, fm_struct)
         CALL cp_fm_set_all(fm_mo_coeff, 0.0_dp)
         IF (.NOT. use_real_wfn) THEN
            CALL cp_fm_create(fm_mo_coeff_im, fm_struct)
            CALL cp_fm_set_all(fm_mo_coeff_im, 0.0_dp)
         END IF
         CALL cp_fm_struct_release(fm_struct)

         DO ispin = 1, nspins
            DO ikp = 1, nkp
               nmo = nmo_spin(ispin)
               col_offset = (ikp - 1)*nmo + (ispin - 1)*nmo_spin(1)*nkp
               IF (ikp >= kp_range(1) .AND. ikp <= kp_range(2)) THEN
                  ikp_loc = ikp - kp_range(1) + 1
                  CALL get_kpoint_env(kp_env(ikp_loc)%kpoint_env, mos=mos_kp)
                  IF (mos_kp(1, ispin)%use_mo_coeff_b) THEN
                     CALL copy_dbcsr_to_fm(mos_kp(1, ispin)%mo_coeff_b, mos_kp(1, ispin)%mo_coeff)
                  END IF
                  CALL cp_fm_to_fm_submat_general(mos_kp(1, ispin)%mo_coeff, fm_mo_coeff, &
                                                  nsgf, nmo, 1, 1, 1, col_offset + 1, blacs_env)
                  mo_energy(1:ao_num, ikp + (ispin - 1)*nkp_mo) = &
                     mos_kp(1, ispin)%eigenvalues(1:ao_num)
                  IF (.NOT. use_real_wfn) THEN
                     IF (mos_kp(2, ispin)%use_mo_coeff_b) THEN
                        CALL copy_dbcsr_to_fm(mos_kp(2, ispin)%mo_coeff_b, mos_kp(2, ispin)%mo_coeff)
                     END IF
                     CALL cp_fm_to_fm_submat_general(mos_kp(2, ispin)%mo_coeff, fm_mo_coeff_im, &
                                                     nsgf, nmo, 1, 1, 1, col_offset + 1, blacs_env)
                  END IF
               ELSE
                  CALL cp_fm_to_fm_submat_general(fm_dummy, fm_mo_coeff, &
                                                  nsgf, nmo, 1, 1, 1, col_offset + 1, blacs_env)
                  IF (.NOT. use_real_wfn) THEN
                     CALL cp_fm_to_fm_submat_general(fm_dummy, fm_mo_coeff_im, &
                                                     nsgf, nmo, 1, 1, 1, col_offset + 1, blacs_env)
                  END IF
               END IF
            END DO
         END DO
         CALL get_kpoint_info(casino_kpoints, para_env_inter_kp=para_env_inter_kp)
         CALL para_env_inter_kp%sum(mo_energy)
      ELSE
         CALL get_qs_env(qs_env, mos=mos)
         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin), nmo=nmo)
            IF (nmo < ao_num) THEN
               CPABORT("CASINO gwfn.data requires a complete MO set. Increase ADDED_MOS.")
            END IF
            nmo_spin(ispin) = nmo
            mo_energy(1:ao_num, 1 + (ispin - 1)*nkp_mo) = mos(ispin)%eigenvalues(1:ao_num)
         END DO
      END IF

      IF (do_kpoints .AND. .NOT. use_real_wfn) THEN
         ALLOCATE (agauge(3*natoms))
         DO iatom = 1, natoms
            scoord(:) = MATMUL(cell%h_inv, particle_set(iatom)%r(1:3))
            agauge(3*(iatom - 1) + 1:3*iatom) = -NINT(scoord(:))
         END DO
      END IF

      IF (ionode) THEN
         IF (npseudo_atoms > 0) THEN
            WRITE (output_unit, "((T2,A,I0,A))") "CASINO| Marked ", npseudo_atoms, &
               " pseudopotential atoms in gwfn.data."
            IF (ngth_pseudo > 0) THEN
               WRITE (output_unit, "((T2,A))") &
                  "CASINO| GTH pseudopotentials require matching external CASINO *_pp.data files."
            END IF
            IF (.NOT. write_pseudos) THEN
               WRITE (output_unit, "((T2,A))") &
                  "CASINO| WRITE_PSEUDOPOTENTIALS is disabled; provide CASINO *_pp.data files manually."
            END IF
         END IF
         WRITE (output_unit, "((T2,A,A))") 'CASINO| Writing gwfn.data file ', TRIM(filename)
         CALL open_file(file_name=filename, file_status="REPLACE", file_action="WRITE", &
                        file_form="FORMATTED", unit_number=iw)
         CALL write_casino_header(iw, periodicity, nspins, e_nn, nel_tot, natoms, coord, &
                                  atomic_number, valence_charge, cell, periodic, nkp_out, nreal_k, &
                                  kvec, shell_num, ao_num, prim_num, shell_ang_mom, shell_type, &
                                  prim_per_shell, first_shell, exponents, coefficients, shell_position)
      END IF

      ALLOCATE (mos_sgf(nsgf, ao_num), mos_sgf_im(nsgf, ao_num))
      mos_sgf(:, :) = 0.0_dp
      mos_sgf_im(:, :) = 0.0_dp

      IF (do_kpoints) THEN
         DO ispin = 1, nspins
            DO ikp_out = 1, nkp_out
               ikp = kp_order(ikp_out)
               col_offset = (ikp - 1)*nmo_spin(ispin) + (ispin - 1)*nmo_spin(1)*nkp
               CALL cp_fm_get_submatrix(fm_mo_coeff, mos_sgf, 1, col_offset + 1, nsgf, ao_num)
               IF (.NOT. use_real_wfn) THEN
                  CALL cp_fm_get_submatrix(fm_mo_coeff_im, mos_sgf_im, 1, col_offset + 1, nsgf, ao_num)
                  DO iao = 1, ao_num
                     iatom = ao_to_atom(iao)
                     kdotg = 2.0_dp*pi*DOT_PRODUCT(xkp(:, ikp), &
                                                   REAL(agauge(3*(iatom - 1) + 1:3*iatom), KIND=dp))
                     cval = COS(kdotg)
                     sval = SIN(kdotg)
                     DO imo = 1, ao_num
                        CALL rotate_complex_pair(mos_sgf(cp2k_to_casino_ao(iao), imo), &
                                                 mos_sgf_im(cp2k_to_casino_ao(iao), imo), cval, sval)
                     END DO
                  END DO
               ELSE
                  mos_sgf_im(:, :) = 0.0_dp
               END IF
               IF (ionode) THEN
                  CALL write_casino_orbitals(iw, mos_sgf, mos_sgf_im, cp2k_to_casino_ao, mo_scale, &
                                             ao_num,.NOT. kp_real(ikp_out))
               END IF
            END DO
         END DO
      ELSE
         DO ispin = 1, nspins
            IF (mos(ispin)%use_mo_coeff_b) THEN
               CALL copy_dbcsr_to_fm(mos(ispin)%mo_coeff_b, mos(ispin)%mo_coeff)
            END IF
            CALL cp_fm_get_submatrix(mos(ispin)%mo_coeff, mos_sgf, 1, 1, nsgf, ao_num)
            mos_sgf_im(:, :) = 0.0_dp
            IF (ionode) THEN
               CALL write_casino_orbitals(iw, mos_sgf, mos_sgf_im, cp2k_to_casino_ao, mo_scale, &
                                          ao_num, .FALSE.)
            END IF
         END DO
      END IF

      IF (ionode) THEN
         WRITE (iw, *)
         IF (periodic) THEN
            WRITE (iw, '(A)') "EIGENVALUES"
            WRITE (iw, '(A)') "-----------"
            DO ikp_out = 1, nkp_out
               ikp = kp_order(ikp_out)
               DO ispin = 1, nspins
                  IF (nspins == 1) THEN
                     WRITE (iw, '(A,I6,3F14.8)') "k", ikp_out, kvec(:, ikp_out)
                  ELSE
                     WRITE (iw, '(A,I3,A,I6,3F14.8)') "spin", ispin, " k", ikp_out, kvec(:, ikp_out)
                  END IF
                  CALL write_real_vector(iw, mo_energy(1:ao_num, ikp + (ispin - 1)*nkp_mo))
               END DO
            END DO
         END IF
         CALL close_file(unit_number=iw)
      END IF

      DEALLOCATE (mos_sgf, mos_sgf_im)
      IF (do_kpoints) THEN
         CALL cp_fm_release(fm_mo_coeff)
         IF (.NOT. use_real_wfn) CALL cp_fm_release(fm_mo_coeff_im)
      END IF
      IF (casino_kpoints_created) CALL kpoint_release(casino_kpoints)
      IF (ALLOCATED(agauge)) DEALLOCATE (agauge)
      DEALLOCATE (ao_to_atom, atomic_number, coefficients, coord, cp2k_to_casino_ao, exponents, &
                  first_shell, kp_order, kp_real, kvec, mo_energy, mo_scale, prim_per_shell, &
                  shell_ang_mom, shell_position, shell_type, valence_charge)

      CALL timestop(handle)
   END SUBROUTINE write_casino

! **************************************************************************************************
!> \brief Write CASINO pseudopotential files for the semilocal ECP kinds.
!> \param kind_set the QS kinds
!> \param particle_set the particle set
!> \param natoms the number of atoms
!> \param gwfn_filename the gwfn.data filename
!> \param output_unit output unit for log messages
! **************************************************************************************************
   SUBROUTINE write_casino_sgp_pseudopotentials(kind_set, particle_set, natoms, gwfn_filename, output_unit)
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: kind_set
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      INTEGER, INTENT(IN)                                :: natoms
      CHARACTER(LEN=*), INTENT(IN)                       :: gwfn_filename
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: iatom, ikind, zatom
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: written
      REAL(KIND=dp)                                      :: zeff
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      NULLIFY (sgp_potential)
      ALLOCATE (written(SIZE(kind_set)))
      written(:) = .FALSE.

      DO iatom = 1, natoms
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, element_symbol=element_symbol, &
                              kind_number=ikind, z=zatom)
         IF (written(ikind)) CYCLE

         CALL get_qs_kind(kind_set(ikind), sgp_potential=sgp_potential, zeff=zeff)
         IF (ASSOCIATED(sgp_potential)) THEN
            IF (ABS(zeff) < 1.0E-8_dp) zeff = REAL(zatom, KIND=dp)
            CALL write_casino_sgp_pseudopotential(sgp_potential, element_symbol, zatom, zeff, &
                                                  gwfn_filename, output_unit)
            written(ikind) = .TRUE.
         END IF
      END DO

      DEALLOCATE (written)
   END SUBROUTINE write_casino_sgp_pseudopotentials

! **************************************************************************************************
!> \brief Write a CASINO tabulated pseudopotential for a CP2K semilocal ECP.
!> \param sgp_potential the CP2K semilocal Gaussian potential
!> \param element_symbol the chemical symbol
!> \param zatom the nuclear charge
!> \param zeff the ECP valence charge
!> \param gwfn_filename the gwfn.data filename
!> \param output_unit output unit for log messages
! **************************************************************************************************
   SUBROUTINE write_casino_sgp_pseudopotential(sgp_potential, element_symbol, zatom, zeff, &
                                               gwfn_filename, output_unit)
      TYPE(sgp_potential_type), INTENT(IN), POINTER      :: sgp_potential
      CHARACTER(LEN=*), INTENT(IN)                       :: element_symbol
      INTEGER, INTENT(IN)                                :: zatom
      REAL(KIND=dp), INTENT(IN)                          :: zeff
      CHARACTER(LEN=*), INTENT(IN)                       :: gwfn_filename
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=default_path_length)                 :: pp_filename
      INTEGER                                            :: igrid, iw, l, local_l, ngrid, nloc, &
                                                            nsemiloc, sl_lmax
      INTEGER, DIMENSION(0:10)                           :: npot
      LOGICAL                                            :: ecp_local, ecp_semi_local, has_nlcc
      REAL(KIND=dp)                                      :: agrid, bgrid, r, rmax, rv_local
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: rgrid
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rpot

      CALL get_potential(potential=sgp_potential, ecp_local=ecp_local, &
                         ecp_semi_local=ecp_semi_local, nloc=nloc, sl_lmax=sl_lmax, &
                         npot=npot, has_nlcc=has_nlcc)
      IF (.NOT. ecp_local .OR. nloc == 0) THEN
         WRITE (output_unit, "((T2,A,A,A))") "CASINO| Cannot write ", TRIM(element_symbol), &
            "_pp.data: only CP2K semilocal ECP potentials are supported."
         RETURN
      END IF

      local_l = MERGE(sl_lmax + 1, 0, ecp_semi_local)
      rmax = 100.0_dp
      agrid = 70.0_dp*EXP(-5.0_dp*LOG(10.0_dp))/REAL(zatom, KIND=dp)
      bgrid = 1.0_dp/70.0_dp

      ngrid = 0
      DO
         r = agrid*(EXP(bgrid*REAL(ngrid, KIND=dp)) - 1.0_dp)
         IF (r > rmax) EXIT
         ngrid = ngrid + 1
      END DO

      ALLOCATE (rgrid(ngrid), rpot(0:local_l, ngrid))
      DO igrid = 1, ngrid
         r = agrid*(EXP(bgrid*REAL(igrid - 1, KIND=dp)) - 1.0_dp)
         rgrid(igrid) = r
         rv_local = casino_sgp_r_times_v(nloc, sgp_potential%nrloc(1:nloc), &
                                         sgp_potential%bloc(1:nloc), &
                                         sgp_potential%aloc(1:nloc), r, zeff, .TRUE.)
         DO l = 0, local_l
            rpot(l, igrid) = rv_local
            IF (l < local_l .AND. ecp_semi_local) THEN
               nsemiloc = npot(l)
               IF (nsemiloc > 0) THEN
                  rpot(l, igrid) = rpot(l, igrid) + &
                                   casino_sgp_r_times_v(nsemiloc, sgp_potential%nrpot(1:nsemiloc, l), &
                                                        sgp_potential%bpot(1:nsemiloc, l), &
                                                        sgp_potential%apot(1:nsemiloc, l), r, zeff, .FALSE.)
               END IF
            END IF
         END DO
      END DO
      rpot(:, :) = 2.0_dp*rpot(:, :)

      CALL casino_pp_filename(gwfn_filename, element_symbol, pp_filename)
      CALL open_file(file_name=pp_filename, file_status="REPLACE", file_action="WRITE", &
                     file_form="FORMATTED", unit_number=iw)
      WRITE (iw, '(A)') "CP2K ECP pseudopotential in real space"
      WRITE (iw, '(A)') "Atomic number and pseudo-charge"
      WRITE (iw, '(I6,1X,F18.10)') zatom, zeff
      WRITE (iw, '(A)') "Energy units (rydberg/hartree/ev):"
      WRITE (iw, '(A)') "rydberg"
      WRITE (iw, '(A)') "Angular momentum of local component (0=s,1=p,2=d..)"
      WRITE (iw, '(I6)') local_l
      WRITE (iw, '(A)') "NLRULE override (1) VMC/DMC (2) config gen (0 ==> input/default value)"
      WRITE (iw, '(2I6)') 0, 0
      WRITE (iw, '(A)') "Number of grid points"
      WRITE (iw, '(I8)') ngrid
      WRITE (iw, '(A)') "R(i) in atomic units"
      DO igrid = 1, ngrid
         WRITE (iw, '(ES20.12)') rgrid(igrid)
      END DO
      DO l = 0, local_l
         WRITE (iw, '(A,I0,A)') "r*potential (L=", l, ") in Ry"
         DO igrid = 1, ngrid
            WRITE (iw, '(ES20.12)') rpot(l, igrid)
         END DO
      END DO
      CALL close_file(unit_number=iw)

      WRITE (output_unit, "((T2,A,A))") "CASINO| Wrote pseudopotential file ", TRIM(pp_filename)
      IF (has_nlcc) THEN
         WRITE (output_unit, "((T2,A,A,A))") "CASINO| NLCC terms for ", TRIM(element_symbol), &
            " are not represented in CASINO *_pp.data."
      END IF

      DEALLOCATE (rgrid, rpot)
   END SUBROUTINE write_casino_sgp_pseudopotential

! **************************************************************************************************
!> \brief Return r times a CP2K semilocal Gaussian ECP channel in Hartree.
!> \param nterm the number of Gaussian terms
!> \param nr the CP2K r**(n-2) exponents
!> \param gaussian_exponent the Gaussian exponents
!> \param coefficient the Gaussian coefficients
!> \param r the radial grid point
!> \param zeff the ECP valence charge
!> \param local_channel true for the local Coulomb-tailed channel
!> \return r times the potential value
! **************************************************************************************************
   FUNCTION casino_sgp_r_times_v(nterm, nr, gaussian_exponent, coefficient, r, zeff, local_channel) RESULT(r_times_v)
      INTEGER, INTENT(IN)                                :: nterm
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nr
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: gaussian_exponent, coefficient
      REAL(KIND=dp), INTENT(IN)                          :: r, zeff
      LOGICAL, INTENT(IN)                                :: local_channel
      REAL(KIND=dp)                                      :: r_times_v

      INTEGER                                            :: iterm

      IF (r == 0.0_dp) THEN
         r_times_v = 0.0_dp
         RETURN
      END IF

      r_times_v = 0.0_dp
      IF (local_channel) r_times_v = -zeff
      DO iterm = 1, nterm
         CPASSERT(nr(iterm) >= 1)
         r_times_v = r_times_v + coefficient(iterm)*r**(nr(iterm) - 1)* &
                     EXP(-gaussian_exponent(iterm)*r*r)
      END DO
   END FUNCTION casino_sgp_r_times_v

! **************************************************************************************************
!> \brief Build the CASINO pseudopotential filename next to gwfn.data.
!> \param gwfn_filename the gwfn.data filename
!> \param element_symbol the chemical symbol
!> \param pp_filename the CASINO pseudopotential filename
! **************************************************************************************************
   SUBROUTINE casino_pp_filename(gwfn_filename, element_symbol, pp_filename)
      CHARACTER(LEN=*), INTENT(IN)                       :: gwfn_filename, element_symbol
      CHARACTER(LEN=*), INTENT(OUT)                      :: pp_filename

      CHARACTER(LEN=2)                                   :: symbol
      INTEGER                                            :: slash

      symbol = ADJUSTL(element_symbol)
      CALL lowercase(symbol)
      slash = INDEX(TRIM(gwfn_filename), "/", BACK=.TRUE.)
      IF (slash > 0) THEN
         pp_filename = gwfn_filename(1:slash)//TRIM(symbol)//"_pp.data"
      ELSE
         pp_filename = TRIM(symbol)//"_pp.data"
      END IF
   END SUBROUTINE casino_pp_filename

! **************************************************************************************************
!> \brief Prepare the k-point object used for CASINO export.
!> \param qs_env the QS environment
!> \param casino_section the CASINO print section
!> \param do_kpoints true when the SCF used k-points
!> \param kpoints_scf the converged SCF k-point object
!> \param kpoints_out the k-point object to write
!> \param created true if kpoints_out must be released by the caller
! **************************************************************************************************
   SUBROUTINE prepare_casino_kpoint_grid(qs_env, casino_section, do_kpoints, kpoints_scf, &
                                         kpoints_out, created)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(section_vals_type), INTENT(IN), POINTER       :: casino_section
      LOGICAL, INTENT(IN)                                :: do_kpoints
      TYPE(kpoint_type), POINTER                         :: kpoints_scf, kpoints_out
      LOGICAL, INTENT(OUT)                               :: created

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_casino_kpoint_grid'

      CHARACTER(LEN=default_string_length)               :: kp_scheme, reuse_reason
      INTEGER                                            :: aligned_blocks, aligned_max_size, &
                                                            handle, nfull, output_unit
      INTEGER, DIMENSION(3)                              :: nkp_grid
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: diis_step, full_grid, full_kpoint_grid, &
                                                            gamma_centered, reuse_scf_mos, &
                                                            reused_scf_mos, symmetry
      REAL(KIND=dp)                                      :: aligned_min_svalue, eps_geo, wsum
      REAL(KIND=dp), DIMENSION(3)                        :: kp_shift
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp_source
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp_source
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      created = .FALSE.
      kpoints_out => kpoints_scf
      NULLIFY (blacs_env, cell, cell_to_index, dft_control, logger, matrix_ks, matrix_s, mos, &
               para_env, particle_set, sab_nl, scf_control, scf_env, wkp_source, xkp_source)

      IF (.NOT. do_kpoints) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CPASSERT(ASSOCIATED(kpoints_scf))

      CALL get_kpoint_info(kpoints_scf, kp_scheme=kp_scheme, symmetry=symmetry, &
                           full_grid=full_grid, nkp_grid=nkp_grid, kp_shift=kp_shift, &
                           gamma_centered=gamma_centered, eps_geo=eps_geo)
      IF (.NOT. symmetry .OR. full_grid) THEN
         CALL timestop(handle)
         RETURN
      END IF

      CALL section_vals_val_get(casino_section, "FULL_KPOINT_GRID", l_val=full_kpoint_grid)
      IF (.NOT. full_kpoint_grid) THEN
         CPABORT("CASINO export requires a full k-point grid. Use PRINT%CASINO%FULL_KPOINT_GRID.")
      END IF

      SELECT CASE (TRIM(kp_scheme))
      CASE ("MONKHORST-PACK", "MACDONALD", "GENERAL")
         ! supported below
      CASE DEFAULT
         CPABORT("CASINO%FULL_KPOINT_GRID supports only MONKHORST-PACK, MACDONALD, and GENERAL k-points.")
      END SELECT

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)
      CALL section_vals_val_get(casino_section, "REUSE_SCF_MOS", l_val=reuse_scf_mos)
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env, cell=cell, &
                      particle_set=particle_set, mos=mos, dft_control=dft_control, &
                      sab_orb=sab_nl, matrix_ks_kp=matrix_ks, matrix_s_kp=matrix_s, &
                      scf_env=scf_env, scf_control=scf_control)
      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(ASSOCIATED(blacs_env))
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(mos))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(sab_nl))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_s))
      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(ASSOCIATED(scf_control))

      NULLIFY (kpoints_out)
      CALL kpoint_create(kpoints_out)
      kpoints_out%kp_scheme = kp_scheme
      kpoints_out%symmetry = .FALSE.
      kpoints_out%full_grid = .TRUE.
      kpoints_out%verbose = .FALSE.
      kpoints_out%use_real_wfn = .FALSE.
      kpoints_out%eps_geo = eps_geo
      kpoints_out%parallel_group_size = para_env%num_pe

      SELECT CASE (TRIM(kp_scheme))
      CASE ("MONKHORST-PACK", "MACDONALD")
         kpoints_out%nkp_grid(1:3) = nkp_grid(1:3)
         kpoints_out%kp_shift(1:3) = kp_shift(1:3)
         kpoints_out%gamma_centered = gamma_centered
         CALL kpoint_initialize(kpoints_out, particle_set, cell)
      CASE ("GENERAL")
         IF (.NOT. ASSOCIATED(kpoints_scf%xkp_input) .OR. &
             .NOT. ASSOCIATED(kpoints_scf%wkp_input)) THEN
            CPABORT("CASINO%FULL_KPOINT_GRID cannot recover the unreduced GENERAL k-point set.")
         END IF
         xkp_source => kpoints_scf%xkp_input
         wkp_source => kpoints_scf%wkp_input
         nfull = SIZE(wkp_source)
         wsum = SUM(wkp_source)
         IF (wsum <= 0.0_dp) CPABORT("CASINO%FULL_KPOINT_GRID found invalid GENERAL k-point weights.")
         kpoints_out%nkp = nfull
         ALLOCATE (kpoints_out%xkp(3, nfull), kpoints_out%wkp(nfull))
         kpoints_out%xkp(1:3, 1:nfull) = xkp_source(1:3, 1:nfull)
         kpoints_out%wkp(1:nfull) = wkp_source(1:nfull)/wsum
      END SELECT

      CALL kpoint_env_initialize(kpoints_out, para_env, blacs_env)
      CALL kpoint_initialize_mos(kpoints_out, mos)
      CALL kpoint_initialize_mo_set(kpoints_out)
      CALL kpoint_init_cell_index(kpoints_out, sab_nl, para_env, dft_control%nimages)

      reused_scf_mos = .FALSE.
      reuse_reason = ""
      aligned_blocks = 0
      aligned_max_size = 0
      aligned_min_svalue = 0.0_dp
      diis_step = .FALSE.
      IF (reuse_scf_mos) THEN
         CALL do_general_diag_kp(matrix_ks, matrix_s, kpoints_scf, scf_env, scf_control, .FALSE., &
                                 diis_step)
         CALL get_kpoint_info(kpoints_out, cell_to_index=cell_to_index)
         CALL prepare_wannier90_scf_mos(kpoints_out, kpoints_scf, matrix_s, matrix_ks, &
                                        cell_to_index, sab_nl, para_env, reused_scf_mos, &
                                        reuse_reason, aligned_blocks, aligned_max_size, &
                                        aligned_min_svalue)
      END IF
      IF (reused_scf_mos) THEN
         IF (output_unit > 0) THEN
            WRITE (output_unit, '(T2,A)') &
               "CASINO| Reused SCF MO coefficients for the full k-point grid."
            IF (aligned_blocks > 0) THEN
               WRITE (output_unit, '(T2,A,I0,A,I0,A,ES10.3)') &
                  "CASINO| Ritz-stabilized ", aligned_blocks, &
                  " degenerate SCF MO subspace(s); largest block has ", aligned_max_size, &
                  " band(s), min metric eigenvalue ", aligned_min_svalue
            END IF
         END IF
      ELSE
         IF (output_unit > 0) THEN
            IF (reuse_scf_mos) THEN
               WRITE (output_unit, '(T2,A,A)') &
                  "CASINO| Could not reuse SCF MOs: ", TRIM(reuse_reason)
            END IF
            WRITE (output_unit, '(T2,A)') &
               "CASINO| Diagonalizing the full k-point grid for export."
         END IF
         diis_step = .FALSE.
         CALL do_general_diag_kp(matrix_ks, matrix_s, kpoints_out, scf_env, scf_control, .FALSE., &
                                 diis_step)
      END IF
      created = .TRUE.

      CALL timestop(handle)
   END SUBROUTINE prepare_casino_kpoint_grid

! **************************************************************************************************
!> \brief Build the CASINO k-point order with all real k-points first.
!> \param cell ...
!> \param periodic ...
!> \param do_kpoints ...
!> \param nkp_total ...
!> \param xkp ...
!> \param eps_kpoint_real ...
!> \param kp_order ...
!> \param kp_real ...
!> \param nkp_out ...
!> \param nreal_k ...
!> \param kvec ...
! **************************************************************************************************
   SUBROUTINE build_kpoint_order(cell, periodic, do_kpoints, nkp_total, xkp, eps_kpoint_real, &
                                 kp_order, kp_real, nkp_out, nreal_k, kvec)
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      LOGICAL, INTENT(IN)                                :: periodic, do_kpoints
      INTEGER, INTENT(IN)                                :: nkp_total
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL, POINTER                               :: xkp
      REAL(KIND=dp), INTENT(IN)                          :: eps_kpoint_real
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: kp_order
      LOGICAL, DIMENSION(:), INTENT(OUT)                 :: kp_real
      INTEGER, INTENT(OUT)                               :: nkp_out, nreal_k
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: kvec

      INTEGER                                            :: ikp, jkp, ncomplex_k
      INTEGER, DIMENSION(nkp_total)                      :: complex_order, real_order
      LOGICAL, DIMENSION(nkp_total)                      :: used

      IF (.NOT. periodic) THEN
         kp_order(1) = 1
         kp_real(1) = .TRUE.
         nkp_out = 1
         nreal_k = 1
         kvec(:, 1) = 0.0_dp
         RETURN
      END IF

      IF (.NOT. do_kpoints) THEN
         kp_order(1) = 1
         kp_real(1) = .TRUE.
         nkp_out = 1
         nreal_k = 1
         kvec(:, 1) = 0.0_dp
         RETURN
      END IF

      used(:) = .FALSE.
      nreal_k = 0
      ncomplex_k = 0

      DO ikp = 1, nkp_total
         IF (is_real_kpoint(cell, xkp(:, ikp), eps_kpoint_real)) THEN
            used(ikp) = .TRUE.
            nreal_k = nreal_k + 1
            real_order(nreal_k) = ikp
         END IF
      END DO

      DO ikp = 1, nkp_total
         IF (.NOT. used(ikp)) THEN
            used(ikp) = .TRUE.
            ncomplex_k = ncomplex_k + 1
            complex_order(ncomplex_k) = ikp
            DO jkp = ikp + 1, nkp_total
               IF (.NOT. used(jkp)) THEN
                  IF (is_conjugate_kpoint(cell, xkp(:, ikp), xkp(:, jkp), eps_kpoint_real)) THEN
                     used(jkp) = .TRUE.
                     EXIT
                  END IF
               END IF
            END DO
         END IF
      END DO

      nkp_out = nreal_k + ncomplex_k
      DO ikp = 1, nreal_k
         kp_order(ikp) = real_order(ikp)
         kp_real(ikp) = .TRUE.
         kvec(:, ikp) = 2.0_dp*pi*MATMUL(TRANSPOSE(cell%h_inv), xkp(:, real_order(ikp)))
      END DO
      DO ikp = 1, ncomplex_k
         jkp = nreal_k + ikp
         kp_order(jkp) = complex_order(ikp)
         kp_real(jkp) = .FALSE.
         kvec(:, jkp) = 2.0_dp*pi*MATMUL(TRANSPOSE(cell%h_inv), xkp(:, complex_order(ikp)))
      END DO
   END SUBROUTINE build_kpoint_order

! **************************************************************************************************
!> \brief Returns true for conjugate k-points modulo reciprocal lattice vectors.
!> \param cell ...
!> \param xk1 ...
!> \param xk2 ...
!> \param eps_kpoint_real ...
!> \return ...
! **************************************************************************************************
   FUNCTION is_conjugate_kpoint(cell, xk1, xk2, eps_kpoint_real) RESULT(is_conjugate)
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xk1, xk2
      REAL(KIND=dp), INTENT(IN)                          :: eps_kpoint_real
      LOGICAL                                            :: is_conjugate

      INTEGER                                            :: idir
      REAL(KIND=dp)                                      :: reduced

      is_conjugate = .TRUE.
      DO idir = 1, 3
         IF (cell%perd(idir) /= 0) THEN
            reduced = xk1(idir) + xk2(idir)
            reduced = reduced - REAL(NINT(reduced), KIND=dp)
            IF (ABS(reduced) > eps_kpoint_real) THEN
               is_conjugate = .FALSE.
               EXIT
            END IF
         END IF
      END DO
   END FUNCTION is_conjugate_kpoint

! **************************************************************************************************
!> \brief Returns true for Gamma/BZ-edge k-points where real Bloch orbitals can be used.
!> \param cell ...
!> \param xk ...
!> \param eps_kpoint_real ...
!> \return ...
! **************************************************************************************************
   FUNCTION is_real_kpoint(cell, xk, eps_kpoint_real) RESULT(is_real)
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xk
      REAL(KIND=dp), INTENT(IN)                          :: eps_kpoint_real
      LOGICAL                                            :: is_real

      INTEGER                                            :: idir
      REAL(KIND=dp)                                      :: reduced

      is_real = .TRUE.
      DO idir = 1, 3
         IF (cell%perd(idir) /= 0) THEN
            reduced = xk(idir) - REAL(NINT(xk(idir)), KIND=dp)
            IF (ABS(reduced) > eps_kpoint_real .AND. &
                ABS(ABS(reduced) - 0.5_dp) > eps_kpoint_real) is_real = .FALSE.
         END IF
      END DO
   END FUNCTION is_real_kpoint

! **************************************************************************************************
!> \brief Write all non-orbital CASINO gwfn.data sections.
!> \param iw ...
!> \param periodicity ...
!> \param nspins ...
!> \param e_nn ...
!> \param nel_tot ...
!> \param natoms ...
!> \param coord ...
!> \param atomic_number ...
!> \param valence_charge ...
!> \param cell ...
!> \param periodic ...
!> \param nkp ...
!> \param nreal_k ...
!> \param kvec ...
!> \param shell_num ...
!> \param ao_num ...
!> \param prim_num ...
!> \param shell_ang_mom ...
!> \param shell_type ...
!> \param prim_per_shell ...
!> \param first_shell ...
!> \param exponents ...
!> \param coefficients ...
!> \param shell_position ...
! **************************************************************************************************
   SUBROUTINE write_casino_header(iw, periodicity, nspins, e_nn, nel_tot, natoms, coord, &
                                  atomic_number, valence_charge, cell, periodic, nkp, nreal_k, kvec, &
                                  shell_num, ao_num, prim_num, shell_ang_mom, shell_type, prim_per_shell, &
                                  first_shell, exponents, coefficients, shell_position)
      INTEGER, INTENT(IN)                                :: iw, periodicity, nspins
      REAL(KIND=dp), INTENT(IN)                          :: e_nn
      INTEGER, INTENT(IN)                                :: nel_tot, natoms
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: coord
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atomic_number
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: valence_charge
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      LOGICAL, INTENT(IN)                                :: periodic
      INTEGER, INTENT(IN)                                :: nkp, nreal_k
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: kvec
      INTEGER, INTENT(IN)                                :: shell_num, ao_num, prim_num
      INTEGER, DIMENSION(:), INTENT(IN)                  :: shell_ang_mom, shell_type, &
                                                            prim_per_shell, first_shell
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: exponents, coefficients
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: shell_position

      INTEGER                                            :: highest_ang_mom, i

      highest_ang_mom = MAXVAL(shell_ang_mom) + 1

      WRITE (iw, '(A)') "CP2K CASINO gwfn.data"
      WRITE (iw, *)
      WRITE (iw, '(A)') "BASIC INFO"
      WRITE (iw, '(A)') "----------"
      WRITE (iw, '(A)') "Generated by:"
      WRITE (iw, '(1X,A)') TRIM(cp2k_version)
      WRITE (iw, '(A)') "Method:"
      WRITE (iw, '(A)') " DFT"
      WRITE (iw, '(A)') "DFT functional:"
      WRITE (iw, '(A)') " CP2K"
      WRITE (iw, '(A)') "Periodicity:"
      WRITE (iw, '(1X,I0)') periodicity
      WRITE (iw, '(A)') "Spin unrestricted:"
      WRITE (iw, '(1X,A)') MERGE(".true. ", ".false.", nspins > 1)
      WRITE (iw, '(A)') "Nuclear repulsion energy (au/atom):"
      WRITE (iw, '(1PE20.13)') e_nn
      WRITE (iw, '(A)') "Number of electrons per primitive cell"
      WRITE (iw, '(1X,I0)') nel_tot
      WRITE (iw, *)

      WRITE (iw, '(A)') "GEOMETRY"
      WRITE (iw, '(A)') "--------"
      WRITE (iw, '(A)') "Number of atoms"
      WRITE (iw, '(1X,I0)') natoms
      WRITE (iw, '(A)') "Atomic positions (au)"
      DO i = 1, natoms
         WRITE (iw, '(3(1PE20.13))') coord(:, i)
      END DO
      WRITE (iw, '(A)') "Atomic numbers for each atom"
      CALL write_integer_vector(iw, atomic_number)
      WRITE (iw, '(A)') "Valence charges for each atom"
      CALL write_real_vector(iw, valence_charge)
      IF (.NOT. periodic) WRITE (iw, *)
      IF (periodic) THEN
         WRITE (iw, '(A)') "Primitive lattice vectors (au)"
         DO i = 1, 3
            WRITE (iw, '(3(1PE20.13))') cell%hmat(:, i)
         END DO
         WRITE (iw, *)

         WRITE (iw, '(A)') "K SPACE NET"
         WRITE (iw, '(A)') "-----------"
         WRITE (iw, '(A)') "Number of k points"
         WRITE (iw, '(1X,I0)') nkp
         WRITE (iw, '(A)') "Number of 'real' k points on BZ edge"
         WRITE (iw, '(1X,I0)') nreal_k
         WRITE (iw, '(A)') "k point coordinates (au)"
         DO i = 1, nkp
            WRITE (iw, '(3(1PE20.13))') kvec(:, i)
         END DO
         WRITE (iw, *)
      END IF

      WRITE (iw, '(A)') "BASIS SET"
      WRITE (iw, '(A)') "---------"
      WRITE (iw, '(A)') "Number of Gaussian centres"
      WRITE (iw, '(1X,I0)') natoms
      WRITE (iw, '(A)') "Number of shells per primitive cell"
      WRITE (iw, '(1X,I0)') shell_num
      WRITE (iw, '(A)') "Number of basis functions ('AO') per primitive cell"
      WRITE (iw, '(1X,I0)') ao_num
      WRITE (iw, '(A)') "Number of Gaussian primitives per primitive cell"
      WRITE (iw, '(1X,I0)') prim_num
      WRITE (iw, '(A)') "Highest shell angular momentum (s/p/d/f/g... 1/2/3/4/5...)"
      WRITE (iw, '(1X,I0)') highest_ang_mom
      WRITE (iw, '(A)') "Code for shell types (s/sp/p/d/f... 1/2/3/4/5...)"
      CALL write_integer_vector(iw, shell_type)
      WRITE (iw, '(A)') "Number of primitive Gaussians in each shell"
      CALL write_integer_vector(iw, prim_per_shell)
      WRITE (iw, '(A)') "Sequence number of first shell on each centre"
      CALL write_integer_vector(iw, first_shell)
      WRITE (iw, '(A)') "Exponents of Gaussian primitives"
      CALL write_real_vector(iw, exponents)
      WRITE (iw, '(A)') "Correctly normalised contraction coefficients"
      CALL write_real_vector(iw, coefficients)
      WRITE (iw, '(A)') "Position of each shell (au)"
      DO i = 1, shell_num
         WRITE (iw, '(3(1PE20.13))') shell_position(:, i)
      END DO
      WRITE (iw, *)

      WRITE (iw, '(A)') "MULTIDETERMINANT INFORMATION"
      WRITE (iw, '(A)') "----------------------------"
      WRITE (iw, '(A)') "GS"
      WRITE (iw, *)
      WRITE (iw, '(A)') "ORBITAL COEFFICIENTS"
      WRITE (iw, '(A)') "---------------------------"
   END SUBROUTINE write_casino_header

! **************************************************************************************************
!> \brief Write one CASINO MO block for a spin/k-point.
!> \param iw ...
!> \param mos_sgf ...
!> \param mos_sgf_im ...
!> \param cp2k_to_casino_ao ...
!> \param mo_scale ...
!> \param ao_num ...
!> \param complex_orbitals ...
! **************************************************************************************************
   SUBROUTINE write_casino_orbitals(iw, mos_sgf, mos_sgf_im, cp2k_to_casino_ao, mo_scale, ao_num, &
                                    complex_orbitals)
      INTEGER, INTENT(IN)                                :: iw
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: mos_sgf, mos_sgf_im
      INTEGER, DIMENSION(:), INTENT(IN)                  :: cp2k_to_casino_ao
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: mo_scale
      INTEGER, INTENT(IN)                                :: ao_num
      LOGICAL, INTENT(IN)                                :: complex_orbitals

      INTEGER                                            :: iao, imo, nbuffer
      REAL(KIND=dp), DIMENSION(4)                        :: buffer

      nbuffer = 0
      buffer(:) = 0.0_dp
      DO imo = 1, ao_num
         DO iao = 1, ao_num
            CALL push_real(iw, buffer, nbuffer, mo_scale(iao)*mos_sgf(cp2k_to_casino_ao(iao), imo))
            IF (complex_orbitals) THEN
               CALL push_real(iw, buffer, nbuffer, mo_scale(iao)*mos_sgf_im(cp2k_to_casino_ao(iao), imo))
            END IF
         END DO
      END DO
      IF (nbuffer > 0) WRITE (iw, '(4(1PE20.13))') buffer(1:nbuffer)
   END SUBROUTINE write_casino_orbitals

! **************************************************************************************************
!> \brief Append one real number to a four-column output buffer.
!> \param iw ...
!> \param buffer ...
!> \param nbuffer ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE push_real(iw, buffer, nbuffer, value)
      INTEGER, INTENT(IN)                                :: iw
      REAL(KIND=dp), DIMENSION(4), INTENT(INOUT)         :: buffer
      INTEGER, INTENT(INOUT)                             :: nbuffer
      REAL(KIND=dp), INTENT(IN)                          :: value

      nbuffer = nbuffer + 1
      buffer(nbuffer) = value
      IF (nbuffer == SIZE(buffer)) THEN
         WRITE (iw, '(4(1PE20.13))') buffer
         nbuffer = 0
      END IF
   END SUBROUTINE push_real

! **************************************************************************************************
!> \brief Write an integer vector in CASINO-friendly fixed-width columns.
!> \param iw ...
!> \param values ...
! **************************************************************************************************
   SUBROUTINE write_integer_vector(iw, values)
      INTEGER, INTENT(IN)                                :: iw
      INTEGER, DIMENSION(:), INTENT(IN)                  :: values

      INTEGER                                            :: i, ilast

      DO i = 1, SIZE(values), 8
         ilast = MIN(i + 7, SIZE(values))
         WRITE (iw, '(8I10)') values(i:ilast)
      END DO
   END SUBROUTINE write_integer_vector

! **************************************************************************************************
!> \brief Write a real vector in CASINO-friendly fixed-width columns.
!> \param iw ...
!> \param values ...
! **************************************************************************************************
   SUBROUTINE write_real_vector(iw, values)
      INTEGER, INTENT(IN)                                :: iw
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: values

      INTEGER                                            :: i, ilast

      DO i = 1, SIZE(values), 4
         ilast = MIN(i + 3, SIZE(values))
         WRITE (iw, '(4(1PE20.13))') values(i:ilast)
      END DO
   END SUBROUTINE write_real_vector

! **************************************************************************************************
!> \brief Rotate a complex value by exp(-i*k.g) using the TREXIO gauge convention.
!> \param re ...
!> \param im ...
!> \param cval ...
!> \param sval ...
! **************************************************************************************************
   SUBROUTINE rotate_complex_pair(re, im, cval, sval)
      REAL(KIND=dp), INTENT(INOUT)                       :: re, im
      REAL(KIND=dp), INTENT(IN)                          :: cval, sval

      REAL(KIND=dp)                                      :: im_old, re_old

      re_old = re
      im_old = im
      re = cval*re_old + sval*im_old
      im = -sval*re_old + cval*im_old
   END SUBROUTINE rotate_complex_pair

! **************************************************************************************************
!> \brief Convert CP2K normalized primitive data to CASINO contraction coefficients.
!> \param l ...
!> \param nprim ...
!> \param zetas ...
!> \param gcc ...
!> \param exponents ...
!> \param coefficients ...
! **************************************************************************************************
   SUBROUTINE casino_shell_coefficients(l, nprim, zetas, gcc, exponents, coefficients)
      INTEGER, INTENT(IN)                                :: l, nprim
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetas, gcc
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: exponents, coefficients

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: contraction_norm, expzet, prefac, &
                                                            prim_cart_fac
      REAL(KIND=dp), DIMENSION(nprim)                    :: raw_coeff

      expzet = 0.25_dp*REAL(2*l + 3, KIND=dp)
      prefac = 2.0_dp**l*(2.0_dp/pi)**0.75_dp
      DO i = 1, nprim
         prim_cart_fac = prefac*zetas(i)**expzet
         raw_coeff(i) = gcc(i)/prim_cart_fac
      END DO
      contraction_norm = casino_contraction_norm(l, nprim, zetas, raw_coeff)
      DO i = 1, nprim
         exponents(i) = zetas(i)
         coefficients(i) = raw_coeff(i)*contraction_norm*casino_primitive_norm(l, zetas(i))
      END DO
   END SUBROUTINE casino_shell_coefficients

! **************************************************************************************************
!> \brief Whole-contraction normalization used by CASINO's molden2qmc converter.
!> \param l ...
!> \param nprim ...
!> \param zetas ...
!> \param raw_coeff ...
!> \return ...
! **************************************************************************************************
   FUNCTION casino_contraction_norm(l, nprim, zetas, raw_coeff) RESULT(norm)
      INTEGER, INTENT(IN)                                :: l, nprim
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetas, raw_coeff
      REAL(KIND=dp)                                      :: norm

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: overlap

      overlap = 0.0_dp
      DO i = 1, nprim
         DO j = 1, nprim
            overlap = overlap + raw_coeff(i)*raw_coeff(j)* &
                      (2.0_dp*SQRT(zetas(i)*zetas(j))/(zetas(i) + zetas(j)))**(l + 1.5_dp)
         END DO
      END DO
      norm = 1.0_dp/SQRT(overlap)
   END FUNCTION casino_contraction_norm

! **************************************************************************************************
!> \brief Primitive m-independent normalization used by CASINO's Gaussian evaluator.
!> \param l ...
!> \param alpha ...
!> \return ...
! **************************************************************************************************
   FUNCTION casino_primitive_norm(l, alpha) RESULT(norm)
      INTEGER, INTENT(IN)                                :: l
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      REAL(KIND=dp)                                      :: norm

      norm = SQRT(2.0_dp**(l + 1.5_dp)*alpha**(l + 1.5_dp))/pi**0.75_dp
      IF (l > 0) norm = norm*SQRT(2.0_dp**l/odd_double_factorial(2*l - 1))
   END FUNCTION casino_primitive_norm

! **************************************************************************************************
!> \brief CASINO shell type code.
!> \param l ...
!> \return ...
! **************************************************************************************************
   FUNCTION casino_shell_type(l) RESULT(shell_type)
      INTEGER, INTENT(IN)                                :: l
      INTEGER                                            :: shell_type

      IF (l == 0) THEN
         shell_type = 1
      ELSE
         shell_type = l + 2
      END IF
   END FUNCTION casino_shell_type

! **************************************************************************************************
!> \brief CP2K AO index for a CASINO/MOLDEN ordered harmonic shell.
!> \param l ...
!> \param k ...
!> \return ...
! **************************************************************************************************
   FUNCTION casino_cp2k_index(l, k) RESULT(idx)
      INTEGER, INTENT(IN)                                :: l, k
      INTEGER                                            :: idx

      INTEGER, DIMENSION(9, 0:max_casino_l), PARAMETER :: map = RESHAPE([1, 0, 0, 0, 0, 0, 0, 0, 0 &
         , 3, 1, 2, 0, 0, 0, 0, 0, 0, 3, 4, 2, 5, 1, 0, 0, 0, 0, 4, 5, 3, 6, 2, 7, 1, 0, 0, 5, 6, 4&
         , 7, 3, 8, 2, 9, 1], [9, max_casino_l + 1])

      idx = map(k, l)
   END FUNCTION casino_cp2k_index

! **************************************************************************************************
!> \brief Scale factors converting MOLDEN harmonic MO coefficients to CASINO conventions.
!> \param l ...
!> \param k ...
!> \return ...
! **************************************************************************************************
   FUNCTION casino_mo_scale(l, k) RESULT(scale)
      INTEGER, INTENT(IN)                                :: l, k
      REAL(KIND=dp)                                      :: scale

      REAL(KIND=dp), DIMENSION(5), PARAMETER :: d_factor = [0.5_dp, 3.0_dp, 3.0_dp, 3.0_dp, 6.0_dp]

      INTEGER                                            :: m

      IF (l <= 1) THEN
         scale = 1.0_dp
      ELSE
         m = casino_m_quantum_number(k)
         scale = casino_m_dependent_factor(l, m)
         IF (l == 2) scale = scale*d_factor(k)
      END IF
   END FUNCTION casino_mo_scale

! **************************************************************************************************
!> \brief m sequence in CASINO/MOLDEN harmonic order: 0,+1,-1,+2,-2,...
!> \param k ...
!> \return ...
! **************************************************************************************************
   FUNCTION casino_m_quantum_number(k) RESULT(m)
      INTEGER, INTENT(IN)                                :: k
      INTEGER                                            :: m

      IF (k == 1) THEN
         m = 0
      ELSE IF (MOD(k, 2) == 0) THEN
         m = k/2
      ELSE
         m = -(k/2)
      END IF
   END FUNCTION casino_m_quantum_number

! **************************************************************************************************
!> \brief CASINO m-dependent normalization factor.
!> \param l ...
!> \param m ...
!> \return ...
! **************************************************************************************************
   FUNCTION casino_m_dependent_factor(l, m) RESULT(factor)
      INTEGER, INTENT(IN)                                :: l, m
      REAL(KIND=dp)                                      :: factor

      INTEGER                                            :: am
      REAL(KIND=dp)                                      :: prefactor

      am = ABS(m)
      prefactor = MERGE(1.0_dp, 2.0_dp, am == 0)
      factor = SQRT(prefactor*factorial(l - am)/factorial(l + am))
   END FUNCTION casino_m_dependent_factor

! **************************************************************************************************
!> \brief Real factorial for small non-negative integers.
!> \param n ...
!> \return ...
! **************************************************************************************************
   FUNCTION factorial(n) RESULT(value)
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp)                                      :: value

      INTEGER                                            :: i

      value = 1.0_dp
      DO i = 2, n
         value = value*REAL(i, KIND=dp)
      END DO
   END FUNCTION factorial

! **************************************************************************************************
!> \brief Odd double factorial.
!> \param n ...
!> \return ...
! **************************************************************************************************
   FUNCTION odd_double_factorial(n) RESULT(value)
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp)                                      :: value

      INTEGER                                            :: i

      value = 1.0_dp
      DO i = MAX(1, n), 1, -2
         value = value*REAL(i, KIND=dp)
      END DO
   END FUNCTION odd_double_factorial

! **************************************************************************************************
!> \brief Computes the nuclear repulsion energy of a molecular system.
!> \param particle_set ...
!> \param kind_set ...
!> \param e_nn ...
! **************************************************************************************************
   SUBROUTINE nuclear_repulsion_energy(particle_set, kind_set, e_nn)
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: kind_set
      REAL(KIND=dp), INTENT(OUT)                         :: e_nn

      INTEGER                                            :: i, ikind, j, jkind, natoms
      REAL(KIND=dp)                                      :: r_ij, zeff_i, zeff_j

      natoms = SIZE(particle_set)
      e_nn = 0.0_dp
      DO i = 1, natoms
         CALL get_atomic_kind(particle_set(i)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(kind_set(ikind), zeff=zeff_i)
         DO j = i + 1, natoms
            r_ij = NORM2(particle_set(i)%r - particle_set(j)%r)
            CALL get_atomic_kind(particle_set(j)%atomic_kind, kind_number=jkind)
            CALL get_qs_kind(kind_set(jkind), zeff=zeff_j)
            e_nn = e_nn + zeff_i*zeff_j/r_ij
         END DO
      END DO
   END SUBROUTINE nuclear_repulsion_energy

! **************************************************************************************************
!> \brief Computes the CASINO-compatible 3D periodic nuclear repulsion energy.
!> \param cell ...
!> \param periodicity ...
!> \param coord ...
!> \param charge ...
!> \param e_nn ...
! **************************************************************************************************
   SUBROUTINE periodic_nuclear_repulsion_energy(cell, periodicity, coord, charge, e_nn)
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      INTEGER, INTENT(IN)                                :: periodicity
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: coord
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: charge
      REAL(KIND=dp), INTENT(OUT)                         :: e_nn

      INTEGER                                            :: gmax, i, ig1, ig2, ig3, j, n1, n2, n3, &
                                                            natoms, nmax
      REAL(KIND=dp) :: alpha, alpha2, cutoff_arg, g_cut, g_sq, min_g, min_h, neut_energy, phase, &
         r, real_cut, real_energy, recip_energy, self_energy, struc_im, struc_re, volume
      REAL(KIND=dp), DIMENSION(3)                        :: delta, g_index, gvec, lattice_shift

      e_nn = 0.0_dp
      IF (periodicity /= 3) RETURN

      volume = ABS(cell%deth)
      IF (volume <= 0.0_dp) CPABORT("CASINO periodic nuclear repulsion requires a non-zero cell volume.")

      natoms = SIZE(charge)
      IF (natoms == 0) RETURN

      min_h = HUGE(1.0_dp)
      min_g = HUGE(1.0_dp)
      DO i = 1, 3
         min_h = MIN(min_h, NORM2(cell%hmat(:, i)))
         g_index = 0.0_dp
         g_index(i) = 1.0_dp
         gvec = 2.0_dp*pi*MATMUL(TRANSPOSE(cell%h_inv), g_index)
         min_g = MIN(min_g, NORM2(gvec))
      END DO
      IF (min_h <= 0.0_dp .OR. min_g <= 0.0_dp) THEN
         CPABORT("CASINO periodic nuclear repulsion requires non-zero lattice vectors.")
      END IF

      cutoff_arg = SQRT(-LOG(1.0E-12_dp))
      alpha = SQRT(pi)*(REAL(natoms, KIND=dp)/volume)**(1.0_dp/3.0_dp)
      alpha2 = alpha*alpha
      real_cut = cutoff_arg/alpha
      g_cut = 2.0_dp*alpha*cutoff_arg
      nmax = MAX(1, CEILING(real_cut/min_h) + 1)
      gmax = MAX(1, CEILING(g_cut/min_g) + 1)

      real_energy = 0.0_dp
      DO i = 1, natoms
         DO j = 1, natoms
            DO n1 = -nmax, nmax
               DO n2 = -nmax, nmax
                  DO n3 = -nmax, nmax
                     IF (i == j .AND. n1 == 0 .AND. n2 == 0 .AND. n3 == 0) CYCLE
                     lattice_shift = REAL(n1, KIND=dp)*cell%hmat(:, 1) + &
                                     REAL(n2, KIND=dp)*cell%hmat(:, 2) + &
                                     REAL(n3, KIND=dp)*cell%hmat(:, 3)
                     delta = coord(:, i) - coord(:, j) + lattice_shift
                     r = NORM2(delta)
                     IF (r <= real_cut) real_energy = real_energy + charge(i)*charge(j)*ERFC(alpha*r)/r
                  END DO
               END DO
            END DO
         END DO
      END DO
      real_energy = 0.5_dp*real_energy

      recip_energy = 0.0_dp
      DO ig1 = -gmax, gmax
         DO ig2 = -gmax, gmax
            DO ig3 = -gmax, gmax
               IF (ig1 == 0 .AND. ig2 == 0 .AND. ig3 == 0) CYCLE
               g_index = [REAL(ig1, KIND=dp), REAL(ig2, KIND=dp), REAL(ig3, KIND=dp)]
               gvec = 2.0_dp*pi*MATMUL(TRANSPOSE(cell%h_inv), g_index)
               g_sq = DOT_PRODUCT(gvec, gvec)
               IF (SQRT(g_sq) > g_cut) CYCLE
               struc_re = 0.0_dp
               struc_im = 0.0_dp
               DO i = 1, natoms
                  phase = DOT_PRODUCT(gvec, coord(:, i))
                  struc_re = struc_re + charge(i)*COS(phase)
                  struc_im = struc_im + charge(i)*SIN(phase)
               END DO
               recip_energy = recip_energy + EXP(-g_sq/(4.0_dp*alpha2))/g_sq* &
                              (struc_re*struc_re + struc_im*struc_im)
            END DO
         END DO
      END DO
      recip_energy = 2.0_dp*pi*recip_energy/volume

      self_energy = -alpha*SUM(charge*charge)/SQRT(pi)
      neut_energy = -pi*SUM(charge)**2/(2.0_dp*alpha2*volume)
      e_nn = real_energy + recip_energy + self_energy + neut_energy
   END SUBROUTINE periodic_nuclear_repulsion_energy

END MODULE casino_utils
