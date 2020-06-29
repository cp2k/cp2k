!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to chargemol
!> \par History
!>      06.2020 created [Hossam Elgabarty]
!> \author Hossam Elgabarty
! **************************************************************************************************
MODULE qs_chargemol
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_files,                        ONLY: open_file
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_get_submatrix
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_generate_filename
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE machine,                         ONLY: m_datum
   USE orbital_pointers,                ONLY: nco,&
                                              nso
   USE orbital_symbols,                 ONLY: cgf_symbol,&
                                              sgf_symbol
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: get_ptable_info
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: mo_set_p_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_subsys_types,                 ONLY: qs_subsys_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_chargemol'
   LOGICAL, PARAMETER                   :: debug_this_module = .FALSE.

   INTEGER, PARAMETER                   :: chargemol_lmax = 5

   PUBLIC :: write_wfx

CONTAINS

   ! **************************************************************************************************
   !> \brief Write chargemol wfx file
   !> \note  The subroutine works only for nspins = 1, execution is aborted otherwise
   !> \note  all nspins = 2 code in only a placeholder for future implementation
   !> \param qs_env ...
   !> \param dft_section ...
   !> \par History
   !>      06.2020 created [Hossam Elgabarty]
   !> \author Hossam Elgabarty
   ! **************************************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dft_section ...
! **************************************************************************************************
   SUBROUTINE write_wfx(qs_env, dft_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_wfx', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: gcc_eps = 1.0E-7_dp

      CHARACTER(LEN=12), DIMENSION(:), POINTER           :: bcgf_symbol
      CHARACTER(len=2)                                   :: asym
      CHARACTER(len=20), ALLOCATABLE, DIMENSION(:)       :: atom_symbols
      CHARACTER(LEN=256)                                 :: datx
      CHARACTER(LEN=6), DIMENSION(:), POINTER            :: bsgf_symbol
      CHARACTER(len=default_string_length)               :: filename
      INTEGER :: handle, i, iatom, icgf, ico, ikind, iloop, imo, ipgf, irow, iset, isgf, ishell, &
         ispin, iwfx, lshell, max_l, ncgf, ncol_global, nelec_alpha, nelec_beta, noccup, &
         nrow_global, nset, nsgf, nspins, num_atoms, output_unit, roks_high, roks_low, tot_npgf = 0
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, atomic_numbers, kind_of
      INTEGER, DIMENSION(:), POINTER                     :: lmax, npgf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: l
      LOGICAL                                            :: do_kpoints, newl, periodic
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atoms_cart, cmatrix, smatrix
      REAL(KIND=dp), DIMENSION(:), POINTER               :: core_charges
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: chargemol_section

!REAL(KIND=dp), PARAMETER                           :: gcc_eps = EPSILON(1.0_dp)

      ! !--------------------------------------------------------------------------------------------!

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(/,T2,A)') &
            '!-----------------------------------------------------------------------------!'
         WRITE (output_unit, '(T32,A)') "Writing Chargemol wfx file..."
         WRITE (output_unit, '(T2,A)') &
            '!-----------------------------------------------------------------------------!'
      END IF

      ! Collect environment info
      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, cell=cell, mos=mos, particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, scf_env=scf_env, scf_control=scf_control, &
                      dft_control=dft_control, energy=energy, force=force, subsys=subsys, &
                      atomic_kind_set=atomic_kind_set, do_kpoints=do_kpoints)

      IF (scf_control%use_ot) THEN
         CPABORT("OT is incompatible with .wfx output. Switch off OT.")
      END IF

      IF (dft_control%lsd .OR. dft_control%do_tddfpt_calculation) THEN
         CALL cp_abort(__LOCATION__, "The output of chargemol . wfx files is currently"// &
                       "implemented only for closed-shell wavefunctions.")
      END IF

      chargemol_section => section_vals_get_subs_vals(dft_section, "PRINT%CHARGEMOL")

      CALL section_vals_val_get(chargemol_section, "PERIODIC", l_val=periodic)

      IF (do_kpoints) THEN
         CALL cp_abort(__LOCATION__, "The output of chargemol . wfx files is currently"// &
                       "not implemented for k-points (only for gamma point).")
      END IF

    !!---------------------------------------------------------------------------------!
      ! output unit
    !!---------------------------------------------------------------------------------!

      filename = cp_print_key_generate_filename(logger, chargemol_section, middle_name="CHARGEMOL", &
                                                extension=".wfx", my_local=.FALSE.)

      CALL open_file(filename, file_status="UNKNOWN", &
                     file_action="WRITE", &
                     unit_number=iwfx)

    !!---------------------------------------------------------------------------------!

      IF (iwfx > 0) THEN

         nspins = SIZE(mos)
         IF (nspins == 1) THEN
            nelec_alpha = mos(1)%mo_set%nelectron/2
            nelec_beta = nelec_alpha
         ELSE
            nelec_alpha = mos(1)%mo_set%nelectron
            nelec_beta = mos(2)%mo_set%nelectron
         END IF

         IF (dft_control%roks) THEN
            IF (mos(1)%mo_set%homo > mos(2)%mo_set%homo) THEN
               roks_high = 1
               roks_low = 2
            ELSE
               roks_high = 2
               roks_low = 1
            END IF
         END IF

       !!---------------------------------------------------------------------------------!
         ! Initial parsing of topology and basis set, check maximum l .le. chargemol_lmax
       !!---------------------------------------------------------------------------------!
         num_atoms = SIZE(particle_set)
         ALLOCATE (atoms_cart(3, num_atoms))
         ALLOCATE (atom_symbols(num_atoms))
         ALLOCATE (atomic_numbers(num_atoms))
         ALLOCATE (core_charges(num_atoms))

         max_l = 0
         DO iatom = 1, num_atoms
            NULLIFY (orb_basis_set)
            atoms_cart(1:3, iatom) = particle_set(iatom)%r(1:3)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, &
                                 element_symbol=asym, kind_number=ikind)
            atom_symbols(iatom) = asym
            CALL get_ptable_info(asym, number=atomic_numbers(iatom))
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, &
                             core_charge=core_charges(iatom))
            IF (ASSOCIATED(orb_basis_set)) THEN
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      nset=nset, lmax=lmax)

               DO iset = 1, nset
                  max_l = MAX(max_l, lmax(iset))
               END DO
            ELSE
               CPABORT("Unknown basis set type. ")
            END IF
            IF (max_l > chargemol_lmax) THEN
               CALL cp_abort(__LOCATION__, "Chargemol supports basis sets with l"// &
                             " up to 5 only (h angular functions).")
            END IF
         END DO

         ! !===========================================================================!
         ! Header comment and Title
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "# Chargemol input file generated by CP2K "
         CALL m_datum(datx)
         WRITE (iwfx, "(A,/)") "# Creation date "//TRIM(datx)

         WRITE (iwfx, "(A)") "<Title>"
         WRITE (iwfx, *) TRIM(logger%iter_info%project_name)
         WRITE (iwfx, "(A,/)") "</Title>"

         ! !===========================================================================!
         ! Keywords
         ! TODO: check: always GTO??
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Keywords>"
         WRITE (iwfx, "(A)") "GTO"
         WRITE (iwfx, "(A,/)") "</Keywords>"

         ! !===========================================================================!
         ! Model
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "#<Model>"
         WRITE (iwfx, "(A)") "#"
         WRITE (iwfx, "(A,/)") "#</Model>"

         ! !===========================================================================!
         ! Number of nuclei
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Number of Nuclei>"
         WRITE (iwfx, "(I6)") num_atoms
         WRITE (iwfx, "(A,/)") "</Number of Nuclei>"

         ! !===========================================================================!
         ! Number of Occupied MOs
         ! Q: Is this really the best way? to rely on homo?
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Number of Occupied Molecular Orbitals>"
         noccup = 0
         IF (dft_control%roks .AND. nspins == 2) THEN
            noccup = MAX(mos(1)%mo_set%homo, mos(2)%mo_set%homo)
         ELSE
            DO ispin = 1, dft_control%nspins
               noccup = noccup + mos(ispin)%mo_set%homo
            END DO
         END IF
         WRITE (iwfx, "(2I6)") noccup
         WRITE (iwfx, "(A,/)") "</Number of Occupied Molecular Orbitals>"

         ! !===========================================================================!
         ! Number of Perturbations
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Number of Perturbations>"
         WRITE (iwfx, "(I6)") 0
         WRITE (iwfx, "(A,/)") "</Number of Perturbations>"

         ! !===========================================================================!
         ! Total charge
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Net Charge>"
         WRITE (iwfx, "(I6)") dft_control%charge
         WRITE (iwfx, "(A,/)") "</Net Charge>"

         ! !===========================================================================!
         ! Number of electrons
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Number of Electrons>"
         WRITE (iwfx, "(I6)") scf_env%nelectron
         WRITE (iwfx, "(A,/)") "</Number of Electrons>"

         ! !===========================================================================!
         ! Number of Alpha electrons
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Number of Alpha Electrons>"
         WRITE (iwfx, "(I6)") nelec_alpha
         WRITE (iwfx, "(A,/)") "</Number of Alpha Electrons>"

         ! !===========================================================================!
         ! Number of Beta electrons
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Number of Beta Electrons>"
         WRITE (iwfx, "(I6)") nelec_beta
         WRITE (iwfx, "(A,/)") "</Number of Beta Electrons>"

         ! !===========================================================================!
         ! Spin multiplicity
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Electron Spin Multiplicity>"
         WRITE (iwfx, "(I4)") dft_control%multiplicity
         WRITE (iwfx, "(A,/)") "</Electron Spin Multiplicity>"

         ! !===========================================================================!
         ! Number of Core Electrons
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Number of Core Electrons>"
         WRITE (iwfx, "(F20.12)") SUM(atomic_numbers) - SUM(core_charges)
         WRITE (iwfx, "(A,/)") "</Number of Core Electrons>"

         ! !===========================================================================!
         ! Nuclear Names
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Nuclear Names>"
         DO iatom = 1, num_atoms
            WRITE (iwfx, "(A4)") atom_symbols(iatom)
         END DO
         WRITE (iwfx, "(A,/)") "</Nuclear Names>"

         ! !===========================================================================!
         ! Atomic Numbers
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Atomic Numbers>"
         DO iatom = 1, num_atoms
            WRITE (iwfx, "(I4)") atomic_numbers(iatom)
         END DO
         WRITE (iwfx, "(A,/)") "</Atomic Numbers>"

         ! !===========================================================================!
         ! Nuclear charges
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Nuclear Charges>"
         DO iatom = 1, num_atoms
            WRITE (iwfx, "(F12.8)") core_charges(iatom)
         END DO
         WRITE (iwfx, "(A,/)") "</Nuclear Charges>"

         ! !===========================================================================!
         ! Nuclear coordinates
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Nuclear Cartesian Coordinates>"
         DO iatom = 1, num_atoms
            WRITE (iwfx, "(3ES26.16E3)") atoms_cart(1:3, iatom)
         END DO
         WRITE (iwfx, "(A,/)") "</Nuclear Cartesian Coordinates>"

         ! !===========================================================================!
         ! Periodic boundary conditions
         ! !===========================================================================!
         IF (periodic) THEN
            CALL get_cell(cell)
            IF (SUM(cell%perd) == 0) THEN
               CALL cp_warn(__LOCATION__, "Writing of periodic cell information"// &
                            " requested in input, but system is not periodic, ignoring request.")
            ELSE
               WRITE (iwfx, "(A)") "<Number of Translation Vectors>"
               WRITE (iwfx, "(I3)") SUM(cell%perd)
               WRITE (iwfx, "(A,/)") "</Number of Translation Vectors>"

               WRITE (iwfx, "(A)") "<Translation Vectors>"
               DO iatom = 1, 3
                  IF (cell%perd(iatom) == 1) THEN
                     WRITE (iwfx, "(3F12.6)") cell%hmat(iatom, 1:3)
                  END IF
               END DO
               WRITE (iwfx, "(A,/)") "</Translation Vectors>"
            END IF
            WRITE (iwfx, "(A)") "<Number of Kpoints>"
            WRITE (iwfx, "(I3)") 1
            WRITE (iwfx, "(A,/)") "</Number of Kpoints>"
            WRITE (iwfx, "(A)") "<Kpoint Weights>"
            WRITE (iwfx, "(I3)") 1
            WRITE (iwfx, "(A,/)") "</Kpoint Weights>"
            WRITE (iwfx, "(A)") "<Kpoint Fractional Coordinates>"
            WRITE (iwfx, "(F8.6)") 0.0
            WRITE (iwfx, "(A,/)") "</Kpoint Fractional Coordinates>"
         END IF

         ! !===========================================================================!
         ! Primitive centers
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Primitive Centers>"
         tot_npgf = 0
         DO iatom = 1, num_atoms
            iloop = 1
            NULLIFY (orb_basis_set)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
            IF (ASSOCIATED(orb_basis_set)) THEN
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      nset=nset, &
                                      npgf=npgf, &
                                      nshell=nshell, &
                                      l=l, &
                                      gcc=gcc)

               DO iset = 1, nset
                  DO ishell = 1, nshell(iset)
                     lshell = l(ishell, iset)
                     DO ico = 1, nco(lshell)
                        DO ipgf = 1, npgf(iset)
                           IF (ABS(gcc(ipgf, ishell, iset)) > gcc_eps) THEN
                              tot_npgf = tot_npgf + 1
                              IF (MOD(iloop + 10, 10) /= 0) THEN
                                 WRITE (iwfx, "(I4)", advance="no") iatom
                                 newl = .TRUE.
                              ELSE
                                 WRITE (iwfx, "(I4)") iatom
                                 newl = .FALSE.
                              END IF
                              iloop = iloop + 1
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
               IF (newl) THEN
                  WRITE (iwfx, *)
               END IF
            ELSE
               CPABORT("Unknown basis set type. ")
            END IF
         END DO
         WRITE (iwfx, "(A,/)") "</Primitive Centers>"

         ! !===========================================================================!
         ! Number of primitives
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Number of Primitives>"
         WRITE (iwfx, "(I8)") tot_npgf
         WRITE (iwfx, "(A,/)") "</Number of Primitives>"

         ! !===========================================================================!
         ! Primitive Types
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Primitive Types>"
         DO iatom = 1, num_atoms
            iloop = 1
            NULLIFY (orb_basis_set)
            NULLIFY (bcgf_symbol)
            NULLIFY (gcc)
            NULLIFY (zet)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
            IF (ASSOCIATED(orb_basis_set)) THEN
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      nset=nset, &
                                      nshell=nshell, &
                                      l=l, &
                                      npgf=npgf, &
                                      ncgf=ncgf, &
                                      zet=zet, &
                                      cgf_symbol=bcgf_symbol, &
                                      sgf_symbol=bsgf_symbol, &
                                      gcc=gcc)

               icgf = 1
               DO iset = 1, nset
                  DO ishell = 1, nshell(iset)
                     lshell = l(ishell, iset)
                     DO ico = 1, nco(lshell)
                        DO ipgf = 1, orb_basis_set%npgf(iset)
                           IF (ABS(gcc(ipgf, ishell, iset)) > gcc_eps) THEN
                              IF (MOD(iloop + 5, 5) /= 0) THEN
                                 WRITE (iwfx, '(I6)', advance="no") &
                                    pgf_type_chargemol(bcgf_symbol(icgf) (3:))
                                 !WRITE (iwfx, '(I2,2x,A6,2x,I2)', advance = "no") &
                                 !     pgf_type_chargemol(bcgf_symbol(icgf)(3:)),&
                                 !bcgf_symbol(icgf), icgf
                                 newl = .TRUE.
                              ELSE
                                 WRITE (iwfx, '(I6)') &
                                    pgf_type_chargemol(bcgf_symbol(icgf) (3:))
                                 !WRITE (iwfx, '(I2,2x,A6,2x,I2)') &
                                 !pgf_type_chargemol(bcgf_symbol(icgf)(3:)),&
                                 !bcgf_symbol(icgf), icgf
                                 newl = .FALSE.
                              END IF
                              iloop = iloop + 1
                           END IF
                        END DO
                        icgf = icgf + 1
                     END DO !ico
                  END DO ! ishell
               END DO ! iset

            ELSE
               CPABORT("Unknown basis set type.")
            END IF
            IF (newl) THEN
               WRITE (iwfx, *)
            END IF
         END DO ! iatom
         WRITE (iwfx, "(A,/)") "</Primitive Types>"

         ! !===========================================================================!
         ! Primitive Exponents
         ! !===========================================================================!
       !! NOTE: CP2K orders the atomic orbitals as follows:
       !! Cartesian orbitals: x (slowest loop), y, z (fastest loop)
       !! Spherical orbitals: m = -l, ..., 0, ..., +l

         WRITE (iwfx, "(A)") "<Primitive Exponents>"
         DO iatom = 1, num_atoms
            NULLIFY (orb_basis_set)
            NULLIFY (gcc)
            NULLIFY (zet)
            iloop = 1

            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
            IF (ASSOCIATED(orb_basis_set)) THEN
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      nset=nset, &
                                      npgf=npgf, &
                                      nshell=nshell, &
                                      l=l, &
                                      zet=zet, &
                                      gcc=gcc)

               DO iset = 1, nset
                  DO ishell = 1, nshell(iset)
                     lshell = l(ishell, iset)
                     DO ico = 1, nco(lshell)
                        DO ipgf = 1, npgf(iset)
                           IF (ABS(gcc(ipgf, ishell, iset)) > gcc_eps) THEN
                              IF (MOD(iloop + 5, 5) /= 0) THEN
                                 WRITE (iwfx, "(ES20.10)", advance="no") zet(ipgf, iset)
                                 newl = .TRUE.
                              ELSE
                                 WRITE (iwfx, "(ES20.10)") zet(ipgf, iset)
                                 newl = .FALSE.
                              END IF
                              iloop = iloop + 1
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
               IF (newl) THEN
                  WRITE (iwfx, *)
               END IF
            ELSE
               CPABORT("Unknown basis set type. ")
            END IF
         END DO
         WRITE (iwfx, "(A,/)") "</Primitive Exponents>"

         ! !===========================================================================!
         !  MO occupation numbers
         !  nonesense for roks at the moment
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Molecular Orbital Occupation Numbers>"
         IF (.NOT. dft_control%roks) THEN
            DO ispin = 1, nspins
               DO imo = 1, mos(ispin)%mo_set%homo
                  WRITE (iwfx, "(f8.6)") &
                     mos(ispin)%mo_set%occupation_numbers(imo)
               END DO
            END DO
         ELSE
            DO imo = 1, mos(roks_low)%mo_set%homo
               WRITE (iwfx, "(*(f8.6,1x))") &
                  mos(roks_low)%mo_set%occupation_numbers(imo) &
                  + mos(roks_low)%mo_set%occupation_numbers(imo)
            END DO
            DO imo = mos(roks_low)%mo_set%homo + 1, mos(roks_high)%mo_set%homo
               WRITE (iwfx, "(f8.6)") &
                  mos(roks_high)%mo_set%occupation_numbers(imo)
            END DO
         END IF
         WRITE (iwfx, "(A,/)") "</Molecular Orbital Occupation Numbers>"

         ! !===========================================================================!
         ! MO energies
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Molecular Orbital Energies>"
         DO ispin = 1, nspins
            DO imo = 1, mos(ispin)%mo_set%homo
               WRITE (iwfx, "(ES20.10)") mos(ispin)%mo_set%eigenvalues(imo)
            END DO
         END DO
         WRITE (iwfx, "(A,/)") "</Molecular Orbital Energies>"

         ! !===========================================================================!
         ! MO Spin types
         ! nonesense for ROKS
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Molecular Orbital Spin Types>"
         DO imo = 1, mos(1)%mo_set%homo
            IF (nspins == 1) THEN
               WRITE (iwfx, '(A15)') "Alpha and Beta"
            ELSEIF (dft_control%uks) THEN
               WRITE (iwfx, '(A6)') "Alpha"
            ELSE
               CALL cp_abort(__LOCATION__, "This calculation method is currently"// &
                             " not supported by the chargemol interface.")
            END IF
         END DO

         IF (nspins == 2) THEN
            DO imo = 1, mos(2)%mo_set%homo
               WRITE (iwfx, '(A5)') "Beta"
            END DO
         END IF

         WRITE (iwfx, "(A,/)") "</Molecular Orbital Spin Types>"

         ! !===========================================================================!
         ! Kpoint index of orbitals
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Kpoint Index for Orbitals>"
         DO ispin = 1, nspins
            DO imo = 1, mos(ispin)%mo_set%homo
               WRITE (iwfx, "(I6)") 1
            END DO
         END DO
         WRITE (iwfx, "(A,/)") "</Kpoint Index for Orbitals>"

         ! !===========================================================================!
         ! MO Primitive Coefficients
         ! broken for ROKS, should i just output mos(roks_high)?
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Molecular Orbital Primitive Coefficients>"
         NULLIFY (orb_basis_set)
         NULLIFY (gcc)

         IF (mos(1)%mo_set%use_mo_coeff_b) THEN
            DO ispin = 1, SIZE(mos)
               IF (.NOT. ASSOCIATED(mos(ispin)%mo_set%mo_coeff_b)) THEN
                  CPASSERT(.FALSE.)
               ENDIF
               CALL copy_dbcsr_to_fm(mos(ispin)%mo_set%mo_coeff_b, &
                                     mos(ispin)%mo_set%mo_coeff)
            ENDDO
         ENDIF

         DO ispin = 1, nspins
            CALL cp_fm_get_info(mos(ispin)%mo_set%mo_coeff, &
                                nrow_global=nrow_global, &
                                ncol_global=ncol_global)

            ALLOCATE (smatrix(nrow_global, ncol_global))

            CALL cp_fm_get_submatrix(mos(ispin)%mo_set%mo_coeff, smatrix)

            CALL get_qs_kind_set(qs_kind_set, ncgf=ncgf, nsgf=nsgf)

            ALLOCATE (cmatrix(ncgf, ncgf))

            cmatrix = 0.0_dp

            ! get MO coefficients in terms of contracted cartesian functions
            icgf = 1
            isgf = 1
            DO iatom = 1, num_atoms
               NULLIFY (orb_basis_set)
               CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), &
                                basis_set=orb_basis_set)
               IF (ASSOCIATED(orb_basis_set)) THEN
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                         nset=nset, &
                                         nshell=nshell, &
                                         l=l)
                  DO iset = 1, nset
                     DO ishell = 1, nshell(iset)
                        lshell = l(ishell, iset)
                        CALL dgemm("T", "N", nco(lshell), mos(ispin)%mo_set%nmo, &
                                   nso(lshell), 1.0_dp, &
                                   orbtramat(lshell)%s2c, nso(lshell), &
                                   smatrix(isgf, 1), nsgf, 0.0_dp, &
                                   cmatrix(icgf, 1), ncgf)
                        icgf = icgf + nco(lshell)
                        isgf = isgf + nso(lshell)
                     END DO
                  END DO
               ELSE
                  CPABORT("Unknown basis set type. ")
               END IF
            END DO ! iatom

            ! Now write MO coefficients in terms of cartesian primitives
            DO imo = 1, mos(ispin)%mo_set%homo

               WRITE (iwfx, "(A)") "<MO Number>"
               IF (nspins == 1 .OR. ispin == 1) THEN
                  WRITE (iwfx, "(I6)") imo
               ELSE
                  WRITE (iwfx, "(I6)") imo + mos(1)%mo_set%homo
               END IF
               WRITE (iwfx, "(A)") "</MO Number>"

               irow = 1
               DO iatom = 1, num_atoms
                  iloop = 1
                  NULLIFY (orb_basis_set)
                  CALL get_atomic_kind(particle_set(iatom)%atomic_kind, &
                                       kind_number=ikind)
                  CALL get_qs_kind(qs_kind_set(ikind), &
                                   basis_set=orb_basis_set)
                  IF (ASSOCIATED(orb_basis_set)) THEN
                     CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                            nset=nset, &
                                            nshell=nshell, &
                                            gcc=gcc, &
                                            l=l)

                     icgf = 1
                     DO iset = 1, nset
                        DO ishell = 1, nshell(iset)
                           lshell = l(ishell, iset)
                           DO ico = 1, nco(lshell)
                              DO ipgf = 1, orb_basis_set%npgf(iset)
                                 IF (ABS(gcc(ipgf, ishell, iset)) > gcc_eps) THEN
                                    IF (MOD(iloop + 5, 5) /= 0) THEN
                                       WRITE (iwfx, '(ES20.10)', advance="no") &
                                          cmatrix(irow, imo)* &
                                          orb_basis_set%norm_cgf(icgf)* &
                                          orb_basis_set%gcc(ipgf, ishell, iset)
                                       newl = .TRUE.
                                    ELSE
                                       WRITE (iwfx, '(ES20.10)') cmatrix(irow, imo)* &
                                          orb_basis_set%norm_cgf(icgf)* &
                                          orb_basis_set%gcc(ipgf, ishell, iset)
                                       newl = .FALSE.
                                    END IF
                                    iloop = iloop + 1
                                 END IF
                              END DO
                              irow = irow + 1
                              !WRITE bcgf_symbol(icgf)
                           END DO !ico
                           icgf = icgf + 1
                        END DO ! ishell
                     END DO ! iset

                  ELSE
                     CPABORT("Unknown basis set type.")
                  END IF
                  IF (newl) THEN
                     WRITE (iwfx, *)
                  END IF
               END DO ! iatom
               !Write (iwfx,*)
            END DO ! imo

            IF (ALLOCATED(cmatrix)) DEALLOCATE (cmatrix)
            IF (ALLOCATED(smatrix)) DEALLOCATE (smatrix)

         END DO ! ispin
         WRITE (iwfx, "(A,/)") "</Molecular Orbital Primitive Coefficients>"

         ! !===========================================================================!
         ! Energy
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Energy>"
         WRITE (iwfx, "(ES26.16E3)") energy%total
         WRITE (iwfx, "(A,/)") "</Energy>"

         ! !===========================================================================!
         ! Virial ratio
         ! !===========================================================================!
         WRITE (iwfx, "(A)") "<Virial Ratio (-V/T)>"
         WRITE (iwfx, "(ES20.12)") - 1*(energy%total - energy%kinetic)/energy%kinetic
         WRITE (iwfx, "(A,/)") "</Virial Ratio (-V/T)>"

         ! !===========================================================================!
         ! Force-related quantities
         ! inactivated for now
         ! !===========================================================================!
         IF (ASSOCIATED(force) .AND. debug_this_module) THEN
            WRITE (iwfx, "(A)") "<Nuclear Cartesian Energy Gradients>"

            CALL get_qs_env(qs_env, force=force, atomic_kind_set=atomic_kind_set)
            ALLOCATE (atom_of_kind(num_atoms))
            ALLOCATE (kind_of(num_atoms))
            CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                     atom_of_kind=atom_of_kind, &
                                     kind_of=kind_of)

            DO iatom = 1, num_atoms
               ikind = kind_of(iatom)
               i = atom_of_kind(iatom)
               WRITE (iwfx, "(3ES20.12E3)") force(ikind)%total(1:3, i)
            END DO

            WRITE (iwfx, "(A,/)") "<Nuclear Cartesian Energy Gradients>"

            ! W
            WRITE (iwfx, "(A)") "<Nuclear Virial of Energy-Gradient-Based Forces on Nuclei, W>"
            WRITE (iwfx, "(A)") ""
            WRITE (iwfx, "(A,/)") "<Nuclear Virial of Energy-Gradient-Based Forces on Nuclei, W>"

            ! Full Virial ratio
            WRITE (iwfx, "(A)") "<Full Virial Ratio, -(V - W)/T)>"
            WRITE (iwfx, "(A)") ""
            WRITE (iwfx, "(A,/)") "</Full Virial Ratio, -(V - W)/T)>"

         END IF

         ! !===========================================================================!
         ! Core-density
         ! !===========================================================================!
         ! Additional Electron Density Function (EDF)
         ! WRITE (iwfx, "(A)") "#<Additional Electron Density Function (EDF)>"
         ! WRITE (iwfx, "(A,/)") "#</Additional Electron Density Function (EDF)>"

         ! !-----------------------------------!
         ! Done
         ! !-----------------------------------!

      ELSE  ! no output unit
         iwfx = -1
      END IF

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(/,T2,A)') &
            '!--------------------------- Chargemol wfx written ---------------------------!'
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_wfx

   ! ! ********************************************************************************
   ! !> \brief Figure out type of cartesian primitive gaussian (S, PX, PY ... HXXXXX)
   ! !> \param  ...
   ! !> \return ...
   ! ! ********************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param l_symb ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION pgf_type_chargemol(l_symb) RESULT(label)
      !INTEGER, INTENT(OUT)                                 :: label
      CHARACTER(len=*), INTENT(IN)                       :: l_symb

      SELECT CASE (l_symb)
      CASE ("s")
         label = 1
      CASE ("px")
         label = 2
      CASE ("py")
         label = 3
      CASE ("pz")
         label = 4
      CASE ("dx2")
         label = 5
      CASE ("dy2")
         label = 6
      CASE ("dz2")
         label = 7
      CASE ("dxy")
         label = 8
      CASE ("dxz")
         label = 9
      CASE ("dyz")
         label = 10
      CASE ("fx3")
         label = 11
      CASE ("fy3")
         label = 12
      CASE ("fz3")
         label = 13
      CASE ("fx2y")
         label = 14
      CASE ("fx2z")
         label = 15
      CASE ("fy2z")
         label = 16
      CASE ("fxy2")
         label = 17
      CASE ("fxz2")
         label = 18
      CASE ("fyz2")
         label = 19
      CASE ("fxyz")
         label = 20
      CASE ("gx4")
         label = 21
      CASE ("gy4")
         label = 22
      CASE ("gz4")
         label = 23
      CASE ("gx3y")
         label = 24
      CASE ("gx3z")
         label = 25
      CASE ("gxy3")
         label = 26
      CASE ("gy3z")
         label = 27
      CASE ("gxz3")
         label = 28
      CASE ("gyz3")
         label = 29
      CASE ("gx2y2")
         label = 30
      CASE ("gx2z2")
         label = 31
      CASE ("gy2z2")
         label = 32
      CASE ("gx2yz")
         label = 33
      CASE ("gxy2z")
         label = 34
      CASE ("gxyz2")
         label = 35
      CASE ("hz5")
         label = 36
      CASE ("hyz4")
         label = 37
      CASE ("hy2z3")
         label = 38
      CASE ("hy3z2")
         label = 39
      CASE ("hy4z")
         label = 40
      CASE ("hy5")
         label = 41
      CASE ("hxz4")
         label = 42
      CASE ("hxyz3")
         label = 43
      CASE ("hxy2z2")
         label = 44
      CASE ("hxy3z")
         label = 45
      CASE ("hxy4")
         label = 46
      CASE ("hx2z3")
         label = 47
      CASE ("hx2yz2")
         label = 48
      CASE ("hx2y2z")
         label = 49
      CASE ("hx2y3")
         label = 50
      CASE ("hx3z2")
         label = 51
      CASE ("hx3yz")
         label = 52
      CASE ("hx3y2")
         label = 53
      CASE ("hx4z")
         label = 54
      CASE ("hx4y")
         label = 55
      CASE ("hx5")
         label = 56
      CASE default
         CPABORT("The chargemol interface supports basis functions up to l=5 only")
         !label = 99
      END SELECT

   END FUNCTION pgf_type_chargemol

END MODULE qs_chargemol
