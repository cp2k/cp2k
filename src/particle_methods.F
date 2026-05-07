!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Define methods related to particle_type
!> \par History
!>            10.2014 Move routines out of particle_types.F [Ole Schuett]
!> \author Ole Schuett
! **************************************************************************************************
MODULE particle_methods
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type
   USE cell_methods,                    ONLY: cell_create,&
                                              set_cell_param
   USE cell_types,                      ONLY: cell_clone,&
                                              cell_release,&
                                              cell_type,&
                                              get_cell,&
                                              pbc,&
                                              real_to_scaled
   USE cp2k_info,                       ONLY: compile_revision,&
                                              cp2k_version,&
                                              r_cwd,&
                                              r_host_name,&
                                              r_user_name
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_unit_nr
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE external_potential_types,        ONLY: fist_potential_type,&
                                              get_potential
   USE input_constants,                 ONLY: dump_atomic,&
                                              dump_dcd,&
                                              dump_dcd_aligned_cell,&
                                              dump_extxyz,&
                                              dump_pdb,&
                                              dump_xmol
   USE input_cp2k_subsys,               ONLY: create_cell_section
   USE input_enumeration_types,         ONLY: enum_i2c,&
                                              enumeration_type
   USE input_keyword_types,             ONLY: keyword_get,&
                                              keyword_type
   USE input_section_types,             ONLY: section_get_keyword,&
                                              section_release,&
                                              section_type,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              sp
   USE machine,                         ONLY: m_timestamp,&
                                              timestamp_length
   USE mathconstants,                   ONLY: degree
   USE mathlib,                         ONLY: angle,&
                                              dihedral_angle,&
                                              gcd
   USE memory_utilities,                ONLY: reallocate
   USE particle_types,                  ONLY: get_particle_pos_or_vel,&
                                              particle_type
   USE periodic_table,                  ONLY: nelem
   USE physcon,                         ONLY: massunit
   USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE shell_potential_types,           ONLY: get_shell,&
                                              shell_kind_type
   USE string_utilities,                ONLY: uppercase
   USE util,                            ONLY: sort,&
                                              sort_unique
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Public subroutines

   PUBLIC :: write_fist_particle_coordinates, &
             write_qs_particle_coordinates, &
             write_particle_distances, &
             write_particle_coordinates, &
             write_structure_data, &
             get_particle_set, &
             write_particle_matrix, &
             write_final_structure

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'particle_methods'

CONTAINS

! **************************************************************************************************
!> \brief   Get the components of a particle set.
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param first_sgf ...
!> \param last_sgf ...
!> \param nsgf ...
!> \param nmao ...
!> \param basis ...
!> \date    14.01.2002
!> \par History
!>      - particle type cleaned (13.10.2003,MK)
!>      - refactoring and add basis set option (17.08.2010,jhu)
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE get_particle_set(particle_set, qs_kind_set, first_sgf, last_sgf, nsgf, &
                               nmao, basis)

      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL     :: first_sgf, last_sgf, nsgf, nmao
      TYPE(gto_basis_set_p_type), DIMENSION(:), OPTIONAL :: basis

      INTEGER                                            :: ikind, iparticle, isgf, nparticle, ns

      CPASSERT(ASSOCIATED(particle_set))

      nparticle = SIZE(particle_set)
      IF (PRESENT(first_sgf)) THEN
         CPASSERT(SIZE(first_sgf) >= nparticle)
      END IF
      IF (PRESENT(last_sgf)) THEN
         CPASSERT(SIZE(last_sgf) >= nparticle)
      END IF
      IF (PRESENT(nsgf)) THEN
         CPASSERT(SIZE(nsgf) >= nparticle)
      END IF
      IF (PRESENT(nmao)) THEN
         CPASSERT(SIZE(nmao) >= nparticle)
      END IF

      IF (PRESENT(first_sgf) .OR. PRESENT(last_sgf) .OR. PRESENT(nsgf)) THEN
         isgf = 0
         DO iparticle = 1, nparticle
            CALL get_atomic_kind(particle_set(iparticle)%atomic_kind, kind_number=ikind)
            IF (PRESENT(basis)) THEN
               IF (ASSOCIATED(basis(ikind)%gto_basis_set)) THEN
                  CALL get_gto_basis_set(gto_basis_set=basis(ikind)%gto_basis_set, nsgf=ns)
               ELSE
                  ns = 0
               END IF
            ELSE
               CALL get_qs_kind(qs_kind_set(ikind), nsgf=ns)
            END IF
            IF (PRESENT(nsgf)) nsgf(iparticle) = ns
            IF (PRESENT(first_sgf)) first_sgf(iparticle) = isgf + 1
            isgf = isgf + ns
            IF (PRESENT(last_sgf)) last_sgf(iparticle) = isgf
         END DO
      END IF

      IF (PRESENT(first_sgf)) THEN
         IF (SIZE(first_sgf) > nparticle) first_sgf(nparticle + 1) = isgf + 1
      END IF

      IF (PRESENT(nmao)) THEN
         DO iparticle = 1, nparticle
            CALL get_atomic_kind(particle_set(iparticle)%atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), mao=ns)
            nmao(iparticle) = ns
         END DO
      END IF

   END SUBROUTINE get_particle_set

! **************************************************************************************************
!> \brief   Should be able to write a few formats e.g. xmol, and some binary
!>          format (dcd) some format can be used for x,v,f
!>
!>          FORMAT   CONTENT                                    UNITS x, v, f
!>          XMOL     POS, VEL, FORCE, POS_VEL, POS_VEL_FORCE    Angstrom, a.u., a.u.
!>
!> \param particle_set ...
!> \param iunit ...
!> \param output_format ...
!> \param content ...
!> \param title ...
!> \param cell ...
!> \param array ...
!> \param unit_conv ...
!> \param charge_occup ...
!> \param charge_beta ...
!> \param charge_extended ...
!> \param print_kind ...
!> \date    14.01.2002
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE write_particle_coordinates(particle_set, iunit, output_format, &
                                         content, title, cell, array, unit_conv, &
                                         charge_occup, charge_beta, &
                                         charge_extended, print_kind)

      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER                                            :: iunit, output_format
      CHARACTER(LEN=*)                                   :: content, title
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: array
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: unit_conv
      LOGICAL, INTENT(IN), OPTIONAL                      :: charge_occup, charge_beta, &
                                                            charge_extended, print_kind

      CHARACTER(len=*), PARAMETER :: routineN = 'write_particle_coordinates'

      CHARACTER(LEN=120)                                 :: line
      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=4)                                   :: name
      CHARACTER(LEN=default_string_length)               :: atm_name, my_format
      INTEGER                                            :: handle, iatom, natom
      LOGICAL                                            :: dummy, my_charge_beta, &
                                                            my_charge_extended, my_charge_occup, &
                                                            my_print_kind
      REAL(KIND=dp)                                      :: angle_alpha, angle_beta, angle_gamma, &
                                                            factor, qeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: arr
      REAL(KIND=dp), DIMENSION(3)                        :: abc, angles, f, r, v
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h
      REAL(KIND=sp), ALLOCATABLE, DIMENSION(:)           :: x4, y4, z4
      TYPE(cell_type), POINTER                           :: cell_dcd
      TYPE(fist_potential_type), POINTER                 :: fist_potential
      TYPE(shell_kind_type), POINTER                     :: shell

      CALL timeset(routineN, handle)

      natom = SIZE(particle_set)
      IF (PRESENT(array)) THEN
         SELECT CASE (TRIM(content))
         CASE ("POS_VEL", "POS_VEL_FORCE")
            CPABORT("Illegal usage")
         END SELECT
      END IF
      factor = 1.0_dp
      IF (PRESENT(unit_conv)) THEN
         factor = unit_conv
      END IF
      SELECT CASE (output_format)
      CASE (dump_xmol, dump_extxyz)
         my_print_kind = .FALSE.
         IF (PRESENT(print_kind)) my_print_kind = print_kind
         WRITE (iunit, "(I8)") natom
         WRITE (iunit, "(A)") TRIM(title)
         DO iatom = 1, natom
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                 element_symbol=element_symbol)
            IF (LEN_TRIM(element_symbol) == 0 .OR. my_print_kind) THEN
               CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                    name=atm_name)
               dummy = qmmm_ff_precond_only_qm(id1=atm_name)
               my_format = "(A,"
               atm_name = TRIM(atm_name)
            ELSE
               my_format = "(T2,A2,"
               atm_name = TRIM(element_symbol)
            END IF
            SELECT CASE (TRIM(content))
            CASE ("POS")
               IF (PRESENT(array)) THEN
                  r(1:3) = get_particle_pos_or_vel(iatom, particle_set, array)
               ELSE
                  r(:) = particle_set(iatom)%r(:)
               END IF
               WRITE (iunit, TRIM(my_format)//"1X,3F20.10)") TRIM(atm_name), r(1:3)*factor
            CASE ("VEL")
               IF (PRESENT(array)) THEN
                  v(1:3) = get_particle_pos_or_vel(iatom, particle_set, array)
               ELSE
                  v(:) = particle_set(iatom)%v(:)
               END IF
               WRITE (iunit, TRIM(my_format)//"1X,3F20.10)") TRIM(atm_name), v(1:3)*factor
            CASE ("FORCE")
               IF (PRESENT(array)) THEN
                  f(:) = array((iatom - 1)*3 + 1:(iatom - 1)*3 + 3)
               ELSE
                  f(:) = particle_set(iatom)%f(:)
               END IF
               WRITE (iunit, TRIM(my_format)//"1X,3F20.10)") TRIM(atm_name), f(1:3)*factor
            CASE ("FORCE_MIXING_LABELS")
               IF (PRESENT(array)) THEN
                  f(:) = array((iatom - 1)*3 + 1:(iatom - 1)*3 + 3)
               ELSE
                  f(:) = particle_set(iatom)%f(:)
               END IF
               WRITE (iunit, TRIM(my_format)//"1X,3F20.10)") TRIM(atm_name), f(1:3)*factor
            END SELECT
         END DO
      CASE (dump_atomic)
         DO iatom = 1, natom
            SELECT CASE (TRIM(content))
            CASE ("POS")
               IF (PRESENT(array)) THEN
                  r(1:3) = get_particle_pos_or_vel(iatom, particle_set, array)
               ELSE
                  r(:) = particle_set(iatom)%r(:)
               END IF
               WRITE (iunit, "(3F20.10)") r(1:3)*factor
            CASE ("VEL")
               IF (PRESENT(array)) THEN
                  v(1:3) = get_particle_pos_or_vel(iatom, particle_set, array)
               ELSE
                  v(:) = particle_set(iatom)%v(:)
               END IF
               WRITE (iunit, "(3F20.10)") v(1:3)*factor
            CASE ("FORCE")
               IF (PRESENT(array)) THEN
                  f(:) = array((iatom - 1)*3 + 1:(iatom - 1)*3 + 3)
               ELSE
                  f(:) = particle_set(iatom)%f(:)
               END IF
               WRITE (iunit, "(3F20.10)") f(1:3)*factor
            CASE ("FORCE_MIXING_LABELS")
               IF (PRESENT(array)) THEN
                  f(:) = array((iatom - 1)*3 + 1:(iatom - 1)*3 + 3)
               ELSE
                  f(:) = particle_set(iatom)%f(:)
               END IF
               WRITE (iunit, "(3F20.10)") f(1:3)*factor
            END SELECT
         END DO
      CASE (dump_dcd, dump_dcd_aligned_cell)
         IF (.NOT. (PRESENT(cell))) THEN
            CPABORT("Cell is not present! Report this bug!")
         END IF
         CALL get_cell(cell, alpha=angle_alpha, beta=angle_beta, gamma=angle_gamma, &
                       abc=abc)
         IF (.NOT. cell%orthorhombic .AND. (output_format == dump_dcd_aligned_cell)) THEN
            ! In the case of a non-orthorhombic cell adopt a common convention
            ! for the orientation of the cell with respect to the Cartesian axes:
            ! Cell vector a is aligned with the x axis and the cell vector b lies
            ! in the xy plane.
            NULLIFY (cell_dcd)
            CALL cell_create(cell_dcd)
            CALL cell_clone(cell, cell_dcd, tag="CELL_DCD")
            angles(1) = angle_alpha/degree
            angles(2) = angle_beta/degree
            angles(3) = angle_gamma/degree
            CALL set_cell_param(cell_dcd, abc, angles, &
                                do_init_cell=.TRUE.)
            h(1:3, 1:3) = MATMUL(cell_dcd%hmat(1:3, 1:3), cell%h_inv(1:3, 1:3))
            CALL cell_release(cell_dcd)
         END IF
         ALLOCATE (arr(3, natom))
         IF (PRESENT(array)) THEN
            arr(1:3, 1:natom) = RESHAPE(array, [3, natom])
         ELSE
            SELECT CASE (TRIM(content))
            CASE ("POS")
               DO iatom = 1, natom
                  arr(1:3, iatom) = particle_set(iatom)%r(1:3)
               END DO
            CASE ("VEL")
               DO iatom = 1, natom
                  arr(1:3, iatom) = particle_set(iatom)%v(1:3)
               END DO
            CASE ("FORCE")
               DO iatom = 1, natom
                  arr(1:3, iatom) = particle_set(iatom)%f(1:3)
               END DO
            CASE DEFAULT
               CPABORT("Illegal DCD dump type")
            END SELECT
         END IF
         ALLOCATE (x4(natom))
         ALLOCATE (y4(natom))
         ALLOCATE (z4(natom))
         IF (.NOT. cell%orthorhombic .AND. (output_format == dump_dcd_aligned_cell)) THEN
            x4(1:natom) = REAL(MATMUL(h(1, 1:3), arr(1:3, 1:natom)), KIND=sp)
            y4(1:natom) = REAL(MATMUL(h(2, 1:3), arr(1:3, 1:natom)), KIND=sp)
            z4(1:natom) = REAL(MATMUL(h(3, 1:3), arr(1:3, 1:natom)), KIND=sp)
         ELSE
            x4(1:natom) = REAL(arr(1, 1:natom), KIND=sp)
            y4(1:natom) = REAL(arr(2, 1:natom), KIND=sp)
            z4(1:natom) = REAL(arr(3, 1:natom), KIND=sp)
         END IF
         WRITE (iunit) abc(1)*factor, angle_gamma, abc(2)*factor, &
            angle_beta, angle_alpha, abc(3)*factor
         WRITE (iunit) x4*REAL(factor, KIND=sp)
         WRITE (iunit) y4*REAL(factor, KIND=sp)
         WRITE (iunit) z4*REAL(factor, KIND=sp)
         ! Release work storage
         DEALLOCATE (arr)
         DEALLOCATE (x4)
         DEALLOCATE (y4)
         DEALLOCATE (z4)
      CASE (dump_pdb)
         my_charge_occup = .FALSE.
         IF (PRESENT(charge_occup)) my_charge_occup = charge_occup
         my_charge_beta = .FALSE.
         IF (PRESENT(charge_beta)) my_charge_beta = charge_beta
         my_charge_extended = .FALSE.
         IF (PRESENT(charge_extended)) my_charge_extended = charge_extended
         IF (LEN_TRIM(title) > 0) THEN
            WRITE (UNIT=iunit, FMT="(A6,T11,A)") &
               "REMARK", TRIM(title)
         END IF
         CALL get_cell(cell, alpha=angle_alpha, beta=angle_beta, gamma=angle_gamma, abc=abc)
         ! COLUMNS       DATA TYPE      CONTENTS
         ! --------------------------------------------------
         !  1 -  6       Record name    "CRYST1"
         !  7 - 15       Real(9.3)      a (Angstroms)
         ! 16 - 24       Real(9.3)      b (Angstroms)
         ! 25 - 33       Real(9.3)      c (Angstroms)
         ! 34 - 40       Real(7.2)      alpha (degrees)
         ! 41 - 47       Real(7.2)      beta (degrees)
         ! 48 - 54       Real(7.2)      gamma (degrees)
         ! 56 - 66       LString        Space group
         ! 67 - 70       Integer        Z value
         WRITE (UNIT=iunit, FMT="(A6,3F9.3,3F7.2)") &
            "CRYST1", abc(1:3)*factor, angle_alpha, angle_beta, angle_gamma
         WRITE (UNIT=line(1:6), FMT="(A6)") "ATOM  "
         DO iatom = 1, natom
            line = ""
            ! COLUMNS        DATA TYPE       CONTENTS
            !  1 -  6        Record name     "ATOM  "
            !  7 - 11        Integer         Atom serial number
            ! 13 - 16        Atom            Atom name
            ! 17             Character       Alternate location indicator
            ! 18 - 20        Residue name    Residue name
            ! 22             Character       Chain identifier
            ! 23 - 26        Integer         Residue sequence number
            ! 27             AChar           Code for insertion of residues
            ! 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstrom
            ! 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstrom
            ! 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstrom
            ! 55 - 60        Real(6.2)       Occupancy
            ! 61 - 66        Real(6.2)       Temperature factor (Default = 0.0)
            ! 73 - 76        LString(4)      Segment identifier, left-justified
            ! 77 - 78        LString(2)      Element symbol, right-justified
            ! 79 - 80        LString(2)      Charge on the atom
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                 element_symbol=element_symbol, name=atm_name, &
                                 fist_potential=fist_potential, shell=shell)
            IF (LEN_TRIM(element_symbol) == 0) THEN
               dummy = qmmm_ff_precond_only_qm(id1=atm_name)
            END IF
            name = TRIM(atm_name)
            IF (ASSOCIATED(fist_potential)) THEN
               CALL get_potential(potential=fist_potential, qeff=qeff)
            ELSE
               qeff = 0.0_dp
            END IF
            IF (ASSOCIATED(shell)) CALL get_shell(shell=shell, charge=qeff)
            WRITE (UNIT=line(1:6), FMT="(A6)") "ATOM  "
            WRITE (UNIT=line(7:11), FMT="(I5)") MODULO(iatom, 100000)
            WRITE (UNIT=line(13:16), FMT="(A4)") ADJUSTL(name)
            ! WRITE (UNIT=line(18:20),FMT="(A3)") TRIM(resname)
            ! WRITE (UNIT=line(23:26),FMT="(I4)") MODULO(idres,10000)
            SELECT CASE (TRIM(content))
            CASE ("POS")
               IF (PRESENT(array)) THEN
                  r(1:3) = get_particle_pos_or_vel(iatom, particle_set, array)
               ELSE
                  r(:) = particle_set(iatom)%r(:)
               END IF
               WRITE (UNIT=line(31:54), FMT="(3F8.3)") r(1:3)*factor
            CASE DEFAULT
               CPABORT("PDB dump only for trajectory available")
            END SELECT
            IF (my_charge_occup) THEN
               WRITE (UNIT=line(55:60), FMT="(F6.2)") qeff
            ELSE
               WRITE (UNIT=line(55:60), FMT="(F6.2)") 0.0_dp
            END IF
            IF (my_charge_beta) THEN
               WRITE (UNIT=line(61:66), FMT="(F6.2)") qeff
            ELSE
               WRITE (UNIT=line(61:66), FMT="(F6.2)") 0.0_dp
            END IF
            ! WRITE (UNIT=line(73:76),FMT="(A4)") ADJUSTL(TRIM(molname))
            WRITE (UNIT=line(77:78), FMT="(A2)") ADJUSTR(TRIM(element_symbol))
            IF (my_charge_extended) THEN
               WRITE (UNIT=line(81:), FMT="(SP,F0.8)") qeff
            END IF
            WRITE (UNIT=iunit, FMT="(A)") TRIM(line)
         END DO
         WRITE (UNIT=iunit, FMT="(A)") "END"
      CASE DEFAULT
         CPABORT("Illegal dump type")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE write_particle_coordinates

! **************************************************************************************************
!> \brief   Write the atomic coordinates to the output unit.
!> \param particle_set ...
!> \param subsys_section ...
!> \param charges ...
!> \date    05.06.2000
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE write_fist_particle_coordinates(particle_set, subsys_section, charges)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: subsys_section
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: charges

      CHARACTER(LEN=default_string_length)               :: name, unit_str
      INTEGER                                            :: iatom, ikind, iw, natom
      REAL(KIND=dp)                                      :: conv, mass, qcore, qeff, qshell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(shell_kind_type), POINTER                     :: shell_kind

      NULLIFY (logger)
      NULLIFY (shell_kind)

      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger, subsys_section, &
                                "PRINT%ATOMIC_COORDINATES", extension=".coordLog")

      CALL section_vals_val_get(subsys_section, "PRINT%ATOMIC_COORDINATES%UNIT", c_val=unit_str)
      conv = cp_unit_from_cp2k(1.0_dp, TRIM(unit_str))
      CALL uppercase(unit_str)
      IF (iw > 0) THEN
         WRITE (UNIT=iw, FMT="(/,/,T2,A)") &
            "MODULE FIST: ATOMIC COORDINATES IN "//TRIM(unit_str)
         WRITE (UNIT=iw, FMT="(/,T4,A,T30,A,T44,A,T58,A,T66,A,T77,A)") &
            "Atom Kind Name", "X", "Y", "Z", "q(eff)", "Mass"
         natom = SIZE(particle_set)
         DO iatom = 1, natom
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                 kind_number=ikind, &
                                 name=name, &
                                 mass=mass, &
                                 qeff=qeff, &
                                 shell=shell_kind)
            IF (PRESENT(charges)) qeff = charges(iatom)
            IF (ASSOCIATED(shell_kind)) THEN
               CALL get_shell(shell=shell_kind, &
                              charge_core=qcore, &
                              charge_shell=qshell)
               qeff = qcore + qshell
            END IF
            WRITE (UNIT=iw, FMT="(T2,I6,1X,I4,1X,A7,3(1X,F13.6),2(1X,F8.4))") &
               iatom, ikind, name, particle_set(iatom)%r(1:3)*conv, qeff, mass/massunit
         END DO
         WRITE (iw, "(A)") ""
      END IF

      CALL cp_print_key_finished_output(iw, logger, subsys_section, &
                                        "PRINT%ATOMIC_COORDINATES")

   END SUBROUTINE write_fist_particle_coordinates

! **************************************************************************************************
!> \brief   Write the atomic coordinates to the output unit.
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param subsys_section ...
!> \param label ...
!> \date    05.06.2000
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE write_qs_particle_coordinates(particle_set, qs_kind_set, subsys_section, label)

      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: subsys_section
      CHARACTER(LEN=*), INTENT(IN)                       :: label

      CHARACTER(len=*), PARAMETER :: routineN = 'write_qs_particle_coordinates'

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: unit_str
      INTEGER                                            :: handle, iatom, ikind, iw, natom, z
      REAL(KIND=dp)                                      :: conv, mass, zeff
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger, subsys_section, &
                                "PRINT%ATOMIC_COORDINATES", extension=".coordLog")

      CALL section_vals_val_get(subsys_section, "PRINT%ATOMIC_COORDINATES%UNIT", c_val=unit_str)
      conv = cp_unit_from_cp2k(1.0_dp, TRIM(unit_str))
      CALL uppercase(unit_str)
      IF (iw > 0) THEN
         WRITE (UNIT=iw, FMT="(/,/,T2,A)") &
            "MODULE "//TRIM(label)//": ATOMIC COORDINATES IN "//TRIM(unit_str)
         WRITE (UNIT=iw, FMT="(/,T4,A,T30,A,T44,A,T58,A,T66,A,T77,A)") &
            "Atom Kind Element", "X", "Y", "Z", "Z(eff)", "Mass"
         natom = SIZE(particle_set)
         DO iatom = 1, natom
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                 kind_number=ikind, &
                                 element_symbol=element_symbol, &
                                 mass=mass, &
                                 z=z)
            CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
            WRITE (UNIT=iw, FMT="(T2,I6,1X,I4,1X,A2,1X,I4,3(1X,F13.6),2(1X,F8.4))") &
               iatom, ikind, element_symbol, z, particle_set(iatom)%r(1:3)*conv, zeff, mass/massunit
         END DO
         WRITE (iw, "(A)") ""
      END IF

      CALL cp_print_key_finished_output(iw, logger, subsys_section, &
                                        "PRINT%ATOMIC_COORDINATES")

      CALL timestop(handle)

   END SUBROUTINE write_qs_particle_coordinates

! **************************************************************************************************
!> \brief   Write the matrix of the particle distances to the output unit.
!> \param particle_set ...
!> \param cell ...
!> \param subsys_section ...
!> \date    06.10.2000
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE write_particle_distances(particle_set, cell, subsys_section)

      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: subsys_section

      CHARACTER(len=*), PARAMETER :: routineN = 'write_particle_distances'

      CHARACTER(LEN=default_string_length)               :: unit_str
      INTEGER                                            :: handle, iatom, iw, jatom, natom
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: explicit
      REAL(KIND=dp)                                      :: conv, dab, dab_abort, dab_min, dab_warn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: distance_matrix
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(subsys_section))

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger, subsys_section, &
                                "PRINT%INTERATOMIC_DISTANCES", extension=".distLog")

      CALL section_vals_val_get(subsys_section, "PRINT%INTERATOMIC_DISTANCES%UNIT", c_val=unit_str)
      conv = cp_unit_from_cp2k(1.0_dp, TRIM(unit_str))
      CALL section_vals_val_get(subsys_section, "PRINT%INTERATOMIC_DISTANCES%CHECK_INTERATOMIC_DISTANCES", &
                                r_val=dab_min, explicit=explicit)

      dab_abort = 0.0_dp
      dab_warn = 0.0_dp
      natom = SIZE(particle_set)

      ! Compute interatomic distances only if their printout or check is explicitly requested
      ! Disable the default check for systems with more than 3000 atoms
      IF (explicit .OR. (iw > 0) .OR. (natom <= 2000)) THEN
         IF (dab_min > 0.0_dp) THEN
            dab_warn = dab_min*conv
         ELSE IF (dab_min < 0.0_dp) THEN
            dab_abort = ABS(dab_min)*conv
         END IF
      END IF

      IF ((iw > 0) .OR. (dab_abort > 0.0_dp) .OR. (dab_warn > 0.0_dp)) THEN
         CALL get_cell(cell=cell, periodic=periodic)
         IF (iw > 0) THEN
            ALLOCATE (distance_matrix(natom, natom))
            distance_matrix(:, :) = 0.0_dp
         END IF
         DO iatom = 1, natom
            DO jatom = iatom + 1, natom
               rab(:) = pbc(particle_set(iatom)%r(:), &
                            particle_set(jatom)%r(:), cell)
               dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))*conv
               IF (dab_abort > 0.0_dp) THEN
                  ! Stop the run for interatomic distances smaller than the requested threshold
                  IF (dab < dab_abort) THEN
                     CALL cp_abort(__LOCATION__, "The distance between the atoms "// &
                                   TRIM(ADJUSTL(cp_to_string(iatom, fmt="(I8)")))//" and "// &
                                   TRIM(ADJUSTL(cp_to_string(jatom, fmt="(I8)")))//" is only "// &
                                   TRIM(ADJUSTL(cp_to_string(dab, fmt="(F6.3)")))//" "// &
                                   TRIM(ADJUSTL(unit_str))//" and thus smaller than the requested threshold of "// &
                                   TRIM(ADJUSTL(cp_to_string(dab_abort, fmt="(F6.3)")))//" "// &
                                   TRIM(ADJUSTL(unit_str)))
                  END IF
               END IF
               IF (dab < dab_warn) THEN
                  ! Print warning for interatomic distances smaller than the requested threshold
                  CALL cp_warn(__LOCATION__, "The distance between the atoms "// &
                               TRIM(ADJUSTL(cp_to_string(iatom, fmt="(I8)")))//" and "// &
                               TRIM(ADJUSTL(cp_to_string(jatom, fmt="(I8)")))//" is only "// &
                               TRIM(ADJUSTL(cp_to_string(dab, fmt="(F6.3)")))//" "// &
                               TRIM(ADJUSTL(unit_str))//" and thus smaller than the threshold of "// &
                               TRIM(ADJUSTL(cp_to_string(dab_warn, fmt="(F6.3)")))//" "// &
                               TRIM(ADJUSTL(unit_str)))
               END IF
               IF (iw > 0) THEN
                  distance_matrix(iatom, jatom) = dab
                  distance_matrix(jatom, iatom) = distance_matrix(iatom, jatom)
               END IF
            END DO
         END DO
         IF (iw > 0) THEN
            ! Print the distance matrix
            WRITE (UNIT=iw, FMT="(/,/,T2,A)") &
               "INTERATOMIC DISTANCES IN "//TRIM(unit_str)
            CALL write_particle_matrix(distance_matrix, particle_set, iw)
            IF (ALLOCATED(distance_matrix)) DEALLOCATE (distance_matrix)
         END IF
         CALL cp_print_key_finished_output(iw, logger, subsys_section, &
                                           "PRINT%INTERATOMIC_DISTANCES")
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_particle_distances

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param particle_set ...
!> \param iw ...
!> \param el_per_part ...
!> \param Ilist ...
!> \param parts_per_line : number of particle columns to be printed in one line
! **************************************************************************************************
   SUBROUTINE write_particle_matrix(matrix, particle_set, iw, el_per_part, Ilist, parts_per_line)
      REAL(KIND=dp), DIMENSION(:, :)                     :: matrix
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: iw
      INTEGER, INTENT(IN), OPTIONAL                      :: el_per_part
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: Ilist
      INTEGER, INTENT(IN), OPTIONAL                      :: parts_per_line

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: fmt_string1, fmt_string2
      INTEGER                                            :: from, i, iatom, icol, jatom, katom, &
                                                            my_el_per_part, my_parts_per_line, &
                                                            natom, to
      INTEGER, DIMENSION(:), POINTER                     :: my_list

      my_el_per_part = 1
      IF (PRESENT(el_per_part)) my_el_per_part = el_per_part
      my_parts_per_line = 5
      IF (PRESENT(parts_per_line)) my_parts_per_line = MAX(parts_per_line, 1)
      WRITE (fmt_string1, FMT='(A,I0,A)') &
         "(/,T2,11X,", my_parts_per_line, "(4X,I5,4X))"
      WRITE (fmt_string2, FMT='(A,I0,A)') &
         "(T2,I5,2X,A2,2X,", my_parts_per_line, "(1X,F12.6))"
      IF (PRESENT(Ilist)) THEN
         natom = SIZE(Ilist)
      ELSE
         natom = SIZE(particle_set)
      END IF
      ALLOCATE (my_list(natom))
      IF (PRESENT(Ilist)) THEN
         my_list = Ilist
      ELSE
         DO i = 1, natom
            my_list(i) = i
         END DO
      END IF
      natom = natom*my_el_per_part
      DO jatom = 1, natom, my_parts_per_line
         from = jatom
         to = MIN(from + my_parts_per_line - 1, natom)
         WRITE (UNIT=iw, FMT=TRIM(fmt_string1)) (icol, icol=from, to)
         DO iatom = 1, natom
            katom = iatom/my_el_per_part
            IF (MOD(iatom, my_el_per_part) /= 0) katom = katom + 1
            CALL get_atomic_kind(atomic_kind=particle_set(my_list(katom))%atomic_kind, &
                                 element_symbol=element_symbol)
            WRITE (UNIT=iw, FMT=TRIM(fmt_string2)) &
               iatom, element_symbol, &
               (matrix(iatom, icol), icol=from, to)
         END DO
      END DO

      DEALLOCATE (my_list)

   END SUBROUTINE write_particle_matrix

! **************************************************************************************************
!> \brief   Write structure data requested by a separate structure data input
!>          section to the output unit.
!>          input_section can be either motion_section or subsys_section.
!>
!> \param particle_set ...
!> \param cell ...
!> \param input_section ...
!> \date    11.03.04
!> \par History
!>          Recovered (23.03.06,MK)
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE write_structure_data(particle_set, cell, input_section)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: input_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_structure_data'

      CHARACTER(LEN=default_string_length)               :: string, unit_str
      INTEGER                                            :: handle, i, i_rep, iw, n, n_rep, n_vals, &
                                                            natom, new_size, old_size, wrk2(2), &
                                                            wrk3(3), wrk4(4)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: work
      INTEGER, DIMENSION(:), POINTER                     :: atomic_indices, index_list
      LOGICAL                                            :: unique
      REAL(KIND=dp)                                      :: conv, dab
      REAL(KIND=dp), DIMENSION(3)                        :: r, rab, rbc, rcd, s
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: section

      CALL timeset(routineN, handle)
      NULLIFY (atomic_indices)
      NULLIFY (index_list)
      NULLIFY (logger)
      NULLIFY (section)
      string = ""

      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger=logger, &
                                basis_section=input_section, &
                                print_key_path="PRINT%STRUCTURE_DATA", &
                                extension=".coordLog")

      CALL section_vals_val_get(input_section, "PRINT%STRUCTURE_DATA%UNIT", c_val=unit_str)
      conv = cp_unit_from_cp2k(1.0_dp, TRIM(unit_str))
      CALL uppercase(unit_str)
      IF (iw > 0) THEN
         natom = SIZE(particle_set)
         section => section_vals_get_subs_vals(section_vals=input_section, &
                                               subsection_name="PRINT%STRUCTURE_DATA")

         WRITE (UNIT=iw, FMT="(/,T2,A)") "REQUESTED STRUCTURE DATA"
         ! Print the requested atomic position vectors
         CALL section_vals_val_get(section_vals=section, &
                                   keyword_name="POSITION", &
                                   n_rep_val=n_rep)
         IF (n_rep > 0) THEN
            WRITE (UNIT=iw, FMT="(/,T3,A,/)") &
               "Position vectors r(i) of the atoms i in "//TRIM(unit_str)
            old_size = 0
            DO i_rep = 1, n_rep
               CALL section_vals_val_get(section_vals=section, &
                                         keyword_name="POSITION", &
                                         i_rep_val=i_rep, &
                                         i_vals=atomic_indices)
               n_vals = SIZE(atomic_indices)
               new_size = old_size + n_vals
               CALL reallocate(index_list, 1, new_size)
               index_list(old_size + 1:new_size) = atomic_indices(1:n_vals)
               old_size = new_size
            END DO
            ALLOCATE (work(new_size))
            CALL sort(index_list, new_size, work)
            DEALLOCATE (work)
            DO i = 1, new_size
               WRITE (UNIT=string, FMT="(A,I0,A)") "(", index_list(i), ")"
               IF ((index_list(i) < 1) .OR. (index_list(i) > natom)) THEN
                  WRITE (UNIT=iw, FMT="(T3,A)") &
                     "Invalid atomic index "//TRIM(string)//" specified. Print request is ignored."
                  CYCLE
               END IF
               IF (i > 1) THEN
                  ! Skip redundant indices
                  IF (index_list(i) == index_list(i - 1)) CYCLE
               END IF
               WRITE (UNIT=iw, FMT="(T3,A,T20,A,3F13.6)") &
                  "r"//TRIM(string), "=", pbc(particle_set(index_list(i))%r(1:3), cell)*conv
            END DO
            DEALLOCATE (index_list)
         END IF

         ! Print the requested atomic position vectors in scaled coordinates
         CALL section_vals_val_get(section_vals=section, &
                                   keyword_name="POSITION_SCALED", &
                                   n_rep_val=n_rep)
         IF (n_rep > 0) THEN
            WRITE (UNIT=iw, FMT="(/,T3,A,/)") &
               "Position vectors s(i) of the atoms i in scaled coordinates"
            old_size = 0
            DO i_rep = 1, n_rep
               CALL section_vals_val_get(section_vals=section, &
                                         keyword_name="POSITION_SCALED", &
                                         i_rep_val=i_rep, &
                                         i_vals=atomic_indices)
               n_vals = SIZE(atomic_indices)
               new_size = old_size + n_vals
               CALL reallocate(index_list, 1, new_size)
               index_list(old_size + 1:new_size) = atomic_indices(1:n_vals)
               old_size = new_size
            END DO
            ALLOCATE (work(new_size))
            CALL sort(index_list, new_size, work)
            DEALLOCATE (work)
            DO i = 1, new_size
               WRITE (UNIT=string, FMT="(A,I0,A)") "(", index_list(i), ")"
               IF ((index_list(i) < 1) .OR. (index_list(i) > natom)) THEN
                  WRITE (UNIT=iw, FMT="(T3,A)") &
                     "Invalid atomic index "//TRIM(string)//" specified. Print request is ignored."
                  CYCLE
               END IF
               IF (i > 1) THEN
                  ! Skip redundant indices
                  IF (index_list(i) == index_list(i - 1)) CYCLE
               END IF
               r(1:3) = pbc(particle_set(index_list(i))%r(1:3), cell)
               CALL real_to_scaled(s, r, cell)
               WRITE (UNIT=iw, FMT="(T3,A,T20,A,3F13.6)") &
                  "s"//TRIM(string), "=", s(1:3)
            END DO
            DEALLOCATE (index_list)
         END IF

         ! Print the requested distances
         CALL section_vals_val_get(section_vals=section, &
                                   keyword_name="DISTANCE", &
                                   n_rep_val=n)
         IF (n > 0) THEN
            WRITE (UNIT=iw, FMT="(/,T3,A,/)") &
               "Distance vector r(i,j) between the atom i and j in "// &
               TRIM(unit_str)
            DO i = 1, n
               CALL section_vals_val_get(section_vals=section, &
                                         keyword_name="DISTANCE", &
                                         i_rep_val=i, &
                                         i_vals=atomic_indices)
               string = ""
               WRITE (UNIT=string, FMT="(A,2(I0,A))") &
                  "(", atomic_indices(1), ",", atomic_indices(2), ")"
               wrk2 = atomic_indices
               CALL sort_unique(wrk2, unique)
               IF (((wrk2(1) >= 1) .AND. (wrk2(SIZE(wrk2)) <= natom)) .AND. unique) THEN
                  rab(:) = pbc(particle_set(atomic_indices(1))%r(:), &
                               particle_set(atomic_indices(2))%r(:), cell)
                  dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
                  WRITE (UNIT=iw, FMT="(T3,A,T20,A,3F13.6,3X,A,F13.6)") &
                     "r"//TRIM(string), "=", rab(:)*conv, &
                     "|r| =", dab*conv
               ELSE
                  WRITE (UNIT=iw, FMT="(T3,A)") &
                     "Invalid atomic indices "//TRIM(string)//" specified. Print request is ignored."
               END IF
            END DO
         END IF

         ! Print the requested angles
         CALL section_vals_val_get(section_vals=section, &
                                   keyword_name="ANGLE", &
                                   n_rep_val=n)
         IF (n > 0) THEN
            WRITE (UNIT=iw, FMT="(/,T3,A,/)") &
               "Angle a(i,j,k) between the atomic distance vectors r(j,i) and "// &
               "r(j,k) in DEGREE"
            DO i = 1, n
               CALL section_vals_val_get(section_vals=section, &
                                         keyword_name="ANGLE", &
                                         i_rep_val=i, &
                                         i_vals=atomic_indices)
               string = ""
               WRITE (UNIT=string, FMT="(A,3(I0,A))") &
                  "(", atomic_indices(1), ",", atomic_indices(2), ",", atomic_indices(3), ")"
               wrk3 = atomic_indices
               CALL sort_unique(wrk3, unique)
               IF (((wrk3(1) >= 1) .AND. (wrk3(SIZE(wrk3)) <= natom)) .AND. unique) THEN
                  rab(:) = pbc(particle_set(atomic_indices(1))%r(:), &
                               particle_set(atomic_indices(2))%r(:), cell)
                  rbc(:) = pbc(particle_set(atomic_indices(2))%r(:), &
                               particle_set(atomic_indices(3))%r(:), cell)
                  WRITE (UNIT=iw, FMT="(T3,A,T26,A,F9.3)") &
                     "a"//TRIM(string), "=", angle(-rab, rbc)*degree
               ELSE
                  WRITE (UNIT=iw, FMT="(T3,A)") &
                     "Invalid atomic indices "//TRIM(string)//" specified. Print request is ignored."
               END IF
            END DO
         END IF

         ! Print the requested dihedral angles
         CALL section_vals_val_get(section_vals=section, &
                                   keyword_name="DIHEDRAL_ANGLE", &
                                   n_rep_val=n)
         IF (n > 0) THEN
            WRITE (UNIT=iw, FMT="(/,T3,A,/)") &
               "Dihedral angle d(i,j,k,l) between the planes (i,j,k) and (j,k,l) "// &
               "in DEGREE"
            DO i = 1, n
               CALL section_vals_val_get(section_vals=section, &
                                         keyword_name="DIHEDRAL_ANGLE", &
                                         i_rep_val=i, &
                                         i_vals=atomic_indices)
               string = ""
               WRITE (UNIT=string, FMT="(A,4(I0,A))") &
                  "(", atomic_indices(1), ",", atomic_indices(2), ",", &
                  atomic_indices(3), ",", atomic_indices(4), ")"
               wrk4 = atomic_indices
               CALL sort_unique(wrk4, unique)
               IF (((wrk4(1) >= 1) .AND. (wrk4(SIZE(wrk4)) <= natom)) .AND. unique) THEN
                  rab(:) = pbc(particle_set(atomic_indices(1))%r(:), &
                               particle_set(atomic_indices(2))%r(:), cell)
                  rbc(:) = pbc(particle_set(atomic_indices(2))%r(:), &
                               particle_set(atomic_indices(3))%r(:), cell)
                  rcd(:) = pbc(particle_set(atomic_indices(3))%r(:), &
                               particle_set(atomic_indices(4))%r(:), cell)
                  WRITE (UNIT=iw, FMT="(T3,A,T26,A,F9.3)") &
                     "d"//TRIM(string), "=", dihedral_angle(rab, rbc, rcd)*degree
               ELSE
                  WRITE (UNIT=iw, FMT="(T3,A)") &
                     "Invalid atomic indices "//TRIM(string)//" specified. Print request is ignored."
               END IF
            END DO
         END IF
      END IF
      CALL cp_print_key_finished_output(iw, logger, input_section, &
                                        "PRINT%STRUCTURE_DATA")

      CALL timestop(handle)

   END SUBROUTINE write_structure_data

! **************************************************************************************************
!> \brief   Write the final geometry and cell information to files
!> \param particle_set pointer to particles with atm_name, element_symbol and position
!> \param cell pointer to cell with abc, angle_alpha, angle_beta, angle_gamma and deth
!> \param input_section pointer to motion_section which has PRINT%FINAL_STRUCTURE
!> \param conv flag for whether convergence is achieved or not in optimization
!> \param keep_angles flag for whether cell optimization keeps initial angles
!> \param keep_symmetry flag for whether cell optimization keeps initial symmetry
!> \param keep_volume flag for whether cell optimization keeps initial volume
!> \param gopt_env_label the geometry optimization label "GEO_OPT", "CELL_OPT", ...
!> \param constraint_label label for directions with constraint in cell optimization
!> \par     Intended to be invoked in gopt_f_methods:write_final_info.
!>          This implementation does not consider higher space groups even if
!>          one is detected, and the chemical formulae are neither written in
!>          the sorted "Hill notation" nor expressed in groups of molecules.
!>          Other potentially useful but yet to be written information includes:
!>          the external pressure from CELL_OPT/EXTERNAL_POTENTIAL and the
!>          stress tensor (virial) for CELL_OPT;
!>          the fixed atoms from MOTION/CONSTRAINT/FIXED_ATOMS for all.
!>
!>          History
!>          04.2026 - Created
!> \author  HE Zilong
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE write_final_structure(particle_set, cell, input_section, conv, &
                                    keep_angles, keep_symmetry, keep_volume, &
                                    gopt_env_label, constraint_label)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      TYPE(section_vals_type), INTENT(IN), POINTER       :: input_section
      LOGICAL, INTENT(IN)                                :: conv, keep_angles, keep_symmetry, &
                                                            keep_volume
      CHARACTER(LEN=default_string_length), INTENT(IN)   :: gopt_env_label
      CHARACTER(LEN=4), INTENT(IN)                       :: constraint_label

      CHARACTER(len=*), PARAMETER :: routineN = 'write_final_structure'

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=2), ALLOCATABLE                      :: element_list(:)
      CHARACTER(LEN=5)                                   :: perd_str
      CHARACTER(LEN=:), ALLOCATABLE                      :: formula_structural, formula_sum
      CHARACTER(LEN=default_path_length)                 :: cell_str, record
      CHARACTER(LEN=default_string_length)               :: atm_name, f_cif, f_cif_label, &
                                                            f_cif_type_symbol, f_xyz
      CHARACTER(LEN=default_string_length), ALLOCATABLE  :: cif_label(:), cif_type_symbol(:), &
                                                            xyz_label(:)
      CHARACTER(LEN=timestamp_length)                    :: timestamp
      INTEGER :: elem_seen, file_unit, gcd_all, handle, i, iatom, ielem, natom, output_unit, &
         symmetry_id, w_cif_label, w_cif_type_symbol, w_xyz_label
      INTEGER, ALLOCATABLE                               :: count_list(:)
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: dummy, elem_in_list, orthorhombic, &
                                                            write_cif, write_xyz
      REAL(KIND=dp)                                      :: angle_alpha, angle_beta, angle_gamma, &
                                                            deth
      REAL(KIND=dp), DIMENSION(3)                        :: abc, r, s
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(enumeration_type), POINTER                    :: enum
      TYPE(keyword_type), POINTER                        :: symmetry_keyword
      TYPE(section_type), POINTER                        :: tmp_cell_section
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      NULLIFY (enum, logger, symmetry_keyword, print_key, tmp_cell_section)
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)
      print_key => section_vals_get_subs_vals(input_section, "PRINT%FINAL_STRUCTURE")

      ! Collect cell information
      perd_str = "F F F"
      CALL get_cell(cell, alpha=angle_alpha, beta=angle_beta, gamma=angle_gamma, &
                    deth=deth, orthorhombic=orthorhombic, abc=abc, periodic=periodic, &
                    h=hmat, symmetry_id=symmetry_id)
      IF (periodic(1) == 1) perd_str(1:1) = "T"
      IF (periodic(2) == 1) perd_str(3:3) = "T"
      IF (periodic(3) == 1) perd_str(5:5) = "T"
      CALL create_cell_section(tmp_cell_section)
      symmetry_keyword => section_get_keyword(tmp_cell_section, "SYMMETRY")
      CALL keyword_get(symmetry_keyword, enum=enum)
      ! cell_str is default_path_length which is longer
      ! than default_string_length and should be enough
      WRITE (UNIT=cell_str, FMT="(9(1X,F19.10))") &
         cp_unit_from_cp2k(hmat(1, 1), "angstrom"), &
         cp_unit_from_cp2k(hmat(2, 1), "angstrom"), &
         cp_unit_from_cp2k(hmat(3, 1), "angstrom"), &
         cp_unit_from_cp2k(hmat(1, 2), "angstrom"), &
         cp_unit_from_cp2k(hmat(2, 2), "angstrom"), &
         cp_unit_from_cp2k(hmat(3, 2), "angstrom"), &
         cp_unit_from_cp2k(hmat(1, 3), "angstrom"), &
         cp_unit_from_cp2k(hmat(2, 3), "angstrom"), &
         cp_unit_from_cp2k(hmat(3, 3), "angstrom")

      ! Collect atom information
      natom = SIZE(particle_set)
      ALLOCATE (element_list(nelem + 1), count_list(nelem + 1))
      count_list(:) = 0
      ALLOCATE (cif_label(natom), cif_type_symbol(natom), xyz_label(natom))
      elem_seen = 0
      w_cif_type_symbol = 0
      w_cif_label = 0
      w_xyz_label = 0
      atom_loop: DO iatom = 1, natom
         elem_in_list = .FALSE.
         CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                              name=atm_name, element_symbol=element_symbol)
         cif_type_symbol(iatom) = TRIM(atm_name)
         ! From write_particle_coordinates above it seems possible
         ! for some atoms to have empty element symbols; whatever
         ! these are, do not count them in the chemical formula
         IF (LEN_TRIM(element_symbol) == 0) THEN
            dummy = qmmm_ff_precond_only_qm(id1=atm_name)
            cif_label(iatom) = TRIM(atm_name)//TRIM(ADJUSTL(cp_to_string(iatom)))
            xyz_label(iatom) = TRIM(atm_name)
         ELSE
            cif_label(iatom) = TRIM(element_symbol)//TRIM(ADJUSTL(cp_to_string(iatom)))
            xyz_label(iatom) = TRIM(element_symbol)
            elem_loop: DO ielem = 1, elem_seen
               IF (element_list(ielem) == element_symbol) THEN
                  elem_in_list = .TRUE.
                  count_list(ielem) = count_list(ielem) + 1
                  EXIT elem_loop
               END IF
            END DO elem_loop
            IF (.NOT. elem_in_list) THEN
               elem_seen = elem_seen + 1
               element_list(elem_seen) = element_symbol
               count_list(elem_seen) = 1
            END IF
         END IF
         IF (LEN_TRIM(cif_type_symbol(iatom)) > w_cif_type_symbol) &
            w_cif_type_symbol = LEN_TRIM(cif_type_symbol(iatom))
         IF (LEN_TRIM(cif_label(iatom)) > w_cif_label) &
            w_cif_label = LEN_TRIM(cif_label(iatom))
         IF (LEN_TRIM(xyz_label(iatom)) > w_xyz_label) &
            w_xyz_label = LEN_TRIM(xyz_label(iatom))
      END DO atom_loop

      ! Determine the format of each line in cif considering width of cif_type_symbol and cif_label
      ! The fields are, in order:
      !  _atom_site_type_symbol, _atom_site_label, _atom_site_symmetry_multiplicity,
      !  _atom_site_fract_x, _atom_site_fract_y, _atom_site_fract_z, _atom_site_occupancy
      ! in which:
      !  _atom_site_type_symbol is taken as atm_name
      !  _atom_site_label is taken as element_symbol//iatom
      !  _atom_site_symmetry_multiplicity and _atom_site_occupancy are always 1
      f_cif_type_symbol = "A"//TRIM(ADJUSTL(cp_to_string(w_cif_type_symbol + 4)))
      f_cif_label = "A"//TRIM(ADJUSTL(cp_to_string(w_cif_label + 4)))
      f_cif = "(T3,"//TRIM(f_cif_type_symbol)//","//TRIM(f_cif_label)//",I4,3F14.8,F8.2)"

      ! Determine the format of each line in xyz
      f_xyz = "(T2,A"//TRIM(ADJUSTL(cp_to_string(w_xyz_label + 4)))//",1X,3F20.10)"

      ! Determine formula_sum
      CPASSERT(elem_seen > 0)
      CPASSERT(count_list(1) > 0)
      formula_sum = "'"
      DO ielem = 1, elem_seen
         formula_sum = formula_sum//TRIM(ADJUSTL(element_list(ielem)))
         formula_sum = formula_sum//TRIM(ADJUSTL(cp_to_string(count_list(ielem))))
         formula_sum = formula_sum//" "
      END DO
      formula_sum = TRIM(ADJUSTL(formula_sum))//"'"

      ! Determine formula_structural and Z
      gcd_all = count_list(1)
      DO ielem = 1, elem_seen
         IF (count_list(ielem) /= 0) THEN
            gcd_all = gcd(gcd_all, count_list(ielem))
         END IF
      END DO
      IF (gcd_all > 1) count_list = count_list/gcd_all
      formula_structural = "'"
      DO ielem = 1, elem_seen
         formula_structural = formula_structural//TRIM(ADJUSTL(element_list(ielem)))
         formula_structural = formula_structural//TRIM(ADJUSTL(cp_to_string(count_list(ielem))))
         formula_structural = formula_structural//" "
      END DO
      formula_structural = TRIM(ADJUSTL(formula_structural))//"'"

      ! Write XYZ
      CALL section_vals_val_get(print_key, "PRINT_XYZ", l_val=write_xyz)
      IF (write_xyz) THEN
         ! Print a message to log
         record = cp_print_key_generate_filename(logger, print_key, &
                                                 extension=".xyz", &
                                                 my_local=.FALSE.)
         IF (output_unit > 0) THEN
            IF (conv) THEN
               WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
                  routineN//": Optimization converged, writing XYZ file gladly:"
            ELSE
               WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
                  routineN//": Optimization not yet converged, writing XYZ file anyway:"
            END IF
            WRITE (UNIT=output_unit, FMT="(T3,A)") TRIM(record)
         END IF

         ! Make timestamp for the file
         CALL m_timestamp(timestamp)

         ! Prepare file unit and write to it
         file_unit = cp_print_key_unit_nr(logger, input_section, "PRINT%FINAL_STRUCTURE", &
                                          file_status="REPLACE", extension=".xyz")
         IF (file_unit > 0) THEN
            WRITE (UNIT=file_unit, FMT="(I8)") &
               natom
            WRITE (UNIT=file_unit, FMT="(A)") &
               'Lattice="'//TRIM(ADJUSTL(cell_str))// &
               '" Properties=species:S:1:pos:R:3 pbc="'// &
               TRIM(perd_str)//'"'
            DO iatom = 1, natom
               DO i = 1, 3
                  r(i) = cp_unit_from_cp2k(particle_set(iatom)%r(i), "angstrom")
               END DO
               WRITE (UNIT=file_unit, FMT=TRIM(f_xyz)) &
                  xyz_label(iatom), r(1:3)
            END DO
         END IF
      END IF

      ! Write CIF
      CALL section_vals_val_get(print_key, "PRINT_CIF", l_val=write_cif)
      IF (write_cif) THEN
         ! Print a message to log
         record = cp_print_key_generate_filename(logger, print_key, &
                                                 extension=".cif", &
                                                 my_local=.FALSE.)
         IF (output_unit > 0) THEN
            IF (conv) THEN
               WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
                  routineN//": Optimization converged, writing CIF file gladly:"
            ELSE
               WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
                  routineN//": Optimization not yet converged, writing CIF file anyway:"
            END IF
            WRITE (UNIT=output_unit, FMT="(T3,A)") TRIM(record)
         END IF

         ! Make timestamp for the file
         CALL m_timestamp(timestamp)

         ! Prepare file unit and write to it
         file_unit = cp_print_key_unit_nr(logger, input_section, "PRINT%FINAL_STRUCTURE", &
                                          file_status="REPLACE", extension=".cif")
         IF (file_unit > 0) THEN
            ! Generic information
            WRITE (UNIT=file_unit, FMT="(A)") &
               "# CIF file created by CP2K "//TRIM(moduleN)//":"//TRIM(routineN)
            WRITE (UNIT=file_unit, FMT="(A)") &
               "data_"//TRIM(logger%iter_info%project_name)
            WRITE (UNIT=file_unit, FMT="(A,T39,A)") &
               "_audit_creation_date", timestamp(:10)
            WRITE (UNIT=file_unit, FMT="(A,/,A,/,A)") &
               "_audit_creation_method", ";", &
               TRIM(cp2k_version)//" (revision "//TRIM(compile_revision)//")"
            WRITE (UNIT=file_unit, FMT="(A,/,A,/,A,/,A)") &
               "Project name "//TRIM(logger%iter_info%project_name), &
               "submitted by "//TRIM(r_user_name)//"@"//TRIM(r_host_name), &
               "processed in "//TRIM(r_cwd), &
               "generated at "//TRIM(timestamp)
            WRITE (UNIT=file_unit, FMT="(T2,A)") &
               "- Optimization type: "//TRIM(gopt_env_label)
            IF (conv) THEN
               WRITE (UNIT=file_unit, FMT="(T2,A)") &
                  "- Optimization converged: TRUE"
            ELSE
               WRITE (UNIT=file_unit, FMT="(T2,A)") &
                  "- Optimization converged: FALSE"
            END IF
            WRITE (UNIT=file_unit, FMT="(T2,A)") &
               "- Requested initial cell symmetry: "//TRIM(enum_i2c(enum, symmetry_id))
            IF (orthorhombic) THEN
               WRITE (UNIT=file_unit, FMT="(T2,A)") &
                  "- Cell is numerically orthorhombic: TRUE"
            ELSE
               WRITE (UNIT=file_unit, FMT="(T2,A)") &
                  "- Cell is numerically orthorhombic: FALSE"
            END IF
            WRITE (UNIT=file_unit, FMT="(T2,A)") &
               "- Periodicity of cell: "//TRIM(perd_str)
            IF (gopt_env_label == "CELL_OPT") THEN
               WRITE (UNIT=file_unit, FMT="(T2,A)") &
                  "- Cell is subject to optimization: TRUE"
               WRITE (UNIT=file_unit, FMT="(T2,A)") &
                  "- Cell has constraint on direction: "//TRIM(ADJUSTL(constraint_label))
               IF (keep_angles) THEN
                  WRITE (UNIT=file_unit, FMT="(T2,A)") &
                     "- Keep angles between the cell vectors during optimization: TRUE"
               ELSE
                  WRITE (UNIT=file_unit, FMT="(T2,A)") &
                     "- Keep angles between the cell vectors during optimization: FALSE"
               END IF
               IF (keep_symmetry) THEN
                  WRITE (UNIT=file_unit, FMT="(T2,A)") &
                     "- Keep initial cell symmetry during optimization: TRUE"
               ELSE
                  WRITE (UNIT=file_unit, FMT="(T2,A)") &
                     "- Keep initial cell symmetry during optimization: FALSE"
               END IF
               IF (keep_volume) THEN
                  WRITE (UNIT=file_unit, FMT="(T2,A)") &
                     "- Keep initial cell volume during optimization: TRUE"
               ELSE
                  WRITE (UNIT=file_unit, FMT="(T2,A)") &
                     "- Keep initial cell volume during optimization: FALSE"
               END IF
            ELSE
               WRITE (UNIT=file_unit, FMT="(T2,A)") &
                  "- Cell is subject to optimization: FALSE"
            END IF
            WRITE (UNIT=file_unit, FMT="(T2,A)") &
               "- Final cell vectors A, B, C by rows [angstrom]:"
            DO i = 1, 3
               WRITE (UNIT=file_unit, FMT="(T3,3(1X,F19.10))") &
                  cp_unit_from_cp2k(hmat(1, i), "angstrom"), &
                  cp_unit_from_cp2k(hmat(2, i), "angstrom"), &
                  cp_unit_from_cp2k(hmat(3, i), "angstrom")
            END DO
            WRITE (UNIT=file_unit, FMT="(A)") ";"
            ! Data of cell and geometry
            WRITE (UNIT=file_unit, FMT="(/,A,T44,A)") &
               "_symmetry_space_group_name_H-M", "'P 1'"
            WRITE (UNIT=file_unit, FMT="(A,T31,F18.8)") &
               "_cell_length_a", cp_unit_from_cp2k(abc(1), "angstrom")
            WRITE (UNIT=file_unit, FMT="(A,T31,F18.8)") &
               "_cell_length_b", cp_unit_from_cp2k(abc(2), "angstrom")
            WRITE (UNIT=file_unit, FMT="(A,T31,F18.8)") &
               "_cell_length_c", cp_unit_from_cp2k(abc(3), "angstrom")
            WRITE (UNIT=file_unit, FMT="(A,T31,F18.8)") &
               "_cell_angle_alpha", angle_alpha
            WRITE (UNIT=file_unit, FMT="(A,T31,F18.8)") &
               "_cell_angle_beta", angle_beta
            WRITE (UNIT=file_unit, FMT="(A,T31,F18.8)") &
               "_cell_angle_gamma", angle_gamma
            WRITE (UNIT=file_unit, FMT="(A,T48,A)") &
               "_symmetry_Int_Tables_number", "1"
            WRITE (UNIT=file_unit, FMT="(A,T36,A)") &
               "_chemical_formula_structural", formula_structural
            WRITE (UNIT=file_unit, FMT="(A,T36,A)") &
               "_chemical_formula_sum", formula_sum
            WRITE (UNIT=file_unit, FMT="(A,T31,F18.8)") &
               "_cell_volume", cp_unit_from_cp2k(ABS(deth), "angstrom^3")
            WRITE (UNIT=file_unit, FMT="(A,T41,I8)") &
               "_cell_formula_units_Z", gcd_all
            WRITE (UNIT=file_unit, FMT="(A,/,T2,A,/,T2,A,/,T3,A)") &
               "loop_", "_symmetry_equiv_pos_site_id", &
               "_symmetry_equiv_pos_as_xyz", "1  'x, y, z'"
            WRITE (UNIT=file_unit, FMT="(A,/,T2,A,/,T2,A,/,T2,A,/,T2,A,/,T2,A,/,T2,A,/,T2,A)") &
               "loop_", "_atom_site_type_symbol", "_atom_site_label", &
               "_atom_site_symmetry_multiplicity", "_atom_site_fract_x", &
               "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy"
            DO iatom = 1, natom
               r(1:3) = pbc(particle_set(iatom)%r(1:3), cell)
               CALL real_to_scaled(s, r, cell)
               DO i = 1, 3
                  s(i) = MODULO(s(i), 1.0_dp)
               END DO
               WRITE (UNIT=file_unit, FMT=TRIM(f_cif)) &
                  cif_type_symbol(iatom), cif_label(iatom), 1, s(1:3), 1.0_dp
            END DO
         END IF
      END IF

      ! Finish
      DEALLOCATE (element_list, count_list, formula_structural, &
                  formula_sum, cif_label, cif_type_symbol, &
                  xyz_label)
      CALL section_release(tmp_cell_section)
      CALL cp_print_key_finished_output(file_unit, logger, input_section, &
                                        "PRINT%FINAL_STRUCTURE")
      IF (output_unit > 0) &
         WRITE (UNIT=output_unit, FMT='(/,T2,A)') &
         routineN//": Done!"

      CALL timestop(handle)

   END SUBROUTINE write_final_structure

END MODULE particle_methods
