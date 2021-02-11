!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Space Group Symmetry Module  (version 1.0, January 16, 2020)
!> \par History
!>      Pierre-André Cazade [pcazade] 01.2020 - University of Limerick
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
MODULE space_groups

   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE bibliography,                    ONLY: Togo2018,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_copy,&
                                              cell_create,&
                                              cell_release,&
                                              cell_type,&
                                              init_cell,&
                                              real_to_scaled,&
                                              scaled_to_real
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE gopt_f_types,                    ONLY: gopt_f_type
   USE input_constants,                 ONLY: default_cell_direct_id,&
                                              default_cell_geo_opt_id,&
                                              default_cell_md_id,&
                                              default_cell_method_id,&
                                              default_minimization_method_id,&
                                              default_ts_method_id
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE machine,                         ONLY: default_output_unit
   USE mathlib,                         ONLY: inv_3x3
   USE particle_list_types,             ONLY: particle_list_type
   USE physcon,                         ONLY: pascal
   USE spglib_f08,                      ONLY: spg_get_international,&
                                              spg_get_multiplicity,&
                                              spg_get_pointgroup,&
                                              spg_get_schoenflies,&
                                              spg_get_symmetry
   USE string_utilities,                ONLY: strlcpy_c2f
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'space_groups'

   TYPE spgr_type
      LOGICAL                                     :: symlib = .FALSE.
      INTEGER                                     :: plevel = -1
      INTEGER                                     :: punit = default_output_unit
      INTEGER                                     :: istriz = -1
      REAL(KIND=dp)                               :: eps_symmetry = 1.0e-4_dp
      INTEGER                                     :: nparticle = 0
      INTEGER                                     :: nparticle_sym = 0
      INTEGER                                     :: n_atom = 0
      INTEGER                                     :: n_core = 0
      INTEGER                                     :: n_shell = 0
      INTEGER                                     :: n_atom_sym = 0
      INTEGER                                     :: n_core_sym = 0
      INTEGER                                     :: n_shell_sym = 0
      INTEGER                                     :: ndigits = 16
      INTEGER, DIMENSION(:), ALLOCATABLE          :: atype
      INTEGER, DIMENSION(:, :), ALLOCATABLE        :: eqatom
      LOGICAL, DIMENSION(:), ALLOCATABLE          :: lop, lat
      REAL(KIND=dp), DIMENSION(3)                 :: pol = 0.0_dp
      !SPGLIB
      CHARACTER(len=11)                           :: international_symbol = ""
      CHARACTER(len=6)                            :: pointgroup_symbol = ""
      CHARACTER(len=7)                            :: schoenflies = ""
      CHARACTER(LEN=default_path_length)          :: filename = ""
      INTEGER                                     :: n_operations = 0
      INTEGER                                     :: n_reduced_operations = 0
      INTEGER                                     :: n_operations_subset = 0
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE    :: rotations
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE    :: rotations_subset
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: translations
      TYPE(cell_type), POINTER                     :: cell_ref
   END TYPE spgr_type

   PUBLIC :: release_spgr_type, spgr_create, identify_space_group, spgr_find_equivalent_atoms
   PUBLIC :: spgr_apply_rotations_coord, spgr_apply_rotations_force, print_spgr, spgr_type
   PUBLIC :: spgr_apply_rotations_stress, spgr_write_stress_tensor

CONTAINS

! **************************************************************************************************
!> \brief Release the SPGR type
!> \param spgr  The SPGR type
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE release_spgr_type(spgr)

      TYPE(spgr_type)                                    :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'release_spgr_type', &
         routineP = moduleN//':'//routineN

      IF (ALLOCATED(spgr%rotations)) THEN
         DEALLOCATE (spgr%rotations)
      END IF
      IF (ALLOCATED(spgr%rotations_subset)) THEN
         DEALLOCATE (spgr%rotations_subset)
      END IF
      IF (ALLOCATED(spgr%translations)) THEN
         DEALLOCATE (spgr%translations)
      END IF
      IF (ALLOCATED(spgr%atype)) THEN
         DEALLOCATE (spgr%atype)
      END IF
      IF (ALLOCATED(spgr%eqatom)) THEN
         DEALLOCATE (spgr%eqatom)
      END IF
      IF (ALLOCATED(spgr%lop)) THEN
         DEALLOCATE (spgr%lop)
      END IF
      IF (ALLOCATED(spgr%lat)) THEN
         DEALLOCATE (spgr%lat)
      END IF

      CALL cell_release(spgr%cell_ref)

   END SUBROUTINE release_spgr_type

! **************************************************************************************************
!> \brief routine creates the space group structure
!> \param spgr ...
!> \param scoor ...
!> \param types ...
!> \param cell ...
!> \param gopt_env ...
!> \param eps_symmetry ...
!> \param pol ...
!> \param ranges ...
!> \param nparticle ...
!> \param n_atom ...
!> \param n_core ...
!> \param n_shell ...
!> \param iounit ...
!> \param plevel ...
!> \param filename ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_create(spgr, scoor, types, cell, gopt_env, eps_symmetry, pol, ranges, &
                          nparticle, n_atom, n_core, n_shell, iounit, plevel, filename)

      TYPE(spgr_type), INTENT(INOUT)                     :: spgr
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: scoor
      INTEGER, DIMENSION(:), INTENT(IN)                  :: types
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      TYPE(gopt_f_type), INTENT(IN), POINTER             :: gopt_env
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_symmetry
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: pol
      INTEGER, DIMENSION(:, :), INTENT(IN), OPTIONAL     :: ranges
      INTEGER, INTENT(IN), OPTIONAL                      :: nparticle, n_atom, n_core, n_shell, &
                                                            iounit, plevel
      CHARACTER(LEN=default_path_length), INTENT(IN), &
         OPTIONAL                                        :: filename

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_create', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ierr, j, n_sr_rep, nop, &
                                                            nchars, tra_mat(3, 3)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: tmp_types
      LOGICAL                                            :: spglib
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: tmp_coor
      CHARACTER(LEN=1000)                                :: buffer

      CALL timeset(routineN, handle)

      !..total number of particles (atoms plus shells)
      IF (PRESENT(nparticle)) THEN
         CPASSERT(nparticle == SIZE(scoor, 2))
         spgr%nparticle = nparticle
      ELSE
         spgr%nparticle = SIZE(scoor, 2)
      ENDIF

      IF (PRESENT(n_atom)) THEN
         spgr%n_atom = n_atom
      ELSE IF (PRESENT(n_core)) THEN
         spgr%n_atom = spgr%nparticle - n_core
      ELSE IF (PRESENT(n_shell)) THEN
         spgr%n_atom = spgr%nparticle - n_shell
      ELSE
         spgr%n_atom = spgr%nparticle
      ENDIF

      IF (PRESENT(n_core)) THEN
         spgr%n_core = n_core
      ELSE IF (PRESENT(n_shell)) THEN
         spgr%n_core = n_shell
      ENDIF

      IF (PRESENT(n_shell)) THEN
         spgr%n_shell = n_shell
      ELSE IF (PRESENT(n_core)) THEN
         spgr%n_shell = n_core
      ENDIF

      IF (.NOT. (spgr%nparticle == (spgr%n_atom + spgr%n_shell))) THEN
         CPABORT("spgr_create: nparticle not equal to natom + nshell.")
      END IF

      spgr%nparticle_sym = spgr%nparticle
      spgr%n_atom_sym = spgr%n_atom
      spgr%n_core_sym = spgr%n_core
      spgr%n_shell_sym = spgr%n_shell

      ! output unit
      IF (PRESENT(filename)) THEN
         spgr%filename = TRIM(filename)
      ELSE IF (PRESENT(iounit)) THEN
         spgr%punit = iounit
      ELSE
         spgr%punit = default_output_unit
      END IF

      ! output unit
      IF (PRESENT(plevel)) THEN
         spgr%plevel = plevel
      END IF

      ! accuracy for symmetry
      IF (PRESENT(eps_symmetry)) THEN
         spgr%eps_symmetry = eps_symmetry
      END IF

      ! vector to test reduced symmetry
      IF (PRESENT(pol)) THEN
         spgr%pol(1) = pol(1)
         spgr%pol(2) = pol(2)
         spgr%pol(3) = pol(3)
      END IF

      ALLOCATE (spgr%lat(spgr%nparticle))
      spgr%lat = .TRUE.

      IF (PRESENT(ranges)) THEN
         n_sr_rep = SIZE(ranges, 2)
         DO i = 1, n_sr_rep
            DO j = ranges(1, i), ranges(2, i)
               spgr%lat(j) = .FALSE.
               spgr%nparticle_sym = spgr%nparticle_sym - 1
               IF (j <= spgr%n_atom) THEN
                  spgr%n_atom_sym = spgr%n_atom_sym - 1
               ELSE IF (j > spgr%n_atom .AND. j <= spgr%nparticle) THEN
                  spgr%n_core_sym = spgr%n_core_sym - 1
                  spgr%n_shell_sym = spgr%n_shell_sym - 1
               ELSE
                  CPABORT("Symmetry exclusion range larger than actual number of particles.")
               ENDIF
            END DO
         END DO
      END IF

      ALLOCATE (tmp_coor(3, spgr%n_atom_sym), tmp_types(spgr%n_atom_sym))

      j = 0
      DO i = 1, spgr%n_atom
         IF (spgr%lat(i)) THEN
            j = j + 1
            tmp_coor(:, j) = scoor(:, i)
            tmp_types(j) = types(i)
         END IF
      END DO

      !..set cell values
      NULLIFY (spgr%cell_ref)
      CALL cell_create(spgr%cell_ref)
      CALL cell_copy(cell, spgr%cell_ref)
      SELECT CASE (gopt_env%type_id)
      CASE (default_minimization_method_id, default_ts_method_id)
         CALL init_cell(spgr%cell_ref, hmat=cell%hmat)
      CASE (default_cell_method_id)
         SELECT CASE (gopt_env%cell_method_id)
         CASE (default_cell_direct_id)
            CALL init_cell(spgr%cell_ref, hmat=gopt_env%h_ref)
         CASE (default_cell_geo_opt_id, default_cell_md_id)
            CPABORT("SPACE_GROUP_SYMMETRY SHOULD NOT BE INVOKED DURING THE CELL STEP.")
         CASE DEFAULT
            CPABORT("SPACE_GROUP_SYMMETRY INVOKED WITH AN UNKNOWN OPTIMIZATION METHOD.")
         END SELECT
      CASE DEFAULT
         CPABORT("SPACE_GROUP_SYMMETRY IS NOT COMPATIBLE WITH MD.")
      END SELECT

      ! atom types
      ALLOCATE (spgr%atype(spgr%nparticle))
      spgr%atype(1:spgr%nparticle) = types(1:spgr%nparticle)

      spgr%n_operations = 0

#ifdef __SPGLIB
      spglib = .TRUE.
      CALL cite_reference(Togo2018)
      ierr = spg_get_international(spgr%international_symbol, TRANSPOSE(cell%hmat), tmp_coor, tmp_types, &
                                   spgr%n_atom_sym, eps_symmetry)
      buffer = ''
      nchars = strlcpy_c2f(buffer, spgr%international_symbol)
      spgr%international_symbol = buffer(1:nchars)
      IF (ierr == 0) THEN
         CPABORT("Symmetry Library SPGLIB failed, most likely due a problem with the coordinates.")
         spglib = .FALSE.
      ELSE
         nop = spg_get_multiplicity(TRANSPOSE(cell%hmat), tmp_coor, tmp_types, &
                                    spgr%n_atom_sym, eps_symmetry)
         ALLOCATE (spgr%rotations(3, 3, nop), spgr%translations(3, nop))
         ALLOCATE (spgr%eqatom(nop, spgr%nparticle))
         ALLOCATE (spgr%lop(nop))
         spgr%n_operations = nop
         spgr%lop = .TRUE.
         ierr = spg_get_symmetry(spgr%rotations, spgr%translations, nop, TRANSPOSE(cell%hmat), &
                                 tmp_coor, tmp_types, spgr%n_atom_sym, eps_symmetry)
         ! Schoenflies Symbol
         ierr = spg_get_schoenflies(spgr%schoenflies, TRANSPOSE(cell%hmat), tmp_coor, tmp_types, &
                                    spgr%n_atom_sym, eps_symmetry)
         buffer = ''
         nchars = strlcpy_c2f(buffer, spgr%schoenflies)
         spgr%schoenflies = buffer(1:nchars)
         ! Point Group
         tra_mat = 0
         ierr = spg_get_pointgroup(spgr%pointgroup_symbol, tra_mat, &
                                   spgr%rotations, spgr%n_operations)
         buffer = ''
         nchars = strlcpy_c2f(buffer, spgr%pointgroup_symbol)
         spgr%pointgroup_symbol = buffer(1:nchars)
      END IF
#else
      CPABORT("Symmetry library SPGLIB not available")
      spglib = .FALSE.
#endif
      spgr%symlib = spglib

      DEALLOCATE (tmp_coor, tmp_types)

      CALL timestop(handle)

   END SUBROUTINE spgr_create

! **************************************************************************************************
!> \brief routine indentifies the space group and finds rotation matrices.
!> \param spgr ...
!> \param subsys ...
!> \param cell ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param iounit ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
!> \note  rotation matrices innclude translations and translation symmetry:
!>        it works with supercells as well.
! **************************************************************************************************
   SUBROUTINE identify_space_group(spgr, subsys, cell, geo_section, gopt_env, iounit)

      TYPE(spgr_type), INTENT(INOUT)                     :: spgr
      TYPE(cp_subsys_type), INTENT(IN), POINTER          :: subsys
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      TYPE(section_vals_type), INTENT(IN), POINTER       :: geo_section
      TYPE(gopt_f_type), INTENT(IN), POINTER             :: gopt_env
      INTEGER                                            :: iounit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'identify_space_group', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: handle, i, k, n_atom, n_core, n_shell, &
                                                            n_sr_rep, nparticle, plevel, &
                                                            shell_index
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atype
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ranges
      INTEGER, DIMENSION(:), POINTER                     :: tmp
      REAL(KIND=dp)                                      :: eps_symmetry
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: scoord
      REAL(KIND=dp), DIMENSION(:), POINTER               :: pol
      TYPE(particle_list_type), POINTER                  :: core_particles, particles, &
                                                            shell_particles

      CALL timeset(routineN, handle)

      n_sr_rep = 0
      nparticle = 0
      n_atom = 0
      n_core = 0
      n_shell = 0

      NULLIFY (particles)
      NULLIFY (core_particles)
      NULLIFY (shell_particles)

      CALL cp_subsys_get(subsys, particles=particles, shell_particles=shell_particles, core_particles=core_particles)

      CPASSERT(ASSOCIATED(particles))
      n_atom = particles%n_els
      ! Check if we have other kinds of particles in this subsystem
      IF (ASSOCIATED(shell_particles)) THEN
         n_shell = shell_particles%n_els
         CPASSERT(ASSOCIATED(core_particles))
         n_core = subsys%core_particles%n_els
         ! The same number of shell and core particles is assumed
         CPASSERT(n_core == n_shell)
      ELSE IF (ASSOCIATED(core_particles)) THEN
         ! This case should not occur at the moment
      ELSE
         n_core = 0
         n_shell = 0
      END IF

      nparticle = n_atom + n_shell
      ALLOCATE (scoord(3, nparticle), atype(nparticle))
      DO i = 1, n_atom
         shell_index = particles%els(i)%shell_index
         IF (shell_index == 0) THEN
            CALL real_to_scaled(scoord(1:3, i), particles%els(i)%r(1:3), cell)
            CALL get_atomic_kind(atomic_kind=particles%els(i)%atomic_kind, kind_number=atype(i))
         ELSE
            CALL real_to_scaled(scoord(1:3, i), core_particles%els(shell_index)%r(1:3), cell)
            CALL get_atomic_kind(atomic_kind=core_particles%els(shell_index)%atomic_kind, kind_number=atype(i))
            k = n_atom + shell_index
            CALL real_to_scaled(scoord(1:3, k), shell_particles%els(shell_index)%r(1:3), cell)
            CALL get_atomic_kind(atomic_kind=shell_particles%els(shell_index)%atomic_kind, kind_number=atype(k))
         END IF
      END DO

      CALL section_vals_val_get(geo_section, "SPGR_PRINT_LEVEL", i_val=plevel)
      CALL section_vals_val_get(geo_section, "SPGR_FILE_NAME", c_val=filename)
      CALL section_vals_val_get(geo_section, "EPS_SYMMETRY", r_val=eps_symmetry)
      CALL section_vals_val_get(geo_section, "SYMM_REDUCTION", r_vals=pol)
      CALL section_vals_val_get(geo_section, "SYMM_EXCLUDE_RANGE", n_rep_val=n_sr_rep)
      IF (n_sr_rep > 0) THEN
         ALLOCATE (ranges(2, n_sr_rep))
         DO i = 1, n_sr_rep
            CALL section_vals_val_get(geo_section, "SYMM_EXCLUDE_RANGE", i_rep_val=i, i_vals=tmp)
            ranges(:, i) = tmp(:)
         END DO
         CALL spgr_create(spgr, scoord, atype, cell, gopt_env, eps_symmetry=eps_symmetry, pol=pol(1:3), &
                          ranges=ranges, nparticle=nparticle, n_atom=n_atom, n_core=n_core, n_shell=n_shell, iounit=iounit, &
                          filename=filename, plevel=plevel)
         DEALLOCATE (ranges)
      ELSE
         CALL spgr_create(spgr, scoord, atype, cell, gopt_env, eps_symmetry=eps_symmetry, pol=pol(1:3), &
                          nparticle=nparticle, n_atom=n_atom, n_core=n_core, n_shell=n_shell, iounit=iounit, &
                          filename=filename, plevel=plevel)
      ENDIF
      CALL spgr_find_equivalent_atoms(spgr, scoord)
      CALL spgr_reduce_symm(spgr)
      CALL spgr_rotations_subset(spgr)
      DEALLOCATE (scoord, atype)

      CALL timestop(handle)

   END SUBROUTINE identify_space_group

! **************************************************************************************************
!> \brief routine indentifies the equivalent atoms for each rotation matrix.
!> \param spgr ...
!> \param scoord ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_find_equivalent_atoms(spgr, scoord)

      TYPE(spgr_type), INTENT(INOUT)                     :: spgr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: scoord

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_find_equivalent_atoms', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ia, ib, ir, j, natom, nop, &
                                                            nshell
      REAL(KIND=dp)                                      :: diff
      REAL(KIND=dp), DIMENSION(3)                        :: rb, ri, ro, tr
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot

      CALL timeset(routineN, handle)

      nop = spgr%n_operations
      natom = spgr%n_atom
      nshell = spgr%n_shell

      IF (.NOT. (spgr%nparticle == (natom + nshell))) THEN
         CPABORT("spgr_find_equivalent_atoms: nparticle not equal to natom + nshell.")
      END IF

      DO ia = 1, spgr%nparticle
         spgr%eqatom(:, ia) = ia
      END DO

      !$OMP PARALLEL DO PRIVATE (ia,ib,ir,ri,rb,ro,rot,tr,diff) SHARED (spgr,scoord,natom,nop) DEFAULT(NONE)
      DO ia = 1, natom
         IF (.NOT. spgr%lat(ia)) CYCLE
         ri(1:3) = scoord(1:3, ia)
         DO ir = 1, nop
            rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
            tr(1:3) = spgr%translations(1:3, ir)
            DO ib = 1, natom
               IF (.NOT. spgr%lat(ib)) CYCLE
               rb(1:3) = scoord(1:3, ib)
               ro(1) = REAL(rot(1, 1), dp)*rb(1) + REAL(rot(2, 1), dp)*rb(2) + REAL(rot(3, 1), dp)*rb(3) + tr(1)
               ro(2) = REAL(rot(1, 2), dp)*rb(1) + REAL(rot(2, 2), dp)*rb(2) + REAL(rot(3, 2), dp)*rb(3) + tr(2)
               ro(3) = REAL(rot(1, 3), dp)*rb(1) + REAL(rot(2, 3), dp)*rb(2) + REAL(rot(3, 3), dp)*rb(3) + tr(3)
               ro(1) = ro(1) - REAL(NINT(ro(1) - ri(1)), dp)
               ro(2) = ro(2) - REAL(NINT(ro(2) - ri(2)), dp)
               ro(3) = ro(3) - REAL(NINT(ro(3) - ri(3)), dp)
               diff = NORM2(ri(:) - ro(:))
               IF ((diff < spgr%eps_symmetry) .AND. (spgr%atype(ia) == spgr%atype(ib))) THEN
                  spgr%eqatom(ir, ia) = ib
                  EXIT
               END IF
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE (i,j,ia,ib,ir,ri,rb,ro,rot,tr,diff) SHARED (spgr,scoord,natom,nshell,nop) DEFAULT(NONE)
      DO i = 1, nshell
         ia = natom + i
         IF (.NOT. spgr%lat(ia)) CYCLE
         ri(1:3) = scoord(1:3, ia)
         DO ir = 1, nop
            rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
            tr(1:3) = spgr%translations(1:3, ir)
            DO j = 1, nshell
               ib = natom + j
               IF (.NOT. spgr%lat(ib)) CYCLE
               rb(1:3) = scoord(1:3, ib)
               ro(1) = REAL(rot(1, 1), dp)*rb(1) + REAL(rot(2, 1), dp)*rb(2) + REAL(rot(3, 1), dp)*rb(3) + tr(1)
               ro(2) = REAL(rot(1, 2), dp)*rb(1) + REAL(rot(2, 2), dp)*rb(2) + REAL(rot(3, 2), dp)*rb(3) + tr(2)
               ro(3) = REAL(rot(1, 3), dp)*rb(1) + REAL(rot(2, 3), dp)*rb(2) + REAL(rot(3, 3), dp)*rb(3) + tr(3)
               ro(1) = ro(1) - REAL(NINT(ro(1) - ri(1)), dp)
               ro(2) = ro(2) - REAL(NINT(ro(2) - ri(2)), dp)
               ro(3) = ro(3) - REAL(NINT(ro(3) - ri(3)), dp)
               diff = NORM2(ri(:) - ro(:))
               IF ((diff < spgr%eps_symmetry) .AND. (spgr%atype(ia) == spgr%atype(ib))) THEN
                  spgr%eqatom(ir, ia) = ib
                  EXIT
               END IF
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE spgr_find_equivalent_atoms

! **************************************************************************************************
!> \brief routine looks for operations compatible with efield
!> \param spgr ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_reduce_symm(spgr)

      TYPE(spgr_type), INTENT(INOUT)                     :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_reduce_symm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ia, ib, ir, ja, jb, nop, nops, &
                                                            nparticle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: x, xold
      REAL(KIND=dp), DIMENSION(3)                        :: ri, ro
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot

      CALL timeset(routineN, handle)

      nop = spgr%n_operations
      nparticle = spgr%nparticle
      ALLOCATE (x(3*nparticle), xold(3*nparticle))
      x = 0.0_dp
      DO ia = 1, nparticle
         ja = 3*(ia - 1)
         x(ja + 1) = x(ja + 1) + spgr%pol(1)
         x(ja + 2) = x(ja + 2) + spgr%pol(2)
         x(ja + 3) = x(ja + 3) + spgr%pol(3)
      END DO
      xold(:) = x(:)

      nops = 0
      DO ir = 1, nop
         x = 0.d0
         spgr%lop(ir) = .TRUE.
         rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
         DO ia = 1, nparticle
            IF (.NOT. spgr%lat(ia)) CYCLE
            ja = 3*(ia - 1)
            ri(1:3) = xold(ja + 1:ja + 3)
            ro(1) = REAL(rot(1, 1), dp)*ri(1) + REAL(rot(2, 1), dp)*ri(2) + REAL(rot(3, 1), dp)*ri(3)
            ro(2) = REAL(rot(1, 2), dp)*ri(1) + REAL(rot(2, 2), dp)*ri(2) + REAL(rot(3, 2), dp)*ri(3)
            ro(3) = REAL(rot(1, 3), dp)*ri(1) + REAL(rot(2, 3), dp)*ri(2) + REAL(rot(3, 3), dp)*ri(3)
            x(ja + 1:ja + 3) = ro(1:3)
         END DO
         DO ia = 1, nparticle
            IF (.NOT. spgr%lat(ia)) CYCLE
            ib = spgr%eqatom(ir, ia)
            ja = 3*(ia - 1)
            jb = 3*(ib - 1)
            ro = x(jb + 1:jb + 3) - xold(ja + 1:ja + 3)
            spgr%lop(ir) = (spgr%lop(ir) .AND. (ABS(ro(1)) < spgr%eps_symmetry) &
                                                .AND. (ABS(ro(2)) < spgr%eps_symmetry) &
                                                .AND. (ABS(ro(3)) < spgr%eps_symmetry))
         END DO
         IF (spgr%lop(ir)) nops = nops + 1
      END DO

      spgr%n_reduced_operations = nops

      DEALLOCATE (x, xold)
      CALL timestop(handle)

   END SUBROUTINE spgr_reduce_symm

! **************************************************************************************************
!> \brief routine looks for unique rotations
!> \param spgr ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************

   SUBROUTINE spgr_rotations_subset(spgr)

      TYPE(spgr_type), INTENT(INOUT)                     :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_rotations_subset', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j
      INTEGER, DIMENSION(3, 3)                           :: d
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: mask

      CALL timeset(routineN, handle)

      ALLOCATE (mask(spgr%n_operations))
      mask = .TRUE.

      DO i = 1, spgr%n_operations
         IF (.NOT. spgr%lop(i)) mask(i) = .FALSE.
      ENDDO

      DO i = 1, spgr%n_operations - 1
         IF (.NOT. mask(i)) CYCLE
         DO j = i + 1, spgr%n_operations
            IF (.NOT. mask(j)) CYCLE
            d(:, :) = spgr%rotations(:, :, j) - spgr%rotations(:, :, i)
            IF (SUM(ABS(d)) == 0) mask(j) = .FALSE.
         END DO
      END DO

      spgr%n_operations_subset = 0
      DO i = 1, spgr%n_operations
         IF (mask(i)) spgr%n_operations_subset = spgr%n_operations_subset + 1
      END DO

      ALLOCATE (spgr%rotations_subset(3, 3, spgr%n_operations_subset))

      j = 0
      DO i = 1, spgr%n_operations
         IF (mask(i)) THEN
            j = j + 1
            spgr%rotations_subset(:, :, j) = spgr%rotations(:, :, i)
         END IF
      END DO

      DEALLOCATE (mask)
      CALL timestop(handle)

   END SUBROUTINE spgr_rotations_subset

! **************************************************************************************************
!> \brief routine applies the rotation matrices to the coordinates.
!> \param spgr ...
!> \param coord ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_apply_rotations_coord(spgr, coord)

      TYPE(spgr_type), INTENT(IN)                        :: spgr
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: coord

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_apply_rotations_coord', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ia, ib, ir, ja, jb, nop, nops, &
                                                            nparticle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cold
      REAL(KIND=dp), DIMENSION(3)                        :: rf, ri, rn, ro, tr
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot

      CALL timeset(routineN, handle)

      ALLOCATE (cold(SIZE(coord)))
      cold(:) = coord(:)

      nop = spgr%n_operations
      nparticle = spgr%nparticle
      nops = spgr%n_reduced_operations

      !$OMP PARALLEL DO PRIVATE (ia,ib,ja,jb,ir,ri,ro,rf,rn,rot,tr) SHARED (spgr,coord,nparticle,nop,nops) DEFAULT(NONE)
      DO ia = 1, nparticle
         IF (.NOT. spgr%lat(ia)) CYCLE
         ja = 3*(ia - 1)
         CALL real_to_scaled(rf(1:3), coord(ja + 1:ja + 3), spgr%cell_ref)
         rn(1:3) = 0.d0
         DO ir = 1, nop
            IF (.NOT. spgr%lop(ir)) CYCLE
            ib = spgr%eqatom(ir, ia)
            rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
            tr(1:3) = spgr%translations(1:3, ir)
            jb = 3*(ib - 1)
            CALL real_to_scaled(ri(1:3), coord(jb + 1:jb + 3), spgr%cell_ref)
            ro(1) = REAL(rot(1, 1), dp)*ri(1) + REAL(rot(2, 1), dp)*ri(2) + REAL(rot(3, 1), dp)*ri(3) + tr(1)
            ro(2) = REAL(rot(1, 2), dp)*ri(1) + REAL(rot(2, 2), dp)*ri(2) + REAL(rot(3, 2), dp)*ri(3) + tr(2)
            ro(3) = REAL(rot(1, 3), dp)*ri(1) + REAL(rot(2, 3), dp)*ri(2) + REAL(rot(3, 3), dp)*ri(3) + tr(3)
            ro(1) = ro(1) - REAL(NINT(ro(1) - rf(1)), dp)
            ro(2) = ro(2) - REAL(NINT(ro(2) - rf(2)), dp)
            ro(3) = ro(3) - REAL(NINT(ro(3) - rf(3)), dp)
            rn(1:3) = rn(1:3) + ro(1:3)
         END DO
         rn = rn/REAL(nops)
         CALL scaled_to_real(coord(ja + 1:ja + 3), rn(1:3), spgr%cell_ref)
      END DO
      !$OMP END PARALLEL DO

      DEALLOCATE (cold)
      CALL timestop(handle)

   END SUBROUTINE spgr_apply_rotations_coord

! **************************************************************************************************
!> \brief routine applies the rotation matrices to the forces.
!> \param spgr ...
!> \param force ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_apply_rotations_force(spgr, force)

      TYPE(spgr_type), INTENT(IN)                        :: spgr
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: force

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_apply_rotations_force', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ia, ib, ir, ja, jb, nop, nops, &
                                                            nparticle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fold
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rn, ro
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot

      CALL timeset(routineN, handle)

      ALLOCATE (fold(SIZE(force)))
      fold(:) = force(:)

      nop = spgr%n_operations
      nparticle = spgr%nparticle
      nops = spgr%n_reduced_operations

      !$OMP PARALLEL DO PRIVATE (ia,ib,ja,jb,ir,ri,ro,rn,rot) SHARED (spgr,force,nparticle,nop,nops) DEFAULT(NONE)
      DO ia = 1, nparticle
         IF (.NOT. spgr%lat(ia)) CYCLE
         ja = 3*(ia - 1)
         rn(1:3) = 0.d0
         DO ir = 1, nop
            IF (.NOT. spgr%lop(ir)) CYCLE
            ib = spgr%eqatom(ir, ia)
            rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
            jb = 3*(ib - 1)
            CALL real_to_scaled(ri(1:3), force(jb + 1:jb + 3), spgr%cell_ref)
            ro(1) = REAL(rot(1, 1), dp)*ri(1) + REAL(rot(2, 1), dp)*ri(2) + REAL(rot(3, 1), dp)*ri(3)
            ro(2) = REAL(rot(1, 2), dp)*ri(1) + REAL(rot(2, 2), dp)*ri(2) + REAL(rot(3, 2), dp)*ri(3)
            ro(3) = REAL(rot(1, 3), dp)*ri(1) + REAL(rot(2, 3), dp)*ri(2) + REAL(rot(3, 3), dp)*ri(3)
            rn(1:3) = rn(1:3) + ro(1:3)
         END DO
         rn = rn/REAL(nops)
         CALL scaled_to_real(force(ja + 1:ja + 3), rn(1:3), spgr%cell_ref)
      END DO
      !$OMP END PARALLEL DO

      DEALLOCATE (fold)
      CALL timestop(handle)

   END SUBROUTINE spgr_apply_rotations_force

! **************************************************************************************************
!> \brief ...
!> \param roti ...
!> \param roto ...
!> \param nop ...
!> \param h1 ...
!> \param h2 ...
! **************************************************************************************************
   SUBROUTINE spgr_change_basis(roti, roto, nop, h1, h2)

      INTEGER, DIMENSION(:, :, :)                        :: roti
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: roto
      INTEGER                                            :: nop
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h1, h2

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_change_basis', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ir
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h1ih2, h2ih1, ih1, ih2, r, s

      CALL timeset(routineN, handle)

      ih1 = inv_3x3(h1)
      ih2 = inv_3x3(h2)
      h2ih1 = MATMUL(h2, ih1)
      h1ih2 = MATMUL(h1, ih2)

      DO ir = 1, nop
         r(:, :) = roti(:, :, ir)
         s = MATMUL(h2ih1, r)
         r = MATMUL(s, h1ih2)
         roto(:, :, ir) = r(:, :)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE spgr_change_basis

! **************************************************************************************************
!> \brief routine applies the rotation matrices to the stress tensor.
!> \param spgr ...
!> \param cell ...
!> \param stress ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_apply_rotations_stress(spgr, cell, stress)

      TYPE(spgr_type), INTENT(IN)                        :: spgr
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: stress

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_apply_rotations_stress', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ir, j, k, l, nop
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: roto
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat1, hmat2, r, stin

      CALL timeset(routineN, handle)

      hmat1 = TRANSPOSE(cell%hmat)

      hmat2 = 0d0
      hmat2(1, 1) = 1.d0
      hmat2(2, 2) = 1.d0
      hmat2(3, 3) = 1.d0

      nop = spgr%n_operations_subset

      ALLOCATE (roto(3, 3, nop))

      CALL spgr_change_basis(spgr%rotations_subset, roto, spgr%n_operations_subset, hmat1, hmat2)

      stin = stress
      stress = 0.d0
      DO ir = 1, nop
         r(:, :) = roto(:, :, ir)
         DO i = 1, 3
            DO j = 1, 3
               DO k = 1, 3
                  DO l = 1, 3
                     stress(i, j) = stress(i, j) + (r(k, i)*r(l, j)*stin(k, l))
                  END DO
               END DO
            END DO
         END DO
      END DO
      stress = stress/FLOAT(nop)

      DEALLOCATE (roto)

      CALL timestop(handle)

   END SUBROUTINE spgr_apply_rotations_stress

! **************************************************************************************************
!> \brief routine prints Space Group Information.
!> \param spgr ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE print_spgr(spgr)

      TYPE(spgr_type), INTENT(IN)                        :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_spgr', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iunit, j

      CALL timeset(routineN, handle)

      IF (LEN_TRIM(spgr%filename) > 0) THEN
         CALL open_file(file_name=spgr%filename, unit_number=iunit, file_status="UNKNOWN", file_action="WRITE")
      ELSE
         iunit = spgr%punit
      END IF

      IF (iunit >= 0) THEN
         WRITE (iunit, '(/,T2,A,A)') "----------------------------------------", &
            "---------------------------------------"
         WRITE (iunit, "(T2,A,T25,A,T77,A)") "----", "SPACE GROUP SYMMETRY INFORMATION", "----"
         WRITE (iunit, '(T2,A,A)') "----------------------------------------", &
            "---------------------------------------"
         IF (spgr%symlib) THEN
            WRITE (iunit, '(T2,A,T55,A16)') "SPGR| INTERNATIONAL SYMBOL:", &
               TRIM(ADJUSTR(spgr%international_symbol))
            WRITE (iunit, '(T2,A,T55,A16)') "SPGR| POINT GROUP SYMBOL:", &
               TRIM(ADJUSTR(spgr%pointgroup_symbol))
            WRITE (iunit, '(T2,A,T55,A16)') "SPGR| SCHOENFLIES SYMBOL:", &
               TRIM(ADJUSTR(spgr%schoenflies))
            WRITE (iunit, '(T2,A,T55,I16)') "SPGR| NUMBER OF SYMMETRY OPERATIONS:", &
               spgr%n_operations
            WRITE (iunit, '(T2,A,T55,I16)') "SPGR| NUMBER OF UNIQUE ROTATIONS:", &
               spgr%n_operations_subset
            WRITE (iunit, '(T2,A,T55,I16)') "SPGR| NUMBER OF REDUCED SYMMETRY OPERATIONS:", &
               spgr%n_reduced_operations
            WRITE (iunit, *) "SPGR| REDUCED SYMMETRY OPERATIONS:", &
               spgr%lop
            WRITE (iunit, '(T2,A,T55,I8,I8)') "SPGR| NUMBER OF PARTICLES AND SYMMETRIZED PARTICLES:", &
               spgr%nparticle, spgr%nparticle_sym
            WRITE (iunit, '(T2,A,T55,I8,I8)') "SPGR| NUMBER OF ATOMS AND SYMMETRIZED ATOMS:", &
               spgr%n_atom, spgr%n_atom_sym
            WRITE (iunit, '(T2,A,T55,I8,I8)') "SPGR| NUMBER OF CORES AND SYMMETRIZED CORES:", &
               spgr%n_core, spgr%n_core_sym
            WRITE (iunit, '(T2,A,T55,I8,I8)') "SPGR| NUMBER OF SHELLS AND SYMMETRIZED SHELLS:", &
               spgr%n_shell, spgr%n_shell_sym
            IF (spgr%plevel > 0) THEN
               WRITE (iunit, '(/,T2,A,A)') "----------------------------------------", &
                  "---------------------------------------"
               WRITE (iunit, '(T2,A,T34,A,T77,A)') "----", "EQUIVALENT ATOMS", "----"
               WRITE (iunit, '(T2,A,A)') "----------------------------------------", &
                  "---------------------------------------"
               DO i = 1, spgr%nparticle
                  DO j = 1, spgr%n_operations
                     WRITE (iunit, '(T2,A,T52,I8,I8,I8)') "SPGR| ATOM | SYMMETRY OPERATION | EQUIVALENT ATOM", &
                        i, j, spgr%eqatom(j, i)
                  END DO
               END DO
               WRITE (iunit, '(T2,A,A)') "----------------------------------------", &
                  "---------------------------------------"
               DO i = 1, spgr%n_operations
                  WRITE (iunit, '(T2,A,T46,i4,T51,3I10,/,T51,3I10,/,T51,3I10)') &
                     "SPGR| SYMMETRY OPERATION #:", i, (spgr%rotations(j, :, i), j=1, 3)
                  WRITE (iunit, '(T51,3F15.7)') spgr%translations(:, i)
               END DO
            END IF
         ELSE
            WRITE (iunit, "(T2,A)") "SPGLIB for Crystal Symmetry Information determination is not availale"
         END IF
      END IF

      IF (LEN_TRIM(spgr%filename) > 0) THEN
         CALL close_file(iunit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_spgr

! **************************************************************************************************
!> \brief Variable precision output of the symmetrized stress tensor
!>
!> \param stress tensor ...
!> \param spgr ...
!> \par History
!>      07.2020 adapted to spgr [pcazade]
!> \author MK (26.08.2010).
! **************************************************************************************************
   SUBROUTINE spgr_write_stress_tensor(stress, spgr)

      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: stress
      TYPE(spgr_type), INTENT(IN)                        :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_write_stress_tensor', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=22)                                  :: fmtstr2
      CHARACTER(LEN=31)                                  :: fmtstr1
      INTEGER                                            :: handle, n
      REAL(KIND=dp), DIMENSION(3, 3)                     :: stress_tensor

      CALL timeset(routineN, handle)

      stress_tensor(:, :) = stress(:, :)*pascal*1.0E-9_dp
      n = MIN(MAX(1, spgr%ndigits), 20)
      fmtstr1 = "(/,T2,A,/,/,T13,A1,2(  X,A1))"
      WRITE (UNIT=fmtstr1(22:23), FMT="(I2)") n + 7
      fmtstr2 = "(T3,A,T5,3(1X,F  .  ))"
      WRITE (UNIT=fmtstr2(16:17), FMT="(I2)") n + 7
      WRITE (UNIT=fmtstr2(19:20), FMT="(I2)") n
      WRITE (UNIT=spgr%punit, FMT=fmtstr1) &
         "SPGR| SYMMETRIZED STRESS TENSOR [GPa]", "X", "Y", "Z"
      WRITE (UNIT=spgr%punit, FMT=fmtstr2) "X", stress_tensor(1, 1:3)
      WRITE (UNIT=spgr%punit, FMT=fmtstr2) "Y", stress_tensor(2, 1:3)
      WRITE (UNIT=spgr%punit, FMT=fmtstr2) "Z", stress_tensor(3, 1:3)

      CALL timestop(handle)

   END SUBROUTINE spgr_write_stress_tensor

END MODULE space_groups
