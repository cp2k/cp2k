!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Space Group Symmetry Module  (version 1.0, January 16, 2020)
!> \par History
!>      Pierre-André Cazade [pcazade] 01.2020 - University of Limerick
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
MODULE space_groups

   USE bibliography,                    ONLY: Togo2018,&
                                              cite_reference
   USE kinds,                           ONLY: dp, sp
   USE spglib_f08,                      ONLY: &
        spg_get_international, spg_get_ir_reciprocal_mesh, spg_get_major_version, &
        spg_get_micro_version, spg_get_minor_version, spg_get_multiplicity, spg_get_pointgroup, &
        spg_get_schoenflies, spg_get_symmetry
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE gopt_f_types,                    ONLY: gopt_f_type
   USE cell_types,                      ONLY: cell_type, cell_copy,&
                                              cell_create, init_cell,&
                                              cell_release,&
                                              real_to_scaled, scaled_to_real
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE input_constants,                 ONLY: default_cell_direct_id,&
                                              default_cell_geo_opt_id,&
                                              default_cell_md_id,&
                                              default_cell_method_id,&
                                              default_minimization_method_id,&
                                              default_shellcore_method_id,&
                                              default_ts_method_id
   USE cp_fm_types,                     ONLY: cp_fm_type, cp_fm_get_info
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE message_passing,                 ONLY: mp_sum
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE physcon,                         ONLY: pascal
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'space_groups'
   
   TYPE spgr_type
      LOGICAL                                     :: symlib = .FALSE.
      INTEGER                                     :: plevel = 0
      INTEGER                                     :: punit = 6
      INTEGER                                     :: istriz = -1
      REAL(KIND=dp)                               :: eps_geo = 1.0e-4_dp
      REAL(KIND=dp)                               :: eps_sym = 1.0e-4_dp
      INTEGER                                     :: nparticle
      INTEGER                                     :: nparticle_sym
      INTEGER                                     :: n_atom
      INTEGER                                     :: n_core
      INTEGER                                     :: n_shell
      INTEGER                                     :: n_atom_sym
      INTEGER                                     :: n_core_sym
      INTEGER                                     :: n_shell_sym
      INTEGER                                     :: ndigits = 16
      INTEGER, DIMENSION(:), ALLOCATABLE          :: atype
      INTEGER, DIMENSION(:,:), ALLOCATABLE        :: eqatom
      LOGICAL, DIMENSION(:), ALLOCATABLE          :: lop,lat
      REAL(KIND=dp), DIMENSION(3)                 :: pol = 0.0_dp
      !SPGLIB
      CHARACTER(len=11)                           :: international_symbol
      CHARACTER(len=6)                            :: pointgroup_symbol
      CHARACTER(len=7)                            :: schoenflies
      INTEGER                                     :: n_operations
      INTEGER                                     :: n_reduced_operations
      INTEGER                                     :: n_operations_subset
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE    :: rotations
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE    :: rotations_subset
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: translations
      TYPE(cell_type),POINTER                     :: cell_ref
   END TYPE spgr_type
   
   PUBLIC :: release_spgr_type, spgr_gen, identify_space_group, spgr_find_equivalent_atoms
   PUBLIC :: spgr_apply_rotations_coord, spgr_apply_rotations_force, print_spgr, spgr_type
   PUBLIC :: spgr_print_cell, spgr_print_xg, spgr_print_hess, spgr_print_hess_sp
   PUBLIC :: spgr_apply_rotations_hess, spgr_apply_rotations_stress, spgr_write_stress_tensor

CONTAINS

! **************************************************************************************************
!> \brief Release the SPGR type
!> \param spgr  The SPGR type
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE release_spgr_type(spgr)
      
      IMPLICIT NONE
      
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
!> \param eps_geo ...
!> \param iounit ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_gen(spgr, scoor, types, cell, gopt_env, eps_geo, eps_sym, pol, ranges, &
                       nparticle, n_atom, n_core, n_shell, iounit)
      
      IMPLICIT NONE
      
      TYPE(spgr_type)                                    :: spgr
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: scoor
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: tmp_coor
      INTEGER, DIMENSION(:), INTENT(IN)                  :: types
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: tmp_types
      INTEGER, DIMENSION(:,:), INTENT(IN), OPTIONAL      :: ranges
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_geo, eps_sym
      REAL(KIND=dp), DIMENSION(3), OPTIONAL              :: pol
      INTEGER, OPTIONAL                                  :: nparticle, n_atom, n_core, n_shell
      INTEGER, INTENT(IN), OPTIONAL                      :: iounit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_gen', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ierr, major, micro, minor, &
                                                            nop, tra_mat(3, 3), i, j, n_sr_rep
      LOGICAL                                            :: spglib

      CALL timeset(routineN, handle)

      !..total number of particles (atoms plus shells)
      IF (PRESENT(nparticle)) THEN
         CPASSERT(nparticle==SIZE(scoor, 2))
         spgr%nparticle = nparticle
      ELSE
         spgr%nparticle = SIZE(scoor, 2)
      ENDIF

      IF (PRESENT(n_atom)) THEN
         spgr%n_atom = n_atom
      ELSE IF (PRESENT(n_core)) THEN
         spgr%n_atom = nparticle - n_core
      ELSE IF (PRESENT(n_shell)) THEN
         spgr%n_atom = nparticle - n_shell
      ELSE
         spgr%n_atom = nparticle
      ENDIF

      IF (PRESENT(n_core)) THEN
         spgr%n_core = n_core
      ELSE IF (PRESENT(n_shell)) THEN
         spgr%n_core = n_shell
      ELSE
         spgr%n_core = 0
      ENDIF

      IF (PRESENT(n_shell)) THEN
         spgr%n_shell = n_shell
      ELSE IF (PRESENT(n_core)) THEN
         spgr%n_shell = n_core
      ELSE
         spgr%n_shell = 0
      ENDIF

      IF ( .NOT.( spgr%nparticle == (spgr%n_atom + spgr%n_shell) ) ) THEN
         CPABORT("spgr_gen: nparticle not equal to natom + nshell.")
      END IF

      spgr%nparticle_sym = spgr%nparticle
      spgr%n_atom_sym = spgr%n_atom
      spgr%n_core_sym = spgr%n_core
      spgr%n_shell_sym = spgr%n_shell

      ! output unit
      IF (PRESENT(iounit)) THEN
         spgr%punit = iounit
      ELSE
         spgr%punit = 6
      END IF

      ! accuracy for symmetry
      IF (PRESENT(eps_geo)) THEN
         spgr%eps_geo = eps_geo
      ELSE
         spgr%eps_geo = 1.e-4_dp
      END IF

      ! accuracy for symmetry
      IF (PRESENT(eps_sym)) THEN
         spgr%eps_sym = eps_sym
      ELSE
         spgr%eps_sym = 1.e-4_dp
      END IF

      ! vector to test reduced symmetry
      IF (PRESENT(pol)) THEN
         spgr%pol(1) = pol(1)
         spgr%pol(2) = pol(2)
         spgr%pol(3) = pol(3)
      ELSE
         spgr%pol = 0.0_dp
      END IF

      ALLOCATE(spgr%lat(nparticle))
      spgr%lat=.TRUE.

      IF (PRESENT(ranges)) THEN
        n_sr_rep = SIZE(ranges, 2)
        DO i = 1, n_sr_rep
          DO j = ranges(1, i), ranges(2, i)
            spgr%lat(j)=.FALSE.
            spgr%nparticle_sym=spgr%nparticle_sym-1
            IF (j <= spgr%n_atom) THEN
               spgr%n_atom_sym = spgr%n_atom_sym -1
            ELSE IF (j > spgr%n_atom .AND. j <= spgr%nparticle) THEN
               spgr%n_core_sym = spgr%n_core_sym -1
               spgr%n_shell_sym = spgr%n_shell_sym -1
            ELSE
               CPABORT("Symmetry exclusion range larger than actual number of particles.")
            ENDIF
          END DO
        END DO
      END IF

      ALLOCATE(tmp_coor(3, spgr%nparticle_sym), tmp_types(spgr%nparticle_sym))

      j = 0
      DO i = 1, nparticle
        IF (spgr%lat(i)) THEN
          j = j+1
          tmp_coor(:, j) = scoor(:, i)
          tmp_types(j) = types(i)
        END IF
      END DO

      !..set cell values
      NULLIFY(spgr%cell_ref)
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
            CPABORT("")
         END SELECT
      CASE DEFAULT
         CPABORT("")
      END SELECT

      ! atom types
      ALLOCATE (spgr%atype(nparticle))
      spgr%atype(1:nparticle) = types(1:nparticle)

      spgr%n_operations = 0

      !..try spglib
      major = spg_get_major_version()
      minor = spg_get_minor_version()
      micro = spg_get_micro_version()
      IF (major == 0) THEN
         CALL cp_warn(__LOCATION__, "Symmetry library SPGLIB not available")
         spglib = .FALSE.
      ELSE
         !IF (major /= 1 .OR. minor /= 9 .OR. micro /= 9) THEN
         !   CALL cp_warn(__LOCATION__, "Version of Symmetry Library SPGLIB not tested")
         !END IF
         spglib = .TRUE.
         CALL cite_reference(Togo2018)
         spgr%international_symbol = '           '
         ierr = spg_get_international(spgr%international_symbol, TRANSPOSE(cell%hmat), tmp_coor, tmp_types, &
                                      spgr%nparticle_sym, eps_sym)
         IF (ierr == 0) THEN
            CALL cp_warn(__LOCATION__, "Symmetry Library SPGLIB failed")
            spglib = .FALSE.
         ELSE
            nop = spg_get_multiplicity(TRANSPOSE(cell%hmat), tmp_coor, tmp_types, &
                                       spgr%nparticle_sym, eps_sym)
            ALLOCATE (spgr%rotations(3, 3, nop), spgr%translations(3, nop))
            ALLOCATE (spgr%eqatom(nop, nparticle))
            ALLOCATE (spgr%lop(nop))
            spgr%n_operations = nop
            spgr%lop = .TRUE.
            ierr = spg_get_symmetry(spgr%rotations, spgr%translations, nop, TRANSPOSE(cell%hmat), &
                                    tmp_coor, tmp_types, spgr%nparticle_sym, eps_sym)
            ! Schoenflies Symbol
            spgr%schoenflies = '       '
            ierr = spg_get_schoenflies(spgr%schoenflies, TRANSPOSE(cell%hmat), tmp_coor, tmp_types, &
                                       spgr%nparticle_sym, eps_sym)
            ! Point Group
            spgr%pointgroup_symbol = '      '
            tra_mat = 0
            ierr = spg_get_pointgroup(spgr%pointgroup_symbol, tra_mat, &
                                      spgr%rotations, spgr%n_operations)
         END IF
      END IF
      spgr%symlib = spglib

      DEALLOCATE(tmp_coor, tmp_types)

      CALL timestop(handle)

   END SUBROUTINE spgr_gen
   
! **************************************************************************************************
!> \brief routine indentifies the space group and finds rotation matrices.
!> \param spgr ...
!> \param particle_set ...
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
   SUBROUTINE identify_space_group(spgr,subsys, cell, geo_section, gopt_env, iounit)
    
      IMPLICIT NONE
      
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      !TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_list_type), POINTER                  :: particles, &
                                                            shell_particles, &
                                                            core_particles

      CHARACTER(LEN=*), PARAMETER :: routineN = 'identify_space_group', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iounit, nparticle
      INTEGER                                            :: n_sr_rep, n_atom, n_core
      INTEGER                                            :: k, shell_index, n_shell
      INTEGER, DIMENSION(:), POINTER                     :: tmp
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atype
      INTEGER, ALLOCATABLE, DIMENSION(:,:)               :: ranges
      REAL(KIND=dp)                                      :: eps_geo, eps_sym
      REAL(KIND=dp), DIMENSION(:), POINTER               :: pol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: scoord
      TYPE(spgr_type)                                    :: spgr

      CALL timeset(routineN, handle)

      n_sr_rep = 0
      nparticle = 0
      n_atom = 0
      n_core = 0
      n_shell = 0

      NULLIFY (particles)
      NULLIFY (core_particles)
      NULLIFY (shell_particles)

      particles => subsys%particles
      shell_particles => subsys%shell_particles
      core_particles => subsys%core_particles

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
         CPASSERT(ASSOCIATED(shell_particles))
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

      CALL section_vals_val_get(geo_section, "EPS_GEO", r_val=eps_geo)
      CALL section_vals_val_get(geo_section, "EPS_SPGR", r_val=eps_sym)
      CALL section_vals_val_get(geo_section, "SYMM_REDUCTION", r_vals=pol)
      CALL section_vals_val_get(geo_section, "SYMM_EXCLUDE_RANGE", n_rep_val=n_sr_rep)
      IF(n_sr_rep > 0)THEN
        ALLOCATE(ranges(2, n_sr_rep))
        DO i = 1, n_sr_rep
          CALL section_vals_val_get(geo_section, "SYMM_EXCLUDE_RANGE", i_rep_val=i, i_vals=tmp)
          ranges(:,i)=tmp(:)
        END DO
        CALL spgr_gen(spgr, scoord, atype, cell, gopt_env, eps_geo=eps_geo, eps_sym=eps_sym, pol=pol(1:3), &
                      ranges=ranges, nparticle=nparticle, n_atom=n_atom, n_core=n_core, n_shell=n_shell, iounit=iounit)
        DEALLOCATE(ranges)
      ELSE
        CALL spgr_gen(spgr, scoord, atype, cell, gopt_env, eps_geo=eps_geo, eps_sym=eps_sym, pol=pol(1:3), &
                      nparticle=nparticle, n_atom=n_atom, n_core=n_core, n_shell=n_shell, iounit=iounit)
      ENDIF
      CALL spgr_find_equivalent_atoms(spgr, scoord)
      CALL spgr_reduce_symm(spgr)
      CALL spgr_rotations_subset(spgr)
      DEALLOCATE(scoord, atype)
      
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
      
      IMPLICIT NONE
      
      TYPE(spgr_type)                                    :: spgr
      
      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_find_equivalent_atoms', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, ia, ib, ir, handle
      INTEGER                                            :: natom, nshell, nop
      REAL(KIND=dp)                                      :: diff
      REAL(KIND=dp), DIMENSION(3)                        :: rb, ri, ro, tr
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: scoord
      
      CALL timeset(routineN, handle)
      
      nop=spgr%n_operations
      natom = spgr%n_atom
      nshell = spgr%n_shell

      IF ( .NOT.( spgr%nparticle == (natom + nshell) ) ) THEN
         CPABORT("spgr_find_equivalent_atoms: nparticle not equal to natom + nshell.")
      END IF

      DO ia = 1, spgr%nparticle
         spgr%eqatom(:,ia) = ia
      END DO

      !$OMP PARALLEL DO PRIVATE (ia,ib,ir,ri,rb,ro,rot,tr,diff) SHARED (spgr,scoord,natom,nop)
      DO ia = 1, natom
         IF (.NOT.spgr%lat(ia)) CYCLE
         ri(1:3) = scoord(1:3, ia)
         DO ir = 1, nop
            rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
            tr(1:3) = spgr%translations(1:3, ir)
            DO ib = 1, natom
               IF (.NOT.spgr%lat(ib)) CYCLE
               rb(1:3) = scoord(1:3, ib)
               ro(1) = REAL(rot(1, 1), dp)*rb(1)+REAL(rot(2, 1), dp)*rb(2)+REAL(rot(3, 1), dp)*rb(3)+tr(1)
               ro(2) = REAL(rot(1, 2), dp)*rb(1)+REAL(rot(2, 2), dp)*rb(2)+REAL(rot(3, 2), dp)*rb(3)+tr(2)
               ro(3) = REAL(rot(1, 3), dp)*rb(1)+REAL(rot(2, 3), dp)*rb(2)+REAL(rot(3, 3), dp)*rb(3)+tr(3)
               ro(1) = ro(1) - REAL(NINT(ro(1)-ri(1)))
               ro(2) = ro(2) - REAL(NINT(ro(2)-ri(2)))
               ro(3) = ro(3) - REAL(NINT(ro(3)-ri(3)))
               diff = SQRT(SUM((ri(:)-ro(:))**2))
               IF ( (diff < spgr%eps_geo) .AND. (spgr%atype(ia)==spgr%atype(ib)) ) THEN
                  spgr%eqatom(ir,ia) = ib
                  EXIT
               END IF
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE (i,j,ia,ib,ir,ri,rb,ro,rot,tr,diff) SHARED (spgr,scoord,natom,nshell,nop)
      DO i = 1, nshell
         ia = natom + i 
         IF (.NOT.spgr%lat(ia)) CYCLE
         ri(1:3) = scoord(1:3, ia)
         DO ir = 1, nop
            rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
            tr(1:3) = spgr%translations(1:3, ir)
            DO j = 1, nshell
               ib = natom + j
               IF (.NOT.spgr%lat(ib)) CYCLE
               rb(1:3) = scoord(1:3, ib)
               ro(1) = REAL(rot(1, 1), dp)*rb(1)+REAL(rot(2, 1), dp)*rb(2)+REAL(rot(3, 1), dp)*rb(3)+tr(1)
               ro(2) = REAL(rot(1, 2), dp)*rb(1)+REAL(rot(2, 2), dp)*rb(2)+REAL(rot(3, 2), dp)*rb(3)+tr(2)
               ro(3) = REAL(rot(1, 3), dp)*rb(1)+REAL(rot(2, 3), dp)*rb(2)+REAL(rot(3, 3), dp)*rb(3)+tr(3)
               ro(1) = ro(1) - REAL(NINT(ro(1)-ri(1)))
               ro(2) = ro(2) - REAL(NINT(ro(2)-ri(2)))
               ro(3) = ro(3) - REAL(NINT(ro(3)-ri(3)))
               diff = SQRT(SUM((ri(:)-ro(:))**2))
               IF ( (diff < spgr%eps_geo) .AND. (spgr%atype(ia)==spgr%atype(ib)) ) THEN
                  spgr%eqatom(ir,ia) = ib
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
!> \param force ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_reduce_symm(spgr)
      
      IMPLICIT NONE
      
      TYPE(spgr_type)                                    :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_apply_rotations_force', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ia, ib, ir, nparticle, nop, handle, &
                                                            ja, jb, nops
      REAL(KIND=dp), DIMENSION(3)                        :: ri, ro
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: x, xold
      
      CALL timeset(routineN, handle)

      nop=spgr%n_operations
      nparticle=spgr%nparticle
      allocate(x(3*nparticle), xold(3*nparticle))
      x=0.0_dp
      DO ia = 1, nparticle
         ja = 3*(ia-1)
         x(ja+1) = x(ja+1) + spgr%pol(1)
         x(ja+2) = x(ja+2) + spgr%pol(2)
         x(ja+3) = x(ja+3) + spgr%pol(3)
      END DO
      xold=x

      nops = 0
      DO ir=1,nop
         x=0.d0
         spgr%lop(ir)=.TRUE.
         rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
         DO ia = 1, nparticle
            IF (.NOT.spgr%lat(ia)) CYCLE
            ja = 3*(ia-1)
            ri(1:3) = xold(ja+1:ja+3)
            ro(1) = REAL(rot(1, 1), dp)*ri(1)+REAL(rot(2, 1), dp)*ri(2)+REAL(rot(3, 1), dp)*ri(3)
            ro(2) = REAL(rot(1, 2), dp)*ri(1)+REAL(rot(2, 2), dp)*ri(2)+REAL(rot(3, 2), dp)*ri(3)
            ro(3) = REAL(rot(1, 3), dp)*ri(1)+REAL(rot(2, 3), dp)*ri(2)+REAL(rot(3, 3), dp)*ri(3)
            x(ja+1:ja+3) = ro(1:3)
         END DO
         DO ia = 1, nparticle
            IF (.NOT.spgr%lat(ia)) CYCLE
            ib = spgr%eqatom(ir, ia)
            ja = 3*(ia-1)
            jb = 3*(ib-1)
            ro=x(jb+1:jb+3)-xold(ja+1:ja+3)
            spgr%lop(ir)=(spgr%lop(ir).AND.(abs(ro(1))<spgr%eps_geo).AND.(abs(ro(2))<spgr%eps_geo).AND.(abs(ro(3))<spgr%eps_geo))
         END DO
         IF(spgr%lop(ir)) nops = nops +1
      END DO

      spgr%n_reduced_operations=nops

      deallocate(x, xold)
      CALL timestop(handle)
   
   END SUBROUTINE spgr_reduce_symm

! **************************************************************************************************
!> \brief routine looks for unique rotations
!> \param spgr ...
!> \param force ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************

   SUBROUTINE spgr_rotations_subset(spgr)

     IMPLICIT NONE

     CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_rotations_subset', &
         routineP = moduleN//':'//routineN
      
     TYPE(spgr_type)                                    :: spgr
     LOGICAL, DIMENSION(:), ALLOCATABLE                 :: mask
     INTEGER                                            :: i, j, handle
     INTEGER, DIMENSION(3, 3)                           :: d

     CALL timeset(routineN, handle)

     ALLOCATE(mask(spgr%n_operations))
     mask=.TRUE.

     DO i = 1, spgr%n_operations
       IF (.NOT. spgr%lop(i)) mask(i) = .FALSE.
     ENDDO
     
     DO i = 1, spgr%n_operations-1
       IF (.NOT. mask(i)) CYCLE
       DO j = i+1, spgr%n_operations
         IF (.NOT. mask(j)) CYCLE
         d(:,:) = spgr%rotations(:, :, j) - spgr%rotations(:, :, i)
         IF (SUM(ABS(d)) == 0) mask(j) = .FALSE.
       END DO
     END DO

     spgr%n_operations_subset=0
     DO i = 1, spgr%n_operations
       IF(mask(i)) spgr%n_operations_subset = spgr%n_operations_subset+1
     END DO

     ALLOCATE(spgr%rotations_subset(3, 3, spgr%n_operations_subset))

     j = 0
     DO i = 1, spgr%n_operations
       IF (mask(i)) THEN
         j = j+1
         spgr%rotations_subset(:, :, j) = spgr%rotations(:, :, i)
       END IF
     END DO

     DEALLOCATE(mask)
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
      
      IMPLICIT NONE
      
      TYPE(spgr_type)                                    :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_apply_rotations_coord', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ia, ib, ir, nparticle, nop, nops, handle, &
                                                            ja, jb
      REAL(KIND=dp), DIMENSION(3)                        :: ri, ro, tr, rf, rn
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot
      REAL(KIND=dp), DIMENSION(:)                        :: coord
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: cold

      INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
      
      CALL timeset(routineN, handle)

      allocate(cold(size(coord)))
      cold=coord
      nop=spgr%n_operations
      nparticle = spgr%nparticle
      nops = spgr%n_reduced_operations
      
      !$OMP PARALLEL DO PRIVATE (ia,ib,ja,jb,ir,ri,ro,rf,rn,rot,tr) SHARED (spgr,coord,nparticle,nop)
      DO ia = 1, nparticle
         IF (.NOT.spgr%lat(ia)) CYCLE
         ja = 3*(ia-1)
         CALL real_to_scaled(rf(1:3), coord(ja+1:ja+3), spgr%cell_ref)
         rn(1:3) = 0.d0
         DO ir = 1, nop
            IF(.NOT.spgr%lop(ir)) CYCLE
            ib = spgr%eqatom(ir, ia)
            rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
            tr(1:3) = spgr%translations(1:3, ir)
            jb = 3*(ib-1)
            CALL real_to_scaled(ri(1:3), coord(jb+1:jb+3), spgr%cell_ref)
            ro(1) = REAL(rot(1, 1), dp)*ri(1)+REAL(rot(2, 1), dp)*ri(2)+REAL(rot(3, 1), dp)*ri(3)+tr(1)
            ro(2) = REAL(rot(1, 2), dp)*ri(1)+REAL(rot(2, 2), dp)*ri(2)+REAL(rot(3, 2), dp)*ri(3)+tr(2)
            ro(3) = REAL(rot(1, 3), dp)*ri(1)+REAL(rot(2, 3), dp)*ri(2)+REAL(rot(3, 3), dp)*ri(3)+tr(3)
            ro(1) = ro(1) - REAL(NINT(ro(1)-rf(1)))
            ro(2) = ro(2) - REAL(NINT(ro(2)-rf(2)))
            ro(3) = ro(3) - REAL(NINT(ro(3)-rf(3)))
            rn(1:3) = rn(1:3) + ro(1:3)
         END DO
         rn = rn / REAL(nops)
         CALL scaled_to_real(coord(ja+1:ja+3), rn(1:3), spgr%cell_ref)
      END DO
      !$OMP END PARALLEL DO

      deallocate(cold)
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
      
      IMPLICIT NONE
      
      TYPE(spgr_type)                                    :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_apply_rotations_force', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ia, ib, ir, nparticle, nop, nops, handle, &
                                                            ja, jb
      REAL(KIND=dp), DIMENSION(3)                        :: ri, ro, rn
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot
      REAL(KIND=dp), DIMENSION(:)                        :: force
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: fold
      
      CALL timeset(routineN, handle)

      allocate(fold(size(force)))
      fold=force
      nop=spgr%n_operations
      nparticle = spgr%nparticle
      nops = spgr%n_reduced_operations

      !$OMP PARALLEL DO PRIVATE (ia,ib,ja,jb,ir,ri,ro,rn,rot) SHARED (spgr,force,nparticle,nop)
      DO ia = 1, nparticle
         IF (.NOT.spgr%lat(ia)) CYCLE
         ja = 3*(ia-1)
         rn(1:3) = 0.d0
         DO ir=1,nop
            IF(.NOT.spgr%lop(ir)) CYCLE
            ib = spgr%eqatom(ir, ia)
            rot(1:3, 1:3) = spgr%rotations(1:3, 1:3, ir)
            jb = 3*(ib-1)
            CALL real_to_scaled(ri(1:3), force(jb+1:jb+3), spgr%cell_ref)
            ro(1) = REAL(rot(1, 1), dp)*ri(1)+REAL(rot(2, 1), dp)*ri(2)+REAL(rot(3, 1), dp)*ri(3)
            ro(2) = REAL(rot(1, 2), dp)*ri(1)+REAL(rot(2, 2), dp)*ri(2)+REAL(rot(3, 2), dp)*ri(3)
            ro(3) = REAL(rot(1, 3), dp)*ri(1)+REAL(rot(2, 3), dp)*ri(2)+REAL(rot(3, 3), dp)*ri(3)
            rn(1:3) = rn(1:3) + ro(1:3)
         END DO
         rn = rn / REAL(nops)
         CALL scaled_to_real(force(ja+1:ja+3), rn(1:3), spgr%cell_ref)
      END DO
      !$OMP END PARALLEL DO

      deallocate(fold)
      CALL timestop(handle)
   
   END SUBROUTINE spgr_apply_rotations_force

   SUBROUTINE spgr_change_basis(roti, roto, nop, h1, h2)

     IMPLICIT NONE

     CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_change_basis', &
         routineP = moduleN//':'//routineN

     INTEGER                           :: ir, nop, handle
     INTEGER, DIMENSION(:, :, :)       :: roti
     REAL(KIND=dp), DIMENSION(3, 3)    :: r, s, h1, h2, ih1, ih2, h1ih2, h2ih1
     REAL(KIND=dp), DIMENSION(:, :, :) :: roto

     CALL timeset(routineN, handle)

     ih1 = inv_3x3(h1)
     ih2 = inv_3x3(h2)
     h2ih1 = MATMUL(h2,ih1)
     h1ih2 = MATMUL(h1,ih2)

     DO ir = 1, nop
       r(:,:) = roti(:, :, ir)
       s = MATMUL(h2ih1, r)
       r = MATMUL(s, h1ih2)
       roto(:, :, ir) = r(:, :)
     ENDDO

     CALL timestop(handle)

   END SUBROUTINE spgr_change_basis

! **************************************************************************************************
!> \brief routine applies the rotation matrices to the stress tensor.
!> \param spgr ...
!> \param stress ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_apply_rotations_stress(spgr, cell, stress)
     
     IMPLICIT NONE
     
     TYPE(spgr_type)                                    :: spgr
     TYPE(cell_type), POINTER                           :: cell

     CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_apply_rotations_force', &
        routineP = moduleN//':'//routineN

     INTEGER                                            :: i, j, k, l, ir, nop, handle
     REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE       :: roto
     REAL(KIND=dp), DIMENSION(3, 3)                     :: r, stress, stin, hmat1, hmat2

     CALL timeset(routineN, handle)

     hmat1 = TRANSPOSE(cell%hmat)

     hmat2 = 0d0
     hmat2(1, 1) = 1.d0
     hmat2(2, 2) = 1.d0
     hmat2(3, 3) = 1.d0

     nop=spgr%n_operations_subset

     ALLOCATE(roto(3, 3, nop))
     
     CALL spgr_change_basis(spgr%rotations_subset, roto, spgr%n_operations_subset, hmat1, hmat2)

     stin = stress
     stress=0.d0
     DO ir = 1, nop
       r(:, :) = roto(:, :, ir)
       DO i = 1, 3
         DO j = 1, 3
           DO k = 1, 3
             DO l = 1, 3
               stress(i, j) = stress(i, j) + ( r(k, i) * r(l, j) * stin(k, l) )
             END DO
           END DO
         END DO
       END DO
     END DO
     stress = stress / FLOAT(nop)
     
     DEALLOCATE(roto)

     CALL timestop(handle)
   
   END SUBROUTINE spgr_apply_rotations_stress

! **************************************************************************************************
!> \brief routine applies the rotation matrices to the hessian.
!> \param ndf ...
!> \param spgr ...
!> \param hess_mat ...
!> \param para_env ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
!> \note  Not currently in use and most likely unnecessary. Will be removed in later version.
! **************************************************************************************************
   SUBROUTINE spgr_apply_rotations_hess(ndf,spgr,hess_mat,para_env)

      IMPLICIT NONE

      INTEGER, INTENT(IN)                                :: ndf
      TYPE(spgr_type)                                    :: spgr
      TYPE(cp_fm_type), POINTER                          :: hess_mat
      TYPE(cp_para_env_type), OPTIONAL, POINTER          :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'pgr_apply_rotations_hess',&
              routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, k, l, ncol_local, &
                                                            nrow_local, ia, ib, nparticle, &
                                                            nop, ja, jb, ka, kb, la, lb, &
                                                            iri, irk
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot, roti, rotk, r, s, t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:)         :: work
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: local_hes


      CALL timeset(routineN,handle)

      ALLOCATE(work(ndf,ndf))

      CALL cp_fm_get_info(hess_mat,row_indices=row_indices,col_indices=col_indices, &
           local_data=local_hes,nrow_local=nrow_local,ncol_local=ncol_local)

      work=0.0_dp
      DO i=1,nrow_local
         j=row_indices(i)
         DO k=1,ncol_local
            l=col_indices(k)
            work(j,l)=local_hes(i,k)
         END DO
      END DO

      CALL mp_sum(work,para_env%group)

      nop=spgr%n_operations
      nparticle = spgr%nparticle
      DO ia = 1, nparticle
         ib = spgr%eqatom(1, ia)
         iri = spgr%eqatom(2, ia)
         roti(1:3, 1:3) = spgr%rotations(1:3, 1:3, iri)
         ja = 3*(ia-1)
         jb = 3*(ib-1)
         DO ka = ia, nparticle
            r=0.0_dp
            kb = spgr%eqatom(1, ka)
            IF (ib==ia .AND. kb==ka) THEN
               CYCLE
            END IF
            irk = spgr%eqatom(2, ka)
            rotk(1:3, 1:3) = spgr%rotations(1:3, 1:3, irk)
            la = 3*(ka-1)
            lb = 3*(kb-1)
            rot=TRANSPOSE(MATMUL(rotk,roti))
            r(:,:) = work(jb+1:jb+3,lb+1:lb+3)
            s = MATMUL(spgr%cell_ref%h_inv,r)
            t = MATMUL(rot,s)
            r = MATMUL(spgr%cell_ref%hmat,t)
            work(ja+1:ja+3,la+1:la+3) = r(:,:)
            work(la+1:la+3,ja+1:ja+3) = TRANSPOSE(r)
         END DO
      END DO
      
      DO i=1,nrow_local
         j=row_indices(i)
         DO k=1,ncol_local
            l=col_indices(k)
            local_hes(i,k) = work(j,l)
         END DO
      END DO

      DEALLOCATE(work)

      CALL timestop(handle)
      
   END SUBROUTINE spgr_apply_rotations_hess

! **************************************************************************************************
!> \brief routine prints Space Group Information.
!> \param spgr ...
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE print_spgr(spgr)
      
      IMPLICIT NONE
      
      TYPE(spgr_type)                                    :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_spgr', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iunit, j, plevel, handle

      CALL timeset(routineN,handle)

      iunit = spgr%punit
      IF (iunit >= 0) THEN
         plevel = spgr%plevel+5
         WRITE (iunit, '(/,T2,A,A)') "----------------------------------------", &
            "---------------------------------------"
         WRITE (iunit, "(T2,A,T25,A,T77,A)") "----", "SPACE GROUP SYMMETRY INFORMATION", "----"
         WRITE (iunit, '(T2,A,A)') "----------------------------------------", &
            "---------------------------------------"
         IF (spgr%symlib) THEN
            CALL spgr_print_cell(spgr%cell_ref,0,iunit)
            WRITE (iunit, *) "SPGR| INTERNATIONAL SYMBOL: ", TRIM(ADJUSTR(spgr%international_symbol))
            WRITE (iunit, *) "SPGR| POINT GROUP SYMBOL  : ", TRIM(ADJUSTR(spgr%pointgroup_symbol))
            WRITE (iunit, *) "SPGR| SCHOENFLIES SYMBOL  : ", TRIM(ADJUSTR(spgr%schoenflies))
            WRITE (iunit, *) "SPGR| NUMBER OF SYMMETRY OPERATIONS:", spgr%n_operations
            WRITE (iunit, *) "SPGR| NUMBER OF UNIQUE ROTATIONS:", spgr%n_operations_subset
            WRITE (iunit, *) "SPGR| NUMBER OF REDUCED SYMMETRY OPERATIONS:", spgr%n_reduced_operations
            WRITE (iunit, *) "SPGR| REDUCED SYMMETRY OPERATIONS:", spgr%lop
            WRITE (iunit, *) "SPGR| NUMBER OF PARTICLES AND SYMMETRIZED PARTICLES:", spgr%nparticle, spgr%nparticle_sym
            WRITE (iunit, *) "SPGR| NUMBER OF ATOMS AND SYMMETRIZED ATOMS:", spgr%n_atom, spgr%n_atom_sym
            WRITE (iunit, *) "SPGR| NUMBER OF CORES AND SYMMETRIZED CORES:", spgr%n_core, spgr%n_core_sym
            WRITE (iunit, *) "SPGR| NUMBER OF SHELLS AND SYMMETRIZED SHELLS:", spgr%n_shell, spgr%n_shell_sym
            IF (plevel > 10) THEN
               WRITE (iunit, '(/,T2,A,A)') "----------------------------------------", &
               "---------------------------------------"
               WRITE (iunit, '(T2,A,T34,A,T77,A)') "----", "EQUIVALENT ATOMS", "----"
               WRITE (iunit, '(T2,A,A)') "----------------------------------------", &
               "---------------------------------------"
               DO i = 1, spgr%nparticle
                  DO j = 1, spgr%n_operations
                     WRITE (iunit, '(T36,I8,I8,I8)') i, j, spgr%eqatom(j,i)
                  END DO
               END DO
               WRITE (iunit, '(T2,A,A)') "----------------------------------------", &
               "---------------------------------------"
               DO i = 1, spgr%n_operations
                  WRITE (iunit, '(A,i4,T51,3I10,/,T51,3I10,/,T51,3I10)') &
                     "           Rotation #: ", i, (spgr%rotations(j, :, i), j=1, 3)
                  WRITE (iunit, '(T36,3F15.7)') spgr%translations(:, i)
               END DO
            END IF
         ELSE
            WRITE (iunit, "(T2,A)") "SPGLIB for Crystal Symmetry Information determination is not availale"
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_spgr

! **************************************************************************************************
!> \brief Variable precision output of the symmetrized stress tensor
!>
!> \param stress tensor ...
!> \param cell ...
!> \param spgr ...
!> \author MK (26.08.2010). pcazade (01.05.2020)
! **************************************************************************************************
   SUBROUTINE spgr_write_stress_tensor(stress, spgr)

     REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: stress
     TYPE(spgr_type)                                    :: spgr
     
     CHARACTER(LEN=*), PARAMETER :: routineN = 'spgr_write_stress_tensor', &
        routineP = moduleN//':'//routineN
     
     CHARACTER(LEN=22)                                  :: fmtstr2
     CHARACTER(LEN=31)                                  :: fmtstr1
     INTEGER                                            :: n, handle
     REAL(KIND=dp), DIMENSION(3, 3)                     :: stress_tensor

     CALL timeset(routineN,handle)
     
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

! **************************************************************************************************
!> \brief routine prints cell details. Only for debugging purpose.
!> \param cell ...
!> \param iter ...
!> \param io output unit
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_print_cell(cell,iter,io)
      IMPLICIT NONE
      
      INTEGER                                            :: iter,io
      TYPE(cell_type)                                    :: cell

      WRITE (io, *) "SPGR| ITERATION: ",iter
      WRITE (io, *) "SPGR| CELL_REF| A: ", cell%hmat(:,1)*0.52917721067
      WRITE (io, *) "SPGR| CELL_REF| B: ", cell%hmat(:,2)*0.52917721067
      WRITE (io, *) "SPGR| CELL_REF| C: ", cell%hmat(:,3)*0.52917721067

   END SUBROUTINE spgr_print_cell

! **************************************************************************************************
!> \brief routine prints coordinates or forces. Only for debugging purpose.
!> \param x ...
!> \param iter ...
!> \param io output unit
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_print_xg(x,iter,io)
      IMPLICIT NONE

      INTEGER                                            :: iter,nparticle,i,ia,io
      REAL(KIND=dp), DIMENSION(:)                        :: x

      nparticle=size(x)/3
      write(io,*) nparticle
      write(io,*) "Iteration: ",iter
      do i = 1,nparticle
         ia = (i-1)*3
         write(io,*) x(ia+1)*0.52917721067,x(ia+2)*0.52917721067,x(ia+3)*0.52917721067
      end do

   END SUBROUTINE spgr_print_xg

! **************************************************************************************************
!> \brief routine prints hessian matrix. Only for debugging purpose.
!> \param hess ...
!> \param iter ...
!> \param io output unit
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_print_hess(hess,iter,io)
      IMPLICIT NONE

      INTEGER                                            :: iter,nl,i,j,io
      REAL(KIND=dp), DIMENSION(:,:)                      :: hess

      nl=size(hess(1,:))
      write(io,*) nl/3
      write(io,*) "Iteration: ",iter
      do i = 1,nl/2
         write(io,*) (hess(j,i),j=1,nl/2)
      end do

      do i = nl/2+1,nl
         write(io+1,*) (hess(j,i),j=nl/2+1,nl)
      end do

      do i = 1,nl/2
         write(io+2,*) (hess(j,i),j=nl/2+1,nl)
      end do

      do i = nl/2+1,nl
         write(io+3,*) (hess(j,i),j=1,nl/2)
      end do

   END SUBROUTINE spgr_print_hess

! **************************************************************************************************
!> \brief routine prints single precision hessian matrix. Only for debugging purpose.
!> \param hess ...
!> \param iter ...
!> \param io output unit
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE spgr_print_hess_sp(hess,iter,io)
      IMPLICIT NONE

      INTEGER                                            :: iter,nl,i,j,io
      REAL(KIND=sp), DIMENSION(:,:)                      :: hess

      nl=size(hess(1,:))
      write(io,*) nl/3
      write(io,*) "Iteration: ",iter
      do i = 1,nl/2
         write(io,*) (hess(j,i),j=1,nl/2)
      end do

      do i = nl/2+1,nl
         write(io+1,*) (hess(j,i),j=nl/2+1,nl)
      end do

      do i = 1,nl/2
         write(io+2,*) (hess(j,i),j=nl/2+1,nl)
      end do

      do i = nl/2+1,nl
         write(io+3,*) (hess(j,i),j=1,nl/2)
      end do

   END SUBROUTINE spgr_print_hess_sp
  
END MODULE space_groups
